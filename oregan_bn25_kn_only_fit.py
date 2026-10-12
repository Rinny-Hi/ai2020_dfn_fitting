"""Fit only negative-electrode kinetics with ORegan2022 and b_n=2.5.

The positive-electrode kinetic prefactor and both solid diffusivities remain at
their EIS/GITT anchors.  Three one-dimensional profile fits are reported:

* dynamic voltage (0.5C, 1C, 2C charge/discharge, central overlap),
* HPPC pulse voltage changes relative to the adjacent rest voltage, and
* an equal-dataset-weight combination of those two objectives.

Absolute HPPC key-point voltage is retained as a validation metric and plot,
but is not used as the kinetic objective because it includes the OCV/initial-
state offset that k_n cannot uniquely correct.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import comprehensive_prefit_cross_cohort as cross
import gitt_ocp_analysis as ga
import kinetics_original_reversion_ablation as kinetics
import prefit_physical_anchor_decision as p0
import priority_initial_state_protocol_recheck as priority

import analyze_oregan_brugg_2p5_anode_potential as anode_analysis
import electrolyte_transport_set_comparison as comparison
import oregan_combined_sensitivity_identifiability as combined


OUT = comparison.ROOT / "results" / "261012_oregan_bn25_kn_only_fit"
BRUGG_N = 2.5
FARADAY = 96485.33212
BASE = dict(combined.INITIAL)
PROFILE_SCALES = np.geomspace(0.35, 2.5, 13)


def apply_kn(base: dict, scale: float) -> dict:
    changed = kinetics.with_prefactors(base, BASE["kn"] * float(scale), BASE["kp"])
    changed = dict(changed)
    changed["params"] = changed["params"].copy()
    changed["params"].update(
        {
            "Negative particle diffusivity [m2.s-1]": BASE["Dsn"],
            "Positive particle diffusivity [m2.s-1]": BASE["Dsp"],
            "Negative electrode Bruggeman coefficient (electrolyte)": BRUGG_N,
        },
        check_already_exists=False,
    )
    return changed


def hppc_pulse_vector(values: np.ndarray) -> np.ndarray:
    """Return discharge drop and charge rise for each of ten SOC levels."""
    matrix = np.asarray(values, dtype=float).reshape(10, len(combined.HPPC_POINT_ORDER))
    discharge_drop = matrix[:, 1] - matrix[:, 0]
    charge_rise = matrix[:, 4] - matrix[:, 3]
    return np.column_stack([discharge_drop, charge_rise]).ravel()


def quadratic_profile_minimum(scales: np.ndarray, values: np.ndarray) -> float:
    """Refine a profile-grid minimum by a local quadratic in ln(scale)."""
    values = np.asarray(values, dtype=float)
    index = int(np.nanargmin(values))
    if index == 0 or index == len(scales) - 1:
        return float(scales[index])
    x = np.log(scales[index - 1 : index + 2])
    y = values[index - 1 : index + 2]
    coefficient = np.polyfit(x, y, 2)
    if coefficient[0] <= 0:
        return float(scales[index])
    x_min = -coefficient[1] / (2.0 * coefficient[0])
    x_min = float(np.clip(x_min, x[0], x[-1]))
    return float(np.exp(x_min))


def interpolate_curve(source: dict, qcell: float, charge: bool, grid: np.ndarray) -> np.ndarray:
    q = priority.transferred_capacity(source, qcell, charge)
    order = np.argsort(q, kind="stable")
    q_unique, index = np.unique(np.asarray(q)[order], return_index=True)
    voltage = np.asarray(source["V"], dtype=float)[order][index]
    return np.interp(grid, q_unique, voltage)


def evaluate_dynamic_metrics(name, runs, experiment, qcell):
    rows = []
    for rate in ga.RATES:
        for charge in (True, False):
            metrics = priority.curve_metrics(
                experiment[(rate, charge)], runs[(rate, charge)], qcell, charge
            )
            rows.append(
                {
                    "Scenario": name,
                    "Trust": "one-parameter profile fit",
                    "Note": "ORegan2022; b_n=2.5; only k_n changed",
                    "Cell": "6-6",
                    "C_rate": rate,
                    "Direction": "Charge" if charge else "Discharge",
                    **metrics,
                    "Capacity_error_mAh": 1000.0
                    * (metrics["Q_end_model_Ah"] - metrics["Q_end_exp_Ah"]),
                }
            )
    return rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    baseline, stage, anode, cathode, experiment, qcell = comparison.current_p0()
    transport = comparison.with_transport_set(baseline, "ORegan2022")

    hppc_targets = combined.load_hppc_targets()
    hppc_times = combined.hppc_target_times()
    hppc_exp = hppc_targets["Voltage_V"].to_numpy(float)
    hppc_pulse_exp = hppc_pulse_vector(hppc_exp)

    baseline_model = apply_kn(transport, 1.0)
    hppc_z0 = priority.monotone_voltage_inverse(
        baseline_model["v_qocv"], baseline_model["soc"], float(hppc_exp[0])
    )
    hppc_stage = priority.stage_at_soc(stage, hppc_z0, False)
    target_times = hppc_times["Time_s"].to_numpy(float)

    cache = {}

    def evaluate(scale: float):
        key = round(float(scale), 9)
        if key in cache:
            return cache[key]
        print("k_n scale {:.6f}".format(scale), flush=True)
        model = apply_kn(transport, scale)
        runs = combined.run_dynamic(model, stage, anode, cathode)
        hppc_model = combined.run_hppc(
            model, hppc_stage, anode, cathode, target_times
        )
        cache[key] = {"model": model, "runs": runs, "hppc": hppc_model}
        return cache[key]

    base_eval = evaluate(1.0)
    grids = combined.make_dynamic_grids(experiment, base_eval["runs"], qcell)
    dynamic_exp = combined.dynamic_vector(experiment, grids, qcell, experiment=True)

    def objective_values(result):
        dynamic_model = combined.dynamic_vector(result["runs"], grids, qcell)
        dynamic_error = 1000.0 * (dynamic_model - dynamic_exp)
        hppc_absolute_error = 1000.0 * (result["hppc"] - hppc_exp)
        hppc_pulse_error = 1000.0 * (
            hppc_pulse_vector(result["hppc"]) - hppc_pulse_exp
        )
        dynamic_mse = float(np.mean(dynamic_error**2))
        pulse_mse = float(np.mean(hppc_pulse_error**2))
        return {
            "Dynamic_RMSE_mV": float(np.sqrt(dynamic_mse)),
            "HPPC_pulse_RMSE_mV": float(np.sqrt(pulse_mse)),
            "HPPC_absolute_RMSE_mV": float(
                np.sqrt(np.mean(hppc_absolute_error**2))
            ),
            "HPPC_absolute_MAE_mV": float(np.mean(np.abs(hppc_absolute_error))),
            "HPPC_absolute_bias_mV": float(np.mean(hppc_absolute_error)),
            "Combined_equal_weight_RMSE_mV": float(
                np.sqrt(0.5 * (dynamic_mse + pulse_mse))
            ),
        }

    profile_rows = []
    for scale in PROFILE_SCALES:
        metrics = objective_values(evaluate(float(scale)))
        profile_rows.append({"kn_scale": float(scale), **metrics})
    profile = pd.DataFrame(profile_rows).sort_values("kn_scale")

    scale_dynamic = quadratic_profile_minimum(
        profile.kn_scale.to_numpy(), profile.Dynamic_RMSE_mV.to_numpy() ** 2
    )
    scale_hppc = quadratic_profile_minimum(
        profile.kn_scale.to_numpy(), profile.HPPC_pulse_RMSE_mV.to_numpy() ** 2
    )
    scale_combined = quadratic_profile_minimum(
        profile.kn_scale.to_numpy(),
        profile.Combined_equal_weight_RMSE_mV.to_numpy() ** 2,
    )
    candidates = {
        "Pre-fit": 1.0,
        "Dynamic-only k_n fit": scale_dynamic,
        "HPPC-pulse-only k_n fit": scale_hppc,
        "Combined k_n fit": scale_combined,
    }

    candidate_rows = []
    dynamic_detail_rows = []
    hppc_rows = []
    candidate_results = {}
    for name, scale in candidates.items():
        result = evaluate(scale)
        candidate_results[name] = result
        metrics = objective_values(result)
        candidate_rows.append(
            {
                "Scenario": name,
                "kn_scale": scale,
                "kn_PyBaMM_prefactor": BASE["kn"] * scale,
                "kn_table_m2p5_mol_m0p5_s": BASE["kn"] * scale / FARADAY,
                "kp_fixed_PyBaMM_prefactor": BASE["kp"],
                "Dsn_fixed_m2_s": BASE["Dsn"],
                "Dsp_fixed_m2_s": BASE["Dsp"],
                "brugg_n_fixed": BRUGG_N,
                **metrics,
            }
        )
        dynamic_detail_rows.extend(
            evaluate_dynamic_metrics(name, result["runs"], experiment, qcell)
        )
        for row_index, target in hppc_targets.iterrows():
            model_voltage = float(result["hppc"][row_index])
            hppc_rows.append(
                {
                    "Scenario": name,
                    "Cycle": int(target.Cycle),
                    "Nominal_SOC": float(target.Nominal_SOC),
                    "Point": str(target.Point),
                    "Experimental_voltage_V": float(target.Voltage_V),
                    "Model_voltage_V": model_voltage,
                    "Error_mV": 1000.0 * (model_voltage - float(target.Voltage_V)),
                }
            )

    candidates_df = pd.DataFrame(candidate_rows)
    detail = pd.DataFrame(dynamic_detail_rows)
    summary = p0.summarize(detail)
    branch = p0.branch_summary(detail, list(candidates))
    hppc_detail = pd.DataFrame(hppc_rows)

    # Negative-electrode potential is a prediction check, not part of the fit.
    anode_rows = []
    for name in ("Pre-fit", "Combined k_n fit"):
        result = candidate_results[name]
        for rate in ga.RATES:
            for charge in (True, False):
                output = result["runs"][(rate, charge)]
                potential = np.asarray(output["potential_V"], dtype=float)
                minimum_index = int(np.nanargmin(potential))
                anode_rows.append(
                    {
                        "Scenario": name,
                        "C_rate": rate,
                        "Direction": "Charge" if charge else "Discharge",
                        "Minimum_anode_potential_mV": 1000.0
                        * float(potential[minimum_index]),
                        "Below_0V_time_fraction_pct": 100.0
                        * anode_analysis.time_fraction(
                            output["time_s"], potential < 0.0
                        ),
                        "Surface_stoichiometry_at_minimum": float(
                            output["surface_sto"][minimum_index]
                        ),
                    }
                )
    anode_summary = pd.DataFrame(anode_rows)

    profile.to_csv(OUT / "kn_profile.csv", index=False, encoding="utf-8-sig")
    candidates_df.to_csv(
        OUT / "candidate_objectives.csv", index=False, encoding="utf-8-sig"
    )
    detail.to_csv(OUT / "dynamic_detail.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(OUT / "dynamic_summary.csv", index=False, encoding="utf-8-sig")
    branch.to_csv(OUT / "dynamic_branch_summary.csv", index=False, encoding="utf-8-sig")
    hppc_detail.to_csv(OUT / "hppc_keypoints.csv", index=False, encoding="utf-8-sig")
    anode_summary.to_csv(OUT / "anode_summary.csv", index=False, encoding="utf-8-sig")

    # Dynamic voltage comparison: pre-fit versus the balanced one-parameter fit.
    fig, axes = plt.subplots(2, 3, figsize=(15.4, 8.4), constrained_layout=True)
    for column, rate in enumerate(ga.RATES):
        for row, charge in enumerate((True, False)):
            ax = axes[row, column]
            exp = experiment[(rate, charge)]
            q_exp = priority.transferred_capacity(exp, qcell, charge)
            ax.plot(q_exp, exp["V"], color="#111111", lw=2.1, label="Experiment")
            for name, color, style in (
                ("Pre-fit", "#E69F00", "--"),
                ("Combined k_n fit", "#0072B2", "-"),
            ):
                run = candidate_results[name]["runs"][(rate, charge)]
                q_model = priority.transferred_capacity(run, qcell, charge)
                ax.plot(q_model, run["V"], color=color, ls=style, lw=1.9, label=name)
            fit_row = branch[
                (branch.Scenario == "Combined k_n fit")
                & (branch.C_rate == rate)
                & (branch.Direction == ("Charge" if charge else "Discharge"))
            ].iloc[0]
            ax.text(
                0.04,
                0.05,
                "fit RMSE: {:.1f} mV\ncapacity error: {:+.2f}%".format(
                    fit_row.Full_RMSE_mV, fit_row.Capacity_error_pct
                ),
                transform=ax.transAxes,
                fontsize=9,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.84),
            )
            ax.set_title("{}C {}".format(rate, "charge" if charge else "discharge"))
            ax.grid(alpha=0.24)
            if row == 1:
                ax.set_xlabel("Transferred capacity (Ah)")
            if column == 0:
                ax.set_ylabel("Terminal voltage (V)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 1.025),
    )
    fig.suptitle(
        "Dynamic voltage: ORegan2022, $b_n=2.5$, only $k_n$ fitted",
        fontsize=15,
        weight="bold",
        y=1.075,
    )
    fig.savefig(OUT / "dynamic_voltage_kn_only_fit.png", dpi=240, bbox_inches="tight")
    plt.close(fig)

    # Absolute HPPC key points remain useful for visual validation of bias/shape.
    point_titles = {
        "rest_pre_end": "Rest before pulse",
        "dchg_end": "End of discharge pulse",
        "rest_30s": "30 s recovery",
        "rest_end": "End of rest",
        "chg_end": "End of charge pulse",
    }
    fig, axes = plt.subplots(2, 3, figsize=(15.4, 8.4), constrained_layout=True)
    axes_flat = axes.ravel()
    for index, point in enumerate(combined.HPPC_POINT_ORDER):
        ax = axes_flat[index]
        target = hppc_targets[hppc_targets.Point.astype(str) == point]
        ax.plot(
            100.0 * target.Nominal_SOC,
            target.Voltage_V,
            "o-",
            color="#111111",
            lw=2.0,
            ms=4.5,
            label="Experiment",
        )
        for name, color, style in (
            ("Pre-fit", "#E69F00", "--"),
            ("Combined k_n fit", "#0072B2", "-"),
        ):
            subset = hppc_detail[
                (hppc_detail.Scenario == name) & (hppc_detail.Point == point)
            ]
            ax.plot(
                100.0 * subset.Nominal_SOC,
                subset.Model_voltage_V,
                color=color,
                ls=style,
                lw=1.9,
                marker="s",
                ms=3.5,
                label=name,
            )
        ax.set_title(point_titles[point])
        ax.set_xlabel("Nominal SOC (%)")
        ax.set_ylabel("Voltage (V)")
        ax.grid(alpha=0.24)
        ax.invert_xaxis()

    ax = axes_flat[5]
    x = np.arange(len(candidates_df))
    width = 0.36
    ax.bar(
        x - width / 2,
        candidates_df.HPPC_pulse_RMSE_mV,
        width,
        color="#0072B2",
        label="Pulse-change RMSE",
    )
    ax.bar(
        x + width / 2,
        candidates_df.HPPC_absolute_RMSE_mV,
        width,
        color="#E69F00",
        label="Absolute key-point RMSE",
    )
    ax.set_xticks(x)
    ax.set_xticklabels(["Pre", "Dynamic", "HPPC", "Combined"])
    ax.set_ylabel("RMSE (mV)")
    ax.set_title("HPPC error comparison")
    ax.grid(axis="y", alpha=0.24)
    ax.legend(frameon=False, fontsize=8)
    axes_flat[0].legend(loc="lower left", frameon=False, fontsize=8)
    fig.suptitle(
        "HPPC validation: measured versus DFN key-point voltage",
        fontsize=15,
        weight="bold",
        y=1.035,
    )
    fig.savefig(OUT / "hppc_kn_only_fit.png", dpi=240, bbox_inches="tight")
    plt.close(fig)

    # Profile plot makes the dataset preference and one-parameter compromise explicit.
    fig, ax = plt.subplots(figsize=(9.6, 6.2), constrained_layout=True)
    ax.semilogx(
        profile.kn_scale, profile.Dynamic_RMSE_mV, "o-", label="Dynamic voltage"
    )
    ax.semilogx(
        profile.kn_scale, profile.HPPC_pulse_RMSE_mV, "s-", label="HPPC pulse change"
    )
    ax.semilogx(
        profile.kn_scale,
        profile.Combined_equal_weight_RMSE_mV,
        "^-",
        label="Equal-weight combined",
    )
    for name, scale in candidates.items():
        if name != "Pre-fit":
            ax.axvline(scale, lw=1.0, ls=":", alpha=0.55)
    ax.set_xlabel("$k_n$ scale relative to EIS anchor")
    ax.set_ylabel("RMSE (mV)")
    ax.set_title("One-dimensional $k_n$ profile")
    ax.grid(alpha=0.24)
    ax.legend(frameon=False)
    fig.savefig(OUT / "kn_objective_profile.png", dpi=240, bbox_inches="tight")
    plt.close(fig)

    metadata = {
        "electrolyte_transport": "ORegan2022",
        "negative_electrode_bruggeman": BRUGG_N,
        "fit_parameter": "negative-electrode kinetic prefactor only",
        "fixed_parameters": {key: BASE[key] for key in ("Dsn", "Dsp", "kp")},
        "hppc_initial_soc": hppc_z0,
        "hppc_fit_metric": "20 pulse voltage changes relative to adjacent rest",
        "dynamic_fit_metric": "central 10-70% capacity-overlap voltage",
        "combined_weighting": "equal dataset weight (mean of two MSE values)",
        "profile_scale_bounds": [float(PROFILE_SCALES[0]), float(PROFILE_SCALES[-1])],
    }
    (OUT / "run_metadata.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print("\n[candidates]\n", candidates_df.to_string(index=False))
    print("\n[dynamic summary]\n", summary.to_string(index=False))
    print("\n[anode]\n", anode_summary.to_string(index=False))
    print("\nSaved:", OUT)


if __name__ == "__main__":
    main()
