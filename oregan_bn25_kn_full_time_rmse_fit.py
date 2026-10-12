"""Fit only k_n by minimizing full time-voltage RMSE over six CC branches.

No HPPC, capacity, endpoint-time, or anode-potential term is included in the
objective.  Each 0.5C/1C/2C charge/discharge branch contributes 400 points on
its common measured/model time interval, so the six branches have equal weight.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

import comprehensive_prefit_cross_cohort as cross
import gitt_ocp_analysis as ga
import priority_initial_state_protocol_recheck as priority

import analyze_oregan_brugg_2p5_anode_potential as anode_analysis
import electrolyte_transport_set_comparison as comparison
import oregan_bn25_kn_only_fit as prior_fit


OUT = Path(__file__).resolve().parent / "results" / "261012_kn_full_time_rmse_fit"
LOWER_SCALE = 0.35
UPPER_SCALE = 2.50


def branch_time_error_mV(experiment: dict, simulation: dict) -> np.ndarray:
    te = np.asarray(experiment["t_min"], dtype=float)
    tm = np.asarray(simulation["time_s"], dtype=float) / 60.0
    te = te - te[0]
    tm = tm - tm[0]
    ve = np.asarray(experiment["V"], dtype=float)
    vm = np.asarray(simulation["V"], dtype=float)
    upper = 0.995 * min(float(te[-1]), float(tm[-1]))
    grid = np.linspace(0.0, upper, 400)
    return 1000.0 * (np.interp(grid, tm, vm) - np.interp(grid, te, ve))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    baseline, stage, anode, cathode, experiment, qcell = comparison.current_p0()
    transport = comparison.with_transport_set(baseline, "ORegan2022")
    cache: dict[float, dict] = {}
    profile_rows: list[dict] = []

    def evaluate(scale: float) -> dict:
        key = round(float(scale), 8)
        if key in cache:
            return cache[key]
        model = prior_fit.apply_kn(transport, float(scale))
        runs = {}
        errors = []
        for rate in ga.RATES:
            for charge in (True, False):
                rest_v, measured_current = comparison.PROTOCOL[(rate, charge)]
                z0 = priority.monotone_voltage_inverse(
                    model["v_qocv"], model["soc"], rest_v
                )
                local_stage = priority.stage_at_soc(stage, z0, charge)
                output = anode_analysis.run_branch(
                    measured_current / cross.NOMINAL_CAPACITY_AH,
                    charge,
                    local_stage,
                    model,
                    anode,
                    cathode,
                )
                runs[(rate, charge)] = output
                errors.append(branch_time_error_mV(experiment[(rate, charge)], output))
        all_errors = np.concatenate(errors)
        result = {
            "scale": float(scale),
            "rmse_mV": float(np.sqrt(np.mean(all_errors**2))),
            "mae_mV": float(np.mean(np.abs(all_errors))),
            "runs": runs,
        }
        cache[key] = result
        profile_rows.append(
            {
                "kn_scale": float(scale),
                "kn_table_value": 1.0e-11 * float(scale),
                "kn_PyBaMM_prefactor": prior_fit.BASE["kn"] * float(scale),
                "Full_time_voltage_RMSE_mV": result["rmse_mV"],
                "Full_time_voltage_MAE_mV": result["mae_mV"],
            }
        )
        print(
            "scale={:.6f}, full-time RMSE={:.4f} mV".format(
                scale, result["rmse_mV"]
            ),
            flush=True,
        )
        return result

    baseline_result = evaluate(1.0)

    optimization = minimize_scalar(
        lambda log_scale: evaluate(float(np.exp(log_scale)))["rmse_mV"] ** 2,
        bounds=(np.log(LOWER_SCALE), np.log(UPPER_SCALE)),
        method="bounded",
        options={"xatol": 1.0e-3, "maxiter": 24},
    )
    best_scale = float(np.exp(optimization.x))
    best = evaluate(best_scale)

    metric_rows = []
    anode_rows = []
    for label, result in (("Retained anchor", baseline_result), ("RMSE fit", best)):
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                output = result["runs"][(rate, charge)]
                exp = experiment[(rate, charge)]
                error = branch_time_error_mV(exp, output)
                capacity = priority.curve_metrics(exp, output, qcell, charge)
                exp_time = np.asarray(exp["t_min"], dtype=float)
                exp_time -= exp_time[0]
                sim_time = np.asarray(output["time_s"], dtype=float) / 60.0
                sim_time -= sim_time[0]
                metric_rows.append(
                    {
                        "Scenario": label,
                        "kn_scale": result["scale"],
                        "C_rate": rate,
                        "Direction": direction,
                        "Full_time_voltage_RMSE_mV": float(
                            np.sqrt(np.mean(error**2))
                        ),
                        "Full_time_voltage_MAE_mV": float(np.mean(np.abs(error))),
                        "Capacity_domain_full_RMSE_mV": capacity["Full_RMSE_mV"],
                        "Center10_70_RMSE_mV": capacity["Center10_70_RMSE_mV"],
                        "Experiment_end_time_min": float(exp_time[-1]),
                        "DFN_end_time_min": float(sim_time[-1]),
                        "End_time_error_min": float(sim_time[-1] - exp_time[-1]),
                        "Capacity_error_pct": capacity["Capacity_error_pct"],
                        "Capacity_error_mAh": 1000.0
                        * (capacity["Q_end_model_Ah"] - capacity["Q_end_exp_Ah"]),
                    }
                )
                potential = np.asarray(output["potential_V"], dtype=float)
                minimum_index = int(np.nanargmin(potential))
                anode_rows.append(
                    {
                        "Scenario": label,
                        "kn_scale": result["scale"],
                        "C_rate": rate,
                        "Direction": direction,
                        "Minimum_anode_potential_mV": 1000.0
                        * float(potential[minimum_index]),
                        "Time_below_0V_pct": 100.0
                        * anode_analysis.time_fraction(
                            output["time_s"], potential < 0.0
                        ),
                        "Surface_stoichiometry_at_minimum": float(
                            output["surface_sto"][minimum_index]
                        ),
                    }
                )

    metrics = pd.DataFrame(metric_rows)
    anode_summary = pd.DataFrame(anode_rows)
    fitted_metrics = metrics[metrics.Scenario == "RMSE fit"]
    capacity_rmse = float(
        np.sqrt(np.mean(fitted_metrics.Capacity_error_pct.to_numpy(float) ** 2))
    )
    max_capacity_error = float(
        np.max(np.abs(fitted_metrics.Capacity_error_mAh.to_numpy(float)))
    )

    profile = pd.DataFrame(profile_rows).drop_duplicates("kn_scale").sort_values("kn_scale")
    profile.to_csv(OUT / "kn_full_time_rmse_profile.csv", index=False, encoding="utf-8-sig")
    metrics.to_csv(OUT / "branch_metrics.csv", index=False, encoding="utf-8-sig")
    anode_summary.to_csv(OUT / "anode_potential_summary.csv", index=False, encoding="utf-8-sig")

    summary = {
        "objective": "full elapsed-time voltage RMSE; six CC branches equal weight",
        "hppc_in_objective": False,
        "capacity_or_endpoint_in_objective": False,
        "best_kn_scale": best_scale,
        "best_kn_table_value": 1.0e-11 * best_scale,
        "best_kn_PyBaMM_prefactor": prior_fit.BASE["kn"] * best_scale,
        "baseline_objective_RMSE_mV": baseline_result["rmse_mV"],
        "best_objective_RMSE_mV": best["rmse_mV"],
        "posthoc_capacity_RMSE_pct": capacity_rmse,
        "posthoc_max_abs_capacity_error_mAh": max_capacity_error,
        "optimizer_success": bool(optimization.success),
    }
    (OUT / "fit_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    fig, axes = plt.subplots(2, 3, figsize=(15.8, 8.6), sharey=True, constrained_layout=True)
    for column, rate in enumerate(ga.RATES):
        for row, charge in enumerate((True, False)):
            direction = "Charge" if charge else "Discharge"
            ax = axes[row, column]
            exp = experiment[(rate, charge)]
            exp_time = np.asarray(exp["t_min"], dtype=float)
            exp_time -= exp_time[0]
            ax.plot(exp_time, exp["V"], color="#202020", lw=2.4, label="Experiment")
            for label, result, color, style in (
                ("Retained anchor", baseline_result, "#E69F00", "--"),
                ("Full-time RMSE fit", best, "#0072B2", "-"),
            ):
                output = result["runs"][(rate, charge)]
                time_min = np.asarray(output["time_s"], dtype=float) / 60.0
                time_min -= time_min[0]
                ax.plot(time_min, output["V"], color=color, ls=style, lw=2.0, label=label)
            row_metric = fitted_metrics[
                (fitted_metrics.C_rate == rate)
                & (fitted_metrics.Direction == direction)
            ].iloc[0]
            ax.text(
                0.03,
                0.04,
                (
                    f"fit RMSE {row_metric.Full_time_voltage_RMSE_mV:.1f} mV\n"
                    f"capacity error {row_metric.Capacity_error_pct:+.2f}%"
                ),
                transform=ax.transAxes,
                fontsize=9.2,
                va="bottom",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.86),
            )
            ax.set_title(f"{rate:g}C {direction}", fontsize=12, weight="bold")
            ax.set_ylim(2.9, 4.3)
            ax.grid(alpha=0.24)
            if row == 1:
                ax.set_xlabel("Elapsed time (min)")
            if column == 0:
                ax.set_ylabel("Cell voltage (V)")

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.suptitle(
        "Cell 6-6: $k_n$ fitted only to full time-voltage RMSE\n"
        "ORegan2022 electrolyte transport + $b_n=2.5$; no HPPC/endpoint term",
        fontsize=15,
        weight="bold",
        y=1.08,
    )
    fig.savefig(OUT / "kn_full_time_rmse_fit_voltage_curves.png", dpi=240, bbox_inches="tight")
    plt.close(fig)

    print("\n[summary]\n", json.dumps(summary, ensure_ascii=False, indent=2))
    print("\n[branch metrics]\n", metrics.to_string(index=False))
    print("\n[anode]\n", anode_summary.to_string(index=False))


if __name__ == "__main__":
    main()
