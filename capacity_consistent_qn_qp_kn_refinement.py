"""Capacity-consistent Qn/Qp-window refinement followed by a constrained kn scan.

The electrode capacities, c_s,max values and stoichiometry widths are changed as
one coupled quantity:

    delta_x = Q_cell / Q_n, delta_y = Q_cell / Q_p

For every capacity candidate only the absolute window positions are re-fitted to
the C/50 qOCV. Dynamic scoring uses the measured branch current and the voltage
at the end of the immediately preceding rest to initialise SOC. The kn scan is
diagnostic only; the current experimental value is retained pending Rct input.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize

import c50_ocp_source_combination_study as c50
import c50_selected_ocp_dynamic_validation as selected
import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import nominal_radius_current_configuration as current
import priority_initial_state_protocol_recheck as priority


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260921_capacity_consistent_qn_qp_kn_refinement"
RESULTS.mkdir(parents=True, exist_ok=True)

NOMINAL_CAPACITY_AH = 2.28

QN_SCALES = (0.98, 1.00, 1.02, 1.04)
QP_SCALES = (1.00, 1.02, 1.04, 49943.0 / 47467.25702278789, 1.06)
KN_VALUES = np.asarray([0.74e-6, 0.85e-6, 0.965e-6, 1.075e-6, 1.20e-6, 1.35e-6])


def fit_capacity_consistent_stage(
    qn: float,
    qp: float,
    qcell: float,
    soc: np.ndarray,
    target: np.ndarray,
    anode: dict,
    cathode: dict,
    initial_stage: dict,
) -> tuple[dict, dict]:
    """Fit absolute positions while Qn/Qp fix the two window widths."""
    dx = qcell / qn
    dy = qcell / qp
    x_min, x_max = float(anode["grid"].min()), float(anode["grid"].max())
    y_min, y_max = float(cathode["grid"].min()), float(cathode["grid"].max())
    bounds = [(x_min, x_max - dx), (y_min, y_max - dy)]
    if bounds[0][1] <= bounds[0][0] or bounds[1][1] <= bounds[1][0]:
        raise ValueError("Capacity candidate produces an OCP window wider than coverage")
    mask = (soc >= 0.02) & (soc <= 0.98)

    def unpack(v):
        x0, y100 = map(float, v)
        return np.asarray([x0, x0 + dx, y100, y100 + dy], dtype=float)

    def objective(v):
        stage = unpack(v)
        voltage = c50.predict(stage, soc, anode, cathode, "equilibrium")
        residual = voltage - target
        endpoint = np.abs(residual[[0, -1]])
        endpoint_excess = np.maximum(endpoint - 0.010, 0.0)
        return float(np.mean(residual[mask] ** 2) + 100.0 * np.mean(endpoint_excess**2))

    x0_seed = np.clip(initial_stage["x0"], *bounds[0])
    y100_seed = np.clip(initial_stage["y100"], *bounds[1])
    starts = [
        np.asarray([x0_seed, y100_seed]),
        np.asarray([bounds[0][0], bounds[1][0]]),
        np.asarray([bounds[0][1], bounds[1][1]]),
        np.asarray([np.mean(bounds[0]), np.mean(bounds[1])]),
    ]
    solutions = [
        minimize(objective, start, method="L-BFGS-B", bounds=bounds, options={"ftol": 1e-15, "maxiter": 2000})
        for start in starts
    ]
    best = min(solutions, key=lambda result: float(result.fun))
    z = unpack(best.x)
    voltage = c50.predict(z, soc, anode, cathode, "equilibrium")
    error_mV = (voltage - target) * 1000.0
    stage = {"x0": z[0], "x100": z[1], "y100": z[2], "y0": z[3]}
    metrics = {
        "Qn_Ah": qn,
        "Qp_Ah": qp,
        "delta_x": dx,
        "delta_y": dy,
        **stage,
        "qOCV_MAE_2_98_mV": float(np.mean(np.abs(error_mV[mask]))),
        "qOCV_RMSE_2_98_mV": float(np.sqrt(np.mean(error_mV[mask] ** 2))),
        "SOC0_error_mV": float(error_mV[0]),
        "SOC100_error_mV": float(error_mV[-1]),
        "max_endpoint_abs_mV": float(np.max(np.abs(error_mV[[0, -1]]))),
        "optimizer_success": bool(best.success),
    }
    return stage, {"voltage": voltage, **metrics}


def model_with_capacities(model: dict, qn: float, qp: float, qn0: float, qp0: float) -> dict:
    changed = dict(model)
    changed["params"] = model["params"].copy()
    changed["csn_max"] = float(model["csn_max"]) * qn / qn0
    changed["csp_max"] = float(model["csp_max"]) * qp / qp0
    changed["qn_cell"] = qn
    changed["qp_cell"] = qp
    changed["params"].update(
        {
            "Maximum concentration in negative electrode [mol.m-3]": changed["csn_max"],
            "Maximum concentration in positive electrode [mol.m-3]": changed["csp_max"],
        },
        check_already_exists=False,
    )
    return changed


def model_with_kn(model: dict, kn: float) -> dict:
    changed = dict(model)
    changed["params"] = model["params"].copy()
    changed["params"].update(
        {
            "Negative electrode exchange-current density [A.m-2]": ga.scale_parameter_function(
                model["params"]["Negative electrode exchange-current density [A.m-2]"],
                float(kn) / current.KN_PREF,
            )
        },
        check_already_exists=False,
    )
    return changed


def run_candidate(
    name: str,
    model: dict,
    stage: dict,
    anode: dict,
    cathode: dict,
    experiment: dict,
    protocol: pd.DataFrame,
    qcell: float,
):
    protocol_index = protocol.set_index(["C_rate", "Direction"])
    rows = []
    runs = {}
    for rate in ga.RATES:
        for charge in (True, False):
            direction = "Charge" if charge else "Discharge"
            rest_v = float(protocol_index.loc[(rate, direction), "Rest_end_V"])
            z_init = priority.monotone_voltage_inverse(
                model["v_qocv"], model["soc"], rest_v
            )
            local_stage = priority.stage_at_soc(stage, z_init, charge)
            measured_current = float(
                np.nanmedian(np.abs(experiment[(rate, charge)]["I_A"]))
            )
            simulation_rate = measured_current / NOMINAL_CAPACITY_AH
            print(f"{name} | {rate:g}C {direction} | z={z_init:.5f}", flush=True)
            simulation = ehq.run_dfn(
                simulation_rate,
                charge,
                local_stage,
                model,
                anode,
                cathode,
                False,
                True,
                positive_initial_h=(-1.0 if charge else 1.0),
            )
            runs[(rate, charge)] = simulation
            metric = priority.curve_metrics(
                experiment[(rate, charge)], simulation, qcell, charge
            )
            rows.append(
                {
                    "Candidate": name,
                    "C_rate": rate,
                    "Direction": direction,
                    "measured_current_A": measured_current,
                    "simulation_rate_C": simulation_rate,
                    "initial_SOC": z_init,
                    **metric,
                    "Capacity_error_mAh": 1000.0
                    * (metric["Q_end_model_Ah"] - metric["Q_end_exp_Ah"]),
                }
            )
    return pd.DataFrame(rows), runs


def summarize(detail: pd.DataFrame) -> pd.DataFrame:
    return detail.groupby(["Candidate", "Direction"], as_index=False).agg(
        mean_center_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
        mean_full_RMSE_mV=("Full_RMSE_mV", "mean"),
        capacity_RMSE_pct=("Capacity_error_pct", lambda x: float(np.sqrt(np.mean(x**2)))),
        mean_capacity_error_mAh=("Capacity_error_mAh", "mean"),
        mean_abs_capacity_error_mAh=("Capacity_error_mAh", lambda x: float(np.mean(np.abs(x)))),
    )


def main() -> None:
    warnings.filterwarnings(
        "ignore", message="The definition of the hysteresis decay rate parameter has changed"
    )
    model0, stage0, anode, cathode, experiment, selected_row, qcell = selected.selected_inputs()
    protocol, _, _ = priority.load_protocol(model0)
    soc = np.asarray(model0["soc"], dtype=float)
    # Use the measured C/50 pseudo-OCV itself. model0["v_qocv"] is the already
    # fitted reconstruction and using it here would leak the old fit into the
    # target and make the baseline qOCV error identically zero.
    measured_c50 = c50.fullcell_c50_target(soc)
    target = np.asarray(measured_c50["qocv"], dtype=float)
    if not np.isclose(float(measured_c50["qcell"]), qcell, rtol=0.0, atol=1e-9):
        raise RuntimeError("C/50 target capacity does not match selected-input capacity")
    qn0, qp0 = float(selected_row.Qn_Ah), float(selected_row.Qp_Ah)

    static_rows = []
    candidates = {}
    for qn_scale in QN_SCALES:
        for qp_scale in QP_SCALES:
            qn, qp = qn0 * qn_scale, qp0 * qp_scale
            stage, fit = fit_capacity_consistent_stage(
                qn, qp, qcell, soc, target, anode, cathode, stage0
            )
            name = f"Qn{qn_scale:.4f}_Qp{qp_scale:.4f}"
            candidate_model = model_with_capacities(model0, qn, qp, qn0, qp0)
            candidate_model["soc"] = soc
            candidate_model["v_qocv"] = fit.pop("voltage")
            fit.update(
                {
                    "Candidate": name,
                    "Qn_scale": qn_scale,
                    "Qp_scale": qp_scale,
                    "csn_max_mol_m3": candidate_model["csn_max"],
                    "csp_max_mol_m3": candidate_model["csp_max"],
                }
            )
            static_rows.append(fit)
            candidates[name] = (candidate_model, stage, fit)

    static = pd.DataFrame(static_rows)
    static["qOCV_admissible"] = (
        (static.qOCV_RMSE_2_98_mV <= 8.0)
        & (static.max_endpoint_abs_mV <= 12.0)
    )
    static.to_csv(RESULTS / "capacity_consistent_static_scan.csv", index=False, encoding="utf-8-sig")

    # Evaluate the full physically coupled grid. This is still small (20 x 6 DFN runs)
    # and avoids choosing Qn/Qp by an arbitrary weighted dynamic objective.
    dynamic_frames = []
    run_cache = {}
    for row in static.itertuples():
        if not row.qOCV_admissible:
            continue
        candidate_model, stage, _ = candidates[row.Candidate]
        detail, runs = run_candidate(
            row.Candidate, candidate_model, stage, anode, cathode, experiment, protocol, qcell
        )
        dynamic_frames.append(detail)
        run_cache[row.Candidate] = runs
    dynamic = pd.concat(dynamic_frames, ignore_index=True)
    dynamic.to_csv(RESULTS / "capacity_consistent_dynamic_detail.csv", index=False, encoding="utf-8-sig")
    dynamic_summary = summarize(dynamic)
    dynamic_summary.to_csv(RESULTS / "capacity_consistent_dynamic_summary.csv", index=False, encoding="utf-8-sig")

    wide = dynamic_summary.pivot(index="Candidate", columns="Direction")
    selection = static.set_index("Candidate").copy()
    for metric in (
        "mean_center_RMSE_mV",
        "capacity_RMSE_pct",
        "mean_abs_capacity_error_mAh",
    ):
        selection[f"charge_{metric}"] = wide[(metric, "Charge")]
        selection[f"discharge_{metric}"] = wide[(metric, "Discharge")]
    selection = selection.reset_index()
    baseline_selection = selection[
        np.isclose(selection.Qn_scale, 1.0) & np.isclose(selection.Qp_scale, 1.0)
    ].iloc[0]
    selection["capacity_change_norm"] = np.sqrt(
        (selection.Qn_scale - 1.0) ** 2 + (selection.Qp_scale - 1.0) ** 2
    )
    selection["dynamic_admissible"] = (
        selection.qOCV_admissible
        & (selection.qOCV_RMSE_2_98_mV <= baseline_selection.qOCV_RMSE_2_98_mV + 0.5)
        & (selection.charge_mean_center_RMSE_mV <= baseline_selection.charge_mean_center_RMSE_mV + 1.0)
        & (selection.discharge_mean_center_RMSE_mV <= baseline_selection.discharge_mean_center_RMSE_mV + 2.0)
        & (selection.charge_capacity_RMSE_pct <= 3.1)
        & (selection.discharge_capacity_RMSE_pct <= 1.0)
    )
    feasible = selection[selection.dynamic_admissible]
    if feasible.empty:
        feasible = selection[selection.qOCV_admissible].copy()
        feasible["fallback_score"] = (
            feasible.charge_capacity_RMSE_pct + feasible.discharge_capacity_RMSE_pct
        )
        capacity_choice = feasible.loc[feasible.fallback_score.idxmin()]
        selection_rule = "minimum charge+discharge capacity RMSE fallback"
    else:
        capacity_choice = feasible.sort_values(
            ["capacity_change_norm", "charge_capacity_RMSE_pct"]
        ).iloc[0]
        selection_rule = (
            "minimum Qn/Qp change subject to: qOCV RMSE <= baseline+0.5 mV, "
            "charge/discharge center RMSE <= baseline+1/+2 mV, and "
            "charge/discharge capacity RMSE <= 3.1/1.0%"
        )
    selection.to_csv(RESULTS / "capacity_candidate_selection.csv", index=False, encoding="utf-8-sig")

    chosen_name = str(capacity_choice.Candidate)
    chosen_model, chosen_stage, chosen_fit = candidates[chosen_name]

    kn_detail_frames = []
    kn_runs = {}
    for kn in KN_VALUES:
        name = f"kn={kn:.4e}"
        local_model = model_with_kn(chosen_model, float(kn))
        detail, runs = run_candidate(
            name, local_model, chosen_stage, anode, cathode, experiment, protocol, qcell
        )
        detail["kn"] = kn
        kn_detail_frames.append(detail)
        kn_runs[float(kn)] = runs
    kn_detail = pd.concat(kn_detail_frames, ignore_index=True)
    kn_detail.to_csv(RESULTS / "kn_refinement_detail.csv", index=False, encoding="utf-8-sig")
    kn_summary = summarize(kn_detail)
    kn_summary["kn"] = kn_summary.Candidate.str.split("=").str[1].astype(float)
    kn_summary.to_csv(RESULTS / "kn_refinement_summary.csv", index=False, encoding="utf-8-sig")

    kn_wide = kn_summary.pivot(index="kn", columns="Direction")
    kn_select = pd.DataFrame(index=kn_wide.index)
    for metric in ("mean_center_RMSE_mV", "capacity_RMSE_pct", "mean_abs_capacity_error_mAh"):
        kn_select[f"charge_{metric}"] = kn_wide[(metric, "Charge")]
        kn_select[f"discharge_{metric}"] = kn_wide[(metric, "Discharge")]
    kn_select = kn_select.reset_index()
    # This scan is diagnostic only.  Do not turn an endpoint-capacity improvement
    # into a new kinetic property; retain the current experimental value until the
    # SOC-resolved Rct analysis supplies an independently identified k_n.
    kn_choice = kn_select.iloc[np.argmin(np.abs(kn_select.kn - current.KN_PREF))]
    kn_rule = "retain current experimental kn; scan is diagnostic pending SOC-resolved Rct"
    kn_select["selected"] = np.isclose(kn_select.kn, float(kn_choice.kn))
    kn_select.to_csv(RESULTS / "kn_candidate_selection.csv", index=False, encoding="utf-8-sig")

    # Trade-off figure
    merged = selection[selection.qOCV_admissible].copy()
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8), constrained_layout=True)
    sc = axes[0].scatter(
        merged.Qp_scale,
        merged.qOCV_RMSE_2_98_mV,
        c=merged.Qn_scale,
        cmap="viridis",
        s=65,
    )
    axes[0].set(xlabel="Qp scale", ylabel="qOCV RMSE [mV]", title="Capacity-consistent qOCV")
    fig.colorbar(sc, ax=axes[0], label="Qn scale")
    axes[1].scatter(
        merged.charge_capacity_RMSE_pct,
        merged.charge_mean_center_RMSE_mV,
        c=merged.Qp_scale,
        cmap="plasma",
        s=65,
    )
    axes[1].set(xlabel="Charge capacity RMSE [%]", ylabel="Charge voltage RMSE [mV]", title="Qn/Qp trade-off")
    axes[2].plot(kn_select.kn * 1e7, kn_select.charge_mean_center_RMSE_mV, "o-", label="charge voltage")
    ax2 = axes[2].twinx()
    ax2.plot(kn_select.kn * 1e7, kn_select.charge_capacity_RMSE_pct, "s--", color="#D55E00", label="charge capacity")
    axes[2].axvline(float(kn_choice.kn) * 1e7, color="black", ls=":")
    axes[2].set(xlabel="kn [1e-7]", ylabel="Voltage RMSE [mV]", title="Constrained kn selection")
    ax2.set_ylabel("Capacity RMSE [%]", color="#D55E00")
    for ax in axes:
        ax.grid(alpha=0.25)
    fig.savefig(RESULTS / "capacity_and_kn_tradeoff.png", dpi=220)
    plt.close(fig)

    # Capacity-feasibility curves at the retained experimental kn
    final_runs = kn_runs[float(kn_choice.kn)]
    fig, axes = plt.subplots(2, 3, figsize=(16, 8.5), constrained_layout=True)
    for col, rate in enumerate(ga.RATES):
        for row_index, charge in enumerate((True, False)):
            ax = axes[row_index, col]
            exp = experiment[(rate, charge)]
            sim = final_runs[(rate, charge)]
            q_exp = priority.transferred_capacity(exp, qcell, charge)
            q_sim = priority.transferred_capacity(sim, qcell, charge)
            ax.plot(q_exp, exp["V"], "k", lw=2.3, label="Experiment")
            ax.plot(q_sim, sim["V"], color="#0072B2", lw=2.0, label="Feasibility model")
            ax.set_title(f"{rate:g}C {'Charge' if charge else 'Discharge'}")
            ax.set_xlabel("Transferred capacity [Ah]")
            ax.set_ylabel("Voltage [V]")
            ax.grid(alpha=0.22)
    axes[0, 0].legend()
    fig.savefig(RESULTS / "refined_charge_discharge_curves.png", dpi=220)
    plt.close(fig)

    report = {
        "method": {
            "capacity_relation": "delta_x=Qcell/Qn, delta_y=Qcell/Qp; c_s_max scales with Qn/Qp at fixed geometry",
            "initial_SOC": "candidate qOCV inverse of the immediately preceding rest-end voltage",
            "current": "measured branch median absolute CC current",
            "hysteresis": "negative off; positive on with preceding-branch history state",
            "capacity_selection_rule": selection_rule,
            "kn_selection_rule": kn_rule,
        },
        "baseline": {
            "Qn_Ah": qn0,
            "Qp_Ah": qp0,
            "csn_max_mol_m3": float(model0["csn_max"]),
            "csp_max_mol_m3": float(model0["csp_max"]),
            **stage0,
        },
        "adoption_decision": "retain baseline; Qn/Qp and kn scans are feasibility diagnostics only",
        "best_capacity_feasibility_candidate": capacity_choice.to_dict(),
        "feasibility_stage": chosen_stage,
        "feasibility_csmax": {
            "negative_mol_m3": float(chosen_model["csn_max"]),
            "positive_mol_m3": float(chosen_model["csp_max"]),
        },
        "retained_kn": kn_choice.to_dict(),
        "qOCV_metrics": {key: value for key, value in chosen_fit.items() if key != "voltage"},
    }
    (RESULTS / "capacity_consistent_refinement_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("\nBEST CAPACITY FEASIBILITY CANDIDATE (NOT ADOPTED)")
    print(capacity_choice.to_string())
    print("\nRETAINED EXPERIMENTAL KN")
    print(kn_choice.to_string())
    print("Saved to", RESULTS)


if __name__ == "__main__":
    main()
