"""Re-identify the OCP window after the measured-loading update, then compare L."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution, minimize

import c50_ocp_source_combination_study as c50
import c50_selected_ocp_dynamic_validation as selected
import comprehensive_prefit_cross_cohort as cross
import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import loading_anchored_thickness_p0_comparison as raw_study
import priority_initial_state_protocol_recheck as priority
import remeasured_electrode_mass_sensitivity as mass_study


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260927_loading_thickness_window_reidentified"
RESULTS.mkdir(parents=True, exist_ok=True)


def fixed_capacity_window(qcell, qn, qp, target, soc, anode, cathode):
    dx, dy = qcell / qn, qcell / qp
    xmin, xmax = float(np.min(anode["grid"])), float(np.max(anode["grid"]))
    ymin, ymax = float(np.min(cathode["grid"])), float(np.max(cathode["grid"]))
    if dx > xmax - xmin or dy > ymax - ymin:
        raise ValueError("Measured-capacity window is outside OCP support")
    mask = (soc >= 0.02) & (soc <= 0.98)

    def stage(v):
        x0, y100 = map(float, v)
        return np.array([x0, x0 + dx, y100, y100 + dy])

    def mse(v):
        err = c50.predict(stage(v), soc, anode, cathode) - target
        return float(np.mean(err[mask] ** 2))

    def endpoint(v):
        pred = c50.predict(stage(v), soc[[0, -1]], anode, cathode)
        return pred - target[[0, -1]]

    def violation(v):
        e = endpoint(v)
        return float(sum(max(0.0, abs(float(x)) - 0.010) ** 2 for x in e))

    def penalized(v):
        return mse(v) + 1e4 * violation(v)

    bounds = [(xmin, xmax - dx), (ymin, ymax - dy)]
    de = differential_evolution(
        penalized, bounds, seed=260927, popsize=20, maxiter=500, tol=1e-11, polish=False
    )
    constraints = [
        {"type": "ineq", "fun": lambda v, i=i: 0.010 - abs(float(endpoint(v)[i]))}
        for i in range(2)
    ]
    local = minimize(
        mse,
        de.x,
        method="SLSQP",
        bounds=bounds,
        constraints=constraints,
        options={"maxiter": 3000, "ftol": 1e-14},
    )
    best = local.x if local.success and violation(local.x) <= violation(de.x) + 1e-12 else de.x
    z = stage(best)
    pred = c50.predict(z, soc, anode, cathode)
    err = 1000.0 * (pred - target)
    return z, pred, {
        "Qn_Ah": qn,
        "Qp_Ah": qp,
        "delta_x": dx,
        "delta_y": dy,
        "qOCV_MAE_mV": float(np.mean(np.abs(err[mask]))),
        "qOCV_RMSE_mV": float(np.sqrt(np.mean(err[mask] ** 2))),
        "endpoint_0_mV": float(err[0]),
        "endpoint_100_mV": float(err[-1]),
        "optimizer_success": bool(local.success),
    }


def main():
    cell_curves, protocol = raw_study.load_july_cells()
    model0, old_stage, anode, cathode, _, selected_row, qcell = selected.selected_inputs()
    target_bundle = c50.fullcell_c50_target(np.asarray(model0["soc"], float))
    target = np.asarray(target_bundle["qocv"], float)

    qn = float(selected_row.Qn_Ah) * mass_study.MASS_RATIO["Negative"]
    qp = float(selected_row.Qp_Ah) * mass_study.MASS_RATIO["Positive"]
    z, qocv, window_metrics = fixed_capacity_window(
        qcell, qn, qp, target, np.asarray(model0["soc"], float), anode, cathode
    )
    stage = {"x0": z[0], "x100": z[1], "y100": z[2], "y0": z[3]}

    models = {
        name: raw_study.build_case(model0, *dims)
        for name, dims in raw_study.CASES.items()
    }
    for model in models.values():
        model["v_qocv"] = qocv.copy()
        model["delta_x"] = window_metrics["delta_x"]
        model["delta_y"] = window_metrics["delta_y"]
        model["qn_cell"] = qn
        model["qp_cell"] = qp

    rows, simulations = [], {}
    for name, model in models.items():
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                p = protocol.loc[("Old July BoL", rate, direction)]
                z0 = priority.monotone_voltage_inverse(
                    model["v_qocv"], model["soc"], float(p.Rest_end_V)
                )
                local_stage = priority.stage_at_soc(stage, z0, charge)
                sim_rate = float(p.Measured_current_A) / cross.NOMINAL_CAPACITY_AH
                print(f"{name} | {rate:g}C {direction}", flush=True)
                sim = ehq.run_dfn(
                    sim_rate, charge, local_stage, model, anode, cathode, False, False
                )
                simulations[(name, rate, charge)] = sim
                for cell, branches in cell_curves.items():
                    metric = priority.curve_metrics(
                        branches[(rate, charge)], sim, qcell, charge
                    )
                    rows.append(
                        {
                            "Case": name,
                            "Cell": cell,
                            "C_rate": rate,
                            "Direction": direction,
                            **metric,
                            "Capacity_error_mAh": 1000.0
                            * (metric["Q_end_model_Ah"] - metric["Q_end_exp_Ah"]),
                        }
                    )
    detail = pd.DataFrame(rows)
    detail.to_csv(RESULTS / "detail.csv", index=False)
    summary = detail.groupby("Case", sort=False, as_index=False).agg(
        Full_RMSE_mV=("Full_RMSE_mV", "mean"),
        Center10_70_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
        Capacity_RMSE_pct=(
            "Capacity_error_pct", lambda x: float(np.sqrt(np.mean(np.asarray(x) ** 2)))
        ),
        Mean_abs_capacity_error_mAh=(
            "Capacity_error_mAh", lambda x: float(np.mean(np.abs(x)))
        ),
        Max_abs_capacity_error_mAh=(
            "Capacity_error_mAh", lambda x: float(np.max(np.abs(x)))
        ),
    )
    summary.to_csv(RESULTS / "summary.csv", index=False)
    branch = detail.groupby(["Case", "C_rate", "Direction"], sort=False, as_index=False).agg(
        Full_RMSE_mV=("Full_RMSE_mV", "mean"),
        Center10_70_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
        Capacity_error_pct=("Capacity_error_pct", "mean"),
        Capacity_error_mAh=("Capacity_error_mAh", "mean"),
    )
    branch.to_csv(RESULTS / "branch_summary.csv", index=False)

    colors = {"Ai2020 paper": "#777777", "Thickness gauge": "#0072B2", "SEM cross-section": "#D55E00"}
    fig, axes = plt.subplots(2, 3, figsize=(16.2, 8.8), constrained_layout=True)
    for row_i, charge in enumerate((True, False)):
        for col_i, rate in enumerate(ga.RATES):
            ax = axes[row_i, col_i]
            direction = "Charge" if charge else "Discharge"
            for branches in cell_curves.values():
                obs = branches[(rate, charge)]
                ax.plot(obs["Q_Ah"], obs["V"], color="black", lw=0.9, alpha=0.27)
            for name in raw_study.CASES:
                sim = simulations[(name, rate, charge)]
                b = branch[(branch.Case == name) & np.isclose(branch.C_rate, rate) & (branch.Direction == direction)].iloc[0]
                ax.plot(sim["Q_Ah"], sim["V"], color=colors[name], lw=1.8,
                        label=f"{name}: {b.Center10_70_RMSE_mV:.1f} mV, {b.Capacity_error_mAh:+.1f} mAh")
            ax.set(title=f"{rate:g}C {direction}", xlabel="Transferred capacity [Ah]", ylabel="Voltage [V]")
            ax.grid(alpha=0.2)
    axes[0, 0].legend(fontsize=6.5)
    fig.suptitle("Measured loading + re-identified OCP window + P0 anchors")
    fig.savefig(RESULTS / "capacity_voltage.png", dpi=220)
    plt.close(fig)

    old_qocv = c50.predict(np.array(list(old_stage.values())), model0["soc"], anode, cathode)
    fig, ax = plt.subplots(figsize=(8.8, 5.2), constrained_layout=True)
    ax.plot(model0["soc"], 1000 * (old_qocv - target), label="P0 old window")
    ax.plot(model0["soc"], 1000 * (qocv - target), label="Loading-fixed widths, re-aligned")
    ax.axhline(0, color="black", lw=0.8)
    ax.set(xlabel="SOC", ylabel="qOCV error [mV]", title="OCP-window consistency check")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(RESULTS / "qocv_window_error.png", dpi=220)
    plt.close(fig)

    report = {
        "old_stage": old_stage,
        "new_stage": {k: float(v) for k, v in stage.items()},
        "window_metrics": window_metrics,
        "mass_ratio": mass_study.MASS_RATIO,
        "method": "Qn/Qp scaled by measured single-face coating-mass ratios; window widths fixed by Qcell/Q; x0 and y100 re-aligned to C/50 qOCV with <=10 mV endpoint constraints",
    }
    (RESULTS / "analysis.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print("\nWINDOW\n", json.dumps(report, ensure_ascii=False, indent=2))
    print("\nSUMMARY\n", summary.to_string(index=False))
    print("\nBRANCH\n", branch.to_string(index=False))


if __name__ == "__main__":
    main()
