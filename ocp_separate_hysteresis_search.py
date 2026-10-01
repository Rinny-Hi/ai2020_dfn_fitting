"""Test whether separate electrode hysteresis scaling is identifiable and useful.

Calibration uses only legacy half-cell OCP and full-cell C/20 curves.  High-rate
0.5C/1C/2C data are held out for validation.  Delta-x and delta-y remain fixed.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize

import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import ocp_grounded_candidate_search as common
import old_notebook_method_comparison as omc


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260920_ocp_separate_hysteresis_search"
RESULTS.mkdir(parents=True, exist_ok=True)


def curves_sep(
    x0: float,
    y100: float,
    negative_scale: float,
    positive_scale: float,
    model: dict[str, Any],
    anode: dict[str, Any],
    cathode: dict[str, Any],
) -> dict[str, np.ndarray | float]:
    soc = model["soc"]
    x = x0 + soc * model["delta_x"]
    y0 = y100 + model["delta_y"]
    y = y0 - soc * model["delta_y"]
    un_eq = np.interp(x, anode["grid"], anode["equilibrium"])
    up_eq = np.interp(y, cathode["grid"], cathode["equilibrium"])
    un_lith = un_eq + negative_scale * (
        np.interp(x, anode["grid"], anode["lithiation"]) - un_eq
    )
    un_delith = un_eq + negative_scale * (
        np.interp(x, anode["grid"], anode["delithiation"]) - un_eq
    )
    up_lith = up_eq + positive_scale * (
        np.interp(y, cathode["grid"], cathode["lithiation"]) - up_eq
    )
    up_delith = up_eq + positive_scale * (
        np.interp(y, cathode["grid"], cathode["delithiation"]) - up_eq
    )
    return {
        "equilibrium": up_eq - un_eq,
        "charge": up_delith - un_lith,
        "discharge": up_lith - un_delith,
        "x0": x0,
        "x100": x0 + model["delta_x"],
        "y100": y100,
        "y0": y0,
    }


def fit_candidate(
    endpoint_tolerance_mV: float,
    regularization: float,
    endpoint_stage: dict[str, float],
    model: dict[str, Any],
    anode: dict[str, Any],
    cathode: dict[str, Any],
    measured_charge: np.ndarray,
    measured_discharge: np.ndarray,
    fit_mask: np.ndarray,
    fixed_window: bool = False,
) -> np.ndarray:
    x_bounds = (
        max(1e-6, float(anode["grid"].min())),
        min(1.0 - model["delta_x"] - 1e-6, float(anode["grid"].max() - model["delta_x"])),
    )
    y_bounds = (
        max(1e-6, float(cathode["grid"].min())),
        min(1.0 - model["delta_y"] - 1e-6, float(cathode["grid"].max() - model["delta_y"])),
    )

    def objective(z: np.ndarray) -> float:
        if fixed_window:
            x0, y100, sn, sp = endpoint_stage["x0"], endpoint_stage["y100"], z[0], z[1]
        else:
            x0, y100, sn, sp = z
        c = curves_sep(x0, y100, sn, sp, model, anode, cathode)
        ec = np.asarray(c["charge"])[fit_mask] - measured_charge[fit_mask]
        ed = np.asarray(c["discharge"])[fit_mask] - measured_discharge[fit_mask]
        return float(0.5 * (np.mean(ec**2) + np.mean(ed**2)) + regularization * (sn - sp) ** 2)

    def endpoint_errors(z: np.ndarray) -> np.ndarray:
        c = curves_sep(z[0], z[1], z[2], z[3], model, anode, cathode)
        return np.array([
            float(np.asarray(c["equilibrium"])[0] - model["v_qocv"][0]),
            float(np.asarray(c["equilibrium"])[-1] - model["v_qocv"][-1]),
        ])

    if fixed_window:
        best = None
        for sn0 in (0.1, 0.4, 0.7, 1.0):
            for sp0 in (0.1, 0.4, 0.7, 1.0):
                result = minimize(objective, [sn0, sp0], method="L-BFGS-B", bounds=[(0, 1), (0, 1)])
                if best is None or result.fun < best.fun:
                    best = result
        return np.array([endpoint_stage["x0"], endpoint_stage["y100"], *best.x], dtype=float)

    tolerance = endpoint_tolerance_mV / 1000.0
    constraints = [
        {"type": "ineq", "fun": lambda z: tolerance - endpoint_errors(z)[0]},
        {"type": "ineq", "fun": lambda z: tolerance + endpoint_errors(z)[0]},
        {"type": "ineq", "fun": lambda z: tolerance - endpoint_errors(z)[1]},
        {"type": "ineq", "fun": lambda z: tolerance + endpoint_errors(z)[1]},
    ]
    best = None
    for sn0 in (0.15, 0.4, 0.7):
        for sp0 in (0.15, 0.4, 0.7):
            result = minimize(
                objective,
                [endpoint_stage["x0"], endpoint_stage["y100"], sn0, sp0],
                method="SLSQP",
                bounds=[x_bounds, y_bounds, (0, 1), (0, 1)],
                constraints=constraints,
                options={"maxiter": 2000, "ftol": 1e-14},
            )
            feasible = np.max(np.abs(endpoint_errors(result.x))) <= tolerance + 1e-7
            if feasible and (best is None or result.fun < best.fun):
                best = result
    if best is None:
        raise RuntimeError("No feasible separate-scale candidate")
    return np.asarray(best.x, dtype=float)


def candidate_metrics(
    name: str,
    z: np.ndarray,
    regularization: float,
    model: dict[str, Any],
    anode: dict[str, Any],
    cathode: dict[str, Any],
    measured_charge: np.ndarray,
    measured_discharge: np.ndarray,
) -> tuple[dict[str, Any], dict[str, np.ndarray | float]]:
    x0, y100, sn, sp = z
    c = curves_sep(x0, y100, sn, sp, model, anode, cathode)
    mask = (model["soc"] >= 0.02) & (model["soc"] <= 0.98)
    eq = (np.asarray(c["equilibrium"]) - model["v_qocv"]) * 1000
    ec = (np.asarray(c["charge"]) - measured_charge) * 1000
    ed = (np.asarray(c["discharge"]) - measured_discharge) * 1000
    return {
        "Candidate": name,
        "x0": x0,
        "x100": c["x100"],
        "y100": y100,
        "y0": c["y0"],
        "negative_scale": sn,
        "positive_scale": sp,
        "regularization": regularization,
        "qOCV_MAE_2_98_mV": float(np.mean(np.abs(eq[mask]))),
        "qOCV_max_endpoint_abs_mV": float(max(abs(eq[0]), abs(eq[-1]))),
        "C20_charge_MAE_2_98_mV": float(np.mean(np.abs(ec[mask]))),
        "C20_discharge_MAE_2_98_mV": float(np.mean(np.abs(ed[mask]))),
        "C20_branch_mean_MAE_2_98_mV": float(0.5 * (np.mean(np.abs(ec[mask])) + np.mean(np.abs(ed[mask])))),
        "inside_halfcell_coverage": bool(
            x0 >= anode["grid"].min()
            and c["x100"] <= anode["grid"].max()
            and y100 >= cathode["grid"].min()
            and c["y0"] <= cathode["grid"].max()
        ),
    }, c


def scale_jacobian_condition(
    z: np.ndarray,
    model: dict[str, Any],
    anode: dict[str, Any],
    cathode: dict[str, Any],
) -> tuple[float, list[float]]:
    mask = (model["soc"] >= 0.02) & (model["soc"] <= 0.98)
    eps = 1e-5
    columns = []
    for index in (2, 3):
        zp, zm = z.copy(), z.copy()
        zp[index] += eps
        zm[index] -= eps
        cp = curves_sep(*zp, model, anode, cathode)
        cm = curves_sep(*zm, model, anode, cathode)
        derivative = np.r_[
            (np.asarray(cp["charge"])[mask] - np.asarray(cm["charge"])[mask]) / (2 * eps),
            (np.asarray(cp["discharge"])[mask] - np.asarray(cm["discharge"])[mask]) / (2 * eps),
        ]
        columns.append(derivative)
    singular = np.linalg.svd(np.column_stack(columns), compute_uv=False)
    return float(singular[0] / singular[-1]), singular.tolist()


def main() -> None:
    paths = ga.get_paths(ROOT)
    bundle = ga.load_legacy_bundle(paths)
    model = ga.build_legacy_model_inputs(bundle)
    anode = ehq.build_anode_detail(bundle, model)
    cathode = ehq.build_cathode_detail(bundle, model)
    model["y_exp"] = cathode["grid"]
    model["up_exp"] = cathode["equilibrium"]
    measured_charge, measured_discharge = common.fullcell_c20_curves(bundle, model["soc"])
    endpoint_stage, _ = omc.endpoint_window(model, {"x": anode["grid"], "equilibrium": anode["equilibrium"]})
    mask = (model["soc"] >= 0.02) & (model["soc"] <= 0.98)

    candidates: list[tuple[str, np.ndarray, float]] = []
    for reg, label in ((0.0, "free"), (1e-5, "weak-reg"), (1e-4, "medium-reg"), (1e-3, "strong-reg")):
        z = fit_candidate(10, reg, endpoint_stage, model, anode, cathode, measured_charge, measured_discharge, mask)
        candidates.append((f"Separate scales, endpoint <=10 mV ({label})", z, reg))
    z_fixed = fit_candidate(0, 0.0, endpoint_stage, model, anode, cathode, measured_charge, measured_discharge, mask, fixed_window=True)
    candidates.append(("Separate scales, exact endpoint window", z_fixed, 0.0))

    # Add the common-scale selected result as a directly comparable control.
    prior = pd.read_csv(ROOT / "results" / "260920_ocp_grounded_candidate_search" / "candidate_integrated_summary.csv")
    control = prior[prior["Candidate"] == "C20-balanced, endpoint <= 10 mV"].iloc[0]
    z_control = np.array([control.x0, control.y100, control.hysteresis_scale, control.hysteresis_scale], dtype=float)
    candidates.append(("Common scale control, endpoint <=10 mV", z_control, np.nan))

    rows, curve_map = [], {}
    for name, z, reg in candidates:
        row, c = candidate_metrics(name, z, reg, model, anode, cathode, measured_charge, measured_discharge)
        condition, singular = scale_jacobian_condition(z, model, anode, cathode)
        row["scale_jacobian_condition_number"] = condition
        row["scale_jacobian_singular_1"] = singular[0]
        row["scale_jacobian_singular_2"] = singular[1]
        rows.append(row)
        curve_map[name] = c

    # Split-SOC stability check: each half independently estimates both scales.
    split_rows = []
    for region, fit_mask in (
        ("SOC 2-50%", (model["soc"] >= 0.02) & (model["soc"] <= 0.50)),
        ("SOC 50-98%", (model["soc"] >= 0.50) & (model["soc"] <= 0.98)),
    ):
        z = fit_candidate(10, 0.0, endpoint_stage, model, anode, cathode, measured_charge, measured_discharge, fit_mask)
        split_rows.append({"Fit_region": region, "x0": z[0], "y100": z[1], "negative_scale": z[2], "positive_scale": z[3]})
    pd.DataFrame(split_rows).to_csv(RESULTS / "separate_scale_split_soc_stability.csv", index=False, encoding="utf-8-sig")

    experiment = omc.load_old_dynamic_data(model["dynamic"], model["q_meas"])
    dynamic_rows, simulations = [], {}
    for name, z, _ in candidates:
        stage = {"x0": z[0], "x100": z[0] + model["delta_x"], "y100": z[1], "y0": z[1] + model["delta_y"]}
        anode_scaled = common.scaled_detail(anode, z[2])
        cathode_scaled = common.scaled_detail(cathode, z[3])
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                print(f"{name} | {rate:g}C | {direction}", flush=True)
                sim = ehq.run_dfn(rate, charge, stage, model, anode_scaled, cathode_scaled, True, True)
                simulations[(name, rate, charge)] = sim
                metric = omc.old_time_metrics(experiment[(rate, charge)], sim)
                metric["capacity_error_Ah"] = metric["end_time_error_min"] * rate * 2.28 / 60.0
                dynamic_rows.append({"Candidate": name, "C_rate": rate, "Direction": direction, **metric})

    ocp = pd.DataFrame(rows)
    dynamic = pd.DataFrame(dynamic_rows)
    summary_dyn = dynamic.groupby("Candidate", as_index=False).agg(
        Dynamic_MAE_SOC10_70_mV=("MAE_SOC10_70_mV", "mean"),
        Mean_abs_capacity_error_Ah=("capacity_error_Ah", lambda x: np.mean(np.abs(x))),
        Max_abs_capacity_error_Ah=("capacity_error_Ah", lambda x: np.max(np.abs(x))),
    )
    by_direction = dynamic.groupby(["Candidate", "Direction"], as_index=False).agg(
        MAE_SOC10_70_mV=("MAE_SOC10_70_mV", "mean"),
        Mean_capacity_error_Ah=("capacity_error_Ah", "mean"),
    )
    charge = by_direction[by_direction.Direction == "Charge"].drop(columns="Direction").rename(columns={"MAE_SOC10_70_mV": "Charge_MAE_SOC10_70_mV", "Mean_capacity_error_Ah": "Charge_mean_capacity_error_Ah"})
    discharge = by_direction[by_direction.Direction == "Discharge"].drop(columns="Direction").rename(columns={"MAE_SOC10_70_mV": "Discharge_MAE_SOC10_70_mV", "Mean_capacity_error_Ah": "Discharge_mean_capacity_error_Ah"})
    summary = ocp.merge(summary_dyn, on="Candidate").merge(charge, on="Candidate").merge(discharge, on="Candidate")
    summary["admissible"] = (
        summary["inside_halfcell_coverage"]
        & (summary["qOCV_max_endpoint_abs_mV"] <= 10.0 + 1e-6)
        & (summary["C20_branch_mean_MAE_2_98_mV"] <= 15.0)
        & (summary["Mean_abs_capacity_error_Ah"] <= 0.05)
    )
    admissible = summary[summary.admissible].sort_values(["Dynamic_MAE_SOC10_70_mV", "C20_branch_mean_MAE_2_98_mV"])
    selected_name = admissible.iloc[0].Candidate if not admissible.empty else summary.sort_values("C20_branch_mean_MAE_2_98_mV").iloc[0].Candidate

    ocp.to_csv(RESULTS / "separate_scale_ocp_metrics.csv", index=False, encoding="utf-8-sig")
    dynamic.to_csv(RESULTS / "separate_scale_dynamic_validation.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(RESULTS / "separate_scale_integrated_summary.csv", index=False, encoding="utf-8-sig")

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
    axes[0].plot(model["soc"], model["v_qocv"], "k", lw=2.4, label="Measured qOCV")
    axes[1].plot(model["soc"], measured_charge, "k", lw=2.4, label="Measured C/20 charge")
    axes[2].plot(model["soc"], measured_discharge, "k", lw=2.4, label="Measured C/20 discharge")
    colors = plt.cm.plasma(np.linspace(0.05, 0.9, len(candidates)))
    color_map = {candidates[i][0]: colors[i] for i in range(len(candidates))}
    for name, *_ in candidates:
        c = curve_map[name]
        label = name + (" (selected)" if name == selected_name else "")
        lw = 2.7 if name == selected_name else 1.2
        axes[0].plot(model["soc"], c["equilibrium"], color=color_map[name], lw=lw, label=label)
        axes[1].plot(model["soc"], c["charge"], color=color_map[name], lw=lw)
        axes[2].plot(model["soc"], c["discharge"], color=color_map[name], lw=lw)
    for ax, title in zip(axes, ("Equilibrium qOCV", "C/20 charge", "C/20 discharge")):
        ax.set(title=title, xlabel="SOC", ylabel="Voltage [V]")
        ax.grid(alpha=0.25)
    axes[0].legend(fontsize=6.5)
    fig.suptitle("Separate electrode hysteresis-scale candidates")
    fig.tight_layout()
    fig.savefig(RESULTS / "separate_scale_ocp_comparison.png", dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.2))
    for col, rate in enumerate(ga.RATES):
        for row_idx, charge_flag in enumerate((True, False)):
            ax = axes[row_idx, col]
            direction = "Charge" if charge_flag else "Discharge"
            metric = dynamic[(dynamic.Candidate == selected_name) & np.isclose(dynamic.C_rate, rate) & (dynamic.Direction == direction)].iloc[0]
            sim = simulations[(selected_name, rate, charge_flag)]
            exp = experiment[(rate, charge_flag)]
            ax.plot(exp["t_min"], exp["V"], "k", lw=2.3, label="Experiment")
            ax.plot(sim["t_min"], sim["V"], color="#d62728", lw=2, label="Selected")
            ax.text(0.04, 0.06, f"MAE 10-70% = {metric.MAE_SOC10_70_mV:.1f} mV\nCapacity error = {metric.capacity_error_Ah:+.3f} Ah", transform=ax.transAxes, fontsize=9, bbox={"facecolor": "white", "alpha": 0.8})
            ax.set(title=f"{rate:g}C {direction}", xlabel="Time [min]", ylabel="Terminal voltage [V]")
            ax.grid(alpha=0.25)
            if row_idx == 0 and col == 0:
                ax.legend(fontsize=8)
    fig.suptitle(f"Held-out dynamic validation: {selected_name}")
    fig.tight_layout()
    fig.savefig(RESULTS / "selected_separate_scale_dynamic_curves.png", dpi=200)
    plt.close(fig)

    report = {
        "selected_candidate": selected_name,
        "new_260918_GITT_used": False,
        "calibration": "legacy half-cell OCP + full-cell C/20 only",
        "validation": "0.5C/1C/2C held out",
        "fixed_delta_x": model["delta_x"],
        "fixed_delta_y": model["delta_y"],
        "admissibility": {"endpoint_mV": 10, "C20_branch_MAE_mV": 15, "mean_capacity_error_Ah": 0.05},
        "split_soc_stability": split_rows,
    }
    (RESULTS / "separate_scale_report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print("\nIntegrated summary")
    print(summary.to_string(index=False))
    print("\nSplit-SOC stability")
    print(pd.DataFrame(split_rows).to_string(index=False))
    print("\nSelected:", selected_name)
    print("Saved to", RESULTS)


if __name__ == "__main__":
    main()
