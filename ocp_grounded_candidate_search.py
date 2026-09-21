"""Search OCP-grounded stoichiometry/hysteresis candidates.

Only legacy half-cell OCP and full-cell C/20 data are used for calibration.
The 0.5C/1C/2C branches are reserved for validation.  Delta x and delta y are
fixed by measured full-cell capacity for every candidate.
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
from scipy.optimize import minimize, minimize_scalar

import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import old_notebook_method_comparison as omc


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260920_ocp_grounded_candidate_search"
RESULTS.mkdir(parents=True, exist_ok=True)
CENTRAL = slice(None)


def fullcell_c20_curves(
    bundle: dict[str, Any], soc: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    full = bundle["full_c20"]
    q_charge, v_charge = ga._extract_qv_branch(full, "source_version", "v2", "charge")
    q_discharge, v_discharge = ga._extract_qv_branch(full, "source_version", "v2", "discharge")
    measured_charge = ga._interp_curve(
        ga._clean_curve(q_charge / q_charge.max(), v_charge), soc
    )
    measured_discharge = ga._interp_curve(
        ga._clean_curve(1.0 - q_discharge / q_discharge.max(), v_discharge), soc
    )
    return measured_charge, measured_discharge


def scaled_detail(detail: dict[str, Any], scale: float) -> dict[str, Any]:
    output = dict(detail)
    equilibrium = np.asarray(detail["equilibrium"], dtype=float)
    output["lithiation"] = equilibrium + scale * (
        np.asarray(detail["lithiation"], dtype=float) - equilibrium
    )
    output["delithiation"] = equilibrium + scale * (
        np.asarray(detail["delithiation"], dtype=float) - equilibrium
    )
    return output


def ocp_curves(
    x0: float,
    y100: float,
    scale: float,
    model: dict[str, Any],
    anode: dict[str, Any],
    cathode: dict[str, Any],
) -> dict[str, np.ndarray | float]:
    soc = model["soc"]
    x100 = x0 + model["delta_x"]
    y0 = y100 + model["delta_y"]
    x = x0 + soc * model["delta_x"]
    y = y0 - soc * model["delta_y"]

    un_eq = np.interp(x, anode["grid"], anode["equilibrium"])
    up_eq = np.interp(y, cathode["grid"], cathode["equilibrium"])
    un_lith = np.interp(x, anode["grid"], anode["lithiation"])
    un_delith = np.interp(x, anode["grid"], anode["delithiation"])
    up_lith = np.interp(y, cathode["grid"], cathode["lithiation"])
    up_delith = np.interp(y, cathode["grid"], cathode["delithiation"])

    un_lith = un_eq + scale * (un_lith - un_eq)
    un_delith = un_eq + scale * (un_delith - un_eq)
    up_lith = up_eq + scale * (up_lith - up_eq)
    up_delith = up_eq + scale * (up_delith - up_eq)

    return {
        "equilibrium": up_eq - un_eq,
        "charge": up_delith - un_lith,
        "discharge": up_lith - un_delith,
        "x0": x0,
        "x100": x100,
        "y100": y100,
        "y0": y0,
    }


def metrics_for_candidate(
    name: str,
    x0: float,
    y100: float,
    scale: float,
    model: dict[str, Any],
    anode: dict[str, Any],
    cathode: dict[str, Any],
    measured_charge: np.ndarray,
    measured_discharge: np.ndarray,
) -> tuple[dict[str, Any], dict[str, np.ndarray | float]]:
    curves = ocp_curves(x0, y100, scale, model, anode, cathode)
    mask = (model["soc"] >= 0.02) & (model["soc"] <= 0.98)
    eq_error = (np.asarray(curves["equilibrium"]) - model["v_qocv"]) * 1000.0
    charge_error = (np.asarray(curves["charge"]) - measured_charge) * 1000.0
    discharge_error = (np.asarray(curves["discharge"]) - measured_discharge) * 1000.0
    measured_gap = (measured_charge - measured_discharge) * 1000.0
    predicted_gap = (
        np.asarray(curves["charge"]) - np.asarray(curves["discharge"])
    ) * 1000.0
    x0_lower = float(anode["grid"].min())
    x0_upper = float(anode["grid"].max() - model["delta_x"])
    y100_lower = float(cathode["grid"].min())
    y100_upper = float(cathode["grid"].max() - model["delta_y"])
    return {
        "Candidate": name,
        "x0": x0,
        "x100": float(curves["x100"]),
        "y100": y100,
        "y0": float(curves["y0"]),
        "hysteresis_scale": scale,
        "qOCV_MAE_2_98_mV": float(np.mean(np.abs(eq_error[mask]))),
        "qOCV_RMSE_2_98_mV": float(np.sqrt(np.mean(eq_error[mask] ** 2))),
        "qOCV_SOC0_error_mV": float(eq_error[0]),
        "qOCV_SOC100_error_mV": float(eq_error[-1]),
        "qOCV_max_endpoint_abs_mV": float(max(abs(eq_error[0]), abs(eq_error[-1]))),
        "C20_charge_MAE_2_98_mV": float(np.mean(np.abs(charge_error[mask]))),
        "C20_discharge_MAE_2_98_mV": float(np.mean(np.abs(discharge_error[mask]))),
        "C20_branch_mean_MAE_2_98_mV": float(
            0.5
            * (
                np.mean(np.abs(charge_error[mask]))
                + np.mean(np.abs(discharge_error[mask]))
            )
        ),
        "measured_C20_hysteresis_MAE_2_98_mV": float(
            np.mean(np.abs(measured_gap[mask]))
        ),
        "predicted_C20_hysteresis_MAE_2_98_mV": float(
            np.mean(np.abs(predicted_gap[mask]))
        ),
        "x0_distance_to_bound": float(min(x0 - x0_lower, x0_upper - x0)),
        "y100_distance_to_bound": float(
            min(y100 - y100_lower, y100_upper - y100)
        ),
        "inside_halfcell_coverage": bool(
            x0 >= x0_lower
            and curves["x100"] <= anode["grid"].max()
            and y100 >= y100_lower
            and curves["y0"] <= cathode["grid"].max()
        ),
    }, curves


def fit_scale_for_fixed_window(
    x0: float,
    y100: float,
    model: dict[str, Any],
    anode: dict[str, Any],
    cathode: dict[str, Any],
    measured_charge: np.ndarray,
    measured_discharge: np.ndarray,
) -> float:
    mask = (model["soc"] >= 0.02) & (model["soc"] <= 0.98)

    def objective(scale: float) -> float:
        curves = ocp_curves(x0, y100, scale, model, anode, cathode)
        e_charge = np.asarray(curves["charge"])[mask] - measured_charge[mask]
        e_discharge = np.asarray(curves["discharge"])[mask] - measured_discharge[mask]
        return float(np.mean(e_charge**2) + np.mean(e_discharge**2))

    result = minimize_scalar(objective, bounds=(0.0, 1.0), method="bounded")
    return float(result.x)


def solve_balanced_candidate(
    endpoint_tolerance_mV: float,
    initial: tuple[float, float, float],
    model: dict[str, Any],
    anode: dict[str, Any],
    cathode: dict[str, Any],
    measured_charge: np.ndarray,
    measured_discharge: np.ndarray,
) -> tuple[np.ndarray, bool, str]:
    mask = (model["soc"] >= 0.02) & (model["soc"] <= 0.98)
    x_bounds = (
        max(1e-6, float(anode["grid"].min())),
        min(1.0 - model["delta_x"] - 1e-6, float(anode["grid"].max() - model["delta_x"])),
    )
    y_bounds = (
        max(1e-6, float(cathode["grid"].min())),
        min(1.0 - model["delta_y"] - 1e-6, float(cathode["grid"].max() - model["delta_y"])),
    )

    def objective(z: np.ndarray) -> float:
        curves = ocp_curves(z[0], z[1], z[2], model, anode, cathode)
        e_charge = np.asarray(curves["charge"])[mask] - measured_charge[mask]
        e_discharge = np.asarray(curves["discharge"])[mask] - measured_discharge[mask]
        return float(0.5 * (np.mean(e_charge**2) + np.mean(e_discharge**2)))

    def endpoint_errors(z: np.ndarray) -> np.ndarray:
        curves = ocp_curves(z[0], z[1], z[2], model, anode, cathode)
        return np.array(
            [
                float(np.asarray(curves["equilibrium"])[0] - model["v_qocv"][0]),
                float(np.asarray(curves["equilibrium"])[-1] - model["v_qocv"][-1]),
            ]
        )

    tolerance = endpoint_tolerance_mV / 1000.0
    if endpoint_tolerance_mV == 0:
        constraints = [
            {"type": "eq", "fun": lambda z: endpoint_errors(z)[0]},
            {"type": "eq", "fun": lambda z: endpoint_errors(z)[1]},
        ]
    else:
        constraints = [
            {"type": "ineq", "fun": lambda z: tolerance - endpoint_errors(z)[0]},
            {"type": "ineq", "fun": lambda z: tolerance + endpoint_errors(z)[0]},
            {"type": "ineq", "fun": lambda z: tolerance - endpoint_errors(z)[1]},
            {"type": "ineq", "fun": lambda z: tolerance + endpoint_errors(z)[1]},
        ]

    best = None
    starts = [initial]
    for scale_start in (0.1, 0.3, 0.5, 0.7, 0.9):
        starts.append((initial[0], initial[1], scale_start))
    for start in starts:
        result = minimize(
            objective,
            x0=np.asarray(start, dtype=float),
            method="SLSQP",
            bounds=[x_bounds, y_bounds, (0.0, 1.0)],
            constraints=constraints,
            options={"maxiter": 1500, "ftol": 1e-14, "disp": False},
        )
        endpoint = endpoint_errors(result.x)
        feasible = (
            np.max(np.abs(endpoint)) <= tolerance + 1e-7
            if endpoint_tolerance_mV > 0
            else np.max(np.abs(endpoint)) <= 1e-7
        )
        if feasible and (best is None or result.fun < best.fun):
            best = result
    if best is None:
        raise RuntimeError(f"No feasible solution for endpoint tolerance {endpoint_tolerance_mV} mV")
    return best.x, bool(best.success), str(best.message)


def main() -> None:
    paths = ga.get_paths(ROOT)
    bundle = ga.load_legacy_bundle(paths)
    model = ga.build_legacy_model_inputs(bundle)
    anode = ehq.build_anode_detail(bundle, model)
    cathode = ehq.build_cathode_detail(bundle, model)
    model["y_exp"] = cathode["grid"]
    model["up_exp"] = cathode["equilibrium"]
    measured_charge, measured_discharge = fullcell_c20_curves(bundle, model["soc"])

    endpoint_stage, _ = omc.endpoint_window(
        model, {"x": anode["grid"], "equilibrium": anode["equilibrium"]}
    )
    endpoint_scale = fit_scale_for_fixed_window(
        endpoint_stage["x0"],
        endpoint_stage["y100"],
        model,
        anode,
        cathode,
        measured_charge,
        measured_discharge,
    )
    old_stage, _ = omc.solve_old_window(
        model, {"x": anode["grid"], "equilibrium": anode["equilibrium"]}
    )
    old_scale = fit_scale_for_fixed_window(
        old_stage["x0"],
        old_stage["y100"],
        model,
        anode,
        cathode,
        measured_charge,
        measured_discharge,
    )

    candidate_parameters: list[tuple[str, float, float, float, str]] = [
        (
            "Endpoint + raw branch gap",
            float(endpoint_stage["x0"]),
            float(endpoint_stage["y100"]),
            1.0,
            "reference",
        ),
        (
            "Endpoint + C20-calibrated gap",
            float(endpoint_stage["x0"]),
            float(endpoint_stage["y100"]),
            endpoint_scale,
            "C20 calibration",
        ),
    ]
    initial = (float(endpoint_stage["x0"]), float(endpoint_stage["y100"]), endpoint_scale)
    solver_notes = {}
    for tolerance in (5.0, 10.0, 20.0):
        solution, success, message = solve_balanced_candidate(
            tolerance,
            initial,
            model,
            anode,
            cathode,
            measured_charge,
            measured_discharge,
        )
        name = f"C20-balanced, endpoint <= {tolerance:g} mV"
        candidate_parameters.append(
            (name, float(solution[0]), float(solution[1]), float(solution[2]), "C20 calibration")
        )
        solver_notes[name] = {"success": success, "message": message}
    candidate_parameters.append(
        (
            "Old qOCV window + C20-calibrated gap",
            float(old_stage["x0"]),
            float(old_stage["y100"]),
            old_scale,
            "reference",
        )
    )

    candidate_rows = []
    candidate_curves = {}
    for name, x0, y100, scale, role in candidate_parameters:
        row, curves = metrics_for_candidate(
            name,
            x0,
            y100,
            scale,
            model,
            anode,
            cathode,
            measured_charge,
            measured_discharge,
        )
        row["calibration_role"] = role
        candidate_rows.append(row)
        candidate_curves[name] = curves
    ocp_table = pd.DataFrame(candidate_rows)
    ocp_table.to_csv(RESULTS / "candidate_ocp_metrics.csv", index=False, encoding="utf-8-sig")

    experiment = omc.load_old_dynamic_data(model["dynamic"], model["q_meas"])
    dynamic_rows = []
    simulations = {}
    for name, x0, y100, scale, _ in candidate_parameters:
        stage = {"x0": x0, "x100": x0 + model["delta_x"], "y100": y100, "y0": y100 + model["delta_y"]}
        anode_scaled = scaled_detail(anode, scale)
        cathode_scaled = scaled_detail(cathode, scale)
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                print(f"{name} | {rate:g}C | {direction}", flush=True)
                sim = ehq.run_dfn(
                    rate,
                    charge,
                    stage,
                    model,
                    anode_scaled,
                    cathode_scaled,
                    negative_hysteresis=True,
                    positive_hysteresis=True,
                )
                simulations[(name, rate, charge)] = sim
                metric = omc.old_time_metrics(experiment[(rate, charge)], sim)
                metric["capacity_error_Ah"] = (
                    metric["end_time_error_min"] * rate * 2.28 / 60.0
                )
                dynamic_rows.append(
                    {
                        "Candidate": name,
                        "C_rate": rate,
                        "Direction": direction,
                        **metric,
                    }
                )
    dynamic = pd.DataFrame(dynamic_rows)
    dynamic.to_csv(RESULTS / "candidate_dynamic_validation.csv", index=False, encoding="utf-8-sig")
    dynamic_summary = (
        dynamic.groupby("Candidate", as_index=False)
        .agg(
            Dynamic_MAE_all_mV=("MAE_all_mV", "mean"),
            Dynamic_MAE_SOC10_70_mV=("MAE_SOC10_70_mV", "mean"),
            Mean_abs_end_time_error_min=("end_time_error_min", lambda x: np.mean(np.abs(x))),
            Mean_abs_capacity_error_Ah=("capacity_error_Ah", lambda x: np.mean(np.abs(x))),
            Max_abs_capacity_error_Ah=("capacity_error_Ah", lambda x: np.max(np.abs(x))),
        )
    )
    charge_summary = (
        dynamic[dynamic["Direction"] == "Charge"]
        .groupby("Candidate", as_index=False)
        .agg(
            Charge_MAE_SOC10_70_mV=("MAE_SOC10_70_mV", "mean"),
            Charge_mean_capacity_error_Ah=("capacity_error_Ah", "mean"),
        )
    )
    discharge_summary = (
        dynamic[dynamic["Direction"] == "Discharge"]
        .groupby("Candidate", as_index=False)
        .agg(
            Discharge_MAE_SOC10_70_mV=("MAE_SOC10_70_mV", "mean"),
            Discharge_mean_capacity_error_Ah=("capacity_error_Ah", "mean"),
        )
    )
    summary = (
        ocp_table.merge(dynamic_summary, on="Candidate")
        .merge(charge_summary, on="Candidate")
        .merge(discharge_summary, on="Candidate")
    )
    summary["endpoint_within_10mV"] = summary["qOCV_max_endpoint_abs_mV"] <= 10.0 + 1e-6
    summary["capacity_mean_abs_within_0p05Ah"] = summary["Mean_abs_capacity_error_Ah"] <= 0.05
    summary["physically_admissible"] = (
        summary["inside_halfcell_coverage"]
        & summary["endpoint_within_10mV"]
        & summary["capacity_mean_abs_within_0p05Ah"]
    )
    summary.to_csv(RESULTS / "candidate_integrated_summary.csv", index=False, encoding="utf-8-sig")

    admissible = summary[summary["physically_admissible"]].copy()
    if admissible.empty:
        selected_name = summary.sort_values(
            ["qOCV_max_endpoint_abs_mV", "C20_branch_mean_MAE_2_98_mV", "Mean_abs_capacity_error_Ah"]
        ).iloc[0]["Candidate"]
        selection_status = "No candidate met all prespecified admissibility checks"
    else:
        selected_name = admissible.sort_values(
            ["C20_branch_mean_MAE_2_98_mV", "Dynamic_MAE_SOC10_70_mV"]
        ).iloc[0]["Candidate"]
        selection_status = "Selected among candidates meeting endpoint, capacity, and coverage checks"

    colors = plt.cm.viridis(np.linspace(0.05, 0.95, len(candidate_parameters)))
    color_map = {item[0]: colors[index] for index, item in enumerate(candidate_parameters)}

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.1))
    axes[0].plot(model["soc"], model["v_qocv"], color="black", lw=2.5, label="Measured qOCV")
    axes[1].plot(model["soc"], measured_charge, color="black", lw=2.5, label="Measured C/20 charge")
    axes[2].plot(model["soc"], measured_discharge, color="black", lw=2.5, label="Measured C/20 discharge")
    for name, *_ in candidate_parameters:
        curves = candidate_curves[name]
        lw = 2.6 if name == selected_name else 1.3
        alpha = 1.0 if name == selected_name else 0.65
        label = f"{name}{' (selected)' if name == selected_name else ''}"
        axes[0].plot(model["soc"], curves["equilibrium"], color=color_map[name], lw=lw, alpha=alpha, label=label)
        axes[1].plot(model["soc"], curves["charge"], color=color_map[name], lw=lw, alpha=alpha, label=label)
        axes[2].plot(model["soc"], curves["discharge"], color=color_map[name], lw=lw, alpha=alpha, label=label)
    for axis, title in zip(axes, ("Equilibrium qOCV", "C/20 charge OCP", "C/20 discharge OCP")):
        axis.set(title=title, xlabel="SOC", ylabel="Cell voltage [V]")
        axis.grid(alpha=0.25)
    axes[0].legend(fontsize=7)
    fig.suptitle("OCP-grounded candidate comparison", fontsize=15)
    fig.tight_layout()
    fig.savefig(RESULTS / "candidate_ocp_comparison.png", dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.1))
    for name, *_ in candidate_parameters:
        row = summary[summary["Candidate"] == name].iloc[0]
        axes[0].scatter(row["qOCV_max_endpoint_abs_mV"], row["C20_branch_mean_MAE_2_98_mV"], color=color_map[name], s=90, label=name)
        axes[1].scatter(row["Mean_abs_capacity_error_Ah"], row["Dynamic_MAE_SOC10_70_mV"], color=color_map[name], s=90, label=name)
        axes[2].scatter(row["qOCV_MAE_2_98_mV"], row["Dynamic_MAE_SOC10_70_mV"], color=color_map[name], s=90, label=name)
    axes[0].axvline(10, color="#d62728", ls="--", lw=1.5, label="10 mV endpoint limit")
    axes[1].axvline(0.05, color="#d62728", ls="--", lw=1.5, label="0.05 Ah capacity limit")
    axes[0].set(xlabel="Maximum qOCV endpoint error [mV]", ylabel="C/20 branch mean MAE [mV]", title="Calibration quality")
    axes[1].set(xlabel="Mean absolute capacity error [Ah]", ylabel="Dynamic MAE 10–70% [mV]", title="Validation quality")
    axes[2].set(xlabel="qOCV MAE 2–98% [mV]", ylabel="Dynamic MAE 10–70% [mV]", title="OCP vs dynamic trade-off")
    for axis in axes:
        axis.grid(alpha=0.25)
    axes[0].legend(fontsize=6.7, loc="best")
    fig.suptitle("Candidate trade-offs: lower-left is better", fontsize=15)
    fig.tight_layout()
    fig.savefig(RESULTS / "candidate_tradeoffs.png", dpi=200)
    plt.close(fig)

    selected = next(item for item in candidate_parameters if item[0] == selected_name)
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.2))
    for col, rate in enumerate(ga.RATES):
        for row_index, charge in enumerate((True, False)):
            axis = axes[row_index, col]
            direction = "Charge" if charge else "Discharge"
            metric = dynamic[
                (dynamic["Candidate"] == selected_name)
                & np.isclose(dynamic["C_rate"], rate)
                & (dynamic["Direction"] == direction)
            ].iloc[0]
            sim = simulations[(selected_name, rate, charge)]
            axis.plot(experiment[(rate, charge)]["t_min"], experiment[(rate, charge)]["V"], color="black", lw=2.3, label="Experiment")
            axis.plot(sim["t_min"], sim["V"], color="#2ca02c", lw=2.0, label="Selected OCP-grounded candidate")
            axis.text(0.04, 0.06, f"MAE 10–70% = {metric['MAE_SOC10_70_mV']:.1f} mV\nCapacity error = {metric['capacity_error_Ah']:+.3f} Ah", transform=axis.transAxes, fontsize=9, bbox={"facecolor": "white", "alpha": 0.82, "edgecolor": "#cccccc"})
            axis.set(title=f"{rate:g}C {direction}", xlabel="Time [min]", ylabel="Terminal voltage [V]")
            axis.grid(alpha=0.25)
            if row_index == 0 and col == 0:
                axis.legend(fontsize=8)
    fig.suptitle(f"Selected candidate validation: {selected_name}", fontsize=14)
    fig.tight_layout()
    fig.savefig(RESULTS / "selected_candidate_dynamic_curves.png", dpi=200)
    plt.close(fig)

    report = {
        "new_260918_GITT_used": False,
        "calibration_data": ["legacy half-cell OCP", "full-cell C/20 charge/discharge"],
        "validation_data": ["0.5C", "1C", "2C", "charge", "discharge"],
        "fixed": {"delta_x": model["delta_x"], "delta_y": model["delta_y"]},
        "prespecified_admissibility": {
            "max_qocv_endpoint_error_mV": 10.0,
            "mean_abs_capacity_error_Ah": 0.05,
            "inside_halfcell_coverage": True,
        },
        "selection_status": selection_status,
        "selected_candidate": selected_name,
        "solver_notes": solver_notes,
    }
    (RESULTS / "selection_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print("\nIntegrated summary")
    print(summary.to_string(index=False))
    print("\nSelected:", selected_name)
    print(selection_status)
    print("Saved to", RESULTS)


if __name__ == "__main__":
    main()
