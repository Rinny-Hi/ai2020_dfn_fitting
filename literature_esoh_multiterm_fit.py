"""Literature-style electrode-capacity / stoichiometry-window identification.

The equilibrium OCP fit uses Qn, Qp and electrode alignment. Measured cell
capacity is enforced as a hard constraint, matching the electrode-SOH
formulation. The objective compares full-cell qOCV, dV/dQ and endpoints.
High-rate curves are held out until after OCP-based candidate selection.
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
from scipy.signal import savgol_filter

import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import ocp_grounded_candidate_search as common
import old_notebook_method_comparison as omc


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260920_literature_esoh_multiterm_fit"
RESULTS.mkdir(parents=True, exist_ok=True)


def window_from_z(z: np.ndarray, q_cell: float) -> dict[str, float]:
    qn, qp, x0, y100 = [float(v) for v in z]
    x100 = x0 + q_cell / qn
    y0 = y100 + q_cell / qp
    q_li = y100 * qp + x100 * qn
    return {
        "Q_n_Ah": qn,
        "Q_p_Ah": qp,
        "Q_Li_Ah": q_li,
        "x0": x0,
        "x100": x100,
        "y100": y100,
        "y0": y0,
        "delta_x": x100 - x0,
        "delta_y": y0 - y100,
    }


def predict(window: dict[str, float], soc: np.ndarray, anode: dict[str, Any], cathode: dict[str, Any]) -> np.ndarray:
    x = window["x0"] + soc * window["delta_x"]
    y = window["y0"] - soc * window["delta_y"]
    return np.interp(y, cathode["grid"], cathode["equilibrium"]) - np.interp(x, anode["grid"], anode["equilibrium"])


def fit_candidate(
    name: str,
    capacity_bound_fraction: float,
    weights: tuple[float, float, float],
    model: dict[str, Any],
    anode: dict[str, Any],
    cathode: dict[str, Any],
    initial: np.ndarray,
) -> tuple[dict[str, Any], np.ndarray]:
    soc = model["soc"]
    q_axis = soc * model["q_meas"]
    measured = model["v_qocv"]
    measured_smooth = savgol_filter(measured, 51, 3)
    measured_dvdq = np.gradient(measured_smooth, q_axis)
    mask_ocv = (soc >= 0.02) & (soc <= 0.98)
    mask_derivative = (soc >= 0.05) & (soc <= 0.95)
    qn0, qp0 = model["qn_cell"], model["qp_cell"]
    qn_bounds = (qn0 * (1 - capacity_bound_fraction), qn0 * (1 + capacity_bound_fraction))
    qp_bounds = (qp0 * (1 - capacity_bound_fraction), qp0 * (1 + capacity_bound_fraction))
    x_bounds = (float(anode["grid"].min()), float(anode["grid"].max()))
    y_bounds = (float(cathode["grid"].min()), float(cathode["grid"].max()))
    w_ocv, w_derivative, w_endpoint = weights

    def components(z: np.ndarray) -> tuple[float, float, float]:
        window = window_from_z(z, model["q_meas"])
        voltage = predict(window, soc, anode, cathode)
        voltage_smooth = savgol_filter(voltage, 51, 3)
        dvdq = np.gradient(voltage_smooth, q_axis)
        j_ocv = float(np.mean(((voltage[mask_ocv] - measured[mask_ocv]) / 0.010) ** 2))
        j_derivative = float(np.mean(((dvdq[mask_derivative] - measured_dvdq[mask_derivative]) / 0.10) ** 2))
        endpoint_error = np.array([voltage[0] - measured[0], voltage[-1] - measured[-1]])
        j_endpoint = float(np.mean((endpoint_error / 0.010) ** 2))
        return j_ocv, j_derivative, j_endpoint

    def objective(z: np.ndarray) -> float:
        j_ocv, j_derivative, j_endpoint = components(z)
        return w_ocv * j_ocv + w_derivative * j_derivative + w_endpoint * j_endpoint

    def coverage(z: np.ndarray) -> np.ndarray:
        w = window_from_z(z, model["q_meas"])
        return np.array([
            w["x0"] - anode["grid"].min(),
            anode["grid"].max() - w["x100"],
            w["y100"] - cathode["grid"].min(),
            cathode["grid"].max() - w["y0"],
        ])

    constraints = [{"type": "ineq", "fun": lambda z, i=i: coverage(z)[i]} for i in range(4)]
    starts = [initial]
    for qn_factor in (1 - capacity_bound_fraction, 1.0, 1 + capacity_bound_fraction):
        for qp_factor in (1 - capacity_bound_fraction, 1.0, 1 + capacity_bound_fraction):
            starts.append(np.array([qn0 * qn_factor, qp0 * qp_factor, initial[2], initial[3]]))
    best = None
    for start in starts:
        result = minimize(
            objective,
            start,
            method="SLSQP",
            bounds=[qn_bounds, qp_bounds, x_bounds, y_bounds],
            constraints=constraints,
            options={"maxiter": 2500, "ftol": 1e-12, "disp": False},
        )
        if np.min(coverage(result.x)) >= -1e-7 and (best is None or result.fun < best.fun):
            best = result
    if best is None:
        raise RuntimeError(f"No feasible solution for {name}")
    window = window_from_z(best.x, model["q_meas"])
    voltage = predict(window, soc, anode, cathode)
    smooth = savgol_filter(voltage, 51, 3)
    dvdq = np.gradient(smooth, q_axis)
    err = (voltage - measured) * 1000
    derr = dvdq - measured_dvdq
    j_ocv, j_derivative, j_endpoint = components(best.x)
    row = {
        "Candidate": name,
        "capacity_bound_percent": capacity_bound_fraction * 100,
        "w_OCV": w_ocv,
        "w_dVdQ": w_derivative,
        "w_endpoint": w_endpoint,
        **window,
        "Q_n_change_percent": 100 * (window["Q_n_Ah"] / qn0 - 1),
        "Q_p_change_percent": 100 * (window["Q_p_Ah"] / qp0 - 1),
        "qOCV_MAE_2_98_mV": float(np.mean(np.abs(err[mask_ocv]))),
        "qOCV_RMSE_2_98_mV": float(np.sqrt(np.mean(err[mask_ocv] ** 2))),
        "qOCV_endpoint_0_error_mV": float(err[0]),
        "qOCV_endpoint_100_error_mV": float(err[-1]),
        "qOCV_max_endpoint_abs_mV": float(max(abs(err[0]), abs(err[-1]))),
        "dVdQ_RMSE_5_95_V_per_Ah": float(np.sqrt(np.mean(derr[mask_derivative] ** 2))),
        "capacity_residual_n_Ah": float(window["Q_n_Ah"] * window["delta_x"] - model["q_meas"]),
        "capacity_residual_p_Ah": float(window["Q_p_Ah"] * window["delta_y"] - model["q_meas"]),
        "J_OCV": j_ocv,
        "J_dVdQ": j_derivative,
        "J_endpoint": j_endpoint,
        "J_total": float(best.fun),
        "optimizer_success": bool(best.success),
        "optimizer_message": str(best.message),
    }
    return row, voltage


def make_model_for_capacity_fit(model: dict[str, Any], row: pd.Series) -> dict[str, Any]:
    output = dict(model)
    output["params"] = model["params"].copy()
    eps_n0 = float(model["params"]["Negative electrode active material volume fraction"])
    eps_p0 = float(model["params"]["Positive electrode active material volume fraction"])
    eps_n = eps_n0 * float(row.Q_n_Ah) / model["qn_cell"]
    eps_p = eps_p0 * float(row.Q_p_Ah) / model["qp_cell"]
    output["params"].update({
        "Negative electrode active material volume fraction": eps_n,
        "Positive electrode active material volume fraction": eps_p,
    })
    output["effective_eps_n"] = eps_n
    output["effective_eps_p"] = eps_p
    output["delta_x"] = float(row.delta_x)
    output["delta_y"] = float(row.delta_y)
    return output


def main() -> None:
    paths = ga.get_paths(ROOT)
    bundle = ga.load_legacy_bundle(paths)
    model = ga.build_legacy_model_inputs(bundle)
    anode = ehq.build_anode_detail(bundle, model)
    cathode = ehq.build_cathode_detail(bundle, model)
    model["y_exp"] = cathode["grid"]
    model["up_exp"] = cathode["equilibrium"]
    measured_charge, measured_discharge = common.fullcell_c20_curves(bundle, model["soc"])
    endpoint, _ = omc.endpoint_window(model, {"x": anode["grid"], "equilibrium": anode["equilibrium"]})
    initial = np.array([model["qn_cell"], model["qp_cell"], endpoint["x0"], endpoint["y100"]], dtype=float)

    specs = [
        ("Voltage+endpoint, capacity fixed (broad)", 0.15, (1.0, 0.0, 1.0)),
        ("OCV+dVdQ+endpoint (tight ±5%)", 0.05, (1.0, 1.0, 1.0)),
        ("OCV+dVdQ+endpoint (broad ±15%)", 0.15, (1.0, 1.0, 1.0)),
        ("Derivative-heavy (broad ±15%)", 0.15, (1.0, 4.0, 1.0)),
        ("Endpoint-heavy (broad ±15%)", 0.15, (1.0, 1.0, 4.0)),
    ]
    rows, curves = [], {}
    for name, bound, weights in specs:
        row, voltage = fit_candidate(name, bound, weights, model, anode, cathode, initial)
        charge_error = (voltage - measured_charge) * 1000
        discharge_error = (voltage - measured_discharge) * 1000
        mask = (model["soc"] >= 0.02) & (model["soc"] <= 0.98)
        row["C20_charge_vs_equilibrium_MAE_2_98_mV"] = float(np.mean(np.abs(charge_error[mask])))
        row["C20_discharge_vs_equilibrium_MAE_2_98_mV"] = float(np.mean(np.abs(discharge_error[mask])))
        rows.append(row)
        curves[name] = voltage

    # Existing fixed-capacity endpoint solution as control.
    control_window = {
        "Q_n_Ah": model["qn_cell"], "Q_p_Ah": model["qp_cell"],
        "x0": endpoint["x0"], "x100": endpoint["x100"],
        "y100": endpoint["y100"], "y0": endpoint["y0"],
        "delta_x": model["delta_x"], "delta_y": model["delta_y"],
    }
    control_voltage = predict(control_window, model["soc"], anode, cathode)
    mask = (model["soc"] >= 0.02) & (model["soc"] <= 0.98)
    control_err = (control_voltage - model["v_qocv"]) * 1000
    control_row = {
        "Candidate": "Existing fixed Qn/Qp endpoint control",
        "capacity_bound_percent": 0.0,
        "w_OCV": np.nan, "w_dVdQ": np.nan, "w_endpoint": np.nan,
        **control_window,
        "Q_Li_Ah": endpoint["y100"] * model["qp_cell"] + endpoint["x100"] * model["qn_cell"],
        "Q_n_change_percent": 0.0, "Q_p_change_percent": 0.0,
        "qOCV_MAE_2_98_mV": float(np.mean(np.abs(control_err[mask]))),
        "qOCV_RMSE_2_98_mV": float(np.sqrt(np.mean(control_err[mask] ** 2))),
        "qOCV_endpoint_0_error_mV": float(control_err[0]),
        "qOCV_endpoint_100_error_mV": float(control_err[-1]),
        "qOCV_max_endpoint_abs_mV": float(max(abs(control_err[0]), abs(control_err[-1]))),
        "dVdQ_RMSE_5_95_V_per_Ah": np.nan,
        "capacity_residual_n_Ah": 0.0, "capacity_residual_p_Ah": 0.0,
        "J_OCV": np.nan, "J_dVdQ": np.nan, "J_endpoint": np.nan, "J_total": np.nan,
        "optimizer_success": True, "optimizer_message": "control",
        "C20_charge_vs_equilibrium_MAE_2_98_mV": float(np.mean(np.abs((control_voltage - measured_charge)[mask]) * 1000)),
        "C20_discharge_vs_equilibrium_MAE_2_98_mV": float(np.mean(np.abs((control_voltage - measured_discharge)[mask]) * 1000)),
    }
    rows.append(control_row)
    curves[control_row["Candidate"]] = control_voltage
    table = pd.DataFrame(rows)

    # Select only among capacity changes that stay within the ±5% measurement prior.
    admissible = table[
        (table.Candidate != "Existing fixed Qn/Qp endpoint control")
        & (table.Q_n_change_percent.abs() <= 5.0 + 1e-6)
        & (table.Q_p_change_percent.abs() <= 5.0 + 1e-6)
        & (table.qOCV_max_endpoint_abs_mV <= 10.0)
    ].copy()
    if admissible.empty:
        selected_name = "Existing fixed Qn/Qp endpoint control"
        selection_status = "No fitted candidate passed ±5% capacity and 10 mV endpoint checks"
    else:
        selected_name = admissible.sort_values(["qOCV_RMSE_2_98_mV", "dVdQ_RMSE_5_95_V_per_Ah"]).iloc[0].Candidate
        selection_status = "Selected among ±5% capacity and 10 mV endpoint-admissible fits"
    selected = table[table.Candidate == selected_name].iloc[0]
    table["selected"] = table.Candidate == selected_name
    table.to_csv(RESULTS / "esoh_multiterm_candidate_metrics.csv", index=False, encoding="utf-8-sig")

    # Held-out dynamic validation uses the same common 0.50 hysteresis scale for
    # both selected and control, so only capacity/window identification changes.
    experiment = omc.load_old_dynamic_data(model["dynamic"], model["q_meas"])
    validation_rows, simulations = [], {}
    validation_names = [
        "Existing fixed Qn/Qp endpoint control",
        selected_name,
        "OCV+dVdQ+endpoint (broad ±15%)",
    ]
    for candidate_name in dict.fromkeys(validation_names):
        row = table[table.Candidate == candidate_name].iloc[0]
        candidate_model = make_model_for_capacity_fit(model, row)
        stage = {"x0": row.x0, "x100": row.x100, "y100": row.y100, "y0": row.y0}
        anode_scaled = common.scaled_detail(anode, 0.50)
        cathode_scaled = common.scaled_detail(cathode, 0.50)
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                print(f"{candidate_name} | {rate:g}C | {direction}", flush=True)
                sim = ehq.run_dfn(rate, charge, stage, candidate_model, anode_scaled, cathode_scaled, True, True)
                simulations[(candidate_name, rate, charge)] = sim
                metric = omc.old_time_metrics(experiment[(rate, charge)], sim)
                metric["capacity_error_Ah"] = metric["end_time_error_min"] * rate * 2.28 / 60.0
                validation_rows.append({"Candidate": candidate_name, "C_rate": rate, "Direction": direction, **metric})
    validation = pd.DataFrame(validation_rows)
    validation.to_csv(RESULTS / "esoh_multiterm_dynamic_validation.csv", index=False, encoding="utf-8-sig")
    summary_dynamic = validation.groupby("Candidate", as_index=False).agg(
        Dynamic_MAE_SOC10_70_mV=("MAE_SOC10_70_mV", "mean"),
        Mean_abs_capacity_error_Ah=("capacity_error_Ah", lambda x: np.mean(np.abs(x))),
        Max_abs_capacity_error_Ah=("capacity_error_Ah", lambda x: np.max(np.abs(x))),
    )
    summary_dynamic.to_csv(RESULTS / "esoh_multiterm_dynamic_summary.csv", index=False, encoding="utf-8-sig")

    # Figure 1: voltage, residual and derivative evidence.
    soc = model["soc"]
    q_axis = soc * model["q_meas"]
    measured_dvdq = np.gradient(savgol_filter(model["v_qocv"], 51, 3), q_axis)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
    axes[0].plot(soc, model["v_qocv"], "k", lw=2.5, label="Measured qOCV")
    colors = {
        "Existing fixed Qn/Qp endpoint control": "#7f7f7f",
        selected_name: "#d62728",
        "OCV+dVdQ+endpoint (broad ±15%)": "#1f77b4",
    }
    for name in dict.fromkeys(validation_names):
        voltage = curves[name]
        label = (
            "Selected Qn/Qp/QLi fit"
            if name == selected_name
            else "Broad-prior diagnostic"
            if "broad ±15%" in name
            else "Fixed-capacity control"
        )
        axes[0].plot(soc, voltage, color=colors[name], lw=2.2, label=label)
        axes[1].plot(soc, (voltage - model["v_qocv"]) * 1000, color=colors[name], lw=1.8, label=label)
        dvdq = np.gradient(savgol_filter(voltage, 51, 3), q_axis)
        axes[2].plot(soc, dvdq, color=colors[name], lw=1.8, label=label)
    axes[2].plot(soc, measured_dvdq, "k", lw=2.2, label="Measured qOCV")
    axes[0].set(title="Equilibrium full-cell OCV", xlabel="SOC", ylabel="Voltage [V]")
    axes[1].axhline(0, color="black", lw=1)
    axes[1].set(title="qOCV residual", xlabel="SOC", ylabel="Model - experiment [mV]")
    axes[2].set(title="Differential-voltage feature", xlabel="SOC", ylabel="dV/dQ [V/Ah]", ylim=(-0.05, 1.0))
    for ax in axes:
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    fig.suptitle("Literature-style Qn/Qp/QLi and stoichiometry-window identification")
    fig.tight_layout()
    fig.savefig(RESULTS / "esoh_multiterm_ocp_comparison.png", dpi=200)
    plt.close(fig)

    # Figure 2: all candidate capacity changes and fit quality.
    fitted = table[table.Candidate != "Existing fixed Qn/Qp endpoint control"].copy()
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
    axes[0].scatter(fitted.Q_n_change_percent, fitted.Q_p_change_percent, s=90)
    for _, row in fitted.iterrows():
        axes[0].annotate(row.Candidate.split(" (")[0], (row.Q_n_change_percent, row.Q_p_change_percent), xytext=(4, 3), textcoords="offset points", fontsize=7)
    axes[0].axvspan(-5, 5, color="#2ca02c", alpha=0.08)
    axes[0].axhspan(-5, 5, color="#2ca02c", alpha=0.08)
    axes[0].set(xlabel="Qn change [%]", ylabel="Qp change [%]", title="Electrode-capacity movement")
    axes[1].scatter(fitted.qOCV_RMSE_2_98_mV, fitted.dVdQ_RMSE_5_95_V_per_Ah, s=90)
    for _, row in fitted.iterrows():
        axes[1].annotate(row.Candidate.split(" (")[0], (row.qOCV_RMSE_2_98_mV, row.dVdQ_RMSE_5_95_V_per_Ah), xytext=(4, 3), textcoords="offset points", fontsize=7)
    axes[1].set(xlabel="qOCV RMSE [mV]", ylabel="dV/dQ RMSE [V/Ah]", title="Curve and feature fit")
    axes[2].scatter(fitted.qOCV_max_endpoint_abs_mV, fitted.qOCV_RMSE_2_98_mV, s=90)
    axes[2].axvline(10, color="red", ls="--", lw=1.3)
    axes[2].set(xlabel="Max endpoint error [mV]", ylabel="qOCV RMSE [mV]", title="Endpoint trade-off")
    for ax in axes:
        ax.grid(alpha=0.25)
    fig.suptitle("Sensitivity to objective weights and capacity priors")
    fig.tight_layout()
    fig.savefig(RESULTS / "esoh_multiterm_candidate_tradeoffs.png", dpi=200)
    plt.close(fig)

    # Figure 3: held-out validation selected versus fixed-capacity control.
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.2))
    for col, rate in enumerate(ga.RATES):
        for row_index, charge in enumerate((True, False)):
            ax = axes[row_index, col]
            direction = "Charge" if charge else "Discharge"
            obs = experiment[(rate, charge)]
            ax.plot(obs["t_min"], obs["V"], "k", lw=2.4, label="Experiment")
            for name in dict.fromkeys(validation_names):
                sim = simulations[(name, rate, charge)]
                metric = validation[(validation.Candidate == name) & np.isclose(validation.C_rate, rate) & (validation.Direction == direction)].iloc[0]
                display = "Selected" if name == selected_name else "Broad" if "broad ±15%" in name else "Control"
                label = f"{display} ({metric.MAE_SOC10_70_mV:.1f} mV)"
                ax.plot(sim["t_min"], sim["V"], color=colors[name], lw=1.8, label=label)
            ax.set(title=f"{rate:g}C {direction}", xlabel="Time [min]", ylabel="Terminal voltage [V]")
            ax.grid(alpha=0.25)
            if row_index == 0 and col == 0:
                ax.legend(fontsize=7.5)
    fig.suptitle("Held-out dynamic validation at common hysteresis scale 0.50")
    fig.tight_layout()
    fig.savefig(RESULTS / "esoh_multiterm_dynamic_validation.png", dpi=200)
    plt.close(fig)

    selected_model = make_model_for_capacity_fit(model, selected)
    report = {
        "selected_candidate": selected_name,
        "selection_status": selection_status,
        "new_260918_GITT_used": False,
        "capacity_treatment": "hard constraint: Q = Qn*delta_x = Qp*delta_y",
        "objective": "J = w_ocv*J_ocv + w_dvdq*J_dvdq + w_endpoint*J_endpoint; J_capacity=0 by construction",
        "selected": {k: (float(v) if isinstance(v, (np.floating, float, int)) and not isinstance(v, bool) else v) for k, v in selected.to_dict().items()},
        "effective_mapping_for_DFN_validation": {
            "c_s_max_changed": False,
            "negative_active_material_fraction": selected_model["effective_eps_n"],
            "positive_active_material_fraction": selected_model["effective_eps_p"],
        },
    }
    (RESULTS / "esoh_multiterm_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print("\nCandidate metrics")
    print(table.to_string(index=False))
    print("\nDynamic summary")
    print(summary_dynamic.to_string(index=False))
    print("\nSelected:", selected_name)
    print(selection_status)
    print("Saved to", RESULTS)


if __name__ == "__main__":
    main()
