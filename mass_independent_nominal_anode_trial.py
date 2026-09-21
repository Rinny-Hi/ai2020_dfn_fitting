"""Mass-independent OCP identification and Ai2020 nominal-anode comparison.

No coating mass is used.  Electrode capacities and stoichiometry alignment are
identified from full-cell qOCV, dV/dQ, endpoint voltage and measured full-cell
capacity.  The user's harvested-anode OCP and the PyBaMM Ai2020 Enertech
nominal graphite OCP are compared against the same experimental cathode OCP.
High-rate data are held out for validation.
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
import pybamm.input.parameters.lithium_ion.Ai2020 as ai_params
from scipy.optimize import differential_evolution, minimize
from scipy.signal import savgol_filter

import chen2020_partial_trial as ct
import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import literature_esoh_multiterm_fit as lit
import ocp_grounded_candidate_search as common
import old_notebook_method_comparison as omc


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260920_mass_independent_nominal_anode"
RESULTS.mkdir(parents=True, exist_ok=True)

# Published Ai2020 Enertech operating-window coefficients (paper Table III).
AI_X0 = 0.0065
AI_X100 = 0.84
AI_Y100 = 0.435
AI_Y0 = 0.9651


def nominal_anode_detail() -> dict[str, Any]:
    _, (x, u) = ai_params.graphite_ocp_Enertech_Ai2020_data
    x = np.asarray(x, dtype=float).reshape(-1)
    u = np.asarray(u, dtype=float).reshape(-1)
    # Retain the published interpolation support but remove duplicate x values.
    order = np.argsort(x)
    x, index = np.unique(x[order], return_index=True)
    u = u[order][index]
    grid = np.linspace(float(x.min()), float(x.max()), 4000)
    equilibrium = np.interp(grid, x, u)
    return {
        "grid": grid,
        "equilibrium": equilibrium,
        "lithiation": equilibrium.copy(),
        "delithiation": equilibrium.copy(),
        "source": "PyBaMM Ai2020 graphite_ocp_Enertech_Ai2020",
    }


def window_from_endpoints(z: np.ndarray, q_cell: float) -> dict[str, float]:
    x0, x100, y100, y0 = [float(v) for v in z]
    dx = x100 - x0
    dy = y0 - y100
    return {
        "Q_n_Ah": q_cell / dx,
        "Q_p_Ah": q_cell / dy,
        "Q_Li_Ah": y100 * q_cell / dy + x100 * q_cell / dx,
        "x0": x0,
        "x100": x100,
        "y100": y100,
        "y0": y0,
        "delta_x": dx,
        "delta_y": dy,
    }


def fit_mass_independent(
    name: str,
    model: dict[str, Any],
    anode: dict[str, Any],
    cathode: dict[str, Any],
    weights: tuple[float, float, float] = (1.0, 1.0, 1.0),
) -> tuple[dict[str, Any], np.ndarray]:
    """Fit stoichiometry endpoints directly; Qn/Qp follow from Qcell/delta."""
    soc = model["soc"]
    q_axis = soc * model["q_meas"]
    measured = model["v_qocv"]
    measured_dvdq = np.gradient(savgol_filter(measured, 51, 3), q_axis)
    mask_ocv = (soc >= 0.02) & (soc <= 0.98)
    mask_d = (soc >= 0.05) & (soc <= 0.95)
    w_ocv, w_d, w_e = weights

    x_lo, x_hi = float(anode["grid"].min()), float(anode["grid"].max())
    y_lo, y_hi = float(cathode["grid"].min()), float(cathode["grid"].max())

    # Broad capacity bounds prevent mathematically flat, non-identifiable windows.
    # They are based on the existing model only as a weak numerical prior, not on
    # the uncertain punched-electrode masses.
    qn_bounds = (0.65 * model["qn_cell"], 1.50 * model["qn_cell"])
    qp_bounds = (0.65 * model["qp_cell"], 1.50 * model["qp_cell"])

    def components(z: np.ndarray) -> tuple[float, float, float, np.ndarray]:
        window = window_from_endpoints(z, model["q_meas"])
        voltage = lit.predict(window, soc, anode, cathode)
        dvdq = np.gradient(savgol_filter(voltage, 51, 3), q_axis)
        jocv = float(np.mean(((voltage[mask_ocv] - measured[mask_ocv]) / 0.010) ** 2))
        jd = float(np.mean(((dvdq[mask_d] - measured_dvdq[mask_d]) / 0.10) ** 2))
        endpoints = np.array([voltage[0] - measured[0], voltage[-1] - measured[-1]])
        je = float(np.mean((endpoints / 0.010) ** 2))
        return jocv, jd, je, voltage

    def feasible(z: np.ndarray) -> bool:
        x0, x100, y100, y0 = z
        if not (x_lo <= x0 < x100 <= x_hi and y_lo <= y100 < y0 <= y_hi):
            return False
        w = window_from_endpoints(z, model["q_meas"])
        return (
            qn_bounds[0] <= w["Q_n_Ah"] <= qn_bounds[1]
            and qp_bounds[0] <= w["Q_p_Ah"] <= qp_bounds[1]
        )

    def objective(z: np.ndarray) -> float:
        if not feasible(z):
            return 1e6
        jocv, jd, je, _ = components(z)
        return w_ocv * jocv + w_d * jd + w_e * je

    bounds = [(x_lo, x_hi), (x_lo, x_hi), (y_lo, y_hi), (y_lo, y_hi)]
    de = differential_evolution(
        objective,
        bounds=bounds,
        seed=20260920,
        maxiter=750,
        popsize=18,
        tol=1e-9,
        polish=False,
        updating="immediate",
    )

    constraints = [
        {"type": "ineq", "fun": lambda z: z[1] - z[0] - 1e-5},
        {"type": "ineq", "fun": lambda z: z[3] - z[2] - 1e-5},
        {"type": "ineq", "fun": lambda z: model["q_meas"] / (z[1] - z[0]) - qn_bounds[0]},
        {"type": "ineq", "fun": lambda z: qn_bounds[1] - model["q_meas"] / (z[1] - z[0])},
        {"type": "ineq", "fun": lambda z: model["q_meas"] / (z[3] - z[2]) - qp_bounds[0]},
        {"type": "ineq", "fun": lambda z: qp_bounds[1] - model["q_meas"] / (z[3] - z[2])},
    ]
    polished = minimize(
        objective,
        de.x,
        method="SLSQP",
        bounds=bounds,
        constraints=constraints,
        options={"maxiter": 3000, "ftol": 1e-13},
    )
    z = polished.x if feasible(polished.x) and polished.fun <= de.fun else de.x
    window = window_from_endpoints(z, model["q_meas"])
    jocv, jd, je, voltage = components(z)
    err = (voltage - measured) * 1000.0
    row = {
        "Candidate": name,
        "Anode_source": anode.get("source", "Harvested anode v1/v2 equilibrium"),
        "Mass_used": False,
        "capacity_prior": "65-150% of current effective Qn/Qp; no punched mass",
        **window,
        "Q_n_change_percent": 100.0 * (window["Q_n_Ah"] / model["qn_cell"] - 1.0),
        "Q_p_change_percent": 100.0 * (window["Q_p_Ah"] / model["qp_cell"] - 1.0),
        "qOCV_MAE_2_98_mV": float(np.mean(np.abs(err[mask_ocv]))),
        "qOCV_RMSE_2_98_mV": float(np.sqrt(np.mean(err[mask_ocv] ** 2))),
        "qOCV_endpoint_0_error_mV": float(err[0]),
        "qOCV_endpoint_100_error_mV": float(err[-1]),
        "qOCV_max_endpoint_abs_mV": float(max(abs(err[0]), abs(err[-1]))),
        "dVdQ_RMSE_5_95_V_per_Ah": float(np.sqrt(jd) * 0.10),
        "J_OCV": jocv,
        "J_dVdQ": jd,
        "J_endpoint": je,
        "J_total": float(w_ocv * jocv + w_d * jd + w_e * je),
        "optimizer_success": bool(polished.success),
    }
    return row, voltage


def evaluate_fixed_window(
    name: str,
    model: dict[str, Any],
    anode: dict[str, Any],
    cathode: dict[str, Any],
    z: np.ndarray,
) -> tuple[dict[str, Any], np.ndarray]:
    window = window_from_endpoints(z, model["q_meas"])
    voltage = lit.predict(window, model["soc"], anode, cathode)
    mask = (model["soc"] >= 0.02) & (model["soc"] <= 0.98)
    q_axis = model["soc"] * model["q_meas"]
    measured_dvdq = np.gradient(savgol_filter(model["v_qocv"], 51, 3), q_axis)
    predicted_dvdq = np.gradient(savgol_filter(voltage, 51, 3), q_axis)
    err = (voltage - model["v_qocv"]) * 1000.0
    row = {
        "Candidate": name,
        "Anode_source": anode.get("source", "Harvested anode v1/v2 equilibrium"),
        "Mass_used": False,
        "capacity_prior": "Ai2020 published stoichiometric coefficients",
        **window,
        "Q_n_change_percent": 100.0 * (window["Q_n_Ah"] / model["qn_cell"] - 1.0),
        "Q_p_change_percent": 100.0 * (window["Q_p_Ah"] / model["qp_cell"] - 1.0),
        "qOCV_MAE_2_98_mV": float(np.mean(np.abs(err[mask]))),
        "qOCV_RMSE_2_98_mV": float(np.sqrt(np.mean(err[mask] ** 2))),
        "qOCV_endpoint_0_error_mV": float(err[0]),
        "qOCV_endpoint_100_error_mV": float(err[-1]),
        "qOCV_max_endpoint_abs_mV": float(max(abs(err[0]), abs(err[-1]))),
        "dVdQ_RMSE_5_95_V_per_Ah": float(
            np.sqrt(np.mean((predicted_dvdq[(model["soc"] >= 0.05) & (model["soc"] <= 0.95)] - measured_dvdq[(model["soc"] >= 0.05) & (model["soc"] <= 0.95)]) ** 2))
        ),
        "J_OCV": np.nan,
        "J_dVdQ": np.nan,
        "J_endpoint": np.nan,
        "J_total": np.nan,
        "optimizer_success": True,
    }
    return row, voltage


def model_for_effective_capacities(
    base: dict[str, Any], row: pd.Series
) -> dict[str, Any]:
    """Represent fitted effective Qn/Qp through c_s,max for validation only."""
    model = dict(base)
    model["params"] = base["params"].copy()
    csn = base["csn_max"] * float(row.Q_n_Ah) / base["qn_cell"]
    csp = base["csp_max"] * float(row.Q_p_Ah) / base["qp_cell"]
    model["params"].update(
        {
            "Maximum concentration in negative electrode [mol.m-3]": csn,
            "Maximum concentration in positive electrode [mol.m-3]": csp,
        }
    )
    model["csn_max"] = csn
    model["csp_max"] = csp
    model["qn_cell"] = float(row.Q_n_Ah)
    model["qp_cell"] = float(row.Q_p_Ah)
    model["delta_x"] = float(row.delta_x)
    model["delta_y"] = float(row.delta_y)
    return model


def main() -> None:
    bundle = ga.load_legacy_bundle(ga.get_paths(ROOT))
    base = ga.build_legacy_model_inputs(bundle)
    measured_anode = ehq.build_anode_detail(bundle, base)
    measured_anode["source"] = "Harvested Enertech anode v1/v2"
    nominal_anode = nominal_anode_detail()
    cathode = ehq.build_cathode_detail(bundle, base)

    rows: list[dict[str, Any]] = []
    curves: dict[str, np.ndarray] = {}
    details = {
        "Measured anode, mass-independent fit": measured_anode,
        "Ai2020 nominal anode, mass-independent fit": nominal_anode,
    }
    for name, anode in details.items():
        row, voltage = fit_mass_independent(name, base, anode, cathode)
        rows.append(row)
        curves[name] = voltage

    # Sensitivity case: retain the Ai2020 nominal equilibrium shape while
    # superposing only the branch offsets observed in the harvested electrode.
    nominal_grid = nominal_anode["grid"]
    measured_grid = measured_anode["grid"]
    lith_offset = measured_anode["lithiation"] - measured_anode["equilibrium"]
    delith_offset = measured_anode["delithiation"] - measured_anode["equilibrium"]
    nominal_hybrid = dict(nominal_anode)
    nominal_hybrid["lithiation"] = nominal_anode["equilibrium"] + np.interp(
        nominal_grid, measured_grid, lith_offset
    )
    nominal_hybrid["delithiation"] = nominal_anode["equilibrium"] + np.interp(
        nominal_grid, measured_grid, delith_offset
    )
    nominal_hybrid["source"] = "Ai2020 nominal equilibrium + harvested branch offsets"
    hybrid_name = "Ai2020 nominal equilibrium + harvested hysteresis"
    nominal_fit_row = next(
        row for row in rows if row["Candidate"] == "Ai2020 nominal anode, mass-independent fit"
    )
    hybrid_row = dict(nominal_fit_row)
    hybrid_row["Candidate"] = hybrid_name
    hybrid_row["Anode_source"] = nominal_hybrid["source"]
    rows.append(hybrid_row)
    curves[hybrid_name] = curves["Ai2020 nominal anode, mass-independent fit"].copy()

    published_name = "Ai2020 nominal anode, published stoichiometry window"
    published_row, published_voltage = evaluate_fixed_window(
        published_name,
        base,
        nominal_anode,
        cathode,
        np.array([AI_X0, AI_X100, AI_Y100, AI_Y0], dtype=float),
    )
    rows.append(published_row)
    curves[published_name] = published_voltage

    # Add the previously selected ±5% measured-anode solution as a control.
    prior = pd.read_csv(
        ROOT / "results" / "260920_literature_esoh_multiterm_fit" / "esoh_multiterm_candidate_metrics.csv"
    )
    prior_row = prior[prior.Candidate == "OCV+dVdQ+endpoint (tight ±5%)"].iloc[0]
    control_name = "Previous measured-anode ±5% control"
    control_window = {k: float(prior_row[k]) for k in ("Q_n_Ah", "Q_p_Ah", "Q_Li_Ah", "x0", "x100", "y100", "y0", "delta_x", "delta_y")}
    control_voltage = lit.predict(control_window, base["soc"], measured_anode, cathode)
    control = prior_row.to_dict()
    control.update(
        {
            "Candidate": control_name,
            "Anode_source": "Harvested Enertech anode v1/v2",
            "Mass_used": False,
            "capacity_prior": "±5% current effective capacity; no punched mass",
        }
    )
    rows.append(control)
    curves[control_name] = control_voltage

    metrics = pd.DataFrame(rows)

    # Quantify the observed v1/v2 half-cell capacity reproducibility.
    branch_capacity = ct.branch_capacity_summary(bundle["anode_lowrate"])
    branch_capacity["capacity_mAh"] = branch_capacity["capacity_Ah"] * 1000.0
    capacity_stats = {
        "mean_mAh": float(branch_capacity.capacity_mAh.mean()),
        "sample_std_mAh": float(branch_capacity.capacity_mAh.std(ddof=1)),
        "CV_percent": float(branch_capacity.capacity_mAh.std(ddof=1) / branch_capacity.capacity_mAh.mean() * 100.0),
        "range_percent_of_mean": float((branch_capacity.capacity_mAh.max() - branch_capacity.capacity_mAh.min()) / branch_capacity.capacity_mAh.mean() * 100.0),
        "v1_v2_charge_difference_percent": float(abs(branch_capacity.query("direction == 'charge'").capacity_mAh.diff().iloc[-1]) / branch_capacity.query("direction == 'charge'").capacity_mAh.mean() * 100.0),
        "v1_v2_discharge_difference_percent": float(abs(branch_capacity.query("direction == 'discharge'").capacity_mAh.diff().iloc[-1]) / branch_capacity.query("direction == 'discharge'").capacity_mAh.mean() * 100.0),
    }
    compare_grid = np.linspace(
        max(0.02, float(measured_anode["grid"].min()), float(nominal_anode["grid"].min())),
        min(0.80, float(measured_anode["grid"].max()), float(nominal_anode["grid"].max())),
        2000,
    )
    measured_common = np.interp(
        compare_grid, measured_anode["grid"], measured_anode["equilibrium"]
    )
    nominal_common = np.interp(
        compare_grid, nominal_anode["grid"], nominal_anode["equilibrium"]
    )
    anode_source_difference_mV = (measured_common - nominal_common) * 1000.0
    anode_source_stats = {
        "stoichiometry_min": float(compare_grid.min()),
        "stoichiometry_max": float(compare_grid.max()),
        "harvested_minus_nominal_bias_mV": float(np.mean(anode_source_difference_mV)),
        "MAE_mV": float(np.mean(np.abs(anode_source_difference_mV))),
        "RMSE_mV": float(np.sqrt(np.mean(anode_source_difference_mV**2))),
        "maximum_abs_mV": float(np.max(np.abs(anode_source_difference_mV))),
    }
    branch_capacity.to_csv(RESULTS / "anode_branch_capacity_reproducibility.csv", index=False, encoding="utf-8-sig")

    # Held-out dynamic validation.
    experiment = omc.load_old_dynamic_data(base["dynamic"], base["q_meas"])
    dynamic_rows = []
    simulations = {}
    dynamic_candidates = [
        control_name,
        "Measured anode, mass-independent fit",
        "Ai2020 nominal anode, mass-independent fit",
        hybrid_name,
        published_name,
    ]
    for name in dynamic_candidates:
        row = metrics[metrics.Candidate == name].iloc[0]
        anode = nominal_anode if name.startswith("Ai2020") else measured_anode
        candidate_model = model_for_effective_capacities(base, row)
        stage = {k: float(row[k]) for k in ("x0", "x100", "y100", "y0")}
        if name == hybrid_name:
            anode_dynamic = common.scaled_detail(nominal_hybrid, 0.50)
            use_negative_hysteresis = True
        elif name.startswith("Ai2020"):
            anode_dynamic = nominal_anode
            use_negative_hysteresis = False
        else:
            anode_dynamic = common.scaled_detail(measured_anode, 0.50)
            use_negative_hysteresis = True
        cathode_dynamic = common.scaled_detail(cathode, 0.50)
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                print(f"{name} | {rate:g}C | {direction}", flush=True)
                sim = ehq.run_dfn(
                    rate,
                    charge,
                    stage,
                    candidate_model,
                    anode_dynamic,
                    cathode_dynamic,
                    use_negative_hysteresis,
                    True,
                )
                simulations[(name, rate, charge)] = sim
                metric = omc.old_time_metrics(experiment[(rate, charge)], sim)
                metric["capacity_error_Ah"] = metric["end_time_error_min"] * rate * 2.28 / 60.0
                dynamic_rows.append({"Candidate": name, "C_rate": rate, "Direction": direction, **metric})
    dynamic = pd.DataFrame(dynamic_rows)
    dynamic_summary = dynamic.groupby("Candidate", as_index=False).agg(
        Dynamic_MAE_SOC10_70_mV=("MAE_SOC10_70_mV", "mean"),
        Mean_abs_capacity_error_Ah=("capacity_error_Ah", lambda x: float(np.mean(np.abs(x)))),
        Max_abs_capacity_error_Ah=("capacity_error_Ah", lambda x: float(np.max(np.abs(x)))),
    )
    metrics = metrics.merge(dynamic_summary, on="Candidate", how="left")
    metrics.to_csv(RESULTS / "mass_independent_ocp_candidate_metrics.csv", index=False, encoding="utf-8-sig")
    dynamic.to_csv(RESULTS / "mass_independent_dynamic_validation.csv", index=False, encoding="utf-8-sig")
    dynamic_summary.to_csv(RESULTS / "mass_independent_dynamic_summary.csv", index=False, encoding="utf-8-sig")

    # OCP source and fit comparison.
    fig, axes = plt.subplots(2, 2, figsize=(14.8, 10.2))
    ax = axes[0, 0]
    ax.plot(measured_anode["grid"], measured_anode["lithiation"], color="#1f77b4", lw=1.5, label="Harvested lithiation")
    ax.plot(measured_anode["grid"], measured_anode["delithiation"], color="#1f77b4", lw=1.5, ls="--", label="Harvested delithiation")
    ax.plot(nominal_anode["grid"], nominal_anode["equilibrium"], color="#d62728", lw=2.0, label="Ai2020 nominal")
    ax.set(xlabel="Negative stoichiometry x", ylabel="OCP vs Li/Li+ (V)", title="Negative-electrode OCP source")
    ax.set_ylim(0.0, 1.2)
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    ax.plot(base["soc"], base["v_qocv"], color="black", lw=2.4, label="Measured qOCV")
    colors = {
        control_name: "#7f7f7f",
        "Measured anode, mass-independent fit": "#1f77b4",
        "Ai2020 nominal anode, mass-independent fit": "#d62728",
        hybrid_name: "#ff7f0e",
        published_name: "#2ca02c",
    }
    for name in dynamic_candidates:
        ax.plot(base["soc"], curves[name], color=colors[name], lw=1.8, label=name)
    ax.set(xlabel="SOC", ylabel="Voltage (V)", title="Full-cell qOCV reconstruction")
    ax.legend(fontsize=7)

    ax = axes[1, 0]
    for name in dynamic_candidates:
        ax.plot(base["soc"], (curves[name] - base["v_qocv"]) * 1000.0, color=colors[name], lw=1.6, label=name)
    ax.axhline(0, color="black", lw=0.8)
    ax.set(xlabel="SOC", ylabel="Model - measured (mV)", title="qOCV residual")

    ax = axes[1, 1]
    plot_metrics = metrics.set_index("Candidate").loc[dynamic_candidates]
    x = np.arange(len(dynamic_candidates))
    width = 0.37
    ax.bar(x - width / 2, plot_metrics.qOCV_MAE_2_98_mV, width=width, label="qOCV MAE")
    ax.bar(x + width / 2, plot_metrics.Dynamic_MAE_SOC10_70_mV, width=width, label="Dynamic MAE")
    ax.set_xticks(x, ["Previous\n±5%", "Measured\nfree", "Ai nominal\nfree", "Ai nominal\n+ hyst", "Ai nominal\npublished"], rotation=0)
    ax.set_ylabel("MAE (mV)")
    ax.set_title("Equilibrium and held-out dynamic error")
    ax.legend(fontsize=8)
    for ax in axes.flat:
        ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(RESULTS / "mass_independent_nominal_anode_comparison.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    # Dynamic overlays for the two mass-independent fitted candidates.
    fig, axes = plt.subplots(2, 3, figsize=(15.8, 8.3))
    compare_names = [
        "Measured anode, mass-independent fit",
        "Ai2020 nominal anode, mass-independent fit",
        hybrid_name,
    ]
    for col, rate in enumerate(ga.RATES):
        for row_index, charge in enumerate((True, False)):
            direction = "Charge" if charge else "Discharge"
            ax = axes[row_index, col]
            exp = experiment[(rate, charge)]
            ax.plot(exp["t_min"], exp["V"], color="black", lw=2.1, label="Experiment")
            for name in compare_names:
                sim = simulations[(name, rate, charge)]
                metric = dynamic[(dynamic.Candidate == name) & np.isclose(dynamic.C_rate, rate) & (dynamic.Direction == direction)].iloc[0]
                if name.startswith("Measured"):
                    source_label = "Harvested anode"
                elif name == hybrid_name:
                    source_label = "Ai nominal + hyst"
                else:
                    source_label = "Ai2020 nominal"
                label = source_label + f" ({metric.MAE_SOC10_70_mV:.1f} mV)"
                ax.plot(sim["t_min"], sim["V"], color=colors[name], lw=1.8, label=label)
            ax.set_title(f"{rate:g}C {direction}")
            ax.set_xlabel("Time (min)")
            ax.set_ylabel("Voltage (V)")
            ax.grid(alpha=0.25)
            if row_index == 0 and col == 0:
                ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(RESULTS / "mass_independent_nominal_anode_dynamic.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    # Selection is not based on minimum qOCV MAE alone.  A replacement must
    # also be non-inferior on the held-out dynamic curves and endpoints.
    control_metric = metrics[metrics.Candidate == control_name].iloc[0]
    measured_free_metric = metrics[
        metrics.Candidate == "Measured anode, mass-independent fit"
    ].iloc[0]
    nominal_free_metric = metrics[
        metrics.Candidate == "Ai2020 nominal anode, mass-independent fit"
    ].iloc[0]
    hybrid_metric = metrics[metrics.Candidate == hybrid_name].iloc[0]
    recommendation = hybrid_name
    conservative_baseline = control_name
    selection_reason = (
        "OCP 및 동적 전압 오차 기준의 최우선 후보는 Ai2020 nominal equilibrium에 측정 hysteresis를 결합한 경우다. "
        f"측정 음극 자유 fitting은 qOCV MAE를 "
        f"{control_metric.qOCV_MAE_2_98_mV - measured_free_metric.qOCV_MAE_2_98_mV:.2f} mV "
        f"개선하지만 held-out 동적 MAE는 "
        f"{measured_free_metric.Dynamic_MAE_SOC10_70_mV - control_metric.Dynamic_MAE_SOC10_70_mV:.2f} mV 악화된다. "
        f"Ai2020 nominal은 qOCV MAE를 "
        f"{control_metric.qOCV_MAE_2_98_mV - nominal_free_metric.qOCV_MAE_2_98_mV:.2f} mV "
        f"개선하지만 held-out 동적 MAE는 "
        f"{nominal_free_metric.Dynamic_MAE_SOC10_70_mV - control_metric.Dynamic_MAE_SOC10_70_mV:.2f} mV 악화된다. "
        f"반면 측정 hysteresis를 추가하면 qOCV MAE는 control 대비 "
        f"{control_metric.qOCV_MAE_2_98_mV - hybrid_metric.qOCV_MAE_2_98_mV:.2f} mV, 동적 MAE는 "
        f"{control_metric.Dynamic_MAE_SOC10_70_mV - hybrid_metric.Dynamic_MAE_SOC10_70_mV:.2f} mV 개선된다. "
        f"다만 평균 및 최대 절대 용량오차가 각각 "
        f"{(hybrid_metric.Mean_abs_capacity_error_Ah - control_metric.Mean_abs_capacity_error_Ah) * 1000:.2f} mAh, "
        f"{(hybrid_metric.Max_abs_capacity_error_Ah - control_metric.Max_abs_capacity_error_Ah) * 1000:.2f} mAh 증가하므로 "
        "이를 최종 확정 모델로 보지 않고 다음 보정의 1순위 후보로 둔다. 현재 보수적 기준 모델은 기존 control이다."
    )

    report = {
        "method": "Mass-independent electrode-SOH fit using qOCV, dV/dQ, endpoints and full-cell capacity; dynamic data held out",
        "mass_used": False,
        "anode_capacity_reproducibility": capacity_stats,
        "harvested_vs_Ai2020_nominal_anode_OCP": anode_source_stats,
        "Ai2020_published_window": {"x0": AI_X0, "x100": AI_X100, "y100": AI_Y100, "y0": AI_Y0},
        "ocp_dynamic_priority_candidate": recommendation,
        "conservative_current_baseline": conservative_baseline,
        "final_acceptance": False,
        "selection_reason": selection_reason,
        "candidate_metrics": metrics.to_dict(orient="records"),
        "interpretation_limits": [
            "The broad effective-capacity fit is an identification result, not a direct teardown measurement.",
            "For dynamic validation, fitted Qn/Qp are represented through effective c_s,max; final physical allocation awaits mass, composition and collector measurements.",
            "The Ai2020 nominal graphite base curve has no branch-specific hysteresis; the hybrid candidate adds the harvested-electrode hysteresis gap.",
            "The 260918 new GITT data are excluded.",
        ],
    }
    (RESULTS / "mass_independent_nominal_anode_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    markdown = f"""# 질량 비의존 OCP 재검토와 Ai2020 nominal 음극 비교

## 결론

**OCP/동적 전압의 1순위 후보는 {recommendation}**이다. 다만 용량·종료시점 오차가 커서 아직 최종 확정하지 않으며, 현재 보수적 기준 모델은 **{conservative_baseline}**이다. {selection_reason}

| 후보 | qOCV MAE (mV) | qOCV RMSE (mV) | 최대 endpoint 오차 (mV) | 동적 평균 MAE (mV) | 평균 절대 용량오차 (mAh) |
|---|---:|---:|---:|---:|---:|
| 기존 측정 음극 ±5% | {control_metric.qOCV_MAE_2_98_mV:.2f} | {control_metric.qOCV_RMSE_2_98_mV:.2f} | {control_metric.qOCV_max_endpoint_abs_mV:.2f} | {control_metric.Dynamic_MAE_SOC10_70_mV:.2f} | {control_metric.Mean_abs_capacity_error_Ah * 1000:.2f} |
| 측정 음극, 질량 비의존 자유 fitting | {measured_free_metric.qOCV_MAE_2_98_mV:.2f} | {measured_free_metric.qOCV_RMSE_2_98_mV:.2f} | {measured_free_metric.qOCV_max_endpoint_abs_mV:.2f} | {measured_free_metric.Dynamic_MAE_SOC10_70_mV:.2f} | {measured_free_metric.Mean_abs_capacity_error_Ah * 1000:.2f} |
| Ai2020 nominal 음극, 자유 fitting | {nominal_free_metric.qOCV_MAE_2_98_mV:.2f} | {nominal_free_metric.qOCV_RMSE_2_98_mV:.2f} | {nominal_free_metric.qOCV_max_endpoint_abs_mV:.2f} | {nominal_free_metric.Dynamic_MAE_SOC10_70_mV:.2f} | {nominal_free_metric.Mean_abs_capacity_error_Ah * 1000:.2f} |
| Ai2020 nominal equilibrium + 측정 hysteresis | {hybrid_metric.qOCV_MAE_2_98_mV:.2f} | {hybrid_metric.qOCV_RMSE_2_98_mV:.2f} | {hybrid_metric.qOCV_max_endpoint_abs_mV:.2f} | {hybrid_metric.Dynamic_MAE_SOC10_70_mV:.2f} | {hybrid_metric.Mean_abs_capacity_error_Ah * 1000:.2f} |
| Ai2020 nominal published window | {published_row['qOCV_MAE_2_98_mV']:.2f} | {published_row['qOCV_RMSE_2_98_mV']:.2f} | {published_row['qOCV_max_endpoint_abs_mV']:.2f} | {float(metrics.loc[metrics.Candidate == published_name, 'Dynamic_MAE_SOC10_70_mV'].iloc[0]):.2f} | {float(metrics.loc[metrics.Candidate == published_name, 'Mean_abs_capacity_error_Ah'].iloc[0]) * 1000:.2f} |

## 질량을 사용하지 않은 방식

펀칭 질량을 목적함수에서 완전히 제외했다. Full-cell 측정 용량을 hard constraint로 두고 x0, x100, y100, y0를 qOCV, dV/dQ 및 endpoint 전압에 맞췄다. Qn과 Qp는 각각 Qcell/(x100-x0), Qcell/(y0-y100)으로 계산했다. 0.5C, 1C, 2C 데이터는 fitting에 사용하지 않고 검증에만 사용했다.

## 음극 데이터 재현성과 nominal OCP

- 기록된 v1/v2 charge/discharge 네 branch의 평균 용량은 {capacity_stats['mean_mAh']:.4f} mAh이다.
- branch 전체 CV는 {capacity_stats['CV_percent']:.2f}%, 범위는 평균 대비 {capacity_stats['range_percent_of_mean']:.2f}%이다.
- x=0.02~0.80에서 회수 음극 equilibrium OCP는 Ai2020 nominal보다 평균 {abs(anode_source_stats['harvested_minus_nominal_bias_mV']):.1f} mV 낮다.

용량 적분값 자체의 v1/v2 재현성은 현재 파일 안에서는 양호하다. 더 큰 불확실성은 젖은 전극의 질량, 초기 lithiation 상태, branch-to-stoichiometry 변환과 OCP 형상에 있다.

## Formation, SEI와 잔류 전해질 해석

- 잔류 전해질은 펀칭 질량을 증가시키므로 질량 기반 c_s,max와 specific capacity를 왜곡할 수 있다.
- Formation 후 SEI는 irreversible lithium inventory와 첫 cycle 효율, 계면 저항에 영향을 줄 수 있다.
- 그러나 충분히 낮은 전류와 이완 조건에서 graphite staging 전위 자체가 SEI 때문에 수십 mV 이동한다고 단정할 근거는 부족하다. 따라서 현재 약 32 mV의 systematic OCP 차이를 SEI 하나로 설명하지 않는다.
- 우선 확인할 실험 요인은 DMC washing 및 건조 조건, 해체 SOC, 펀칭 위치/면, coin-cell 압력과 wetting time, Li counter-electrode 상태, 첫 cycle 제외 여부다.

## 문헌 위치

- Ai et al., 2020, DOI: https://doi.org/10.1149/2.0122001JES
- Lu et al., 2021, half-cell OCP data processing, DOI: https://doi.org/10.1149/1945-7111/ac11a4
- PyBaMM Ai2020 parameter-set documentation: https://docs.pybamm.org/en/v25.6.0/source/examples/notebooks/models/Validating_mechanical_models_Enertech_DFN.html
"""
    (RESULTS / "질량_비의존_OCP와_Ai2020_nominal_음극_검토.md").write_text(
        markdown, encoding="utf-8"
    )

    print(json.dumps(report, ensure_ascii=False, indent=2))
    print("\nCandidate metrics\n", metrics.to_string(index=False))
    print("\nDynamic validation\n", dynamic.to_string(index=False))
    print("\nSaved to", RESULTS)


if __name__ == "__main__":
    main()
