"""Chen et al. (2020)-style physical parameterization trial.

The negative maximum concentration is calculated from the measured practical
half-cell areal capacity, electrode thickness and active-material fraction,
following Eq. 7 of Chen et al. The positive maximum concentration is calculated
from the measured one-sided coating mass, LiCoO2 molar mass, electrode thickness
and active-material fraction, following Eq. 5. Qn and Qp are then fixed; only
electrode alignment is fitted to full-cell qOCV with OCV, dV/dQ and endpoint
terms.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.signal import savgol_filter

import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import literature_esoh_multiterm_fit as lit
import ocp_grounded_candidate_search as common
import old_notebook_method_comparison as omc


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260920_chen2020_partial_trial"
RESULTS.mkdir(parents=True, exist_ok=True)

# User-provided one-sided coating masses for a 10 mm punched disk.
POSITIVE_COATING_MASS_KG = 17.25e-6
NEGATIVE_COATING_MASS_KG = 10.40e-6
LCO_MOLAR_MASS_KG_PER_MOL = 97.873e-3
POSITIVE_LITHIUM_PER_FORMULA = 1.0


def branch_capacity_summary(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for source in ("v1", "v2"):
        for direction in ("charge", "discharge"):
            q, _ = ga._extract_qv_branch(frame, "source_version", source, direction)
            rows.append({"source": source, "direction": direction, "capacity_Ah": float(q.max())})
    return pd.DataFrame(rows)


def fit_fixed_capacities(model, anode, cathode, qn, qp):
    soc = model["soc"]
    q_axis = soc * model["q_meas"]
    measured = model["v_qocv"]
    measured_dvdq = np.gradient(savgol_filter(measured, 51, 3), q_axis)
    mask_ocv = (soc >= 0.02) & (soc <= 0.98)
    mask_d = (soc >= 0.05) & (soc <= 0.95)
    dx = model["q_meas"] / qn
    dy = model["q_meas"] / qp

    def window(z):
        x0, y100 = [float(v) for v in z]
        return {
            "Q_n_Ah": qn,
            "Q_p_Ah": qp,
            "x0": x0,
            "x100": x0 + dx,
            "y100": y100,
            "y0": y100 + dy,
            "delta_x": dx,
            "delta_y": dy,
        }

    def components(z):
        w = window(z)
        voltage = lit.predict(w, soc, anode, cathode)
        dvdq = np.gradient(savgol_filter(voltage, 51, 3), q_axis)
        jocv = np.mean(((voltage[mask_ocv] - measured[mask_ocv]) / 0.010) ** 2)
        jd = np.mean(((dvdq[mask_d] - measured_dvdq[mask_d]) / 0.10) ** 2)
        je = np.mean((np.array([voltage[0] - measured[0], voltage[-1] - measured[-1]]) / 0.010) ** 2)
        return float(jocv), float(jd), float(je), voltage

    def objective(z):
        a, b, c, _ = components(z)
        return a + b + c

    x_bounds = (float(anode["grid"].min()), float(anode["grid"].max() - dx))
    y_bounds = (float(cathode["grid"].min()), float(cathode["grid"].max() - dy))
    if x_bounds[1] < x_bounds[0] or y_bounds[1] < y_bounds[0]:
        raise RuntimeError("Chen fixed-capacity window does not fit inside measured electrode OCP support")
    starts = [
        np.array([x_bounds[0], y_bounds[0]]),
        np.array([(x_bounds[0] + x_bounds[1]) / 2, (y_bounds[0] + y_bounds[1]) / 2]),
        np.array([x_bounds[1], y_bounds[1]]),
    ]
    best = None
    for start in starts:
        result = minimize(objective, start, method="L-BFGS-B", bounds=[x_bounds, y_bounds])
        if best is None or result.fun < best.fun:
            best = result
    w = window(best.x)
    jocv, jd, je, voltage = components(best.x)
    error = (voltage - measured) * 1000
    w["Q_Li_Ah"] = w["y100"] * qp + w["x100"] * qn
    w.update(
        {
            "qOCV_MAE_2_98_mV": float(np.mean(np.abs(error[mask_ocv]))),
            "qOCV_RMSE_2_98_mV": float(np.sqrt(np.mean(error[mask_ocv] ** 2))),
            "endpoint_0_error_mV": float(error[0]),
            "endpoint_100_error_mV": float(error[-1]),
            "max_endpoint_abs_mV": float(max(abs(error[0]), abs(error[-1]))),
            "dVdQ_RMSE_5_95_V_per_Ah": float(np.sqrt(jd) * 0.10),
            "J_OCV": jocv,
            "J_dVdQ": jd,
            "J_endpoint": je,
            "J_total": float(best.fun),
            "optimizer_success": bool(best.success),
        }
    )
    return w, voltage


def main():
    bundle = ga.load_legacy_bundle(ga.get_paths(ROOT))
    base = ga.build_legacy_model_inputs(bundle)
    params = base["params"]
    disk_diameter = float(bundle["prior"]("disk_diameter_m"))
    disk_area = np.pi * (disk_diameter / 2) ** 2
    cell_area = (
        float(params["Electrode height [m]"])
        * float(params["Electrode width [m]"])
        * float(params["Number of electrodes connected in parallel to make a cell"])
    )
    capacities = branch_capacity_summary(bundle["anode_lowrate"])
    qn_disk = float(capacities.capacity_Ah.mean())
    eps_n = float(params["Negative electrode active material volume fraction"])
    length_n = float(params["Negative electrode thickness [m]"])
    csn_chen = qn_disk * 3600 / (ga.FARADAY_CONSTANT * eps_n * length_n * disk_area)
    qn_chen = qn_disk * cell_area / disk_area
    eps_p = float(params["Positive electrode active material volume fraction"])
    length_p = float(params["Positive electrode thickness [m]"])
    positive_coating_mass_areal = POSITIVE_COATING_MASS_KG / disk_area
    csp_chen = (
        positive_coating_mass_areal
        * POSITIVE_LITHIUM_PER_FORMULA
        / (LCO_MOLAR_MASS_KG_PER_MOL * length_p * eps_p)
    )
    qp_chen = ga.FARADAY_CONSTANT * csp_chen * eps_p * length_p * cell_area / 3600

    chen = dict(base)
    chen["params"] = params.copy()
    chen["params"].update(
        {
            "Maximum concentration in negative electrode [mol.m-3]": csn_chen,
            "Maximum concentration in positive electrode [mol.m-3]": csp_chen,
        }
    )
    chen["csn_max"] = csn_chen
    chen["csp_max"] = csp_chen
    chen["qn_cell"] = qn_chen
    chen["qp_cell"] = qp_chen
    chen["delta_x"] = chen["q_meas"] / qn_chen
    chen["delta_y"] = chen["q_meas"] / chen["qp_cell"]
    anode_chen = ehq.build_anode_detail(bundle, chen)
    cathode = ehq.build_cathode_detail(bundle, chen)
    chen["y_exp"] = cathode["grid"]
    chen["up_exp"] = cathode["equilibrium"]

    fit, chen_voltage = fit_fixed_capacities(
        chen, anode_chen, cathode, chen["qn_cell"], chen["qp_cell"]
    )

    # Reconstruct previous control and selected OCP curves for a direct comparison.
    previous = pd.read_csv(
        ROOT / "results" / "260920_literature_esoh_multiterm_fit" / "esoh_multiterm_candidate_metrics.csv"
    )
    base_anode = ehq.build_anode_detail(bundle, base)
    base_cathode = ehq.build_cathode_detail(bundle, base)
    comparison_rows = []
    curves = {"Chen 2020 partial": chen_voltage}
    for label, candidate in (
        ("Existing fixed", "Existing fixed Qn/Qp endpoint control"),
        ("Previous ±5%", "OCV+dVdQ+endpoint (tight ±5%)"),
    ):
        row = previous[previous.Candidate == candidate].iloc[0]
        window = {k: float(row[k]) for k in ("Q_n_Ah", "Q_p_Ah", "x0", "x100", "y100", "y0", "delta_x", "delta_y")}
        curves[label] = lit.predict(window, base["soc"], base_anode, base_cathode)
        comparison_rows.append(
            {
                "Candidate": label,
                "Q_n_Ah": float(row.Q_n_Ah),
                "Q_p_Ah": float(row.Q_p_Ah),
                "qOCV_MAE_2_98_mV": float(row.qOCV_MAE_2_98_mV),
                "qOCV_RMSE_2_98_mV": float(row.qOCV_RMSE_2_98_mV),
                "max_endpoint_abs_mV": float(row.qOCV_max_endpoint_abs_mV),
            }
        )
    comparison_rows.append({"Candidate": "Chen 2020 partial", **{k: fit[k] for k in ("Q_n_Ah", "Q_p_Ah", "qOCV_MAE_2_98_mV", "qOCV_RMSE_2_98_mV", "max_endpoint_abs_mV")}})
    comparison = pd.DataFrame(comparison_rows)

    # Held-out dynamic validation for the Chen candidate.
    experiment = omc.load_old_dynamic_data(chen["dynamic"], chen["q_meas"])
    validation_rows = []
    simulations = {}
    stage = {k: fit[k] for k in ("x0", "x100", "y100", "y0")}
    anode_scaled = common.scaled_detail(anode_chen, 0.50)
    cathode_scaled = common.scaled_detail(cathode, 0.50)
    for rate in ga.RATES:
        for charge in (True, False):
            direction = "Charge" if charge else "Discharge"
            print(f"Chen partial | {rate:g}C | {direction}", flush=True)
            sim = ehq.run_dfn(rate, charge, stage, chen, anode_scaled, cathode_scaled, True, True)
            simulations[(rate, charge)] = sim
            metric = omc.old_time_metrics(experiment[(rate, charge)], sim)
            metric["capacity_error_Ah"] = metric["end_time_error_min"] * rate * 2.28 / 60.0
            validation_rows.append({"C_rate": rate, "Direction": direction, **metric})
    validation = pd.DataFrame(validation_rows)
    validation.to_csv(RESULTS / "chen2020_partial_dynamic_validation.csv", index=False, encoding="utf-8-sig")
    dynamic_summary = {
        "Dynamic_MAE_SOC10_70_mV": float(validation.MAE_SOC10_70_mV.mean()),
        "Mean_abs_capacity_error_Ah": float(validation.capacity_error_Ah.abs().mean()),
        "Max_abs_capacity_error_Ah": float(validation.capacity_error_Ah.abs().max()),
    }

    comparison["Dynamic_MAE_SOC10_70_mV"] = np.nan
    previous_dynamic = pd.read_csv(
        ROOT / "results" / "260920_literature_esoh_multiterm_fit" / "esoh_multiterm_dynamic_summary.csv"
    )
    comparison.loc[comparison.Candidate == "Existing fixed", "Dynamic_MAE_SOC10_70_mV"] = float(
        previous_dynamic[previous_dynamic.Candidate == "Existing fixed Qn/Qp endpoint control"].Dynamic_MAE_SOC10_70_mV.iloc[0]
    )
    comparison.loc[comparison.Candidate == "Previous ±5%", "Dynamic_MAE_SOC10_70_mV"] = float(
        previous_dynamic[previous_dynamic.Candidate == "OCV+dVdQ+endpoint (tight ±5%)"].Dynamic_MAE_SOC10_70_mV.iloc[0]
    )
    comparison.loc[comparison.Candidate == "Chen 2020 partial", "Dynamic_MAE_SOC10_70_mV"] = dynamic_summary["Dynamic_MAE_SOC10_70_mV"]
    comparison.to_csv(RESULTS / "chen2020_partial_comparison.csv", index=False, encoding="utf-8-sig")

    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.2))
    for col, rate in enumerate(ga.RATES):
        for row_index, charge in enumerate((True, False)):
            ax = axes[row_index, col]
            exp = experiment[(rate, charge)]
            sim = simulations[(rate, charge)]
            direction = "Charge" if charge else "Discharge"
            metric = validation[(validation.C_rate == rate) & (validation.Direction == direction)].iloc[0]
            ax.plot(exp["t_min"], exp["V"], color="k", lw=2, label="Experiment")
            ax.plot(sim["t_min"], sim["V"], color="#2ca02c", lw=1.8, label="Chen 2020 trial")
            ax.set_title(f"{rate:g}C {direction}\nMAE(10–70%)={metric.MAE_SOC10_70_mV:.1f} mV")
            ax.set_xlabel("Time (min)")
            ax.set_ylabel("Voltage (V)")
            ax.grid(alpha=0.25)
            if row_index == 0 and col == 0:
                ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(RESULTS / "chen2020_partial_dynamic_validation.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    soc = base["soc"]
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
    axes[0].plot(soc, base["v_qocv"], "k", lw=2.5, label="Measured qOCV")
    colors = {"Existing fixed": "#7f7f7f", "Previous ±5%": "#d62728", "Chen 2020 partial": "#2ca02c"}
    for name, voltage in curves.items():
        axes[0].plot(soc, voltage, lw=2, color=colors[name], label=name)
        axes[1].plot(soc, (voltage - base["v_qocv"]) * 1000, lw=1.8, color=colors[name], label=name)
    axes[0].set(xlabel="SOC", ylabel="Voltage (V)", title="qOCV reconstruction")
    axes[1].axhline(0, color="k", lw=0.8)
    axes[1].set(xlabel="SOC", ylabel="Model - measured (mV)", title="qOCV residual")
    axes[0].legend(fontsize=8)
    names = comparison.Candidate.tolist()
    x = np.arange(len(names))
    axes[2].bar(x - 0.18, comparison.qOCV_MAE_2_98_mV, width=0.36, label="qOCV MAE")
    axes[2].bar(x + 0.18, comparison.Dynamic_MAE_SOC10_70_mV, width=0.36, label="Dynamic MAE")
    axes[2].set_xticks(x, names, rotation=20, ha="right")
    axes[2].set_ylabel("MAE (mV)")
    axes[2].set_title("Equilibrium vs held-out dynamic error")
    axes[2].legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(RESULTS / "chen2020_partial_comparison.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    report = {
        "method": "Chen 2020 Eq. 7 for negative c_s_max and Eq. 5 for positive c_s_max",
        "assumption": "17.25 mg positive and 10.40 mg negative are one-sided coating masses for a 10 mm disk; positive coating is treated as 100% LiCoO2, matching Chen's primary calculation convention",
        "anode_halfcell_capacity_Ah_per_10mm_disk": qn_disk,
        "anode_areal_capacity_mAh_cm2": qn_disk * 1000 / (disk_area * 1e4),
        "positive_coating_mass_mg_per_10mm_disk": POSITIVE_COATING_MASS_KG * 1e6,
        "negative_coating_mass_mg_per_10mm_disk": NEGATIVE_COATING_MASS_KG * 1e6,
        "negative_specific_capacity_mAh_g_coating": qn_disk * 1000 / (NEGATIVE_COATING_MASS_KG * 1000),
        "positive_coating_mass_mg_cm2": positive_coating_mass_areal * 100,
        "cell_effective_area_cm2": cell_area * 1e4,
        "negative_c_s_max_original_mol_m3": float(base["csn_max"]),
        "negative_c_s_max_chen_mol_m3": csn_chen,
        "positive_c_s_max_original_mol_m3": float(base["csp_max"]),
        "positive_c_s_max_chen_mol_m3": csp_chen,
        "fit": fit,
        "dynamic_summary": dynamic_summary,
        "limitations": [
            "Positive coating is assumed to be 100% LiCoO2; binder/conductive-additive fraction is not measured",
            "No three-electrode full-cell OCV; terminal qOCV difference fitting is used instead",
            "Cell effective area and layer count remain inherited from the current model",
        ],
    }
    (RESULTS / "chen2020_partial_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    print("\nComparison\n", comparison.to_string(index=False))
    print("\nDynamic validation\n", validation.to_string(index=False))
    print("Saved to", RESULTS)


if __name__ == "__main__":
    main()
