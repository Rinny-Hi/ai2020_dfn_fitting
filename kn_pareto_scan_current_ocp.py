"""Scan kn with the current OCP and report a non-weighted Pareto comparison."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import c50_selected_ocp_dynamic_validation as selected
import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import nominal_radius_current_configuration as current
import priority_initial_state_protocol_recheck as priority


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260923_kn_pareto_current_ocp"
RESULTS.mkdir(parents=True, exist_ok=True)
NOMINAL_CAPACITY_AH = 2.28
KN_VALUES = np.array([0.40, 0.50, 0.60, 0.70, 0.74, 0.80, 0.90, 1.00, 1.20]) * 1e-6


def with_kn(model0: dict, kn: float) -> dict:
    changed = dict(model0)
    changed["params"] = model0["params"].copy()
    key = "Negative electrode exchange-current density [A.m-2]"
    changed["params"].update(
        {key: ga.scale_parameter_function(model0["params"][key], kn / current.KN_PREF)},
        check_already_exists=False,
    )
    return changed


def pareto_mask(frame: pd.DataFrame, columns: list[str]) -> np.ndarray:
    values = frame[columns].to_numpy(float)
    keep = np.ones(len(frame), dtype=bool)
    for i in range(len(frame)):
        for j in range(len(frame)):
            if i == j:
                continue
            if np.all(values[j] <= values[i]) and np.any(values[j] < values[i]):
                keep[i] = False
                break
    return keep


def main() -> None:
    model0, stage, anode, cathode, experiment, _, qcell = selected.selected_inputs()
    protocol, _, _ = priority.load_protocol(model0)
    pidx = protocol.set_index(["C_rate", "Direction"])
    rows = []
    for kn in KN_VALUES:
        model = with_kn(model0, float(kn))
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                rest_v = float(pidx.loc[(rate, direction), "Rest_end_V"])
                initial_soc = priority.monotone_voltage_inverse(model["v_qocv"], model["soc"], rest_v)
                local_stage = priority.stage_at_soc(stage, initial_soc, charge)
                measured_current = float(np.nanmedian(np.abs(experiment[(rate, charge)]["I_A"])))
                sim_rate = measured_current / NOMINAL_CAPACITY_AH
                print(f"kn={kn:.2e} | {rate:g}C {direction}", flush=True)
                sim = ehq.run_dfn(
                    sim_rate,
                    charge,
                    local_stage,
                    model,
                    anode,
                    cathode,
                    False,
                    True,
                    positive_initial_h=(-1.0 if charge else 1.0),
                )
                metric = priority.curve_metrics(experiment[(rate, charge)], sim, qcell, charge)
                rows.append(
                    {
                        "kn": kn,
                        "C_rate": rate,
                        "Direction": direction,
                        "initial_SOC": initial_soc,
                        "measured_current_A": measured_current,
                        **metric,
                        "Capacity_error_mAh": 1000 * (metric["Q_end_model_Ah"] - metric["Q_end_exp_Ah"]),
                    }
                )
    detail = pd.DataFrame(rows)
    detail.to_csv(RESULTS / "kn_scan_detail.csv", index=False, encoding="utf-8-sig")
    summary = detail.groupby(["kn", "Direction"], sort=True, as_index=False).agg(
        Mean_full_RMSE_mV=("Full_RMSE_mV", "mean"),
        Mean_center_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
        Capacity_RMSE_pct=("Capacity_error_pct", lambda x: float(np.sqrt(np.mean(x**2)))),
        Mean_capacity_error_mAh=("Capacity_error_mAh", "mean"),
        Mean_abs_capacity_error_mAh=("Capacity_error_mAh", lambda x: float(np.mean(np.abs(x)))),
    )
    summary.to_csv(RESULTS / "kn_scan_summary.csv", index=False, encoding="utf-8-sig")

    charge = summary[summary.Direction == "Charge"].copy().reset_index(drop=True)
    objectives = ["Mean_full_RMSE_mV", "Mean_center_RMSE_mV", "Capacity_RMSE_pct"]
    charge["Charge_Pareto"] = pareto_mask(charge, objectives)
    discharge = summary[summary.Direction == "Discharge"].drop(columns="Direction")
    comparison = charge.merge(discharge, on="kn", suffixes=("_Charge", "_Discharge"))
    current_index = int(np.argmin(np.abs(comparison.kn.to_numpy(float) - current.KN_PREF)))
    current_row = comparison.iloc[current_index]
    all_objectives = [
        "Mean_full_RMSE_mV_Charge",
        "Mean_center_RMSE_mV_Charge",
        "Capacity_RMSE_pct_Charge",
        "Mean_full_RMSE_mV_Discharge",
        "Mean_center_RMSE_mV_Discharge",
        "Capacity_RMSE_pct_Discharge",
    ]
    comparison["Dominates_current_EIS_all_metrics"] = [
        bool(
            np.all(row[all_objectives].to_numpy(float) <= current_row[all_objectives].to_numpy(float))
            and np.any(row[all_objectives].to_numpy(float) < current_row[all_objectives].to_numpy(float))
        )
        for _, row in comparison.iterrows()
    ]
    comparison.to_csv(RESULTS / "kn_pareto_comparison.csv", index=False, encoding="utf-8-sig")

    fig, axes = plt.subplots(2, 3, figsize=(15.8, 8.0), constrained_layout=True)
    metrics = [
        ("Mean_full_RMSE_mV", "Full-range voltage RMSE [mV]"),
        ("Mean_center_RMSE_mV", "10-70% voltage RMSE [mV]"),
        ("Capacity_RMSE_pct", "Capacity RMSE [%]"),
    ]
    for row_index, direction in enumerate(("Charge", "Discharge")):
        part = summary[summary.Direction == direction]
        for col, (metric, label) in enumerate(metrics):
            ax = axes[row_index, col]
            ax.plot(part.kn * 1e6, part[metric], marker="o", color="#0072B2")
            ax.axvline(current.KN_PREF * 1e6, color="#D55E00", ls="--", label="Current EIS kn")
            ax.set(title=f"{direction}: {label}", xlabel="kn prefactor [1e-6]", ylabel=label)
            ax.grid(alpha=0.22)
            if row_index == 0 and col == 0:
                ax.legend()
    fig.suptitle("kn scan with current OCP, rest-based initial SOC, and measured current")
    fig.savefig(RESULTS / "kn_pareto_scan.png", dpi=220)
    plt.close(fig)

    report = {
        "kn_values": KN_VALUES.tolist(),
        "kp_fixed": current.KP_PREF,
        "selection_rule": "No weighted sum. Charge Pareto on full RMSE, 10-70% RMSE, and capacity RMSE; discharge is cross-validation.",
        "current_EIS_kn": current.KN_PREF,
        "charge_pareto_kn": charge.loc[charge.Charge_Pareto, "kn"].tolist(),
        "kn_dominating_current_on_all_charge_and_discharge_metrics": comparison.loc[
            comparison.Dominates_current_EIS_all_metrics, "kn"
        ].tolist(),
        "comparison": comparison.to_dict(orient="records"),
    }
    (RESULTS / "kn_pareto_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print("\nSUMMARY")
    print(summary.to_string(index=False))
    print("\nPARETO")
    print(comparison.to_string(index=False))


if __name__ == "__main__":
    main()
