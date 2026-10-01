"""Map the OCP/dynamic/capacity Pareto trade-off for one common hysteresis scale."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import ocp_grounded_candidate_search as common
import old_notebook_method_comparison as omc


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260920_ocp_common_scale_pareto"
RESULTS.mkdir(parents=True, exist_ok=True)
SCALES = [0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.80, 1.00]


def pareto_mask(values: np.ndarray) -> np.ndarray:
    keep = np.ones(len(values), dtype=bool)
    for i, row in enumerate(values):
        dominated = np.any(np.all(values <= row, axis=1) & np.any(values < row, axis=1))
        keep[i] = not dominated
    return keep


def main() -> None:
    paths = ga.get_paths(ROOT)
    bundle = ga.load_legacy_bundle(paths)
    model = ga.build_legacy_model_inputs(bundle)
    anode = ehq.build_anode_detail(bundle, model)
    cathode = ehq.build_cathode_detail(bundle, model)
    measured_charge, measured_discharge = common.fullcell_c20_curves(bundle, model["soc"])
    prior = pd.read_csv(ROOT / "results" / "260920_ocp_grounded_candidate_search" / "candidate_integrated_summary.csv")
    base = prior[prior.Candidate == "C20-balanced, endpoint <= 10 mV"].iloc[0]
    stage = {"x0": float(base.x0), "x100": float(base.x100), "y100": float(base.y100), "y0": float(base.y0)}
    exp = omc.load_old_dynamic_data(model["dynamic"], model["q_meas"])

    ocp_rows: list[dict[str, Any]] = []
    dynamic_rows: list[dict[str, Any]] = []
    simulations: dict[tuple[float, float, bool], dict[str, np.ndarray]] = {}
    for scale in SCALES:
        name = f"Common scale {scale:.2f}"
        row, _ = common.metrics_for_candidate(name, stage["x0"], stage["y100"], scale, model, anode, cathode, measured_charge, measured_discharge)
        ocp_rows.append(row)
        anode_scaled = common.scaled_detail(anode, scale)
        cathode_scaled = common.scaled_detail(cathode, scale)
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                print(f"scale={scale:.2f} | {rate:g}C | {direction}", flush=True)
                sim = ehq.run_dfn(rate, charge, stage, model, anode_scaled, cathode_scaled, True, True)
                simulations[(scale, rate, charge)] = sim
                metric = omc.old_time_metrics(exp[(rate, charge)], sim)
                metric["capacity_error_Ah"] = metric["end_time_error_min"] * rate * 2.28 / 60.0
                dynamic_rows.append({"Candidate": name, "scale": scale, "C_rate": rate, "Direction": direction, **metric})

    ocp = pd.DataFrame(ocp_rows)
    dynamic = pd.DataFrame(dynamic_rows)
    overall = dynamic.groupby(["Candidate", "scale"], as_index=False).agg(
        Dynamic_MAE_SOC10_70_mV=("MAE_SOC10_70_mV", "mean"),
        Mean_abs_capacity_error_Ah=("capacity_error_Ah", lambda x: np.mean(np.abs(x))),
        Max_abs_capacity_error_Ah=("capacity_error_Ah", lambda x: np.max(np.abs(x))),
    )
    dirs = dynamic.groupby(["Candidate", "Direction"], as_index=False).agg(
        MAE_SOC10_70_mV=("MAE_SOC10_70_mV", "mean"),
        Mean_capacity_error_Ah=("capacity_error_Ah", "mean"),
    )
    for label in ("Charge", "Discharge"):
        part = dirs[dirs.Direction == label].drop(columns="Direction").rename(columns={
            "MAE_SOC10_70_mV": f"{label}_MAE_SOC10_70_mV",
            "Mean_capacity_error_Ah": f"{label}_mean_capacity_error_Ah",
        })
        overall = overall.merge(part, on="Candidate")
    summary = ocp.merge(overall, on="Candidate")
    summary["admissible"] = (
        (summary.qOCV_max_endpoint_abs_mV <= 10.0 + 1e-6)
        & (summary.C20_branch_mean_MAE_2_98_mV <= 20.0)
        & (summary.Mean_abs_capacity_error_Ah <= 0.05)
        & summary.inside_halfcell_coverage
    )
    criteria = summary[["C20_branch_mean_MAE_2_98_mV", "Dynamic_MAE_SOC10_70_mV", "Mean_abs_capacity_error_Ah"]].to_numpy(float)
    summary["pareto_optimal"] = pareto_mask(criteria)
    admissible = summary[summary.admissible].copy()
    if admissible.empty:
        selected = summary.sort_values(["C20_branch_mean_MAE_2_98_mV", "Mean_abs_capacity_error_Ah"]).iloc[0]
        selection_status = "No scale passed all OCP/capacity gates"
    else:
        # All admitted candidates already satisfy absolute OCP/capacity gates;
        # choose the lowest held-out voltage error, then the lower C/20 error.
        selected = admissible.sort_values(["Dynamic_MAE_SOC10_70_mV", "C20_branch_mean_MAE_2_98_mV"]).iloc[0]
        selection_status = "Lowest held-out dynamic MAE among OCP/capacity-admissible scales"
    selected_scale = float(selected.scale)
    summary.to_csv(RESULTS / "common_scale_pareto_summary.csv", index=False, encoding="utf-8-sig")
    dynamic.to_csv(RESULTS / "common_scale_pareto_dynamic_detail.csv", index=False, encoding="utf-8-sig")

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
    axes[0].plot(summary.scale, summary.C20_branch_mean_MAE_2_98_mV, "o-", label="C/20 branch MAE")
    axes[0].plot(summary.scale, summary.Dynamic_MAE_SOC10_70_mV, "o-", label="0.5C/1C/2C MAE")
    axes[0].axhline(20, color="red", ls="--", lw=1.3, label="C/20 gate")
    axes[0].set(xlabel="Common hysteresis scale", ylabel="MAE [mV]", title="OCP vs dynamic error")
    axes[0].legend(fontsize=8)
    axes[1].plot(summary.scale, summary.Mean_abs_capacity_error_Ah, "o-", label="Mean absolute")
    axes[1].plot(summary.scale, summary.Max_abs_capacity_error_Ah, "o-", label="Maximum absolute")
    axes[1].axhline(0.05, color="red", ls="--", lw=1.3, label="Mean capacity gate")
    axes[1].set(xlabel="Common hysteresis scale", ylabel="Capacity error [Ah]", title="Cutoff capacity error")
    axes[1].legend(fontsize=8)
    for _, row in summary.iterrows():
        marker = "*" if np.isclose(row.scale, selected_scale) else "o"
        size = 180 if marker == "*" else 70
        axes[2].scatter(row.C20_branch_mean_MAE_2_98_mV, row.Dynamic_MAE_SOC10_70_mV, s=size, marker=marker, c=[row.scale], cmap="viridis", vmin=min(SCALES), vmax=max(SCALES))
        axes[2].annotate(f"{row.scale:.2f}", (row.C20_branch_mean_MAE_2_98_mV, row.Dynamic_MAE_SOC10_70_mV), xytext=(4, 3), textcoords="offset points", fontsize=8)
    axes[2].set(xlabel="C/20 branch MAE [mV]", ylabel="Held-out dynamic MAE [mV]", title="Pareto map (* selected)")
    for ax in axes:
        ax.grid(alpha=0.25)
    fig.suptitle("Common hysteresis-scale trade-off with fixed OCP window")
    fig.tight_layout()
    fig.savefig(RESULTS / "common_scale_pareto_tradeoff.png", dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.2))
    for col, rate in enumerate(ga.RATES):
        for row_idx, charge in enumerate((True, False)):
            ax = axes[row_idx, col]
            label = "Charge" if charge else "Discharge"
            metric = dynamic[np.isclose(dynamic.scale, selected_scale) & np.isclose(dynamic.C_rate, rate) & (dynamic.Direction == label)].iloc[0]
            sim = simulations[(selected_scale, rate, charge)]
            obs = exp[(rate, charge)]
            ax.plot(obs["t_min"], obs["V"], "k", lw=2.3, label="Experiment")
            ax.plot(sim["t_min"], sim["V"], color="#9467bd", lw=2, label=f"Selected scale {selected_scale:.2f}")
            ax.text(0.04, 0.06, f"MAE 10-70% = {metric.MAE_SOC10_70_mV:.1f} mV\nCapacity error = {metric.capacity_error_Ah:+.3f} Ah", transform=ax.transAxes, fontsize=9, bbox={"facecolor": "white", "alpha": 0.82})
            ax.set(title=f"{rate:g}C {label}", xlabel="Time [min]", ylabel="Terminal voltage [V]")
            ax.grid(alpha=0.25)
            if row_idx == 0 and col == 0:
                ax.legend(fontsize=8)
    fig.suptitle(f"Selected OCP-admissible compromise: common scale {selected_scale:.2f}")
    fig.tight_layout()
    fig.savefig(RESULTS / "selected_common_scale_dynamic_curves.png", dpi=200)
    plt.close(fig)

    report = {
        "selected_scale": selected_scale,
        "selection_status": selection_status,
        "new_260918_GITT_used": False,
        "calibration": "fixed endpoint<=10 mV qOCV window and legacy C/20",
        "validation": "0.5C/1C/2C held out until ranking",
        "fixed_window": stage,
        "gates": {"qOCV_endpoint_mV": 10, "C20_branch_MAE_mV": 20, "mean_capacity_error_Ah": 0.05},
    }
    (RESULTS / "common_scale_pareto_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print("\nSummary")
    print(summary.to_string(index=False))
    print("\nSelected scale:", selected_scale)
    print(selection_status)
    print("Saved to", RESULTS)


if __name__ == "__main__":
    main()
