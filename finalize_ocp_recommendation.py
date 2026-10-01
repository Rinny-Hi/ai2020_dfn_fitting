"""Create the final OCP-first recommendation with capacity-normalized gates."""

from __future__ import annotations

import json
from pathlib import Path

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
RESULTS = ROOT / "results" / "260920_final_ocp_recommendation"
RESULTS.mkdir(parents=True, exist_ok=True)


def main() -> None:
    paths = ga.get_paths(ROOT)
    bundle = ga.load_legacy_bundle(paths)
    model = ga.build_legacy_model_inputs(bundle)
    anode = ehq.build_anode_detail(bundle, model)
    cathode = ehq.build_cathode_detail(bundle, model)
    measured_charge, measured_discharge = common.fullcell_c20_curves(bundle, model["soc"])
    sweep = pd.read_csv(ROOT / "results" / "260920_ocp_common_scale_pareto" / "common_scale_pareto_summary.csv")
    detail = pd.read_csv(ROOT / "results" / "260920_ocp_common_scale_pareto" / "common_scale_pareto_dynamic_detail.csv")

    mean_capacity_gate = 0.02 * model["q_meas"]
    max_capacity_gate = 0.03 * model["q_meas"]
    sweep["final_admissible"] = (
        (sweep.qOCV_max_endpoint_abs_mV <= 10.0 + 1e-6)
        & (sweep.C20_branch_mean_MAE_2_98_mV <= 20.0)
        & (sweep.Mean_abs_capacity_error_Ah <= mean_capacity_gate)
        & (sweep.Max_abs_capacity_error_Ah <= max_capacity_gate)
        & sweep.inside_halfcell_coverage
    )
    selected = sweep[sweep.final_admissible].sort_values(["Dynamic_MAE_SOC10_70_mV", "C20_branch_mean_MAE_2_98_mV"]).iloc[0]
    scale = float(selected.scale)
    stage = {"x0": float(selected.x0), "x100": float(selected.x100), "y100": float(selected.y100), "y0": float(selected.y0)}
    curves = common.ocp_curves(stage["x0"], stage["y100"], scale, model, anode, cathode)

    exp = omc.load_old_dynamic_data(model["dynamic"], model["q_meas"])
    simulations = {}
    final_rows = []
    for rate in ga.RATES:
        for charge in (True, False):
            direction = "Charge" if charge else "Discharge"
            print(f"Final selected scale={scale:.2f} | {rate:g}C | {direction}", flush=True)
            sim = ehq.run_dfn(
                rate,
                charge,
                stage,
                model,
                common.scaled_detail(anode, scale),
                common.scaled_detail(cathode, scale),
                True,
                True,
            )
            simulations[(rate, charge)] = sim
            metric = omc.old_time_metrics(exp[(rate, charge)], sim)
            metric["capacity_error_Ah"] = metric["end_time_error_min"] * rate * 2.28 / 60.0
            metric["capacity_error_percent_Q"] = 100 * metric["capacity_error_Ah"] / model["q_meas"]
            metric["negative_surface_sto_min"] = sim["negative_surface_sto_min"]
            metric["negative_surface_sto_max"] = sim["negative_surface_sto_max"]
            metric["positive_surface_sto_min"] = sim["positive_surface_sto_min"]
            metric["positive_surface_sto_max"] = sim["positive_surface_sto_max"]
            metric["negative_OCP_extrapolation_max"] = max(
                0.0,
                sim["negative_surface_sto_max"] - float(anode["grid"].max()),
                float(anode["grid"].min()) - sim["negative_surface_sto_min"],
            )
            metric["positive_OCP_extrapolation_max"] = max(
                0.0,
                sim["positive_surface_sto_max"] - float(cathode["grid"].max()),
                float(cathode["grid"].min()) - sim["positive_surface_sto_min"],
            )
            final_rows.append({"C_rate": rate, "Direction": direction, **metric})
    final_detail = pd.DataFrame(final_rows)
    final_detail.to_csv(RESULTS / "final_selected_dynamic_metrics.csv", index=False, encoding="utf-8-sig")
    sweep.to_csv(RESULTS / "final_candidate_decision_table.csv", index=False, encoding="utf-8-sig")

    # OCP evidence figure.
    mask = (model["soc"] >= 0.02) & (model["soc"] <= 0.98)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
    observed = [model["v_qocv"], measured_charge, measured_discharge]
    predicted = [curves["equilibrium"], curves["charge"], curves["discharge"]]
    titles = [
        f"Equilibrium qOCV\nMAE={selected.qOCV_MAE_2_98_mV:.1f} mV, endpoint max={selected.qOCV_max_endpoint_abs_mV:.1f} mV",
        f"C/20 charge OCP\nMAE={selected.C20_charge_MAE_2_98_mV:.1f} mV",
        f"C/20 discharge OCP\nMAE={selected.C20_discharge_MAE_2_98_mV:.1f} mV",
    ]
    for ax, obs, pred, title in zip(axes, observed, predicted, titles):
        ax.plot(model["soc"], obs, "k", lw=2.5, label="Measured")
        ax.plot(model["soc"], pred, color="#9467bd", lw=2.2, label=f"Selected scale {scale:.2f}")
        ax.axvspan(0.02, 0.98, color="#d9d9d9", alpha=0.15)
        ax.set(title=title, xlabel="SOC", ylabel="Cell voltage [V]")
        ax.grid(alpha=0.25)
    axes[0].legend(fontsize=8)
    fig.suptitle("Final OCP-first candidate: equilibrium and directional C/20 evidence")
    fig.tight_layout()
    fig.savefig(RESULTS / "final_selected_ocp_curves.png", dpi=200)
    plt.close(fig)

    # Held-out validation figure.
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.2))
    for col, rate in enumerate(ga.RATES):
        for row_idx, charge in enumerate((True, False)):
            ax = axes[row_idx, col]
            direction = "Charge" if charge else "Discharge"
            metric = final_detail[np.isclose(final_detail.C_rate, rate) & (final_detail.Direction == direction)].iloc[0]
            sim = simulations[(rate, charge)]
            obs = exp[(rate, charge)]
            ax.plot(obs["t_min"], obs["V"], "k", lw=2.3, label="Experiment")
            ax.plot(sim["t_min"], sim["V"], color="#9467bd", lw=2, label=f"Selected scale {scale:.2f}")
            ax.text(
                0.04,
                0.06,
                f"MAE 10-70% = {metric.MAE_SOC10_70_mV:.1f} mV\nEnd-time error = {metric.end_time_error_min:+.2f} min\nCapacity error = {metric.capacity_error_Ah:+.3f} Ah ({metric.capacity_error_percent_Q:+.1f}%)",
                transform=ax.transAxes,
                fontsize=8.5,
                bbox={"facecolor": "white", "alpha": 0.84, "edgecolor": "#cccccc"},
            )
            ax.set(title=f"{rate:g}C {direction}", xlabel="Time [min]", ylabel="Terminal voltage [V]")
            ax.grid(alpha=0.25)
            if row_idx == 0 and col == 0:
                ax.legend(fontsize=8)
    fig.suptitle("Held-out 0.5C/1C/2C validation of final OCP-first candidate")
    fig.tight_layout()
    fig.savefig(RESULTS / "final_selected_dynamic_curves.png", dpi=200)
    plt.close(fig)

    # Decision figure showing why the raw-gap and scale 0.70 candidates are rejected.
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
    axes[0].plot(sweep.scale, sweep.C20_branch_mean_MAE_2_98_mV, "o-", color="#1f77b4")
    axes[0].axhline(20, color="red", ls="--", lw=1.3, label="20 mV gate")
    axes[0].scatter([scale], [selected.C20_branch_mean_MAE_2_98_mV], marker="*", s=220, color="#9467bd", zorder=5, label="Selected")
    axes[0].set(xlabel="Common hysteresis scale", ylabel="C/20 branch MAE [mV]", title="Directional OCP fit")
    axes[1].plot(sweep.scale, sweep.Dynamic_MAE_SOC10_70_mV, "o-", color="#2ca02c")
    axes[1].scatter([scale], [selected.Dynamic_MAE_SOC10_70_mV], marker="*", s=220, color="#9467bd", zorder=5)
    axes[1].set(xlabel="Common hysteresis scale", ylabel="Held-out MAE [mV]", title="0.5C/1C/2C voltage")
    axes[2].plot(sweep.scale, 100 * sweep.Mean_abs_capacity_error_Ah / model["q_meas"], "o-", label="Mean absolute")
    axes[2].plot(sweep.scale, 100 * sweep.Max_abs_capacity_error_Ah / model["q_meas"], "o-", label="Maximum absolute")
    axes[2].axhline(2, color="#ff7f0e", ls="--", lw=1.3, label="2% mean gate")
    axes[2].axhline(3, color="red", ls="--", lw=1.3, label="3% max gate")
    axes[2].scatter([scale], [100 * selected.Max_abs_capacity_error_Ah / model["q_meas"]], marker="*", s=220, color="#9467bd", zorder=5)
    axes[2].set(xlabel="Common hysteresis scale", ylabel="Capacity error [% of Q]", title="Cutoff-capacity control")
    for ax in axes:
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    fig.suptitle("Selection logic: low voltage MAE is accepted only inside OCP and capacity gates")
    fig.tight_layout()
    fig.savefig(RESULTS / "final_selection_tradeoffs.png", dpi=200)
    plt.close(fig)

    report = {
        "recommendation": "provisional OCP-first candidate",
        "selected_common_hysteresis_scale": scale,
        "new_260918_GITT_used": False,
        "fixed_window": stage,
        "fixed_delta_x": model["delta_x"],
        "fixed_delta_y": model["delta_y"],
        "calibration_data": ["legacy half-cell OCP", "full-cell C/20 charge/discharge"],
        "held_out_validation": ["0.5C", "1C", "2C", "charge", "discharge"],
        "gates": {
            "qOCV_endpoint_max_mV": 10.0,
            "C20_branch_mean_MAE_mV": 20.0,
            "mean_abs_capacity_error_percent_Q": 2.0,
            "max_abs_capacity_error_percent_Q": 3.0,
        },
        "selected_metrics": {
            "qOCV_MAE_2_98_mV": float(selected.qOCV_MAE_2_98_mV),
            "qOCV_max_endpoint_abs_mV": float(selected.qOCV_max_endpoint_abs_mV),
            "C20_branch_mean_MAE_2_98_mV": float(selected.C20_branch_mean_MAE_2_98_mV),
            "Dynamic_MAE_SOC10_70_mV": float(selected.Dynamic_MAE_SOC10_70_mV),
            "Mean_abs_capacity_error_Ah": float(selected.Mean_abs_capacity_error_Ah),
            "Max_abs_capacity_error_Ah": float(selected.Max_abs_capacity_error_Ah),
            "Max_negative_OCP_stoichiometry_extrapolation": float(final_detail.negative_OCP_extrapolation_max.max()),
            "Max_positive_OCP_stoichiometry_extrapolation": float(final_detail.positive_OCP_extrapolation_max.max()),
        },
        "limitation": "Held-out discharge MAE remains high; do not treat as a final validated DFN parameter set.",
    }
    (RESULTS / "final_recommendation_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print("\nSelected scale", scale)
    print(pd.DataFrame([selected]).to_string(index=False))
    print("\nPer-condition metrics")
    print(final_detail.to_string(index=False))
    print("Saved to", RESULTS)


if __name__ == "__main__":
    main()
