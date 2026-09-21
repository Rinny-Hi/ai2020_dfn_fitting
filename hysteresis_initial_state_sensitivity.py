"""Sensitivity of dual-electrode hysteresis initial state to cutoff timing."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import old_notebook_method_comparison as omc


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260920_hysteresis_initial_state_sensitivity"
RESULTS.mkdir(parents=True, exist_ok=True)
ALPHAS = [-1.0, -0.5, 0.0, 0.5, 1.0]


def main() -> None:
    paths = ga.get_paths(ROOT)
    bundle = ga.load_legacy_bundle(paths)
    model = ga.build_legacy_model_inputs(bundle)
    anode = ehq.build_anode_detail(bundle, model)
    cathode = ehq.build_cathode_detail(bundle, model)
    model["y_exp"] = cathode["grid"]
    model["up_exp"] = cathode["equilibrium"]
    stage, _ = omc.endpoint_window(
        model, {"x": anode["grid"], "equilibrium": anode["equilibrium"]}
    )
    exp = omc.load_old_dynamic_data(model["dynamic"], model["q_meas"])

    rows = []
    sims = {}
    for alpha in ALPHAS:
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                target_n = -1.0 if charge else 1.0
                target_p = 1.0 if charge else -1.0
                print(f"alpha={alpha:+.1f} | {rate:g}C | {direction}", flush=True)
                sim = ehq.run_dfn(
                    rate,
                    charge,
                    stage,
                    model,
                    anode,
                    cathode,
                    negative_hysteresis=True,
                    positive_hysteresis=True,
                    negative_initial_h=alpha * target_n,
                    positive_initial_h=alpha * target_p,
                )
                sims[(alpha, rate, charge)] = sim
                metric = omc.old_time_metrics(exp[(rate, charge)], sim)
                current_a = rate * 2.28
                metric["capacity_error_Ah"] = metric["end_time_error_min"] * current_a / 60.0
                rows.append(
                    {
                        "initial_state_alpha": alpha,
                        "initial_state_interpretation": (
                            "history/opposite branch" if alpha == -1
                            else "neutral" if alpha == 0
                            else "target branch" if alpha == 1
                            else "partial"
                        ),
                        "C_rate": rate,
                        "Direction": direction,
                        **metric,
                    }
                )
    metrics = pd.DataFrame(rows)
    metrics.to_csv(RESULTS / "initial_state_sensitivity_detail.csv", index=False, encoding="utf-8-sig")
    summary = (
        metrics.groupby("initial_state_alpha", as_index=False)
        .agg(
            Mean_MAE_all_mV=("MAE_all_mV", "mean"),
            Mean_MAE_SOC10_70_mV=("MAE_SOC10_70_mV", "mean"),
            Mean_abs_end_time_error_min=("end_time_error_min", lambda x: np.mean(np.abs(x))),
            Mean_abs_capacity_error_Ah=("capacity_error_Ah", lambda x: np.mean(np.abs(x))),
        )
    )
    charge_summary = (
        metrics[metrics["Direction"] == "Charge"]
        .groupby("initial_state_alpha", as_index=False)
        .agg(
            Charge_MAE_SOC10_70_mV=("MAE_SOC10_70_mV", "mean"),
            Charge_mean_end_time_error_min=("end_time_error_min", "mean"),
            Charge_mean_capacity_error_Ah=("capacity_error_Ah", "mean"),
        )
    )
    discharge_summary = (
        metrics[metrics["Direction"] == "Discharge"]
        .groupby("initial_state_alpha", as_index=False)
        .agg(
            Discharge_MAE_SOC10_70_mV=("MAE_SOC10_70_mV", "mean"),
            Discharge_mean_end_time_error_min=("end_time_error_min", "mean"),
            Discharge_mean_capacity_error_Ah=("capacity_error_Ah", "mean"),
        )
    )
    summary = summary.merge(charge_summary, on="initial_state_alpha").merge(
        discharge_summary, on="initial_state_alpha"
    )
    summary.to_csv(RESULTS / "initial_state_sensitivity_summary.csv", index=False, encoding="utf-8-sig")

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
    axes[0].plot(summary["initial_state_alpha"], summary["Mean_MAE_SOC10_70_mV"], "o-", lw=2, label="All")
    axes[0].plot(summary["initial_state_alpha"], summary["Charge_MAE_SOC10_70_mV"], "o-", lw=2, label="Charge")
    axes[0].plot(summary["initial_state_alpha"], summary["Discharge_MAE_SOC10_70_mV"], "o-", lw=2, label="Discharge")
    axes[0].set(title="Voltage error", xlabel="Initial-state alpha", ylabel="Mean MAE 10–70% [mV]")
    axes[0].legend()
    axes[1].plot(summary["initial_state_alpha"], summary["Charge_mean_end_time_error_min"], "o-", lw=2, label="Charge")
    axes[1].plot(summary["initial_state_alpha"], summary["Discharge_mean_end_time_error_min"], "o-", lw=2, label="Discharge")
    axes[1].axhline(0, color="black", lw=1)
    axes[1].set(title="Cutoff timing", xlabel="Initial-state alpha", ylabel="Mean end-time error [min]")
    axes[1].legend()
    axes[2].plot(summary["initial_state_alpha"], summary["Charge_mean_capacity_error_Ah"], "o-", lw=2, label="Charge")
    axes[2].plot(summary["initial_state_alpha"], summary["Discharge_mean_capacity_error_Ah"], "o-", lw=2, label="Discharge")
    axes[2].axhline(0, color="black", lw=1)
    axes[2].set(title="Equivalent capacity error", xlabel="Initial-state alpha", ylabel="Mean capacity error [Ah]")
    axes[2].legend()
    for axis in axes:
        axis.grid(alpha=0.25)
        axis.set_xticks(ALPHAS)
    fig.suptitle("Dual-electrode hysteresis initial-state sensitivity\n−1: previous/history branch, 0: neutral, +1: immediate target branch", fontsize=14)
    fig.tight_layout()
    fig.savefig(RESULTS / "initial_state_sensitivity.png", dpi=200)
    plt.close(fig)

    selected_alpha = float(
        summary.loc[summary["Mean_abs_end_time_error_min"].idxmin(), "initial_state_alpha"]
    )
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.2))
    for col, rate in enumerate(ga.RATES):
        for row, charge in enumerate((True, False)):
            axis = axes[row, col]
            axis.plot(exp[(rate, charge)]["t_min"], exp[(rate, charge)]["V"], color="black", lw=2.3, label="Experiment")
            for alpha, color, label in (
                (1.0, "#2ca02c", "Immediate target branch"),
                (0.0, "#ff7f0e", "Neutral state"),
                (-1.0, "#1f77b4", "History-consistent state"),
            ):
                axis.plot(sims[(alpha, rate, charge)]["t_min"], sims[(alpha, rate, charge)]["V"], color=color, lw=1.7, label=label)
            axis.set(title=f"{rate:g}C {'Charge' if charge else 'Discharge'}", xlabel="Time [min]", ylabel="Terminal voltage [V]")
            axis.grid(alpha=0.25)
            if row == 0 and col == 0:
                axis.legend(fontsize=7.5)
    fig.suptitle("Effect of hysteresis initial state on voltage and cutoff", fontsize=15)
    fig.tight_layout()
    fig.savefig(RESULTS / "initial_state_voltage_curves.png", dpi=200)
    plt.close(fig)

    print("\nSummary")
    print(summary.to_string(index=False))
    print("\nMinimum mean absolute end-time error alpha:", selected_alpha)
    print("Saved to", RESULTS)


if __name__ == "__main__":
    main()
