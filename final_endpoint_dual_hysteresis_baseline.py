"""Final legacy-data baseline: endpoint window and dual-electrode hysteresis."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybamm

import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import old_notebook_method_comparison as omc


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260920_final_endpoint_dual_hysteresis_baseline"
RESULTS.mkdir(parents=True, exist_ok=True)


def soc_residual(
    experiment: dict[str, np.ndarray], simulation: dict[str, np.ndarray]
) -> tuple[np.ndarray, np.ndarray]:
    exp_soc, exp_voltage = ga._clean_curve(experiment["SOC"], experiment["V"])
    sim_soc, sim_voltage = ga._clean_curve(simulation["SOC"], simulation["V"])
    overlap = (exp_soc >= sim_soc.min()) & (exp_soc <= sim_soc.max())
    residual = (
        np.interp(exp_soc[overlap], sim_soc, sim_voltage) - exp_voltage[overlap]
    ) * 1000.0
    return exp_soc[overlap], residual


def main() -> None:
    paths = ga.get_paths(ROOT)
    bundle = ga.load_legacy_bundle(paths)
    model = ga.build_legacy_model_inputs(bundle)
    anode = ehq.build_anode_detail(bundle, model)
    cathode = ehq.build_cathode_detail(bundle, model)

    model["y_exp"] = cathode["grid"]
    model["up_exp"] = cathode["equilibrium"]
    anode_for_window = {
        "x": anode["grid"],
        "equilibrium": anode["equilibrium"],
    }
    stage, qocv_model = omc.endpoint_window(model, anode_for_window)
    exp = omc.load_old_dynamic_data(model["dynamic"], model["q_meas"])

    simulations = {}
    rows = []
    for rate in ga.RATES:
        for charge in (True, False):
            direction = "Charge" if charge else "Discharge"
            print(f"Final baseline | {rate:g}C | {direction}", flush=True)
            sim = ehq.run_dfn(
                rate,
                charge,
                stage,
                model,
                anode,
                cathode,
                negative_hysteresis=True,
                positive_hysteresis=True,
            )
            simulations[(rate, charge)] = sim
            rows.append(
                {
                    "C_rate": rate,
                    "Direction": direction,
                    **omc.old_time_metrics(exp[(rate, charge)], sim),
                }
            )

    metrics = pd.DataFrame(rows)
    metrics.to_csv(RESULTS / "final_dynamic_metrics.csv", index=False, encoding="utf-8-sig")
    summary = pd.DataFrame(
        [
            {
                "Mean_MAE_all_mV": metrics["MAE_all_mV"].mean(),
                "Mean_MAE_SOC10_70_mV": metrics["MAE_SOC10_70_mV"].mean(),
                "Mean_abs_end_time_error_min": metrics["end_time_error_min"].abs().mean(),
                "qOCV_RMSE_SOC2_98_mV": stage["qOCV_RMSE_2_98_mV"],
                "qOCV_MAE_SOC2_98_mV": stage["qOCV_MAE_2_98_mV"],
                "qOCV_SOC0_error_mV": stage["SOC_0_endpoint_error_mV"],
                "qOCV_SOC100_error_mV": stage["SOC_100_endpoint_error_mV"],
            }
        ]
    )
    summary.to_csv(RESULTS / "final_baseline_summary.csv", index=False, encoding="utf-8-sig")

    fig, ax = plt.subplots(figsize=(8.8, 5.2))
    ax.plot(model["soc"], model["v_qocv"], color="black", lw=2.6, label="Measured C/20 qOCV")
    ax.plot(model["soc"], qocv_model, color="#2ca02c", lw=2.2, label="Endpoint-window equilibrium OCP")
    ax.scatter([0, 1], [qocv_model[0], qocv_model[-1]], color="#d62728", s=45, zorder=4, label="Constrained endpoints")
    ax.text(0.03, 0.05, f"SOC 2–98% MAE = {stage['qOCV_MAE_2_98_mV']:.1f} mV\nEndpoint errors ≈ 0 mV", transform=ax.transAxes, fontsize=10, bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "#cccccc"})
    ax.set(xlabel="SOC", ylabel="Cell OCV [V]", title="Final baseline qOCV: endpoint-constrained window")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(RESULTS / "final_qocv_endpoint_window.png", dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.2))
    for col, rate in enumerate(ga.RATES):
        for row, charge in enumerate((True, False)):
            axis = axes[row, col]
            direction = "Charge" if charge else "Discharge"
            case = metrics[
                np.isclose(metrics["C_rate"], rate)
                & (metrics["Direction"] == direction)
            ].iloc[0]
            axis.plot(exp[(rate, charge)]["t_min"], exp[(rate, charge)]["V"], color="black", lw=2.3, label="Experiment")
            axis.plot(simulations[(rate, charge)]["t_min"], simulations[(rate, charge)]["V"], color="#2ca02c", lw=2.0, label="Final baseline")
            axis.text(0.04, 0.06, f"MAE 10–70% = {case['MAE_SOC10_70_mV']:.1f} mV\nFull MAE = {case['MAE_all_mV']:.1f} mV\nEnd-time error = {case['end_time_error_min']:+.1f} min", transform=axis.transAxes, fontsize=9, bbox={"facecolor": "white", "alpha": 0.82, "edgecolor": "#cccccc"})
            axis.set(title=f"{rate:g}C {direction}", xlabel="Time [min]", ylabel="Terminal voltage [V]")
            axis.grid(alpha=0.25)
            if row == 0 and col == 0:
                axis.legend(fontsize=9)
    fig.suptitle("Final DFN baseline: endpoint window + anode/cathode hysteresis", fontsize=15)
    fig.tight_layout()
    fig.savefig(RESULTS / "final_dynamic_voltage_curves.png", dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(15.5, 7.8), sharey=True)
    for col, rate in enumerate(ga.RATES):
        for row, charge in enumerate((True, False)):
            axis = axes[row, col]
            soc, residual = soc_residual(exp[(rate, charge)], simulations[(rate, charge)])
            axis.axhspan(-20, 20, color="#2ca02c", alpha=0.10, label="±20 mV")
            axis.axhline(0, color="black", lw=1)
            axis.plot(soc, residual, color="#9467bd", lw=1.8)
            axis.set(title=f"{rate:g}C {'Charge' if charge else 'Discharge'}", xlabel="Experimental SOC", ylabel="Model − experiment [mV]")
            axis.grid(alpha=0.25)
            if row == 0 and col == 0:
                axis.legend(fontsize=9)
    fig.suptitle("Final baseline voltage residuals", fontsize=15)
    fig.tight_layout()
    fig.savefig(RESULTS / "final_dynamic_voltage_residuals.png", dpi=200)
    plt.close(fig)

    manifest = {
        "new_260918_GITT_used": False,
        "window": "endpoint constrained; delta x and delta y fixed by measured C/20 capacity",
        "negative_ocp": "legacy low-rate v1/v2 directional branches with one-state hysteresis",
        "positive_ocp": "legacy GITT 7-7/7-8 rest-end directional branches with one-state hysteresis",
        "kinetics": "Ai2020 nominal",
        "pybamm_version": pybamm.__version__,
        "stage": {key: float(value) if isinstance(value, (float, np.floating)) else value for key, value in stage.items()},
    }
    (RESULTS / "final_baseline_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print("\nFinal dynamic metrics")
    print(metrics.to_string(index=False))
    print("\nSummary")
    print(summary.to_string(index=False))
    print("\nSaved to", RESULTS)


if __name__ == "__main__":
    main()
