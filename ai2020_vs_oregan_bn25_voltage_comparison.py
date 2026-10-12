"""Voltage validation of Ai2020 vs ORegan2022 electrolyte at fixed b_n=2.5.

All non-electrolyte parameters are identical.  The comparison uses the same
Cell 6-6 July-BOL 0.5C/1C/2C charge/discharge branches and independent rest-
voltage initialization as the retained baseline analysis.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import comprehensive_prefit_cross_cohort as cross
import gitt_ocp_analysis as ga
import priority_initial_state_protocol_recheck as priority

import analyze_oregan_brugg_2p5_anode_potential as anode_analysis
import electrolyte_transport_set_comparison as comparison
import plot_oregan_brugg_2p5_time_voltage as time_plot


OUT = Path(__file__).resolve().parent / "results" / "261012_ai2020_vs_oregan_bn25_voltage"
SETS = ("ORegan2022", "Ai2020 corrected")
COLORS = {"ORegan2022": "#D55E00", "Ai2020 corrected": "#0072B2"}
STYLES = {"ORegan2022": "--", "Ai2020 corrected": ":"}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    baseline, stage, anode, cathode, experiment, qcell = comparison.current_p0()
    metric_rows: list[dict] = []
    curve_rows: list[dict] = []
    runs: dict[tuple, dict] = {}

    for set_name in SETS:
        model = comparison.with_transport_set(baseline, set_name)
        model["params"].update(
            {"Negative electrode Bruggeman coefficient (electrolyte)": 2.5},
            check_already_exists=False,
        )
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                rest_v, measured_current = comparison.PROTOCOL[(rate, charge)]
                z0 = priority.monotone_voltage_inverse(model["v_qocv"], model["soc"], rest_v)
                local_stage = priority.stage_at_soc(stage, z0, charge)
                print(f"{set_name}: {rate:g}C {direction}", flush=True)
                output = anode_analysis.run_branch(
                    measured_current / cross.NOMINAL_CAPACITY_AH,
                    charge,
                    local_stage,
                    model,
                    anode,
                    cathode,
                )
                runs[(set_name, rate, charge)] = output
                time_metrics = time_plot.time_voltage_metrics(
                    experiment[(rate, charge)], output, qcell, measured_current
                )
                capacity_metrics = priority.curve_metrics(
                    experiment[(rate, charge)], output, qcell, charge
                )
                metric_rows.append(
                    {
                        "Electrolyte_set": set_name,
                        "C_rate": rate,
                        "Direction": direction,
                        "Initial_SOC_pct": 100.0 * z0,
                        **time_metrics,
                        "Capacity_domain_full_RMSE_mV": capacity_metrics["Full_RMSE_mV"],
                        "Center10_70_RMSE_mV": capacity_metrics["Center10_70_RMSE_mV"],
                        "Capacity_error_pct": capacity_metrics["Capacity_error_pct"],
                        "Capacity_error_mAh": 1000.0
                        * (capacity_metrics["Q_end_model_Ah"] - capacity_metrics["Q_end_exp_Ah"]),
                    }
                )

    metrics = pd.DataFrame(metric_rows)
    summary_rows = []
    for set_name in SETS:
        subset = metrics[metrics.Electrolyte_set == set_name]
        summary_rows.append(
            {
                "Electrolyte_set": set_name,
                "Overall_time_voltage_RMSE_mV": float(
                    np.sqrt(np.mean(subset.Full_time_RMSE_mV.to_numpy(float) ** 2))
                ),
                "Overall_center10_70_RMSE_mV": float(
                    np.sqrt(np.mean(subset.Center10_70_RMSE_mV.to_numpy(float) ** 2))
                ),
                "Capacity_RMSE_pct": float(
                    np.sqrt(np.mean(subset.Capacity_error_pct.to_numpy(float) ** 2))
                ),
                "Maximum_abs_capacity_error_mAh": float(
                    np.max(np.abs(subset.Capacity_error_mAh.to_numpy(float)))
                ),
            }
        )
    summary = pd.DataFrame(summary_rows)
    metrics.to_csv(OUT / "branch_metrics.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(OUT / "overall_summary.csv", index=False, encoding="utf-8-sig")

    fig, axes = plt.subplots(2, 3, figsize=(15.8, 8.6), sharey=True, constrained_layout=True)
    for column, rate in enumerate(ga.RATES):
        for row, charge in enumerate((True, False)):
            direction = "Charge" if charge else "Discharge"
            ax = axes[row, column]
            exp = experiment[(rate, charge)]
            exp_t = np.array(exp["t_min"], dtype=float, copy=True)
            exp_t -= exp_t[0]
            ax.plot(exp_t, exp["V"], color="#202020", lw=2.5, label="Experiment")
            annotation = []
            for set_name in SETS:
                output = runs[(set_name, rate, charge)]
                time_min = np.asarray(output["time_s"], dtype=float) / 60.0
                time_min -= time_min[0]
                ax.plot(
                    time_min,
                    output["V"],
                    color=COLORS[set_name],
                    ls=STYLES[set_name],
                    lw=2.1,
                    label=set_name,
                )
                row_metric = metrics[
                    (metrics.Electrolyte_set == set_name)
                    & (metrics.C_rate == rate)
                    & (metrics.Direction == direction)
                ].iloc[0]
                short = "ORegan" if set_name == "ORegan2022" else "Ai2020"
                annotation.append(
                    f"{short}: {row_metric.Full_time_RMSE_mV:.1f} mV"
                )
                for i in range(len(time_min)):
                    curve_rows.append(
                        {
                            "Electrolyte_set": set_name,
                            "C_rate": rate,
                            "Direction": direction,
                            "Elapsed_time_min": float(time_min[i]),
                            "Voltage_V": float(output["V"][i]),
                        }
                    )
            ax.text(
                0.03,
                0.04,
                "\n".join(annotation),
                transform=ax.transAxes,
                fontsize=9.2,
                va="bottom",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.86),
            )
            ax.set_title(f"{rate:g}C {direction}", fontsize=12, weight="bold")
            ax.set_ylim(2.9, 4.3)
            ax.grid(alpha=0.24)
            if row == 1:
                ax.set_xlabel("Elapsed time (min)")
            if column == 0:
                ax.set_ylabel("Cell voltage (V)")

    pd.DataFrame(curve_rows).to_csv(
        OUT / "model_voltage_curves.csv", index=False, encoding="utf-8-sig"
    )
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.suptitle(
        "Cell 6-6 voltage validation: electrolyte transport only\n"
        "ORegan2022 vs Ai2020 corrected; all other parameters fixed, $b_n=2.5$",
        fontsize=15,
        weight="bold",
        y=1.08,
    )
    figure = OUT / "ai2020_vs_oregan_bn25_voltage_curves.png"
    fig.savefig(figure, dpi=240, bbox_inches="tight")
    plt.close(fig)
    print("\nBRANCH METRICS")
    print(metrics.to_string(index=False))
    print("\nOVERALL")
    print(summary.to_string(index=False))
    print(figure)


if __name__ == "__main__":
    main()
