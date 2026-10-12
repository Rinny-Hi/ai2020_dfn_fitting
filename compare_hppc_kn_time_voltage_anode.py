"""Compare the retained k_n anchor with the HPPC-pulse-only preferred k_n."""

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
import oregan_bn25_kn_only_fit as kn_fit
from plot_oregan_brugg_2p5_time_voltage import time_voltage_metrics


OUT = Path(__file__).resolve().parent / "results" / "261012_hppc_kn_time_voltage_anode"
HPPC_SCALE = 2.256441945848353


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    baseline, stage, anode, cathode, experiment, qcell = comparison.current_p0()
    transport = comparison.with_transport_set(baseline, "ORegan2022")
    scenarios = {
        "Retained anchor ($k_n=1.00\\times10^{-11}$)": kn_fit.apply_kn(transport, 1.0),
        "HPPC preferred ($k_n=2.256\\times10^{-11}$)": kn_fit.apply_kn(
            transport, HPPC_SCALE
        ),
    }

    fig, axes = plt.subplots(2, 3, figsize=(15.8, 8.6), sharey=True, constrained_layout=True)
    metric_rows: list[dict] = []
    anode_rows: list[dict] = []

    for column, rate in enumerate(ga.RATES):
        for row, charge in enumerate((True, False)):
            direction = "Charge" if charge else "Discharge"
            rest_v, measured_current = comparison.PROTOCOL[(rate, charge)]
            exp = experiment[(rate, charge)]
            exp_time = np.asarray(exp["t_min"], dtype=float)
            exp_time -= exp_time[0]
            ax = axes[row, column]
            ax.plot(exp_time, exp["V"], color="#202020", linewidth=2.4, label="Experiment")

            scenario_results = {}
            for name, model in scenarios.items():
                z0 = priority.monotone_voltage_inverse(
                    model["v_qocv"], model["soc"], rest_v
                )
                local_stage = priority.stage_at_soc(stage, z0, charge)
                output = anode_analysis.run_branch(
                    measured_current / cross.NOMINAL_CAPACITY_AH,
                    charge,
                    local_stage,
                    model,
                    anode,
                    cathode,
                )
                scenario_results[name] = output
                time_min = np.asarray(output["time_s"], dtype=float) / 60.0
                time_min -= time_min[0]
                style = "--" if name.startswith("Retained") else "-"
                color = "#E69F00" if name.startswith("Retained") else "#0072B2"
                ax.plot(time_min, output["V"], color=color, linestyle=style, linewidth=2.0, label=name)

                time_metrics = time_voltage_metrics(exp, output, qcell, measured_current)
                capacity_metrics = priority.curve_metrics(exp, output, qcell, charge)
                metric_rows.append(
                    {
                        "Scenario": name,
                        "C_rate": rate,
                        "Direction": direction,
                        **time_metrics,
                        "Capacity_error_pct": capacity_metrics["Capacity_error_pct"],
                        "Capacity_error_mAh": 1000.0
                        * (
                            capacity_metrics["Q_end_model_Ah"]
                            - capacity_metrics["Q_end_exp_Ah"]
                        ),
                    }
                )

                potential = np.asarray(output["potential_V"], dtype=float)
                minimum_index = int(np.nanargmin(potential))
                anode_rows.append(
                    {
                        "Scenario": name,
                        "C_rate": rate,
                        "Direction": direction,
                        "Minimum_anode_potential_mV": 1000.0
                        * float(potential[minimum_index]),
                        "Time_at_minimum_min": float(time_min[minimum_index]),
                        "Surface_stoichiometry_at_minimum": float(
                            output["surface_sto"][minimum_index]
                        ),
                        "Time_below_0V_pct": 100.0
                        * anode_analysis.time_fraction(
                            np.asarray(output["time_s"], dtype=float),
                            potential < 0.0,
                        ),
                    }
                )

            preferred_name = next(name for name in scenarios if name.startswith("HPPC"))
            preferred_metrics = next(
                item
                for item in metric_rows
                if item["Scenario"] == preferred_name
                and item["C_rate"] == rate
                and item["Direction"] == direction
            )
            ax.text(
                0.03,
                0.04,
                (
                    "HPPC-$k_n$ result\n"
                    f"time RMSE {preferred_metrics['Full_time_RMSE_mV']:.1f} mV\n"
                    f"capacity error {preferred_metrics['Capacity_error_pct']:+.2f}%"
                ),
                transform=ax.transAxes,
                fontsize=9.1,
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

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.suptitle(
        "Cell 6-6: effect of HPPC-preferred negative-electrode kinetics\n"
        "ORegan2022 electrolyte transport + $b_n=2.5$; all other parameters fixed",
        fontsize=15,
        weight="bold",
        y=1.08,
    )
    fig.savefig(OUT / "hppc_preferred_kn_time_voltage_comparison.png", dpi=240, bbox_inches="tight")
    plt.close(fig)

    metrics = pd.DataFrame(metric_rows)
    anode_summary = pd.DataFrame(anode_rows)
    metrics.to_csv(OUT / "time_voltage_comparison.csv", index=False, encoding="utf-8-sig")
    anode_summary.to_csv(OUT / "anode_potential_comparison.csv", index=False, encoding="utf-8-sig")
    print("\n[time voltage]\n", metrics.to_string(index=False))
    print("\n[anode potential]\n", anode_summary.to_string(index=False))


if __name__ == "__main__":
    main()
