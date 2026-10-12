"""Plot measured and no-fit DFN curves for ORegan2022 with b_n = 2.5."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from typing import Dict, List, Union

import comprehensive_prefit_cross_cohort as cross
import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import priority_initial_state_protocol_recheck as priority

import electrolyte_transport_set_comparison as comparison


SCENARIO = "ORegan2022, b_n = 2.5 (no fit)"


def main() -> None:
    comparison.OUT.mkdir(parents=True, exist_ok=True)
    baseline, stage, anode, cathode, experiment, qcell = comparison.current_p0()
    model = comparison.with_transport_set(baseline, "ORegan2022")
    model["params"].update(
        {"Negative electrode Bruggeman coefficient (electrolyte)": 2.5}
    )

    fig, axes = plt.subplots(
        2,
        3,
        figsize=(15.6, 8.5),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    curve_rows: List[Dict[str, Union[float, str]]] = []
    metric_rows: List[Dict[str, Union[float, str]]] = []

    for column, rate in enumerate(ga.RATES):
        for row, charge in enumerate((True, False)):
            direction = "Charge" if charge else "Discharge"
            rest_v, measured_current = comparison.PROTOCOL[(rate, charge)]
            z0 = priority.monotone_voltage_inverse(
                model["v_qocv"], model["soc"], rest_v
            )
            local_stage = priority.stage_at_soc(stage, z0, charge)
            simulation = ehq.run_dfn(
                measured_current / cross.NOMINAL_CAPACITY_AH,
                charge,
                local_stage,
                model,
                anode,
                cathode,
                False,
                False,
            )
            exp = experiment[(rate, charge)]
            q_exp = priority.transferred_capacity(exp, qcell, charge)
            q_sim = priority.transferred_capacity(simulation, qcell, charge)
            metrics = priority.curve_metrics(exp, simulation, qcell, charge)
            metric_rows.append(
                {
                    "C_rate": rate,
                    "Direction": direction,
                    "Initial_SOC_pct": 100.0 * z0,
                    **metrics,
                    "Capacity_error_mAh": 1000.0
                    * (metrics["Q_end_model_Ah"] - metrics["Q_end_exp_Ah"]),
                }
            )

            ax = axes[row, column]
            ax.plot(
                q_exp,
                exp["V"],
                color="#202020",
                linewidth=2.3,
                label="Experiment",
                zorder=2,
            )
            ax.plot(
                q_sim,
                simulation["V"],
                color="#D55E00",
                linewidth=2.0,
                linestyle="--",
                label="DFN: ORegan2022, $b_n=2.5$",
                zorder=3,
            )
            ax.set_title(f"{rate:g}C {direction}", fontsize=12, weight="bold")
            ax.set_xlim(0.0, 2.42)
            ax.set_ylim(2.9, 4.3)
            ax.grid(True, color="#B0B0B0", alpha=0.25, linewidth=0.7)
            ax.text(
                0.03,
                0.04,
                (
                    f"Full RMSE {metrics['Full_RMSE_mV']:.1f} mV\n"
                    f"Center 10-70% RMSE {metrics['Center10_70_RMSE_mV']:.1f} mV\n"
                    f"capacity error {metrics['Capacity_error_pct']:+.2f}%"
                ),
                transform=ax.transAxes,
                fontsize=9.5,
                va="bottom",
                ha="left",
                bbox={
                    "boxstyle": "round,pad=0.3",
                    "facecolor": "white",
                    "edgecolor": "#B0B0B0",
                    "alpha": 0.86,
                },
            )
            if row == 1:
                ax.set_xlabel("Transferred capacity (Ah)")
            if column == 0:
                ax.set_ylabel("Cell voltage (V)")

            for source, capacity, voltage in (
                ("Experiment", q_exp, np.asarray(exp["V"], dtype=float)),
                ("DFN", q_sim, np.asarray(simulation["V"], dtype=float)),
            ):
                for q_value, v_value in zip(capacity, voltage):
                    curve_rows.append(
                        {
                            "C_rate": rate,
                            "Direction": direction,
                            "Source": source,
                            "Transferred_capacity_Ah": float(q_value),
                            "Voltage_V": float(v_value),
                        }
                    )

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=2,
        frameon=False,
        bbox_to_anchor=(0.5, 1.02),
    )
    fig.suptitle(
        "Cell 6-6: measured vs DFN voltage curves\n"
        "ORegan2022 electrolyte transport + negative-electrode $b_n=2.5$",
        fontsize=15,
        weight="bold",
        y=1.08,
    )

    figure_path = comparison.OUT / "oregan2022_brugg_n_2p5_voltage_curves.png"
    curve_path = comparison.OUT / "oregan2022_brugg_n_2p5_curve_data.csv"
    metric_path = comparison.OUT / "oregan2022_brugg_n_2p5_curve_metrics.csv"
    fig.savefig(figure_path, dpi=240, bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame(curve_rows).to_csv(curve_path, index=False, encoding="utf-8-sig")
    pd.DataFrame(metric_rows).to_csv(metric_path, index=False, encoding="utf-8-sig")
    print(figure_path)
    print(pd.DataFrame(metric_rows).to_string(index=False))


if __name__ == "__main__":
    main()
