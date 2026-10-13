"""Plot the current no-fit ORegan2022/b_n=2.5 result against elapsed time.

The six branches are the Cell 6-6 July-BOL constant-current measurements.
Charge and discharge are initialized independently from each measured rest
voltage, then simulated to the measured protocol cutoffs (4.2/3.0 V).
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


OUT = Path(__file__).resolve().parent / "results" / "261012_oregan_bn25_time_voltage"


def time_voltage_metrics(
    experiment: dict,
    simulation: dict,
    qcell_ah: float,
    measured_current_a: float,
) -> dict[str, float]:
    """Return voltage errors on a common elapsed-time grid.

    The central interval is the same physical interval used in the earlier
    capacity-domain audit: 10--70% of Qcell.  Under constant current this maps
    uniquely to elapsed time.
    """

    te = np.asarray(experiment["t_min"], dtype=float)
    tm = np.asarray(simulation["time_s"], dtype=float) / 60.0
    te = te - te[0]
    tm = tm - tm[0]
    ve = np.asarray(experiment["V"], dtype=float)
    vm = np.asarray(simulation["V"], dtype=float)

    upper = 0.995 * min(float(te[-1]), float(tm[-1]))
    full_grid = np.linspace(0.0, upper, 400)
    center_lo = 0.10 * qcell_ah / measured_current_a * 60.0
    center_hi = min(0.70 * qcell_ah / measured_current_a * 60.0, upper)
    center_grid = np.linspace(center_lo, center_hi, 160)

    def error(grid: np.ndarray) -> np.ndarray:
        return (np.interp(grid, tm, vm) - np.interp(grid, te, ve)) * 1000.0

    full_error = error(full_grid)
    center_error = error(center_grid)
    return {
        "Full_time_RMSE_mV": float(np.sqrt(np.mean(full_error**2))),
        "Full_time_MAE_mV": float(np.mean(np.abs(full_error))),
        "Center10_70_Qcell_time_RMSE_mV": float(
            np.sqrt(np.mean(center_error**2))
        ),
        "Center10_70_Qcell_time_MAE_mV": float(np.mean(np.abs(center_error))),
        "Experiment_end_time_min": float(te[-1]),
        "DFN_end_time_min": float(tm[-1]),
        "End_time_error_min": float(tm[-1] - te[-1]),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    baseline, stage, anode, cathode, experiment, qcell = comparison.current_p0()
    model = comparison.with_transport_set(baseline, "ORegan2022")
    model["params"].update(
        {"Negative electrode Bruggeman coefficient (electrolyte)": 2.5}
    )

    fig, axes = plt.subplots(
        2,
        3,
        figsize=(15.6, 8.5),
        sharey=True,
        constrained_layout=True,
    )
    voltage_rows: list[dict] = []
    metric_rows: list[dict] = []
    anode_rows: list[dict] = []

    for column, rate in enumerate(ga.RATES):
        for row, charge in enumerate((True, False)):
            direction = "Charge" if charge else "Discharge"
            rest_v, measured_current = comparison.PROTOCOL[(rate, charge)]
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
            exp = experiment[(rate, charge)]
            exp_time = np.asarray(exp["t_min"], dtype=float)
            exp_time = exp_time - exp_time[0]
            sim_time = np.asarray(output["time_s"], dtype=float) / 60.0
            sim_time = sim_time - sim_time[0]

            voltage_metrics = time_voltage_metrics(
                exp, output, qcell, measured_current
            )
            capacity_metrics = priority.curve_metrics(
                exp, output, qcell, charge
            )
            metric_rows.append(
                {
                    "C_rate": rate,
                    "Direction": direction,
                    "Measured_current_A": measured_current,
                    "Initial_rest_voltage_V": rest_v,
                    "Initial_SOC_pct": 100.0 * z0,
                    **voltage_metrics,
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
            below_zero_fraction = anode_analysis.time_fraction(
                np.asarray(output["time_s"], dtype=float), potential < 0.0
            )
            anode_rows.append(
                {
                    "C_rate": rate,
                    "Direction": direction,
                    "Initial_SOC_pct": 100.0 * z0,
                    "Minimum_anode_potential_mV": 1000.0
                    * float(potential[minimum_index]),
                    "Time_at_minimum_min": float(sim_time[minimum_index]),
                    "Transferred_capacity_at_minimum_Ah": float(
                        output["capacity_Ah"][minimum_index]
                    ),
                    "Surface_stoichiometry_at_minimum": float(
                        output["surface_sto"][minimum_index]
                    ),
                    "Time_below_0V_pct": 100.0 * below_zero_fraction,
                }
            )

            ax = axes[row, column]
            ax.plot(
                exp_time,
                np.asarray(exp["V"], dtype=float),
                color="#202020",
                linewidth=2.3,
                label="Experiment",
            )
            ax.plot(
                sim_time,
                np.asarray(output["V"], dtype=float),
                color="#D55E00",
                linewidth=2.0,
                linestyle="--",
                label="DFN: ORegan2022, $b_n=2.5$",
            )
            ax.set_title(f"{rate:g}C {direction}", fontsize=12, weight="bold")
            ax.set_ylim(2.9, 4.3)
            ax.grid(True, color="#B0B0B0", alpha=0.25, linewidth=0.7)
            ax.text(
                0.03,
                0.04,
                (
                    f"Time RMSE {voltage_metrics['Full_time_RMSE_mV']:.1f} mV\n"
                    f"end-time error {voltage_metrics['End_time_error_min']:+.2f} min\n"
                    f"capacity error {capacity_metrics['Capacity_error_pct']:+.2f}%"
                ),
                transform=ax.transAxes,
                fontsize=9.2,
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
                ax.set_xlabel("Elapsed time (min)")
            if column == 0:
                ax.set_ylabel("Cell voltage (V)")

            for source, time_min, voltage in (
                ("Experiment", exp_time, np.asarray(exp["V"], dtype=float)),
                ("DFN", sim_time, np.asarray(output["V"], dtype=float)),
            ):
                for time_value, voltage_value in zip(time_min, voltage):
                    voltage_rows.append(
                        {
                            "C_rate": rate,
                            "Direction": direction,
                            "Source": source,
                            "Elapsed_time_min": float(time_value),
                            "Voltage_V": float(voltage_value),
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
        "Cell 6-6 measured vs DFN voltage-time curves\n"
        "ORegan2022 electrolyte transport + negative-electrode $b_n=2.5$ (no fit)",
        fontsize=15,
        weight="bold",
        y=1.08,
    )

    figure_path = OUT / "oregan2022_brugg_n_2p5_time_voltage_curves.png"
    fig.savefig(figure_path, dpi=240, bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame(voltage_rows).to_csv(
        OUT / "time_voltage_curve_data.csv", index=False, encoding="utf-8-sig"
    )
    pd.DataFrame(metric_rows).to_csv(
        OUT / "time_voltage_metrics.csv", index=False, encoding="utf-8-sig"
    )
    pd.DataFrame(anode_rows).to_csv(
        OUT / "anode_potential_summary.csv", index=False, encoding="utf-8-sig"
    )
    print(figure_path)
    print("\n[time-voltage metrics]\n", pd.DataFrame(metric_rows).to_string(index=False))
    print("\n[anode potential]\n", pd.DataFrame(anode_rows).to_string(index=False))


if __name__ == "__main__":
    main()
