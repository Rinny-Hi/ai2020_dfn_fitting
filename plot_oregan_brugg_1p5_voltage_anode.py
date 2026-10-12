"""Compare ORegan2022 b_n=1.5 against the current b_n=2.5 result."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import comprehensive_prefit_cross_cohort as cross
import gitt_ocp_analysis as ga
import priority_initial_state_protocol_recheck as priority

import analyze_oregan_brugg_2p5_anode_potential as anode_analysis
import electrolyte_transport_set_comparison as comparison


BRUGG_N = 1.5


def main() -> None:
    comparison.OUT.mkdir(parents=True, exist_ok=True)
    baseline, stage, anode, cathode, experiment, qcell = comparison.current_p0()
    model = comparison.with_transport_set(baseline, "ORegan2022")
    model["params"].update(
        {"Negative electrode Bruggeman coefficient (electrolyte)": BRUGG_N}
    )

    previous_voltage = pd.read_csv(
        comparison.OUT / "oregan2022_brugg_n_2p5_curve_data.csv"
    )
    previous_anode = pd.read_csv(
        comparison.OUT / "oregan2022_brugg_n_2p5_anode_potential_curves.csv"
    )

    voltage_figure, voltage_axes = plt.subplots(
        2,
        3,
        figsize=(15.6, 8.5),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    anode_figure, anode_axes = plt.subplots(
        2,
        3,
        figsize=(15.6, 8.5),
        sharex=True,
        sharey="row",
        constrained_layout=True,
    )
    summary_rows = []
    curve_rows = []

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
            q_exp = priority.transferred_capacity(exp, qcell, charge)
            metrics = priority.curve_metrics(exp, output, qcell, charge)

            potential = output["potential_V"]
            minimum_index = int(np.nanargmin(potential))
            below_zero_fraction = anode_analysis.time_fraction(
                output["time_s"], potential < 0.0
            )
            summary_rows.append(
                {
                    "C_rate": rate,
                    "Direction": direction,
                    "Initial_SOC_pct": 100.0 * z0,
                    **metrics,
                    "Capacity_error_mAh": 1000.0
                    * (metrics["Q_end_model_Ah"] - metrics["Q_end_exp_Ah"]),
                    "Minimum_anode_potential_V": float(potential[minimum_index]),
                    "Surface_stoichiometry_at_minimum": float(
                        output["surface_sto"][minimum_index]
                    ),
                    "Below_0V_time_fraction_pct": 100.0 * below_zero_fraction,
                }
            )

            prior_v = previous_voltage[
                np.isclose(previous_voltage["C_rate"], rate)
                & (previous_voltage["Direction"] == direction)
                & (previous_voltage["Source"] == "DFN")
            ]
            vax = voltage_axes[row, column]
            vax.plot(
                q_exp,
                exp["V"],
                color="#202020",
                linewidth=2.3,
                label="Experiment",
            )
            vax.plot(
                prior_v["Transferred_capacity_Ah"],
                prior_v["Voltage_V"],
                color="#D55E00",
                linewidth=1.8,
                linestyle="-.",
                label="$b_n=2.5$",
            )
            vax.plot(
                output["capacity_Ah"],
                output["V"],
                color="#0072B2",
                linewidth=2.0,
                linestyle="--",
                label="$b_n=1.5$",
            )
            vax.set_title("{}C {}".format(rate, direction), fontsize=12, weight="bold")
            vax.set_xlim(0.0, 2.42)
            vax.set_ylim(2.9, 4.3)
            vax.grid(True, color="#B0B0B0", alpha=0.25, linewidth=0.7)
            vax.text(
                0.03,
                0.04,
                (
                    "$b_n=1.5$\n"
                    "Full RMSE {:.1f} mV\n"
                    "Center RMSE {:.1f} mV\n"
                    "capacity error {:+.2f}%"
                ).format(
                    metrics["Full_RMSE_mV"],
                    metrics["Center10_70_RMSE_mV"],
                    metrics["Capacity_error_pct"],
                ),
                transform=vax.transAxes,
                fontsize=9.2,
                va="bottom",
                bbox={
                    "boxstyle": "round,pad=0.3",
                    "facecolor": "white",
                    "edgecolor": "#B0B0B0",
                    "alpha": 0.86,
                },
            )
            if row == 1:
                vax.set_xlabel("Transferred capacity (Ah)")
            if column == 0:
                vax.set_ylabel("Cell voltage (V)")

            prior_a = previous_anode[
                np.isclose(previous_anode["C_rate"], rate)
                & (previous_anode["Direction"] == direction)
            ]
            aax = anode_axes[row, column]
            aax.plot(
                prior_a["Transferred_capacity_Ah"],
                prior_a["Anode_potential_V"],
                color="#D55E00",
                linewidth=1.9,
                linestyle="-.",
                label="$b_n=2.5$ anode potential",
            )
            aax.plot(
                output["capacity_Ah"],
                potential,
                color="#0072B2",
                linewidth=2.2,
                label="$b_n=1.5$ anode potential",
            )
            aax.plot(
                output["capacity_Ah"],
                output["ocp_V"],
                color="#555555",
                linewidth=1.6,
                linestyle="--",
                label="$b_n=1.5$ local OCP",
            )
            aax.axhline(0.0, color="#202020", linewidth=1.1, linestyle=":")
            aax.scatter(
                [output["capacity_Ah"][minimum_index]],
                [potential[minimum_index]],
                color="#0072B2",
                s=28,
                zorder=4,
            )
            aax.set_title("{}C {}".format(rate, direction), fontsize=12, weight="bold")
            aax.set_xlim(0.0, 2.42)
            aax.set_ylim((-0.12, 0.35) if charge else (-0.03, 0.85))
            aax.grid(True, color="#B0B0B0", alpha=0.25, linewidth=0.7)
            aax.text(
                0.03,
                0.05,
                "$b_n=1.5$: min {:+.1f} mV\n< 0 V: {:.1f}% of time".format(
                    1000.0 * potential[minimum_index],
                    100.0 * below_zero_fraction,
                ),
                transform=aax.transAxes,
                fontsize=9.2,
                va="bottom",
                bbox={
                    "boxstyle": "round,pad=0.3",
                    "facecolor": "white",
                    "edgecolor": "#B0B0B0",
                    "alpha": 0.86,
                },
            )
            if row == 1:
                aax.set_xlabel("Transferred capacity (Ah)")
            if column == 0:
                aax.set_ylabel("Anode potential (V vs. Li/Li+)")

            for index in range(len(output["capacity_Ah"])):
                curve_rows.append(
                    {
                        "C_rate": rate,
                        "Direction": direction,
                        "Transferred_capacity_Ah": float(output["capacity_Ah"][index]),
                        "Terminal_voltage_V": float(output["V"][index]),
                        "Anode_potential_V": float(potential[index]),
                        "Local_equilibrium_OCP_V": float(output["ocp_V"][index]),
                        "Negative_surface_stoichiometry": float(
                            output["surface_sto"][index]
                        ),
                    }
                )

    voltage_handles, voltage_labels = voltage_axes[0, 0].get_legend_handles_labels()
    voltage_figure.legend(
        voltage_handles,
        voltage_labels,
        loc="upper center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 1.02),
    )
    voltage_figure.suptitle(
        "Cell 6-6 terminal-voltage comparison\n"
        "ORegan2022: $b_n=1.5$ versus current $b_n=2.5$",
        fontsize=15,
        weight="bold",
        y=1.08,
    )
    anode_handles, anode_labels = anode_axes[0, 0].get_legend_handles_labels()
    anode_figure.legend(
        anode_handles,
        anode_labels,
        loc="upper center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 1.02),
    )
    anode_figure.suptitle(
        "Cell 6-6 anode-potential comparison\n"
        "ORegan2022: $b_n=1.5$ versus current $b_n=2.5$",
        fontsize=15,
        weight="bold",
        y=1.08,
    )

    voltage_path = comparison.OUT / "oregan2022_brugg_n_1p5_voltage_comparison.png"
    anode_path = comparison.OUT / "oregan2022_brugg_n_1p5_anode_comparison.png"
    summary_path = comparison.OUT / "oregan2022_brugg_n_1p5_summary.csv"
    curves_path = comparison.OUT / "oregan2022_brugg_n_1p5_curves.csv"
    voltage_figure.savefig(voltage_path, dpi=240, bbox_inches="tight")
    anode_figure.savefig(anode_path, dpi=240, bbox_inches="tight")
    plt.close(voltage_figure)
    plt.close(anode_figure)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(summary_path, index=False, encoding="utf-8-sig")
    pd.DataFrame(curve_rows).to_csv(curves_path, index=False, encoding="utf-8-sig")
    print(voltage_path)
    print(anode_path)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
