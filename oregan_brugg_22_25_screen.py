"""Fine screen of ORegan2022 negative-electrode Bruggeman coefficient."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import comprehensive_prefit_cross_cohort as cross
import gitt_ocp_analysis as ga
import prefit_physical_anchor_decision as p0
import priority_initial_state_protocol_recheck as priority

import analyze_oregan_brugg_2p5_anode_potential as anode_analysis
import electrolyte_transport_set_comparison as comparison


BRUGG_VALUES = (2.20, 2.25, 2.30, 2.35, 2.40, 2.45, 2.50)
OUT = comparison.ROOT / "results" / "261012_oregan_brugg_22_25_screen"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    baseline, stage, anode, cathode, experiment, qcell = comparison.current_p0()
    detail_rows = []

    for brugg_n in BRUGG_VALUES:
        model = comparison.with_transport_set(baseline, "ORegan2022")
        model["params"].update(
            {"Negative electrode Bruggeman coefficient (electrolyte)": brugg_n}
        )
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                rest_v, measured_current = comparison.PROTOCOL[(rate, charge)]
                z0 = priority.monotone_voltage_inverse(
                    model["v_qocv"], model["soc"], rest_v
                )
                local_stage = priority.stage_at_soc(stage, z0, charge)
                print(
                    "b_n={:.2f} | {}C {}".format(brugg_n, rate, direction),
                    flush=True,
                )
                output = anode_analysis.run_branch(
                    measured_current / cross.NOMINAL_CAPACITY_AH,
                    charge,
                    local_stage,
                    model,
                    anode,
                    cathode,
                )
                metrics = priority.curve_metrics(
                    experiment[(rate, charge)], output, qcell, charge
                )
                potential = output["potential_V"]
                minimum_index = int(np.nanargmin(potential))
                below_zero_fraction = anode_analysis.time_fraction(
                    output["time_s"], potential < 0.0
                )
                detail_rows.append(
                    {
                        "Scenario": "ORegan2022 b_n={:.2f}".format(brugg_n),
                        "Trust": "diagnostic fine screen",
                        "Note": "ORegan2022; only brugg_n varied; no fitting",
                        "brugg_n": brugg_n,
                        "Cell": "6-6",
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

    detail = pd.DataFrame(detail_rows)
    summary = p0.summarize(detail)
    scenario_to_brugg = detail[["Scenario", "brugg_n"]].drop_duplicates()
    summary = summary.merge(scenario_to_brugg, on="Scenario", how="left")
    charge = detail[detail["Direction"] == "Charge"].copy()
    charge_risk = (
        charge.groupby("brugg_n", as_index=False)
        .agg(
            Worst_charge_anode_potential_mV=(
                "Minimum_anode_potential_V",
                lambda x: 1000.0 * float(np.min(x)),
            ),
            Mean_charge_below_0V_time_pct=("Below_0V_time_fraction_pct", "mean"),
            Worst_charge_below_0V_time_pct=("Below_0V_time_fraction_pct", "max"),
            Max_charge_surface_stoichiometry=(
                "Surface_stoichiometry_at_minimum",
                "max",
            ),
        )
    )
    summary = summary.merge(charge_risk, on="brugg_n", how="left")
    detail.to_csv(OUT / "brugg_22_25_detail.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(OUT / "brugg_22_25_summary.csv", index=False, encoding="utf-8-sig")

    plot_summary = summary.sort_values("brugg_n")
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 8.6), constrained_layout=True)
    axes[0, 0].plot(plot_summary["brugg_n"], plot_summary["Full_RMSE_mV"], "o-", label="Full")
    axes[0, 0].plot(
        plot_summary["brugg_n"], plot_summary["Center10_70_RMSE_mV"], "s--", label="Center 10-70%"
    )
    axes[0, 0].set(ylabel="Voltage RMSE (mV)", title="Voltage fit")
    axes[0, 0].legend()
    axes[0, 1].plot(plot_summary["brugg_n"], plot_summary["Capacity_RMSE_pct"], "o-")
    axes[0, 1].set(ylabel="Capacity RMSE (%)", title="Cutoff capacity")
    axes[1, 0].plot(
        plot_summary["brugg_n"], plot_summary["Worst_charge_anode_potential_mV"], "o-"
    )
    axes[1, 0].axhline(0.0, color="black", linestyle=":", linewidth=1.0)
    axes[1, 0].set(
        xlabel="$b_n$",
        ylabel="Worst charge anode potential (mV)",
        title="Anode-potential risk",
    )
    axes[1, 1].plot(
        plot_summary["brugg_n"], plot_summary["Worst_charge_below_0V_time_pct"], "o-"
    )
    axes[1, 1].set(
        xlabel="$b_n$",
        ylabel="Worst charge time below 0 V (%)",
        title="Sub-zero duration",
    )
    for ax in axes.flat:
        ax.grid(alpha=0.25)
    fig.suptitle("ORegan2022 fine screen: negative-electrode Bruggeman coefficient")
    fig.savefig(OUT / "brugg_22_25_tradeoff.png", dpi=240)
    plt.close(fig)
    print("\n", summary.to_string(index=False))


if __name__ == "__main__":
    main()
