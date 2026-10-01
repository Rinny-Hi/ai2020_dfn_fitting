"""Plot the current P0 curves and test whether an extra series R0 is justified."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import c50_selected_ocp_dynamic_validation as selected
import comprehensive_prefit_cross_cohort as cross
import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import prefit_physical_anchor_decision as p0
import priority_initial_state_protocol_recheck as priority


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260927_current_p0_curve_and_r0_diagnostic"
RESULTS.mkdir(parents=True, exist_ok=True)


def main():
    curves, qc_rows = {}, []
    for cell, path in cross.OLD_FILES.items():
        cell_curves, rows = cross.parse_cell(
            "Old July BoL", cell, path, cross.OLD_BRANCH_STEPS
        )
        curves[cell] = cell_curves
        qc_rows.extend(rows)
    qc = pd.DataFrame(qc_rows)
    protocol = cross.cohort_protocol_means(qc)

    model0, stage, anode, cathode, _, _, qcell = selected.selected_inputs()
    scenario = p0.Scenario("P0 anchored reference")
    model = p0.build_model(model0, scenario)

    simulations, rows = {}, []
    for rate in ga.RATES:
        for charge in (True, False):
            direction = "Charge" if charge else "Discharge"
            pr = protocol.loc[("Old July BoL", rate, direction)]
            rest_v = float(pr.Rest_end_V)
            z0 = priority.monotone_voltage_inverse(model["v_qocv"], model["soc"], rest_v)
            local_stage = priority.stage_at_soc(stage, z0, charge)
            sim_rate = float(pr.Measured_current_A) / cross.NOMINAL_CAPACITY_AH
            print(f"P0 | {rate:g}C {direction}", flush=True)
            sim = ehq.run_dfn(
                sim_rate, charge, local_stage, model, anode, cathode, False, False
            )
            simulations[(rate, charge)] = sim
            for cell, branches in curves.items():
                obs = branches[(rate, charge)]
                m = priority.curve_metrics(obs, sim, qcell, charge)
                rows.append(
                    {
                        "Cell": cell,
                        "C_rate": rate,
                        "Direction": direction,
                        **m,
                        "Capacity_error_mAh": 1000.0
                        * (m["Q_end_model_Ah"] - m["Q_end_exp_Ah"]),
                    }
                )
    detail = pd.DataFrame(rows)
    detail.to_csv(RESULTS / "p0_curve_metrics_by_cell.csv", index=False)
    summary = detail.groupby(["C_rate", "Direction"], as_index=False).agg(
        Full_RMSE_mV=("Full_RMSE_mV", "mean"),
        Center10_70_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
        Capacity_error_pct=("Capacity_error_pct", "mean"),
        Capacity_error_mAh=("Capacity_error_mAh", "mean"),
    )
    summary.to_csv(RESULTS / "p0_curve_summary.csv", index=False)

    colors = {0.5: "#0072B2", 1.0: "#009E73", 2.0: "#D55E00"}
    for axis_kind, x_label, file_name in (
        ("t_min", "Time [min]", "p0_time_voltage.png"),
        ("Q_Ah", "Transferred capacity [Ah]", "p0_capacity_voltage.png"),
    ):
        fig, axes = plt.subplots(2, 3, figsize=(16.2, 8.8), constrained_layout=True)
        for row_i, charge in enumerate((True, False)):
            direction = "Charge" if charge else "Discharge"
            for col_i, rate in enumerate(ga.RATES):
                ax = axes[row_i, col_i]
                for cell, branches in curves.items():
                    obs = branches[(rate, charge)]
                    ax.plot(obs[axis_kind], obs["V"], color="black", alpha=0.28, lw=1.0)
                sim = simulations[(rate, charge)]
                b = summary[
                    np.isclose(summary.C_rate, rate) & (summary.Direction == direction)
                ].iloc[0]
                ax.plot(
                    sim[axis_kind], sim["V"], color=colors[rate], lw=2.1,
                    label=(
                        f"P0: center RMSE {b.Center10_70_RMSE_mV:.1f} mV\n"
                        f"capacity error {b.Capacity_error_mAh:+.1f} mAh"
                    ),
                )
                ax.set(title=f"{rate:g}C {direction}", xlabel=x_label, ylabel="Voltage [V]")
                ax.grid(alpha=0.2)
                ax.legend(fontsize=8)
        axes[0, 0].plot([], [], color="black", alpha=0.45, label="July BoL cells (n=3)")
        axes[0, 0].legend(fontsize=8)
        fig.suptitle("Current P0 before dynamic fitting | CC-only")
        fig.savefig(RESULTS / file_name, dpi=220)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(10.2, 5.8), constrained_layout=True)
    x = np.arange(len(ga.RATES))
    width = 0.36
    ch = summary[summary.Direction == "Charge"].set_index("C_rate").loc[list(ga.RATES)]
    dis = summary[summary.Direction == "Discharge"].set_index("C_rate").loc[list(ga.RATES)]
    ax.bar(x - width / 2, ch.Capacity_error_mAh, width, label="Charge")
    ax.bar(x + width / 2, dis.Capacity_error_mAh, width, label="Discharge")
    ax.axhline(0, color="black", lw=0.8)
    ax.set(xticks=x, xticklabels=[f"{r:g}C" for r in ga.RATES], ylabel="Model - experiment [mAh]", title="P0 CC cutoff capacity error")
    ax.grid(axis="y", alpha=0.2)
    ax.legend()
    fig.savefig(RESULTS / "p0_capacity_error.png", dpi=220)
    plt.close(fig)

    # Extra-R0 diagnostic: compare the measured and DFN-predicted initial loaded
    # voltage jump from the immediately preceding rest voltage.
    r_rows = []
    for rate in ga.RATES:
        for charge in (True, False):
            direction = "Charge" if charge else "Discharge"
            sim = simulations[(rate, charge)]
            sub = qc[np.isclose(qc.C_rate, rate) & (qc.Direction == direction)]
            model_rest = float(protocol.loc[("Old July BoL", rate, direction), "Rest_end_V"])
            model_loaded = float(sim["V"][0])
            model_jump = (model_loaded - model_rest) if charge else (model_rest - model_loaded)
            for _, row in sub.iterrows():
                obs = curves[row.Cell][(rate, charge)]
                exp_loaded = float(obs["V"][0])
                exp_jump = (exp_loaded - row.Rest_end_V) if charge else (row.Rest_end_V - exp_loaded)
                current = float(row.Median_current_A)
                r_rows.append(
                    {
                        "Cell": row.Cell,
                        "C_rate": rate,
                        "Direction": direction,
                        "Current_A": current,
                        "Experimental_jump_mV": 1000.0 * exp_jump,
                        "DFN_jump_mV": 1000.0 * model_jump,
                        "Unexplained_jump_mV": 1000.0 * (exp_jump - model_jump),
                        "Apparent_extra_R0_mOhm": 1000.0 * (exp_jump - model_jump) / current,
                    }
                )
    r0 = pd.DataFrame(r_rows)
    r0.to_csv(RESULTS / "r0_initial_jump_diagnostic.csv", index=False)
    r0_summary = r0.groupby(["C_rate", "Direction"], as_index=False).agg(
        Experimental_jump_mV=("Experimental_jump_mV", "mean"),
        DFN_jump_mV=("DFN_jump_mV", "mean"),
        Unexplained_jump_mV=("Unexplained_jump_mV", "mean"),
        Apparent_extra_R0_mOhm=("Apparent_extra_R0_mOhm", "mean"),
        Extra_R0_SD_mOhm=("Apparent_extra_R0_mOhm", "std"),
    )
    r0_summary.to_csv(RESULTS / "r0_initial_jump_summary.csv", index=False)

    manifest = {
        "baseline": "P0 anchored reference",
        "scope": "July BoL n=3, CC-only, no dynamic fitting",
        "r0_interpretation": "diagnostic only; an additional R0 should be added only if unexplained jump/I is positive and consistent across cells, rates, and directions",
    }
    (RESULTS / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("\nSUMMARY\n", summary.to_string(index=False))
    print("\nR0 DIAGNOSTIC\n", r0_summary.to_string(index=False))


if __name__ == "__main__":
    main()
