"""Screen Rn/Rp combinations and re-check OCP candidates at the best pair.

Fixed conditions follow the user-requested transport case:
  Dsn=2.1e-14 m2/s, Dsp=4.4e-14 m2/s,
  electrolyte Bruggeman b_n=b_p=b_s=1.5.

The full 3x3 radius screen is run with the best prior OCP candidate
(Ai2020 nominal equilibrium + harvested hysteresis).  The existing Rp=3 um
slice is reused.  All OCP candidates are then validated at the best radius
pair to check whether the OCP ranking changes.
"""

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
import mass_independent_nominal_anode_trial as mit
import ocp_grounded_candidate_search as common
import old_notebook_method_comparison as omc
import ocp_rn_sweep_fixed_transport as rs


ROOT = Path(__file__).resolve().parent
PRIOR = ROOT / "results" / "260920_ocp_rn_sweep_fixed_transport"
OCP_SOURCE = ROOT / "results" / "260920_mass_independent_nominal_anode"
RESULTS = ROOT / "results" / "260921_ocp_rn_rp_combination"
RESULTS.mkdir(parents=True, exist_ok=True)

RN_VALUES_UM = (2.0, 3.0, 4.0)
RP_VALUES_UM = (2.0, 3.0, 4.0)
HYBRID = "Ai2020 nominal equilibrium + harvested hysteresis"
NAMES = [
    "Previous measured-anode ±5% control",
    "Measured anode, mass-independent fit",
    "Ai2020 nominal anode, mass-independent fit",
    HYBRID,
    "Ai2020 nominal anode, published stoichiometry window",
]


def simulate_case(
    name,
    rn_um,
    rp_um,
    metrics,
    base,
    experiment,
    measured,
    nominal,
    hybrid,
    cathode,
    simulations,
):
    metric_row = metrics[metrics.Candidate == name].iloc[0]
    effective_model = mit.model_for_effective_capacities(base, metric_row)
    model = rs.update_transport(effective_model, rn_um, rp_um * 1e-6)
    stage = {key: float(metric_row[key]) for key in ("x0", "x100", "y100", "y0")}
    if name == HYBRID:
        anode = common.scaled_detail(hybrid, 0.50)
        negative_hysteresis = True
    elif name.startswith("Ai2020"):
        anode = nominal
        negative_hysteresis = False
    else:
        anode = common.scaled_detail(measured, 0.50)
        negative_hysteresis = True
    cathode_dynamic = common.scaled_detail(cathode, 0.50)
    rows = []
    for rate in ga.RATES:
        for charge in (True, False):
            direction = "Charge" if charge else "Discharge"
            print(
                f"{name} | Rn={rn_um:g} um | Rp={rp_um:g} um | {rate:g}C | {direction}",
                flush=True,
            )
            sim = ehq.run_dfn(
                rate,
                charge,
                stage,
                model,
                anode,
                cathode_dynamic,
                negative_hysteresis,
                True,
            )
            simulations[(name, rn_um, rp_um, rate, charge)] = sim
            result = omc.old_time_metrics(experiment[(rate, charge)], sim)
            capacity_error = result["end_time_error_min"] * rate * 2.28 / 60.0
            rows.append(
                {
                    "Candidate": name,
                    "Rn_um": rn_um,
                    "Rp_um": rp_um,
                    "Dsn_m2_s": rs.DSN,
                    "Dsp_m2_s": rs.DSP,
                    "Bruggeman_electrolyte": rs.BRUGGEMAN,
                    "C_rate": rate,
                    "Direction": direction,
                    **result,
                    "capacity_error_Ah": capacity_error,
                    "status": "ok",
                }
            )
    return rows


def summarize(frame):
    return frame.groupby(["Candidate", "Rn_um", "Rp_um"], as_index=False).agg(
        Dynamic_MAE_SOC10_70_mV=("MAE_SOC10_70_mV", "mean"),
        Charge_MAE_SOC10_70_mV=("MAE_SOC10_70_mV", lambda x: float(x[frame.loc[x.index, "Direction"] == "Charge"].mean())),
        Discharge_MAE_SOC10_70_mV=("MAE_SOC10_70_mV", lambda x: float(x[frame.loc[x.index, "Direction"] == "Discharge"].mean())),
        Mean_abs_capacity_error_mAh=("capacity_error_Ah", lambda x: float(np.mean(np.abs(x))) * 1000.0),
        Max_abs_capacity_error_mAh=("capacity_error_Ah", lambda x: float(np.max(np.abs(x))) * 1000.0),
    )


def main():
    bundle = ga.load_legacy_bundle(ga.get_paths(ROOT))
    base = ga.build_legacy_model_inputs(bundle)
    measured, nominal, hybrid, cathode = rs.candidate_details(bundle, base)
    experiment = omc.load_old_dynamic_data(base["dynamic"], base["q_meas"])
    metrics = pd.read_csv(OCP_SOURCE / "mass_independent_ocp_candidate_metrics.csv")
    simulations = {}

    # Reuse the already-computed Rp=3 um hybrid slice.
    prior = pd.read_csv(PRIOR / "ocp_rn_dynamic_validation.csv")
    grid_rows = prior[(prior.Candidate == HYBRID) & (prior.Rp_um == 3.0)].copy().to_dict(orient="records")
    for rp_um in (2.0, 4.0):
        for rn_um in RN_VALUES_UM:
            grid_rows.extend(
                simulate_case(
                    HYBRID,
                    rn_um,
                    rp_um,
                    metrics,
                    base,
                    experiment,
                    measured,
                    nominal,
                    hybrid,
                    cathode,
                    simulations,
                )
            )
    grid_validation = pd.DataFrame(grid_rows)
    grid_summary = summarize(grid_validation)
    grid_summary["tau_n_s"] = (grid_summary.Rn_um * 1e-6) ** 2 / rs.DSN
    grid_summary["tau_p_s"] = (grid_summary.Rp_um * 1e-6) ** 2 / rs.DSP
    best_geometry = grid_summary.sort_values(
        ["Dynamic_MAE_SOC10_70_mV", "Mean_abs_capacity_error_mAh"]
    ).iloc[0]
    best_rn = float(best_geometry.Rn_um)
    best_rp = float(best_geometry.Rp_um)

    # Re-evaluate all OCP sources at the selected geometry.
    ocp_rows = []
    for name in NAMES:
        ocp_rows.extend(
            simulate_case(
                name,
                best_rn,
                best_rp,
                metrics,
                base,
                experiment,
                measured,
                nominal,
                hybrid,
                cathode,
                simulations,
            )
        )
    ocp_validation = pd.DataFrame(ocp_rows)
    ocp_summary = summarize(ocp_validation)
    qcols = ["Candidate", "qOCV_MAE_2_98_mV", "qOCV_RMSE_2_98_mV", "qOCV_max_endpoint_abs_mV"]
    ocp_summary = ocp_summary.merge(metrics[qcols], on="Candidate", how="left")

    grid_validation.to_csv(RESULTS / "rn_rp_grid_dynamic_validation.csv", index=False, encoding="utf-8-sig")
    grid_summary.to_csv(RESULTS / "rn_rp_grid_summary.csv", index=False, encoding="utf-8-sig")
    ocp_validation.to_csv(RESULTS / "best_geometry_ocp_validation.csv", index=False, encoding="utf-8-sig")
    ocp_summary.to_csv(RESULTS / "best_geometry_ocp_summary.csv", index=False, encoding="utf-8-sig")

    pivot_mae = grid_summary.pivot(index="Rn_um", columns="Rp_um", values="Dynamic_MAE_SOC10_70_mV")
    pivot_cap = grid_summary.pivot(index="Rn_um", columns="Rp_um", values="Mean_abs_capacity_error_mAh")
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.5), constrained_layout=True)
    for ax, pivot, title, cbar_label in (
        (axes[0], pivot_mae, "Dynamic voltage MAE", "mV"),
        (axes[1], pivot_cap, "Mean absolute capacity error", "mAh"),
    ):
        im = ax.imshow(pivot.to_numpy(), origin="lower", cmap="viridis")
        ax.set_xticks(range(len(pivot.columns)), [f"{v:g}" for v in pivot.columns])
        ax.set_yticks(range(len(pivot.index)), [f"{v:g}" for v in pivot.index])
        ax.set_xlabel("Rp (um)")
        ax.set_ylabel("Rn (um)")
        ax.set_title(title)
        for i in range(len(pivot.index)):
            for j in range(len(pivot.columns)):
                value = pivot.iloc[i, j]
                ax.text(j, i, f"{value:.1f}", ha="center", va="center", color="white" if value > pivot.to_numpy().mean() else "black")
        fig.colorbar(im, ax=ax, label=cbar_label)
    fig.suptitle("Ai2020 nominal equilibrium + measured hysteresis | Rn/Rp screen")
    fig.savefig(RESULTS / "rn_rp_grid_heatmap.png", dpi=220)
    plt.close(fig)

    ranked = ocp_summary.sort_values("Dynamic_MAE_SOC10_70_mV")
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.5), constrained_layout=True)
    short = {
        NAMES[0]: "Measured ±5%",
        NAMES[1]: "Measured free",
        NAMES[2]: "Ai nominal",
        NAMES[3]: "Ai nominal+hyst",
        NAMES[4]: "Ai published",
    }
    labels = [short[v] for v in ranked.Candidate]
    axes[0].bar(labels, ranked.Dynamic_MAE_SOC10_70_mV, color="#4c78a8")
    axes[0].set_ylabel("Dynamic MAE (mV)")
    axes[0].set_title(f"OCP comparison at Rn={best_rn:g}, Rp={best_rp:g} um")
    axes[1].bar(labels, ranked.Mean_abs_capacity_error_mAh, color="#f58518")
    axes[1].set_ylabel("Mean absolute capacity error (mAh)")
    axes[1].set_title("Termination-capacity error")
    for ax in axes:
        ax.tick_params(axis="x", rotation=22)
        ax.grid(axis="y", alpha=0.25)
    fig.savefig(RESULTS / "best_geometry_ocp_comparison.png", dpi=220)
    plt.close(fig)

    best_name = str(ranked.iloc[0].Candidate)
    fig, axes = plt.subplots(2, 3, figsize=(16, 8.5), constrained_layout=True)
    for col, rate in enumerate(ga.RATES):
        for row_index, charge in enumerate((True, False)):
            ax = axes[row_index, col]
            exp = experiment[(rate, charge)]
            sim = simulations[(best_name, best_rn, best_rp, rate, charge)]
            ax.plot(exp["t_min"], exp["V"], color="black", lw=2.3, label="Experiment")
            ax.plot(sim["t_min"], sim["V"], color="#ff7f0e", lw=2.0, label="Model")
            ax.set_title(f"{rate:g}C {'Charge' if charge else 'Discharge'}")
            ax.set_xlabel("Time (min)")
            ax.set_ylabel("Voltage (V)")
            ax.grid(alpha=0.2)
    axes[0, 0].legend()
    fig.suptitle(f"{best_name} | Rn={best_rn:g} um, Rp={best_rp:g} um")
    fig.savefig(RESULTS / "best_geometry_dynamic_curves.png", dpi=220)
    plt.close(fig)

    report = {
        "fixed": {
            "Dsn_m2_s": rs.DSN,
            "Dsp_m2_s": rs.DSP,
            "electrolyte_Bruggeman": rs.BRUGGEMAN,
            "solid_phase_Bruggeman": 0.0,
        },
        "radius_input_conversion": "Values labelled um are multiplied by 1e-6 before PyBaMM input in metres.",
        "best_geometry_from_hybrid_screen": best_geometry.to_dict(),
        "ocp_ranking_at_best_geometry": ranked.to_dict(orient="records"),
        "method": "3x3 radius screen with hybrid OCP; then all OCP candidates revalidated at selected geometry",
    }
    (RESULTS / "rn_rp_combination_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print("\nRADIUS GRID")
    print(grid_summary[["Rn_um", "Rp_um", "Dynamic_MAE_SOC10_70_mV", "Mean_abs_capacity_error_mAh", "Max_abs_capacity_error_mAh"]].sort_values(["Rn_um", "Rp_um"]).to_string(index=False))
    print("\nBEST GEOMETRY", best_rn, best_rp)
    print("\nOCP RANKING")
    print(ranked[["Candidate", "Dynamic_MAE_SOC10_70_mV", "Mean_abs_capacity_error_mAh", "Max_abs_capacity_error_mAh", "qOCV_MAE_2_98_mV"]].to_string(index=False))
    print("\nSaved to", RESULTS)


if __name__ == "__main__":
    main()
