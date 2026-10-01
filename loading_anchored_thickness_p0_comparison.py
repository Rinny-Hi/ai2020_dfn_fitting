"""Compare three thickness measurements after applying the remeasured loading.

The supplied 17.35 mg (positive) and 10.10 mg (negative) values are already
single-face coating masses: (whole punched electrode - collector) / 2.  Since
the active-material mass fractions remain confidential, the absolute split of
solid volume between active material/additives cannot be reconstructed.  The
defensible model update is therefore the measured new/old coating-mass ratio,
followed by a same-inventory thickness comparison in which porosity closes the
volume balance.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import comprehensive_prefit_cross_cohort as cross
import c50_selected_ocp_dynamic_validation as selected
import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import prefit_physical_anchor_decision as p0
import priority_initial_state_protocol_recheck as priority
import remeasured_electrode_mass_sensitivity as mass_study


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260927_loading_anchored_thickness_p0_comparison"
RESULTS.mkdir(parents=True, exist_ok=True)

CASES = {
    "Ai2020 paper": (76.5, 68.0),
    "Thickness gauge": (73.0, 62.5),
    "SEM cross-section": (81.52, 62.965),
}


def load_july_cells():
    curves = {}
    qc_rows = []
    for cell, path in cross.OLD_FILES.items():
        cell_curves, rows = cross.parse_cell(
            "Old July BoL", cell, path, cross.OLD_BRANCH_STEPS
        )
        curves[cell] = cell_curves
        qc_rows.extend(rows)
    qc = pd.DataFrame(qc_rows)
    return curves, cross.cohort_protocol_means(qc)


def build_case(model0: dict, ln_um: float, lp_um: float) -> dict:
    scenario = p0.Scenario(
        "loading-thickness case",
        ln_um=ln_um,
        lp_um=lp_um,
        rn_um=3.7542,
        rp_um=4.2387,
        kn=p0.AI_K,
        kp=4.18e-7,
        positive_hysteresis=False,
        initialization="rest",
    )
    anchored = p0.build_model(model0, scenario)
    return mass_study.apply_mass_ratio(anchored)


def main():
    cell_curves, protocol = load_july_cells()
    model0, stage, anode, cathode, _, _, qcell = selected.selected_inputs()
    models = {name: build_case(model0, *dims) for name, dims in CASES.items()}

    rows = []
    simulations = {}
    for name, model in models.items():
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                p = protocol.loc[("Old July BoL", rate, direction)]
                z0 = priority.monotone_voltage_inverse(
                    model["v_qocv"], model["soc"], float(p.Rest_end_V)
                )
                local_stage = priority.stage_at_soc(stage, z0, charge)
                sim_rate = float(p.Measured_current_A) / cross.NOMINAL_CAPACITY_AH
                print(f"{name} | {rate:g}C {direction}", flush=True)
                sim = ehq.run_dfn(
                    sim_rate,
                    charge,
                    local_stage,
                    model,
                    anode,
                    cathode,
                    False,
                    False,
                )
                simulations[(name, rate, charge)] = sim
                for cell, branches in cell_curves.items():
                    metric = priority.curve_metrics(
                        branches[(rate, charge)], sim, qcell, charge
                    )
                    rows.append(
                        {
                            "Case": name,
                            "Cell": cell,
                            "C_rate": rate,
                            "Direction": direction,
                            "Ln_um": CASES[name][0],
                            "Lp_um": CASES[name][1],
                            "eps_s_n": float(model["params"]["Negative electrode active material volume fraction"]),
                            "eps_e_n": float(model["params"]["Negative electrode porosity"]),
                            "eps_s_p": float(model["params"]["Positive electrode active material volume fraction"]),
                            "eps_e_p": float(model["params"]["Positive electrode porosity"]),
                            **metric,
                            "Capacity_error_mAh": 1000.0
                            * (metric["Q_end_model_Ah"] - metric["Q_end_exp_Ah"]),
                        }
                    )

    detail = pd.DataFrame(rows)
    detail.to_csv(RESULTS / "loading_thickness_detail.csv", index=False)

    summary = (
        detail.groupby("Case", sort=False, as_index=False)
        .agg(
            Full_RMSE_mV=("Full_RMSE_mV", "mean"),
            Center10_70_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
            Capacity_RMSE_pct=(
                "Capacity_error_pct",
                lambda x: float(np.sqrt(np.mean(np.asarray(x, float) ** 2))),
            ),
            Mean_abs_capacity_error_mAh=(
                "Capacity_error_mAh", lambda x: float(np.mean(np.abs(x)))
            ),
            Max_abs_capacity_error_mAh=(
                "Capacity_error_mAh", lambda x: float(np.max(np.abs(x)))
            ),
        )
    )
    summary.to_csv(RESULTS / "loading_thickness_summary.csv", index=False)

    branch = (
        detail.groupby(["Case", "C_rate", "Direction"], sort=False, as_index=False)
        .agg(
            Full_RMSE_mV=("Full_RMSE_mV", "mean"),
            Center10_70_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
            Capacity_error_pct=("Capacity_error_pct", "mean"),
            Capacity_error_mAh=("Capacity_error_mAh", "mean"),
        )
    )
    branch.to_csv(RESULTS / "loading_thickness_branch_summary.csv", index=False)

    geometry = []
    for name, model in models.items():
        geometry.append(
            {
                "Case": name,
                "Ln_um": CASES[name][0],
                "Lp_um": CASES[name][1],
                "eps_s_n": float(model["params"]["Negative electrode active material volume fraction"]),
                "eps_e_n": float(model["params"]["Negative electrode porosity"]),
                "eps_s_p": float(model["params"]["Positive electrode active material volume fraction"]),
                "eps_e_p": float(model["params"]["Positive electrode porosity"]),
            }
        )
    pd.DataFrame(geometry).to_csv(RESULTS / "loading_thickness_geometry.csv", index=False)

    colors = {
        "Ai2020 paper": "#777777",
        "Thickness gauge": "#0072B2",
        "SEM cross-section": "#D55E00",
    }
    fig, axes = plt.subplots(2, 3, figsize=(16.2, 8.8), constrained_layout=True)
    for row_i, charge in enumerate((True, False)):
        for col_i, rate in enumerate(ga.RATES):
            ax = axes[row_i, col_i]
            for _, branches in cell_curves.items():
                obs = branches[(rate, charge)]
                ax.plot(obs["Q_Ah"], obs["V"], color="black", lw=0.9, alpha=0.27)
            direction = "Charge" if charge else "Discharge"
            for name in CASES:
                sim = simulations[(name, rate, charge)]
                b = branch[
                    (branch.Case == name)
                    & np.isclose(branch.C_rate, rate)
                    & (branch.Direction == direction)
                ].iloc[0]
                ax.plot(
                    sim["Q_Ah"],
                    sim["V"],
                    color=colors[name],
                    lw=1.8,
                    label=(
                        f"{name}: center {b.Center10_70_RMSE_mV:.1f} mV, "
                        f"dQ {b.Capacity_error_mAh:+.1f} mAh"
                    ),
                )
            ax.set(
                title=f"{rate:g}C {direction}",
                xlabel="Transferred capacity [Ah]",
                ylabel="Voltage [V]",
            )
            ax.grid(alpha=0.2)
    axes[0, 0].legend(fontsize=6.5)
    fig.suptitle("Measured loading + P0 anchors: three thickness cases")
    fig.savefig(RESULTS / "loading_thickness_capacity_voltage.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.5, 6.3), constrained_layout=True)
    for _, r in summary.iterrows():
        ax.scatter(r.Capacity_RMSE_pct, r.Center10_70_RMSE_mV, s=90, color=colors[r.Case])
        ax.annotate(r.Case, (r.Capacity_RMSE_pct, r.Center10_70_RMSE_mV), xytext=(6, 5), textcoords="offset points")
    ax.set(xlabel="Capacity RMSE [%]", ylabel="10-70% voltage RMSE [mV]", title="Thickness decision trade-off")
    ax.grid(alpha=0.25)
    fig.savefig(RESULTS / "loading_thickness_tradeoff.png", dpi=220)
    plt.close(fig)

    manifest = {
        "coating_mass_definition": "single-face coating mass = (punched electrode - collector) / 2",
        "new_mass_mg": mass_study.NEW_MASS_MG,
        "old_mass_mg": mass_study.OLD_MASS_MG,
        "mass_ratio": mass_study.MASS_RATIO,
        "punch_diameter_mm": 10.0,
        "loading_mg_cm2": {
            key: value / (np.pi * 0.5**2)
            for key, value in mass_study.NEW_MASS_MG.items()
        },
        "fixed": {
            "OCP": "C/50 selected nominal anode + old cathode",
            "Rn_Rp_um": [3.7542, 4.2387],
            "Dsn_Dsp_m2_s": [2.1e-14, 4.4e-14],
            "kn_kp": [p0.AI_K, 4.18e-7],
            "Bruggeman_n_p_s": [2.914, 1.83, 1.5],
            "hysteresis": "negative off, positive off",
            "initialization": "preceding-rest OCP inversion",
        },
        "limitation": "active/additive mass fractions unavailable; absolute porosity is represented by scaling the established effective solid inventory by the measured new/old coating-mass ratio",
    }
    (RESULTS / "analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print("\nGEOMETRY\n", pd.DataFrame(geometry).to_string(index=False))
    print("\nSUMMARY\n", summary.to_string(index=False))
    print("\nBRANCH\n", branch.to_string(index=False))


if __name__ == "__main__":
    main()
