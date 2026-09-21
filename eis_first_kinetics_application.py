"""Apply first-pass EIS kinetic prefactors to the current nominal-radius DFN."""

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
import mass_independent_nominal_anode_trial as mit
import nominal_radius_current_configuration as current
import ocp_grounded_candidate_search as common
import old_notebook_method_comparison as omc
import ocp_rn_sweep_fixed_transport as rs
import restored_brugg_de_csmax_audit as audit


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "results" / "260920_mass_independent_nominal_anode"
RESULTS = ROOT / "results" / "260921_eis_first_kinetics"
RESULTS.mkdir(parents=True, exist_ok=True)

# Values read from the user's EIS first-pass result images.  They have the same
# dimensions as the Dualfoil/Ai2020 m_ref prefactor, not a dimensionless factor.
KN_EIS = 2.48e-7
KN_EIS_N4_ONLY = 7.40e-7
KP_EIS = 3.12e-7
FARADAY = 96485.33212
AI2020_MREF = 1.0e-11 * FARADAY

CASES = {
    "Current baseline": (1.0, 1.0),
    "EIS kn only": (KN_EIS / AI2020_MREF, 1.0),
    "EIS kp only": (1.0, KP_EIS / AI2020_MREF),
    "EIS kn + kp": (KN_EIS / AI2020_MREF, KP_EIS / AI2020_MREF),
    "EIS kn N4 only": (KN_EIS_N4_ONLY / AI2020_MREF, 1.0),
    "EIS kn N4 + kp": (KN_EIS_N4_ONLY / AI2020_MREF, KP_EIS / AI2020_MREF),
}


def with_kinetics(model, kn_multiplier, kp_multiplier):
    changed = dict(model)
    changed["params"] = model["params"].copy()
    raw = pybamm.ParameterValues("Ai2020")
    changed["params"].update(
        {
            "Negative electrode exchange-current density [A.m-2]": ga.scale_parameter_function(
                raw["Negative electrode exchange-current density [A.m-2]"], kn_multiplier
            ),
            "Positive electrode exchange-current density [A.m-2]": ga.scale_parameter_function(
                raw["Positive electrode exchange-current density [A.m-2]"], kp_multiplier
            ),
        },
        check_already_exists=False,
    )
    return changed


def main():
    bundle = ga.load_legacy_bundle(ga.get_paths(ROOT))
    base = audit.build_base(bundle, 1e-4, audit.CSN_CORRECTED)
    _, _, hybrid, cathode = rs.candidate_details(bundle, base)
    experiment = omc.load_old_dynamic_data(base["dynamic"], base["q_meas"])
    candidates = pd.read_csv(SOURCE / "mass_independent_ocp_candidate_metrics.csv")
    ocp_row = candidates[candidates.Candidate == current.NAME].iloc[0]
    stage = {key: float(ocp_row[key]) for key in ("x0", "x100", "y100", "y0")}
    effective = mit.model_for_effective_capacities(base, ocp_row)
    nominal = audit.set_dynamic_parameters(effective, current.RN_UM, current.RP_UM, 1e-4)
    anode = common.scaled_detail(hybrid, 0.50)
    cathode_dynamic = common.scaled_detail(cathode, 0.50)

    rows = []
    simulations = {}
    for case, (kn_multiplier, kp_multiplier) in CASES.items():
        model = with_kinetics(nominal, kn_multiplier, kp_multiplier)
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                print(f"{case} | {rate:g}C {direction}", flush=True)
                sim = ehq.run_dfn(
                    rate,
                    charge,
                    stage,
                    model,
                    anode,
                    cathode_dynamic,
                    True,
                    True,
                )
                simulations[(case, rate, charge)] = sim
                result = omc.old_time_metrics(experiment[(rate, charge)], sim)
                capacity_error = result["end_time_error_min"] * rate * 2.28 / 60.0
                rows.append(
                    {
                        "Case": case,
                        "kn_prefactor": AI2020_MREF * kn_multiplier,
                        "kp_prefactor": AI2020_MREF * kp_multiplier,
                        "kn_multiplier_vs_Ai2020": kn_multiplier,
                        "kp_multiplier_vs_Ai2020": kp_multiplier,
                        "C_rate": rate,
                        "Direction": direction,
                        **result,
                        "capacity_error_Ah": capacity_error,
                        "capacity_error_mAh": capacity_error * 1000.0,
                    }
                )

    branches = pd.DataFrame(rows)
    branches.to_csv(RESULTS / "eis_kinetics_branch_metrics.csv", index=False, encoding="utf-8-sig")
    summary = branches.groupby("Case", as_index=False, sort=False).agg(
        mean_MAE_SOC10_70_mV=("MAE_SOC10_70_mV", "mean"),
        mean_abs_capacity_error_mAh=("capacity_error_mAh", lambda s: float(np.mean(np.abs(s)))),
        max_abs_capacity_error_mAh=("capacity_error_mAh", lambda s: float(np.max(np.abs(s)))),
    )
    charge_summary = branches[branches.Direction == "Charge"].groupby("Case", sort=False)["MAE_SOC10_70_mV"].mean()
    discharge_summary = branches[branches.Direction == "Discharge"].groupby("Case", sort=False)["MAE_SOC10_70_mV"].mean()
    summary["charge_mean_MAE_SOC10_70_mV"] = summary.Case.map(charge_summary)
    summary["discharge_mean_MAE_SOC10_70_mV"] = summary.Case.map(discharge_summary)
    summary.to_csv(RESULTS / "eis_kinetics_summary.csv", index=False, encoding="utf-8-sig")

    colors = {
        "Current baseline": "#ff7f0e",
        "EIS kn + kp": "#999999",
        "EIS kn N4 only": "#009E73",
        "EIS kn N4 + kp": "#0072B2",
    }
    fig, axes = plt.subplots(2, 3, figsize=(16, 8.5), constrained_layout=True)
    for col, rate in enumerate(ga.RATES):
        for row_index, charge in enumerate((True, False)):
            ax = axes[row_index, col]
            exp = experiment[(rate, charge)]
            ax.plot(exp["t_min"], exp["V"], color="black", lw=2.4, label="Experiment")
            for case in ("Current baseline", "EIS kn + kp", "EIS kn N4 only", "EIS kn N4 + kp"):
                sim = simulations[(case, rate, charge)]
                metric = branches[
                    (branches.Case == case)
                    & (branches.C_rate == rate)
                    & (branches.Direction == ("Charge" if charge else "Discharge"))
                ].iloc[0]
                ax.plot(
                    sim["t_min"], sim["V"], color=colors[case], lw=1.9,
                    label=f"{case} ({metric.MAE_SOC10_70_mV:.1f} mV)",
                )
            ax.set_title(f"{rate:g}C {'Charge' if charge else 'Discharge'}")
            ax.set_xlabel("Time (min)")
            ax.set_ylabel("Voltage (V)")
            ax.grid(alpha=0.2)
            ax.legend(fontsize=8)
    fig.suptitle("First EIS kinetics applied to current nominal-radius configuration")
    fig.savefig(RESULTS / "eis_kinetics_charge_discharge_comparison.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), constrained_layout=True)
    x = np.arange(len(summary))
    width = 0.36
    axes[0].bar(x - width / 2, summary.charge_mean_MAE_SOC10_70_mV, width, label="Charge", color="#D55E00")
    axes[0].bar(x + width / 2, summary.discharge_mean_MAE_SOC10_70_mV, width, label="Discharge", color="#0072B2")
    axes[0].set_ylabel("Mean MAE, 10-70% SOC (mV)")
    axes[0].set_xticks(x, summary.Case, rotation=18, ha="right")
    axes[0].legend()
    axes[0].grid(axis="y", alpha=0.2)
    axes[1].bar(x, summary.mean_abs_capacity_error_mAh, color="#009E73")
    axes[1].set_ylabel("Mean absolute capacity error (mAh)")
    axes[1].set_xticks(x, summary.Case, rotation=18, ha="right")
    axes[1].grid(axis="y", alpha=0.2)
    fig.suptitle("Kinetic-prefactor ablation: which electrode drives the change?")
    fig.savefig(RESULTS / "eis_kinetics_ablation_summary.png", dpi=220)
    plt.close(fig)

    report = {
        "unit_mapping": {
            "Ai2020_m_ref": AI2020_MREF,
            "kn_EIS": KN_EIS,
            "kp_EIS": KP_EIS,
            "kn_EIS_N4_only": KN_EIS_N4_ONLY,
            "kn_multiplier_vs_Ai2020": KN_EIS / AI2020_MREF,
            "kp_multiplier_vs_Ai2020": KP_EIS / AI2020_MREF,
            "kn_N4_only_multiplier_vs_Ai2020": KN_EIS_N4_ONLY / AI2020_MREF,
            "interpretation": "EIS values replace the Ai2020/Dualfoil m_ref prefactors.",
        },
        "unchanged_configuration": {
            "Rn_um": current.RN_UM,
            "Rp_um": current.RP_UM,
            "Dsn_m2_s": rs.DSN,
            "Dsp_m2_s": rs.DSP,
            "Bruggeman_n_p_s": [2.914, 1.83, 1.5],
            "De_multiplier_vs_raw_Ai2020_function": 1e-4,
            "OCP_candidate": current.NAME,
            "qOCV_MAE_mV": float(ocp_row.qOCV_MAE_2_98_mV),
        },
        "summary": summary.to_dict(orient="records"),
        "branches": branches.to_dict(orient="records"),
    }
    (RESULTS / "eis_first_kinetics_application.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("\nUnit mapping")
    print(json.dumps(report["unit_mapping"], indent=2))
    print("\nSummary")
    print(summary.to_string(index=False))
    print("\nBranches: current baseline vs EIS combinations")
    print(
        branches[branches.Case.isin(["Current baseline", "EIS kn + kp", "EIS kn N4 + kp"])][
            ["Case", "C_rate", "Direction", "MAE_SOC10_70_mV", "capacity_error_mAh"]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
