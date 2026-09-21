"""Apply provisional teardown thicknesses to the current EIS/Bruggeman model."""

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
import nominal_radius_current_configuration as current
import ocp_grounded_candidate_search as common
import old_notebook_method_comparison as omc
import ocp_rn_sweep_fixed_transport as rs
import restored_brugg_de_csmax_audit as audit


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "results" / "260920_mass_independent_nominal_anode"
RESULTS = ROOT / "results" / "260921_current_thickness_estimate"
RESULTS.mkdir(parents=True, exist_ok=True)

LN_ESTIMATE = 72.5e-6
LP_ESTIMATE = 61.5e-6


def set_thickness_preserve_capacity(model, ln, lp):
    changed = dict(model)
    changed["params"] = model["params"].copy()
    old_ln = float(model["params"]["Negative electrode thickness [m]"])
    old_lp = float(model["params"]["Positive electrode thickness [m]"])
    old_csn = float(model["params"]["Maximum concentration in negative electrode [mol.m-3]"])
    old_csp = float(model["params"]["Maximum concentration in positive electrode [mol.m-3]"])
    new_csn = old_csn * old_ln / ln
    new_csp = old_csp * old_lp / lp
    changed["params"].update(
        {
            "Negative electrode thickness [m]": ln,
            "Positive electrode thickness [m]": lp,
            "Maximum concentration in negative electrode [mol.m-3]": new_csn,
            "Maximum concentration in positive electrode [mol.m-3]": new_csp,
        },
        check_already_exists=False,
    )
    changed["csn_max"] = new_csn
    changed["csp_max"] = new_csp
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
    nominal = current.apply_current_kinetics_and_bruggeman(nominal)
    estimated = set_thickness_preserve_capacity(nominal, LN_ESTIMATE, LP_ESTIMATE)
    models = {
        "Ai2020 nominal thickness": nominal,
        "Estimated teardown thickness": estimated,
    }
    anode = common.scaled_detail(hybrid, 0.50)
    cathode = common.scaled_detail(cathode, 0.50)

    rows = []
    simulations = {}
    for scenario, model in models.items():
        print(scenario, flush=True)
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                sim = ehq.run_dfn(rate, charge, stage, model, anode, cathode, True, True)
                simulations[(scenario, rate, charge)] = sim
                metric = omc.old_time_metrics(experiment[(rate, charge)], sim)
                capacity_error = metric["end_time_error_min"] * rate * 2.28 / 60.0
                rows.append(
                    {
                        "Scenario": scenario,
                        "Ln_um": float(model["params"]["Negative electrode thickness [m]"]) * 1e6,
                        "Lp_um": float(model["params"]["Positive electrode thickness [m]"]) * 1e6,
                        "csn_max_mol_m3": float(model["params"]["Maximum concentration in negative electrode [mol.m-3]"]),
                        "csp_max_mol_m3": float(model["params"]["Maximum concentration in positive electrode [mol.m-3]"]),
                        "C_rate": rate,
                        "Direction": direction,
                        **metric,
                        "capacity_error_mAh": capacity_error * 1000.0,
                    }
                )
    branches = pd.DataFrame(rows)
    branches.to_csv(RESULTS / "thickness_branch_metrics.csv", index=False, encoding="utf-8-sig")
    summary = branches.groupby(
        ["Scenario", "Ln_um", "Lp_um", "csn_max_mol_m3", "csp_max_mol_m3"],
        as_index=False,
        sort=False,
    ).agg(
        mean_MAE_SOC10_70_mV=("MAE_SOC10_70_mV", "mean"),
        mean_abs_capacity_error_mAh=("capacity_error_mAh", lambda s: float(np.mean(np.abs(s)))),
        max_abs_capacity_error_mAh=("capacity_error_mAh", lambda s: float(np.max(np.abs(s)))),
    )
    summary.to_csv(RESULTS / "thickness_summary.csv", index=False, encoding="utf-8-sig")

    colors = {"Ai2020 nominal thickness": "#ff7f0e", "Estimated teardown thickness": "#0072B2"}
    fig, axes = plt.subplots(2, 3, figsize=(16, 8.5), constrained_layout=True)
    for col, rate in enumerate(ga.RATES):
        for row_idx, charge in enumerate((True, False)):
            ax = axes[row_idx, col]
            exp = experiment[(rate, charge)]
            ax.plot(exp["t_min"], exp["V"], color="black", lw=2.4, label="Experiment")
            for scenario in models:
                sim = simulations[(scenario, rate, charge)]
                metric = branches[
                    (branches.Scenario == scenario)
                    & (branches.C_rate == rate)
                    & (branches.Direction == ("Charge" if charge else "Discharge"))
                ].iloc[0]
                ax.plot(sim["t_min"], sim["V"], color=colors[scenario], lw=1.8, label=f"{scenario} ({metric.MAE_SOC10_70_mV:.1f} mV)")
            ax.set_title(f"{rate:g}C {'Charge' if charge else 'Discharge'}")
            ax.set_xlabel("Time (min)")
            ax.set_ylabel("Voltage (V)")
            ax.grid(alpha=0.2)
            ax.legend(fontsize=7.5)
    fig.suptitle("Provisional coating thicknesses with Qn/Qp preserved")
    fig.savefig(RESULTS / "thickness_dynamic_comparison.png", dpi=220)
    plt.close(fig)

    report = {
        "assumptions": {
            "positive_double_sided_total_um": 138.0,
            "negative_double_sided_total_um": 155.0,
            "assumed_Al_collector_um": 15.0,
            "assumed_Cu_collector_um": 10.0,
            "formula": "single coating thickness = (double-sided total - collector) / 2",
            "capacity_constraint": "csmax rescaled inversely with thickness to preserve fitted Qn and Qp",
        },
        "fixed_current_configuration": {
            "kn": current.KN_PREF,
            "kp": current.KP_PREF,
            "Bruggeman": [current.BRUGG_N, current.BRUGG_P, current.BRUGG_S],
            "qOCV_MAE_mV": float(ocp_row.qOCV_MAE_2_98_mV),
        },
        "summary": summary.to_dict(orient="records"),
    }
    (RESULTS / "current_thickness_estimate.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("\nSummary")
    print(summary.to_string(index=False))
    print("\nEstimated branches")
    print(branches[branches.Scenario == "Estimated teardown thickness"][["C_rate", "Direction", "MAE_SOC10_70_mV", "capacity_error_mAh"]].to_string(index=False))


if __name__ == "__main__":
    main()
