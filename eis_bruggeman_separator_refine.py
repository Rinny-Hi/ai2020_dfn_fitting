"""Refine separator Bruggeman for the Pareto electrode-Bruggeman candidates."""

from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import electrode_hysteresis_quantification as ehq
import eis_bruggeman_screen as screen
import eis_first_kinetics_application as eis
import gitt_ocp_analysis as ga
import mass_independent_nominal_anode_trial as mit
import nominal_radius_current_configuration as current
import ocp_grounded_candidate_search as common
import old_notebook_method_comparison as omc
import ocp_rn_sweep_fixed_transport as rs
import restored_brugg_de_csmax_audit as audit


ROOT = screen.ROOT
RESULTS = screen.RESULTS
ELECTRODE_CASES = ((2.5, 1.5), (2.5, 1.83), (2.5, 2.0))
BS_VALUES = (1.0, 1.25, 1.5, 1.75, 2.0)


def main():
    bundle = ga.load_legacy_bundle(ga.get_paths(ROOT))
    base = audit.build_base(bundle, 1e-4, audit.CSN_CORRECTED)
    _, _, hybrid, cathode = rs.candidate_details(bundle, base)
    experiment = omc.load_old_dynamic_data(base["dynamic"], base["q_meas"])
    candidates = pd.read_csv(screen.SOURCE / "mass_independent_ocp_candidate_metrics.csv")
    ocp_row = candidates[candidates.Candidate == current.NAME].iloc[0]
    stage = {key: float(ocp_row[key]) for key in ("x0", "x100", "y100", "y0")}
    effective = mit.model_for_effective_capacities(base, ocp_row)
    nominal = audit.set_dynamic_parameters(effective, current.RN_UM, current.RP_UM, 1e-4)
    eis_fixed = eis.with_kinetics(
        nominal,
        eis.KN_EIS_N4_ONLY / eis.AI2020_MREF,
        eis.KP_EIS / eis.AI2020_MREF,
    )
    anode = common.scaled_detail(hybrid, 0.50)
    cathode = common.scaled_detail(cathode, 0.50)
    rows = []
    for bn, bp in ELECTRODE_CASES:
        for bs in BS_VALUES:
            case = f"bn={bn:g}, bp={bp:g}, bs={bs:g}"
            print(case, flush=True)
            model = screen.set_bruggeman(eis_fixed, bn, bp, bs)
            for rate in ga.RATES:
                for charge in (True, False):
                    direction = "Charge" if charge else "Discharge"
                    sim = ehq.run_dfn(rate, charge, stage, model, anode, cathode, True, True)
                    result = omc.old_time_metrics(experiment[(rate, charge)], sim)
                    capacity_error = result["end_time_error_min"] * rate * 2.28 / 60.0
                    rows.append(
                        {
                            "Case": case,
                            "bn": bn,
                            "bp": bp,
                            "bs": bs,
                            "C_rate": rate,
                            "Direction": direction,
                            **result,
                            "capacity_error_mAh": capacity_error * 1000.0,
                        }
                    )
    branches = pd.DataFrame(rows)
    branches.to_csv(RESULTS / "eis_bruggeman_separator_branch_metrics.csv", index=False, encoding="utf-8-sig")
    summary = branches.groupby(["Case", "bn", "bp", "bs"], as_index=False, sort=False).agg(
        mean_MAE_SOC10_70_mV=("MAE_SOC10_70_mV", "mean"),
        mean_abs_capacity_error_mAh=("capacity_error_mAh", lambda s: float(np.mean(np.abs(s)))),
        max_abs_capacity_error_mAh=("capacity_error_mAh", lambda s: float(np.max(np.abs(s)))),
    )
    charge = branches[branches.Direction == "Charge"].groupby("Case")["MAE_SOC10_70_mV"].mean()
    discharge = branches[branches.Direction == "Discharge"].groupby("Case")["MAE_SOC10_70_mV"].mean()
    summary["charge_mean_MAE_mV"] = summary.Case.map(charge)
    summary["discharge_mean_MAE_mV"] = summary.Case.map(discharge)
    summary["Pareto"] = screen.pareto_mask(summary)
    summary.to_csv(RESULTS / "eis_bruggeman_separator_summary.csv", index=False, encoding="utf-8-sig")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    for (bn, bp), group in summary.groupby(["bn", "bp"], sort=False):
        label = f"bn={bn:g}, bp={bp:g}"
        group = group.sort_values("bs")
        axes[0].plot(group.bs, group.mean_MAE_SOC10_70_mV, marker="o", label=label)
        axes[1].plot(group.bs, group.mean_abs_capacity_error_mAh, marker="o", label=label)
    axes[0].set_ylabel("Mean voltage MAE (mV)")
    axes[1].set_ylabel("Mean absolute capacity error (mAh)")
    for ax in axes:
        ax.set_xlabel("Separator Bruggeman bs")
        ax.grid(alpha=0.2)
        ax.legend()
    fig.suptitle("Separator Bruggeman sensitivity with both EIS kinetics fixed")
    fig.savefig(RESULTS / "eis_bruggeman_separator_refine.png", dpi=220)
    plt.close(fig)

    pareto = summary[summary.Pareto].sort_values("mean_abs_capacity_error_mAh")
    report = {
        "fixed_kinetics": {"kn": eis.KN_EIS_N4_ONLY, "kp": eis.KP_EIS},
        "pareto": pareto.to_dict(orient="records"),
    }
    (RESULTS / "eis_bruggeman_separator_refine.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("\nPareto")
    print(pareto[["bn", "bp", "bs", "mean_MAE_SOC10_70_mV", "mean_abs_capacity_error_mAh", "max_abs_capacity_error_mAh"]].to_string(index=False))


if __name__ == "__main__":
    main()
