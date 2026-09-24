"""Current recommended configuration using Ai2020 nominal particle radii."""

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
import ocp_grounded_candidate_search as common
import old_notebook_method_comparison as omc
import ocp_rn_sweep_fixed_transport as rs
import restored_brugg_de_csmax_audit as audit


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "results" / "260920_mass_independent_nominal_anode"
RESULTS = ROOT / "results" / "260921_nominal_radius_current_configuration"
RESULTS.mkdir(parents=True, exist_ok=True)

NAME = "Ai2020 nominal equilibrium + harvested hysteresis"
RN_UM = 5.0
RP_UM = 3.0
LN_UM = 76.5
LP_UM = 68.0
# Simulation baselines, not equally trusted physical measurements.
# The anode half-cell Rct route was rejected on 2026-09-25 because the harvested
# anode was not reproducible and retained visible residue after extended washing.
# KN_PREF therefore remains only as a legacy comparison value until later fitting.
KN_PREF = 7.40e-7
# KP_PREF is the legacy default. Cathode Rct re-evaluation produced provisional
# alternatives 4.18e-7 (physical c_s,max representation) and 4.32e-7 (the same
# Rct re-expressed with the current model's effective c_s,max); see
# results/260924_kp_eis_recalculation/README_KO.md.
KP_PREF = 3.12e-7
AI2020_MREF = 1.0e-11 * 96485.33212
BRUGG_N = 2.914
BRUGG_P = 1.83
BRUGG_S = 1.5


def apply_provisional_thickness_preserve_capacity(model):
    """Apply teardown coating thicknesses while preserving fitted Qn and Qp."""
    changed = dict(model)
    changed["params"] = model["params"].copy()
    old_ln = float(model["params"]["Negative electrode thickness [m]"])
    old_lp = float(model["params"]["Positive electrode thickness [m]"])
    old_csn = float(model["params"]["Maximum concentration in negative electrode [mol.m-3]"])
    old_csp = float(model["params"]["Maximum concentration in positive electrode [mol.m-3]"])
    new_ln = LN_UM * 1e-6
    new_lp = LP_UM * 1e-6
    new_csn = old_csn * old_ln / new_ln
    new_csp = old_csp * old_lp / new_lp
    changed["params"].update(
        {
            "Negative electrode thickness [m]": new_ln,
            "Positive electrode thickness [m]": new_lp,
            "Maximum concentration in negative electrode [mol.m-3]": new_csn,
            "Maximum concentration in positive electrode [mol.m-3]": new_csp,
        },
        check_already_exists=False,
    )
    changed["csn_max"] = new_csn
    changed["csp_max"] = new_csp
    return changed


def apply_current_kinetics_and_bruggeman(model):
    """Apply the selected experimental kinetics and grounded Bruggeman candidate."""
    changed = dict(model)
    changed["params"] = model["params"].copy()
    raw = pybamm.ParameterValues("Ai2020")
    changed["params"].update(
        {
            "Negative electrode exchange-current density [A.m-2]": ga.scale_parameter_function(
                raw["Negative electrode exchange-current density [A.m-2]"], KN_PREF / AI2020_MREF
            ),
            "Positive electrode exchange-current density [A.m-2]": ga.scale_parameter_function(
                raw["Positive electrode exchange-current density [A.m-2]"], KP_PREF / AI2020_MREF
            ),
            "Negative electrode Bruggeman coefficient (electrolyte)": BRUGG_N,
            "Positive electrode Bruggeman coefficient (electrolyte)": BRUGG_P,
            "Separator Bruggeman coefficient (electrolyte)": BRUGG_S,
        },
        check_already_exists=False,
    )
    return changed


def main():
    bundle = ga.load_legacy_bundle(ga.get_paths(ROOT))
    base = audit.build_base(bundle, 1e-4, audit.CSN_CORRECTED)
    _, _, hybrid, cathode = rs.candidate_details(bundle, base)
    experiment = omc.load_old_dynamic_data(base["dynamic"], base["q_meas"])
    metrics = pd.read_csv(SOURCE / "mass_independent_ocp_candidate_metrics.csv")
    row = metrics[metrics.Candidate == NAME].iloc[0]
    stage = {key: float(row[key]) for key in ("x0", "x100", "y100", "y0")}
    effective = mit.model_for_effective_capacities(base, row)
    model = audit.set_dynamic_parameters(effective, RN_UM, RP_UM, 1e-4)
    model = apply_provisional_thickness_preserve_capacity(model)
    model = apply_current_kinetics_and_bruggeman(model)
    anode = common.scaled_detail(hybrid, 0.50)
    cathode_dynamic = common.scaled_detail(cathode, 0.50)

    validation_rows = []
    simulations = {}
    for rate in ga.RATES:
        for charge in (True, False):
            direction = "Charge" if charge else "Discharge"
            print(f"{rate:g}C {direction}", flush=True)
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
            simulations[(rate, charge)] = sim
            result = omc.old_time_metrics(experiment[(rate, charge)], sim)
            capacity_error = result["end_time_error_min"] * rate * 2.28 / 60.0
            validation_rows.append(
                {
                    "C_rate": rate,
                    "Direction": direction,
                    **result,
                    "capacity_error_Ah": capacity_error,
                    "capacity_error_mAh": capacity_error * 1000.0,
                }
            )
    validation = pd.DataFrame(validation_rows)
    validation.to_csv(RESULTS / "nominal_radius_branch_metrics.csv", index=False, encoding="utf-8-sig")

    fig, axes = plt.subplots(2, 3, figsize=(16, 8.5), constrained_layout=True)
    for col, rate in enumerate(ga.RATES):
        for row_index, charge in enumerate((True, False)):
            ax = axes[row_index, col]
            exp = experiment[(rate, charge)]
            sim = simulations[(rate, charge)]
            branch = validation[(validation.C_rate == rate) & (validation.Direction == ("Charge" if charge else "Discharge"))].iloc[0]
            ax.plot(exp["t_min"], exp["V"], color="black", lw=2.4, label="Experiment")
            ax.plot(sim["t_min"], sim["V"], color="#ff7f0e", lw=2.0, label="Model")
            ax.set_title(
                f"{rate:g}C {'Charge' if charge else 'Discharge'}\n"
                f"MAE(10-70%)={branch.MAE_SOC10_70_mV:.1f} mV, capacity error={branch.capacity_error_mAh:+.1f} mAh"
            )
            ax.set_xlabel("Time (min)")
            ax.set_ylabel("Voltage (V)")
            ax.grid(alpha=0.2)
    axes[0, 0].legend()
    fig.suptitle(
        "Current configuration | experimental kn/kp + provisional Ln/Lp"
    )
    fig.savefig(RESULTS / "nominal_radius_charge_discharge_curves.png", dpi=220)
    plt.close(fig)

    report = {
        "configuration": {
            "OCP": NAME,
            "negative_OCP_charge_branch": "lithiation",
            "negative_OCP_discharge_branch": "delithiation",
            "positive_OCP_charge_branch": "delithiation",
            "positive_OCP_discharge_branch": "lithiation",
            "Rn_um": RN_UM,
            "Rp_um": RP_UM,
            "Ln_um_provisional": LN_UM,
            "Lp_um_provisional": LP_UM,
            "thickness_assumption": "(double-sided total - assumed collector thickness) / 2",
            "assumed_Cu_collector_um": 10.0,
            "assumed_Al_collector_um": 15.0,
            "Dsn_m2_s": rs.DSN,
            "Dsp_m2_s": rs.DSP,
            "kn_prefactor": KN_PREF,
            "kp_prefactor": KP_PREF,
            "kn_multiplier_vs_Ai2020": KN_PREF / AI2020_MREF,
            "kp_multiplier_vs_Ai2020": KP_PREF / AI2020_MREF,
            "Bruggeman_n_p_s": [BRUGG_N, BRUGG_P, BRUGG_S],
            "De_multiplier_vs_raw_Ai2020_function": 1e-4,
            "Stage0_csn_max_mol_m3": audit.CSN_CORRECTED,
            "effective_csn_max_mol_m3": float(model["csn_max"]),
            "effective_csp_max_mol_m3": float(model["csp_max"]),
            "Qn_Ah": float(row.Q_n_Ah),
            "Qp_Ah": float(row.Q_p_Ah),
            **stage,
        },
        "qOCV": {
            "MAE_mV": float(row.qOCV_MAE_2_98_mV),
            "RMSE_mV": float(row.qOCV_RMSE_2_98_mV),
            "max_endpoint_abs_mV": float(row.qOCV_max_endpoint_abs_mV),
        },
        "dynamic_summary": {
            "mean_MAE_SOC10_70_mV": float(validation.MAE_SOC10_70_mV.mean()),
            "charge_mean_MAE_SOC10_70_mV": float(validation[validation.Direction == "Charge"].MAE_SOC10_70_mV.mean()),
            "discharge_mean_MAE_SOC10_70_mV": float(validation[validation.Direction == "Discharge"].MAE_SOC10_70_mV.mean()),
            "mean_abs_capacity_error_mAh": float(validation.capacity_error_mAh.abs().mean()),
            "max_abs_capacity_error_mAh": float(validation.capacity_error_mAh.abs().max()),
        },
        "branches": validation.to_dict(orient="records"),
    }
    (RESULTS / "nominal_radius_current_configuration.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(validation[["C_rate", "Direction", "MAE_SOC10_70_mV", "capacity_error_mAh"]].to_string(index=False))
    print(json.dumps(report["dynamic_summary"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
