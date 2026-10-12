"""Compare electrolyte transport sets on the current 2--98% pre-fit P0.

Only electrolyte bulk-transport inputs are replaced.  Electrode OCP, geometry,
solid diffusion, kinetics, stoichiometry window, and active inventory are held
at the current Enertech P0 values.  The Ai2020 diffusivity unit correction is
already present in the baseline; no 1e-4 multiplier is applied to the other
PyBaMM parameter sets, whose diffusivities are already in SI units.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pybamm

import c50_selected_ocp_dynamic_validation as selected
import comprehensive_prefit_cross_cohort as cross
import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import paper_teardown_geometry as geometry
import prefit_physical_anchor_decision as p0
import priority_initial_state_protocol_recheck as priority


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "261012_electrolyte_transport_set_comparison"

SETS = ("Ai2020 corrected", "Chen2020", "OKane2022", "ORegan2022")
TRANSPORT_KEYS = (
    "Electrolyte diffusivity [m2.s-1]",
    "Electrolyte conductivity [S.m-1]",
    "Cation transference number",
    "Thermodynamic factor",
    "Initial concentration in electrolyte [mol.m-3]",
)

# Cell 6-6 July BoL protocol values recorded in the executed final notebook.
PROTOCOL = {
    (0.5, True): (3.0601, 1.139948),
    (0.5, False): (4.1783, 1.139974),
    (1.0, True): (3.2073, 2.279965),
    (1.0, False): (4.1780, 2.279931),
    (2.0, True): (3.3263, 4.559894),
    (2.0, False): (4.1787, 4.560041),
}


def current_p0():
    model0, stage, anode, cathode, experiment, selected_ocp, qcell = selected.selected_inputs()
    scenario = p0.Scenario(
        name="P0 anchored reference",
        ln_um=76.5,
        lp_um=68.0,
        rn_um=3.7542,
        rp_um=4.2387,
        kn=9.648533212e-7,
        kp=5.899519906708763e-7,
        dsn="measured",
        dsp="measured",
        positive_hysteresis=False,
        initialization="rest",
    )
    model = geometry.apply_overlap_area(p0.build_model(model0, scenario))
    model["params"].update(
        {
            "Negative particle diffusivity [m2.s-1]": 4.1707741459682473e-14,
            "Positive particle diffusivity [m2.s-1]": 6.888175323807748e-14,
            "Negative electrode Bruggeman coefficient (electrolyte)": 2.914,
            "Positive electrode Bruggeman coefficient (electrolyte)": 1.83,
            "Separator Bruggeman coefficient (electrolyte)": 1.5,
            "Negative electrode active material volume fraction": 0.61,
            "Positive electrode active material volume fraction": 0.62,
            "Maximum concentration in negative electrode [mol.m-3]": 25182.812491955036,
            "Maximum concentration in positive electrode [mol.m-3]": 49140.111895837676,
        },
        check_already_exists=False,
    )
    model["csn_max"] = 25182.812491955036
    model["csp_max"] = 49140.111895837676
    model["qn_cell"] = float(selected_ocp.Qn_Ah)
    model["qp_cell"] = float(selected_ocp.Qp_Ah)
    model["q_meas"] = float(qcell)
    return model, stage, anode, cathode, experiment, float(qcell)


def with_transport_set(model, name):
    changed = dict(model)
    changed["params"] = model["params"].copy()
    if name == "Ai2020 corrected":
        return changed
    source = pybamm.ParameterValues(name)
    updates = {key: source[key] for key in TRANSPORT_KEYS if key in source}
    changed["params"].update(updates, check_already_exists=False)
    return changed


def scalar_value(parameter_values, value, concentration=1000.0, temperature=298.15):
    if not callable(value):
        return float(value)
    for args in ((concentration, temperature), (concentration,), (temperature,)):
        try:
            result = value(*args)
            if hasattr(result, "evaluate"):
                result = parameter_values.evaluate(result)
            return float(np.asarray(result).squeeze())
        except (TypeError, ValueError):
            continue
    return np.nan


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    baseline, stage, anode, cathode, experiment, qcell = current_p0()
    detail_rows = []
    property_rows = []

    for set_name in SETS:
        model = with_transport_set(baseline, set_name)
        params = model["params"]
        property_rows.append(
            {
                "Electrolyte_set": set_name,
                "De_at_1M_298K_m2_s": scalar_value(params, params[TRANSPORT_KEYS[0]]),
                "Conductivity_at_1M_298K_S_m": scalar_value(params, params[TRANSPORT_KEYS[1]]),
                "t_plus_at_1M_298K": scalar_value(params, params[TRANSPORT_KEYS[2]]),
                "Thermodynamic_factor_at_1M_298K": scalar_value(params, params[TRANSPORT_KEYS[3]]),
                "Initial_electrolyte_concentration_mol_m3": float(params[TRANSPORT_KEYS[4]]),
            }
        )
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                rest_v, measured_current = PROTOCOL[(rate, charge)]
                z0 = priority.monotone_voltage_inverse(
                    model["v_qocv"], model["soc"], rest_v
                )
                local_stage = priority.stage_at_soc(stage, z0, charge)
                simulation_rate = measured_current / cross.NOMINAL_CAPACITY_AH
                print(f"{set_name} | {rate:g}C {direction}", flush=True)
                try:
                    simulation = ehq.run_dfn(
                        simulation_rate,
                        charge,
                        local_stage,
                        model,
                        anode,
                        cathode,
                        False,
                        False,
                    )
                    metrics = priority.curve_metrics(
                        experiment[(rate, charge)], simulation, qcell, charge
                    )
                    detail_rows.append(
                        {
                            "Scenario": set_name,
                            "Trust": "transport sensitivity",
                            "Note": "cell 6-6 integrated dynamic data; no fitting",
                            "Cell": "6-6",
                            "C_rate": rate,
                            "Direction": direction,
                            "Initial_SOC_pct": 100.0 * z0,
                            **metrics,
                            "Capacity_error_mAh": 1000.0
                            * (metrics["Q_end_model_Ah"] - metrics["Q_end_exp_Ah"]),
                            "Status": "ok",
                        }
                    )
                except Exception as exc:
                    detail_rows.append(
                        {
                            "Scenario": set_name,
                            "Trust": "transport sensitivity",
                            "Note": "cell 6-6 integrated dynamic data; no fitting",
                            "Cell": "6-6",
                            "C_rate": rate,
                            "Direction": direction,
                            "Initial_SOC_pct": 100.0 * z0,
                            "Status": f"failed: {type(exc).__name__}: {exc}",
                        }
                    )

    detail = pd.DataFrame(detail_rows)
    properties = pd.DataFrame(property_rows)
    valid = detail[detail["Status"] == "ok"].copy()
    summary = p0.summarize(valid) if not valid.empty else pd.DataFrame()
    branch = p0.branch_summary(valid, list(SETS)) if not valid.empty else pd.DataFrame()

    detail.to_csv(OUT / "detail.csv", index=False, encoding="utf-8-sig")
    properties.to_csv(OUT / "electrolyte_properties.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(OUT / "summary.csv", index=False, encoding="utf-8-sig")
    branch.to_csv(OUT / "branch_summary.csv", index=False, encoding="utf-8-sig")
    print("\n[properties]\n", properties.to_string(index=False))
    print("\n[summary]\n", summary.to_string(index=False))
    print("\n[branch]\n", branch.to_string(index=False))
    failures = detail[detail["Status"] != "ok"]
    if not failures.empty:
        print("\n[failures]\n", failures[["Scenario", "C_rate", "Direction", "Status"]].to_string(index=False))


if __name__ == "__main__":
    main()
