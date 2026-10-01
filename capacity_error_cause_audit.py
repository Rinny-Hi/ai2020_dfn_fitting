"""Ablate actual test current and uniform electrode-capacity scaling."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import c50_selected_ocp_dynamic_validation as selected
import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import priority_initial_state_protocol_recheck as priority


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260921_capacity_error_cause_audit"
RESULTS.mkdir(parents=True, exist_ok=True)
NOMINAL_CAPACITY_AH = 2.28
ACTUAL_MEAN_CURRENT_A = {
    (0.5, True): 1.134967,
    (0.5, False): 1.130667,
    (1.0, True): 2.258038,
    (1.0, False): 2.242548,
    (2.0, True): 4.451286,
    (2.0, False): 4.412905,
}


def scale_electrode_capacity(
    model: dict, negative_scale: float, positive_scale: float
) -> dict:
    changed = dict(model)
    changed["params"] = model["params"].copy()
    changed["csn_max"] = float(model["csn_max"]) * negative_scale
    changed["csp_max"] = float(model["csp_max"]) * positive_scale
    changed["params"].update(
        {
            "Maximum concentration in negative electrode [mol.m-3]": changed["csn_max"],
            "Maximum concentration in positive electrode [mol.m-3]": changed["csp_max"],
        },
        check_already_exists=False,
    )
    return changed


def main() -> None:
    model, stage, anode, cathode, experiment, _, qcell = selected.selected_inputs()
    nominal_csn_scale = 28700.0 / float(model["csn_max"])
    nominal_csp_scale = 49943.0 / float(model["csp_max"])
    scenarios = [
        ("baseline", 1.00, 1.00, False),
        ("actual current", 1.00, 1.00, True),
        ("both capacity +2%", 1.02, 1.02, False),
        ("both capacity +4%", 1.04, 1.04, False),
        ("negative capacity +4%", 1.04, 1.00, False),
        ("positive capacity +4%", 1.00, 1.04, False),
        ("both capacity +4% + actual current", 1.04, 1.04, True),
        ("Ai2020 nominal csn only", nominal_csn_scale, 1.00, False),
        ("Ai2020 nominal csp only", 1.00, nominal_csp_scale, False),
        (
            "Ai2020 nominal csn+csp",
            nominal_csn_scale,
            nominal_csp_scale,
            False,
        ),
    ]
    rows = []
    for name, negative_scale, positive_scale, use_actual_current in scenarios:
        local_model = scale_electrode_capacity(model, negative_scale, positive_scale)
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                simulation_rate = (
                    ACTUAL_MEAN_CURRENT_A[(rate, charge)] / NOMINAL_CAPACITY_AH
                    if use_actual_current
                    else rate
                )
                print(f"{name} | {rate:g}C {direction}", flush=True)
                sim = ehq.run_dfn(
                    simulation_rate,
                    charge,
                    stage,
                    local_model,
                    anode,
                    cathode,
                    False,
                    True,
                )
                metric = priority.curve_metrics(
                    experiment[(rate, charge)], sim, qcell, charge
                )
                rows.append(
                    {
                        "Scenario": name,
                        "negative_capacity_scale": negative_scale,
                        "positive_capacity_scale": positive_scale,
                        "C_rate_label": rate,
                        "simulation_rate_C": simulation_rate,
                        "Direction": direction,
                        **metric,
                        "Capacity_error_mAh": 1000
                        * (metric["Q_end_model_Ah"] - metric["Q_end_exp_Ah"]),
                    }
                )
    detail = pd.DataFrame(rows)
    summary = detail.groupby(["Scenario", "Direction"], as_index=False, sort=False).agg(
        mean_center_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
        mean_capacity_error_mAh=("Capacity_error_mAh", "mean"),
        mean_abs_capacity_error_mAh=("Capacity_error_mAh", lambda x: float(np.mean(np.abs(x)))),
    )
    detail.to_csv(RESULTS / "capacity_cause_detail.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(RESULTS / "capacity_cause_summary.csv", index=False, encoding="utf-8-sig")
    print("\nSUMMARY")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
