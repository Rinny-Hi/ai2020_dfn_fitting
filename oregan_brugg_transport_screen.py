"""Diagnostic ORegan2022 screen for negative-electrode transport resistance."""

from __future__ import annotations

import numpy as np
import pandas as pd

import comprehensive_prefit_cross_cohort as cross
import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import prefit_physical_anchor_decision as p0
import priority_initial_state_protocol_recheck as priority

import electrolyte_transport_set_comparison as comparison


def main():
    baseline, stage, anode, cathode, experiment, qcell = comparison.current_p0()
    oregan = comparison.with_transport_set(baseline, "ORegan2022")
    rows = []
    for brugg_n in (2.914, 2.5, 2.0, 1.5):
        model = dict(oregan)
        model["params"] = oregan["params"].copy()
        model["params"].update(
            {"Negative electrode Bruggeman coefficient (electrolyte)": brugg_n}
        )
        name = f"ORegan2022 brugg_n={brugg_n:g}"
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                rest_v, measured_current = comparison.PROTOCOL[(rate, charge)]
                z0 = priority.monotone_voltage_inverse(
                    model["v_qocv"], model["soc"], rest_v
                )
                local_stage = priority.stage_at_soc(stage, z0, charge)
                sim = ehq.run_dfn(
                    measured_current / cross.NOMINAL_CAPACITY_AH,
                    charge,
                    local_stage,
                    model,
                    anode,
                    cathode,
                    False,
                    False,
                )
                metrics = priority.curve_metrics(
                    experiment[(rate, charge)], sim, qcell, charge
                )
                rows.append(
                    {
                        "Scenario": name,
                        "Trust": "diagnostic only",
                        "Note": "ORegan2022; only brugg_n changed; no fitting",
                        "Cell": "6-6",
                        "C_rate": rate,
                        "Direction": direction,
                        **metrics,
                        "Capacity_error_mAh": 1000
                        * (metrics["Q_end_model_Ah"] - metrics["Q_end_exp_Ah"]),
                    }
                )
    detail = pd.DataFrame(rows)
    summary = p0.summarize(detail)
    branch = p0.branch_summary(detail, detail["Scenario"].drop_duplicates().tolist())
    out = comparison.OUT
    detail.to_csv(out / "oregan_brugg_screen_detail.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(out / "oregan_brugg_screen_summary.csv", index=False, encoding="utf-8-sig")
    branch.to_csv(out / "oregan_brugg_screen_branch.csv", index=False, encoding="utf-8-sig")
    print("\n[summary]\n", summary.to_string(index=False))
    print("\n[branch]\n", branch.to_string(index=False))


if __name__ == "__main__":
    main()
