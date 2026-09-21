"""Compare micrometer, Ai2020, and SEM-ratio proxy coating thicknesses."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import c50_selected_ocp_dynamic_validation as selected
import current_thickness_estimate_application as thickness
import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import priority_initial_state_protocol_recheck as priority


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260921_thickness_measurement_method_sensitivity"
RESULTS.mkdir(parents=True, exist_ok=True)

POSITIVE_SEM_RATIO = 64.3 / 58.5
NEGATIVE_SEM_RATIO = 96.0 / 79.5

SCENARIOS = {
    "Enertech micrometer": (74.0, 66.0),
    "Ai2020 nominal": (76.5, 68.0),
    "Total-minus-collector": (73.0, 62.5),
    "Hyundai SEM-ratio proxy": (
        74.0 * NEGATIVE_SEM_RATIO,
        66.0 * POSITIVE_SEM_RATIO,
    ),
}


def main() -> None:
    base, stage, anode, cathode, experiment, _, qcell = selected.selected_inputs()
    rows = []
    for name, (ln_um, lp_um) in SCENARIOS.items():
        model = thickness.set_thickness_preserve_capacity(
            base, ln_um * 1e-6, lp_um * 1e-6
        )
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                print(f"{name} | {rate:g}C {direction}", flush=True)
                sim = ehq.run_dfn(
                    rate, charge, stage, model, anode, cathode, False, True
                )
                metric = priority.curve_metrics(
                    experiment[(rate, charge)], sim, qcell, charge
                )
                rows.append(
                    {
                        "Scenario": name,
                        "Ln_um": ln_um,
                        "Lp_um": lp_um,
                        "C_rate": rate,
                        "Direction": direction,
                        **metric,
                        "Capacity_error_mAh": 1000
                        * (metric["Q_end_model_Ah"] - metric["Q_end_exp_Ah"]),
                    }
                )
    detail = pd.DataFrame(rows)
    summary = detail.groupby(
        ["Scenario", "Ln_um", "Lp_um"], as_index=False, sort=False
    ).agg(
        mean_center_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
        mean_center_MAE_mV=("Center10_70_MAE_mV", "mean"),
        mean_abs_capacity_error_mAh=("Capacity_error_mAh", lambda x: float(np.mean(np.abs(x)))),
        charge_mean_capacity_error_mAh=(
            "Capacity_error_mAh",
            lambda x: float(np.mean(x[detail.loc[x.index, "Direction"] == "Charge"])),
        ),
        discharge_mean_abs_capacity_error_mAh=(
            "Capacity_error_mAh",
            lambda x: float(np.mean(np.abs(x[detail.loc[x.index, "Direction"] == "Discharge"]))),
        ),
    )
    detail.to_csv(RESULTS / "thickness_method_detail.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(RESULTS / "thickness_method_summary.csv", index=False, encoding="utf-8-sig")
    print("\nSUMMARY")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
