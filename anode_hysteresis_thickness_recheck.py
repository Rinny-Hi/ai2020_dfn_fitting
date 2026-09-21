"""Recheck harvested-anode hysteresis after applying measured collector thicknesses."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

import c50_ocp_source_combination_study as c50
import c50_selected_ocp_dynamic_validation as validation
import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import nominal_radius_current_configuration as current
import ocp_rn_sweep_fixed_transport as radius
import priority_initial_state_protocol_recheck as priority


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260921_anode_hysteresis_thickness_recheck"
RESULTS.mkdir(parents=True, exist_ok=True)
SCALES = (0.0, 0.25, 0.50, 0.75, 1.0)


def make_hybrid(nominal: dict, measured: dict, scale: float) -> dict:
    grid = np.asarray(nominal["grid"], dtype=float)
    measured_grid = np.asarray(measured["grid"], dtype=float)
    equilibrium = np.asarray(nominal["equilibrium"], dtype=float)
    result = dict(nominal)
    for branch in ("lithiation", "delithiation"):
        offset = np.asarray(measured[branch]) - np.asarray(measured["equilibrium"])
        result[branch] = equilibrium + scale * np.interp(grid, measured_grid, offset)
    return result


def main() -> None:
    model, stage, nominal, cathode, experiment, selected, qcell = validation.selected_inputs()
    bundle = ga.load_legacy_bundle(ga.get_paths(ROOT))
    measured = ehq.build_anode_detail(bundle, model)

    target_soc = np.linspace(0.0, 1.0, 1001)
    target = c50.fullcell_c50_target(target_soc)
    mask = (target_soc >= 0.10) & (target_soc <= 0.90)

    dynamic_rows = []
    low_rate_rows = []
    for scale in SCALES:
        anode = nominal if scale == 0 else make_hybrid(nominal, measured, scale)
        negative_hysteresis = scale > 0

        equilibrium = c50.predict(
            np.array([stage["x0"], stage["x100"], stage["y100"], stage["y0"]]),
            target_soc,
            nominal,
            cathode,
            "equilibrium",
        )
        charge = c50.predict(
            np.array([stage["x0"], stage["x100"], stage["y100"], stage["y0"]]),
            target_soc,
            anode,
            cathode,
            "charge",
        )
        discharge = c50.predict(
            np.array([stage["x0"], stage["x100"], stage["y100"], stage["y0"]]),
            target_soc,
            anode,
            cathode,
            "discharge",
        )
        low_rate_rows.append(
            {
                "anode_hysteresis_scale": scale,
                "qOCV_MAE_10_90_mV": float(np.mean(np.abs(equilibrium[mask] - target["qocv"][mask])) * 1000),
                "C50_charge_MAE_10_90_mV": float(np.mean(np.abs(charge[mask] - target["charge"][mask])) * 1000),
                "C50_discharge_MAE_10_90_mV": float(np.mean(np.abs(discharge[mask] - target["discharge"][mask])) * 1000),
            }
        )

        for rate in ga.RATES:
            for charge_flag in (True, False):
                direction = "Charge" if charge_flag else "Discharge"
                print(f"scale={scale:.2f} | {rate:g}C {direction}", flush=True)
                sim = ehq.run_dfn(
                    rate,
                    charge_flag,
                    stage,
                    model,
                    anode,
                    cathode,
                    negative_hysteresis,
                    True,
                )
                metrics = priority.curve_metrics(
                    experiment[(rate, charge_flag)], sim, qcell, charge_flag
                )
                dynamic_rows.append(
                    {
                        "anode_hysteresis_scale": scale,
                        "C_rate": rate,
                        "Direction": direction,
                        **metrics,
                        "Capacity_error_mAh": 1000
                        * (metrics["Q_end_model_Ah"] - metrics["Q_end_exp_Ah"]),
                    }
                )

    low_rate = pd.DataFrame(low_rate_rows)
    low_rate["C50_branch_mean_MAE_10_90_mV"] = 0.5 * (
        low_rate["C50_charge_MAE_10_90_mV"]
        + low_rate["C50_discharge_MAE_10_90_mV"]
    )
    dynamic = pd.DataFrame(dynamic_rows)
    summary = dynamic.groupby(["anode_hysteresis_scale", "Direction"], as_index=False).agg(
        mean_center_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
        mean_center_MAE_mV=("Center10_70_MAE_mV", "mean"),
        mean_abs_capacity_error_mAh=("Capacity_error_mAh", lambda x: float(np.mean(np.abs(x)))),
    )
    overall = dynamic.groupby("anode_hysteresis_scale", as_index=False).agg(
        mean_center_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
        mean_center_MAE_mV=("Center10_70_MAE_mV", "mean"),
        mean_abs_capacity_error_mAh=("Capacity_error_mAh", lambda x: float(np.mean(np.abs(x)))),
    )

    low_rate.to_csv(RESULTS / "c50_hysteresis_scale_metrics.csv", index=False, encoding="utf-8-sig")
    dynamic.to_csv(RESULTS / "dynamic_hysteresis_scale_metrics.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(RESULTS / "dynamic_hysteresis_scale_summary.csv", index=False, encoding="utf-8-sig")
    overall.to_csv(RESULTS / "dynamic_hysteresis_scale_overall.csv", index=False, encoding="utf-8-sig")

    params = model["params"]
    report = {
        "applied_coating_thickness_um": {
            "Ln": current.LN_UM,
            "Lp": current.LP_UM,
        },
        "effective_csmax_mol_m3": {
            "negative": float(params["Maximum concentration in negative electrode [mol.m-3]"]),
            "positive": float(params["Maximum concentration in positive electrode [mol.m-3]"]),
        },
        "low_rate": low_rate.to_dict(orient="records"),
        "dynamic_by_direction": summary.to_dict(orient="records"),
        "dynamic_overall": overall.to_dict(orient="records"),
    }
    (RESULTS / "decision_data.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("\nC50")
    print(low_rate.to_string(index=False))
    print("\nDYNAMIC BY DIRECTION")
    print(summary.to_string(index=False))
    print("\nDYNAMIC OVERALL")
    print(overall.to_string(index=False))


if __name__ == "__main__":
    main()
