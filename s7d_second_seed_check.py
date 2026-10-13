"""One independent PSO seed for the four-parameter S7d stability check."""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd

import cumulative_six_stage_model_comparison as cumulative
import electrolyte_transport_set_comparison as transport
import s7_fit_subset_comparison as comparison
import s7_pso_voltage_fitting as s7


OUT = Path(__file__).resolve().parent / "results" / "261013_s7_subset_fit_comparison"


def main() -> None:
    _, _, _, _, experiment, qcell = transport.current_p0()
    entry = cumulative.build_stages(qcell)["S6"]
    grids, targets = s7.make_experimental_grids(experiment)
    cache = {}

    def evaluate(values):
        key = tuple(round(float(values[name]), 18) for name in comparison.ALL_PARAMETERS)
        if key in cache:
            return cache[key]
        started = time.perf_counter()
        try:
            model = comparison.apply_values(entry["model"], values)
            runs = s7.run_model(model, entry["stage"], entry["anode"], entry["cathode"])
            residual, feasibility = s7.residual_vector(runs, grids, targets)
            status = "ok"
        except Exception as exc:
            runs = None
            residual = np.full(6 * s7.N_GRID, 1000.0)
            feasibility = {"Minimum_time_margin_min": np.nan, "Early_cutoff_branches": 6}
            status = f"failed:{type(exc).__name__}"
        result = {
            "runs": runs,
            "residual": residual,
            "Objective_RMSE_mV": float(np.sqrt(np.mean(residual**2))),
            "Objective_MAE_mV": float(np.mean(np.abs(residual))),
            "Status": status,
            "Elapsed_s": time.perf_counter() - started,
            **feasibility,
        }
        cache[key] = result
        return result

    values, result, history = comparison.pso(
        "S7d_seed1409",
        comparison.SUBSETS["S7d_bn_kp_kn_bp"],
        evaluate,
        1409,
    )
    pd.DataFrame([{**values, "Objective_RMSE_mV": result["Objective_RMSE_mV"], "Objective_MAE_mV": result["Objective_MAE_mV"]}]).to_csv(
        OUT / "s7d_seed1409_result.csv", index=False, encoding="utf-8-sig"
    )
    history.to_csv(OUT / "s7d_seed1409_history.csv", index=False, encoding="utf-8-sig")
    print(values, result["Objective_RMSE_mV"], result["Objective_MAE_mV"])


if __name__ == "__main__":
    main()
