"""Four-parameter S7d fit initialized with Ai2020 kn=kp."""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd

import cumulative_six_stage_model_comparison as cumulative
import electrolyte_transport_set_comparison as transport
import kinetics_original_reversion_ablation as kinetics
import s7_fit_subset_comparison as comparison
import s7_pso_voltage_fitting as s7


OUT = Path(__file__).resolve().parent / "results" / "261013_s7_subset_fit_comparison"
LABEL = "S7e_Ai2020_kn_kp_start"
NAMES = comparison.SUBSETS["S7d_bn_kp_kn_bp"]
AI_PREF = kinetics.AI2020_PREF
START = {**comparison.BASE, "kn": AI_PREF, "kp": AI_PREF}


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

    start_vector = np.array([comparison.encode(name, START[name]) for name in NAMES])
    start_result = evaluate(START)
    rng = np.random.default_rng(2020)
    positions = np.empty((comparison.N_PARTICLES, len(NAMES)))
    positions[0] = start_vector
    positions[1:] = np.clip(
        start_vector + rng.normal(0.0, 0.18, size=(comparison.N_PARTICLES - 1, len(NAMES))),
        0.0,
        1.0,
    )
    velocities = rng.normal(0.0, 0.04, size=positions.shape)
    personal_positions = positions.copy()
    personal_scores = np.full(comparison.N_PARTICLES, np.inf)
    global_position = start_vector.copy()
    global_score = np.inf
    history = []

    for iteration in range(comparison.N_ITERATIONS):
        for particle in range(comparison.N_PARTICLES):
            values = comparison.values_from_vector(positions[particle], NAMES)
            result = evaluate(values)
            score = result["Objective_RMSE_mV"]
            if score < personal_scores[particle]:
                personal_scores[particle] = score
                personal_positions[particle] = positions[particle].copy()
            if score < global_score:
                global_score = score
                global_position = positions[particle].copy()
            history.append({"Scenario": LABEL, "Iteration": iteration + 1, "Particle": particle + 1, "Objective_RMSE_mV": score, **values})
        print(f"{LABEL} iteration {iteration + 1}/{comparison.N_ITERATIONS}: {global_score:.4f} mV", flush=True)
        r1 = rng.random(size=positions.shape)
        r2 = rng.random(size=positions.shape)
        velocities = 0.72 * velocities + 1.49 * r1 * (personal_positions - positions) + 1.49 * r2 * (global_position - positions)
        positions = np.clip(positions + np.clip(velocities, -0.22, 0.22), 0.0, 1.0)

    values = comparison.values_from_vector(global_position, NAMES)
    result = evaluate(values)
    detail, anode = s7.scenario_metrics(LABEL, result, experiment, qcell)
    detail_frame = pd.DataFrame(detail)
    summary = s7.aggregate_metrics(detail_frame)
    result_row = {
        "Scenario": LABEL,
        "Start_kn": START["kn"],
        "Start_kp": START["kp"],
        "Start_objective_RMSE_mV": start_result["Objective_RMSE_mV"],
        **values,
        "Objective_RMSE_mV": result["Objective_RMSE_mV"],
        "Objective_MAE_mV": result["Objective_MAE_mV"],
    }
    pd.DataFrame([result_row]).to_csv(OUT / "s7e_ai2020_start_parameters.csv", index=False, encoding="utf-8-sig")
    detail_frame.to_csv(OUT / "s7e_ai2020_start_branch_metrics.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(anode).to_csv(OUT / "s7e_ai2020_start_anode.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(OUT / "s7e_ai2020_start_summary.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(history).to_csv(OUT / "s7e_ai2020_start_history.csv", index=False, encoding="utf-8-sig")
    print("START", START, start_result["Objective_RMSE_mV"])
    print("FINAL", result_row)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
