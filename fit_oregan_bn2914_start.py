"""Fit k_n and b_n for ORegan2022 starting from the Ai2020 b_n=2.914.

This is a controlled start-point check against the existing S7a calculation.
It reuses the same data, bounds, fixed-time voltage objective, early-cutoff
penalty, PSO settings, and post-fit diagnostics as s7_pso_voltage_fitting.py.
"""

from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np
import pandas as pd

import cumulative_six_stage_model_comparison as cumulative
import electrolyte_transport_set_comparison as transport
import s7_pso_voltage_fitting as s7


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "261020_oregan_bn2914_start_fit"
START = {
    "Dsn": 4.1707741459682473e-14,
    "kn": 9.648533212e-7,
    "brugg_n": 2.914,
}
FITTED = ("kn", "brugg_n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    # Reuse all numerical definitions from the existing S7 workflow, changing
    # only the official starting point for b_n.
    s7.BASE = dict(START)

    _, _, _, _, experiment, qcell = transport.current_p0()
    s6_entry = cumulative.build_stages(qcell)["S6"]
    s6_model = s6_entry["model"]
    stage = s6_entry["stage"]
    anode = s6_entry["anode"]
    cathode = s6_entry["cathode"]
    grids, targets = s7.make_experimental_grids(experiment)

    cache: dict[tuple[float, float, float], dict] = {}
    evaluations: list[dict] = []

    def evaluate(values: dict[str, float]) -> dict:
        key = tuple(round(float(values[name]), 18) for name in s7.PARAMETERS)
        if key in cache:
            return cache[key]
        started = time.perf_counter()
        status = "ok"
        try:
            model = s7.apply_values(s6_model, values)
            runs = s7.run_model(model, stage, anode, cathode)
            residual, feasibility = s7.residual_vector(runs, grids, targets)
        except Exception as exc:
            runs = None
            residual = np.full(6 * s7.N_GRID, 1000.0)
            feasibility = {
                "Minimum_time_margin_min": np.nan,
                "Early_cutoff_branches": 6,
            }
            status = f"failed:{type(exc).__name__}"
        result = {
            "Objective_RMSE_mV": float(np.sqrt(np.mean(residual**2))),
            "Objective_MAE_mV": float(np.mean(np.abs(residual))),
            "residual": residual,
            "runs": runs,
            "Status": status,
            **feasibility,
        }
        cache[key] = result
        evaluations.append(
            {
                "Evaluation": len(evaluations) + 1,
                **values,
                "Objective_RMSE_mV": result["Objective_RMSE_mV"],
                "Objective_MAE_mV": result["Objective_MAE_mV"],
                "Minimum_time_margin_min": result["Minimum_time_margin_min"],
                "Early_cutoff_branches": result["Early_cutoff_branches"],
                "Status": status,
                "Elapsed_s": time.perf_counter() - started,
            }
        )
        if len(evaluations) == 1 or len(evaluations) % 10 == 0:
            print(
                f"evaluation {len(evaluations)}: "
                f"{result['Objective_RMSE_mV']:.3f} mV",
                flush=True,
            )
            pd.DataFrame(evaluations).to_csv(
                OUT / "all_model_evaluations.csv",
                index=False,
                encoding="utf-8-sig",
            )
        return result

    print("ORegan2022 baseline with b_n=2.914", flush=True)
    start_result = evaluate(START)

    print("PSO fit: k_n + b_n", flush=True)
    fitted_result, history = s7.pso_fit(
        "ORegan b_n=2.914 start",
        FITTED,
        evaluate,
        informed=True,
        seed=261013,
    )

    results = {
        "Start b_n=2.914": {**START, **start_result},
        "Fitted from b_n=2.914": fitted_result,
    }
    detail_rows: list[dict] = []
    anode_rows: list[dict] = []
    for scenario, result in results.items():
        rows, potentials = s7.scenario_metrics(
            scenario, result, experiment, qcell
        )
        detail_rows.extend(rows)
        anode_rows.extend(potentials)

    detail = pd.DataFrame(detail_rows)
    summary = s7.aggregate_metrics(detail)
    parameters = pd.DataFrame(
        [
            {
                "Scenario": scenario,
                "Dsn_m2_s": result["Dsn"],
                "kn_prefactor": result["kn"],
                "brugg_n": result["brugg_n"],
                "Objective_RMSE_mV": result["Objective_RMSE_mV"],
                "Objective_MAE_mV": result["Objective_MAE_mV"],
                "Minimum_time_margin_min": result["Minimum_time_margin_min"],
                "Early_cutoff_branches": result["Early_cutoff_branches"],
            }
            for scenario, result in results.items()
        ]
    )

    detail.to_csv(OUT / "branch_metrics.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(OUT / "summary.csv", index=False, encoding="utf-8-sig")
    parameters.to_csv(OUT / "parameters.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(anode_rows).to_csv(
        OUT / "anode_potential.csv", index=False, encoding="utf-8-sig"
    )
    history.to_csv(OUT / "pso_history.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(evaluations).to_csv(
        OUT / "all_model_evaluations.csv", index=False, encoding="utf-8-sig"
    )

    manifest = {
        "electrolyte": "ORegan2022",
        "official_start": START,
        "fitted_parameters": FITTED,
        "bounds": s7.BOUNDS,
        "objective": (
            "fixed physical-time terminal-voltage RMSE over six Cell 6-6 "
            "charge/discharge branches"
        ),
        "capacity_in_objective": False,
        "anode_potential_in_objective": False,
        "pso_particles": s7.N_PARTICLES,
        "pso_iterations": s7.N_ITERATIONS,
        "seed": 261013,
    }
    (OUT / "fit_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print("\nParameters", flush=True)
    print(parameters.to_string(index=False), flush=True)
    print("\nSummary", flush=True)
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
