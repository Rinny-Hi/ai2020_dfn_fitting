"""Compare bn+kn, bn+kn+bp, and bn+kn+kp+bp voltage fits.

All cases use the S7 fixed-time six-branch voltage objective.  Capacity and
anode potential are validation metrics only.  S7a is evaluated at its saved
optimum; the larger subsets use the same 7-particle, 7-iteration PSO budget.
"""

from __future__ import annotations

from pathlib import Path
import time

import numpy as np
import pandas as pd

import cumulative_six_stage_model_comparison as cumulative
import electrolyte_transport_set_comparison as transport
import kinetics_original_reversion_ablation as kinetics
import s7_pso_voltage_fitting as s7


OUT = Path(__file__).resolve().parent / "results" / "261013_s7_subset_fit_comparison"
ALL_PARAMETERS = ("kn", "kp", "brugg_n", "brugg_p")
BASE = {
    "kn": 9.648533212e-7,
    "kp": 5.899519906708763e-7,
    "brugg_n": 2.5,
    "brugg_p": 1.83,
}
BOUNDS = {
    "kn": (3.0e-7, 3.0e-6),
    "kp": (0.3 * BASE["kp"], 3.0 * BASE["kp"]),
    "brugg_n": (1.5, 3.5),
    "brugg_p": (1.5, 3.5),
}
SUBSETS = {
    "S7a_bn_kn": ("brugg_n", "kn"),
    "S7c_bn_kn_bp": ("brugg_n", "kn", "brugg_p"),
    "S7d_bn_kp_kn_bp": ("brugg_n", "kp", "kn", "brugg_p"),
}
S7A_VALUES = {
    **BASE,
    "kn": 7.746509582521757e-7,
    "brugg_n": 2.4949377341788863,
}
N_PARTICLES = 7
N_ITERATIONS = 7


def encode(name: str, value: float) -> float:
    lo, hi = BOUNDS[name]
    if name in {"kn", "kp"}:
        return float((np.log(value) - np.log(lo)) / (np.log(hi) - np.log(lo)))
    return float((value - lo) / (hi - lo))


def decode(name: str, coordinate: float) -> float:
    lo, hi = BOUNDS[name]
    z = float(np.clip(coordinate, 0.0, 1.0))
    if name in {"kn", "kp"}:
        return float(np.exp(np.log(lo) + z * (np.log(hi) - np.log(lo))))
    return float(lo + z * (hi - lo))


def values_from_vector(vector: np.ndarray, names: tuple[str, ...]) -> dict[str, float]:
    values = dict(BASE)
    for name, coordinate in zip(names, vector):
        values[name] = decode(name, float(coordinate))
    return values


def apply_values(model: dict, values: dict[str, float]) -> dict:
    changed = kinetics.with_prefactors(model, values["kn"], values["kp"])
    changed = dict(changed)
    changed["params"] = changed["params"].copy()
    changed["params"].update(
        {
            "Negative particle diffusivity [m2.s-1]": 4.1707741459682473e-14,
            "Positive particle diffusivity [m2.s-1]": 6.888175323807748e-14,
            "Negative electrode Bruggeman coefficient (electrolyte)": values["brugg_n"],
            "Positive electrode Bruggeman coefficient (electrolyte)": values["brugg_p"],
            "Separator Bruggeman coefficient (electrolyte)": 1.5,
        },
        check_already_exists=False,
    )
    return changed


def pso(label, names, evaluate, seed):
    rng = np.random.default_rng(seed)
    base = np.array([encode(name, BASE[name]) for name in names])
    positions = np.empty((N_PARTICLES, len(names)))
    positions[0] = base
    positions[1:] = np.clip(base + rng.normal(0.0, 0.18, size=(N_PARTICLES - 1, len(names))), 0.0, 1.0)
    velocities = rng.normal(0.0, 0.04, size=positions.shape)
    personal_positions = positions.copy()
    personal_scores = np.full(N_PARTICLES, np.inf)
    global_position = base.copy()
    global_score = np.inf
    history = []
    for iteration in range(N_ITERATIONS):
        for particle in range(N_PARTICLES):
            values = values_from_vector(positions[particle], names)
            result = evaluate(values)
            score = result["Objective_RMSE_mV"]
            if score < personal_scores[particle]:
                personal_scores[particle] = score
                personal_positions[particle] = positions[particle].copy()
            if score < global_score:
                global_score = score
                global_position = positions[particle].copy()
            history.append({"Scenario": label, "Iteration": iteration + 1, "Particle": particle + 1, "Objective_RMSE_mV": score, **values})
        print(f"{label} iteration {iteration + 1}/{N_ITERATIONS}: {global_score:.4f} mV", flush=True)
        r1 = rng.random(size=positions.shape)
        r2 = rng.random(size=positions.shape)
        velocities = 0.72 * velocities + 1.49 * r1 * (personal_positions - positions) + 1.49 * r2 * (global_position - positions)
        positions = np.clip(positions + np.clip(velocities, -0.22, 0.22), 0.0, 1.0)
    values = values_from_vector(global_position, names)
    return values, evaluate(values), pd.DataFrame(history)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    _, _, _, _, experiment, qcell = transport.current_p0()
    entry = cumulative.build_stages(qcell)["S6"]
    grids, targets = s7.make_experimental_grids(experiment)
    cache = {}
    evaluations = []

    def evaluate(values):
        key = tuple(round(float(values[name]), 18) for name in ALL_PARAMETERS)
        if key in cache:
            return cache[key]
        started = time.perf_counter()
        model = apply_values(entry["model"], values)
        try:
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
            **feasibility,
        }
        cache[key] = result
        evaluations.append({"Evaluation": len(evaluations) + 1, **values, "Objective_RMSE_mV": result["Objective_RMSE_mV"], "Objective_MAE_mV": result["Objective_MAE_mV"], "Status": status, "Elapsed_s": time.perf_counter() - started})
        return result

    results = {}
    histories = []
    results["S7a_bn_kn"] = (S7A_VALUES, evaluate(S7A_VALUES))
    for label, seed in (("S7c_bn_kn_bp", 731), ("S7d_bn_kp_kn_bp", 947)):
        values, result, history = pso(label, SUBSETS[label], evaluate, seed)
        results[label] = (values, result)
        histories.append(history)

    parameter_rows = []
    detail_rows = []
    anode_rows = []
    for label, (values, result) in results.items():
        parameter_rows.append({"Scenario": label, "Fitted_parameters": "+".join(SUBSETS[label]), **values, "Objective_RMSE_mV": result["Objective_RMSE_mV"], "Objective_MAE_mV": result["Objective_MAE_mV"]})
        detail, anode = s7.scenario_metrics(label, result, experiment, qcell)
        detail_rows.extend(detail)
        anode_rows.extend(anode)

    parameters = pd.DataFrame(parameter_rows)
    detail = pd.DataFrame(detail_rows)
    summary = s7.aggregate_metrics(detail).merge(parameters, on="Scenario", how="left")
    parameters.to_csv(OUT / "fitted_parameters.csv", index=False, encoding="utf-8-sig")
    detail.to_csv(OUT / "branch_metrics.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(anode_rows).to_csv(OUT / "anode_potential.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(OUT / "summary.csv", index=False, encoding="utf-8-sig")
    pd.concat(histories, ignore_index=True).to_csv(OUT / "pso_history.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(evaluations).to_csv(OUT / "all_evaluations.csv", index=False, encoding="utf-8-sig")
    print("\nPARAMETERS\n", parameters.to_string(index=False))
    print("\nSUMMARY\n", summary.to_string(index=False))


if __name__ == "__main__":
    main()
