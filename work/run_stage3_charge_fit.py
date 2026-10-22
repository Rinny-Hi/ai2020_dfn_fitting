"""Charge-only two-parameter fit, independent discharge validation, and potential audit.

Four branches: two OCP choices × Ai2020-corrected or ORegan2022 electrolyte.
Only the five electrolyte bulk properties are transplanted from ORegan2022.
"""

from __future__ import annotations

import json
import math
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pybamm
from scipy.optimize import least_squares

from run_stage2_two_ocp import cases, kinetic_constants, scope, stage_with_values


OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(exist_ok=True)
TRANSPORT = (
    "Electrolyte diffusivity [m2.s-1]",
    "Electrolyte conductivity [S.m-1]",
    "Cation transference number",
    "Thermodynamic factor",
    "Initial concentration in electrolyte [mol.m-3]",
)
OREGAN = pybamm.ParameterValues("ORegan2022")
BN_KEY = "Negative electrode Bruggeman coefficient (electrolyte)"
KN_BOUNDS = (3.0e-7, 3.0e-6)
BN_BOUNDS = (1.5, 3.5)
CHARGE_RATES = (0.5, 1.0, 2.0)
CONDITIONS = tuple((rate, charge) for rate in CHARGE_RATES for charge in (True, False))
N_TRAIN = 140
SHORTFALL_WEIGHT_MV = 300.0
LOG_PATH = OUT / "Stage3_charge_fit_evaluations.csv"


def model_stage(case, electrolyte, kn, bn):
    values = kinetic_constants(case)
    values["kn"] = float(kn)
    stage = stage_with_values(case, values)
    parameters = stage["parameters"].copy()
    updates = {BN_KEY: float(bn)}
    if electrolyte == "ORegan2022":
        updates.update({key: OREGAN[key] for key in TRANSPORT})
    parameters.update(updates)
    stage["parameters"] = parameters
    return stage


def train_data():
    blocks = {}
    for rate in CHARGE_RATES:
        exp = scope["experimental_curves"][(rate, True)]
        t = np.array(exp["time_min"], dtype=float, copy=True)
        t -= t[0]
        grid = np.linspace(0, 0.98 * t[-1], N_TRAIN)
        blocks[rate] = {"time": grid, "voltage": np.interp(grid, t, exp["voltage_v"]),
                        "experimental_end_min": float(t[-1])}
    return blocks


TRAIN = train_data()


def coordinate_values(z):
    return float(np.exp(z[0])), float(z[1])


def simulate_charge(case, electrolyte, kn, bn):
    stage = model_stage(case, electrolyte, kn, bn)
    return {rate: scope["run_dfn_branch"](stage, rate, True) for rate in CHARGE_RATES}


def residuals_from_runs(runs):
    blocks, penalties = [], []
    for rate in CHARGE_RATES:
        data = TRAIN[rate]
        run = runs[rate]
        t = np.asarray(run["time_min"], float)
        v = np.asarray(run["voltage_v"], float)
        prediction = np.interp(data["time"], t, v)
        residual = 1000.0 * (prediction - data["voltage"])
        shortfall = max(0.0, data["time"][-1] - t[-1]) / data["time"][-1]
        blocks.append(residual)
        penalties.extend([SHORTFALL_WEIGHT_MV * shortfall] * 20)
    return np.concatenate([*blocks, np.asarray(penalties)])


def fit_one(case_name, case, electrolyte):
    initial = kinetic_constants(case)
    initial_kn = initial["kn"]
    seeds = (
        (initial_kn, 2.914),
        (max(KN_BOUNDS[0], min(KN_BOUNDS[1], 1.35 * initial_kn)), 2.1),
        (max(KN_BOUNDS[0], min(KN_BOUNDS[1], 0.75 * initial_kn)), 3.3),
    )
    lower = (math.log(KN_BOUNDS[0]), BN_BOUNDS[0])
    upper = (math.log(KN_BOUNDS[1]), BN_BOUNDS[1])
    evaluations = []
    cache = {}

    def objective(z):
        kn, bn = coordinate_values(z)
        key = (round(kn, 15), round(bn, 10))
        if key not in cache:
            try:
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", category=pybamm.SolverWarning)
                    runs = simulate_charge(case, electrolyte, kn, bn)
                residual = residuals_from_runs(runs)
                status = "ok"
            except Exception as exc:
                residual = np.full(3 * N_TRAIN + 60, 1000.0)
                status = f"failed:{type(exc).__name__}"
            cache[key] = residual
            voltage_rms = float(np.sqrt(np.mean(residual[:3 * N_TRAIN] ** 2)))
            record = {"OCP_case": case_name, "Electrolyte": electrolyte,
                      "Eval": len(evaluations) + 1, "kn": kn, "brugg_n": bn,
                      "Train_voltage_RMSE_mV": voltage_rms,
                      "Objective_RMS_mV": float(np.sqrt(np.mean(residual**2))),
                      "Status": status}
            evaluations.append(record)
            if len(evaluations) % 5 == 0:
                print("EVAL", case_name, electrolyte, len(evaluations),
                      "best", round(min(x["Objective_RMS_mV"] for x in evaluations), 3), flush=True)
                pd.DataFrame(evaluations).to_csv(LOG_PATH, mode="a", header=not LOG_PATH.exists(),
                                                 index=False, encoding="utf-8-sig")
                evaluations.clear()
        return cache[key]

    results = []
    for start, (kn0, bn0) in enumerate(seeds, 1):
        print("FIT", case_name, electrolyte, "start", start, kn0, bn0, flush=True)
        result = least_squares(objective, [math.log(kn0), bn0], bounds=(lower, upper),
                               max_nfev=20, ftol=1e-4, xtol=1e-4, gtol=1e-4,
                               diff_step=1e-3)
        kn, bn = coordinate_values(result.x)
        full_residual = objective(result.x)
        results.append({"OCP_case": case_name, "Electrolyte": electrolyte,
                        "Start": start, "kn": kn, "brugg_n": bn,
                        "Train_voltage_RMSE_mV": float(np.sqrt(np.mean(full_residual[:3 * N_TRAIN] ** 2))),
                        "Objective_RMS_mV": float(np.sqrt(np.mean(full_residual**2))),
                        "N_function_evals": int(result.nfev),
                        "Optimizer_success": bool(result.success),
                        "Optimizer_message": str(result.message)})
        print("FITRESULT", results[-1], flush=True)
    if evaluations:
        pd.DataFrame(evaluations).to_csv(LOG_PATH, mode="a", header=not LOG_PATH.exists(),
                                         index=False, encoding="utf-8-sig")
    best = min(results, key=lambda row: row["Objective_RMS_mV"])
    return best, results


def potential_metrics(run):
    pot = np.asarray(run["anode_potential_v"], float)
    soc = np.asarray(run["full_cell_soc_pct"], float)
    time_min = np.asarray(run["time_min"], float)
    window = (soc >= 10) & (soc <= 70)
    return {
        "Anode_min_full_mV": 1000 * float(np.min(pot)),
        "Anode_min_SOC10_70_mV": 1000 * float(np.min(pot[window])) if window.any() else np.nan,
        "Anode_final_mV": 1000 * float(pot[-1]),
        "Time_anode_below_0V_min": float(np.trapezoid((pot < 0).astype(float), time_min)),
    }


def evaluate_scenario(case_name, case, electrolyte, phase, kn, bn):
    stage = model_stage(case, electrolyte, kn, bn)
    metrics, curves = [], []
    for rate, charge in CONDITIONS:
        run = scope["run_dfn_branch"](stage, rate, charge)
        exp = scope["experimental_curves"][(rate, charge)]
        score = scope["compare_with_experiment"](exp, run)
        direction = "Charge" if charge else "Discharge"
        metrics.append({"OCP_case": case_name, "Electrolyte": electrolyte,
                        "Phase": phase, "C_rate": rate, "Direction": direction,
                        "kn": kn, "brugg_n": bn, **score, **potential_metrics(run),
                        "Sim_duration_min": float(run["time_min"][-1]),
                        "Experiment_duration_min": float(exp["time_min"][-1] - exp["time_min"][0])})
        for i in range(len(run["time_min"])):
            curves.append({"OCP_case": case_name, "Electrolyte": electrolyte,
                           "Phase": phase, "C_rate": rate, "Direction": direction,
                           "Series": "DFN", "Time_min": float(run["time_min"][i]),
                           "Voltage_V": float(run["voltage_v"][i]),
                           "SOC_pct": float(run["full_cell_soc_pct"][i]),
                           "Anode_potential_mV": 1000 * float(run["anode_potential_v"][i])})
        for t, v in zip(exp["time_min"] - exp["time_min"][0], exp["voltage_v"]):
            curves.append({"OCP_case": case_name, "Electrolyte": electrolyte,
                           "Phase": phase, "C_rate": rate, "Direction": direction,
                           "Series": "Experiment", "Time_min": float(t), "Voltage_V": float(v),
                           "SOC_pct": np.nan, "Anode_potential_mV": np.nan})
    return metrics, curves


def write_checkpoint(best_rows, start_rows, metric_rows, curve_rows):
    pd.DataFrame(best_rows).to_csv(OUT / "Stage3_charge_fit_best.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(start_rows).to_csv(OUT / "Stage3_charge_fit_multistart.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(metric_rows).to_csv(OUT / "Stage3_charge_fit_branch_metrics.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(curve_rows).to_csv(OUT / "Stage3_charge_fit_curves.csv", index=False, encoding="utf-8-sig")


def main():
    start = time.monotonic()
    best_rows, start_rows, metric_rows, curve_rows = [], [], [], []
    for case_name, case in cases().items():
        for electrolyte in ("Ai2020_corrected", "ORegan2022"):
            initial_kn = kinetic_constants(case)["kn"]
            print("BASE", case_name, electrolyte, flush=True)
            metrics, curves = evaluate_scenario(case_name, case, electrolyte,
                                                "experimental_initial", initial_kn, 2.914)
            metric_rows.extend(metrics)
            curve_rows.extend(curves)
            best, starts = fit_one(case_name, case, electrolyte)
            best_rows.append(best)
            start_rows.extend(starts)
            print("VALIDATE", case_name, electrolyte, best["kn"], best["brugg_n"], flush=True)
            metrics, curves = evaluate_scenario(case_name, case, electrolyte,
                                                "charge_fit_discharge_validation",
                                                best["kn"], best["brugg_n"])
            metric_rows.extend(metrics)
            curve_rows.extend(curves)
            write_checkpoint(best_rows, start_rows, metric_rows, curve_rows)
    branches = pd.DataFrame(metric_rows)
    grouped = branches.groupby(["OCP_case", "Electrolyte", "Phase", "Direction"], sort=False).agg(
        Full_RMSE_mV=("full_RMSE_mV", lambda s: float(np.sqrt(np.mean(np.square(s))))),
        Full_MAE_mV=("full_MAE_mV", "mean"),
        Center_10_70_RMSE_mV=("center_10_70_RMSE_mV", lambda s: float(np.sqrt(np.mean(np.square(s))))),
        Center_10_70_MAE_mV=("center_10_70_MAE_mV", "mean"),
        Capacity_error_RMS_pct=("capacity_error_pct", lambda s: float(np.sqrt(np.mean(np.square(s))))),
        Min_anode_full_mV=("Anode_min_full_mV", "min"),
        Min_anode_SOC10_70_mV=("Anode_min_SOC10_70_mV", "min"),
    ).reset_index()
    grouped.to_csv(OUT / "Stage3_charge_fit_direction_summary.csv", index=False, encoding="utf-8-sig")
    print("BEST\n", pd.DataFrame(best_rows).to_string(index=False), flush=True)
    print("SUMMARY\n", grouped.to_string(index=False), flush=True)
    print("ELAPSED", round(time.monotonic() - start, 1), "s", flush=True)


if __name__ == "__main__":
    main()
