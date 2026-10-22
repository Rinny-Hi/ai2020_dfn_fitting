"""ORegan: discharge fit -> charge validation, then b_n=2.5 k_n-only cross fits."""

from __future__ import annotations

import json
import math
import time
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybamm
from scipy.optimize import least_squares, minimize_scalar

from run_stage2_two_ocp import cases, kinetic_constants, scope
from run_stage3_charge_fit import model_stage, potential_metrics

OUT = Path(__file__).resolve().parents[1] / "outputs"
KN_BOUNDS = (3e-7, 3e-6)
BN_BOUNDS = (1.5, 3.5)
RATES = (0.5, 1.0, 2.0)
N_GRID = 140
PENALTY_MV = 300.0
Q_CELL = float(scope["CELL_CAPACITY_AH"])
TRACE = OUT / "ORegan2022_교차피팅_탐색기록.csv"
FIT = OUT / "ORegan2022_교차피팅_파라미터.csv"
METRIC = OUT / "ORegan2022_교차피팅_조건별결과.csv"
CURVE = OUT / "ORegan2022_교차피팅_곡선.csv"


def train_blocks(charge):
    blocks = {}
    for rate in RATES:
        exp = scope["experimental_curves"][(rate, charge)]
        t = np.asarray(exp["time_min"], float).copy()
        t -= t[0]
        grid = np.linspace(0, .98 * t[-1], N_GRID)
        blocks[rate] = {"t": grid, "v": np.interp(grid, t, exp["voltage_v"]),
                        "end": float(t[-1])}
    return blocks


BLOCKS = {True: train_blocks(True), False: train_blocks(False)}


def residuals(runs, charge):
    pieces = []
    penalties = []
    for rate in RATES:
        target = BLOCKS[charge][rate]
        run = runs[rate]
        pred = np.interp(target["t"], run["time_min"], run["voltage_v"])
        pieces.append(1000 * (pred - target["v"]))
        short = max(0.0, target["t"][-1] - run["time_min"][-1]) / target["t"][-1]
        penalties.extend([PENALTY_MV * short] * 20)
    return np.concatenate([*pieces, np.asarray(penalties)])


def fit_one(case_name, case, charge, fixed_bn=None):
    label = "charge" if charge else "discharge"
    initial = kinetic_constants(case)
    cache = {}
    trace = []
    def evaluate(kn, bn):
        key = (round(float(kn), 15), round(float(bn), 9))
        if key not in cache:
            try:
                stage = model_stage(case, "ORegan2022", kn, bn)
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", category=pybamm.SolverWarning)
                    runs = {rate: scope["run_dfn_branch"](stage, rate, charge) for rate in RATES}
                rr = residuals(runs, charge)
                status = "ok"
            except Exception as exc:
                rr = np.full(3 * N_GRID + 60, 1000.0)
                status = f"{type(exc).__name__}: {str(exc)[:120]}"
            cache[key] = rr
            row = {"OCP_case": case_name, "fit_direction": label,
                   "fixed_bn": fixed_bn, "kn": kn, "brugg_n": bn,
                   "Train_voltage_RMSE_mV": float(np.sqrt(np.mean(rr[:3*N_GRID]**2))),
                   "Objective_RMS_mV": float(np.sqrt(np.mean(rr**2))),
                   "Status": status}
            trace.append(row)
            pd.DataFrame([row]).to_csv(TRACE, mode="a", header=not TRACE.exists(),
                                        index=False, encoding="utf-8-sig")
            print("EVAL", case_name, label, fixed_bn, len(trace),
                  round(kn, 10), round(bn, 5), round(row["Objective_RMS_mV"], 3),
                  status, flush=True)
        return cache[key]

    if fixed_bn is None:
        starts = [(initial["kn"], 2.914), (8.726541e-7, 2.533583)]
        results = []
        for i, (kn0, bn0) in enumerate(starts, 1):
            print("START", case_name, label, i, kn0, bn0, flush=True)
            x0 = [math.log(max(KN_BOUNDS[0], min(KN_BOUNDS[1], kn0))), bn0]
            res = least_squares(
                lambda z: evaluate(float(np.exp(z[0])), float(z[1])), x0,
                bounds=([math.log(KN_BOUNDS[0]), BN_BOUNDS[0]],
                        [math.log(KN_BOUNDS[1]), BN_BOUNDS[1]]),
                diff_step=1e-3, max_nfev=18, ftol=1e-4, xtol=1e-4, gtol=1e-4,
            )
            kn, bn = float(np.exp(res.x[0])), float(res.x[1])
            rr = evaluate(kn, bn)
            results.append({"start": i, "kn": kn, "brugg_n": bn,
                            "Objective_RMS_mV": float(np.sqrt(np.mean(rr**2))),
                            "Train_voltage_RMSE_mV": float(np.sqrt(np.mean(rr[:3*N_GRID]**2))),
                            "optimizer_success": bool(res.success), "nfev": res.nfev})
        best = min(results, key=lambda x: x["Objective_RMS_mV"])
    else:
        # One-dimensional bounded search, with endpoints and EIS starting value
        # included explicitly so boundary optima are not hidden.
        res = minimize_scalar(
            lambda x: float(np.mean(evaluate(float(np.exp(x)), fixed_bn)**2)),
            bounds=(math.log(KN_BOUNDS[0]), math.log(KN_BOUNDS[1])),
            method="bounded", options={"xatol": 0.005, "maxiter": 18},
        )
        candidates = [float(np.exp(res.x)), *KN_BOUNDS, initial["kn"]]
        choices = []
        for kn in candidates:
            rr = evaluate(kn, fixed_bn)
            choices.append((float(np.sqrt(np.mean(rr**2))), kn,
                            float(np.sqrt(np.mean(rr[:3*N_GRID]**2)))))
        obj, kn, v_rmse = min(choices)
        best = {"start": "bounded_scalar", "kn": kn, "brugg_n": fixed_bn,
                "Objective_RMS_mV": obj, "Train_voltage_RMSE_mV": v_rmse,
                "optimizer_success": bool(res.success), "nfev": res.nfev}
    return {"OCP_case": case_name, "fit_direction": label,
            "fit_type": "kn+bn" if fixed_bn is None else "kn_only_bn2.5", **best}


def true_soc_metric(exp, run, charge):
    exp_t = np.asarray(exp["time_min"], float) - float(exp["time_min"][0])
    exp_q = np.asarray(exp["capacity_ah"], float) - float(exp["capacity_ah"][0])
    initial = float(run["full_cell_soc_pct"][0])
    soc = initial + (1 if charge else -1) * 100 * exp_q / Q_CELL
    end = .995 * min(float(exp_t[-1]), float(run["time_min"][-1]))
    keep = (soc >= 10) & (soc <= 70) & (exp_t <= end)
    if keep.sum() < 2:
        return {"True_SOC10_70_RMSE_mV": np.nan, "True_SOC10_70_MAE_mV": np.nan,
                "SOC_window_low_pct": np.nan, "SOC_window_high_pct": np.nan}
    grid = np.linspace(float(exp_t[keep][0]), float(exp_t[keep][-1]), 160)
    err = 1000 * (np.interp(grid, run["time_min"], run["voltage_v"])
                  - np.interp(grid, exp_t, exp["voltage_v"]))
    return {"True_SOC10_70_RMSE_mV": float(np.sqrt(np.mean(err**2))),
            "True_SOC10_70_MAE_mV": float(np.mean(np.abs(err))),
            "SOC_window_low_pct": float(min(soc[keep][0], soc[keep][-1])),
            "SOC_window_high_pct": float(max(soc[keep][0], soc[keep][-1]))}


def evaluate_both(case_name, case, fit):
    stage = model_stage(case, "ORegan2022", fit["kn"], fit["brugg_n"])
    rows, curves = [], []
    for rate in RATES:
        for charge in (True, False):
            direction = "Charge" if charge else "Discharge"
            run = scope["run_dfn_branch"](stage, rate, charge)
            exp = scope["experimental_curves"][(rate, charge)]
            m = scope["compare_with_experiment"](exp, run)
            m.update(true_soc_metric(exp, run, charge))
            m.update(potential_metrics(run))
            rows.append({"OCP_case": case_name, "fit_type": fit["fit_type"],
                         "fit_direction": fit["fit_direction"], "C_rate": rate,
                         "Direction": direction,
                         "Data_role": "Fit" if fit["fit_direction"] == direction.lower() else "Validation",
                         "kn": fit["kn"], "brugg_n": fit["brugg_n"], **m,
                         "Sim_duration_min": float(run["time_min"][-1]),
                         "Experiment_duration_min": float(exp["time_min"][-1]-exp["time_min"][0])})
            for series, t, v, soc, pot in (
                ("DFN", run["time_min"], run["voltage_v"], run["full_cell_soc_pct"],
                 1000 * run["anode_potential_v"]),
                ("Experiment", exp["time_min"]-exp["time_min"][0], exp["voltage_v"],
                 np.full(len(exp["voltage_v"]), np.nan), np.full(len(exp["voltage_v"]), np.nan)),
            ):
                for vals in zip(t, v, soc, pot):
                    curves.append({"OCP_case": case_name, "fit_type": fit["fit_type"],
                                   "fit_direction": fit["fit_direction"], "C_rate": rate,
                                   "Direction": direction, "Series": series,
                                   "Time_min": float(vals[0]), "Voltage_V": float(vals[1]),
                                   "SOC_pct": float(vals[2]), "Anode_potential_mV": float(vals[3])})
            print("VALIDATE", case_name, fit["fit_type"], fit["fit_direction"],
                  rate, direction, round(m["full_RMSE_mV"], 2),
                  round(m["capacity_error_pct"], 2), flush=True)
    return rows, curves


def write_checkpoint(fits, metrics, curves):
    pd.DataFrame(fits).to_csv(FIT, index=False, encoding="utf-8-sig")
    pd.DataFrame(metrics).to_csv(METRIC, index=False, encoding="utf-8-sig")
    pd.DataFrame(curves).to_csv(CURVE, index=False, encoding="utf-8-sig")


def main():
    began = time.monotonic()
    fits, metrics, curves = [], [], []
    for case_name, case in cases().items():
        # Discharge two-parameter fit first; the opposite direction is untouched validation.
        for charge, bn in ((False, None), (True, 2.5), (False, 2.5)):
            fit = fit_one(case_name, case, charge, bn)
            fits.append(fit)
            print("FITRESULT", fit, flush=True)
            new_m, new_c = evaluate_both(case_name, case, fit)
            metrics.extend(new_m)
            curves.extend(new_c)
            write_checkpoint(fits, metrics, curves)
    print("TOTAL_SECONDS", round(time.monotonic() - began, 1), flush=True)


if __name__ == "__main__":
    main()
