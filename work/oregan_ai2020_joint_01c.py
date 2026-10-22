"""Fit Ai2020 OCP ORegan k_n and b_n on 0.5/1C both directions; hold out 2C."""
from __future__ import annotations

import math
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pybamm
from scipy.optimize import least_squares, minimize_scalar

from oregan_cross_direction_fit import BLOCKS, KN_BOUNDS, BN_BOUNDS, N_GRID, PENALTY_MV
from oregan_cross_direction_fit import evaluate_both
from run_stage2_two_ocp import cases
from run_stage3_charge_fit import model_stage
from stage2_setup_probe import scope

OUT = Path(__file__).resolve().parents[1] / "outputs"
CASE_NAME = "Ai2020_anode__measured_cathode"
CASE = cases()[CASE_NAME]
RATES = (0.5, 1.0)
TRACE = []
CACHE = {}


def residual(kn: float, bn: float) -> np.ndarray:
    key = (round(kn, 15), round(bn, 9))
    if key not in CACHE:
        try:
            stage = model_stage(CASE, "ORegan2022", kn, bn)
            pieces = []
            penalties = []
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=pybamm.SolverWarning)
                for charge in (True, False):
                    for rate in RATES:
                        run = scope["run_dfn_branch"](stage, rate, charge)
                        target = BLOCKS[charge][rate]
                        pred = np.interp(target["t"], run["time_min"], run["voltage_v"])
                        pieces.append(1000 * (pred - target["v"]))
                        short = max(0.0, target["t"][-1] - run["time_min"][-1]) / target["t"][-1]
                        penalties.extend([PENALTY_MV * short] * 20)
            rr = np.concatenate([*pieces, np.asarray(penalties)])
            status = "ok"
        except Exception as exc:
            rr = np.full(4 * N_GRID + 80, 1000.0)
            status = f"{type(exc).__name__}: {str(exc)[:120]}"
        CACHE[key] = rr
        row = {"kn": kn, "brugg_n": bn, "Objective_RMS_mV": float(np.sqrt(np.mean(rr**2))), "status": status}
        TRACE.append(row)
        pd.DataFrame(TRACE).to_csv(OUT / "ORegan2022_Ai2020_1C이하_공동피팅_탐색.csv", index=False, encoding="utf-8-sig")
        print("EVAL", len(TRACE), kn, bn, row["Objective_RMS_mV"], status, flush=True)
    return CACHE[key]


def fit_joint():
    options = []
    for kn, bn in ((8.726541e-7, 2.533583), (7.439385e-7, 2.44688)):
        z0 = [math.log(kn), bn]
        opt = least_squares(
            lambda z: residual(float(np.exp(z[0])), float(z[1])), z0,
            bounds=([math.log(KN_BOUNDS[0]), BN_BOUNDS[0]], [math.log(KN_BOUNDS[1]), BN_BOUNDS[1]]),
            diff_step=1e-3, max_nfev=20, ftol=1e-4, xtol=1e-4, gtol=1e-4,
        )
        k, b = float(np.exp(opt.x[0])), float(opt.x[1])
        options.append((float(np.sqrt(np.mean(residual(k, b)**2))), k, b))
    obj, kn, bn = min(options)
    return {"OCP_case": CASE_NAME, "fit_type": "joint_0p5_1C_kn_bn", "fit_direction": "both", "kn": kn, "brugg_n": bn, "Objective_RMS_mV": obj}


def fit_fixed():
    res = minimize_scalar(
        lambda logkn: float(np.mean(residual(float(np.exp(logkn)), 2.5)**2)),
        bounds=(math.log(KN_BOUNDS[0]), math.log(KN_BOUNDS[1])), method="bounded",
        options={"xatol": 0.005, "maxiter": 18},
    )
    options = [(float(np.sqrt(np.mean(residual(kn, 2.5)**2))), kn) for kn in (float(np.exp(res.x)), *KN_BOUNDS)]
    obj, kn = min(options)
    return {"OCP_case": CASE_NAME, "fit_type": "joint_0p5_1C_kn_only_bn2p5", "fit_direction": "both", "kn": kn, "brugg_n": 2.5, "Objective_RMS_mV": obj}


def main():
    fits = [fit_joint(), fit_fixed()]
    metric_rows, curve_rows = [], []
    for fit in fits:
        print("FIT", fit, flush=True)
        rows, curves = evaluate_both(CASE_NAME, CASE, fit)
        for row in rows:
            row["Data_role"] = "Fit" if row["C_rate"] in RATES else "Validation"
        metric_rows.extend(rows)
        curve_rows.extend(curves)
    pd.DataFrame(fits).to_csv(OUT / "ORegan2022_Ai2020_1C이하_공동피팅_파라미터.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(metric_rows).to_csv(OUT / "ORegan2022_Ai2020_1C이하_공동피팅_조건별결과.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(curve_rows).to_csv(OUT / "ORegan2022_Ai2020_1C이하_공동피팅_곡선.csv", index=False, encoding="utf-8-sig")


if __name__ == "__main__":
    main()
