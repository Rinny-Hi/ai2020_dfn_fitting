"""Stage 2 pre-fit and k_n-only charge fit from specified k_n=5.54e-7.

The model remains Ai2020-corrected electrolyte with Ai2020 anode OCP and
measured cathode GITT OCP. No final electrolyte change is made here.
"""
from __future__ import annotations

import math
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pybamm
from scipy.optimize import least_squares

from run_stage3_charge_fit import (
    KN_BOUNDS, N_TRAIN, SHORTFALL_WEIGHT_MV,
    model_stage, cases, scope, residuals_from_runs, potential_metrics,
)

OUT = Path(__file__).resolve().parents[1] / 'outputs'
CASE_NAME = 'Ai2020_anode__measured_cathode'
CASE = cases()[CASE_NAME]
ELECTROLYTE = 'Ai2020_corrected'
KN_INITIAL = 5.54e-7  # explicitly chosen from the user's anode Rct slide
BN_INITIAL = 2.914
KP_FIXED = 4.768215585284518e-7
RATES = (0.5, 1.0, 2.0)
CONDITIONS = tuple((rate, charge) for rate in RATES for charge in (True, False))
EVAL_PATH = OUT / 'Stage2_user_kn554_charge_fit_evaluations.csv'


def charge_runs(kn, bn):
    stage = model_stage(CASE, ELECTROLYTE, kn, bn)
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=pybamm.SolverWarning)
        return {rate: scope['run_dfn_branch'](stage, rate, True) for rate in RATES}


def evaluate(kn, bn, phase):
    stage = model_stage(CASE, ELECTROLYTE, kn, bn)
    rows = []
    for rate, charge in CONDITIONS:
        print('EVALUATE', phase, f'{rate:g}C', 'Charge' if charge else 'Discharge', flush=True)
        run = scope['run_dfn_branch'](stage, rate, charge)
        exp = scope['experimental_curves'][(rate, charge)]
        scores = scope['compare_with_experiment'](exp, run)
        rows.append({'Phase': phase, 'C_rate': rate,
                     'Direction': 'Charge' if charge else 'Discharge',
                     'Role': 'Fit' if charge else 'Independent validation',
                     'kn': kn, 'brugg_n': bn, 'kp': KP_FIXED,
                     **scores, **potential_metrics(run)})
        pd.DataFrame(rows).to_csv(OUT / f'Stage2_user_kn554_{phase}_branch_metrics.csv',
                                  index=False, encoding='utf-8-sig')
    return rows


def fit():
    cache = {}
    eval_rows = []

    def objective(z):
        kn, bn = float(np.exp(z[0])), BN_INITIAL
        key = round(kn, 15)
        if key not in cache:
            try:
                residual = residuals_from_runs(charge_runs(kn, bn))
                status = 'ok'
            except Exception as exc:
                residual = np.full(3 * N_TRAIN + 60, 1000.0)
                status = f'failed:{type(exc).__name__}'
            cache[key] = residual
            record = {'eval': len(eval_rows)+1, 'kn': kn, 'brugg_n': bn,
                      'train_voltage_RMSE_mV': float(np.sqrt(np.mean(residual[:3*N_TRAIN]**2))),
                      'objective_RMS_mV': float(np.sqrt(np.mean(residual**2))),
                      'status': status}
            eval_rows.append(record)
            pd.DataFrame(eval_rows).to_csv(EVAL_PATH, index=False, encoding='utf-8-sig')
            print('FIT EVAL', record, flush=True)
        return cache[key]

    lo = [math.log(KN_BOUNDS[0])]
    hi = [math.log(KN_BOUNDS[1])]
    seed = [math.log(KN_INITIAL)]
    result = least_squares(objective, seed, bounds=(lo, hi), max_nfev=20,
                           ftol=1e-4, xtol=1e-4, gtol=1e-4, diff_step=1e-3)
    kn, bn = float(np.exp(result.x[0])), BN_INITIAL
    residual = objective(result.x)
    fit_record = {'fit_parameters': 'kn only', 'start_kn': KN_INITIAL, 'start_brugg_n': BN_INITIAL,
                  'fitted_kn': kn, 'fitted_brugg_n': bn, 'kp_fixed': KP_FIXED,
                  'train_voltage_RMSE_mV': float(np.sqrt(np.mean(residual[:3*N_TRAIN]**2))),
                  'objective_RMS_mV': float(np.sqrt(np.mean(residual**2))),
                  'optimizer_success': bool(result.success), 'message': str(result.message),
                  'n_function_evals': int(result.nfev),
                  'shortfall_weight_mV': SHORTFALL_WEIGHT_MV}
    pd.DataFrame([fit_record]).to_csv(OUT / 'Stage2_user_kn554_charge_fit_result.csv',
                                       index=False, encoding='utf-8-sig')
    print('FIT RESULT', fit_record, flush=True)
    return kn, bn


def main():
    started = time.monotonic()
    initial = evaluate(KN_INITIAL, BN_INITIAL, 'Initial_kn554')
    kn, bn = fit()
    fitted = evaluate(kn, bn, 'Charge_fit_discharge_validation')
    metrics = pd.DataFrame(initial + fitted)
    metrics.to_csv(OUT / 'Stage2_user_kn554_initial_vs_fitted_branch_metrics.csv',
                   index=False, encoding='utf-8-sig')
    summary = metrics.groupby(['Phase', 'Direction'], sort=False).agg(
        Full_RMSE_mV=('full_RMSE_mV', lambda a: float(np.sqrt(np.mean(np.square(a))))),
        Full_MAE_mV=('full_MAE_mV', 'mean'),
        Central_RMSE_mV=('center_10_70_RMSE_mV', lambda a: float(np.sqrt(np.mean(np.square(a))))),
        Central_MAE_mV=('center_10_70_MAE_mV', 'mean'),
        Capacity_error_RMS_pct=('capacity_error_pct', lambda a: float(np.sqrt(np.mean(np.square(a))))),
        Min_anode_mV=('Anode_min_full_mV', 'min'),
    ).reset_index()
    summary.to_csv(OUT / 'Stage2_user_kn554_initial_vs_fitted_direction_summary.csv',
                   index=False, encoding='utf-8-sig')
    print('SUMMARY', summary.to_string(index=False), flush=True)
    print('ELAPSED_S', time.monotonic()-started, flush=True)


if __name__ == '__main__':
    main()
