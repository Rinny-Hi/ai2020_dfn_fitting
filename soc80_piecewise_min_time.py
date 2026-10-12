"""Estimate the 0--80% minimum charge time with a 0--6C current bound.

Two estimates are produced:
1. Voltage-only: 6C CC followed by 4.2 V CV (the time-optimal boundary policy
   when only current and upper-voltage constraints are considered).
2. Voltage + anode-potential: a practical greedy boundary policy with eight
   10%-SOC constant-current blocks.  Each block uses the largest current that
   keeps terminal voltage <= 4.2 V and separator-side anode potential >= 0 V.

This is a screening calculation, not a formally certified optimal-control
solution.  It is intentionally reproducible and records the applied schedule.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pybamm
from scipy.integrate import cumulative_trapezoid

import comprehensive_prefit_cross_cohort as cross
import electrode_hysteresis_quantification as ehq
import priority_initial_state_protocol_recheck as priority

import analyze_oregan_brugg_2p5_anode_potential as anode_analysis
import electrolyte_transport_set_comparison as comparison
import oregan_bn25_kn_only_fit as kn_fit


TARGET_SOC = 0.80
MAX_C = 6.0
SOC_BLOCK = float(os.environ.get("SOC80_BLOCK", "0.10"))
KN_SCALE = float(os.environ.get("SOC80_KN_SCALE", "1.0"))
ELECTROLYTE_SET = os.environ.get("SOC80_ELECTROLYTE_SET", "ORegan2022")
BLOCK_LABEL = int(round(100.0 * SOC_BLOCK))
KN_LABEL = "anchor_kn" if abs(KN_SCALE - 1.0) < 1e-12 else f"kn_scale_{KN_SCALE:.6f}"
SET_LABEL = ELECTROLYTE_SET.lower().replace(" ", "_")
OUT = (
    Path(__file__).resolve().parent
    / "results"
    / f"261012_soc80_min_time_{BLOCK_LABEL}pct_blocks_{SET_LABEL}_{KN_LABEL}"
)
V_MAX = 4.2
ANODE_MIN = 0.0


def baseline_inputs():
    baseline, stage, anode, cathode, _, qcell = comparison.current_p0()
    model_inputs = comparison.with_transport_set(baseline, ELECTROLYTE_SET)
    model_inputs = kn_fit.apply_kn(model_inputs, KN_SCALE)
    model_inputs["params"].update(
        {"Negative electrode Bruggeman coefficient (electrolyte)": 2.5},
        check_already_exists=False,
    )
    initial_stage = priority.stage_at_soc(stage, 0.0, True)
    values = ehq.make_parameter_values(
        model_inputs, initial_stage, anode, cathode, True, False, False
    )
    model = pybamm.lithium_ion.DFN(
        {"open-circuit potential": ("single", "single"), "thermal": "isothermal"}
    )
    return model, values, qcell


def solve_steps(model, values, steps, period="5 seconds"):
    sim = pybamm.Simulation(
        model,
        parameter_values=values,
        experiment=pybamm.Experiment([tuple(steps)], period=period),
        var_pts={"x_n": 15, "x_s": 15, "x_p": 15, "r_n": 15, "r_p": 15},
        solver=pybamm.IDAKLUSolver(rtol=1e-6, atol=1e-8),
    )
    return sim.solve()


def arrays(solution):
    t = np.asarray(solution["Time [s]"].entries, dtype=float)
    current = np.abs(np.asarray(solution["Current [A]"].entries, dtype=float))
    q = cumulative_trapezoid(current, t, initial=0.0) / 3600.0
    v = np.asarray(solution["Terminal voltage [V]"].entries, dtype=float)
    p = np.squeeze(
        np.asarray(
            solution[
                "Negative electrode surface potential difference at separator interface [V]"
            ].entries,
            dtype=float,
        )
    )
    return t, q, v, p, current


def voltage_only(model, values, target_ah):
    solution = solve_steps(
        model,
        values,
        ["Charge at 6C until 4.2 V", "Hold at 4.2 V until C/100"],
        period="2 seconds",
    )
    t, q, v, p, current = arrays(solution)
    if q[-1] < target_ah:
        raise RuntimeError("CC-CV simulation did not reach 80% SOC")
    target_t = float(np.interp(target_ah, q, t))
    mask = t <= target_t
    p_target = float(np.interp(target_t, t, p))
    p_window = np.r_[p[mask], p_target]
    return {
        "time_min": target_t / 60.0,
        "minimum_anode_potential_mV": 1000.0 * float(np.min(p_window)),
        "current_at_80pct_C": float(np.interp(target_t, t, current) / cross.NOMINAL_CAPACITY_AH),
        "cc_end_time_min": float(solution.cycles[0].steps[0]["Time [s]"].entries[-1] / 60.0),
        "cc_end_SOC_pct": 100.0
        * float(
            np.interp(
                solution.cycles[0].steps[0]["Time [s]"].entries[-1], t, q
            )
            / (target_ah / TARGET_SOC)
        ),
    }


def step_text(rate, delta_q):
    seconds = 3600.0 * delta_q / (rate * cross.NOMINAL_CAPACITY_AH)
    return f"Charge at {rate * cross.NOMINAL_CAPACITY_AH:.9f} A for {seconds:.9f} seconds"


def evaluate_schedule(model, values, rates, delta_q, expected_q):
    steps = [step_text(rate, delta_q) for rate in rates]
    try:
        solution = solve_steps(model, values, steps)
        t, q, v, p, current = arrays(solution)
        reached = bool(q[-1] >= expected_q - 2e-4)
        return {
            "feasible": reached and float(np.max(v)) <= V_MAX + 2e-4 and float(np.min(p)) >= ANODE_MIN - 2e-5,
            "reached": reached,
            "end_q_Ah": float(q[-1]),
            "time_min": float(t[-1] / 60.0),
            "max_voltage_V": float(np.max(v)),
            "min_anode_V": float(np.min(p)),
            "solution": solution,
        }
    except Exception as exc:
        return {"feasible": False, "reached": False, "error": repr(exc)}


def greedy_anode_safe(model, values, qcell):
    delta_q = SOC_BLOCK * qcell
    n_blocks = int(round(TARGET_SOC / SOC_BLOCK))
    rates = []
    rows = []
    for block in range(n_blocks):
        expected_q = (block + 1) * delta_q
        lo, hi = 0.05, MAX_C
        best = None
        # A low rate should always be feasible; retain it explicitly.
        low_result = evaluate_schedule(model, values, rates + [lo], delta_q, expected_q)
        if not low_result["feasible"]:
            raise RuntimeError(f"No feasible current in block {block}: {low_result}")
        best = low_result
        for _ in range(11):
            mid = 0.5 * (lo + hi)
            result = evaluate_schedule(model, values, rates + [mid], delta_q, expected_q)
            if result["feasible"]:
                lo, best = mid, result
            else:
                hi = mid
        rates.append(lo)
        rows.append(
            {
                "SOC_start_pct": 100.0 * block * SOC_BLOCK,
                "SOC_end_pct": 100.0 * (block + 1) * SOC_BLOCK,
                "C_rate": lo,
                "Current_A": lo * cross.NOMINAL_CAPACITY_AH,
                "Block_time_min": 60.0 * delta_q / (lo * cross.NOMINAL_CAPACITY_AH),
                "Cumulative_time_min": best["time_min"],
                "Maximum_voltage_to_date_V": best["max_voltage_V"],
                "Minimum_anode_potential_to_date_mV": 1000.0 * best["min_anode_V"],
            }
        )
        print(f"block {block}: {lo:.4f}C, cumulative {best['time_min']:.3f} min", flush=True)
    final = evaluate_schedule(model, values, rates, delta_q, TARGET_SOC * qcell)
    return pd.DataFrame(rows), final


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    model, values, qcell = baseline_inputs()
    target_ah = TARGET_SOC * qcell
    voltage = voltage_only(model, values, target_ah)
    schedule, safe = greedy_anode_safe(model, values, qcell)
    schedule.to_csv(
        OUT / f"anode_safe_{BLOCK_LABEL}pct_schedule.csv",
        index=False,
        encoding="utf-8-sig",
    )
    summary = {
        "model": f"{ELECTROLYTE_SET} transport; b_n=2.5",
        "electrolyte_set": ELECTROLYTE_SET,
        "kn_scale_vs_retained_anchor": KN_SCALE,
        "kn_table_value": 1.0e-11 * KN_SCALE,
        "C_rate_basis_Ah": cross.NOMINAL_CAPACITY_AH,
        "C50_capacity_Ah": qcell,
        "target_throughput_Ah": target_ah,
        "ideal_6C_no_constraint_time_min": 60.0 * target_ah / (MAX_C * cross.NOMINAL_CAPACITY_AH),
        "voltage_only_6C_CC_CV": voltage,
        f"voltage_and_anode_safe_piecewise_{BLOCK_LABEL}pct": {
            "time_min": safe["time_min"],
            "max_voltage_V": safe["max_voltage_V"],
            "minimum_anode_potential_mV": 1000.0 * safe["min_anode_V"],
            "rates_C": schedule["C_rate"].tolist(),
        },
        "qualification": f"piecewise result is a greedy {BLOCK_LABEL}%-SOC screening estimate, not a certified global optimum",
    }
    (OUT / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(schedule.to_string(index=False))


if __name__ == "__main__":
    main()
