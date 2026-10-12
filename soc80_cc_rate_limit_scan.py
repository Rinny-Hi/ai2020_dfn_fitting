"""Scan constant-current charging from 0 to 6C for 0--80% full-cell SOC.

The current P0 baseline is used: ORegan2022 electrolyte transport, negative
electrode Bruggeman exponent 2.5, and the retained (Ai2020/teardown-v2) k_n.
The target is defined from the measured C/50 full-cell capacity and the C-rate
uses the project's 2.28 Ah nominal-capacity convention.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

import comprehensive_prefit_cross_cohort as cross
import priority_initial_state_protocol_recheck as priority

import analyze_oregan_brugg_2p5_anode_potential as anode_analysis
import electrolyte_transport_set_comparison as comparison


OUT = Path(__file__).resolve().parent / "results" / "261012_soc80_cc_rate_scan"
TARGET_SOC = 0.80
MAX_C_RATE = 6.0


def simulate(rate: float, stage: dict, model: dict, anode: dict, cathode: dict) -> dict:
    local_stage = priority.stage_at_soc(stage, 0.0, True)
    return anode_analysis.run_branch(rate, True, local_stage, model, anode, cathode)


def metrics(rate: float, output: dict, target_ah: float) -> dict:
    q = np.asarray(output["capacity_Ah"], dtype=float)
    time_min = np.asarray(output["time_s"], dtype=float) / 60.0
    voltage = np.asarray(output["V"], dtype=float)
    potential = np.asarray(output["potential_V"], dtype=float)
    reached = bool(q[-1] >= target_ah - 1e-5)
    if reached:
        target_time = float(np.interp(target_ah, q, time_min))
        target_voltage = float(np.interp(target_ah, q, voltage))
        target_potential = float(np.interp(target_ah, q, potential))
        keep = q <= target_ah
        if not np.any(np.isclose(q, target_ah)):
            potential_window = np.r_[potential[keep], target_potential]
        else:
            potential_window = potential[keep]
        minimum_potential = float(np.min(potential_window))
    else:
        target_time = np.nan
        target_voltage = np.nan
        target_potential = np.nan
        minimum_potential = float(np.min(potential))
    return {
        "C_rate": float(rate),
        "Current_A": float(rate * cross.NOMINAL_CAPACITY_AH),
        "Reached_80pct_before_4p2V": reached,
        "Time_to_80pct_min": target_time,
        "Voltage_at_80pct_V": target_voltage,
        "Anode_potential_at_80pct_mV": 1000.0 * target_potential,
        "Minimum_anode_potential_to_80pct_mV": 1000.0 * minimum_potential,
        "Cutoff_capacity_Ah": float(q[-1]),
        "Cutoff_SOC_pct": 100.0 * float(q[-1] / target_ah * TARGET_SOC),
        "Cutoff_time_min": float(time_min[-1]),
        "Cutoff_voltage_V": float(voltage[-1]),
    }


def evaluate(rate: float, stage: dict, model: dict, anode: dict, cathode: dict, target_ah: float) -> dict:
    print(f"{rate:.5f}C", flush=True)
    try:
        return metrics(rate, simulate(rate, stage, model, anode, cathode), target_ah)
    except Exception as exc:
        return {
            "C_rate": float(rate),
            "Current_A": float(rate * cross.NOMINAL_CAPACITY_AH),
            "Reached_80pct_before_4p2V": False,
            "Time_to_80pct_min": np.nan,
            "Voltage_at_80pct_V": np.nan,
            "Anode_potential_at_80pct_mV": np.nan,
            "Minimum_anode_potential_to_80pct_mV": np.nan,
            "Cutoff_capacity_Ah": np.nan,
            "Cutoff_SOC_pct": np.nan,
            "Cutoff_time_min": np.nan,
            "Cutoff_voltage_V": np.nan,
            "Error": repr(exc),
        }


def bisect_boundary(predicate, lo: float, hi: float, evaluator, iterations: int = 9) -> tuple[float, dict]:
    best_rate = lo
    best = evaluator(lo)
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        row = evaluator(mid)
        if predicate(row):
            best_rate, best, lo = mid, row, mid
        else:
            hi = mid
    return best_rate, best


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    baseline, stage, anode, cathode, _, qcell = comparison.current_p0()
    model = comparison.with_transport_set(baseline, "ORegan2022")
    model["params"].update(
        {"Negative electrode Bruggeman coefficient (electrolyte)": 2.5},
        check_already_exists=False,
    )
    target_ah = TARGET_SOC * qcell
    cache: dict[float, dict] = {}

    def evaluator(rate: float) -> dict:
        key = round(float(rate), 7)
        if key not in cache:
            cache[key] = evaluate(key, stage, model, anode, cathode, target_ah)
        return cache[key]

    coarse_rates = np.arange(0.25, MAX_C_RATE + 0.001, 0.25)
    coarse = [evaluator(float(rate)) for rate in coarse_rates]

    voltage_ok = lambda row: bool(row["Reached_80pct_before_4p2V"])
    plating_ok = lambda row: voltage_ok(row) and row["Minimum_anode_potential_to_80pct_mV"] >= 0.0

    def bracket(predicate):
        passing = [row["C_rate"] for row in coarse if predicate(row)]
        if not passing:
            return None
        lo = max(passing)
        failing = [row["C_rate"] for row in coarse if row["C_rate"] > lo and not predicate(row)]
        return (lo, min(failing)) if failing else (lo, MAX_C_RATE)

    results = {}
    for name, predicate in (("voltage_only", voltage_ok), ("voltage_and_anode", plating_ok)):
        pair = bracket(predicate)
        if pair is None:
            results[name] = None
            continue
        lo, hi = pair
        rate, row = (lo, evaluator(lo)) if hi == lo else bisect_boundary(predicate, lo, hi, evaluator)
        results[name] = {"max_feasible_C_rate": rate, **row}

    all_rows = pd.DataFrame(cache.values()).sort_values("C_rate")
    all_rows.to_csv(OUT / "soc80_cc_rate_scan.csv", index=False, encoding="utf-8-sig")
    summary = {
        "model": "ORegan2022 transport; b_n=2.5; retained k_n=1e-11 table value",
        "initial_SOC_pct": 0.0,
        "target_SOC_pct": 80.0,
        "C_rate_basis_capacity_Ah": cross.NOMINAL_CAPACITY_AH,
        "C50_full_cell_capacity_Ah": qcell,
        "target_throughput_Ah": target_ah,
        "ideal_6C_time_min": 60.0 * target_ah / (MAX_C_RATE * cross.NOMINAL_CAPACITY_AH),
        "upper_voltage_limit_V": 4.2,
        "anode_constraint_V_vs_Li": 0.0,
        **results,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print("\nSUMMARY")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print("\nCOARSE")
    print(pd.DataFrame(coarse).to_string(index=False))


if __name__ == "__main__":
    main()
