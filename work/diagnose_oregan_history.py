"""Controlled 2C comparison of the October ORegan parameter revisions.

All rows use the same current measured cathode OCP, window, experimental curve,
and ORegan bulk-electrolyte functions. Only the named parameter is changed.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_stage2_two_ocp import cases, kinetic_constants, scope, stage_with_values
from run_stage3_charge_fit import OREGAN, TRANSPORT, BN_KEY

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "ORegan2022_과거현재_통제비교.csv"
CASE = cases()["Ai2020_anode__measured_cathode"]
BASE = kinetic_constants(CASE)
OLD_KN = 9.648533212e-7
OLD_KP = 5.899519906708763e-7

VARIANTS = [
    ("현재 실험 초기값", BASE["kn"], BASE["kp"], 2.914, 3.692, 4.239),
    ("옛 kn만 복원", OLD_KN, BASE["kp"], 2.914, 3.692, 4.239),
    ("옛 kp만 복원", BASE["kn"], OLD_KP, 2.914, 3.692, 4.239),
    ("옛 kn+kp 복원", OLD_KN, OLD_KP, 2.914, 3.692, 4.239),
    ("옛 kn+kp+입자반경 복원", OLD_KN, OLD_KP, 2.914, 3.7542, 4.2387),
    ("현재 EIS+k에 옛 공동피팅 b", BASE["kn"], BASE["kp"], 2.533583, 3.692, 4.239),
    ("옛 kn+kp에 공동피팅 b", OLD_KN, OLD_KP, 2.533583, 3.7542, 4.2387),
    ("10월20일 공동피팅값", 8.726541e-7, OLD_KP, 2.533583, 3.7542, 4.2387),
    ("이번 충전전용 피팅값", 3e-6, BASE["kp"], 3.05, 3.692, 4.239),
]


def run_variant(name, kn, kp, bn, rn, rp, charge):
    values = dict(BASE)
    values.update(kn=kn, kp=kp)
    stage = stage_with_values(CASE, values)
    params = stage["parameters"].copy()
    params.update({key: OREGAN[key] for key in TRANSPORT})
    params.update({BN_KEY: bn,
                   "Negative particle radius [m]": rn * 1e-6,
                   "Positive particle radius [m]": rp * 1e-6})
    stage["parameters"] = params
    sim = scope["run_dfn_branch"](stage, 2.0, charge)
    exp = scope["experimental_curves"][(2.0, charge)]
    m = scope["compare_with_experiment"](exp, sim)
    qe = exp["capacity_ah"] - exp["capacity_ah"][0]
    qm = sim["capacity_ah"] - sim["capacity_ah"][0]
    qend = 0.995 * min(float(qe[-1]), float(qm[-1]))
    qgrid = np.linspace(0, qend, 400)
    ev = np.interp(qgrid, qe, exp["voltage_v"])
    mv = np.interp(qgrid, qm, sim["voltage_v"])
    cap_rmse = float(np.sqrt(np.mean((1000 * (mv - ev)) ** 2)))
    return {"Variant": name, "Direction": "Charge" if charge else "Discharge",
            "kn": kn, "kp": kp, "brugg_n": bn, "rn_um": rn, "rp_um": rp,
            "Time_RMSE_mV": m["full_RMSE_mV"],
            "Capacity_aligned_RMSE_mV": cap_rmse,
            "Cutoff_capacity_error_pct": m["capacity_error_pct"],
            "Simulation_duration_min": float(sim["time_min"][-1]),
            "Experimental_duration_min": float(exp["time_min"][-1] - exp["time_min"][0])}


def main():
    rows = pd.read_csv(OUT, encoding="utf-8-sig").to_dict("records") if OUT.exists() else []
    for variant in VARIANTS:
        if len(sys.argv) > 1 and sys.argv[1] not in variant[0]:
            continue
        for charge in (True, False):
            row = run_variant(*variant, charge)
            rows = [old for old in rows if (old["Variant"], old["Direction"]) !=
                    (row["Variant"], row["Direction"])]
            rows.append(row)
            print(row, flush=True)
            pd.DataFrame(rows).to_csv(OUT, index=False, encoding="utf-8-sig")


if __name__ == "__main__":
    main()
