"""One-at-a-time anode-potential ablation from historical fit to present fit."""
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import pybamm

from run_stage2_two_ocp import cases, kinetic_constants, scope, stage_with_values
from run_stage3_charge_fit import OREGAN, TRANSPORT, BN_KEY

OUT = Path(__file__).resolve().parents[1] / "outputs"
CASE = cases()["Ai2020_anode__measured_cathode"]
OLD = {"kn": 8.726541e-7, "kp": 5.899519906708763e-7,
       "bn": 2.533583, "rn_um": 3.7542, "rp_um": 4.2387}
NEW = {"kn": 9.70279348377558e-7, "kp": 4.768215585284518e-7,
       "bn": 2.5, "rn_um": 3.692, "rp_um": 4.239}
VARIANTS = [("Old fit restored", OLD)] + [
    (f"Only {name} updated", {**OLD, name: NEW[name]})
    for name in ("kn", "kp", "bn", "rn_um")
] + [("Current selected", NEW)]
ROWS = []
for label, vals in VARIANTS:
    kinetic = kinetic_constants(CASE)
    kinetic.update(kn=vals["kn"], kp=vals["kp"])
    stage = stage_with_values(CASE, kinetic)
    params = stage["parameters"].copy()
    params.update({key: OREGAN[key] for key in TRANSPORT})
    params.update({BN_KEY: vals["bn"],
                   "Negative particle radius [m]": vals["rn_um"] * 1e-6,
                   "Positive particle radius [m]": vals["rp_um"] * 1e-6})
    stage["parameters"] = params
    for rate in (1.0, 2.0):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=pybamm.SolverWarning)
            run = scope["run_dfn_branch"](stage, rate, True)
        soc = np.asarray(run["full_cell_soc_pct"], float)
        pot = 1000 * np.asarray(run["anode_potential_v"], float)
        row = {"Variant": label, "C_rate": rate, **vals,
               "Min_anode_mV": float(np.min(pot)), "End_SOC_pct": float(soc[-1]),
               "SOC50_mV": float(np.interp(50, soc, pot)),
               "SOC70_mV": float(np.interp(70, soc, pot)) if soc[-1] >= 70 else np.nan}
        ROWS.append(row)
        print(row, flush=True)
    pd.DataFrame(ROWS).to_csv(OUT / "ORegan2022_과거현재_음극전위_파라미터기여.csv", index=False, encoding="utf-8-sig")
