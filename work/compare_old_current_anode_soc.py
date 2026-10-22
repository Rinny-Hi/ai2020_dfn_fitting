"""Compare 10/20 joint fit and present selected fit at matched full-cell SoC."""
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import pybamm

from run_stage2_two_ocp import cases, kinetic_constants, scope, stage_with_values
from run_stage3_charge_fit import OREGAN, TRANSPORT, BN_KEY

OUT = Path(__file__).resolve().parents[1] / "outputs"
CASE = cases()["Ai2020_anode__measured_cathode"]
OLD_KN = 8.726541e-7
OLD_KP = 5.899519906708763e-7
OLD_BN = 2.533583
RATES = (0.5, 1.0, 2.0)
SOC_POINTS = (30.0, 50.0, 70.0, 75.0, 78.0)

values = kinetic_constants(CASE)
values.update(kn=OLD_KN, kp=OLD_KP)
stage = stage_with_values(CASE, values)
params = stage["parameters"].copy()
params.update({key: OREGAN[key] for key in TRANSPORT})
params.update({
    BN_KEY: OLD_BN,
    "Negative particle radius [m]": 3.7542e-6,
    "Positive particle radius [m]": 4.2387e-6,
})
stage["parameters"] = params

new = pd.read_csv(OUT / "ORegan2022_교차피팅_곡선.csv")
new = new.loc[
    (new.OCP_case == "Ai2020_anode__measured_cathode")
    & (new.fit_type == "kn_only_bn2.5")
    & (new.fit_direction == "charge")
    & (new.Direction == "Charge")
    & (new.Series == "DFN")
]
rows, curves = [], []
for rate in RATES:
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=pybamm.SolverWarning)
        old = scope["run_dfn_branch"](stage, rate, True)
    old_soc = np.asarray(old["full_cell_soc_pct"], float)
    old_pot = 1000 * np.asarray(old["anode_potential_v"], float)
    present = new.loc[new.C_rate == rate].sort_values("SOC_pct")
    new_soc = present.SOC_pct.to_numpy(float)
    new_pot = present.Anode_potential_mV.to_numpy(float)
    for label, soc, pot in (("2026-10-20_joint_fit", old_soc, old_pot),
                            ("current_bn2p5_charge_fit", new_soc, new_pot)):
        below = soc[pot < 0]
        data = {"fit": label, "C_rate": rate,
                "Min_mV": float(np.min(pot)), "End_SOC_pct": float(soc[-1]),
                "Negative_SOC_start_pct": float(below.min()) if len(below) else np.nan,
                "Negative_SOC_end_pct": float(below.max()) if len(below) else np.nan}
        for target in SOC_POINTS:
            data[f"SOC{target:g}_mV"] = float(np.interp(target, soc, pot)) if soc[0] <= target <= soc[-1] else np.nan
        rows.append(data)
        curves.extend({"fit": label, "C_rate": rate, "SOC_pct": float(s),
                       "Anode_potential_mV": float(p)} for s, p in zip(soc, pot))
    print(rate, rows[-2], rows[-1], flush=True)
pd.DataFrame(rows).to_csv(OUT / "ORegan2022_과거공동피팅_vs_현재선정_음극전위_요약.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(curves).to_csv(OUT / "ORegan2022_과거공동피팅_vs_현재선정_음극전위_곡선.csv", index=False, encoding="utf-8-sig")
