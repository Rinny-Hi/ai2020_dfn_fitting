"""Mesh refinement for selected b_n=2.5 charge-fit cross-validation candidates."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from run_stage2_two_ocp import cases, scope
from run_stage3_charge_fit import model_stage

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
nb = json.loads((OUT / "261020_Ai2020_1단계_OCP_Stoichiometry_Measured_Anode_비교.ipynb").read_text(encoding="utf-8"))
src = "".join(nb["cells"][13]["source"])
start = src.index("def run_dfn_branch(")
end = src.index("def compare_with_experiment(", start)
fn = src[start:end]
needle = 'var_pts={"x_n": 15, "x_s": 15, "x_p": 15, "r_n": 15, "r_p": 15}'
assert needle in fn
fn = fn.replace(needle, 'var_pts={"x_n": 25, "x_s": 25, "x_p": 25, "r_n": 25, "r_p": 25}')
env = dict(scope)
exec(compile(fn, "refined_DFN", "exec"), env)
run25 = env["run_dfn_branch"]

fits = pd.read_csv(OUT / "ORegan2022_교차피팅_파라미터.csv")
curves = pd.read_csv(OUT / "ORegan2022_교차피팅_곡선.csv")
rows = []
for ocp, case in cases().items():
    fit = fits.loc[(fits.OCP_case == ocp) & (fits.fit_type == "kn_only_bn2.5") &
                   (fits.fit_direction == "charge")].iloc[0]
    stage = model_stage(case, "ORegan2022", fit.kn, 2.5)
    for rate, charge in ((1., True), (2., True), (2., False)):
        direction = "Charge" if charge else "Discharge"
        run = run25(stage, rate, charge)
        ref = curves.loc[(curves.OCP_case == ocp) & (curves.fit_type == "kn_only_bn2.5") &
                         (curves.fit_direction == "charge") & (curves.C_rate == rate) &
                         (curves.Direction == direction) & (curves.Series == "DFN")]
        t = np.linspace(0, min(run["time_min"][-1], ref.Time_min.iloc[-1]) * .995, 400)
        dv = 1000 * (np.interp(t, run["time_min"], run["voltage_v"])
                     - np.interp(t, ref.Time_min, ref.Voltage_V))
        dp = 1000 * (np.interp(t, run["time_min"], run["anode_potential_v"])
                     - np.interp(t, ref.Time_min, ref.Anode_potential_mV / 1000))
        exp = scope["experimental_curves"][(rate, charge)]
        score = scope["compare_with_experiment"](exp, run)
        rows.append({"OCP_case": ocp, "C_rate": rate, "Direction": direction,
                     "Voltage_RMS_mesh_difference_mV": float(np.sqrt(np.mean(dv**2))),
                     "Potential_RMS_mesh_difference_mV": float(np.sqrt(np.mean(dp**2))),
                     "Mesh25_full_RMSE_mV": score["full_RMSE_mV"],
                     "Mesh25_capacity_error_pct": score["capacity_error_pct"],
                     "Mesh25_min_anode_mV": 1000 * float(np.min(run["anode_potential_v"])),
                     "Mesh15_duration_min": float(ref.Time_min.iloc[-1]),
                     "Mesh25_duration_min": float(run["time_min"][-1])})
        pd.DataFrame(rows).to_csv(OUT / "ORegan2022_bn2.5_mesh15_vs25.csv", index=False,
                                  encoding="utf-8-sig")
        print(rows[-1], flush=True)
