"""Audit model errors on an actual full-cell SoC 10–70% time window."""
from pathlib import Path
import numpy as np
import pandas as pd

from run_stage2_two_ocp import scope

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
Q = float(scope["CELL_CAPACITY_AH"])


def score_one(model, exp, rate, direction):
    sim = model.loc[model.Series.eq("DFN")].sort_values("Time_min")
    exp_t = np.asarray(exp["time_min"], float) - float(exp["time_min"][0])
    exp_q = np.asarray(exp["capacity_ah"], float)
    exp_q -= exp_q[0]
    sign = 1 if direction == "Charge" else -1
    initial_soc = float(sim.SOC_pct.iloc[0])
    exp_soc = initial_soc + sign * 100 * exp_q / Q
    sim_end = float(sim.Time_min.iloc[-1])
    common_end = .995 * min(exp_t[-1], sim_end)
    full_grid = np.linspace(0, common_end, 400)
    def error(grid):
        pred = np.interp(grid, sim.Time_min, sim.Voltage_V)
        obs = np.interp(grid, exp_t, exp["voltage_v"])
        return 1000 * (pred - obs)
    full_e = error(full_grid)
    # The full-cell SoC window is defined by experimental coulomb counting
    # from the branch-specific rest-voltage-derived initial SoC.
    keep = (exp_soc >= 10) & (exp_soc <= 70) & (exp_t <= common_end)
    if np.count_nonzero(keep) < 2:
        raise ValueError((rate, direction, exp_soc.min(), exp_soc.max(), common_end))
    t0, t1 = float(exp_t[keep][0]), float(exp_t[keep][-1])
    soc_grid = np.linspace(t0, t1, 160)
    mid_e = error(soc_grid)
    return {
        "Initial_SOC_pct": initial_soc,
        "Full_RMSE_mV": float(np.sqrt(np.mean(full_e**2))),
        "Full_MAE_mV": float(np.mean(abs(full_e))),
        "SOC10_70_RMSE_mV": float(np.sqrt(np.mean(mid_e**2))),
        "SOC10_70_MAE_mV": float(np.mean(abs(mid_e))),
        "SOC_window_actual_low_pct": float(np.interp(t0, exp_t, exp_soc)),
        "SOC_window_actual_high_pct": float(np.interp(t1, exp_t, exp_soc)),
        "SOC_window_start_min": t0,
        "SOC_window_end_min": t1,
        "Experiment_SOC_at_common_end_pct": float(np.interp(common_end, exp_t, exp_soc)),
        "Sim_end_min": sim_end,
        "Experiment_end_min": float(exp_t[-1]),
        "N_points": len(soc_grid),
    }


rows = []
for stem, phase in (("Stage3_prefit_two_OCP_curves.csv", "experimental_initial"),
                    ("Stage3_charge_fit_curves.csv", None)):
    curves = pd.read_csv(OUT / stem)
    fields = ["OCP_case", "C_rate", "Direction"] if phase else ["OCP_case", "Electrolyte", "Phase", "C_rate", "Direction"]
    for keys, group in curves.groupby(fields, sort=False):
        identity = dict(zip(fields, keys))
        rate = float(identity["C_rate"])
        direction = identity["Direction"]
        exp = scope["experimental_curves"][(rate, direction == "Charge")]
        row = {**identity, "Electrolyte": identity.get("Electrolyte", "Ai2020_corrected"),
               "Phase": identity.get("Phase", phase)}
        row.update(score_one(group, exp, rate, direction))
        rows.append(row)

frame = pd.DataFrame(rows)
frame = frame.drop_duplicates(["OCP_case", "Electrolyte", "Phase", "C_rate", "Direction"], keep="last")
frame.to_csv(OUT / "Stage3_true_fullcell_SOC10_70_metrics.csv", index=False, encoding="utf-8-sig")
print(frame.to_string(index=False))
