"""Build clear Ai2020 OCP ORegan2022 fit/validation tables and figures."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parents[1] / "outputs"
CASE = "Ai2020_anode__measured_cathode"
RATES = (0.5, 1.0, 2.0)
DIRECTIONS = ("Charge", "Discharge")

initial = pd.read_csv(OUT / "Stage3_true_fullcell_SOC10_70_metrics.csv")
initial = initial.loc[(initial.OCP_case == CASE) &
                      (initial.Electrolyte == "ORegan2022") &
                      (initial.Phase == "experimental_initial")].copy()
initial_extra = pd.read_csv(OUT / "Stage3_charge_fit_branch_metrics.csv")
initial_extra = initial_extra.loc[(initial_extra.OCP_case == CASE) &
                                  (initial_extra.Electrolyte == "ORegan2022") &
                                  (initial_extra.Phase == "experimental_initial")]
initial = initial.merge(initial_extra[["C_rate", "Direction", "capacity_error_pct", "Anode_min_full_mV"]],
                        on=["C_rate", "Direction"], validate="one_to_one")
before = pd.read_csv(OUT / "ORegan2022_Ai2020_bn2.5_피팅전_조건별결과.csv")
after = pd.read_csv(OUT / "ORegan2022_교차피팅_조건별결과.csv")
after = after.loc[(after.OCP_case == CASE) & (after.fit_type == "kn_only_bn2.5") &
                  (after.fit_direction == "charge")].copy()

rows = []
for _, r in initial.iterrows():
    rows.append({"Stage": "A: EIS kn, original bn=2.914", "C_rate": r.C_rate,
                 "Direction": r.Direction, "Data_role": "Before fit",
                 "kn": 6.307602286116924e-7, "brugg_n": 2.914,
                 "kp": 4.768215585284518e-7,
                 "Full_RMSE_mV": r.Full_RMSE_mV, "Full_MAE_mV": r.Full_MAE_mV,
                 "SOC10_70_RMSE_mV": r.SOC10_70_RMSE_mV,
                 "SOC10_70_MAE_mV": r.SOC10_70_MAE_mV,
                 "SOC_window_low_pct": r.SOC_window_actual_low_pct,
                 "SOC_window_high_pct": r.SOC_window_actual_high_pct,
                 "Capacity_error_pct": r.capacity_error_pct,
                 "Anode_min_mV": r.Anode_min_full_mV})
for df, label, role in ((before, "B: EIS kn, bn=2.5", "Before k fit"),
                        (after, "C: charge-fitted kn, bn=2.5", None)):
    for _, r in df.iterrows():
        rows.append({"Stage": label, "C_rate": r.C_rate,
                     "Direction": r.Direction,
                     "Data_role": role or ("Fit" if r.Direction == "Charge" else "Validation"),
                     "kn": r.kn, "brugg_n": r.brugg_n, "kp": 4.768215585284518e-7,
                     "Full_RMSE_mV": r.full_RMSE_mV, "Full_MAE_mV": r.full_MAE_mV,
                     "SOC10_70_RMSE_mV": r.True_SOC10_70_RMSE_mV,
                     "SOC10_70_MAE_mV": r.True_SOC10_70_MAE_mV,
                     "SOC_window_low_pct": r.SOC_window_low_pct,
                     "SOC_window_high_pct": r.SOC_window_high_pct,
                     "Capacity_error_pct": r.capacity_error_pct,
                     "Anode_min_mV": r.Anode_min_full_mV})
metrics = pd.DataFrame(rows)
metrics.to_csv(OUT / "ORegan2022_Ai2020_최종_3단계_조건별_RMSE_MAE.csv", index=False, encoding="utf-8-sig")

summaries = []
for (stage, direction), g in metrics.groupby(["Stage", "Direction"], sort=False):
    summaries.append({"Stage": stage, "Direction": direction,
                      "Full_RMSE_3rates_mV": float(np.sqrt(np.mean(g.Full_RMSE_mV**2))),
                      "Full_MAE_3rates_mV": float(np.mean(g.Full_MAE_mV)),
                      "SOC10_70_RMSE_3rates_mV": float(np.sqrt(np.mean(g.SOC10_70_RMSE_mV**2))),
                      "SOC10_70_MAE_3rates_mV": float(np.mean(g.SOC10_70_MAE_mV)),
                      "Capacity_RMS_3rates_pct": float(np.sqrt(np.mean(g.Capacity_error_pct**2)))})
pd.DataFrame(summaries).to_csv(OUT / "ORegan2022_Ai2020_최종_3단계_방향별요약.csv", index=False, encoding="utf-8-sig")

old_curve = pd.read_csv(OUT / "Stage3_charge_fit_curves.csv")
old_curve = old_curve.loc[(old_curve.OCP_case == CASE) & (old_curve.Electrolyte == "ORegan2022") &
                          (old_curve.Phase == "experimental_initial")]
before_curve = pd.read_csv(OUT / "ORegan2022_Ai2020_bn2.5_피팅전_곡선.csv")
after_curve = pd.read_csv(OUT / "ORegan2022_교차피팅_곡선.csv")
after_curve = after_curve.loc[(after_curve.OCP_case == CASE) &
                              (after_curve.fit_type == "kn_only_bn2.5") &
                              (after_curve.fit_direction == "charge")]

plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
styles = [("A: EIS k, b=2.914", old_curve, "#7B61A8", ":"),
          ("B: EIS k, b=2.5", before_curve, "#0072B2", "--"),
          ("C: charge-fitted k, b=2.5", after_curve, "#D55E00", "-")]
for final_only, filename in ((False, "ORegan2022_Ai2020_피팅전후_충방전6곡선.png"),
                             (True, "ORegan2022_Ai2020_최종_충방전6곡선.png")):
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True, sharey=True)
    for col, rate in enumerate(RATES):
        for row, direction in enumerate(DIRECTIONS):
            ax = axes[row, col]
            exp = after_curve.loc[(after_curve.C_rate == rate) & (after_curve.Direction == direction) &
                                  (after_curve.Series == "Experiment")]
            ax.plot(exp.Time_min, exp.Voltage_V, color="#222222", lw=2.0, label="Experiment")
            for label, source, color, ls in (styles[-1:] if final_only else styles):
                sim = source.loc[(source.C_rate == rate) & (source.Direction == direction) &
                                 (source.Series == "DFN")]
                ax.plot(sim.Time_min, sim.Voltage_V, color=color, lw=2.0 if final_only else 1.7,
                        ls=ls, label=label)
            m = metrics.loc[(metrics.Stage == "C: charge-fitted kn, bn=2.5") &
                            (metrics.C_rate == rate) & (metrics.Direction == direction)].iloc[0]
            tag = "FIT" if direction == "Charge" else "VALIDATION"
            ax.set_title(f"{rate:g}C {direction} — {tag} | final RMSE {m.Full_RMSE_mV:.1f} mV")
            ax.set_xlabel("Elapsed time (min)")
            if col == 0:
                ax.set_ylabel("Cell voltage (V)")
            ax.grid(alpha=.2)
    axes[0, 0].legend(loc="lower right", fontsize=8, frameon=True, framealpha=.9)
    fig.suptitle("Ai2020 graphite OCP + measured cathode OCP + ORegan2022", fontsize=15, fontweight="bold")
    fig.savefig(OUT / filename, dpi=180)
    plt.close(fig)

final_charge = after_curve.loc[(after_curve.Direction == "Charge") & (after_curve.Series == "DFN")]
soc_rows = []
fig, ax = plt.subplots(figsize=(9.4, 5.4), constrained_layout=True)
colors = {0.5: "#0072B2", 1.0: "#009E73", 2.0: "#D55E00"}
for rate, g in final_charge.groupby("C_rate"):
    g = g.sort_values("SOC_pct")
    soc = g.SOC_pct.to_numpy(float)
    pot = g.Anode_potential_mV.to_numpy(float)
    ax.plot(soc, pot, color=colors[rate], lw=2.3, label=f"{rate:g}C charge")
    row = {"C_rate": rate, "Min_anode_mV": float(np.min(pot)),
           "End_SOC_pct": float(soc[-1])}
    below = soc[pot < 0]
    row["Below_zero_SOC_start_pct"] = float(below.min()) if len(below) else np.nan
    row["Below_zero_SOC_end_pct"] = float(below.max()) if len(below) else np.nan
    for target in (10, 30, 50, 70, 75, 78):
        row[f"SOC{target}_mV"] = float(np.interp(target, soc, pot)) if soc[0] <= target <= soc[-1] else np.nan
    soc_rows.append(row)
ax.axhline(0, color="#333333", lw=1.1)
ax.axvline(70, color="#888888", ls=":", lw=1)
ax.set_xlim(0, 90)
ax.set_ylim(-100, 180)
ax.set_xlabel("Full-cell SoC (%)")
ax.set_ylabel("Predicted anode potential at separator (mV)")
ax.set_title("Selected model: SoC-resolved charging anode potential")
ax.legend(frameon=False)
ax.grid(alpha=.18)
fig.savefig(OUT / "ORegan2022_Ai2020_최종_SoC별음극전위.png", dpi=190)
plt.close(fig)
pd.DataFrame(soc_rows).to_csv(OUT / "ORegan2022_Ai2020_최종_SoC별음극전위.csv", index=False, encoding="utf-8-sig")
print(metrics[["Stage", "C_rate", "Direction", "Full_RMSE_mV", "SOC10_70_RMSE_mV"]].round(2).to_string(index=False))
