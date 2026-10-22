"""Presentation-ready voltage and negative-potential figures from the audit CSVs."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parents[1] / "outputs"
curves = pd.read_csv(OUT / "Stage3_charge_fit_curves.csv")
metrics = pd.read_csv(OUT / "Stage3_true_fullcell_SOC10_70_metrics.csv")
pots = pd.read_csv(OUT / "Stage3_charge_fit_branch_metrics.csv")
names = {
    "Ai2020_anode__measured_cathode": "Ai2020 anode OCP",
    "measured_anode__measured_cathode": "Measured anode OCP",
}
styles = {"experimental_initial": ("#147da5", ":"),
          "charge_fit_discharge_validation": ("#c44e32", "--")}

for (ocp, electrolyte), block in curves.groupby(["OCP_case", "Electrolyte"], sort=False):
    fig, axes = plt.subplots(2, 3, figsize=(14.4, 7.3), sharey=True, layout="constrained")
    for col, rate in enumerate((0.5, 1., 2.)):
        for row, direction in enumerate(("Charge", "Discharge")):
            ax = axes[row, col]
            sub = block.loc[block.C_rate.eq(rate) & block.Direction.eq(direction)]
            exp = sub.loc[sub.Series.eq("Experiment") & sub.Phase.eq("experimental_initial")]
            ax.plot(exp.Time_min, exp.Voltage_V, color="#252525", lw=2, label="Experiment")
            for phase, (color, ls) in styles.items():
                m = sub.loc[sub.Series.eq("DFN") & sub.Phase.eq(phase)]
                ax.plot(m.Time_min, m.Voltage_V, color=color, ls=ls, lw=2, label="Initial" if phase == "experimental_initial" else "Charge-fit")
            met = metrics.loc[metrics.OCP_case.eq(ocp) & metrics.Electrolyte.eq(electrolyte)
                              & metrics.Phase.eq("charge_fit_discharge_validation") & metrics.C_rate.eq(rate)
                              & metrics.Direction.eq(direction)].iloc[0]
            ax.set_title(f"{rate:g}C {direction}  |  fit RMSE {met.Full_RMSE_mV:.1f} mV", fontsize=10.5)
            ax.grid(alpha=.2)
            if row == 1: ax.set_xlabel("Time (min)")
            if col == 0: ax.set_ylabel("Cell voltage (V)")
    axes[0, 0].legend(frameon=False, fontsize=9)
    fig.suptitle(f"{names[ocp]} | {electrolyte} electrolyte | charge fit / discharge validation", fontsize=14)
    fig.savefig(OUT / f"Stage3_fit_{ocp}_{electrolyte}_six_curves.png", dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9.6, 5.4), layout="constrained")
    for phase, (color, ls) in styles.items():
        for rate, alpha in ((0.5, .6), (1., .8), (2., 1.)):
            sub = block.loc[block.C_rate.eq(rate) & block.Direction.eq("Charge")
                            & block.Phase.eq(phase) & block.Series.eq("DFN")]
            ax.plot(sub.SOC_pct, sub.Anode_potential_mV, color=color, ls=ls, alpha=alpha,
                    lw=1.8 + .3 * rate, label=f"{'Initial' if phase == 'experimental_initial' else 'Fit'} {rate:g}C")
    ax.axhline(0, color="black", lw=1)
    ax.axhline(20, color="#777", ls="-.", lw=.9)
    ax.set(xlabel="Full-cell SOC (%)", ylabel="Anode surface potential vs electrolyte (mV)",
           title=f"{names[ocp]} | {electrolyte} electrolyte | charge")
    ax.grid(alpha=.2)
    ax.legend(ncol=3, fontsize=8, frameon=False)
    fig.savefig(OUT / f"Stage3_fit_{ocp}_{electrolyte}_anode_potential.png", dpi=170)
    plt.close(fig)

# One slide-size comparison of fitted charge/discharge voltage for the two electrolyte sets.
for ocp, block in curves.groupby("OCP_case", sort=False):
    fig, axes = plt.subplots(2, 3, figsize=(14.4, 7.3), sharey=True, layout="constrained")
    for col, rate in enumerate((.5, 1., 2.)):
        for row, direction in enumerate(("Charge", "Discharge")):
            ax = axes[row, col]
            sub = block.loc[block.C_rate.eq(rate) & block.Direction.eq(direction)
                            & block.Phase.eq("charge_fit_discharge_validation")]
            exp = sub.loc[sub.Series.eq("Experiment") & sub.Electrolyte.eq("Ai2020_corrected")]
            ax.plot(exp.Time_min, exp.Voltage_V, color="#222", lw=2, label="Experiment")
            for electrolyte, color in (("Ai2020_corrected", "#147da5"), ("ORegan2022", "#c44e32")):
                dat = sub.loc[sub.Series.eq("DFN") & sub.Electrolyte.eq(electrolyte)]
                ax.plot(dat.Time_min, dat.Voltage_V, ls="--", color=color, lw=2, label=electrolyte)
            ax.set_title(f"{rate:g}C {direction}", fontsize=10.5)
            ax.grid(alpha=.2)
            if row == 1: ax.set_xlabel("Time (min)")
            if col == 0: ax.set_ylabel("Cell voltage (V)")
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.suptitle(f"{names[ocp]} | charge-fitted electrolyte comparison", fontsize=14)
    fig.savefig(OUT / f"Stage3_fit_{ocp}_electrolyte_comparison.png", dpi=170)
    plt.close(fig)

print("Created final figures")
