"""Six-condition scientific voltage plot for both OCP alternatives."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


OUT = Path(__file__).resolve().parents[1] / "outputs"
curves = pd.read_csv(OUT / "Stage3_prefit_two_OCP_curves.csv")
fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.0), sharey=True, layout="constrained")
for col, rate in enumerate((0.5, 1.0, 2.0)):
    for row, direction in enumerate(("Charge", "Discharge")):
        ax = axes[row, col]
        base = curves.loc[(curves.C_rate == rate) & (curves.Direction == direction)]
        exp = base.loc[(base.Series == "Experiment") & (base.OCP_case == "Ai2020_anode__measured_cathode")]
        ax.plot(exp.Time_min, exp.Voltage_V, color="#202020", lw=2.5, label="Experiment")
        for name, title, color in (
            ("Ai2020_anode__measured_cathode", "Ai2020 anode", "#2458a8"),
            ("measured_anode__measured_cathode", "Measured anode", "#d05c28"),
        ):
            data = base.loc[(base.Series == "DFN") & (base.OCP_case == name)]
            ax.plot(data.Time_min, data.Voltage_V, color=color, lw=2, label=title)
        ax.set_title(f"{rate:g}C {direction}")
        ax.grid(alpha=.24)
        if row == 1:
            ax.set_xlabel("Time (min)")
        if col == 0:
            ax.set_ylabel("Cell voltage (V)")
axes[0, 0].legend(frameon=False, fontsize=9)
fig.suptitle("Experimental Ds and kinetics | two OCP alternatives", fontsize=16)
fig.savefig(OUT / "Stage3_prefit_two_OCP_overlay_six_curves.png", dpi=190)
plt.close(fig)
