"""Plot SoC-resolved anode potential for the selected and rejected Ai2020 fits."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

OUT = Path(__file__).resolve().parents[1] / "outputs"
old = pd.read_csv(OUT / "ORegan2022_교차피팅_곡선.csv")
old = old.loc[
    (old.OCP_case == "Ai2020_anode__measured_cathode")
    & (old.fit_type == "kn_only_bn2.5")
    & (old.fit_direction == "charge")
    & (old.Direction == "Charge")
    & (old.Series == "DFN")
]
new = pd.read_csv(OUT / "ORegan2022_Ai2020_1C이하_공동피팅_곡선.csv")
new = new.loc[
    (new.fit_type == "joint_0p5_1C_kn_bn")
    & (new.Direction == "Charge")
    & (new.Series == "DFN")
    & (new.C_rate == 1.0)
]
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(figsize=(9.0, 5.4), constrained_layout=True)
colors = {0.5: "#0072B2", 1.0: "#009E73", 2.0: "#D55E00"}
for rate, rows in old.groupby("C_rate"):
    ax.plot(rows.SOC_pct, rows.Anode_potential_mV, color=colors[rate], lw=2.5,
            label=f"{rate:g}C, selected $b_n=2.5$")
ax.plot(new.SOC_pct, new.Anode_potential_mV, color="#6A3D9A", lw=2.0, ls="--",
        label="1C, free $k_n+b_n$ fit")
ax.axhline(0, color="#343434", lw=1.2)
ax.axvline(70, color="#777777", lw=1.0, ls=":")
ax.axvspan(79.0, 80.1, color="#009E73", alpha=.12)
ax.set_xlim(0, 90)
ax.set_ylim(-100, 180)
ax.set_xlabel("Full-cell SoC (%)")
ax.set_ylabel("Predicted anode potential at separator (mV)")
ax.set_title("Ai2020 graphite OCP + ORegan2022: charging anode potential")
ax.legend(loc="upper right", fontsize=9, frameon=False)
ax.annotate("Selected 1C: -3.3 mV only near 79-80% SoC",
            xy=(79.5, -3.3), xytext=(39, -45), fontsize=9,
            arrowprops={"arrowstyle": "->", "color": "#009E73"})
fig.savefig(OUT / "ORegan2022_Ai2020_SoC별_음극전위_선정비교.png", dpi=190)
