"""Visual comparison of past joint fit and current selected anode potential."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

OUT = Path(__file__).resolve().parents[1] / "outputs"
d = pd.read_csv(OUT / "ORegan2022_과거공동피팅_vs_현재선정_음극전위_곡선.csv")
fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True, constrained_layout=True)
colors = {"2026-10-20_joint_fit": "#6A3D9A", "current_bn2p5_charge_fit": "#009E73"}
labels = {"2026-10-20_joint_fit": "Previous joint fit (10/20)",
          "current_bn2p5_charge_fit": "Current selected fit"}
for ax, (rate, group) in zip(axes, d.groupby("C_rate")):
    for name, sub in group.groupby("fit"):
        ax.plot(sub.SOC_pct, sub.Anode_potential_mV, lw=2.1, color=colors[name],
                ls="--" if name == "2026-10-20_joint_fit" else "-", label=labels[name])
    ax.axhline(0, color="#444444", lw=1)
    ax.set_xlim(0, 90)
    ax.set_ylim(-80, 150)
    ax.set_title(f"{rate:g}C charge")
    ax.set_xlabel("Full-cell SoC (%)")
    ax.grid(alpha=.18)
axes[0].set_ylabel("Predicted anode potential (mV)")
axes[0].legend(loc="upper right", fontsize=8, frameon=False)
fig.suptitle("Ai2020 graphite OCP + ORegan2022: previous vs current", fontsize=14, fontweight="bold")
fig.savefig(OUT / "ORegan2022_과거공동피팅_vs_현재선정_음극전위.png", dpi=190)
