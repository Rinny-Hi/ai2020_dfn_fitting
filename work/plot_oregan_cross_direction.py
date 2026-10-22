"""Plot six-condition voltage checks and the charge/discharge fit tradeoff."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parents[1] / "outputs"
new = pd.read_csv(OUT / "ORegan2022_교차피팅_곡선.csv")
old = pd.read_csv(OUT / "Stage3_charge_fit_curves.csv")
old = old.loc[(old.Electrolyte == "ORegan2022") &
              (old.Phase == "charge_fit_discharge_validation")].copy()
old_metric = pd.read_csv(OUT / "Stage3_charge_fit_branch_metrics.csv")
old_metric = old_metric.loc[(old_metric.Electrolyte == "ORegan2022") &
                            (old_metric.Phase == "charge_fit_discharge_validation")].copy()
new_metric = pd.read_csv(OUT / "ORegan2022_교차피팅_조건별결과.csv")

STYLE = {
    "free_charge": ("Free b: charge fit", "#D1495B", ":"),
    "free_discharge": ("Free b: discharge fit", "#8556A3", "--"),
    "fixed_charge": ("b=2.5: charge fit", "#2166AC", "-"),
    "fixed_discharge": ("b=2.5: discharge fit", "#2A9D65", "-."),
}

def key(row):
    if row.get("fit_type") == "kn+bn":
        return "free_discharge"
    return "fixed_" + row["fit_direction"]

for case in new.OCP_case.unique():
    fig, axes = plt.subplots(2, 3, figsize=(14.8, 7.8), sharey=True)
    for col, rate in enumerate((0.5, 1.0, 2.0)):
        for ir, direction in enumerate(("Charge", "Discharge")):
            ax = axes[ir, col]
            g = new.loc[(new.OCP_case == case) & (new.C_rate == rate) &
                        (new.Direction == direction)]
            exp = g.loc[g.Series == "Experiment"].iloc[:]
            exp = exp.loc[(exp.fit_type == "kn+bn") & (exp.fit_direction == "discharge")]
            ax.plot(exp.Time_min, exp.Voltage_V, color="black", linewidth=2.1,
                    label="Experiment")
            prev = old.loc[(old.OCP_case == case) & (old.C_rate == rate) &
                           (old.Direction == direction) & (old.Series == "DFN")]
            name, color, ls = STYLE["free_charge"]
            ax.plot(prev.Time_min, prev.Voltage_V, color=color, linestyle=ls,
                    linewidth=1.65, label=name)
            for (kind, fitdir), sim in g.loc[g.Series == "DFN"].groupby(
                    ["fit_type", "fit_direction"], sort=False):
                label = "free_discharge" if kind == "kn+bn" else "fixed_" + fitdir
                name, color, ls = STYLE[label]
                ax.plot(sim.Time_min, sim.Voltage_V, color=color, linestyle=ls,
                        linewidth=1.75, label=name)
            ax.set_title(f"{rate:g}C {direction}")
            ax.set_xlabel("Elapsed time (min)")
            if col == 0:
                ax.set_ylabel("Terminal voltage (V)")
            ax.grid(alpha=.25)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5,
               bbox_to_anchor=(.5, -.01), fontsize=9)
    title = "Ai2020 graphite + measured cathode" if case.startswith("Ai2020") else "Measured graphite + measured cathode"
    fig.suptitle("ORegan2022 electrolyte: " + title, fontsize=15)
    fig.tight_layout(rect=(0, .065, 1, .955))
    suffix = "Ai2020_OCP" if case.startswith("Ai2020") else "Measured_OCP"
    fig.savefig(OUT / f"ORegan2022_교차피팅_6곡선_{suffix}.png", dpi=180,
                bbox_inches="tight")
    plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.4), sharex=True, sharey=True)
for ax, case in zip(axes, new_metric.OCP_case.unique()):
    rows = []
    for (kind, fitdir), group in new_metric.loc[new_metric.OCP_case == case].groupby(
            ["fit_type", "fit_direction"]):
        score = {d: float(np.sqrt(np.mean(group.loc[group.Direction == d, "full_RMSE_mV"]**2)))
                 for d in ("Charge", "Discharge")}
        rows.append(("free_discharge" if kind == "kn+bn" else "fixed_" + fitdir, score))
    g = old_metric.loc[old_metric.OCP_case == case]
    rows.append(("free_charge", {d: float(np.sqrt(np.mean(g.loc[g.Direction == d, "full_RMSE_mV"]**2)))
                                 for d in ("Charge", "Discharge")}))
    for label, score in rows:
        name, color, _ = STYLE[label]
        ax.scatter(score["Charge"], score["Discharge"], color=color, s=90,
                   label=name, zorder=3)
    ax.set_title("Ai2020 graphite" if case.startswith("Ai2020") else "Measured graphite")
    ax.set_xlabel("Charge RMSE across 0.5/1/2C (mV)")
    ax.grid(alpha=.25)
axes[0].set_ylabel("Discharge RMSE across 0.5/1/2C (mV)")
fig.suptitle("ORegan2022: opposite-direction validation tradeoff", fontsize=14)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=9,
           bbox_to_anchor=(.5, -.005))
fig.tight_layout(rect=(0, .07, 1, .95))
fig.savefig(OUT / "ORegan2022_교차피팅_충방전_절충.png", dpi=180,
            bbox_inches="tight")
plt.close(fig)
