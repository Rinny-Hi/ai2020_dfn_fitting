from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parents[1] / "outputs"
ranking = pd.read_csv(OUT / "Stage2_two_OCP_sensitivity.csv")
case_names = list(ranking.OCP_case.unique())
labels = ["Ai2020 anode + measured cathode", "Measured anode + cathode"]
parameters = ["Dsn", "Dsp", "kn", "kp"]

fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.3), layout="constrained")
for col, (case, label) in enumerate(zip(case_names, labels)):
    data = ranking[ranking.OCP_case == case].set_index("Parameter").loc[parameters]
    axes[0, col].barh(parameters, data.RMS_mV_per_log_parameter, color="#2458b4")
    axes[0, col].set_title(label)
    axes[0, col].set_xlabel("RMS voltage sensitivity (mV / ln parameter)")
    axes[0, col].grid(axis="x", alpha=0.2)
    corr = pd.read_csv(OUT / f"Stage2_{case}_Pearson_correlation.csv", index_col=0).loc[parameters, parameters]
    im = axes[1, col].imshow(corr, vmin=-1, vmax=1, cmap="coolwarm")
    axes[1, col].set_xticks(range(4), parameters)
    axes[1, col].set_yticks(range(4), parameters)
    axes[1, col].set_title("Pearson correlation of normalized voltage sensitivities")
    for i in range(4):
        for j in range(4):
            axes[1, col].text(j, i, f"{corr.iloc[i, j]:.2f}", ha="center", va="center", fontsize=9)
fig.colorbar(im, ax=axes[1, :], fraction=0.025, pad=0.03)
fig.savefig(OUT / "Stage2_two_OCP_sensitivity_identifiability.png", dpi=180)
plt.close(fig)
