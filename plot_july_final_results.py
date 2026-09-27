from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
REPO = Path(
    r"C:\Users\user\Documents\Codex\2026-09-18"
    r"\https-github-com-rinny-hi-ai2020\work\ai2020_dfn_fitting"
)
SOURCE = REPO / "results" / "260927_charge_subset_fit_validation" / "conservative"
OUT = PROJECT / "outputs" / "260928_july_final_result_graph"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = pd.read_csv(SOURCE / "best_postfit_summary.csv")
    result = json.loads((SOURCE / "best_result.json").read_text(encoding="utf-8"))

    rates = np.array([0.5, 1.0, 2.0])
    charge = summary[summary.Direction == "Charge"].set_index("C_rate").loc[rates]
    discharge = summary[summary.Direction == "Discharge"].set_index("C_rate").loc[rates]
    x = np.arange(len(rates))
    width = 0.34

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(16.0, 5.2),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [1.2, 1.2, 1.0]},
    )
    axes[0].bar(x - width / 2, charge.Full_RMSE_mV, width, color="#0072B2", label="Charge fit")
    axes[0].bar(x + width / 2, discharge.Full_RMSE_mV, width, color="#D55E00", label="Discharge validation")
    axes[0].set_xticks(x, [f"{rate:g}C" for rate in rates])
    axes[0].set(ylabel="Full-curve RMSE [mV]", title="Voltage accuracy by rate")
    axes[0].grid(axis="y", alpha=0.25)
    axes[0].legend(fontsize=9)
    for container in axes[0].containers:
        axes[0].bar_label(container, fmt="%.1f", padding=3, fontsize=9)

    axes[1].bar(x - width / 2, charge.Capacity_RMSE_pct, width, color="#0072B2", label="Charge fit")
    axes[1].bar(x + width / 2, discharge.Capacity_RMSE_pct, width, color="#D55E00", label="Discharge validation")
    axes[1].set_xticks(x, [f"{rate:g}C" for rate in rates])
    axes[1].set(ylabel="Capacity RMSE [%]", title="Post-hoc cutoff-capacity error")
    axes[1].grid(axis="y", alpha=0.25)
    for container in axes[1].containers:
        axes[1].bar_label(container, fmt="%.2f", padding=3, fontsize=9)

    axes[2].axis("off")
    lines = [
        "Selected July model",
        "",
        "Fitted: Dsn + kn",
        f"Dsn  {result['Dsn']:.4e} m2/s",
        f"kn    {result['kn']:.4e}",
        f"brugg_n  {result['brugg_n']:.3f} (fixed)",
        "",
        f"Charge objective  {result['Charge_objective_RMSE_mV']:.2f} mV",
        f"Charge full RMSE  {result['Charge_full_RMSE_mV']:.2f} mV",
        f"Discharge full RMSE  {result['Discharge_full_RMSE_mV']:.2f} mV",
        "",
        f"Charge capacity RMSE  {result['Charge_capacity_RMSE_pct']:.2f}%",
        f"Discharge capacity RMSE  {result['Discharge_capacity_RMSE_pct']:.2f}%",
        "",
        "Capacity was not used in fitting.",
    ]
    axes[2].text(
        0.03,
        0.97,
        "\n".join(lines),
        va="top",
        ha="left",
        fontsize=11,
        linespacing=1.35,
        bbox=dict(boxstyle="round,pad=0.7", facecolor="#F3F6FA", edgecolor="#8A9AA9"),
    )
    fig.suptitle("July BoL final conservative model", fontsize=16)
    fig.savefig(OUT / "july_final_rate_metrics.png", dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    main()
