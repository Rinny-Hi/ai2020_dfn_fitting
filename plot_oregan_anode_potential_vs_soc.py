"""Plot ORegan2022 anode-potential predictions against full-cell SOC."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "results" / "261012_anode_soc_crate_electrolyte"
OUT = ROOT / "results" / "261012_oregan_anode_potential_vs_soc"
RATES = (0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    curves = pd.read_csv(SOURCE / "anode_potential_curves.csv")
    curves = curves[
        (curves.Electrolyte_set == "ORegan2022")
        & curves.C_rate.isin(RATES)
    ].copy()
    cmap = plt.get_cmap("viridis")
    colors = {rate: cmap(i / (len(RATES) - 1)) for i, rate in enumerate(RATES)}
    summary_rows = []

    fig, axes = plt.subplots(1, 2, figsize=(15.2, 6.0), constrained_layout=True)
    for rate in RATES:
        data = curves[curves.C_rate == rate].sort_values("SOC_pct")
        soc = data.SOC_pct.to_numpy(float)
        potential = data.Anode_potential_mV.to_numpy(float)
        label = f"{rate:g}C"
        for ax in axes:
            ax.plot(soc, potential, lw=2.3, color=colors[rate], label=label)
            ax.scatter(
                [soc[-1]], [potential[-1]], s=38, color=colors[rate],
                edgecolor="white", linewidth=0.7, zorder=4,
            )
        below = np.flatnonzero(potential < 0.0)
        if len(below):
            index = int(below[0])
            if index == 0:
                crossing_soc = float(soc[0])
            else:
                crossing_soc = float(
                    np.interp(
                        0.0,
                        [potential[index], potential[index - 1]],
                        [soc[index], soc[index - 1]],
                    )
                )
        else:
            crossing_soc = np.nan
        summary_rows.append(
            {
                "C_rate": rate,
                "SOC_at_4p2V_cutoff_pct": float(soc[-1]),
                "Minimum_anode_potential_mV": float(np.min(potential)),
                "First_SOC_at_0V_pct": crossing_soc,
            }
        )

    for ax in axes:
        ax.axhspan(-120.0, 0.0, color="#D55E00", alpha=0.09, label="Predicted < 0 V")
        ax.axhline(0.0, color="#202020", lw=1.5, ls="--")
        ax.set_xlim(0.0, 100.0)
        ax.set_xlabel("Full-cell SOC (%)")
        ax.set_ylabel("Anode potential (mV vs. Li/Li$^+$)")
        ax.grid(alpha=0.23)
    axes[0].set_ylim(-100.0, 750.0)
    axes[0].set_title("Full range", weight="bold")
    axes[1].set_ylim(-100.0, 150.0)
    axes[1].set_title("Control-relevant zoom near 0 V", weight="bold")
    handles, labels = axes[1].get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    fig.legend(
        unique.values(), unique.keys(), loc="lower center", ncol=8,
        frameon=False, bbox_to_anchor=(0.5, -0.07),
    )
    fig.suptitle(
        "DFN-predicted anode potential vs full-cell SOC during CC charge\n"
        "ORegan2022 electrolyte, $b_n=2.5$, retained $k_n$; curves end at 4.2 V cutoff",
        fontsize=15, weight="bold", y=1.07,
    )
    figure = OUT / "oregan2022_anode_potential_vs_soc_by_crate.png"
    fig.savefig(figure, dpi=240, bbox_inches="tight")
    plt.close(fig)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT / "oregan2022_anode_potential_vs_soc_summary.csv", index=False, encoding="utf-8-sig")
    print(summary.to_string(index=False))
    print(figure)


if __name__ == "__main__":
    main()
