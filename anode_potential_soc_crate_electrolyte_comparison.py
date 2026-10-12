"""Compare charge anode-potential maps for ORegan2022 and Ai2020 electrolyte.

All non-electrolyte parameters are fixed at the retained P0 baseline, including
b_n=2.5 and the retained k_n.  Each constant-current charge begins at 0% SOC
and terminates at the 4.2 V cell-voltage limit.  Values beyond each cutoff are
left missing rather than extrapolated.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import priority_initial_state_protocol_recheck as priority

import analyze_oregan_brugg_2p5_anode_potential as anode_analysis
import electrolyte_transport_set_comparison as comparison


OUT = Path(__file__).resolve().parent / "results" / "261012_anode_soc_crate_electrolyte"
RATES = np.arange(0.25, 6.0001, 0.25)
SOC_GRID = np.linspace(0.0, 80.0, 321)
SETS = ("ORegan2022", "Ai2020 corrected")


def model_for_set(baseline: dict, set_name: str) -> dict:
    model = comparison.with_transport_set(baseline, set_name)
    model["params"].update(
        {"Negative electrode Bruggeman coefficient (electrolyte)": 2.5},
        check_already_exists=False,
    )
    return model


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    baseline, stage, anode, cathode, _, qcell = comparison.current_p0()
    initial_stage = priority.stage_at_soc(stage, 0.0, True)
    maps: dict[str, np.ndarray] = {}
    curve_rows: list[dict] = []
    summary_rows: list[dict] = []

    for set_name in SETS:
        model = model_for_set(baseline, set_name)
        potential_map = np.full((len(RATES), len(SOC_GRID)), np.nan)
        for row, rate in enumerate(RATES):
            print(f"{set_name}: {rate:.2f}C", flush=True)
            output = anode_analysis.run_branch(
                float(rate), True, initial_stage, model, anode, cathode
            )
            soc = 100.0 * np.asarray(output["capacity_Ah"], dtype=float) / qcell
            potential_mV = 1000.0 * np.asarray(output["potential_V"], dtype=float)
            voltage = np.asarray(output["V"], dtype=float)
            time_min = np.asarray(output["time_s"], dtype=float) / 60.0
            max_soc = float(soc[-1])
            mask = SOC_GRID <= max_soc + 1e-9
            potential_map[row, mask] = np.interp(SOC_GRID[mask], soc, potential_mV)
            crossing = np.flatnonzero(potential_mV < 0.0)
            first_negative_soc = float(soc[crossing[0]]) if len(crossing) else np.nan
            summary_rows.append(
                {
                    "Electrolyte_set": set_name,
                    "C_rate": float(rate),
                    "Cutoff_SOC_pct": max_soc,
                    "Cutoff_time_min": float(time_min[-1]),
                    "Minimum_anode_potential_mV": float(np.min(potential_mV)),
                    "First_SOC_below_0V_pct": first_negative_soc,
                    "End_voltage_V": float(voltage[-1]),
                }
            )
            for i in range(len(soc)):
                curve_rows.append(
                    {
                        "Electrolyte_set": set_name,
                        "C_rate": float(rate),
                        "SOC_pct": float(soc[i]),
                        "Elapsed_time_min": float(time_min[i]),
                        "Cell_voltage_V": float(voltage[i]),
                        "Anode_potential_mV": float(potential_mV[i]),
                    }
                )
        maps[set_name] = potential_map

    difference = maps["ORegan2022"] - maps["Ai2020 corrected"]
    pd.DataFrame(curve_rows).to_csv(
        OUT / "anode_potential_curves.csv", index=False, encoding="utf-8-sig"
    )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT / "cutoff_and_minimum_summary.csv", index=False, encoding="utf-8-sig")

    for set_name, values in maps.items():
        pd.DataFrame(values, index=RATES, columns=SOC_GRID).to_csv(
            OUT / f"{set_name.lower().replace(' ', '_')}_map_mV.csv",
            encoding="utf-8-sig",
        )

    fig, axes = plt.subplots(1, 3, figsize=(18.0, 5.7), constrained_layout=True)
    main_levels = np.linspace(-100.0, 100.0, 41)
    diff_levels = np.linspace(-40.0, 40.0, 33)
    titles = {
        "ORegan2022": "Current model: ORegan2022 electrolyte",
        "Ai2020 corrected": "Ai2020 electrolyte",
    }
    for ax, set_name in zip(axes[:2], SETS):
        values = maps[set_name]
        contour = ax.contourf(
            SOC_GRID, RATES, values, levels=main_levels, cmap="RdBu_r", extend="both"
        )
        ax.contour(
            SOC_GRID,
            RATES,
            values,
            levels=[0.0],
            colors="black",
            linewidths=2.0,
        )
        selected = summary[summary.Electrolyte_set == set_name]
        ax.plot(
            selected.Cutoff_SOC_pct,
            selected.C_rate,
            color="#F0C808",
            lw=2.3,
            label="4.2 V cutoff",
        )
        ax.set_title(titles[set_name], weight="bold")
        ax.set_xlabel("Full-cell SOC (%)")
        ax.set_ylabel("Charge C-rate")
        ax.set_xlim(0.0, 80.0)
        ax.set_ylim(0.25, 6.0)
        ax.legend(loc="upper right", framealpha=0.9)
        fig.colorbar(contour, ax=ax, label="Anode potential (mV vs. Li/Li$^+$)")

    ax = axes[2]
    contour = ax.contourf(
        SOC_GRID, RATES, difference, levels=diff_levels, cmap="PiYG", extend="both"
    )
    ax.contour(
        SOC_GRID,
        RATES,
        difference,
        levels=[0.0],
        colors="black",
        linewidths=1.2,
    )
    ax.set_title("Difference: ORegan2022 - Ai2020", weight="bold")
    ax.set_xlabel("Full-cell SOC (%)")
    ax.set_ylabel("Charge C-rate")
    ax.set_xlim(0.0, 80.0)
    ax.set_ylim(0.25, 6.0)
    fig.colorbar(contour, ax=ax, label="Anode-potential difference (mV)")
    fig.suptitle(
        "DFN-predicted separator-side anode potential during CC charge\n"
        "All parameters fixed except electrolyte transport; $b_n=2.5$, retained $k_n$",
        fontsize=15,
        weight="bold",
    )
    figure = OUT / "anode_potential_soc_crate_electrolyte_comparison.png"
    fig.savefig(figure, dpi=240, bbox_inches="tight")
    plt.close(fig)

    pivot = summary[summary.C_rate.isin([0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0])]
    print("\nSELECTED SUMMARY")
    print(pivot.to_string(index=False))
    print(figure)


if __name__ == "__main__":
    main()
