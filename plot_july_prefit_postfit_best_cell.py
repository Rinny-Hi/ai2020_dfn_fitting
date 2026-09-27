from __future__ import annotations

import sys
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
OUT = PROJECT / "outputs" / "260928_july_final_result_graph"
sys.path.insert(0, str(REPO))

import charge_physical_time_common as common
import gitt_ocp_analysis as ga
import priority_initial_state_protocol_recheck as priority


BEFORE = {
    "Dsn": 4.170774145968247e-14,
    "kn": 9.648533212e-7,
    "brugg_n": 2.914,
}
AFTER = {
    "Dsn": 7.466200473885764e-14,
    "kn": 1.4502357402082455e-6,
    "brugg_n": 2.914,
}


def physical_time_metric(experiment: dict, simulation: dict) -> tuple[float, bool, float]:
    view = common._cc_view(experiment)
    grid = np.linspace(0.0, float(view["t_min"][-1]), common.N_POINTS)
    sim_t = np.asarray(simulation["t_min"], float)
    sim_v = np.asarray(simulation["V"], float)
    margin = float(sim_t[-1] - grid[-1])
    if margin < -1e-9:
        return np.nan, False, margin
    observed = np.interp(grid, view["t_min"], view["V"])
    predicted = np.interp(grid, sim_t, sim_v)
    return float(np.sqrt(np.mean((1000.0 * (predicted - observed)) ** 2))), True, margin


def choose_best_cell(context, after_charge, after_discharge) -> tuple[str, pd.DataFrame]:
    rows = []
    for cell in context.curves:
        values = []
        feasible = True
        for charge, runs in ((True, after_charge), (False, after_discharge)):
            for rate in ga.RATES:
                rmse, ok, margin = physical_time_metric(
                    context.curves[cell][(rate, charge)], runs[rate]
                )
                feasible = feasible and ok
                if ok:
                    values.append(rmse)
                rows.append(
                    {
                        "Cell": cell,
                        "C_rate": rate,
                        "Direction": "Charge" if charge else "Discharge",
                        "Physical_time_RMSE_mV": rmse,
                        "Full_time_feasible": ok,
                        "Time_margin_min": margin,
                    }
                )
        rows.append(
            {
                "Cell": cell,
                "C_rate": "ALL",
                "Direction": "ALL",
                "Physical_time_RMSE_mV": float(np.sqrt(np.mean(np.asarray(values) ** 2)))
                if feasible
                else np.nan,
                "Full_time_feasible": feasible,
                "Time_margin_min": np.nan,
            }
        )
    detail = pd.DataFrame(rows)
    aggregate = detail[(detail.C_rate == "ALL") & detail.Full_time_feasible].sort_values(
        "Physical_time_RMSE_mV"
    )
    if aggregate.empty:
        raise RuntimeError("No cell is full-time feasible in all six after-fit branches")
    return str(aggregate.iloc[0].Cell), detail


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    context = common.build_context()
    before_charge = common.simulate_direction(context, BEFORE, True)
    before_discharge = common.simulate_direction(context, BEFORE, False)
    after_charge = common.simulate_direction(context, AFTER, True)
    after_discharge = common.simulate_direction(context, AFTER, False)
    best_cell, ranking = choose_best_cell(context, after_charge, after_discharge)
    ranking.to_csv(OUT / "july_experiment_cell_selection.csv", index=False, encoding="utf-8-sig")

    fig, axes = plt.subplots(2, 3, figsize=(16.2, 9.5), constrained_layout=True)
    rows = []
    for row, (charge, before_runs, after_runs) in enumerate(
        (
            (True, before_charge, after_charge),
            (False, before_discharge, after_discharge),
        )
    ):
        direction = "Charge" if charge else "Discharge"
        for col, rate in enumerate(ga.RATES):
            ax = axes[row, col]
            experiment = context.curves[best_cell][(rate, charge)]
            experiment_view = common._cc_view(experiment)
            before = before_runs[rate]
            after = after_runs[rate]
            before_rmse, before_feasible, before_margin = physical_time_metric(experiment, before)
            after_rmse, after_feasible, after_margin = physical_time_metric(experiment, after)
            before_text = f"{before_rmse:.1f}" if before_feasible else "early cutoff"
            after_text = f"{after_rmse:.1f}" if after_feasible else "early cutoff"
            rows.append(
                {
                    "Cell": best_cell,
                    "C_rate": rate,
                    "Direction": direction,
                    "Before_physical_time_RMSE_mV": before_rmse,
                    "After_physical_time_RMSE_mV": after_rmse,
                    "Before_full_time_feasible": before_feasible,
                    "After_full_time_feasible": after_feasible,
                    "Before_time_margin_min": before_margin,
                    "After_time_margin_min": after_margin,
                    "RMSE_improvement_mV": before_rmse - after_rmse
                    if before_feasible and after_feasible
                    else np.nan,
                }
            )
            ax.plot(
                experiment_view["t_min"],
                experiment_view["V"],
                color="black",
                lw=2.4,
                label=f"Experiment ({best_cell})",
                zorder=4,
            )
            ax.plot(
                before["t_min"],
                before["V"],
                color="#999999",
                lw=2.0,
                ls="--",
                label="Before fitting",
            )
            ax.plot(
                after["t_min"],
                after["V"],
                color="#0072B2",
                lw=2.3,
                label="After fitting",
            )
            ax.set(
                title=(
                    f"{rate:g}C {direction}\n"
                    f"Physical-time RMSE: {before_text} -> {after_text} mV"
                ),
                xlabel="Time [min]",
                ylabel="Voltage [V]",
            )
            ax.grid(alpha=0.23)
            if row == 0 and col == 0:
                ax.legend(fontsize=9)

    fig.suptitle(
        "July BoL: before and after fitting, best-matching experiment cell only",
        fontsize=15,
    )
    fig.text(
        0.5,
        0.002,
        "Before: Dsn=4.171e-14, kn=9.649e-7   |   "
        "After: Dsn=7.466e-14, kn=1.450e-6   |   brugg_n=2.914 fixed",
        ha="center",
        fontsize=10,
    )
    fig.savefig(OUT / "july_prefit_postfit_best_cell.png", dpi=220)
    plt.close(fig)
    pd.DataFrame(rows).to_csv(
        OUT / "july_prefit_postfit_best_cell_metrics.csv",
        index=False,
        encoding="utf-8-sig",
    )


if __name__ == "__main__":
    main()
