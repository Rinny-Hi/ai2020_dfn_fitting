from __future__ import annotations

import json
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
OUT = PROJECT / "outputs" / "260928_july_kn_2p59e7_check"
sys.path.insert(0, str(REPO))

import charge_physical_time_common as common
import gitt_ocp_analysis as ga
import priority_initial_state_protocol_recheck as priority


SCENARIOS = {
    "July final": {
        "Dsn": 7.466200473885764e-14,
        "kn": 1.4502357402082455e-6,
        "brugg_n": 2.914,
    },
    "kn=2.59e-7": {
        "Dsn": 7.466200473885764e-14,
        "kn": 2.59e-7,
        "brugg_n": 2.914,
    },
}
BEST_CELL = "6-8"


def physical_time_metric(experiment: dict, simulation: dict) -> dict:
    view = common._cc_view(experiment)
    grid = np.linspace(0.0, float(view["t_min"][-1]), common.N_POINTS)
    sim_t = np.asarray(simulation["t_min"], float)
    sim_v = np.asarray(simulation["V"], float)
    margin = float(sim_t[-1] - grid[-1])
    overlap_end = min(float(grid[-1]), float(sim_t[-1]))
    overlap_grid = np.linspace(0.0, overlap_end, common.N_POINTS)
    observed_overlap = np.interp(overlap_grid, view["t_min"], view["V"])
    predicted_overlap = np.interp(overlap_grid, sim_t, sim_v)
    overlap_error = 1000.0 * (predicted_overlap - observed_overlap)
    out = {
        "Full_time_feasible": bool(margin >= -1e-9),
        "Time_margin_min": margin,
        "Overlap_time_RMSE_mV": float(np.sqrt(np.mean(overlap_error**2))),
        "Overlap_time_bias_mV": float(np.mean(overlap_error)),
    }
    if margin >= -1e-9:
        observed = np.interp(grid, view["t_min"], view["V"])
        predicted = np.interp(grid, sim_t, sim_v)
        error = 1000.0 * (predicted - observed)
        out.update(
            {
                "Physical_time_RMSE_mV": float(np.sqrt(np.mean(error**2))),
                "Physical_time_bias_mV": float(np.mean(error)),
            }
        )
    else:
        out.update({"Physical_time_RMSE_mV": np.nan, "Physical_time_bias_mV": np.nan})
    return out


def summarize(detail: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (scenario, direction), part in detail.groupby(["Scenario", "Direction"]):
        all_feasible = bool(part.Full_time_feasible.all())
        rows.append(
            {
                "Scenario": scenario,
                "Direction": direction,
                "All_full_time_feasible": all_feasible,
                "Minimum_time_margin_min": float(part.Time_margin_min.min()),
                "Physical_time_RMSE_mV": float(
                    np.sqrt(np.mean(part.Physical_time_RMSE_mV.to_numpy(float) ** 2))
                )
                if all_feasible
                else np.nan,
                "Overlap_time_RMSE_mV": float(
                    np.sqrt(np.mean(part.Overlap_time_RMSE_mV.to_numpy(float) ** 2))
                ),
                "Capacity_RMSE_pct": float(
                    np.sqrt(np.mean(part.Capacity_error_pct.to_numpy(float) ** 2))
                ),
                "Capacity_MAE_mAh": float(np.mean(np.abs(part.Capacity_error_mAh))),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    context = common.build_context()
    run_cache = {}
    rows = []
    for scenario, values in SCENARIOS.items():
        for charge in (True, False):
            direction = "Charge" if charge else "Discharge"
            runs = common.simulate_direction(context, values, charge)
            run_cache[(scenario, charge)] = runs
            for rate in ga.RATES:
                for cell, branches in context.curves.items():
                    experiment = branches[(rate, charge)]
                    simulation = runs[rate]
                    time_metric = physical_time_metric(experiment, simulation)
                    capacity_metric = priority.curve_metrics(
                        experiment, simulation, context.qcell, charge
                    )
                    rows.append(
                        {
                            "Scenario": scenario,
                            "Cell": cell,
                            "C_rate": rate,
                            "Direction": direction,
                            **time_metric,
                            "Capacity_aligned_full_RMSE_mV": capacity_metric["Full_RMSE_mV"],
                            "Experimental_cutoff_capacity_Ah": capacity_metric["Q_end_exp_Ah"],
                            "Model_cutoff_capacity_Ah": capacity_metric["Q_end_model_Ah"],
                            "Capacity_error_pct": capacity_metric["Capacity_error_pct"],
                            "Capacity_error_mAh": 1000.0
                            * (capacity_metric["Q_end_model_Ah"] - capacity_metric["Q_end_exp_Ah"]),
                        }
                    )
    detail = pd.DataFrame(rows)
    summary = summarize(detail)
    detail.to_csv(OUT / "kn_2p59e7_detail.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(OUT / "kn_2p59e7_summary.csv", index=False, encoding="utf-8-sig")

    fig, axes = plt.subplots(2, 3, figsize=(16.2, 9.5), constrained_layout=True)
    for row, charge in enumerate((True, False)):
        direction = "Charge" if charge else "Discharge"
        for col, rate in enumerate(ga.RATES):
            ax = axes[row, col]
            experiment = common._cc_view(context.curves[BEST_CELL][(rate, charge)])
            base_run = run_cache[("July final", charge)][rate]
            test_run = run_cache[("kn=2.59e-7", charge)][rate]
            test_row = detail[
                (detail.Scenario == "kn=2.59e-7")
                & (detail.Cell == BEST_CELL)
                & (detail.C_rate == rate)
                & (detail.Direction == direction)
            ].iloc[0]
            status = (
                f"RMSE {test_row.Physical_time_RMSE_mV:.1f} mV"
                if bool(test_row.Full_time_feasible)
                else f"early cutoff {abs(test_row.Time_margin_min):.2f} min"
            )
            ax.plot(experiment["t_min"], experiment["V"], color="black", lw=2.4, label=f"Experiment ({BEST_CELL})")
            ax.plot(base_run["t_min"], base_run["V"], color="#0072B2", lw=2.2, label="July final kn")
            ax.plot(test_run["t_min"], test_run["V"], color="#D55E00", lw=2.2, ls="--", label="kn=2.59e-7")
            ax.set(
                title=f"{rate:g}C {direction}\nkn=2.59e-7: {status}",
                xlabel="Time from 95% current [min]",
                ylabel="Voltage [V]",
            )
            ax.grid(alpha=0.24)
            if row == 0 and col == 0:
                ax.legend(fontsize=8)
    fig.suptitle("July BoL sensitivity check: lowering kn to 2.59e-7", fontsize=15)
    fig.savefig(OUT / "july_kn_2p59e7_comparison.png", dpi=220)
    plt.close(fig)

    manifest = {
        "fixed": {"Dsn": SCENARIOS["July final"]["Dsn"], "brugg_n": 2.914},
        "reference_kn": SCENARIOS["July final"]["kn"],
        "test_kn": SCENARIOS["kn=2.59e-7"]["kn"],
        "original_kn_lower_bound": 3.0e-7,
        "capacity_in_objective": False,
        "summary": summary.to_dict(orient="records"),
    }
    (OUT / "analysis_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
