"""Compare every previously discussed negative-electrode kinetic value.

September-specific fitted values are intentionally excluded.  The two Rct/EIS
values are kept visibly separate from an unverified historical check value and
from model-fitted/literature reference values.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import gitt_ocp_analysis as ga
import prefit_physical_anchor_decision as p0
import priority_initial_state_protocol_recheck as priority

import analyze_oregan_brugg_2p5_anode_potential as anode_analysis
import electrolyte_transport_set_comparison as comparison
import oregan_bn25_kn_only_fit as fit
import oregan_combined_sensitivity_identifiability as combined


OUT = comparison.ROOT / "results" / "261012_oregan_bn25_historical_kn_comparison"

CANDIDATES = {
    "EIS N3+N4": {
        "kn": 2.48e-7,
        "status": "past experimental interpretation; rejected",
        "source": "anode half-cell Rct/EIS",
    },
    "Past 2.59e-7 check": {
        "kn": 2.59e-7,
        "status": "historical check; experimental provenance unconfirmed",
        "source": "past standalone comparison",
    },
    "EIS N4-only": {
        "kn": 7.40e-7,
        "status": "past experimental interpretation; rejected",
        "source": "anode half-cell Rct/EIS",
    },
    "Ai2020 anchor": {
        "kn": fit.BASE["kn"],
        "status": "current retained reference",
        "source": "Ai2020/Dualfoil nominal",
    },
    "Current combined fit": {
        "kn": fit.BASE["kn"] * 1.0992678196191192,
        "status": "diagnostic only; not adopted",
        "source": "current Dynamic+HPPC profile fit",
    },
    "July prior fit": {
        "kn": 1.4502357402082455e-6,
        "status": "past fitted effective value; not experimental",
        "source": "July charge fit with prior model configuration",
    },
}


def objective_metrics(result, grids, dynamic_exp, hppc_exp):
    dynamic_model = combined.dynamic_vector(result["runs"], grids, result["qcell"])
    dynamic_error = 1000.0 * (dynamic_model - dynamic_exp)
    pulse_exp = fit.hppc_pulse_vector(hppc_exp)
    pulse_model = fit.hppc_pulse_vector(result["hppc"])
    pulse_error = 1000.0 * (pulse_model - pulse_exp)
    absolute_error = 1000.0 * (result["hppc"] - hppc_exp)
    return {
        "Dynamic_objective_RMSE_mV": float(np.sqrt(np.mean(dynamic_error**2))),
        "HPPC_pulse_RMSE_mV": float(np.sqrt(np.mean(pulse_error**2))),
        "HPPC_absolute_RMSE_mV": float(np.sqrt(np.mean(absolute_error**2))),
        "HPPC_absolute_bias_mV": float(np.mean(absolute_error)),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    baseline, stage, anode, cathode, experiment, qcell = comparison.current_p0()
    transport = comparison.with_transport_set(baseline, "ORegan2022")
    hppc_targets = combined.load_hppc_targets()
    hppc_times = combined.hppc_target_times()
    hppc_exp = hppc_targets.Voltage_V.to_numpy(float)
    hppc_z0 = priority.monotone_voltage_inverse(
        transport["v_qocv"], transport["soc"], float(hppc_exp[0])
    )
    hppc_stage = priority.stage_at_soc(stage, hppc_z0, False)

    results = {}
    for name, metadata in CANDIDATES.items():
        print(name, metadata["kn"], flush=True)
        model = fit.apply_kn(transport, metadata["kn"] / fit.BASE["kn"])
        results[name] = {
            "model": model,
            "runs": combined.run_dynamic(model, stage, anode, cathode),
            "hppc": combined.run_hppc(
                model,
                hppc_stage,
                anode,
                cathode,
                hppc_times.Time_s.to_numpy(float),
            ),
            "qcell": qcell,
        }

    grids = combined.make_dynamic_grids(
        experiment, results["Ai2020 anchor"]["runs"], qcell
    )
    dynamic_exp = combined.dynamic_vector(experiment, grids, qcell, experiment=True)

    score_rows = []
    dynamic_rows = []
    hppc_rows = []
    anode_rows = []
    for name, metadata in CANDIDATES.items():
        result = results[name]
        score_rows.append(
            {
                "Scenario": name,
                "kn_PyBaMM_prefactor": metadata["kn"],
                "kn_table_m2p5_mol_m0p5_s": metadata["kn"] / fit.FARADAY,
                "Relative_to_Ai2020": metadata["kn"] / fit.BASE["kn"],
                "Status": metadata["status"],
                "Source": metadata["source"],
                **objective_metrics(result, grids, dynamic_exp, hppc_exp),
            }
        )
        dynamic_rows.extend(
            fit.evaluate_dynamic_metrics(name, result["runs"], experiment, qcell)
        )
        for index, target in hppc_targets.iterrows():
            model_voltage = float(result["hppc"][index])
            hppc_rows.append(
                {
                    "Scenario": name,
                    "Cycle": int(target.Cycle),
                    "Nominal_SOC": float(target.Nominal_SOC),
                    "Point": str(target.Point),
                    "Experimental_voltage_V": float(target.Voltage_V),
                    "Model_voltage_V": model_voltage,
                    "Error_mV": 1000.0 * (model_voltage - float(target.Voltage_V)),
                }
            )
        for rate in ga.RATES:
            for charge in (True, False):
                run = result["runs"][(rate, charge)]
                potential = np.asarray(run["potential_V"], dtype=float)
                index = int(np.nanargmin(potential))
                anode_rows.append(
                    {
                        "Scenario": name,
                        "C_rate": rate,
                        "Direction": "Charge" if charge else "Discharge",
                        "Minimum_anode_potential_mV": 1000.0 * float(potential[index]),
                        "Below_0V_time_fraction_pct": 100.0
                        * anode_analysis.time_fraction(run["time_s"], potential < 0.0),
                    }
                )

    scores = pd.DataFrame(score_rows)
    detail = pd.DataFrame(dynamic_rows)
    summary = p0.summarize(detail)
    branch = p0.branch_summary(detail, list(CANDIDATES))
    hppc_detail = pd.DataFrame(hppc_rows)
    anode_summary = pd.DataFrame(anode_rows)
    scores = scores.merge(
        summary[
            [
                "Scenario",
                "Full_RMSE_mV",
                "Center10_70_RMSE_mV",
                "Capacity_RMSE_pct",
            ]
        ],
        on="Scenario",
        how="left",
    )

    scores.to_csv(OUT / "historical_kn_scores.csv", index=False, encoding="utf-8-sig")
    detail.to_csv(OUT / "dynamic_detail.csv", index=False, encoding="utf-8-sig")
    branch.to_csv(OUT / "dynamic_branch_summary.csv", index=False, encoding="utf-8-sig")
    hppc_detail.to_csv(OUT / "hppc_keypoints.csv", index=False, encoding="utf-8-sig")
    anode_summary.to_csv(OUT / "anode_summary.csv", index=False, encoding="utf-8-sig")

    order = list(CANDIDATES)
    plot = scores.set_index("Scenario").loc[order]
    colors = ["#D55E00", "#CC79A7", "#009E73", "#E69F00", "#0072B2", "#6A3D9A"]
    fig, axes = plt.subplots(2, 2, figsize=(14.8, 9.0), constrained_layout=True)
    panels = (
        ("Dynamic_objective_RMSE_mV", "Dynamic fit-window RMSE (mV)"),
        ("Full_RMSE_mV", "Dynamic full-curve RMSE (mV)"),
        ("Capacity_RMSE_pct", "Dynamic capacity RMSE (%)"),
        ("HPPC_pulse_RMSE_mV", "HPPC pulse-change RMSE (mV)"),
    )
    for ax, (column, title) in zip(axes.ravel(), panels):
        bars = ax.bar(np.arange(len(order)), plot[column], color=colors)
        ax.bar_label(bars, fmt="%.1f", padding=3, fontsize=8.5)
        ax.set_xticks(np.arange(len(order)))
        ax.set_xticklabels(
            ["EIS\nN3+N4", "Past\n2.59", "EIS\nN4", "Ai2020", "Combined", "July\nfit"]
        )
        ax.set_title(title)
        ax.grid(axis="y", alpha=0.24)
    fig.suptitle(
        "Historical $k_n$ comparison: ORegan2022, $b_n=2.5$",
        fontsize=15,
        weight="bold",
    )
    fig.savefig(OUT / "historical_kn_metric_comparison.png", dpi=240, bbox_inches="tight")
    plt.close(fig)

    # HPPC pulse voltage changes, the portion most directly connected to kinetics.
    exp_matrix = hppc_exp.reshape(10, len(combined.HPPC_POINT_ORDER))
    exp_drop = 1000.0 * (exp_matrix[:, 1] - exp_matrix[:, 0])
    exp_rise = 1000.0 * (exp_matrix[:, 4] - exp_matrix[:, 3])
    soc = 100.0 * hppc_targets.drop_duplicates("Cycle").Nominal_SOC.to_numpy(float)
    fig, axes = plt.subplots(1, 2, figsize=(14.0, 5.6), constrained_layout=True)
    axes[0].plot(soc, exp_drop, "o-", color="black", lw=2.4, label="Experiment")
    axes[1].plot(soc, exp_rise, "o-", color="black", lw=2.4, label="Experiment")
    for (name, result), color in zip(results.items(), colors):
        matrix = np.asarray(result["hppc"]).reshape(10, len(combined.HPPC_POINT_ORDER))
        axes[0].plot(soc, 1000.0 * (matrix[:, 1] - matrix[:, 0]), "s--", color=color, label=name)
        axes[1].plot(soc, 1000.0 * (matrix[:, 4] - matrix[:, 3]), "s--", color=color, label=name)
    for ax, title in zip(axes, ("30 s discharge pulse drop", "10 s charge pulse rise")):
        ax.set_title(title)
        ax.set_xlabel("Nominal SOC (%)")
        ax.set_ylabel("Voltage change (mV)")
        ax.invert_xaxis()
        ax.grid(alpha=0.24)
    axes[1].legend(frameon=False, fontsize=8, ncol=2)
    fig.suptitle("HPPC pulse response by historical $k_n$", fontsize=15, weight="bold")
    fig.savefig(OUT / "historical_kn_hppc_pulse_comparison.png", dpi=240, bbox_inches="tight")
    plt.close(fig)

    print("\n", scores.to_string(index=False))
    print("\nSaved:", OUT)


if __name__ == "__main__":
    main()
