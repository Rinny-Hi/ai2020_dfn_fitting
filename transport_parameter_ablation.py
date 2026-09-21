"""Factorial ablation of Dsn, Dsp and electrolyte Bruggeman changes.

The OCP, stoichiometry window, effective capacities, and particle radii are
held fixed.  All 2^3 combinations are simulated against held-out 0.5C/1C/2C
charge and discharge data.  Exact Shapley allocations are computed so that
nonlinear interactions are assigned transparently.
"""

from __future__ import annotations

import itertools
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import mass_independent_nominal_anode_trial as mit
import ocp_grounded_candidate_search as common
import old_notebook_method_comparison as omc
import ocp_rn_sweep_fixed_transport as rs


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "results" / "260920_mass_independent_nominal_anode"
RESULTS = ROOT / "results" / "260921_transport_parameter_ablation"
RESULTS.mkdir(parents=True, exist_ok=True)

NAME = "Ai2020 nominal equilibrium + harvested hysteresis"
FACTORS = ("Dsn", "Dsp", "Bruggeman")
RN_UM = 5.0
RP_UM = 3.0


def label(enabled):
    if not enabled:
        return "Baseline Ai2020 transport"
    return "+".join(factor for factor in FACTORS if factor in enabled)


def update_case(model, enabled):
    changed = dict(model)
    changed["params"] = model["params"].copy()
    updates = {
        "Negative particle radius [m]": RN_UM * 1e-6,
        "Positive particle radius [m]": RP_UM * 1e-6,
    }
    if "Dsn" in enabled:
        updates["Negative particle diffusivity [m2.s-1]"] = rs.DSN
    if "Dsp" in enabled:
        updates["Positive particle diffusivity [m2.s-1]"] = rs.DSP
    if "Bruggeman" in enabled:
        updates.update({key: rs.BRUGGEMAN for key in rs.BRUGGEMAN_KEYS})
    changed["params"].update(updates, check_already_exists=False)
    return changed


def exact_shapley(values):
    """Return exact three-factor Shapley allocation for a scalar metric."""
    n = len(FACTORS)
    result = {}
    universe = set(FACTORS)
    for factor in FACTORS:
        contribution = 0.0
        others = universe - {factor}
        for size in range(n):
            for subset_tuple in itertools.combinations(sorted(others), size):
                subset = frozenset(subset_tuple)
                weight = math.factorial(size) * math.factorial(n - size - 1) / math.factorial(n)
                contribution += weight * (values[subset | {factor}] - values[subset])
        result[factor] = contribution
    return result


def main():
    bundle = ga.load_legacy_bundle(ga.get_paths(ROOT))
    base = ga.build_legacy_model_inputs(bundle)
    measured, nominal, hybrid, cathode = rs.candidate_details(bundle, base)
    experiment = omc.load_old_dynamic_data(base["dynamic"], base["q_meas"])
    metrics = pd.read_csv(SOURCE / "mass_independent_ocp_candidate_metrics.csv")
    metric_row = metrics[metrics.Candidate == NAME].iloc[0]
    effective = mit.model_for_effective_capacities(base, metric_row)
    stage = {key: float(metric_row[key]) for key in ("x0", "x100", "y100", "y0")}
    anode = common.scaled_detail(hybrid, 0.50)
    cathode_dynamic = common.scaled_detail(cathode, 0.50)

    subsets = []
    for size in range(len(FACTORS) + 1):
        subsets.extend(frozenset(v) for v in itertools.combinations(FACTORS, size))

    rows = []
    simulations = {}
    for enabled in subsets:
        case = label(enabled)
        model = update_case(effective, enabled)
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                print(f"{case} | {rate:g}C | {direction}", flush=True)
                sim = ehq.run_dfn(
                    rate,
                    charge,
                    stage,
                    model,
                    anode,
                    cathode_dynamic,
                    True,
                    True,
                )
                simulations[(case, rate, charge)] = sim
                result = omc.old_time_metrics(experiment[(rate, charge)], sim)
                capacity_error = result["end_time_error_min"] * rate * 2.28 / 60.0
                rows.append(
                    {
                        "Case": case,
                        "Dsn_changed": "Dsn" in enabled,
                        "Dsp_changed": "Dsp" in enabled,
                        "Bruggeman_changed": "Bruggeman" in enabled,
                        "C_rate": rate,
                        "Direction": direction,
                        **result,
                        "capacity_error_Ah": capacity_error,
                    }
                )
    validation = pd.DataFrame(rows)
    summary = validation.groupby(
        ["Case", "Dsn_changed", "Dsp_changed", "Bruggeman_changed"], as_index=False
    ).agg(
        Dynamic_MAE_SOC10_70_mV=("MAE_SOC10_70_mV", "mean"),
        Charge_MAE_SOC10_70_mV=("MAE_SOC10_70_mV", lambda x: float(x[validation.loc[x.index, "Direction"] == "Charge"].mean())),
        Discharge_MAE_SOC10_70_mV=("MAE_SOC10_70_mV", lambda x: float(x[validation.loc[x.index, "Direction"] == "Discharge"].mean())),
        Mean_abs_capacity_error_mAh=("capacity_error_Ah", lambda x: float(np.mean(np.abs(x))) * 1000.0),
        Max_abs_capacity_error_mAh=("capacity_error_Ah", lambda x: float(np.max(np.abs(x))) * 1000.0),
    )
    baseline = summary[summary.Case == label(frozenset())].iloc[0]
    for metric in (
        "Dynamic_MAE_SOC10_70_mV",
        "Charge_MAE_SOC10_70_mV",
        "Discharge_MAE_SOC10_70_mV",
        "Mean_abs_capacity_error_mAh",
        "Max_abs_capacity_error_mAh",
    ):
        summary[f"Delta_{metric}"] = summary[metric] - float(baseline[metric])

    key_to_subset = {
        row.Case: frozenset(
            factor
            for factor, flag in zip(
                FACTORS,
                (row.Dsn_changed, row.Dsp_changed, row.Bruggeman_changed),
            )
            if bool(flag)
        )
        for _, row in summary.iterrows()
    }
    shapley_rows = []
    for metric in (
        "Dynamic_MAE_SOC10_70_mV",
        "Charge_MAE_SOC10_70_mV",
        "Discharge_MAE_SOC10_70_mV",
        "Mean_abs_capacity_error_mAh",
    ):
        values = {
            key_to_subset[row.Case]: float(row[metric])
            for _, row in summary.iterrows()
        }
        allocations = exact_shapley(values)
        for factor, value in allocations.items():
            shapley_rows.append(
                {
                    "Metric": metric,
                    "Factor": factor,
                    "Shapley_contribution": value,
                    "Baseline": values[frozenset()],
                    "All_changed": values[frozenset(FACTORS)],
                    "Total_change": values[frozenset(FACTORS)] - values[frozenset()],
                }
            )
    shapley = pd.DataFrame(shapley_rows)

    validation.to_csv(RESULTS / "transport_ablation_validation.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(RESULTS / "transport_ablation_summary.csv", index=False, encoding="utf-8-sig")
    shapley.to_csv(RESULTS / "transport_ablation_shapley.csv", index=False, encoding="utf-8-sig")

    order = [
        "Baseline Ai2020 transport",
        "Dsn",
        "Dsp",
        "Bruggeman",
        "Dsn+Dsp",
        "Dsn+Bruggeman",
        "Dsp+Bruggeman",
        "Dsn+Dsp+Bruggeman",
    ]
    plot = summary.set_index("Case").loc[order].reset_index()
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.8), constrained_layout=True)
    axes[0].bar(plot.Case, plot.Dynamic_MAE_SOC10_70_mV, color="#4c78a8")
    axes[0].axhline(float(baseline.Dynamic_MAE_SOC10_70_mV), color="black", ls="--", lw=1)
    axes[0].set_ylabel("Dynamic MAE, SOC 10-70% (mV)")
    axes[0].set_title("All 2^3 transport-parameter combinations")
    axes[1].bar(plot.Case, plot.Mean_abs_capacity_error_mAh, color="#f58518")
    axes[1].axhline(float(baseline.Mean_abs_capacity_error_mAh), color="black", ls="--", lw=1)
    axes[1].set_ylabel("Mean absolute capacity error (mAh)")
    axes[1].set_title("Termination-capacity error")
    for ax in axes:
        ax.tick_params(axis="x", rotation=28)
        ax.grid(axis="y", alpha=0.2)
    fig.suptitle("Transport ablation | hybrid OCP | Rn=5 um, Rp=3 um")
    fig.savefig(RESULTS / "transport_ablation_summary.png", dpi=220)
    plt.close(fig)

    sh = shapley[shapley.Metric == "Dynamic_MAE_SOC10_70_mV"].copy()
    fig, ax = plt.subplots(figsize=(8.5, 5.2), constrained_layout=True)
    colors = ["#d62728" if value > 0 else "#2ca02c" for value in sh.Shapley_contribution]
    ax.bar(sh.Factor, sh.Shapley_contribution, color=colors)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_ylabel("Allocated change in dynamic MAE (mV)")
    ax.set_title("Exact Shapley allocation including interactions")
    ax.grid(axis="y", alpha=0.2)
    fig.savefig(RESULTS / "transport_ablation_shapley.png", dpi=220)
    plt.close(fig)

    curve_cases = ["Baseline Ai2020 transport", "Dsn", "Dsp", "Bruggeman", "Dsn+Dsp+Bruggeman"]
    colors = {
        curve_cases[0]: "#000000",
        curve_cases[1]: "#1f77b4",
        curve_cases[2]: "#d62728",
        curve_cases[3]: "#2ca02c",
        curve_cases[4]: "#ff7f0e",
    }
    fig, axes = plt.subplots(2, 3, figsize=(16, 8.6), constrained_layout=True)
    for col, rate in enumerate(ga.RATES):
        for row_index, charge in enumerate((True, False)):
            ax = axes[row_index, col]
            exp = experiment[(rate, charge)]
            ax.plot(exp["t_min"], exp["V"], color="#555555", lw=3.0, label="Experiment")
            for case in curve_cases:
                sim = simulations[(case, rate, charge)]
                ax.plot(sim["t_min"], sim["V"], color=colors[case], lw=1.5, label=case)
            ax.set_title(f"{rate:g}C {'Charge' if charge else 'Discharge'}")
            ax.set_xlabel("Time (min)")
            ax.set_ylabel("Voltage (V)")
            ax.grid(alpha=0.2)
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("One-at-a-time transport changes and all-changed case")
    fig.savefig(RESULTS / "transport_ablation_dynamic_curves.png", dpi=220)
    plt.close(fig)

    report = {
        "fixed": {
            "OCP": NAME,
            "Rn_um": RN_UM,
            "Rp_um": RP_UM,
            "stoichiometry_window": {key: float(metric_row[key]) for key in ("x0", "x100", "y100", "y0")},
        },
        "old_values": {
            "Dsn_m2_s": 3.9e-14,
            "Dsp_m2_s": 5.387e-15,
            "bruggeman_n_p_s": [2.914, 1.83, 1.5],
        },
        "new_values": {
            "Dsn_m2_s": rs.DSN,
            "Dsp_m2_s": rs.DSP,
            "bruggeman_n_p_s": [1.5, 1.5, 1.5],
        },
        "summary": summary.to_dict(orient="records"),
        "shapley": shapley.to_dict(orient="records"),
    }
    (RESULTS / "transport_ablation_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print("\nSUMMARY")
    print(plot[["Case", "Dynamic_MAE_SOC10_70_mV", "Charge_MAE_SOC10_70_mV", "Discharge_MAE_SOC10_70_mV", "Mean_abs_capacity_error_mAh", "Delta_Dynamic_MAE_SOC10_70_mV"]].to_string(index=False))
    print("\nSHAPLEY DYNAMIC MAE")
    print(sh[["Factor", "Shapley_contribution", "Total_change"]].to_string(index=False))
    print("\nSaved to", RESULTS)


if __name__ == "__main__":
    main()
