"""Calibrate hysteresis magnitude and transition rate from full-cell C/20 history.

The equilibrium stoichiometry window is fixed to the endpoint<=10 mV solution.
Only legacy half-cell OCP and full-cell C/20 curves are used in the grid search.
The 0.5C/1C/2C curves remain held-out validation data.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import ocp_grounded_candidate_search as common
import old_notebook_method_comparison as omc


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260920_ocp_history_transition_search"
RESULTS.mkdir(parents=True, exist_ok=True)
SCALES = [0.35, 0.50, 0.65, 0.80, 1.00]
DECAYS = [0.5, 1.0, 2.0, 3.2, 5.0, 8.0, 12.0]


def history_initial_states(charge: bool) -> tuple[float, float]:
    """Start on the branch established by the preceding opposite direction."""
    target_n = -1.0 if charge else 1.0
    target_p = 1.0 if charge else -1.0
    return -target_n, -target_p


def interpolate_sim(sim: dict[str, np.ndarray], soc: np.ndarray) -> np.ndarray:
    order = np.argsort(sim["SOC"])
    x = np.asarray(sim["SOC"])[order]
    y = np.asarray(sim["V"])[order]
    unique, idx = np.unique(x, return_index=True)
    return np.interp(soc, unique, y[idx])


def c20_metrics(
    charge_sim: dict[str, np.ndarray],
    discharge_sim: dict[str, np.ndarray],
    measured_charge: np.ndarray,
    measured_discharge: np.ndarray,
    soc: np.ndarray,
    q_meas: float,
) -> dict[str, float | bool]:
    mask = (soc >= 0.02) & (soc <= 0.98)
    charge_pred = interpolate_sim(charge_sim, soc)
    discharge_pred = interpolate_sim(discharge_sim, soc)
    charge_error = (charge_pred - measured_charge) * 1000
    discharge_error = (discharge_pred - measured_discharge) * 1000
    coverage = (
        charge_sim["SOC"].min() <= 0.02
        and charge_sim["SOC"].max() >= 0.98
        and discharge_sim["SOC"].min() <= 0.02
        and discharge_sim["SOC"].max() >= 0.98
    )
    charge_capacity = float((charge_sim["SOC"].max() - charge_sim["SOC"].min()) * q_meas)
    discharge_capacity = float((discharge_sim["SOC"].max() - discharge_sim["SOC"].min()) * q_meas)
    return {
        "C20_charge_MAE_2_98_mV": float(np.mean(np.abs(charge_error[mask]))),
        "C20_discharge_MAE_2_98_mV": float(np.mean(np.abs(discharge_error[mask]))),
        "C20_branch_mean_MAE_2_98_mV": float(0.5 * (np.mean(np.abs(charge_error[mask])) + np.mean(np.abs(discharge_error[mask])))),
        "C20_charge_capacity_Ah": charge_capacity,
        "C20_discharge_capacity_Ah": discharge_capacity,
        "C20_full_SOC_coverage": bool(coverage),
    }


def run_case(
    rate: float,
    charge: bool,
    stage: dict[str, float],
    model: dict[str, Any],
    anode: dict[str, Any],
    cathode: dict[str, Any],
    scale: float,
    decay: float,
    history_start: bool,
) -> dict[str, np.ndarray]:
    n0, p0 = history_initial_states(charge) if history_start else (None, None)
    return ehq.run_dfn(
        rate,
        charge,
        stage,
        model,
        common.scaled_detail(anode, scale),
        common.scaled_detail(cathode, scale),
        True,
        True,
        negative_initial_h=n0,
        positive_initial_h=p0,
        negative_decay_rate=decay,
        positive_decay_rate=decay,
    )


def main() -> None:
    paths = ga.get_paths(ROOT)
    bundle = ga.load_legacy_bundle(paths)
    model = ga.build_legacy_model_inputs(bundle)
    anode = ehq.build_anode_detail(bundle, model)
    cathode = ehq.build_cathode_detail(bundle, model)
    measured_charge, measured_discharge = common.fullcell_c20_curves(bundle, model["soc"])

    prior = pd.read_csv(ROOT / "results" / "260920_ocp_grounded_candidate_search" / "candidate_integrated_summary.csv")
    base = prior[prior.Candidate == "C20-balanced, endpoint <= 10 mV"].iloc[0]
    stage = {"x0": float(base.x0), "x100": float(base.x100), "y100": float(base.y100), "y0": float(base.y0)}

    grid_rows: list[dict[str, Any]] = []
    c20_sims: dict[tuple[float, float, bool], dict[str, np.ndarray]] = {}
    for scale in SCALES:
        for decay in DECAYS:
            print(f"C/20 calibration | scale={scale:.2f} | decay={decay:g}", flush=True)
            charge_sim = run_case(0.05, True, stage, model, anode, cathode, scale, decay, True)
            discharge_sim = run_case(0.05, False, stage, model, anode, cathode, scale, decay, True)
            c20_sims[(scale, decay, True)] = charge_sim
            c20_sims[(scale, decay, False)] = discharge_sim
            row = {"scale": scale, "decay_rate": decay, **c20_metrics(charge_sim, discharge_sim, measured_charge, measured_discharge, model["soc"], model["q_meas"])}
            grid_rows.append(row)
    grid = pd.DataFrame(grid_rows)
    grid.to_csv(RESULTS / "history_transition_c20_grid.csv", index=False, encoding="utf-8-sig")

    eligible = grid[grid.C20_full_SOC_coverage].sort_values("C20_branch_mean_MAE_2_98_mV")
    # Validate the best three distinct C/20-calibrated candidates plus controls.
    selected_pairs: list[tuple[float, float, str, bool]] = []
    for _, row in eligible.head(3).iterrows():
        selected_pairs.append((float(row.scale), float(row.decay_rate), f"History C20 fit s={row.scale:.2f}, gamma={row.decay_rate:g}", True))
    selected_pairs.extend([
        (float(base.hysteresis_scale), 3.2, "Immediate-target common-scale control", False),
        (1.0, 3.2, "Immediate-target raw-gap control", False),
    ])

    exp = omc.load_old_dynamic_data(model["dynamic"], model["q_meas"])
    dynamic_rows: list[dict[str, Any]] = []
    validation_sims: dict[tuple[str, float, bool], dict[str, np.ndarray]] = {}
    for scale, decay, name, history_start in selected_pairs:
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                print(f"Validation | {name} | {rate:g}C | {direction}", flush=True)
                sim = run_case(rate, charge, stage, model, anode, cathode, scale, decay, history_start)
                validation_sims[(name, rate, charge)] = sim
                metric = omc.old_time_metrics(exp[(rate, charge)], sim)
                metric["capacity_error_Ah"] = metric["end_time_error_min"] * rate * 2.28 / 60.0
                dynamic_rows.append({"Candidate": name, "scale": scale, "decay_rate": decay, "history_start": history_start, "C_rate": rate, "Direction": direction, **metric})
    dynamic = pd.DataFrame(dynamic_rows)
    dynamic.to_csv(RESULTS / "history_transition_dynamic_validation.csv", index=False, encoding="utf-8-sig")

    dyn_summary = dynamic.groupby(["Candidate", "scale", "decay_rate", "history_start"], as_index=False).agg(
        Dynamic_MAE_SOC10_70_mV=("MAE_SOC10_70_mV", "mean"),
        Mean_abs_capacity_error_Ah=("capacity_error_Ah", lambda x: np.mean(np.abs(x))),
        Max_abs_capacity_error_Ah=("capacity_error_Ah", lambda x: np.max(np.abs(x))),
    )
    direction = dynamic.groupby(["Candidate", "Direction"], as_index=False).agg(
        MAE_SOC10_70_mV=("MAE_SOC10_70_mV", "mean"),
        Mean_capacity_error_Ah=("capacity_error_Ah", "mean"),
    )
    for label in ("Charge", "Discharge"):
        part = direction[direction.Direction == label].drop(columns="Direction").rename(columns={
            "MAE_SOC10_70_mV": f"{label}_MAE_SOC10_70_mV",
            "Mean_capacity_error_Ah": f"{label}_mean_capacity_error_Ah",
        })
        dyn_summary = dyn_summary.merge(part, on="Candidate")

    c20_lookup = grid.set_index(["scale", "decay_rate"])
    c20_values = []
    for _, row in dyn_summary.iterrows():
        if bool(row.history_start) and (row.scale, row.decay_rate) in c20_lookup.index:
            c20_values.append(float(c20_lookup.loc[(row.scale, row.decay_rate), "C20_branch_mean_MAE_2_98_mV"]))
        elif "common-scale" in row.Candidate:
            c20_values.append(float(base.C20_branch_mean_MAE_2_98_mV))
        else:
            raw = prior[prior.Candidate == "Endpoint + raw branch gap"].iloc[0]
            c20_values.append(float(raw.C20_branch_mean_MAE_2_98_mV))
    dyn_summary["C20_branch_mean_MAE_2_98_mV"] = c20_values
    dyn_summary["qOCV_MAE_2_98_mV"] = float(base.qOCV_MAE_2_98_mV)
    dyn_summary["qOCV_max_endpoint_abs_mV"] = float(base.qOCV_max_endpoint_abs_mV)
    dyn_summary["admissible"] = (
        (dyn_summary.C20_branch_mean_MAE_2_98_mV <= 20.0)
        & (dyn_summary.Mean_abs_capacity_error_Ah <= 0.05)
        & (dyn_summary.qOCV_max_endpoint_abs_mV <= 10.0 + 1e-6)
    )
    admitted = dyn_summary[dyn_summary.admissible].sort_values(["Dynamic_MAE_SOC10_70_mV", "C20_branch_mean_MAE_2_98_mV"])
    selected_name = admitted.iloc[0].Candidate if not admitted.empty else dyn_summary.sort_values("C20_branch_mean_MAE_2_98_mV").iloc[0].Candidate
    dyn_summary.to_csv(RESULTS / "history_transition_integrated_summary.csv", index=False, encoding="utf-8-sig")

    pivot = grid.pivot(index="decay_rate", columns="scale", values="C20_branch_mean_MAE_2_98_mV")
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
    image = axes[0].imshow(pivot.values, origin="lower", aspect="auto", cmap="viridis_r")
    axes[0].set_xticks(range(len(pivot.columns)), [f"{x:.2f}" for x in pivot.columns])
    axes[0].set_yticks(range(len(pivot.index)), [f"{x:g}" for x in pivot.index])
    axes[0].set(xlabel="Hysteresis scale", ylabel="Decay rate gamma", title="C/20 branch MAE [mV]")
    fig.colorbar(image, ax=axes[0], label="mV")
    for _, row in dyn_summary.iterrows():
        axes[1].scatter(row.C20_branch_mean_MAE_2_98_mV, row.Dynamic_MAE_SOC10_70_mV, s=90, label=row.Candidate)
        axes[2].scatter(row.Mean_abs_capacity_error_Ah, row.Dynamic_MAE_SOC10_70_mV, s=90, label=row.Candidate)
    axes[1].axvline(20, color="red", ls="--", lw=1.4)
    axes[2].axvline(0.05, color="red", ls="--", lw=1.4)
    axes[1].set(xlabel="C/20 branch MAE [mV]", ylabel="Held-out dynamic MAE [mV]", title="OCP calibration vs validation")
    axes[2].set(xlabel="Mean absolute capacity error [Ah]", ylabel="Held-out dynamic MAE [mV]", title="Capacity vs voltage validation")
    for ax in axes[1:]:
        ax.grid(alpha=0.25)
    axes[1].legend(fontsize=6.5)
    fig.suptitle("History-aware hysteresis calibration and held-out validation")
    fig.tight_layout()
    fig.savefig(RESULTS / "history_transition_tradeoffs.png", dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.2))
    for col, rate in enumerate(ga.RATES):
        for row_idx, charge in enumerate((True, False)):
            ax = axes[row_idx, col]
            direction_label = "Charge" if charge else "Discharge"
            metric = dynamic[(dynamic.Candidate == selected_name) & np.isclose(dynamic.C_rate, rate) & (dynamic.Direction == direction_label)].iloc[0]
            sim = validation_sims[(selected_name, rate, charge)]
            obs = exp[(rate, charge)]
            ax.plot(obs["t_min"], obs["V"], "k", lw=2.3, label="Experiment")
            ax.plot(sim["t_min"], sim["V"], color="#1f77b4", lw=2, label="Selected")
            ax.text(0.04, 0.06, f"MAE 10-70% = {metric.MAE_SOC10_70_mV:.1f} mV\nCapacity error = {metric.capacity_error_Ah:+.3f} Ah", transform=ax.transAxes, fontsize=9, bbox={"facecolor": "white", "alpha": 0.82})
            ax.set(title=f"{rate:g}C {direction_label}", xlabel="Time [min]", ylabel="Terminal voltage [V]")
            ax.grid(alpha=0.25)
            if row_idx == 0 and col == 0:
                ax.legend(fontsize=8)
    fig.suptitle(f"Held-out validation: {selected_name}")
    fig.tight_layout()
    fig.savefig(RESULTS / "selected_history_transition_dynamic_curves.png", dpi=200)
    plt.close(fig)

    report = {
        "selected_candidate": selected_name,
        "new_260918_GITT_used": False,
        "calibration": "legacy half-cell OCP + full-cell C/20 history",
        "validation": "0.5C/1C/2C held out",
        "fixed_window": stage,
        "fixed_delta_x": model["delta_x"],
        "fixed_delta_y": model["delta_y"],
        "admissibility": {"qOCV_endpoint_mV": 10, "C20_branch_MAE_mV": 20, "mean_capacity_error_Ah": 0.05},
    }
    (RESULTS / "history_transition_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print("\nBest C/20 grid points")
    print(eligible.head(10).to_string(index=False))
    print("\nIntegrated summary")
    print(dyn_summary.to_string(index=False))
    print("\nSelected:", selected_name)
    print("Saved to", RESULTS)


if __name__ == "__main__":
    main()
