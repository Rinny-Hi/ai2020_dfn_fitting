"""Screen electrode Bruggeman coefficients with both experimental EIS kinetics fixed."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import electrode_hysteresis_quantification as ehq
import eis_first_kinetics_application as eis
import gitt_ocp_analysis as ga
import mass_independent_nominal_anode_trial as mit
import nominal_radius_current_configuration as current
import ocp_grounded_candidate_search as common
import old_notebook_method_comparison as omc
import ocp_rn_sweep_fixed_transport as rs
import restored_brugg_de_csmax_audit as audit


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "results" / "260920_mass_independent_nominal_anode"
RESULTS = ROOT / "results" / "260921_eis_bruggeman_screen"
RESULTS.mkdir(parents=True, exist_ok=True)

BN_VALUES = (1.5, 2.0, 2.5, 2.914, 3.5)
BP_VALUES = (1.5, 1.83, 2.0, 2.5, 3.0)
BS_FIXED = 1.5


def set_bruggeman(model, bn, bp, bs=BS_FIXED):
    changed = dict(model)
    changed["params"] = model["params"].copy()
    changed["params"].update(
        {
            "Negative electrode Bruggeman coefficient (electrolyte)": bn,
            "Positive electrode Bruggeman coefficient (electrolyte)": bp,
            "Separator Bruggeman coefficient (electrolyte)": bs,
        },
        check_already_exists=False,
    )
    return changed


def pareto_mask(frame):
    values = frame[["mean_MAE_SOC10_70_mV", "mean_abs_capacity_error_mAh"]].to_numpy()
    keep = np.ones(len(values), dtype=bool)
    for i, point in enumerate(values):
        dominates = np.all(values <= point, axis=1) & np.any(values < point, axis=1)
        if np.any(dominates):
            keep[i] = False
    return keep


def main():
    bundle = ga.load_legacy_bundle(ga.get_paths(ROOT))
    base = audit.build_base(bundle, 1e-4, audit.CSN_CORRECTED)
    _, _, hybrid, cathode = rs.candidate_details(bundle, base)
    experiment = omc.load_old_dynamic_data(base["dynamic"], base["q_meas"])
    candidates = pd.read_csv(SOURCE / "mass_independent_ocp_candidate_metrics.csv")
    ocp_row = candidates[candidates.Candidate == current.NAME].iloc[0]
    stage = {key: float(ocp_row[key]) for key in ("x0", "x100", "y100", "y0")}
    effective = mit.model_for_effective_capacities(base, ocp_row)
    nominal = audit.set_dynamic_parameters(effective, current.RN_UM, current.RP_UM, 1e-4)
    eis_fixed = eis.with_kinetics(
        nominal,
        eis.KN_EIS_N4_ONLY / eis.AI2020_MREF,
        eis.KP_EIS / eis.AI2020_MREF,
    )
    anode = common.scaled_detail(hybrid, 0.50)
    cathode_dynamic = common.scaled_detail(cathode, 0.50)

    rows = []
    for bn in BN_VALUES:
        for bp in BP_VALUES:
            model = set_bruggeman(eis_fixed, bn, bp)
            case = f"bn={bn:g}, bp={bp:g}, bs={BS_FIXED:g}"
            print(case, flush=True)
            for rate in ga.RATES:
                for charge in (True, False):
                    direction = "Charge" if charge else "Discharge"
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
                    result = omc.old_time_metrics(experiment[(rate, charge)], sim)
                    capacity_error = result["end_time_error_min"] * rate * 2.28 / 60.0
                    rows.append(
                        {
                            "Case": case,
                            "bn": bn,
                            "bp": bp,
                            "bs": BS_FIXED,
                            "C_rate": rate,
                            "Direction": direction,
                            **result,
                            "capacity_error_mAh": capacity_error * 1000.0,
                        }
                    )

    branches = pd.DataFrame(rows)
    branches.to_csv(RESULTS / "eis_bruggeman_branch_metrics.csv", index=False, encoding="utf-8-sig")
    summary = branches.groupby(["Case", "bn", "bp", "bs"], as_index=False, sort=False).agg(
        mean_MAE_SOC10_70_mV=("MAE_SOC10_70_mV", "mean"),
        mean_abs_capacity_error_mAh=("capacity_error_mAh", lambda s: float(np.mean(np.abs(s)))),
        max_abs_capacity_error_mAh=("capacity_error_mAh", lambda s: float(np.max(np.abs(s)))),
    )
    charge = branches[branches.Direction == "Charge"].groupby("Case")["MAE_SOC10_70_mV"].mean()
    discharge = branches[branches.Direction == "Discharge"].groupby("Case")["MAE_SOC10_70_mV"].mean()
    summary["charge_mean_MAE_mV"] = summary.Case.map(charge)
    summary["discharge_mean_MAE_mV"] = summary.Case.map(discharge)
    summary["Pareto"] = pareto_mask(summary)
    summary.to_csv(RESULTS / "eis_bruggeman_summary.csv", index=False, encoding="utf-8-sig")

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2), constrained_layout=True)
    for ax, column, title, cmap in (
        (axes[0], "mean_MAE_SOC10_70_mV", "Mean voltage MAE (mV)", "viridis_r"),
        (axes[1], "mean_abs_capacity_error_mAh", "Mean absolute capacity error (mAh)", "magma_r"),
    ):
        pivot = summary.pivot(index="bn", columns="bp", values=column).sort_index(ascending=False)
        im = ax.imshow(pivot.to_numpy(), aspect="auto", cmap=cmap)
        ax.set_xticks(range(len(pivot.columns)), [f"{v:g}" for v in pivot.columns])
        ax.set_yticks(range(len(pivot.index)), [f"{v:g}" for v in pivot.index])
        ax.set_xlabel("Positive electrode Bruggeman bp")
        ax.set_ylabel("Negative electrode Bruggeman bn")
        ax.set_title(title)
        for i in range(len(pivot.index)):
            for j in range(len(pivot.columns)):
                ax.text(j, i, f"{pivot.iloc[i, j]:.1f}", ha="center", va="center", fontsize=8)
        fig.colorbar(im, ax=ax, shrink=0.82)
    fig.suptitle("Bruggeman screen with fixed EIS kn=7.40e-7, kp=3.12e-7, bs=1.5")
    fig.savefig(RESULTS / "eis_bruggeman_heatmaps.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    ax.scatter(
        summary.mean_abs_capacity_error_mAh,
        summary.mean_MAE_SOC10_70_mV,
        c=summary.bn,
        s=45 + 35 * summary.bp,
        cmap="viridis",
        alpha=0.8,
    )
    pareto = summary[summary.Pareto]
    ax.scatter(
        pareto.mean_abs_capacity_error_mAh,
        pareto.mean_MAE_SOC10_70_mV,
        facecolors="none",
        edgecolors="red",
        s=180,
        linewidths=1.8,
        label="Pareto candidates",
    )
    for row in pareto.itertuples():
        ax.annotate(f"bn={row.bn:g}, bp={row.bp:g}", (row.mean_abs_capacity_error_mAh, row.mean_MAE_SOC10_70_mV), xytext=(5, 5), textcoords="offset points", fontsize=8)
    ax.set_xlabel("Mean absolute capacity error (mAh)")
    ax.set_ylabel("Mean voltage MAE, 10-70% SOC (mV)")
    ax.grid(alpha=0.2)
    ax.legend()
    ax.set_title("Voltage/capacity trade-off with both EIS kinetics fixed")
    fig.savefig(RESULTS / "eis_bruggeman_pareto.png", dpi=220)
    plt.close(fig)

    report = {
        "fixed": {
            "kn": eis.KN_EIS_N4_ONLY,
            "kp": eis.KP_EIS,
            "bs": BS_FIXED,
            "Rn_um": current.RN_UM,
            "Rp_um": current.RP_UM,
            "Dsn_m2_s": rs.DSN,
            "Dsp_m2_s": rs.DSP,
            "De_multiplier": 1e-4,
            "qOCV_MAE_mV": float(ocp_row.qOCV_MAE_2_98_mV),
        },
        "pareto": pareto.to_dict(orient="records"),
        "best_voltage": summary.nsmallest(5, "mean_MAE_SOC10_70_mV").to_dict(orient="records"),
        "best_capacity": summary.nsmallest(5, "mean_abs_capacity_error_mAh").to_dict(orient="records"),
    }
    (RESULTS / "eis_bruggeman_screen.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("\nPareto candidates")
    print(
        pareto[["bn", "bp", "mean_MAE_SOC10_70_mV", "mean_abs_capacity_error_mAh", "max_abs_capacity_error_mAh"]]
        .sort_values("mean_abs_capacity_error_mAh")
        .to_string(index=False)
    )
    print("\nBest voltage fits")
    print(summary.nsmallest(8, "mean_MAE_SOC10_70_mV")[["bn", "bp", "mean_MAE_SOC10_70_mV", "mean_abs_capacity_error_mAh"]].to_string(index=False))


if __name__ == "__main__":
    main()
