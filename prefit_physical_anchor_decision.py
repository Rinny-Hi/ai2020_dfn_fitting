"""Decide which independently anchored parameters improve the July BoL pre-fit.

This is a no-optimization ablation.  Every scenario uses the same C/50-selected
OCP and stoichiometry window, and is evaluated on the absolute-Ah axis against
all three July BoL cells.  The script separates physically defensible anchors
from legacy numerical comparators so a low error cannot silently become a
material-property claim.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybamm

import c50_selected_ocp_dynamic_validation as selected
import comprehensive_prefit_cross_cohort as cross
import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import kinetics_original_reversion_ablation as kinetics
import nominal_radius_current_configuration as current
import priority_initial_state_protocol_recheck as priority
import thickness_radius_six_case_comparison as radius_study
import three_thickness_measurement_comparison as thickness_study


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260927_prefit_physical_anchor_decision"
RESULTS.mkdir(parents=True, exist_ok=True)
AI_K = 1.0e-11 * 96485.33212


@dataclass(frozen=True)
class Scenario:
    name: str
    ln_um: float = 76.5
    lp_um: float = 68.0
    rn_um: float = 3.7542
    rp_um: float = 4.2387
    kn: float = AI_K
    kp: float = 4.18e-7
    dsn: str = "measured"
    dsp: str = "measured"
    positive_hysteresis: bool = False
    initialization: str = "rest"
    trust: str = "anchored"
    note: str = ""


def scenarios() -> list[Scenario]:
    return [
        Scenario("P0 anchored reference"),
        Scenario("P1 anchored + positive hysteresis", positive_hysteresis=True),
        Scenario("P2 anchored + endpoint start", initialization="endpoint"),
        Scenario("R0 Ai2020 radii", rn_um=5.0, rp_um=3.0),
        Scenario("L0 gauge thickness", ln_um=73.0, lp_um=62.5),
        Scenario("L1 SEM thickness", ln_um=81.52, lp_um=62.965),
        Scenario("L2 consensus thickness", ln_um=76.5, lp_um=62.7325),
        Scenario("D0 nominal Dsn", dsn="nominal"),
        Scenario("D1 nominal Dsp", dsp="nominal"),
        Scenario("D2 nominal Dsn+Dsp", dsn="nominal", dsp="nominal"),
        Scenario(
            "K0 legacy kp",
            kp=current.KP_PREF,
            trust="sensitivity-only",
            note="legacy kp is not the preferred cathode-EIS conversion",
        ),
        Scenario(
            "K1 Ai2020 kp",
            kp=AI_K,
            trust="literature comparator",
            note="Ai2020 cathode kinetics instead of cell-specific EIS",
        ),
        Scenario(
            "K2 legacy kn+kp",
            kn=current.KN_PREF,
            kp=current.KP_PREF,
            positive_hysteresis=True,
            trust="invalidated numerical comparator",
            note="anode-Rct kn was rejected because harvested anode was not reproducible",
        ),
        Scenario(
            "OLD old-notebook-like",
            rn_um=5.0,
            rp_um=3.0,
            kn=current.KN_PREF,
            kp=current.KP_PREF,
            positive_hysteresis=True,
            initialization="endpoint",
            trust="method comparator",
            note="approximates the earlier notebook choices; not a final physical set",
        ),
    ]


def build_model(model0: dict, s: Scenario) -> dict:
    model = thickness_study.same_loading_geometry(model0, s.ln_um, s.lp_um)
    model = radius_study.with_radius(model, s.rn_um, s.rp_um)
    model = kinetics.with_prefactors(model, s.kn, s.kp)
    if s.dsn == "nominal" or s.dsp == "nominal":
        raw = pybamm.ParameterValues("Ai2020")
        changed = dict(model)
        changed["params"] = model["params"].copy()
        updates = {}
        if s.dsn == "nominal":
            updates["Negative particle diffusivity [m2.s-1]"] = raw[
                "Negative particle diffusivity [m2.s-1]"
            ]
        if s.dsp == "nominal":
            updates["Positive particle diffusivity [m2.s-1]"] = raw[
                "Positive particle diffusivity [m2.s-1]"
            ]
        changed["params"].update(updates, check_already_exists=False)
        model = changed
    return model


def evaluate():
    old_curves = {}
    qc_rows = []
    for cell, path in cross.OLD_FILES.items():
        cell_curves, cell_rows = cross.parse_cell(
            "Old July BoL", cell, path, cross.OLD_BRANCH_STEPS
        )
        old_curves[cell] = cell_curves
        qc_rows.extend(cell_rows)
    qc = pd.DataFrame(qc_rows)
    protocol = cross.cohort_protocol_means(qc)
    model0, stage, anode, cathode, _, _, qcell = selected.selected_inputs()
    rows: list[dict] = []
    simulations = {}

    for s in scenarios():
        model = build_model(model0, s)
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                p = protocol.loc[("Old July BoL", rate, direction)]
                if s.initialization == "rest":
                    z0 = priority.monotone_voltage_inverse(
                        model["v_qocv"], model["soc"], float(p.Rest_end_V)
                    )
                    local_stage = priority.stage_at_soc(stage, z0, charge)
                else:
                    z0 = 0.0 if charge else 1.0
                    local_stage = stage
                sim_rate = float(p.Measured_current_A) / cross.NOMINAL_CAPACITY_AH
                print(f"{s.name} | {rate:g}C {direction}", flush=True)
                sim = ehq.run_dfn(
                    sim_rate,
                    charge,
                    local_stage,
                    model,
                    anode,
                    cathode,
                    False,
                    s.positive_hysteresis,
                    positive_initial_h=(-1.0 if charge else 1.0),
                )
                simulations[(s.name, rate, charge)] = sim
                for cell, branches in old_curves.items():
                    m = priority.curve_metrics(branches[(rate, charge)], sim, qcell, charge)
                    rows.append(
                        {
                            "Scenario": s.name,
                            "Trust": s.trust,
                            "Note": s.note,
                            "Cell": cell,
                            "C_rate": rate,
                            "Direction": direction,
                            "Ln_um": s.ln_um,
                            "Lp_um": s.lp_um,
                            "Rn_um": s.rn_um,
                            "Rp_um": s.rp_um,
                            "kn": s.kn,
                            "kp": s.kp,
                            "Dsn_source": s.dsn,
                            "Dsp_source": s.dsp,
                            "Positive_hysteresis": s.positive_hysteresis,
                            "Initialization": s.initialization,
                            **m,
                            "Capacity_error_mAh": 1000.0
                            * (m["Q_end_model_Ah"] - m["Q_end_exp_Ah"]),
                        }
                    )
    return pd.DataFrame(rows), simulations, old_curves


def summarize(detail: pd.DataFrame) -> pd.DataFrame:
    summary = (
        detail.groupby(["Scenario", "Trust", "Note"], as_index=False)
        .agg(
            Full_MAE_mV=("Full_MAE_mV", "mean"),
            Full_RMSE_mV=("Full_RMSE_mV", "mean"),
            Center10_70_MAE_mV=("Center10_70_MAE_mV", "mean"),
            Center10_70_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
            Capacity_RMSE_pct=(
                "Capacity_error_pct",
                lambda x: float(np.sqrt(np.mean(np.asarray(x, float) ** 2))),
            ),
            Mean_abs_capacity_error_mAh=(
                "Capacity_error_mAh", lambda x: float(np.mean(np.abs(x)))
            ),
            Max_abs_capacity_error_mAh=(
                "Capacity_error_mAh", lambda x: float(np.max(np.abs(x)))
            ),
        )
        .sort_values(["Capacity_RMSE_pct", "Center10_70_RMSE_mV"])
        .reset_index(drop=True)
    )
    center = summary["Center10_70_RMSE_mV"].to_numpy(float)
    capacity = summary["Capacity_RMSE_pct"].to_numpy(float)
    summary["Pareto_center_capacity"] = [
        not np.any(
            (center <= center[i])
            & (capacity <= capacity[i])
            & ((center < center[i]) | (capacity < capacity[i]))
        )
        for i in range(len(summary))
    ]
    return summary


def branch_summary(detail: pd.DataFrame, names: list[str]) -> pd.DataFrame:
    subset = detail[detail.Scenario.isin(names)]
    return (
        subset.groupby(["Scenario", "C_rate", "Direction"], as_index=False)
        .agg(
            Full_RMSE_mV=("Full_RMSE_mV", "mean"),
            Center10_70_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
            Center10_70_MAE_mV=("Center10_70_MAE_mV", "mean"),
            Capacity_error_pct=("Capacity_error_pct", "mean"),
            Capacity_error_mAh=("Capacity_error_mAh", "mean"),
        )
        .sort_values(["Scenario", "C_rate", "Direction"])
    )


def plot_tradeoff(summary: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(10.5, 7.0), constrained_layout=True)
    trust_colors = {
        "anchored": "#1f77b4",
        "literature comparator": "#2ca02c",
        "sensitivity-only": "#ff7f0e",
        "invalidated numerical comparator": "#d62728",
        "method comparator": "#9467bd",
    }
    for _, row in summary.iterrows():
        ax.scatter(
            row.Capacity_RMSE_pct,
            row.Center10_70_RMSE_mV,
            s=95 if row.Pareto_center_capacity else 55,
            color=trust_colors.get(row.Trust, "#777777"),
            marker="o" if row.Pareto_center_capacity else "x",
            zorder=3,
        )
        ax.annotate(
            row.Scenario.split()[0],
            (row.Capacity_RMSE_pct, row.Center10_70_RMSE_mV),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
        )
    ax.set_xlabel("Capacity RMSE [%]")
    ax.set_ylabel("10-70% voltage RMSE [mV]")
    ax.set_title("July BoL pre-fit: physical anchors and legacy comparators")
    ax.grid(alpha=0.25)
    fig.savefig(RESULTS / "prefit_anchor_tradeoff.png", dpi=220)
    plt.close(fig)


def plot_selected_curves(old_curves, simulations, names: list[str]):
    colors = ["#1f77b4", "#d62728", "#2ca02c"]
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.5), constrained_layout=True)
    for row, charge in enumerate((True, False)):
        for col, rate in enumerate(ga.RATES):
            ax = axes[row, col]
            for cell, branches in old_curves.items():
                curve = branches[(rate, charge)]
                ax.plot(curve["Q_Ah"], curve["V"], color="#333333", alpha=0.28, lw=1.0)
            for color, name in zip(colors, names):
                sim = simulations[(name, rate, charge)]
                ax.plot(sim["Q_Ah"], sim["V"], color=color, lw=1.7, label=name)
            ax.set_title(f"{rate:g}C {'Charge' if charge else 'Discharge'}")
            ax.set_xlabel("Transferred capacity [Ah]")
            ax.set_ylabel("Voltage [V]")
            ax.grid(alpha=0.2)
    axes[0, 0].plot([], [], color="#333333", alpha=0.6, label="July cells")
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("Best defensible pre-fit candidates on the absolute-Ah axis")
    fig.savefig(RESULTS / "selected_prefit_capacity_voltage.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.5), constrained_layout=True)
    for row, charge in enumerate((True, False)):
        for col, rate in enumerate(ga.RATES):
            ax = axes[row, col]
            for cell, branches in old_curves.items():
                curve = branches[(rate, charge)]
                ax.plot(curve["t_min"], curve["V"], color="#333333", alpha=0.28, lw=1.0)
            for color, name in zip(colors, names):
                sim = simulations[(name, rate, charge)]
                ax.plot(sim["t_min"], sim["V"], color=color, lw=1.7, label=name)
            ax.set_title(f"{rate:g}C {'Charge' if charge else 'Discharge'}")
            ax.set_xlabel("Time [min]")
            ax.set_ylabel("Voltage [V]")
            ax.grid(alpha=0.2)
    axes[0, 0].plot([], [], color="#333333", alpha=0.6, label="July cells")
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("July BoL pre-fit time-voltage comparison")
    fig.savefig(RESULTS / "selected_prefit_time_voltage.png", dpi=220)
    plt.close(fig)


def main():
    detail, simulations, old_curves = evaluate()
    summary = summarize(detail)
    selected_names = [
        "P0 anchored reference",
        "P2 anchored + endpoint start",
        "OLD old-notebook-like",
    ]
    branches = branch_summary(detail, selected_names)
    detail.to_csv(RESULTS / "prefit_anchor_detail.csv", index=False)
    summary.to_csv(RESULTS / "prefit_anchor_summary.csv", index=False)
    branches.to_csv(RESULTS / "selected_branch_summary.csv", index=False)
    plot_tradeoff(summary)
    plot_selected_curves(old_curves, simulations, selected_names)
    manifest = {
        "scope": "July BoL three-cell pre-fit; no dynamic optimization",
        "ocp": "C/50-selected nominal anode + old cathode",
        "axis": "absolute transferred capacity",
        "scenarios": [s.__dict__ for s in scenarios()],
    }
    (RESULTS / "analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print("\nSUMMARY\n", summary.to_string(index=False))
    print("\nSELECTED BRANCHES\n", branches.to_string(index=False))


if __name__ == "__main__":
    main()
