"""Build the final integrated S1--S7 result package from saved analyses."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "261013_final_integrated_s1_s7"
PREFIT = ROOT / "results" / "261013_cumulative_seven_stage_model_comparison"
S7 = ROOT / "results" / "261013_s7_pso_voltage_fitting"
SUBSET = ROOT / "results" / "261013_s7_subset_fit_comparison"

S7_LABELS = {
    "S7a": "kn + bn fit (recommended)",
    "S7b": "Dsn + kn + bn fit",
    "S7c": "kn + bn + bp fit",
    "S7d": "kn + kp + bn + bp fit",
    "S7e": "four-parameter fit, Ai/Ai start",
}
SELECTION = {
    "S1": "construction path",
    "S2": "construction path",
    "S2b": "construction path",
    "S3": "construction path",
    "S4": "construction path",
    "S5": "construction path",
    "S6": "prefit baseline",
    "S7a": "recommended final model",
    "S7b": "not selected: Dsn low sensitivity",
    "S7c": "not selected: marginal voltage gain",
    "S7d": "exploratory: kn-kp non-identifiable",
    "S7e": "exploratory: initialization dependent",
}


def read_summary() -> pd.DataFrame:
    prefit = pd.read_csv(PREFIT / "seven_stage_summary.csv")
    prefit["Selection"] = prefit.Stage.map(SELECTION)
    prefit["Objective_RMSE_mV"] = np.nan

    s7_base = pd.read_csv(S7 / "s6_s7_summary.csv")
    s7_base = s7_base[s7_base.Scenario.isin(["S7a", "S7b"])].copy()
    s7_base["Stage"] = s7_base.Scenario
    s7_base["Label"] = s7_base.Stage.map(S7_LABELS)
    saved_objectives = pd.read_csv(S7 / "s6_s7_parameters.csv").set_index("Scenario")
    s7_base["Objective_RMSE_mV"] = s7_base.Stage.map(
        saved_objectives.Objective_RMSE_mV.to_dict()
    )

    subset = pd.read_csv(SUBSET / "summary.csv")
    subset = subset[subset.Scenario.isin(["S7c_bn_kn_bp", "S7d_bn_kp_kn_bp"])].copy()
    subset["Stage"] = subset.Scenario.map(
        {"S7c_bn_kn_bp": "S7c", "S7d_bn_kp_kn_bp": "S7d"}
    )
    subset["Label"] = subset.Stage.map(S7_LABELS)

    s7e = pd.read_csv(SUBSET / "s7e_ai2020_start_summary.csv")
    s7e["Stage"] = "S7e"
    s7e["Label"] = S7_LABELS["S7e"]
    s7e_params = pd.read_csv(SUBSET / "s7e_ai2020_start_parameters.csv").iloc[0]
    s7e["Objective_RMSE_mV"] = float(s7e_params.Objective_RMSE_mV)

    s7_frames = []
    for frame in (s7_base, subset, s7e):
        local = frame.copy()
        local["Successful_branches"] = 6
        local["Failed_branches"] = 0
        local["Selection"] = local.Stage.map(SELECTION)
        if "Objective_RMSE_mV" not in local:
            local["Objective_RMSE_mV"] = np.nan
        s7_frames.append(local)

    # Add mean absolute endpoint-capacity error from the integrated branch table later.
    columns = [
        "Stage",
        "Label",
        "Successful_branches",
        "Failed_branches",
        "Overall_full_time_RMSE_mV",
        "Overall_full_time_MAE_mV",
        "Overall_center_time_RMSE_mV",
        "Overall_center_time_MAE_mV",
        "Capacity_RMSE_pct",
        "Max_abs_capacity_error_mAh",
        "Objective_RMSE_mV",
        "Selection",
    ]
    combined = pd.concat(
        [prefit.reindex(columns=columns), *[frame.reindex(columns=columns) for frame in s7_frames]],
        ignore_index=True,
    )
    order = ["S1", "S2", "S2b", "S3", "S4", "S5", "S6", "S7a", "S7b", "S7c", "S7d", "S7e"]
    combined["_order"] = combined.Stage.map({name: i for i, name in enumerate(order)})
    return combined.sort_values("_order").drop(columns="_order").reset_index(drop=True)


def read_branches() -> pd.DataFrame:
    columns = [
        "Stage",
        "Label",
        "C_rate",
        "Direction",
        "Full_time_RMSE_mV",
        "Full_time_MAE_mV",
        "Center10_70_Qcell_time_RMSE_mV",
        "Center10_70_Qcell_time_MAE_mV",
        "Experiment_end_time_min",
        "DFN_end_time_min",
        "End_time_error_min",
        "Capacity_error_pct",
        "Capacity_error_mAh",
    ]
    prefit = pd.read_csv(PREFIT / "seven_stage_detail.csv").reindex(columns=columns)

    base = pd.read_csv(S7 / "s6_s7_branch_metrics.csv")
    base = base[base.Scenario.isin(["S7a", "S7b"])].copy()
    base["Stage"] = base.Scenario
    base["Label"] = base.Stage.map(S7_LABELS)

    subset = pd.read_csv(SUBSET / "branch_metrics.csv")
    subset = subset[subset.Scenario.isin(["S7c_bn_kn_bp", "S7d_bn_kp_kn_bp"])].copy()
    subset["Stage"] = subset.Scenario.map(
        {"S7c_bn_kn_bp": "S7c", "S7d_bn_kp_kn_bp": "S7d"}
    )
    subset["Label"] = subset.Stage.map(S7_LABELS)

    s7e = pd.read_csv(SUBSET / "s7e_ai2020_start_branch_metrics.csv")
    s7e["Stage"] = "S7e"
    s7e["Label"] = S7_LABELS["S7e"]
    return pd.concat(
        [prefit, base.reindex(columns=columns), subset.reindex(columns=columns), s7e.reindex(columns=columns)],
        ignore_index=True,
    )


def parameter_table() -> pd.DataFrame:
    dsn = 4.1707741459682473e-14
    dsp = 6.888175323807748e-14
    kn = 9.648533212e-7
    kp = 5.899519906708763e-7
    rows = [
        {"Stage": "S6", "Fitted_parameters": "none", "Start_basis": "kn Ai2020; kp EIS", "Dsn": dsn, "Dsp": dsp, "kn": kn, "kp": kp, "brugg_n": 2.5, "brugg_p": 1.83, "brugg_s": 1.5},
    ]
    saved = pd.read_csv(S7 / "s6_s7_parameters.csv").set_index("Scenario")
    rows.extend(
        [
            {"Stage": "S7a", "Fitted_parameters": "kn+brugg_n", "Start_basis": "kn Ai2020; kp EIS", "Dsn": float(saved.loc["S7a", "Dsn_m2_s"]), "Dsp": dsp, "kn": float(saved.loc["S7a", "kn_prefactor"]), "kp": kp, "brugg_n": float(saved.loc["S7a", "brugg_n"]), "brugg_p": 1.83, "brugg_s": 1.5},
            {"Stage": "S7b", "Fitted_parameters": "Dsn+kn+brugg_n", "Start_basis": "kn Ai2020; kp EIS", "Dsn": float(saved.loc["S7b", "Dsn_m2_s"]), "Dsp": dsp, "kn": float(saved.loc["S7b", "kn_prefactor"]), "kp": kp, "brugg_n": float(saved.loc["S7b", "brugg_n"]), "brugg_p": 1.83, "brugg_s": 1.5},
        ]
    )
    subset = pd.read_csv(SUBSET / "fitted_parameters.csv").set_index("Scenario")
    for source, stage in (("S7c_bn_kn_bp", "S7c"), ("S7d_bn_kp_kn_bp", "S7d")):
        row = subset.loc[source]
        rows.append(
            {"Stage": stage, "Fitted_parameters": str(row.Fitted_parameters), "Start_basis": "kn Ai2020; kp EIS", "Dsn": dsn, "Dsp": dsp, "kn": float(row.kn), "kp": float(row.kp), "brugg_n": float(row.brugg_n), "brugg_p": float(row.brugg_p), "brugg_s": 1.5}
        )
    e = pd.read_csv(SUBSET / "s7e_ai2020_start_parameters.csv").iloc[0]
    rows.append(
        {"Stage": "S7e", "Fitted_parameters": "kn+kp+brugg_n+brugg_p", "Start_basis": "kn Ai2020; kp Ai2020", "Dsn": dsn, "Dsp": dsp, "kn": float(e.kn), "kp": float(e.kp), "brugg_n": float(e.brugg_n), "brugg_p": float(e.brugg_p), "brugg_s": 1.5}
    )
    frame = pd.DataFrame(rows)
    frame["Selected"] = frame.Stage.eq("S7a")
    return frame


def make_plots(summary: pd.DataFrame, parameters: pd.DataFrame) -> None:
    prefit = summary[summary.Stage.isin(["S1", "S2", "S2b", "S3", "S4", "S5", "S6"])]
    fits = summary[summary.Stage.isin(["S6", "S7a", "S7b", "S7c", "S7d", "S7e"])]
    fig, axes = plt.subplots(2, 2, figsize=(14.2, 8.4), constrained_layout=True)
    axes[0, 0].plot(prefit.Stage, prefit.Overall_full_time_RMSE_mV, "o-", label="Full RMSE")
    axes[0, 0].plot(prefit.Stage, prefit.Overall_center_time_RMSE_mV, "s-", label="Center 10-70%")
    axes[0, 0].set_title("S1-S6 construction path: voltage")
    axes[0, 0].set_ylabel("RMSE (mV)")
    axes[0, 0].legend(frameon=False)
    axes[1, 0].plot(prefit.Stage, prefit.Capacity_RMSE_pct, "o-", color="#D55E00")
    axes[1, 0].set_title("S1-S6 construction path: capacity")
    axes[1, 0].set_ylabel("Capacity RMSE (%)")
    axes[1, 0].set_xlabel("Stage")
    axes[0, 1].bar(np.arange(len(fits)) - 0.18, fits.Overall_full_time_RMSE_mV, width=0.36, label="Full RMSE")
    axes[0, 1].bar(np.arange(len(fits)) + 0.18, fits.Overall_center_time_RMSE_mV, width=0.36, label="Center 10-70%")
    axes[0, 1].set_xticks(np.arange(len(fits)), fits.Stage)
    axes[0, 1].set_title("S6 and S7 alternatives: voltage")
    axes[0, 1].set_ylabel("RMSE (mV)")
    axes[0, 1].legend(frameon=False)
    colors = ["#999999" if stage == "S6" else "#0072B2" if stage == "S7a" else "#CC79A7" for stage in fits.Stage]
    axes[1, 1].bar(fits.Stage, fits.Capacity_RMSE_pct, color=colors)
    axes[1, 1].set_title("S6 and S7 alternatives: capacity")
    axes[1, 1].set_ylabel("Capacity RMSE (%)")
    axes[1, 1].set_xlabel("Stage")
    for ax in axes.flat:
        ax.grid(axis="y", alpha=0.22)
    fig.suptitle("Final integrated S1-S7 comparison", fontsize=15, weight="bold")
    fig.savefig(OUT / "final_integrated_error_comparison.png", dpi=240)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.8), constrained_layout=True)
    zoom_stages = {"S6", "S7a", "S7b", "S7c", "S7d", "S7e"}
    zoom_offsets = {
        "S6": (5, 4),
        "S7a": (5, -15),
        "S7b": (5, -3),
        "S7c": (5, 8),
        "S7d": (-24, 8),
        "S7e": (-25, -14),
    }
    for ax, zoom in zip(axes, (False, True)):
        source = (
            summary[summary.Stage.isin(zoom_stages)]
            if zoom
            else summary[summary.Stage.isin(["S1", "S2", "S2b", "S3", "S4", "S5", "S6"])]
        )
        for _, row in source.iterrows():
            selected = row.Stage == "S7a"
            ax.scatter(row.Overall_full_time_RMSE_mV, row.Capacity_RMSE_pct, s=100 if selected else 55, color="#0072B2" if selected else "#777777", zorder=3)
            offset = zoom_offsets.get(row.Stage, (5, 5)) if zoom else (5, 5)
            ax.annotate(row.Stage, (row.Overall_full_time_RMSE_mV, row.Capacity_RMSE_pct), xytext=offset, textcoords="offset points", weight="bold" if selected else "normal")
        ax.set_xlabel("Overall voltage RMSE (mV)")
        ax.grid(alpha=0.22)
    axes[0].set_ylabel("Capacity RMSE (%)")
    axes[0].set_title("S1-S6 construction path")
    axes[1].set_title("S6-S7 zoom (S7a highlighted)")
    axes[1].set_xlim(36.4, 39.0)
    axes[1].set_ylim(0.75, 3.2)
    fig.suptitle("Voltage-capacity trade-off", fontsize=15, weight="bold")
    fig.savefig(OUT / "final_voltage_capacity_tradeoff.png", dpi=240)
    plt.close(fig)

    normalized = parameters.copy()
    base = normalized[normalized.Stage == "S6"].iloc[0]
    for key in ("kn", "kp", "brugg_n", "brugg_p"):
        normalized[key] = normalized[key] / float(base[key])
    fig, ax = plt.subplots(figsize=(10.5, 5.8), constrained_layout=True)
    x = np.arange(len(normalized))
    width = 0.19
    for offset, key in zip((-1.5, -0.5, 0.5, 1.5), ("kn", "kp", "brugg_n", "brugg_p")):
        ax.bar(x + offset * width, normalized[key], width=width, label=key)
    ax.axhline(1.0, color="black", lw=1, ls="--")
    ax.set_xticks(x, normalized.Stage)
    ax.set_ylabel("Value / S6 value")
    ax.set_title("S7 fitted-parameter movement")
    ax.legend(frameon=False, ncol=4)
    ax.grid(axis="y", alpha=0.22)
    fig.savefig(OUT / "final_s7_parameter_movement.png", dpi=240)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    branches = read_branches()
    summary = read_summary()
    mean_abs = branches.groupby("Stage").Capacity_error_mAh.apply(lambda x: float(np.mean(np.abs(x)))).rename("Mean_abs_capacity_error_mAh")
    summary = summary.drop(columns=["Mean_abs_capacity_error_mAh"], errors="ignore").merge(mean_abs, on="Stage", how="left")
    parameters = parameter_table()
    summary.to_csv(OUT / "final_integrated_summary.csv", index=False, encoding="utf-8-sig")
    branches.to_csv(OUT / "final_integrated_branch_metrics.csv", index=False, encoding="utf-8-sig")
    parameters.to_csv(OUT / "final_s7_parameter_comparison.csv", index=False, encoding="utf-8-sig")
    make_plots(summary, parameters)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
