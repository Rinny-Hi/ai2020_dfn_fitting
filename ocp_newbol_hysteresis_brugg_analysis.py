from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybamm


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
REPO = Path(
    r"C:\Users\user\Documents\Codex\2026-09-18\https-github-com-rinny-hi-ai2020\work\ai2020_dfn_fitting"
)
OUT = PROJECT / "outputs" / "260927_ocp_newbol_hysteresis_brugg_comparison"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REPO))

import c50_ocp_source_combination_study as c50
import charge_physical_time_common as common
import comprehensive_prefit_cross_cohort as cross
import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import p0_five_parameter_sensitivity_area as legacy
import priority_initial_state_protocol_recheck as priority


FINAL = {
    "Dsn": 7.466200473885764e-14,
    "kn": 1.4502357402082455e-6,
    "brugg_n": 2.914,
}
SCENARIOS = {
    "Baseline": {"brugg_n": 2.914, "positive_hysteresis": False},
    "+Cathode hysteresis": {"brugg_n": 2.914, "positive_hysteresis": True},
    "b_n = 1.5": {"brugg_n": 1.5, "positive_hysteresis": False},
    "Hysteresis + b_n = 1.5": {"brugg_n": 1.5, "positive_hysteresis": True},
}
COLORS = {
    "Baseline": "#222222",
    "+Cathode hysteresis": "#0072B2",
    "b_n = 1.5": "#E69F00",
    "Hysteresis + b_n = 1.5": "#CC79A7",
}


def clean_qv(data: dict, qcell: float, charge: bool) -> tuple[np.ndarray, np.ndarray]:
    q = priority.transferred_capacity(data, qcell, charge)
    v = np.asarray(data["V"], float)
    order = np.argsort(q, kind="stable")
    q, v = q[order], v[order]
    q, index = np.unique(q, return_index=True)
    return q, v[index]


def make_ocp_figure() -> pd.DataFrame:
    source = REPO / "results" / "260921_c50_ocp_source_combinations"
    summary = pd.read_csv(source / "ocp_source_combination_summary.csv")
    summary.to_csv(OUT / "ocp_selection_metrics.csv", index=False, encoding="utf-8-sig")

    bundle = ga.load_legacy_bundle(ga.get_paths(REPO))
    model = ga.build_legacy_model_inputs(bundle)
    old_anode = c50.old_detail(ehq.build_anode_detail(bundle, model))
    old_cathode = c50.old_detail(ehq.build_cathode_detail(bundle, model))
    new_anode = c50.new_anode_detail()
    new_cathode, _ = c50.extract_new_cathode_delith()
    params = pybamm.ParameterValues("Ai2020")
    nominal_anode = c50.nominal_detail(params["Negative electrode OCP [V]"])
    nominal_cathode = c50.nominal_detail(params["Positive electrode OCP [V]"])
    anodes = {"Old measured": old_anode, "New measured": new_anode, "Ai2020 nominal": nominal_anode}
    cathodes = {"Old measured": old_cathode, "New measured": new_cathode, "Ai2020 nominal": nominal_cathode}

    selected = summary[(summary.anode_source == "nominal") & (summary.cathode_source == "old")].iloc[0]
    soc = np.linspace(0.0, 1.0, 1001)
    target = c50.fullcell_c50_target(soc)
    keys = [("nominal", "old"), ("old", "old"), ("nominal", "new"), ("nominal", "nominal")]
    raw_anodes = {"old": old_anode, "new": new_anode, "nominal": nominal_anode}
    raw_cathodes = {"old": old_cathode, "new": new_cathode, "nominal": nominal_cathode}

    fig, axes = plt.subplots(2, 2, figsize=(14.8, 10.2), constrained_layout=True)
    ax = axes[0, 0]
    for label, detail in anodes.items():
        chosen = label == "Ai2020 nominal"
        ax.plot(detail["grid"], detail["equilibrium"], lw=3.0 if chosen else 1.7,
                color="#0072B2" if chosen else None, label=label + ("  SELECTED" if chosen else ""))
    ax.axvspan(float(selected.x0), float(selected.x100), color="#0072B2", alpha=0.10,
               label="Selected operating window")
    ax.set(xlabel="Graphite stoichiometry x", ylabel="Anode OCP [V vs Li/Li+]",
           title="A. Anode OCP sources")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)
    ax.text(0.02, 0.03, "Nominal chosen: absolute stoichiometry anchor\n+ best valid C/50 full-cell pairing",
            transform=ax.transAxes, fontsize=9, va="bottom",
            bbox=dict(facecolor="white", alpha=0.85, edgecolor="none"))

    ax = axes[0, 1]
    for label, detail in cathodes.items():
        chosen = label == "Old measured"
        ax.plot(detail["grid"], detail["equilibrium"], lw=3.0 if chosen else 1.7,
                color="#D55E00" if chosen else None, label=label + ("  SELECTED" if chosen else ""))
    ax.axvspan(float(selected.y100), float(selected.y0), color="#D55E00", alpha=0.10,
               label="Selected operating window")
    ax.set(xlabel="Cathode stoichiometry y", ylabel="Cathode OCP [V vs Li/Li+]",
           title="B. Cathode OCP sources")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)
    ax.text(0.02, 0.03, "Old measured chosen: absolute anchor + both directions\nNew cathode reverse branch only 73.6% complete",
            transform=ax.transAxes, fontsize=9, va="bottom",
            bbox=dict(facecolor="white", alpha=0.85, edgecolor="none"))

    ax = axes[1, 0]
    ax.plot(soc * 100, target["qocv"], color="black", lw=3.0, label="Measured C/50 qOCV")
    styles = {
        ("nominal", "old"): ("#D55E00", "-", 2.6, "nominal / old  SELECTED"),
        ("old", "old"): ("#0072B2", "--", 1.8, "old / old"),
        ("nominal", "new"): ("#CC79A7", ":", 2.2, "nominal / new  INVALID anchor"),
        ("nominal", "nominal"): ("#009E73", "-.", 1.8, "nominal / nominal"),
    }
    for key in keys:
        row = summary[(summary.anode_source == key[0]) & (summary.cathode_source == key[1])].iloc[0]
        z = np.array([row.x0, row.x100, row.y100, row.y0], float)
        pred = c50.predict(z, soc, raw_anodes[key[0]], raw_cathodes[key[1]], "equilibrium")
        color, ls, lw, label = styles[key]
        ax.plot(soc * 100, pred, color=color, ls=ls, lw=lw, label=label)
    ax.set(xlabel="Full-cell SOC [%]", ylabel="Voltage [V]", title="C. C/50 qOCV reconstruction")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    eligible = summary.primary_eligible.astype(bool)
    ax.scatter(summary.loc[eligible, "qOCV_RMSE_2_98_mV"],
               summary.loc[eligible, "C50_branch_mean_MAE_10_90_mV"],
               s=75, color="#0072B2", alpha=0.8, label="Eligible")
    ax.scatter(summary.loc[~eligible, "qOCV_RMSE_2_98_mV"],
               summary.loc[~eligible, "C50_branch_mean_MAE_10_90_mV"],
               s=75, facecolor="none", edgecolor="#CC79A7", linewidth=1.8, label="Not eligible")
    short_name = {"old": "Old", "new": "New", "nominal": "Nom"}
    for _, row in summary.iterrows():
        label = short_name[str(row.anode_source)] + "/" + short_name[str(row.cathode_source)]
        ax.annotate(label, (row.qOCV_RMSE_2_98_mV, row.C50_branch_mean_MAE_10_90_mV),
                    xytext=(4, 3), textcoords="offset points", fontsize=8)
    ax.scatter([selected.qOCV_RMSE_2_98_mV], [selected.C50_branch_mean_MAE_10_90_mV],
               marker="*", s=280, color="#D55E00", edgecolor="black", zorder=5, label="Selected Nom/Old")
    ax.set(xlabel="qOCV RMSE, 2-98% [mV]", ylabel="C/50 charge/discharge mean MAE [mV]",
           title="D. Accuracy vs data validity")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.suptitle("OCP source decision: accuracy alone is not enough", fontsize=15)
    fig.savefig(OUT / "ocp_selection_evidence.png", dpi=220)
    plt.close(fig)
    return summary


def rest_analysis(context: common.Context) -> tuple[pd.DataFrame, pd.DataFrame]:
    qc_path = REPO / "results" / "260927_comprehensive_prefit_cross_cohort" / "branch_qc_all_cells.csv"
    qc = pd.read_csv(qc_path)
    qc["Abs_tail_slope_mV_per_min"] = qc.Rest_tail_slope_mV_per_min.abs()
    qc["Onset_jump_mV"] = 1000.0 * (qc.Onset_V - qc.Rest_end_V)
    qc["Abs_onset_jump_mV"] = qc.Onset_jump_mV.abs()
    qc["Initial_SOC_pct"] = [
        100.0 * priority.monotone_voltage_inverse(
            context.model_geometry["v_qocv"], context.model_geometry["soc"], float(v)
        )
        for v in qc.Rest_end_V
    ]
    qc.to_csv(OUT / "bol_rest_branch_audit_with_inferred_soc.csv", index=False, encoding="utf-8-sig")
    grouped = (
        qc.groupby(["Cohort", "C_rate", "Direction"], as_index=False)
        .agg(
            n_cells=("Cell", "count"),
            Rest_duration_min=("Rest_duration_min", "mean"),
            Mean_abs_tail_slope_mV_per_min=("Abs_tail_slope_mV_per_min", "mean"),
            SD_abs_tail_slope_mV_per_min=("Abs_tail_slope_mV_per_min", "std"),
            Mean_abs_onset_jump_mV=("Abs_onset_jump_mV", "mean"),
            SD_abs_onset_jump_mV=("Abs_onset_jump_mV", "std"),
            Mean_initial_SOC_pct=("Initial_SOC_pct", "mean"),
            SD_initial_SOC_pct=("Initial_SOC_pct", "std"),
        )
    )
    old = grouped[grouped.Cohort == "Old July BoL"].set_index(["C_rate", "Direction"])
    new = grouped[grouped.Cohort == "New September"].set_index(["C_rate", "Direction"])
    ratio_rows = []
    for key in old.index.intersection(new.index):
        ratio_rows.append({
            "C_rate": key[0],
            "Direction": key[1],
            "Tail_slope_reduction_factor_July_over_September":
                old.loc[key, "Mean_abs_tail_slope_mV_per_min"] / new.loc[key, "Mean_abs_tail_slope_mV_per_min"],
        })
    ratios = pd.DataFrame(ratio_rows)
    grouped = grouped.merge(ratios, on=["C_rate", "Direction"], how="left")
    grouped.to_csv(OUT / "bol_rest_cohort_summary.csv", index=False, encoding="utf-8-sig")

    labels = [f"{r:g}C\n{d[:3]}" for r in (0.5, 1.0, 2.0) for d in ("Charge", "Discharge")]
    order = [(r, d) for r in (0.5, 1.0, 2.0) for d in ("Charge", "Discharge")]
    x = np.arange(len(order))
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.3), constrained_layout=True)
    width = 0.36
    for offset, cohort, color, name in [(-width/2, "Old July BoL", "#999999", "July: 10 min"),
                                         (width/2, "New September", "#009E73", "September: 120 min")]:
        part = grouped[grouped.Cohort == cohort].set_index(["C_rate", "Direction"])
        slope = np.array([part.loc[k, "Mean_abs_tail_slope_mV_per_min"] for k in order])
        jump = np.array([part.loc[k, "Mean_abs_onset_jump_mV"] for k in order])
        socv = np.array([part.loc[k, "Mean_initial_SOC_pct"] for k in order])
        axes[0].bar(x + offset, slope, width, color=color, label=name)
        axes[1].bar(x + offset, jump, width, color=color, label=name)
        axes[2].bar(x + offset, socv, width, color=color, label=name)
    axes[0].set_yscale("log")
    axes[0].axhline(0.1, color="#D55E00", ls="--", lw=1.2, label="0.1 mV/min guide")
    axes[0].set(ylabel="|rest-tail slope| [mV/min]", title="A. Remaining relaxation at rest end")
    axes[1].set(ylabel="|first-load voltage - rest voltage| [mV]", title="B. Observed current-step jump")
    axes[2].set(ylabel="qOCV-inferred initial SOC [%]", title="C. Initial state implied by rest voltage")
    for ax in axes:
        ax.set_xticks(x, labels)
        ax.grid(axis="y", alpha=0.25)
    axes[0].legend(fontsize=8)
    axes[2].text(0.02, 0.97, "Different starting SOC by rate is real protocol state,\nnot removed by longer rest.",
                 transform=axes[2].transAxes, va="top", fontsize=9,
                 bbox=dict(facecolor="white", alpha=0.85, edgecolor="none"))
    fig.suptitle("Does the 2-hour rest improve initial-state definition?", fontsize=15)
    fig.savefig(OUT / "new_bol_rest_equilibrium_check.png", dpi=220)
    plt.close(fig)
    return qc, grouped


def simulate_case(
    context: common.Context,
    scenario: str,
    spec: dict,
    protocol_row: pd.Series,
    rate: float,
    charge: bool,
) -> dict:
    values = dict(FINAL)
    values["brugg_n"] = float(spec["brugg_n"])
    model = legacy.apply_parameters(context.model_geometry, common.complete_values(values, context))
    z0 = priority.monotone_voltage_inverse(model["v_qocv"], model["soc"], float(protocol_row.Rest_end_V))
    local_stage = priority.stage_at_soc(context.stage, z0, charge)
    sim_rate = float(protocol_row.Measured_current_A) / cross.NOMINAL_CAPACITY_AH
    print(f"{scenario} | {protocol_row.name if hasattr(protocol_row, 'name') else ''} | {rate:g}C | {'Charge' if charge else 'Discharge'}", flush=True)
    return ehq.run_dfn(
        sim_rate,
        charge,
        local_stage,
        model,
        context.anode,
        context.cathode,
        False,
        bool(spec["positive_hysteresis"]),
        positive_initial_h=(-1.0 if charge else 1.0),
    )


def physical_time_metric(experiment: dict, simulation: dict) -> tuple[float, bool, float]:
    view = common._cc_view(experiment)
    grid = np.linspace(0.0, float(view["t_min"][-1]), 100)
    sim_t = np.asarray(simulation["t_min"], float)
    sim_v = np.asarray(simulation["V"], float)
    margin = float(sim_t[-1] - grid[-1])
    if margin < -1e-9:
        return np.nan, False, margin
    exp_v = np.interp(grid, view["t_min"], view["V"])
    mod_v = np.interp(grid, sim_t, sim_v)
    return float(np.sqrt(np.mean((1000.0 * (mod_v - exp_v)) ** 2))), True, margin


def scenario_comparison(context: common.Context, rest_qc: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    july_protocol = context.protocol
    new_protocol = (
        rest_qc[rest_qc.Cohort == "New September"]
        .groupby(["C_rate", "Direction"], as_index=True)
        .agg(Rest_end_V=("Rest_end_V", "mean"), Measured_current_A=("Median_current_A", "mean"))
    )
    runs: dict[tuple[str, str, float, bool], dict] = {}
    detail_rows = []
    new_rows = []
    for scenario, spec in SCENARIOS.items():
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                prow = july_protocol.loc[("Old July BoL", rate, direction)]
                sim = simulate_case(context, scenario, spec, prow, rate, charge)
                runs[(scenario, "July", rate, charge)] = sim
                for cell, branches in context.curves.items():
                    metric = priority.curve_metrics(branches[(rate, charge)], sim, context.qcell, charge)
                    time_rmse, feasible, margin = physical_time_metric(branches[(rate, charge)], sim)
                    detail_rows.append({
                        "Scenario": scenario,
                        "Cohort": "Old July BoL",
                        "Cell": cell,
                        "C_rate": rate,
                        "Direction": direction,
                        "brugg_n": spec["brugg_n"],
                        "Positive_hysteresis": spec["positive_hysteresis"],
                        **metric,
                        "Physical_time_RMSE_mV": time_rmse,
                        "Physical_time_feasible": feasible,
                        "Time_margin_min": margin,
                    })

                nrow = new_protocol.loc[(rate, direction)]
                nsim = simulate_case(context, scenario, spec, nrow, rate, charge)
                runs[(scenario, "September", rate, charge)] = nsim
                observed = rest_qc[(rest_qc.Cohort == "New September") &
                                   (rest_qc.C_rate == rate) & (rest_qc.Direction == direction)]
                q_model = float(np.asarray(nsim["Q_Ah"], float)[-1])
                t_model = float(np.asarray(nsim["t_min"], float)[-1])
                for _, row in observed.iterrows():
                    q_exp = float(row.CC_capacity_Ah)
                    new_rows.append({
                        "Scenario": scenario,
                        "Cohort": "New September",
                        "Cell": row.Cell,
                        "C_rate": rate,
                        "Direction": direction,
                        "brugg_n": spec["brugg_n"],
                        "Positive_hysteresis": spec["positive_hysteresis"],
                        "Q_end_exp_Ah": q_exp,
                        "Q_end_model_Ah": q_model,
                        "Capacity_error_pct": 100.0 * (q_model - q_exp) / q_exp,
                        "CC_duration_exp_min": float(row.CC_duration_min),
                        "CC_duration_model_min": t_model,
                        "Duration_error_min": t_model - float(row.CC_duration_min),
                    })
    detail = pd.DataFrame(detail_rows)
    new_detail = pd.DataFrame(new_rows)
    detail.to_csv(OUT / "july_hysteresis_brugg_detailed_metrics.csv", index=False, encoding="utf-8-sig")
    new_detail.to_csv(OUT / "september_hysteresis_brugg_endpoint_metrics.csv", index=False, encoding="utf-8-sig")

    july_summary = (
        detail.groupby(["Scenario", "Direction"], sort=False, as_index=False)
        .agg(
            Full_RMSE_mV=("Full_RMSE_mV", "mean"),
            Center_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
            Mean_abs_initial_error_mV=("Initial_voltage_error_mV", lambda x: float(np.mean(np.abs(x)))),
            Capacity_RMSE_pct=("Capacity_error_pct", lambda x: float(np.sqrt(np.mean(np.asarray(x, float) ** 2)))),
            Physical_time_RMSE_mV=("Physical_time_RMSE_mV", lambda x: float(np.sqrt(np.nanmean(np.asarray(x, float) ** 2))) if np.any(np.isfinite(x)) else np.nan),
            Physical_time_all_feasible=("Physical_time_feasible", "all"),
            Minimum_time_margin_min=("Time_margin_min", "min"),
        )
    )
    new_summary = (
        new_detail.groupby(["Scenario", "Direction"], sort=False, as_index=False)
        .agg(
            Capacity_RMSE_pct=("Capacity_error_pct", lambda x: float(np.sqrt(np.mean(np.asarray(x, float) ** 2)))),
            Mean_abs_capacity_error_mAh=("Q_end_model_Ah", lambda x: np.nan),
            Duration_RMSE_min=("Duration_error_min", lambda x: float(np.sqrt(np.mean(np.asarray(x, float) ** 2)))),
        )
    )
    # Correct the mean absolute capacity error with paired model/experimental values.
    mae = (new_detail.assign(abs_err_mAh=1000.0 * (new_detail.Q_end_model_Ah - new_detail.Q_end_exp_Ah).abs())
           .groupby(["Scenario", "Direction"], sort=False).abs_err_mAh.mean())
    new_summary["Mean_abs_capacity_error_mAh"] = [mae.loc[(r.Scenario, r.Direction)] for r in new_summary.itertuples()]
    july_summary["Cohort"] = "Old July BoL"
    new_summary["Cohort"] = "New September (endpoint only)"
    july_summary.to_csv(OUT / "july_hysteresis_brugg_summary.csv", index=False, encoding="utf-8-sig")
    new_summary.to_csv(OUT / "september_hysteresis_brugg_summary.csv", index=False, encoding="utf-8-sig")

    # Curves on transferred-capacity axes: July has full raw curves and supports shape scoring.
    fig, axes = plt.subplots(3, 2, figsize=(14.5, 14.2), constrained_layout=True, sharex=False)
    for i, rate in enumerate(ga.RATES):
        for j, charge in enumerate((True, False)):
            ax = axes[i, j]
            for cell, branches in context.curves.items():
                qe, ve = clean_qv(branches[(rate, charge)], context.qcell, charge)
                ax.plot(qe, ve, color="#B5B5B5", lw=1.0, alpha=0.65)
            for scenario in SCENARIOS:
                qm, vm = clean_qv(runs[(scenario, "July", rate, charge)], context.qcell, charge)
                ax.plot(qm, vm, color=COLORS[scenario], lw=2.0, label=scenario)
            ax.set(xlabel="Transferred capacity [Ah]", ylabel="Voltage [V]",
                   title=f"{rate:g}C {'charge' if charge else 'discharge'}")
            ax.grid(alpha=0.25)
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("Hysteresis and negative-electrode Bruggeman comparison (July full curves)", fontsize=15)
    fig.savefig(OUT / "hysteresis_brugg_curve_comparison.png", dpi=220)
    plt.close(fig)

    # Summary metrics for July full curves and September endpoint capacity.
    fig, axes = plt.subplots(1, 3, figsize=(17.2, 5.5), constrained_layout=True)
    labels = list(SCENARIOS)
    x = np.arange(len(labels))
    width = 0.36
    for offset, direction, color in [(-width/2, "Charge", "#D55E00"), (width/2, "Discharge", "#0072B2")]:
        part = july_summary[july_summary.Direction == direction].set_index("Scenario")
        axes[0].bar(x + offset, [part.loc[s, "Full_RMSE_mV"] for s in labels], width, color=color, label=direction)
        axes[1].bar(x + offset, [part.loc[s, "Capacity_RMSE_pct"] for s in labels], width, color=color, label=direction)
        npart = new_summary[new_summary.Direction == direction].set_index("Scenario")
        axes[2].bar(x + offset, [npart.loc[s, "Capacity_RMSE_pct"] for s in labels], width, color=color, label=direction)
    axes[0].set(ylabel="Mean full-curve RMSE [mV]", title="A. July voltage shape")
    axes[1].set(ylabel="Cutoff capacity RMSE [%]", title="B. July endpoint (post hoc)")
    axes[2].set(ylabel="Cutoff capacity RMSE [%]", title="C. September endpoint only")
    short = ["Baseline", "+Hyst", "b=1.5", "Both"]
    for ax in axes:
        ax.set_xticks(x, short)
        ax.grid(axis="y", alpha=0.25)
    axes[0].legend(fontsize=8)
    fig.suptitle("Effect of cathode hysteresis and b_n=1.5", fontsize=15)
    fig.savefig(OUT / "hysteresis_brugg_metric_comparison.png", dpi=220)
    plt.close(fig)
    return july_summary, new_summary, runs


def write_report(
    ocp: pd.DataFrame,
    rest_summary: pd.DataFrame,
    july: pd.DataFrame,
    september: pd.DataFrame,
) -> None:
    selected = ocp[(ocp.anode_source == "nominal") & (ocp.cathode_source == "old")].iloc[0]
    oldold = ocp[(ocp.anode_source == "old") & (ocp.cathode_source == "old")].iloc[0]
    nomnew = ocp[(ocp.anode_source == "nominal") & (ocp.cathode_source == "new")].iloc[0]
    lines = [
        "# OCP 선택, 2시간 rest, hysteresis + Bruggeman 비교",
        "",
        "## 결론",
        "",
        "1. OCP는 **Ai2020 nominal graphite 음극 + 기존 측정 cathode 평형 OCP**를 유지한다. "
        f"유효 후보 중 qOCV RMSE가 {selected.qOCV_RMSE_2_98_mV:.2f} mV이고 C/50 branch MAE가 "
        f"{selected.C50_branch_mean_MAE_10_90_mV:.2f} mV다. 기존/기존 조합은 "
        f"{oldold.qOCV_RMSE_2_98_mV:.2f}/{oldold.C50_branch_mean_MAE_10_90_mV:.2f} mV로 더 나쁘다.",
        "2. nominal/new cathode 조합의 qOCV RMSE는 "
        f"{nomnew.qOCV_RMSE_2_98_mV:.2f} mV로 더 작지만, new cathode는 reverse branch가 73.6%만 완성되어 "
        "absolute/directional anchor가 없다. 그래서 수치만 보고 선택하지 않았다.",
        "3. September BoL의 120분 rest는 July 10분 rest보다 초기 SOC/OCP 정의를 **확실히 개선**한다. "
        "다만 0.5C charge 직전 저전압 rest에는 약 0.30 mV/min의 잔류 relaxation이 있어 완전 평형으로 간주하면 안 된다.",
        "4. `brugg_n=1.5`는 negative electrode electrolyte Bruggeman exponent만 2.914에서 1.5로 바꾼 경우다. "
        "positive=1.83, separator=1.5, fitted Dsn/kn은 고정했다. nominal anode에는 방향별 hysteresis 데이터가 없으므로 "
        "hysteresis 비교는 측정 방향성이 있는 cathode에만 적용했다.",
        "",
        "## 1. OCP 선택 근거",
        "",
        "- 음극: Ai2020 nominal은 절대 stoichiometry anchor가 있고 선택된 full-cell window와 직접 결합할 수 있다. "
        "새 음극 측정치는 양방향 형상은 있으나 절대 anchor가 없어 주 OCP로 쓰지 않았다.",
        "- 양극: 기존 측정 cathode는 절대 anchor와 lithiation/delithiation 양방향이 모두 있다. "
        "새 cathode는 역방향 미완성이라 diagnostic 후보로만 두었다.",
        "- 최종 선택은 단일 qOCV 오차가 아니라 C/50 charge/discharge branch 재현성과 데이터 완결성을 함께 본 결과다.",
        "",
        "## 2. 2시간 rest가 초기 상태 문제를 얼마나 줄였나",
        "",
        "| Rate | Direction | July slope | September slope | 감소 배수 | September inferred SOC |",
        "|---:|---|---:|---:|---:|---:|",
    ]
    rs = rest_summary.set_index(["Cohort", "C_rate", "Direction"])
    for rate in ga.RATES:
        for direction in ("Charge", "Discharge"):
            old = rs.loc[("Old July BoL", rate, direction)]
            new = rs.loc[("New September", rate, direction)]
            lines.append(
                f"| {rate:g}C | {direction} | {old.Mean_abs_tail_slope_mV_per_min:.3f} mV/min | "
                f"{new.Mean_abs_tail_slope_mV_per_min:.3f} mV/min | "
                f"{old.Mean_abs_tail_slope_mV_per_min/new.Mean_abs_tail_slope_mV_per_min:.1f}x | "
                f"{new.Mean_initial_SOC_pct:.2f}% |"
            )
    lines += [
        "",
        "- Charge 직전 rest 기준 기울기 감소는 0.5/1/2C에서 각각 약 8/29/90배다. "
        "Discharge 직전 고전압 rest도 약 30배 전후로 평탄해졌다.",
        "- 따라서 July의 큰 initial voltage error를 OCP fitting으로 흡수할 위험은 크게 줄었다. "
        "하지만 긴 rest가 cell inventory나 protocol history 차이를 없애지는 않는다. September는 rate별 rest 전압과 inferred SOC가 다르며, "
        "July보다 capacity도 작으므로 cohort-specific initial SOC/QLi/window가 여전히 필요하다.",
        "",
        "## 3. Hysteresis와 brugg_n=1.5 비교",
        "",
        "July는 원본 full voltage curve에 대해 RMSE를 계산했고, September는 현재 원본 파일이 연결되지 않아 저장된 read-only audit의 "
        "rest/current/cutoff 정보로 endpoint capacity만 비교했다. Capacity는 fitting 목적함수가 아니라 여기서도 사후 지표다.",
        "",
        "### July full-curve 결과",
        "",
        "| Scenario | Direction | Full RMSE | Center RMSE | Initial | Capacity RMSE | Physical-time feasible |",
        "|---|---|---:|---:|---:|---:|---|",
    ]
    for scenario in SCENARIOS:
        for direction in ("Charge", "Discharge"):
            row = july[(july.Scenario == scenario) & (july.Direction == direction)].iloc[0]
            lines.append(
                f"| {scenario} | {direction} | {row.Full_RMSE_mV:.2f} mV | {row.Center_RMSE_mV:.2f} mV | "
                f"{row.Mean_abs_initial_error_mV:.2f} mV | {row.Capacity_RMSE_pct:.2f}% | "
                f"{'yes' if row.Physical_time_all_feasible else 'no'} |"
            )
    lines += [
        "",
        "### September endpoint-only 결과",
        "",
        "| Scenario | Direction | Capacity RMSE | Mean abs capacity error | Duration RMSE |",
        "|---|---|---:|---:|---:|",
    ]
    for scenario in SCENARIOS:
        for direction in ("Charge", "Discharge"):
            row = september[(september.Scenario == scenario) & (september.Direction == direction)].iloc[0]
            lines.append(
                f"| {scenario} | {direction} | {row.Capacity_RMSE_pct:.2f}% | "
                f"{row.Mean_abs_capacity_error_mAh:.1f} mAh | {row.Duration_RMSE_min:.2f} min |"
            )
    b0c = july[(july.Scenario == "Baseline") & (july.Direction == "Charge")].iloc[0]
    b0d = july[(july.Scenario == "Baseline") & (july.Direction == "Discharge")].iloc[0]
    bothc = july[(july.Scenario == "Hysteresis + b_n = 1.5") & (july.Direction == "Charge")].iloc[0]
    bothd = july[(july.Scenario == "Hysteresis + b_n = 1.5") & (july.Direction == "Discharge")].iloc[0]
    lines += [
        "",
        "### 판단",
        "",
        f"- 결합 변경은 July charge full RMSE를 {b0c.Full_RMSE_mV:.2f}→{bothc.Full_RMSE_mV:.2f} mV, "
        f"discharge를 {b0d.Full_RMSE_mV:.2f}→{bothd.Full_RMSE_mV:.2f} mV로 바꾼다.",
        "- hysteresis는 rest 직후 offset을 일부 이동시키지만 전체 형상/중간 SOC가 같이 좋아지는지는 별개다. "
        "양극 old-GITT hysteresis를 1:1로 적용하는 것은 transfer test이지 확정 calibration이 아니다.",
        "- brugg_n=1.5는 electrolyte effective transport를 크게 높이는 구조적 변경이다. 기존 fitted Dsn/kn과 함께 바꾸면 "
        "파라미터 보상이 깨지므로, 성능이 좋아 보여도 새 bound 고정값으로 바로 채택하지 말고 재식별해야 한다.",
        "",
        "## 제한",
        "",
        "- September 원본 `E:\\예린\\260927 6-1/6-2 ...xlsx`가 현재 세션에서 마운트되지 않았다. "
        "2시간 rest 진단은 프로젝트에 이미 저장된 record-level QC 결과를 사용했다.",
        "- September의 scenario 결과는 full voltage RMSE가 아니라 cutoff endpoint 비교다. 원본이 다시 연결되면 동일한 physical-time grid로 full rerun해야 한다.",
        "- 모든 scenario는 current selected OCP와 final conservative Dsn/kn을 고정한 국소 비교다.",
    ]
    (OUT / "README_KO.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    print("Building context", flush=True)
    context = common.build_context()
    print("OCP comparison", flush=True)
    ocp = make_ocp_figure()
    print("Rest analysis", flush=True)
    rest_qc, rest_summary = rest_analysis(context)
    print("Scenario simulations", flush=True)
    july, september, _ = scenario_comparison(context, rest_qc)
    write_report(ocp, rest_summary, july, september)
    manifest = {
        "selected_ocp": "Ai2020 nominal graphite equilibrium + old measured cathode equilibrium",
        "final_parameters": FINAL,
        "scenario_definition": SCENARIOS,
        "new_bol_raw_available": False,
        "new_bol_rest_source": str(REPO / "results" / "260927_comprehensive_prefit_cross_cohort" / "branch_qc_all_cells.csv"),
    }
    (OUT / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Done: {OUT}", flush=True)


if __name__ == "__main__":
    main()
