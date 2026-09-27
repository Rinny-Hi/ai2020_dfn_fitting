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


REPO = Path(
    r"C:\Users\user\Documents\Codex\2026-09-18\https-github-com-rinny-hi-ai2020"
    r"\work\ai2020_dfn_fitting"
)
OUT = Path(
    r"C:\Users\user\Documents\Codex\2026-09-27"
    r"\github-260921-ocp-eis-current-progress\outputs\260927_ocp_initial_state_deep_dive"
)
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


def rmse(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(values**2)))


def charge_coulomb_start_soc(qc: pd.DataFrame, qcell: float) -> dict[tuple[float, str], float]:
    starts: dict[tuple[float, str], float] = {}
    for cell in sorted(qc.Cell.unique()):
        inv = 0.0
        starts[(0.5, cell)] = 0.0
        for rate, next_rate in ((0.5, 1.0), (1.0, 2.0)):
            charge = qc[
                (qc.Cell == cell)
                & np.isclose(qc.C_rate, rate)
                & (qc.Direction == "Charge")
            ].iloc[0]
            discharge = qc[
                (qc.Cell == cell)
                & np.isclose(qc.C_rate, rate)
                & (qc.Direction == "Discharge")
            ].iloc[0]
            inv += float(charge.CC_plus_CV_Ah) - float(discharge.CC_capacity_Ah)
            starts[(next_rate, cell)] = float(np.clip(inv / qcell, 0.0, 1.0))
    return starts


def source_details():
    paths = ga.get_paths(REPO)
    bundle = ga.load_legacy_bundle(paths)
    raw_model = ga.build_legacy_model_inputs(bundle)
    old_anode = c50.old_detail(ehq.build_anode_detail(bundle, raw_model))
    old_cathode = c50.old_detail(ehq.build_cathode_detail(bundle, raw_model))
    new_anode = c50.new_anode_detail()
    new_cathode, _ = c50.extract_new_cathode_delith()
    params = pybamm.ParameterValues("Ai2020")
    nominal_anode = c50.nominal_detail(params["Negative electrode OCP [V]"])
    nominal_cathode = c50.nominal_detail(params["Positive electrode OCP [V]"])
    return (
        {"old": old_anode, "new": new_anode, "nominal": nominal_anode},
        {"old": old_cathode, "new": new_cathode, "nominal": nominal_cathode},
    )


def candidate_initial_state_table(qc: pd.DataFrame, qcell: float, coulomb_soc: dict):
    summary = pd.read_csv(
        REPO
        / "results"
        / "260921_c50_ocp_source_combinations"
        / "ocp_source_combination_summary.csv"
    )
    anodes, cathodes = source_details()
    soc_grid = np.linspace(0.0, 1.0, 4001)
    detail_rows = []
    summary_rows = []
    charge_qc = qc[qc.Direction == "Charge"].copy()
    for row in summary.itertuples(index=False):
        key = f"{row.anode_source}/{row.cathode_source}"
        z = np.array([row.x0, row.x100, row.y100, row.y0], dtype=float)
        qocv = c50.predict(z, soc_grid, anodes[row.anode_source], cathodes[row.cathode_source])
        per_candidate = []
        for obs in charge_qc.itertuples(index=False):
            z_rest = priority.monotone_voltage_inverse(qocv, soc_grid, float(obs.Rest_end_V))
            z_coul = float(coulomb_soc[(float(obs.C_rate), str(obs.Cell))])
            v_coul = float(np.interp(z_coul, soc_grid, qocv))
            rec = {
                "candidate": key,
                "anode_source": row.anode_source,
                "cathode_source": row.cathode_source,
                "cell": obs.Cell,
                "C_rate": float(obs.C_rate),
                "rest_end_V": float(obs.Rest_end_V),
                "rest_tail_slope_mV_per_min": float(obs.Rest_tail_slope_mV_per_min),
                "rest_inverted_SOC_pct": 100.0 * z_rest,
                "coulomb_bookkeeping_SOC_pct": 100.0 * z_coul,
                "rest_minus_coulomb_SOC_pp": 100.0 * (z_rest - z_coul),
                "qOCV_at_coulomb_SOC_V": v_coul,
                "rest_minus_qOCV_at_coulomb_mV": 1000.0 * (float(obs.Rest_end_V) - v_coul),
                "onset_minus_rest_mV": 1000.0 * (float(obs.Onset_V) - float(obs.Rest_end_V)),
            }
            per_candidate.append(rec)
            detail_rows.append(rec)
        frame = pd.DataFrame(per_candidate)
        summary_rows.append(
            {
                "candidate": key,
                "anode_source": row.anode_source,
                "cathode_source": row.cathode_source,
                "qOCV_RMSE_2_98_mV": float(row.qOCV_RMSE_2_98_mV),
                "C50_branch_mean_MAE_10_90_mV": float(row.C50_branch_mean_MAE_10_90_mV),
                "primary_eligible": bool(row.primary_eligible),
                "anode_absolute_anchor": bool(row.anode_absolute_anchor),
                "cathode_absolute_anchor": bool(row.cathode_absolute_anchor),
                "anode_directional_complete": bool(row.anode_directional_complete),
                "cathode_directional_complete": bool(row.cathode_directional_complete),
                "mean_abs_rest_minus_coulomb_SOC_pp": float(
                    np.mean(np.abs(frame.rest_minus_coulomb_SOC_pp))
                ),
                "mean_abs_rest_minus_qOCV_at_coulomb_mV": float(
                    np.mean(np.abs(frame.rest_minus_qOCV_at_coulomb_mV))
                ),
            }
        )
    return pd.DataFrame(detail_rows), pd.DataFrame(summary_rows)


def initial_soc_lookup(
    scenario: str,
    rate: float,
    cell: str,
    qc_index: pd.DataFrame,
    protocol: pd.DataFrame,
    model: dict,
    coulomb_soc: dict,
) -> float:
    if scenario == "rest_mean":
        voltage = float(protocol.loc[("Old July BoL", rate, "Charge")].Rest_end_V)
        return priority.monotone_voltage_inverse(model["v_qocv"], model["soc"], voltage)
    if scenario in {"rest_cell", "rest_mean_cathode_hysteresis"}:
        voltage = float(qc_index.loc[(cell, rate, "Charge")].Rest_end_V)
        return priority.monotone_voltage_inverse(model["v_qocv"], model["soc"], voltage)
    if scenario == "endpoint":
        return 0.0
    if scenario == "coulomb_cell":
        return float(coulomb_soc[(rate, cell)])
    raise KeyError(scenario)


def simulate_initialization_scenarios(context: common.Context, qc: pd.DataFrame, coulomb_soc: dict):
    model = legacy.apply_parameters(context.model_geometry, common.complete_values(FINAL, context))
    qc_index = qc.set_index(["Cell", "C_rate", "Direction"])
    scenarios = [
        "rest_mean",
        "rest_cell",
        "endpoint",
        "coulomb_cell",
        "rest_mean_cathode_hysteresis",
    ]
    rows = []
    simulations = {}
    for scenario in scenarios:
        for rate in ga.RATES:
            for cell in sorted(context.curves):
                z0 = initial_soc_lookup(
                    scenario, float(rate), cell, qc_index, context.protocol, model, coulomb_soc
                )
                local_stage = priority.stage_at_soc(context.stage, z0, True)
                p = qc_index.loc[(cell, float(rate), "Charge")]
                sim_rate = float(p.Median_current_A) / cross.NOMINAL_CAPACITY_AH
                positive_hysteresis = scenario == "rest_mean_cathode_hysteresis"
                sim = ehq.run_dfn(
                    sim_rate,
                    True,
                    local_stage,
                    model,
                    context.anode,
                    context.cathode,
                    False,
                    positive_hysteresis,
                    positive_initial_h=(-1.0 if positive_hysteresis else None),
                )
                simulations[(scenario, float(rate), cell)] = sim
                grid = context.grids[(float(rate), cell)]
                exp = context.experimental_blocks[(float(rate), cell)]
                sim_t = np.asarray(sim["t_min"], float)
                sim_v = np.asarray(sim["V"], float)
                covered = float(min(1.0, sim_t[-1] / grid[-1]))
                pred = np.interp(np.minimum(grid, sim_t[-1]), sim_t, sim_v)
                err = 1000.0 * (pred - exp)
                n = len(err)
                early = err[: max(1, int(np.ceil(0.10 * n)))]
                middle = err[int(np.floor(0.10 * n)) : int(np.ceil(0.70 * n))]
                late = err[int(np.floor(0.70 * n)) :]
                qocv0 = float(np.interp(z0, model["soc"], model["v_qocv"]))
                rows.append(
                    {
                        "scenario": scenario,
                        "cell": cell,
                        "C_rate": float(rate),
                        "initial_SOC_pct": 100.0 * z0,
                        "rest_end_V": float(p.Rest_end_V),
                        "qOCV_initial_V": qocv0,
                        "experiment_onset_V": float(exp[0]),
                        "model_onset_V": float(pred[0]),
                        "qOCV_minus_rest_mV": 1000.0 * (qocv0 - float(p.Rest_end_V)),
                        "experiment_jump_from_rest_mV": 1000.0 * (float(exp[0]) - float(p.Rest_end_V)),
                        "model_jump_from_qOCV_mV": 1000.0 * (float(pred[0]) - qocv0),
                        "initial_error_mV": float(err[0]),
                        "full_RMSE_mV": rmse(err),
                        "full_bias_mV": float(np.mean(err)),
                        "early_0_10_RMSE_mV": rmse(early),
                        "early_0_10_bias_mV": float(np.mean(early)),
                        "middle_10_70_RMSE_mV": rmse(middle),
                        "middle_10_70_bias_mV": float(np.mean(middle)),
                        "late_70_100_RMSE_mV": rmse(late),
                        "late_70_100_bias_mV": float(np.mean(late)),
                        "time_coverage_fraction": covered,
                        "model_time_margin_min": float(sim_t[-1] - grid[-1]),
                    }
                )
    return pd.DataFrame(rows), simulations


def summaries(detail: pd.DataFrame):
    by_rate = (
        detail.groupby(["scenario", "C_rate"], as_index=False)
        .agg(
            initial_SOC_pct=("initial_SOC_pct", "mean"),
            initial_error_mV=("initial_error_mV", "mean"),
            full_RMSE_mV=("full_RMSE_mV", lambda x: rmse(np.asarray(x))),
            full_bias_mV=("full_bias_mV", "mean"),
            early_0_10_RMSE_mV=("early_0_10_RMSE_mV", lambda x: rmse(np.asarray(x))),
            middle_10_70_RMSE_mV=("middle_10_70_RMSE_mV", lambda x: rmse(np.asarray(x))),
            late_70_100_RMSE_mV=("late_70_100_RMSE_mV", lambda x: rmse(np.asarray(x))),
            min_time_coverage=("time_coverage_fraction", "min"),
            min_time_margin_min=("model_time_margin_min", "min"),
        )
    )
    overall = (
        detail.groupby("scenario", as_index=False)
        .agg(
            initial_error_RMSE_mV=("initial_error_mV", lambda x: rmse(np.asarray(x))),
            full_RMSE_mV=("full_RMSE_mV", lambda x: rmse(np.asarray(x))),
            early_0_10_RMSE_mV=("early_0_10_RMSE_mV", lambda x: rmse(np.asarray(x))),
            middle_10_70_RMSE_mV=("middle_10_70_RMSE_mV", lambda x: rmse(np.asarray(x))),
            late_70_100_RMSE_mV=("late_70_100_RMSE_mV", lambda x: rmse(np.asarray(x))),
            min_time_coverage=("time_coverage_fraction", "min"),
        )
    )
    return by_rate, overall


def make_plots(qc: pd.DataFrame, candidate_summary: pd.DataFrame, dynamic_by_rate: pd.DataFrame):
    charge = qc[qc.Direction == "Charge"].copy()
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), constrained_layout=True)
    for cell, part in charge.groupby("Cell"):
        axes[0].plot(part.C_rate, part.Rest_tail_slope_mV_per_min, "o-", label=cell)
        axes[1].plot(part.C_rate, 1000.0 * (part.Onset_V - part.Rest_end_V), "o-", label=cell)
    axes[0].axhline(0.1, color="black", lw=1, ls="--", label="0.1 mV/min reference")
    axes[0].set(xlabel="Charge rate", ylabel="Last-rest slope [mV/min]", title="10-min rest is still relaxing")
    axes[1].set(xlabel="Charge rate", ylabel="Onset minus rest [mV]", title="Measured current-step voltage rise")
    for ax in axes:
        ax.grid(alpha=0.25)
    axes[0].legend(fontsize=8)
    fig.savefig(OUT / "rest_non_equilibrium_diagnostics.png", dpi=220)
    plt.close(fig)

    eligible = candidate_summary[candidate_summary.primary_eligible].sort_values("qOCV_RMSE_2_98_mV")
    fig, ax = plt.subplots(figsize=(10.5, 5.2), constrained_layout=True)
    x = np.arange(len(eligible))
    ax.bar(x, eligible.qOCV_RMSE_2_98_mV, color="#4472C4")
    ax.set_xticks(x, eligible.candidate, rotation=30, ha="right")
    ax.set_ylabel("C/50 qOCV RMSE [mV]")
    ax.set_title("OCP source combinations with complete cathode coverage")
    ax.grid(axis="y", alpha=0.25)
    fig.savefig(OUT / "ocp_candidate_qocv_comparison.png", dpi=220)
    plt.close(fig)

    selected = dynamic_by_rate[
        dynamic_by_rate.scenario.isin(["rest_mean", "rest_cell", "endpoint", "coulomb_cell"])
    ]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.0), constrained_layout=True)
    names = list(dict.fromkeys(selected.scenario))
    colors = ["#4472C4", "#70AD47", "#A5A5A5", "#ED7D31"]
    width = 0.18
    rates = np.array(sorted(selected.C_rate.unique()), dtype=float)
    for offset, name, color in zip(np.linspace(-1.5, 1.5, len(names)) * width, names, colors):
        part = selected[selected.scenario == name].set_index("C_rate").loc[rates]
        axes[0].bar(rates + offset, part.initial_error_mV, width, label=name, color=color)
        axes[1].bar(rates + offset, part.full_RMSE_mV, width, label=name, color=color)
    axes[0].set(xlabel="Charge rate", ylabel="Model - experiment [mV]", title="First-point error")
    axes[1].set(xlabel="Charge rate", ylabel="RMSE [mV]", title="Full physical-time voltage error")
    for ax in axes:
        ax.set_xticks(rates, [f"{r:g}C" for r in rates])
        ax.axhline(0, color="black", lw=0.8)
        ax.grid(axis="y", alpha=0.25)
    axes[0].legend(fontsize=8)
    fig.savefig(OUT / "initialization_scenario_comparison.png", dpi=220)
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    context = common.build_context()
    qc = context.qc.copy()
    coulomb_soc = charge_coulomb_start_soc(qc, context.qcell)

    candidate_detail, candidate_summary = candidate_initial_state_table(
        qc, context.qcell, coulomb_soc
    )
    dynamic_detail, _ = simulate_initialization_scenarios(context, qc, coulomb_soc)
    dynamic_by_rate, dynamic_overall = summaries(dynamic_detail)

    candidate_detail.to_csv(OUT / "ocp_candidate_initial_soc_detail.csv", index=False, encoding="utf-8-sig")
    candidate_summary.to_csv(OUT / "ocp_candidate_summary.csv", index=False, encoding="utf-8-sig")
    dynamic_detail.to_csv(OUT / "initialization_dynamic_detail.csv", index=False, encoding="utf-8-sig")
    dynamic_by_rate.to_csv(OUT / "initialization_dynamic_by_rate.csv", index=False, encoding="utf-8-sig")
    dynamic_overall.to_csv(OUT / "initialization_dynamic_overall.csv", index=False, encoding="utf-8-sig")
    qc.to_csv(OUT / "july_protocol_charge_rest_audit.csv", index=False, encoding="utf-8-sig")
    make_plots(qc, candidate_summary, dynamic_by_rate)

    selected_detail = candidate_detail[candidate_detail.candidate == "nominal/old"]
    manifest = {
        "selected_model_parameters": FINAL,
        "qcell_Ah": context.qcell,
        "coulomb_start_soc_definition": (
            "0.5C charge starts at zero; later charge starts accumulate preceding "
            "CC+CV charge minus CC discharge, divided by C/50 Qcell"
        ),
        "selected_ocp_rest_inversion_by_rate": selected_detail.groupby("C_rate")[
            "rest_inverted_SOC_pct"
        ].mean().to_dict(),
        "coulomb_bookkeeping_soc_by_rate": selected_detail.groupby("C_rate")[
            "coulomb_bookkeeping_SOC_pct"
        ].mean().to_dict(),
        "dynamic_overall": dynamic_overall.to_dict(orient="records"),
    }
    (OUT / "analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print("\nCANDIDATE SUMMARY")
    print(candidate_summary.sort_values("qOCV_RMSE_2_98_mV").to_string(index=False))
    print("\nDYNAMIC BY RATE")
    print(dynamic_by_rate.to_string(index=False))
    print("\nDYNAMIC OVERALL")
    print(dynamic_overall.to_string(index=False))
    print("\nOUTPUT", OUT)


if __name__ == "__main__":
    main()
