"""Compare Ai2020 and measured low-rate anode OCP through S3-S7a.

The cathode OCP, geometry, particle-radius changes, measured kinetic/solid-
diffusion anchors, electrolyte transport, fitting objective, parameter bounds,
PSO settings, and experimental branches are held identical.  Each anode OCP
uses its own C/50-fitted stoichiometry window and capacity-consistent inventory.
"""

from __future__ import annotations

from itertools import combinations
from pathlib import Path
import sys
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import cumulative_six_stage_model_comparison as cumulative
import electrolyte_transport_set_comparison as transport
import gitt_ocp_analysis as ga
import ocp_dynamic_combination_comparison as ocp_compare
import s7_pso_voltage_fitting as s7


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "261019_lowrate_anode_full_pipeline_comparison"
STAGES = ("S3", "S4", "S5", "S6")
ANODE_LABELS = {
    "Ai2020": "Ai2020 anode",
    "Low-rate": "Measured low-rate anode",
}


def markdown_table(frame: pd.DataFrame) -> str:
    """Render a compact Markdown table without the optional tabulate package."""

    printable = frame.copy()
    columns = [str(column) for column in printable.columns]
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in printable.itertuples(index=False, name=None):
        values = [str(value).replace("|", "\\|") for value in row]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def write_readme(
    pipeline_summary: pd.DataFrame,
    fit_summary: pd.DataFrame,
    fit_parameters: pd.DataFrame,
    selected_parameters: tuple[str, ...],
) -> None:
    ai_s6 = fit_summary[fit_summary.Scenario == "Ai2020 S6"].iloc[0]
    low_s6 = fit_summary[fit_summary.Scenario == "Low-rate S6"].iloc[0]
    ai_s7 = fit_summary[fit_summary.Scenario == "Ai2020 S7a"].iloc[0]
    low_s7 = fit_summary[fit_summary.Scenario == "Low-rate S7a"].iloc[0]
    parameter_table = fit_parameters[
        [
            "Scenario",
            "Fitted_parameters",
            "Dsn_m2_s",
            "kn_prefactor",
            "brugg_n",
            "Objective_RMSE_mV",
        ]
    ].copy()
    parameter_table["Dsn_m2_s"] = parameter_table.Dsn_m2_s.map(
        lambda value: f"{float(value):.6e}"
    )
    parameter_table["kn_prefactor"] = parameter_table.kn_prefactor.map(
        lambda value: f"{float(value):.6e}"
    )
    parameter_table["brugg_n"] = parameter_table.brugg_n.map(
        lambda value: f"{float(value):.6f}"
    )
    parameter_table["Objective_RMSE_mV"] = parameter_table.Objective_RMSE_mV.map(
        lambda value: f"{float(value):.3f}"
    )
    readme = [
        "# Ai2020 vs measured low-rate anode OCP through S3-S7a",
        "",
        "양극 OCP, geometry, Rn/Rp, kp/Dsn/Dsp, ORegan2022 electrolyte, b_n=2.5, 실험 데이터와 오차 정의를 동일하게 유지했다. 각 음극 OCP에는 C/50에서 독립 선정한 window와 Qn/Qp, capacity-consistent cs,max를 사용했다.",
        "",
        "## S3-S6 동일 변경 적용 결과",
        "",
        markdown_table(pipeline_summary.round(4)),
        "",
        "## S7a 동일 민감도/식별도 및 PSO fitting",
        "",
        f"- 저율 OCP 기준으로 선별된 fitting subset: `{'+'.join(selected_parameters)}`",
        f"- Ai2020 S6→S7a 전체 RMSE: {ai_s6.Overall_full_time_RMSE_mV:.2f}→{ai_s7.Overall_full_time_RMSE_mV:.2f} mV",
        f"- 저율 S6→S7a 전체 RMSE: {low_s6.Overall_full_time_RMSE_mV:.2f}→{low_s7.Overall_full_time_RMSE_mV:.2f} mV",
        "",
        markdown_table(fit_summary.round(4)),
        "",
        markdown_table(parameter_table),
        "",
        "## 해석",
        "",
        "- 저율 음극 OCP에도 이후 변경과 fitting을 동일하게 적용할 수 있다.",
        "- 다만 terminal-voltage fitting이 OCP 차이를 kinetic/transport parameter로 보상할 수 있으므로, fitted RMSE만으로 음극 OCP의 물리적 우수성을 판단하지 않는다.",
        "- 저율 OCP의 absolute x 환산은 Ai2020 c_s,max, active fraction, thickness에 의존하므로 독립적인 실험 anchor가 아니다.",
        "- 일부 고율 분기에서 저율 OCP 보간 상한의 미세 외삽 경고가 발생했다.",
        "- 음극 전위 절대값은 3전극 검증 전까지 후보 간 경향으로만 해석한다.",
    ]
    (OUT / "README_KO.md").write_text("\n".join(readme) + "\n", encoding="utf-8")


def aggregate_with_anode_potential(
    detail: pd.DataFrame, anode_detail: pd.DataFrame
) -> pd.DataFrame:
    summary = s7.aggregate_metrics(detail)
    charge_minimum = (
        anode_detail[anode_detail.Direction == "Charge"]
        .groupby("Scenario", as_index=False)
        .agg(Minimum_charge_anode_potential_mV=("Minimum_anode_potential_mV", "min"))
    )
    return summary.merge(charge_minimum, on="Scenario", how="left")


def evaluate_stage(
    scenario: str,
    model: dict,
    stage: dict,
    anode: dict,
    cathode: dict,
    experiment: dict,
    qcell: float,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    runs = s7.run_model(model, stage, anode, cathode)
    result = {**s7.BASE, "runs": runs}
    branch_rows, anode_rows = s7.scenario_metrics(
        scenario, result, experiment, qcell
    )
    return pd.DataFrame(branch_rows), pd.DataFrame(anode_rows), result


def sensitivity_and_selection(evaluate) -> tuple[pd.DataFrame, pd.DataFrame, tuple[str, ...]]:
    base_coordinate = np.array(
        [s7.encode_value(name, s7.BASE[name]) for name in s7.PARAMETERS]
    )
    derivative_columns = []
    sensitivity_rows = []

    for index, name in enumerate(s7.PARAMETERS):
        lower = base_coordinate.copy()
        upper = base_coordinate.copy()
        lower[index] = max(0.0, lower[index] - s7.SENSITIVITY_STEP)
        upper[index] = min(1.0, upper[index] + s7.SENSITIVITY_STEP)
        lower_values = s7.vector_to_values(lower, s7.PARAMETERS)
        upper_values = s7.vector_to_values(upper, s7.PARAMETERS)
        lower_residual = evaluate(lower_values)["residual"]
        upper_residual = evaluate(upper_values)["residual"]
        derivative = (upper_residual - lower_residual) / (upper[index] - lower[index])
        derivative_columns.append(derivative)
        sensitivity_rows.append(
            {
                "Parameter": name,
                "Baseline": s7.BASE[name],
                "RMS_mV_per_normalized_coordinate": float(
                    np.sqrt(np.mean(derivative**2))
                ),
                "Max_abs_mV_per_normalized_coordinate": float(
                    np.max(np.abs(derivative))
                ),
            }
        )

    jacobian = np.column_stack(derivative_columns)
    sensitivity = pd.DataFrame(sensitivity_rows)
    sensitivity["Relative_sensitivity"] = (
        sensitivity.RMS_mV_per_normalized_coordinate
        / sensitivity.RMS_mV_per_normalized_coordinate.max()
    )

    subset_rows = []
    for size in range(1, len(s7.PARAMETERS) + 1):
        for subset in combinations(s7.PARAMETERS, size):
            indexes = [s7.PARAMETERS.index(name) for name in subset]
            subset_rows.append(s7.identifiability_report(jacobian[:, indexes], subset))
    subsets = pd.DataFrame(subset_rows)
    eligible_parameters = set(
        sensitivity.loc[
            sensitivity.Relative_sensitivity >= s7.MIN_RELATIVE_SENSITIVITY,
            "Parameter",
        ].astype(str)
    )
    eligible_subset = subsets.Parameters.apply(
        lambda value: set(str(value).split("+")).issubset(eligible_parameters)
    )
    passing = subsets[subsets.Pass & eligible_subset].sort_values(
        ["N_parameters", "Normalized_condition_number"], ascending=[False, True]
    )
    if passing.empty:
        selected = (
            str(sensitivity.sort_values("Relative_sensitivity").iloc[-1].Parameter),
        )
    else:
        selected = tuple(str(passing.iloc[0].Parameters).split("+"))
    return sensitivity, subsets, selected


def make_plot(pipeline_summary: pd.DataFrame, fit_summary: pd.DataFrame) -> None:
    figure, axes = plt.subplots(1, 3, figsize=(15.5, 4.8), constrained_layout=True)
    colors = {"Ai2020": "#D55E00", "Low-rate": "#0072B2"}

    for anode_name in ("Ai2020", "Low-rate"):
        part = pipeline_summary[pipeline_summary.Anode_OCP == anode_name]
        axes[0].plot(
            part.Stage,
            part.Overall_full_time_RMSE_mV,
            "o-",
            lw=2.1,
            color=colors[anode_name],
            label=ANODE_LABELS[anode_name],
        )
        axes[1].plot(
            part.Stage,
            part.Capacity_RMSE_pct,
            "o-",
            lw=2.1,
            color=colors[anode_name],
            label=ANODE_LABELS[anode_name],
        )

    positions = np.arange(len(fit_summary))
    axes[2].bar(
        positions,
        fit_summary.Overall_full_time_RMSE_mV,
        color=["#D55E00", "#0072B2", "#E69F00", "#56B4E9"],
    )
    axes[2].set_xticks(positions, fit_summary.Scenario, rotation=20, ha="right")

    axes[0].set(title="S3-S6 voltage", ylabel="Overall time RMSE [mV]")
    axes[1].set(title="S3-S6 capacity", ylabel="Capacity RMSE [%]")
    axes[2].set(title="S6 vs identifiability-screened S7a", ylabel="Overall time RMSE [mV]")
    for axis in axes:
        axis.grid(axis="y", alpha=0.25)
    axes[0].legend(frameon=False)
    figure.suptitle(
        "Ai2020 vs measured low-rate anode OCP under the same downstream pipeline",
        weight="bold",
    )
    figure.savefig(OUT / "full_pipeline_comparison.png", dpi=240, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    if "--report-only" in sys.argv:
        pipeline_summary = pd.read_csv(OUT / "pipeline_summary.csv")
        fit_summary = pd.read_csv(OUT / "s6_s7a_summary.csv")
        fit_parameters = pd.read_csv(OUT / "s6_s7a_parameters.csv")
        selected = tuple(
            str(
                fit_parameters.loc[
                    fit_parameters.Scenario == "Low-rate S7a", "Fitted_parameters"
                ].iloc[0]
            ).split("+")
        )
        write_readme(pipeline_summary, fit_summary, fit_parameters, selected)
        print(pipeline_summary.to_string(index=False), flush=True)
        print(fit_summary.to_string(index=False), flush=True)
        print(fit_parameters.to_string(index=False), flush=True)
        return
    _, _, _, _, experiment, qcell = transport.current_p0()
    stage_entries = cumulative.build_stages(qcell)

    source = pd.read_csv(
        ocp_compare.SOURCE / "ocp_source_combination_summary.csv"
    )
    selected_row = source[
        (source.anode_source == "nominal") & (source.cathode_source == "old")
    ].iloc[0]
    lowrate_row = source[
        (source.anode_source == "old") & (source.cathode_source == "old")
    ].iloc[0]
    anodes, cathodes = ocp_compare.load_ocp_sources(stage_entries["S6"]["model"])

    pipeline_detail_parts = []
    pipeline_anode_parts = []
    prepared = {}
    for stage_name in STAGES:
        entry = stage_entries[stage_name]
        candidates = {
            "Ai2020": (
                entry["model"],
                entry["stage"],
                entry["anode"],
                entry["cathode"],
            )
        }
        lowrate_model, lowrate_stage = ocp_compare.candidate_inputs(
            entry["model"],
            selected_row,
            lowrate_row,
            anodes["old"],
            cathodes["old"],
        )
        candidates["Low-rate"] = (
            lowrate_model,
            lowrate_stage,
            anodes["old"],
            cathodes["old"],
        )

        for anode_name, inputs in candidates.items():
            scenario = f"{stage_name} {anode_name}"
            print(scenario, flush=True)
            branch, anode_branch, result = evaluate_stage(
                scenario, *inputs, experiment, qcell
            )
            pipeline_detail_parts.append(branch)
            pipeline_anode_parts.append(anode_branch)
            prepared[(stage_name, anode_name)] = (*inputs, result)

    pipeline_detail = pd.concat(pipeline_detail_parts, ignore_index=True)
    pipeline_anode = pd.concat(pipeline_anode_parts, ignore_index=True)
    pipeline_summary = aggregate_with_anode_potential(
        pipeline_detail, pipeline_anode
    )
    pipeline_summary[["Stage", "Anode_OCP"]] = pipeline_summary.Scenario.str.split(
        " ", n=1, expand=True
    )
    pipeline_summary["qOCV_RMSE_2_98_mV"] = pipeline_summary.Anode_OCP.map(
        {
            "Ai2020": float(selected_row.qOCV_RMSE_2_98_mV),
            "Low-rate": float(lowrate_row.qOCV_RMSE_2_98_mV),
        }
    )

    low_model, low_stage, low_anode, low_cathode, low_s6_result = prepared[
        ("S6", "Low-rate")
    ]
    grids, targets = s7.make_experimental_grids(experiment)
    cache: dict[tuple[float, float, float], dict] = {}
    evaluations = []

    def evaluate(values: dict[str, float]) -> dict:
        key = tuple(round(float(values[name]), 18) for name in s7.PARAMETERS)
        if key in cache:
            return cache[key]
        started = time.perf_counter()
        try:
            model = s7.apply_values(low_model, values)
            runs = s7.run_model(model, low_stage, low_anode, low_cathode)
            residual, feasibility = s7.residual_vector(runs, grids, targets)
            status = "ok"
        except Exception as error:
            runs = None
            residual = np.full(6 * s7.N_GRID, 1000.0)
            feasibility = {
                "Minimum_time_margin_min": np.nan,
                "Early_cutoff_branches": 6,
            }
            status = f"failed:{type(error).__name__}"
        result = {
            "Objective_RMSE_mV": float(np.sqrt(np.mean(residual**2))),
            "Objective_MAE_mV": float(np.mean(np.abs(residual))),
            "residual": residual,
            "runs": runs,
            "Status": status,
            **feasibility,
        }
        cache[key] = result
        evaluations.append(
            {
                "Evaluation": len(evaluations) + 1,
                **values,
                "Objective_RMSE_mV": result["Objective_RMSE_mV"],
                "Objective_MAE_mV": result["Objective_MAE_mV"],
                "Status": status,
                "Elapsed_s": time.perf_counter() - started,
            }
        )
        return result

    base_low_result = evaluate(s7.BASE)
    sensitivity, subsets, selected_parameters = sensitivity_and_selection(evaluate)
    print(f"Low-rate selected parameters: {selected_parameters}", flush=True)
    low_fit_result, history = s7.pso_fit(
        "Low-rate S7a",
        selected_parameters,
        evaluate,
        informed=True,
        seed=261013,
    )

    low_s6_detail, low_s6_anode = s7.scenario_metrics(
        "Low-rate S6", base_low_result, experiment, qcell
    )
    low_s7_detail, low_s7_anode = s7.scenario_metrics(
        "Low-rate S7a", low_fit_result, experiment, qcell
    )

    ai_fit_detail = pd.read_csv(
        ROOT
        / "results"
        / "261013_s7_pso_voltage_fitting"
        / "s6_s7_branch_metrics.csv"
    )
    ai_fit_anode = pd.read_csv(
        ROOT
        / "results"
        / "261013_s7_pso_voltage_fitting"
        / "s6_s7_anode_potential.csv"
    )
    ai_fit_detail = ai_fit_detail[ai_fit_detail.Scenario.isin(["S6", "S7a"])].copy()
    ai_fit_anode = ai_fit_anode[ai_fit_anode.Scenario.isin(["S6", "S7a"])].copy()
    ai_fit_detail["Scenario"] = ai_fit_detail.Scenario.map(
        {"S6": "Ai2020 S6", "S7a": "Ai2020 S7a"}
    )
    ai_fit_anode["Scenario"] = ai_fit_anode.Scenario.map(
        {"S6": "Ai2020 S6", "S7a": "Ai2020 S7a"}
    )

    fit_detail = pd.concat(
        [ai_fit_detail, pd.DataFrame(low_s6_detail), pd.DataFrame(low_s7_detail)],
        ignore_index=True,
    )
    fit_anode = pd.concat(
        [ai_fit_anode, pd.DataFrame(low_s6_anode), pd.DataFrame(low_s7_anode)],
        ignore_index=True,
    )
    fit_summary = aggregate_with_anode_potential(fit_detail, fit_anode)

    ai_parameters = pd.read_csv(
        ROOT
        / "results"
        / "261013_s7_pso_voltage_fitting"
        / "s6_s7_parameters.csv"
    )
    ai_parameters = ai_parameters[ai_parameters.Scenario.isin(["S6", "S7a"])].copy()
    ai_parameters["Scenario"] = ai_parameters.Scenario.map(
        {"S6": "Ai2020 S6", "S7a": "Ai2020 S7a"}
    )
    low_parameters = pd.DataFrame(
        [
            {
                "Scenario": "Low-rate S6",
                "Fitted_parameters": "none",
                "Dsn_m2_s": s7.BASE["Dsn"],
                "kn_prefactor": s7.BASE["kn"],
                "brugg_n": s7.BASE["brugg_n"],
                "Objective_RMSE_mV": base_low_result["Objective_RMSE_mV"],
                "Objective_MAE_mV": base_low_result["Objective_MAE_mV"],
                "Early_cutoff_branches": base_low_result["Early_cutoff_branches"],
            },
            {
                "Scenario": "Low-rate S7a",
                "Fitted_parameters": "+".join(selected_parameters),
                "Dsn_m2_s": low_fit_result["Dsn"],
                "kn_prefactor": low_fit_result["kn"],
                "brugg_n": low_fit_result["brugg_n"],
                "Objective_RMSE_mV": low_fit_result["Objective_RMSE_mV"],
                "Objective_MAE_mV": low_fit_result["Objective_MAE_mV"],
                "Early_cutoff_branches": low_fit_result["Early_cutoff_branches"],
            },
        ]
    )
    fit_parameters = pd.concat(
        [ai_parameters, low_parameters], ignore_index=True, sort=False
    )

    pipeline_detail.to_csv(OUT / "pipeline_branch_metrics.csv", index=False, encoding="utf-8-sig")
    pipeline_anode.to_csv(OUT / "pipeline_anode_potential.csv", index=False, encoding="utf-8-sig")
    pipeline_summary.to_csv(OUT / "pipeline_summary.csv", index=False, encoding="utf-8-sig")
    sensitivity.to_csv(OUT / "lowrate_sensitivity.csv", index=False, encoding="utf-8-sig")
    subsets.to_csv(OUT / "lowrate_identifiability.csv", index=False, encoding="utf-8-sig")
    history.to_csv(OUT / "lowrate_pso_history.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(evaluations).to_csv(OUT / "lowrate_all_evaluations.csv", index=False, encoding="utf-8-sig")
    fit_detail.to_csv(OUT / "s6_s7a_branch_metrics.csv", index=False, encoding="utf-8-sig")
    fit_anode.to_csv(OUT / "s6_s7a_anode_potential.csv", index=False, encoding="utf-8-sig")
    fit_summary.to_csv(OUT / "s6_s7a_summary.csv", index=False, encoding="utf-8-sig")
    fit_parameters.to_csv(OUT / "s6_s7a_parameters.csv", index=False, encoding="utf-8-sig")
    make_plot(pipeline_summary, fit_summary)

    write_readme(pipeline_summary, fit_summary, fit_parameters, selected_parameters)
    print(pipeline_summary.to_string(index=False), flush=True)
    print(fit_summary.to_string(index=False), flush=True)
    print(fit_parameters.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
