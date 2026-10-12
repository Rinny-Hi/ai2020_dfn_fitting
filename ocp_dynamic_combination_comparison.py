"""Cross-validate C/50-fitted OCP combinations on Cell 6-6 CC branches.

Each OCP pair uses its own C/50-fitted stoichiometry window and therefore its
own capacity-consistent Qn/Qp and cs,max.  Geometry, solid diffusion, kinetics,
electrolyte transport (ORegan2022), and b_n=2.5 are held fixed.  No dynamic
parameter is re-fitted for an OCP candidate.
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybamm
import pybamm.input.parameters.lithium_ion.Ai2020 as ai_params

import analyze_oregan_brugg_2p5_anode_potential as anode_analysis
import c50_ocp_source_combination_study as c50
import comprehensive_prefit_cross_cohort as cross
import electrode_hysteresis_quantification as ehq
import electrolyte_transport_set_comparison as transport
import gitt_ocp_analysis as ga
import plot_oregan_brugg_2p5_time_voltage as time_plot
import priority_initial_state_protocol_recheck as priority


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "results" / "260921_c50_ocp_source_combinations"
OUT = ROOT / "results" / "261012_ocp_dynamic_combination_comparison"

# These are the six combinations used in the OCP selection discussion.  The
# three omitted new/new cross-combinations can be added later, but do not add
# decision value while the new cathode reverse branch and new anode absolute
# stoichiometry are unresolved.
CANDIDATES = (
    ("nominal", "old"),
    ("old", "old"),
    ("old", "nominal"),
    ("nominal", "nominal"),
    ("new", "old"),
    ("nominal", "new"),
)

LABELS = {
    ("nominal", "old"): "Ai2020 anode + existing GITT cathode",
    ("old", "old"): "Existing low-rate anode + existing GITT cathode",
    ("old", "nominal"): "Existing low-rate anode + Ai2020 cathode",
    ("nominal", "nominal"): "Ai2020 anode + Ai2020 cathode",
    ("new", "old"): "New GITT anode + existing GITT cathode",
    ("nominal", "new"): "Ai2020 anode + new GITT cathode",
}

SHORT_LABELS = {
    ("nominal", "old"): "Ai anode + old cathode",
    ("old", "old"): "old anode + old cathode",
    ("old", "nominal"): "old anode + Ai cathode",
    ("nominal", "nominal"): "Ai anode + Ai cathode",
    ("new", "old"): "new anode + old cathode",
    ("nominal", "new"): "Ai anode + new cathode",
}

COLORS = {
    ("nominal", "old"): "#D55E00",
    ("old", "old"): "#0072B2",
    ("old", "nominal"): "#009E73",
    ("nominal", "nominal"): "#CC79A7",
    ("new", "old"): "#E69F00",
    ("nominal", "new"): "#56B4E9",
}


def new_anode_detail_with_input_fallback() -> dict:
    """Build the new-anode OCP without silently substituting another dataset."""

    paths = ga.get_paths(ROOT)
    configured = os.environ.get("AI2020_NEW_GITT_BOOK")
    candidates = [
        paths.gitt_book,
        Path(configured) if configured else None,
        Path(
            r"C:\Users\user\Documents\Codex\2026-09-18\https-github-com-rinny-hi-ai2020\work\ai2020_dfn_fitting\data\raw\260918_Anode_GITT\260918 Anode GITT v2-v3.xlsx"
        ),
    ]
    resolved = next((path for path in candidates if path and path.exists()), None)
    if resolved is None:
        raise FileNotFoundError(
            "260918 Anode GITT v2-v3.xlsx was not found. Set "
            "AI2020_NEW_GITT_BOOK to the original combined workbook."
        )
    paths.gitt_book = resolved
    gitt = ga.load_new_gitt(paths)
    built = ga.build_new_ocp(gitt)
    curves = built["curves"]

    def mean_branch(branch: str) -> tuple[np.ndarray, np.ndarray]:
        return ga._mean_overlap(
            [curves["v2"][branch], curves["v3"][branch]], n=3001
        )

    lith = mean_branch("lithiation")
    delith = mean_branch("delithiation")
    lo = max(lith[0].min(), delith[0].min())
    hi = min(lith[0].max(), delith[0].max())
    grid = np.linspace(lo, hi, 3001)
    lith_v = c50.interp_curve(lith, grid)
    delith_v = c50.interp_curve(delith, grid)
    return {
        "grid": grid,
        "equilibrium": 0.5 * (lith_v + delith_v),
        "lithiation": lith_v,
        "delithiation": delith_v,
        "directional_complete": True,
        "absolute_anchor": False,
        "source_book": str(resolved),
    }


def nominal_detail_from_data(data: tuple) -> dict:
    """Create a scalar equilibrium OCP from the packaged Ai2020 tabulation."""

    _, (stoichiometry, voltage) = data
    grid_raw = np.asarray(stoichiometry, dtype=float).reshape(-1)
    voltage_raw = np.asarray(voltage, dtype=float).reshape(-1)
    order = np.argsort(grid_raw)
    grid_raw, unique_index = np.unique(grid_raw[order], return_index=True)
    voltage_raw = voltage_raw[order][unique_index]
    grid = np.linspace(float(grid_raw.min()), float(grid_raw.max()), 3001)
    values = np.interp(grid, grid_raw, voltage_raw)
    return {
        "grid": grid,
        "equilibrium": values,
        "lithiation": values.copy(),
        "delithiation": values.copy(),
        "directional_complete": False,
        "absolute_anchor": True,
    }


def load_ocp_sources(base: dict) -> tuple[dict, dict]:
    """Build old/new/nominal OCP details exactly as in the C/50 study."""

    bundle = ga.load_legacy_bundle(ga.get_paths(ROOT))
    new_cathode, _ = c50.extract_new_cathode_delith()
    anodes = {
        "old": c50.old_detail(ehq.build_anode_detail(bundle, base)),
        "new": new_anode_detail_with_input_fallback(),
        "nominal": nominal_detail_from_data(
            ai_params.graphite_ocp_Enertech_Ai2020_data
        ),
    }
    cathodes = {
        "old": c50.old_detail(ehq.build_cathode_detail(bundle, base)),
        "new": new_cathode,
        "nominal": nominal_detail_from_data(ai_params.lico2_ocp_Ai2020_data),
    }
    return anodes, cathodes


def candidate_inputs(
    baseline: dict,
    selected_row: pd.Series,
    row: pd.Series,
    anode: dict,
    cathode: dict,
) -> tuple[dict, dict]:
    """Apply a candidate's C/50 window and capacity-consistent inventories."""

    model = dict(baseline)
    model["params"] = baseline["params"].copy()

    qn = float(row.Qn_Ah)
    qp = float(row.Qp_Ah)
    csn = float(baseline["csn_max"]) * qn / float(selected_row.Qn_Ah)
    csp = float(baseline["csp_max"]) * qp / float(selected_row.Qp_Ah)
    model["params"].update(
        {
            "Maximum concentration in negative electrode [mol.m-3]": csn,
            "Maximum concentration in positive electrode [mol.m-3]": csp,
        },
        check_already_exists=False,
    )
    model["csn_max"] = csn
    model["csp_max"] = csp
    model["qn_cell"] = qn
    model["qp_cell"] = qp

    stage = {
        "x0": float(row.x0),
        "x100": float(row.x100),
        "y0": float(row.y0),
        "y100": float(row.y100),
    }
    soc = np.linspace(0.0, 1.0, 1001)
    model["soc"] = soc
    model["v_qocv"] = c50.predict(
        np.array([stage["x0"], stage["x100"], stage["y100"], stage["y0"]]),
        soc,
        anode,
        cathode,
        "equilibrium",
    )
    return model, stage


def safe_time_metrics(
    experiment: dict,
    simulation: dict,
    qcell: float,
    measured_current: float,
) -> dict[str, float]:
    """Time-domain metrics with an explicit central-interval coverage flag."""

    result = time_plot.time_voltage_metrics(
        experiment, simulation, qcell, measured_current
    )
    center_start = 0.10 * qcell / measured_current * 60.0
    result["Center10_70_time_covered"] = bool(
        result["DFN_end_time_min"] >= center_start
    )
    if not result["Center10_70_time_covered"]:
        result["Center10_70_Qcell_time_RMSE_mV"] = np.nan
    return result


def combined_summary(detail: pd.DataFrame) -> pd.DataFrame:
    valid = detail[detail.Status == "ok"].copy()
    rows = []
    for scenario, part in valid.groupby("Scenario", sort=False):
        failures = int(
            (detail.loc[detail.Scenario == scenario, "Status"] != "ok").sum()
        )
        rows.append(
            {
                "Scenario": scenario,
                "Trust": part.Trust.iloc[0],
                "Successful_branches": len(part),
                "Failed_branches": failures,
                "Overall_full_time_RMSE_mV": float(
                    np.sqrt(np.mean(part.Full_time_RMSE_mV.to_numpy(float) ** 2))
                ),
                "Overall_center_time_RMSE_mV": float(
                    np.sqrt(
                        np.nanmean(
                            part.Center10_70_Qcell_time_RMSE_mV.to_numpy(float) ** 2
                        )
                    )
                ),
                "Capacity_RMSE_pct": float(
                    np.sqrt(np.mean(part.Capacity_error_pct.to_numpy(float) ** 2))
                ),
                "Mean_abs_capacity_error_mAh": float(
                    np.mean(np.abs(part.Capacity_error_mAh.to_numpy(float)))
                ),
                "Max_abs_capacity_error_mAh": float(
                    np.max(np.abs(part.Capacity_error_mAh.to_numpy(float)))
                ),
                "Minimum_charge_anode_potential_mV": float(
                    part.loc[
                        part.Direction == "Charge", "Minimum_anode_potential_mV"
                    ].min()
                ),
                "Maximum_charge_time_below_0V_pct": float(
                    part.loc[
                        part.Direction == "Charge", "Time_below_0V_pct"
                    ].max()
                ),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["Failed_branches", "Overall_full_time_RMSE_mV"]
    )


def branch_table(detail: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "Scenario",
        "C_rate",
        "Direction",
        "Full_time_RMSE_mV",
        "Center10_70_Qcell_time_RMSE_mV",
        "Capacity_error_pct",
        "Capacity_error_mAh",
        "Minimum_anode_potential_mV",
        "Time_below_0V_pct",
        "Status",
    ]
    return detail[columns].sort_values(["Scenario", "C_rate", "Direction"])


def write_readme(summary: pd.DataFrame, source_rows: pd.DataFrame) -> None:
    lines = [
        "# OCP 조합별 동적 교차검증",
        "",
        "## 고정 조건",
        "",
        "- Cell 6-6 July-BOL 0.5C/1C/2C 충·방전 6개 CC branch를 사용했다.",
        "- ORegan2022 전해질 수송식과 음극 Bruggeman 계수 `b_n=2.5`를 사용했다.",
        "- geometry, Dsn/Dsp, kn/kp는 현재 no-fit 기준값으로 고정했다.",
        "- 각 OCP 조합에는 해당 조합의 C/50 fitted stoichiometry window, Qn/Qp, cs,max를 적용했다.",
        "- 동적 파라미터는 OCP 조합별로 다시 fitting하지 않았다.",
        "- 전압 RMSE는 공통 경과시간 구간에서 계산하고 용량 오차는 cutoff 종료용량으로 별도 평가했다.",
        "",
        "## 종합 결과",
        "",
        "| OCP 조합 | 신뢰 구분 | 성공/6 | 전체 시간 RMSE | 중앙 시간 RMSE | 용량 RMSE | 최대 용량 오차 | 충전 최저 음극전위 |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for _, row in summary.iterrows():
        lines.append(
            "| {} | {} | {}/6 | {:.2f} mV | {:.2f} mV | {:.2f}% | {:.1f} mAh | {:+.1f} mV |".format(
                row.Scenario,
                row.Trust,
                int(row.Successful_branches),
                row.Overall_full_time_RMSE_mV,
                row.Overall_center_time_RMSE_mV,
                row.Capacity_RMSE_pct,
                row.Max_abs_capacity_error_mAh,
                row.Minimum_charge_anode_potential_mV,
            )
        )
    lines.extend(
        [
            "",
            "## 해석 제한",
            "",
            "- `new GITT anode`는 절대 stoichiometry anchor가 확립되지 않아 sensitivity-only이다.",
            "- `new GITT cathode`는 reverse lithiation branch가 미완결이므로 sensitivity-only이다.",
            "- 동적 전압 RMSE가 낮더라도 측정 완결성과 절대 stoichiometry가 부족한 조합을 최종 OCP로 선정하지 않는다.",
            "- OCP 변경에 따른 음극전위 차이는 모델 예측 민감도이며 3전극 측정에 의한 절대 검증값이 아니다.",
            "",
            "## C/50 입력 window",
            "",
            "| OCP 조합 | x0 | x100 | y0 | y100 | Qn | Qp | qOCV RMSE | primary eligible |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---|",
        ]
    )
    for _, row in source_rows.iterrows():
        key = (row.anode_source, row.cathode_source)
        lines.append(
            "| {} | {:.6f} | {:.6f} | {:.6f} | {:.6f} | {:.4f} Ah | {:.4f} Ah | {:.2f} mV | {} |".format(
                LABELS[key],
                row.x0,
                row.x100,
                row.y0,
                row.y100,
                row.Qn_Ah,
                row.Qp_Ah,
                row.qOCV_RMSE_2_98_mV,
                bool(row.primary_eligible),
            )
        )
    (OUT / "README_KO.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    baseline, _, _, _, experiment, qcell = transport.current_p0()
    baseline = transport.with_transport_set(baseline, "ORegan2022")
    baseline["params"].update(
        {"Negative electrode Bruggeman coefficient (electrolyte)": 2.5}
    )

    source = pd.read_csv(SOURCE / "ocp_source_combination_summary.csv")
    selected_row = source[
        (source.anode_source == "nominal") & (source.cathode_source == "old")
    ].iloc[0]
    anodes, cathodes = load_ocp_sources(baseline)

    fig_v, axes_v = plt.subplots(
        2, 3, figsize=(16.5, 9.2), sharey=True, constrained_layout=True
    )
    fig_a, axes_a = plt.subplots(
        2, 3, figsize=(16.5, 9.2), sharey="row", constrained_layout=True
    )
    detail_rows: list[dict] = []
    curve_rows: list[dict] = []

    # Draw experiment once per panel.
    for column, rate in enumerate(ga.RATES):
        for row_index, charge in enumerate((True, False)):
            exp = experiment[(rate, charge)]
            exp_time = np.asarray(exp["t_min"], dtype=float)
            exp_time = exp_time - exp_time[0]
            axes_v[row_index, column].plot(
                exp_time,
                np.asarray(exp["V"], dtype=float),
                color="#202020",
                linewidth=2.7,
                label="Experiment",
                zorder=10,
            )

    for anode_name, cathode_name in CANDIDATES:
        key = (anode_name, cathode_name)
        label = LABELS[key]
        short = SHORT_LABELS[key]
        row = source[
            (source.anode_source == anode_name)
            & (source.cathode_source == cathode_name)
        ].iloc[0]
        trust = "primary eligible" if bool(row.primary_eligible) else "sensitivity-only"
        anode = anodes[anode_name]
        cathode = cathodes[cathode_name]
        model, stage = candidate_inputs(
            baseline, selected_row, row, anode, cathode
        )

        for column, rate in enumerate(ga.RATES):
            for row_index, charge in enumerate((True, False)):
                direction = "Charge" if charge else "Discharge"
                rest_v, measured_current = transport.PROTOCOL[(rate, charge)]
                z0 = priority.monotone_voltage_inverse(
                    model["v_qocv"], model["soc"], rest_v
                )
                local_stage = priority.stage_at_soc(stage, z0, charge)
                print(f"{label} | {rate:g}C {direction}", flush=True)
                common = {
                    "Scenario": label,
                    "Short_label": short,
                    "Anode_source": anode_name,
                    "Cathode_source": cathode_name,
                    "Trust": trust,
                    "C_rate": rate,
                    "Direction": direction,
                    "Initial_rest_voltage_V": rest_v,
                    "Initial_SOC_pct": 100.0 * z0,
                    "x0": stage["x0"],
                    "x100": stage["x100"],
                    "y0": stage["y0"],
                    "y100": stage["y100"],
                    "Qn_Ah": float(row.Qn_Ah),
                    "Qp_Ah": float(row.Qp_Ah),
                    "csn_max_mol_m3": float(model["csn_max"]),
                    "csp_max_mol_m3": float(model["csp_max"]),
                }
                try:
                    output = anode_analysis.run_branch(
                        measured_current / cross.NOMINAL_CAPACITY_AH,
                        charge,
                        local_stage,
                        model,
                        anode,
                        cathode,
                    )
                    exp = experiment[(rate, charge)]
                    time_metrics = safe_time_metrics(
                        exp, output, qcell, measured_current
                    )
                    capacity_metrics = priority.curve_metrics(
                        exp, output, qcell, charge
                    )
                    potential = np.asarray(output["potential_V"], dtype=float)
                    minimum_index = int(np.nanargmin(potential))
                    below_zero = anode_analysis.time_fraction(
                        np.asarray(output["time_s"], dtype=float), potential < 0.0
                    )
                    detail_rows.append(
                        {
                            **common,
                            **time_metrics,
                            "Capacity_error_pct": capacity_metrics[
                                "Capacity_error_pct"
                            ],
                            "Capacity_error_mAh": 1000.0
                            * (
                                capacity_metrics["Q_end_model_Ah"]
                                - capacity_metrics["Q_end_exp_Ah"]
                            ),
                            "Minimum_anode_potential_mV": 1000.0
                            * float(potential[minimum_index]),
                            "Time_at_minimum_min": float(
                                output["time_s"][minimum_index] / 60.0
                            ),
                            "Surface_stoichiometry_at_minimum": float(
                                output["surface_sto"][minimum_index]
                            ),
                            "Time_below_0V_pct": 100.0 * below_zero,
                            "Status": "ok",
                        }
                    )

                    sim_time = np.asarray(output["time_s"], dtype=float) / 60.0
                    sim_time = sim_time - sim_time[0]
                    axes_v[row_index, column].plot(
                        sim_time,
                        np.asarray(output["V"], dtype=float),
                        color=COLORS[key],
                        linewidth=1.55,
                        alpha=0.92,
                        label=short,
                    )
                    axes_a[row_index, column].plot(
                        np.asarray(output["capacity_Ah"], dtype=float),
                        potential,
                        color=COLORS[key],
                        linewidth=1.7,
                        alpha=0.95,
                        label=short,
                    )
                    for time_min, voltage, anode_v, capacity in zip(
                        sim_time,
                        np.asarray(output["V"], dtype=float),
                        potential,
                        np.asarray(output["capacity_Ah"], dtype=float),
                    ):
                        curve_rows.append(
                            {
                                "Scenario": label,
                                "C_rate": rate,
                                "Direction": direction,
                                "Elapsed_time_min": float(time_min),
                                "Transferred_capacity_Ah": float(capacity),
                                "Cell_voltage_V": float(voltage),
                                "Anode_potential_V": float(anode_v),
                            }
                        )
                except Exception as exc:
                    detail_rows.append(
                        {
                            **common,
                            "Status": f"failed: {type(exc).__name__}: {exc}",
                        }
                    )
                    print(f"FAILED: {type(exc).__name__}: {exc}", flush=True)

    for column, rate in enumerate(ga.RATES):
        for row_index, charge in enumerate((True, False)):
            direction = "Charge" if charge else "Discharge"
            ax_v = axes_v[row_index, column]
            ax_v.set_title(f"{rate:g}C {direction}", weight="bold")
            ax_v.set_ylim(2.9, 4.3)
            ax_v.grid(True, alpha=0.25)
            if row_index == 1:
                ax_v.set_xlabel("Elapsed time (min)")
            if column == 0:
                ax_v.set_ylabel("Cell voltage (V)")

            ax_a = axes_a[row_index, column]
            ax_a.axhline(0.0, color="#202020", linestyle=":", linewidth=1.1)
            ax_a.set_title(f"{rate:g}C {direction}", weight="bold")
            ax_a.grid(True, alpha=0.25)
            if row_index == 1:
                ax_a.set_xlabel("Transferred capacity (Ah)")
            if column == 0:
                ax_a.set_ylabel("Anode interfacial potential (V)")

    handles_v, labels_v = axes_v[0, 0].get_legend_handles_labels()
    fig_v.legend(
        handles_v,
        labels_v,
        loc="upper center",
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.5, 1.04),
    )
    fig_v.suptitle(
        "Cell 6-6 voltage-time cross-validation of C/50-fitted OCP pairs\n"
        "ORegan2022 electrolyte, $b_n=2.5$; no dynamic re-fitting",
        fontsize=15,
        weight="bold",
        y=1.09,
    )
    fig_v.savefig(OUT / "ocp_pair_time_voltage_overlay.png", dpi=240, bbox_inches="tight")
    plt.close(fig_v)

    handles_a, labels_a = axes_a[0, 0].get_legend_handles_labels()
    fig_a.legend(
        handles_a,
        labels_a,
        loc="upper center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 1.04),
    )
    fig_a.suptitle(
        "Predicted anode potential sensitivity to the OCP pair\n"
        "ORegan2022 electrolyte, $b_n=2.5$; no three-electrode validation",
        fontsize=15,
        weight="bold",
        y=1.09,
    )
    fig_a.savefig(OUT / "ocp_pair_anode_potential_overlay.png", dpi=240, bbox_inches="tight")
    plt.close(fig_a)

    detail = pd.DataFrame(detail_rows)
    summary = combined_summary(detail)
    branches = branch_table(detail)
    selected_source_rows = pd.concat(
        [
            source[
                (source.anode_source == anode_name)
                & (source.cathode_source == cathode_name)
            ]
            for anode_name, cathode_name in CANDIDATES
        ],
        ignore_index=True,
    )
    detail.to_csv(OUT / "dynamic_detail.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(OUT / "dynamic_summary.csv", index=False, encoding="utf-8-sig")
    branches.to_csv(OUT / "branch_summary.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(curve_rows).to_csv(
        OUT / "curve_data.csv", index=False, encoding="utf-8-sig"
    )
    selected_source_rows.to_csv(
        OUT / "c50_window_inputs.csv", index=False, encoding="utf-8-sig"
    )
    write_readme(summary, selected_source_rows)

    print("\n[summary]\n", summary.to_string(index=False))
    failures = detail[detail.Status != "ok"]
    if not failures.empty:
        print(
            "\n[failures]\n",
            failures[["Scenario", "C_rate", "Direction", "Status"]].to_string(
                index=False
            ),
        )


if __name__ == "__main__":
    main()
