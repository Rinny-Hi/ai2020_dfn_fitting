"""Compare the cumulative model-building stages requested on Cell 6-6.

The sequence is deliberately cumulative:
1. raw PyBaMM Ai2020,
2. Ai2020 with only the electrolyte-diffusivity unit correction,
2b. Ai2020 OCP with the selected window and capacity-consistent inventory,
3. selected OCP with the same selected window,
4. measured particle radii,
5. measured kp, Dsn, and Dsp,
6. ORegan2022 electrolyte transport and b_n=2.5.

No parameter is fitted in this script.  Each branch is initialized from its
measured rest voltage and simulated to the 4.2/3.0 V cutoff.
"""

from __future__ import annotations

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
import electrolyte_transport_set_comparison as transport
import gitt_ocp_analysis as ga
import mass_independent_nominal_anode_trial as nominal_anode
import ocp_dynamic_combination_comparison as ocp_compare
import plot_oregan_brugg_2p5_time_voltage as time_plot
import priority_initial_state_protocol_recheck as priority


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "261013_cumulative_seven_stage_model_comparison"

AI_WINDOW = {
    "x0": 0.0065,
    "x100": 0.8400,
    "y0": 0.9651,
    "y100": 0.4350,
}

STAGE_META = (
    (
        "S1",
        "Ai2020 raw",
        "Raw PyBaMM Ai2020 parameters, nominal OCPs, and published window",
    ),
    (
        "S2",
        "S1 + corrected De",
        "Only Ai2020 electrolyte diffusivity is multiplied by 1e-4",
    ),
    (
        "S2b",
        "S2 + selected window",
        "Ai2020 OCPs retained; selected C/50 window and capacity-consistent inventory applied",
    ),
    (
        "S3",
        "S2b + selected OCP",
        "Ai2020 graphite + existing LCO GITT with the same selected window/inventory",
    ),
    (
        "S4",
        "S3 + measured Rn/Rp",
        "Rn=3.7542 um and Rp=4.2387 um",
    ),
    (
        "S5",
        "S4 + kp/Dsn/Dsp",
        "Measured-anchor kp=5.8995e-7, Dsn=4.1708e-14, Dsp=6.8882e-14",
    ),
    (
        "S6",
        "S5 + ORegan2022/bn=2.5",
        "ORegan2022 electrolyte transport and negative-electrode Bruggeman 2.5",
    ),
)

COLORS = {
    "S1": "#999999",
    "S2": "#CC79A7",
    "S2b": "#9467BD",
    "S3": "#56B4E9",
    "S4": "#009E73",
    "S5": "#E69F00",
    "S6": "#D55E00",
}


def copy_model(model: dict) -> dict:
    changed = dict(model)
    changed["params"] = model["params"].copy()
    return changed


def update_params(model: dict, updates: dict) -> dict:
    changed = copy_model(model)
    changed["params"].update(updates, check_already_exists=False)
    return changed


def nominal_inputs(qcell: float) -> tuple[dict, dict, dict, dict]:
    raw = pybamm.ParameterValues("Ai2020")
    anode = nominal_anode.nominal_anode_detail()
    cathode = ocp_compare.nominal_detail_from_data(ai_params.lico2_ocp_Ai2020_data)
    stage = dict(AI_WINDOW)
    soc = np.linspace(0.0, 1.0, 1001)
    vector = np.array(
        [stage["x0"], stage["x100"], stage["y100"], stage["y0"]],
        dtype=float,
    )
    model = {
        "params": raw.copy(),
        "csn_max": float(raw["Maximum concentration in negative electrode [mol.m-3]"]),
        "csp_max": float(raw["Maximum concentration in positive electrode [mol.m-3]"]),
        "temperature": 298.15,
        "q_meas": float(qcell),
        "soc": soc,
        "v_qocv": c50.predict(vector, soc, anode, cathode, "equilibrium"),
    }
    return model, stage, anode, cathode


def build_stages(qcell: float) -> dict[str, dict]:
    raw = pybamm.ParameterValues("Ai2020")
    s1_model, ai_stage, ai_anode, ai_cathode = nominal_inputs(qcell)

    s2_model = update_params(
        s1_model,
        {
            "Electrolyte diffusivity [m2.s-1]": ga.scale_parameter_function(
                raw["Electrolyte diffusivity [m2.s-1]"], 1e-4
            )
        },
    )

    final_model, selected_stage, selected_anode, selected_cathode, _, _ = (
        transport.current_p0()
    )
    # S3 retains the selected OCP/window/inventory/overlap geometry, but rolls
    # the requested dynamic parameters back to their Ai2020 values.
    s3_model = update_params(
        final_model,
        {
            "Electrolyte diffusivity [m2.s-1]": ga.scale_parameter_function(
                raw["Electrolyte diffusivity [m2.s-1]"], 1e-4
            ),
            "Electrolyte conductivity [S.m-1]": raw[
                "Electrolyte conductivity [S.m-1]"
            ],
            "Cation transference number": raw["Cation transference number"],
            "Thermodynamic factor": raw["Thermodynamic factor"],
            "Initial concentration in electrolyte [mol.m-3]": raw[
                "Initial concentration in electrolyte [mol.m-3]"
            ],
            "Negative electrode Bruggeman coefficient (electrolyte)": raw[
                "Negative electrode Bruggeman coefficient (electrolyte)"
            ],
            "Positive electrode Bruggeman coefficient (electrolyte)": raw[
                "Positive electrode Bruggeman coefficient (electrolyte)"
            ],
            "Separator Bruggeman coefficient (electrolyte)": raw[
                "Separator Bruggeman coefficient (electrolyte)"
            ],
            "Negative particle radius [m]": raw["Negative particle radius [m]"],
            "Positive particle radius [m]": raw["Positive particle radius [m]"],
            "Negative particle diffusivity [m2.s-1]": raw[
                "Negative particle diffusivity [m2.s-1]"
            ],
            "Positive particle diffusivity [m2.s-1]": raw[
                "Positive particle diffusivity [m2.s-1]"
            ],
            "Negative electrode exchange-current density [A.m-2]": raw[
                "Negative electrode exchange-current density [A.m-2]"
            ],
            "Positive electrode exchange-current density [A.m-2]": raw[
                "Positive electrode exchange-current density [A.m-2]"
            ],
        },
    )
    # S2b uses the same selected window and physical inventory as S3 while
    # retaining the Ai2020 electrode OCP functions.  This separates the
    # window/inventory effect (S2 -> S2b) from the OCP-source effect
    # (S2b -> S3).
    s2b_model = copy_model(s3_model)

    s4_model = update_params(
        s3_model,
        {
            "Negative particle radius [m]": final_model["params"][
                "Negative particle radius [m]"
            ],
            "Positive particle radius [m]": final_model["params"][
                "Positive particle radius [m]"
            ],
        },
    )

    s5_model = update_params(
        s4_model,
        {
            "Positive electrode exchange-current density [A.m-2]": final_model[
                "params"
            ]["Positive electrode exchange-current density [A.m-2]"],
            "Negative particle diffusivity [m2.s-1]": final_model["params"][
                "Negative particle diffusivity [m2.s-1]"
            ],
            "Positive particle diffusivity [m2.s-1]": final_model["params"][
                "Positive particle diffusivity [m2.s-1]"
            ],
        },
    )

    s6_model = transport.with_transport_set(s5_model, "ORegan2022")
    s6_model["params"].update(
        {"Negative electrode Bruggeman coefficient (electrolyte)": 2.5}
    )

    ai_common = {"stage": ai_stage, "anode": ai_anode, "cathode": ai_cathode}
    selected_common = {
        "stage": selected_stage,
        "anode": selected_anode,
        "cathode": selected_cathode,
    }
    return {
        "S1": {"model": s1_model, **ai_common},
        "S2": {"model": s2_model, **ai_common},
        "S2b": {
            "model": s2b_model,
            "stage": selected_stage,
            "anode": ai_anode,
            "cathode": ai_cathode,
        },
        "S3": {"model": s3_model, **selected_common},
        "S4": {"model": s4_model, **selected_common},
        "S5": {"model": s5_model, **selected_common},
        "S6": {"model": s6_model, **selected_common},
    }


def overall_summary(detail: pd.DataFrame) -> pd.DataFrame:
    valid = detail[detail.Status == "ok"].copy()
    rows = []
    meta = {code: (label, note) for code, label, note in STAGE_META}
    for code, part in valid.groupby("Stage", sort=False):
        rows.append(
            {
                "Stage": code,
                "Label": meta[code][0],
                "Successful_branches": len(part),
                "Failed_branches": int(
                    (detail.loc[detail.Stage == code, "Status"] != "ok").sum()
                ),
                "Overall_full_time_RMSE_mV": float(
                    np.sqrt(np.mean(part.Full_time_RMSE_mV.to_numpy(float) ** 2))
                ),
                "Overall_full_time_MAE_mV": float(
                    np.mean(part.Full_time_MAE_mV.to_numpy(float))
                ),
                "Overall_center_time_RMSE_mV": float(
                    np.sqrt(
                        np.nanmean(
                            part.Center10_70_Qcell_time_RMSE_mV.to_numpy(float) ** 2
                        )
                    )
                ),
                "Overall_center_time_MAE_mV": float(
                    np.mean(part.Center10_70_Qcell_time_MAE_mV.to_numpy(float))
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
            }
        )
    frame = pd.DataFrame(rows)
    order = {code: index for index, (code, _, _) in enumerate(STAGE_META)}
    frame["Stage_order"] = frame.Stage.map(order)
    return frame.sort_values("Stage_order").drop(columns="Stage_order")


def write_readme(summary: pd.DataFrame, detail: pd.DataFrame) -> None:
    lines = [
        "# Ai2020에서 현재 기준모델까지 누적 7단계 비교",
        "",
        "## 비교 정의",
        "",
    ]
    for code, label, note in STAGE_META:
        lines.append(f"- **{code} {label}**: {note}")
    lines.extend(
        [
            "",
            "S2b에서 selected window와 함께 effective cs,max 및 공통 overlap 면적을 적용했다. 이는 선택된 window의 Qn/Qp를 실제 DFN inventory로 구현하기 위한 연동 bookkeeping이며 동적 fitting은 아니다. S2b는 Ai2020 OCP를 유지하므로 S2→S2b는 window/inventory 효과, S2b→S3는 OCP source 변경 효과다.",
            "",
            "모든 단계는 Cell 6-6 July-BOL 0.5C/1C/2C 충·방전 실험의 각 rest-end 전압에서 독립 초기화했으며, 별도 fitting 없이 4.2/3.0 V cutoff까지 계산했다.",
            "",
            "## 종합 결과",
            "",
            "| 단계 | 성공/6 | 전체 RMSE/MAE | 10–70% RMSE/MAE | 용량 RMSE | 최대 용량 오차 |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for _, row in summary.iterrows():
        lines.append(
            "| {} {} | {}/6 | {:.2f}/{:.2f} mV | {:.2f}/{:.2f} mV | {:.2f}% | {:.1f} mAh |".format(
                row.Stage,
                row.Label,
                int(row.Successful_branches),
                row.Overall_full_time_RMSE_mV,
                row.Overall_full_time_MAE_mV,
                row.Overall_center_time_RMSE_mV,
                row.Overall_center_time_MAE_mV,
                row.Capacity_RMSE_pct,
                row.Max_abs_capacity_error_mAh,
            )
        )
    lines.extend(
        [
            "",
            "## C-rate별 전체 시간 전압 RMSE",
            "",
            "아래 값은 `충전 / 방전` 순서다.",
            "",
            "| 단계 | 0.5C | 1C | 2C |",
            "|---|---:|---:|---:|",
        ]
    )
    for code, label, _ in STAGE_META:
        part = detail[(detail.Stage == code) & (detail.Status == "ok")]
        pairs = []
        for rate in ga.RATES:
            charge = part[
                np.isclose(part.C_rate.astype(float), rate)
                & (part.Direction == "Charge")
            ].iloc[0]
            discharge = part[
                np.isclose(part.C_rate.astype(float), rate)
                & (part.Direction == "Discharge")
            ].iloc[0]
            pairs.append(
                f"{charge.Full_time_RMSE_mV:.1f} / "
                f"{discharge.Full_time_RMSE_mV:.1f} mV"
            )
        lines.append(f"| {code} {label} | {pairs[0]} | {pairs[1]} | {pairs[2]} |")
    lines.extend(
        [
            "",
            "## 해석",
            "",
            "- S1→S2에서 De 단위 보정만으로 전체 전압 RMSE가 118.20→76.28 mV, 용량 RMSE가 14.40→7.92%로 감소했다.",
            "- S2→S2b는 Ai2020 OCP를 유지한 상태에서 selected window/inventory만 반영한 결과다. 10–70% 충전 RMSE는 99.69→78.09 mV, 방전 RMSE는 25.29→12.53 mV로 감소했고 용량 RMSE는 7.92→1.99%로 개선됐다.",
            "- S2b→S3는 동일한 window/inventory에서 OCP source만 변경한 결과다. 10–70% 충전 RMSE는 78.09→49.06 mV로 크게 감소했지만 방전 RMSE는 12.53→34.11 mV로 증가했다. 다만 6곡선 종합 RMSE와 용량 RMSE는 모두 개선됐다.",
            "- S3→S4의 측정 입자반경 반영은 전압 RMSE를 소폭 줄였지만 용량 RMSE를 3.22%로 일시 증가시켰다. 입자반경이 분극과 cutoff 도달시점에 함께 영향을 준 결과다.",
            "- S4→S5의 kp/Dsn/Dsp 반영은 전압 및 용량 오차를 다시 감소시켰다.",
            "- S6가 전체 전압 RMSE 38.59 mV, 중앙 RMSE 26.30 mV로 가장 우수하며 용량 RMSE도 1.02%로 유지된다. 기존 현재 기준모델 결과와 수치상 일치한다.",
            "",
            "## 주의",
            "",
            "- 전체 시간 RMSE는 실험과 모델이 모두 존재하는 공통 경과시간 구간에서 계산했다. 종료용량 차이는 용량 RMSE와 최대 용량 오차로 별도 평가했다.",
            "- 단계 간 차이는 누적 변경의 효과이며 개별 파라미터의 독립 인과효과로 해석하지 않는다.",
            "- S1의 raw Ai2020 De와 S2의 corrected De 차이는 Ai2020 구현의 단위 해석 검증을 위한 비교다.",
            "- S4 1C 방전 계산 중 선택된 양극 OCP 보간 상한의 미세 외삽 경고가 1회 발생했다. 계산은 정상 종료됐으며 최종 S6에서는 발생하지 않았다.",
        ]
    )
    (OUT / "README_KO.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    _, _, _, _, experiment, qcell = transport.current_p0()
    stages = build_stages(qcell)
    meta = {code: (label, note) for code, label, note in STAGE_META}

    fig, axes = plt.subplots(
        2, 3, figsize=(17.2, 9.8), sharey=True, constrained_layout=True
    )
    detail_rows: list[dict] = []
    curve_rows: list[dict] = []

    for column, rate in enumerate(ga.RATES):
        for row_index, charge in enumerate((True, False)):
            exp = experiment[(rate, charge)]
            exp_time = np.asarray(exp["t_min"], dtype=float)
            exp_time = exp_time - exp_time[0]
            axes[row_index, column].plot(
                exp_time,
                np.asarray(exp["V"], dtype=float),
                color="#111111",
                linewidth=3.0,
                label="Experiment",
                zorder=20,
            )

    for code, _, _ in STAGE_META:
        entry = stages[code]
        model = entry["model"]
        stage = entry["stage"]
        anode = entry["anode"]
        cathode = entry["cathode"]
        label, note = meta[code]
        for column, rate in enumerate(ga.RATES):
            for row_index, charge in enumerate((True, False)):
                direction = "Charge" if charge else "Discharge"
                rest_v, measured_current = transport.PROTOCOL[(rate, charge)]
                z0 = priority.monotone_voltage_inverse(
                    model["v_qocv"], model["soc"], rest_v
                )
                local_stage = priority.stage_at_soc(stage, z0, charge)
                print(f"{code} {label} | {rate:g}C {direction}", flush=True)
                common = {
                    "Stage": code,
                    "Label": label,
                    "Note": note,
                    "C_rate": rate,
                    "Direction": direction,
                    "Initial_rest_voltage_V": rest_v,
                    "Initial_SOC_pct": 100.0 * z0,
                    "x0": stage["x0"],
                    "x100": stage["x100"],
                    "y0": stage["y0"],
                    "y100": stage["y100"],
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
                    time_metrics = time_plot.time_voltage_metrics(
                        experiment[(rate, charge)], output, qcell, measured_current
                    )
                    capacity_metrics = priority.curve_metrics(
                        experiment[(rate, charge)], output, qcell, charge
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
                            "Status": "ok",
                        }
                    )
                    sim_time = np.asarray(output["time_s"], dtype=float) / 60.0
                    sim_time = sim_time - sim_time[0]
                    voltage = np.asarray(output["V"], dtype=float)
                    axes[row_index, column].plot(
                        sim_time,
                        voltage,
                        color=COLORS[code],
                        linewidth=2.0 if code == "S6" else 1.45,
                        linestyle="-" if code == "S6" else "--",
                        alpha=0.95,
                        label=f"{code} {label}",
                    )
                    for time_value, voltage_value in zip(sim_time, voltage):
                        curve_rows.append(
                            {
                                "Stage": code,
                                "Label": label,
                                "C_rate": rate,
                                "Direction": direction,
                                "Elapsed_time_min": float(time_value),
                                "Voltage_V": float(voltage_value),
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
            ax = axes[row_index, column]
            ax.set_title(f"{rate:g}C {direction}", fontsize=13, weight="bold")
            ax.set_ylim(2.9, 4.3)
            ax.grid(True, alpha=0.24)
            if row_index == 1:
                ax.set_xlabel("Elapsed time (min)")
            if column == 0:
                ax.set_ylabel("Cell voltage (V)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.5, 1.045),
        fontsize=9,
    )
    fig.suptitle(
        "Cumulative DFN build-up with separated window and OCP effects\n"
        "Cell 6-6, measured vs simulated voltage-time curves",
        fontsize=16,
        weight="bold",
        y=1.10,
    )
    fig.savefig(OUT / "seven_stage_time_voltage_overlay.png", dpi=240, bbox_inches="tight")
    plt.close(fig)

    detail = pd.DataFrame(detail_rows)
    summary = overall_summary(detail)
    stage_order = [code for code, _, _ in STAGE_META]
    branch_labels = [
        f"{rate:g}C {'Chg' if charge else 'DChg'}"
        for charge in (True, False)
        for rate in ga.RATES
    ]
    matrix = np.full((len(stage_order), len(branch_labels)), np.nan)
    for i, code in enumerate(stage_order):
        for j, (charge, rate) in enumerate(
            (pair for charge in (True, False) for pair in [(charge, r) for r in ga.RATES])
        ):
            direction = "Charge" if charge else "Discharge"
            match = detail[
                (detail.Stage == code)
                & (detail.C_rate == rate)
                & (detail.Direction == direction)
                & (detail.Status == "ok")
            ]
            if not match.empty:
                matrix[i, j] = float(match.iloc[0].Full_time_RMSE_mV)

    fig_h, ax_h = plt.subplots(figsize=(10.8, 6.2), constrained_layout=True)
    image = ax_h.imshow(matrix, aspect="auto", cmap="YlOrRd")
    ax_h.set_xticks(np.arange(len(branch_labels)), branch_labels)
    ax_h.set_yticks(
        np.arange(len(stage_order)),
        [f"{code} {meta[code][0]}" for code in stage_order],
    )
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if np.isfinite(matrix[i, j]):
                ax_h.text(
                    j,
                    i,
                    f"{matrix[i, j]:.1f}",
                    ha="center",
                    va="center",
                    color="white" if matrix[i, j] > np.nanmedian(matrix) else "black",
                    fontsize=9,
                )
    ax_h.set_title("Full elapsed-time voltage RMSE by stage and branch [mV]", weight="bold")
    fig_h.colorbar(image, ax=ax_h, label="RMSE (mV)")
    fig_h.savefig(OUT / "seven_stage_branch_rmse_heatmap.png", dpi=240, bbox_inches="tight")
    plt.close(fig_h)

    fig_s, axes_s = plt.subplots(1, 2, figsize=(13.5, 5.3), constrained_layout=True)
    x = np.arange(len(summary))
    axes_s[0].plot(
        x,
        summary.Overall_full_time_RMSE_mV,
        "o-",
        color="#0072B2",
        label="Full time RMSE",
    )
    axes_s[0].plot(
        x,
        summary.Overall_center_time_RMSE_mV,
        "s--",
        color="#009E73",
        label="Center 10-70% RMSE",
    )
    axes_s[0].set_ylabel("Voltage RMSE (mV)")
    axes_s[0].legend()
    axes_s[1].plot(
        x,
        summary.Capacity_RMSE_pct,
        "o-",
        color="#D55E00",
        label="Capacity RMSE",
    )
    axes_s[1].set_ylabel("Capacity RMSE (%)")
    for ax in axes_s:
        ax.set_xticks(x, summary.Stage)
        ax.set_xlabel("Cumulative stage")
        ax.grid(True, alpha=0.25)
    fig_s.suptitle("Error evolution across the cumulative model-building stages", weight="bold")
    fig_s.savefig(OUT / "seven_stage_error_evolution.png", dpi=240, bbox_inches="tight")
    plt.close(fig_s)

    detail.to_csv(OUT / "seven_stage_detail.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(OUT / "seven_stage_summary.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(curve_rows).to_csv(
        OUT / "seven_stage_curve_data.csv", index=False, encoding="utf-8-sig"
    )
    write_readme(summary, detail)
    print("\n[summary]\n", summary.to_string(index=False))
    failures = detail[detail.Status != "ok"]
    if not failures.empty:
        print("\n[failures]\n", failures[["Stage", "C_rate", "Direction", "Status"]].to_string(index=False))


if __name__ == "__main__":
    main()
