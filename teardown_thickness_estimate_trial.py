"""Apply provisional coating thicknesses inferred from double-sided electrodes.

The measured total thickness is assumed to include one current collector and
two equal coating layers.  Current-collector thicknesses are provisionally
taken from the current Ai2020-based parameter set.  Maximum solid
concentrations are recomputed from the same measured half-cell capacity and
coating mass used in ``chen2020_partial_trial.py`` so that geometry and active
material inventory remain physically consistent.
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

import chen2020_partial_trial as ct
import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import ocp_grounded_candidate_search as common
import old_notebook_method_comparison as omc


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260920_teardown_thickness_estimate"
RESULTS.mkdir(parents=True, exist_ok=True)

POSITIVE_TOTAL_DOUBLE_SIDED_M = 138e-6
NEGATIVE_TOTAL_DOUBLE_SIDED_M = 155e-6


def make_consistent_model(
    base: dict[str, Any],
    bundle: dict[str, Any],
    length_n: float,
    length_p: float,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any], np.ndarray]:
    """Update L and c_s,max together while preserving measured inventory."""
    params = base["params"].copy()
    disk_diameter = float(bundle["prior"]("disk_diameter_m"))
    disk_area = np.pi * (disk_diameter / 2.0) ** 2
    cell_area = (
        float(params["Electrode height [m]"])
        * float(params["Electrode width [m]"])
        * float(params["Number of electrodes connected in parallel to make a cell"])
    )
    capacities = ct.branch_capacity_summary(bundle["anode_lowrate"])
    qn_disk = float(capacities.capacity_Ah.mean())
    eps_n = float(params["Negative electrode active material volume fraction"])
    eps_p = float(params["Positive electrode active material volume fraction"])

    csn = qn_disk * 3600.0 / (
        ga.FARADAY_CONSTANT * eps_n * length_n * disk_area
    )
    positive_mass_areal = ct.POSITIVE_COATING_MASS_KG / disk_area
    csp = positive_mass_areal / (
        ct.LCO_MOLAR_MASS_KG_PER_MOL * length_p * eps_p
    )
    qn = ga.FARADAY_CONSTANT * csn * eps_n * length_n * cell_area / 3600.0
    qp = ga.FARADAY_CONSTANT * csp * eps_p * length_p * cell_area / 3600.0

    params.update(
        {
            "Negative electrode thickness [m]": length_n,
            "Positive electrode thickness [m]": length_p,
            "Maximum concentration in negative electrode [mol.m-3]": csn,
            "Maximum concentration in positive electrode [mol.m-3]": csp,
        }
    )
    model = dict(base)
    model["params"] = params
    model["csn_max"] = csn
    model["csp_max"] = csp
    model["qn_cell"] = qn
    model["qp_cell"] = qp
    model["delta_x"] = model["q_meas"] / qn
    model["delta_y"] = model["q_meas"] / qp

    anode = ehq.build_anode_detail(bundle, model)
    cathode = ehq.build_cathode_detail(bundle, model)
    model["y_exp"] = cathode["grid"]
    model["up_exp"] = cathode["equilibrium"]
    fit, qocv = ct.fit_fixed_capacities(model, anode, cathode, qn, qp)
    return model, anode, cathode, fit, qocv


def run_dynamic(
    scenario: str,
    model: dict[str, Any],
    anode: dict[str, Any],
    cathode: dict[str, Any],
    fit: dict[str, Any],
    experiment: dict[tuple[float, bool], dict[str, np.ndarray]],
) -> tuple[pd.DataFrame, dict[tuple[float, bool], dict[str, np.ndarray]]]:
    stage = {k: fit[k] for k in ("x0", "x100", "y100", "y0")}
    anode_scaled = common.scaled_detail(anode, 0.50)
    cathode_scaled = common.scaled_detail(cathode, 0.50)
    rows = []
    simulations = {}
    for rate in ga.RATES:
        for charge in (True, False):
            direction = "Charge" if charge else "Discharge"
            print(f"{scenario} | {rate:g}C | {direction}", flush=True)
            sim = ehq.run_dfn(
                rate,
                charge,
                stage,
                model,
                anode_scaled,
                cathode_scaled,
                True,
                True,
            )
            simulations[(rate, charge)] = sim
            metric = omc.old_time_metrics(experiment[(rate, charge)], sim)
            metric["capacity_error_Ah"] = (
                metric["end_time_error_min"] * rate * 2.28 / 60.0
            )
            rows.append(
                {
                    "Scenario": scenario,
                    "C_rate": rate,
                    "Direction": direction,
                    **metric,
                }
            )
    return pd.DataFrame(rows), simulations


def main() -> None:
    bundle = ga.load_legacy_bundle(ga.get_paths(ROOT))
    base = ga.build_legacy_model_inputs(bundle)
    params = base["params"]

    cc_p = float(params["Positive current collector thickness [m]"])
    cc_n = float(params["Negative current collector thickness [m]"])
    old_lp = float(params["Positive electrode thickness [m]"])
    old_ln = float(params["Negative electrode thickness [m]"])
    new_lp = (POSITIVE_TOTAL_DOUBLE_SIDED_M - cc_p) / 2.0
    new_ln = (NEGATIVE_TOTAL_DOUBLE_SIDED_M - cc_n) / 2.0
    if new_lp <= 0 or new_ln <= 0:
        raise ValueError("Inferred coating thickness must be positive")

    scenarios = {}
    for name, ln, lp in (
        ("Previous Chen-consistent thickness", old_ln, old_lp),
        ("Estimated teardown thickness", new_ln, new_lp),
    ):
        scenarios[name] = make_consistent_model(base, bundle, ln, lp)

    assumption_rows = [
        {
            "Electrode": "Positive",
            "Measured_double_sided_total_um": POSITIVE_TOTAL_DOUBLE_SIDED_M * 1e6,
            "Assumed_current_collector_um": cc_p * 1e6,
            "Previous_single_side_coating_um": old_lp * 1e6,
            "Estimated_single_side_coating_um": new_lp * 1e6,
            "Thickness_change_percent": (new_lp / old_lp - 1.0) * 100.0,
            "Formula": "(double-sided total - current collector) / 2",
        },
        {
            "Electrode": "Negative",
            "Measured_double_sided_total_um": NEGATIVE_TOTAL_DOUBLE_SIDED_M * 1e6,
            "Assumed_current_collector_um": cc_n * 1e6,
            "Previous_single_side_coating_um": old_ln * 1e6,
            "Estimated_single_side_coating_um": new_ln * 1e6,
            "Thickness_change_percent": (new_ln / old_ln - 1.0) * 100.0,
            "Formula": "(double-sided total - current collector) / 2",
        },
    ]
    assumptions = pd.DataFrame(assumption_rows)
    assumptions.to_csv(
        RESULTS / "teardown_thickness_assumptions.csv", index=False, encoding="utf-8-sig"
    )

    experiment = omc.load_old_dynamic_data(base["dynamic"], base["q_meas"])
    all_dynamic = []
    all_simulations = {}
    summary_rows = []
    for name, (model, anode, cathode, fit, qocv) in scenarios.items():
        dynamic, simulations = run_dynamic(
            name, model, anode, cathode, fit, experiment
        )
        all_dynamic.append(dynamic)
        all_simulations[name] = simulations
        summary_rows.append(
            {
                "Scenario": name,
                "Negative_thickness_um": float(model["params"]["Negative electrode thickness [m]"]) * 1e6,
                "Positive_thickness_um": float(model["params"]["Positive electrode thickness [m]"]) * 1e6,
                "Negative_c_s_max_mol_m3": model["csn_max"],
                "Positive_c_s_max_mol_m3": model["csp_max"],
                "Q_n_Ah": fit["Q_n_Ah"],
                "Q_p_Ah": fit["Q_p_Ah"],
                "x0": fit["x0"],
                "x100": fit["x100"],
                "y100": fit["y100"],
                "y0": fit["y0"],
                "qOCV_MAE_2_98_mV": fit["qOCV_MAE_2_98_mV"],
                "qOCV_RMSE_2_98_mV": fit["qOCV_RMSE_2_98_mV"],
                "max_endpoint_abs_mV": fit["max_endpoint_abs_mV"],
                "Dynamic_MAE_SOC10_70_mV": float(dynamic["MAE_SOC10_70_mV"].mean()),
                "Mean_abs_capacity_error_Ah": float(dynamic["capacity_error_Ah"].abs().mean()),
                "Max_abs_capacity_error_Ah": float(dynamic["capacity_error_Ah"].abs().max()),
            }
        )
    dynamic_all = pd.concat(all_dynamic, ignore_index=True)
    summary = pd.DataFrame(summary_rows)
    dynamic_all.to_csv(
        RESULTS / "teardown_thickness_dynamic_validation.csv",
        index=False,
        encoding="utf-8-sig",
    )
    summary.to_csv(
        RESULTS / "teardown_thickness_comparison_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )

    # Dynamic overlays: experiment, previous thickness, provisional teardown thickness.
    fig, axes = plt.subplots(2, 3, figsize=(15.8, 8.3))
    old_name = "Previous Chen-consistent thickness"
    new_name = "Estimated teardown thickness"
    for col, rate in enumerate(ga.RATES):
        for row_index, charge in enumerate((True, False)):
            ax = axes[row_index, col]
            direction = "Charge" if charge else "Discharge"
            exp = experiment[(rate, charge)]
            old_sim = all_simulations[old_name][(rate, charge)]
            new_sim = all_simulations[new_name][(rate, charge)]
            old_metric = dynamic_all[
                (dynamic_all.Scenario == old_name)
                & np.isclose(dynamic_all.C_rate, rate)
                & (dynamic_all.Direction == direction)
            ].iloc[0]
            new_metric = dynamic_all[
                (dynamic_all.Scenario == new_name)
                & np.isclose(dynamic_all.C_rate, rate)
                & (dynamic_all.Direction == direction)
            ].iloc[0]
            ax.plot(exp["t_min"], exp["V"], color="black", lw=2.0, label="Experiment")
            ax.plot(
                old_sim["t_min"], old_sim["V"], color="#7f7f7f", lw=1.6,
                ls="--", label="Previous thickness"
            )
            ax.plot(
                new_sim["t_min"], new_sim["V"], color="#1f77b4", lw=1.8,
                label="Estimated thickness"
            )
            ax.set_title(
                f"{rate:g}C {direction}\nMAE old {old_metric.MAE_SOC10_70_mV:.1f} → new {new_metric.MAE_SOC10_70_mV:.1f} mV"
            )
            ax.set_xlabel("Time (min)")
            ax.set_ylabel("Voltage (V)")
            ax.grid(alpha=0.25)
            if row_index == 0 and col == 0:
                ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(
        RESULTS / "teardown_thickness_dynamic_comparison.png",
        dpi=220,
        bbox_inches="tight",
    )
    plt.close(fig)

    # Show the invariant equilibrium fit and the dynamic/capacity trade-off.
    old_model, _, _, old_fit, old_qocv = scenarios[old_name]
    new_model, _, _, new_fit, new_qocv = scenarios[new_name]
    fig, axes = plt.subplots(1, 3, figsize=(16.2, 5.2))
    axes[0].plot(base["soc"], base["v_qocv"], color="black", lw=2.3, label="Measured qOCV")
    axes[0].plot(base["soc"], old_qocv, color="#7f7f7f", lw=2.0, ls="--", label="Previous thickness")
    axes[0].plot(base["soc"], new_qocv, color="#1f77b4", lw=1.6, label="Estimated thickness")
    axes[0].set(xlabel="SOC", ylabel="Voltage (V)", title="qOCV reconstruction")
    axes[0].legend(fontsize=8)

    order = [(r, c) for r in ga.RATES for c in (True, False)]
    labels = [f"{r:g}C\n{'Chg' if c else 'Dch'}" for r, c in order]
    x = np.arange(len(order))
    width = 0.36
    for i, (name, color) in enumerate(((old_name, "#7f7f7f"), (new_name, "#1f77b4"))):
        vals = []
        caps = []
        for rate, charge in order:
            direction = "Charge" if charge else "Discharge"
            row = dynamic_all[
                (dynamic_all.Scenario == name)
                & np.isclose(dynamic_all.C_rate, rate)
                & (dynamic_all.Direction == direction)
            ].iloc[0]
            vals.append(row.MAE_SOC10_70_mV)
            caps.append(abs(row.capacity_error_Ah) * 1000.0)
        offset = (-0.5 if i == 0 else 0.5) * width
        axes[1].bar(x + offset, vals, width=width, color=color, label="Previous" if i == 0 else "Estimated")
        axes[2].bar(x + offset, caps, width=width, color=color, label="Previous" if i == 0 else "Estimated")
    axes[1].set(title="Dynamic voltage error", ylabel="MAE, SOC 10–70% (mV)")
    axes[2].set(title="End-capacity error magnitude", ylabel="|capacity error| (mAh)")
    for ax in axes[1:]:
        ax.set_xticks(x, labels)
        ax.legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(
        RESULTS / "teardown_thickness_effect_summary.png",
        dpi=220,
        bbox_inches="tight",
    )
    plt.close(fig)

    old_row = summary[summary.Scenario == old_name].iloc[0]
    new_row = summary[summary.Scenario == new_name].iloc[0]
    report = {
        "assumption": "Double-sided total thickness includes one current collector and two equal coating layers.",
        "collector_thickness_source": "Current Ai2020-based model values; replace after direct measurement.",
        "positive": {
            "double_sided_total_um": POSITIVE_TOTAL_DOUBLE_SIDED_M * 1e6,
            "assumed_current_collector_um": cc_p * 1e6,
            "estimated_single_side_coating_um": new_lp * 1e6,
        },
        "negative": {
            "double_sided_total_um": NEGATIVE_TOTAL_DOUBLE_SIDED_M * 1e6,
            "assumed_current_collector_um": cc_n * 1e6,
            "estimated_single_side_coating_um": new_ln * 1e6,
        },
        "method": "Thickness and c_s_max updated together; Qn and Qp fixed by measured half-cell capacity/coating mass.",
        "old": old_row.to_dict(),
        "estimated": new_row.to_dict(),
        "dynamic_delta_estimated_minus_old": {
            "MAE_SOC10_70_mV": float(new_row.Dynamic_MAE_SOC10_70_mV - old_row.Dynamic_MAE_SOC10_70_mV),
            "mean_abs_capacity_error_Ah": float(new_row.Mean_abs_capacity_error_Ah - old_row.Mean_abs_capacity_error_Ah),
            "max_abs_capacity_error_Ah": float(new_row.Max_abs_capacity_error_Ah - old_row.Max_abs_capacity_error_Ah),
        },
        "important_limitations": [
            "Current collector thickness is provisional until Monday measurement.",
            "Equal coating thickness on both sides is assumed.",
            "Porosity, particle radius, diffusion and kinetic parameters are not changed in this trial.",
            "Positive coating mass is treated as LiCoO2 mass as in the prior Chen-style trial.",
        ],
    }
    (RESULTS / "teardown_thickness_estimate_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    markdown = f"""# 전극 두께 1차 반영 결과

## 입력과 환산

- 양극 양면 총두께 138 µm, 모델 Al 집전체 15 µm 가정 → 단면 코팅 두께 **{new_lp * 1e6:.1f} µm**
- 음극 양면 총두께 155 µm, 모델 Cu 집전체 10 µm 가정 → 단면 코팅 두께 **{new_ln * 1e6:.1f} µm**
- 식: 단면 코팅 두께 = (양면 총두께 - 집전체 두께) / 2

## 물리적으로 일관된 반영

두께만 낮추면 전극 용량이 인위적으로 줄어든다. 따라서 기존 10 mm 펀칭 시편의 코팅 질량과 음극 반쪽전지 용량을 그대로 만족하도록 c_s,max를 역보정했다. 이 경우 L·ε·c_s,max가 보존되어 Qn/Qp, stoichiometry window, qOCV 재구성은 거의 변하지 않고 DFN의 두께 방향 수송만 달라진다.

| 항목 | 기존 | 추정 두께 반영 |
|---|---:|---:|
| 양극 단면 두께 (µm) | {old_lp * 1e6:.1f} | {new_lp * 1e6:.1f} |
| 음극 단면 두께 (µm) | {old_ln * 1e6:.1f} | {new_ln * 1e6:.1f} |
| 양극 c_s,max (mol/m³) | {old_row.Positive_c_s_max_mol_m3:.1f} | {new_row.Positive_c_s_max_mol_m3:.1f} |
| 음극 c_s,max (mol/m³) | {old_row.Negative_c_s_max_mol_m3:.1f} | {new_row.Negative_c_s_max_mol_m3:.1f} |
| Qp (Ah) | {old_row.Q_p_Ah:.6f} | {new_row.Q_p_Ah:.6f} |
| Qn (Ah) | {old_row.Q_n_Ah:.6f} | {new_row.Q_n_Ah:.6f} |
| qOCV MAE (mV) | {old_row.qOCV_MAE_2_98_mV:.3f} | {new_row.qOCV_MAE_2_98_mV:.3f} |
| 동적 평균 MAE (mV) | {old_row.Dynamic_MAE_SOC10_70_mV:.3f} | {new_row.Dynamic_MAE_SOC10_70_mV:.3f} |
| 평균 절대 용량 오차 (mAh) | {old_row.Mean_abs_capacity_error_Ah * 1000:.2f} | {new_row.Mean_abs_capacity_error_Ah * 1000:.2f} |

## 판단

집전체 실측 전까지는 이 값을 **임시 teardown geometry**로 사용할 수 있다. 다만 qOCV 개선을 위한 조정값은 아니다. qOCV는 전극별 총 활성물질 재고가 같으면 변하지 않는 것이 정상이다. 월요일 집전체 두께 측정 후 위 식에 실측값만 대입하고 동일 검증을 다시 수행해야 한다.
"""
    (RESULTS / "전극_두께_1차_반영_결과.md").write_text(markdown, encoding="utf-8")

    print(json.dumps(report, ensure_ascii=False, indent=2))
    print("\nDynamic validation\n", dynamic_all.to_string(index=False))
    print("\nSaved to", RESULTS)


if __name__ == "__main__":
    main()
