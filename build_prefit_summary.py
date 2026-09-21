"""Build the consolidated baseline summary before dynamic parameter fitting."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260921_prefit_summary"
RESULTS.mkdir(parents=True, exist_ok=True)

DYNAMIC_DIR = ROOT / "results" / "260921_c50_selected_ocp_dynamic_validation"
OCP_DIR = ROOT / "results" / "260921_c50_ocp_source_combinations"


def main() -> None:
    dynamic_json = json.loads(
        (DYNAMIC_DIR / "dynamic_validation.json").read_text(encoding="utf-8")
    )
    detail = pd.read_csv(DYNAMIC_DIR / "dynamic_condition_metrics.csv")
    detail = detail[detail["Scenario"] == "C50 OCP endpoint start"].copy()
    ocp = pd.read_csv(OCP_DIR / "ocp_source_combination_summary.csv")
    ocp = ocp[(ocp.anode_source == "nominal") & (ocp.cathode_source == "old")].iloc[0]
    fixed = dynamic_json["fixed_dynamic_parameters"]
    stage = dynamic_json["stage"]

    parameters = pd.DataFrame(
        [
            ("OCP", "Negative OCP", "Ai2020 nominal graphite", "fixed", "No anode hysteresis"),
            ("OCP", "Positive OCP", "Old cathode GITT", "fixed", "Cathode directional branches used"),
            ("Stoichiometry", "x0", stage["x0"], "identified", "C/50 qOCV"),
            ("Stoichiometry", "x100", stage["x100"], "identified", "C/50 qOCV"),
            ("Stoichiometry", "y100", stage["y100"], "identified", "C/50 qOCV"),
            ("Stoichiometry", "y0", stage["y0"], "identified", "C/50 qOCV"),
            ("Capacity", "Qcell", dynamic_json["qcell_C50_Ah"], "measured", "C/50 full-cell mean [Ah]"),
            ("Capacity", "Qn", dynamic_json["Qn_Ah"], "identified", "Qcell / delta x [Ah]"),
            ("Capacity", "Qp", dynamic_json["Qp_Ah"], "identified", "Qcell / delta y [Ah]"),
            ("Geometry", "Ln", fixed["Ln_um"], "fixed", "Ai2020 nominal [um]"),
            ("Geometry", "Lp", fixed["Lp_um"], "fixed", "Ai2020 nominal [um]"),
            ("Geometry", "Rn", fixed["Rn_um"], "fixed", "Ai2020 nominal [um]"),
            ("Geometry", "Rp", fixed["Rp_um"], "fixed", "Ai2020 nominal [um]"),
            ("Solid diffusion", "Dsn", fixed["Dsn_m2_s"], "fixed", "Specified [m2/s]"),
            ("Solid diffusion", "Dsp", fixed["Dsp_m2_s"], "fixed", "Specified [m2/s]"),
            ("Kinetics", "kn prefactor", fixed["kn_prefactor"], "fixed", "EIS N4-only"),
            ("Kinetics", "kp prefactor", fixed["kp_prefactor"], "fixed", "EIS first estimate"),
            ("Electrolyte", "Bruggeman n", fixed["Bruggeman_n_p_s"][0], "fixed", "Ai2020 nominal"),
            ("Electrolyte", "Bruggeman p", fixed["Bruggeman_n_p_s"][1], "fixed", "Ai2020 nominal"),
            ("Electrolyte", "Bruggeman s", fixed["Bruggeman_n_p_s"][2], "fixed", "Ai2020 nominal"),
            ("Electrolyte", "De multiplier", 1e-4, "fixed", "Ai2020 raw function unit conversion"),
            ("Effective inventory", "csn,max", fixed["effective_csn_max_mol_m3"], "derived", "Qn preserved [mol/m3]"),
            ("Effective inventory", "csp,max", fixed["effective_csp_max_mol_m3"], "derived", "Qp preserved [mol/m3]"),
        ],
        columns=["Group", "Parameter", "Value", "Status", "Basis"],
    )
    parameters.to_csv(RESULTS / "prefit_parameter_table.csv", index=False, encoding="utf-8-sig")

    dynamic = detail[
        [
            "C_rate",
            "Direction",
            "Center10_70_MAE_mV",
            "Center10_70_RMSE_mV",
            "Full_RMSE_mV",
            "Initial_voltage_error_mV",
            "Q_end_exp_Ah",
            "Q_end_model_Ah",
            "Capacity_error_mAh",
        ]
    ].copy()
    dynamic.to_csv(RESULTS / "prefit_dynamic_metrics.csv", index=False, encoding="utf-8-sig")

    labels = [f"{r:g}C\n{d}" for r, d in zip(dynamic.C_rate, dynamic.Direction)]
    colors = ["#D55E00" if d == "Charge" else "#0072B2" for d in dynamic.Direction]
    x = np.arange(len(dynamic))

    fig, axes = plt.subplots(2, 2, figsize=(13.2, 8.6), constrained_layout=True)

    ax = axes[0, 0]
    bars = ax.bar(x, dynamic.Center10_70_RMSE_mV, color=colors)
    ax.bar_label(bars, fmt="%.1f", padding=3, fontsize=9)
    ax.set_xticks(x, labels)
    ax.set_ylabel("RMSE [mV]")
    ax.set_title("Dynamic voltage error, SOC 10–70%")
    ax.grid(axis="y", alpha=0.25)

    ax = axes[0, 1]
    bars = ax.bar(x, dynamic.Capacity_error_mAh, color=colors)
    ax.bar_label(bars, fmt="%+.0f", padding=3, fontsize=9)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(x, labels)
    ax.set_ylabel("Model − experiment [mAh]")
    ax.set_title("CC cutoff-capacity error")
    ax.grid(axis="y", alpha=0.25)

    ax = axes[1, 0]
    names = ["qOCV\nMAE", "qOCV\nRMSE", "C/50 charge\nMAE", "C/50 discharge\nMAE"]
    values = [
        ocp.qOCV_MAE_2_98_mV,
        ocp.qOCV_RMSE_2_98_mV,
        ocp.C50_charge_MAE_10_90_mV,
        ocp.C50_discharge_MAE_10_90_mV,
    ]
    bars = ax.bar(names, values, color=["#4C78A8", "#72B7B2", "#D55E00", "#0072B2"])
    ax.bar_label(bars, fmt="%.1f", padding=3, fontsize=9)
    ax.set_ylabel("Error [mV]")
    ax.set_title("Low-rate OCP baseline")
    ax.grid(axis="y", alpha=0.25)

    ax = axes[1, 1]
    ax.hlines(1.0, stage["x0"], stage["x100"], color="#444444", lw=9)
    ax.scatter([stage["x0"], stage["x100"]], [1, 1], color="#111111", zorder=3)
    ax.hlines(0.0, stage["y100"], stage["y0"], color="#4C78A8", lw=9)
    ax.scatter([stage["y100"], stage["y0"]], [0, 0], color="#1F4E79", zorder=3)
    ax.text(stage["x0"], 1.15, f"x0={stage['x0']:.4f}", ha="left", fontsize=9)
    ax.text(stage["x100"], 1.15, f"x100={stage['x100']:.4f}", ha="right", fontsize=9)
    ax.text(stage["y100"], 0.15, f"y100={stage['y100']:.4f}", ha="left", fontsize=9)
    ax.text(stage["y0"], 0.15, f"y0={stage['y0']:.4f}", ha="right", fontsize=9)
    ax.set_yticks([0, 1], ["Positive", "Negative"])
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.35, 1.45)
    ax.set_xlabel("Stoichiometry")
    ax.set_title("Identified operating windows")
    ax.grid(axis="x", alpha=0.25)

    fig.suptitle("Pre-fitting baseline | OCP fixed, dynamic parameters not optimized", fontsize=15)
    fig.savefig(RESULTS / "prefit_summary_graph.png", dpi=220)
    plt.close(fig)

    dynamic_headers = list(dynamic.columns)
    dynamic_md = "| " + " | ".join(dynamic_headers) + " |\n"
    dynamic_md += "|" + "|".join(["---"] * len(dynamic_headers)) + "|\n"
    for row in dynamic.itertuples(index=False, name=None):
        cells = [f"{value:.3f}" if isinstance(value, (float, np.floating)) else str(value) for value in row]
        dynamic_md += "| " + " | ".join(cells) + " |\n"

    report = f"""# 동적 fitting 전 baseline 정리

## 범위

- OCP 및 stoichiometry window는 C/50 데이터로 선정된 상태이다.
- 아래 0.5C/1C/2C 결과에는 동적 데이터에 대한 parameter optimization을 아직 적용하지 않았다.
- 두께, 반경, Bruggeman은 Ai2020 nominal을 사용한다.
- Dsn/Dsp와 kn/kp는 지정값/EIS 값을 고정했다.

## OCP 및 용량

| 항목 | 값 |
|---|---:|
| OCP | Ai2020 nominal anode + old cathode GITT |
| Qcell | {dynamic_json['qcell_C50_Ah']:.5f} Ah |
| Qn | {dynamic_json['Qn_Ah']:.5f} Ah |
| Qp | {dynamic_json['Qp_Ah']:.5f} Ah |
| x0 / x100 | {stage['x0']:.6f} / {stage['x100']:.6f} |
| y100 / y0 | {stage['y100']:.6f} / {stage['y0']:.6f} |
| qOCV MAE / RMSE | {ocp.qOCV_MAE_2_98_mV:.2f} / {ocp.qOCV_RMSE_2_98_mV:.2f} mV |
| C/50 charge / discharge MAE | {ocp.C50_charge_MAE_10_90_mV:.2f} / {ocp.C50_discharge_MAE_10_90_mV:.2f} mV |

## 동적 사전 검증

{dynamic_md}

## 해석

- 중앙 SOC 전압 형상은 6조건 모두 약 7–18 mV RMSE 수준이다.
- 충전 CC cutoff 용량은 모든 C-rate에서 모델이 작으며, 동적 fitting의 핵심 잔차이다.
- 방전 용량오차는 충전보다 작으므로 Qn/Qp 전체를 동시에 확대하는 방식은 우선하지 않는다.
- 다음 fitting에서는 OCP/Qn/Qp를 우선 고정하고 actual current, cathode-side capacity/kinetics, 초기 이력을 분리한다.
"""
    (RESULTS / "README_KO.md").write_text(report, encoding="utf-8")
    print(parameters.to_string(index=False))
    print("\nDYNAMIC")
    print(dynamic.to_string(index=False))
    print("\nSaved", RESULTS)


if __name__ == "__main__":
    main()
