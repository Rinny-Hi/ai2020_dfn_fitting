"""Build the compact, Stage-1-only Colab notebook."""

from __future__ import annotations

import json
from pathlib import Path
from textwrap import dedent


ROOT = Path(__file__).resolve().parent
SOURCE_NOTEBOOK = ROOT / "261019_Ai2020_S1_S3_현대차제출용_간결본.ipynb"
OUTPUT = ROOT / "261020_Ai2020_1단계_OCP_Stoichiometry_간결본.ipynb"


def lines(text: str) -> list[str]:
    return (dedent(text).strip("\n") + "\n").splitlines(keepends=True)


def markdown(text: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": lines(text)}


def code(text: str) -> dict:
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": lines(text),
    }


old = json.loads(SOURCE_NOTEBOOK.read_text(encoding="utf-8"))


def old_source(index: int) -> str:
    return "".join(old["cells"][index]["source"])


# Reuse only the already validated raw-data parsing and DFN evaluation helpers.
data_loading = old_source(2)
dynamic_helpers = old_source(12)

ocp_source = old_source(6)
duration_helpers = ocp_source[
    ocp_source.index("def duration_to_seconds") : ocp_source.index(
        "def extract_long_half_cell_branch"
    )
]
cathode_builder = ocp_source[
    ocp_source.index("def build_cathode_gitt_ocp") : ocp_source.index(
        "ANODE_LOW_RATE_OCP"
    )
]

fit_source = old_source(9)
fit_helpers = fit_source[: fit_source.index("s2_vector =")]


cells = [
    markdown(
        r"""
        # Ai2020 기반 DFN 모델 — 1단계 간결본

        **목적:** 공통 물리 보정 후, 양극 OCP 교체와 stoichiometry window 재최적화의 효과만 검증한다.

        1. 공통 기반: Ai2020 `D_e` 단위 보정 + 측정 `R_n/R_p`
        2. 기준: Ai2020 음극/양극 OCP + Ai2020 공개 window
        3. 개선: Ai2020 음극 OCP + 실험 양극 GITT OCP + window 재최적화
        4. 검증: C/50 qOCV 및 0.5C/1C/2C 충·방전 전압·cut-off 용량

        `D_s`, `k`, Bruggeman 계수, 전해질 parameter set 및 fitting은 이 노트북의 범위가 아니다.
        """
    ),
    markdown(
        r"""
        ## 0. 실행 환경과 입력 파일

        아래 5개 원본 파일을 노트북과 같은 폴더에 둔다.

        - `260610 저율데이터 1.xlsx`
        - `260610 저율데이터 2.xlsx`
        - `260808 Eneterch셀 열화데이터 v1 6-6.xlsx`
        - `260816 Positive GITT 7-7.xlsx`
        - `260816 Positive GITT 7-8.xlsx`

        음극 OCP는 Ai2020을 사용하므로 음극 half-cell 파일은 필요하지 않다.
        """
    ),
    code(
        r'''
        # Colab 최초 실행 시 주석을 해제한다.
        # PyBaMM 26.7.1의 Ai2020 D_e를 명시적으로 보정해 기존 분석을 재현한다.
        # %pip install -q pybamm==26.7.1.0 openpyxl

        from pathlib import Path
        import warnings

        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
        import pybamm
        from IPython.display import display
        from scipy.integrate import cumulative_trapezoid
        from scipy.optimize import differential_evolution, minimize

        warnings.filterwarnings("ignore", message="Workbook contains no default style")
        REQUIRED_PYBAMM = "26.7.1.0"
        if pybamm.__version__ != REQUIRED_PYBAMM:
            raise RuntimeError(
                f"PyBaMM {REQUIRED_PYBAMM}이 필요합니다. 현재 버전: {pybamm.__version__}"
            )

        FULL_CELL_LOW_RATE_FILES = [
            "260610 저율데이터 1.xlsx",
            "260610 저율데이터 2.xlsx",
        ]
        DYNAMIC_FILE = "260808 Eneterch셀 열화데이터 v1 6-6.xlsx"
        CATHODE_GITT_FILES = [
            "260816 Positive GITT 7-7.xlsx",
            "260816 Positive GITT 7-8.xlsx",
        ]
        required_files = FULL_CELL_LOW_RATE_FILES + [DYNAMIC_FILE] + CATHODE_GITT_FILES
        missing_files = [name for name in required_files if not Path(name).exists()]
        if missing_files:
            raise FileNotFoundError("업로드가 필요한 파일: " + ", ".join(missing_files))

        print("PyBaMM:", pybamm.__version__)
        print("입력 파일 확인:", len(required_files), "개")
        '''
    ),
    code(data_loading),
    markdown(
        r"""
        ## 1. 공통 기반 고정

        OCP 비교 전에 두 모델에 동일하게 적용한다.

        - Ai2020 electrolyte diffusivity: 원 구현의 단위 해석에 따라 `× 10⁻⁴`
        - 음극 평균 입자반경: `3.7542 μm`
        - 양극 평균 입자반경: `4.2387 μm`
        - 전극 겹침 면적: `45.4 mm × 51.0 mm × 34면`

        따라서 뒤의 차이는 입자반경 변경이 아니라 OCP/window 변경에서 발생한다.
        """
    ),
    code(
        r'''
        FARADAY_CONSTANT = 96485.33212
        NEGATIVE_RADIUS_M = 3.7542e-6
        POSITIVE_RADIUS_M = 4.2387e-6
        OVERLAP_AREA_M2 = 45.4e-3 * 51.0e-3 * 34.0

        raw_parameters = pybamm.ParameterValues("Ai2020")
        raw_electrolyte_diffusivity = raw_parameters[
            "Electrolyte diffusivity [m2.s-1]"
        ]


        def corrected_electrolyte_diffusivity(concentration, temperature):
            return 1e-4 * raw_electrolyte_diffusivity(concentration, temperature)


        common_parameters = raw_parameters.copy()
        common_parameters.update(
            {
                "Electrolyte diffusivity [m2.s-1]": corrected_electrolyte_diffusivity,
                "Negative particle radius [m]": NEGATIVE_RADIUS_M,
                "Positive particle radius [m]": POSITIVE_RADIUS_M,
                "Electrode height [m]": 51.0e-3,
                "Electrode width [m]": 45.4e-3,
                "Number of electrodes connected in parallel to make a cell": 34.0,
            }
        )

        common_changes = pd.DataFrame(
            [
                ["Electrolyte diffusivity", "Ai2020 function", "Ai2020 function × 1e-4", "unit correction"],
                ["Negative particle radius", f"{1e6 * float(raw_parameters['Negative particle radius [m]']):.4f} μm", "3.7542 μm", "SEM/ImageJ mean"],
                ["Positive particle radius", f"{1e6 * float(raw_parameters['Positive particle radius [m]']):.4f} μm", "4.2387 μm", "SEM/ImageJ mean"],
            ],
            columns=["Parameter", "Ai2020", "Stage-1 common value", "Basis"],
        )
        display(common_changes)
        '''
    ),
    markdown(
        r"""
        ## 2. OCP와 stoichiometry window

        음극은 absolute stoichiometry 축이 정의된 Ai2020 graphite OCP를 유지한다.
        양극은 두 GITT 반복시험의 장시간 완화전압을 평균해 사용한다.

        Window 최적화 목적함수는 C/50 충·방전 평균 qOCV의 2–98% SOC RMSE이며,
        SOC 0/100% 끝점 오차는 각각 10 mV 이내로 제한한다.

        양극 OCP 비교 범위는 최적화된 full-cell 3.0–4.2 V 운전 window의
        `y100–y0`로 제한한다. 4.2 V는 full-cell 단자전압이므로 양극 OCP를
        `4.2 V vs Li/Li+`에서 직접 절단하지 않는다.
        """
    ),
    code(
        r'''
        PUBLISHED_WINDOW = {"x0": 0.0065, "x100": 0.8400, "y0": 0.9651, "y100": 0.4350}


        def evaluate_parameter_function(parameter_values, name, grid):
            function = parameter_values[name]
            return np.array(
                [
                    float(
                        np.asarray(
                            function(pybamm.Scalar(float(value))).evaluate(), dtype=float
                        ).reshape(-1)[0]
                    )
                    for value in grid
                ]
            )


        ocp_grid = np.linspace(0.001, 0.999, 3001)
        AI2020_OCP = {
            "negative": (
                ocp_grid,
                evaluate_parameter_function(common_parameters, "Negative electrode OCP [V]", ocp_grid),
            ),
            "positive": (
                ocp_grid,
                evaluate_parameter_function(common_parameters, "Positive electrode OCP [V]", ocp_grid),
            ),
        }
        '''
    ),
    code(duration_helpers + "\n" + cathode_builder + dedent(r'''

        CATHODE_GITT_OCP = build_cathode_gitt_ocp(
            CATHODE_GITT_FILES, raw_parameters
        )
        ''')),
    code(fit_helpers + dedent(r'''

        def vector_to_window(vector):
            x0, x100, y100, y0 = map(float, vector)
            return {"x0": x0, "x100": x100, "y0": y0, "y100": y100}


        published_vector = np.array([
            PUBLISHED_WINDOW["x0"], PUBLISHED_WINDOW["x100"],
            PUBLISHED_WINDOW["y100"], PUBLISHED_WINDOW["y0"],
        ])
        ai_refit_vector = fit_stoichiometry_window(
            AI2020_OCP["negative"], AI2020_OCP["positive"], random_seed=260924
        )
        improved_vector = fit_stoichiometry_window(
            AI2020_OCP["negative"], CATHODE_GITT_OCP, random_seed=260921
        )
        AI_REFIT_WINDOW = vector_to_window(ai_refit_vector)
        IMPROVED_WINDOW = vector_to_window(improved_vector)
        ''')),
    code(
        r'''
        # 양극 OCP는 full-cell 3.0–4.2 V에 대응하는 최종 stoichiometry
        # window에서 비교한다. 4.2 V는 full-cell 단자전압이므로 양극 OCP를
        # 4.2 V vs Li/Li+에서 수평 절단하지 않는다.
        y_charge_side = IMPROVED_WINDOW["y100"]
        y_discharge_side = IMPROVED_WINDOW["y0"]
        comparison_y = np.linspace(y_charge_side, y_discharge_side, 1000)
        ai2020_positive_in_window = np.interp(
            comparison_y, *AI2020_OCP["positive"]
        )
        measured_positive_in_window = np.interp(
            comparison_y, *CATHODE_GITT_OCP
        )

        figure, axis = plt.subplots(figsize=(8.2, 4.8), constrained_layout=True)
        # 범위 밖 곡선은 출처와 측정 범위 확인용으로만 흐리게 표시한다.
        axis.plot(
            *AI2020_OCP["positive"], color="#777777", lw=1.2, alpha=0.25
        )
        axis.plot(
            *CATHODE_GITT_OCP, color="#D55E00", lw=1.2, alpha=0.25
        )
        axis.axvspan(
            y_charge_side,
            y_discharge_side,
            color="#56B4E9",
            alpha=0.09,
            label="Full-cell 3.0–4.2 V stoichiometry window",
        )
        axis.plot(
            comparison_y,
            ai2020_positive_in_window,
            color="#666666",
            lw=2.3,
            label="Ai2020 cathode OCP",
        )
        axis.plot(
            comparison_y,
            measured_positive_in_window,
            color="#D55E00",
            lw=2.5,
            label="Measured cathode GITT OCP",
        )
        axis.axvline(
            y_charge_side,
            color="#CC3311",
            ls="--",
            lw=1.5,
            label=f"Full-cell 4.2 V side: y100={y_charge_side:.4f}",
        )
        axis.axvline(
            y_discharge_side,
            color="#0072B2",
            ls="--",
            lw=1.5,
            label=f"Full-cell 3.0 V side: y0={y_discharge_side:.4f}",
        )
        margin = 0.03
        axis.set_xlim(
            max(float(CATHODE_GITT_OCP[0].min()), y_charge_side - margin),
            min(float(CATHODE_GITT_OCP[0].max()), y_discharge_side + margin),
        )
        axis.set(
            xlabel="Cathode stoichiometry y",
            ylabel="OCP vs Li/Li+ [V]",
            title="Cathode OCP over the full-cell operating window",
        )
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
        plt.show()

        # qOCV: 공개 window → window만 재최적화 → 실험 양극 OCP+window 재최적화
        qocv_cases = [
            ("Ai2020 published", published_vector, AI2020_OCP["negative"], AI2020_OCP["positive"], "#999999", "--"),
            ("Ai2020 OCP + refit window", ai_refit_vector, AI2020_OCP["negative"], AI2020_OCP["positive"], "#0072B2", "-"),
            ("Measured cathode OCP + refit window", improved_vector, AI2020_OCP["negative"], CATHODE_GITT_OCP, "#D55E00", "-"),
        ]
        fitting_region = (SOC_GRID >= 0.02) & (SOC_GRID <= 0.98)
        qocv_rows, qocv_predictions = [], {}

        for name, vector, anode_ocp, cathode_ocp, color, style in qocv_cases:
            prediction = predict_full_cell_ocv(vector, anode_ocp, cathode_ocp)
            qocv_predictions[name] = prediction
            error_mv = 1000.0 * (prediction - low_rate_data["qocv"])
            x0, x100, y100, y0 = vector
            qocv_rows.append(
                {
                    "Case": name,
                    "x0": x0, "x100": x100, "y0": y0, "y100": y100,
                    "Qn [Ah]": CELL_CAPACITY_AH / (x100 - x0),
                    "Qp [Ah]": CELL_CAPACITY_AH / (y0 - y100),
                    "qOCV RMSE 2–98% [mV]": np.sqrt(np.mean(error_mv[fitting_region] ** 2)),
                    "SOC0 error [mV]": error_mv[0],
                    "SOC100 error [mV]": error_mv[-1],
                }
            )

        window_results = pd.DataFrame(qocv_rows)
        display(window_results.round(6))

        figure, axes = plt.subplots(1, 2, figsize=(13, 4.5), constrained_layout=True)
        axes[0].plot(SOC_GRID * 100, low_rate_data["qocv"], color="black", lw=2.7, label="Measured C/50 qOCV")
        axes[1].axhline(0.0, color="black", lw=0.8)
        for name, _, _, _, color, style in qocv_cases:
            prediction = qocv_predictions[name]
            axes[0].plot(SOC_GRID * 100, prediction, color=color, ls=style, lw=2.0, label=name)
            axes[1].plot(SOC_GRID * 100, 1000 * (prediction - low_rate_data["qocv"]), color=color, ls=style, lw=1.8, label=name)
        axes[0].set(title="qOCV reconstruction", xlabel="Full-cell SOC [%]", ylabel="Voltage [V]")
        axes[1].set(title="qOCV error", xlabel="Full-cell SOC [%]", ylabel="Model − experiment [mV]")
        for axis in axes:
            axis.grid(alpha=0.25)
            axis.legend(fontsize=8)
        plt.show()
        '''
    ),
    markdown(
        r"""
        ## 3. Capacity-consistent inventory와 DFN 검증

        개선 window에서 계산된 `Q_n`, `Q_p`를 같은 전극 면적·두께·활물질 분율로 구현하도록
        `c_s,max,n/p`를 역산한다. 이 값은 독립 fitting 변수가 아니라 window와 용량을 연결하는 bookkeeping이다.
        """
    ),
    code(
        r'''
        def build_improved_parameters(base_parameters, window, positive_ocp):
            parameters = base_parameters.copy()
            qn_ah = CELL_CAPACITY_AH / (window["x100"] - window["x0"])
            qp_ah = CELL_CAPACITY_AH / (window["y0"] - window["y100"])
            csmax_n = 3600.0 * qn_ah / (
                FARADAY_CONSTANT * OVERLAP_AREA_M2
                * float(parameters["Negative electrode thickness [m]"])
                * float(parameters["Negative electrode active material volume fraction"])
            )
            csmax_p = 3600.0 * qp_ah / (
                FARADAY_CONSTANT * OVERLAP_AREA_M2
                * float(parameters["Positive electrode thickness [m]"])
                * float(parameters["Positive electrode active material volume fraction"])
            )
            positive_grid, positive_voltage = positive_ocp

            def measured_positive_ocp(stoichiometry):
                return pybamm.Interpolant(
                    positive_grid, positive_voltage, stoichiometry,
                    name="Measured cathode GITT OCP", interpolator="linear"
                )

            parameters.update(
                {
                    "Maximum concentration in negative electrode [mol.m-3]": csmax_n,
                    "Maximum concentration in positive electrode [mol.m-3]": csmax_p,
                    "Positive electrode OCP [V]": measured_positive_ocp,
                }
            )
            return parameters, csmax_n, csmax_p


        improved_parameters, improved_csmax_n, improved_csmax_p = build_improved_parameters(
            common_parameters, IMPROVED_WINDOW, CATHODE_GITT_OCP
        )
        baseline_csmax_n = float(common_parameters["Maximum concentration in negative electrode [mol.m-3]"])
        baseline_csmax_p = float(common_parameters["Maximum concentration in positive electrode [mol.m-3]"])

        parameter_summary = pd.DataFrame(
            [
                ["Negative OCP", "Ai2020", "Ai2020", "absolute x-axis retained"],
                ["Positive OCP", "Ai2020", "Measured GITT mean", "cell-specific relaxed voltage"],
                ["x0 / x100", f"{PUBLISHED_WINDOW['x0']:.6f} / {PUBLISHED_WINDOW['x100']:.6f}", f"{IMPROVED_WINDOW['x0']:.6f} / {IMPROVED_WINDOW['x100']:.6f}", "C/50 qOCV fit"],
                ["y0 / y100", f"{PUBLISHED_WINDOW['y0']:.6f} / {PUBLISHED_WINDOW['y100']:.6f}", f"{IMPROVED_WINDOW['y0']:.6f} / {IMPROVED_WINDOW['y100']:.6f}", "C/50 qOCV fit"],
                ["c_s,max,n [mol m-3]", f"{baseline_csmax_n:.1f}", f"{improved_csmax_n:.1f}", "Qn/window consistency"],
                ["c_s,max,p [mol m-3]", f"{baseline_csmax_p:.1f}", f"{improved_csmax_p:.1f}", "Qp/window consistency"],
            ],
            columns=["Item", "Baseline", "Improved", "Basis"],
        )
        display(parameter_summary)
        '''
    ),
    code(dynamic_helpers),
    code(
        r'''
        stages = {
            "Baseline: Ai2020 OCP/window": {
                "parameters": common_parameters,
                "window": PUBLISHED_WINDOW,
                "c_s_max_n": baseline_csmax_n,
                "c_s_max_p": baseline_csmax_p,
                "anode_ocp": AI2020_OCP["negative"],
                "cathode_ocp": AI2020_OCP["positive"],
            },
            "Improved: measured cathode OCP/refit window": {
                "parameters": improved_parameters,
                "window": IMPROVED_WINDOW,
                "c_s_max_n": improved_csmax_n,
                "c_s_max_p": improved_csmax_p,
                "anode_ocp": AI2020_OCP["negative"],
                "cathode_ocp": CATHODE_GITT_OCP,
            },
        }

        simulations, metric_rows = {}, []
        for stage_name, stage in stages.items():
            for c_rate in TEST_C_RATES:
                for is_charge in (True, False):
                    direction = "Charge" if is_charge else "Discharge"
                    print(f"{stage_name} | {c_rate:g}C {direction}")
                    simulation = run_dfn_branch(stage, c_rate, is_charge)
                    simulations[(stage_name, c_rate, is_charge)] = simulation
                    metric_rows.append(
                        {
                            "Model": stage_name,
                            "C-rate": c_rate,
                            "Direction": direction,
                            **compare_with_experiment(
                                experimental_curves[(c_rate, is_charge)], simulation
                            ),
                        }
                    )

        metrics = pd.DataFrame(metric_rows)
        display(metrics.round(3))
        '''
    ),
    code(
        r'''
        def rms(values):
            values = np.asarray(values, dtype=float)
            return float(np.sqrt(np.mean(values**2)))


        model_summary = (
            metrics.groupby("Model", as_index=False)
            .agg(
                Voltage_full_RMSE_mV=("full_RMSE_mV", rms),
                Voltage_10_70_RMSE_mV=("center_10_70_RMSE_mV", rms),
                Cutoff_capacity_RMSE_pct=("capacity_error_pct", rms),
                Max_abs_cutoff_capacity_error_pct=("capacity_error_pct", lambda x: float(np.max(np.abs(x)))),
            )
        )
        display(model_summary.round(3))

        colors = {
            "Baseline: Ai2020 OCP/window": "#888888",
            "Improved: measured cathode OCP/refit window": "#D55E00",
        }
        figure, axes = plt.subplots(2, 3, figsize=(15.5, 8.0), sharey=True, constrained_layout=True)
        for column, c_rate in enumerate(TEST_C_RATES):
            for row, is_charge in enumerate((True, False)):
                direction = "Charge" if is_charge else "Discharge"
                axis = axes[row, column]
                experiment = experimental_curves[(c_rate, is_charge)]
                axis.plot(experiment["time_min"], experiment["voltage_v"], color="black", lw=2.7, label="Experiment")
                for stage_name in stages:
                    simulation = simulations[(stage_name, c_rate, is_charge)]
                    branch = metrics[
                        (metrics["Model"] == stage_name)
                        & np.isclose(metrics["C-rate"], c_rate)
                        & (metrics["Direction"] == direction)
                    ].iloc[0]
                    axis.plot(
                        simulation["time_min"], simulation["voltage_v"],
                        color=colors[stage_name], lw=2.0,
                        ls="--" if stage_name.startswith("Baseline") else "-",
                        label=f"{stage_name.split(':')[0]} ({branch.full_RMSE_mV:.1f} mV)",
                    )
                axis.set_title(f"{c_rate:g}C {direction}", weight="bold")
                axis.set_ylim(2.9, 4.3)
                axis.grid(alpha=0.25)
                if row == 1:
                    axis.set_xlabel("Elapsed time [min]")
                if column == 0:
                    axis.set_ylabel("Cell voltage [V]")
        handles, labels = axes[0, 0].get_legend_handles_labels()
        figure.legend(handles, labels, loc="outside upper center", ncol=3, frameon=False)
        figure.suptitle("Stage 1: baseline vs OCP/window improvement", weight="bold")
        plt.show()
        '''
    ),
    markdown(
        r"""
        ## 해석 원칙

        - qOCV 표는 **window 재최적화 효과**와 **양극 OCP 교체의 추가 효과**를 분리해 보여준다.
        - 동적 표의 전압 오차는 실험·모델 공통 시간 구간의 RMSE이다.
        - `cut-off capacity error`는 4.2/3.0 V 도달 용량의 상대오차이며, 독립적인 용량 예측 오차로 보고한다.
        - 개선 모델에서 양극 표면 stoichiometry가 GITT OCP 상한을 미세하게 넘으면 PyBaMM 경고가 발생할 수 있다. 경고를 숨기지 않고 OCP 측정 범위의 한계로 남긴다.
        - 이 단계에서 `D_s`, `k`, Bruggeman 계수는 변경하지 않았으므로 이후 단계의 fitting 효과와 섞이지 않는다.
        """
    ),
]

for index, cell in enumerate(cells):
    cell["id"] = f"stage1-compact-{index:02d}"

notebook = {
    "cells": cells,
    "metadata": {
        "kernelspec": {
            "display_name": "Python 3 (ipykernel)",
            "language": "python",
            "name": "python3",
        },
        "language_info": {"name": "python", "version": "3.12"},
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}

OUTPUT.write_text(json.dumps(notebook, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
print(OUTPUT)
