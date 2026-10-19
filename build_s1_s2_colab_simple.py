"""Build a compact, Colab-friendly S1/S2 DFN notebook."""

from __future__ import annotations

import json
from pathlib import Path
from textwrap import dedent


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "261019_Ai2020_S1_S2_현대차제출용_간결본.ipynb"


def source_lines(text: str) -> list[str]:
    return (dedent(text).strip("\n") + "\n").splitlines(keepends=True)


def markdown(text: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": source_lines(text)}


def code(text: str) -> dict:
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": source_lines(text),
    }


cells = [
    markdown(
        r"""
        # Ai2020 기반 DFN 모델: S1 → S2

        - **[0]** PyBaMM 실행 및 실험 데이터 불러오기
        - **[1] S1** Ai2020 원형 parameter set
        - **[2] S2** 전해질 확산계수 (D_e) 단위 보정 + 확정 stoichiometry window 적용

        Colab의 파일 창에 아래 Excel 파일 3개를 업로드한 후 셀을 순서대로 실행한다.
        """
    ),
    code(
        r'''
        # [0-1] Colab 최초 실행 시 아래 한 줄의 주석을 해제한다.
        # %pip install -q pybamm openpyxl

        from pathlib import Path
        import warnings

        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
        import pybamm
        import pybamm.input.parameters.lithium_ion.Ai2020 as ai2020_data
        from IPython.display import display
        from scipy.integrate import cumulative_trapezoid

        warnings.filterwarnings("ignore", message="Workbook contains no default style")
        print("PyBaMM version:", pybamm.__version__)

        # Colab에 업로드할 파일 이름
        LOW_RATE_FILES = [
            "260610 저율데이터 1.xlsx",
            "260610 저율데이터 2.xlsx",
        ]
        DYNAMIC_DATA_FILE = "260902 Enertech셀 데이터 모음.xlsx"

        required_files = [*LOW_RATE_FILES, DYNAMIC_DATA_FILE]
        missing_files = [name for name in required_files if not Path(name).exists()]
        if missing_files:
            raise FileNotFoundError("먼저 Colab에 업로드하세요: " + ", ".join(missing_files))
        '''
    ),
    code(
        r'''
        # [0-2] 저율 qOCV와 0.5C/1C/2C 실험 데이터 불러오기
        TEST_C_RATES = (0.5, 1.0, 2.0)
        SOC_GRID = np.linspace(0.0, 1.0, 1001)


        def clean_curve(horizontal, vertical):
            horizontal = np.asarray(horizontal, dtype=float).reshape(-1)
            vertical = np.asarray(vertical, dtype=float).reshape(-1)
            valid = np.isfinite(horizontal) & np.isfinite(vertical)
            order = np.argsort(horizontal[valid], kind="stable")
            horizontal = horizontal[valid][order]
            vertical = vertical[valid][order]
            horizontal, unique_index = np.unique(horizontal, return_index=True)
            return horizontal, vertical[unique_index]


        def load_low_rate_data(file_names, soc_grid):
            charge_curves, discharge_curves = [], []
            charge_capacities, discharge_capacities = [], []

            for file_name in file_names:
                data = pd.read_excel(file_name, sheet_name="record")
                contiguous_step = data["Step Type"].astype(str).ne(
                    data["Step Type"].astype(str).shift()
                ).cumsum()
                branches = {"CC Chg": [], "CC DChg": []}

                for _, branch in data.groupby(contiguous_step):
                    step_type = str(branch["Step Type"].iloc[0])
                    if step_type not in branches:
                        continue
                    capacity = pd.to_numeric(
                        branch["Capacity(Ah)"], errors="coerce"
                    ).to_numpy(float)
                    voltage = pd.to_numeric(
                        branch["Voltage(V)"], errors="coerce"
                    ).to_numpy(float)
                    current = float(
                        pd.to_numeric(branch["Current(A)"], errors="coerce").median()
                    )
                    valid = np.isfinite(capacity) & np.isfinite(voltage)
                    capacity, voltage = capacity[valid], voltage[valid]
                    if len(capacity) < 50 or capacity.max() < 2.2 or abs(current) > 0.06:
                        continue

                    end_capacity = float(capacity.max())
                    branch_soc = (
                        capacity / end_capacity
                        if step_type == "CC Chg"
                        else 1.0 - capacity / end_capacity
                    )
                    branch_soc, voltage = clean_curve(branch_soc, voltage)
                    branches[step_type].append(
                        (end_capacity, np.interp(soc_grid, branch_soc, voltage))
                    )

                # 첫 방전은 부분 초기충전 직후이므로 제외: 각 셀 complete cycles 2~5 사용
                for capacity, voltage in branches["CC Chg"]:
                    charge_capacities.append(capacity)
                    charge_curves.append(voltage)
                for capacity, voltage in branches["CC DChg"][1:]:
                    discharge_capacities.append(capacity)
                    discharge_curves.append(voltage)

            mean_charge = np.mean(charge_curves, axis=0)
            mean_discharge = np.mean(discharge_curves, axis=0)
            cell_capacity = 0.5 * (
                np.mean(charge_capacities) + np.mean(discharge_capacities)
            )
            return {
                "charge": mean_charge,
                "discharge": mean_discharge,
                "qocv": 0.5 * (mean_charge + mean_discharge),
                "capacity_ah": float(cell_capacity),
                "number_of_charge_curves": len(charge_curves),
                "number_of_discharge_curves": len(discharge_curves),
            }


        def load_dynamic_data(file_name, cell_capacity_ah):
            data = pd.read_excel(file_name, sheet_name="DYNAMIC", header=8)
            experiment = {}
            for c_rate in TEST_C_RATES:
                for is_charge in (True, False):
                    direction = "charge" if is_charge else "discharge"
                    branch = data[
                        np.isclose(pd.to_numeric(data["rate_C"], errors="coerce"), c_rate)
                        & (data["direction"].astype(str).str.lower() == direction)
                    ].copy()
                    capacity = pd.to_numeric(
                        branch["Q_Ah"], errors="coerce"
                    ).to_numpy(float)
                    capacity = capacity - capacity[0]
                    experiment[(c_rate, is_charge)] = {
                        "time_min": pd.to_numeric(
                            branch["t_min"], errors="coerce"
                        ).to_numpy(float),
                        "voltage_v": pd.to_numeric(
                            branch["V_V"], errors="coerce"
                        ).to_numpy(float),
                        "capacity_ah": capacity,
                    }
            return experiment


        low_rate_data = load_low_rate_data(LOW_RATE_FILES, SOC_GRID)
        CELL_CAPACITY_AH = low_rate_data["capacity_ah"]
        experimental_curves = load_dynamic_data(DYNAMIC_DATA_FILE, CELL_CAPACITY_AH)

        print(f"C/50 charge curves    : {low_rate_data['number_of_charge_curves']}")
        print(f"C/50 discharge curves : {low_rate_data['number_of_discharge_curves']}")
        print(f"Cell capacity          : {CELL_CAPACITY_AH:.6f} Ah")
        '''
    ),
    code(
        r'''
        # [0-3] DFN 실행에 공통으로 사용하는 함수
        FARADAY_CONSTANT = 96485.33212
        NOMINAL_CAPACITY_AH = 2.28

        MEASURED_PROTOCOL = {
            (0.5, True):  {"rest_voltage_v": 3.0601, "current_a": 1.139948},
            (0.5, False): {"rest_voltage_v": 4.1783, "current_a": 1.139974},
            (1.0, True):  {"rest_voltage_v": 3.2073, "current_a": 2.279965},
            (1.0, False): {"rest_voltage_v": 4.1780, "current_a": 2.279931},
            (2.0, True):  {"rest_voltage_v": 3.3263, "current_a": 4.559894},
            (2.0, False): {"rest_voltage_v": 4.1787, "current_a": 4.560041},
        }


        def ai2020_ocp_arrays():
            _, (negative_stoichiometry, negative_voltage) = (
                ai2020_data.graphite_ocp_Enertech_Ai2020_data
            )
            _, (positive_stoichiometry, positive_voltage) = (
                ai2020_data.lico2_ocp_Ai2020_data
            )
            negative_stoichiometry, negative_voltage = clean_curve(
                negative_stoichiometry, negative_voltage
            )
            positive_stoichiometry, positive_voltage = clean_curve(
                positive_stoichiometry, positive_voltage
            )
            return (
                negative_stoichiometry,
                negative_voltage,
                positive_stoichiometry,
                positive_voltage,
            )


        def full_cell_ocv(window):
            negative_axis, negative_voltage, positive_axis, positive_voltage = (
                ai2020_ocp_arrays()
            )
            negative_stoichiometry = window["x0"] + SOC_GRID * (
                window["x100"] - window["x0"]
            )
            positive_stoichiometry = window["y0"] - SOC_GRID * (
                window["y0"] - window["y100"]
            )
            return np.interp(
                positive_stoichiometry, positive_axis, positive_voltage
            ) - np.interp(
                negative_stoichiometry, negative_axis, negative_voltage
            )


        def initial_soc_from_rest_voltage(window, rest_voltage_v):
            qocv = full_cell_ocv(window)
            order = np.argsort(qocv, kind="stable")
            unique_voltage, unique_index = np.unique(qocv[order], return_index=True)
            return float(np.interp(rest_voltage_v, unique_voltage, SOC_GRID[order][unique_index]))


        def run_dfn_branch(parameters, window, c_s_max_n, c_s_max_p, c_rate, is_charge):
            protocol = MEASURED_PROTOCOL[(c_rate, is_charge)]
            initial_soc = initial_soc_from_rest_voltage(
                window, protocol["rest_voltage_v"]
            )
            initial_x = window["x0"] + initial_soc * (window["x100"] - window["x0"])
            initial_y = window["y0"] - initial_soc * (window["y0"] - window["y100"])

            branch_parameters = parameters.copy()
            branch_parameters.update(
                {
                    "Initial concentration in negative electrode [mol.m-3]": initial_x * c_s_max_n,
                    "Initial concentration in positive electrode [mol.m-3]": initial_y * c_s_max_p,
                }
            )

            direction = "Charge" if is_charge else "Discharge"
            cutoff_voltage = 4.2 if is_charge else 3.0
            actual_c_rate = protocol["current_a"] / NOMINAL_CAPACITY_AH

            model = pybamm.lithium_ion.DFN({"thermal": "isothermal"})
            experiment = pybamm.Experiment(
                [f"{direction} at {actual_c_rate}C until {cutoff_voltage} V"],
                period="10 seconds",
            )
            simulation = pybamm.Simulation(
                model,
                parameter_values=branch_parameters,
                experiment=experiment,
                var_pts={"x_n": 15, "x_s": 15, "x_p": 15, "r_n": 15, "r_p": 15},
                solver=pybamm.IDAKLUSolver(rtol=1e-6, atol=1e-8),
            )
            solution = simulation.solve()
            time_s = np.asarray(solution["Time [s]"].entries, dtype=float)
            voltage_v = np.asarray(solution["Terminal voltage [V]"].entries, dtype=float)
            current_a = np.asarray(solution["Current [A]"].entries, dtype=float)
            capacity_ah = cumulative_trapezoid(
                np.abs(current_a), time_s, initial=0.0
            ) / 3600.0
            return {
                "time_min": time_s / 60.0,
                "voltage_v": voltage_v,
                "capacity_ah": capacity_ah,
                "initial_soc_pct": 100.0 * initial_soc,
            }


        def compare_with_experiment(experiment, simulation, measured_current_a):
            experiment_time = experiment["time_min"] - experiment["time_min"][0]
            simulation_time = simulation["time_min"] - simulation["time_min"][0]
            common_end = 0.995 * min(experiment_time[-1], simulation_time[-1])
            full_grid = np.linspace(0.0, common_end, 400)
            center_start = 0.10 * CELL_CAPACITY_AH / measured_current_a * 60.0
            center_end = min(
                0.70 * CELL_CAPACITY_AH / measured_current_a * 60.0,
                common_end,
            )
            center_grid = np.linspace(center_start, center_end, 160)

            def error_mv(time_grid):
                return 1000.0 * (
                    np.interp(time_grid, simulation_time, simulation["voltage_v"])
                    - np.interp(time_grid, experiment_time, experiment["voltage_v"])
                )

            full_error = error_mv(full_grid)
            center_error = error_mv(center_grid)
            experiment_capacity = experiment["capacity_ah"][-1]
            simulation_capacity = simulation["capacity_ah"][-1]
            return {
                "full_RMSE_mV": float(np.sqrt(np.mean(full_error**2))),
                "center_10_70_RMSE_mV": float(np.sqrt(np.mean(center_error**2))),
                "capacity_error_pct": float(
                    100.0 * (simulation_capacity - experiment_capacity) / experiment_capacity
                ),
            }
        '''
    ),
    markdown(
        r"""
        ## [1] S1 — Ai2020 원형
        """
    ),
    code(
        r'''
        S1_PARAMETERS = pybamm.ParameterValues("Ai2020")
        S1_WINDOW = {"x0": 0.0065, "x100": 0.8400, "y0": 0.9651, "y100": 0.4350}
        S1_C_S_MAX_N = float(
            S1_PARAMETERS["Maximum concentration in negative electrode [mol.m-3]"]
        )
        S1_C_S_MAX_P = float(
            S1_PARAMETERS["Maximum concentration in positive electrode [mol.m-3]"]
        )

        print("S1 window:", S1_WINDOW)
        print(f"c_s,max,n = {S1_C_S_MAX_N:.1f} mol m^-3")
        print(f"c_s,max,p = {S1_C_S_MAX_P:.1f} mol m^-3")
        '''
    ),
    markdown(
        r"""
        ## [2] S2 — (D_e) 단위 보정 + 확정 stoichiometry window

        아래 window는 앞선 C/50 qOCV fitting에서 확정한 값으로, 이 노트북에서는 다시 fitting하지 않고 그대로 사용한다. (Q_n,Q_p)와 (c_{s,\max})는 코드에서 다시 계산한다.
        """
    ),
    code(
        r'''
        S2_WINDOW = {
            "x0": 0.0038062237,
            "x100": 0.9542212582,
            "y0": 0.9707457852,
            "y100": 0.4316424031,
        }

        S2_PARAMETERS = pybamm.ParameterValues("Ai2020")
        original_electrolyte_diffusivity = S2_PARAMETERS[
            "Electrolyte diffusivity [m2.s-1]"
        ]


        def corrected_electrolyte_diffusivity(concentration, temperature):
            return 1e-4 * original_electrolyte_diffusivity(concentration, temperature)


        negative_capacity_ah = CELL_CAPACITY_AH / (
            S2_WINDOW["x100"] - S2_WINDOW["x0"]
        )
        positive_capacity_ah = CELL_CAPACITY_AH / (
            S2_WINDOW["y0"] - S2_WINDOW["y100"]
        )

        overlap_area_m2 = 45.4e-3 * 51.0e-3 * 34.0
        negative_thickness_m = float(S2_PARAMETERS["Negative electrode thickness [m]"])
        positive_thickness_m = float(S2_PARAMETERS["Positive electrode thickness [m]"])
        negative_active_fraction = float(
            S2_PARAMETERS["Negative electrode active material volume fraction"]
        )
        positive_active_fraction = float(
            S2_PARAMETERS["Positive electrode active material volume fraction"]
        )

        S2_C_S_MAX_N = 3600.0 * negative_capacity_ah / (
            FARADAY_CONSTANT
            * overlap_area_m2
            * negative_thickness_m
            * negative_active_fraction
        )
        S2_C_S_MAX_P = 3600.0 * positive_capacity_ah / (
            FARADAY_CONSTANT
            * overlap_area_m2
            * positive_thickness_m
            * positive_active_fraction
        )

        S2_PARAMETERS.update(
            {
                "Electrolyte diffusivity [m2.s-1]": corrected_electrolyte_diffusivity,
                "Electrode height [m]": 51.0e-3,
                "Electrode width [m]": 45.4e-3,
                "Number of electrodes connected in parallel to make a cell": 34.0,
                "Maximum concentration in negative electrode [mol.m-3]": S2_C_S_MAX_N,
                "Maximum concentration in positive electrode [mol.m-3]": S2_C_S_MAX_P,
            }
        )

        print("S2 window:", S2_WINDOW)
        print(f"Qn = {negative_capacity_ah:.6f} Ah, Qp = {positive_capacity_ah:.6f} Ah")
        print(f"c_s,max,n = {S2_C_S_MAX_N:.1f} mol m^-3")
        print(f"c_s,max,p = {S2_C_S_MAX_P:.1f} mol m^-3")
        '''
    ),
    code(
        r'''
        # S1/S2의 0.5C, 1C, 2C 충·방전 실행 및 비교
        stage_inputs = {
            "S1 Ai2020 raw": (S1_PARAMETERS, S1_WINDOW, S1_C_S_MAX_N, S1_C_S_MAX_P),
            "S2 corrected De + window": (
                S2_PARAMETERS,
                S2_WINDOW,
                S2_C_S_MAX_N,
                S2_C_S_MAX_P,
            ),
        }

        simulations = {}
        result_rows = []

        for stage_name, (parameters, window, c_s_max_n, c_s_max_p) in stage_inputs.items():
            for c_rate in TEST_C_RATES:
                for is_charge in (True, False):
                    direction = "Charge" if is_charge else "Discharge"
                    print(f"{stage_name} | {c_rate:g}C {direction}")
                    simulation = run_dfn_branch(
                        parameters,
                        window,
                        c_s_max_n,
                        c_s_max_p,
                        c_rate,
                        is_charge,
                    )
                    simulations[(stage_name, c_rate, is_charge)] = simulation
                    metrics = compare_with_experiment(
                        experimental_curves[(c_rate, is_charge)],
                        simulation,
                        MEASURED_PROTOCOL[(c_rate, is_charge)]["current_a"],
                    )
                    result_rows.append(
                        {
                            "Stage": stage_name,
                            "C-rate": c_rate,
                            "Direction": direction,
                            **metrics,
                        }
                    )

        results = pd.DataFrame(result_rows)
        display(results.round(3))
        '''
    ),
    code(
        r'''
        # 전압-시간 곡선
        colors = {"S1 Ai2020 raw": "#888888", "S2 corrected De + window": "#0072B2"}
        line_styles = {"S1 Ai2020 raw": "--", "S2 corrected De + window": "-"}

        figure, axes = plt.subplots(2, 3, figsize=(16, 8.5), sharey=True, constrained_layout=True)
        for column, c_rate in enumerate(TEST_C_RATES):
            for row, is_charge in enumerate((True, False)):
                direction = "Charge" if is_charge else "Discharge"
                axis = axes[row, column]
                experiment = experimental_curves[(c_rate, is_charge)]
                experiment_time = experiment["time_min"] - experiment["time_min"][0]
                axis.plot(experiment_time, experiment["voltage_v"], "k", lw=2.7, label="Experiment")

                for stage_name in stage_inputs:
                    simulation = simulations[(stage_name, c_rate, is_charge)]
                    metric = results[
                        (results["Stage"] == stage_name)
                        & np.isclose(results["C-rate"], c_rate)
                        & (results["Direction"] == direction)
                    ].iloc[0]
                    axis.plot(
                        simulation["time_min"],
                        simulation["voltage_v"],
                        color=colors[stage_name],
                        ls=line_styles[stage_name],
                        lw=2.0,
                        label=f"{stage_name}: {metric.full_RMSE_mV:.1f} mV",
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
        plt.show()

        # 방향별 핵심 결과
        def rms(values):
            values = np.asarray(values, dtype=float)
            return float(np.sqrt(np.mean(values**2)))

        summary = (
            results.groupby(["Stage", "Direction"], as_index=False)
            .agg(
                Full_RMSE_mV=("full_RMSE_mV", rms),
                Center_10_70_RMSE_mV=("center_10_70_RMSE_mV", rms),
                Capacity_RMSE_pct=("capacity_error_pct", rms),
            )
        )
        display(summary.round(3))
        '''
    ),
]

for index, cell in enumerate(cells):
    cell["id"] = f"simple-s1-s2-{index:02d}"

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
