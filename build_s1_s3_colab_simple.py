"""Build the compact Colab notebook for the S1-S3 OCP construction flow."""

from __future__ import annotations

import json
from pathlib import Path
from textwrap import dedent


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "261019_Ai2020_S1_S3_현대차제출용_간결본.ipynb"


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
        # Ai2020 기반 DFN 모델: S1 → S2 → S3

        - **[0]** PyBaMM 실행 및 원본 실험 데이터 불러오기
        - **[1] S1** Ai2020 원형
        - **[2] S2** 전해질 확산계수 단위 보정 + Ai2020 OCP에서 stoichiometry window 선정
        - **[3] S3** 전극 OCP 비교 + 양극 GITT OCP 반영 + stoichiometry window 재선정

        아래 파일을 Colab의 같은 폴더에 업로드한 후 셀을 순서대로 실행한다.
        """
    ),
    code(
        r'''
        # [0-1] Colab 최초 실행 시 아래 줄의 주석을 해제한다.
        # PyBaMM 26.10.1부터 Ai2020 D_e가 내부에서 보정되므로,
        # S1(raw) → S2(D_e correction) 비교 재현을 위해 26.7.1.0을 고정한다.
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
        print("PyBaMM version:", pybamm.__version__)
        if pybamm.__version__ != "26.7.1.0":
            raise RuntimeError("PyBaMM 26.7.1.0으로 실행하세요. 설치 후 런타임을 재시작해야 할 수 있습니다.")

        # Colab에 업로드할 실제 원본 파일 이름
        FULL_CELL_LOW_RATE_FILES = [
            "260610 저율데이터 1.xlsx",
            "260610 저율데이터 2.xlsx",
        ]
        DYNAMIC_FILE = "260808 Eneterch셀 열화데이터 v1 6-6.xlsx"
        ANODE_LOW_RATE_FILES = [
            "260816 Anode GITT + OCP v1.xlsx",
            "260816 Anode GITT + OCP v2.xlsx",
        ]
        ANODE_GITT_FILES = [
            # v1은 v2/v3보다 branch 용량이 작아 중앙 형상 비교에서 제외한다.
            "260918 Anode GITT v2 1-6.xlsx",
            "260918 Anode GITT v3 1-7.xlsx",
        ]
        CATHODE_GITT_FILES = [
            "260816 Positive GITT 7-7.xlsx",
            "260816 Positive GITT 7-8.xlsx",
        ]

        required_files = (
            FULL_CELL_LOW_RATE_FILES
            + [DYNAMIC_FILE]
            + ANODE_LOW_RATE_FILES
            + ANODE_GITT_FILES
            + CATHODE_GITT_FILES
        )
        missing_files = [name for name in required_files if not Path(name).exists()]
        if missing_files:
            raise FileNotFoundError("먼저 Colab에 업로드하세요: " + ", ".join(missing_files))
        '''
    ),
    code(
        r'''
        # [0-2] C/50 qOCV와 Cell 6-6의 0.5C/1C/2C 원본 데이터 불러오기
        TEST_C_RATES = (0.5, 1.0, 2.0)
        SOC_GRID = np.linspace(0.0, 1.0, 1001)
        DYNAMIC_STEP_INDEX = {
            (0.5, True): 7,
            (0.5, False): 10,
            (1.0, True): 12,
            (1.0, False): 15,
            (2.0, True): 17,
            (2.0, False): 20,
        }


        def clean_curve(horizontal, vertical):
            horizontal = np.asarray(horizontal, dtype=float).reshape(-1)
            vertical = np.asarray(vertical, dtype=float).reshape(-1)
            valid = np.isfinite(horizontal) & np.isfinite(vertical)
            order = np.argsort(horizontal[valid], kind="stable")
            horizontal = horizontal[valid][order]
            vertical = vertical[valid][order]
            unique_horizontal, inverse = np.unique(horizontal, return_inverse=True)
            if len(unique_horizontal) == len(horizontal):
                return horizontal, vertical
            return unique_horizontal, np.bincount(inverse, weights=vertical) / np.bincount(inverse)


        def mean_overlap(curves, number_of_points=3001):
            lower = max(curve[0].min() for curve in curves)
            upper = min(curve[0].max() for curve in curves)
            grid = np.linspace(lower, upper, number_of_points)
            values = np.vstack([np.interp(grid, *curve) for curve in curves])
            return grid, values.mean(axis=0)


        def load_full_cell_low_rate(file_names):
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
                    capacity = pd.to_numeric(branch["Capacity(Ah)"], errors="coerce").to_numpy(float)
                    voltage = pd.to_numeric(branch["Voltage(V)"], errors="coerce").to_numpy(float)
                    current = float(pd.to_numeric(branch["Current(A)"], errors="coerce").median())
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
                        (end_capacity, np.interp(SOC_GRID, branch_soc, voltage))
                    )

                # 첫 방전은 부분 초기충전 직후이므로 제외한다.
                for capacity, voltage in branches["CC Chg"]:
                    charge_capacities.append(capacity)
                    charge_curves.append(voltage)
                for capacity, voltage in branches["CC DChg"][1:]:
                    discharge_capacities.append(capacity)
                    discharge_curves.append(voltage)

            mean_charge = np.mean(charge_curves, axis=0)
            mean_discharge = np.mean(discharge_curves, axis=0)
            return {
                "charge": mean_charge,
                "discharge": mean_discharge,
                "qocv": 0.5 * (mean_charge + mean_discharge),
                "capacity_ah": float(
                    0.5 * (np.mean(charge_capacities) + np.mean(discharge_capacities))
                ),
                "number_of_charge_curves": len(charge_curves),
                "number_of_discharge_curves": len(discharge_curves),
            }


        def load_dynamic_data(file_name):
            step = pd.read_excel(file_name, sheet_name="step")
            record = pd.read_excel(file_name, sheet_name="record")
            record["Date"] = pd.to_datetime(record["Date"])
            experiment = {}

            for key, step_index in DYNAMIC_STEP_INDEX.items():
                branch_row = step[
                    (step["Cycle Index"] == 1) & (step["Step Index"] == step_index)
                ].iloc[0]
                previous_row = step[
                    (step["Cycle Index"] == 1)
                    & (step["Step Number"] == int(branch_row["Step Number"]) - 1)
                ].iloc[0]

                branch = record[
                    (record["Date"] >= pd.to_datetime(branch_row["Oneset Date"]))
                    & (record["Date"] <= pd.to_datetime(branch_row["End Date"]))
                    & (record["Step Type"].astype(str) == str(branch_row["Step Type"]))
                ].copy()
                previous_rest = record[
                    (record["Date"] >= pd.to_datetime(previous_row["Oneset Date"]))
                    & (record["Date"] <= pd.to_datetime(previous_row["End Date"]))
                    & (record["Step Type"].astype(str) == "Rest")
                ].copy()

                time_min = (
                    branch["Date"] - branch["Date"].iloc[0]
                ).dt.total_seconds().to_numpy(float) / 60.0
                capacity = pd.to_numeric(branch["Capacity(Ah)"], errors="coerce").to_numpy(float)
                capacity = capacity - capacity[0]
                experiment[key] = {
                    "time_min": time_min,
                    "voltage_v": pd.to_numeric(branch["Voltage(V)"], errors="coerce").to_numpy(float),
                    "capacity_ah": capacity,
                    "current_a": abs(float(pd.to_numeric(branch["Current(A)"], errors="coerce").median())),
                    "rest_voltage_v": float(pd.to_numeric(previous_rest["Voltage(V)"], errors="coerce").iloc[-1]),
                }
            return experiment


        low_rate_data = load_full_cell_low_rate(FULL_CELL_LOW_RATE_FILES)
        CELL_CAPACITY_AH = low_rate_data["capacity_ah"]
        experimental_curves = load_dynamic_data(DYNAMIC_FILE)

        print(f"C/50 charge curves    : {low_rate_data['number_of_charge_curves']}")
        print(f"C/50 discharge curves : {low_rate_data['number_of_discharge_curves']}")
        print(f"Cell capacity          : {CELL_CAPACITY_AH:.6f} Ah")
        '''
    ),
    markdown(
        r"""
        ## [1] S1 — Ai2020 원형
        """
    ),
    code(
        r'''
        FARADAY_CONSTANT = 96485.33212
        S1_PARAMETERS = pybamm.ParameterValues("Ai2020")
        S1_WINDOW = {"x0": 0.0065, "x100": 0.8400, "y0": 0.9651, "y100": 0.4350}


        def ai2020_ocp_curves():
            grid = np.linspace(0.001, 0.999, 3001)

            def evaluate_ocp(parameter_name):
                function = S1_PARAMETERS[parameter_name]
                return np.array(
                    [
                        float(
                            np.asarray(
                                function(pybamm.Scalar(float(stoichiometry))).evaluate(),
                                dtype=float,
                            ).reshape(-1)[0]
                        )
                        for stoichiometry in grid
                    ],
                    dtype=float,
                )

            return {
                "negative": (grid, evaluate_ocp("Negative electrode OCP [V]")),
                "positive": (grid, evaluate_ocp("Positive electrode OCP [V]")),
            }


        AI2020_OCP = ai2020_ocp_curves()
        print("Published Ai2020 window:", S1_WINDOW)
        '''
    ),
    markdown(
        r"""
        ## [2] S2 — D_e 단위 보정 + Ai2020 OCP stoichiometry window 선정

        C/50 충·방전 평균으로 qOCV를 만들고, 2–98% SOC 전압 오차를 최소화한다. 양 끝점 오차는 각각 10 mV 이내로 제한한다.
        """
    ),
    code(
        r'''
        # [2-1] S2/S3에 사용할 전극 OCP 원본 준비
        def duration_to_seconds(value):
            if pd.isna(value):
                return np.nan
            if hasattr(value, "total_seconds"):
                return float(value.total_seconds())
            if isinstance(value, (int, float, np.integer, np.floating)):
                return float(value) * 86400.0
            parts = str(value).strip().split(":")
            if len(parts) != 3:
                return np.nan
            hours, minutes, seconds = map(float, parts)
            return 3600.0 * hours + 60.0 * minutes + seconds


        def duration_seconds(values):
            return values.map(duration_to_seconds)


        def extract_long_half_cell_branch(file_name, step_type):
            step = pd.read_excel(file_name, sheet_name="step")
            record = pd.read_excel(file_name, sheet_name="record")
            record["Date"] = pd.to_datetime(record["Date"])
            candidates = step[(step["Cycle Index"] == 2) & (step["Step Type"] == step_type)]
            selected = candidates.loc[candidates["Capacity(Ah)"].astype(float).idxmax()]
            branch = record[
                (record["Date"] >= pd.to_datetime(selected["Oneset Date"]))
                & (record["Date"] <= pd.to_datetime(selected["End Date"]))
                & (record["Step Type"].astype(str) == step_type)
            ]
            capacity = pd.to_numeric(branch["Capacity(Ah)"], errors="coerce").to_numpy(float)
            voltage = pd.to_numeric(branch["Voltage(V)"], errors="coerce").to_numpy(float)
            capacity = capacity - capacity[0]
            return capacity, voltage


        def build_anode_low_rate_candidate(file_names, parameters):
            disk_area_m2 = np.pi * (0.01 / 2.0) ** 2
            disk_capacity_ah = (
                FARADAY_CONSTANT
                * float(parameters["Maximum concentration in negative electrode [mol.m-3]"])
                * float(parameters["Negative electrode active material volume fraction"])
                * float(parameters["Negative electrode thickness [m]"])
                * disk_area_m2
                / 3600.0
            )
            curves = []
            for file_name in file_names:
                charge_q, charge_v = extract_long_half_cell_branch(file_name, "CC Chg")
                discharge_q, discharge_v = extract_long_half_cell_branch(file_name, "CC DChg")
                charge_span = charge_q.max() / disk_capacity_ah
                curves.append(clean_curve(charge_span - charge_q / disk_capacity_ah, charge_v))
                curves.append(clean_curve(discharge_q / disk_capacity_ah, discharge_v))
            return mean_overlap(curves)


        def build_anode_gitt_candidate(file_names):
            replicate_curves = []
            for file_name in file_names:
                step = pd.read_excel(file_name, sheet_name="step")
                step["duration_s"] = duration_seconds(step["Step Time"])
                cycle = step[step["Cycle Index"] == 2].reset_index(drop=True)
                first_rest = cycle[
                    (cycle["Step Type"] == "Rest") & cycle["duration_s"].between(7100, 7300)
                ].index[0]
                gitt = cycle.loc[first_rest:].reset_index(drop=True)

                cumulative = {"CC DChg": 0.0, "CC Chg": 0.0}
                points = {"CC DChg": [], "CC Chg": []}
                last_pulse = None
                for row_index, row in gitt.iterrows():
                    step_type = str(row["Step Type"])
                    if step_type in cumulative:
                        cumulative[step_type] += float(row["Capacity(Ah)"])
                        last_pulse = step_type
                    elif step_type == "Rest" and 7100 <= row["duration_s"] <= 7300:
                        if row_index == 0:
                            points["CC DChg"].append((0.0, float(row["End Voltage(V)"])))
                        elif last_pulse is not None:
                            points[last_pulse].append(
                                (cumulative[last_pulse], float(row["End Voltage(V)"]))
                            )
                            last_pulse = None

                lith = np.asarray(points["CC DChg"], dtype=float)
                delith = np.asarray(points["CC Chg"], dtype=float)
                delith = np.vstack([[0.0, lith[-1, 1]], delith])
                lith_curve = clean_curve(lith[:, 0] / lith[:, 0].max(), lith[:, 1])
                delith_curve = clean_curve(1.0 - delith[:, 0] / delith[:, 0].max(), delith[:, 1])
                replicate_curves.append(mean_overlap([lith_curve, delith_curve]))
            return mean_overlap(replicate_curves)


        def build_cathode_gitt_ocp(file_names, parameters):
            disk_area_m2 = np.pi * (0.01 / 2.0) ** 2
            disk_capacity_ah = (
                FARADAY_CONSTANT
                * float(parameters["Maximum concentration in positive electrode [mol.m-3]"])
                * float(parameters["Positive electrode active material volume fraction"])
                * float(parameters["Positive electrode thickness [m]"])
                * disk_area_m2
                / 3600.0
            )
            lithiation_curves, delithiation_curves = [], []

            for file_name in file_names:
                step = pd.read_excel(file_name, sheet_name="step")
                record = pd.read_excel(file_name, sheet_name="record")
                step["duration_s"] = duration_seconds(step["Step Time"])
                cycle = step[step["Cycle Index"] == 2].reset_index(drop=True)
                record["Date"] = pd.to_datetime(record["Date"])
                record = record.sort_values("Date").reset_index(drop=True)
                record_time_ns = record["Date"].to_numpy(dtype="datetime64[ns]").astype(np.int64)
                cumulative = {"CC Chg": 0.0, "CC DChg": 0.0}
                points = {"CC Chg": [], "CC DChg": []}

                for row_index, row in cycle.iterrows():
                    step_type = str(row["Step Type"])
                    if step_type not in cumulative or not (0 < row["duration_s"] <= 700):
                        continue
                    rest_index = row_index + 1
                    following_rests = []
                    while rest_index < len(cycle) and str(cycle.iloc[rest_index]["Step Type"]) == "Rest":
                        if cycle.iloc[rest_index]["duration_s"] >= 3500:
                            following_rests.append(cycle.iloc[rest_index])
                        rest_index += 1
                    if not following_rests:
                        continue
                    relaxed_state = following_rests[-1]

                    onset_ns = pd.Timestamp(row["Oneset Date"]).value
                    end_ns = pd.Timestamp(row["End Date"]).value
                    left = int(np.searchsorted(record_time_ns, onset_ns, side="left"))
                    right = int(np.searchsorted(record_time_ns, end_ns, side="right"))
                    pulse_record = record.iloc[left:right]
                    pulse_record = pulse_record[pulse_record["Step Type"].astype(str) == step_type]
                    pulse_capacity = pd.to_numeric(
                        pulse_record["Capacity(Ah)"], errors="coerce"
                    ).max()
                    cumulative[step_type] += float(pulse_capacity)
                    points[step_type].append(
                        (cumulative[step_type], float(relaxed_state["End Voltage(V)"]))
                    )

                charge = np.asarray(points["CC Chg"], dtype=float)
                discharge = np.asarray(points["CC DChg"], dtype=float)
                delithiation_curves.append(
                    clean_curve(1.0 - charge[:, 0] / disk_capacity_ah, charge[:, 1])
                )
                discharge_span = discharge[:, 0].max() / disk_capacity_ah
                lithiation_curves.append(
                    clean_curve((1.0 - discharge_span) + discharge[:, 0] / disk_capacity_ah, discharge[:, 1])
                )

            mean_lithiation = mean_overlap(lithiation_curves)
            mean_delithiation = mean_overlap(delithiation_curves)
            return mean_overlap([mean_lithiation, mean_delithiation])


        ANODE_LOW_RATE_OCP = build_anode_low_rate_candidate(ANODE_LOW_RATE_FILES, S1_PARAMETERS)
        ANODE_GITT_OCP = build_anode_gitt_candidate(ANODE_GITT_FILES)
        CATHODE_GITT_OCP = build_cathode_gitt_ocp(CATHODE_GITT_FILES, S1_PARAMETERS)
        '''
    ),
    markdown(
        r"""
        ## [3] S3 — 전극 OCP 선정 및 실험 양극 GITT 반영

        양극은 GITT 완화전압을 사용한다. 음극은 Ai2020, 저율 실험, 신규 GITT를 동일 조건에서 비교한 뒤
        qOCV, 동적 전압, 용량 및 absolute stoichiometry 축의 신뢰성을 함께 고려하여 기준모델을 선정한다.
        """
    ),
    code(
        r'''
        # [3-1] OCP 후보 비교: PPT의 전극 OCP 선정 페이지용
        figure, axes = plt.subplots(1, 3, figsize=(16, 4.6), constrained_layout=True)

        axes[0].plot(*AI2020_OCP["positive"], color="#777777", lw=2.0, label="Ai2020 LCO")
        axes[0].plot(*CATHODE_GITT_OCP, color="#D55E00", lw=2.3, label="Measured GITT mean")
        axes[0].set(title="Cathode OCP", xlabel="Cathode stoichiometry y", ylabel="OCP vs Li/Li+ [V]")

        axes[1].plot(*AI2020_OCP["negative"], color="#0072B2", lw=2.3, label="Ai2020 graphite (reference)")
        axes[1].plot(*ANODE_LOW_RATE_OCP, color="#009E73", lw=1.9, label="Measured low-rate OCP")
        axes[1].set(title="Anode OCP with absolute axis", xlabel="Anode stoichiometry x", ylabel="OCP vs Li/Li+ [V]")

        axes[2].plot(*ANODE_GITT_OCP, color="#E69F00", lw=2.2, label="Measured GITT v2/v3")
        axes[2].set(
            title="Anode GITT candidate",
            xlabel="Normalized GITT capacity (no absolute x anchor)",
            ylabel="Relaxed voltage vs Li/Li+ [V]",
        )

        for axis in axes:
            axis.grid(alpha=0.25)
            axis.legend(fontsize=8)
        plt.show()

        ocp_selection = pd.DataFrame(
            [
                ["Cathode", "Ai2020 LCO", "Reference", "Original S1/S2 baseline"],
                ["Cathode", "260816 Positive GITT 7-7, 7-8", "Selected", "Relaxed-voltage mean with capacity-based y axis"],
                ["Anode", "Ai2020 graphite", "Current reference", "Published absolute x axis; must be compared with measured candidates"],
                ["Anode", "260816 Anode GITT + OCP v1, v2", "Candidate", "Capacity-mapped x depends on Ai2020 c_s,max, epsilon_s and L"],
                ["Anode", "260918 Anode GITT v2 1-6, v3 1-7", "Candidate", "Shape comparison only; absolute x anchor is not yet established"],
            ],
            columns=["Electrode", "OCP source", "Use", "Reason"],
        )
        display(ocp_selection)
        '''
    ),
    code(
        r'''
        # [3-2] 각 OCP 조합에서 stoichiometry window를 독립적으로 선정한다.
        def predict_full_cell_ocv(window_vector, anode_curve, cathode_curve, soc=SOC_GRID):
            x0, x100, y100, y0 = map(float, window_vector)
            negative_stoichiometry = x0 + soc * (x100 - x0)
            positive_stoichiometry = y0 - soc * (y0 - y100)
            negative_ocp = np.interp(negative_stoichiometry, *anode_curve)
            positive_ocp = np.interp(positive_stoichiometry, *cathode_curve)
            return positive_ocp - negative_ocp


        def fit_stoichiometry_window(anode_curve, cathode_curve, random_seed):
            anode_min, anode_max = float(anode_curve[0].min()), float(anode_curve[0].max())
            cathode_min, cathode_max = float(cathode_curve[0].min()), float(cathode_curve[0].max())
            minimum_width = CELL_CAPACITY_AH / 6.0
            maximum_width = CELL_CAPACITY_AH / 2.4
            fitting_region = (SOC_GRID >= 0.02) & (SOC_GRID <= 0.98)

            def voltage_mse(window_vector):
                error = (
                    predict_full_cell_ocv(window_vector, anode_curve, cathode_curve)
                    - low_rate_data["qocv"]
                )
                return float(np.mean(error[fitting_region] ** 2))

            def constraint_violation(window_vector):
                x0, x100, y100, y0 = window_vector
                delta_x, delta_y = x100 - x0, y0 - y100
                endpoint_error = (
                    predict_full_cell_ocv(
                        window_vector, anode_curve, cathode_curve, SOC_GRID[[0, -1]]
                    )
                    - low_rate_data["qocv"][[0, -1]]
                )
                violations = [
                    max(0.0, minimum_width - delta_x),
                    max(0.0, delta_x - maximum_width),
                    max(0.0, minimum_width - delta_y),
                    max(0.0, delta_y - maximum_width),
                    max(0.0, abs(endpoint_error[0]) - 0.010),
                    max(0.0, abs(endpoint_error[1]) - 0.010),
                ]
                return float(np.sum(np.square(violations)))

            def penalized_objective(window_vector):
                return voltage_mse(window_vector) + 1e4 * constraint_violation(window_vector)

            bounds = [
                (anode_min, anode_max),
                (anode_min, anode_max),
                (cathode_min, cathode_max),
                (cathode_min, cathode_max),
            ]
            global_result = differential_evolution(
                penalized_objective,
                bounds=bounds,
                seed=random_seed,
                popsize=16,
                maxiter=350,
                tol=1e-10,
                polish=False,
                workers=1,
            )
            constraints = [
                {"type": "ineq", "fun": lambda z: (z[1] - z[0]) - minimum_width},
                {"type": "ineq", "fun": lambda z: maximum_width - (z[1] - z[0])},
                {"type": "ineq", "fun": lambda z: (z[3] - z[2]) - minimum_width},
                {"type": "ineq", "fun": lambda z: maximum_width - (z[3] - z[2])},
                {
                    "type": "ineq",
                    "fun": lambda z: 0.010
                    - abs(
                        predict_full_cell_ocv(z, anode_curve, cathode_curve, SOC_GRID[[0]])[0]
                        - low_rate_data["qocv"][0]
                    ),
                },
                {
                    "type": "ineq",
                    "fun": lambda z: 0.010
                    - abs(
                        predict_full_cell_ocv(z, anode_curve, cathode_curve, SOC_GRID[[-1]])[0]
                        - low_rate_data["qocv"][-1]
                    ),
                },
            ]
            local_result = minimize(
                voltage_mse,
                x0=global_result.x,
                method="SLSQP",
                bounds=bounds,
                constraints=constraints,
                options={"maxiter": 3000, "ftol": 1e-14, "disp": False},
            )
            if constraint_violation(local_result.x) <= constraint_violation(global_result.x) + 1e-12:
                return np.asarray(local_result.x, dtype=float)
            return np.asarray(global_result.x, dtype=float)


        s2_vector = fit_stoichiometry_window(
            AI2020_OCP["negative"], AI2020_OCP["positive"], random_seed=260924
        )
        s3_ai2020_anode_vector = fit_stoichiometry_window(
            AI2020_OCP["negative"], CATHODE_GITT_OCP, random_seed=260921
        )
        s3_low_rate_anode_vector = fit_stoichiometry_window(
            ANODE_LOW_RATE_OCP, CATHODE_GITT_OCP, random_seed=260922
        )
        s3_gitt_anode_vector = fit_stoichiometry_window(
            ANODE_GITT_OCP, CATHODE_GITT_OCP, random_seed=260923
        )

        def vector_to_window(vector):
            x0, x100, y100, y0 = map(float, vector)
            return {"x0": x0, "x100": x100, "y0": y0, "y100": y100}

        S2_WINDOW = vector_to_window(s2_vector)
        S3_AI2020_ANODE_WINDOW = vector_to_window(s3_ai2020_anode_vector)
        S3_LOW_RATE_ANODE_WINDOW = vector_to_window(s3_low_rate_anode_vector)
        S3_GITT_ANODE_WINDOW = vector_to_window(s3_gitt_anode_vector)
        '''
    ),
    code(
        r'''
        # [3-3] Ai2020 OCP 조합(S2)과 실험 양극 GITT 조합(S3)의 qOCV 비교
        window_rows = []
        qocv_predictions = {}
        for stage_name, vector, anode_curve, cathode_curve in [
            ("S2 Ai2020/Ai2020 OCP", s2_vector, AI2020_OCP["negative"], AI2020_OCP["positive"]),
            ("S3a Ai2020 anode", s3_ai2020_anode_vector, AI2020_OCP["negative"], CATHODE_GITT_OCP),
            ("S3b low-rate anode", s3_low_rate_anode_vector, ANODE_LOW_RATE_OCP, CATHODE_GITT_OCP),
            ("S3c new GITT anode", s3_gitt_anode_vector, ANODE_GITT_OCP, CATHODE_GITT_OCP),
        ]:
            prediction = predict_full_cell_ocv(vector, anode_curve, cathode_curve)
            qocv_predictions[stage_name] = prediction
            error_mv = 1000.0 * (prediction - low_rate_data["qocv"])
            fitting_region = (SOC_GRID >= 0.02) & (SOC_GRID <= 0.98)
            x0, x100, y100, y0 = vector
            window_rows.append(
                {
                    "Stage": stage_name,
                    "x0": x0,
                    "x100": x100,
                    "y0": y0,
                    "y100": y100,
                    "Qn_Ah": CELL_CAPACITY_AH / (x100 - x0),
                    "Qp_Ah": CELL_CAPACITY_AH / (y0 - y100),
                    "qOCV_RMSE_2_98_mV": np.sqrt(np.mean(error_mv[fitting_region] ** 2)),
                    "SOC0_error_mV": error_mv[0],
                    "SOC100_error_mV": error_mv[-1],
                }
            )

        window_results = pd.DataFrame(window_rows)
        display(window_results.round(6))

        figure, axes = plt.subplots(1, 2, figsize=(13, 4.7), constrained_layout=True)
        axes[0].plot(SOC_GRID * 100, low_rate_data["qocv"], color="black", lw=2.6, label="Measured C/50 qOCV")
        axes[1].axhline(0.0, color="black", lw=0.8)
        colors = {
            "S2 Ai2020/Ai2020 OCP": "#777777",
            "S3a Ai2020 anode": "#D55E00",
            "S3b low-rate anode": "#009E73",
            "S3c new GITT anode": "#CC79A7",
        }
        for stage_name, prediction in qocv_predictions.items():
            axes[0].plot(SOC_GRID * 100, prediction, color=colors[stage_name], lw=2.0, label=stage_name)
            axes[1].plot(
                SOC_GRID * 100,
                1000.0 * (prediction - low_rate_data["qocv"]),
                color=colors[stage_name],
                lw=1.8,
                label=stage_name,
            )
        axes[0].set(title="Full-cell qOCV reconstruction", xlabel="Full-cell SOC [%]", ylabel="Voltage [V]")
        axes[1].set(title="qOCV error", xlabel="Full-cell SOC [%]", ylabel="Model - experiment [mV]")
        for axis in axes:
            axis.axvspan(0, 2, color="#dddddd", alpha=0.45)
            axis.axvspan(98, 100, color="#dddddd", alpha=0.45)
            axis.grid(alpha=0.25)
            axis.legend(fontsize=8)
        plt.show()
        '''
    ),
    code(
        r'''
        # [3-4] window에서 Qn/Qp와 c_s,max를 일관되게 계산한다.
        def capacity_consistent_parameters(
            base_parameters, window, negative_ocp=None, positive_ocp=None
        ):
            parameters = base_parameters.copy()
            original_diffusivity = parameters["Electrolyte diffusivity [m2.s-1]"]

            def corrected_electrolyte_diffusivity(concentration, temperature):
                return 1e-4 * original_diffusivity(concentration, temperature)

            negative_capacity_ah = CELL_CAPACITY_AH / (window["x100"] - window["x0"])
            positive_capacity_ah = CELL_CAPACITY_AH / (window["y0"] - window["y100"])
            overlap_area_m2 = 45.4e-3 * 51.0e-3 * 34.0
            negative_c_s_max = 3600.0 * negative_capacity_ah / (
                FARADAY_CONSTANT
                * overlap_area_m2
                * float(parameters["Negative electrode thickness [m]"])
                * float(parameters["Negative electrode active material volume fraction"])
            )
            positive_c_s_max = 3600.0 * positive_capacity_ah / (
                FARADAY_CONSTANT
                * overlap_area_m2
                * float(parameters["Positive electrode thickness [m]"])
                * float(parameters["Positive electrode active material volume fraction"])
            )
            updates = {
                "Electrolyte diffusivity [m2.s-1]": corrected_electrolyte_diffusivity,
                "Electrode height [m]": 51.0e-3,
                "Electrode width [m]": 45.4e-3,
                "Number of electrodes connected in parallel to make a cell": 34.0,
                "Maximum concentration in negative electrode [mol.m-3]": negative_c_s_max,
                "Maximum concentration in positive electrode [mol.m-3]": positive_c_s_max,
            }
            if negative_ocp is not None:
                negative_grid, negative_voltage = negative_ocp

                def measured_negative_ocp(stoichiometry):
                    return pybamm.Interpolant(
                        negative_grid,
                        negative_voltage,
                        stoichiometry,
                        name="Measured graphite OCP",
                        interpolator="linear",
                    )

                updates["Negative electrode OCP [V]"] = measured_negative_ocp
            if positive_ocp is not None:
                positive_grid, positive_voltage = positive_ocp

                def measured_positive_ocp(stoichiometry):
                    return pybamm.Interpolant(
                        positive_grid,
                        positive_voltage,
                        stoichiometry,
                        name="Measured LCO GITT OCP",
                        interpolator="linear",
                    )

                updates["Positive electrode OCP [V]"] = measured_positive_ocp
            parameters.update(updates, check_already_exists=False)
            return parameters, negative_c_s_max, positive_c_s_max


        S2_PARAMETERS, S2_C_S_MAX_N, S2_C_S_MAX_P = capacity_consistent_parameters(
            pybamm.ParameterValues("Ai2020"), S2_WINDOW
        )
        S3_AI2020_ANODE_PARAMETERS, S3_AI2020_C_S_MAX_N, S3_AI2020_C_S_MAX_P = capacity_consistent_parameters(
            pybamm.ParameterValues("Ai2020"),
            S3_AI2020_ANODE_WINDOW,
            positive_ocp=CATHODE_GITT_OCP,
        )
        S3_LOW_RATE_ANODE_PARAMETERS, S3_LOW_RATE_C_S_MAX_N, S3_LOW_RATE_C_S_MAX_P = capacity_consistent_parameters(
            pybamm.ParameterValues("Ai2020"),
            S3_LOW_RATE_ANODE_WINDOW,
            negative_ocp=ANODE_LOW_RATE_OCP,
            positive_ocp=CATHODE_GITT_OCP,
        )
        S3_GITT_ANODE_PARAMETERS, S3_GITT_C_S_MAX_N, S3_GITT_C_S_MAX_P = capacity_consistent_parameters(
            pybamm.ParameterValues("Ai2020"),
            S3_GITT_ANODE_WINDOW,
            negative_ocp=ANODE_GITT_OCP,
            positive_ocp=CATHODE_GITT_OCP,
        )
        S1_C_S_MAX_N = float(S1_PARAMETERS["Maximum concentration in negative electrode [mol.m-3]"])
        S1_C_S_MAX_P = float(S1_PARAMETERS["Maximum concentration in positive electrode [mol.m-3]"])

        stage_summary = pd.DataFrame(
            [
                ["S1", "Ai2020 raw", "Ai2020 / Ai2020", "Published"],
                ["S2", "D_e correction + window", "Ai2020 / Ai2020", "C/50 fitted"],
                ["S3a", "Measured cathode + anode comparison", "Ai2020 / measured GITT", "C/50 re-fitted"],
                ["S3b", "Measured cathode + anode comparison", "Low-rate / measured GITT", "C/50 re-fitted"],
                ["S3c", "Measured cathode + anode comparison", "New GITT / measured GITT", "C/50 re-fitted"],
            ],
            columns=["Stage", "Change", "Anode / Cathode OCP", "Stoichiometry window"],
        )
        display(stage_summary)
        '''
    ),
    code(
        r'''
        # [결과-1] DFN 실행 및 실험 대비 오차 계산
        def stage_ocv(stage):
            window = stage["window"]
            return predict_full_cell_ocv(
                [window["x0"], window["x100"], window["y100"], window["y0"]],
                stage["anode_ocp"],
                stage["cathode_ocp"],
            )


        def initial_soc_from_rest_voltage(stage, rest_voltage_v):
            qocv = stage_ocv(stage)
            order = np.argsort(qocv, kind="stable")
            unique_voltage, unique_index = np.unique(qocv[order], return_index=True)
            return float(np.interp(rest_voltage_v, unique_voltage, SOC_GRID[order][unique_index]))


        def run_dfn_branch(stage, c_rate, is_charge):
            experiment = experimental_curves[(c_rate, is_charge)]
            initial_soc = initial_soc_from_rest_voltage(stage, experiment["rest_voltage_v"])
            window = stage["window"]
            initial_x = window["x0"] + initial_soc * (window["x100"] - window["x0"])
            initial_y = window["y0"] - initial_soc * (window["y0"] - window["y100"])
            branch_parameters = stage["parameters"].copy()
            branch_parameters.update(
                {
                    "Initial concentration in negative electrode [mol.m-3]": initial_x * stage["c_s_max_n"],
                    "Initial concentration in positive electrode [mol.m-3]": initial_y * stage["c_s_max_p"],
                }
            )

            direction = "Charge" if is_charge else "Discharge"
            cutoff_voltage = 4.2 if is_charge else 3.0
            protocol = f"{direction} at {experiment['current_a']} A until {cutoff_voltage} V"
            simulation = pybamm.Simulation(
                pybamm.lithium_ion.DFN({"thermal": "isothermal"}),
                parameter_values=branch_parameters,
                experiment=pybamm.Experiment([protocol], period="10 seconds"),
                var_pts={"x_n": 15, "x_s": 15, "x_p": 15, "r_n": 15, "r_p": 15},
                solver=pybamm.IDAKLUSolver(rtol=1e-6, atol=1e-8),
            )
            solution = simulation.solve()
            time_s = np.asarray(solution["Time [s]"].entries, dtype=float)
            voltage_v = np.asarray(solution["Terminal voltage [V]"].entries, dtype=float)
            current_a = np.asarray(solution["Current [A]"].entries, dtype=float)
            anode_potential_v = np.asarray(
                solution["Negative electrode surface potential difference at separator interface [V]"].entries,
                dtype=float,
            ).squeeze()
            capacity_ah = cumulative_trapezoid(np.abs(current_a), time_s, initial=0.0) / 3600.0
            if is_charge:
                full_cell_soc_pct = 100.0 * (initial_soc + capacity_ah / CELL_CAPACITY_AH)
            else:
                full_cell_soc_pct = 100.0 * (initial_soc - capacity_ah / CELL_CAPACITY_AH)
            return {
                "time_min": time_s / 60.0,
                "voltage_v": voltage_v,
                "anode_potential_v": anode_potential_v,
                "capacity_ah": capacity_ah,
                "full_cell_soc_pct": np.clip(full_cell_soc_pct, 0.0, 100.0),
                "initial_soc_pct": 100.0 * initial_soc,
            }


        def compare_with_experiment(experiment, simulation):
            experiment_time = experiment["time_min"] - experiment["time_min"][0]
            simulation_time = simulation["time_min"] - simulation["time_min"][0]
            common_end = 0.995 * min(experiment_time[-1], simulation_time[-1])
            full_grid = np.linspace(0.0, common_end, 400)
            center_start = 0.10 * CELL_CAPACITY_AH / experiment["current_a"] * 60.0
            center_end = min(0.70 * CELL_CAPACITY_AH / experiment["current_a"] * 60.0, common_end)
            if center_end <= center_start:
                center_grid = full_grid[full_grid >= 0.10 * common_end]
            else:
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
                "full_MAE_mV": float(np.mean(np.abs(full_error))),
                "center_10_70_RMSE_mV": float(np.sqrt(np.mean(center_error**2))),
                "center_10_70_MAE_mV": float(np.mean(np.abs(center_error))),
                "capacity_error_pct": float(100.0 * (simulation_capacity - experiment_capacity) / experiment_capacity),
            }
        '''
    ),
    code(
        r'''
        # [결과-2] S1/S2/S3의 0.5C, 1C, 2C 충·방전 실행
        stages = {
            "S1 Ai2020 raw": {
                "parameters": S1_PARAMETERS,
                "window": S1_WINDOW,
                "c_s_max_n": S1_C_S_MAX_N,
                "c_s_max_p": S1_C_S_MAX_P,
                "anode_ocp": AI2020_OCP["negative"],
                "cathode_ocp": AI2020_OCP["positive"],
            },
            "S2 corrected D_e + Ai2020 window": {
                "parameters": S2_PARAMETERS,
                "window": S2_WINDOW,
                "c_s_max_n": S2_C_S_MAX_N,
                "c_s_max_p": S2_C_S_MAX_P,
                "anode_ocp": AI2020_OCP["negative"],
                "cathode_ocp": AI2020_OCP["positive"],
            },
            "S3a Ai2020 anode": {
                "parameters": S3_AI2020_ANODE_PARAMETERS,
                "window": S3_AI2020_ANODE_WINDOW,
                "c_s_max_n": S3_AI2020_C_S_MAX_N,
                "c_s_max_p": S3_AI2020_C_S_MAX_P,
                "anode_ocp": AI2020_OCP["negative"],
                "cathode_ocp": CATHODE_GITT_OCP,
            },
            "S3b low-rate anode": {
                "parameters": S3_LOW_RATE_ANODE_PARAMETERS,
                "window": S3_LOW_RATE_ANODE_WINDOW,
                "c_s_max_n": S3_LOW_RATE_C_S_MAX_N,
                "c_s_max_p": S3_LOW_RATE_C_S_MAX_P,
                "anode_ocp": ANODE_LOW_RATE_OCP,
                "cathode_ocp": CATHODE_GITT_OCP,
            },
            "S3c new GITT anode": {
                "parameters": S3_GITT_ANODE_PARAMETERS,
                "window": S3_GITT_ANODE_WINDOW,
                "c_s_max_n": S3_GITT_C_S_MAX_N,
                "c_s_max_p": S3_GITT_C_S_MAX_P,
                "anode_ocp": ANODE_GITT_OCP,
                "cathode_ocp": CATHODE_GITT_OCP,
            },
        }

        simulations = {}
        result_rows = []
        for stage_name, stage in stages.items():
            for c_rate in TEST_C_RATES:
                for is_charge in (True, False):
                    direction = "Charge" if is_charge else "Discharge"
                    print(f"{stage_name} | {c_rate:g}C {direction}")
                    simulation = run_dfn_branch(stage, c_rate, is_charge)
                    simulations[(stage_name, c_rate, is_charge)] = simulation
                    result_rows.append(
                        {
                            "Stage": stage_name,
                            "C-rate": c_rate,
                            "Direction": direction,
                            "minimum_anode_potential_mV": (
                                1000.0 * float(np.min(simulation["anode_potential_v"]))
                                if is_charge
                                else np.nan
                            ),
                            **compare_with_experiment(experimental_curves[(c_rate, is_charge)], simulation),
                        }
                    )

        results = pd.DataFrame(result_rows)
        display(results.round(3))
        '''
    ),
    code(
        r'''
        # [결과-3] 단계별 핵심 결과
        def root_mean_square(values):
            values = np.asarray(values, dtype=float)
            return float(np.sqrt(np.mean(values**2)))


        summary = (
            results.groupby(["Stage", "Direction"], as_index=False)
            .agg(
                Full_RMSE_mV=("full_RMSE_mV", root_mean_square),
                Full_MAE_mV=("full_MAE_mV", "mean"),
                Center_10_70_RMSE_mV=("center_10_70_RMSE_mV", root_mean_square),
                Center_10_70_MAE_mV=("center_10_70_MAE_mV", "mean"),
                Capacity_RMSE_pct=("capacity_error_pct", root_mean_square),
            )
        )
        display(summary.round(3))

        s3_comparison = (
            results[results["Stage"].str.startswith("S3")]
            .groupby("Stage", as_index=False)
            .agg(
                Dynamic_full_RMSE_mV=("full_RMSE_mV", root_mean_square),
                Dynamic_center_RMSE_mV=("center_10_70_RMSE_mV", root_mean_square),
                Capacity_RMSE_pct=("capacity_error_pct", root_mean_square),
                Maximum_capacity_error_pct=("capacity_error_pct", lambda x: float(np.max(np.abs(x)))),
                Minimum_charge_anode_potential_mV=("minimum_anode_potential_mV", "min"),
            )
        )
        s3_qocv = window_results[window_results["Stage"].str.startswith("S3")][
            ["Stage", "qOCV_RMSE_2_98_mV", "Qn_Ah", "Qp_Ah"]
        ]
        s3_comparison = s3_comparison.merge(s3_qocv, on="Stage", how="left")
        s3_comparison["Physical interpretation"] = s3_comparison["Stage"].map(
            {
                "S3a Ai2020 anode": "공개된 absolute x 축을 가진 현재 기준 후보",
                "S3b low-rate anode": "조건부 실험 후보: x 환산에 Ai2020 c_s,max, epsilon_s, L 사용",
                "S3c new GITT anode": "민감도 후보: 독립적인 absolute x anchor 미확보",
            }
        )
        s3_comparison["Current decision"] = s3_comparison["Stage"].map(
            {
                "S3a Ai2020 anode": "기준모델 유지",
                "S3b low-rate anode": "보류",
                "S3c new GITT anode": "absolute x 확보 후 재검토",
            }
        )
        display(s3_comparison.round(3))

        print(
            "현재 판단: Ai2020 음극 OCP는 측정 후보를 배제해서가 아니라, "
            "세 후보 중 absolute x 축과 qOCV·용량·동적 전압의 종합 균형을 함께 만족하므로 유지한다."
        )
        print(
            "단, 신규 GITT에 독립적인 x anchor를 부여하고 3전극 음극 전위로 검증하면 "
            "실험 OCP로 교체할 수 있다."
        )
        '''
    ),
    code(
        r'''
        # [결과-4] PPT p2용: Ai2020 원형의 실험 대비 전압-시간 곡선
        def plot_voltage_time(stage_names, title):
            colors = {
                "S1 Ai2020 raw": "#888888",
                "S2 corrected D_e + Ai2020 window": "#0072B2",
                "S3a Ai2020 anode": "#D55E00",
                "S3b low-rate anode": "#009E73",
                "S3c new GITT anode": "#CC79A7",
            }
            line_styles = {
                "S1 Ai2020 raw": "--",
                "S2 corrected D_e + Ai2020 window": "-",
                "S3a Ai2020 anode": "-",
                "S3b low-rate anode": "--",
                "S3c new GITT anode": ":",
            }
            figure, axes = plt.subplots(2, 3, figsize=(16, 8.5), sharey=True, constrained_layout=True)
            for column, c_rate in enumerate(TEST_C_RATES):
                for row, is_charge in enumerate((True, False)):
                    direction = "Charge" if is_charge else "Discharge"
                    axis = axes[row, column]
                    experiment = experimental_curves[(c_rate, is_charge)]
                    axis.plot(experiment["time_min"], experiment["voltage_v"], "k", lw=2.7, label="Experiment")
                    for stage_name in stage_names:
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
            figure.legend(handles, labels, loc="outside upper center", ncol=len(handles), frameon=False)
            figure.suptitle(title, weight="bold", fontsize=15)
            plt.show()


        plot_voltage_time(["S1 Ai2020 raw"], "S1: original Ai2020 model")
        '''
    ),
    code(
        r'''
        # [결과-5] PPT p3용: S2와 음극 OCP 후보 3개를 비교
        plot_voltage_time(
            [
                "S2 corrected D_e + Ai2020 window",
                "S3a Ai2020 anode",
                "S3b low-rate anode",
                "S3c new GITT anode",
            ],
            "S2 vs S3: cathode GITT with three anode OCP candidates",
        )
        '''
    ),
    code(
        r'''
        # [결과-6] 음극 OCP 후보별 충전 중 음극 전위 예측
        # 3전극 실험으로 검증하기 전까지는 절대값 확정치가 아니라 후보 간 경향 비교로 해석한다.
        anode_candidates = [
            "S3a Ai2020 anode",
            "S3b low-rate anode",
            "S3c new GITT anode",
        ]
        colors = {
            "S3a Ai2020 anode": "#D55E00",
            "S3b low-rate anode": "#009E73",
            "S3c new GITT anode": "#CC79A7",
        }
        labels = {
            "S3a Ai2020 anode": "Ai2020 anode",
            "S3b low-rate anode": "Measured low-rate anode",
            "S3c new GITT anode": "Measured GITT anode",
        }

        figure, axes = plt.subplots(1, 3, figsize=(15, 4.6), sharey=True, constrained_layout=True)
        for axis, c_rate in zip(axes, TEST_C_RATES):
            for stage_name in anode_candidates:
                simulation = simulations[(stage_name, c_rate, True)]
                axis.plot(
                    simulation["full_cell_soc_pct"],
                    1000.0 * simulation["anode_potential_v"],
                    lw=2.2,
                    color=colors[stage_name],
                    label=labels[stage_name],
                )
            axis.axhline(0.0, color="black", ls="--", lw=1.2, label="0 mV vs Li/Li+" if c_rate == TEST_C_RATES[0] else None)
            axis.set_title(f"{c_rate:g}C Charge", weight="bold")
            axis.set_xlabel("Full-cell SOC [%]")
            axis.grid(alpha=0.25)
        axes[0].set_ylabel("Predicted anode potential [mV vs Li/Li+]")
        handles, legend_labels = axes[0].get_legend_handles_labels()
        figure.legend(handles, legend_labels, loc="outside upper center", ncol=4, frameon=False)
        figure.suptitle("S3 anode OCP candidates: DFN-predicted anode potential (unvalidated)", weight="bold")
        plt.show()
        '''
    ),
]

for index, cell in enumerate(cells):
    cell["id"] = f"simple-s1-s3-{index:02d}"

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
