import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
SOURCE = Path(r"C:\Users\user\Downloads\261020_Ai2020_1단계_OCP_Stoichiometry.ipynb")
TARGET = ROOT / "outputs" / "261020_Ai2020_1단계_OCP_Stoichiometry_Measured_Anode_비교.ipynb"

notebook = json.loads(SOURCE.read_text(encoding="utf-8"))
cells = notebook["cells"]


def source(index):
    return "".join(cells[index]["source"])


def set_source(index, value):
    cells[index]["source"] = value.strip("\n").splitlines(keepends=True)


def replace_once(index, old, new):
    value = source(index)
    assert value.count(old) == 1, (index, old[:80], value.count(old))
    set_source(index, value.replace(old, new))


set_source(0, """
# Ai2020 기반 DFN 모델 — 1단계 OCP 비교

**목적:** 공통 물리 보정 후 Ai2020 음극, 측정 저율 음극, 측정 GITT 양극 OCP와 각 stoichiometry window의 영향을 같은 데이터로 비교한다.

1. Ai2020 음극·양극 OCP와 공개 window
2. Ai2020 음극·양극 OCP와 재최적화 window
3. Ai2020 음극 + 측정 양극 GITT OCP와 재최적화 window
4. **측정 저율 음극 + 같은 측정 양극 GITT OCP와 재최적화 window**

3과 4는 양극 OCP, geometry, 수송·반응 파라미터, 전류, 평가 구간을 공유한다. 음극 OCP와 그에 따른 window·inventory 차이를 함께 평가한다. `D_s`, `k`, Bruggeman 및 전해질 parameter set의 fitting은 범위 밖이다.
""")
set_source(1, """
## 0. 실행 환경과 입력 파일

아래 6개 파일을 노트북과 같은 폴더에 둔다. 마지막 파일은 [저장소의 통합 Excel](https://github.com/Rinny-Hi/ai2020_dfn_fitting/blob/main/260902%20Enertech%EC%85%80%20%EB%8D%B0%EC%9D%B4%ED%84%B0%20%EB%AA%A8%EC%9D%8C.xlsx)이며 `ANODE_LOWRATE`(v1/v2 충·방전)와 `PARAMETER_PRIOR` 시트를 사용한다.

- `260610 저율데이터 1.xlsx`
- `260610 저율데이터 2.xlsx`
- `260808 Eneterch셀 열화데이터 v1 6-6.xlsx`
- `260816 Positive GITT 7-7.xlsx`
- `260816 Positive GITT 7-8.xlsx`
- `260902 Enertech셀 데이터 모음.xlsx`

측정 음극은 **저율 half-cell OCP**를 사용한다. 신규 음극 GITT OCP는 절대 stoichiometry 기준이 확정되지 않아 이 비교에 포함하지 않는다.
""")
replace_once(2, 'CATHODE_GITT_FILES = [\n    "260816 Positive GITT 7-7.xlsx",\n    "260816 Positive GITT 7-8.xlsx",\n]\nrequired_files = FULL_CELL_LOW_RATE_FILES + [DYNAMIC_FILE] + CATHODE_GITT_FILES', 'CATHODE_GITT_FILES = [\n    "260816 Positive GITT 7-7.xlsx",\n    "260816 Positive GITT 7-8.xlsx",\n]\nMEASURED_ANODE_BOOK = "260902 Enertech셀 데이터 모음.xlsx"\nrequired_files = FULL_CELL_LOW_RATE_FILES + [DYNAMIC_FILE] + CATHODE_GITT_FILES + [MEASURED_ANODE_BOOK]')
set_source(4, """
## 1. 공통 기반 고정

모든 OCP 후보에 동일하게 적용한다.

- Ai2020 electrolyte diffusivity: 원 구현의 단위 해석에 따라 `× 10⁻⁴`
- 음극 평균 입자반경: **3.6920 μm**
- 양극 평균 입자반경: **4.2390 μm**
- 전극 겹침 면적: `45.4 mm × 51.0 mm × 34면`

위 입자반경은 **첨부 노트북의 실행 코드값**이다. 브랜치 간결본의 3.7542/4.2387 μm와 다르므로, 브랜치에 저장된 수치를 이번 실행 결과로 인용하지 않는다.
""")
set_source(6, """
## 2. OCP와 stoichiometry window

양극은 두 GITT 반복시험의 장시간 완화전압을 평균해 사용한다. 측정 음극은 통합 Excel의 저율 half-cell v1/v2 충·방전을 같은 용량 기준으로 변환·평균한다. 이때 절대 x축은 Ai2020 `c_s,max`, 활물질 분율·두께 및 10 mm disk 면적에 의존하므로 **독립적인 실측 절대 stoichiometry 기준이 아니다**.

각 OCP 조합의 window는 같은 C/50 충·방전 평균 qOCV의 2–98% SOC RMSE를 최소화한다. SOC 0/100% 끝점 오차는 각각 10 mV 이내로 제한한다. 측정 OCP의 보간 범위와 선택 window의 관계를 함께 표시한다.

4.2 V는 full-cell 단자전압이므로 양극 OCP 자체를 `4.2 V vs Li/Li+`에서 절단하지 않는다.
""")
anode_code = '''

def build_measured_anode_lowrate_ocp(book_name, parameters):
    frame = pd.read_excel(book_name, sheet_name="ANODE_LOWRATE", header=8)
    prior = pd.read_excel(book_name, sheet_name="PARAMETER_PRIOR", header=7)
    diameter_row = prior.loc[prior["parameter"].astype(str).eq("disk_diameter_m"), "value"]
    if diameter_row.empty:
        raise ValueError("PARAMETER_PRIOR에 disk_diameter_m이 없습니다")
    disk_diameter_m = float(diameter_row.iloc[0])
    disk_area_m2 = np.pi * (disk_diameter_m / 2.0) ** 2
    disk_capacity_ah = (
        FARADAY_CONSTANT
        * float(parameters["Maximum concentration in negative electrode [mol.m-3]"])
        * float(parameters["Negative electrode active material volume fraction"])
        * float(parameters["Negative electrode thickness [m]"])
        * disk_area_m2 / 3600.0
    )
    if disk_capacity_ah <= 0:
        raise ValueError("음극 disk 기준 용량이 양수가 아닙니다")

    lithiation_curves, delithiation_curves = [], []
    for version in ("v1", "v2"):
        branches = {}
        for direction in ("charge", "discharge"):
            selected = frame.loc[
                frame["source_version"].astype(str).eq(version)
                & frame["direction"].astype(str).eq(direction)
            ]
            capacity = pd.to_numeric(selected["Q_Ah"], errors="coerce").to_numpy(float)
            voltage = pd.to_numeric(selected["V_V"], errors="coerce").to_numpy(float)
            finite = np.isfinite(capacity) & np.isfinite(voltage)
            capacity, voltage = capacity[finite], voltage[finite]
            if capacity.size < 2:
                raise ValueError(f"ANODE_LOWRATE {version}/{direction} 데이터가 부족합니다")
            branches[direction] = (capacity, voltage)
        q_discharge, v_discharge = branches["discharge"]
        q_charge, v_charge = branches["charge"]
        lithiation_curves.append(clean_curve(q_discharge / disk_capacity_ah, v_discharge))
        delithiation_curves.append(
            clean_curve(q_charge.max() / disk_capacity_ah - q_charge / disk_capacity_ah, v_charge)
        )

    mean_lithiation = mean_overlap(lithiation_curves)
    mean_delithiation = mean_overlap(delithiation_curves)
    return mean_overlap([mean_lithiation, mean_delithiation]), disk_capacity_ah


MEASURED_ANODE_OCP, ANODE_DISK_CAPACITY_AH = build_measured_anode_lowrate_ocp(
    MEASURED_ANODE_BOOK, raw_parameters
)
print(
    f"Measured anode x support: {MEASURED_ANODE_OCP[0].min():.6f}–"
    f"{MEASURED_ANODE_OCP[0].max():.6f}; assumed disk Q={ANODE_DISK_CAPACITY_AH:.6f} Ah"
)
'''
set_source(8, source(8) + anode_code)
replace_once(9, 'AI_REFIT_WINDOW = vector_to_window(ai_refit_vector)\nIMPROVED_WINDOW = vector_to_window(improved_vector)', 'measured_anode_vector = fit_stoichiometry_window(\n    MEASURED_ANODE_OCP, CATHODE_GITT_OCP, random_seed=261020\n)\nAI_REFIT_WINDOW = vector_to_window(ai_refit_vector)\nIMPROVED_WINDOW = vector_to_window(improved_vector)\nMEASURED_ANODE_WINDOW = vector_to_window(measured_anode_vector)')
replace_once(10, '    ("Measured cathode OCP + refit window", improved_vector, AI2020_OCP["negative"], CATHODE_GITT_OCP, "#D55E00", "-"),', '    ("Measured cathode OCP + refit window", improved_vector, AI2020_OCP["negative"], CATHODE_GITT_OCP, "#D55E00", "-"),\n    ("Measured anode + cathode OCP + refit window", measured_anode_vector, MEASURED_ANODE_OCP, CATHODE_GITT_OCP, "#009E73", "-"),')
replace_once(10, 'window_results = pd.DataFrame(qocv_rows)\ndisplay(window_results.round(6))', 'window_results = pd.DataFrame(qocv_rows)\nwindow_results["Anode support x_min"] = [np.nan, np.nan, np.nan, MEASURED_ANODE_OCP[0].min()]\nwindow_results["Anode support x_max"] = [np.nan, np.nan, np.nan, MEASURED_ANODE_OCP[0].max()]\ndisplay(window_results.round(6))')
set_source(11, """
## 3. Capacity-consistent inventory와 DFN 검증

각 재최적화 window에서 `Q_n`, `Q_p`를 계산하고, 동일한 면적·두께·활물질 분율을 쓰도록 `c_s,max,n/p`를 역산한다. 공통 OCP 쌍에 대해 **Ai2020 음극 → 측정 음극** 변화가 qOCV와 0.5C/1C/2C 충·방전 전압·cut-off 용량에 미치는 영향을 비교한다. `c_s,max`는 window에 따른 bookkeeping 값이며 독립 측정값이 아니다.
""")
set_source(12, '''
def build_candidate_parameters(base_parameters, window, positive_ocp, negative_ocp=None):
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

    updates = {
        "Maximum concentration in negative electrode [mol.m-3]": csmax_n,
        "Maximum concentration in positive electrode [mol.m-3]": csmax_p,
        "Positive electrode OCP [V]": measured_positive_ocp,
    }
    if negative_ocp is not None:
        negative_grid, negative_voltage = negative_ocp

        def measured_negative_ocp(stoichiometry):
            return pybamm.Interpolant(
                negative_grid, negative_voltage, stoichiometry,
                name="Measured anode low-rate OCP", interpolator="linear"
            )

        updates["Negative electrode OCP [V]"] = measured_negative_ocp
    parameters.update(updates)
    return parameters, csmax_n, csmax_p


improved_parameters, improved_csmax_n, improved_csmax_p = build_candidate_parameters(
    common_parameters, IMPROVED_WINDOW, CATHODE_GITT_OCP
)
measured_anode_parameters, measured_anode_csmax_n, measured_anode_csmax_p = build_candidate_parameters(
    common_parameters, MEASURED_ANODE_WINDOW, CATHODE_GITT_OCP, MEASURED_ANODE_OCP
)
baseline_csmax_n = float(common_parameters["Maximum concentration in negative electrode [mol.m-3]"])
baseline_csmax_p = float(common_parameters["Maximum concentration in positive electrode [mol.m-3]"])

parameter_summary = pd.DataFrame(
    [
        ["Negative OCP", "Ai2020", "Ai2020", "Measured low-rate v1/v2"],
        ["Positive OCP", "Ai2020", "Measured GITT mean", "Measured GITT mean"],
        ["x0 / x100", f"{PUBLISHED_WINDOW['x0']:.6f} / {PUBLISHED_WINDOW['x100']:.6f}", f"{IMPROVED_WINDOW['x0']:.6f} / {IMPROVED_WINDOW['x100']:.6f}", f"{MEASURED_ANODE_WINDOW['x0']:.6f} / {MEASURED_ANODE_WINDOW['x100']:.6f}"],
        ["y0 / y100", f"{PUBLISHED_WINDOW['y0']:.6f} / {PUBLISHED_WINDOW['y100']:.6f}", f"{IMPROVED_WINDOW['y0']:.6f} / {IMPROVED_WINDOW['y100']:.6f}" , f"{MEASURED_ANODE_WINDOW['y0']:.6f} / {MEASURED_ANODE_WINDOW['y100']:.6f}"],
        ["c_s,max,n [mol m-3]", f"{baseline_csmax_n:.1f}", f"{improved_csmax_n:.1f}", f"{measured_anode_csmax_n:.1f}"],
        ["c_s,max,p [mol m-3]", f"{baseline_csmax_p:.1f}", f"{improved_csmax_p:.1f}", f"{measured_anode_csmax_p:.1f}"],
    ],
    columns=["Item", "Baseline", "Measured cathode", "Measured anode + cathode"],
)
display(parameter_summary)
''')
replace_once(14, '    "Improved: measured cathode OCP/refit window": {\n        "parameters": improved_parameters,\n        "window": IMPROVED_WINDOW,\n        "c_s_max_n": improved_csmax_n,\n        "c_s_max_p": improved_csmax_p,\n        "anode_ocp": AI2020_OCP["negative"],\n        "cathode_ocp": CATHODE_GITT_OCP,\n    },', '    "Improved: measured cathode OCP/refit window": {\n        "parameters": improved_parameters,\n        "window": IMPROVED_WINDOW,\n        "c_s_max_n": improved_csmax_n,\n        "c_s_max_p": improved_csmax_p,\n        "anode_ocp": AI2020_OCP["negative"],\n        "cathode_ocp": CATHODE_GITT_OCP,\n    },\n    "Measured: anode + cathode OCP/refit window": {\n        "parameters": measured_anode_parameters,\n        "window": MEASURED_ANODE_WINDOW,\n        "c_s_max_n": measured_anode_csmax_n,\n        "c_s_max_p": measured_anode_csmax_p,\n        "anode_ocp": MEASURED_ANODE_OCP,\n        "cathode_ocp": CATHODE_GITT_OCP,\n    },')
replace_once(15, '    "Improved: measured cathode OCP/refit window": "#D55E00",', '    "Improved: measured cathode OCP/refit window": "#D55E00",\n    "Measured: anode + cathode OCP/refit window": "#009E73",')
replace_once(15, 'figure.suptitle("Stage 1: baseline vs OCP/window improvement", weight="bold")', 'figure.suptitle("Stage 1: measured anode OCP comparison", weight="bold")')
set_source(16, """
## 해석 원칙

- **측정 양극만** 적용한 경우와 **측정 음극·양극 모두** 적용한 경우를 비교하면, 같은 양극 아래 음극 OCP와 재최적화 window의 영향을 볼 수 있다. 이는 음극 OCP만의 고정-window 효과와 동일하지 않다.
- qOCV 표에는 window, `Q_n/Q_p`, 측정 음극 OCP 지원 x범위를 표시한다. 고율 DFN에서 OCP 보간범위 외삽 경고가 나오면 그 branch를 별도로 해석한다.
- 전압 RMSE는 실험·모델 공통 시간 구간만 계산한다. 조기 cutoff가 있다면 그 뒤의 누락 구간이 RMSE에 포함되지 않으므로, 해당 RMSE만으로 전체 구간 달성을 주장하지 않는다.
- `cut-off capacity error`는 4.2/3.0 V 도달 용량의 상대오차다. `D_s`, `k`, Bruggeman은 재피팅하지 않았다.
- 측정 저율 음극의 절대 x축은 Ai2020 inventory 가정에 의존한다. 음극 전위 예측은 기준전극 실측 검증이 아니다.
""")

for cell in cells:
    if cell["cell_type"] == "code":
        cell["outputs"] = []
        cell["execution_count"] = None

TARGET.parent.mkdir(parents=True, exist_ok=True)
TARGET.write_text(json.dumps(notebook, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
print(TARGET)
