"""Append the local-minimum Rct calculation to the executed Stage-1 notebook."""

from pathlib import Path
import csv
import json


OUT = Path(__file__).resolve().parents[1] / "outputs"
notebook_path = next(OUT.glob("*Measured_Anode*ipynb"))
notebook = json.loads(notebook_path.read_text(encoding="utf-8"))
input_description = "".join(notebook["cells"][1]["source"])
drt_bullet = "- `Ai2020_Anode_DRT_Frequency_Origin_Paste.xlsx` (4절 N3/N4 적분)"
if drt_bullet not in input_description:
    input_description = input_description.replace(
        "- `260902 Enertech셀 데이터 모음.xlsx`",
        "- `260902 Enertech셀 데이터 모음.xlsx`\n" + drt_bullet,
    )
notebook["cells"][1]["source"] = input_description.splitlines(keepends=True)
notebook["cells"] = [
    cell for cell in notebook["cells"]
    if "## 4. 음극 DRT local-minimum" not in "".join(cell.get("source", []))
    and "# [4-1] N3+N4 local-minimum" not in "".join(cell.get("source", []))
]

markdown = """## 4. 음극 DRT local-minimum 적분으로 Rct·j0·k_n 초기값 계산

입력: `Ai2020_Anode_DRT_Frequency_Origin_Paste.xlsx`의 `Origin_Paste` 시트. γ는 10 mm 전극 면적을 반영한 Ω·cm² 단위이다. 각 SoC 곡선에서 N4/LF 골(0.15–0.45 Hz), N3/N4 골(3–8 Hz), N2/N3 골(50–120 Hz)을 Savitzky–Golay 평활화 후 찾는다. **저주파 바깥 골부터 고주파 바깥 골까지 원본 γ를 `d ln f`에 대해 적분**한다. 중간 골은 N3와 N4 기여도를 분리하는 데에만 사용한다. `f=1/(2πτ)`이므로 적분은 `d ln τ`와 같은 크기이다.

`R_ct,A = ∫γ d ln f`, `a_s=3ε_s/r_n`, `j₀=RT/(F R_ct,A a_s L_n)`, `k_n=j₀/[c_s,max,n √(c_e x_n(1−x_n))]`을 사용한다. `R_ct,A`는 계산 전에 Ω·cm²에서 Ω·m²로 바꾼다. 25 °C, `c_e=1000 mol m⁻³`, `ε_s,n=0.610`, `L_n=76.5 µm`, 최신 노트북 값 `r_n=3.692 µm`를 사용한다. **A2 슬라이드의 Ai2020 음극 + 측정 양극 OCP 조합**에는 `Q_cell=2.356528 Ah`, `x₀=0.003806`, `x₁₀₀=0.954221`을 대입하여 `c_s,max,n=25,182.8 mol m⁻³`를 사용한다. 현재 노트북의 같은 조합 재최적화 값은 `25,184.2 mol m⁻³`이며, 다른 **Measured anode + cathode OCP** 조합의 `32,447.1 mol m⁻³`와 구분한다. 표에 단일 값을 제시할 때는 **50% SoC**를 대표로 삼는다.

옛날 DRT 면적표에는 >10 kHz 접촉 피크가 따로 세어지지 않아, ~20 Hz와 ~1 Hz charge-transfer 피크가 **옛 N2/N3**이었다. 이후 SoC별 그림은 접촉 피크를 N1로 세어 같은 주파수대의 피크를 **현재 N3/N4**로 표기한다. 옛 N2+N3 면적은 `23.74489+6.95225=30.69714 Ω·cm²`이며, 옛 N3+N4 `12.41252 Ω·cm²`을 이력값으로 인용한 것은 번호 혼동이다. 아래의 새 Rct는 별도의 SoC별 DRT 원본을 local-minimum 방식으로 다시 적분한 값이다.

이 환산은 N3/N4를 모두 음극 charge transfer에 배정하는 가정이 있으므로, `k_n`은 DFN fitting **초기값**으로 사용한다. 원시 EIS/DRT와 전극 오염의 영향을 고려해 독립적인 확정 물성치로 취급하지 않는다.
"""

code = """# [4-1] N3+N4 local-minimum integration, A2 slide window
import openpyxl
from scipy.signal import savgol_filter

drt_candidates = [
    Path('Ai2020_Anode_DRT_Frequency_Origin_Paste.xlsx'),
    Path('outputs/Ai2020_Anode_DRT_Frequency_Origin_Paste.xlsx'),
    Path.home() / 'Downloads' / 'Ai2020_Anode_DRT_Frequency_Origin_Paste.xlsx',
]
drt_path = next((p for p in drt_candidates if p.exists()), None)
if drt_path is None:
    raise FileNotFoundError('Place Ai2020_Anode_DRT_Frequency_Origin_Paste.xlsx beside the notebook.')
drt_rows = list(openpyxl.load_workbook(drt_path, read_only=True, data_only=True)['Origin_Paste'].values)
slide_x0, slide_x100 = 0.003806, 0.954221
slide_cell_capacity_ah = 2.356528
slide_csmax_n = 3600.0 * slide_cell_capacity_ah / (
    (slide_x100 - slide_x0) * FARADAY_CONSTANT * OVERLAP_AREA_M2
    * 76.5e-6 * 0.610
)

def minimum_in_band(frequency, smoothed, lower_hz, upper_hz):
    indexes = np.flatnonzero((frequency >= lower_hz) & (frequency <= upper_hz))
    if len(indexes) < 3:
        raise ValueError(f'No local-minimum search points: {lower_hz}–{upper_hz} Hz')
    location = int(indexes[np.argmin(smoothed[indexes])])
    if location in (indexes[0], indexes[-1]):
        raise ValueError(f'Minimum on search-window edge: {frequency[location]:g} Hz')
    return location

drt_results = []
for column in range(0, 18, 2):
    soc_pct = int(str(drt_rows[0][column]).split('%')[0])
    pairs = np.array([
        (row[column], row[column + 1]) for row in drt_rows[3:]
        if isinstance(row[column], (int, float)) and isinstance(row[column + 1], (int, float))
    ], dtype=float)
    pairs = pairs[np.argsort(pairs[:, 0])]
    frequency, gamma = pairs.T
    smoothed = savgol_filter(gamma, 21, 3)
    low = minimum_in_band(frequency, smoothed, 0.15, 0.45)
    mid = minimum_in_band(frequency, smoothed, 3.0, 8.0)
    high = minimum_in_band(frequency, smoothed, 50.0, 120.0)
    r_n4 = np.trapezoid(gamma[low:mid + 1], np.log(frequency[low:mid + 1]))
    r_n3 = np.trapezoid(gamma[mid:high + 1], np.log(frequency[mid:high + 1]))
    r_ct = r_n3 + r_n4
    a_s_n = 3.0 * 0.610 / NEGATIVE_RADIUS_M
    j0_n = 8.314462618 * 298.15 / (
        FARADAY_CONSTANT * (r_ct * 1e-4) * a_s_n * 76.5e-6
    )
    x_n = slide_x0 + soc_pct / 100.0 * (slide_x100 - slide_x0)
    k_n = j0_n / (slide_csmax_n * np.sqrt(1000.0 * x_n * (1.0 - x_n)))
    measured_x_n = MEASURED_ANODE_WINDOW['x0'] + soc_pct / 100.0 * (
        MEASURED_ANODE_WINDOW['x100'] - MEASURED_ANODE_WINDOW['x0']
    )
    measured_k_n = j0_n / (
        measured_anode_csmax_n * np.sqrt(1000.0 * measured_x_n * (1.0 - measured_x_n))
    )
    drt_results.append({
        'SoC [%]': soc_pct,
        'f_low [Hz]': frequency[low], 'f_mid [Hz]': frequency[mid],
        'f_high [Hz]': frequency[high],
        'R_N4 [Ω cm²]': r_n4, 'R_N3 [Ω cm²]': r_n3,
        'Rct [Ω cm²]': r_ct, 'x_n': x_n,
        'j0 [A m⁻²]': j0_n, 'k_n A2 [A m⁻²(m³ mol⁻¹)^1.5]': k_n,
        'k_n measured anode [A m⁻²(m³ mol⁻¹)^1.5]': measured_k_n,
    })

drt_table = pd.DataFrame(drt_results).sort_values('SoC [%]', ascending=False)
display(drt_table.round({'f_low [Hz]': 4, 'f_mid [Hz]': 3, 'f_high [Hz]': 2,
                         'R_N4 [Ω cm²]': 4, 'R_N3 [Ω cm²]': 4, 'Rct [Ω cm²]': 4,
                         'x_n': 4, 'j0 [A m⁻²]': 4}))
representative = drt_table.loc[drt_table['SoC [%]'] == 50].iloc[0]
print(f"A2 c_s,max,n={slide_csmax_n:.1f} mol m⁻³")
print(f"50% SoC: Rct={representative['Rct [Ω cm²]']:.4f} Ω cm², "
      f"j0={representative['j0 [A m⁻²]']:.4f} A m⁻², "
      f"k_n(A2)={representative['k_n A2 [A m⁻²(m³ mol⁻¹)^1.5]']:.4e}")
"""

rows = list(csv.DictReader((OUT / "Anode_DRT_local_minimum_Rct_results.csv").open(encoding="utf-8-sig")))
lines = [" SoC  f_low   f_mid  f_high     R_N4     R_N3      Rct       j0           k_n\n"]
for row in rows:
    lines.append(
        f"{int(row['soc_pct']):4d}  {float(row['f_low_min_Hz']):6.4f}  "
        f"{float(row['f_mid_min_Hz']):6.3f}  {float(row['f_high_min_Hz']):6.2f}  "
        f"{float(row['R_N4_ohm_cm2']):8.4f}  {float(row['R_N3_ohm_cm2']):8.4f}  "
        f"{float(row['R_ct_N3_N4_ohm_cm2']):8.4f}  {float(row['j0_A_m2']):7.4f}  "
        f"{float(row['k_n_A_m-2_(m3_mol-1)^1.5']):.4e}\n"
    )
row50 = next(row for row in rows if int(row["soc_pct"]) == 50)
lines.append("A2 c_s,max,n=25182.8 mol m⁻³\n")
lines.append(
    f"50% SoC: Rct={float(row50['R_ct_N3_N4_ohm_cm2']):.4f} Ω cm², "
    f"j0={float(row50['j0_A_m2']):.4f} A m⁻², "
    f"k_n(A2)={float(row50['k_n_A_m-2_(m3_mol-1)^1.5']):.4e}\n"
)

notebook["cells"].extend([
    {
        "cell_type": "markdown", "metadata": {},
        "source": markdown.splitlines(keepends=True),
    },
    {
        "cell_type": "code", "execution_count": 12, "metadata": {},
        "source": code.splitlines(keepends=True),
        "outputs": [{"name": "stdout", "output_type": "stream", "text": lines}],
    },
])
notebook_path.write_text(json.dumps(notebook, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
print(notebook_path)
