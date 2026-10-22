"""Append latest k_p correction, 3-stage audit, figures and method review."""
import base64
from pathlib import Path

import nbformat

OUT = Path(__file__).resolve().parents[1] / "outputs"
PATH = OUT / "261021_Ai2020_DFN_ORegan_양방향검증_간결본.ipynb"
nb = nbformat.read(PATH, as_version=4)
marker = "## [최신 정리] k_p 정정, 충전 피팅·방전 검증, 목적함수 점검"
if not any(marker in cell.source for cell in nb.cells):
    nb.cells.append(nbformat.v4.new_markdown_cell(f"""{marker}

**현재 모든 재계산의 `k_p=4.768215585×10⁻⁷`**이다. 이전 비교에 등장한 `5.899520×10⁻⁷`은 10/20 과거 공동 피팅 파라미터를 복원할 때만 썼다. 최종 조합은 Ai2020 음극 OCP·측정 GITT 양극 OCP·ORegan2022 전해질, `b_n=2.5`, 충전 0.5/1/2C로 피팅한 `k_n=9.702793×10⁻⁷`이다. 방전 세 율은 독립 검증이다.

| 조건 | 역할 | 전체 RMSE/MAE (mV) | 실제 SoC 10–70% RMSE/MAE (mV) | 용량 오차 |
| --- | --- | ---: | ---: | ---: |
| 0.5C 충전 | 피팅 | 35.94 / 28.97 | 37.60 / 35.16 | −1.84% |
| 0.5C 방전 | 검증 | 35.60 / 27.26 | 30.01 / 26.28 | +0.46% |
| 1C 충전 | 피팅 | 42.18 / 28.57 | 26.21 / 24.00 | −1.78% |
| 1C 방전 | 검증 | 36.01 / 27.68 | 29.52 / 25.07 | +0.91% |
| 2C 충전 | 피팅 | 46.74 / 23.13 | 10.69 / 9.67 | −3.40% |
| 2C 방전 | 검증 | 21.41 / 11.30 | 8.31 / 4.95 | +0.03% |

실험 초기값 A(`k_n=6.3076×10⁻⁷,b_n=2.914`)→`b_n`만 고정한 B(`k_n=6.3076×10⁻⁷,b_n=2.5`)→최종 충전 피팅 C의 여섯 조건 전구간 RMSE는 **66.69→39.70→37.14 mV**이다. C는 B보다 2C/종료 용량·1C 음극 전위가 좋지만, **0.5C 두 방향과 일부 중앙 SoC 오차는 나쁘다.** 중앙 구간의 가중치를 정하기 전까지 절대 최적 물성이라고 해석하지 않는다. 2C 충전의 중앙 계산은 실제 도달 가능한 약 SoC 10.8–59.2%다.

목적함수 조사: 현재 1변수 bounded scalar 탐색 자체는 안정적이지만, 조기 컷오프 이후 전압을 평탄하게 외삽하고 임의 크기의 시간 부족 벌점을 넣는 방법은 개선 여지가 있다. 충전/방전 분할을 유지하고 공통 관측구간 전압·종료 용량·실제 SoC 중앙 구간을 분리 평가한 뒤, 반복 측정 불확도에 근거해 비용 가중치를 정하는 방안을 권고한다. PyBOP의 데이터·시뮬레이터·비용함수·최적화 분리 구조를 참고할 수 있다. 실제 PyBOP 포팅과 새 비용함수 재피팅은 수행하지 않았다.

상세 방법 및 연구 코드 근거: `ORegan2022_Ai2020_최종결과_방법검토.md`. 표 원자료: `ORegan2022_Ai2020_최종_3단계_조건별_RMSE_MAE.csv`.
"""))
    images = [
        ("피팅 전후 충방전 6곡선", "ORegan2022_Ai2020_피팅전후_충방전6곡선.png", "before_after_six.png"),
        ("최종 충방전 6곡선", "ORegan2022_Ai2020_최종_충방전6곡선.png", "final_six.png"),
        ("최종 SoC별 음극 전위", "ORegan2022_Ai2020_최종_SoC별음극전위.png", "final_anode_soc.png"),
    ]
    for title, name, attachment in images:
        cell = nbformat.v4.new_markdown_cell(f"### {title}\n\n![{title}](attachment:{attachment})")
        cell.attachments = {attachment: {"image/png": base64.b64encode((OUT / name).read_bytes()).decode()}}
        nb.cells.append(cell)
    nb.cells.append(nbformat.v4.new_code_cell(
        """from pathlib import Path
import pandas as pd

here = Path.cwd()
out = here if (here / 'ORegan2022_Ai2020_최종_3단계_조건별_RMSE_MAE.csv').exists() else here / 'outputs'
metrics = pd.read_csv(out / 'ORegan2022_Ai2020_최종_3단계_조건별_RMSE_MAE.csv')
display(metrics[['Stage', 'C_rate', 'Direction', 'Data_role', 'Full_RMSE_mV', 'Full_MAE_mV',
                 'SOC10_70_RMSE_mV', 'SOC10_70_MAE_mV', 'SOC_window_low_pct',
                 'SOC_window_high_pct', 'Capacity_error_pct', 'Anode_min_mV']].round(3))
soc = pd.read_csv(out / 'ORegan2022_Ai2020_최종_SoC별음극전위.csv')
display(soc.round(3))
"""
    ))
    nbformat.validate(nb)
    nbformat.write(nb, PATH)
print(len(nb.cells))
