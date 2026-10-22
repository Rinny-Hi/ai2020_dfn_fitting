"""Append the selected-model audit to the user-facing compact notebook."""
import base64
from pathlib import Path

import nbformat

OUT = Path(__file__).resolve().parents[1] / "outputs"
PATH = OUT / "261021_Ai2020_DFN_ORegan_양방향검증_간결본.ipynb"
nb = nbformat.read(PATH, as_version=4)
marker = "## [최종 판정] Ai2020 OCP·ORegan2022의 양방향성과 SoC별 음극 전위"
if not any(marker in cell.source for cell in nb.cells):
    nb.cells.append(nbformat.v4.new_markdown_cell(
        f"""{marker}

**권고:** Ai2020 음극 OCP + 측정 GITT 양극 OCP + ORegan2022 전해질에서 음극 Bruggeman 계수 `b_n=2.5`를 고정하고, `k_n=9.7028×10⁻⁷`을 0.5/1/2C 충전으로 피팅한다. 방전은 독립 검증이다. 전체 여섯 조건 RMSE는 37.14 mV, 방전 검증 31.74 mV, 2C 방전 21.41 mV다.

OCP별 stoichiometry window와 `c_s,max,n/p`를 다시 적용했고, EIS `j₀`에서 환산한 초기 `k_n/k_p`에도 각 window의 50% SoC 조성을 반영했다. Ai2020 분기의 `c_s,max,n/p=25,184.17/49,143.14 mol m⁻³`, 초기 EIS `k_n/k_p=6.3076/4.7682×10⁻⁷`. 피팅 후 `k_n`은 전압으로 재추정한 값이며 `k_p`는 EIS 환산값으로 유지한다.

두 변수를 0.5/1C 충·방전으로 공동 피팅하면 `k_n=1.599×10⁻⁶, b_n=2.983`이나, 보류한 2C 방전 RMSE가 112.66 mV로 커지고 1C 음극 전위가 SoC 약 69%부터 음수가 된다. 이 결과는 채택하지 않았다. `b_n=2.5` 선택은 음극 전위를 0 mV 위로 강제로 끌어올린 피팅이 아니다.

권고 조합의 1C 음극 분리막 계면 전위는 SoC 50/70/75/78%에서 +44.74/+19.54/+8.94/+1.87 mV이고, 약 79–80%에서 최저 −3.33 mV이다. 격자 15→25점 검증에서 이 최저값은 −3.334→−3.332 mV였다. **1C에서 0 mV 이상이라는 엄격한 조건은 만족하지 않지만**, SoC에 따라 매끄럽게 내려가고 작은 음수는 마지막 약 1% SoC에 국한된다. 2C 최저는 −49.55 mV로 고율 충전 용도로 별도 평가가 필요하다. 이 전위는 DFN 예측이며 도금 실측값이 아니다.

아래 표는 데이터 파일에서 다시 계산할 수 있다. 자세한 근거: `ORegan2022_Ai2020_최종선정_및_SoC음극전위.md`.
"""
    ))
    nb.cells.append(nbformat.v4.new_code_cell(
        """from pathlib import Path
import pandas as pd

out = Path.cwd()  # 이 노트북을 outputs 폴더에서 열었을 때
cross = pd.read_csv(out / 'ORegan2022_교차피팅_조건별결과.csv')
cross = cross[(cross.OCP_case == 'Ai2020_anode__measured_cathode') &
              (cross.fit_type == 'kn_only_bn2.5') &
              (cross.fit_direction == 'charge')]
display(cross[['C_rate', 'Direction', 'Data_role', 'full_RMSE_mV',
               'full_MAE_mV', 'capacity_error_pct', 'Anode_min_full_mV',
               'Anode_min_SOC10_70_mV']].round(3))

joint = pd.read_csv(out / 'ORegan2022_Ai2020_1C이하_공동피팅_조건별결과.csv')
display(joint[['fit_type', 'C_rate', 'Direction', 'Data_role',
               'full_RMSE_mV', 'capacity_error_pct', 'Anode_min_full_mV']].round(3))
"""
    ))
    raw = (OUT / "ORegan2022_Ai2020_SoC별_음극전위_선정비교.png").read_bytes()
    cell = nbformat.v4.new_markdown_cell("![SoC별 음극 전위](attachment:anode_soc_audit.png)")
    cell.attachments = {"anode_soc_audit.png": {"image/png": base64.b64encode(raw).decode()}}
    nb.cells.append(cell)
    nbformat.validate(nb)
    nbformat.write(nb, PATH)
print(len(nb.cells), PATH)
