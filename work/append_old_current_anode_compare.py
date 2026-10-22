"""Add old-vs-current anode potential comparison to the compact notebook."""
import base64
from pathlib import Path
import nbformat

OUT = Path(__file__).resolve().parents[1] / "outputs"
path = OUT / "261021_Ai2020_DFN_ORegan_양방향검증_간결본.ipynb"
nb = nbformat.read(path, as_version=4)
marker = "## [추가 비교] 10/20 공동 피팅 대비 음극 전위"
if not any(marker in c.source for c in nb.cells):
    nb.cells.append(nbformat.v4.new_markdown_cell(f"""{marker}

동일한 현재 Ai2020 음극 OCP/window, 측정 양극 OCP, ORegan2022 전해질 및 Cell 6-6 실험 조건에서 과거 10/20 공동 피팅값(`k_n=8.726541×10⁻⁷`, `k_p=5.899520×10⁻⁷`, `b_n=2.533583`, `r_n=3.7542 µm`, `r_p=4.2387 µm`)을 복원해 다시 계산했다. 이번 선정값은 `k_n=9.702793×10⁻⁷`, `k_p=4.768216×10⁻⁷`, `b_n=2.5`, `r_n=3.692 µm`, `r_p=4.239 µm`이다. 두 조합은 여러 파라미터가 함께 다르므로 아래 변화는 `b_n` 하나의 효과가 아니다.

| 충전율 | 과거 최저 | 현재 최저 | 현재−과거 |
| --- | ---: | ---: | ---: |
| 0.5C | +31.43 mV | +36.72 mV | **+5.30 mV** |
| 1C | −11.20 mV | −3.33 mV | **+7.87 mV** |
| 2C | −59.74 mV | −49.55 mV | **+10.19 mV** |

1C의 SoC 70/75/78%에서 과거는 +12.74/+1.65/−5.75 mV, 현재는 +19.54/+8.94/+1.87 mV이다. 0 mV 아래로 내려가는 시작 SoC가 과거 약 75.8%에서 현재 약 79.0%로 늦어진다. 2C의 시작 SoC는 약 29.7%→34.5%다. 0.5C는 두 조합 모두 양수다. 이는 DFN 예측의 상대 비교이며 도금 실측이 아니다.

원자료: `ORegan2022_과거공동피팅_vs_현재선정_음극전위_요약.csv`, `..._곡선.csv`.
"""))
    pic = (OUT / "ORegan2022_과거공동피팅_vs_현재선정_음극전위.png").read_bytes()
    cell = nbformat.v4.new_markdown_cell("![과거 대비 현재 음극 전위](attachment:old_current_anode_soc.png)")
    cell.attachments = {"old_current_anode_soc.png": {"image/png": base64.b64encode(pic).decode()}}
    nb.cells.append(cell)
    nbformat.validate(nb)
    nbformat.write(nb, path)
print(len(nb.cells))
