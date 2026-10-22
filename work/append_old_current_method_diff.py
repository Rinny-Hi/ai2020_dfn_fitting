"""Append method and parameter differences to compact notebook."""
from pathlib import Path
import nbformat

OUT = Path(__file__).resolve().parents[1] / "outputs"
path = OUT / "261021_Ai2020_DFN_ORegan_양방향검증_간결본.ipynb"
nb = nbformat.read(path, as_version=4)
marker = "### 과거·현재의 피팅 방식과 파라미터 차이"
if not any(marker in cell.source for cell in nb.cells):
    nb.cells.append(nbformat.v4.new_markdown_cell(f"""{marker}

| | 과거 10/20 | 현재 선정 |
| --- | --- | --- |
| 피팅/검증 | 충·방전 6곡선 공동 피팅, 독립 방전 검증 없음 | 충전 3곡선 피팅, 방전 3곡선 보류 검증 |
| 피팅 변수 | `k_n` + `b_n` | `b_n=2.5` 고정, `k_n`만 |
| 최적화 | PSO 7입자×7반복 | 1차원 bounded scalar |
| 전압 목적함수 | 120점/곡선, 실험시간 99.5%, 조기 컷오프 잔차 페널티 500 mV 비율 | 140점/곡선, 실험시간 98%, 조기 컷오프 별도 20항/곡선·300 mV 비율 |

| 파라미터 | 과거 | 현재 |
| --- | ---: | ---: |
| `k_n` | `8.726541×10⁻⁷` | `9.702793×10⁻⁷` |
| `b_n` | 2.533583 | 2.500000 |
| `k_p` | `5.899520×10⁻⁷` | `4.768216×10⁻⁷` |
| `r_n` | 3.7542 µm | 3.6920 µm |
| `c_s,max,n` | 25,182.81 | 25,184.17 mol m⁻³ |

현재 모델에서 값을 한 번에 하나씩 바꾼 통제 계산의 **동일 SoC 50%** 전위 차이: `k_n` 변경은 1C/2C에서 +4.23/+5.19 mV, `b_n` 변경은 +1.51/+2.83 mV, `r_n` 변경은 +0.70/+0.88 mV이다. `k_p` 변경은 같은 SoC 전위에 거의 영향이 없으나 종료 SoC를 바꿔 최저값에 영향을 준다. 이는 국소적인 하나씩 변경 계산이며 효과는 단순 합산할 수 없다. 상세: `ORegan2022_과거현재_방식_파라미터_차이.md`, 원자료 `ORegan2022_과거현재_음극전위_파라미터기여.csv`.
"""))
    nbformat.validate(nb)
    nbformat.write(nb, path)
print(len(nb.cells))
