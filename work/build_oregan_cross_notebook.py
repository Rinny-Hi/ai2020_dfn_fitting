"""Append audited ORegan cross-direction fits to the existing reproducible notebook."""
from __future__ import annotations

import base64
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
SOURCE = OUT / "261021_Ai2020_DFN_최종_간결본.ipynb"
TARGET = OUT / "261021_Ai2020_DFN_ORegan_양방향검증_간결본.ipynb"
nb = json.loads(SOURCE.read_text(encoding="utf-8"))


def markdown(s):
    return {"cell_type": "markdown", "metadata": {}, "source": s.splitlines(True)}


def code(s, outputs=None):
    return {"cell_type": "code", "metadata": {}, "execution_count": None,
            "outputs": outputs or [], "source": s.splitlines(True)}


def table(frame):
    lines = ["| " + " | ".join(str(c) for c in frame.columns) + " |",
             "| " + " | ".join("---" for _ in frame.columns) + " |"]
    for row in frame.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(str(v) for v in row) + " |")
    return "\n".join(lines)


def embedded_image(path, caption):
    content = base64.b64encode(path.read_bytes()).decode("ascii")
    src = f'from IPython.display import Image, display\ndisplay(Image(filename="{path.name}"))\n'
    output = {"output_type": "display_data", "data": {
        "image/png": content, "text/plain": [f"<Figure: {caption}>"]}, "metadata": {}}
    return code(src, [output])


summary = pd.read_csv(OUT / "ORegan2022_교차피팅_방향별요약.csv")
summary["OCP"] = summary.OCP_case.map(lambda s: "Ai2020 음극" if s.startswith("Ai2020") else "측정 음극")
summary["설정"] = summary.apply(lambda x:
    "자유 b / 충전 피팅" if x.fit_type == "kn+bn" and x.fit_direction == "charge" else
    "자유 b / 방전 피팅" if x.fit_type == "kn+bn" else
    "b=2.5 / 충전 피팅" if x.fit_direction == "charge" else "b=2.5 / 방전 피팅", axis=1)
summary["k_n"] = summary.kn.map(lambda x: f"{x:.3g}")
summary["b_n"] = summary.brugg_n.map(lambda x: f"{x:.3f}")
for col in ("Charge_RMSE_3rates_mV", "Discharge_RMSE_3rates_mV",
            "Charge_2C_full_RMSE_mV", "Discharge_2C_full_RMSE_mV",
            "Six_branch_RMSE_mV", "Charge_2C_min_anode_mV"):
    summary[col] = summary[col].map(lambda x: f"{x:.1f}")
summary = summary[["OCP", "설정", "k_n", "b_n", "Charge_RMSE_3rates_mV",
                   "Discharge_RMSE_3rates_mV", "Charge_2C_full_RMSE_mV",
                   "Discharge_2C_full_RMSE_mV", "Charge_2C_min_anode_mV"]]
summary.columns = ["OCP", "설정", "k_n", "b_n", "충전 RMSE", "방전 RMSE",
                   "2C 충전", "2C 방전", "2C 최저 음극전위"]

rows = pd.read_csv(OUT / "ORegan2022_교차피팅_조건별결과.csv")
rows = rows.loc[(rows.fit_type == "kn_only_bn2.5") & (rows.fit_direction == "charge")].copy()
rows["OCP"] = rows.OCP_case.map(lambda s: "Ai2020" if s.startswith("Ai2020") else "Measured")
rows["조건"] = rows.C_rate.map(lambda x: f"{x:g}C") + " " + rows.Direction.map({"Charge":"충전", "Discharge":"방전"})
rows["전체 RMSE/MAE"] = rows.apply(lambda x: f"{x.full_RMSE_mV:.1f}/{x.full_MAE_mV:.1f}", axis=1)
rows["SoC10–70 RMSE/MAE"] = rows.apply(lambda x: f"{x.True_SOC10_70_RMSE_mV:.1f}/{x.True_SOC10_70_MAE_mV:.1f}", axis=1)
rows["용량 오차"] = rows.capacity_error_pct.map(lambda x: f"{x:+.2f}%")
rows = rows[["OCP", "조건", "전체 RMSE/MAE", "SoC10–70 RMSE/MAE", "용량 오차"]]

cells = nb["cells"]
cells.append(markdown("""
# 업데이트: ORegan2022 충·방전 교차 검증과 `b_n=2.5` 고정

앞 절의 자유 `b_n` **충전 피팅→방전 검증**에 이어, 여기서는 **방전 피팅→충전 검증**을 추가한다. 그 후 `b_n=2.5`를 고정하고 `k_n`만 충전 기준과 방전 기준으로 각각 피팅해 반대 방향을 검증한다. 두 OCP, 전해질, 실험 `D_s` 및 양극 `k_p`는 앞 절과 같다. 아래 결과는 저장된 실행값이며 뒤의 코드 셀로 재계산할 수 있다.

`b_n=2.5`는 문헌에서 보고된 음극 Bruggeman 계수 범위에 속하는 **모델 비교값**이다. 현재 Enertech 전극에서 직접 측정한 값은 아니다. 기존 방전 결과를 보고 이 고정값을 채택했으므로, 반대 방향 검증은 `k_n` 피팅에서 제외한 방향의 검증이지 완전한 블라인드 모델 선정 실험은 아니다.

목적함수: 대상 방향의 0.5C/1C/2C 실험 곡선에서 각각 140개 시간점(실험 종료시간의 98%까지)의 전압 잔차(mV)와 조기 cutoff 페널티를 결합한 RMS. `k_n∈[3×10⁻⁷,3×10⁻⁶]`; 자유 피팅의 `b_n∈[1.5,3.5]`. 반대 방향은 피팅 목적함수에 넣지 않았다. 전 구간 오차는 실험·모델이 공통으로 존재하는 시간 구간에서 계산하고, 종료 용량 오차는 따로 기록한다.
"""))
cells.append(markdown("## 방향별 결과 (RMSE 및 전위 단위 mV)\n\n" + table(summary) + "\n"))
cells.append(markdown("""
**판정:** Ai2020 음극 OCP에서 `b_n=2.5`·충전 `k_n` 피팅이 충전 41.9 mV, 반대 방향 방전 31.7 mV의 세 율 통합 RMSE와 2C 방전 21.4 mV를 보였다. 측정 음극은 2C 방전 17.8 mV지만 0.5C/1C 방전 오차가 커 전체 방전 37.2 mV다. 두 OCP 모두 2C 충전 최저 음극 전위는 음수이므로 전압 정합을 도금 안전성 확보로 해석하지 않는다. 기존 Ai2020 보정 전해질·실험 초기값의 여섯 조건 RMSE(37.11 mV)와 Ai2020 OCP의 ORegan 고정 `b_n` 결과(37.14 mV)는 사실상 같다. 따라서 전해질 교체가 **단자전압 전체 오차를 개선했다는 주장도 하지 않는다.**
"""))
cells.append(markdown("## `b_n=2.5` 충전 피팅의 조건별 결과 (mV)\n\n" + table(rows) + "\n"))

# The earlier notebook cells define cases, kinetic_constants, model_stage,
# potential_metrics, and the experimental data. The appended code reuses them.
fit_src = (ROOT / "work" / "oregan_cross_direction_fit.py").read_text(encoding="utf-8")
fit_src = fit_src.replace("from __future__ import annotations\n", "")
fit_src = fit_src.replace("from run_stage2_two_ocp import cases, kinetic_constants, scope\n", "")
fit_src = fit_src.replace("from run_stage3_charge_fit import model_stage, potential_metrics\n", "")
fit_src = fit_src.replace('OUT = Path(__file__).resolve().parents[1] / "outputs"',
                          'scope = globals()\nOUT = Path.cwd()')
fit_src = fit_src.replace('if __name__ == "__main__":\n    main()\n', '')
cells.append(markdown("## 재실행 코드\n\n앞 절의 계산 셀을 실행한 뒤 아래 셀을 실행한다. 여섯 피팅과 36개 조건 재검증이므로 시간이 걸리며, 탐색 기록과 결과 CSV가 노트북 작업 폴더에 저장된다.\n"))
cells.append(code(fit_src))
cells.append(code("main()\n"))

plot_src = (ROOT / "work" / "plot_oregan_cross_direction.py").read_text(encoding="utf-8")
plot_src = plot_src.replace('OUT = Path(__file__).resolve().parents[1] / "outputs"',
                            'OUT = Path.cwd()')
cells.append(code(plot_src))
cells.append(markdown("## 저장된 전압 곡선과 방향별 절충\n\n검정은 실험, 색상은 파라미터 설정이다. 빨간 점선은 앞 절의 자유 `b_n` 충전 피팅이며, 청색은 `b_n=2.5` 충전 피팅이다.\n"))
for filename, caption in (
    ("ORegan2022_교차피팅_6곡선_Ai2020_OCP.png", "Ai2020 음극 OCP 여섯 곡선"),
    ("ORegan2022_교차피팅_6곡선_Measured_OCP.png", "측정 음극 OCP 여섯 곡선"),
    ("ORegan2022_교차피팅_충방전_절충.png", "충·방전 교차 검증 절충"),
):
    cells.append(embedded_image(OUT / filename, caption))

mesh = pd.read_csv(OUT / "ORegan2022_bn2.5_mesh15_vs25.csv")
mesh["OCP"] = mesh.OCP_case.map(lambda s: "Ai2020" if s.startswith("Ai2020") else "Measured")
mesh["조건"] = mesh.C_rate.map(lambda x: f"{x:g}C") + " " + mesh.Direction
mesh["전압 격자차 RMS"] = mesh.Voltage_RMS_mesh_difference_mV.map(lambda x: f"{x:.2f}")
mesh["25점 RMSE"] = mesh.Mesh25_full_RMSE_mV.map(lambda x: f"{x:.2f}")
mesh["25점 최저 음극전위"] = mesh.Mesh25_min_anode_mV.map(lambda x: f"{x:.2f}")
mesh = mesh[["OCP", "조건", "전압 격자차 RMS", "25점 RMSE", "25점 최저 음극전위"]]
cells.append(markdown("""
## 격자 정련 및 해석 한계

15→25점 정련 시 핵심 2C 전압 차이는 1 mV 미만이다. 측정 음극 OCP의 1C 최저 전위는 +0.41→−1.07 mV로 부호가 바뀌고, 고율 충전에서 측정 OCP 보간 상한 외삽 경고가 발생한다. 따라서 이 경우의 1C 전위가 안전하다고 주장하지 않는다. 2C에서는 두 OCP 모두 음수다. 음극 전위는 3전극 실측이 아닌 DFN 예측이다.

참고: [음극 Bruggeman 계수 범위를 보고한 재파라미터화 연구](https://doi.org/10.1016/j.etran.2023.100231), [Chen et al.의 전극별 Bruggeman 측정 사례](https://wrap.warwick.ac.uk/id/eprint/136795/2/WRAP-development-experimental-techniques-parameterization-multi-scale-lithium-ion-battery-models-Chen-2020.pdf).
""" + "\n" + table(mesh) + "\n"))

for index, cell in enumerate(cells):
    cell.setdefault("id", hashlib.sha1(f"{index}:".encode() +
                                   "".join(cell["source"]).encode("utf-8")).hexdigest()[:12])
TARGET.write_text(json.dumps(nb, ensure_ascii=False, indent=1), encoding="utf-8")
print(TARGET, TARGET.stat().st_size, len(cells))
