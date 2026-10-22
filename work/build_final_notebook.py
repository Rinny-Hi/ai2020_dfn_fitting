"""Package the verified workflow and embedded result figures in a concise notebook."""
from __future__ import annotations
import base64
import json
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
source = json.loads((OUT / "261021_Ai2020_간결본_1-2단계_두OCP.ipynb").read_text(encoding="utf-8"))
cells = source["cells"]

def md_table(frame):
    cols=[str(c) for c in frame.columns]
    lines=['| '+' | '.join(cols)+' |','| '+' | '.join(['---']*len(cols))+' |']
    lines += ['| '+' | '.join(str(v).replace('|','/') for v in row)+' |' for row in frame.itertuples(index=False,name=None)]
    return '\n'.join(lines)

def markdown(s):
    return {"cell_type":"markdown","metadata":{},"source":s.splitlines(True)}

def code(s, outputs=None):
    return {"cell_type":"code","metadata":{},"execution_count":None,
            "outputs":outputs or [],"source":s.splitlines(True)}

def image_cell(path, caption):
    text = f'from IPython.display import Image, display\ndisplay(Image(filename="{path.name}"))\n'
    data = base64.b64encode(path.read_bytes()).decode("ascii")
    output = {"output_type":"display_data","data":{"image/png":data,"text/plain":[f"<Figure: {caption}>"]},"metadata":{}}
    return code(text,[output])

selected = [2,3,5,7,8,9,10,12,13]
new=[]
new.append(markdown("""# Ai2020 DFN 최종 간결본 — 두 OCP, 실험값, 충전 피팅·방전 검증

이 노트북은 0.5C/1C/2C 충·방전 실험을 같은 시간 축에서 비교한다. **전체 오차**와 **실제 풀셀 SoC 10–70% 오차**를 별도로 계산한다. 2C 충전은 70%에 도달하지 않으면 도달한 SoC까지만 점수를 낸다. 음극 전위는 분리막 계면에서의 `φ_s,n−φ_e,n`이며 측정 전위가 아닌 DFN 예측이다.

실행 순서: 입력 파일 확인 → OCP/window와 실험 파라미터 → 7변수 민감도/식별도 → Ai2020·ORegan2022 전해질 각각 충전 피팅 → 방전 검증 → 전압/음극 전위 그래프. 저장된 결과와 그림은 아래에서 즉시 확인할 수 있다.
"""))
new.append(markdown("""## 입력 파일

같은 폴더에 두거나 환경변수 `AI2020_INPUT_DIR`로 위치를 지정한다. 이 작업공간에서는 `../work/stage2_inputs` 복사본을 자동 인식한다.

- `260610 저율데이터 1.xlsx`, `260610 저율데이터 2.xlsx`
- `260808 Eneterch셀 열화데이터 v1 6-6.xlsx`
- `260902 Enertech셀 데이터 모음.xlsx` (`CATHODE_GITT`, `ANODE_LOWRATE` 시트)

양극 GITT rest-end는 통합 Excel에서 읽어 큰 원본 기록을 재파싱하지 않는다. 음극 DRT 적분의 `R_ct`와 `D_s` 결과는 아래에 명시한 검증값을 사용한다. 세부 산출 근거는 앞 단계 노트북과 함께 제공된 CSV에 있다.
"""))

for i in selected:
    c = json.loads(json.dumps(cells[i]))
    if i == 2:
        s=''.join(c['source'])
        s=s.replace('pybamm==26.7.1.0','pybamm==26.8.0.0').replace('"26.7.1.0"','"26.8.0.0"')
        s=s.replace('from pathlib import Path\n', 'from pathlib import Path\nimport os\n')
        target='required_files = FULL_CELL_LOW_RATE_FILES + [DYNAMIC_FILE] + CATHODE_GITT_FILES + [MEASURED_ANODE_BOOK]'
        new_target='''# 작업 폴더의 복사본, 노트북 폴더, 또는 AI2020_INPUT_DIR의 원본을 읽는다.
INPUT_DIR = Path(os.environ.get("AI2020_INPUT_DIR", "")).expanduser() if os.environ.get("AI2020_INPUT_DIR") else None
if INPUT_DIR is None:
    candidates = [Path.cwd(), Path.cwd() / "work" / "stage2_inputs", Path.cwd().parent / "work" / "stage2_inputs"]
    INPUT_DIR = next((p for p in candidates if (p / DYNAMIC_FILE).exists()), Path.cwd())
FULL_CELL_LOW_RATE_FILES = [str(INPUT_DIR / p) for p in FULL_CELL_LOW_RATE_FILES]
DYNAMIC_FILE = str(INPUT_DIR / DYNAMIC_FILE)
CATHODE_GITT_FILES = [str(INPUT_DIR / p) for p in CATHODE_GITT_FILES]
MEASURED_ANODE_BOOK = str(INPUT_DIR / MEASURED_ANODE_BOOK)
required_files = FULL_CELL_LOW_RATE_FILES + [DYNAMIC_FILE] + [MEASURED_ANODE_BOOK]'''
        s=s.replace(target,new_target)
        c['source']=s.splitlines(True)
        c['outputs']=[]; c['execution_count']=None
    elif i == 8:
        s=''.join(c['source'])
        compact='''def build_cathode_gitt_ocp_compact(book_name, parameters):
    frame = pd.read_excel(book_name, sheet_name="CATHODE_GITT", header=8)
    disk_area_m2 = np.pi * (0.01 / 2.0) ** 2
    disk_capacity_ah = (FARADAY_CONSTANT
        * float(parameters["Maximum concentration in positive electrode [mol.m-3]"])
        * float(parameters["Positive electrode active material volume fraction"])
        * float(parameters["Positive electrode thickness [m]"])
        * disk_area_m2 / 3600.0)
    lithiation_curves, delithiation_curves = [], []
    for channel in ("7-7", "7-8"):
        block = frame.loc[frame["source_channel"].astype(str).eq(channel)]
        charge = block.loc[block["pulse_direction"].eq("charge_pulse"), ["Qcum_Ah", "Vrest_V"]].to_numpy(float)
        discharge = block.loc[block["pulse_direction"].eq("discharge_pulse"), ["Qcum_Ah", "Vrest_V"]].to_numpy(float)
        delithiation_curves.append(clean_curve(1.0 - charge[:, 0] / disk_capacity_ah, charge[:, 1]))
        discharge_span = discharge[:, 0].max() / disk_capacity_ah
        lithiation_curves.append(clean_curve((1.0 - discharge_span) + discharge[:, 0] / disk_capacity_ah, discharge[:, 1]))
    return mean_overlap([mean_overlap(lithiation_curves), mean_overlap(delithiation_curves)])

CATHODE_GITT_OCP = build_cathode_gitt_ocp_compact(MEASURED_ANODE_BOOK, raw_parameters)
'''
        s=compact+s[s.index('def build_measured_anode_lowrate_ocp('):]
        c['source']=s.splitlines(True); c['outputs']=[]; c['execution_count']=None
    elif c['cell_type']=='code':
        # Static result pages below are the audited outputs. Computational cells
        # are deliberately clear so re-running cannot silently retain old numbers.
        c['outputs']=[]; c['execution_count']=None
    new.append(c)

new.append(markdown("""## 실험 입력과 파라미터

`D_sn=4.1708×10⁻¹⁴`, `D_sp=6.8882×10⁻¹⁴ m²/s`; 음극 N3+N4 local-minimum 적분 `R_ct=27.0009 Ω·cm²`, 양극 `25.387 Ω·cm²`. 입자 반경은 음극 `3.692 µm`, 양극 `4.239 µm`. OCP마다 window와 `c_s,max`를 다시 계산하므로 서로 섞지 않는다. 음극 OCP 두 버전의 비교는 실험값 적용과 피팅 후에도 유지한다.
"""))

stage2=''.join(cells[20]['source'])
stage2=stage2.replace('main()\n','')
new.append(code(stage2))
new.append(markdown("""## 7변수 민감도와 식별도

`D_sn,D_sp,k_n,k_p,b_n,b_s,b_p`에 각각 5% 로그 교란을 주고 0.5C/1C/2C 충·방전 전압의 국소 미분을 계산한다. Jacobian 열 상관계수로 식별도를 본다. `b_n`은 가장 민감하지만 porosity와 구조적으로 결합되어 **고유 물성값으로 해석하지 않는다**. `k_n`과 `k_p`, `b_s`와 `b_p`는 서로 높은 상관을 보인다. 이에 실험 `D_s`와 `k_p`, `b_s`,`b_p`는 고정하고 `k_n,b_n`만 제한적으로 피팅한다.
"""))
lsa=''.join(cells[25]['source'])
new.append(code(lsa))

fit=(ROOT/'work'/'run_stage3_charge_fit.py').read_text(encoding='utf-8')
fit=fit.replace('from run_stage2_two_ocp import cases, kinetic_constants, scope, stage_with_values\n','')
fit=fit.replace('OUT = Path(__file__).resolve().parents[1] / "outputs"','scope = globals()\nOUT = Path.cwd()')
fit=fit.replace('if __name__ == "__main__":\n    main()\n','')
new.append(markdown("""## 충전 피팅 → 방전 검증

두 OCP × 두 전해질을 각각 충전 0.5C/1C/2C 전압으로만 least-squares 피팅한다. 목적함수는 각 율의 140개 공통 시간점 전압 잔차(mV)와 조기 cutoff 페널티를 합친 RMS다. `k_n∈[3×10⁻⁷,3×10⁻⁶]`, `b_n∈[1.5,3.5]`; EIS 초기값과 통상 Bruggeman 1.5를 참고했다. `b_n=3.5`는 현재 Ai2020 기준값 2.914보다 넓게 둔 탐색 여유이며 특정 논문의 측정치가 아니다. 경계에 붙은 결과는 고유 식별로 주장하지 않는다. 3개 시작점을 비교한다. ORegan2022에서는 전해질 확산계수·전도도·양이온 전달수·열역학 계수·초기 전해질 농도만 이식한다. 모든 방전 곡선은 피팅에 쓰지 않고 검증한다.

근거: [COBRAPRO / JES 2024](https://pangea.stanford.edu/ERE/pdf/OnoriPDF/Journals/96.pdf), [ORegan et al. 2022](https://www.sciencedirect.com/science/article/abs/pii/S0013468622008593).
"""))
new.append(code(fit))
new.append(code('main()\n'))

soc=(ROOT/'work'/'recompute_true_soc_metrics.py').read_text(encoding='utf-8')
soc=soc.replace('from run_stage2_two_ocp import scope\n','scope = globals()\n')
soc=soc.replace('ROOT = Path(__file__).resolve().parents[1]\nOUT = ROOT / "outputs"','OUT = Path.cwd()')
soc=soc.replace('for stem, phase in (("Stage3_prefit_two_OCP_curves.csv", "experimental_initial"),\n                    ("Stage3_charge_fit_curves.csv", None)):','for stem, phase in (("Stage3_charge_fit_curves.csv", None),):')
new.append(markdown('## 실제 풀셀 SoC 10–70% 점수와 경계 확인\n\nSoC는 각 branch의 초기 rest 전압으로 정한 초기 SoC에 실험 전류 적산 용량을 더하거나 빼서 계산한다. 시간 축이 겹치는 구간에서만 160점 비교하며, 실제 비교 SoC 하한/상한을 CSV에 기록한다.\n'))
new.append(code(soc))

figcode=(ROOT/'work'/'plot_stage3_final.py').read_text(encoding='utf-8')
figcode=figcode.replace('OUT = Path(__file__).resolve().parents[1] / "outputs"','OUT = Path.cwd()')
new.append(code(figcode))
diagnostic=(ROOT/'work'/'oregan_brugg_diagnostic.py').read_text(encoding='utf-8')
diagnostic=diagnostic.replace('from run_stage3_charge_fit import cases, kinetic_constants, model_stage, scope\n','scope = globals()\n')
diagnostic=diagnostic.replace("OUT=Path(__file__).resolve().parents[1]/'outputs'","OUT=Path.cwd()")
new.append(markdown('## 선택 후 진단 격자 (독립 검증 모델 선정에는 사용하지 않음)\n\nORegan2022 측정 음극 OCP에서 `k_n` 두 값과 `b_n` 네 값을 훑어 충전과 2C 방전 간 trade-off를 확인한다. 재실행 시 비용이 들지만 모든 수치가 코드에 남아 있다.\n'))
new.append(code(diagnostic))

# Verified, inspectable results with embedded images. These remain visible when
# the notebook is opened without re-running the long optimization.
metrics=pd.read_csv(OUT/'Stage3_true_fullcell_SOC10_70_metrics.csv')
best=pd.read_csv(OUT/'Stage3_charge_fit_best.csv')
branch=pd.read_csv(OUT/'Stage3_charge_fit_branch_metrics.csv')
summary=[]
for (ocp,el,phase,direction),g in metrics.groupby(['OCP_case','Electrolyte','Phase','Direction'],sort=False):
    summary.append({'OCP': 'Ai2020' if ocp.startswith('Ai2020') else 'Measured',
                    'Electrolyte':el,'Phase':'Initial' if phase=='experimental_initial' else 'Charge-fit',
                    'Direction':direction,
                    'Full RMSE':round(float((g.Full_RMSE_mV.pow(2).mean())**.5),1),
                    'Full MAE':round(float(g.Full_MAE_mV.mean()),1),
                    'SoC10–70 RMSE':round(float((g.SOC10_70_RMSE_mV.pow(2).mean())**.5),1),
                    'SoC10–70 MAE':round(float(g.SOC10_70_MAE_mV.mean()),1)})
summary=pd.DataFrame(summary)
new.append(markdown('## 재계산된 핵심 결과 (단위 mV)\n\n'+md_table(summary)+'\n\n각 행은 0.5C/1C/2C 세 branch의 집계다. 조건별 전체·SoC10–70 RMSE/MAE와 실제 SoC 범위는 `Stage3_true_fullcell_SOC10_70_metrics.csv`를 재실행하여 확인한다.\n'))
def condition_view(frame):
    f=frame.copy()
    f['OCP']=f.OCP_case.map(lambda x:'Ai2020' if x.startswith('Ai2020') else 'Measured')
    f['조건']=f.C_rate.map(lambda x:f'{x:g}C')+' '+f.Direction.map({'Charge':'충전','Discharge':'방전'})
    f['전체 RMSE/MAE']=f.apply(lambda r:f'{r.Full_RMSE_mV:.1f}/{r.Full_MAE_mV:.1f}',axis=1)
    f['SoC RMSE/MAE']=f.apply(lambda r:f'{r.SOC10_70_RMSE_mV:.1f}/{r.SOC10_70_MAE_mV:.1f}',axis=1)
    f['실제 SoC 범위']=f.apply(lambda r:f'{r.SOC_window_actual_low_pct:.1f}→{r.SOC_window_actual_high_pct:.1f}',axis=1)
    return f[['OCP','조건','전체 RMSE/MAE','SoC RMSE/MAE','실제 SoC 범위']]
prefit_rows=metrics.loc[metrics.Electrolyte.eq('Ai2020_corrected')&metrics.Phase.eq('experimental_initial')]
oregan_rows=metrics.loc[metrics.Electrolyte.eq('ORegan2022')&metrics.Phase.eq('charge_fit_discharge_validation')]
new.append(markdown('### 실험값 반영: 두 OCP × 0.5/1/2C 조건별 전체·SoC 오차 (mV)\n\n'+md_table(condition_view(prefit_rows))+'\n'))
new.append(markdown('### 전해질 변경을 검토한 이유: 피팅 전 음극 전위\n\nAi2020 전해질·실험 파라미터 단계의 충전 최소값은 Ai2020 음극 OCP에서 0.5C `+23.8`, 1C `−21.7`, 2C `−77.8 mV`; 측정 음극 OCP에서 0.5C `+17.7`, 1C `−21.8`, 2C `−78.5 mV`였다. 두 OCP 모두 1–2C의 음극 전위가 0 mV를 하회하여 충전 프로토콜용 예측 모델로 그대로 채택하기 어려웠다. 이에 전해질 수송 parameter set 변경을 검토했다.\n'))
new.append(markdown('### ORegan2022 충전 피팅 후: 독립 방전 포함 조건별 오차 (mV)\n\n'+md_table(condition_view(oregan_rows))+'\n'))
best_display=best[['OCP_case','Electrolyte','kn','brugg_n','Train_voltage_RMSE_mV','Objective_RMS_mV']].copy()
best_display['kn']=best_display['kn'].map(lambda x:f'{x:.3g}')
for col in ['brugg_n','Train_voltage_RMSE_mV','Objective_RMS_mV']:
    best_display[col]=best_display[col].map(lambda x:f'{x:.2f}')
new.append(markdown('## 충전 피팅 파라미터 및 목적함수\n\n'+md_table(best_display)+'\n\n경계 도달과 방전 검증 실패는 결과 해석에서 별도 표기한다.\n'))

pot=branch.loc[branch.Direction.eq('Charge') & branch.Phase.eq('charge_fit_discharge_validation'),
               ['OCP_case','Electrolyte','C_rate','Anode_min_full_mV','Time_anode_below_0V_min']].copy()
pot['Anode_min_full_mV']=pot['Anode_min_full_mV'].map(lambda x:f'{x:.1f}')
pot['Time_anode_below_0V_min']=pot['Time_anode_below_0V_min'].map(lambda x:f'{x:.1f}')
new.append(markdown('## 충전 중 음극 전위 (모델 예측)\n\n'+md_table(pot)+'\n\n0 mV는 열역학적 참조선이며 도금의 실험적 확인값은 아니다. 음극 전위를 직접 측정한 기준 자료가 없어, 전압 적합만으로 도금 안전성을 검증했다고 주장할 수 없다.\n'))
mesh_path=OUT/'Stage3_ORegan_mesh15_vs25.csv'
if mesh_path.exists():
    mesh=pd.read_csv(mesh_path)
    mm=mesh[['OCP_case','C_rate','Direction','Voltage_RMS_mesh_difference_mV','Potential_RMS_mesh_difference_mV']].copy()
    mm['OCP_case']=mm.OCP_case.map(lambda x:'Ai2020' if x.startswith('Ai2020') else 'Measured')
    for c in ['Voltage_RMS_mesh_difference_mV','Potential_RMS_mesh_difference_mV']:
        mm[c]=mm[c].map(lambda x:f'{x:.2f}')
    new.append(markdown('## 격자 정련 재검증: 15점 → 25점\n\n'+md_table(mm)+'\n\n격자 차이보다 2C 방전 오차와 음극 전위 0 mV 하회가 훨씬 크므로 결론은 유지된다.\n'))
new.append(markdown('## 최종 판정\n\nORegan2022 전해질 후보는 충전 전압의 학습 오차를 낮췄지만, 2C 방전에서 큰 오차와 조기 cutoff가 남는다. 모델 음극 전위는 1–2C에서 0 mV를 하회한다. `k_n`도 탐색 상한에 도달했다. 따라서 **현재 실험과 목적함수만으로 최종 목표 달성 또는 도금 안전성 검증을 선언할 수 없다.** 2C 방전 잔차 원인과 reference-electrode 음극 전위를 추가 확인해야 한다.\n'))
new.append(markdown('### 추가 진단: 두 파라미터의 충전·방전 trade-off\n\n측정 음극 OCP·ORegan2022 전해질에서 `k_n=5.06×10⁻⁷`, `b_n=2.5`는 2C 방전 RMSE `22.7 mV`를 보였으나 2C 충전 RMSE `60.7 mV`, 1C 최소 음극 전위 `−19.7 mV`였다. 충전 최적해 `k_n≈3×10⁻⁶`, `b_n≈3.13`은 2C 충전 RMSE `38.3 mV`이나 2C 방전 `124.6 mV`였다. 이 격자 진단은 방전 데이터를 본 뒤 수행했으므로 **독립 검증을 유지하기 위해 모델 선정에는 사용하지 않았다**. 재현 데이터는 `Stage3_ORegan_Brugg_diagnostic.csv`에 있다.\n'))

figure_names=[
    ('Stage3_prefit_two_OCP_overlay_six_curves.png','두 OCP 실험 파라미터 반영'),
    ('Stage3_fit_Ai2020_anode__measured_cathode_electrolyte_comparison.png','Ai2020 OCP 전해질 비교'),
    ('Stage3_fit_measured_anode__measured_cathode_electrolyte_comparison.png','Measured OCP 전해질 비교'),
    ('Stage3_fit_Ai2020_anode__measured_cathode_ORegan2022_anode_potential.png','ORegan2022 음극 전위'),
]
new.append(markdown('## 검증 그래프\n'))
for name,caption in figure_names:
    path=OUT/name
    if path.exists():
        new.append(markdown(f'### {caption}\n'))
        new.append(image_cell(path,caption))

notebook={'cells':new,'metadata':source['metadata'],'nbformat':4,'nbformat_minor':5}
target=OUT/'261021_Ai2020_DFN_최종_간결본.ipynb'
target.write_text(json.dumps(notebook,ensure_ascii=False,indent=1),encoding='utf-8')
print(target,len(new),'cells',target.stat().st_size,'bytes')
