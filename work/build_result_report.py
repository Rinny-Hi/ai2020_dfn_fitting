"""Write a Korean audit report from the verified tabular outputs."""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs'

def md_table(frame):
    cols=[str(c) for c in frame.columns]
    lines=['| '+' | '.join(cols)+' |','| '+' | '.join(['---']*len(cols))+' |']
    lines += ['| '+' | '.join(str(v).replace('|','/') for v in row)+' |' for row in frame.itertuples(index=False,name=None)]
    return '\n'.join(lines)
met=pd.read_csv(OUT/'Stage3_true_fullcell_SOC10_70_metrics.csv')
best=pd.read_csv(OUT/'Stage3_charge_fit_best.csv')
branch=pd.read_csv(OUT/'Stage3_charge_fit_branch_metrics.csv')
mesh_path=OUT/'Stage3_ORegan_mesh15_vs25.csv'
mesh=pd.read_csv(mesh_path) if mesh_path.exists() else None
pre=met.loc[met.Electrolyte.eq('Ai2020_corrected') & met.Phase.eq('experimental_initial')].copy()

def short(s): return 'Ai2020 음극' if str(s).startswith('Ai2020') else '측정 음극'
def condition_table(df):
    f=df.copy()
    f['OCP']=f.OCP_case.map(short)
    f['전해질']=f.Electrolyte.map({'Ai2020_corrected':'Ai2020','ORegan2022':'ORegan2022'})
    f['조건']=f.C_rate.map(lambda x:f'{x:g}C')+' '+f.Direction.map({'Charge':'충전','Discharge':'방전'})
    cols={'Full_RMSE_mV':'전체 RMSE','Full_MAE_mV':'전체 MAE',
      'SOC10_70_RMSE_mV':'SoC10–70 RMSE','SOC10_70_MAE_mV':'SoC10–70 MAE',
      'SOC_window_actual_low_pct':'실제 SoC 시작','SOC_window_actual_high_pct':'실제 SoC 끝'}
    f=f.rename(columns=cols)
    use=['OCP','전해질','조건',*cols.values()]
    for c in cols.values():f[c]=f[c].map(lambda x:f'{x:.1f}')
    return md_table(f[use])

line=[]
line.append('# 261021 Ai2020 DFN 검증 결과\n')
line.append('**판정:** 충전 전압만으로 맞춘 ORegan2022 전해질 후보는 2C 방전과 음극 전위 검증에서 최종 목표를 충족했다고 볼 수 없다. 이 결과를 최종 선정 성공으로 표시하지 않았다. 두 음극 OCP 모두 비교했고 기존 PPT의 비어 있던 영역에 결과와 한계를 추가했다.\n')
line.append('## 1. 실험값 반영, 피팅 전 (Ai2020 전해질)\n')
line.append('GITT `D_sn=4.1708×10⁻¹⁴`, `D_sp=6.8882×10⁻¹⁴ m²/s`; 음극 DRT N3+N4 local-minimum `R_ct=27.0009 Ω·cm²`, 양극 `25.387 Ω·cm²`. 음극 `r=3.692 µm`, 양극 `4.239 µm`. OCP마다 window와 `c_s,max`를 별도 계산했다.\n')
line.append('오차는 각 조건에서 실험·DFN 공통 시간 구간의 전압 차이(mV)다. 실제 풀셀 SoC는 branch 초기 rest 전압에서 얻은 초기 SoC와 실험 적산용량으로 계산했다. 2C 충전은 전압 cutoff로 SoC 70% 미도달이며 표의 실제 끝값까지만 평가한다.\n')
line.append(condition_table(pre)+'\n')
line.append('**전해질 변경을 검토한 음극 전위 근거:** Ai2020 전해질·실험 파라미터 단계에서 충전 최소 음극 전위는 Ai2020 음극 OCP 기준 0.5C `+23.8 mV`, 1C `−21.7 mV`, 2C `−77.8 mV`이고, 측정 음극 OCP 기준 0.5C `+17.7 mV`, 1C `−21.8 mV`, 2C `−78.5 mV`였다. 두 OCP 모두 1–2C에서 0 mV를 하회하여 충전 프로토콜용 음극 전위 예측 모델로 그대로 채택하기 어려웠다. 따라서 전해질 수송 parameter set 변경을 검토했다. ORegan2022 변경 뒤에도 하회가 남으므로 변경만으로 문제 해결을 주장할 수 없다.\n')
line.append('## 2. 민감도·식별도 및 충전 피팅\n')
line.append('7개 후보 `D_sn,D_sp,k_n,k_p,b_n,b_s,b_p`의 5% 로그 교란을 0.5C/1C/2C 충·방전 전압에 적용했다. `b_n`의 RMS 전압변화는 약 9 mV, `k_n`·`k_p`는 약 2 mV, `D_s`는 0.12 mV 이하였다. `k_n–k_p`와 `b_s–b_p` Pearson 상관은 약 0.99–1.00이다. 따라서 `D_s,k_p,b_s,b_p`는 측정/기존값으로 고정하고 `k_n,b_n`만 피팅했다. `b_n`은 porosity와 구조적으로 중복되어 유효계수로만 해석한다.\n')
f=best[['OCP_case','Electrolyte','kn','brugg_n','Train_voltage_RMSE_mV','Objective_RMS_mV']].copy()
f['OCP_case']=f.OCP_case.map(short)
f['kn']=f.kn.map(lambda x:f'{x:.3g}')
for c in ('brugg_n','Train_voltage_RMSE_mV','Objective_RMS_mV'):f[c]=f[c].map(lambda x:f'{x:.2f}')
line.append(md_table(f)+'\n')
line.append('충전 세 율에서 각 140점의 전압 잔차와 조기 cutoff 페널티를 결합한 RMS 목적함수다. `k_n` 범위 `[3×10⁻⁷,3×10⁻⁶]`는 EIS 초기값 약 `5–6×10⁻⁷`을 포괄하는 탐색 범위다. `b_n` 범위 `[1.5,3.5]`의 하한은 통상 Bruggeman 1.5, 상한은 현재 Ai2020 기준값 2.914보다 넓게 잡은 **탐색 여유**이며 3.5가 특정 논문의 측정값이라는 뜻은 아니다. 각 조합 3개 시작점으로 반복했고, 상한에 닿는 `k_n`/`b_n`은 잘 식별된 최적값으로 주장할 수 없다.\n')
line.append('## 3. 전해질 변경 후 독립 방전 검증\n')
line.append('ORegan2022에서는 전해질 확산계수, 전도도, 양이온 전달수, 열역학 계수, 초기 농도만 가져왔다. 전극 OCP/기하/활물질 농도와 실험 `D_s`·`k` 출발값은 현재 셀에 맞춘 상태로 유지했다. 아래는 **피팅에 쓰지 않은 방전** 결과다.\n')
dis=met.loc[met.Phase.eq('charge_fit_discharge_validation')&met.Direction.eq('Discharge')].copy()
line.append(condition_table(dis)+'\n')
cap=branch.loc[branch.Phase.eq('charge_fit_discharge_validation')&branch.Direction.eq('Discharge')&branch.C_rate.eq(2),
               ['OCP_case','Electrolyte','full_RMSE_mV','capacity_error_pct']].copy()
cap['OCP_case']=cap.OCP_case.map(short)
cap['full_RMSE_mV']=cap.full_RMSE_mV.map(lambda x:f'{x:.1f}')
cap['capacity_error_pct']=cap.capacity_error_pct.map(lambda x:f'{x:.1f}')
line.append('2C 방전의 전체 오차·용량 차이:\n\n'+md_table(cap)+'\n')
line.append('## 4. 충전 중 음극 전위\n')
pot=branch.loc[branch.Phase.eq('charge_fit_discharge_validation')&branch.Direction.eq('Charge'),
               ['OCP_case','Electrolyte','C_rate','Anode_min_full_mV','Time_anode_below_0V_min']].copy()
pot['OCP_case']=pot.OCP_case.map(short)
for c in ('Anode_min_full_mV','Time_anode_below_0V_min'):pot[c]=pot[c].map(lambda x:f'{x:.1f}')
line.append(md_table(pot)+'\n')
line.append('음극 전위는 `Negative electrode surface potential difference at separator interface`, 즉 `φ_s,n−φ_e,n`의 DFN 예측이다. 실측 reference-electrode 전위가 없으며 0 mV 하회가 곧 도금 실험 확인을 뜻하지는 않는다. 다만 1–2C에서 음수로 내려가는 후보를 현재 목적에 맞는 안전 예측 모델이라고 검증할 수는 없다.\n')
line.append('**추가 진단(선정에 사용하지 않음):** 측정 음극 OCP·ORegan 전해질에서 `k_n=5.06×10⁻⁷`, `b_n=2.5`로 고정하면 2C 방전 RMSE가 `22.7 mV`로 낮아지지만 2C 충전 RMSE는 `60.7 mV`, 1C 충전 최소 음극 전위는 `−19.7 mV`다. 반대로 충전 최적해(`k_n≈3×10⁻⁶`, `b_n≈3.13`)는 2C 충전 RMSE `38.3 mV`이나 2C 방전 RMSE `124.6 mV`다. 이 8점 격자 진단은 충전·방전 간 trade-off를 보여주지만, 방전 결과를 본 뒤 파라미터를 선택하면 독립 검증이 아니므로 최종 모델 재선정에는 사용하지 않았다.\n')
if mesh is not None:
    line.append('## 5. 격자 정련 확인\n')
    m=mesh[['OCP_case','C_rate','Direction','Voltage_RMS_mesh_difference_mV','Potential_RMS_mesh_difference_mV','Mesh25_full_RMSE_mV']].copy()
    m['OCP_case']=m.OCP_case.map(short)
    line.append(md_table(m.round(2))+'\n')
line.append('## 출처 및 재현 파일\n')
line.append('- COBRAPRO의 상관 기반 식별도·파라미터 분리: https://pangea.stanford.edu/ERE/pdf/OnoriPDF/Journals/96.pdf\n')
line.append('- ORegan2022 원논문: https://www.sciencedirect.com/science/article/abs/pii/S0013468622008593\n')
line.append('- 계산 코드·입력 확인·그림·전체 CSV: `261021_Ai2020_DFN_최종_간결본.ipynb` 및 함께 저장된 `Stage3_*.csv`/PNG.\n')
target=OUT/'261021_DFN_검증결과_요약.md'
target.write_text('\n'.join(line),encoding='utf-8')
print(target)
