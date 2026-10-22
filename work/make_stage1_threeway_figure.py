"""Recreate the three comparable Stage-1 qOCV curves from the compact notebook."""
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = next((ROOT / 'outputs').glob('261020_Ai2020_1*.ipynb'))
INPUTS = ROOT / 'work' / 'run_inputs'
OUT = ROOT / 'outputs'

font_path = Path('C:/Windows/Fonts/NotoSansKR-VF.ttf')
if font_path.exists():
    font_manager.fontManager.addfont(str(font_path))
    plt.rcParams['font.family'] = 'Noto Sans KR'
plt.rcParams['axes.unicode_minus'] = False

nb = json.loads(NOTEBOOK.read_text(encoding='utf-8'))
namespace = {'__name__': '__main__'}
os.chdir(INPUTS)
for index in [2, 3, 5, 7, 8]:
    code = ''.join(nb['cells'][index]['source'])
    code = '\n'.join(line for line in code.splitlines() if not line.lstrip().startswith('%'))
    if index == 3:
        code = code.replace('experimental_curves = load_dynamic_data(DYNAMIC_FILE)', 'experimental_curves = {}  # dynamic data not needed for qOCV figure')
    if index == 8:
        code = code.split('def build_measured_anode_lowrate_ocp', 1)[0]
    exec(compile(code, f'{NOTEBOOK.name}:cell{index}', 'exec'), namespace)

soc = namespace['SOC_GRID']
observed = np.asarray(namespace['low_rate_data']['qocv'], float)
neg = namespace['AI2020_OCP']['negative']
pos_ai = namespace['AI2020_OCP']['positive']
pos_meas = namespace['CATHODE_GITT_OCP']
stages = [
    ('A  공통 물리 보정 + 공개 window', np.array([0.0065, 0.8400, 0.4350, 0.9651]), pos_ai, '#808A97', '--'),
    ('B  Ai2020 OCP + window 재최적화', np.array([0.003180, 0.985066, 0.441824, 0.867768]), pos_ai, '#2371B6', '-'),
    ('C  측정 양극 OCP + window 재최적화', np.array([0.003803, 0.954167, 0.431646, 0.970716]), pos_meas, '#D66422', '-'),
]

def full_cell_ocv(z, ncurve, pcurve):
    x0,x100,y100,y0 = z
    n = np.interp(x0 + soc*(x100-x0), *ncurve)
    p = np.interp(y0 - soc*(y0-y100), *pcurve)
    return p-n

fig, axes = plt.subplots(1, 2, figsize=(14.4, 6.4), dpi=200, gridspec_kw={'width_ratios':[1.08,1]})
fig.subplots_adjust(left=.07,right=.985,top=.85,bottom=.19,wspace=.19)
ax, er = axes
ax.plot(100*soc, observed, color='#20262D', linewidth=2.8, label='실측 C/50 qOCV', zorder=5)
region=(soc>=.02)&(soc<=.98)
metrics=[]
for label,z,pcurve,color,ls in stages:
    v=full_cell_ocv(z,neg,pcurve)
    rmse=float(np.sqrt(np.mean((1000*(v[region]-observed[region]))**2)))
    metrics.append({'stage':label,'rmse_2_98_mV':rmse,'window':{'x0':z[0],'x100':z[1],'y100':z[2],'y0':z[3]}})
    ax.plot(100*soc,v,color=color,linestyle=ls,linewidth=2.25,label=f'{label}  ({rmse:.2f} mV)')
    er.plot(100*soc,1000*(v-observed),color=color,linestyle=ls,linewidth=2.25,label=label)

ax.set_title('C/50 qOCV 재구성',fontsize=18,weight='bold',pad=12)
er.set_title('실측과의 전압 차이',fontsize=18,weight='bold',pad=12)
ax.set(xlim=(0,100),ylim=(2.98,4.26),xlabel='풀셀 SoC (%)',ylabel='전압 (V)')
er.set(xlim=(0,100),ylim=(-95,60),xlabel='풀셀 SoC (%)',ylabel='모델 − 실측 (mV)')
er.axhline(0,color='#2E3338',linewidth=1.1)
for a in axes:
    a.grid(color='#DADEE3',alpha=.65)
    a.tick_params(labelsize=12)
    a.set_axisbelow(True)
ax.legend(loc='lower right',fontsize=11,framealpha=.94)
fig.text(.07,.945,'공통 입력  ·  Ai2020 Dₑ 함수 ×10⁻⁴  |  rₙ 3.692 µm  |  rₚ 4.239 µm',fontsize=15,color='#156082',weight='bold')
fig.text(.07,.072,'qOCV RMSE (SoC 2–98%):  A 39.55  →  B 16.35  →  C 6.24 mV',fontsize=17,color='#156082',weight='bold')
fig.text(.07,.025,'A/B는 Ai2020 양극 OCP 유지. C에서만 측정 양극 GITT OCP를 적용하고 window를 다시 최적화.',fontsize=11.5,color='#4A5560')

png=OUT/'Stage1_De_RpRn_OCP_window_양극OCP_3단계_qOCV_비교.png'
svg=OUT/'Stage1_De_RpRn_OCP_window_양극OCP_3단계_qOCV_비교.svg'
fig.savefig(png,dpi=220,facecolor='white')
fig.savefig(svg,facecolor='white')
plt.close(fig)
(OUT/'Stage1_De_RpRn_OCP_window_양극OCP_3단계_qOCV_수치.json').write_text(json.dumps(metrics,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'png':str(png),'svg':str(svg),'metrics':metrics},ensure_ascii=False,indent=2))
