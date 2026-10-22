"""Run the compact notebook's DFN settings for three Stage-1 dynamic variants."""
import json, os
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
NB=next((ROOT/'outputs').glob('261020_Ai2020_1*.ipynb'))
OUT=ROOT/'outputs'
INPUTS=ROOT/'work'/'run_inputs'
font_manager.fontManager.addfont('C:/Windows/Fonts/malgun.ttf')
font_manager.fontManager.addfont('C:/Windows/Fonts/malgunbd.ttf')
plt.rcParams.update({'font.family':'Malgun Gothic','axes.unicode_minus':False})

nb=json.loads(NB.read_text(encoding='utf-8'))
ns={'__name__':'__main__'}
os.chdir(INPUTS)
for i in [2,3,5,7,8,9,12,13]:
    code=''.join(nb['cells'][i]['source'])
    code='\n'.join(line for line in code.splitlines() if not line.lstrip().startswith('%'))
    if i==8:code=code.split('def build_measured_anode_lowrate_ocp',1)[0]
    if i==9:code=code.split('published_vector =',1)[0]
    if i==12:code=code.split('improved_parameters, improved_csmax_n',1)[0]
    print('Loading notebook cell',i,flush=True)
    exec(compile(code,f'{NB.name}:cell{i}','exec'),ns)

pv=ns['common_parameters']; cap=ns['CELL_CAPACITY_AH']; faraday=ns['FARADAY_CONSTANT']; area=ns['OVERLAP_AREA_M2']
neg=ns['AI2020_OCP']['negative']; pos=ns['AI2020_OCP']['positive']; measured_pos=ns['CATHODE_GITT_OCP']
wa=ns['PUBLISHED_WINDOW']
wb={'x0':0.0031795298633619,'x100':0.9850660533269036,'y0':0.8677683066214189,'y100':0.441823602664651}
wc={'x0':0.003803,'x100':0.954167,'y0':0.970716,'y100':0.431646}
bp,bn,bpmax=ns['build_candidate_parameters'](pv,wb,pos)
cp,cn,cpmax=ns['build_candidate_parameters'](pv,wc,measured_pos)
stages={
 'A · 공통 물리 보정':{'parameters':pv,'window':wa,
   'c_s_max_n':float(pv['Maximum concentration in negative electrode [mol.m-3]']),
   'c_s_max_p':float(pv['Maximum concentration in positive electrode [mol.m-3]']),
   'anode_ocp':neg,'cathode_ocp':pos},
 'B · Ai2020 OCP + window 재최적화':{'parameters':bp,'window':wb,'c_s_max_n':bn,'c_s_max_p':bpmax,'anode_ocp':neg,'cathode_ocp':pos},
 'C · 측정 양극 OCP + window 재최적화':{'parameters':cp,'window':wc,'c_s_max_n':cn,'c_s_max_p':cpmax,'anode_ocp':neg,'cathode_ocp':measured_pos},
}
colors={'A · 공통 물리 보정':'#87909B','B · Ai2020 OCP + window 재최적화':'#2371B6','C · 측정 양극 OCP + window 재최적화':'#D66422'}
styles={'A · 공통 물리 보정':'--','B · Ai2020 OCP + window 재최적화':'-.','C · 측정 양극 OCP + window 재최적화':'-'}
sim={};rows=[];curve_rows=[]
for name,stage in stages.items():
    for rate in (0.5,1.0,2.0):
        for charge in (True,False):
            print('Simulating',name,rate,'Charge' if charge else 'Discharge',flush=True)
            s=ns['run_dfn_branch'](stage,rate,charge)
            sim[(name,rate,charge)]=s
            m=ns['compare_with_experiment'](ns['experimental_curves'][(rate,charge)],s)
            rows.append({'stage':name,'C_rate':rate,'direction':'Charge' if charge else 'Discharge',**m})
            for t,v in zip(s['time_min'],s['voltage_v']):
                curve_rows.append({'stage':name,'C_rate':rate,'direction':'Charge' if charge else 'Discharge','series':'DFN','time_min':t,'voltage_V':v})
for rate in (0.5,1.0,2.0):
    for charge in (True,False):
        e=ns['experimental_curves'][(rate,charge)]
        for t,v in zip(e['time_min'],e['voltage_v']):
            curve_rows.append({'stage':'Experiment','C_rate':rate,'direction':'Charge' if charge else 'Discharge','series':'Experiment','time_min':t,'voltage_V':v})

metrics=pd.DataFrame(rows)
metrics.to_csv(OUT/'Stage1_De_RpRn_OCP_window_양극OCP_3단계_동적전압_RMSE_MAE.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(curve_rows).to_csv(OUT/'Stage1_De_RpRn_OCP_window_양극OCP_3단계_동적전압_곡선.csv',index=False,encoding='utf-8-sig')
agg=metrics.groupby('stage')['full_RMSE_mV'].apply(lambda a:float(np.sqrt(np.mean(np.asarray(a)**2))))
print('Six-curve full RMSE:',agg.to_dict(),flush=True)

fig,axes=plt.subplots(2,3,figsize=(15.5,8.8),sharey=True,dpi=160)
fig.subplots_adjust(left=.065,right=.985,top=.85,bottom=.17,wspace=.14,hspace=.24)
for col,rate in enumerate((0.5,1.0,2.0)):
  for row,charge in enumerate((True,False)):
    ax=axes[row,col];e=ns['experimental_curves'][(rate,charge)]
    ax.plot(e['time_min'],e['voltage_v'],color='#20262D',lw=2.65,label='실험',zorder=5)
    for name in stages:
      d=sim[(name,rate,charge)]
      ax.plot(d['time_min'],d['voltage_v'],color=colors[name],ls=styles[name],lw=2.25,label=name)
    ax.set_title(f'{rate:g}C {"충전" if charge else "방전"}',fontsize=17,weight='bold',pad=9)
    ax.set_ylim(2.95,4.27);ax.grid(alpha=.24);ax.set_axisbelow(True)
    ax.tick_params(labelsize=11)
    if row==1:ax.set_xlabel('경과시간 (분)',fontsize=13)
    if col==0:ax.set_ylabel('셀 전압 (V)',fontsize=13)
handles,labels=axes[0,0].get_legend_handles_labels()
fig.legend(handles,labels,loc='upper center',ncol=4,bbox_to_anchor=(.51,.945),frameon=False,fontsize=12)
fig.text(.065,.975,'공통 입력: Ai2020 De×1e-4 · rn 3.692 µm · rp 4.239 µm · Ds/k/Bruggeman 미피팅',fontsize=15,color='#156082',weight='bold',va='top')
fig.text(.065,.087,'6곡선 전체 RMSE (mV)   '+'   →   '.join(f'{k[0]} {agg[k]:.2f}' for k in stages),fontsize=17,color='#156082',weight='bold')
fig.text(.065,.04,'B는 Ai2020 양극 OCP를 유지하고 window/cs,max만 재계산. C에서 측정 양극 GITT OCP를 반영해 window/cs,max를 다시 계산.',fontsize=11.5,color='#4A5560')
png=OUT/'Stage1_De_RpRn_OCP_window_양극OCP_3단계_동적전압_6곡선.png'
svg=OUT/'Stage1_De_RpRn_OCP_window_양극OCP_3단계_동적전압_6곡선.svg'
fig.savefig(png,dpi=200,facecolor='white');fig.savefig(svg,facecolor='white');plt.close(fig)
print('Saved',png,svg,flush=True)
