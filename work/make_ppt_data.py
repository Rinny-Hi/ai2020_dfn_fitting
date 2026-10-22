"""Small verified data payload for editable PPT tables and captions."""
import json
from pathlib import Path
import pandas as pd
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs'
met=pd.read_csv(OUT/'Stage3_true_fullcell_SOC10_70_metrics.csv')
best=pd.read_csv(OUT/'Stage3_charge_fit_best.csv')
branch=pd.read_csv(OUT/'Stage3_charge_fit_branch_metrics.csv')
par=pd.read_csv(OUT/'Stage2_two_OCP_parameters.csv')
sens=pd.read_csv(OUT/'Stage2_seven_parameter_sensitivity.csv')

def short(x): return 'Ai2020' if str(x).startswith('Ai2020') else 'Measured'
payload={'prefit':[],'scenario':[],'fit':[],'potential':[],'window':[],'sensitivity':[],'branch':[]}
for _,r in met.loc[met.Electrolyte.eq('Ai2020_corrected') & met.Phase.eq('experimental_initial')].iterrows():
    payload['prefit'].append({'ocp':short(r.OCP_case),'rate':float(r.C_rate),'direction':r.Direction,
       'full_rmse':round(r.Full_RMSE_mV,1),'full_mae':round(r.Full_MAE_mV,1),
       'soc_rmse':round(r.SOC10_70_RMSE_mV,1),'soc_mae':round(r.SOC10_70_MAE_mV,1),
       'soc_low':round(r.SOC_window_actual_low_pct,1),'soc_high':round(r.SOC_window_actual_high_pct,1)})
for (ocp,el,phase,direction),g in met.groupby(['OCP_case','Electrolyte','Phase','Direction'],sort=False):
    payload['scenario'].append({'ocp':short(ocp),'electrolyte':el,'phase':phase,'direction':direction,
      'full_rmse':round(float(np.sqrt(np.mean(g.Full_RMSE_mV**2))),1),
      'full_mae':round(float(np.mean(g.Full_MAE_mV)),1),
      'soc_rmse':round(float(np.sqrt(np.mean(g.SOC10_70_RMSE_mV**2))),1),
      'soc_mae':round(float(np.mean(g.SOC10_70_MAE_mV)),1)})
for _,r in best.iterrows():
    payload['fit'].append({'ocp':short(r.OCP_case),'electrolyte':r.Electrolyte,
      'kn':f'{r.kn:.3g}','bn':round(r.brugg_n,3),
      'train_rmse':round(r.Train_voltage_RMSE_mV,1),'objective':round(r.Objective_RMS_mV,1),
      'boundary': bool(r.kn>2.999e-6 or r.brugg_n>3.499)})
for _,r in branch.loc[branch.Direction.eq('Charge')].iterrows():
    payload['potential'].append({'ocp':short(r.OCP_case),'electrolyte':r.Electrolyte,
       'phase':r.Phase,'rate':float(r.C_rate),'min_mv':round(r.Anode_min_full_mV,1),
       'time_below_zero':round(r.Time_anode_below_0V_min,1)})
for _,r in branch.iterrows():
    payload['branch'].append({'ocp':short(r.OCP_case),'electrolyte':r.Electrolyte,
       'phase':r.Phase,'rate':float(r.C_rate),'direction':r.Direction,
       'full_rmse':round(r.full_RMSE_mV,1),'capacity_error':round(r.capacity_error_pct,1)})
for _,r in par.iterrows():
    payload['window'].append({'ocp':short(r.OCP_case),'x0':round(r.x0,6),'x100':round(r.x100,6),
      'y0':round(r.y0,6),'y100':round(r.y100,6),
      'csn':round(r.c_s_max_n_mol_m3),'csp':round(r.c_s_max_p_mol_m3)})
for _,r in sens.iterrows():
    if r.Parameter in ('brugg_n','kn','kp','brugg_p','brugg_s','Dsp','Dsn'):
        payload['sensitivity'].append({'ocp':short(r.OCP_case),'parameter':r.Parameter,
          'mV_5pct':round(r.Voltage_change_for_5pct_mV_RMS,2)})
(ROOT/'work'/'ppt_data.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2),encoding='utf-8')
print(ROOT/'work'/'ppt_data.json')
