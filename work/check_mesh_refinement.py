"""Recompute key ORegan fitted branches on a 25-point DFN mesh."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from run_stage3_charge_fit import cases, model_stage, scope

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs'
nb=json.loads((OUT/'261020_Ai2020_1단계_OCP_Stoichiometry_Measured_Anode_비교.ipynb').read_text(encoding='utf-8'))
src=''.join(nb['cells'][13]['source'])
start=src.index('def run_dfn_branch(')
end=src.index('def compare_with_experiment(',start)
fn=src[start:end]
needle='var_pts={"x_n": 15, "x_s": 15, "x_p": 15, "r_n": 15, "r_p": 15}'
assert needle in fn
fn=fn.replace(needle,'var_pts={"x_n": 25, "x_s": 25, "x_p": 25, "r_n": 25, "r_p": 25}')
env=dict(scope)
exec(compile(fn,'refined_DFN','exec'),env)
run25=env['run_dfn_branch']
best=pd.read_csv(OUT/'Stage3_charge_fit_best.csv')
curves=pd.read_csv(OUT/'Stage3_charge_fit_curves.csv')
rows=[]
for ocp_name,case in cases().items():
    fit=best.loc[best.OCP_case.eq(ocp_name)&best.Electrolyte.eq('ORegan2022')].iloc[0]
    stage=model_stage(case,'ORegan2022',fit.kn,fit.brugg_n)
    for rate in (1.,2.):
        for charge in (True,False):
            direction='Charge' if charge else 'Discharge'
            r=run25(stage,rate,charge)
            ref=curves.loc[curves.OCP_case.eq(ocp_name)&curves.Electrolyte.eq('ORegan2022')
                           &curves.Phase.eq('charge_fit_discharge_validation')
                           &curves.C_rate.eq(rate)&curves.Direction.eq(direction)&curves.Series.eq('DFN')]
            t=np.linspace(0,min(r['time_min'][-1],ref.Time_min.iloc[-1])*.995,400)
            dv=1000*(np.interp(t,r['time_min'],r['voltage_v'])-np.interp(t,ref.Time_min,ref.Voltage_V))
            dp=1000*(np.interp(t,r['time_min'],r['anode_potential_v'])-np.interp(t,ref.Time_min,ref.Anode_potential_mV/1000))
            exp=scope['experimental_curves'][(rate,charge)]
            score=scope['compare_with_experiment'](exp,r)
            rows.append({'OCP_case':ocp_name,'C_rate':rate,'Direction':direction,
                'Mesh15_duration_min':float(ref.Time_min.iloc[-1]),'Mesh25_duration_min':float(r['time_min'][-1]),
                'Voltage_RMS_mesh_difference_mV':float(np.sqrt(np.mean(dv**2))),
                'Voltage_max_mesh_difference_mV':float(np.max(abs(dv))),
                'Potential_RMS_mesh_difference_mV':float(np.sqrt(np.mean(dp**2))),
                'Potential_max_mesh_difference_mV':float(np.max(abs(dp))),
                'Mesh25_full_RMSE_mV':score['full_RMSE_mV'],
                'Mesh25_min_anode_mV':1000*float(np.min(r['anode_potential_v']))})
            pd.DataFrame(rows).to_csv(OUT/'Stage3_ORegan_mesh15_vs25.csv',index=False,encoding='utf-8-sig')
            print(rows[-1],flush=True)
