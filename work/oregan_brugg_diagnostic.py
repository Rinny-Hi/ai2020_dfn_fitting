"""Diagnostic ORegan sweep; no discharge residual is used to tune the fit."""
from pathlib import Path
import numpy as np
import pandas as pd
from run_stage3_charge_fit import cases, kinetic_constants, model_stage, scope

OUT=Path(__file__).resolve().parents[1]/'outputs'
rows=[]
for ocp_name, case in cases().items():
    if not ocp_name.startswith('measured_anode'): continue
    k0=kinetic_constants(case)['kn']
    for kn in (k0, 3e-6):
        for bn in (1.5,2.0,2.5,3.13):
            stage=model_stage(case,'ORegan2022',kn,bn)
            for rate,charge in ((1.0,True),(2.0,True),(2.0,False)):
                try:
                    run=scope['run_dfn_branch'](stage,rate,charge)
                    exp=scope['experimental_curves'][(rate,charge)]
                    met=scope['compare_with_experiment'](exp,run)
                    rows.append({'OCP_case':ocp_name,'kn':kn,'brugg_n':bn,'C_rate':rate,
                                 'Direction':'Charge' if charge else 'Discharge',
                                 'Full_RMSE_mV':met['full_RMSE_mV'],
                                 'Capacity_error_pct':met['capacity_error_pct'],
                                 'Min_anode_mV':1000*float(np.min(run['anode_potential_v']))})
                except Exception as exc:
                    rows.append({'OCP_case':ocp_name,'kn':kn,'brugg_n':bn,'C_rate':rate,
                                 'Direction':'Charge' if charge else 'Discharge',
                                 'error':repr(exc)})
            pd.DataFrame(rows).to_csv(OUT/'Stage3_ORegan_Brugg_diagnostic.csv',index=False,encoding='utf-8-sig')
            print(ocp_name,kn,bn,rows[-3:],flush=True)
