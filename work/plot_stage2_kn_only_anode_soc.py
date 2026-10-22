"""Compare Stage 2 separator-side anode potential vs estimated full-cell SoC.

Both versions use the same Ai2020 electrolyte, OCP, geometry, GITT D_s,
cathode k_p and fixed negative electrolyte Bruggeman coefficient.
"""
from pathlib import Path
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybamm

from run_stage3_charge_fit import cases, model_stage, scope

OUT = Path(__file__).resolve().parents[1] / 'outputs'
CASE = cases()['Ai2020_anode__measured_cathode']
fit = pd.read_csv(OUT / 'Stage2_user_kn554_charge_fit_result.csv').iloc[0]
assert fit.fit_parameters == 'kn only'
assert np.isclose(fit.fitted_brugg_n, fit.start_brugg_n)
BN = 2.914
VERSIONS = [('Initial', 5.54e-7), ('k_n fitted', float(fit.fitted_kn))]
COLORS = {'Initial': '#87909B', 'k_n fitted': '#D66422'}

rows = []
summary = []
for label, kn in VERSIONS:
    stage = model_stage(CASE, 'Ai2020_corrected', kn, BN)
    for rate in (0.5, 1.0, 2.0):
        print('SIMULATE', label, rate, flush=True)
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=pybamm.SolverWarning)
            run = scope['run_dfn_branch'](stage, rate, True)
        soc = np.asarray(run['full_cell_soc_pct'], float)
        pot = 1000 * np.asarray(run['anode_potential_v'], float)
        time = np.asarray(run['time_min'], float)
        for s, p, t in zip(soc, pot, time):
            rows.append({'Version': label, 'k_n': kn, 'brugg_n': BN, 'C_rate': rate,
                         'Estimated_full_cell_SoC_pct': s,
                         'Anode_potential_separator_mV': p, 'Time_min': t})
        neg = pot < 0
        summary.append({'Version': label, 'k_n': kn, 'brugg_n': BN, 'C_rate': rate,
                        'Start_SoC_pct': float(soc[0]), 'End_SoC_pct': float(soc[-1]),
                        'Min_potential_mV': float(np.min(pot)),
                        'SoC_at_min_pct': float(soc[np.argmin(pot)]),
                        'First_SoC_below_0_pct': float(soc[np.flatnonzero(neg)[0]]) if neg.any() else np.nan,
                        'Time_below_0_min': float(np.trapezoid(neg.astype(float), time))})

data = pd.DataFrame(rows)
data.to_csv(OUT / 'Stage2_kn554_initial_vs_kn_only_fit_anode_potential_SoC_curves.csv', index=False)
pd.DataFrame(summary).to_csv(OUT / 'Stage2_kn554_initial_vs_kn_only_fit_anode_potential_summary.csv', index=False)

fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.6), dpi=160, sharex='col')
fig.subplots_adjust(left=.075, right=.985, top=.88, bottom=.10, wspace=.17, hspace=.16)
max_pot = max(300.0, float(np.ceil(data.Anode_potential_separator_mV.max()/100)*100))
min_pot = min(-110.0, float(np.floor(data.Anode_potential_separator_mV.min()/20)*20))
for col, rate in enumerate((.5, 1., 2.)):
    for row in (0, 1):
        ax = axes[row, col]
        for label, kn in VERSIONS:
            d = data.loc[(data.Version == label) & (data.C_rate == rate)]
            ax.plot(d.Estimated_full_cell_SoC_pct, d.Anode_potential_separator_mV,
                    color=COLORS[label], lw=2.5, ls='--' if label == 'Initial' else '-',
                    label=f'{label}  (kₙ = {kn*1e7:.2f} × 10⁻⁷)')
        if row == 0:
            ax.set_ylim(min_pot, max_pot)
            ax.set_title(f'{rate:g}C charge', fontsize=16, weight='bold')
        else:
            ax.set_ylim(-120, 200)
            ax.axhspan(-120, 0, color='#F7E4E0', alpha=.55, zorder=0)
            ax.axhline(0, color='#B24035', lw=1.15)
            ax.set_xlabel('Estimated full-cell SoC (%)', fontsize=12)
        ax.set_xlim(0, 100)
        ax.grid(alpha=.22)
        ax.tick_params(labelsize=10)
        if col == 0:
            ax.set_ylabel('Anode potential (mV)', fontsize=12)
axes[0, 0].text(.02, .94, 'Full potential range', transform=axes[0, 0].transAxes,
                ha='left', va='top', color='#344251', fontsize=10.5)
axes[1, 0].text(.02, .94, 'Near-zero detail', transform=axes[1, 0].transAxes,
                ha='left', va='top', color='#344251', fontsize=10.5)
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', ncol=2, bbox_to_anchor=(.5, .95),
           frameon=False, fontsize=12)
fig.text(.075, .987, 'Stage 2 intermediate model: estimated anode potential by SoC',
         fontsize=17, color='#156082', weight='bold', va='top')
fig.text(.075, .035, 'Separator-side negative electrode surface potential difference; Ai2020 electrolyte and bₙ = 2.914 fixed.',
         fontsize=10.8, color='#4A5560')
stem = OUT / 'Stage2_kn554_initial_vs_kn_only_fit_anode_potential_SoC'
fig.savefig(stem.with_suffix('.png'), dpi=200, facecolor='white')
fig.savefig(stem.with_suffix('.svg'), facecolor='white')
plt.close(fig)
print(pd.DataFrame(summary).to_string(index=False), flush=True)
print(stem.with_suffix('.png'), flush=True)
