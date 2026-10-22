"""Render English versions of the existing Stage 1 DFN figure and error tables.

Uses the saved curve and metric CSVs, so no simulation or recalculation is needed.
"""
from pathlib import Path
import glob

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs'
curves = pd.read_csv(glob.glob(str(OUT / 'Stage1*3단계*동적전압_곡선.csv'))[0])
metrics = pd.read_csv(glob.glob(str(OUT / 'Stage1*3단계*동적전압_RMSE_MAE.csv'))[0])
stages = list(metrics.stage.drop_duplicates())
labels = {
    stages[0]: 'A  Ai2020 baseline window',
    stages[1]: 'B  Ai2020 OCP + refitted window',
    stages[2]: 'C  Measured cathode OCP + refitted window',
}
colors = {stages[0]: '#87909B', stages[1]: '#2371B6', stages[2]: '#D66422'}
styles = {stages[0]: '--', stages[1]: '-.', stages[2]: '-'}
plt.rcParams.update({'font.family': 'DejaVu Sans', 'axes.unicode_minus': False})


def make_figure():
    agg = {s: float(np.sqrt(np.mean(metrics.loc[metrics.stage == s, 'full_RMSE_mV'] ** 2))) for s in stages}
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.8), sharey=True, dpi=160)
    fig.subplots_adjust(left=.065, right=.985, top=.85, bottom=.17, wspace=.14, hspace=.24)
    for col, rate in enumerate((0.5, 1.0, 2.0)):
        for row, direction in enumerate(('Charge', 'Discharge')):
            ax = axes[row, col]
            exp = curves[(curves.stage == 'Experiment') & (curves.C_rate == rate) & (curves.direction == direction)]
            ax.plot(exp.time_min, exp.voltage_V, color='#20262D', lw=2.65, label='Experiment', zorder=5)
            for stage in stages:
                d = curves[(curves.stage == stage) & (curves.C_rate == rate) & (curves.direction == direction)]
                ax.plot(d.time_min, d.voltage_V, color=colors[stage], ls=styles[stage], lw=2.25, label=labels[stage])
            ax.set_title(f'{rate:g}C {direction}', fontsize=17, weight='bold', pad=9)
            ax.set_ylim(2.95, 4.27)
            ax.grid(alpha=.24)
            ax.set_axisbelow(True)
            ax.tick_params(labelsize=11)
            if row == 1:
                ax.set_xlabel('Elapsed time (min)', fontsize=13)
            if col == 0:
                ax.set_ylabel('Cell voltage (V)', fontsize=13)
    handles, legends = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, legends, loc='upper center', ncol=4, bbox_to_anchor=(.51, .945), frameon=False, fontsize=11)
    fig.text(.065, .975, r'Common inputs: Ai2020 $D_e\times10^{-4}$ · $r_n$ 3.692 µm · $r_p$ 4.239 µm · $D_s$, $k$, Bruggeman not fitted',
             fontsize=14, color='#156082', weight='bold', va='top')
    fig.text(.065, .087, 'Six-curve full-range RMSE (mV)   ' + '   →   '.join(f'{chr(65+i)} {agg[s]:.2f}' for i, s in enumerate(stages)),
             fontsize=17, color='#156082', weight='bold')
    fig.text(.065, .04, 'B retains the Ai2020 cathode OCP; C uses the measured cathode GITT OCP. The window and cₛ,ₘₐₓ are recalculated in each.',
             fontsize=11.5, color='#4A5560')
    stem = OUT / 'Stage1_dynamic_voltage_0p5C_1C_2C_English'
    for ext in ('.png', '.svg'):
        fig.savefig(stem.with_suffix(ext), dpi=200, facecolor='white')
        print(stem.with_suffix(ext))
    plt.close(fig)


def make_table(kind):
    if kind == 'full':
        title = 'Stage 1 dynamic voltage error · full common time interval'
        rmse, mae = 'full_RMSE_mV', 'full_MAE_mV'
        foot = '400 common elapsed-time points per condition. Final row: root mean square of six condition RMSEs.'
        out = OUT / 'Stage1_dynamic_voltage_full_error_table_English.png'
    else:
        title = 'Stage 1 dynamic voltage error · central 10–70% interval'
        rmse, mae = 'center_10_70_RMSE_mV', 'center_10_70_MAE_mV'
        foot = 'Central interval: 10–70% of transferred capacity in the compact notebook, not full-cell SOC. Final row: RMS of six RMSEs.'
        out = OUT / 'Stage1_dynamic_voltage_central_error_table_English.png'
    conditions = [(.5, 'Charge'), (.5, 'Discharge'), (1., 'Charge'), (1., 'Discharge'), (2., 'Charge'), (2., 'Discharge')]
    rows = []
    for rate, direction in conditions:
        vals = []
        for stage in stages:
            r = metrics[(metrics.stage == stage) & (metrics.C_rate == rate) & (metrics.direction == direction)].iloc[0]
            vals.append(f'{r[rmse]:.2f} / {r[mae]:.2f}')
        rows.append([f'{rate:g}C {direction}', *vals])
    overall = [float(np.sqrt(np.mean(metrics.loc[metrics.stage == stage, rmse] ** 2))) for stage in stages]
    rows.append(['Six-curve RMSE', *[f'{v:.2f}' for v in overall]])
    fig, ax = plt.subplots(figsize=(14.8, 6.5), dpi=180)
    fig.patch.set_facecolor('white')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    ax.text(.035, .925, title, fontsize=21, fontweight='bold', color='#156082', va='center')
    ax.text(.035, .872, 'Unit: mV  |  Values in each cell: RMSE / MAE', fontsize=15, color='#344251', va='center')
    x = [.035, .235, .49, .745, .965]
    top, rowh = .805, .091
    headers = ['Condition', 'A  Ai2020\nbaseline window', 'B  Ai2020 OCP\nrefitted window', 'C  Measured cathode OCP\nrefitted window']
    for j in range(4):
        ax.add_patch(Rectangle((x[j], top-rowh), x[j+1]-x[j], rowh, facecolor='#156082', edgecolor='white', linewidth=1))
        ax.text((x[j]+x[j+1])/2, top-rowh/2, headers[j], ha='center', va='center', color='white', fontsize=11.5, fontweight='bold', linespacing=1.15)
    for i, row in enumerate(rows):
        y = top-(i+2)*rowh
        fill = '#C4DAEB' if i == len(rows)-1 else ('#DCE7F2' if i % 2 == 0 else '#F3F6FA')
        for j, val in enumerate(row):
            ax.add_patch(Rectangle((x[j], y), x[j+1]-x[j], rowh, facecolor=fill, edgecolor='white', linewidth=1))
            ax.text((x[j]+x[j+1])/2, y+rowh/2, val, ha='center', va='center', fontsize=15 if j else 13.5,
                    color='#202A35', fontweight='bold' if i == len(rows)-1 else 'normal')
    ax.text(.035, .047, foot, fontsize=11.5, color='#4A5560', va='center')
    fig.savefig(out, dpi=200, facecolor='white', bbox_inches='tight', pad_inches=.12)
    plt.close(fig)
    print(out)


make_figure()
make_table('full')
make_table('central')
