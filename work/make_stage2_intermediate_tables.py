"""Slide-ready English error tables for the pre-fit Stage 2 model."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs'
raw = pd.read_csv(OUT / 'Stage2_two_OCP_branch_metrics.csv')
case = raw.loc[raw.OCP_case == 'Ai2020_anode__measured_cathode'].copy()
phases = [
    ('OCP_only', 'Stage 1 C: OCP + window'),
    ('plus_Ds_GITT', '+ GITT Dₛ'),
    ('plus_kp_EIS', '+ cathode EIS kₚ'),
    ('plus_kn_DRT_initial', '+ anode DRT kₙ (initial)'),
]
conditions = [(0.5, 'charge'), (0.5, 'discharge'), (1.0, 'charge'),
              (1.0, 'discharge'), (2.0, 'charge'), (2.0, 'discharge')]
plt.rcParams.update({'font.family': 'DejaVu Sans', 'axes.unicode_minus': False})


def rms(frame, col):
    return float(np.sqrt(np.mean(frame[col].to_numpy() ** 2)))


def draw_table(title, subtitle, headers, rows, out, foot, widths, height=5.3, rowh=.095):
    fig, ax = plt.subplots(figsize=(14.8, height), dpi=180)
    fig.patch.set_facecolor('white')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    ax.text(.035, .93, title, fontsize=21, weight='bold', color='#156082', va='center')
    ax.text(.035, .865, subtitle, fontsize=12.6, color='#344251', va='center')
    x = [.035]
    for w in widths:
        x.append(x[-1] + w)
    top = .77
    for j, header in enumerate(headers):
        ax.add_patch(Rectangle((x[j], top-rowh), x[j+1]-x[j], rowh, facecolor='#156082', edgecolor='white', linewidth=1))
        ax.text((x[j]+x[j+1])/2, top-rowh/2, header, ha='center', va='center', color='white', fontsize=11.7,
                weight='bold', linespacing=1.15)
    for i, row in enumerate(rows):
        y = top - (i+2)*rowh
        fill = '#C4DAEB' if i == len(rows)-1 else ('#DCE7F2' if i % 2 == 0 else '#F3F6FA')
        for j, val in enumerate(row):
            ax.add_patch(Rectangle((x[j], y), x[j+1]-x[j], rowh, facecolor=fill, edgecolor='white', linewidth=1))
            ax.text(x[j]+.015 if j == 0 else (x[j]+x[j+1])/2, y+rowh/2, val,
                    ha='left' if j == 0 else 'center', va='center', fontsize=12.5 if j == 0 else 14.5,
                    color='#202A35', weight='bold' if i == len(rows)-1 else 'normal')
    ax.text(.035, .055, foot, fontsize=11, color='#4A5560', va='center')
    fig.savefig(out, dpi=200, facecolor='white', bbox_inches='tight', pad_inches=.12)
    plt.close(fig)
    print(out)


summary = []
for phase, label in phases:
    frame = case.loc[case.Phase == phase]
    summary.append([label, f'{rms(frame, "full_RMSE_mV"):.2f}',
                    f'{rms(frame, "center_10_70_RMSE_mV"):.2f}'])
draw_table(
    'Stage 2 intermediate model: experimental Dₛ and k',
    'Ai2020 anode OCP + measured cathode OCP; stoichiometry window and cₛ,ₘₐₓ recalculated',
    ['Sequential parameter application', 'Full-range RMSE (mV)', 'Central 10–70% RMSE (mV)'],
    summary,
    OUT / 'Stage2_intermediate_sequential_error_table_English.png',
    'Six-curve RMS across 0.5C, 1C and 2C charge/discharge. kₙ is a provisional pre-fit input.',
    [.42, .25, .26], height=3.6, rowh=.13)

detail = []
selected = case.loc[case.Phase == 'plus_kn_DRT_initial'].copy()
for rate, direction in conditions:
    r = selected.loc[(selected.C_rate == rate) & (selected.Direction == direction)].iloc[0]
    detail.append([f'{rate:g}C {direction.title()}',
                   f'{r.full_RMSE_mV:.2f}', f'{r.full_MAE_mV:.2f}',
                   f'{r.center_10_70_RMSE_mV:.2f}', f'{r.center_10_70_MAE_mV:.2f}'])
detail.append(['Six-curve RMSE', f'{rms(selected, "full_RMSE_mV"):.2f}', '—',
               f'{rms(selected, "center_10_70_RMSE_mV"):.2f}', '—'])
draw_table(
    'Stage 2 intermediate model: voltage error by condition',
    'GITT Dₛ + cathode EIS kₚ + anode DRT kₙ (initial); before voltage fitting and electrolyte update',
    ['Condition', 'Full RMSE\n(mV)', 'Full MAE\n(mV)', 'Central RMSE\n(mV)', 'Central MAE\n(mV)'],
    detail,
    OUT / 'Stage2_intermediate_six_curve_error_table_English.png',
    'Central interval: 10–70% of transferred capacity in the compact notebook. Final-row RMSE: RMS of six conditions.',
    [.25, .17, .17, .17, .17], height=6.5, rowh=.078)

export = selected[['C_rate', 'Direction', 'full_RMSE_mV', 'full_MAE_mV',
                   'center_10_70_RMSE_mV', 'center_10_70_MAE_mV', 'capacity_error_pct']].copy()
export.to_csv(OUT / 'Stage2_intermediate_six_curve_error_table_English.csv', index=False)
