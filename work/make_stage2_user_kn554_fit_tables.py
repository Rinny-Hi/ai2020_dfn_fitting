"""English tables for the user-selected Stage 2 k_n initial value and charge fit."""
from pathlib import Path
import shutil

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs'
earlier = pd.read_csv(OUT / 'Stage2_two_OCP_branch_metrics.csv')
earlier = earlier.loc[earlier.OCP_case == 'Ai2020_anode__measured_cathode']
new = pd.read_csv(OUT / 'Stage2_user_kn554_initial_vs_fitted_branch_metrics.csv')
fit = pd.read_csv(OUT / 'Stage2_user_kn554_charge_fit_result.csv').iloc[0]
plt.rcParams.update({'font.family': 'DejaVu Sans', 'axes.unicode_minus': False})


def rms(df, col):
    return float(np.sqrt(np.mean(np.square(df[col].to_numpy()))))


def table(title, subtitle, headers, rows, widths, foot, path, height, rowh, top=.78):
    fig, ax = plt.subplots(figsize=(14.8, height), dpi=180)
    fig.patch.set_facecolor('white')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    ax.text(.035, .93, title, fontsize=20.5, weight='bold', color='#156082', va='center')
    ax.text(.035, .865, subtitle, fontsize=12.5, color='#344251', va='center')
    x = [.035]
    for width in widths:
        x.append(x[-1] + width)
    for j, h in enumerate(headers):
        ax.add_patch(Rectangle((x[j], top-rowh), widths[j], rowh, facecolor='#156082', edgecolor='white'))
        ax.text((x[j]+x[j+1])/2, top-rowh/2, h, ha='center', va='center', color='white',
                fontsize=11.4, weight='bold', linespacing=1.1)
    for i, row in enumerate(rows):
        y = top - (i+2)*rowh
        fill = '#DCE7F2' if i % 2 == 0 else '#F3F6FA'
        for j, val in enumerate(row):
            ax.add_patch(Rectangle((x[j], y), widths[j], rowh, facecolor=fill, edgecolor='white'))
            ax.text(x[j]+.015 if j == 0 else (x[j]+x[j+1])/2, y+rowh/2, str(val),
                    ha='left' if j == 0 else 'center', va='center',
                    fontsize=12.5 if j == 0 else 14, color='#202A35',
                    weight='bold' if i == len(rows)-1 else 'normal')
    ax.text(.035, .048, foot, fontsize=10.8, color='#4A5560', va='center')
    fig.savefig(path, dpi=200, facecolor='white', bbox_inches='tight', pad_inches=.12)
    plt.close(fig)
    print(path)


initial = new.loc[new.Phase == 'Initial_kn554']
fitted = new.loc[new.Phase == 'Charge_fit_discharge_validation']
sequence = []
for phase, label in [('OCP_only', 'Stage 1 C: OCP + window'),
                     ('plus_Ds_GITT', '+ GITT Dₛ'),
                     ('plus_kp_EIS', '+ cathode EIS kₚ')]:
    d = earlier.loc[earlier.Phase == phase]
    sequence.append((label, d, 'Direction', 'charge', 'discharge'))
sequence += [('+ kₙ = 5.54 × 10⁻⁷ (initial)', initial, 'Direction', 'Charge', 'Discharge'),
             ('+ charge fit: kₙ only', fitted, 'Direction', 'Charge', 'Discharge')]

summary_rows, summary_csv = [], []
for label, df, direction_col, charge_label, discharge_label in sequence:
    ch = df.loc[df[direction_col] == charge_label]
    dis = df.loc[df[direction_col] == discharge_label]
    vals = [rms(ch, 'full_RMSE_mV'), rms(dis, 'full_RMSE_mV'),
            rms(df, 'full_RMSE_mV'), rms(df, 'center_10_70_RMSE_mV')]
    summary_rows.append([label, *[f'{v:.2f}' for v in vals]])
    summary_csv.append(dict(Stage=label, Charge_full_RMSE_mV=vals[0],
                            Discharge_full_RMSE_mV=vals[1],
                            Six_curve_full_RMSE_mV=vals[2], Six_curve_central_RMSE_mV=vals[3]))
pd.DataFrame(summary_csv).to_csv(OUT / 'Stage2_intermediate_kn554_to_charge_fit_summary_English.csv', index=False)

table(
    'Stage 2 intermediate model: parameter application and fitting',
    'Ai2020 anode OCP + measured cathode OCP; Ai2020 electrolyte retained',
    ['Model stage', 'Charge RMSE\n(mV)', 'Discharge RMSE\n(mV)', 'Six-curve full\nRMSE (mV)', 'Six-curve central\nRMSE (mV)'],
    summary_rows, [.35, .145, .145, .145, .145],
    'Only kₙ fitted to charge; bₙ = 2.914 fixed. Discharge independently validated. Central = transferred capacity 10–70%.',
    OUT / 'Stage2_intermediate_kn554_to_charge_fit_table_English.png', 5.0, .105)

def detailed(df, title, subtitle, path, foot, show_roles):
    rows = []
    for rate in (.5, 1., 2.):
        for direction in ('Charge', 'Discharge'):
            r = df.loc[(df.C_rate == rate) & (df.Direction == direction)].iloc[0]
            role = (' (fit)' if direction == 'Charge' else ' (validation)') if show_roles else ''
            rows.append([f'{rate:g}C {direction}{role}', f'{r.full_RMSE_mV:.2f}',
                         f'{r.full_MAE_mV:.2f}', f'{r.center_10_70_RMSE_mV:.2f}',
                         f'{r.center_10_70_MAE_mV:.2f}'])
    rows.append(['Six-curve RMSE', f'{rms(df, "full_RMSE_mV"):.2f}', '—',
                 f'{rms(df, "center_10_70_RMSE_mV"):.2f}', '—'])
    table(title, subtitle,
          ['Condition', 'Full RMSE\n(mV)', 'Full MAE\n(mV)', 'Central RMSE\n(mV)', 'Central MAE\n(mV)'],
          rows, [.29, .16, .16, .16, .16], foot, path, 6.5, .078)


detailed(initial, 'Stage 2 intermediate model: kₙ initial value',
         'GITT Dₛ + cathode EIS kₚ = 4.768 × 10⁻⁷ + anode kₙ = 5.54 × 10⁻⁷; bₙ = 2.914',
         OUT / 'Stage2_intermediate_kn554_initial_error_table_English.png',
         'Central = 10–70% of transferred capacity. No voltage fitting or electrolyte update at this stage.', False)

detailed(fitted, 'Stage 2 intermediate model: kₙ fit and discharge validation',
         f'Fitted kₙ = {fit.fitted_kn*1e7:.3f} × 10⁻⁷; bₙ = {fit.fitted_brugg_n:.3f} fixed; cathode kₚ = 4.768 × 10⁻⁷ fixed',
         OUT / 'Stage2_intermediate_kn554_charge_fit_error_table_English.png',
         '0.5C, 1C and 2C charge used for fitting; discharge is independent validation. Central = transferred capacity 10–70%.', True)

# Replace the previously shared Stage 2 initial tables so their values no longer
# conflict with the user's revised initial k_n.
shutil.copyfile(OUT / 'Stage2_intermediate_kn554_initial_error_table_English.png',
                OUT / 'Stage2_intermediate_six_curve_error_table_English.png')
shutil.copyfile(OUT / 'Stage2_intermediate_kn554_to_charge_fit_table_English.png',
                OUT / 'Stage2_intermediate_sequential_error_table_English.png')
initial[['C_rate', 'Direction', 'full_RMSE_mV', 'full_MAE_mV',
         'center_10_70_RMSE_mV', 'center_10_70_MAE_mV', 'capacity_error_pct']].to_csv(
             OUT / 'Stage2_intermediate_six_curve_error_table_English.csv', index=False)
