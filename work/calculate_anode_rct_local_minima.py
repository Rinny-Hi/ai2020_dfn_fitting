"""Calculate anode N3+N4 DRT resistance from adjacent local minima.

Source: Ai2020_Anode_DRT_Frequency_Origin_Paste.xlsx, Origin_Paste.
The workbook ordinates are already normalized to the 10 mm disk area (ohm cm2).
"""

from pathlib import Path
import csv
import math

import numpy as np
import openpyxl
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt


SOURCE = Path(r"C:\Users\user\Downloads\Ai2020_Anode_DRT_Frequency_Origin_Paste.xlsx")
OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(exist_ok=True)

R_GAS = 8.314462618
FARADAY = 96485.33212
TEMP_K = 298.15
EPS_S_N = 0.610
L_N_M = 76.5e-6
R_N_M = 3.692e-6
C_E = 1000.0
X_0 = 0.003806
X_100 = 0.954221
Q_CELL_AH = (2.355758 + 2.357298) / 2.0
A_M2 = 787.236e-4
C_S_MAX_N = 3600 * Q_CELL_AH / (
    (X_100 - X_0) * FARADAY * A_M2 * L_N_M * EPS_S_N
)  # A2 slide: Ai2020 anode + measured cathode OCP
MEASURED_C_S_MAX_N = 32447.1  # separate measured-anode OCP case
MEASURED_X_0 = 0.002150
MEASURED_X_100 = 0.739784


def local_minimum(frequency, smoothed_gamma, left, right):
    candidates = np.flatnonzero((frequency >= left) & (frequency <= right))
    if not len(candidates):
        raise ValueError(f"No points within {left}–{right} Hz")
    index = candidates[np.argmin(smoothed_gamma[candidates])]
    if index in (candidates[0], candidates[-1]):
        raise ValueError(f"Window-edge minimum at {frequency[index]:g} Hz")
    return int(index)


workbook = openpyxl.load_workbook(SOURCE, read_only=True, data_only=True)
rows = list(workbook["Origin_Paste"].values)
results = []
curves = {}
for column in range(0, 18, 2):
    soc = int(str(rows[0][column]).split("%")[0])
    data = np.array(
        [
            (row[column], row[column + 1])
            for row in rows[3:]
            if isinstance(row[column], (float, int))
            and isinstance(row[column + 1], (float, int))
        ],
        dtype=float,
    )
    data = data[np.argsort(data[:, 0])]
    frequency, gamma = data.T
    smoothed = savgol_filter(gamma, 21, 3)
    low = local_minimum(frequency, smoothed, 0.15, 0.45)  # N4 / LF valley
    middle = local_minimum(frequency, smoothed, 3.0, 8.0)  # N4 / N3 valley
    high = local_minimum(frequency, smoothed, 50.0, 120.0)  # N3 / N2 valley
    if not low < middle < high:
        raise ValueError(f"Unexpected valley order at {soc}% SoC")

    # f = 1/(2*pi*tau), so |d ln(tau)| = d ln(f).
    r_n4 = float(np.trapezoid(gamma[low : middle + 1], np.log(frequency[low : middle + 1])))
    r_n3 = float(np.trapezoid(gamma[middle : high + 1], np.log(frequency[middle : high + 1])))
    r_ct = r_n3 + r_n4
    a_s = 3 * EPS_S_N / R_N_M
    j_0 = R_GAS * TEMP_K / (FARADAY * (r_ct * 1e-4) * a_s * L_N_M)
    x = X_0 + (X_100 - X_0) * soc / 100.0
    k_n = j_0 / (C_S_MAX_N * math.sqrt(C_E * x * (1 - x)))
    measured_x = MEASURED_X_0 + (MEASURED_X_100 - MEASURED_X_0) * soc / 100.0
    measured_k = j_0 / (
        MEASURED_C_S_MAX_N * math.sqrt(C_E * measured_x * (1 - measured_x))
    )
    results.append(
        {
            "soc_pct": soc,
            "f_low_min_Hz": frequency[low],
            "f_mid_min_Hz": frequency[middle],
            "f_high_min_Hz": frequency[high],
            "R_N4_ohm_cm2": r_n4,
            "R_N3_ohm_cm2": r_n3,
            "R_ct_N3_N4_ohm_cm2": r_ct,
            "x_n": x,
            "j0_A_m2": j_0,
            "k_n_A_m-2_(m3_mol-1)^1.5": k_n,
            "measured_anode_x_n": measured_x,
            "measured_anode_k_n_A_m-2_(m3_mol-1)^1.5": measured_k,
        }
    )
    curves[soc] = (frequency, gamma, smoothed, (low, middle, high))

csv_path = OUT / "Anode_DRT_local_minimum_Rct_results.csv"
with csv_path.open("w", newline="", encoding="utf-8-sig") as file:
    writer = csv.DictWriter(file, fieldnames=list(results[0]))
    writer.writeheader()
    writer.writerows(results)

for row in results:
    print(
        f"SoC {row['soc_pct']:2d}%: "
        f"min {row['f_low_min_Hz']:.4f}/{row['f_mid_min_Hz']:.3f}/{row['f_high_min_Hz']:.2f} Hz, "
        f"N4 {row['R_N4_ohm_cm2']:.4f}, N3 {row['R_N3_ohm_cm2']:.4f}, "
        f"Rct {row['R_ct_N3_N4_ohm_cm2']:.4f} ohm cm2, "
        f"j0 {row['j0_A_m2']:.4f} A/m2, k {row['k_n_A_m-2_(m3_mol-1)^1.5']:.4e}"
    )

print(f"Saved {csv_path}")

frequency, gamma, smoothed, (low, middle, high) = curves[50]
fig, ax = plt.subplots(figsize=(8, 4.5), layout="constrained")
ax.semilogx(frequency, gamma, color="#2458b4", lw=1.5, label="DRT, 50% SoC")
ax.fill_between(
    frequency[low : high + 1], gamma[low : high + 1],
    color="#80b4ee", alpha=0.35, label="N3 + N4 integral"
)
for index, label in ((low, "N4/LF min"), (middle, "N3/N4 min"), (high, "N2/N3 min")):
    ax.plot(frequency[index], gamma[index], "o", color="#c43132", ms=5)
    ax.annotate(label, (frequency[index], gamma[index]), xytext=(4, 8),
                textcoords="offset points", fontsize=8)
ax.set_xlim(0.08, 300)
ax.set_ylim(bottom=0)
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel(r"DRT $\gamma$ ($\Omega$ cm$^2$)")
ax.set_title("Anode N3 + N4 local-minimum integration, 50% SoC")
ax.grid(True, which="both", alpha=0.2)
ax.legend(loc="upper left", fontsize=8)
plot_path = OUT / "Anode_DRT_local_minimum_50SOC.png"
fig.savefig(plot_path, dpi=180)
plt.close(fig)
print(f"Saved {plot_path}")
