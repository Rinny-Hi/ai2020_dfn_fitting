"""Experimental-parameter DFN voltage and negative-interface potential evidence."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from run_stage2_two_ocp import cases, kinetic_constants, scope, stage_with_values


OUT = Path(__file__).resolve().parents[1] / "outputs"
CONDITIONS = tuple((rate, charge) for rate in (0.5, 1.0, 2.0) for charge in (True, False))
rows, potential_rows, curves = [], [], []
for case_name, case in cases().items():
    stage = stage_with_values(case, kinetic_constants(case))
    for rate, charge in CONDITIONS:
        run = scope["run_dfn_branch"](stage, rate, charge)
        exp = scope["experimental_curves"][(rate, charge)]
        metrics = scope["compare_with_experiment"](exp, run)
        direction = "Charge" if charge else "Discharge"
        rows.append({"OCP_case": case_name, "C_rate": rate, "Direction": direction, **metrics})
        time = np.asarray(run["time_min"], float)
        voltage = np.asarray(run["voltage_v"], float)
        potential = np.asarray(run["anode_potential_v"], float)
        soc = np.asarray(run["full_cell_soc_pct"], float)
        center = (soc >= 10) & (soc <= 70)
        at70 = float(np.interp(70, soc, potential)) if charge and soc[0] <= 70 <= soc[-1] else np.nan
        potential_rows.append({
            "OCP_case": case_name, "C_rate": rate, "Direction": direction,
            "Initial_SOC_pct": float(run["initial_soc_pct"]),
            "Potential_min_full_mV": 1000 * float(np.min(potential)),
            "Potential_min_SOC10_70_mV": 1000 * float(np.min(potential[center])) if center.any() else np.nan,
            "Potential_at_SOC70_mV": 1000 * at70 if np.isfinite(at70) else np.nan,
            "Potential_final_mV": 1000 * float(potential[-1]),
            "Time_below_0V_min": float(np.trapezoid((potential < 0).astype(float), time)),
        })
        for i in range(len(time)):
            curves.append({"OCP_case": case_name, "C_rate": rate, "Direction": direction,
                           "Series": "DFN", "Time_min": float(time[i]), "Voltage_V": float(voltage[i]),
                           "SOC_pct": float(soc[i]), "Anode_potential_V": float(potential[i])})
        for t, v in zip(exp["time_min"] - exp["time_min"][0], exp["voltage_v"]):
            curves.append({"OCP_case": case_name, "C_rate": rate, "Direction": direction,
                           "Series": "Experiment", "Time_min": float(t), "Voltage_V": float(v),
                           "SOC_pct": np.nan, "Anode_potential_V": np.nan})

metrics = pd.DataFrame(rows)
potentials = pd.DataFrame(potential_rows)
curves = pd.DataFrame(curves)
metrics.to_csv(OUT / "Stage3_prefit_two_OCP_voltage_metrics.csv", index=False, encoding="utf-8-sig")
potentials.to_csv(OUT / "Stage3_prefit_two_OCP_anode_potential.csv", index=False, encoding="utf-8-sig")
curves.to_csv(OUT / "Stage3_prefit_two_OCP_curves.csv", index=False, encoding="utf-8-sig")

for case_name, title in (
    ("Ai2020_anode__measured_cathode", "Ai2020 anode + measured cathode"),
    ("measured_anode__measured_cathode", "Measured anode + measured cathode"),
):
    fig, axes = plt.subplots(2, 3, figsize=(15.2, 8.0), sharey=True, layout="constrained")
    for col, rate in enumerate((0.5, 1.0, 2.0)):
        for row, direction in enumerate(("Charge", "Discharge")):
            ax = axes[row, col]
            subset = curves.loc[(curves.OCP_case == case_name) & (curves.C_rate == rate) & (curves.Direction == direction)]
            for series, color, style in (("Experiment", "#202020", "-"), ("DFN", "#2458a8", "--")):
                data = subset.loc[subset.Series == series]
                ax.plot(data.Time_min, data.Voltage_V, color=color, linestyle=style, lw=2, label=series)
            met = metrics.loc[(metrics.OCP_case == case_name) & (metrics.C_rate == rate) & (metrics.Direction == direction)].iloc[0]
            ax.set_title(f"{rate:g}C {direction}  |  RMSE {met.full_RMSE_mV:.1f} mV", fontsize=11)
            ax.grid(alpha=.23)
            if row == 1:
                ax.set_xlabel("Time (min)")
            if col == 0:
                ax.set_ylabel("Voltage (V)")
    axes[0, 0].legend(frameon=False)
    fig.suptitle(title + " | experimental Ds and kinetics, before voltage fitting", fontsize=15)
    fig.savefig(OUT / f"Stage3_prefit_{case_name}_six_curves.png", dpi=180)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(9, 5.2), layout="constrained")
    for rate, color in ((0.5, "#2458a8"), (1.0, "#cf6028"), (2.0, "#6b49a9")):
        data = curves.loc[(curves.OCP_case == case_name) & (curves.C_rate == rate) & (curves.Direction == "Charge") & (curves.Series == "DFN")]
        ax.plot(data.SOC_pct, 1000 * data.Anode_potential_V, color=color, lw=2, label=f"{rate:g}C")
    ax.axhline(0, color="black", ls="--", lw=1)
    ax.set(xlabel="Full-cell SOC (%)", ylabel="Negative electrode surface potential at separator interface (mV)",
           title=title + " | charge potential before voltage fitting")
    ax.grid(alpha=.23)
    ax.legend(frameon=False)
    fig.savefig(OUT / f"Stage3_prefit_{case_name}_anode_potential.png", dpi=180)
    plt.close(fig)

print(metrics.to_string(index=False))
print(potentials.to_string(index=False))
