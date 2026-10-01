"""Analyze two 2026-09-27 RPT workbooks and validate six geometry cases."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import c50_selected_ocp_dynamic_validation as selected
import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import kp_eis_recalculation_ablation as kpstudy
import priority_initial_state_protocol_recheck as priority
import thickness_radius_six_case_comparison as six
import three_thickness_measurement_comparison as thickness3


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260927_two_cell_rpt_validation"
RESULTS.mkdir(parents=True, exist_ok=True)

FILES = {
    "6-1": Path(r"E:\예린\260927 6-1 0.5C 1C 2C.xlsx"),
    "6-2": Path(r"E:\예린\260927 6-2 0.5C 1C 2C.xlsx"),
}

# The first partial charge is conditioning.  The paired full branches are below.
BRANCH_STEPS = {
    (0.5, False): 5,
    (0.5, True): 7,
    (1.0, False): 10,
    (1.0, True): 12,
    (2.0, False): 15,
    (2.0, True): 17,
}

NOMINAL_CAPACITY_AH = 2.28


def parse_workbook(cell: str, path: Path):
    step = pd.read_excel(path, sheet_name="step")
    record = pd.read_excel(path, sheet_name="record")
    for col in ("Oneset Date", "End Date"):
        step[col] = pd.to_datetime(step[col])
    record["Date"] = pd.to_datetime(record["Date"])
    branches = {}
    rows = []
    step_by_number = step.set_index("Step Number")
    for (rate, charge), step_number in BRANCH_STEPS.items():
        direction = "Charge" if charge else "Discharge"
        branch = step_by_number.loc[step_number]
        previous = step_by_number.loc[step_number - 1]
        seg = record[
            (record.Date >= branch["Oneset Date"])
            & (record.Date <= branch["End Date"])
            & (record["Step Type"] == branch["Step Type"])
        ].copy()
        if seg.empty:
            raise RuntimeError(f"No record data for {cell}, step {step_number}")
        seg["t_min"] = (seg.Date - seg.Date.iloc[0]).dt.total_seconds() / 60.0
        q = np.asarray(seg["Capacity(Ah)"], float)
        q = np.abs(q - q[0])
        seg["Q_Ah"] = q
        seg = seg.rename(columns={"Voltage(V)": "V", "Current(A)": "I_A"})
        branches[(rate, charge)] = seg[["t_min", "V", "I_A", "Q_Ah"]].copy()

        rest = record[
            (record.Date >= previous["Oneset Date"])
            & (record.Date <= previous["End Date"])
            & (record["Step Type"] == "Rest")
        ].copy()
        rest["t_min"] = (rest.Date - rest.Date.iloc[0]).dt.total_seconds() / 60.0
        tail = rest[rest.t_min >= rest.t_min.max() - 10.0]
        rest_slope = float(np.polyfit(tail.t_min, tail["Voltage(V)"], 1)[0] * 1000.0)
        current = float(np.nanmedian(np.abs(seg.I_A)))
        rows.append(
            {
                "Cell": cell,
                "C_rate_label": rate,
                "Direction": direction,
                "Step_number": step_number,
                "Step_type": str(branch["Step Type"]),
                "Median_current_A": current,
                "Actual_C_rate_vs_2.28Ah": current / NOMINAL_CAPACITY_AH,
                "CC_duration_min": float(seg.t_min.iloc[-1]),
                "CC_capacity_Ah_record": float(branch["Capacity(Ah)"]),
                "CC_capacity_Ah_extracted": float(seg.Q_Ah.iloc[-1]),
                "Onset_V": float(branch["Oneset Volt.(V)"]),
                "End_V": float(branch["End Voltage(V)"]),
                "Preceding_rest_duration_min": float((previous["End Date"] - previous["Oneset Date"]).total_seconds() / 60.0),
                "Rest_start_V": float(rest["Voltage(V)"].iloc[0]),
                "Rest_end_V": float(rest["Voltage(V)"].iloc[-1]),
                "Rest_relaxation_mV": 1000.0 * float(rest["Voltage(V)"].iloc[-1] - rest["Voltage(V)"].iloc[0]),
                "Rest_last10min_slope_mV_per_min": rest_slope,
            }
        )
    return branches, pd.DataFrame(rows), step


def capacity_repeatability(qc: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (rate, direction), group in qc.groupby(["C_rate_label", "Direction"], sort=True):
        q = np.asarray(group.CC_capacity_Ah_extracted, float)
        rows.append(
            {
                "C_rate": rate,
                "Direction": direction,
                "Mean_CC_capacity_Ah": float(np.mean(q)),
                "Difference_mAh_6_2_minus_6_1": float(1000.0 * (q[1] - q[0])),
                "Range_mAh": float(1000.0 * np.ptp(q)),
                "Range_pct_of_mean": float(100.0 * np.ptp(q) / np.mean(q)),
            }
        )
    return pd.DataFrame(rows)


def build_models(model0: dict) -> dict:
    models = {}
    for thickness_case, lv in thickness3.CASES_UM.items():
        geometry = thickness3.same_loading_geometry(model0, lv["Ln"], lv["Lp"])
        for radius_case, rv in six.RADIUS_CASES_UM.items():
            models[(thickness_case, radius_case)] = six.with_radius(
                geometry, rv["Rn"], rv["Rp"]
            )
    return models


def diffusion_timescales(qc: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    dsn, dsp = 2.1e-14, 4.4e-14
    rows = []
    for radius_case, rv in six.RADIUS_CASES_UM.items():
        for electrode, r_um, ds in (
            ("Negative", rv["Rn"], dsn),
            ("Positive", rv["Rp"], dsp),
        ):
            tau = (r_um * 1e-6) ** 2 / ds
            rows.append(
                {
                    "Radius_case": radius_case,
                    "Electrode": electrode,
                    "Radius_um": r_um,
                    "Ds_m2_s": ds,
                    "R2_over_D_s": tau,
                    "R2_over_D_min": tau / 60.0,
                    "first_mode_R2_over_pi2D_s": tau / np.pi**2,
                    "first_mode_R2_over_pi2D_min": tau / np.pi**2 / 60.0,
                    "Fourier_number_for_10min_GITT_pulse": 600.0 / tau,
                }
            )
    times = pd.DataFrame(rows)
    duration = qc.groupby(["C_rate_label", "Direction"], as_index=False).CC_duration_min.mean()
    comparisons = []
    for _, branch in duration.iterrows():
        for _, ts in times.iterrows():
            comparisons.append(
                {
                    "C_rate": branch.C_rate_label,
                    "Direction": branch.Direction,
                    "Mean_CC_duration_min": branch.CC_duration_min,
                    "Radius_case": ts.Radius_case,
                    "Electrode": ts.Electrode,
                    "R2_over_D_min": ts.R2_over_D_min,
                    "CC_duration_over_R2D": branch.CC_duration_min / ts.R2_over_D_min,
                }
            )
    return times, pd.DataFrame(comparisons)


def main() -> None:
    all_branches = {}
    qc_frames = []
    step_tables = {}
    for cell, path in FILES.items():
        branches, qc, step = parse_workbook(cell, path)
        all_branches[cell] = branches
        qc_frames.append(qc)
        step_tables[cell] = step
    qc = pd.concat(qc_frames, ignore_index=True)
    qc.to_csv(RESULTS / "new_rpt_branch_qc.csv", index=False, encoding="utf-8-sig")
    repeat = capacity_repeatability(qc)
    repeat.to_csv(RESULTS / "two_cell_repeatability.csv", index=False, encoding="utf-8-sig")
    timescales, time_comparison = diffusion_timescales(qc)
    timescales.to_csv(RESULTS / "diffusion_timescales.csv", index=False, encoding="utf-8-sig")
    time_comparison.to_csv(RESULTS / "branch_duration_vs_diffusion_time.csv", index=False, encoding="utf-8-sig")

    model0, stage, anode, cathode, _, _, qcell = selected.selected_inputs()
    model0 = kpstudy.with_kp(model0, kpstudy.KP_EIS)
    models = build_models(model0)

    score_rows = []
    simulations = {}
    for (thickness_case, radius_case), model in models.items():
        for cell, branches in all_branches.items():
            local_qc = qc[qc.Cell == cell].set_index(["C_rate_label", "Direction"])
            for rate in ga.RATES:
                for charge in (True, False):
                    direction = "Charge" if charge else "Discharge"
                    exp = branches[(rate, charge)]
                    rest_v = float(local_qc.loc[(rate, direction), "Rest_end_V"])
                    z0 = priority.monotone_voltage_inverse(
                        model["v_qocv"], model["soc"], rest_v
                    )
                    local_stage = priority.stage_at_soc(stage, z0, charge)
                    current = float(np.nanmedian(np.abs(exp.I_A)))
                    sim_rate = current / NOMINAL_CAPACITY_AH
                    print(
                        f"{thickness_case} | {radius_case} | {cell} | {rate:g}C {direction}",
                        flush=True,
                    )
                    sim = ehq.run_dfn(
                        sim_rate,
                        charge,
                        local_stage,
                        model,
                        anode,
                        cathode,
                        False,
                        True,
                        positive_initial_h=(-1.0 if charge else 1.0),
                    )
                    simulations[(thickness_case, radius_case, cell, rate, charge)] = sim
                    metric = priority.curve_metrics(exp, sim, qcell, charge)
                    score_rows.append(
                        {
                            "Thickness_case": thickness_case,
                            "Radius_case": radius_case,
                            "Cell": cell,
                            "C_rate": rate,
                            "Direction": direction,
                            "Measured_current_A": current,
                            "Initial_SOC_from_rest": z0,
                            **metric,
                            "Capacity_error_mAh": 1000.0
                            * (metric["Q_end_model_Ah"] - metric["Q_end_exp_Ah"]),
                        }
                    )
    scores = pd.DataFrame(score_rows)
    scores.to_csv(RESULTS / "six_case_new_data_detail.csv", index=False, encoding="utf-8-sig")
    summary = scores.groupby(
        ["Thickness_case", "Radius_case", "Direction"], sort=False, as_index=False
    ).agg(
        Mean_full_RMSE_mV=("Full_RMSE_mV", "mean"),
        Mean_center_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
        Capacity_RMSE_pct=("Capacity_error_pct", lambda x: float(np.sqrt(np.mean(x**2)))),
        Mean_capacity_error_mAh=("Capacity_error_mAh", "mean"),
        Mean_abs_capacity_error_mAh=("Capacity_error_mAh", lambda x: float(np.mean(np.abs(x)))),
        Max_abs_capacity_error_mAh=("Capacity_error_mAh", lambda x: float(np.max(np.abs(x)))),
    )
    summary.to_csv(RESULTS / "six_case_new_data_summary.csv", index=False, encoding="utf-8-sig")
    overall = scores.groupby(
        ["Thickness_case", "Radius_case"], sort=False, as_index=False
    ).agg(
        Mean_full_RMSE_mV=("Full_RMSE_mV", "mean"),
        Mean_center_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
        Capacity_RMSE_pct=("Capacity_error_pct", lambda x: float(np.sqrt(np.mean(x**2)))),
        Mean_abs_capacity_error_mAh=("Capacity_error_mAh", lambda x: float(np.mean(np.abs(x)))),
        Max_abs_capacity_error_mAh=("Capacity_error_mAh", lambda x: float(np.max(np.abs(x)))),
    ).sort_values(["Mean_full_RMSE_mV", "Capacity_RMSE_pct"])
    overall.to_csv(RESULTS / "six_case_new_data_overall.csv", index=False, encoding="utf-8-sig")

    # Plot both cells and the two Pareto candidates retained from the old-data study.
    candidates = [
        ("Ai2020 paper", "SEM ImageJ R"),
        ("Thickness gauge", "SEM ImageJ R"),
    ]
    colors = {"6-1": "#111111", "6-2": "#777777"}
    model_colors = {
        ("Ai2020 paper", "SEM ImageJ R"): "#0072B2",
        ("Thickness gauge", "SEM ImageJ R"): "#D55E00",
    }
    fig, axes = plt.subplots(2, 3, figsize=(16.8, 9.0), constrained_layout=True)
    for col, rate in enumerate(ga.RATES):
        for row_idx, charge in enumerate((True, False)):
            ax = axes[row_idx, col]
            direction = "Charge" if charge else "Discharge"
            for cell in FILES:
                exp = all_branches[cell][(rate, charge)]
                ax.plot(exp.Q_Ah, exp.V, color=colors[cell], lw=2.2, label=f"Experiment {cell}")
            for candidate in candidates:
                # Average the two simulations on a common capacity grid for display.
                sims = [simulations[(*candidate, cell, rate, charge)] for cell in FILES]
                qmax = min(
                    float(np.asarray(s["Q_Ah"], dtype=float)[-1])
                    if "Q_Ah" in s
                    else float(priority.transferred_capacity(s, qcell, charge)[-1])
                    for s in sims
                )
                grid = np.linspace(0.0, qmax, 500)
                values = []
                for sim in sims:
                    qsim = priority.transferred_capacity(sim, qcell, charge)
                    values.append(np.interp(grid, qsim, np.asarray(sim["V"], dtype=float)))
                metric = scores[
                    (scores.Thickness_case == candidate[0])
                    & (scores.Radius_case == candidate[1])
                    & np.isclose(scores.C_rate, rate)
                    & (scores.Direction == direction)
                ]
                ax.plot(
                    grid,
                    np.mean(values, axis=0),
                    color=model_colors[candidate],
                    lw=1.9,
                    label=f"{candidate[0]} + SEM R: {metric.Full_RMSE_mV.mean():.1f} mV",
                )
            ax.set(title=f"{rate:g}C {direction}", xlabel="CC capacity [Ah]", ylabel="Voltage [V]")
            ax.grid(alpha=0.22)
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("Independent validation on cells 6-1 and 6-2")
    fig.savefig(RESULTS / "new_rpt_pareto_candidate_curves.png", dpi=220)
    plt.close(fig)

    # Raw experimental repeatability plot.
    fig, axes = plt.subplots(2, 3, figsize=(16.5, 8.7), constrained_layout=True)
    for col, rate in enumerate(ga.RATES):
        for row_idx, charge in enumerate((True, False)):
            ax = axes[row_idx, col]
            for cell in FILES:
                exp = all_branches[cell][(rate, charge)]
                ax.plot(exp.Q_Ah, exp.V, lw=2.0, label=cell)
            ax.set(title=f"{rate:g}C {'Charge' if charge else 'Discharge'}", xlabel="CC capacity [Ah]", ylabel="Voltage [V]")
            ax.grid(alpha=0.22)
    axes[0, 0].legend()
    fig.suptitle("Two-cell RPT repeatability")
    fig.savefig(RESULTS / "new_rpt_two_cell_repeatability.png", dpi=220)
    plt.close(fig)

    report = {
        "files": {key: str(value) for key, value in FILES.items()},
        "branch_steps": {f"{rate:g}C {'Charge' if charge else 'Discharge'}": step for (rate, charge), step in BRANCH_STEPS.items()},
        "protocol_note": "Step 2-3 is a partial conditioning charge and is not used as the 0.5C validation charge.",
        "fixed_model": {
            "Dsn_m2_s": 2.1e-14,
            "Dsp_m2_s": 4.4e-14,
            "kp": kpstudy.KP_EIS,
            "kn": "legacy 7.40e-7 placeholder",
            "OCP": "current C/50-selected OCP/window",
        },
        "repeatability": repeat.to_dict(orient="records"),
        "diffusion_timescales": timescales.to_dict(orient="records"),
        "model_overall": overall.to_dict(orient="records"),
        "model_direction_summary": summary.to_dict(orient="records"),
    }
    (RESULTS / "new_rpt_validation_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("\nQC")
    print(qc.to_string(index=False))
    print("\nREPEATABILITY")
    print(repeat.to_string(index=False))
    print("\nDIFFUSION TIMESCALES")
    print(timescales.to_string(index=False))
    print("\nMODEL OVERALL")
    print(overall.to_string(index=False))
    print("\nMODEL BY DIRECTION")
    print(summary.to_string(index=False))
    print("\nSaved to", RESULTS)


if __name__ == "__main__":
    main()
