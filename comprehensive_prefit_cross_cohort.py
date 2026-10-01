"""Cross-cohort pre-fit audit for the Enertech 0.5C/1C/2C data.

The script keeps source workbooks read-only.  It compares the three-cell July
BoL RPT with the two-cell September RPT, audits normalization, and evaluates a
small set of physically interpretable fixed-parameter DFN candidates before
any new optimization is performed.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import c50_selected_ocp_dynamic_validation as selected
import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import kinetics_original_reversion_ablation as kinetics
import kp_eis_recalculation_ablation as kpstudy
import nominal_radius_current_configuration as current
import priority_initial_state_protocol_recheck as priority
import thickness_radius_six_case_comparison as radius_study
import three_thickness_measurement_comparison as thickness_study


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260927_comprehensive_prefit_cross_cohort"
RESULTS.mkdir(parents=True, exist_ok=True)

DATA_ROOT = Path(
    r"C:\Users\user\OneDrive\01. 적응형 충전 프로토콜"
    r"\4. Enertech 파우치셀 실험데이터"
)

OLD_FILES = {
    "6-5": DATA_ROOT / "260703 실험데이터" / "260703 BoL 데이터 6-5.xlsx",
    "6-6": DATA_ROOT / "260703 실험데이터" / "260703 BoL 데이터 6-6.xlsx",
    "6-8": DATA_ROOT / "260703 실험데이터" / "260703 BoL 데이터 6-8.xlsx",
}
NEW_FILES = {
    "6-1": Path(r"E:\예린\260927 6-1 0.5C 1C 2C.xlsx"),
    "6-2": Path(r"E:\예린\260927 6-2 0.5C 1C 2C.xlsx"),
}

OLD_BRANCH_STEPS = {
    (0.5, True): 7,
    (0.5, False): 10,
    (1.0, True): 12,
    (1.0, False): 15,
    (2.0, True): 17,
    (2.0, False): 20,
}
NEW_BRANCH_STEPS = {
    (0.5, False): 5,
    (0.5, True): 7,
    (1.0, False): 10,
    (1.0, True): 12,
    (2.0, False): 15,
    (2.0, True): 17,
}
NOMINAL_CAPACITY_AH = 2.28
AI2020_PREF = 1.0e-11 * 96485.33212


@dataclass(frozen=True)
class Candidate:
    name: str
    ln_um: float
    lp_um: float
    rn_um: float
    rp_um: float
    kn: float
    kp: float
    positive_hysteresis: bool
    initialization: str = "rest"
    physical_note: str = ""


def _duration_min(row: pd.Series) -> float:
    return float((pd.to_datetime(row["End Date"]) - pd.to_datetime(row["Oneset Date"])).total_seconds() / 60.0)


def parse_cell(
    cohort: str,
    cell: str,
    path: Path,
    branch_steps: dict[tuple[float, bool], int],
) -> tuple[dict[tuple[float, bool], dict[str, np.ndarray]], list[dict]]:
    step = pd.read_excel(path, sheet_name="step")
    record = pd.read_excel(path, sheet_name="record")
    for col in ("Oneset Date", "End Date"):
        step[col] = pd.to_datetime(step[col])
    record["Date"] = pd.to_datetime(record["Date"])
    by_number = step.set_index("Step Number")
    curves: dict[tuple[float, bool], dict[str, np.ndarray]] = {}
    rows: list[dict] = []
    for (rate, charge), step_no in branch_steps.items():
        branch = by_number.loc[step_no]
        previous = by_number.loc[step_no - 1]
        seg = record[
            (record.Date >= branch["Oneset Date"])
            & (record.Date <= branch["End Date"])
            & (record["Step Type"].astype(str) == str(branch["Step Type"]))
        ].copy()
        if seg.empty:
            raise RuntimeError(f"Missing branch {cohort} {cell} step {step_no}")
        seg["t_min"] = (seg.Date - seg.Date.iloc[0]).dt.total_seconds() / 60.0
        q = np.abs(pd.to_numeric(seg["Capacity(Ah)"], errors="coerce").to_numpy(float))
        q = np.abs(q - q[0])
        v = pd.to_numeric(seg["Voltage(V)"], errors="coerce").to_numpy(float)
        i = pd.to_numeric(seg["Current(A)"], errors="coerce").to_numpy(float)
        valid = np.isfinite(q) & np.isfinite(v) & np.isfinite(i)
        curves[(rate, charge)] = {
            "t_min": seg.loc[valid, "t_min"].to_numpy(float),
            "V": v[valid],
            "I_A": i[valid],
            "Q_Ah": q[valid],
        }

        rest = record[
            (record.Date >= previous["Oneset Date"])
            & (record.Date <= previous["End Date"])
            & (record["Step Type"].astype(str) == "Rest")
        ].copy()
        rest["t_min"] = (rest.Date - rest.Date.iloc[0]).dt.total_seconds() / 60.0
        tail_width = min(10.0, max(2.0, float(rest.t_min.max()) * 0.25))
        tail = rest[rest.t_min >= rest.t_min.max() - tail_width]
        slope = float(np.polyfit(tail.t_min, tail["Voltage(V)"], 1)[0] * 1000.0)

        cv_capacity = np.nan
        if charge and step_no + 1 in by_number.index:
            following = by_number.loc[step_no + 1]
            if str(following["Step Type"]) == "CV Chg":
                cv_capacity = float(following["Capacity(Ah)"])
        rows.append(
            {
                "Cohort": cohort,
                "Cell": cell,
                "File": path.name,
                "C_rate": rate,
                "Direction": "Charge" if charge else "Discharge",
                "Step_number": step_no,
                "Rest_duration_min": _duration_min(previous),
                "Rest_end_V": float(rest["Voltage(V)"].iloc[-1]),
                "Rest_tail_slope_mV_per_min": slope,
                "Onset_V": float(branch["Oneset Volt.(V)"]),
                "End_V": float(branch["End Voltage(V)"]),
                "Median_current_A": float(np.nanmedian(np.abs(i))),
                "CC_duration_min": float(curves[(rate, charge)]["t_min"][-1]),
                "CC_capacity_Ah": float(curves[(rate, charge)]["Q_Ah"][-1]),
                "CV_capacity_Ah": cv_capacity,
                "CC_plus_CV_Ah": float(curves[(rate, charge)]["Q_Ah"][-1])
                + (0.0 if not np.isfinite(cv_capacity) else cv_capacity),
            }
        )
    return curves, rows


def load_cohorts():
    curves: dict[str, dict[str, dict[tuple[float, bool], dict[str, np.ndarray]]]] = {
        "Old July BoL": {},
        "New September": {},
    }
    rows: list[dict] = []
    for cell, path in OLD_FILES.items():
        cell_curves, cell_rows = parse_cell("Old July BoL", cell, path, OLD_BRANCH_STEPS)
        curves["Old July BoL"][cell] = cell_curves
        rows.extend(cell_rows)
    for cell, path in NEW_FILES.items():
        cell_curves, cell_rows = parse_cell("New September", cell, path, NEW_BRANCH_STEPS)
        curves["New September"][cell] = cell_curves
        rows.extend(cell_rows)
    return curves, pd.DataFrame(rows)


def cohort_summary(qc: pd.DataFrame) -> pd.DataFrame:
    return (
        qc.groupby(["Cohort", "C_rate", "Direction"], as_index=False)
        .agg(
            N_cells=("Cell", "nunique"),
            Mean_CC_capacity_Ah=("CC_capacity_Ah", "mean"),
            SD_CC_capacity_mAh=("CC_capacity_Ah", lambda x: float(np.std(x, ddof=1) * 1000.0)),
            Range_CC_capacity_mAh=("CC_capacity_Ah", lambda x: float(np.ptp(x) * 1000.0)),
            Mean_total_charge_or_CC_Ah=("CC_plus_CV_Ah", "mean"),
            Mean_duration_min=("CC_duration_min", "mean"),
            Mean_onset_V=("Onset_V", "mean"),
            Mean_rest_end_V=("Rest_end_V", "mean"),
            Mean_rest_duration_min=("Rest_duration_min", "mean"),
            Mean_rest_tail_slope_mV_per_min=("Rest_tail_slope_mV_per_min", "mean"),
        )
        .sort_values(["C_rate", "Direction", "Cohort"])
        .reset_index(drop=True)
    )


def between_cohort_capacity(summary: pd.DataFrame) -> pd.DataFrame:
    old = summary[summary.Cohort == "Old July BoL"].set_index(["C_rate", "Direction"])
    new = summary[summary.Cohort == "New September"].set_index(["C_rate", "Direction"])
    rows = []
    for key in old.index:
        rate, direction = key
        old_cc = float(old.loc[key, "Mean_CC_capacity_Ah"])
        new_cc = float(new.loc[key, "Mean_CC_capacity_Ah"])
        old_total = float(old.loc[key, "Mean_total_charge_or_CC_Ah"])
        new_total = float(new.loc[key, "Mean_total_charge_or_CC_Ah"])
        rows.append(
            {
                "C_rate": rate,
                "Direction": direction,
                "Old_CC_Ah": old_cc,
                "New_CC_Ah": new_cc,
                "New_minus_old_CC_mAh": 1000.0 * (new_cc - old_cc),
                "New_minus_old_CC_pct": 100.0 * (new_cc - old_cc) / old_cc,
                "Old_total_charge_or_CC_Ah": old_total,
                "New_total_charge_or_CC_Ah": new_total,
                "New_minus_old_total_mAh": 1000.0 * (new_total - old_total),
                "New_minus_old_total_pct": 100.0 * (new_total - old_total) / old_total,
            }
        )
    return pd.DataFrame(rows).sort_values(["C_rate", "Direction"])


def mean_curve_by_progress(cell_curves, key, grid):
    values = []
    for branches in cell_curves.values():
        curve = branches[key]
        q = np.asarray(curve["Q_Ah"], float)
        progress = q / q[-1]
        values.append(np.interp(grid, progress, curve["V"]))
    array = np.vstack(values)
    return np.mean(array, axis=0), np.std(array, axis=0, ddof=1) if len(array) > 1 else np.zeros_like(grid)


def experimental_curve_differences(curves) -> pd.DataFrame:
    grid = np.linspace(0.0, 1.0, 1001)
    rows = []
    for rate in ga.RATES:
        for charge in (True, False):
            old, _ = mean_curve_by_progress(curves["Old July BoL"], (rate, charge), grid)
            new, _ = mean_curve_by_progress(curves["New September"], (rate, charge), grid)
            err = (new - old) * 1000.0
            center = (grid >= 0.10) & (grid <= 0.90)
            rows.append(
                {
                    "C_rate": rate,
                    "Direction": "Charge" if charge else "Discharge",
                    "Normalized_progress_full_RMSE_mV": float(np.sqrt(np.mean(err**2))),
                    "Normalized_progress_10_90_RMSE_mV": float(np.sqrt(np.mean(err[center] ** 2))),
                    "Normalized_progress_10_90_bias_mV_new_minus_old": float(np.mean(err[center])),
                    "Initial_voltage_delta_mV": float(err[0]),
                    "Midpoint_voltage_delta_mV": float(err[len(grid) // 2]),
                }
            )
    return pd.DataFrame(rows)


def plot_experiment_comparison(curves, qc):
    colors = {"Old July BoL": "#0072B2", "New September": "#D55E00"}
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.7), constrained_layout=True)
    for row, charge in enumerate((True, False)):
        for col, rate in enumerate(ga.RATES):
            ax = axes[row, col]
            for cohort, cells in curves.items():
                for cell, branches in cells.items():
                    branch = branches[(rate, charge)]
                    ax.plot(
                        branch["t_min"],
                        branch["V"],
                        color=colors[cohort],
                        alpha=0.34,
                        lw=1.0,
                    )
                # Add the cohort mean on a normalized-time representation scaled by mean duration.
                grid = np.linspace(0.0, 1.0, 1001)
                mean_v, _ = mean_curve_by_progress(cells, (rate, charge), grid)
                mean_t = float(
                    qc[(qc.Cohort == cohort) & (qc.C_rate == rate) & (qc.Direction == ("Charge" if charge else "Discharge"))]["CC_duration_min"].mean()
                )
                ax.plot(grid * mean_t, mean_v, color=colors[cohort], lw=2.2, label=cohort)
            ax.set_title(f"{rate:g}C {'Charge' if charge else 'Discharge'}")
            ax.set_xlabel("Time [min]")
            ax.set_ylabel("Voltage [V]")
            ax.grid(alpha=0.22)
    axes[0, 0].legend()
    fig.suptitle("Old and new RPT time-voltage curves")
    fig.savefig(RESULTS / "old_vs_new_time_voltage.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.7), constrained_layout=True)
    for row, charge in enumerate((True, False)):
        for col, rate in enumerate(ga.RATES):
            ax = axes[row, col]
            for cohort, cells in curves.items():
                for cell, branches in cells.items():
                    branch = branches[(rate, charge)]
                    ax.plot(branch["Q_Ah"], branch["V"], color=colors[cohort], alpha=0.34, lw=1.0)
            ax.set_title(f"{rate:g}C {'Charge' if charge else 'Discharge'}")
            ax.set_xlabel("Transferred capacity [Ah]")
            ax.set_ylabel("Voltage [V]")
            ax.grid(alpha=0.22)
    axes[0, 0].plot([], [], color=colors["Old July BoL"], label="Old July BoL")
    axes[0, 0].plot([], [], color=colors["New September"], label="New September")
    axes[0, 0].legend()
    fig.suptitle("Old and new RPT capacity-voltage curves")
    fig.savefig(RESULTS / "old_vs_new_capacity_voltage.png", dpi=220)
    plt.close(fig)


def candidate_list() -> list[Candidate]:
    sem_rn, sem_rp = 3.7542, 4.2387
    nominal_rn, nominal_rp = 5.0, 3.0
    gauge_ln, gauge_lp = 73.0, 62.5
    sem_ln, sem_lp = 81.52, 62.965
    consensus_ln, consensus_lp = 76.5, 0.5 * (gauge_lp + sem_lp)
    return [
        Candidate(
            "A Literature kinetics + Ai2020 L + SEM R",
            76.5,
            68.0,
            sem_rn,
            sem_rp,
            AI2020_PREF,
            kpstudy.KP_EIS,
            True,
            physical_note="literature kn; measured R; provisional measured kp",
        ),
        Candidate(
            "B Legacy kn + Ai2020 L + SEM R",
            76.5,
            68.0,
            sem_rn,
            sem_rp,
            current.KN_PREF,
            kpstudy.KP_EIS,
            True,
            physical_note="invalidated anode-Rct kn retained only as numerical comparator",
        ),
        Candidate(
            "C Literature kinetics + gauge L + SEM R",
            gauge_ln,
            gauge_lp,
            sem_rn,
            sem_rp,
            AI2020_PREF,
            kpstudy.KP_EIS,
            True,
            physical_note="direct gauge thickness; measured R",
        ),
        Candidate(
            "D Legacy kn + gauge L + SEM R",
            gauge_ln,
            gauge_lp,
            sem_rn,
            sem_rp,
            current.KN_PREF,
            kpstudy.KP_EIS,
            True,
            physical_note="best-new-data geometry with non-anchored kn",
        ),
        Candidate(
            "E Literature kinetics + SEM L + SEM R",
            sem_ln,
            sem_lp,
            sem_rn,
            sem_rp,
            AI2020_PREF,
            kpstudy.KP_EIS,
            True,
            physical_note="cross-section thickness; measured R",
        ),
        Candidate(
            "F Legacy kn + SEM L + SEM R",
            sem_ln,
            sem_lp,
            sem_rn,
            sem_rp,
            current.KN_PREF,
            kpstudy.KP_EIS,
            True,
            physical_note="best-old-capacity geometry with non-anchored kn",
        ),
        Candidate(
            "G Literature kinetics + consensus L + SEM R",
            consensus_ln,
            consensus_lp,
            sem_rn,
            sem_rp,
            AI2020_PREF,
            kpstudy.KP_EIS,
            True,
            physical_note="Ln bracket center; Lp agreement of gauge and SEM; measured R",
        ),
        Candidate(
            "H Legacy kn + consensus L + SEM R",
            consensus_ln,
            consensus_lp,
            sem_rn,
            sem_rp,
            current.KN_PREF,
            kpstudy.KP_EIS,
            True,
            physical_note="consensus geometry with non-anchored kn",
        ),
        Candidate(
            "I Literature kinetics + Ai2020 L/R",
            76.5,
            68.0,
            nominal_rn,
            nominal_rp,
            AI2020_PREF,
            kpstudy.KP_EIS,
            True,
            physical_note="literature geometry and kinetics comparator",
        ),
        Candidate(
            "J Literature kinetics + consensus L + SEM R, no p hysteresis",
            consensus_ln,
            consensus_lp,
            sem_rn,
            sem_rp,
            AI2020_PREF,
            kpstudy.KP_EIS,
            False,
            physical_note="tests whether old cathode hysteresis transfers to the new cohort",
        ),
        Candidate(
            "K Literature kinetics + Ai2020 L + SEM R, no p hysteresis",
            76.5,
            68.0,
            sem_rn,
            sem_rp,
            AI2020_PREF,
            kpstudy.KP_EIS,
            False,
            physical_note="no-positive-hysteresis comparator",
        ),
        Candidate(
            "L Literature kinetics + consensus L + SEM R, endpoint init",
            consensus_ln,
            consensus_lp,
            sem_rn,
            sem_rp,
            AI2020_PREF,
            kpstudy.KP_EIS,
            True,
            initialization="endpoint",
            physical_note="protocol-inconsistent initialization control",
        ),
        Candidate(
            "M Literature kinetics + gauge L + SEM R, no p hysteresis",
            gauge_ln,
            gauge_lp,
            sem_rn,
            sem_rp,
            AI2020_PREF,
            kpstudy.KP_EIS,
            False,
            physical_note="physical gauge case without positive hysteresis",
        ),
        Candidate(
            "N Legacy kn + gauge L + SEM R, no p hysteresis",
            gauge_ln,
            gauge_lp,
            sem_rn,
            sem_rp,
            current.KN_PREF,
            kpstudy.KP_EIS,
            False,
            physical_note="numerical gauge comparator without positive hysteresis",
        ),
        Candidate(
            "O Literature kinetics + gauge L + SEM R, legacy kp",
            gauge_ln,
            gauge_lp,
            sem_rn,
            sem_rp,
            AI2020_PREF,
            current.KP_PREF,
            True,
            physical_note="tests legacy kp against the cathode-Rct provisional value",
        ),
        Candidate(
            "P Legacy kn/kp + gauge L + SEM R",
            gauge_ln,
            gauge_lp,
            sem_rn,
            sem_rp,
            current.KN_PREF,
            current.KP_PREF,
            True,
            physical_note="fully legacy kinetics numerical comparator",
        ),
    ]


def build_candidate_model(model0: dict, candidate: Candidate) -> dict:
    model = thickness_study.same_loading_geometry(model0, candidate.ln_um, candidate.lp_um)
    model = radius_study.with_radius(model, candidate.rn_um, candidate.rp_um)
    model = kinetics.with_prefactors(model, candidate.kn, candidate.kp)
    return model


def cohort_protocol_means(qc: pd.DataFrame) -> pd.DataFrame:
    return (
        qc.groupby(["Cohort", "C_rate", "Direction"], as_index=False)
        .agg(
            Rest_end_V=("Rest_end_V", "mean"),
            Measured_current_A=("Median_current_A", "mean"),
        )
        .set_index(["Cohort", "C_rate", "Direction"])
    )


def simulate_candidates(curves, qc):
    model0, stage, anode, cathode, _, _, qcell = selected.selected_inputs()
    protocol = cohort_protocol_means(qc)
    rows = []
    simulations = {}
    for candidate in candidate_list():
        model = build_candidate_model(model0, candidate)
        for cohort, cell_curves in curves.items():
            for rate in ga.RATES:
                for charge in (True, False):
                    direction = "Charge" if charge else "Discharge"
                    p = protocol.loc[(cohort, rate, direction)]
                    if candidate.initialization == "rest":
                        z0 = priority.monotone_voltage_inverse(
                            model["v_qocv"], model["soc"], float(p.Rest_end_V)
                        )
                        local_stage = priority.stage_at_soc(stage, z0, charge)
                    else:
                        z0 = 0.0 if charge else 1.0
                        local_stage = stage
                    sim_rate = float(p.Measured_current_A) / NOMINAL_CAPACITY_AH
                    print(
                        f"{candidate.name} | {cohort} | {rate:g}C {direction}",
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
                        candidate.positive_hysteresis,
                        positive_initial_h=(-1.0 if charge else 1.0),
                    )
                    simulations[(candidate.name, cohort, rate, charge)] = sim
                    for cell, branches in cell_curves.items():
                        metric = priority.curve_metrics(
                            branches[(rate, charge)], sim, qcell, charge
                        )
                        rows.append(
                            {
                                "Candidate": candidate.name,
                                "Cohort": cohort,
                                "Cell": cell,
                                "C_rate": rate,
                                "Direction": direction,
                                "Ln_um": candidate.ln_um,
                                "Lp_um": candidate.lp_um,
                                "Rn_um": candidate.rn_um,
                                "Rp_um": candidate.rp_um,
                                "kn": candidate.kn,
                                "kp": candidate.kp,
                                "Positive_hysteresis": candidate.positive_hysteresis,
                                "Initialization": candidate.initialization,
                                "Initial_SOC": z0,
                                "Physical_note": candidate.physical_note,
                                **metric,
                                "Capacity_error_mAh": 1000.0
                                * (metric["Q_end_model_Ah"] - metric["Q_end_exp_Ah"]),
                            }
                        )
    return pd.DataFrame(rows), simulations


def summarize_candidates(detail: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    by_cohort = (
        detail.groupby(["Candidate", "Cohort"], sort=False, as_index=False)
        .agg(
            Mean_full_MAE_mV=("Full_MAE_mV", "mean"),
            Mean_full_RMSE_mV=("Full_RMSE_mV", "mean"),
            Mean_center_MAE_mV=("Center10_70_MAE_mV", "mean"),
            Mean_center_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
            Capacity_RMSE_pct=("Capacity_error_pct", lambda x: float(np.sqrt(np.mean(x**2)))),
            Mean_abs_capacity_error_mAh=("Capacity_error_mAh", lambda x: float(np.mean(np.abs(x)))),
            Max_abs_capacity_error_mAh=("Capacity_error_mAh", lambda x: float(np.max(np.abs(x)))),
        )
    )
    pivot = by_cohort.pivot(index="Candidate", columns="Cohort")
    robust_rows = []
    for candidate in by_cohort.Candidate.unique():
        group = by_cohort[by_cohort.Candidate == candidate].set_index("Cohort")
        old = group.loc["Old July BoL"]
        new = group.loc["New September"]
        sample = detail[detail.Candidate == candidate].iloc[0]
        robust_rows.append(
            {
                "Candidate": candidate,
                "Ln_um": sample.Ln_um,
                "Lp_um": sample.Lp_um,
                "Rn_um": sample.Rn_um,
                "Rp_um": sample.Rp_um,
                "kn": sample.kn,
                "kp": sample.kp,
                "Positive_hysteresis": sample.Positive_hysteresis,
                "Initialization": sample.Initialization,
                "Physical_note": sample.Physical_note,
                "Old_full_RMSE_mV": old.Mean_full_RMSE_mV,
                "New_full_RMSE_mV": new.Mean_full_RMSE_mV,
                "Old_center_RMSE_mV": old.Mean_center_RMSE_mV,
                "New_center_RMSE_mV": new.Mean_center_RMSE_mV,
                "Old_capacity_RMSE_pct": old.Capacity_RMSE_pct,
                "New_capacity_RMSE_pct": new.Capacity_RMSE_pct,
                "Old_mean_abs_capacity_error_mAh": old.Mean_abs_capacity_error_mAh,
                "New_mean_abs_capacity_error_mAh": new.Mean_abs_capacity_error_mAh,
                "Worst_cohort_full_RMSE_mV": max(old.Mean_full_RMSE_mV, new.Mean_full_RMSE_mV),
                "Worst_cohort_center_RMSE_mV": max(old.Mean_center_RMSE_mV, new.Mean_center_RMSE_mV),
                "Worst_cohort_capacity_RMSE_pct": max(old.Capacity_RMSE_pct, new.Capacity_RMSE_pct),
                "Mean_two_cohort_full_RMSE_mV": 0.5 * (old.Mean_full_RMSE_mV + new.Mean_full_RMSE_mV),
                "Mean_two_cohort_capacity_RMSE_pct": 0.5 * (old.Capacity_RMSE_pct + new.Capacity_RMSE_pct),
            }
        )
    robust = pd.DataFrame(robust_rows)
    metrics = [
        "Old_center_RMSE_mV",
        "New_center_RMSE_mV",
        "Old_capacity_RMSE_pct",
        "New_capacity_RMSE_pct",
    ]
    values = robust[metrics].to_numpy(float)
    pareto = np.ones(len(robust), dtype=bool)
    for i in range(len(robust)):
        for j in range(len(robust)):
            if i == j:
                continue
            if np.all(values[j] <= values[i]) and np.any(values[j] < values[i]):
                pareto[i] = False
                break
    robust["Pareto_4metric"] = pareto
    robust = robust.sort_values(
        ["Pareto_4metric", "Worst_cohort_full_RMSE_mV", "Worst_cohort_capacity_RMSE_pct"],
        ascending=[False, True, True],
    ).reset_index(drop=True)
    return by_cohort, robust


def plot_candidate_tradeoff(robust: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(9.5, 6.8), constrained_layout=True)
    colors = np.where(robust.Pareto_4metric, "#D55E00", "#8A8A8A")
    ax.scatter(
        robust.Worst_cohort_center_RMSE_mV,
        robust.Worst_cohort_capacity_RMSE_pct,
        c=colors,
        s=85,
        alpha=0.9,
    )
    for _, row in robust.iterrows():
        ax.annotate(
            str(row.Candidate).split()[0],
            (row.Worst_cohort_center_RMSE_mV, row.Worst_cohort_capacity_RMSE_pct),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=9,
        )
    ax.set_xlabel("Worst-cohort 10-70% voltage RMSE [mV]")
    ax.set_ylabel("Worst-cohort capacity RMSE [%]")
    ax.grid(alpha=0.25)
    ax.set_title("Fixed-parameter pre-fit candidates across old and new cohorts")
    fig.savefig(RESULTS / "cross_cohort_candidate_tradeoff.png", dpi=220)
    plt.close(fig)


def main():
    curves, qc = load_cohorts()
    qc.to_csv(RESULTS / "branch_qc_all_cells.csv", index=False, encoding="utf-8-sig")
    summary = cohort_summary(qc)
    summary.to_csv(RESULTS / "cohort_branch_summary.csv", index=False, encoding="utf-8-sig")
    capacity = between_cohort_capacity(summary)
    capacity.to_csv(RESULTS / "old_vs_new_capacity_shift.csv", index=False, encoding="utf-8-sig")
    curve_diff = experimental_curve_differences(curves)
    curve_diff.to_csv(RESULTS / "old_vs_new_normalized_shape_difference.csv", index=False, encoding="utf-8-sig")
    plot_experiment_comparison(curves, qc)

    detail, _ = simulate_candidates(curves, qc)
    detail.to_csv(RESULTS / "candidate_detail.csv", index=False, encoding="utf-8-sig")
    by_cohort, robust = summarize_candidates(detail)
    by_cohort.to_csv(RESULTS / "candidate_by_cohort.csv", index=False, encoding="utf-8-sig")
    robust.to_csv(RESULTS / "candidate_robust_comparison.csv", index=False, encoding="utf-8-sig")
    plot_candidate_tradeoff(robust)

    report = {
        "source_files": {
            "old": [str(path) for path in OLD_FILES.values()],
            "new": [str(path) for path in NEW_FILES.values()],
        },
        "definitions": {
            "curve_error": "absolute transferred-capacity domain over 99.5% of common support",
            "center_error": "absolute 0.10-0.70 times 2.356528 Ah target capacity",
            "capacity_error": "model CC endpoint minus experimental CC endpoint",
            "candidate_selection": "four-metric Pareto only; no weighted scalar score",
        },
        "fixed_parameters": {
            "Dsn_m2_s": 2.1e-14,
            "Dsp_m2_s": 4.4e-14,
            "brugg_n": 2.914,
            "brugg_p": 1.83,
            "brugg_s": 1.5,
            "kp_prefactor": kpstudy.KP_EIS,
            "OCP": "Ai2020 nominal anode + old cathode; negative hysteresis off",
        },
        "counts": {
            "experimental_branches": int(len(qc)),
            "candidate_rows": int(len(detail)),
            "candidates": int(robust.Candidate.nunique()),
        },
    }
    (RESULTS / "analysis_manifest.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("\nCAPACITY SHIFT\n", capacity.to_string(index=False))
    print("\nCURVE DIFFERENCE\n", curve_diff.to_string(index=False))
    print("\nROBUST CANDIDATES\n", robust.to_string(index=False))


if __name__ == "__main__":
    main()
