"""260918 anode GITT OCP candidate analysis.

This module intentionally preserves the Stage 0/Stage 1 logic and DFN setup from
``260908 Charge Fitting Discharge Validation.ipynb``.  The only comparison
variable is the negative-electrode OCP source.
"""

from __future__ import annotations

import json
import os
import platform
import time
import warnings
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any

os.environ.setdefault("PYBAMM_DISABLE_TELEMETRY", "true")
os.environ.setdefault("PYBAMM_TELEMETRY", "false")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import openpyxl
import pandas as pd
import pybamm
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import CubicSpline
from scipy.optimize import least_squares


warnings.filterwarnings(
    "ignore",
    message="Workbook contains no default style",
)

RANDOM_SEED = 260918
np.random.seed(RANDOM_SEED)

FARADAY_CONSTANT = 96485.33212
DE_MULTIPLIER = 1e-4
CSN_MAX_CORRECTED = 29700.0
RATES = (0.5, 1.0, 2.0)


@dataclass
class AnalysisPaths:
    root: Path
    legacy_book: Path
    gitt_book: Path
    results: Path


def get_paths(root: Path | str | None = None) -> AnalysisPaths:
    root_path = Path(root or Path.cwd()).resolve()
    return AnalysisPaths(
        root=root_path,
        legacy_book=root_path / "260902 Enertech셀 데이터 모음.xlsx",
        gitt_book=(
            root_path
            / "data"
            / "raw"
            / "260918_Anode_GITT"
            / "260918 Anode GITT v2-v3.xlsx"
        ),
        results=root_path / "results" / "260918_gitt_ocp_comparison",
    )


def _parse_datetime(value: Any) -> datetime:
    if isinstance(value, datetime):
        return value
    return datetime.fromisoformat(str(value))


def _duration_seconds(value: Any) -> float:
    if isinstance(value, timedelta):
        return value.total_seconds()
    if isinstance(value, (int, float)):
        return float(value) * 86400.0
    pieces = str(value).split(":")
    if len(pieces) != 3:
        raise ValueError(f"Unsupported duration: {value!r}")
    hours, minutes, seconds = (float(piece) for piece in pieces)
    return hours * 3600.0 + minutes * 60.0 + seconds


def _stream_sheet_rows(workbook: openpyxl.Workbook, sheet_name: str):
    worksheet = workbook[sheet_name]
    # The source workbook advertises A1:A1 dimensions despite having many
    # populated columns/rows. Resetting dimensions is required for read-only
    # streaming and does not modify the source file.
    worksheet.reset_dimensions()
    yield from worksheet.iter_rows(values_only=True)


def _read_small_sheet(workbook: openpyxl.Workbook, sheet_name: str) -> pd.DataFrame:
    rows = list(_stream_sheet_rows(workbook, sheet_name))
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows[1:], columns=rows[0])


def _prepare_step_table(frame: pd.DataFrame, replicate: str) -> pd.DataFrame:
    out = frame.copy()
    out["replicate"] = replicate
    out["Step Number"] = pd.to_numeric(out["Step Number"], errors="coerce").astype(int)
    out["Cycle Index"] = pd.to_numeric(out["Cycle Index"], errors="coerce").astype(int)
    out["Capacity(Ah)"] = pd.to_numeric(out["Capacity(Ah)"], errors="coerce")
    out["Oneset Date"] = pd.to_datetime(out["Oneset Date"])
    out["End Date"] = pd.to_datetime(out["End Date"])
    out["step_duration_s"] = out["Step Time"].map(_duration_seconds)
    out["actual_duration_s"] = (
        out["End Date"] - out["Oneset Date"]
    ).dt.total_seconds()
    return out


def _annotate_gitt_steps(steps: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    cycle2 = steps[steps["Cycle Index"] == 2].copy().reset_index(drop=True)
    rest_2h = (
        cycle2["Step Type"].eq("Rest")
        & cycle2["step_duration_s"].between(7100, 7300)
    )
    if not rest_2h.any():
        raise RuntimeError(f"No 2 h GITT rest found for {steps['replicate'].iloc[0]}")
    first_index = int(np.flatnonzero(rest_2h.to_numpy())[0])
    gitt = cycle2.iloc[first_index:].copy().reset_index(drop=True)

    lith_capacity = 0.0
    delith_capacity = 0.0
    lith_pulse = 0
    delith_pulse = 0
    last_branch: str | None = None
    pulse_rows: list[dict[str, Any]] = []
    rest_rows: list[dict[str, Any]] = []

    for row_index, row in gitt.iterrows():
        step_type = str(row["Step Type"])
        if step_type == "CC DChg":
            lith_pulse += 1
            lith_capacity += float(row["Capacity(Ah)"])
            last_branch = "lithiation"
            pulse_rows.append(
                {
                    **row.to_dict(),
                    "branch": last_branch,
                    "pulse_no": lith_pulse,
                    "cumulative_capacity_Ah": lith_capacity,
                    "cutoff_pulse": bool(row["step_duration_s"] < 590),
                }
            )
        elif step_type == "CC Chg":
            delith_pulse += 1
            delith_capacity += float(row["Capacity(Ah)"])
            last_branch = "delithiation"
            pulse_rows.append(
                {
                    **row.to_dict(),
                    "branch": last_branch,
                    "pulse_no": delith_pulse,
                    "cumulative_capacity_Ah": delith_capacity,
                    "cutoff_pulse": bool(row["step_duration_s"] < 590),
                }
            )
        elif step_type == "Rest" and 7100 <= float(row["step_duration_s"]) <= 7300:
            is_initial = row_index == 0
            branch = "lithiation" if is_initial else last_branch
            if branch is None:
                continue
            pulse_no = 0 if is_initial else (lith_pulse if branch == "lithiation" else delith_pulse)
            cumulative = (
                0.0
                if is_initial
                else (lith_capacity if branch == "lithiation" else delith_capacity)
            )
            rest_rows.append(
                {
                    **row.to_dict(),
                    "branch": branch,
                    "pulse_no": pulse_no,
                    "cumulative_capacity_Ah": cumulative,
                    "is_initial_gitt_rest": is_initial,
                }
            )

    pulse_table = pd.DataFrame(pulse_rows)
    rest_table = pd.DataFrame(rest_rows)
    if pulse_table.empty or rest_table.empty:
        raise RuntimeError(f"Could not identify GITT pulses/rests for {steps['replicate'].iloc[0]}")
    return pulse_table, rest_table


def _stream_rest_measurements(
    workbook: openpyxl.Workbook,
    sheet_name: str,
    rest_table: pd.DataFrame,
    replicate: str,
) -> pd.DataFrame:
    intervals: list[dict[str, Any]] = []
    for _, row in rest_table.sort_values("Oneset Date").iterrows():
        intervals.append(
            {
                "metadata": row.to_dict(),
                "start": row["Oneset Date"].to_pydatetime(),
                "end": row["End Date"].to_pydatetime(),
                "first_time": None,
                "first_voltage": None,
                "last_time": None,
                "last_voltage": None,
                "tail_time_min": [],
                "tail_voltage": [],
                "record_count": 0,
            }
        )

    row_iter = _stream_sheet_rows(workbook, sheet_name)
    header = next(row_iter)
    index = {str(name): position for position, name in enumerate(header)}
    required = ("Step Type", "Voltage(V)", "Date")
    missing = [name for name in required if name not in index]
    if missing:
        raise RuntimeError(f"Missing {missing} in {sheet_name}")

    interval_index = 0
    total_rows = 0
    for values in row_iter:
        total_rows += 1
        if total_rows % 250000 == 0:
            print(f"{replicate}: streamed {total_rows:,} record rows", flush=True)
        if str(values[index["Step Type"]]) != "Rest":
            continue
        raw_date = values[index["Date"]]
        raw_voltage = values[index["Voltage(V)"]]
        if raw_date is None or raw_voltage is None:
            continue
        timestamp = _parse_datetime(raw_date)

        while interval_index < len(intervals) and timestamp > intervals[interval_index]["end"]:
            interval_index += 1
        if interval_index >= len(intervals):
            break
        interval = intervals[interval_index]
        if timestamp < interval["start"]:
            continue
        voltage = float(raw_voltage)
        if interval["first_time"] is None:
            interval["first_time"] = timestamp
            interval["first_voltage"] = voltage
        interval["last_time"] = timestamp
        interval["last_voltage"] = voltage
        interval["record_count"] += 1
        if timestamp >= interval["end"] - timedelta(minutes=10):
            interval["tail_time_min"].append(
                (timestamp - (interval["end"] - timedelta(minutes=10))).total_seconds() / 60.0
            )
            interval["tail_voltage"].append(voltage)

    output_rows: list[dict[str, Any]] = []
    for interval in intervals:
        metadata = interval["metadata"]
        tail_t = np.asarray(interval["tail_time_min"], dtype=float)
        tail_v = np.asarray(interval["tail_voltage"], dtype=float)
        if len(tail_t) >= 2 and np.ptp(tail_t) > 0:
            slope_v_per_min = float(np.polyfit(tail_t, tail_v, 1)[0])
        else:
            slope_v_per_min = np.nan
        measured_duration = (
            (interval["last_time"] - interval["first_time"]).total_seconds()
            if interval["first_time"] is not None and interval["last_time"] is not None
            else np.nan
        )
        output_rows.append(
            {
                "replicate": replicate,
                "step_number": int(metadata["Step Number"]),
                "branch": metadata["branch"],
                "pulse_no": int(metadata["pulse_no"]),
                "cumulative_capacity_Ah": float(metadata["cumulative_capacity_Ah"]),
                "rest_start": metadata["Oneset Date"],
                "rest_end": metadata["End Date"],
                "programmed_rest_duration_min": float(metadata["step_duration_s"]) / 60.0,
                "measured_rest_duration_min": measured_duration / 60.0,
                "rest_start_voltage_V": interval["first_voltage"],
                "rest_end_voltage_V": interval["last_voltage"],
                "last_10_min_slope_mV_per_min": slope_v_per_min * 1000.0,
                "last_10_min_abs_dVdt_mV_per_min": abs(slope_v_per_min) * 1000.0,
                "record_count": int(interval["record_count"]),
                "tail_record_count": int(len(tail_t)),
                "is_initial_gitt_rest": bool(metadata["is_initial_gitt_rest"]),
            }
        )
    return pd.DataFrame(output_rows)


def load_new_gitt(paths: AnalysisPaths) -> dict[str, Any]:
    if not paths.gitt_book.exists():
        raise FileNotFoundError(paths.gitt_book)
    started = time.time()
    workbook = openpyxl.load_workbook(paths.gitt_book, read_only=True, data_only=True)
    info_rows = list(_stream_sheet_rows(workbook, "INFO"))
    info = {str(key): value for key, value, *_ in info_rows}
    replicates: dict[str, dict[str, pd.DataFrame]] = {}
    all_rests: list[pd.DataFrame] = []
    cached_rest_path = paths.results / "gitt_rest_equilibrium.csv"
    cached_rests = None
    if (
        cached_rest_path.exists()
        and cached_rest_path.stat().st_mtime >= paths.gitt_book.stat().st_mtime
        and cached_rest_path.stat().st_size > 100
    ):
        cached_rests = pd.read_csv(
            cached_rest_path,
            parse_dates=["rest_start", "rest_end"],
        )
        print("Using cached rest-end extraction; raw workbook is unchanged", flush=True)

    for replicate in ("v2", "v3"):
        steps = _prepare_step_table(_read_small_sheet(workbook, f"{replicate}_step"), replicate)
        pulses, rests = _annotate_gitt_steps(steps)
        branch_capacity = pulses.groupby("branch")["Capacity(Ah)"].sum().to_dict()
        if cached_rests is None:
            rest_measurements = _stream_rest_measurements(
                workbook,
                f"{replicate}_record",
                rests,
                replicate,
            )
            rest_measurements["branch_capacity_Ah"] = rest_measurements["branch"].map(branch_capacity)
            rest_measurements["q_norm"] = (
                rest_measurements["cumulative_capacity_Ah"]
                / rest_measurements["branch_capacity_Ah"]
            )
        else:
            rest_measurements = cached_rests[cached_rests["replicate"] == replicate].copy()
        replicates[replicate] = {
            "steps": steps,
            "pulses": pulses,
            "rests": rest_measurements,
        }
        all_rests.append(rest_measurements)

    workbook.close()
    print(f"Raw GITT streaming completed in {(time.time() - started) / 60:.1f} min", flush=True)
    return {
        "info": info,
        "replicates": replicates,
        "rests": pd.concat(all_rests, ignore_index=True),
    }


def _clean_curve(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    valid = np.isfinite(x) & np.isfinite(y)
    x = np.asarray(x, dtype=float)[valid]
    y = np.asarray(y, dtype=float)[valid]
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    unique_x, inverse = np.unique(x, return_inverse=True)
    if len(unique_x) == len(x):
        return x, y
    sums = np.bincount(inverse, weights=y)
    counts = np.bincount(inverse)
    return unique_x, sums / counts


def _interp_curve(curve: tuple[np.ndarray, np.ndarray], grid: np.ndarray) -> np.ndarray:
    return np.interp(grid, curve[0], curve[1])


def _mean_overlap(curves: list[tuple[np.ndarray, np.ndarray]], n: int = 3001):
    lower = max(float(curve[0].min()) for curve in curves)
    upper = min(float(curve[0].max()) for curve in curves)
    if upper <= lower:
        raise RuntimeError("OCP curves have no common overlap")
    grid = np.linspace(lower, upper, n)
    values = np.vstack([_interp_curve(curve, grid) for curve in curves])
    return grid, values.mean(axis=0)


def build_new_ocp(gitt: dict[str, Any]) -> dict[str, Any]:
    curves: dict[str, dict[str, tuple[np.ndarray, np.ndarray]]] = {}
    for replicate, tables in gitt["replicates"].items():
        rests = tables["rests"]
        lith = rests[rests["branch"] == "lithiation"].copy()
        delith = rests[rests["branch"] == "delithiation"].copy()
        if lith.empty or delith.empty:
            raise RuntimeError(f"Both GITT branches are required for {replicate}")

        # The last lithiation rest is also the q=0 equilibrium point for the
        # subsequent delithiation branch.
        shared_endpoint = lith.sort_values("pulse_no").iloc[-1]
        delith_zero = pd.DataFrame(
            [
                {
                    "q_norm": 0.0,
                    "rest_end_voltage_V": shared_endpoint["rest_end_voltage_V"],
                }
            ]
        )
        delith_points = pd.concat(
            [delith_zero, delith[["q_norm", "rest_end_voltage_V"]]],
            ignore_index=True,
        )

        lith_curve = _clean_curve(
            lith["q_norm"].to_numpy(dtype=float),
            lith["rest_end_voltage_V"].to_numpy(dtype=float),
        )
        delith_curve = _clean_curve(
            1.0 - delith_points["q_norm"].to_numpy(dtype=float),
            delith_points["rest_end_voltage_V"].to_numpy(dtype=float),
        )
        mean_curve = _mean_overlap([lith_curve, delith_curve])
        curves[replicate] = {
            "lithiation": lith_curve,
            "delithiation": delith_curve,
            "mean": mean_curve,
        }

    representative = _mean_overlap([curves["v2"]["mean"], curves["v3"]["mean"]])
    return {"curves": curves, "representative": representative}


def reproducibility_table(gitt: dict[str, Any], new_ocp: dict[str, Any]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for branch in ("lithiation", "delithiation", "mean"):
        v2_curve = new_ocp["curves"]["v2"][branch]
        v3_curve = new_ocp["curves"]["v3"][branch]
        lower = max(v2_curve[0].min(), v3_curve[0].min())
        upper = min(v2_curve[0].max(), v3_curve[0].max())
        grid = np.linspace(lower, upper, 1001)
        difference_mV = (
            _interp_curve(v2_curve, grid) - _interp_curve(v3_curve, grid)
        ) * 1000.0
        for region, mask in (
            ("full_overlap", np.ones_like(grid, dtype=bool)),
            ("normalized_capacity_10_90", (grid >= 0.10) & (grid <= 0.90)),
        ):
            values = difference_mV[mask]
            v2_capacity = np.nan
            v3_capacity = np.nan
            capacity_difference_pct = np.nan
            if branch != "mean":
                v2_capacity = float(
                    gitt["replicates"]["v2"]["pulses"]
                    .query("branch == @branch")["Capacity(Ah)"]
                    .sum()
                )
                v3_capacity = float(
                    gitt["replicates"]["v3"]["pulses"]
                    .query("branch == @branch")["Capacity(Ah)"]
                    .sum()
                )
                capacity_difference_pct = (
                    abs(v2_capacity - v3_capacity)
                    / np.mean([v2_capacity, v3_capacity])
                    * 100.0
                )
            rows.append(
                {
                    "branch": branch,
                    "region": region,
                    "q_norm_min": float(grid[mask].min()),
                    "q_norm_max": float(grid[mask].max()),
                    "n_grid": int(mask.sum()),
                    "MAE_mV": float(np.mean(np.abs(values))),
                    "RMSE_mV": float(np.sqrt(np.mean(values**2))),
                    "median_abs_difference_mV": float(np.median(np.abs(values))),
                    "maximum_abs_difference_mV": float(np.max(np.abs(values))),
                    "v2_branch_capacity_Ah": v2_capacity,
                    "v3_branch_capacity_Ah": v3_capacity,
                    "capacity_difference_pct": capacity_difference_pct,
                }
            )
    return pd.DataFrame(rows)


def branch_hysteresis_table(new_ocp: dict[str, Any]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for replicate in ("v2", "v3"):
        lith = new_ocp["curves"][replicate]["lithiation"]
        delith = new_ocp["curves"][replicate]["delithiation"]
        lower = max(lith[0].min(), delith[0].min())
        upper = min(lith[0].max(), delith[0].max())
        grid = np.linspace(lower, upper, 1001)
        delta_mV = (_interp_curve(lith, grid) - _interp_curve(delith, grid)) * 1000.0
        for region, mask in (
            ("full_overlap", np.ones_like(grid, dtype=bool)),
            ("normalized_capacity_10_90", (grid >= 0.10) & (grid <= 0.90)),
        ):
            values = delta_mV[mask]
            rows.append(
                {
                    "replicate": replicate,
                    "region": region,
                    "mean_lithiation_minus_delithiation_mV": float(np.mean(values)),
                    "MAE_hysteresis_mV": float(np.mean(np.abs(values))),
                    "RMSE_hysteresis_mV": float(np.sqrt(np.mean(values**2))),
                    "median_abs_hysteresis_mV": float(np.median(np.abs(values))),
                    "maximum_abs_hysteresis_mV": float(np.max(np.abs(values))),
                }
            )
    return pd.DataFrame(rows)


def data_qc_tables(gitt: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_rows: list[dict[str, Any]] = []
    capacity_rows: list[dict[str, Any]] = []
    rests = gitt["rests"]
    for replicate, tables in gitt["replicates"].items():
        steps = tables["steps"]
        pulses = tables["pulses"]
        rep_rests = rests[rests["replicate"] == replicate]
        experiment_start = steps["Oneset Date"].min()
        experiment_end = steps["End Date"].max()
        gitt_start = rep_rests["rest_start"].min()
        gitt_end = rep_rests["rest_end"].max()
        for branch in ("lithiation", "delithiation"):
            branch_pulses = pulses[pulses["branch"] == branch]
            branch_rests = rep_rests[rep_rests["branch"] == branch]
            slopes = branch_rests["last_10_min_abs_dVdt_mV_per_min"].dropna()
            summary_rows.append(
                {
                    "replicate": replicate,
                    "branch": branch,
                    "experiment_start": experiment_start,
                    "experiment_end": experiment_end,
                    "experiment_duration_h": (experiment_end - experiment_start).total_seconds() / 3600.0,
                    "gitt_start": gitt_start,
                    "gitt_end": gitt_end,
                    "gitt_duration_h": (gitt_end - gitt_start).total_seconds() / 3600.0,
                    "pulse_count": int(len(branch_pulses)),
                    "rest_count": int(len(branch_rests)),
                    "median_abs_dVdt_mV_per_min": float(slopes.median()),
                    "mean_abs_dVdt_mV_per_min": float(slopes.mean()),
                    "P90_abs_dVdt_mV_per_min": float(slopes.quantile(0.90)),
                    "maximum_abs_dVdt_mV_per_min": float(slopes.max()),
                    "normal_pulse_duration_min_median": float(
                        branch_pulses.loc[~branch_pulses["cutoff_pulse"], "step_duration_s"].median()
                        / 60.0
                    ),
                    "rest_duration_min_median": float(branch_rests["measured_rest_duration_min"].median()),
                    "cutoff_reached": bool(branch_pulses["cutoff_pulse"].any()),
                    "cutoff_pulse_duration_min": float(
                        branch_pulses.loc[branch_pulses["cutoff_pulse"], "step_duration_s"].min()
                        / 60.0
                    ) if branch_pulses["cutoff_pulse"].any() else np.nan,
                }
            )
            capacity_rows.append(
                {
                    "replicate": replicate,
                    "capacity_type": f"gitt_{branch}",
                    "step_number": np.nan,
                    "capacity_Ah": float(branch_pulses["Capacity(Ah)"].sum()),
                    "start_voltage_V": float(branch_pulses.iloc[0]["Oneset Volt.(V)"]),
                    "end_voltage_V": float(branch_pulses.iloc[-1]["End Voltage(V)"]),
                }
            )
        all_slopes = rep_rests["last_10_min_abs_dVdt_mV_per_min"].dropna()
        summary_rows.append(
            {
                "replicate": replicate,
                "branch": "all_GITT_rests",
                "experiment_start": experiment_start,
                "experiment_end": experiment_end,
                "experiment_duration_h": (experiment_end - experiment_start).total_seconds() / 3600.0,
                "gitt_start": gitt_start,
                "gitt_end": gitt_end,
                "gitt_duration_h": (gitt_end - gitt_start).total_seconds() / 3600.0,
                "pulse_count": int(len(pulses)),
                "rest_count": int(len(rep_rests)),
                "median_abs_dVdt_mV_per_min": float(all_slopes.median()),
                "mean_abs_dVdt_mV_per_min": float(all_slopes.mean()),
                "P90_abs_dVdt_mV_per_min": float(all_slopes.quantile(0.90)),
                "maximum_abs_dVdt_mV_per_min": float(all_slopes.max()),
                "normal_pulse_duration_min_median": float(
                    pulses.loc[~pulses["cutoff_pulse"], "step_duration_s"].median() / 60.0
                ),
                "rest_duration_min_median": float(rep_rests["measured_rest_duration_min"].median()),
                "cutoff_reached": bool(pulses["cutoff_pulse"].any()),
                "cutoff_pulse_duration_min": float(
                    pulses.loc[pulses["cutoff_pulse"], "step_duration_s"].min() / 60.0
                ),
            }
        )
        initial_full = steps[
            (steps["Cycle Index"].isin([1, 2]))
            & steps["Step Type"].isin(["CC DChg", "CC Chg"])
        ].head(4)
        for _, row in initial_full.iterrows():
            capacity_rows.append(
                {
                    "replicate": replicate,
                    "capacity_type": (
                        "initial_full_lithiation" if row["Step Type"] == "CC DChg" else "initial_full_delithiation"
                    ),
                    "step_number": int(row["Step Number"]),
                    "capacity_Ah": float(row["Capacity(Ah)"]),
                    "start_voltage_V": float(row["Oneset Volt.(V)"]),
                    "end_voltage_V": float(row["End Voltage(V)"]),
                }
            )
    return pd.DataFrame(summary_rows), pd.DataFrame(capacity_rows)


def load_legacy_bundle(paths: AnalysisPaths) -> dict[str, Any]:
    if not paths.legacy_book.exists():
        raise FileNotFoundError(paths.legacy_book)
    anode_lowrate = pd.read_excel(paths.legacy_book, sheet_name="ANODE_LOWRATE", header=8)
    cathode_gitt = pd.read_excel(paths.legacy_book, sheet_name="CATHODE_GITT", header=8)
    full_c20 = pd.read_excel(paths.legacy_book, sheet_name="FULLCELL_LOWRATE", header=8)
    dynamic = pd.read_excel(paths.legacy_book, sheet_name="DYNAMIC", header=8)
    prior_frame = pd.read_excel(paths.legacy_book, sheet_name="PARAMETER_PRIOR", header=7)

    def prior(key: str, default=None, cast=float):
        selected = prior_frame[prior_frame["parameter"].astype(str) == str(key)]
        if selected.empty:
            return default
        value = selected.iloc[0]["value"]
        if cast is None:
            return value
        try:
            return cast(value)
        except Exception:
            return default

    return {
        "anode_lowrate": anode_lowrate,
        "cathode_gitt": cathode_gitt,
        "full_c20": full_c20,
        "dynamic": dynamic,
        "prior": prior,
        "prior_frame": prior_frame,
    }


def scale_parameter_function(function, factor):
    def wrapped(*args, **kwargs):
        return factor * function(*args, **kwargs)

    return wrapped


def make_source_corrected_ai2020():
    base = pybamm.ParameterValues("Ai2020")
    values = base.copy()
    values.update(
        {
            "Electrolyte diffusivity [m2.s-1]": scale_parameter_function(
                base["Electrolyte diffusivity [m2.s-1]"], DE_MULTIPLIER
            ),
            "Maximum concentration in negative electrode [mol.m-3]": CSN_MAX_CORRECTED,
        }
    )
    return values


def _extract_qv_branch(frame, source_col, source, direction):
    selected = frame[(frame[source_col] == source) & (frame["direction"] == direction)]
    return (
        selected["Q_Ah"].to_numpy(dtype=float),
        selected["V_V"].to_numpy(dtype=float),
    )


def _extract_gitt_branch(frame, source_col, source, pulse_direction):
    selected = frame[
        (frame[source_col] == source) & (frame["pulse_direction"] == pulse_direction)
    ]
    return np.c_[
        selected["Qcum_Ah"].to_numpy(dtype=float),
        selected["Vrest_V"].to_numpy(dtype=float),
    ]


def build_legacy_model_inputs(bundle: dict[str, Any]) -> dict[str, Any]:
    prior = bundle["prior"]
    params = make_source_corrected_ai2020()
    csn_max = float(params["Maximum concentration in negative electrode [mol.m-3]"])
    csp_max = float(params["Maximum concentration in positive electrode [mol.m-3]"])
    eps_n = float(params["Negative electrode active material volume fraction"])
    eps_p = float(params["Positive electrode active material volume fraction"])
    length_n = float(params["Negative electrode thickness [m]"])
    length_p = float(params["Positive electrode thickness [m]"])
    disk_diameter = prior("disk_diameter_m")
    temperature = prior("temperature_K")
    disk_area = np.pi * (disk_diameter / 2.0) ** 2
    qn_disk = FARADAY_CONSTANT * csn_max * eps_n * length_n * disk_area / 3600.0
    qp_disk = FARADAY_CONSTANT * csp_max * eps_p * length_p * disk_area / 3600.0

    anode = bundle["anode_lowrate"]
    anode_curves: list[tuple[np.ndarray, np.ndarray]] = []
    for source in ("v1", "v2"):
        charge = _extract_qv_branch(anode, "source_version", source, "charge")
        discharge = _extract_qv_branch(anode, "source_version", source, "discharge")
        x_lith = discharge[0] / qn_disk
        lith = _clean_curve(x_lith, discharge[1])
        span = charge[0].max() / qn_disk
        x_delith = span - charge[0] / qn_disk
        delith = _clean_curve(x_delith, charge[1])
        anode_curves.extend([lith, delith])
    x_legacy, un_legacy = _mean_overlap(anode_curves)

    cathode = bundle["cathode_gitt"]
    cathode_curves: list[tuple[np.ndarray, np.ndarray]] = []
    for source in ("7-7", "7-8"):
        for pulse_direction, direction in (
            ("discharge_pulse", "lith"),
            ("charge_pulse", "delith"),
        ):
            branch = _extract_gitt_branch(cathode, "source_channel", source, pulse_direction)
            if direction == "lith":
                span = branch[:, 0].max() / qp_disk
                y = (1.0 - span) + branch[:, 0] / qp_disk
            else:
                y = 1.0 - branch[:, 0] / qp_disk
            cathode_curves.append(_clean_curve(y, branch[:, 1]))
    y_exp, up_exp = _mean_overlap(cathode_curves)

    full = bundle["full_c20"]
    capacity_charge, voltage_charge = _extract_qv_branch(full, "source_version", "v2", "charge")
    capacity_discharge, voltage_discharge = _extract_qv_branch(full, "source_version", "v2", "discharge")
    soc = np.linspace(0.0, 1.0, 1001)
    v_charge = _interp_curve(
        _clean_curve(capacity_charge / capacity_charge.max(), voltage_charge), soc
    )
    v_discharge = _interp_curve(
        _clean_curve(1.0 - capacity_discharge / capacity_discharge.max(), voltage_discharge), soc
    )
    v_qocv = 0.5 * (v_charge + v_discharge)
    q_meas = 0.5 * (capacity_charge.max() + capacity_discharge.max())

    cell_area = (
        float(params["Electrode height [m]"])
        * float(params["Electrode width [m]"])
        * float(params["Number of electrodes connected in parallel to make a cell"])
    )
    qn_cell = FARADAY_CONSTANT * csn_max * eps_n * length_n * cell_area / 3600.0
    qp_cell = FARADAY_CONSTANT * csp_max * eps_p * length_p * cell_area / 3600.0
    delta_x = q_meas / qn_cell
    delta_y = q_meas / qp_cell

    return {
        "params": params,
        "csn_max": csn_max,
        "csp_max": csp_max,
        "temperature": temperature,
        "x_legacy": x_legacy,
        "un_legacy": un_legacy,
        "y_exp": y_exp,
        "up_exp": up_exp,
        "soc": soc,
        "v_qocv": v_qocv,
        "q_meas": q_meas,
        "qn_cell": qn_cell,
        "qp_cell": qp_cell,
        "delta_x": delta_x,
        "delta_y": delta_y,
        "dynamic": bundle["dynamic"],
        "q_meas_bundle": prior("q_meas_Ah"),
    }


def solve_stage1_case(
    case_name: str,
    x_ocp: np.ndarray,
    un_ocp: np.ndarray,
    model_inputs: dict[str, Any],
) -> tuple[dict[str, Any], np.ndarray]:
    soc = model_inputs["soc"]
    v_qocv = model_inputs["v_qocv"]
    y_exp = model_inputs["y_exp"]
    up_exp = model_inputs["up_exp"]
    delta_x = float(model_inputs["delta_x"])
    delta_y = float(model_inputs["delta_y"])

    def un(x):
        return np.interp(x, x_ocp, un_ocp)

    def up(y):
        return np.interp(y, y_exp, up_exp)

    x0_lower = max(1e-6, float(x_ocp.min()))
    x0_upper = min(1.0 - delta_x - 1e-6, float(x_ocp.max()) - delta_x)
    y100_lower = max(1e-6, float(y_exp.min()))
    y100_upper = min(1.0 - delta_y - 1e-6, float(y_exp.max()) - delta_y)
    if x0_upper <= x0_lower or y100_upper <= y100_lower:
        raise RuntimeError(f"{case_name}: experimental OCP coverage cannot cover Stage 1 window")

    def residuals(offsets):
        x0, y100 = offsets
        x100 = x0 + delta_x
        y0 = y100 + delta_y
        return np.array(
            [
                up(y0) - un(x0) - v_qocv[0],
                up(y100) - un(x100) - v_qocv[-1],
            ]
        )

    initial = np.array([(x0_lower + x0_upper) / 2.0, (y100_lower + y100_upper) / 2.0])
    result = least_squares(
        residuals,
        x0=initial,
        bounds=([x0_lower, y100_lower], [x0_upper, y100_upper]),
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
    )
    x0, y100 = result.x
    x100 = x0 + delta_x
    y0 = y100 + delta_y
    voltage = up(y0 - soc * (y0 - y100)) - un(x0 + soc * (x100 - x0))
    error_mV = (voltage - v_qocv) * 1000.0
    central = (soc >= 0.02) & (soc <= 0.98)
    endpoint_error = residuals(result.x) * 1000.0
    tolerance = 1e-5
    bound_names: list[str] = []
    for name, value, lower, upper in (
        ("x0", x0, x0_lower, x0_upper),
        ("y100", y100, y100_lower, y100_upper),
    ):
        if value - lower <= tolerance:
            bound_names.append(f"{name}_lower")
        if upper - value <= tolerance:
            bound_names.append(f"{name}_upper")

    model_span = x100 - x0
    outside_span = max(0.0, float(x_ocp.min()) - x0) + max(0.0, x100 - float(x_ocp.max()))
    row = {
        "Case": case_name,
        "x0": float(x0),
        "x100": float(x100),
        "y100": float(y100),
        "y0": float(y0),
        "Delta_x": delta_x,
        "Delta_y": delta_y,
        "Q_MEAS_Ah": float(model_inputs["q_meas"]),
        "Qn_Ah": float(model_inputs["qn_cell"]),
        "Qp_Ah": float(model_inputs["qp_cell"]),
        "SOC_0_endpoint_error_mV": float(endpoint_error[0]),
        "SOC_100_endpoint_error_mV": float(endpoint_error[1]),
        "qOCV_RMSE_full_mV": float(np.sqrt(np.mean(error_mV**2))),
        "qOCV_MAE_full_mV": float(np.mean(np.abs(error_mV))),
        "qOCV_RMSE_2_98_mV": float(np.sqrt(np.mean(error_mV[central] ** 2))),
        "qOCV_MAE_2_98_mV": float(np.mean(np.abs(error_mV[central]))),
        "negative_OCP_measured_min": float(x_ocp.min()),
        "negative_OCP_measured_max": float(x_ocp.max()),
        "negative_OCP_model_required_min": float(x0),
        "negative_OCP_model_required_max": float(x100),
        "stage1_extrapolated_fraction": float(outside_span / model_span),
        "stage1_extrapolation_side": (
            "low-x" if x0 < x_ocp.min() else "high-x" if x100 > x_ocp.max() else "none"
        ),
        "physical_bounds_satisfied": bool(0 <= x0 < x100 <= 1 and 0 <= y100 < y0 <= 1),
        "optimizer_boundary_hit": bool(bound_names),
        "optimizer_boundary": ";".join(bound_names) if bound_names else "none",
        "Solver_success": bool(result.success),
        "Solver_message": str(result.message),
    }
    return row, voltage


def pybamm_interp(grid, values, name):
    grid = np.asarray(grid, dtype=float)
    values = np.asarray(values, dtype=float)

    def ocp(sto):
        return pybamm.Interpolant(
            grid,
            values,
            [sto],
            name=name,
            interpolator="linear",
            extrapolate=True,
        )

    return ocp


def load_dynamic_data(frame: pd.DataFrame) -> dict[tuple[float, bool], dict[str, np.ndarray]]:
    output: dict[tuple[float, bool], dict[str, np.ndarray]] = {}
    for rate in RATES:
        for charge in (True, False):
            direction = "charge" if charge else "discharge"
            selected = frame[
                np.isclose(pd.to_numeric(frame["rate_C"], errors="coerce"), rate)
                & (frame["direction"].astype(str).str.lower() == direction)
            ].copy()
            if len(selected) <= 5:
                raise RuntimeError(f"Dynamic data missing: {rate}C {direction}")
            output[(rate, charge)] = {
                "t_min": pd.to_numeric(selected["t_min"], errors="coerce").to_numpy(dtype=float),
                "V": pd.to_numeric(selected["V_V"], errors="coerce").to_numpy(dtype=float),
                "SOC": pd.to_numeric(selected["SOC"], errors="coerce").to_numpy(dtype=float),
            }
    return output


def make_dfn_parameter_values(
    model_inputs: dict[str, Any],
    x_init: float,
    y_init: float,
    x_ocp: np.ndarray,
    un_ocp: np.ndarray,
):
    base = model_inputs["params"].copy()
    values = base.copy()
    values.update(
        {
            "Initial concentration in negative electrode [mol.m-3]": float(x_init)
            * model_inputs["csn_max"],
            "Initial concentration in positive electrode [mol.m-3]": float(y_init)
            * model_inputs["csp_max"],
            "Negative electrode OCP [V]": pybamm_interp(
                x_ocp, un_ocp, "Experimental negative OCP"
            ),
            "Positive electrode OCP [V]": pybamm_interp(
                model_inputs["y_exp"], model_inputs["up_exp"], "Experimental positive OCP"
            ),
            "Initial temperature [K]": model_inputs["temperature"],
            "Ambient temperature [K]": model_inputs["temperature"],
        },
        check_already_exists=False,
    )
    return values


def _surface_stoichiometry_at_separator(solution, n_time: int) -> np.ndarray:
    entries = np.asarray(solution["Negative particle surface stoichiometry"].entries, dtype=float)
    entries = np.squeeze(entries)
    if entries.ndim == 1:
        return entries
    if entries.shape[-1] == n_time:
        return entries[-1, :]
    if entries.shape[0] == n_time:
        return entries[:, -1]
    raise RuntimeError(f"Unexpected surface stoichiometry shape: {entries.shape}")


def run_cc_dfn(
    rate: float,
    charge: bool,
    stage1: dict[str, Any],
    x_ocp: np.ndarray,
    un_ocp: np.ndarray,
    model_inputs: dict[str, Any],
) -> dict[str, np.ndarray]:
    x_init = stage1["x0"] if charge else stage1["x100"]
    y_init = stage1["y0"] if charge else stage1["y100"]
    model = pybamm.lithium_ion.DFN({"thermal": "isothermal"})
    model.variables["Anode potential [V]"] = model.variables[
        "Negative electrode surface potential difference at separator interface [V]"
    ]
    step = (
        f"Charge at {rate}C until 4.2 V"
        if charge
        else f"Discharge at {rate}C until 3.0 V"
    )
    simulation = pybamm.Simulation(
        model,
        parameter_values=make_dfn_parameter_values(
            model_inputs, x_init, y_init, x_ocp, un_ocp
        ),
        experiment=pybamm.Experiment([step], period="10 seconds"),
        var_pts={"x_n": 15, "x_s": 15, "x_p": 15, "r_n": 15, "r_p": 15},
        solver=pybamm.IDAKLUSolver(rtol=1e-6, atol=1e-8),
    )
    solution = simulation.solve()
    time_s = np.asarray(solution["Time [s]"].entries, dtype=float)
    voltage = np.asarray(solution["Terminal voltage [V]"].entries, dtype=float)
    current = np.asarray(solution["Current [A]"].entries, dtype=float)
    anode_potential = np.asarray(solution["Anode potential [V]"].entries, dtype=float).squeeze()
    capacity_abs = cumulative_trapezoid(np.abs(current), time_s, initial=0.0) / 3600.0
    soc = capacity_abs / model_inputs["q_meas"] if charge else 1.0 - capacity_abs / model_inputs["q_meas"]
    x_surface = _surface_stoichiometry_at_separator(solution, len(time_s))
    return {
        "t_min": time_s / 60.0,
        "V": voltage,
        "SOC": np.clip(soc, 0.0, 1.0),
        "Anode_potential_V": anode_potential,
        "x_surface_separator": x_surface,
    }


def _time_fraction(t_min: np.ndarray, condition: np.ndarray) -> float:
    if len(t_min) < 2 or t_min[-1] <= t_min[0]:
        return np.nan
    weights = np.diff(t_min)
    segment_fraction = 0.5 * (condition[:-1].astype(float) + condition[1:].astype(float))
    return float(np.sum(weights * segment_fraction) / (t_min[-1] - t_min[0]))


def _dynamic_error_metrics(experiment: dict[str, np.ndarray], simulation: dict[str, np.ndarray]):
    exp_soc, exp_voltage = _clean_curve(experiment["SOC"], experiment["V"])
    sim_soc, sim_voltage = _clean_curve(simulation["SOC"], simulation["V"])
    overlap = (exp_soc >= sim_soc.min()) & (exp_soc <= sim_soc.max())
    model_on_exp = np.interp(exp_soc[overlap], sim_soc, sim_voltage)
    error_mV = (model_on_exp - exp_voltage[overlap]) * 1000.0
    soc_overlap = exp_soc[overlap]
    central = (soc_overlap >= 0.10) & (soc_overlap <= 0.90)
    if not central.any():
        central_error = np.array([np.nan])
    else:
        central_error = error_mV[central]
    return {
        "overlap_SOC_min": float(soc_overlap.min()),
        "overlap_SOC_max": float(soc_overlap.max()),
        "overlap_point_count": int(len(error_mV)),
        "RMSE_full_overlap_mV": float(np.sqrt(np.mean(error_mV**2))),
        "MAE_full_overlap_mV": float(np.mean(np.abs(error_mV))),
        "RMSE_SOC_10_90_mV": float(np.sqrt(np.nanmean(central_error**2))),
        "MAE_SOC_10_90_mV": float(np.nanmean(np.abs(central_error))),
        "end_time_error_min": float(simulation["t_min"][-1] - experiment["t_min"][-1]),
    }


def run_dynamic_comparison(
    cases: dict[str, dict[str, Any]],
    model_inputs: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[tuple[str, float, bool], dict[str, np.ndarray]]]:
    experiment = load_dynamic_data(model_inputs["dynamic"])
    dynamic_rows: list[dict[str, Any]] = []
    anode_rows: list[dict[str, Any]] = []
    coverage_rows: list[dict[str, Any]] = []
    simulations: dict[tuple[str, float, bool], dict[str, np.ndarray]] = {}

    for case_name, case in cases.items():
        for rate in RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                print(f"DFN {case_name}: {rate}C {direction}", flush=True)
                base = {
                    "Case": case_name,
                    "C_rate": rate,
                    "direction": direction,
                }
                try:
                    simulation = run_cc_dfn(
                        rate,
                        charge,
                        case["stage1"],
                        case["x_ocp"],
                        case["un_ocp"],
                        model_inputs,
                    )
                    simulations[(case_name, rate, charge)] = simulation
                    metrics = _dynamic_error_metrics(experiment[(rate, charge)], simulation)
                    dynamic_rows.append({**base, **metrics, "simulation_success": True, "failure_message": ""})

                    x_surface = simulation["x_surface_separator"]
                    x_min = float(np.min(x_surface))
                    x_max = float(np.max(x_surface))
                    measured_min = float(case["x_ocp"].min())
                    measured_max = float(case["x_ocp"].max())
                    outside = (x_surface < measured_min) | (x_surface > measured_max)
                    coverage_rows.append(
                        {
                            **base,
                            "x_surf_min": x_min,
                            "x_surf_max": x_max,
                            "measured_OCP_min": measured_min,
                            "measured_OCP_max": measured_max,
                            "inside_measured_OCP_coverage": bool(not outside.any()),
                            "measured_upper_limit_excess": max(0.0, x_max - measured_max),
                            "measured_lower_limit_shortfall": max(0.0, measured_min - x_min),
                            "out_of_range_time_fraction": _time_fraction(simulation["t_min"], outside),
                            "simulation_success": True,
                            "failure_message": "",
                        }
                    )

                    if charge:
                        potential = simulation["Anode_potential_V"]
                        minimum_index = int(np.argmin(potential))
                        below_zero = potential < 0.0
                        anode_rows.append(
                            {
                                **base,
                                "minimum_anode_potential_V": float(potential[minimum_index]),
                                "minimum_time_min": float(simulation["t_min"][minimum_index]),
                                "minimum_SOC": float(simulation["SOC"][minimum_index]),
                                "negative_surface_stoichiometry_at_minimum": float(x_surface[minimum_index]),
                                "anode_potential_below_0V": bool(below_zero.any()),
                                "below_0V_time_fraction": _time_fraction(simulation["t_min"], below_zero),
                                "simulation_success": True,
                                "failure_message": "",
                            }
                        )
                except Exception as error:
                    message = f"{type(error).__name__}: {error}"
                    print(f"FAILED {case_name} {rate}C {direction}: {message}", flush=True)
                    dynamic_rows.append(
                        {
                            **base,
                            "overlap_SOC_min": np.nan,
                            "overlap_SOC_max": np.nan,
                            "overlap_point_count": 0,
                            "RMSE_full_overlap_mV": np.nan,
                            "MAE_full_overlap_mV": np.nan,
                            "RMSE_SOC_10_90_mV": np.nan,
                            "MAE_SOC_10_90_mV": np.nan,
                            "end_time_error_min": np.nan,
                            "simulation_success": False,
                            "failure_message": message,
                        }
                    )
                    coverage_rows.append(
                        {
                            **base,
                            "x_surf_min": np.nan,
                            "x_surf_max": np.nan,
                            "measured_OCP_min": float(case["x_ocp"].min()),
                            "measured_OCP_max": float(case["x_ocp"].max()),
                            "inside_measured_OCP_coverage": False,
                            "measured_upper_limit_excess": np.nan,
                            "measured_lower_limit_shortfall": np.nan,
                            "out_of_range_time_fraction": np.nan,
                            "simulation_success": False,
                            "failure_message": message,
                        }
                    )
                    if charge:
                        anode_rows.append(
                            {
                                **base,
                                "minimum_anode_potential_V": np.nan,
                                "minimum_time_min": np.nan,
                                "minimum_SOC": np.nan,
                                "negative_surface_stoichiometry_at_minimum": np.nan,
                                "anode_potential_below_0V": False,
                                "below_0V_time_fraction": np.nan,
                                "simulation_success": False,
                                "failure_message": message,
                            }
                        )
    return (
        pd.DataFrame(dynamic_rows),
        pd.DataFrame(anode_rows),
        pd.DataFrame(coverage_rows),
        simulations,
    )


def create_ocp_comparison_table(
    legacy: dict[str, Any], new_ocp: dict[str, Any]
) -> pd.DataFrame:
    grid = np.linspace(0.0, 1.0, 3001)
    legacy_normalized_x = (
        (legacy["x_legacy"] - legacy["x_legacy"].min())
        / (legacy["x_legacy"].max() - legacy["x_legacy"].min())
    )
    output = {
        "normalized_capacity": grid,
        "legacy_low_rate_OCP_V": np.interp(grid, legacy_normalized_x, legacy["un_legacy"]),
    }
    for replicate in ("v2", "v3"):
        for branch in ("lithiation", "delithiation", "mean"):
            curve = new_ocp["curves"][replicate][branch]
            output[f"{replicate}_{branch}_OCP_V"] = np.interp(grid, curve[0], curve[1])
    representative = new_ocp["representative"]
    output["new_GITT_representative_OCP_V"] = np.interp(grid, representative[0], representative[1])
    output["v2_minus_v3_mean_mV"] = (
        output["v2_mean_OCP_V"] - output["v3_mean_OCP_V"]
    ) * 1000.0
    output["new_minus_legacy_mV"] = (
        output["new_GITT_representative_OCP_V"] - output["legacy_low_rate_OCP_V"]
    ) * 1000.0
    return pd.DataFrame(output)


def save_figures(
    paths: AnalysisPaths,
    gitt: dict[str, Any],
    new_ocp: dict[str, Any],
    legacy: dict[str, Any],
    stage1_voltage: dict[str, np.ndarray],
    simulations: dict[tuple[str, float, bool], dict[str, np.ndarray]],
):
    plt.rcParams.update({"figure.figsize": (8.0, 5.0), "axes.grid": True, "grid.alpha": 0.25})
    results = paths.results

    fig, axis = plt.subplots(figsize=(9.2, 5.8))
    legacy_x_norm = (legacy["x_legacy"] - legacy["x_legacy"].min()) / (
        legacy["x_legacy"].max() - legacy["x_legacy"].min()
    )
    axis.plot(legacy_x_norm, legacy["un_legacy"], color="black", lw=2.2, label="Legacy low-rate OCP")
    colors = {"v2": "#1f77b4", "v3": "#ff7f0e"}
    for replicate in ("v2", "v3"):
        for branch, linestyle in (("lithiation", "--"), ("delithiation", ":")):
            curve = new_ocp["curves"][replicate][branch]
            axis.plot(
                curve[0],
                curve[1],
                color=colors[replicate],
                linestyle=linestyle,
                alpha=0.65,
                label=f"{replicate} {branch}",
            )
        curve = new_ocp["curves"][replicate]["mean"]
        axis.plot(curve[0], curve[1], color=colors[replicate], lw=1.8, label=f"{replicate} branch mean")
    axis.plot(
        new_ocp["representative"][0],
        new_ocp["representative"][1],
        color="#d62728",
        lw=2.5,
        label="New GITT v2/v3 mean",
    )
    axis.set(xlabel="Normalized capacity / candidate x", ylabel="Anode OCP [V vs. Li/Li+]", title="Legacy low-rate Anode OCP vs New GITT Anode OCP")
    axis.legend(ncol=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(results / "anode_ocp_legacy_vs_new.png", dpi=180)
    plt.close(fig)


def write_summary_report(
    paths: AnalysisPaths,
    qc_summary: pd.DataFrame,
    capacity_summary: pd.DataFrame,
    reproducibility: pd.DataFrame,
    branch_hysteresis: pd.DataFrame,
    stage1: pd.DataFrame,
    dynamic: pd.DataFrame,
    anode: pd.DataFrame,
    coverage: pd.DataFrame,
    gitt: dict[str, Any],
    new_ocp: dict[str, Any],
    legacy: dict[str, Any],
    stage1_voltage: dict[str, np.ndarray],
    simulations: dict[tuple[str, float, bool], dict[str, np.ndarray]],
) -> None:
    results = paths.results
    def qc_row(replicate):
        return qc_summary[
            (qc_summary["replicate"] == replicate)
            & (qc_summary["branch"] == "all_GITT_rests")
        ].iloc[0]

    def capacity(replicate, kind):
        return float(
            capacity_summary[
                (capacity_summary["replicate"] == replicate)
                & (capacity_summary["capacity_type"] == kind)
            ].iloc[0]["capacity_Ah"]
        )

    v2_lith = capacity("v2", "gitt_lithiation")
    v3_lith = capacity("v3", "gitt_lithiation")
    v2_delith = capacity("v2", "gitt_delithiation")
    v3_delith = capacity("v3", "gitt_delithiation")
    lith_diff = abs(v2_lith - v3_lith) / np.mean([v2_lith, v3_lith]) * 100.0
    delith_diff = abs(v2_delith - v3_delith) / np.mean([v2_delith, v3_delith]) * 100.0
    rep_full = reproducibility[
        (reproducibility["branch"] == "mean")
        & (reproducibility["region"] == "full_overlap")
    ].iloc[0]
    rep_mid = reproducibility[
        (reproducibility["branch"] == "mean")
        & (reproducibility["region"] == "normalized_capacity_10_90")
    ].iloc[0]
    legacy_stage = stage1[stage1["Case"] == "Legacy"].iloc[0]
    new_stage = stage1[stage1["Case"] == "New GITT"].iloc[0]
    legacy_success = dynamic[dynamic["Case"] == "Legacy"]["simulation_success"].sum()
    new_charge = dynamic[(dynamic["Case"] == "New GITT") & (dynamic["direction"] == "Charge")]
    new_discharge = dynamic[(dynamic["Case"] == "New GITT") & (dynamic["direction"] == "Discharge")]

    lines = [
        "# 260918 GITT OCP Stage 1 comparison",
        "",
        "## 1. Data QC",
        "",
    ]
    for replicate in ("v2", "v3"):
        row = qc_row(replicate)
        lines.append(
            f"- {replicate}: full experiment {row['experiment_duration_h']:.2f} h; "
            f"GITT section {row['gitt_duration_h']:.2f} h; "
            f"{int(row['pulse_count'])} pulses and {int(row['rest_count'])} two-hour rests. "
            f"Last-10-min |dV/dt| median/mean/P90/max = "
            f"{row['median_abs_dVdt_mV_per_min']:.4f}/"
            f"{row['mean_abs_dVdt_mV_per_min']:.4f}/"
            f"{row['P90_abs_dVdt_mV_per_min']:.4f}/"
            f"{row['maximum_abs_dVdt_mV_per_min']:.4f} mV/min."
        )
    lines.extend(
        [
            f"- GITT usable capacity: v2 lithiation/delithiation = {v2_lith*1000:.3f}/{v2_delith*1000:.3f} mAh; "
            f"v3 = {v3_lith*1000:.3f}/{v3_delith*1000:.3f} mAh. "
            f"Replicate differences are {lith_diff:.2f}% and {delith_diff:.2f}%.",
            f"- Replicate mean-OCP difference: full-range MAE/RMSE = {rep_full['MAE_mV']:.2f}/{rep_full['RMSE_mV']:.2f} mV; "
            f"10–90% MAE/RMSE = {rep_mid['MAE_mV']:.2f}/{rep_mid['RMSE_mV']:.2f} mV. "
            f"The full-range maximum is {rep_full['maximum_abs_difference_mV']:.1f} mV and is endpoint-driven.",
            "- The legacy bundle contains only aggregated 1 h GITT rest-end points, not the underlying final-10-min record data. A like-for-like 1 h vs 2 h |dV/dt| comparison is therefore not available.",
            "",
            "## 2. OCP comparison",
            "",
        ]
    )
    for replicate in ("v2", "v3"):
        h = branch_hysteresis[
            (branch_hysteresis["replicate"] == replicate)
            & (branch_hysteresis["region"] == "normalized_capacity_10_90")
        ].iloc[0]
        lines.append(
            f"- {replicate} lithiation-delithiation hysteresis over 10–90%: "
            f"MAE {h['MAE_hysteresis_mV']:.2f} mV, RMSE {h['RMSE_hysteresis_mV']:.2f} mV, "
            f"maximum {h['maximum_abs_hysteresis_mV']:.2f} mV."
        )
    lines.extend(
        [
            "- v2 and v3 were normalized by their own usable branch capacity before interpolation and averaging. Raw GITT capacity was not converted directly to absolute stoichiometry.",
            "- For this first-pass candidate, the normalized usable range was mapped to candidate x=0..1. This improves nominal measured coverage but does not establish the absolute stoichiometry anchor.",
            "",
            "## 3. Stage 1 comparison",
            "",
            f"- Legacy: x0/x100 = {legacy_stage['x0']:.6f}/{legacy_stage['x100']:.6f}; "
            f"y100/y0 = {legacy_stage['y100']:.6f}/{legacy_stage['y0']:.6f}; "
            f"qOCV RMSE full/2–98% = {legacy_stage['qOCV_RMSE_full_mV']:.2f}/{legacy_stage['qOCV_RMSE_2_98_mV']:.2f} mV.",
            f"- New GITT: x0/x100 = {new_stage['x0']:.6f}/{new_stage['x100']:.6f}; "
            f"y100/y0 = {new_stage['y100']:.6f}/{new_stage['y0']:.6f}; "
            f"qOCV RMSE full/2–98% = {new_stage['qOCV_RMSE_full_mV']:.2f}/{new_stage['qOCV_RMSE_2_98_mV']:.2f} mV.",
            f"- New GITT selected the x0 lower bound and left endpoint errors of "
            f"{new_stage['SOC_0_endpoint_error_mV']:.2f} mV (SOC 0) and "
            f"{new_stage['SOC_100_endpoint_error_mV']:.2f} mV (SOC 100).",
            "",
            "## 4. Dynamic baseline",
            "",
            f"- Legacy completed {int(legacy_success)}/6 nominal DFN cases.",
            f"- New GITT completed {int(new_discharge['simulation_success'].sum())}/3 discharge cases and "
            f"{int(new_charge['simulation_success'].sum())}/3 charge cases. All New GITT charge cases failed at initialization with x0 near zero.",
        ]
    )
    for rate in RATES:
        old = dynamic[
            (dynamic["Case"] == "Legacy")
            & (dynamic["direction"] == "Discharge")
            & np.isclose(dynamic["C_rate"], rate)
        ].iloc[0]
        new = dynamic[
            (dynamic["Case"] == "New GITT")
            & (dynamic["direction"] == "Discharge")
            & np.isclose(dynamic["C_rate"], rate)
        ].iloc[0]
        lines.append(
            f"- {rate:g}C discharge full-overlap RMSE: Legacy {old['RMSE_full_overlap_mV']:.2f} mV; "
            f"New GITT {new['RMSE_full_overlap_mV']:.2f} mV."
        )
    lines.extend(["", "## 5. Anode potential", ""])
    for _, row in anode[anode["Case"] == "Legacy"].iterrows():
        lines.append(
            f"- Legacy {row['C_rate']:g}C charge: minimum {row['minimum_anode_potential_V']:.4f} V at "
            f"SOC {row['minimum_SOC']:.3f} and {row['minimum_time_min']:.2f} min; "
            f"time fraction below 0 V = {row['below_0V_time_fraction']:.3f}."
        )
    lines.extend(
        [
            "- New GITT charge anode-potential metrics are unavailable because all charge simulations failed at initialization.",
            "- A modeled anode potential below 0 V is treated only as a DFN electrochemical-state/risk indicator, not proof of lithium plating.",
            "",
            "## 6. Recommendation",
            "",
            "**Do not adopt yet**",
            "",
            "The 2 h GITT experiment has good central-range replicate reproducibility and the New GITT OCP improves the three discharge baseline RMSE values. However, it worsens Stage 1 qOCV error, drives x0 to the lower bound, fails to satisfy the SOC 0 endpoint closely, and makes all nominal charge DFN simulations fail at initialization. The absolute stoichiometry anchor of the normalized capacity axis also remains unverified.",
            "",
            "Before Stage 2–4 is considered:",
            "",
            "1. establish an independent absolute stoichiometry/capacity anchor for the new half-cells;",
            "2. investigate the low-x endpoint mismatch and endpoint replicate divergence without tuning dynamic parameters;",
            "3. resolve the x0 boundary solution and reproduce all three charge baselines;",
            "4. then reassess surface-stoichiometry coverage and anode potential for New GITT.",
            "",
            "Stage 2–4 fitting was not run.",
        ]
    )
    (paths.results / "comparison_report.md").write_text("\n".join(lines), encoding="utf-8")

    fig, axis = plt.subplots(figsize=(8.5, 5.2))
    for branch, linestyle in (("lithiation", "--"), ("delithiation", ":"), ("mean", "-")):
        v2_curve = new_ocp["curves"]["v2"][branch]
        v3_curve = new_ocp["curves"]["v3"][branch]
        lower = max(v2_curve[0].min(), v3_curve[0].min())
        upper = min(v2_curve[0].max(), v3_curve[0].max())
        grid = np.linspace(lower, upper, 1001)
        delta = (_interp_curve(v2_curve, grid) - _interp_curve(v3_curve, grid)) * 1000.0
        axis.plot(grid, delta, linestyle=linestyle, label=branch)
    axis.axhline(0.0, color="black", lw=0.8)
    axis.axvspan(0.0, 0.1, color="#eeeeee", alpha=0.7)
    axis.axvspan(0.9, 1.0, color="#eeeeee", alpha=0.7)
    axis.set(xlabel="Normalized capacity / candidate x", ylabel="v2 - v3 [mV]", title="New GITT v2-v3 difference")
    axis.legend()
    fig.tight_layout()
    fig.savefig(results / "anode_gitt_v2_v3_reproducibility.png", dpi=180)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(9.0, 5.4))
    for (replicate, branch), frame in gitt["rests"].groupby(["replicate", "branch"]):
        plot_values = np.maximum(
            frame["last_10_min_abs_dVdt_mV_per_min"].to_numpy(dtype=float),
            1e-4,
        )
        axis.plot(
            frame["q_norm"],
            plot_values,
            marker="o",
            ms=3,
            lw=1,
            label=f"{replicate} {branch}",
        )
    axis.set_yscale("log")
    axis.set_ylim(1e-4, 1.0)
    axis.set(
        xlabel="Normalized branch capacity",
        ylabel="Last-10-min |dV/dt| [mV/min] (values <1e-4 clipped)",
        title="Rest equilibrium quality",
    )
    axis.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(results / "gitt_rest_equilibrium.png", dpi=180)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(8.5, 5.2))
    axis.plot(legacy["soc"], legacy["v_qocv"], color="black", lw=2.5, label="Measured C/20 qOCV")
    axis.plot(legacy["soc"], stage1_voltage["Legacy"], lw=1.8, label="Legacy OCP")
    axis.plot(legacy["soc"], stage1_voltage["New GITT"], lw=1.8, label="New GITT OCP")
    axis.set(xlabel="SOC", ylabel="Cell OCV [V]", title="Stage 1 qOCV comparison")
    axis.legend()
    fig.tight_layout()
    fig.savefig(results / "qocv_stage1_comparison.png", dpi=180)
    plt.close(fig)

    experiment = load_dynamic_data(legacy["dynamic"])
    fig, axes = plt.subplots(2, 3, figsize=(14.0, 8.0), sharex=True)
    for row, charge in enumerate((True, False)):
        for column, rate in enumerate(RATES):
            axis = axes[row, column]
            exp = experiment[(rate, charge)]
            axis.plot(exp["SOC"], exp["V"], color="black", lw=2, label="Experiment")
            for case_name, color in (("Legacy", "#1f77b4"), ("New GITT", "#d62728")):
                simulation = simulations.get((case_name, rate, charge))
                if simulation is not None:
                    axis.plot(simulation["SOC"], simulation["V"], color=color, lw=1.4, label=case_name)
                elif case_name == "New GITT":
                    axis.text(
                        0.04,
                        0.06,
                        "New GITT: solver failure",
                        color=color,
                        fontsize=8,
                        transform=axis.transAxes,
                    )
            axis.set_title(f"{rate}C {'Charge' if charge else 'Discharge'}")
            axis.set_xlabel("SOC")
            axis.set_ylabel("Terminal voltage [V]")
    axes[0, 0].legend(loc="best", fontsize=8)
    fig.suptitle("Dynamic voltage baseline comparison", y=0.985)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(results / "dynamic_voltage_baseline_comparison.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.5), sharey=True)
    for column, rate in enumerate(RATES):
        axis = axes[column]
        for case_name, color in (("Legacy", "#1f77b4"), ("New GITT", "#d62728")):
            simulation = simulations.get((case_name, rate, True))
            if simulation is not None:
                axis.plot(simulation["SOC"], simulation["Anode_potential_V"], color=color, lw=1.6, label=case_name)
            elif case_name == "New GITT":
                axis.text(
                    0.04,
                    0.06,
                    "New GITT: solver failure",
                    color=color,
                    fontsize=8,
                    transform=axis.transAxes,
                )
        axis.axhline(0.0, color="black", lw=0.9, linestyle="--")
        axis.set_title(f"{rate}C Charge")
        axis.set_xlabel("SOC")
        axis.set_ylabel("Anode potential [V vs. Li/Li+]")
    axes[0].legend(loc="best", fontsize=8)
    fig.suptitle("Anode potential comparison", y=0.985)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(results / "anode_potential_comparison.png", dpi=180)
    plt.close(fig)


def run_analysis(root: Path | str | None = None, run_dfn: bool = True) -> dict[str, Any]:
    started = time.time()
    paths = get_paths(root)
    paths.results.mkdir(parents=True, exist_ok=True)
    print(f"Python {platform.python_version()} | PyBaMM {pybamm.__version__}", flush=True)
    print(f"Repository: {paths.root}", flush=True)
    print(f"Results: {paths.results}", flush=True)

    gitt = load_new_gitt(paths)
    new_ocp = build_new_ocp(gitt)
    reproducibility = reproducibility_table(gitt, new_ocp)
    branch_hysteresis = branch_hysteresis_table(new_ocp)
    qc_summary, capacity_summary = data_qc_tables(gitt)
    bundle = load_legacy_bundle(paths)
    legacy = build_legacy_model_inputs(bundle)
    ocp_comparison = create_ocp_comparison_table(legacy, new_ocp)

    stage1_rows: list[dict[str, Any]] = []
    stage1_voltage: dict[str, np.ndarray] = {}
    cases: dict[str, dict[str, Any]] = {}
    for case_name, x_ocp, un_ocp in (
        ("Legacy", legacy["x_legacy"], legacy["un_legacy"]),
        ("New GITT", new_ocp["representative"][0], new_ocp["representative"][1]),
    ):
        row, voltage = solve_stage1_case(case_name, x_ocp, un_ocp, legacy)
        stage1_rows.append(row)
        stage1_voltage[case_name] = voltage
        cases[case_name] = {"stage1": row, "x_ocp": x_ocp, "un_ocp": un_ocp}
    stage1 = pd.DataFrame(stage1_rows)

    reproducibility.to_csv(paths.results / "gitt_reproducibility.csv", index=False, encoding="utf-8-sig")
    branch_hysteresis.to_csv(paths.results / "gitt_branch_hysteresis.csv", index=False, encoding="utf-8-sig")
    gitt["rests"].to_csv(paths.results / "gitt_rest_equilibrium.csv", index=False, encoding="utf-8-sig")
    qc_summary.to_csv(paths.results / "gitt_data_qc_summary.csv", index=False, encoding="utf-8-sig")
    capacity_summary.to_csv(paths.results / "gitt_capacity_summary.csv", index=False, encoding="utf-8-sig")
    ocp_comparison.to_csv(paths.results / "legacy_vs_new_ocp.csv", index=False, encoding="utf-8-sig")
    stage1.to_csv(paths.results / "stage1_stoichiometry_comparison.csv", index=False, encoding="utf-8-sig")

    dynamic = pd.DataFrame()
    anode = pd.DataFrame()
    coverage = pd.DataFrame()
    simulations: dict[tuple[str, float, bool], dict[str, np.ndarray]] = {}
    if run_dfn:
        dynamic, anode, coverage, simulations = run_dynamic_comparison(cases, legacy)
        dynamic.to_csv(paths.results / "dynamic_baseline_comparison.csv", index=False, encoding="utf-8-sig")
        anode.to_csv(paths.results / "anode_potential_comparison.csv", index=False, encoding="utf-8-sig")
        coverage.to_csv(paths.results / "surface_stoichiometry_coverage.csv", index=False, encoding="utf-8-sig")
    else:
        for filename, columns in (
            ("dynamic_baseline_comparison.csv", []),
            ("anode_potential_comparison.csv", []),
            ("surface_stoichiometry_coverage.csv", []),
        ):
            pd.DataFrame(columns=columns).to_csv(paths.results / filename, index=False, encoding="utf-8-sig")

    save_figures(paths, gitt, new_ocp, legacy, stage1_voltage, simulations)
    if run_dfn:
        write_summary_report(
            paths,
            qc_summary,
            capacity_summary,
            reproducibility,
            branch_hysteresis,
            stage1,
            dynamic,
            anode,
            coverage,
            gitt,
            new_ocp,
            legacy,
            stage1_voltage,
            simulations,
        )
    manifest = {
        "random_seed": RANDOM_SEED,
        "python_version": platform.python_version(),
        "pybamm_version": pybamm.__version__,
        "source_branch": "260918-gitt-ocp-test",
        "legacy_notebook": "260908 Charge Fitting Discharge Validation.ipynb",
        "negative_ocp_cases": ["Legacy", "New GITT"],
        "fixed_items": [
            "Full-cell C/20 qOCV",
            "Q_MEAS",
            "Delta x and Delta y",
            "1A endpoint equations",
            "Stage 0 source corrections",
            "DFN solver and mesh",
            "temperature",
            "nominal dynamic parameters",
        ],
        "new_gitt_absolute_stoichiometry_note": (
            "New GITT capacity is normalized independently by branch and replicate. "
            "The candidate OCP maps the full normalized usable range to x=0..1; "
            "the absolute stoichiometry anchor remains unverified."
        ),
        "elapsed_min": (time.time() - started) / 60.0,
    }
    (paths.results / "analysis_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"Analysis completed in {manifest['elapsed_min']:.1f} min", flush=True)
    return {
        "paths": paths,
        "gitt": gitt,
        "new_ocp": new_ocp,
        "legacy": legacy,
        "reproducibility": reproducibility,
        "branch_hysteresis": branch_hysteresis,
        "qc_summary": qc_summary,
        "capacity_summary": capacity_summary,
        "stage1": stage1,
        "dynamic": dynamic,
        "anode": anode,
        "coverage": coverage,
        "simulations": simulations,
        "manifest": manifest,
    }


if __name__ == "__main__":
    run_analysis()
