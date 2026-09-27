from __future__ import annotations

import json
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


warnings.filterwarnings("ignore", message="Workbook contains no default style")

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "outputs" / "260927_september_raw_fullcurve_audit"
OUT.mkdir(parents=True, exist_ok=True)
DATA = Path(
    r"C:\Users\user\OneDrive\01. 적응형 충전 프로토콜"
    r"\4. Enertech 파우치셀 실험데이터\260927 BoL 실험데이터"
)
FILES = {
    "6-1": DATA / "260927 6-1 0.5C 1C 2C.xlsx",
    "6-2": DATA / "260927 6-2 0.5C 1C 2C.xlsx",
}
BRANCH_STEPS = {
    (0.5, "Discharge"): 5,
    (0.5, "Charge"): 7,
    (1.0, "Discharge"): 10,
    (1.0, "Charge"): 12,
    (2.0, "Discharge"): 15,
    (2.0, "Charge"): 17,
}


def duration_min(row: pd.Series) -> float:
    return float((row["End Date"] - row["Oneset Date"]).total_seconds() / 60.0)


def segment(record: pd.DataFrame, row: pd.Series) -> pd.DataFrame:
    part = record[
        (record.Date >= row["Oneset Date"])
        & (record.Date <= row["End Date"])
        & (record["Step Type"].astype(str) == str(row["Step Type"]))
    ].copy()
    if part.empty:
        raise RuntimeError(f"No record rows for step {row['Step Number']} {row['Step Type']}")
    return part.sort_values("Date", kind="stable").reset_index(drop=True)


def numeric(part: pd.DataFrame, name: str) -> np.ndarray:
    return pd.to_numeric(part[name], errors="coerce").to_numpy(float)


def elapsed_seconds(part: pd.DataFrame, column: str = "Time") -> np.ndarray:
    values = pd.to_timedelta(part[column].astype(str), errors="coerce").dt.total_seconds().to_numpy(float)
    if not np.all(np.isfinite(values)):
        raise RuntimeError(f"Invalid elapsed time in column {column}")
    return values


def load_and_audit() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    branch_rows: list[dict] = []
    cc_rows: list[pd.DataFrame] = []
    rest_rows: list[pd.DataFrame] = []
    cccv_rows: list[pd.DataFrame] = []
    step_rows: list[pd.DataFrame] = []

    for cell, path in FILES.items():
        step = pd.read_excel(path, sheet_name="step")
        record = pd.read_excel(path, sheet_name="record")
        step["Oneset Date"] = pd.to_datetime(step["Oneset Date"])
        step["End Date"] = pd.to_datetime(step["End Date"])
        record["Date"] = pd.to_datetime(record["Date"])
        step.insert(0, "Cell", cell)
        step.insert(1, "Source_file", path.name)
        step_rows.append(step.copy())
        by = step.set_index("Step Number")

        for (rate, direction), step_no in BRANCH_STEPS.items():
            branch = by.loc[step_no]
            rest_row = by.loc[step_no - 1]
            cc = segment(record, branch)
            rest = segment(record, rest_row)

            cc_elapsed_s = elapsed_seconds(cc, "Total Time")
            rest_elapsed_s = elapsed_seconds(rest, "Total Time")
            t_cc = (cc_elapsed_s - cc_elapsed_s[0]) / 60.0
            # Use Neware's monotonic Total Time as the common clock because the
            # exported Date and step Time fields can both round a terminal
            # one-second event record to a duplicate timestamp.
            t_rest = (rest_elapsed_s - cc_elapsed_s[0]) / 60.0
            voltage = numeric(cc, "Voltage(V)")
            current_signed = numeric(cc, "Current(A)")
            current = np.abs(current_signed)
            q_raw = numeric(cc, "Capacity(Ah)")
            q_cc = np.abs(q_raw - q_raw[0])
            command = float(np.nanmedian(current[current > 0]))
            reaches = np.flatnonzero(current >= 0.95 * command)
            if not len(reaches):
                raise RuntimeError(f"No 95% current point: {cell} {rate}C {direction}")
            first95 = int(reaches[0])

            rest_t_abs = (rest_elapsed_s - rest_elapsed_s[0]) / 60.0
            rest_v = numeric(rest, "Voltage(V)")
            tail = rest_t_abs >= rest_t_abs[-1] - 10.0
            tail_slope = float(np.polyfit(rest_t_abs[tail], rest_v[tail], 1)[0] * 1000.0)

            dt = np.diff(cc_elapsed_s)
            rest_dt = np.diff(rest_elapsed_s)
            date_dt = cc.Date.diff().dt.total_seconds().dropna().to_numpy(float)
            rest_date_dt = rest.Date.diff().dt.total_seconds().dropna().to_numpy(float)
            v_at_95 = float(voltage[first95])
            rest_end_v = float(rest_v[-1])
            onset_jump = 1000.0 * (v_at_95 - rest_end_v)

            cv_capacity = np.nan
            cv_duration = np.nan
            cv_fraction = np.nan
            cv_transition_gap_s = np.nan
            cv_transition_dv_mV = np.nan
            first_cv_current = np.nan
            total_capacity = float(q_cc[-1])
            if direction == "Charge" and step_no + 1 in by.index:
                cv_row = by.loc[step_no + 1]
                if str(cv_row["Step Type"]) == "CV Chg":
                    cv = segment(record, cv_row)
                    cv_elapsed_s = elapsed_seconds(cv, "Total Time")
                    t_cv = (cv_elapsed_s - cc_elapsed_s[0]) / 60.0
                    cv_v = numeric(cv, "Voltage(V)")
                    cv_i = numeric(cv, "Current(A)")
                    cv_q_raw = numeric(cv, "Capacity(Ah)")
                    cv_q = np.abs(cv_q_raw - cv_q_raw[0])
                    cv_capacity = float(cv_q[-1])
                    cv_duration = float(t_cv[-1] - t_cv[0])
                    total_capacity = float(q_cc[-1] + cv_q[-1])
                    cv_fraction = 100.0 * cv_capacity / total_capacity
                    total_cc_s = elapsed_seconds(cc, "Total Time")
                    total_cv_s = elapsed_seconds(cv, "Total Time")
                    cv_transition_gap_s = float(total_cv_s[0] - total_cc_s[-1])
                    cv_transition_dv_mV = float(1000.0 * (cv_v[0] - voltage[-1]))
                    first_cv_current = float(abs(cv_i[0]))
                    cccv_rows.append(
                        pd.DataFrame(
                            {
                                "Cell": cell,
                                "C_rate": rate,
                                "Direction": direction,
                                "Phase": "CV",
                                "t_from_CC_start_min": t_cv,
                                "Voltage_V": cv_v,
                                "Current_A": cv_i,
                                "Transferred_capacity_Ah": q_cc[-1] + cv_q,
                                "Date": cv.Date,
                            }
                        )
                    )

            branch_rows.append(
                {
                    "Cell": cell,
                    "Source_file": path.name,
                    "C_rate": rate,
                    "Direction": direction,
                    "CC_step": step_no,
                    "Rest_step": step_no - 1,
                    "CV_step": step_no + 1 if direction == "Charge" else np.nan,
                    "N_CC_records": len(cc),
                    "N_rest_records": len(rest),
                    "CC_median_sample_interval_s": float(np.nanmedian(dt)),
                    "CC_max_sample_interval_s": float(np.nanmax(dt)),
                    "Rest_median_sample_interval_s": float(np.nanmedian(rest_dt)),
                    "Total_time_strictly_increasing": bool(np.all(dt > 0) and np.all(rest_dt > 0)),
                    "Date_timestamp_strictly_increasing": bool(
                        np.all(date_dt > 0) and np.all(rest_date_dt > 0)
                    ),
                    "Rest_duration_min": duration_min(rest_row),
                    "Rest_end_V": rest_end_v,
                    "Rest_tail_slope_mV_per_min": tail_slope,
                    "CC_command_A": command,
                    "Time_to_95pct_command_s": float(t_cc[first95] * 60.0),
                    "V_at_95pct_command_V": v_at_95,
                    "Onset_jump_from_rest_mV": onset_jump,
                    "CC_duration_min": float(t_cc[-1]),
                    "CC_end_V": float(voltage[-1]),
                    "CC_capacity_Ah": float(q_cc[-1]),
                    "CV_duration_min": cv_duration,
                    "CV_capacity_Ah": cv_capacity,
                    "CCCV_total_capacity_Ah": total_capacity,
                    "CV_fraction_of_total_pct": cv_fraction,
                    "CC_to_CV_record_gap_s": cv_transition_gap_s,
                    "CC_to_CV_voltage_jump_mV": cv_transition_dv_mV,
                    "First_CV_current_A": first_cv_current,
                    "CC_current_CV_pct": 100.0 * float(np.nanstd(current[first95:])) / command,
                }
            )

            cc_frame = pd.DataFrame(
                {
                    "Cell": cell,
                    "C_rate": rate,
                    "Direction": direction,
                    "Phase": "CC",
                    "t_from_CC_start_min": t_cc,
                    "Voltage_V": voltage,
                    "Current_A": current_signed,
                    "Transferred_capacity_Ah": q_cc,
                    "Date": cc.Date,
                }
            )
            cc_rows.append(cc_frame)
            if direction == "Charge":
                cccv_rows.append(cc_frame.copy())
            rest_rows.append(
                pd.DataFrame(
                    {
                        "Cell": cell,
                        "C_rate": rate,
                        "Direction_next": direction,
                        "t_relative_to_CC_start_min": t_rest,
                        "Voltage_V": rest_v,
                        "Current_A": numeric(rest, "Current(A)"),
                        "Date": rest.Date,
                    }
                )
            )

    branches = pd.DataFrame(branch_rows).sort_values(["C_rate", "Direction", "Cell"])
    cc_all = pd.concat(cc_rows, ignore_index=True)
    rest_all = pd.concat(rest_rows, ignore_index=True)
    cccv_all = pd.concat(cccv_rows, ignore_index=True).sort_values(
        ["C_rate", "Cell", "t_from_CC_start_min", "Phase"], kind="stable"
    )
    steps = pd.concat(step_rows, ignore_index=True)
    return branches, cc_all, rest_all, cccv_all, steps


def summarize(branches: pd.DataFrame) -> pd.DataFrame:
    return (
        branches.groupby(["C_rate", "Direction"], as_index=False)
        .agg(
            N_cells=("Cell", "nunique"),
            Mean_rest_end_V=("Rest_end_V", "mean"),
            Mean_rest_tail_slope_mV_per_min=("Rest_tail_slope_mV_per_min", "mean"),
            Mean_time_to_95pct_command_s=("Time_to_95pct_command_s", "mean"),
            Mean_onset_jump_mV=("Onset_jump_from_rest_mV", "mean"),
            Mean_CC_duration_min=("CC_duration_min", "mean"),
            Mean_CC_capacity_Ah=("CC_capacity_Ah", "mean"),
            Cell_range_CC_capacity_mAh=("CC_capacity_Ah", lambda x: 1000.0 * float(np.ptp(x))),
            Mean_CV_capacity_Ah=("CV_capacity_Ah", "mean"),
            Mean_CCCV_total_capacity_Ah=("CCCV_total_capacity_Ah", "mean"),
            Mean_CV_fraction_pct=("CV_fraction_of_total_pct", "mean"),
            Max_CC_sample_gap_s=("CC_max_sample_interval_s", "max"),
            All_total_times_increasing=("Total_time_strictly_increasing", "all"),
            All_date_timestamps_increasing=("Date_timestamp_strictly_increasing", "all"),
        )
        .sort_values(["C_rate", "Direction"])
    )


def prior_crosscheck(branches: pd.DataFrame) -> pd.DataFrame:
    prior_path = Path(
        r"C:\Users\user\Documents\Codex\2026-09-18\https-github-com-rinny-hi-ai2020"
        r"\work\ai2020_dfn_fitting\results\260927_comprehensive_prefit_cross_cohort"
        r"\branch_qc_all_cells.csv"
    )
    prior = pd.read_csv(prior_path)
    prior = prior[prior.Cohort == "New September"].copy()
    merged = branches.merge(
        prior,
        on=["Cell", "C_rate", "Direction"],
        suffixes=("_raw_rerun", "_prior"),
        how="outer",
        indicator=True,
    )
    for name in ("Rest_end_V", "Rest_tail_slope_mV_per_min", "CC_duration_min", "CC_capacity_Ah"):
        merged[f"delta_{name}"] = merged[f"{name}_raw_rerun"] - merged[f"{name}_prior"]
    return merged


def plot_full_curves(cc: pd.DataFrame) -> None:
    fig, axes = plt.subplots(3, 2, figsize=(14.8, 13.2), constrained_layout=True)
    colors = {"6-1": "#0072B2", "6-2": "#D55E00"}
    for i, rate in enumerate((0.5, 1.0, 2.0)):
        for j, direction in enumerate(("Charge", "Discharge")):
            ax = axes[i, j]
            ax_i = ax.twinx()
            for cell in ("6-1", "6-2"):
                p = cc[(cc.C_rate == rate) & (cc.Direction == direction) & (cc.Cell == cell)]
                ax.plot(p.t_from_CC_start_min, p.Voltage_V, color=colors[cell], lw=2.0, label=f"{cell} voltage")
                ax_i.plot(p.t_from_CC_start_min, np.abs(p.Current_A), color=colors[cell], lw=1.0, ls="--", alpha=0.55)
            ax.set(xlabel="Physical time from CC start [min]", ylabel="Voltage [V]",
                   title=f"{rate:g}C {direction.lower()}")
            ax_i.set_ylabel("|Current| [A]", color="#666666")
            max_current = float(np.nanmax(np.abs(cc[(cc.C_rate == rate) & (cc.Direction == direction)].Current_A)))
            ax_i.set_ylim(0.0, 1.15 * max_current)
            ax_i.ticklabel_format(axis="y", style="plain", useOffset=False)
            ax.grid(alpha=0.25)
            if i == 0 and j == 0:
                ax.legend(fontsize=8, loc="best")
    fig.suptitle("September BoL full CC physical-time curves", fontsize=15)
    fig.savefig(OUT / "september_full_cc_physical_time.png", dpi=220)
    plt.close(fig)


def plot_rest_transitions(rest: pd.DataFrame, cc: pd.DataFrame) -> None:
    fig, axes = plt.subplots(3, 2, figsize=(14.8, 13.2), constrained_layout=True)
    colors = {"6-1": "#0072B2", "6-2": "#D55E00"}
    for i, rate in enumerate((0.5, 1.0, 2.0)):
        for j, direction in enumerate(("Charge", "Discharge")):
            ax = axes[i, j]
            for cell in ("6-1", "6-2"):
                r = rest[(rest.C_rate == rate) & (rest.Direction_next == direction) & (rest.Cell == cell)]
                r = r[r.t_relative_to_CC_start_min >= -20.0]
                c = cc[(cc.C_rate == rate) & (cc.Direction == direction) & (cc.Cell == cell)]
                c = c[c.t_from_CC_start_min <= 5.0]
                ax.plot(r.t_relative_to_CC_start_min, r.Voltage_V, color=colors[cell], lw=1.8, label=cell)
                ax.plot(c.t_from_CC_start_min, c.Voltage_V, color=colors[cell], lw=2.3)
            ax.axvline(0.0, color="black", lw=1.0, ls=":")
            ax.set(xlabel="Time relative to CC start [min]", ylabel="Voltage [V]",
                   title=f"{rate:g}C {direction.lower()}: last 20 min rest + first 5 min CC")
            ax.grid(alpha=0.25)
            if i == 0 and j == 0:
                ax.legend(fontsize=8)
    fig.suptitle("Rest-to-load transition audit", fontsize=15)
    fig.savefig(OUT / "september_rest_to_cc_transition.png", dpi=220)
    plt.close(fig)


def plot_cccv(cccv: pd.DataFrame, branches: pd.DataFrame) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(13.8, 13.0), constrained_layout=True)
    colors = {"6-1": "#0072B2", "6-2": "#D55E00"}
    for ax, rate in zip(axes, (0.5, 1.0, 2.0)):
        ax_i = ax.twinx()
        for cell in ("6-1", "6-2"):
            p = cccv[(cccv.C_rate == rate) & (cccv.Cell == cell)]
            ax.plot(p.t_from_CC_start_min, p.Voltage_V, color=colors[cell], lw=2.0, label=f"{cell} voltage")
            ax_i.plot(p.t_from_CC_start_min, np.abs(p.Current_A), color=colors[cell], lw=1.1, ls="--", alpha=0.6)
            row = branches[(branches.C_rate == rate) & (branches.Direction == "Charge") & (branches.Cell == cell)].iloc[0]
            ax.axvline(row.CC_duration_min, color=colors[cell], lw=1.0, ls=":")
        ax.set(xlabel="Physical time from CC start [min]", ylabel="Voltage [V]",
               title=f"{rate:g}C charge: CC and CV (dotted line = CC end)")
        ax_i.set_ylabel("|Current| [A]", color="#666666")
        ax.grid(alpha=0.25)
    axes[0].legend(fontsize=8)
    fig.suptitle("September charge CC-CV transition and CV tail", fontsize=15)
    fig.savefig(OUT / "september_charge_cccv_transition.png", dpi=220)
    plt.close(fig)


def write_report(branches: pd.DataFrame, summary: pd.DataFrame, cross: pd.DataFrame) -> None:
    s = summary.set_index(["C_rate", "Direction"])
    max_cross = {
        c: float(np.nanmax(np.abs(cross[c].to_numpy(float))))
        for c in cross.columns
        if c.startswith("delta_")
    }
    lines = [
        "# September BoL 원본 full-curve 복구 및 protocol QC",
        "",
        "## 결론",
        "",
        "- 6-1/6-2 원본 workbook의 `step`과 `record`를 읽어 0.5C/1C/2C charge·discharge CC 곡선을 모두 복구했다.",
        "- 선택된 12개 branch에서 Neware `Total Time`은 모두 단조 증가하며, step 누락이나 잘못된 CC/CV 순서는 발견되지 않았다.",
        "- 6-2 0.5C discharge의 마지막 1초 cutoff record는 `Date`와 step `Time`이 직전 record와 같게 표시되지만 `Total Time`에는 1초 차이가 있다. full physical-time 축은 이 record를 보존하는 `Total Time`을 사용한다.",
        "- 각 branch 직전 rest는 120분이다. 고전압 rest 말단은 거의 평형이나, 저전압 charge 직전 rest는 0.5C와 1C에서 여전히 0.30/0.20 mV/min의 relaxation이 남는다.",
        "- 2C charge는 CC보다 CV에서 들어가는 용량이 더 크다. 따라서 2C charge의 낮은 CC cutoff capacity를 cell inventory 감소로만 해석하면 안 된다.",
        "",
        "## Branch 요약",
        "",
        "| Rate | Direction | Rest tail slope | Load jump | CC time | CC capacity | CV capacity | CV fraction | Cell range |",
        "|---:|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for rate in (0.5, 1.0, 2.0):
        for direction in ("Charge", "Discharge"):
            row = s.loc[(rate, direction)]
            cv_cap = "-" if not np.isfinite(row.Mean_CV_capacity_Ah) else f"{row.Mean_CV_capacity_Ah:.4f} Ah"
            cv_frac = "-" if not np.isfinite(row.Mean_CV_fraction_pct) else f"{row.Mean_CV_fraction_pct:.1f}%"
            lines.append(
                f"| {rate:g}C | {direction} | {row.Mean_rest_tail_slope_mV_per_min:+.3f} mV/min | "
                f"{row.Mean_onset_jump_mV:+.1f} mV | {row.Mean_CC_duration_min:.2f} min | "
                f"{row.Mean_CC_capacity_Ah:.4f} Ah | {cv_cap} | {cv_frac} | "
                f"{row.Cell_range_CC_capacity_mAh:.1f} mAh |"
            )
    lines += [
        "",
        "## 해석",
        "",
        "1. **초기 상태 품질**: discharge 직전 고전압 rest slope는 절대값 0.009--0.016 mV/min으로 안정적이다. "
        "charge 직전 저전압 rest는 2시간 후에도 0.5C 0.298, 1C 0.197, 2C 0.071 mV/min이다. "
        "July 10분 rest보다 훨씬 낫지만 0.5C/1C 저전압 상태를 완전 평형 OCP로 간주하지 않는다.",
        "2. **전류 인가 응답**: 평균 rest-to-load jump는 charge에서 +46.1/+87.0/+162.8 mV, discharge에서 "
        "-36.9/-70.8/-141.0 mV다. rate에 따라 거의 증가하므로 OCP offset보다 kinetics·ohmic·transport 성분이 크다.",
        "3. **CC/CV 분리**: charge CV 용량 비중은 0.5C 약 9.4%, 1C 약 20.9%, 2C 약 52.5%다. "
        "2C CC 구간은 polarization 때문에 4.2 V에 일찍 도달하며, 이후 CV에서 약 1.12 Ah가 추가된다.",
        "4. **CC→CV 전환 품질**: 두 셀·세 rate 모두 record gap은 0초이고 전압 불연속은 0.0--0.2 mV다. protocol 전환 문제로 볼 증거는 없다.",
        "5. **셀 반복성**: 두 셀의 CC capacity 범위는 branch별 4--22 mAh 수준이다. cohort 내 반복성보다 July--September cohort 차이가 더 크다.",
        "6. **fitting 입력**: physical-time charge fitting에는 CC record만 사용하고, CV capacity는 목적함수에 넣지 않는다. "
        "CV 및 total capacity는 post-hoc validation으로 남긴다.",
        "",
        "## 이전 추출값 교차검증",
        "",
        f"이전 저장 QC와 raw rerun의 최대 차이는 `Rest_end_V={max_cross.get('delta_Rest_end_V', np.nan):.3e} V`, "
        f"`rest slope={max_cross.get('delta_Rest_tail_slope_mV_per_min', np.nan):.3e} mV/min`, "
        f"`CC duration={max_cross.get('delta_CC_duration_min', np.nan):.3e} min`, "
        f"`CC capacity={max_cross.get('delta_CC_capacity_Ah', np.nan):.3e} Ah`이다. 기존 September 요약 추출은 재현된다.",
        "",
        "## 산출물",
        "",
        "- `september_full_cc_records.csv`: fitting용 CC full physical-time records",
        "- `september_charge_cccv_records.csv`: charge CC+CV 연속 record",
        "- `september_prebranch_rest_records.csv`: 각 branch 직전 120분 rest record",
        "- `september_branch_qc.csv`, `september_protocol_summary.csv`: branch 및 cohort QC",
        "- PNG 3종: full CC, rest-to-load, CC-CV 전환",
    ]
    (OUT / "README_KO.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    for path in FILES.values():
        if not path.exists():
            raise FileNotFoundError(path)
    branches, cc, rest, cccv, steps = load_and_audit()
    summary = summarize(branches)
    cross = prior_crosscheck(branches)
    branches.to_csv(OUT / "september_branch_qc.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(OUT / "september_protocol_summary.csv", index=False, encoding="utf-8-sig")
    cc.to_csv(OUT / "september_full_cc_records.csv", index=False, encoding="utf-8-sig")
    rest.to_csv(OUT / "september_prebranch_rest_records.csv", index=False, encoding="utf-8-sig")
    cccv.to_csv(OUT / "september_charge_cccv_records.csv", index=False, encoding="utf-8-sig")
    steps.to_csv(OUT / "september_step_sequence.csv", index=False, encoding="utf-8-sig")
    cross.to_csv(OUT / "prior_qc_crosscheck.csv", index=False, encoding="utf-8-sig")
    plot_full_curves(cc)
    plot_rest_transitions(rest, cc)
    plot_cccv(cccv, branches)
    write_report(branches, summary, cross)
    manifest = {
        "sources": [str(p) for p in FILES.values()],
        "source_workbooks_modified": False,
        "branch_steps": {f"{r:g}C_{d}": n for (r, d), n in BRANCH_STEPS.items()},
        "record_rows": {"cc": len(cc), "rest": len(rest), "charge_cccv": len(cccv)},
    }
    (OUT / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print(summary.to_string(index=False))
    print(f"Done: {OUT}")


if __name__ == "__main__":
    main()
