"""September BoL conservative Dsn+kn fit on fixed physical-time grids.

The two source workbooks are never modified.  This script consumes the
read-only full-curve CSV extraction produced by september_raw_fullcurve_audit.py.
Capacity and the CV phase are excluded from the optimizer residual and are
reported only after the voltage fit.
"""

from __future__ import annotations

import json
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution, least_squares
from scipy.stats import qmc


PROJECT = Path(__file__).resolve().parents[1]
REPO = Path(
    r"C:\Users\user\Documents\Codex\2026-09-18"
    r"\https-github-com-rinny-hi-ai2020\work\ai2020_dfn_fitting"
)
SOURCE = PROJECT / "outputs" / "260927_september_raw_fullcurve_audit"
OUT = PROJECT / "outputs" / "260927_september_recommended_dsn_kn_fit"
sys.path.insert(0, str(REPO))

import c50_selected_ocp_dynamic_validation as selected
import charge_physical_time_common as july_common
import comprehensive_prefit_cross_cohort as cross
import electrode_hysteresis_quantification as ehq
import paper_geometry_area_porosity_ablation as area_study
import p0_five_parameter_sensitivity_area as legacy
import prefit_physical_anchor_decision as p0
import priority_initial_state_protocol_recheck as priority


RATES = (0.5, 1.0, 2.0)
CELLS = ("6-1", "6-2")
N_POINTS = 100
PARAMETERS = ("Dsn", "kn")
BOUNDS = {
    "Dsn": (1.0e-14, 8.0e-14),
    "kn": (3.0e-7, 3.0e-6),
}
GITT_START = {"Dsn": 4.170774145968247e-14, "kn": 9.648533212e-7}
JULY_FINAL = {"Dsn": 7.466200473885764e-14, "kn": 1.4502357402082455e-6}
FIXED = {
    "Dsp": 6.888175323807748e-14,
    "brugg_n": 2.914,
    "brugg_p": 1.83,
    "brugg_s": 1.5,
    "hysteresis": False,
}


@dataclass
class FitContext:
    curves: dict
    qc: pd.DataFrame
    model_geometry: dict
    stage: dict
    anode: dict
    cathode: dict
    qcell: float
    kp: float
    grids: dict
    observed: dict
    initial_states: pd.DataFrame


def normalized_to_values(z) -> dict[str, float]:
    z = np.asarray(z, float)
    out = {}
    for coordinate, name in zip(z, PARAMETERS):
        lo, hi = BOUNDS[name]
        out[name] = float(np.exp(np.log(lo) + coordinate * (np.log(hi) - np.log(lo))))
    return out


def values_to_normalized(values: dict[str, float]) -> np.ndarray:
    return np.asarray(
        [
            (np.log(float(values[name])) - np.log(BOUNDS[name][0]))
            / (np.log(BOUNDS[name][1]) - np.log(BOUNDS[name][0]))
            for name in PARAMETERS
        ],
        float,
    )


def prepare_branch(part: pd.DataFrame) -> dict[str, np.ndarray]:
    part = part.sort_values("t_from_CC_start_min", kind="stable")
    t = part.t_from_CC_start_min.to_numpy(float)
    v = part.Voltage_V.to_numpy(float)
    i_signed = part.Current_A.to_numpy(float)
    q = part.Transferred_capacity_Ah.to_numpy(float)
    current = np.abs(i_signed)
    command = float(np.nanmedian(current[current > 0]))
    matches = np.flatnonzero(current >= 0.95 * command)
    if not len(matches):
        raise RuntimeError("No point reaches 95% of the CC current")
    first = int(matches[0])
    return {
        "t_min": t[first:] - t[first],
        "V": v[first:],
        "I_A": i_signed[first:],
        "Q_Ah": q[first:] - q[first],
        "command_A": command,
    }


def build_context() -> FitContext:
    records = pd.read_csv(SOURCE / "september_full_cc_records.csv")
    qc = pd.read_csv(SOURCE / "september_branch_qc.csv")
    curves = {}
    for (cell, rate, direction), part in records.groupby(
        ["Cell", "C_rate", "Direction"], sort=False
    ):
        curves[(str(cell), float(rate), str(direction))] = prepare_branch(part)

    model0, stage, anode, cathode, _, _, qcell = selected.selected_inputs()
    anchored = p0.build_model(model0, p0.Scenario("P0 anchored reference"))
    geometry = area_study.apply_paper_overlap(anchored, preserve_capacity=True)
    kp = float(legacy.kp_from_area_specific_rct(geometry)["kp"])

    state_rows = []
    grids = {}
    observed = {}
    for _, row in qc.iterrows():
        cell = str(row.Cell)
        rate = float(row.C_rate)
        direction = str(row.Direction)
        z0 = priority.monotone_voltage_inverse(
            geometry["v_qocv"], geometry["soc"], float(row.Rest_end_V)
        )
        x0 = float(stage["x0"] + z0 * (stage["x100"] - stage["x0"]))
        y0 = float(stage["y0"] - z0 * (stage["y0"] - stage["y100"]))
        state_rows.append(
            {
                "Cell": cell,
                "C_rate": rate,
                "Direction": direction,
                "Rest_end_V": float(row.Rest_end_V),
                "Rest_tail_slope_mV_per_min": float(row.Rest_tail_slope_mV_per_min),
                "Initial_SOC_pct": 100.0 * float(z0),
                "Initial_negative_stoichiometry_x": x0,
                "Initial_positive_stoichiometry_y": y0,
                "Measured_current_A": float(row.CC_command_A),
            }
        )
        curve = curves[(cell, rate, direction)]
        grid = np.linspace(0.0, float(curve["t_min"][-1]), N_POINTS)
        grids[(cell, rate, direction)] = grid
        observed[(cell, rate, direction)] = np.interp(
            grid, curve["t_min"], curve["V"]
        )
    return FitContext(
        curves=curves,
        qc=qc,
        model_geometry=geometry,
        stage=stage,
        anode=anode,
        cathode=cathode,
        qcell=float(qcell),
        kp=kp,
        grids=grids,
        observed=observed,
        initial_states=pd.DataFrame(state_rows).sort_values(
            ["Direction", "C_rate", "Cell"]
        ),
    )


def model_from_values(context: FitContext, values: dict[str, float]) -> dict:
    return legacy.apply_parameters(
        context.model_geometry,
        {
            "Dsn": float(values["Dsn"]),
            "Dsp": FIXED["Dsp"],
            "kn": float(values["kn"]),
            "kp": context.kp,
            "brugg_n": FIXED["brugg_n"],
        },
    )


def simulate_direction(
    context: FitContext, values: dict[str, float], direction: str
) -> dict:
    charge = direction == "Charge"
    model = model_from_values(context, values)
    runs = {}
    states = context.initial_states[context.initial_states.Direction == direction]
    for _, row in states.iterrows():
        local_stage = priority.stage_at_soc(
            context.stage, float(row.Initial_SOC_pct) / 100.0, charge
        )
        rate_for_model = float(row.Measured_current_A) / cross.NOMINAL_CAPACITY_AH
        key = (str(row.Cell), float(row.C_rate), direction)
        runs[key] = ehq.run_dfn(
            rate_for_model,
            charge,
            local_stage,
            model,
            context.anode,
            context.cathode,
            False,
            False,
        )
    return runs


def physical_time_residual(
    context: FitContext,
    runs: dict,
    direction: str,
    penalty_mV: float | None = None,
) -> tuple[np.ndarray | None, dict]:
    blocks = []
    margins = {}
    missing = False
    for rate in RATES:
        for cell in CELLS:
            key = (cell, rate, direction)
            grid = context.grids[key]
            sim_t = np.asarray(runs[key]["t_min"], float)
            sim_v = np.asarray(runs[key]["V"], float)
            margin = float(sim_t[-1] - grid[-1])
            margins[key] = margin
            if margin < -1.0e-9:
                missing = True
            prediction = np.interp(np.minimum(grid, sim_t[-1]), sim_t, sim_v)
            residual = 1000.0 * (prediction - context.observed[key])
            if margin < -1.0e-9 and penalty_mV is not None:
                overrun = np.maximum(grid - sim_t[-1], 0.0) / max(grid[-1], 1e-12)
                residual = residual + penalty_mV * (1.0 + overrun)
            blocks.append(residual)
    vector = np.concatenate(blocks)
    if missing and penalty_mV is None:
        vector = None
    return vector, {
        "feasible": not missing,
        "minimum_time_margin_min": float(min(margins.values())),
        "time_margins_min": {
            f"{cell}_{rate:g}C_{direction}": value
            for (cell, rate, _), value in margins.items()
        },
    }


class Evaluator:
    def __init__(self, context: FitContext):
        self.context = context
        self.cache = {}
        self.history = []

    def evaluate(self, z):
        z = np.clip(np.asarray(z, float), 0.0, 1.0)
        cache_key = tuple(np.round(z, 12))
        if cache_key in self.cache:
            return self.cache[cache_key]
        values = normalized_to_values(z)
        started = time.perf_counter()
        status = "ok"
        runs = None
        try:
            runs = simulate_direction(self.context, values, "Charge")
            exact, feasibility = physical_time_residual(
                self.context, runs, "Charge", penalty_mV=None
            )
            if exact is None:
                residual, _ = physical_time_residual(
                    self.context, runs, "Charge", penalty_mV=1000.0
                )
                status = "invalid_early_cutoff"
            else:
                residual = exact
        except Exception as exc:
            residual = np.full(len(CELLS) * len(RATES) * N_POINTS, 2000.0)
            feasibility = {
                "feasible": False,
                "minimum_time_margin_min": np.nan,
                "time_margins_min": {},
            }
            status = f"failed:{type(exc).__name__}"
        rmse = float(np.sqrt(np.mean(residual**2)))
        self.history.append(
            {
                "Evaluation": len(self.history) + 1,
                "Objective_RMSE_mV": rmse,
                "Status": status,
                "Feasible": bool(feasibility["feasible"]),
                "Minimum_time_margin_min": feasibility["minimum_time_margin_min"],
                "Elapsed_s": time.perf_counter() - started,
                **values,
            }
        )
        if len(self.history) == 1 or len(self.history) % 10 == 0:
            print(
                f"evaluation {len(self.history)}: {rmse:.3f} mV ({status})",
                flush=True,
            )
            self.save_history()
        result = residual, runs, feasibility, values
        self.cache[cache_key] = result
        return result

    def residual(self, z):
        return self.evaluate(z)[0]

    def scalar(self, z):
        r = self.residual(z)
        return float(np.mean(r**2))

    def save_history(self):
        pd.DataFrame(self.history).to_csv(
            OUT / "optimizer_evaluation_history.csv",
            index=False,
            encoding="utf-8-sig",
        )


def finite_difference_identifiability(
    evaluator: Evaluator, center: np.ndarray, label: str, step: float = 0.025
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    center = np.asarray(center, float)
    center_residual, _, center_feas, _ = evaluator.evaluate(center)
    if not center_feas["feasible"]:
        return (
            pd.DataFrame(),
            pd.DataFrame(),
            {
                "Label": label,
                "Pass": False,
                "Reason": "anchor is not full-grid feasible",
            },
        )
    columns = []
    rows = []
    for index, name in enumerate(PARAMETERS):
        lo = center.copy()
        hi = center.copy()
        lo[index] = max(0.0, center[index] - step)
        hi[index] = min(1.0, center[index] + step)
        lo_r, _, lo_f, _ = evaluator.evaluate(lo)
        hi_r, _, hi_f, _ = evaluator.evaluate(hi)
        if lo_f["feasible"] and hi_f["feasible"]:
            derivative = (hi_r - lo_r) / (hi[index] - lo[index])
            scheme = "central"
        elif hi_f["feasible"]:
            derivative = (hi_r - center_residual) / (hi[index] - center[index])
            scheme = "forward"
        elif lo_f["feasible"]:
            derivative = (center_residual - lo_r) / (center[index] - lo[index])
            scheme = "backward"
        else:
            return (
                pd.DataFrame(),
                pd.DataFrame(),
                {"Label": label, "Pass": False, "Reason": f"no feasible side for {name}"},
            )
        columns.append(derivative)
        rows.append(
            {
                "Label": label,
                "Parameter": name,
                "Difference_scheme": scheme,
                "RMS_mV_per_normalized_coordinate": float(
                    np.sqrt(np.mean(derivative**2))
                ),
            }
        )
    jacobian = np.column_stack(columns)
    sens = pd.DataFrame(rows)
    sens["Relative_sensitivity"] = (
        sens.RMS_mV_per_normalized_coordinate
        / sens.RMS_mV_per_normalized_coordinate.max()
    )
    unit = jacobian / np.maximum(
        np.linalg.norm(jacobian, axis=0), np.finfo(float).eps
    )
    cosine = unit.T @ unit
    singular = np.linalg.svd(unit, compute_uv=False)
    condition = float(singular[0] / singular[-1])
    max_pair = float(abs(cosine[0, 1]))
    rank = int(np.linalg.matrix_rank(jacobian))
    passed = bool(
        rank == len(PARAMETERS)
        and float(sens.Relative_sensitivity.min()) >= 0.10
        and max_pair <= 0.95
        and condition <= 20.0
    )
    corr = pd.DataFrame(cosine, index=PARAMETERS, columns=PARAMETERS)
    result = {
        "Label": label,
        "Pass": passed,
        "Matrix_rank": rank,
        "Min_relative_sensitivity": float(sens.Relative_sensitivity.min()),
        "Max_abs_pairwise_cosine": max_pair,
        "Normalized_condition_number": condition,
        "Normalized_step": step,
    }
    return sens, corr, result


def branch_metrics(context: FitContext, runs: dict, direction: str) -> pd.DataFrame:
    rows = []
    for rate in RATES:
        for cell in CELLS:
            key = (cell, rate, direction)
            exp = context.curves[key]
            sim = runs[key]
            grid = context.grids[key]
            sim_t = np.asarray(sim["t_min"], float)
            sim_v = np.asarray(sim["V"], float)
            feasible = bool(sim_t[-1] >= grid[-1] - 1e-9)
            pred = np.interp(np.minimum(grid, sim_t[-1]), sim_t, sim_v)
            err = 1000.0 * (pred - context.observed[key])
            q_exp = float(np.asarray(exp["Q_Ah"], float)[-1])
            q_model = float(np.asarray(sim["Q_Ah"], float)[-1])
            rows.append(
                {
                    "Cell": cell,
                    "C_rate": rate,
                    "Direction": direction,
                    "Physical_time_feasible": feasible,
                    "Time_margin_min": float(sim_t[-1] - grid[-1]),
                    "Physical_time_RMSE_mV": float(np.sqrt(np.mean(err**2)))
                    if feasible
                    else np.nan,
                    "Physical_time_MAE_mV": float(np.mean(np.abs(err)))
                    if feasible
                    else np.nan,
                    "Physical_time_bias_mV": float(np.mean(err)) if feasible else np.nan,
                    "Experimental_CC_capacity_Ah": q_exp,
                    "Model_cutoff_capacity_Ah": q_model,
                    "Capacity_error_mAh": 1000.0 * (q_model - q_exp),
                    "Capacity_error_pct": 100.0 * (q_model - q_exp) / q_exp,
                }
            )
    return pd.DataFrame(rows)


def plot_fit(context: FitContext, charge_runs: dict, discharge_runs: dict) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(16.0, 9.3), constrained_layout=True)
    colors = {"6-1": "#0072B2", "6-2": "#D55E00"}
    for row_index, (direction, runs) in enumerate(
        (("Charge", charge_runs), ("Discharge", discharge_runs))
    ):
        for col_index, rate in enumerate(RATES):
            ax = axes[row_index, col_index]
            for cell in CELLS:
                key = (cell, rate, direction)
                exp = context.curves[key]
                sim = runs[key]
                ax.plot(
                    exp["t_min"],
                    exp["V"],
                    color=colors[cell],
                    lw=2.0,
                    label=f"{cell} experiment",
                )
                ax.plot(
                    sim["t_min"],
                    sim["V"],
                    color=colors[cell],
                    lw=1.6,
                    ls="--",
                    label=f"{cell} DFN",
                )
            ax.set(
                title=f"{rate:g}C {direction}",
                xlabel="Physical time from 95% current [min]",
                ylabel="Voltage [V]",
            )
            ax.grid(alpha=0.25)
            if row_index == 0 and col_index == 0:
                ax.legend(fontsize=8)
    fig.suptitle(
        "September conservative fit: charge training and untouched discharge validation",
        fontsize=15,
    )
    fig.savefig(OUT / "charge_fit_discharge_validation.png", dpi=220)
    plt.close(fig)


def plot_state_audit(context: FitContext) -> None:
    data = context.initial_states.copy()
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.9), constrained_layout=True)
    colors = {"Charge": "#0072B2", "Discharge": "#D55E00"}
    offsets = {"6-1": -0.04, "6-2": 0.04}
    for direction in ("Charge", "Discharge"):
        for cell in CELLS:
            p = data[(data.Direction == direction) & (data.Cell == cell)]
            axes[0].plot(
                p.C_rate + offsets[cell],
                p.Initial_SOC_pct,
                "o-",
                color=colors[direction],
                alpha=0.65 if cell == "6-2" else 1.0,
                label=f"{direction}, {cell}",
            )
    axes[0].set(
        xlabel="C-rate",
        ylabel="Rest-voltage-inferred initial SOC [%]",
        title="Per-cell, per-branch initialization",
    )
    axes[0].grid(alpha=0.25)
    axes[0].legend(fontsize=8)

    stage = context.stage
    labels = ["x0", "x100", "y100", "y0"]
    values = [stage[name] for name in labels]
    axes[1].bar(labels, values, color=["#4472C4", "#4472C4", "#ED7D31", "#ED7D31"])
    axes[1].set(
        ylabel="Electrode stoichiometry",
        title="Retained identifiable OCP window",
        ylim=(0, 1.05),
    )
    axes[1].grid(axis="y", alpha=0.25)
    axes[1].text(
        0.02,
        0.98,
        "September has dynamic branches and rest anchors,\nnot an independent full qOCV curve.\nWindow refit would be underdetermined.",
        transform=axes[1].transAxes,
        va="top",
        fontsize=9,
        bbox=dict(facecolor="white", alpha=0.9, edgecolor="none"),
    )
    fig.savefig(OUT / "initial_state_and_window_audit.png", dpi=220)
    plt.close(fig)


def write_report(
    context: FitContext,
    candidates: pd.DataFrame,
    best: dict,
    identifiability: pd.DataFrame,
    metrics: pd.DataFrame,
) -> None:
    summary = (
        metrics.groupby("Direction", as_index=False)
        .agg(
            Physical_time_RMSE_mV=(
                "Physical_time_RMSE_mV",
                lambda x: float(np.sqrt(np.nanmean(np.asarray(x, float) ** 2))),
            ),
            Capacity_RMSE_pct=(
                "Capacity_error_pct",
                lambda x: float(np.sqrt(np.mean(np.asarray(x, float) ** 2))),
            ),
            Capacity_MAE_mAh=(
                "Capacity_error_mAh", lambda x: float(np.mean(np.abs(x)))
            ),
        )
    )
    charge = summary[summary.Direction == "Charge"].iloc[0]
    discharge = summary[summary.Direction == "Discharge"].iloc[0]
    state = context.initial_states
    low = state[state.Direction == "Charge"]
    high = state[state.Direction == "Discharge"]
    ident = identifiability.set_index("Label")
    stage = context.stage
    lines = [
        "# September BoL 권장 설정 기반 Dsn+kn fitting",
        "",
        "## 결론",
        "",
        f"- 9월 6-1/6-2의 각 rate·cell 직전 2시간 rest 종점전압으로 초기 SOC를 개별 설정했다. "
        f"charge 시작 SOC는 {low.Initial_SOC_pct.min():.2f}--{low.Initial_SOC_pct.max():.2f}%, "
        f"discharge 시작 SOC는 {high.Initial_SOC_pct.min():.2f}--{high.Initial_SOC_pct.max():.2f}%다.",
        "- OCP는 `Ai2020 nominal graphite equilibrium + old measured cathode equilibrium`, hysteresis는 OFF, "
        "`brugg_n/p/s=2.914/1.83/1.5`를 유지했다.",
        "- 9월 데이터에는 독립적인 full qOCV 곡선이 없어서 4개 OCP-window endpoint와 QLi를 새로 동시에 식별할 수 없다. "
        "따라서 근거 없는 window refit은 하지 않고, 기존 식별 window를 유지한 채 branch별 초기 상태만 9월 rest 전압으로 갱신했다.",
        f"- 선택된 보수 모델은 `Dsn={best['Dsn']:.6e} m2/s`, `kn={best['kn']:.6e}`이며 "
        f"charge physical-time objective는 {best['Objective_RMSE_mV']:.2f} mV다.",
        f"- untouched discharge physical-time RMSE는 {discharge.Physical_time_RMSE_mV:.2f} mV다. "
        f"사후 capacity RMSE는 charge {charge.Capacity_RMSE_pct:.2f}%, discharge {discharge.Capacity_RMSE_pct:.2f}%다.",
        "",
        "## 초기상태와 OCP window 결정",
        "",
        f"- 유지한 window: `x0={stage['x0']:.7f}`, `x100={stage['x100']:.7f}`, "
        f"`y100={stage['y100']:.7f}`, `y0={stage['y0']:.7f}`, `Qcell={context.qcell:.7f} Ah`.",
        "- 고SOC rest에서 셀·rate 간 초기 SOC가 매우 좁게 모이는지, 저SOC rest에서 discharge rate에 따른 cutoff-state 차이가 "
        "합리적으로 남는지를 확인했다. 이 구조는 모든 branch를 0%/100%로 강제하는 것보다 측정 protocol을 잘 반영한다.",
        "- 저SOC rest에는 2시간 후에도 양의 relaxation slope가 남으므로, rest 종점 기반 SOC는 완전 평형값이 아니라 "
        "현재 데이터에서 가장 직접적인 초기상태 anchor다.",
        "",
        "## 민감도·식별도",
        "",
        "| 위치 | 통과 | 최소 상대민감도 | max |cosine| | condition |",
        "|---|---:|---:|---:|---:|",
    ]
    for label in ("pre-fit July-selected anchor", "post-fit optimum"):
        row = ident.loc[label]
        lines.append(
            f"| {label} | {bool(row.Pass)} | {row.Min_relative_sensitivity:.3f} | "
            f"{row.Max_abs_pairwise_cosine:.3f} | {row.Normalized_condition_number:.2f} |"
        )
    lines += [
        "",
        "통과 기준은 full rank, 최소 상대민감도 0.10 이상, max |cosine| 0.95 이하, normalized condition 20 이하이다.",
        "",
        "## 최적화 후보",
        "",
        "| 후보 | charge RMSE | Dsn | kn | full-time feasible |",
        "|---|---:|---:|---:|---:|",
    ]
    for _, row in candidates.sort_values("Objective_RMSE_mV").iterrows():
        lines.append(
            f"| {row.Candidate} | {row.Objective_RMSE_mV:.3f} mV | {row.Dsn:.4e} | "
            f"{row.kn:.4e} | {bool(row.Feasible)} |"
        )
    lines += [
        "",
        "## 목적함수와 사후지표 분리",
        "",
        "- 목적함수: charge CC 0.5C/1C/2C, 2셀, 각 100개 physical-time voltage point의 ordinary RMSE.",
        "- 제외: capacity, CV 구간, slope/shape penalty, prior/regularization.",
        "- 사후평가: charge/discharge CC cutoff capacity error와 untouched discharge time-voltage RMSE.",
        "",
        "## 제한",
        "",
        "- 이 fitting은 2개의 새 셀과 6개 charge branch에 기반한다. Dsn은 상수 effective parameter이며 intrinsic graphite 물성 확정값이 아니다.",
        "- 9월 전용 OCP window를 독립적으로 재식별하려면 같은 셀의 저율 full-cell qOCV 또는 양·음극 half-cell OCP가 필요하다.",
        "- 2C charge의 CV 용량 비중이 크므로 CC cutoff capacity만으로 lithium inventory를 역산하지 않았다.",
    ]
    (OUT / "README_KO.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    context = build_context()
    context.initial_states.to_csv(
        OUT / "september_initial_states.csv", index=False, encoding="utf-8-sig"
    )
    plot_state_audit(context)

    evaluator = Evaluator(context)
    pre_sens, pre_corr, pre_ident = finite_difference_identifiability(
        evaluator,
        values_to_normalized(JULY_FINAL),
        "pre-fit July-selected anchor",
    )
    if not pre_ident.get("Pass", False):
        raise RuntimeError(f"Pre-fit Dsn+kn subset failed identifiability: {pre_ident}")

    starts = [
        ("TRF GITT nominal", values_to_normalized(GITT_START)),
        ("TRF July final", values_to_normalized(JULY_FINAL)),
        ("TRF bound center", np.full(2, 0.5)),
    ]
    lhs = qmc.LatinHypercube(d=2, seed=260927).random(1)[0]
    starts.append(("TRF log-space LHS", lhs))
    candidate_rows = []
    candidate_vectors = {}
    for label, start in starts:
        print(label, flush=True)
        begun = time.perf_counter()
        result = least_squares(
            evaluator.residual,
            start,
            bounds=(np.zeros(2), np.ones(2)),
            method="trf",
            x_scale="jac",
            diff_step=1e-3,
            max_nfev=60,
            ftol=1e-4,
            xtol=1e-4,
            gtol=1e-4,
        )
        residual, _, feas, values = evaluator.evaluate(result.x)
        candidate_vectors[label] = np.asarray(result.x, float)
        candidate_rows.append(
            {
                "Candidate": label,
                "Algorithm": "TRF",
                "Success": bool(result.success),
                "Message": str(result.message),
                "NFEV_reported": int(result.nfev),
                "Wall_s": time.perf_counter() - begun,
                "Objective_RMSE_mV": float(np.sqrt(np.mean(residual**2))),
                "Feasible": bool(feas["feasible"]),
                "Minimum_time_margin_min": feas["minimum_time_margin_min"],
                **values,
            }
        )

    print("Differential evolution -> TRF", flush=True)
    begun = time.perf_counter()
    de = differential_evolution(
        evaluator.scalar,
        [(0.0, 1.0), (0.0, 1.0)],
        popsize=5,
        maxiter=6,
        seed=260927,
        polish=False,
        workers=1,
        updating="immediate",
        tol=0.01,
    )
    polished = least_squares(
        evaluator.residual,
        de.x,
        bounds=(np.zeros(2), np.ones(2)),
        method="trf",
        x_scale="jac",
        diff_step=1e-3,
        max_nfev=60,
        ftol=1e-4,
        xtol=1e-4,
        gtol=1e-4,
    )
    residual, _, feas, values = evaluator.evaluate(polished.x)
    label = "DE -> TRF"
    candidate_vectors[label] = np.asarray(polished.x, float)
    candidate_rows.append(
        {
            "Candidate": label,
            "Algorithm": "DE + TRF",
            "Success": bool(polished.success),
            "Message": str(polished.message),
            "NFEV_reported": int(polished.nfev),
            "Wall_s": time.perf_counter() - begun,
            "Objective_RMSE_mV": float(np.sqrt(np.mean(residual**2))),
            "Feasible": bool(feas["feasible"]),
            "Minimum_time_margin_min": feas["minimum_time_margin_min"],
            **values,
        }
    )

    candidates = pd.DataFrame(candidate_rows)
    candidates.to_csv(
        OUT / "optimizer_candidates.csv", index=False, encoding="utf-8-sig"
    )
    evaluator.save_history()
    feasible = candidates[candidates.Feasible.astype(bool)].sort_values(
        "Objective_RMSE_mV"
    )
    if feasible.empty:
        raise RuntimeError("No full-time feasible optimizer candidate")
    best_row = feasible.iloc[0]
    best_label = str(best_row.Candidate)
    best_vector = candidate_vectors[best_label]
    best_values = normalized_to_values(best_vector)
    best = {
        "Candidate": best_label,
        "Objective_RMSE_mV": float(best_row.Objective_RMSE_mV),
        **best_values,
    }

    post_sens, post_corr, post_ident = finite_difference_identifiability(
        evaluator, best_vector, "post-fit optimum"
    )
    sensitivity = pd.concat([pre_sens, post_sens], ignore_index=True)
    sensitivity.to_csv(
        OUT / "sensitivity_ranking.csv", index=False, encoding="utf-8-sig"
    )
    pre_corr.to_csv(
        OUT / "sensitivity_cosine_prefit.csv", encoding="utf-8-sig"
    )
    post_corr.to_csv(
        OUT / "sensitivity_cosine_postfit.csv", encoding="utf-8-sig"
    )
    identifiability = pd.DataFrame([pre_ident, post_ident])
    identifiability.to_csv(
        OUT / "identifiability_summary.csv", index=False, encoding="utf-8-sig"
    )

    charge_runs = simulate_direction(context, best_values, "Charge")
    discharge_runs = simulate_direction(context, best_values, "Discharge")
    metrics = pd.concat(
        [
            branch_metrics(context, charge_runs, "Charge"),
            branch_metrics(context, discharge_runs, "Discharge"),
        ],
        ignore_index=True,
    )
    metrics.to_csv(OUT / "branch_metrics.csv", index=False, encoding="utf-8-sig")
    summary = (
        metrics.groupby(["Direction", "C_rate"], as_index=False)
        .agg(
            N_cells=("Cell", "count"),
            All_physical_time_feasible=("Physical_time_feasible", "all"),
            Mean_physical_time_RMSE_mV=("Physical_time_RMSE_mV", "mean"),
            Mean_physical_time_MAE_mV=("Physical_time_MAE_mV", "mean"),
            Mean_physical_time_bias_mV=("Physical_time_bias_mV", "mean"),
            Capacity_RMSE_pct=(
                "Capacity_error_pct",
                lambda x: float(np.sqrt(np.mean(np.asarray(x, float) ** 2))),
            ),
            Capacity_MAE_mAh=(
                "Capacity_error_mAh", lambda x: float(np.mean(np.abs(x)))
            ),
        )
    )
    summary.to_csv(OUT / "summary_by_rate.csv", index=False, encoding="utf-8-sig")
    plot_fit(context, charge_runs, discharge_runs)

    manifest = {
        "selected_ocp": "Ai2020 nominal graphite equilibrium + old measured cathode equilibrium",
        "window_policy": "retain previously identified window; September branch rest voltages define per-cell/per-rate initial SOC because an independent September full qOCV curve is unavailable",
        "stage": {key: float(context.stage[key]) for key in ("x0", "x100", "y100", "y0")},
        "qcell_Ah": context.qcell,
        "objective": "ordinary charge CC voltage RMSE on 100-point fixed physical-time grids, equal branch weight",
        "capacity_in_objective": False,
        "cv_in_objective": False,
        "fixed": {**FIXED, "kp": context.kp},
        "bounds": BOUNDS,
        "best": best,
        "prefit_identifiability": pre_ident,
        "postfit_identifiability": post_ident,
    }
    (OUT / "analysis_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_report(context, candidates, best, identifiability, metrics)
    print(json.dumps(best, ensure_ascii=False, indent=2), flush=True)
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
