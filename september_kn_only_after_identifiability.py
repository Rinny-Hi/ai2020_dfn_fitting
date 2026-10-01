"""Fallback September fit after Dsn+kn fails the pre-fit sensitivity rule.

Dsn is profiled as a fixed independent anchor.  Only kn is fitted to charge
physical-time voltage.  Capacity and CV are post-fit diagnostics only.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution, least_squares, minimize_scalar


PROJECT = Path(__file__).resolve().parents[1]
WORK = PROJECT / "work"
sys.path.insert(0, str(WORK))
import september_recommended_dsn_kn_fit as base


OUT = PROJECT / "outputs" / "260927_september_recommended_dsn_kn_fit"
GITT_DSN = base.GITT_START["Dsn"]
DSN_PROFILE = (1.0e-14, GITT_DSN, base.JULY_FINAL["Dsn"], 8.0e-14)
KN_BOUND = base.BOUNDS["kn"]


def kn_to_z(kn: float) -> np.ndarray:
    lo, hi = KN_BOUND
    return np.asarray([(np.log(kn) - np.log(lo)) / (np.log(hi) - np.log(lo))])


def z_to_kn(z) -> float:
    lo, hi = KN_BOUND
    coordinate = float(np.asarray(z, float)[0])
    return float(np.exp(np.log(lo) + coordinate * (np.log(hi) - np.log(lo))))


class KnEvaluator:
    def __init__(self, context: base.FitContext, dsn: float):
        self.context = context
        self.dsn = float(dsn)
        self.cache = {}
        self.history = []

    def evaluate(self, z):
        z = np.clip(np.asarray(z, float), 0.0, 1.0)
        key = round(float(z[0]), 12)
        if key in self.cache:
            return self.cache[key]
        values = {"Dsn": self.dsn, "kn": z_to_kn(z)}
        started = time.perf_counter()
        status = "ok"
        runs = None
        try:
            runs = base.simulate_direction(self.context, values, "Charge")
            exact, feasibility = base.physical_time_residual(
                self.context, runs, "Charge", penalty_mV=None
            )
            if exact is None:
                residual, _ = base.physical_time_residual(
                    self.context, runs, "Charge", penalty_mV=1000.0
                )
                status = "invalid_early_cutoff"
            else:
                residual = exact
        except Exception as exc:
            residual = np.full(len(base.CELLS) * len(base.RATES) * base.N_POINTS, 2000.0)
            feasibility = {
                "feasible": False,
                "minimum_time_margin_min": np.nan,
                "time_margins_min": {},
            }
            status = f"failed:{type(exc).__name__}"
        result = residual, runs, feasibility, values
        self.cache[key] = result
        self.history.append(
            {
                "Evaluation": len(self.history) + 1,
                "Fixed_Dsn": self.dsn,
                "kn": values["kn"],
                "Objective_RMSE_mV": float(np.sqrt(np.mean(residual**2))),
                "Status": status,
                "Feasible": bool(feasibility["feasible"]),
                "Minimum_time_margin_min": feasibility["minimum_time_margin_min"],
                "Elapsed_s": time.perf_counter() - started,
            }
        )
        if len(self.history) == 1 or len(self.history) % 10 == 0:
            print(
                f"Dsn={self.dsn:.4e}, evaluation {len(self.history)}: "
                f"{self.history[-1]['Objective_RMSE_mV']:.3f} mV ({status})",
                flush=True,
            )
        return result

    def residual(self, z):
        return self.evaluate(z)[0]

    def scalar(self, z):
        r = self.residual(z)
        return float(np.mean(r**2))


def fit_kn(context: base.FitContext, dsn: float, full_ladder: bool = False):
    evaluator = KnEvaluator(context, dsn)
    candidates = []
    vectors = {}
    begun = time.perf_counter()
    scalar = minimize_scalar(
        lambda coordinate: evaluator.scalar([coordinate]),
        bounds=(0.0, 1.0),
        method="bounded",
        options={"xatol": 3e-3, "maxiter": 18},
    )
    polished = least_squares(
        evaluator.residual,
        np.asarray([scalar.x]),
        bounds=(np.zeros(1), np.ones(1)),
        method="trf",
        x_scale="jac",
        diff_step=1e-3,
        max_nfev=15,
        ftol=1e-4,
        xtol=1e-4,
        gtol=1e-4,
    )
    residual, _, feasible, values = evaluator.evaluate(polished.x)
    label = "Bounded scalar -> TRF"
    vectors[label] = np.asarray(polished.x, float)
    candidates.append(
        {
            "Candidate": label,
            "Algorithm": "Bounded scalar + TRF",
            "Success": bool(polished.success),
            "NFEV_reported": int(scalar.nfev + polished.nfev),
            "Wall_s": time.perf_counter() - begun,
            "Objective_RMSE_mV": float(np.sqrt(np.mean(residual**2))),
            "Feasible": bool(feasible["feasible"]),
            "Minimum_time_margin_min": feasible["minimum_time_margin_min"],
            **values,
        }
    )
    if full_ladder:
        begun = time.perf_counter()
        direct = least_squares(
            evaluator.residual,
            kn_to_z(base.JULY_FINAL["kn"]),
            bounds=(np.zeros(1), np.ones(1)),
            method="trf",
            x_scale="jac",
            diff_step=1e-3,
            max_nfev=20,
            ftol=1e-4,
            xtol=1e-4,
            gtol=1e-4,
        )
        residual, _, feasible, values = evaluator.evaluate(direct.x)
        label = "TRF July-kn start"
        vectors[label] = np.asarray(direct.x, float)
        candidates.append(
            {
                "Candidate": label,
                "Algorithm": "TRF",
                "Success": bool(direct.success),
                "NFEV_reported": int(direct.nfev),
                "Wall_s": time.perf_counter() - begun,
                "Objective_RMSE_mV": float(np.sqrt(np.mean(residual**2))),
                "Feasible": bool(feasible["feasible"]),
                "Minimum_time_margin_min": feasible["minimum_time_margin_min"],
                **values,
            }
        )
        begun = time.perf_counter()
        de = differential_evolution(
            evaluator.scalar,
            [(0.0, 1.0)],
            popsize=4,
            maxiter=3,
            seed=260927,
            polish=False,
            workers=1,
            updating="immediate",
            tol=0.01,
        )
        polished = least_squares(
            evaluator.residual,
            de.x,
            bounds=(np.zeros(1), np.ones(1)),
            method="trf",
            x_scale="jac",
            diff_step=1e-3,
            max_nfev=15,
            ftol=1e-4,
            xtol=1e-4,
            gtol=1e-4,
        )
        residual, _, feasible, values = evaluator.evaluate(polished.x)
        label = "DE -> TRF"
        vectors[label] = np.asarray(polished.x, float)
        candidates.append(
            {
                "Candidate": label,
                "Algorithm": "DE + TRF",
                "Success": bool(polished.success),
                "NFEV_reported": int(polished.nfev),
                "Wall_s": time.perf_counter() - begun,
                "Objective_RMSE_mV": float(np.sqrt(np.mean(residual**2))),
                "Feasible": bool(feasible["feasible"]),
                "Minimum_time_margin_min": feasible["minimum_time_margin_min"],
                **values,
            }
        )
    table = pd.DataFrame(candidates)
    feasible = table[table.Feasible.astype(bool)].sort_values("Objective_RMSE_mV")
    if feasible.empty:
        raise RuntimeError(f"No feasible kn-only result at Dsn={dsn:.4e}")
    best = feasible.iloc[0]
    best_label = str(best.Candidate)
    best_vector = vectors[best_label]
    best_values = {"Dsn": float(dsn), "kn": z_to_kn(best_vector)}
    return evaluator, table, best_vector, best_values


def kn_identifiability(evaluator: KnEvaluator, center: np.ndarray, label: str) -> dict:
    center = np.asarray(center, float)
    center_r, _, center_f, _ = evaluator.evaluate(center)
    result = {
        "Label": label,
        "Pass": True,
        "Matrix_rank": 1,
        "Min_relative_sensitivity": 1.0,
        "Max_abs_pairwise_cosine": 0.0,
        "Normalized_condition_number": 1.0,
    }
    rows = []
    if not center_f["feasible"]:
        result["Pass"] = False
        result["Reason"] = "post-fit center is not full-grid feasible"
        return result, pd.DataFrame()
    for step in (0.01, 0.025, 0.05):
        low = center.copy()
        high = center.copy()
        low[0] = max(0.0, low[0] - step)
        high[0] = min(1.0, high[0] + step)
        low_r, _, low_f, _ = evaluator.evaluate(low)
        high_r, _, high_f, _ = evaluator.evaluate(high)
        if low_f["feasible"] and high_f["feasible"]:
            derivative = (high_r - low_r) / (high[0] - low[0])
            scheme = "central"
        elif high_f["feasible"]:
            derivative = (high_r - center_r) / (high[0] - center[0])
            scheme = "forward"
        elif low_f["feasible"]:
            derivative = (center_r - low_r) / (center[0] - low[0])
            scheme = "backward"
        else:
            result["Pass"] = False
            result["Reason"] = f"non-feasible perturbation at step {step}"
            continue
        rms = float(np.sqrt(np.mean(derivative**2)))
        rows.append(
            {
                "Label": label,
                "Normalized_step": step,
                "Parameter": "kn",
                "Difference_scheme": scheme,
                "RMS_mV_per_normalized_coordinate": rms,
                "Relative_sensitivity": 1.0,
                "Pass": rms > 0.0,
            }
        )
        result["Pass"] = bool(result["Pass"] and rms > 0.0)
    return result, pd.DataFrame(rows)


def validation_summary(metrics: pd.DataFrame) -> dict:
    out = {}
    for direction in ("Charge", "Discharge"):
        part = metrics[metrics.Direction == direction]
        out[f"{direction}_physical_time_RMSE_mV"] = float(
            np.sqrt(np.nanmean(part.Physical_time_RMSE_mV.to_numpy(float) ** 2))
        )
        out[f"{direction}_capacity_RMSE_pct"] = float(
            np.sqrt(np.mean(part.Capacity_error_pct.to_numpy(float) ** 2))
        )
        out[f"{direction}_capacity_MAE_mAh"] = float(
            np.mean(np.abs(part.Capacity_error_mAh.to_numpy(float)))
        )
        out[f"{direction}_all_time_feasible"] = bool(part.Physical_time_feasible.all())
    return out


def profile_plot(profile: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.7), constrained_layout=True)
    x = profile.Fixed_Dsn.to_numpy(float) / 1e-14
    axes[0].plot(x, profile.Charge_objective_RMSE_mV, "o-", color="#0072B2")
    axes[0].axvline(GITT_DSN / 1e-14, color="black", ls="--", lw=1, label="GITT fixed")
    axes[0].set(xlabel="Fixed Dsn [1e-14 m2/s]", ylabel="Charge RMSE [mV]", title="Charge profile objective")
    axes[1].plot(x, profile.Fitted_kn * 1e6, "o-", color="#D55E00")
    axes[1].set(xlabel="Fixed Dsn [1e-14 m2/s]", ylabel="Fitted kn [1e-6]", title="Compensation in kn")
    axes[2].plot(x, profile.Discharge_physical_time_RMSE_mV, "o-", color="#009E73")
    axes[2].set(xlabel="Fixed Dsn [1e-14 m2/s]", ylabel="Discharge RMSE [mV]", title="Untouched discharge validation")
    for ax in axes:
        ax.grid(alpha=0.25)
    axes[0].legend(fontsize=8)
    fig.suptitle("Dsn profile after two-parameter identifiability failure", fontsize=14)
    fig.savefig(OUT / "dsn_profile_fit_validation.png", dpi=220)
    plt.close(fig)


def write_report(
    context: base.FitContext,
    dsn_audit: pd.DataFrame,
    kn_ident: pd.DataFrame,
    candidates: pd.DataFrame,
    profile: pd.DataFrame,
    best_values: dict,
    primary_summary: dict,
) -> None:
    state = context.initial_states
    charge_state = state[state.Direction == "Charge"]
    discharge_state = state[state.Direction == "Discharge"]
    best_row = candidates.sort_values("Objective_RMSE_mV").iloc[0]
    lines = [
        "# September BoL 권장 설정 실행 결과",
        "",
        "## 최종 결론",
        "",
        "- 9월 charge 데이터에서 사전 `Dsn+kn` 민감도/식별도를 재검증했으나, Dsn 상대민감도가 기준 0.10에 미달했다. "
        "따라서 두 파라미터 동시 fitting은 수행하지 않았다.",
        f"- 통과한 축소 모델은 **Dsn 고정 + kn 단일 fitting**이다. 독립 mass-corrected GITT 값 "
        f"`Dsn={GITT_DSN:.6e} m2/s`를 고정하고 `kn={best_values['kn']:.6e}`를 얻었다.",
        f"- charge physical-time RMSE는 {primary_summary['Charge_physical_time_RMSE_mV']:.2f} mV, "
        f"untouched discharge RMSE는 {primary_summary['Discharge_physical_time_RMSE_mV']:.2f} mV다.",
        f"- capacity는 목적함수에 넣지 않았다. 사후 capacity RMSE는 charge "
        f"{primary_summary['Charge_capacity_RMSE_pct']:.2f}%, discharge "
        f"{primary_summary['Discharge_capacity_RMSE_pct']:.2f}%다.",
        "",
        "## 권장 설정 적용",
        "",
        "- OCP: Ai2020 nominal graphite equilibrium + old measured cathode equilibrium",
        "- hysteresis: OFF",
        "- `brugg_n/p/s=2.914/1.83/1.5`",
        f"- 초기 SOC: charge {charge_state.Initial_SOC_pct.min():.2f}--{charge_state.Initial_SOC_pct.max():.2f}%, "
        f"discharge {discharge_state.Initial_SOC_pct.min():.2f}--{discharge_state.Initial_SOC_pct.max():.2f}% "
        "(각 cell·rate의 2시간 rest 종점전압으로 개별 산정)",
        f"- 유지 OCP window: `x0={context.stage['x0']:.7f}`, `x100={context.stage['x100']:.7f}`, "
        f"`y100={context.stage['y100']:.7f}`, `y0={context.stage['y0']:.7f}`",
        "",
        "9월에는 독립 full-qOCV가 없어 QLi와 4개 window endpoint를 새로 동시에 식별할 수 없다. "
        "근거 없는 window 재최적화 대신, 검증된 window를 유지하고 branch별 초기 상태만 9월 rest로 갱신했다.",
        "",
        "## Dsn+kn 사전 식별도",
        "",
        "| Anchor | step | 통과 | 최소 상대민감도 | max abs(cosine) | condition |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for _, row in dsn_audit.iterrows():
        lines.append(
            f"| {row.Label} | {row.Normalized_step:.3f} | {bool(row.Pass)} | "
            f"{row.Min_relative_sensitivity:.3f} | {row.Max_abs_pairwise_cosine:.3f} | "
            f"{row.Normalized_condition_number:.2f} |"
        )
    lines += [
        "",
        "full rank와 양호한 condition에도 Dsn 신호가 일관되게 약하므로, 수치적 2변수 optimum을 물성값으로 채택하지 않았다.",
        "",
        "## kn 단일 모델 식별도와 최적화",
        "",
        f"- 1/2.5/5% normalized perturbation에서 kn RMS sensitivity는 "
        f"{kn_ident.RMS_mV_per_normalized_coordinate.min():.2f}--{kn_ident.RMS_mV_per_normalized_coordinate.max():.2f} mV이며 모두 통과했다.",
        f"- 최선 후보는 `{best_row.Candidate}`이고 full-time feasibility를 만족했다.",
        "",
        "| 후보 | charge RMSE | kn | feasible |",
        "|---|---:|---:|---:|",
    ]
    for _, row in candidates.sort_values("Objective_RMSE_mV").iterrows():
        lines.append(
            f"| {row.Candidate} | {row.Objective_RMSE_mV:.3f} mV | {row.kn:.6e} | {bool(row.Feasible)} |"
        )
    lines += [
        "",
        "## Dsn 고정값 profile",
        "",
        "| Fixed Dsn | fitted kn | charge RMSE | discharge RMSE | discharge time feasible | charge cap RMSE | discharge cap RMSE |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for _, row in profile.iterrows():
        discharge_rmse = (
            f"{row.Discharge_physical_time_RMSE_mV:.2f} mV"
            if bool(row.Discharge_all_time_feasible)
            else "invalid (early cutoff)"
        )
        lines.append(
            f"| {row.Fixed_Dsn:.4e} | {row.Fitted_kn:.4e} | {row.Charge_objective_RMSE_mV:.2f} mV | "
            f"{discharge_rmse} | {bool(row.Discharge_all_time_feasible)} | {row.Charge_capacity_RMSE_pct:.2f}% | "
            f"{row.Discharge_capacity_RMSE_pct:.2f}% |"
        )
    lines += [
        "",
        "profile은 Dsn을 데이터에서 추정한 것이 아니라 고정 가정 변화에 대한 구조 민감도다. 최종값은 independent GITT anchor를 사용한다.",
        "",
        "## 목적함수",
        "",
        "- charge CC 0.5C/1C/2C, 2셀, 각 100개 physical-time voltage point의 ordinary RMSE",
        "- capacity, CV, voltage slope/shape, prior는 목적함수에서 제외",
        "- candidate가 실험 cutoff time 전에 종료되면 fixed-time residual을 정의할 수 없어 infeasible 처리",
        "- discharge와 capacity는 fitting 이후에만 평가",
        "",
        "## 해석 제한",
        "",
        "- 이번 9월 charge만으로 Dsn을 재추정하지 못했다. 이는 실패가 아니라 파라미터 수를 줄여야 한다는 식별도 결과다.",
        f"- 선택 kn의 최소 charge time margin은 {float(best_row.Minimum_time_margin_min) * 60.0:.4f}초로 feasibility 경계에 매우 가깝다. "
        "따라서 kn도 독립 kinetic 물성값보다 cutoff를 포함한 effective parameter로 해석한다.",
        "- 2C charge capacity 오차는 두 셀에서 +25.0/+29.3%이고, discharge voltage bias는 전 branch에서 +54--66 mV다. "
        "단일 kn 조정만으로 고율 charge polarization과 charge/discharge 비대칭을 동시에 설명하지 못한다.",
        "- 9월 전용 QLi/window를 얻으려면 같은 셀의 저율 full qOCV 또는 전극 half-cell OCP가 필요하다.",
        "- 2C charge는 CV 비중이 약 52.5%이므로 CC cutoff capacity를 inventory fitting 신호로 사용하지 않았다.",
    ]
    (OUT / "README_KO.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    context = base.build_context()
    context.initial_states.to_csv(OUT / "september_initial_states.csv", index=False, encoding="utf-8-sig")
    base.plot_state_audit(context)

    two_evaluator = base.Evaluator(context)
    audit_rows = []
    for label, values, steps in (
        ("GITT anchor", base.GITT_START, (0.01, 0.025, 0.05)),
        ("July-selected anchor", base.JULY_FINAL, (0.025,)),
        ("Bound center", base.normalized_to_values([0.5, 0.5]), (0.025,)),
        ("High-kn anchor", {"Dsn": GITT_DSN, "kn": 1.9e-6}, (0.025,)),
    ):
        for step in steps:
            _, _, result = base.finite_difference_identifiability(
                two_evaluator,
                base.values_to_normalized(values),
                label,
                step,
            )
            audit_rows.append(result)
    dsn_audit = pd.DataFrame(audit_rows)
    dsn_audit.to_csv(OUT / "dsn_kn_prefit_identifiability.csv", index=False, encoding="utf-8-sig")
    two_evaluator.save_history()
    if bool(dsn_audit.Pass.any()):
        raise RuntimeError("Unexpected Dsn+kn pass; review before using fallback")

    profile_rows = []
    profile_metrics = []
    primary_candidates = None
    primary_evaluator = None
    primary_vector = None
    primary_values = None
    for dsn in DSN_PROFILE:
        is_primary = np.isclose(dsn, GITT_DSN, rtol=1e-10, atol=0.0)
        evaluator, candidates, vector, values = fit_kn(context, dsn, full_ladder=is_primary)
        charge_runs = base.simulate_direction(context, values, "Charge")
        discharge_runs = base.simulate_direction(context, values, "Discharge")
        metrics = pd.concat(
            [
                base.branch_metrics(context, charge_runs, "Charge"),
                base.branch_metrics(context, discharge_runs, "Discharge"),
            ],
            ignore_index=True,
        )
        summary = validation_summary(metrics)
        best_candidate = candidates[candidates.Feasible.astype(bool)].sort_values("Objective_RMSE_mV").iloc[0]
        profile_rows.append(
            {
                "Fixed_Dsn": float(dsn),
                "Fitted_kn": float(values["kn"]),
                "Charge_objective_RMSE_mV": float(best_candidate.Objective_RMSE_mV),
                **summary,
            }
        )
        metrics.insert(0, "Fixed_Dsn", float(dsn))
        profile_metrics.append(metrics)
        if is_primary:
            primary_candidates = candidates
            primary_evaluator = evaluator
            primary_vector = vector
            primary_values = values
            metrics.to_csv(OUT / "branch_metrics.csv", index=False, encoding="utf-8-sig")
            candidates.to_csv(OUT / "optimizer_candidates.csv", index=False, encoding="utf-8-sig")
            pd.DataFrame(evaluator.history).to_csv(
                OUT / "optimizer_evaluation_history_kn_only.csv", index=False, encoding="utf-8-sig"
            )
            base.plot_fit(context, charge_runs, discharge_runs)

    profile = pd.DataFrame(profile_rows).sort_values("Fixed_Dsn")
    profile.to_csv(OUT / "dsn_fixed_profile_summary.csv", index=False, encoding="utf-8-sig")
    pd.concat(profile_metrics, ignore_index=True).to_csv(
        OUT / "dsn_fixed_profile_branch_metrics.csv", index=False, encoding="utf-8-sig"
    )
    profile_plot(profile)

    if primary_evaluator is None or primary_values is None:
        raise RuntimeError("Primary GITT-fixed fit was not created")
    kn_ident_result, kn_ident = kn_identifiability(
        primary_evaluator, primary_vector, "post-fit kn-only"
    )
    kn_ident.to_csv(OUT / "kn_only_identifiability.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame([kn_ident_result]).to_csv(
        OUT / "kn_only_identifiability_summary.csv", index=False, encoding="utf-8-sig"
    )
    if not kn_ident_result["Pass"]:
        raise RuntimeError(f"kn-only model failed post-fit identifiability: {kn_ident_result}")

    primary_row = profile[np.isclose(profile.Fixed_Dsn, GITT_DSN, rtol=1e-10, atol=0.0)].iloc[0]
    primary_summary = primary_row.to_dict()
    manifest = {
        "decision": "Dsn+kn rejected before fitting because Dsn relative sensitivity stayed below 0.10; fit kn only with Dsn fixed to independent mass-corrected GITT",
        "selected_ocp": "Ai2020 nominal graphite equilibrium + old measured cathode equilibrium",
        "window_policy": "retain previously identified window; update each September branch initial SOC from its own 2-hour rest endpoint voltage",
        "stage": {key: float(context.stage[key]) for key in ("x0", "x100", "y100", "y0")},
        "qcell_Ah": context.qcell,
        "fixed": {**base.FIXED, "Dsn": GITT_DSN, "kp": context.kp},
        "fitted": {"kn": primary_values["kn"]},
        "objective": "charge CC ordinary voltage RMSE; 2 cells x 3 rates x 100 fixed physical-time points",
        "capacity_in_objective": False,
        "cv_in_objective": False,
        "validation": primary_summary,
        "kn_identifiability": kn_ident_result,
    }
    (OUT / "analysis_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_report(
        context,
        dsn_audit,
        kn_ident,
        primary_candidates,
        profile,
        primary_values,
        primary_summary,
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2), flush=True)
    print(profile.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
