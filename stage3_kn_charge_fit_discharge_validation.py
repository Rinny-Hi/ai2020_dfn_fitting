"""Stage-3 first fit: fit kn on charge, validate without refitting on discharge.

The primary objective is the standard equal-curve voltage least-squares loss.
Charge termination capacity is handled independently by (i) reporting it,
(ii) Pareto analysis, and (iii) a sweep of explicit capacity-RMSE constraints.
No Bayesian prior or ad-hoc voltage/capacity weighted sum is used.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import nominal_radius_current_configuration as current
import stage2_seven_parameter_sensitivity as stage2


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260921_stage3_kn_charge_fit"
RESULTS.mkdir(parents=True, exist_ok=True)

KN_MIN = 1.0e-7
KN_MAX = 3.0e-6
N_SCAN = 25
N_VOLTAGE_POINTS = 120
CAPACITY_CONSTRAINTS_PCT = (0.5, 1.0, 2.0, 3.0, 5.0)


def model_with_kn(baseline, kn_value):
    changed = dict(baseline)
    changed["params"] = baseline["params"].copy()
    nominal = baseline["params"][stage2.PARAMETER_KEYS["kn"]]
    scaled = ga.scale_parameter_function(nominal, float(kn_value) / current.KN_PREF)
    changed["params"].update(
        {stage2.PARAMETER_KEYS["kn"]: scaled}, check_already_exists=False
    )
    return changed


def q_axis(data, q_meas, charge):
    soc = np.asarray(data["SOC"], dtype=float)
    return soc * q_meas if charge else (1.0 - soc) * q_meas


def end_capacity(data, q_meas, charge):
    return float(q_axis(data, q_meas, charge)[-1])


def fixed_voltage_grids(experiment, baseline_runs, q_meas, charge):
    grids = {}
    for rate in ga.RATES:
        exp_q = q_axis(experiment[(rate, charge)], q_meas, charge)
        sim_q = q_axis(baseline_runs[(rate, charge)], q_meas, charge)
        lower = 0.10 * q_meas
        upper = min(0.70 * q_meas, 0.90 * float(sim_q[-1]), float(exp_q[-1]))
        if upper <= lower:
            raise RuntimeError(f"No fixed voltage grid for {rate:g}C")
        grids[rate] = np.linspace(lower, upper, N_VOLTAGE_POINTS)
    return grids


def interpolate_on_q(data, grid, q_meas, charge):
    q = q_axis(data, q_meas, charge)
    v = np.asarray(data["V"], dtype=float)
    order = np.argsort(q, kind="stable")
    q_unique, index = np.unique(q[order], return_index=True)
    if grid[-1] > q_unique[-1] + 1e-10:
        raise RuntimeError("Simulation ended before the fixed voltage grid")
    return np.interp(grid, q_unique, v[order][index])


def run_direction(kn_value, charge, baseline, stage, anode, cathode):
    model = model_with_kn(baseline, kn_value)
    runs = {}
    for rate in ga.RATES:
        print(
            f"  kn={kn_value:.5e} | {rate:g}C {'charge' if charge else 'discharge'}",
            flush=True,
        )
        runs[(rate, charge)] = ehq.run_dfn(
            rate, charge, stage, model, anode, cathode, True, True
        )
    return runs


def direction_metrics(kn_value, charge, runs, experiment, grids, q_meas):
    point_rows = []
    condition_rows = []
    for rate in ga.RATES:
        grid = grids[rate]
        exp = experiment[(rate, charge)]
        sim = runs[(rate, charge)]
        exp_v = interpolate_on_q(exp, grid, q_meas, charge)
        sim_v = interpolate_on_q(sim, grid, q_meas, charge)
        residual_v = sim_v - exp_v
        q_exp = end_capacity(exp, q_meas, charge)
        q_sim = end_capacity(sim, q_meas, charge)
        rel_capacity_error = (q_sim - q_exp) / q_exp
        condition_rows.append(
            {
                "kn": kn_value,
                "C_rate": rate,
                "Direction": "Charge" if charge else "Discharge",
                "MSE_V2": float(np.mean(residual_v**2)),
                "RMSE_mV": float(np.sqrt(np.mean(residual_v**2)) * 1000.0),
                "MAE_mV": float(np.mean(np.abs(residual_v)) * 1000.0),
                "Bias_mV": float(np.mean(residual_v) * 1000.0),
                "MaxAbs_mV": float(np.max(np.abs(residual_v)) * 1000.0),
                "Q_end_exp_Ah": q_exp,
                "Q_end_model_Ah": q_sim,
                "Capacity_error_mAh": (q_sim - q_exp) * 1000.0,
                "Capacity_error_pct": rel_capacity_error * 100.0,
            }
        )
        for q, exp_value, sim_value, residual in zip(grid, exp_v, sim_v, residual_v):
            point_rows.append(
                {
                    "kn": kn_value,
                    "C_rate": rate,
                    "Direction": "Charge" if charge else "Discharge",
                    "Q_Ah": q,
                    "Experiment_V": exp_value,
                    "Model_V": sim_value,
                    "Residual_mV": residual * 1000.0,
                }
            )
    conditions = pd.DataFrame(condition_rows)
    summary = {
        "kn": kn_value,
        # Equal weight for each C-rate curve, independent of sample count.
        "Voltage_MSE_V2": float(conditions.MSE_V2.mean()),
        "Voltage_RMSE_mV": float(np.sqrt(conditions.MSE_V2.mean()) * 1000.0),
        "Voltage_MAE_mV": float(conditions.MAE_mV.mean()),
        "Voltage_abs_bias_mV": float(conditions.Bias_mV.abs().mean()),
        "Capacity_RMSE_pct": float(np.sqrt(np.mean(conditions.Capacity_error_pct**2))),
        "Capacity_MAE_pct": float(conditions.Capacity_error_pct.abs().mean()),
        "Capacity_MAE_mAh": float(conditions.Capacity_error_mAh.abs().mean()),
        "Capacity_MaxAbs_mAh": float(conditions.Capacity_error_mAh.abs().max()),
    }
    return summary, conditions, pd.DataFrame(point_rows)


def pareto_mask(frame):
    values = frame[["Voltage_MSE_V2", "Capacity_RMSE_pct"]].to_numpy(dtype=float)
    mask = np.ones(len(values), dtype=bool)
    for i, value in enumerate(values):
        dominated = np.any(
            np.all(values <= value, axis=1) & np.any(values < value, axis=1)
        )
        mask[i] = not dominated
    return mask


def choose_methods(scan):
    methods = []
    voltage = scan.loc[scan.Voltage_MSE_V2.idxmin()]
    capacity = scan.loc[scan.Capacity_RMSE_pct.idxmin()]
    methods.append(("Voltage-only LS", voltage))
    methods.append(("Capacity-only reference", capacity))

    pareto = scan[scan.Pareto].sort_values("Capacity_RMSE_pct").copy()
    for column in ("Voltage_MSE_V2", "Capacity_RMSE_pct"):
        span = float(pareto[column].max() - pareto[column].min())
        pareto[f"norm_{column}"] = (
            (pareto[column] - pareto[column].min()) / span if span else 0.0
        )
    pareto["Utopia_distance"] = np.sqrt(
        pareto["norm_Voltage_MSE_V2"] ** 2
        + pareto["norm_Capacity_RMSE_pct"] ** 2
    )
    knee = pareto.loc[pareto.Utopia_distance.idxmin()]
    methods.append(("Pareto knee (descriptive)", knee))

    feasibility = []
    for tolerance in CAPACITY_CONSTRAINTS_PCT:
        feasible = scan[scan.Capacity_RMSE_pct <= tolerance]
        if feasible.empty:
            feasibility.append(
                {
                    "Capacity_constraint_pct": tolerance,
                    "Feasible": False,
                    "kn": np.nan,
                    "Voltage_RMSE_mV": np.nan,
                    "Capacity_RMSE_pct": np.nan,
                }
            )
            continue
        selected = feasible.loc[feasible.Voltage_MSE_V2.idxmin()]
        label = f"Capacity-constrained {tolerance:g}%"
        methods.append((label, selected))
        feasibility.append(
            {
                "Capacity_constraint_pct": tolerance,
                "Feasible": True,
                "kn": float(selected.kn),
                "Voltage_RMSE_mV": float(selected.Voltage_RMSE_mV),
                "Capacity_RMSE_pct": float(selected.Capacity_RMSE_pct),
            }
        )

    method_rows = []
    for label, row in methods:
        method_rows.append(
            {
                "Method": label,
                **{
                    key: row[key]
                    for key in (
                        "kn",
                        "Voltage_MSE_V2",
                        "Voltage_RMSE_mV",
                        "Voltage_MAE_mV",
                        "Capacity_RMSE_pct",
                        "Capacity_MAE_pct",
                        "Capacity_MAE_mAh",
                        "Capacity_MaxAbs_mAh",
                    )
                },
            }
        )
    return pd.DataFrame(method_rows).drop_duplicates("Method"), pd.DataFrame(feasibility), pareto


def plot_tradeoff(scan, methods):
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.0), constrained_layout=True)
    ax = axes[0]
    ax.semilogx(scan.kn, scan.Voltage_RMSE_mV, "o-", label="Charge voltage RMSE")
    ax.set_xlabel("kn")
    ax.set_ylabel("Charge voltage RMSE (mV)")
    ax.grid(alpha=0.25)
    ax2 = ax.twinx()
    ax2.semilogx(scan.kn, scan.Capacity_RMSE_pct, "s--", color="#D55E00", label="Charge capacity RMSE")
    ax2.set_ylabel("Charge capacity RMSE (%)", color="#D55E00")
    ax.set_title("Objective scan")

    ax = axes[1]
    ax.scatter(scan.Capacity_RMSE_pct, scan.Voltage_RMSE_mV, c="#999999", s=30, label="Dominated")
    pareto = scan[scan.Pareto].sort_values("Capacity_RMSE_pct")
    ax.plot(pareto.Capacity_RMSE_pct, pareto.Voltage_RMSE_mV, "o-", color="#0072B2", label="Pareto front")
    for _, row in methods.iterrows():
        if row.Method in ("Voltage-only LS", "Capacity-only reference", "Pareto knee (descriptive)"):
            ax.annotate(
                row.Method.replace(" (descriptive)", ""),
                (row.Capacity_RMSE_pct, row.Voltage_RMSE_mV),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize=8,
            )
    ax.set_xlabel("Charge capacity RMSE (%)")
    ax.set_ylabel("Charge voltage RMSE (mV)")
    ax.set_title("Voltage-capacity Pareto trade-off")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(RESULTS / "kn_charge_objective_and_pareto.png", dpi=220)
    plt.close(fig)


def plot_curves(experiment, runs_by_method, condition_metrics, q_meas, charge, filename):
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)
    colors = {
        "Baseline EIS kn": "#666666",
        "Voltage-only LS": "#0072B2",
        "Pareto knee (descriptive)": "#D55E00",
    }
    preferred = [name for name in colors if name in runs_by_method]
    for ax, rate in zip(axes, ga.RATES):
        exp = experiment[(rate, charge)]
        ax.plot(q_axis(exp, q_meas, charge), exp["V"], color="black", lw=2.4, label="Experiment")
        for method in preferred:
            sim = runs_by_method[method][(rate, charge)]
            metric = condition_metrics[
                (condition_metrics.Method == method)
                & np.isclose(condition_metrics.C_rate, rate)
            ].iloc[0]
            ax.plot(
                q_axis(sim, q_meas, charge),
                sim["V"],
                color=colors[method],
                lw=1.7,
                label=f"{method}: {metric.RMSE_mV:.1f} mV",
            )
        ax.set_title(f"{rate:g}C {'charge fit' if charge else 'discharge validation'}")
        ax.set_xlabel("Transferred capacity (Ah)")
        ax.set_ylabel("Voltage (V)")
        ax.grid(alpha=0.2)
    axes[0].legend(fontsize=8)
    fig.savefig(RESULTS / filename, dpi=220)
    plt.close(fig)


def main():
    warnings.filterwarnings(
        "ignore", message="The definition of the hysteresis decay rate parameter has changed"
    )
    _, baseline, stage, anode, cathode, experiment, ocp_row = stage2.build_current_inputs()
    q_meas = float(baseline["q_meas"])

    print("Baseline charge/discharge for fixed grids", flush=True)
    baseline_charge = run_direction(current.KN_PREF, True, baseline, stage, anode, cathode)
    baseline_discharge = run_direction(current.KN_PREF, False, baseline, stage, anode, cathode)
    charge_grids = fixed_voltage_grids(experiment, baseline_charge, q_meas, True)
    discharge_grids = fixed_voltage_grids(experiment, baseline_discharge, q_meas, False)

    kn_values = np.unique(
        np.concatenate(
            [
                np.geomspace(KN_MIN, KN_MAX, N_SCAN),
                # Dense refinement around the coarse voltage/capacity optima.
                np.linspace(0.80e-6, 2.10e-6, 34),
                np.asarray([current.KN_PREF]),
            ]
        )
    )
    scan_rows = []
    charge_condition_frames = []
    charge_point_frames = []
    charge_cache = {float(current.KN_PREF): baseline_charge}
    for kn_value in kn_values:
        try:
            runs = charge_cache.get(float(kn_value))
            if runs is None:
                runs = run_direction(kn_value, True, baseline, stage, anode, cathode)
                charge_cache[float(kn_value)] = runs
            summary, conditions, points = direction_metrics(
                kn_value, True, runs, experiment, charge_grids, q_meas
            )
            summary["Valid"] = True
            scan_rows.append(summary)
            charge_condition_frames.append(conditions)
            charge_point_frames.append(points)
        except Exception as exc:
            print(f"  invalid kn={kn_value:.5e}: {exc}", flush=True)
            scan_rows.append({"kn": kn_value, "Valid": False, "Error": str(exc)})

    scan = pd.DataFrame(scan_rows)
    valid = scan[scan.Valid].copy().reset_index(drop=True)
    valid["Pareto"] = pareto_mask(valid)
    scan = scan.merge(valid[["kn", "Pareto"]], on="kn", how="left")
    scan["Pareto"] = scan.Pareto.fillna(False)
    scan.to_csv(RESULTS / "kn_charge_scan.csv", index=False, encoding="utf-8-sig")

    methods, constraints, pareto = choose_methods(valid)
    baseline_summary, baseline_conditions, baseline_points = direction_metrics(
        current.KN_PREF,
        True,
        baseline_charge,
        experiment,
        charge_grids,
        q_meas,
    )
    baseline_row = {"Method": "Baseline EIS kn", **baseline_summary}
    methods = pd.concat([pd.DataFrame([baseline_row]), methods], ignore_index=True)
    methods.to_csv(RESULTS / "charge_selected_methods.csv", index=False, encoding="utf-8-sig")
    constraints.to_csv(RESULTS / "capacity_constraint_sweep.csv", index=False, encoding="utf-8-sig")
    pareto.to_csv(RESULTS / "kn_charge_pareto_front.csv", index=False, encoding="utf-8-sig")

    charge_conditions_all = pd.concat(
        [baseline_conditions.assign(Method="Baseline EIS kn")]
        + [
            pd.concat(charge_condition_frames, ignore_index=True)
            .loc[
                lambda x, k=float(row.kn): np.isclose(
                    x.kn, k, rtol=1e-12, atol=0.0
                )
            ]
            .assign(Method=row.Method)
            for _, row in methods[methods.Method != "Baseline EIS kn"].iterrows()
        ],
        ignore_index=True,
    )
    charge_conditions_all.to_csv(
        RESULTS / "charge_selected_condition_metrics.csv", index=False, encoding="utf-8-sig"
    )

    print("Discharge validation for selected unique kn values", flush=True)
    discharge_cache = {float(current.KN_PREF): baseline_discharge}
    discharge_summaries = []
    discharge_conditions = []
    discharge_points = []
    runs_by_method = {"Baseline EIS kn": baseline_discharge}
    charge_runs_by_method = {"Baseline EIS kn": baseline_charge}
    for _, row in methods.iterrows():
        method = str(row.Method)
        kn_value = float(row.kn)
        charge_runs_by_method[method] = charge_cache[kn_value]
        runs = discharge_cache.get(kn_value)
        if runs is None:
            runs = run_direction(kn_value, False, baseline, stage, anode, cathode)
            discharge_cache[kn_value] = runs
        runs_by_method[method] = runs
        summary, conditions, points = direction_metrics(
            kn_value, False, runs, experiment, discharge_grids, q_meas
        )
        discharge_summaries.append({"Method": method, **summary})
        discharge_conditions.append(conditions.assign(Method=method))
        discharge_points.append(points.assign(Method=method))

    discharge_summary = pd.DataFrame(discharge_summaries).drop_duplicates("Method")
    discharge_condition = pd.concat(discharge_conditions, ignore_index=True)
    discharge_point = pd.concat(discharge_points, ignore_index=True)
    discharge_summary.to_csv(
        RESULTS / "discharge_validation_summary.csv", index=False, encoding="utf-8-sig"
    )
    discharge_condition.to_csv(
        RESULTS / "discharge_validation_condition_metrics.csv", index=False, encoding="utf-8-sig"
    )
    discharge_point.to_csv(
        RESULTS / "discharge_validation_points.csv", index=False, encoding="utf-8-sig"
    )

    combined = methods.merge(
        discharge_summary,
        on=["Method", "kn"],
        how="left",
        suffixes=("_ChargeFit", "_DischargeValidation"),
    )
    combined.to_csv(RESULTS / "fit_validation_method_comparison.csv", index=False, encoding="utf-8-sig")

    plot_tradeoff(valid, methods[methods.Method != "Baseline EIS kn"])
    plot_curves(
        experiment,
        charge_runs_by_method,
        charge_conditions_all,
        q_meas,
        True,
        "charge_fit_curves.png",
    )
    plot_curves(
        experiment,
        runs_by_method,
        discharge_condition,
        q_meas,
        False,
        "discharge_validation_curves.png",
    )

    report = {
        "fixed_configuration": {
            "qOCV_MAE_mV": float(ocp_row.qOCV_MAE_2_98_mV),
            "kn_EIS": current.KN_PREF,
            "kp_fixed": current.KP_PREF,
            "Dsn_fixed_m2_s": 2.1e-14,
            "Dsp_fixed_m2_s": 4.4e-14,
            "brugg_n_p_s_fixed": [current.BRUGG_N, current.BRUGG_P, current.BRUGG_S],
        },
        "objective": {
            "primary": "equal-curve ordinary least squares of charge voltage on fixed transferred-capacity grids",
            "equation": "mean_over_rates(mean_over_Q((V_model - V_exp)^2))",
            "capacity": "reported separately; Pareto front; capacity-RMSE constraint sweep",
            "prior": "none",
            "kn_bounds": [KN_MIN, KN_MAX],
            "capacity_constraints_pct": list(CAPACITY_CONSTRAINTS_PCT),
        },
        "selected_methods": methods.to_dict(orient="records"),
        "discharge_validation": discharge_summary.to_dict(orient="records"),
    }
    (RESULTS / "stage3_kn_fit_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print("\nCharge selected methods")
    print(
        methods[
            ["Method", "kn", "Voltage_RMSE_mV", "Capacity_RMSE_pct", "Capacity_MAE_mAh"]
        ].to_string(index=False)
    )
    print("\nDischarge validation")
    print(
        discharge_summary[
            ["Method", "kn", "Voltage_RMSE_mV", "Capacity_RMSE_pct", "Capacity_MAE_mAh"]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
