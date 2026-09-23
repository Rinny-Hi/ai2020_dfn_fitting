"""Absolute-Ah C/50 reconstruction, eSOH fitting, and held-out DFN validation.

The previous C/50 target normalized every charge/discharge branch to its own
capacity before averaging.  This audit keeps the measured Ah coordinate, fits
Qn/Qp and electrode alignment without a weighted multi-term objective, and
uses the 0.5C/1C/2C curves only for held-out validation.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution, minimize
from scipy.signal import savgol_filter

import c50_ocp_source_combination_study as c50
import c50_selected_ocp_dynamic_validation as selected
import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import nominal_radius_current_configuration as current
import priority_initial_state_protocol_recheck as priority


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260923_absolute_capacity_esoh_recheck"
RESULTS.mkdir(parents=True, exist_ok=True)
NOMINAL_CAPACITY_AH = 2.28


def load_absolute_branches() -> tuple[list[dict[str, Any]], pd.DataFrame]:
    curves: list[dict[str, Any]] = []
    rows: list[dict[str, Any]] = []
    for file_name in c50.JUNE_FILES:
        frame = pd.read_excel(file_name, sheet_name="record")
        step_type = frame["Step Type"].astype(str)
        group = step_type.ne(step_type.shift()).cumsum()
        found: dict[str, list[dict[str, Any]]] = {"CC Chg": [], "CC DChg": []}
        for _, part in frame.groupby(group):
            kind = str(part["Step Type"].iloc[0])
            if kind not in found:
                continue
            q = pd.to_numeric(part["Capacity(Ah)"], errors="coerce").to_numpy(float)
            v = pd.to_numeric(part["Voltage(V)"], errors="coerce").to_numpy(float)
            current_a = float(pd.to_numeric(part["Current(A)"], errors="coerce").median())
            valid = np.isfinite(q) & np.isfinite(v)
            q, v = q[valid], v[valid]
            if len(q) < 50 or float(np.max(q)) < 2.2 or abs(current_a) > 0.06:
                continue
            q = q - q[0]
            qmax = float(np.max(q))
            # Both directions use capacity measured from the empty state.
            q_empty = q if kind == "CC Chg" else qmax - q
            order = np.argsort(q_empty, kind="stable")
            q_unique, index = np.unique(q_empty[order], return_index=True)
            found[kind].append(
                {
                    "q": q_unique,
                    "v": v[order][index],
                    "capacity_Ah": qmax,
                    "current_A": current_a,
                }
            )

        selected_groups = {
            "Charge": found["CC Chg"],
            # The first discharge follows a partial initial charge.
            "Discharge": found["CC DChg"][1:],
        }
        for direction, entries in selected_groups.items():
            for cycle_index, entry in enumerate(entries, start=1):
                entry = dict(entry)
                entry.update(
                    {
                        "file": file_name.name,
                        "direction": direction,
                        "selected_cycle_index": cycle_index,
                    }
                )
                curves.append(entry)
                rows.append(
                    {
                        "file": file_name.name,
                        "direction": direction,
                        "selected_cycle_index": cycle_index,
                        "capacity_Ah": entry["capacity_Ah"],
                        "current_A": entry["current_A"],
                    }
                )
    return curves, pd.DataFrame(rows)


def absolute_target(curves: list[dict[str, Any]], q_ref: float) -> dict[str, Any]:
    q = np.linspace(0.0, q_ref, 1201)
    by_direction: dict[str, np.ndarray] = {}
    counts: dict[str, np.ndarray] = {}
    for direction in ("Charge", "Discharge"):
        matrix = []
        for curve in curves:
            if curve["direction"] != direction:
                continue
            values = np.interp(q, curve["q"], curve["v"])
            values[q > float(curve["capacity_Ah"])] = np.nan
            matrix.append(values)
        array = np.vstack(matrix)
        by_direction[direction] = np.nanmean(array, axis=0)
        counts[direction] = np.sum(np.isfinite(array), axis=0)
    qocv = 0.5 * (by_direction["Charge"] + by_direction["Discharge"])
    return {
        "q_Ah": q,
        "charge_V": by_direction["Charge"],
        "discharge_V": by_direction["Discharge"],
        "qocv_V": qocv,
        "charge_count": counts["Charge"],
        "discharge_count": counts["Discharge"],
    }


def normalized_target_from_same_branches(
    curves: list[dict[str, Any]], soc: np.ndarray
) -> np.ndarray:
    """Reproduce the previous per-branch normalization on the same branches."""
    directional = []
    for direction in ("Charge", "Discharge"):
        matrix = []
        for curve in curves:
            if curve["direction"] != direction:
                continue
            branch_soc = np.asarray(curve["q"], float) / float(curve["capacity_Ah"])
            matrix.append(np.interp(soc, branch_soc, curve["v"]))
        directional.append(np.mean(np.vstack(matrix), axis=0))
    return 0.5 * (directional[0] + directional[1])


def predict_absolute(
    q: np.ndarray,
    qn: float,
    qp: float,
    x0: float,
    y0: float,
    anode: dict[str, Any],
    cathode: dict[str, Any],
) -> np.ndarray:
    x = x0 + q / qn
    y = y0 - q / qp
    un = np.interp(x, anode["grid"], anode["equilibrium"])
    up = np.interp(y, cathode["grid"], cathode["equilibrium"])
    return up - un


def fit_case(
    name: str,
    bound_fraction: float,
    target: dict[str, Any],
    anode: dict[str, Any],
    cathode: dict[str, Any],
    baseline: pd.Series,
    seed: int,
) -> tuple[dict[str, Any], np.ndarray]:
    q = target["q_Ah"]
    measured = target["qocv_V"]
    q_ref = float(q[-1])
    mask = np.isfinite(measured) & (q >= 0.02 * q_ref) & (q <= 0.98 * q_ref)
    derivative_mask = np.isfinite(measured) & (q >= 0.05 * q_ref) & (q <= 0.95 * q_ref)
    xmin, xmax = float(np.min(anode["grid"])), float(np.max(anode["grid"]))
    ymin, ymax = float(np.min(cathode["grid"])), float(np.max(cathode["grid"]))
    qn0, qp0 = float(baseline.Qn_Ah), float(baseline.Qp_Ah)
    if bound_fraction == 0:
        bounds = [(xmin, xmax), (ymin, ymax)]

        def decode(z: np.ndarray) -> tuple[float, float, float, float]:
            return qn0, qp0, float(z[0]), float(z[1])

    else:
        bounds = [
            (qn0 * (1 - bound_fraction), qn0 * (1 + bound_fraction)),
            (qp0 * (1 - bound_fraction), qp0 * (1 + bound_fraction)),
            (xmin, xmax),
            (ymin, ymax),
        ]

        def decode(z: np.ndarray) -> tuple[float, float, float, float]:
            return tuple(map(float, z))  # type: ignore[return-value]

    def coverage(z: np.ndarray) -> np.ndarray:
        qn, qp, x0, y0 = decode(z)
        return np.array([x0 - xmin, xmax - (x0 + q_ref / qn), ymin * -1 + (y0 - q_ref / qp), ymax - y0])

    def voltage_mse(z: np.ndarray) -> float:
        qn, qp, x0, y0 = decode(z)
        voltage = predict_absolute(q, qn, qp, x0, y0, anode, cathode)
        return float(np.mean((voltage[mask] - measured[mask]) ** 2))

    def endpoint_error(z: np.ndarray) -> np.ndarray:
        qn, qp, x0, y0 = decode(z)
        voltage = predict_absolute(q[[0, -1]], qn, qp, x0, y0, anode, cathode)
        return voltage - measured[[0, -1]]

    def violation(z: np.ndarray) -> float:
        terms = [max(0.0, -float(v)) for v in coverage(z)]
        terms += [max(0.0, abs(float(v)) - 0.010) for v in endpoint_error(z)]
        return float(np.sum(np.square(terms)))

    def penalized(z: np.ndarray) -> float:
        return voltage_mse(z) + 1e5 * violation(z)

    de = differential_evolution(
        penalized,
        bounds=bounds,
        seed=seed,
        popsize=15,
        maxiter=240,
        tol=1e-10,
        polish=False,
        workers=1,
    )
    constraints = [
        {"type": "ineq", "fun": lambda z, i=i: coverage(z)[i]}
        for i in range(4)
    ] + [
        {"type": "ineq", "fun": lambda z, i=i: 0.010 - abs(float(endpoint_error(z)[i]))}
        for i in range(2)
    ]
    local = minimize(
        voltage_mse,
        de.x,
        method="SLSQP",
        bounds=bounds,
        constraints=constraints,
        options={"maxiter": 3000, "ftol": 1e-14, "disp": False},
    )
    z = local.x if violation(local.x) <= violation(de.x) + 1e-12 else de.x
    qn, qp, x0, y0 = decode(z)
    x100 = x0 + q_ref / qn
    y100 = y0 - q_ref / qp
    voltage = predict_absolute(q, qn, qp, x0, y0, anode, cathode)
    err_mv = (voltage - measured) * 1000.0
    measured_smooth = savgol_filter(measured, 51, 3)
    voltage_smooth = savgol_filter(voltage, 51, 3)
    dvdq_err = np.gradient(voltage_smooth, q) - np.gradient(measured_smooth, q)
    row = {
        "Candidate": name,
        "capacity_bound_percent": 100 * bound_fraction,
        "Q_n_Ah": qn,
        "Q_p_Ah": qp,
        "Q_Li_proxy_Ah": qn * x0 + qp * y0,
        "Q_n_change_percent": 100 * (qn / qn0 - 1),
        "Q_p_change_percent": 100 * (qp / qp0 - 1),
        "x0": x0,
        "x100": x100,
        "y100": y100,
        "y0": y0,
        "delta_x": x100 - x0,
        "delta_y": y0 - y100,
        "qOCV_MAE_2_98_mV": float(np.mean(np.abs(err_mv[mask]))),
        "qOCV_RMSE_2_98_mV": float(np.sqrt(np.mean(err_mv[mask] ** 2))),
        "qOCV_endpoint_empty_error_mV": float(err_mv[0]),
        "qOCV_endpoint_full_error_mV": float(err_mv[-1]),
        "qOCV_max_endpoint_abs_mV": float(max(abs(err_mv[0]), abs(err_mv[-1]))),
        "dVdQ_RMSE_5_95_V_per_Ah": float(np.sqrt(np.mean(dvdq_err[derivative_mask] ** 2))),
        "coverage_min": float(np.min(coverage(z))),
        "constraint_violation": violation(z),
        "optimizer_success": bool(local.success),
        "optimizer_message": str(local.message),
    }
    return row, voltage


def baseline_on_absolute_grid(
    target: dict[str, Any], anode: dict[str, Any], cathode: dict[str, Any], baseline: pd.Series
) -> tuple[dict[str, Any], np.ndarray]:
    q = target["q_Ah"]
    q_ref = float(q[-1])
    x0, y0 = float(baseline.x0), float(baseline.y0)
    qn, qp = float(baseline.Qn_Ah), float(baseline.Qp_Ah)
    voltage = predict_absolute(q, qn, qp, x0, y0, anode, cathode)
    measured = target["qocv_V"]
    mask = np.isfinite(measured) & (q >= 0.02 * q_ref) & (q <= 0.98 * q_ref)
    dmask = np.isfinite(measured) & (q >= 0.05 * q_ref) & (q <= 0.95 * q_ref)
    err = (voltage - measured) * 1000
    dvdq_err = np.gradient(savgol_filter(voltage, 51, 3), q) - np.gradient(savgol_filter(measured, 51, 3), q)
    row = {
        "Candidate": "Previous normalized-SOC solution",
        "capacity_bound_percent": np.nan,
        "Q_n_Ah": qn,
        "Q_p_Ah": qp,
        "Q_Li_proxy_Ah": qn * x0 + qp * y0,
        "Q_n_change_percent": 0.0,
        "Q_p_change_percent": 0.0,
        "x0": x0,
        "x100": x0 + q_ref / qn,
        "y100": y0 - q_ref / qp,
        "y0": y0,
        "delta_x": q_ref / qn,
        "delta_y": q_ref / qp,
        "qOCV_MAE_2_98_mV": float(np.mean(np.abs(err[mask]))),
        "qOCV_RMSE_2_98_mV": float(np.sqrt(np.mean(err[mask] ** 2))),
        "qOCV_endpoint_empty_error_mV": float(err[0]),
        "qOCV_endpoint_full_error_mV": float(err[-1]),
        "qOCV_max_endpoint_abs_mV": float(max(abs(err[0]), abs(err[-1]))),
        "dVdQ_RMSE_5_95_V_per_Ah": float(np.sqrt(np.mean(dvdq_err[dmask] ** 2))),
        "coverage_min": np.nan,
        "constraint_violation": np.nan,
        "optimizer_success": True,
        "optimizer_message": "existing control",
    }
    return row, voltage


def model_for_candidate(model0: dict[str, Any], baseline: pd.Series, row: pd.Series, q_ref: float) -> dict[str, Any]:
    model = dict(model0)
    model["params"] = model0["params"].copy()
    csn0 = float(model0["params"]["Maximum concentration in negative electrode [mol.m-3]"])
    csp0 = float(model0["params"]["Maximum concentration in positive electrode [mol.m-3]"])
    model["params"].update(
        {
            "Maximum concentration in negative electrode [mol.m-3]": csn0 * float(row.Q_n_Ah) / float(baseline.Qn_Ah),
            "Maximum concentration in positive electrode [mol.m-3]": csp0 * float(row.Q_p_Ah) / float(baseline.Qp_Ah),
        },
        check_already_exists=False,
    )
    model["q_meas"] = q_ref
    model["qn_cell"] = float(row.Q_n_Ah)
    model["qp_cell"] = float(row.Q_p_Ah)
    model["delta_x"] = float(row.delta_x)
    model["delta_y"] = float(row.delta_y)
    soc = np.linspace(0, 1, 1001)
    model["soc"] = soc
    model["v_qocv"] = predict_absolute(
        soc * q_ref,
        float(row.Q_n_Ah),
        float(row.Q_p_Ah),
        float(row.x0),
        float(row.y0),
        model["anode_detail"],
        model["cathode_detail"],
    )
    return model


def equilibrium_from_transferred(
    q: np.ndarray,
    charge: bool,
    initial_soc: float,
    row: pd.Series,
    q_ref: float,
    anode: dict[str, Any],
    cathode: dict[str, Any],
) -> np.ndarray:
    q0 = initial_soc * q_ref
    absolute_q = q0 + q if charge else q0 - q
    absolute_q = np.clip(absolute_q, 0, q_ref)
    return predict_absolute(
        absolute_q,
        float(row.Q_n_Ah),
        float(row.Q_p_Ah),
        float(row.x0),
        float(row.y0),
        anode,
        cathode,
    )


def validate_candidates(
    candidates: pd.DataFrame,
    curves: dict[str, np.ndarray],
    baseline: pd.Series,
    model0: dict[str, Any],
    stage0: dict[str, float],
    anode: dict[str, Any],
    cathode: dict[str, Any],
    experiment: dict[tuple[float, bool], dict[str, np.ndarray]],
    q_ref: float,
) -> tuple[pd.DataFrame, dict[tuple[str, float, bool], dict[str, np.ndarray]], pd.DataFrame]:
    protocol, _, _ = priority.load_protocol(model0)
    pidx = protocol.set_index(["C_rate", "Direction"])
    rows: list[dict[str, Any]] = []
    decomposition: list[dict[str, Any]] = []
    simulations: dict[tuple[str, float, bool], dict[str, np.ndarray]] = {}
    model0 = dict(model0)
    model0["anode_detail"] = anode
    model0["cathode_detail"] = cathode

    for candidate_name in candidates.Candidate:
        row = candidates[candidates.Candidate == candidate_name].iloc[0]
        model = model_for_candidate(model0, baseline, row, q_ref)
        stage = {key: float(row[key]) for key in ("x0", "x100", "y100", "y0")}
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                rest_v = float(pidx.loc[(rate, direction), "Rest_end_V"])
                initial_soc = priority.monotone_voltage_inverse(model["v_qocv"], model["soc"], rest_v)
                local_stage = priority.stage_at_soc(stage, initial_soc, charge)
                measured_current = float(np.nanmedian(np.abs(experiment[(rate, charge)]["I_A"])))
                sim_rate = measured_current / NOMINAL_CAPACITY_AH
                print(f"{candidate_name} | {rate:g}C {direction}", flush=True)
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
                simulations[(candidate_name, rate, charge)] = sim
                metric = priority.curve_metrics(experiment[(rate, charge)], sim, q_ref, charge)
                rows.append(
                    {
                        "Candidate": candidate_name,
                        "C_rate": rate,
                        "Direction": direction,
                        "initial_SOC": initial_soc,
                        "measured_current_A": measured_current,
                        **metric,
                        "Capacity_error_mAh": 1000 * (metric["Q_end_model_Ah"] - metric["Q_end_exp_Ah"]),
                    }
                )

                qe = priority.transferred_capacity(experiment[(rate, charge)], q_ref, charge)
                qm = priority.transferred_capacity(sim, q_ref, charge)
                ve = np.asarray(experiment[(rate, charge)]["V"], float)
                vm = np.asarray(sim["V"], float)
                upper = 0.995 * min(float(np.max(qe)), float(np.max(qm)))
                grid = np.linspace(0, upper, 301)
                exp_v = np.interp(grid, qe, ve)
                mod_v = np.interp(grid, qm, vm)
                eq_v = equilibrium_from_transferred(grid, charge, initial_soc, row, q_ref, anode, cathode)
                exp_pol = exp_v - eq_v
                mod_pol = mod_v - eq_v
                pol_err = (mod_pol - exp_pol) * 1000
                for fraction in (0.1, 0.5, 0.9):
                    index = int(round(fraction * (len(grid) - 1)))
                    decomposition.append(
                        {
                            "Candidate": candidate_name,
                            "C_rate": rate,
                            "Direction": direction,
                            "fraction_of_common_range": fraction,
                            "capacity_Ah": grid[index],
                            "experiment_polarization_mV": exp_pol[index] * 1000,
                            "model_polarization_mV": mod_pol[index] * 1000,
                            "polarization_error_mV": pol_err[index],
                        }
                    )
                rows[-1]["Polarization_RMSE_mV"] = float(np.sqrt(np.mean(pol_err**2)))
                rows[-1]["Polarization_bias_mV"] = float(np.mean(pol_err))
    return pd.DataFrame(rows), simulations, pd.DataFrame(decomposition)


def main() -> None:
    model0, stage0, anode, cathode, experiment, baseline, qcell = selected.selected_inputs()
    absolute_curves, branch_table = load_absolute_branches()
    q_ref = float(branch_table.capacity_Ah.mean())
    target = absolute_target(absolute_curves, q_ref)
    branch_table.to_csv(RESULTS / "absolute_capacity_branch_inventory.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(
        {
            "capacity_from_empty_Ah": target["q_Ah"],
            "charge_mean_V": target["charge_V"],
            "discharge_mean_V": target["discharge_V"],
            "qOCV_mean_V": target["qocv_V"],
            "charge_branch_count": target["charge_count"],
            "discharge_branch_count": target["discharge_count"],
        }
    ).to_csv(RESULTS / "absolute_capacity_c50_target.csv", index=False, encoding="utf-8-sig")

    rows: list[dict[str, Any]] = []
    ocv_curves: dict[str, np.ndarray] = {}
    row, curve = baseline_on_absolute_grid(target, anode, cathode, baseline)
    rows.append(row)
    ocv_curves[row["Candidate"]] = curve
    for i, (name, fraction) in enumerate(
        [
            ("Absolute Ah, Qn/Qp fixed", 0.0),
            ("Absolute Ah, Qn/Qp ±5%", 0.05),
            ("Absolute Ah, Qn/Qp ±15% diagnostic", 0.15),
        ]
    ):
        row, curve = fit_case(name, fraction, target, anode, cathode, baseline, 260923 + i)
        rows.append(row)
        ocv_curves[name] = curve
    fit_table = pd.DataFrame(rows)
    fit_table.to_csv(RESULTS / "absolute_capacity_esoh_candidates.csv", index=False, encoding="utf-8-sig")

    # Validate control, the no-capacity-change correction, and the physically
    # provisional ±5% eSOH candidate.  The ±15% case remains a diagnostic only.
    validation_names = [
        "Previous normalized-SOC solution",
        "Absolute Ah, Qn/Qp fixed",
        "Absolute Ah, Qn/Qp ±5%",
    ]
    validation_candidates = fit_table[fit_table.Candidate.isin(validation_names)].copy()
    dynamic, simulations, decomposition = validate_candidates(
        validation_candidates,
        ocv_curves,
        baseline,
        model0,
        stage0,
        anode,
        cathode,
        experiment,
        q_ref,
    )
    dynamic.to_csv(RESULTS / "heldout_dynamic_metrics.csv", index=False, encoding="utf-8-sig")
    decomposition.to_csv(RESULTS / "polarization_decomposition_points.csv", index=False, encoding="utf-8-sig")
    summary = dynamic.groupby(["Candidate", "Direction"], sort=False, as_index=False).agg(
        Mean_full_RMSE_mV=("Full_RMSE_mV", "mean"),
        Mean_center_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
        Mean_polarization_RMSE_mV=("Polarization_RMSE_mV", "mean"),
        Capacity_RMSE_pct=("Capacity_error_pct", lambda x: float(np.sqrt(np.mean(x**2)))),
        Mean_capacity_error_mAh=("Capacity_error_mAh", "mean"),
        Mean_abs_capacity_error_mAh=("Capacity_error_mAh", lambda x: float(np.mean(np.abs(x)))),
    )
    summary.to_csv(RESULTS / "heldout_dynamic_summary.csv", index=False, encoding="utf-8-sig")

    colors = {
        "Previous normalized-SOC solution": "#777777",
        "Absolute Ah, Qn/Qp fixed": "#0072B2",
        "Absolute Ah, Qn/Qp ±5%": "#D55E00",
        "Absolute Ah, Qn/Qp ±15% diagnostic": "#CC79A7",
    }
    q = target["q_Ah"]
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.0), constrained_layout=True)
    axes[0].plot(q, target["charge_V"], color="#E69F00", lw=1.4, label="C/50 charge mean")
    axes[0].plot(q, target["discharge_V"], color="#56B4E9", lw=1.4, label="C/50 discharge mean")
    axes[0].plot(q, target["qocv_V"], "k", lw=2.4, label="C/50 qOCV")
    for name, voltage in ocv_curves.items():
        axes[0].plot(q, voltage, color=colors[name], lw=1.8, label=name)
        axes[1].plot(q, (voltage - target["qocv_V"]) * 1000, color=colors[name], lw=1.7, label=name)
    measured_dvdq = np.gradient(savgol_filter(target["qocv_V"], 51, 3), q)
    axes[2].plot(q, measured_dvdq, "k", lw=2.3, label="C/50 qOCV")
    for name, voltage in ocv_curves.items():
        axes[2].plot(q, np.gradient(savgol_filter(voltage, 51, 3), q), color=colors[name], lw=1.6, label=name)
    axes[0].set(title="Absolute-capacity OCP fit", xlabel="Capacity from empty [Ah]", ylabel="Voltage [V]")
    axes[1].axhline(0, color="black", lw=0.8)
    axes[1].set(title="qOCV residual", xlabel="Capacity from empty [Ah]", ylabel="Model - experiment [mV]")
    axes[2].set(title="Differential-voltage check", xlabel="Capacity from empty [Ah]", ylabel="dV/dQ [V/Ah]", ylim=(-0.05, 1.0))
    for ax in axes:
        ax.grid(alpha=0.22)
    axes[0].legend(fontsize=7)
    fig.suptitle("C/50 reconstruction without per-branch SOC normalization")
    fig.savefig(RESULTS / "absolute_capacity_ocp_comparison.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(16.5, 8.8), constrained_layout=True)
    for col, rate in enumerate(ga.RATES):
        for r, charge in enumerate((True, False)):
            ax = axes[r, col]
            direction = "Charge" if charge else "Discharge"
            obs = experiment[(rate, charge)]
            q_exp = priority.transferred_capacity(obs, q_ref, charge)
            ax.plot(q_exp, obs["V"], "k", lw=2.5, label="Experiment")
            for name in validation_names:
                sim = simulations[(name, rate, charge)]
                metric = dynamic[(dynamic.Candidate == name) & np.isclose(dynamic.C_rate, rate) & (dynamic.Direction == direction)].iloc[0]
                q_sim = priority.transferred_capacity(sim, q_ref, charge)
                ax.plot(q_sim, sim["V"], color=colors[name], lw=1.7, label=f"{name}: {metric.Full_RMSE_mV:.1f} mV")
            ax.set(title=f"{rate:g}C {direction}", xlabel="Transferred capacity [Ah]", ylabel="Voltage [V]")
            ax.grid(alpha=0.22)
    axes[0, 0].legend(fontsize=6.6)
    fig.suptitle("Held-out 0.5C / 1C / 2C validation after absolute-Ah eSOH fit")
    fig.savefig(RESULTS / "heldout_dynamic_curves.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(15.8, 8.2), constrained_layout=True)
    for col, rate in enumerate(ga.RATES):
        for r, direction in enumerate(("Charge", "Discharge")):
            ax = axes[r, col]
            for name in validation_names:
                part = decomposition[
                    (decomposition.Candidate == name)
                    & np.isclose(decomposition.C_rate, rate)
                    & (decomposition.Direction == direction)
                ]
                ax.plot(
                    100 * part.fraction_of_common_range,
                    part.polarization_error_mV,
                    marker="o",
                    lw=1.8,
                    color=colors[name],
                    label=name,
                )
            ax.axhline(0, color="black", lw=0.8)
            ax.set(
                title=f"{rate:g}C {direction}",
                xlabel="Fraction of common capacity range [%]",
                ylabel="Model - measured polarization [mV]",
            )
            ax.grid(alpha=0.22)
    axes[0, 0].legend(fontsize=6.6)
    fig.suptitle("Polarization-residual decomposition at 10%, 50%, and 90%")
    fig.savefig(RESULTS / "polarization_residual_decomposition.png", dpi=220)
    plt.close(fig)

    template = pd.DataFrame(
        {
            "electrode": ["negative", "positive"],
            "SOC_or_stoichiometry": [np.nan, np.nan],
            "temperature_C": [25.0, 25.0],
            "Rct_ohm_cm2": [np.nan, np.nan],
            "electrolyte_concentration_mol_m3": [1000.0, 1000.0],
            "replicate_count": [np.nan, np.nan],
            "status": ["awaiting measured Rct", "awaiting measured Rct"],
        }
    )
    template.to_csv(RESULTS / "rct_measurement_input_template.csv", index=False, encoding="utf-8-sig")

    branch_stats = branch_table.groupby("direction").capacity_Ah.agg(["count", "mean", "std", "min", "max"]).reset_index()
    branch_stats.to_csv(RESULTS / "absolute_capacity_branch_statistics.csv", index=False, encoding="utf-8-sig")
    normalized_qocv = normalized_target_from_same_branches(absolute_curves, target["q_Ah"] / q_ref)
    normalization_error_mv = (normalized_qocv - target["qocv_V"]) * 1000
    normalization_mask = (target["q_Ah"] >= 0.02 * q_ref) & (target["q_Ah"] <= 0.98 * q_ref)
    report = {
        "absolute_capacity_reference_Ah": q_ref,
        "previous_normalized_capacity_Ah": float(qcell),
        "branch_capacity_range_mAh": 1000 * float(branch_table.capacity_Ah.max() - branch_table.capacity_Ah.min()),
        "per_branch_normalization_effect": {
            "MAE_2_98_mV": float(np.mean(np.abs(normalization_error_mv[normalization_mask]))),
            "RMSE_2_98_mV": float(np.sqrt(np.mean(normalization_error_mv[normalization_mask] ** 2))),
            "max_abs_2_98_mV": float(np.max(np.abs(normalization_error_mv[normalization_mask]))),
        },
        "fit_objective": "unweighted qOCV voltage MSE over 2-98% absolute Ah",
        "constraints": "electrode OCP-domain coverage and <=10 mV error at empty/full endpoints",
        "dVdQ_use": "diagnostic only; not in objective",
        "high_rate_use": "held-out validation only",
        "capacity_mapping_for_validation": "Qn/Qp changes represented by proportional effective c_s,max scaling; diagnostic until loading is measured",
        "fit_candidates": fit_table.to_dict(orient="records"),
        "dynamic_summary": summary.to_dict(orient="records"),
    }
    (RESULTS / "absolute_capacity_esoh_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print("\nBRANCH STATS")
    print(branch_stats.to_string(index=False))
    print("\nFIT CANDIDATES")
    print(fit_table.to_string(index=False))
    print("\nHELD-OUT SUMMARY")
    print(summary.to_string(index=False))
    print("Saved to", RESULTS)


if __name__ == "__main__":
    main()
