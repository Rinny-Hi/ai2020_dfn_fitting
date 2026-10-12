"""Fit identifiable transport combinations to dynamic and HPPC voltage data.

The local analysis shows kn and kp sensitivities are nearly collinear.  Their
measured baseline ratio is therefore preserved and one common kinetic scale is
fitted together with Dsn and Dsp.
"""

from __future__ import annotations

import json
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

import gitt_ocp_analysis as ga
import prefit_physical_anchor_decision as p0
import priority_initial_state_protocol_recheck as priority

import oregan_combined_sensitivity_identifiability as analysis
import electrolyte_transport_set_comparison as comparison


OUT = comparison.ROOT / "results" / "261012_oregan_identifiable_combined_fit"
FIT_NAMES = ("Dsn_scale", "Dsp_scale", "common_k_scale")
LOWER = np.log([0.5, 0.5, 0.3])
UPPER = np.log([2.0, 2.0, 3.0])


def values_from_log_scale(log_scale):
    scale = np.exp(np.asarray(log_scale, dtype=float))
    return {
        "Dsn": analysis.INITIAL["Dsn"] * scale[0],
        "Dsp": analysis.INITIAL["Dsp"] * scale[1],
        "kn": analysis.INITIAL["kn"] * scale[2],
        "kp": analysis.INITIAL["kp"] * scale[2],
        "Dsn_scale": scale[0],
        "Dsp_scale": scale[1],
        "common_k_scale": scale[2],
    }


def dynamic_prediction_and_missing(runs, grids, qcell):
    blocks = []
    missing_blocks = []
    for rate in ga.RATES:
        for charge in (True, False):
            source = runs[(rate, charge)]
            q = priority.transferred_capacity(source, qcell, charge)
            order = np.argsort(q, kind="stable")
            q_unique, index = np.unique(np.asarray(q)[order], return_index=True)
            voltage = np.asarray(source["V"], dtype=float)[order][index]
            grid = grids[(rate, charge)]
            blocks.append(np.interp(grid, q_unique, voltage))
            missing_blocks.append(np.maximum(grid - q_unique[-1], 0.0))
    return np.concatenate(blocks), np.concatenate(missing_blocks)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    hppc_targets = analysis.load_hppc_targets()
    hppc_times = analysis.hppc_target_times()
    hppc_exp = hppc_targets["Voltage_V"].to_numpy(float)
    baseline, stage, anode, cathode, experiment, qcell = comparison.current_p0()
    base_transport = comparison.with_transport_set(baseline, "ORegan2022")
    base_model = analysis.apply_parameters(base_transport, analysis.INITIAL)
    first_rest_voltage = float(hppc_targets.iloc[0]["Voltage_V"])
    hppc_z0 = priority.monotone_voltage_inverse(
        base_model["v_qocv"], base_model["soc"], first_rest_voltage
    )
    hppc_stage = priority.stage_at_soc(stage, hppc_z0, False)

    print("Preparing common objective", flush=True)
    baseline_dynamic_runs = analysis.run_dynamic(base_model, stage, anode, cathode)
    grids = analysis.make_dynamic_grids(experiment, baseline_dynamic_runs, qcell)
    dynamic_exp = analysis.dynamic_vector(experiment, grids, qcell, experiment=True)
    n_dynamic = len(dynamic_exp)
    n_hppc = len(hppc_exp)
    cache = {}
    history = []

    def evaluate(log_scale, keep_runs=False):
        clipped = np.clip(np.asarray(log_scale, dtype=float), LOWER, UPPER)
        key = tuple(np.round(clipped, 10))
        if key in cache and not keep_runs:
            return cache[key]
        values = values_from_log_scale(clipped)
        model_values = {name: values[name] for name in analysis.PARAMETERS}
        model = analysis.apply_parameters(base_transport, model_values)
        started = time.perf_counter()
        status = "ok"
        dynamic_runs = None
        try:
            dynamic_runs = analysis.run_dynamic(model, stage, anode, cathode)
            dynamic_model, missing = dynamic_prediction_and_missing(
                dynamic_runs, grids, qcell
            )
            dynamic_error = 1000.0 * (dynamic_model - dynamic_exp)
            if np.any(missing > 0.0):
                dynamic_error += 500.0 * missing / qcell
                status = "early_cutoff_penalty"
            hppc_model = analysis.run_hppc(
                model, hppc_stage, anode, cathode, hppc_times["Time_s"]
            )
            hppc_error = 1000.0 * (hppc_model - hppc_exp)
            residual = np.concatenate(
                [
                    dynamic_error / np.sqrt(n_dynamic),
                    hppc_error / np.sqrt(n_hppc),
                ]
            )
        except Exception as exc:
            residual = np.full(n_dynamic + n_hppc, 100.0)
            dynamic_error = np.full(n_dynamic, 1000.0)
            hppc_error = np.full(n_hppc, 1000.0)
            hppc_model = np.full(n_hppc, np.nan)
            status = "failed:{}".format(type(exc).__name__)
        record = {
            "Evaluation": len(history) + 1,
            "Status": status,
            "Dynamic_RMSE_mV": float(np.sqrt(np.mean(dynamic_error**2))),
            "HPPC_RMSE_mV": float(np.sqrt(np.mean(hppc_error**2))),
            "Equal_dataset_objective": float(np.sum(residual**2)),
            "Elapsed_s": time.perf_counter() - started,
            **values,
        }
        history.append(record)
        print(
            "eval {Evaluation}: dyn={Dynamic_RMSE_mV:.2f} mV, "
            "HPPC={HPPC_RMSE_mV:.2f} mV, scale=({Dsn_scale:.3f},"
            "{Dsp_scale:.3f},{common_k_scale:.3f}) [{Status}]".format(**record),
            flush=True,
        )
        pd.DataFrame(history).to_csv(
            OUT / "optimizer_history.csv", index=False, encoding="utf-8-sig"
        )
        result = (residual, dynamic_runs, hppc_model, record)
        cache[key] = result
        return result

    def finite_difference_jacobian(log_scale):
        center = np.asarray(log_scale, dtype=float)
        columns = []
        for index in range(len(center)):
            lower = center.copy()
            upper = center.copy()
            lower[index] = max(LOWER[index], center[index] - analysis.SENSITIVITY_STEP)
            upper[index] = min(UPPER[index], center[index] + analysis.SENSITIVITY_STEP)
            lower_residual = evaluate(lower)[0]
            upper_residual = evaluate(upper)[0]
            columns.append(
                (upper_residual - lower_residual) / (upper[index] - lower[index])
            )
        return np.column_stack(columns)

    result = least_squares(
        lambda x: evaluate(x)[0],
        np.zeros(3),
        jac=finite_difference_jacobian,
        bounds=(LOWER, UPPER),
        method="trf",
        x_scale="jac",
        max_nfev=10,
        ftol=1.0e-3,
        xtol=1.0e-3,
        gtol=1.0e-3,
        verbose=1,
    )

    candidates = {
        "Pre-fit b_n=2.45": np.zeros(3),
        "Identifiable combined fit": result.x,
    }
    detail_rows = []
    anode_rows = []
    candidate_runs = {}
    candidate_hppc = {}
    parameter_rows = []
    for name, vector in candidates.items():
        residual, dynamic_runs, hppc_model, record = evaluate(vector, keep_runs=True)
        candidate_runs[name] = dynamic_runs
        candidate_hppc[name] = hppc_model
        values = values_from_log_scale(vector)
        parameter_rows.append({"Scenario": name, **values})
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                sim = dynamic_runs[(rate, charge)]
                metrics = priority.curve_metrics(
                    experiment[(rate, charge)], sim, qcell, charge
                )
                potential = np.asarray(sim["potential_V"], dtype=float)
                minimum_index = int(np.nanargmin(potential))
                below_zero = analysis.anode_analysis.time_fraction(
                    sim["time_s"], potential < 0.0
                )
                detail_rows.append(
                    {
                        "Scenario": name,
                        "Trust": "combined dynamic+HPPC fit",
                        "Note": "Dsn,Dsp and common kn/kp scale; b_n=2.45",
                        "Cell": "6-6",
                        "C_rate": rate,
                        "Direction": direction,
                        **metrics,
                        "Capacity_error_mAh": 1000.0
                        * (metrics["Q_end_model_Ah"] - metrics["Q_end_exp_Ah"]),
                    }
                )
                anode_rows.append(
                    {
                        "Scenario": name,
                        "C_rate": rate,
                        "Direction": direction,
                        "Minimum_anode_potential_V": float(potential[minimum_index]),
                        "Below_0V_time_fraction_pct": 100.0 * below_zero,
                        "Surface_stoichiometry_at_minimum": float(
                            sim["surface_sto"][minimum_index]
                        ),
                    }
                )

    detail = pd.DataFrame(detail_rows)
    summary = p0.summarize(detail)
    anode_table = pd.DataFrame(anode_rows)
    hppc_rows = []
    for name, model_values in candidate_hppc.items():
        errors = 1000.0 * (model_values - hppc_exp)
        hppc_rows.append(
            {
                "Scenario": name,
                "HPPC_RMSE_mV": float(np.sqrt(np.mean(errors**2))),
                "HPPC_MAE_mV": float(np.mean(np.abs(errors))),
                "HPPC_Bias_mV": float(np.mean(errors)),
            }
        )
    hppc_summary = pd.DataFrame(hppc_rows)
    detail.to_csv(OUT / "dynamic_detail.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(OUT / "dynamic_summary.csv", index=False, encoding="utf-8-sig")
    anode_table.to_csv(OUT / "anode_potential_summary.csv", index=False, encoding="utf-8-sig")
    hppc_summary.to_csv(OUT / "hppc_summary.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(parameter_rows).to_csv(
        OUT / "fitted_parameters.csv", index=False, encoding="utf-8-sig"
    )

    colors = {"Pre-fit b_n=2.45": "#D55E00", "Identifiable combined fit": "#0072B2"}
    fig, axes = plt.subplots(2, 3, figsize=(15.6, 8.5), sharex=True, sharey=True, constrained_layout=True)
    for column, rate in enumerate(ga.RATES):
        for row, charge in enumerate((True, False)):
            direction = "Charge" if charge else "Discharge"
            ax = axes[row, column]
            exp = experiment[(rate, charge)]
            q_exp = priority.transferred_capacity(exp, qcell, charge)
            ax.plot(q_exp, exp["V"], color="#202020", lw=2.3, label="Experiment")
            for name in candidates:
                sim = candidate_runs[name][(rate, charge)]
                ax.plot(
                    sim["capacity_Ah"],
                    sim["V"],
                    color=colors[name],
                    lw=2.0,
                    linestyle="-." if name.startswith("Pre-fit") else "-",
                    label=name,
                )
            ax.set_title("{}C {}".format(rate, direction), weight="bold")
            ax.set_xlim(0.0, 2.42)
            ax.set_ylim(2.9, 4.3)
            ax.grid(alpha=0.25)
            if row == 1:
                ax.set_xlabel("Transferred capacity (Ah)")
            if column == 0:
                ax.set_ylabel("Cell voltage (V)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.suptitle("ORegan2022, $b_n=2.45$: before/after combined dynamic+HPPC fit", y=1.06, weight="bold")
    fig.savefig(OUT / "voltage_before_after.png", dpi=240, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(15.6, 8.5), sharex=True, sharey="row", constrained_layout=True)
    for column, rate in enumerate(ga.RATES):
        for row, charge in enumerate((True, False)):
            direction = "Charge" if charge else "Discharge"
            ax = axes[row, column]
            for name in candidates:
                sim = candidate_runs[name][(rate, charge)]
                ax.plot(
                    sim["capacity_Ah"],
                    sim["potential_V"],
                    color=colors[name],
                    lw=2.0,
                    linestyle="-." if name.startswith("Pre-fit") else "-",
                    label=name,
                )
            ax.axhline(0.0, color="#202020", lw=1.0, linestyle=":")
            ax.set_title("{}C {}".format(rate, direction), weight="bold")
            ax.set_xlim(0.0, 2.42)
            ax.set_ylim((-0.12, 0.35) if charge else (-0.03, 0.85))
            ax.grid(alpha=0.25)
            if row == 1:
                ax.set_xlabel("Transferred capacity (Ah)")
            if column == 0:
                ax.set_ylabel("Anode potential (V vs. Li/Li+)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.suptitle("Anode-potential prediction before/after combined fit", y=1.06, weight="bold")
    fig.savefig(OUT / "anode_potential_before_after.png", dpi=240, bbox_inches="tight")
    plt.close(fig)

    manifest = {
        "reason_for_common_kinetics_scale": "kn/kp normalized sensitivity cosine = 0.993; their ratio is not separately identifiable",
        "fixed_brugg_n": analysis.BRUGG_N,
        "fit_coordinates": FIT_NAMES,
        "bounds_scale": {
            "Dsn": [0.5, 2.0],
            "Dsp": [0.5, 2.0],
            "common_kn_kp": [0.3, 3.0],
        },
        "objective": "equal dataset-level least-squares weight for 480 dynamic center-overlap voltage points and 50 HPPC key-point voltages",
        "capacity_in_objective": False,
        "optimization_success": bool(result.success),
        "optimization_message": str(result.message),
        "nfev": int(result.nfev),
    }
    (OUT / "fit_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("\nPARAMETERS\n", pd.DataFrame(parameter_rows).to_string(index=False))
    print("\nDYNAMIC SUMMARY\n", summary.to_string(index=False))
    print("\nHPPC SUMMARY\n", hppc_summary.to_string(index=False))
    print("\nANODE\n", anode_table.to_string(index=False))


if __name__ == "__main__":
    main()
