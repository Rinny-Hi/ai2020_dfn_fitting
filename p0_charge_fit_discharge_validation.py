"""Five-parameter charge fit followed by untouched discharge validation.

Four optimizers see the same bounded log-parameter problem.  The objective is
terminal-voltage residual only, equally sampled for each July BoL cell and for
0.5C/1C/2C.  CC cutoff capacity is calculated only after fitting.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution, dual_annealing, least_squares, minimize

import c50_selected_ocp_dynamic_validation as selected
import comprehensive_prefit_cross_cohort as cross
import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import paper_geometry_area_porosity_ablation as area_study
import prefit_physical_anchor_decision as p0
import priority_initial_state_protocol_recheck as priority
import p0_five_parameter_sensitivity_area as sens


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260927_p0_charge_fit_discharge_validation"
RESULTS.mkdir(parents=True, exist_ok=True)
PARAMETERS = sens.PARAMETERS


def decode(log_values: np.ndarray) -> dict[str, float]:
    return {name: float(np.exp(value)) for name, value in zip(PARAMETERS, log_values)}


def silent_simulate(model, stage, anode, cathode, protocol, charge):
    runs = {}
    direction = "Charge" if charge else "Discharge"
    for rate in ga.RATES:
        row = protocol.loc[("Old July BoL", rate, direction)]
        z0 = priority.monotone_voltage_inverse(
            model["v_qocv"], model["soc"], float(row.Rest_end_V)
        )
        local_stage = priority.stage_at_soc(stage, z0, charge)
        sim_rate = float(row.Measured_current_A) / cross.NOMINAL_CAPACITY_AH
        runs[rate] = ehq.run_dfn(
            sim_rate, charge, local_stage, model, anode, cathode, False, False
        )
    return runs


def main():
    curves, qc, protocol = sens.load_july()
    model0, stage, anode, cathode, _, _, qcell = selected.selected_inputs()
    p0_base = p0.build_model(model0, p0.Scenario("P0 anchored reference"))
    geometry = area_study.apply_paper_overlap(p0_base, preserve_capacity=True)
    kp_meta = sens.kp_from_area_specific_rct(geometry)
    initial = {**sens.INITIAL_COMMON, "kp": kp_meta["kp"]}
    x0 = np.log([initial[name] for name in PARAMETERS])
    lower = np.log([sens.BOUNDS[name][0] for name in PARAMETERS])
    upper = np.log([sens.BOUNDS[name][1] for name in PARAMETERS])
    bounds_pairs = list(zip(lower, upper))

    base_model = sens.apply_parameters(geometry, initial)
    base_charge = silent_simulate(base_model, stage, anode, cathode, protocol, True)
    grids = sens.objective_grids(curves, qcell, base_charge)
    experimental = sens.voltage_vector(None, curves, grids, experiment=True)
    n_residual = experimental.size

    cache: dict[tuple[float, ...], tuple[np.ndarray, dict | None]] = {}
    history: list[dict] = []
    evaluation_count = 0

    def evaluate(log_values):
        nonlocal evaluation_count
        clipped = np.clip(np.asarray(log_values, float), lower, upper)
        key = tuple(np.round(clipped, 12))
        if key in cache:
            return cache[key]
        evaluation_count += 1
        values = decode(clipped)
        start = time.perf_counter()
        status = "ok"
        runs = None
        try:
            model = sens.apply_parameters(geometry, values)
            runs = silent_simulate(model, stage, anode, cathode, protocol, True)
            blocks = []
            for rate in ga.RATES:
                sim = runs[rate]
                sim_t = np.asarray(sim["t_min"], float)
                sim_v = np.asarray(sim["V"], float)
                for cell in curves:
                    grid = grids[(rate, cell)]
                    prediction = np.interp(np.minimum(grid, sim_t[-1]), sim_t, sim_v)
                    block = prediction - np.interp(
                        grid,
                        curves[cell][(rate, True)]["t_min"],
                        curves[cell][(rate, True)]["V"],
                    )
                    # Missing experimental-time coverage is a time-voltage
                    # feasibility penalty, not a capacity-error objective.
                    missing = np.maximum(grid - sim_t[-1], 0.0)
                    if np.any(missing > 0):
                        block = block + 0.50 * missing / max(grid[-1], 1.0e-9)
                        status = "early_cutoff_penalty"
                    blocks.append(block)
            residual = 1000.0 * np.concatenate(blocks)
        except Exception as exc:  # optimizer must receive a fixed-length vector
            residual = np.full(n_residual, 1000.0)
            status = f"failed:{type(exc).__name__}"
        elapsed = time.perf_counter() - start
        rmse = float(np.sqrt(np.mean(residual**2)))
        history.append(
            {
                "Evaluation": evaluation_count,
                "RMSE_mV": rmse,
                "Status": status,
                "Elapsed_s": elapsed,
                **values,
            }
        )
        if evaluation_count == 1 or evaluation_count % 10 == 0:
            print(f"evaluation {evaluation_count}: {rmse:.3f} mV ({status})", flush=True)
            pd.DataFrame(history).to_csv(
                RESULTS / "optimizer_evaluation_history.csv",
                index=False,
                encoding="utf-8-sig",
            )
        cache[key] = (residual, runs)
        return residual, runs

    def residual_function(log_values):
        return evaluate(log_values)[0]

    def scalar_function(log_values):
        residual = residual_function(log_values)
        return float(np.mean(residual**2))

    algorithm_rows = []
    fitted = {}

    start = time.perf_counter()
    result = least_squares(
        residual_function,
        x0,
        bounds=(lower, upper),
        method="trf",
        x_scale="jac",
        max_nfev=28,
        ftol=2.0e-4,
        xtol=2.0e-4,
        gtol=2.0e-4,
    )
    fitted["TRF least_squares"] = result.x
    algorithm_rows.append(
        {
            "Algorithm": "TRF least_squares",
            "Success": bool(result.success),
            "Message": str(result.message),
            "NFEV_reported": int(result.nfev),
            "Wall_s": time.perf_counter() - start,
        }
    )

    start = time.perf_counter()
    result = minimize(
        scalar_function,
        x0,
        method="L-BFGS-B",
        bounds=bounds_pairs,
        options={"maxiter": 12, "maxfun": 75, "ftol": 1.0e-7},
    )
    fitted["L-BFGS-B"] = result.x
    algorithm_rows.append(
        {
            "Algorithm": "L-BFGS-B",
            "Success": bool(result.success),
            "Message": str(result.message),
            "NFEV_reported": int(result.nfev),
            "Wall_s": time.perf_counter() - start,
        }
    )

    start = time.perf_counter()
    result = differential_evolution(
        scalar_function,
        bounds_pairs,
        maxiter=2,
        popsize=3,
        seed=260927,
        polish=False,
        updating="immediate",
        workers=1,
        tol=0.01,
    )
    fitted["Differential evolution"] = result.x
    algorithm_rows.append(
        {
            "Algorithm": "Differential evolution",
            "Success": bool(result.success),
            "Message": str(result.message),
            "NFEV_reported": int(result.nfev),
            "Wall_s": time.perf_counter() - start,
        }
    )

    start = time.perf_counter()
    result = dual_annealing(
        scalar_function,
        bounds_pairs,
        maxfun=40,
        seed=260927,
        no_local_search=True,
        x0=x0,
    )
    fitted["Dual annealing"] = result.x
    algorithm_rows.append(
        {
            "Algorithm": "Dual annealing",
            "Success": bool(result.success),
            "Message": str(result.message),
            "NFEV_reported": int(result.nfev),
            "Wall_s": time.perf_counter() - start,
        }
    )

    # Add the pre-fit point for a transparent before/after comparison.
    candidates = {"Pre-fit initial": x0, **fitted}
    algorithm_metrics = []
    all_detail = []
    all_runs = {}
    for name, vector in candidates.items():
        values = decode(vector)
        model = sens.apply_parameters(geometry, values)
        charge_runs = silent_simulate(model, stage, anode, cathode, protocol, True)
        discharge_runs = silent_simulate(model, stage, anode, cathode, protocol, False)
        all_runs[name] = {True: charge_runs, False: discharge_runs}
        objective_residual = evaluate(vector)[0]
        algorithm_metrics.append(
            {
                "Algorithm": name,
                "Charge_objective_RMSE_mV": float(
                    np.sqrt(np.mean(objective_residual**2))
                ),
                **values,
            }
        )
        detail = sens.summarize_runs(
            name, {True: charge_runs, False: discharge_runs}, curves, qcell
        )
        detail = detail.rename(columns={"Area_policy": "Algorithm"})
        all_detail.append(detail)

    metrics = pd.DataFrame(algorithm_metrics).sort_values("Charge_objective_RMSE_mV")
    detail = pd.concat(all_detail, ignore_index=True)
    post = (
        detail.groupby(["Algorithm", "Direction"], as_index=False)
        .agg(
            Full_RMSE_mV=("Full_RMSE_mV", "mean"),
            Center10_70_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
            Capacity_RMSE_pct=(
                "Capacity_error_pct",
                lambda x: float(np.sqrt(np.mean(np.asarray(x, float) ** 2))),
            ),
            Mean_capacity_error_mAh=("Capacity_error_mAh", "mean"),
        )
        .merge(metrics, on="Algorithm", how="left")
        .sort_values(["Direction", "Center10_70_RMSE_mV"])
    )
    algorithms = pd.DataFrame(algorithm_rows).merge(metrics, on="Algorithm", how="left")
    algorithms = algorithms.sort_values("Charge_objective_RMSE_mV")
    algorithms.to_csv(RESULTS / "algorithm_comparison.csv", index=False, encoding="utf-8-sig")
    detail.to_csv(RESULTS / "postfit_detail_by_cell.csv", index=False, encoding="utf-8-sig")
    post.to_csv(RESULTS / "postfit_summary.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(history).to_csv(
        RESULTS / "optimizer_evaluation_history.csv", index=False, encoding="utf-8-sig"
    )

    best_name = str(metrics.iloc[0].Algorithm)
    colors = {
        "Pre-fit initial": "#7f7f7f",
        "TRF least_squares": "#0072B2",
        "L-BFGS-B": "#009E73",
        "Differential evolution": "#D55E00",
        "Dual annealing": "#CC79A7",
    }
    for charge, filename in (
        (True, "charge_fitting_time_voltage.png"),
        (False, "discharge_validation_time_voltage.png"),
    ):
        direction = "Charge" if charge else "Discharge"
        fig, axes = plt.subplots(1, 3, figsize=(16.2, 5.2), constrained_layout=True)
        for ax, rate in zip(axes, ga.RATES):
            for idx, (cell, branches) in enumerate(curves.items()):
                obs = branches[(rate, charge)]
                ax.plot(
                    obs["t_min"],
                    obs["V"],
                    color="black",
                    alpha=0.28,
                    lw=1.0,
                    label="Experiment (July BoL, n=3)" if idx == 0 else None,
                )
            for name in candidates:
                sim = all_runs[name][charge][rate]
                style = "--" if name == "Pre-fit initial" else "-"
                width = 2.5 if name == best_name else 1.35
                ax.plot(
                    sim["t_min"],
                    sim["V"],
                    color=colors[name],
                    ls=style,
                    lw=width,
                    alpha=0.92,
                    label=name,
                )
            ax.set(
                title=f"{rate:g}C {direction}",
                xlabel="Time [min]",
                ylabel="Voltage [V]",
            )
            ax.grid(alpha=0.2)
            ax.legend(fontsize=7.2)
        fig.suptitle(
            "Charge voltage fitting" if charge else "Discharge validation (not fitted)"
        )
        fig.savefig(RESULTS / filename, dpi=220)
        plt.close(fig)

    identifiability = []
    fit_only = metrics[metrics.Algorithm != "Pre-fit initial"]
    for parameter in PARAMETERS:
        values = fit_only[parameter].to_numpy(float)
        identifiability.append(
            {
                "Parameter": parameter,
                "Min_across_optimizers": float(np.min(values)),
                "Max_across_optimizers": float(np.max(values)),
                "Max_over_min": float(np.max(values) / np.min(values)),
                "Initial": initial[parameter],
                "Lower_bound": sens.BOUNDS[parameter][0],
                "Upper_bound": sens.BOUNDS[parameter][1],
            }
        )
    pd.DataFrame(identifiability).to_csv(
        RESULTS / "cross_optimizer_parameter_spread.csv",
        index=False,
        encoding="utf-8-sig",
    )

    manifest = {
        "cohort": "July BoL cells 6-5/6-6/6-8",
        "geometry": "teardown cathode overlap area, Qn/Qp preserved",
        "objective": (
            "charge terminal-voltage residual only on fixed experimental-time grids; "
            "0.5C/1C/2C and all cells equally sampled; capacity excluded"
        ),
        "objective_range": "10-55% of q_meas, truncated to 95% of pre-fit cutoff time",
        "postfit_metrics": "full and 10-70% voltage RMSE plus CC cutoff capacity error",
        "validation": "discharge curves were never shown to the charge optimizers",
        "initial_parameters": initial,
        "bounds": sens.BOUNDS,
        "kp_recalculation": kp_meta,
        "best_charge_objective_algorithm": best_name,
        "warning": (
            "kn and kp are nearly collinear in local charge sensitivity; optimizer "
            "agreement and validation, not the smallest training RMSE alone, must decide."
        ),
    }
    (RESULTS / "analysis_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("\nALGORITHMS\n", algorithms.to_string(index=False))
    print("\nPOST-FIT\n", post.to_string(index=False))
    print(f"\nBest charge objective: {best_name}")


if __name__ == "__main__":
    main()
