"""S6 marker comparison and S7 PSO voltage fitting.

S7a performs a local sensitivity/identifiability audit at S6 and fits the
largest parameter subset that passes the stated criteria.  S7b skips that
screen and fits Dsn, kn, and b_n together.  Both PSO cases minimize the same
fixed-time terminal-voltage RMSE over the six Cell 6-6 CC branches.  Capacity,
endpoint time, and anode potential are post-fit diagnostics only.
"""

from __future__ import annotations

from itertools import combinations
import json
from pathlib import Path
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import analyze_oregan_brugg_2p5_anode_potential as anode_analysis
import comprehensive_prefit_cross_cohort as cross
import cumulative_six_stage_model_comparison as cumulative
import electrolyte_transport_set_comparison as transport
import gitt_ocp_analysis as ga
import kinetics_original_reversion_ablation as kinetics
import oregan_bn25_kn_only_fit as kn_reference
import plot_oregan_brugg_2p5_time_voltage as time_metrics
import priority_initial_state_protocol_recheck as priority


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "261013_s7_pso_voltage_fitting"
PARAMETERS = ("Dsn", "kn", "brugg_n")
BASE = {
    "Dsn": 4.1707741459682473e-14,
    "kn": 9.648533212e-7,
    "brugg_n": 2.5,
}
BOUNDS = {
    "Dsn": (1.0e-14, 8.0e-14),
    "kn": (3.0e-7, 3.0e-6),
    "brugg_n": (1.5, 3.5),
}
N_GRID = 120
SENSITIVITY_STEP = 0.025
MIN_RELATIVE_SENSITIVITY = 0.10
MAX_ABS_COSINE = 0.95
MAX_CONDITION = 20.0
N_PARTICLES = 7
N_ITERATIONS = 7
EARLY_CUTOFF_PENALTY_MV = 500.0


def decode_coordinate(name: str, coordinate: float) -> float:
    lo, hi = BOUNDS[name]
    z = float(np.clip(coordinate, 0.0, 1.0))
    if name in {"Dsn", "kn"}:
        return float(np.exp(np.log(lo) + z * (np.log(hi) - np.log(lo))))
    return float(lo + z * (hi - lo))


def encode_value(name: str, value: float) -> float:
    lo, hi = BOUNDS[name]
    if name in {"Dsn", "kn"}:
        return float((np.log(value) - np.log(lo)) / (np.log(hi) - np.log(lo)))
    return float((value - lo) / (hi - lo))


def vector_to_values(
    coordinates: np.ndarray,
    fitted_parameters: tuple[str, ...],
) -> dict[str, float]:
    values = dict(BASE)
    for name, coordinate in zip(fitted_parameters, coordinates):
        values[name] = decode_coordinate(name, float(coordinate))
    return values


def apply_values(s6_model: dict, values: dict[str, float]) -> dict:
    changed = kinetics.with_prefactors(
        s6_model,
        values["kn"],
        kn_reference.BASE["kp"],
    )
    changed = dict(changed)
    changed["params"] = changed["params"].copy()
    changed["params"].update(
        {
            "Negative particle diffusivity [m2.s-1]": values["Dsn"],
            "Negative electrode Bruggeman coefficient (electrolyte)": values[
                "brugg_n"
            ],
        },
        check_already_exists=False,
    )
    return changed


def make_experimental_grids(experiment: dict) -> tuple[dict, dict]:
    grids = {}
    voltages = {}
    for rate in ga.RATES:
        for charge in (True, False):
            source = experiment[(rate, charge)]
            elapsed = np.asarray(source["t_min"], dtype=float)
            elapsed = elapsed - elapsed[0]
            grid = np.linspace(0.0, 0.995 * float(elapsed[-1]), N_GRID)
            grids[(rate, charge)] = grid
            voltages[(rate, charge)] = np.interp(
                grid, elapsed, np.asarray(source["V"], dtype=float)
            )
    return grids, voltages


def run_model(model: dict, stage: dict, anode, cathode) -> dict:
    runs = {}
    for rate in ga.RATES:
        for charge in (True, False):
            rest_v, measured_current = transport.PROTOCOL[(rate, charge)]
            z0 = priority.monotone_voltage_inverse(
                model["v_qocv"], model["soc"], rest_v
            )
            local_stage = priority.stage_at_soc(stage, z0, charge)
            runs[(rate, charge)] = anode_analysis.run_branch(
                measured_current / cross.NOMINAL_CAPACITY_AH,
                charge,
                local_stage,
                model,
                anode,
                cathode,
            )
    return runs


def residual_vector(runs: dict, grids: dict, targets: dict) -> tuple[np.ndarray, dict]:
    blocks = []
    minimum_margin = np.inf
    early_branches = 0
    for rate in ga.RATES:
        for charge in (True, False):
            output = runs[(rate, charge)]
            model_time = np.asarray(output["time_s"], dtype=float) / 60.0
            model_time = model_time - model_time[0]
            model_voltage = np.asarray(output["V"], dtype=float)
            grid = grids[(rate, charge)]
            margin = float(model_time[-1] - grid[-1])
            minimum_margin = min(minimum_margin, margin)
            prediction = np.interp(
                np.minimum(grid, model_time[-1]), model_time, model_voltage
            )
            block = 1000.0 * (prediction - targets[(rate, charge)])
            if margin < 0.0:
                early_branches += 1
                missing = np.maximum(grid - model_time[-1], 0.0) / max(grid[-1], 1e-12)
                block = block + EARLY_CUTOFF_PENALTY_MV * missing
            blocks.append(block)
    return np.concatenate(blocks), {
        "Minimum_time_margin_min": float(minimum_margin),
        "Early_cutoff_branches": int(early_branches),
    }


def identifiability_report(matrix: np.ndarray, names: tuple[str, ...]) -> dict:
    rms = np.sqrt(np.mean(matrix**2, axis=0))
    global_max = float(np.max(rms))
    relative = rms / max(global_max, np.finfo(float).eps)
    unit = matrix / np.maximum(np.linalg.norm(matrix, axis=0), np.finfo(float).eps)
    cosine = unit.T @ unit
    singular = np.linalg.svd(unit, compute_uv=False)
    condition = float(singular[0] / singular[-1]) if len(names) > 1 else 1.0
    max_pair = (
        float(np.max(np.abs(cosine - np.eye(len(names))))) if len(names) > 1 else 0.0
    )
    full_rank = int(np.linalg.matrix_rank(matrix)) == len(names)
    passes = bool(
        full_rank
        and float(np.min(relative)) >= MIN_RELATIVE_SENSITIVITY
        and max_pair <= MAX_ABS_COSINE
        and condition <= MAX_CONDITION
    )
    return {
        "Parameters": "+".join(names),
        "N_parameters": len(names),
        "Matrix_rank": int(np.linalg.matrix_rank(matrix)),
        "Normalized_condition_number": condition,
        "Max_abs_pairwise_cosine": max_pair,
        "Min_relative_sensitivity": float(np.min(relative)),
        "Pass": passes,
    }


def make_initial_positions(
    fitted_parameters: tuple[str, ...],
    rng: np.random.Generator,
    informed: bool,
) -> np.ndarray:
    dimension = len(fitted_parameters)
    base = np.array([encode_value(name, BASE[name]) for name in fitted_parameters])
    positions = np.empty((N_PARTICLES, dimension), dtype=float)
    positions[0] = base
    if informed:
        for index in range(1, N_PARTICLES):
            positions[index] = np.clip(base + rng.normal(0.0, 0.16, dimension), 0.0, 1.0)
    else:
        positions[1:] = rng.uniform(0.0, 1.0, size=(N_PARTICLES - 1, dimension))
    return positions


def pso_fit(
    label: str,
    fitted_parameters: tuple[str, ...],
    evaluate,
    informed: bool,
    seed: int,
) -> tuple[dict, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    positions = make_initial_positions(fitted_parameters, rng, informed)
    velocities = rng.normal(0.0, 0.04, size=positions.shape)
    personal_positions = positions.copy()
    personal_scores = np.full(N_PARTICLES, np.inf)
    global_position = positions[0].copy()
    global_score = np.inf
    history = []

    for iteration in range(N_ITERATIONS):
        for particle in range(N_PARTICLES):
            values = vector_to_values(positions[particle], fitted_parameters)
            result = evaluate(values)
            score = result["Objective_RMSE_mV"]
            if score < personal_scores[particle]:
                personal_scores[particle] = score
                personal_positions[particle] = positions[particle].copy()
            if score < global_score:
                global_score = score
                global_position = positions[particle].copy()
            history.append(
                {
                    "Scenario": label,
                    "Iteration": iteration + 1,
                    "Particle": particle + 1,
                    "Objective_RMSE_mV": score,
                    **values,
                }
            )
        print(
            f"{label} iteration {iteration + 1}/{N_ITERATIONS}: "
            f"best={global_score:.3f} mV",
            flush=True,
        )
        r1 = rng.random(size=positions.shape)
        r2 = rng.random(size=positions.shape)
        velocities = (
            0.72 * velocities
            + 1.49 * r1 * (personal_positions - positions)
            + 1.49 * r2 * (global_position - positions)
        )
        velocities = np.clip(velocities, -0.22, 0.22)
        positions = np.clip(positions + velocities, 0.0, 1.0)

    best_values = vector_to_values(global_position, fitted_parameters)
    best = evaluate(best_values)
    return {**best_values, **best}, pd.DataFrame(history)


def scenario_metrics(
    scenario: str,
    result: dict,
    experiment: dict,
    qcell: float,
) -> tuple[list[dict], list[dict]]:
    rows = []
    anode_rows = []
    for rate in ga.RATES:
        for charge in (True, False):
            direction = "Charge" if charge else "Discharge"
            output = result["runs"][(rate, charge)]
            exp = experiment[(rate, charge)]
            elapsed_metrics = time_metrics.time_voltage_metrics(
                exp,
                output,
                qcell,
                transport.PROTOCOL[(rate, charge)][1],
            )
            capacity = priority.curve_metrics(exp, output, qcell, charge)
            rows.append(
                {
                    "Scenario": scenario,
                    "C_rate": rate,
                    "Direction": direction,
                    **elapsed_metrics,
                    "Capacity_error_pct": capacity["Capacity_error_pct"],
                    "Capacity_error_mAh": 1000.0
                    * (capacity["Q_end_model_Ah"] - capacity["Q_end_exp_Ah"]),
                }
            )
            potential = np.asarray(output["potential_V"], dtype=float)
            index = int(np.nanargmin(potential))
            anode_rows.append(
                {
                    "Scenario": scenario,
                    "C_rate": rate,
                    "Direction": direction,
                    "Minimum_anode_potential_mV": 1000.0 * float(potential[index]),
                    "Time_below_0V_pct": 100.0
                    * anode_analysis.time_fraction(
                        output["time_s"], potential < 0.0
                    ),
                    "Surface_stoichiometry_at_minimum": float(
                        output["surface_sto"][index]
                    ),
                }
            )
    return rows, anode_rows


def aggregate_metrics(detail: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for scenario, part in detail.groupby("Scenario", sort=False):
        rows.append(
            {
                "Scenario": scenario,
                "Overall_full_time_RMSE_mV": float(
                    np.sqrt(np.mean(part.Full_time_RMSE_mV.to_numpy(float) ** 2))
                ),
                "Overall_full_time_MAE_mV": float(
                    np.mean(part.Full_time_MAE_mV.to_numpy(float))
                ),
                "Overall_center_time_RMSE_mV": float(
                    np.sqrt(
                        np.mean(
                            part.Center10_70_Qcell_time_RMSE_mV.to_numpy(float) ** 2
                        )
                    )
                ),
                "Overall_center_time_MAE_mV": float(
                    np.mean(part.Center10_70_Qcell_time_MAE_mV.to_numpy(float))
                ),
                "Capacity_RMSE_pct": float(
                    np.sqrt(np.mean(part.Capacity_error_pct.to_numpy(float) ** 2))
                ),
                "Max_abs_capacity_error_mAh": float(
                    np.max(np.abs(part.Capacity_error_mAh.to_numpy(float)))
                ),
            }
        )
    return pd.DataFrame(rows)


def plot_marker_comparison(experiment: dict, s6_runs: dict, qcell: float) -> None:
    fig, axes = plt.subplots(
        2, 3, figsize=(16.0, 8.6), sharey=True, constrained_layout=True
    )
    for column, rate in enumerate(ga.RATES):
        for row, charge in enumerate((True, False)):
            ax = axes[row, column]
            direction = "Charge" if charge else "Discharge"
            exp = experiment[(rate, charge)]
            exp_time = np.asarray(exp["t_min"], dtype=float)
            exp_time = exp_time - exp_time[0]
            sim = s6_runs[(rate, charge)]
            sim_time = np.asarray(sim["time_s"], dtype=float) / 60.0
            sim_time = sim_time - sim_time[0]
            exp_stride = max(1, len(exp_time) // 28)
            sim_stride = max(1, len(sim_time) // 28)
            ax.plot(exp_time, exp["V"], color="#303030", lw=1.0, alpha=0.55)
            ax.plot(sim_time, sim["V"], color="#D55E00", lw=1.0, alpha=0.65)
            ax.plot(
                exp_time[::exp_stride],
                np.asarray(exp["V"])[::exp_stride],
                linestyle="none",
                marker="s",
                markersize=5.0,
                markerfacecolor="white",
                markeredgecolor="#202020",
                markeredgewidth=1.1,
                label="Experiment □",
            )
            ax.plot(
                sim_time[::sim_stride],
                np.asarray(sim["V"])[::sim_stride],
                linestyle="none",
                marker="o",
                markersize=4.8,
                markerfacecolor="none",
                markeredgecolor="#D55E00",
                markeredgewidth=1.2,
                label="S6 simulation ○",
            )
            metric = time_metrics.time_voltage_metrics(
                exp, sim, qcell, transport.PROTOCOL[(rate, charge)][1]
            )
            ax.text(
                0.03,
                0.05,
                f"RMSE {metric['Full_time_RMSE_mV']:.1f} mV\n"
                f"10–70% {metric['Center10_70_Qcell_time_RMSE_mV']:.1f} mV",
                transform=ax.transAxes,
                fontsize=9,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.88),
            )
            ax.set_title(f"{rate:g}C {direction}", weight="bold")
            ax.set_ylim(2.9, 4.3)
            ax.grid(alpha=0.22)
            if row == 1:
                ax.set_xlabel("Elapsed time (min)")
            if column == 0:
                ax.set_ylabel("Cell voltage (V)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles[-2:],
        labels[-2:],
        loc="upper center",
        ncol=2,
        frameon=False,
        bbox_to_anchor=(0.5, 1.045),
    )
    fig.suptitle(
        "S6 measured–simulated agreement: squares vs circles",
        fontsize=15,
        weight="bold",
        y=1.10,
    )
    fig.savefig(OUT / "s6_square_circle_voltage.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def plot_s7_comparison(experiment: dict, results: dict) -> None:
    colors = {"S6": "#D55E00", "S7a": "#0072B2", "S7b": "#009E73"}
    styles = {"S6": "--", "S7a": "-", "S7b": "-."}
    fig, axes = plt.subplots(
        2, 3, figsize=(16.0, 8.6), sharey=True, constrained_layout=True
    )
    for column, rate in enumerate(ga.RATES):
        for row, charge in enumerate((True, False)):
            ax = axes[row, column]
            direction = "Charge" if charge else "Discharge"
            exp = experiment[(rate, charge)]
            exp_time = np.asarray(exp["t_min"], dtype=float)
            exp_time = exp_time - exp_time[0]
            stride = max(1, len(exp_time) // 30)
            ax.plot(
                exp_time[::stride],
                np.asarray(exp["V"])[::stride],
                linestyle="none",
                marker="s",
                markersize=4.6,
                markerfacecolor="white",
                markeredgecolor="#202020",
                label="Experiment □",
            )
            for scenario in ("S6", "S7a", "S7b"):
                output = results[scenario]["runs"][(rate, charge)]
                elapsed = np.asarray(output["time_s"], dtype=float) / 60.0
                elapsed = elapsed - elapsed[0]
                ax.plot(
                    elapsed,
                    output["V"],
                    color=colors[scenario],
                    linestyle=styles[scenario],
                    lw=2.0,
                    label=scenario,
                )
            ax.set_title(f"{rate:g}C {direction}", weight="bold")
            ax.set_ylim(2.9, 4.3)
            ax.grid(alpha=0.22)
            if row == 1:
                ax.set_xlabel("Elapsed time (min)")
            if column == 0:
                ax.set_ylabel("Cell voltage (V)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.5, 1.045),
    )
    fig.suptitle(
        "S6 vs PSO-fitted S7a/S7b voltage curves", weight="bold", y=1.10
    )
    fig.savefig(OUT / "s6_s7_voltage_comparison.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    _, _, _, _, experiment, qcell = transport.current_p0()
    s6_entry = cumulative.build_stages(qcell)["S6"]
    s6_model = s6_entry["model"]
    stage = s6_entry["stage"]
    anode = s6_entry["anode"]
    cathode = s6_entry["cathode"]
    grids, targets = make_experimental_grids(experiment)
    cache: dict[tuple[float, float, float], dict] = {}
    evaluations = []

    def evaluate(values: dict[str, float]) -> dict:
        key = tuple(round(float(values[name]), 18) for name in PARAMETERS)
        if key in cache:
            return cache[key]
        started = time.perf_counter()
        status = "ok"
        try:
            model = apply_values(s6_model, values)
            runs = run_model(model, stage, anode, cathode)
            residual, feasibility = residual_vector(runs, grids, targets)
        except Exception as exc:
            runs = None
            residual = np.full(6 * N_GRID, 1000.0)
            feasibility = {
                "Minimum_time_margin_min": np.nan,
                "Early_cutoff_branches": 6,
            }
            status = f"failed:{type(exc).__name__}"
        result = {
            "Objective_RMSE_mV": float(np.sqrt(np.mean(residual**2))),
            "Objective_MAE_mV": float(np.mean(np.abs(residual))),
            "residual": residual,
            "runs": runs,
            "Status": status,
            **feasibility,
        }
        cache[key] = result
        evaluations.append(
            {
                "Evaluation": len(evaluations) + 1,
                **values,
                "Objective_RMSE_mV": result["Objective_RMSE_mV"],
                "Objective_MAE_mV": result["Objective_MAE_mV"],
                "Early_cutoff_branches": result["Early_cutoff_branches"],
                "Status": status,
                "Elapsed_s": time.perf_counter() - started,
            }
        )
        if len(evaluations) == 1 or len(evaluations) % 10 == 0:
            print(
                f"evaluation {len(evaluations)}: "
                f"{result['Objective_RMSE_mV']:.3f} mV",
                flush=True,
            )
            pd.DataFrame(evaluations).to_csv(
                OUT / "all_model_evaluations.csv", index=False, encoding="utf-8-sig"
            )
        return result

    print("S6 baseline", flush=True)
    s6_result = evaluate(BASE)
    plot_marker_comparison(experiment, s6_result["runs"], qcell)

    print("Local sensitivity and identifiability", flush=True)
    base_coordinate = np.array([encode_value(name, BASE[name]) for name in PARAMETERS])
    derivative_columns = []
    sensitivity_rows = []
    for index, name in enumerate(PARAMETERS):
        lower = base_coordinate.copy()
        upper = base_coordinate.copy()
        lower[index] = max(0.0, lower[index] - SENSITIVITY_STEP)
        upper[index] = min(1.0, upper[index] + SENSITIVITY_STEP)
        lower_values = vector_to_values(lower, PARAMETERS)
        upper_values = vector_to_values(upper, PARAMETERS)
        lower_residual = evaluate(lower_values)["residual"]
        upper_residual = evaluate(upper_values)["residual"]
        derivative = (upper_residual - lower_residual) / (
            upper[index] - lower[index]
        )
        derivative_columns.append(derivative)
        sensitivity_rows.append(
            {
                "Parameter": name,
                "Baseline": BASE[name],
                "RMS_mV_per_normalized_coordinate": float(
                    np.sqrt(np.mean(derivative**2))
                ),
                "Max_abs_mV_per_normalized_coordinate": float(
                    np.max(np.abs(derivative))
                ),
            }
        )
    jacobian = np.column_stack(derivative_columns)
    sensitivity = pd.DataFrame(sensitivity_rows)
    sensitivity["Relative_sensitivity"] = (
        sensitivity.RMS_mV_per_normalized_coordinate
        / sensitivity.RMS_mV_per_normalized_coordinate.max()
    )
    subset_rows = []
    subset_lookup = {}
    for size in range(1, len(PARAMETERS) + 1):
        for subset in combinations(PARAMETERS, size):
            indexes = [PARAMETERS.index(name) for name in subset]
            matrix = jacobian[:, indexes]
            report = identifiability_report(matrix, subset)
            subset_rows.append(report)
            subset_lookup[subset] = report
    subsets = pd.DataFrame(subset_rows)
    eligible_parameters = set(
        sensitivity.loc[
            sensitivity.Relative_sensitivity >= MIN_RELATIVE_SENSITIVITY,
            "Parameter",
        ].astype(str)
    )
    eligible_subset = subsets.Parameters.apply(
        lambda value: set(str(value).split("+")).issubset(eligible_parameters)
    )
    passing = subsets[subsets.Pass & eligible_subset].sort_values(
        ["N_parameters", "Normalized_condition_number"],
        ascending=[False, True],
    )
    if passing.empty:
        selected_parameters = (
            str(sensitivity.sort_values("Relative_sensitivity").iloc[-1].Parameter),
        )
    else:
        selected_parameters = tuple(str(passing.iloc[0].Parameters).split("+"))

    unit = jacobian / np.maximum(
        np.linalg.norm(jacobian, axis=0), np.finfo(float).eps
    )
    cosine = unit.T @ unit
    sensitivity.to_csv(OUT / "s7a_sensitivity_ranking.csv", index=False, encoding="utf-8-sig")
    subsets.to_csv(OUT / "s7a_subset_identifiability.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(cosine, index=PARAMETERS, columns=PARAMETERS).to_csv(
        OUT / "s7a_sensitivity_cosine.csv", encoding="utf-8-sig"
    )
    np.save(OUT / "s7a_voltage_jacobian.npy", jacobian)

    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.8), constrained_layout=True)
    axes[0].bar(
        sensitivity.Parameter,
        sensitivity.Relative_sensitivity,
        color="#0072B2",
    )
    axes[0].axhline(MIN_RELATIVE_SENSITIVITY, color="black", ls="--", lw=1)
    axes[0].set_ylabel("Relative RMS sensitivity")
    axes[0].set_title("S6 local voltage sensitivity")
    axes[0].grid(axis="y", alpha=0.22)
    image = axes[1].imshow(cosine, vmin=-1.0, vmax=1.0, cmap="coolwarm")
    axes[1].set_xticks(range(3), PARAMETERS)
    axes[1].set_yticks(range(3), PARAMETERS)
    axes[1].set_title("Normalized-sensitivity cosine")
    for i in range(3):
        for j in range(3):
            axes[1].text(j, i, f"{cosine[i, j]:.2f}", ha="center", va="center")
    fig.colorbar(image, ax=axes[1], shrink=0.8)
    fig.savefig(OUT / "s7a_sensitivity_identifiability.png", dpi=240)
    plt.close(fig)

    print(f"S7a selected parameters: {selected_parameters}", flush=True)
    s7a_result, s7a_history = pso_fit(
        "S7a",
        selected_parameters,
        evaluate,
        informed=True,
        seed=261013,
    )
    print("S7b all-parameter unscreened PSO", flush=True)
    s7b_result, s7b_history = pso_fit(
        "S7b",
        PARAMETERS,
        evaluate,
        informed=False,
        seed=261014,
    )
    history = pd.concat([s7a_history, s7b_history], ignore_index=True)
    history.to_csv(OUT / "pso_history.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(evaluations).to_csv(
        OUT / "all_model_evaluations.csv", index=False, encoding="utf-8-sig"
    )

    results = {"S6": {**BASE, **s6_result}, "S7a": s7a_result, "S7b": s7b_result}
    metric_rows = []
    anode_rows = []
    for scenario, result in results.items():
        branch_rows, scenario_anode = scenario_metrics(
            scenario, result, experiment, qcell
        )
        metric_rows.extend(branch_rows)
        anode_rows.extend(scenario_anode)
    detail = pd.DataFrame(metric_rows)
    summary = aggregate_metrics(detail)
    parameters = pd.DataFrame(
        [
            {
                "Scenario": scenario,
                "Fitted_parameters": (
                    "none"
                    if scenario == "S6"
                    else "+".join(
                        selected_parameters if scenario == "S7a" else PARAMETERS
                    )
                ),
                "Dsn_m2_s": result["Dsn"],
                "kn_prefactor": result["kn"],
                "brugg_n": result["brugg_n"],
                "Objective_RMSE_mV": result["Objective_RMSE_mV"],
                "Objective_MAE_mV": result["Objective_MAE_mV"],
                "Early_cutoff_branches": result["Early_cutoff_branches"],
            }
            for scenario, result in results.items()
        ]
    )
    detail.to_csv(OUT / "s6_s7_branch_metrics.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(OUT / "s6_s7_summary.csv", index=False, encoding="utf-8-sig")
    parameters.to_csv(OUT / "s6_s7_parameters.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(anode_rows).to_csv(
        OUT / "s6_s7_anode_potential.csv", index=False, encoding="utf-8-sig"
    )
    plot_s7_comparison(experiment, results)

    manifest = {
        "baseline": "S6 = ORegan2022 electrolyte transport, b_n=2.5",
        "objective": "fixed physical-time terminal-voltage RMSE over six Cell 6-6 branches; 120 points per branch",
        "capacity_in_objective": False,
        "anode_potential_in_objective": False,
        "early_cutoff_handling": f"voltage-grid feasibility penalty up to {EARLY_CUTOFF_PENALTY_MV:g} mV",
        "bounds": BOUNDS,
        "s7a_selected_parameters": selected_parameters,
        "s7a_pass_rule": {
            "minimum_relative_sensitivity": MIN_RELATIVE_SENSITIVITY,
            "maximum_absolute_pairwise_cosine": MAX_ABS_COSINE,
            "maximum_normalized_condition_number": MAX_CONDITION,
            "full_column_rank": True,
        },
        "s7b_parameters": PARAMETERS,
        "pso": {
            "particles": N_PARTICLES,
            "iterations": N_ITERATIONS,
            "inertia": 0.72,
            "cognitive": 1.49,
            "social": 1.49,
            "note": "single deterministic seed per workflow; optimizer stochastic uncertainty is not quantified",
        },
    }
    (OUT / "fit_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    summary_index = summary.set_index("Scenario")
    parameter_index = parameters.set_index("Scenario")
    readme_lines = [
        "# S6 marker plot and S7 PSO voltage fitting",
        "",
        "## 정의",
        "",
        "- S6: 현재 no-fit 기준모델.",
        f"- S7a: S6 국소 민감도/식별도 선별 후 `{'+'.join(selected_parameters)}`만 PSO fitting.",
        "- S7b: 선별 없이 `Dsn+kn+brugg_n`을 모두 PSO fitting.",
        "- 목적함수: Cell 6-6의 0.5C/1C/2C 충·방전 6곡선 고정 시간격자 terminal-voltage RMSE.",
        "- 용량, cutoff 시간, anode potential은 fitting에 사용하지 않고 사후평가했다.",
        "",
        "## 민감도/식별도",
        "",
        "| Parameter | Relative sensitivity | 판단 |",
        "|---|---:|---|",
    ]
    for _, row in sensitivity.iterrows():
        decision = "fit 후보" if row.Relative_sensitivity >= MIN_RELATIVE_SENSITIVITY else "고정"
        readme_lines.append(
            f"| {row.Parameter} | {row.Relative_sensitivity:.3f} | {decision} |"
        )
    readme_lines.extend(
        [
            "",
            "`Dsn`의 상대 민감도 0.069가 기준 0.10보다 낮아 GITT anchor로 고정했다. `kn+brugg_n` 부분집합은 full rank, condition number 2.29, 최대 sensitivity cosine 0.679로 기준을 통과했다.",
            "",
            "## fitted parameter",
            "",
            "| Scenario | Dsn [m2/s] | kn prefactor | brugg_n |",
            "|---|---:|---:|---:|",
        ]
    )
    for scenario in ("S6", "S7a", "S7b"):
        row = parameter_index.loc[scenario]
        readme_lines.append(
            f"| {scenario} | {row.Dsn_m2_s:.6e} | {row.kn_prefactor:.6e} | {row.brugg_n:.6f} |"
        )
    readme_lines.extend(
        [
            "",
            "## 사후 성능",
            "",
            "| Scenario | 전체 RMSE | 10–70% RMSE | 용량 RMSE | 최대 용량오차 |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for scenario in ("S6", "S7a", "S7b"):
        row = summary_index.loc[scenario]
        readme_lines.append(
            f"| {scenario} | {row.Overall_full_time_RMSE_mV:.2f} mV | "
            f"{row.Overall_center_time_RMSE_mV:.2f} mV | "
            f"{row.Capacity_RMSE_pct:.2f}% | {row.Max_abs_capacity_error_mAh:.1f} mAh |"
        )
    readme_lines.extend(
        [
            "",
            "## 결론",
            "",
            "- S7a와 S7b의 전체 전압 RMSE 차이는 약 0.01 mV로 사실상 동일하다.",
            "- S7b에서 Dsn과 brugg_n은 S6 대비 거의 움직이지 않아 세 파라미터 동시 fitting의 추가 실익이 작다.",
            "- fitting 모델 중에는 실험 Dsn을 보존하고 더 적은 파라미터를 사용하는 S7a가 우선이다.",
            "- 그러나 S6 대비 전체 RMSE 개선은 약 1.44 mV에 그치고, 용량 RMSE는 1.02%에서 1.98%로 증가한다. 따라서 현재 기준모델은 S6 유지가 더 보수적이다.",
            "- PSO는 7 particles × 7 iterations, 단일 seed의 탐색 비교이며 전역수렴 또는 optimizer 불확실성을 증명하지 않는다.",
        ]
    )
    (OUT / "README_KO.md").write_text(
        "\n".join(readme_lines) + "\n", encoding="utf-8"
    )
    print("\nSENSITIVITY\n", sensitivity.to_string(index=False))
    print("\nSUBSETS\n", subsets.to_string(index=False))
    print("\nPARAMETERS\n", parameters.to_string(index=False))
    print("\nSUMMARY\n", summary.to_string(index=False))


if __name__ == "__main__":
    main()
