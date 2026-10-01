"""Physical-time sensitivity and local identifiability for the chosen subsets."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import charge_physical_time_common as common


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "260927_charge_physical_time_lsa"
ANCHOR = {**common.INITIAL, "kn": 1.90e-6}
STEPS = (0.01, 0.025, 0.05)
MAIN_STEP = 0.025


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    context = common.build_context()
    simulation_cache: dict[tuple[float, float, float], tuple[dict, np.ndarray | None, dict]] = {}

    def evaluate(values: dict[str, float]):
        key = tuple(float(values[name]) for name in ("Dsn", "kn", "brugg_n"))
        if key not in simulation_cache:
            print("simulate", {name: f"{values[name]:.6g}" for name in common.INITIAL}, flush=True)
            runs = common.simulate_direction(context, values, True)
            residual, feasibility = common.charge_residual(context, runs)
            simulation_cache[key] = (runs, residual, feasibility)
        return simulation_cache[key]

    baseline_runs = common.simulate_direction(context, common.INITIAL, True)
    baseline_residual, baseline_feasibility = common.charge_residual(context, baseline_runs)
    anchor_runs, anchor_residual, anchor_feasibility = evaluate(ANCHOR)
    if anchor_residual is None:
        raise RuntimeError("Chosen physical-time LSA anchor is not full-grid feasible")

    sensitivity_rows = []
    correlation_rows = []
    subset_rows = []
    matrices = {}
    for step in STEPS:
        derivatives = {}
        for parameter in common.MODEL_PARAMETERS["main"]:
            center = common.values_to_normalized(ANCHOR, (parameter,))[0]
            lo_z, hi_z = center - step, center + step
            if lo_z < 0 or hi_z > 1:
                raise RuntimeError(f"Finite-difference step leaves bound for {parameter}")
            lower = common.normalized_to_values([lo_z], (parameter,))
            upper = common.normalized_to_values([hi_z], (parameter,))
            lower.update({k: v for k, v in ANCHOR.items() if k != parameter})
            upper.update({k: v for k, v in ANCHOR.items() if k != parameter})
            _, lower_residual, lower_feasible = evaluate(lower)
            _, upper_residual, upper_feasible = evaluate(upper)
            if lower_residual is None or upper_residual is None:
                raise RuntimeError(
                    f"Non-feasible central difference for {parameter}, step={step}: "
                    f"{lower_feasible}, {upper_feasible}"
                )
            # residual = model - experiment; experimental voltage cancels.
            derivatives[parameter] = (upper_residual - lower_residual) / (2.0 * step)

        for model_name, parameters in common.MODEL_PARAMETERS.items():
            matrix = np.column_stack([derivatives[name] for name in parameters])
            matrices[(model_name, step)] = matrix
            rms = np.sqrt(np.mean(matrix**2, axis=0))
            relative = rms / np.max(rms)
            unit = matrix / np.maximum(np.linalg.norm(matrix, axis=0), np.finfo(float).eps)
            cosine = unit.T @ unit
            singular = np.linalg.svd(unit, compute_uv=False)
            condition = float(singular[0] / singular[-1])
            max_pair = float(np.max(np.abs(cosine - np.eye(len(parameters)))))
            full_rank = int(np.linalg.matrix_rank(matrix)) == len(parameters)
            passes = bool(
                full_rank
                and float(np.min(relative)) >= 0.10
                and max_pair <= 0.95
                and condition <= 20.0
            )
            subset_rows.append(
                {
                    "Model": model_name,
                    "Normalized_step": step,
                    "N_parameters": len(parameters),
                    "Matrix_rank": int(np.linalg.matrix_rank(matrix)),
                    "Normalized_condition_number": condition,
                    "Max_abs_pairwise_cosine": max_pair,
                    "Min_relative_sensitivity": float(np.min(relative)),
                    "Pass": passes,
                }
            )
            for name, sensitivity, rel in zip(parameters, rms, relative):
                sensitivity_rows.append(
                    {
                        "Model": model_name,
                        "Normalized_step": step,
                        "Parameter": name,
                        "RMS_mV_per_normalized_coordinate": float(sensitivity),
                        "Relative_sensitivity": float(rel),
                    }
                )
            for i, name_i in enumerate(parameters):
                for j, name_j in enumerate(parameters):
                    correlation_rows.append(
                        {
                            "Model": model_name,
                            "Normalized_step": step,
                            "Parameter_i": name_i,
                            "Parameter_j": name_j,
                            "Sensitivity_cosine": float(cosine[i, j]),
                        }
                    )

    sensitivity = pd.DataFrame(sensitivity_rows)
    correlations = pd.DataFrame(correlation_rows)
    subsets = pd.DataFrame(subset_rows)
    sensitivity.to_csv(OUT / "sensitivity_ranking.csv", index=False, encoding="utf-8-sig")
    correlations.to_csv(OUT / "sensitivity_cosine.csv", index=False, encoding="utf-8-sig")
    subsets.to_csv(OUT / "subset_identifiability.csv", index=False, encoding="utf-8-sig")
    for (model_name, step), matrix in matrices.items():
        np.save(OUT / f"jacobian_{model_name}_step_{str(step).replace('.', 'p')}.npy", matrix)

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.2), constrained_layout=True)
    main_sens = sensitivity[np.isclose(sensitivity.Normalized_step, MAIN_STEP)]
    pivot = main_sens.pivot(index="Parameter", columns="Model", values="Relative_sensitivity")
    pivot.plot(kind="bar", ax=axes[0], color=["#0072B2", "#D55E00"])
    axes[0].axhline(0.10, color="black", ls="--", lw=1)
    axes[0].set(title="Physical-time local sensitivity", ylabel="Relative RMS sensitivity")
    main_corr = correlations[
        (correlations.Model == "main")
        & np.isclose(correlations.Normalized_step, MAIN_STEP)
    ].pivot(index="Parameter_i", columns="Parameter_j", values="Sensitivity_cosine")
    main_corr = main_corr.loc[list(common.MODEL_PARAMETERS["main"]), list(common.MODEL_PARAMETERS["main"])]
    image = axes[1].imshow(main_corr.to_numpy(float), vmin=-1, vmax=1, cmap="coolwarm")
    axes[1].set_xticks(range(3), main_corr.columns, rotation=30)
    axes[1].set_yticks(range(3), main_corr.index)
    axes[1].set_title("Main-model sensitivity cosine")
    for i in range(3):
        for j in range(3):
            axes[1].text(j, i, f"{main_corr.iloc[i, j]:.2f}", ha="center", va="center")
    fig.colorbar(image, ax=axes[1], shrink=0.82)
    fig.savefig(OUT / "physical_time_identifiability.png", dpi=220)
    plt.close(fig)

    manifest = {
        "objective": "ordinary voltage residual on 100-point fixed physical-time grids for each July BoL cell and 0.5C/1C/2C charge branch",
        "capacity_in_objective": False,
        "grid_start": "first CC record reaching 95% of median command current",
        "grid_end": "experimental CC cutoff",
        "baseline": common.INITIAL,
        "fixed": {**common.FIXED, "kp": context.kp},
        "baseline_full_grid_feasibility": baseline_feasibility,
        "baseline_full_grid_RMSE_mV": None
        if baseline_residual is None
        else float(np.sqrt(np.mean(baseline_residual**2))),
        "lsa_anchor": ANCHOR,
        "lsa_anchor_reason": "nearest simple common anchor tested with baseline brugg_n and increased kn that leaves central finite differences feasible on every experimental-time grid",
        "lsa_anchor_feasibility": anchor_feasibility,
        "lsa_anchor_RMSE_mV": float(np.sqrt(np.mean(anchor_residual**2))),
        "parameter_coordinate": "Dsn and kn log-bound normalized to [0,1]; brugg_n linear-bound normalized to [0,1]",
        "finite_difference_steps": STEPS,
        "pass_rule": {
            "full_column_rank": True,
            "minimum_relative_sensitivity": 0.10,
            "maximum_absolute_pairwise_sensitivity_cosine": 0.95,
            "maximum_normalized_condition_number": 20.0,
            "must_be_stable_across_all_steps": True,
        },
    }
    (OUT / "analysis_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("\nSUBSETS\n", subsets.to_string(index=False))
    print("\nSENSITIVITY\n", sensitivity.to_string(index=False))


if __name__ == "__main__":
    main()
