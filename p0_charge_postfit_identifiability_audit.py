"""Audit post-fit Jacobians without differentiating through cutoff penalties."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

import charge_physical_time_common as common


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "results" / "260927_charge_subset_fit_validation"
OUT = ROOT / "results" / "260927_charge_postfit_identifiability_audit"
STEPS = (0.005, 0.01, 0.02)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    context = common.build_context()
    summary_rows = []
    derivative_rows = []
    for model_name, parameters in common.MODEL_PARAMETERS.items():
        best = json.loads((SOURCE / model_name / "best_result.json").read_text(encoding="utf-8"))
        values = {name: float(best[name]) for name in common.INITIAL}
        center = common.values_to_normalized(values, parameters)
        center_runs = common.simulate_direction(context, values, True)
        center_residual, center_feasibility = common.charge_residual(context, center_runs)
        if center_residual is None:
            raise RuntimeError(f"Stored optimum is not feasible: {model_name}")
        for step in STEPS:
            columns = []
            modes = []
            for index, parameter in enumerate(parameters):
                lower = center.copy()
                upper = center.copy()
                lower[index] = max(0.0, center[index] - step)
                upper[index] = min(1.0, center[index] + step)
                lower_values = common.normalized_to_values(lower, parameters)
                upper_values = common.normalized_to_values(upper, parameters)
                lower_residual, lower_feasibility = common.charge_residual(
                    context, common.simulate_direction(context, lower_values, True)
                )
                upper_residual, upper_feasibility = common.charge_residual(
                    context, common.simulate_direction(context, upper_values, True)
                )
                if lower_residual is not None and upper_residual is not None:
                    derivative = (upper_residual - lower_residual) / (
                        upper[index] - lower[index]
                    )
                    mode = "central"
                elif upper_residual is not None:
                    derivative = (upper_residual - center_residual) / (
                        upper[index] - center[index]
                    )
                    mode = "forward_feasible"
                elif lower_residual is not None:
                    derivative = (center_residual - lower_residual) / (
                        center[index] - lower[index]
                    )
                    mode = "backward_feasible"
                else:
                    raise RuntimeError(
                        f"No feasible finite-difference side: {model_name} {parameter} {step}"
                    )
                columns.append(derivative)
                modes.append(mode)
                derivative_rows.append(
                    {
                        "Model": model_name,
                        "Step": step,
                        "Parameter": parameter,
                        "Difference_mode": mode,
                        "Lower_feasible": bool(lower_feasibility["feasible"]),
                        "Upper_feasible": bool(upper_feasibility["feasible"]),
                        "RMS_mV_per_normalized_coordinate": float(
                            np.sqrt(np.mean(derivative**2))
                        ),
                    }
                )
            matrix = np.column_stack(columns)
            rms = np.sqrt(np.mean(matrix**2, axis=0))
            relative = rms / np.max(rms)
            unit = matrix / np.maximum(np.linalg.norm(matrix, axis=0), np.finfo(float).eps)
            cosine = unit.T @ unit
            singular = np.linalg.svd(unit, compute_uv=False)
            condition = float(singular[0] / singular[-1])
            max_pair = float(np.max(np.abs(cosine - np.eye(len(parameters)))))
            summary_rows.append(
                {
                    "Model": model_name,
                    "Step": step,
                    "Condition_number": condition,
                    "Max_abs_pairwise_cosine": max_pair,
                    "Min_relative_sensitivity": float(np.min(relative)),
                    "Difference_modes": ",".join(modes),
                    "Pass": bool(
                        np.linalg.matrix_rank(matrix) == len(parameters)
                        and condition <= 20
                        and max_pair <= 0.95
                        and np.min(relative) >= 0.10
                    ),
                }
            )
            np.save(OUT / f"jacobian_{model_name}_step_{str(step).replace('.', 'p')}.npy", matrix)
            pd.DataFrame(cosine, index=parameters, columns=parameters).to_csv(
                OUT / f"cosine_{model_name}_step_{str(step).replace('.', 'p')}.csv",
                encoding="utf-8-sig",
            )
    summary = pd.DataFrame(summary_rows)
    detail = pd.DataFrame(derivative_rows)
    summary.to_csv(OUT / "postfit_identifiability.csv", index=False, encoding="utf-8-sig")
    detail.to_csv(OUT / "postfit_derivative_detail.csv", index=False, encoding="utf-8-sig")
    (OUT / "analysis_manifest.json").write_text(
        json.dumps(
            {
                "reason": "best fits lie on the full-time feasibility boundary, so central differences that cross early cutoff differentiate the penalty rather than voltage",
                "method": "central difference when both sides are feasible; otherwise use the feasible one-sided voltage derivative",
                "steps": STEPS,
                "capacity_in_analysis": False,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    print(summary.to_string(index=False))
    print(detail.to_string(index=False))


if __name__ == "__main__":
    main()
