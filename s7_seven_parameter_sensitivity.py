"""Seven-parameter local voltage sensitivity at the current S6 baseline.

The perturbation coordinate is a +/-5% physical effect: k and D are scaled
multiplicatively, while each Bruggeman exponent is shifted so eps**b changes
by the same multiplicative amount.  The voltage residual definition is the
same fixed-time six-branch objective used by ``s7_pso_voltage_fitting.py``.
"""

from __future__ import annotations

from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

import cumulative_six_stage_model_comparison as cumulative
import electrolyte_transport_set_comparison as transport
import kinetics_original_reversion_ablation as kinetics
import s7_pso_voltage_fitting as s7


OUT = Path(__file__).resolve().parent / "results" / "261013_s7_seven_parameter_sensitivity"
PARAMETERS = ("kn", "kp", "Dsn", "Dsp", "brugg_n", "brugg_p", "brugg_s")
BASE = {
    "kn": 9.648533212e-7,
    "kp": 5.899519906708763e-7,
    "Dsn": 4.1707741459682473e-14,
    "Dsp": 6.888175323807748e-14,
    "brugg_n": 2.5,
    "brugg_p": 1.83,
    "brugg_s": 1.5,
}
PARAMETER_KEYS = {
    "Dsn": "Negative particle diffusivity [m2.s-1]",
    "Dsp": "Positive particle diffusivity [m2.s-1]",
    "brugg_n": "Negative electrode Bruggeman coefficient (electrolyte)",
    "brugg_p": "Positive electrode Bruggeman coefficient (electrolyte)",
    "brugg_s": "Separator Bruggeman coefficient (electrolyte)",
}
POROSITY_KEYS = {
    "brugg_n": "Negative electrode porosity",
    "brugg_p": "Positive electrode porosity",
    "brugg_s": "Separator porosity",
}
STEP = float(np.log(1.05))


def apply_values(model: dict, values: dict[str, float]) -> dict:
    changed = kinetics.with_prefactors(model, values["kn"], values["kp"])
    changed = dict(changed)
    changed["params"] = changed["params"].copy()
    changed["params"].update(
        {
            PARAMETER_KEYS["Dsn"]: values["Dsn"],
            PARAMETER_KEYS["Dsp"]: values["Dsp"],
            PARAMETER_KEYS["brugg_n"]: values["brugg_n"],
            PARAMETER_KEYS["brugg_p"]: values["brugg_p"],
            PARAMETER_KEYS["brugg_s"]: values["brugg_s"],
        },
        check_already_exists=False,
    )
    return changed


def perturbed_values(model: dict, name: str, direction: int) -> dict[str, float]:
    values = dict(BASE)
    if name.startswith("brugg_"):
        porosity = float(model["params"][POROSITY_KEYS[name]])
        values[name] = BASE[name] + direction * STEP / np.log(porosity)
    else:
        values[name] = BASE[name] * np.exp(direction * STEP)
    return values


def subset_report(matrix: np.ndarray, names: tuple[str, ...], relative: dict) -> dict:
    unit = matrix / np.maximum(np.linalg.norm(matrix, axis=0), np.finfo(float).eps)
    cosine = unit.T @ unit
    singular = np.linalg.svd(unit, compute_uv=False)
    condition = float(singular[0] / singular[-1]) if len(names) > 1 else 1.0
    max_pair = float(np.max(np.abs(cosine - np.eye(len(names))))) if len(names) > 1 else 0.0
    minimum_sensitivity = min(relative[name] for name in names)
    full_rank = int(np.linalg.matrix_rank(matrix)) == len(names)
    return {
        "Parameters": "+".join(names),
        "N_parameters": len(names),
        "Matrix_rank": int(np.linalg.matrix_rank(matrix)),
        "Normalized_condition_number": condition,
        "Max_abs_pairwise_cosine": max_pair,
        "Min_relative_sensitivity": minimum_sensitivity,
        "Pass": bool(full_rank and condition <= 20.0 and max_pair <= 0.95 and minimum_sensitivity >= 0.10),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    _, _, _, _, experiment, qcell = transport.current_p0()
    entry = cumulative.build_stages(qcell)["S6"]
    base_model = entry["model"]
    grids, targets = s7.make_experimental_grids(experiment)

    columns = []
    rows = []
    for name in PARAMETERS:
        residuals = []
        for direction in (-1, 1):
            values = perturbed_values(base_model, name, direction)
            model = apply_values(base_model, values)
            runs = s7.run_model(model, entry["stage"], entry["anode"], entry["cathode"])
            residual, _ = s7.residual_vector(runs, grids, targets)
            residuals.append(residual)
        derivative = (residuals[1] - residuals[0]) / (2.0 * STEP)
        columns.append(derivative)
        rows.append(
            {
                "Parameter": name,
                "Baseline": BASE[name],
                "RMS_mV_per_log_effect": float(np.sqrt(np.mean(derivative**2))),
                "Max_abs_mV_per_log_effect": float(np.max(np.abs(derivative))),
            }
        )
        print(name, rows[-1]["RMS_mV_per_log_effect"], flush=True)

    jacobian = np.column_stack(columns)
    ranking = pd.DataFrame(rows)
    ranking["Relative_sensitivity"] = ranking.RMS_mV_per_log_effect / ranking.RMS_mV_per_log_effect.max()
    relative = dict(zip(ranking.Parameter, ranking.Relative_sensitivity))

    unit = jacobian / np.maximum(np.linalg.norm(jacobian, axis=0), np.finfo(float).eps)
    cosine = unit.T @ unit
    subset_rows = []
    for size in (1, 2, 3):
        for names in combinations(PARAMETERS, size):
            indexes = [PARAMETERS.index(name) for name in names]
            subset_rows.append(subset_report(jacobian[:, indexes], names, relative))

    ranking.sort_values("Relative_sensitivity", ascending=False).to_csv(
        OUT / "sensitivity_ranking.csv", index=False, encoding="utf-8-sig"
    )
    pd.DataFrame(cosine, index=PARAMETERS, columns=PARAMETERS).to_csv(
        OUT / "sensitivity_cosine.csv", encoding="utf-8-sig"
    )
    subsets = pd.DataFrame(subset_rows)
    subsets.to_csv(OUT / "subset_identifiability_up_to_3.csv", index=False, encoding="utf-8-sig")
    np.save(OUT / "voltage_jacobian.npy", jacobian)
    print("\nRANKING\n", ranking.sort_values("Relative_sensitivity", ascending=False).to_string(index=False))
    print("\nPASSING THREE-PARAMETER SUBSETS\n", subsets[(subsets.N_parameters == 3) & subsets.Pass].sort_values("Normalized_condition_number").to_string(index=False))


if __name__ == "__main__":
    main()
