"""One-shot Dsn upper-bound expansion diagnostic required by the fit design."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

import charge_physical_time_common as common


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "results" / "260927_charge_subset_fit_validation"
OUT = ROOT / "results" / "260927_charge_dsn_boundary_test"
EXPANDED_DSN = (1.0e-14, 1.2e-13)


def encode(values: dict[str, float], parameters: tuple[str, ...]) -> np.ndarray:
    result = []
    for name in parameters:
        lo, hi = EXPANDED_DSN if name == "Dsn" else common.BOUNDS[name]
        if name in {"Dsn", "kn"}:
            result.append((np.log(values[name]) - np.log(lo)) / (np.log(hi) - np.log(lo)))
        else:
            result.append((values[name] - lo) / (hi - lo))
    return np.asarray(result, float)


def decode(coordinates, parameters: tuple[str, ...]) -> dict[str, float]:
    values = dict(common.INITIAL)
    for z, name in zip(coordinates, parameters):
        lo, hi = EXPANDED_DSN if name == "Dsn" else common.BOUNDS[name]
        if name in {"Dsn", "kn"}:
            values[name] = float(np.exp(np.log(lo) + float(z) * (np.log(hi) - np.log(lo))))
        else:
            values[name] = float(lo + float(z) * (hi - lo))
    return values


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    context = common.build_context()
    rows = []
    for model_name, parameters in common.MODEL_PARAMETERS.items():
        original = json.loads((SOURCE / model_name / "best_result.json").read_text(encoding="utf-8"))
        start_values = {name: float(original[name]) for name in common.INITIAL}
        cache = {}

        def evaluate(coordinates):
            key = tuple(np.round(np.clip(coordinates, 0, 1), 12))
            if key not in cache:
                values = decode(key, parameters)
                runs = common.simulate_direction(context, values, True)
                exact, feasibility = common.charge_residual(context, runs)
                if exact is None:
                    exact, _ = common.charge_residual(context, runs, invalid_penalty_mV=1000.0)
                cache[key] = exact, feasibility, values
                print(
                    model_name,
                    len(cache),
                    f"{np.sqrt(np.mean(exact**2)):.4f} mV",
                    "ok" if feasibility["feasible"] else "invalid",
                    flush=True,
                )
            return cache[key]

        x0 = encode(start_values, parameters)
        result = least_squares(
            lambda x: evaluate(x)[0],
            x0,
            bounds=(np.zeros(len(parameters)), np.ones(len(parameters))),
            method="trf",
            x_scale="jac",
            diff_step=1.0e-3,
            max_nfev=80 if model_name == "main" else 60,
            ftol=1.0e-4,
            xtol=1.0e-4,
            gtol=1.0e-4,
        )
        residual, feasibility, values = evaluate(result.x)
        distance_to_upper = np.log(EXPANDED_DSN[1] / values["Dsn"]) / np.log(
            EXPANDED_DSN[1] / EXPANDED_DSN[0]
        )
        rows.append(
            {
                "Model": model_name,
                "Original_Dsn_upper": common.BOUNDS["Dsn"][1],
                "Expanded_Dsn_upper": EXPANDED_DSN[1],
                "Original_RMSE_mV": original["Charge_objective_RMSE_mV"],
                "Expanded_RMSE_mV": float(np.sqrt(np.mean(residual**2))),
                "RMSE_improvement_mV": float(original["Charge_objective_RMSE_mV"] - np.sqrt(np.mean(residual**2))),
                "Feasible": bool(feasibility["feasible"]),
                "Minimum_time_margin_min": feasibility["minimum_time_margin_min"],
                "Dsn": values["Dsn"],
                "kn": values["kn"],
                "brugg_n": values["brugg_n"],
                "Dsn_log_distance_to_expanded_upper_fraction": float(distance_to_upper),
                "Touches_expanded_upper_5pct": bool(distance_to_upper <= 0.05),
                "NFEV": int(result.nfev),
                "Success": bool(result.success),
                "Message": str(result.message),
            }
        )
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "boundary_test.csv", index=False, encoding="utf-8-sig")
    (OUT / "analysis_manifest.json").write_text(
        json.dumps(
            {
                "purpose": "diagnose whether the first-pass Dsn upper bound truncates the voltage optimum",
                "method": "single TRF restart from each first-pass best point with only Dsn upper bound expanded",
                "capacity_in_objective": False,
                "expanded_Dsn_bound_m2_s": EXPANDED_DSN,
                "interpretation": "diagnostic range, not automatic physical acceptance range",
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    print(frame.to_string(index=False))


if __name__ == "__main__":
    main()
