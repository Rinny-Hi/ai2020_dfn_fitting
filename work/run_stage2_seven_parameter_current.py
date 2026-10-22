"""Current two-OCP, seven-parameter voltage LSA including electrolyte Bruggeman exponents."""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from run_stage2_two_ocp import cases, kinetic_constants, scope, stage_with_values


OUT = Path(__file__).resolve().parents[1] / "outputs"
OUT.mkdir(exist_ok=True)
PARAMETERS = ("Dsn", "Dsp", "kn", "kp", "brugg_n", "brugg_s", "brugg_p")
BRUGG_KEYS = {
    "brugg_n": "Negative electrode Bruggeman coefficient (electrolyte)",
    "brugg_s": "Separator Bruggeman coefficient (electrolyte)",
    "brugg_p": "Positive electrode Bruggeman coefficient (electrolyte)",
}
CONDITIONS = tuple((rate, charge) for rate in (0.5, 1.0, 2.0) for charge in (True, False))
STEP = math.log(1.05)


def stage_with_seven(case, values):
    stage = stage_with_values(case, values)
    parameters = stage["parameters"].copy()
    parameters.update({key: values[name] for name, key in BRUGG_KEYS.items()})
    stage["parameters"] = parameters
    return stage


def run_six(stage):
    return {condition: scope["run_dfn_branch"](stage, *condition) for condition in CONDITIONS}


def metrics(case_name, phase, runs):
    rows = []
    for (rate, charge), run in runs.items():
        metric = scope["compare_with_experiment"](scope["experimental_curves"][(rate, charge)], run)
        rows.append({"OCP_case": case_name, "Phase": phase, "C_rate": rate,
                     "Direction": "charge" if charge else "discharge", **metric})
    return rows


def grids_for(runs):
    grids = {}
    for condition in CONDITIONS:
        exp = scope["experimental_curves"][condition]
        first_end = min(run[condition]["time_min"][-1] for run in runs)
        t0 = 0.10 * scope["CELL_CAPACITY_AH"] / exp["current_a"] * 60
        t1 = min(0.70 * scope["CELL_CAPACITY_AH"] / exp["current_a"] * 60,
                 0.98 * first_end, 0.98 * exp["time_min"][-1])
        if t1 <= t0:
            raise RuntimeError(f"No common sensitivity interval for {condition}")
        grids[condition] = np.linspace(t0, t1, 80)
    return grids


def sensitivity(case_name, case, values, base_runs):
    variants = {}
    for parameter in PARAMETERS:
        for sign in (-1, +1):
            perturbed = dict(values)
            perturbed[parameter] *= math.exp(sign * STEP)
            print("SENS", case_name, parameter, sign, flush=True)
            variants[(parameter, sign)] = run_six(stage_with_seven(case, perturbed))
    grids = grids_for([base_runs, *variants.values()])
    nominal = np.concatenate([
        np.interp(grids[condition], base_runs[condition]["time_min"], base_runs[condition]["voltage_v"])
        for condition in CONDITIONS
    ])
    columns = []
    rows = []
    condition_rows = []
    for parameter in PARAMETERS:
        derivative_blocks = []
        nonlinear_blocks = []
        difference_blocks = []
        for condition in CONDITIONS:
            grid = grids[condition]
            plus = variants[(parameter, +1)][condition]
            minus = variants[(parameter, -1)][condition]
            base = base_runs[condition]
            vp = np.interp(grid, plus["time_min"], plus["voltage_v"])
            vm = np.interp(grid, minus["time_min"], minus["voltage_v"])
            v0 = np.interp(grid, base["time_min"], base["voltage_v"])
            derivative = (vp - vm) / (2 * STEP)
            derivative_blocks.append(derivative)
            nonlinear_blocks.append(vp + vm - 2 * v0)
            difference_blocks.append(vp - vm)
            condition_rows.append({
                "OCP_case": case_name, "Parameter": parameter,
                "C_rate": condition[0], "Direction": "charge" if condition[1] else "discharge",
                "RMS_mV_per_log_parameter": 1000 * float(np.sqrt(np.mean(derivative**2))),
                "Window_start_min": float(grid[0]), "Window_end_min": float(grid[-1]),
            })
        derivative = np.concatenate(derivative_blocks)
        normalized = derivative / nominal
        columns.append(normalized)
        rows.append({
            "OCP_case": case_name, "Parameter": parameter, "Baseline": values[parameter],
            "RMS_normalized_sensitivity": float(np.sqrt(np.mean(normalized**2))),
            "RMS_mV_per_log_parameter": 1000 * float(np.sqrt(np.mean(derivative**2))),
            "Voltage_change_for_5pct_mV_RMS": 1000 * STEP * float(np.sqrt(np.mean(derivative**2))),
            "Central_nonlinearity_ratio": float(
                np.linalg.norm(np.concatenate(nonlinear_blocks))
                / np.linalg.norm(np.concatenate(difference_blocks))
            ),
        })
    matrix = np.column_stack(columns)
    corr = pd.DataFrame(np.corrcoef(matrix, rowvar=False), index=PARAMETERS, columns=PARAMETERS)
    rank = pd.DataFrame(rows).sort_values("RMS_normalized_sensitivity", ascending=False)
    rank["Relative_sensitivity"] = rank.RMS_normalized_sensitivity / rank.RMS_normalized_sensitivity.max()
    subset = {}
    for threshold in (.8, .95, .99):
        chosen = []
        for candidate in rank.Parameter:
            if all(abs(corr.loc[candidate, selected]) < threshold for selected in chosen):
                chosen.append(candidate)
        subset[str(threshold)] = chosen
    singular = np.linalg.svd(matrix / np.linalg.norm(matrix, axis=0), compute_uv=False)
    diagnostic = {
        "OCP_case": case_name,
        "Observation_count": int(len(matrix)),
        "Candidate_parameters": list(PARAMETERS),
        "Singular_values_unit_columns": [float(x) for x in singular],
        "Condition_number_unit_columns": float(singular[0] / singular[-1]),
        "Correlation_greedy_subsets": subset,
    }
    return rank, corr, condition_rows, diagnostic


def main():
    all_metrics = []
    all_sensitivities = []
    all_condition_sensitivities = []
    all_diagnostics = []
    baseline = {
        "parameters": scope["common_parameters"],
        "window": scope["PUBLISHED_WINDOW"],
        "c_s_max_n": scope["baseline_csmax_n"],
        "c_s_max_p": scope["baseline_csmax_p"],
        "anode_ocp": scope["AI2020_OCP"]["negative"],
        "cathode_ocp": scope["AI2020_OCP"]["positive"],
    }
    print("RUN common geometry + Ai2020 OCP/nominal transport/kinetics baseline", flush=True)
    all_metrics.extend(metrics("Ai2020_common_baseline", "baseline", run_six(baseline)))
    for case_name, case in cases().items():
        values = kinetic_constants(case)
        values.update({name: float(case["parameters"][key]) for name, key in BRUGG_KEYS.items()})
        print("RUN final", case_name, values, flush=True)
        runs = run_six(stage_with_seven(case, values))
        all_metrics.extend(metrics(case_name, "experimental_all", runs))
        rank, corr, detail, diagnostic = sensitivity(case_name, case, values, runs)
        all_sensitivities.append(rank)
        all_condition_sensitivities.extend(detail)
        all_diagnostics.append(diagnostic)
        corr.to_csv(OUT / f"Stage2_{case_name}_seven_parameter_Pearson.csv", encoding="utf-8-sig")
    branches = pd.DataFrame(all_metrics)
    summary = branches.groupby(["OCP_case", "Phase"], sort=False).agg(
        Full_RMSE_mV=("full_RMSE_mV", lambda x: float(np.sqrt(np.mean(np.square(x))))),
        Center_10_70_RMSE_mV=("center_10_70_RMSE_mV", lambda x: float(np.sqrt(np.mean(np.square(x))))),
        Capacity_error_RMSE_pct=("capacity_error_pct", lambda x: float(np.sqrt(np.mean(np.square(x))))),
    ).reset_index()
    branches.to_csv(OUT / "Stage2_full_experimental_branch_comparison.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(OUT / "Stage2_full_experimental_summary.csv", index=False, encoding="utf-8-sig")
    pd.concat(all_sensitivities, ignore_index=True).to_csv(
        OUT / "Stage2_seven_parameter_sensitivity.csv", index=False, encoding="utf-8-sig"
    )
    pd.DataFrame(all_condition_sensitivities).to_csv(
        OUT / "Stage2_seven_parameter_condition_sensitivity.csv", index=False, encoding="utf-8-sig"
    )
    (OUT / "Stage2_seven_parameter_identifiability.json").write_text(
        json.dumps(all_diagnostics, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("\nSUMMARY\n", summary.to_string(index=False), flush=True)
    print("\nSENSITIVITY\n", pd.concat(all_sensitivities).to_string(index=False), flush=True)
    print("\nIDENTIFIABILITY\n", json.dumps(all_diagnostics, indent=2), flush=True)


if __name__ == "__main__":
    main()
