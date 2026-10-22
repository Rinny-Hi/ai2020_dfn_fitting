"""Two-OCP Stage 2: experimental Ds/kinetics and voltage identifiability.

Uses the executed Stage-1 notebook definitions, with the equivalent compact
GITT rest-end table to avoid reparsing two very large raw record sheets.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from time import monotonic

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from stage2_setup_probe import scope


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)

FARADAY = 96485.33212
R_GAS = 8.314462618
T_K = 298.15
C_E = 1000.0
D_SN = 4.1707741459682473e-14
D_SP = 6.888175323807748e-14
R_CT_N = 27.000934087550903  # current 50% SoC, N3+N4 local-minimum integral
R_CT_P = 25.387  # measured cathode EIS-DRT slide
L_N, EPS_N, R_N = 76.5e-6, 0.610, 3.692e-6
L_P, EPS_P, R_P = 68.0e-6, 0.620, 4.239e-6
AI2020_K_REF = 1e-11 * FARADAY
PARAMETERS = ("Dsn", "Dsp", "kn", "kp")
RATES = (0.5, 1.0, 2.0)
CONDITIONS = tuple((rate, charge) for rate in RATES for charge in (True, False))
LOG_STEP = math.log(1.05)


def j0_from_rct(rct_ohm_cm2: float, eps_s: float, radius_m: float, thickness_m: float) -> float:
    area_resistance = rct_ohm_cm2 * 1e-4
    surface_area = 3.0 * eps_s / radius_m
    return R_GAS * T_K / (FARADAY * area_resistance * surface_area * thickness_m)


J0_N = j0_from_rct(R_CT_N, EPS_N, R_N, L_N)
J0_P = j0_from_rct(R_CT_P, EPS_P, R_P, L_P)


def cases():
    return {
        "Ai2020_anode__measured_cathode": {
            "parameters": scope["improved_parameters"],
            "window": scope["IMPROVED_WINDOW"],
            "c_s_max_n": scope["improved_csmax_n"],
            "c_s_max_p": scope["improved_csmax_p"],
            "anode_ocp": scope["AI2020_OCP"]["negative"],
            "cathode_ocp": scope["CATHODE_GITT_OCP"],
        },
        "measured_anode__measured_cathode": {
            "parameters": scope["measured_anode_parameters"],
            "window": scope["MEASURED_ANODE_WINDOW"],
            "c_s_max_n": scope["measured_anode_csmax_n"],
            "c_s_max_p": scope["measured_anode_csmax_p"],
            "anode_ocp": scope["MEASURED_ANODE_OCP"],
            "cathode_ocp": scope["CATHODE_GITT_OCP"],
        },
    }


def kinetic_constants(case):
    window = case["window"]
    x50 = (window["x0"] + window["x100"]) / 2.0
    y50 = (window["y0"] + window["y100"]) / 2.0
    kn = J0_N / (case["c_s_max_n"] * math.sqrt(C_E * x50 * (1 - x50)))
    kp = J0_P / (case["c_s_max_p"] * math.sqrt(C_E * y50 * (1 - y50)))
    return {"Dsn": D_SN, "Dsp": D_SP, "kn": kn, "kp": kp, "x50": x50, "y50": y50}


def stage_with_values(case, values, use_ds=True, use_kp=True, use_kn=True):
    changed = dict(case)
    params = case["parameters"].copy()
    updates = {}
    if use_ds:
        updates.update({
            "Negative particle diffusivity [m2.s-1]": values["Dsn"],
            "Positive particle diffusivity [m2.s-1]": values["Dsp"],
        })
    for use, label, parameter_name in (
        (use_kp, "kp", "Positive electrode exchange-current density [A.m-2]"),
        (use_kn, "kn", "Negative electrode exchange-current density [A.m-2]"),
    ):
        if use:
            original = scope["raw_parameters"][parameter_name]
            factor = values[label] / AI2020_K_REF

            def scaled(c_e, c_s_surf, c_s_max, temperature, original=original, factor=factor):
                return factor * original(c_e, c_s_surf, c_s_max, temperature)

            updates[parameter_name] = scaled
    params.update(updates)
    changed["parameters"] = params
    return changed


def run_conditions(stage, label):
    runs = {}
    for rate, charge in CONDITIONS:
        direction = "charge" if charge else "discharge"
        print("RUN", label, rate, direction, flush=True)
        runs[(rate, charge)] = scope["run_dfn_branch"](stage, rate, charge)
    return runs


def summarize_runs(case_name, phase, runs):
    rows = []
    for condition, run in runs.items():
        rate, charge = condition
        metrics = scope["compare_with_experiment"](scope["experimental_curves"][condition], run)
        rows.append({"OCP_case": case_name, "Phase": phase, "C_rate": rate,
                     "Direction": "charge" if charge else "discharge", **metrics})
    return rows


def build_sensitivity(case_name, case, values, base_runs):
    variants = {}
    for parameter in PARAMETERS:
        for sign in (-1, 1):
            changed = dict(values)
            changed[parameter] *= math.exp(sign * LOG_STEP)
            stage = stage_with_values(case, changed)
            variants[(parameter, sign)] = run_conditions(stage, f"{case_name}/{parameter}/{sign:+d}")

    # A common fixed-time window per condition is shared by all Jacobian columns.
    grids = {}
    for condition in CONDITIONS:
        exp = scope["experimental_curves"][condition]
        rate, charge = condition
        earliest_end = min(
            [base_runs[condition]["time_min"][-1]]
            + [runs[condition]["time_min"][-1] for runs in variants.values()]
        )
        t0 = 0.10 * scope["CELL_CAPACITY_AH"] / exp["current_a"] * 60.0
        t1 = min(0.70 * scope["CELL_CAPACITY_AH"] / exp["current_a"] * 60.0,
                 0.98 * earliest_end, 0.98 * exp["time_min"][-1])
        if t1 <= t0:
            raise RuntimeError(f"No common 10–70% capacity window: {case_name} {condition}")
        grids[condition] = np.linspace(t0, t1, 80)

    columns = []
    rows = []
    condition_rows = []
    baseline_voltage = []
    for condition in CONDITIONS:
        run = base_runs[condition]
        baseline_voltage.append(np.interp(grids[condition], run["time_min"], run["voltage_v"]))
    nominal = np.concatenate(baseline_voltage)

    for parameter in PARAMETERS:
        raw_blocks = []
        nonlin_blocks = []
        for condition in CONDITIONS:
            grid = grids[condition]
            plus = variants[(parameter, +1)][condition]
            minus = variants[(parameter, -1)][condition]
            v_plus = np.interp(grid, plus["time_min"], plus["voltage_v"])
            v_minus = np.interp(grid, minus["time_min"], minus["voltage_v"])
            base = np.interp(grid, base_runs[condition]["time_min"], base_runs[condition]["voltage_v"])
            derivative = (v_plus - v_minus) / (2.0 * LOG_STEP)
            raw_blocks.append(derivative)
            nonlin_blocks.append((v_plus + v_minus - 2 * base, v_plus - v_minus))
            rate, charge = condition
            condition_rows.append({
                "OCP_case": case_name, "C_rate": rate,
                "Direction": "charge" if charge else "discharge",
                "Parameter": parameter,
                "RMS_mV_per_log_parameter": 1000 * float(np.sqrt(np.mean(derivative**2))),
                "MaxAbs_mV_per_log_parameter": 1000 * float(np.max(np.abs(derivative))),
                "Window_start_min": float(grid[0]), "Window_end_min": float(grid[-1]),
            })
        derivative = np.concatenate(raw_blocks)
        normalized = derivative / nominal
        columns.append(normalized)
        numerator = np.linalg.norm(np.concatenate([x for x, _ in nonlin_blocks]))
        denominator = np.linalg.norm(np.concatenate([y for _, y in nonlin_blocks]))
        rows.append({
            "OCP_case": case_name, "Parameter": parameter,
            "Baseline": values[parameter],
            "RMS_normalized_sensitivity": float(np.sqrt(np.mean(normalized**2))),
            "RMS_mV_per_log_parameter": 1000 * float(np.sqrt(np.mean(derivative**2))),
            "Voltage_change_for_5pct_mV_RMS": 1000 * LOG_STEP * float(np.sqrt(np.mean(derivative**2))),
            "Central_nonlinearity_ratio": float(numerator / denominator),
        })

    matrix = np.column_stack(columns)
    corr = np.corrcoef(matrix, rowvar=False)
    rank = pd.DataFrame(rows).sort_values("RMS_normalized_sensitivity", ascending=False)
    rank["Relative_sensitivity"] = rank["RMS_normalized_sensitivity"] / rank["RMS_normalized_sensitivity"].max()
    selected = {}
    for beta in (0.8, 0.95, 0.99):
        chosen = []
        for name in rank.Parameter:
            index = PARAMETERS.index(name)
            if all(abs(corr[index, PARAMETERS.index(other)]) < beta for other in chosen):
                chosen.append(name)
        selected[str(beta)] = chosen

    unit_columns = matrix / np.linalg.norm(matrix, axis=0)
    singular = np.linalg.svd(unit_columns, compute_uv=False)
    diagnostics = {
        "OCP_case": case_name,
        "N_observations": int(len(matrix)),
        "N_parameters": len(PARAMETERS),
        "Pearson_abs_max_offdiagonal": float(np.max(np.abs(corr - np.eye(len(PARAMETERS))))),
        "Normalized_column_singular_values": [float(x) for x in singular],
        "Normalized_column_condition_number": float(singular[0] / singular[-1]),
        "Correlation_greedy_subsets": selected,
    }
    return rank, pd.DataFrame(corr, index=PARAMETERS, columns=PARAMETERS), pd.DataFrame(condition_rows), diagnostics


def main():
    start = monotonic()
    cases_dict = cases()
    parameter_rows = []
    branch_rows = []
    sensitivities = []
    correlations = {}
    condition_sensitivities = []
    diagnostics = []
    for case_name, case in cases_dict.items():
        values = kinetic_constants(case)
        parameter_rows.append({
            "OCP_case": case_name,
            "x0": case["window"]["x0"], "x100": case["window"]["x100"],
            "y0": case["window"]["y0"], "y100": case["window"]["y100"],
            "c_s_max_n_mol_m3": case["c_s_max_n"], "c_s_max_p_mol_m3": case["c_s_max_p"],
            "r_n_um": R_N * 1e6, "r_p_um": R_P * 1e6,
            "Dsn_m2_s": D_SN, "Dsp_m2_s": D_SP,
            "Rct_n_ohm_cm2": R_CT_N, "Rct_p_ohm_cm2": R_CT_P,
            "j0_n_A_m2": J0_N, "j0_p_A_m2": J0_P,
            "x_n_at_50pct": values["x50"], "y_p_at_50pct": values["y50"],
            "kn": values["kn"], "kp": values["kp"],
        })
        print("CASE", case_name, "k", values["kn"], values["kp"], flush=True)
        phases = (
            ("OCP_only", False, False, False),
            ("plus_Ds_GITT", True, False, False),
            ("plus_kp_EIS", True, True, False),
            ("plus_kn_DRT_initial", True, True, True),
        )
        final_runs = None
        for label, use_ds, use_kp, use_kn in phases:
            stage = stage_with_values(case, values, use_ds, use_kp, use_kn)
            runs = run_conditions(stage, case_name + "/" + label)
            branch_rows.extend(summarize_runs(case_name, label, runs))
            if use_kn:
                final_runs = runs
        assert final_runs is not None
        rank, corr, detail, diagnostic = build_sensitivity(case_name, case, values, final_runs)
        sensitivities.append(rank)
        correlations[case_name] = corr
        condition_sensitivities.append(detail)
        diagnostics.append(diagnostic)

    parameters = pd.DataFrame(parameter_rows)
    branches = pd.DataFrame(branch_rows)
    aggregate = branches.groupby(["OCP_case", "Phase"], sort=False).agg(
        Full_RMSE_mV=("full_RMSE_mV", lambda x: float(np.sqrt(np.mean(np.square(x))))),
        Center_10_70_RMSE_mV=("center_10_70_RMSE_mV", lambda x: float(np.sqrt(np.mean(np.square(x))))),
        Mean_abs_capacity_error_pct=("capacity_error_pct", lambda x: float(np.mean(np.abs(x)))),
    ).reset_index()
    ranking = pd.concat(sensitivities, ignore_index=True)
    detail = pd.concat(condition_sensitivities, ignore_index=True)

    for name, frame in (
        ("Stage2_two_OCP_parameters.csv", parameters),
        ("Stage2_two_OCP_branch_metrics.csv", branches),
        ("Stage2_two_OCP_phase_summary.csv", aggregate),
        ("Stage2_two_OCP_sensitivity.csv", ranking),
        ("Stage2_two_OCP_condition_sensitivity.csv", detail),
    ):
        frame.to_csv(OUT / name, index=False, encoding="utf-8-sig")
    for case_name, matrix in correlations.items():
        matrix.to_csv(OUT / f"Stage2_{case_name}_Pearson_correlation.csv", encoding="utf-8-sig")
    (OUT / "Stage2_two_OCP_identifiability.json").write_text(
        json.dumps(diagnostics, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), layout="constrained")
    for ax, (case_name, group) in zip(axes, aggregate.groupby("OCP_case", sort=False)):
        x = np.arange(len(group))
        ax.plot(x, group["Center_10_70_RMSE_mV"], "o-", color="#2458b4", lw=2)
        ax.set_xticks(x, ["OCP", "+Ds", "+kp", "+kn"])
        ax.set_title(case_name.replace("__", " + ").replace("_", " "), fontsize=10)
        ax.set_ylabel("10–70% window RMSE (mV)")
        ax.grid(alpha=0.25)
    fig.savefig(OUT / "Stage2_two_OCP_phase_RMSE.png", dpi=180)
    plt.close(fig)

    print("\nPARAMETERS\n", parameters.to_string(index=False), flush=True)
    print("\nAGGREGATE\n", aggregate.to_string(index=False), flush=True)
    print("\nSENSITIVITY\n", ranking.to_string(index=False), flush=True)
    print("\nIDENTIFIABILITY\n", json.dumps(diagnostics, ensure_ascii=False, indent=2), flush=True)
    print("ELAPSED", round(monotonic() - start, 1), "seconds", flush=True)


if __name__ == "__main__":
    main()
