"""Stage 2 sensitivity/identifiability for kinetics, diffusion, and all Bruggeman coefficients.

The current OCP, stoichiometry window, geometry, EIS kinetics and nominal
Bruggeman settings are used as the baseline.  C-rate and HPPC sensitivities are
computed independently so that either dataset can remain held out later.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybamm

import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import mass_independent_nominal_anode_trial as mit
import nominal_radius_current_configuration as current
import ocp_grounded_candidate_search as common
import old_notebook_method_comparison as omc
import ocp_rn_sweep_fixed_transport as rs
import restored_brugg_de_csmax_audit as audit


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "results" / "260920_mass_independent_nominal_anode"
RESULTS = ROOT / "results" / "260921_stage2_seven_parameter_sensitivity"
RESULTS.mkdir(parents=True, exist_ok=True)

HPPC_FILE = Path(
    r"C:\Users\user\OneDrive\01. 적응형 충전 프로토콜\4. Enertech 파우치셀 실험데이터"
    r"\260808 Enertech셀 열화데이터\260808 Eneterch셀 열화데이터 v1 6-6.xlsx"
)

PARAMETERS = ("kn", "kp", "Dsn", "Dsp", "brugg_n", "brugg_p", "brugg_s")
PARAMETER_KEYS = {
    "kn": "Negative electrode exchange-current density [A.m-2]",
    "kp": "Positive electrode exchange-current density [A.m-2]",
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
EFFECT_FRACTION = 0.05
LOG_EFFECT_STEP = float(np.log1p(EFFECT_FRACTION))
SENSITIVITY_FLOOR = 0.10
COSINE_THRESHOLD = 0.95
N_POINTS_PER_CRATE_CONDITION = 100


def build_current_inputs():
    bundle = ga.load_legacy_bundle(ga.get_paths(ROOT))
    base = audit.build_base(bundle, 1e-4, audit.CSN_CORRECTED)
    _, _, hybrid, cathode = rs.candidate_details(bundle, base)
    experiment = omc.load_old_dynamic_data(base["dynamic"], base["q_meas"])
    candidates = pd.read_csv(SOURCE / "mass_independent_ocp_candidate_metrics.csv")
    ocp_row = candidates[candidates.Candidate == current.NAME].iloc[0]
    stage = {key: float(ocp_row[key]) for key in ("x0", "x100", "y100", "y0")}
    effective = mit.model_for_effective_capacities(base, ocp_row)
    model = audit.set_dynamic_parameters(effective, current.RN_UM, current.RP_UM, 1e-4)
    model = current.apply_provisional_thickness_preserve_capacity(model)
    model = current.apply_current_kinetics_and_bruggeman(model)
    anode = common.scaled_detail(hybrid, 0.50)
    cathode = common.scaled_detail(cathode, 0.50)
    return bundle, model, stage, anode, cathode, experiment, ocp_row


def perturb_model(model_inputs, parameter, direction):
    changed = dict(model_inputs)
    changed["params"] = model_inputs["params"].copy()
    key = PARAMETER_KEYS[parameter]
    nominal = model_inputs["params"][key]
    if parameter in POROSITY_KEYS:
        porosity = float(model_inputs["params"][POROSITY_KEYS[parameter]])
        # Coordinate z = log(eps**b / eps**b0).  A step in z therefore has
        # the same 5% effective-transport scale as the k and D perturbations.
        value = float(nominal) + direction * LOG_EFFECT_STEP / np.log(porosity)
    else:
        value = ga.scale_parameter_function(nominal, np.exp(direction * LOG_EFFECT_STEP)) if callable(nominal) else float(nominal) * np.exp(direction * LOG_EFFECT_STEP)
    changed["params"].update({key: value}, check_already_exists=False)
    return changed


def c_rate_grids(experiment, baseline_runs):
    """Build a fixed pre-cutoff comparison grid for every C-rate condition.

    The 2C charge simulation reaches its voltage cutoff before the experiment
    reaches 70% SOC.  Use the common portion of the requested experimental
    10--70% SOC window and the nominal model trajectory.  The separate
    endpoint analysis below retains the cutoff/capacity information.
    """
    grids = {}
    for rate in ga.RATES:
        for charge in (True, False):
            exp = experiment[(rate, charge)]
            mask = (exp["SOC"] >= 0.10) & (exp["SOC"] <= 0.70)
            lower = max(
                float(np.min(exp["t_min"][mask])),
                float(np.min(baseline_runs[(rate, charge)]["t_min"])),
            )
            # Keep a 10% guard band from the nominal cutoff so the same grid
            # remains valid under the local +/-5% perturbations.
            upper = min(
                float(np.max(exp["t_min"][mask])),
                0.90 * float(np.max(baseline_runs[(rate, charge)]["t_min"])),
            )
            if upper <= lower:
                raise RuntimeError(
                    f"No common pre-cutoff C-rate window: {rate:g}C "
                    f"{'charge' if charge else 'discharge'}"
                )
            grids[(rate, charge)] = np.linspace(
                lower,
                upper,
                N_POINTS_PER_CRATE_CONDITION,
            )
    return grids


def run_c_rate_set(model_inputs, stage, anode, cathode):
    runs = {}
    for rate in ga.RATES:
        for charge in (True, False):
            print(f"  C-rate {rate:g}C {'charge' if charge else 'discharge'}", flush=True)
            runs[(rate, charge)] = ehq.run_dfn(
                rate, charge, stage, model_inputs, anode, cathode, True, True
            )
    return runs


def vectorize_c_rate(runs, grids, direction_filter=None):
    vectors = []
    condition_vectors = {}
    for rate in ga.RATES:
        for charge in (True, False):
            if direction_filter is not None and charge != direction_filter:
                continue
            sim = runs[(rate, charge)]
            grid = grids[(rate, charge)]
            if grid[-1] > sim["t_min"][-1] + 1e-8:
                raise RuntimeError(
                    f"Simulation ended before the 10-70% grid: {rate:g}C "
                    f"{'charge' if charge else 'discharge'}"
                )
            values = np.interp(grid, sim["t_min"], sim["V"])
            condition_vectors[(rate, charge)] = values
            vectors.append(values)
    return np.concatenate(vectors), condition_vectors


def load_hppc_targets():
    if not HPPC_FILE.exists():
        raise FileNotFoundError(HPPC_FILE)
    step = pd.read_excel(HPPC_FILE, sheet_name="step")
    record = pd.read_excel(HPPC_FILE, sheet_name="record")
    for col in ("Oneset Date", "End Date"):
        step[col] = pd.to_datetime(step[col])
    record["Date"] = pd.to_datetime(record["Date"])
    step["Cycle Index"] = pd.to_numeric(step["Cycle Index"], errors="coerce")
    step["Step Index"] = pd.to_numeric(step["Step Index"], errors="coerce")

    def get_step(cycle, step_index):
        selected = step[(step["Cycle Index"] == cycle) & (step["Step Index"] == step_index)]
        if len(selected) != 1:
            raise RuntimeError(f"HPPC step not unique: cycle={cycle}, step={step_index}")
        return selected.iloc[0]

    def segment(row):
        selected = record[
            (record["Date"] >= row["Oneset Date"])
            & (record["Date"] <= row["End Date"])
            & (record["Step Type"].astype(str) == str(row["Step Type"]))
        ].sort_values("Date")
        if selected.empty:
            raise RuntimeError(f"Missing HPPC record segment: {row['Oneset Date']}")
        return selected

    def closest_voltage(frame, time):
        idx = (frame["Date"] - time).abs().idxmin()
        return float(frame.loc[idx, "Voltage(V)"])

    rows = []
    for cycle in range(1, 11):
        rest_pre = get_step(cycle, 24)
        dchg = get_step(cycle, 25)
        rest = get_step(cycle, 26)
        chg = get_step(cycle, 27)
        rest_pre_data = segment(rest_pre)
        dchg_data = segment(dchg)
        rest_data = segment(rest)
        chg_data = segment(chg)
        values = {
            "rest_pre_end": float(rest_pre_data["Voltage(V)"].iloc[-1]),
            "dchg_end": float(dchg_data["Voltage(V)"].iloc[-1]),
            "rest_30s": closest_voltage(rest_data, rest["Oneset Date"] + pd.Timedelta(seconds=30)),
            "rest_end": float(rest_data["Voltage(V)"].iloc[-1]),
            "chg_end": float(chg_data["Voltage(V)"].iloc[-1]),
        }
        for point, voltage in values.items():
            rows.append(
                {
                    "Cycle": cycle,
                    "Nominal_SOC": 1.0 - 0.1 * (cycle - 1),
                    "Point": point,
                    "Voltage_V": voltage,
                }
            )
    targets = pd.DataFrame(rows)
    order = ["rest_pre_end", "dchg_end", "rest_30s", "rest_end", "chg_end"]
    targets["Point"] = pd.Categorical(targets.Point, categories=order, ordered=True)
    return targets.sort_values(["Cycle", "Point"]).reset_index(drop=True)


def make_hppc_experiment_and_times():
    steps = []
    rows = []
    time_s = 0.0
    for cycle in range(1, 11):
        steps.extend(
            [
                "Rest for 10 minutes",
                "Discharge at 2.28 A for 30 seconds",
                "Rest for 40 seconds",
                "Charge at 1.71 A for 10 seconds",
            ]
        )
        time_s += 600.0
        rows.append({"Cycle": cycle, "Point": "rest_pre_end", "Time_s": time_s})
        time_s += 30.0
        rows.append({"Cycle": cycle, "Point": "dchg_end", "Time_s": time_s})
        time_s += 30.0
        rows.append({"Cycle": cycle, "Point": "rest_30s", "Time_s": time_s})
        time_s += 10.0
        rows.append({"Cycle": cycle, "Point": "rest_end", "Time_s": time_s})
        time_s += 10.0
        rows.append({"Cycle": cycle, "Point": "chg_end", "Time_s": time_s})
        if cycle < 10:
            steps.append("Discharge at 0.76 A for 18 minutes")
            time_s += 1080.0
    return pybamm.Experiment(steps, period="10 seconds"), pd.DataFrame(rows)


def interp_first_duplicate(time, values, targets):
    order = np.argsort(time, kind="stable")
    unique_time, index = np.unique(np.asarray(time)[order], return_index=True)
    unique_values = np.asarray(values)[order][index]
    return np.interp(np.asarray(targets, dtype=float), unique_time, unique_values)


def run_hppc(model_inputs, stage, anode, cathode, experiment, target_times):
    model = pybamm.lithium_ion.DFN(
        {"open-circuit potential": ("one-state hysteresis", "one-state hysteresis"), "thermal": "isothermal"}
    )
    values = ehq.make_parameter_values(
        model_inputs, stage, anode, cathode, False, True, True
    )
    simulation = pybamm.Simulation(
        model,
        parameter_values=values,
        experiment=experiment,
        var_pts={"x_n": 15, "x_s": 15, "x_p": 15, "r_n": 15, "r_p": 15},
        solver=pybamm.IDAKLUSolver(rtol=1e-6, atol=1e-8),
    )
    solution = simulation.solve()
    time_s = np.asarray(solution["Time [s]"].entries, dtype=float)
    voltage = np.asarray(solution["Terminal voltage [V]"].entries, dtype=float)
    return interp_first_duplicate(time_s, voltage, target_times.Time_s)


def summarize_matrix(dataset, matrix):
    rms = np.sqrt(np.mean(matrix**2, axis=0)) * 1000.0
    relative = rms / np.max(rms)
    norms = np.linalg.norm(matrix, axis=0)
    normalized = matrix / np.where(norms > 0, norms, 1.0)
    cosine = normalized.T @ normalized
    singular = np.linalg.svd(normalized, compute_uv=False)
    rel_singular = singular / singular[0]
    condition = np.inf if singular[-1] == 0 else singular[0] / singular[-1]

    order = np.argsort(-rms)
    selected_indices = []
    for idx in order:
        if relative[idx] < SENSITIVITY_FLOOR:
            continue
        if selected_indices and max(abs(cosine[idx, j]) for j in selected_indices) >= COSINE_THRESHOLD:
            continue
        selected_indices.append(idx)
    selected = [PARAMETERS[idx] for idx in selected_indices]
    ranking = pd.DataFrame(
        {
            "Dataset": dataset,
            "Parameter": PARAMETERS,
            "RMS_mV_per_log_effect": rms,
            "Relative_sensitivity": relative,
            "Selected": [p in selected for p in PARAMETERS],
        }
    ).sort_values("RMS_mV_per_log_effect", ascending=False)
    correlation = pd.DataFrame(cosine, index=PARAMETERS, columns=PARAMETERS)
    singular_table = pd.DataFrame(
        {
            "Dataset": dataset,
            "Index": np.arange(1, len(singular) + 1),
            "Singular_value": singular,
            "Relative_singular_value": rel_singular,
            "Condition_number": condition,
            "Rank_relative_1e-2": int(np.sum(rel_singular >= 1e-2)),
            "Rank_relative_1e-3": int(np.sum(rel_singular >= 1e-3)),
        }
    )
    return ranking, correlation, singular_table, selected


def main():
    bundle, baseline, stage, anode, cathode, experiment, ocp_row = build_current_inputs()
    hppc_targets = load_hppc_targets()
    hppc_experiment, hppc_times = make_hppc_experiment_and_times()
    hppc_targets.to_csv(RESULTS / "hppc_experimental_targets.csv", index=False, encoding="utf-8-sig")

    print("Baseline", flush=True)
    base_c_runs = run_c_rate_set(baseline, stage, anode, cathode)
    grids = c_rate_grids(experiment, base_c_runs)
    base_vectors = {
        "C-rate combined": vectorize_c_rate(base_c_runs, grids)[0],
        "C-rate charge": vectorize_c_rate(base_c_runs, grids, True)[0],
        "C-rate discharge": vectorize_c_rate(base_c_runs, grids, False)[0],
    }
    print("  HPPC", flush=True)
    base_hppc = run_hppc(baseline, stage, anode, cathode, hppc_experiment, hppc_times)

    c_columns = {key: [] for key in base_vectors}
    hppc_columns = []
    condition_rows = []
    endpoint_rows = []
    linearity_rows = []

    for parameter in PARAMETERS:
        print(f"Parameter {parameter}", flush=True)
        plus = perturb_model(baseline, parameter, +1)
        minus = perturb_model(baseline, parameter, -1)
        plus_runs = run_c_rate_set(plus, stage, anode, cathode)
        minus_runs = run_c_rate_set(minus, stage, anode, cathode)
        print("  HPPC plus/minus", flush=True)
        plus_hppc = run_hppc(plus, stage, anode, cathode, hppc_experiment, hppc_times)
        minus_hppc = run_hppc(minus, stage, anode, cathode, hppc_experiment, hppc_times)

        for dataset, direction_filter in (
            ("C-rate combined", None),
            ("C-rate charge", True),
            ("C-rate discharge", False),
        ):
            plus_vector = vectorize_c_rate(plus_runs, grids, direction_filter)[0]
            minus_vector = vectorize_c_rate(minus_runs, grids, direction_filter)[0]
            sensitivity = (plus_vector - minus_vector) / (2.0 * LOG_EFFECT_STEP)
            c_columns[dataset].append(sensitivity)
            numerator = np.linalg.norm(plus_vector + minus_vector - 2.0 * base_vectors[dataset])
            denominator = np.linalg.norm(plus_vector - minus_vector)
            linearity_rows.append(
                {
                    "Dataset": dataset,
                    "Parameter": parameter,
                    "Central_nonlinearity_ratio": float(numerator / denominator) if denominator else np.nan,
                }
            )

        hppc_sensitivity = (plus_hppc - minus_hppc) / (2.0 * LOG_EFFECT_STEP)
        hppc_columns.append(hppc_sensitivity)
        numerator = np.linalg.norm(plus_hppc + minus_hppc - 2.0 * base_hppc)
        denominator = np.linalg.norm(plus_hppc - minus_hppc)
        linearity_rows.append(
            {
                "Dataset": "HPPC",
                "Parameter": parameter,
                "Central_nonlinearity_ratio": float(numerator / denominator) if denominator else np.nan,
            }
        )

        for rate in ga.RATES:
            for charge in (True, False):
                pv = vectorize_c_rate(plus_runs, grids)[1][(rate, charge)]
                mv = vectorize_c_rate(minus_runs, grids)[1][(rate, charge)]
                sens = (pv - mv) / (2.0 * LOG_EFFECT_STEP)
                condition_rows.append(
                    {
                        "Dataset": "C-rate",
                        "C_rate": rate,
                        "Direction": "Charge" if charge else "Discharge",
                        "Parameter": parameter,
                        "RMS_mV_per_log_effect": float(np.sqrt(np.mean(sens**2)) * 1000.0),
                        "MaxAbs_mV_per_log_effect": float(np.max(np.abs(sens)) * 1000.0),
                    }
                )
                end_time_sensitivity = (
                    plus_runs[(rate, charge)]["t_min"][-1]
                    - minus_runs[(rate, charge)]["t_min"][-1]
                ) / (2.0 * LOG_EFFECT_STEP)
                endpoint_rows.append(
                    {
                        "C_rate": rate,
                        "Direction": "Charge" if charge else "Discharge",
                        "Parameter": parameter,
                        "End_time_min_per_log_effect": float(end_time_sensitivity),
                        "Capacity_mAh_per_log_effect": float(end_time_sensitivity * rate * 2.28 / 60.0 * 1000.0),
                    }
                )

        for point in hppc_targets.Point.astype(str).unique():
            mask = (hppc_targets.Point.astype(str) == point).to_numpy()
            sens = hppc_sensitivity[mask]
            condition_rows.append(
                {
                    "Dataset": "HPPC",
                    "C_rate": np.nan,
                    "Direction": point,
                    "Parameter": parameter,
                    "RMS_mV_per_log_effect": float(np.sqrt(np.mean(sens**2)) * 1000.0),
                    "MaxAbs_mV_per_log_effect": float(np.max(np.abs(sens)) * 1000.0),
                }
            )

    matrices = {name: np.column_stack(columns) for name, columns in c_columns.items()}
    matrices["HPPC"] = np.column_stack(hppc_columns)
    rankings = []
    correlations = {}
    singular_tables = []
    selected_subsets = {}
    for dataset, matrix in matrices.items():
        ranking, correlation, singular, selected = summarize_matrix(dataset, matrix)
        rankings.append(ranking)
        correlations[dataset] = correlation
        singular_tables.append(singular)
        selected_subsets[dataset] = selected
        safe_name = dataset.lower().replace("-", "_").replace(" ", "_")
        correlation.to_csv(RESULTS / f"{safe_name}_cosine_correlation.csv", encoding="utf-8-sig")
        np.save(RESULTS / f"{safe_name}_sensitivity_matrix.npy", matrix)

    ranking_all = pd.concat(rankings, ignore_index=True)
    singular_all = pd.concat(singular_tables, ignore_index=True)
    condition_detail = pd.DataFrame(condition_rows)
    endpoint_detail = pd.DataFrame(endpoint_rows)
    linearity = pd.DataFrame(linearity_rows)
    ranking_all.to_csv(RESULTS / "sensitivity_ranking.csv", index=False, encoding="utf-8-sig")
    singular_all.to_csv(RESULTS / "identifiability_singular_values.csv", index=False, encoding="utf-8-sig")
    condition_detail.to_csv(RESULTS / "condition_sensitivity_detail.csv", index=False, encoding="utf-8-sig")
    endpoint_detail.to_csv(RESULTS / "endpoint_sensitivity.csv", index=False, encoding="utf-8-sig")
    linearity.to_csv(RESULTS / "central_linearity_check.csv", index=False, encoding="utf-8-sig")

    datasets_to_plot = ("C-rate combined", "HPPC")
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), constrained_layout=True)
    for ax, dataset in zip(axes, datasets_to_plot):
        frame = ranking_all[ranking_all.Dataset == dataset].sort_values("RMS_mV_per_log_effect", ascending=True)
        colors = ["#0072B2" if value else "#999999" for value in frame.Selected]
        ax.barh(frame.Parameter, frame.RMS_mV_per_log_effect, color=colors)
        ax.set_xlabel("RMS sensitivity (mV per log-effect)")
        ax.set_title(dataset)
        ax.grid(axis="x", alpha=0.2)
    fig.suptitle("Seven-parameter sensitivity at current configuration (5% effect step)")
    fig.savefig(RESULTS / "sensitivity_ranking.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2), constrained_layout=True)
    for ax, dataset in zip(axes, datasets_to_plot):
        corr = correlations[dataset]
        image = ax.imshow(corr.to_numpy(), vmin=-1, vmax=1, cmap="coolwarm")
        ax.set_xticks(range(len(PARAMETERS)), PARAMETERS, rotation=35, ha="right")
        ax.set_yticks(range(len(PARAMETERS)), PARAMETERS)
        ax.set_title(f"{dataset} cosine correlation")
        for i in range(len(PARAMETERS)):
            for j in range(len(PARAMETERS)):
                ax.text(j, i, f"{corr.iloc[i, j]:.2f}", ha="center", va="center", fontsize=8)
        fig.colorbar(image, ax=ax, shrink=0.8)
    fig.savefig(RESULTS / "sensitivity_cosine_correlation.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2), constrained_layout=True)
    crate_pivot = condition_detail[condition_detail.Dataset == "C-rate"].pivot_table(
        index=["C_rate", "Direction"], columns="Parameter", values="RMS_mV_per_log_effect"
    )
    hppc_pivot = condition_detail[condition_detail.Dataset == "HPPC"].pivot_table(
        index="Direction", columns="Parameter", values="RMS_mV_per_log_effect"
    )
    for ax, pivot, title in (
        (axes[0], crate_pivot, "C-rate condition sensitivity"),
        (axes[1], hppc_pivot, "HPPC point-type sensitivity"),
    ):
        pivot = pivot.loc[:, list(PARAMETERS)]
        image = ax.imshow(pivot.to_numpy(), aspect="auto", cmap="viridis")
        ax.set_xticks(range(len(PARAMETERS)), PARAMETERS, rotation=35, ha="right")
        ax.set_yticks(range(len(pivot.index)), [" / ".join(map(str, idx)) if isinstance(idx, tuple) else str(idx) for idx in pivot.index])
        ax.set_title(title)
        for i in range(len(pivot.index)):
            for j in range(len(PARAMETERS)):
                ax.text(j, i, f"{pivot.iloc[i, j]:.1f}", ha="center", va="center", fontsize=7, color="white" if pivot.iloc[i, j] > np.nanmedian(pivot.to_numpy()) else "black")
        fig.colorbar(image, ax=ax, shrink=0.8, label="mV per log-effect")
    fig.savefig(RESULTS / "condition_sensitivity_heatmaps.png", dpi=220)
    plt.close(fig)

    report = {
        "baseline": {
            "kn": current.KN_PREF,
            "kp": current.KP_PREF,
            "Dsn_m2_s": rs.DSN,
            "Dsp_m2_s": rs.DSP,
            "brugg_n": current.BRUGG_N,
            "brugg_p": current.BRUGG_P,
            "brugg_s": current.BRUGG_S,
            "Ln_um": current.LN_UM,
            "Lp_um": current.LP_UM,
            "Rn_um": current.RN_UM,
            "Rp_um": current.RP_UM,
            "qOCV_MAE_mV": float(ocp_row.qOCV_MAE_2_98_mV),
        },
        "method": {
            "effect_fraction": EFFECT_FRACTION,
            "log_effect_step": LOG_EFFECT_STEP,
            "positive_parameters": "central derivative with respect to log(parameter)",
            "Bruggeman_coefficients": "central derivatives with respect to log(eps_i**b_i / eps_i**b_i0), i=n,p,s",
            "C_rate_voltage_window": "common pre-cutoff portion of experimental 10-70% SOC window; upper bound=min(experimental 70% SOC time, 90% nominal-model cutoff time); 100 points per condition",
            "HPPC_points": int(len(hppc_targets)),
            "sensitivity_floor": SENSITIVITY_FLOOR,
            "cosine_threshold": COSINE_THRESHOLD,
        },
        "selected_subsets": selected_subsets,
        "rankings": ranking_all.to_dict(orient="records"),
        "singular_values": singular_all.to_dict(orient="records"),
    }
    (RESULTS / "stage2_sensitivity_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    lines = [
        "# Stage 2 seven-parameter sensitivity and identifiability",
        "",
        "Current OCP/stoichiometry/geometry were fixed. The candidate set was `kn, kp, Dsn, Dsp, brugg_n, brugg_p, brugg_s`.",
        "Positive parameters were perturbed by ±5% in log-effect coordinates. Each Bruggeman coefficient was perturbed to create the same ±5% change in its local `eps_i**b_i`, making RMS effects comparable across all candidates.",
        "",
        "## Selected subsets",
        "",
    ]
    for dataset, selected in selected_subsets.items():
        lines.append(f"- {dataset}: `{', '.join(selected) if selected else 'none'}`")
    lines.extend(["", "## Ranking", ""])
    for dataset in datasets_to_plot:
        lines.extend([f"### {dataset}", "", "| Parameter | RMS (mV/log-effect) | Relative | Selected |", "|---|---:|---:|---|"])
        frame = ranking_all[ranking_all.Dataset == dataset]
        for row in frame.itertuples():
            lines.append(f"| {row.Parameter} | {row.RMS_mV_per_log_effect:.3f} | {row.Relative_sensitivity:.3f} | {'yes' if row.Selected else 'no'} |")
        lines.append("")
    lines.extend(
        [
            "## Interpretation guardrails",
            "",
            "- Selection is local to the current parameter point and the 5% effect scale.",
            "- HPPC raw pulse records contain only pulse start/end points; diffusion identifiability must therefore be interpreted conservatively.",
            "- A parameter excluded from one dataset can still be retained in the other dataset's independent fitting path.",
            "- Stage 3 should fit only the subset selected from its own training dataset to prevent validation leakage.",
        ]
    )
    (RESULTS / "STAGE2_SENSITIVITY_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\nSelected subsets")
    print(json.dumps(selected_subsets, indent=2))
    print("\nRanking")
    print(ranking_all[ranking_all.Dataset.isin(datasets_to_plot)].to_string(index=False))
    print("\nSingular values")
    print(singular_all[singular_all.Dataset.isin(datasets_to_plot)].to_string(index=False))


if __name__ == "__main__":
    main()
