"""Dynamic + HPPC sensitivity and identifiability for four DFN parameters."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybamm

import comprehensive_prefit_cross_cohort as cross
import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import kinetics_original_reversion_ablation as kinetics
import priority_initial_state_protocol_recheck as priority

import analyze_oregan_brugg_2p5_anode_potential as anode_analysis
import electrolyte_transport_set_comparison as comparison


OUT = comparison.ROOT / "results" / "261012_oregan_combined_sensitivity_identifiability"
BRUGG_N = 2.45
PARAMETERS = ("Dsn", "Dsp", "kn", "kp")
INITIAL = {
    "Dsn": 4.1707741459682473e-14,
    "Dsp": 6.888175323807748e-14,
    "kn": 9.648533212e-7,
    "kp": 5.899519906708763e-7,
}
HPPC_FILE = Path(
    r"C:\Users\user\OneDrive\01. 적응형 충전 프로토콜\4. Enertech 파우치셀 실험데이터"
    r"\260808 Enertech셀 열화데이터\260808 Eneterch셀 열화데이터 v1 6-6.xlsx"
)
HPPC_POINT_ORDER = ("rest_pre_end", "dchg_end", "rest_30s", "rest_end", "chg_end")
HPPC_DCHG_CURRENT_A = 2.28
HPPC_CHG_CURRENT_A = 1.71
HPPC_SHIFT_CURRENT_A = 0.76
SENSITIVITY_STEP = np.log(1.05)


def apply_parameters(base: dict, values: dict) -> dict:
    changed = kinetics.with_prefactors(base, values["kn"], values["kp"])
    changed = dict(changed)
    changed["params"] = changed["params"].copy()
    changed["params"].update(
        {
            "Negative particle diffusivity [m2.s-1]": values["Dsn"],
            "Positive particle diffusivity [m2.s-1]": values["Dsp"],
            "Negative electrode Bruggeman coefficient (electrolyte)": BRUGG_N,
        },
        check_already_exists=False,
    )
    return changed


def load_hppc_targets():
    if not HPPC_FILE.exists():
        raise FileNotFoundError(HPPC_FILE)
    step = pd.read_excel(HPPC_FILE, sheet_name="step")
    record = pd.read_excel(HPPC_FILE, sheet_name="record")
    for column in ("Oneset Date", "End Date"):
        step[column] = pd.to_datetime(step[column])
    record["Date"] = pd.to_datetime(record["Date"])
    step["Cycle Index"] = pd.to_numeric(step["Cycle Index"], errors="coerce")
    step["Step Index"] = pd.to_numeric(step["Step Index"], errors="coerce")

    def get_step(cycle, step_index):
        selected = step[
            (step["Cycle Index"] == cycle) & (step["Step Index"] == step_index)
        ]
        if len(selected) != 1:
            raise RuntimeError(
                "HPPC step not unique: cycle={}, step={}".format(cycle, step_index)
            )
        return selected.iloc[0]

    def segment(step_row):
        selected = record[
            (record["Date"] >= step_row["Oneset Date"])
            & (record["Date"] <= step_row["End Date"])
            & (record["Step Type"].astype(str) == str(step_row["Step Type"]))
        ].copy()
        if selected.empty:
            raise RuntimeError("Empty HPPC record segment")
        selected = selected.sort_values("Date")
        selected["Voltage(V)"] = pd.to_numeric(
            selected["Voltage(V)"], errors="coerce"
        )
        return selected

    def closest_voltage(data, target):
        index = (data["Date"] - target).abs().idxmin()
        return float(data.loc[index, "Voltage(V)"])

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
            "rest_30s": closest_voltage(
                rest_data, rest["Oneset Date"] + pd.Timedelta(seconds=30)
            ),
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
    targets["Point"] = pd.Categorical(
        targets["Point"], categories=HPPC_POINT_ORDER, ordered=True
    )
    targets = targets.sort_values(["Cycle", "Point"]).reset_index(drop=True)
    return targets


def make_hppc_experiment():
    steps = []
    for cycle in range(1, 11):
        steps.extend(
            [
                "Rest for 10 minutes",
                "Discharge at {} A for 30 seconds".format(HPPC_DCHG_CURRENT_A),
                "Rest for 40 seconds",
                "Charge at {} A for 10 seconds".format(HPPC_CHG_CURRENT_A),
            ]
        )
        if cycle < 10:
            steps.append(
                "Discharge at {} A for 18 minutes".format(HPPC_SHIFT_CURRENT_A)
            )
    return pybamm.Experiment(steps, period="10 seconds")


def hppc_target_times():
    rows = []
    time_s = 0.0
    for cycle in range(1, 11):
        time_s += 600.0
        rows.append((cycle, "rest_pre_end", time_s))
        time_s += 30.0
        rows.append((cycle, "dchg_end", time_s))
        time_s += 30.0
        rows.append((cycle, "rest_30s", time_s))
        time_s += 10.0
        rows.append((cycle, "rest_end", time_s))
        time_s += 10.0
        rows.append((cycle, "chg_end", time_s))
        if cycle < 10:
            time_s += 1080.0
    return pd.DataFrame(rows, columns=["Cycle", "Point", "Time_s"])


def interp_first_duplicate(time, values, targets):
    time = np.asarray(time, dtype=float)
    values = np.asarray(values, dtype=float)
    order = np.argsort(time, kind="stable")
    unique_time, index = np.unique(time[order], return_index=True)
    unique_values = values[order][index]
    return np.interp(np.asarray(targets, dtype=float), unique_time, unique_values)


def run_hppc(model_inputs, initial_stage, anode, cathode, target_times):
    model = pybamm.lithium_ion.DFN(
        {"open-circuit potential": ("single", "single"), "thermal": "isothermal"}
    )
    simulation = pybamm.Simulation(
        model,
        parameter_values=ehq.make_parameter_values(
            model_inputs, initial_stage, anode, cathode, False, False, False
        ),
        experiment=make_hppc_experiment(),
        var_pts={"x_n": 15, "x_s": 15, "x_p": 15, "r_n": 15, "r_p": 15},
        solver=pybamm.IDAKLUSolver(rtol=1e-6, atol=1e-8),
    )
    solution = simulation.solve()
    return interp_first_duplicate(
        solution["Time [s]"].entries,
        solution["Terminal voltage [V]"].entries,
        target_times,
    )


def run_dynamic(model, stage, anode, cathode):
    runs = {}
    for rate in ga.RATES:
        for charge in (True, False):
            rest_v, current = comparison.PROTOCOL[(rate, charge)]
            z0 = priority.monotone_voltage_inverse(model["v_qocv"], model["soc"], rest_v)
            local_stage = priority.stage_at_soc(stage, z0, charge)
            runs[(rate, charge)] = anode_analysis.run_branch(
                current / cross.NOMINAL_CAPACITY_AH,
                charge,
                local_stage,
                model,
                anode,
                cathode,
            )
    return runs


def make_dynamic_grids(experiment, baseline_runs, qcell, n_points=80):
    grids = {}
    for rate in ga.RATES:
        for charge in (True, False):
            exp = experiment[(rate, charge)]
            q_exp = priority.transferred_capacity(exp, qcell, charge)
            q_mod = priority.transferred_capacity(baseline_runs[(rate, charge)], qcell, charge)
            lo = 0.10 * qcell
            hi = min(0.70 * qcell, 0.98 * float(q_exp[-1]), 0.98 * float(q_mod[-1]))
            if hi <= lo:
                raise RuntimeError("Insufficient dynamic overlap")
            grids[(rate, charge)] = np.linspace(lo, hi, n_points)
    return grids


def dynamic_vector(data, grids, qcell, experiment=False):
    blocks = []
    for rate in ga.RATES:
        for charge in (True, False):
            source = data[(rate, charge)]
            q = priority.transferred_capacity(source, qcell, charge)
            order = np.argsort(q, kind="stable")
            q_unique, index = np.unique(np.asarray(q)[order], return_index=True)
            voltage = np.asarray(source["V"], dtype=float)[order][index]
            blocks.append(np.interp(grids[(rate, charge)], q_unique, voltage))
    return np.concatenate(blocks)


def matrix_report(name, matrix):
    rms = np.sqrt(np.mean(matrix**2, axis=0))
    unit = matrix / np.maximum(np.linalg.norm(matrix, axis=0), np.finfo(float).eps)
    cosine = unit.T @ unit
    singular = np.linalg.svd(unit, compute_uv=False)
    condition = float(singular[0] / singular[-1])
    pd.DataFrame(cosine, index=PARAMETERS, columns=PARAMETERS).to_csv(
        OUT / "{}_cosine.csv".format(name), encoding="utf-8-sig"
    )
    pd.DataFrame(
        {
            "Index": np.arange(1, len(singular) + 1),
            "Singular_value": singular,
            "Relative_sigma": singular / singular[0],
        }
    ).to_csv(OUT / "{}_singular_values.csv".format(name), index=False, encoding="utf-8-sig")
    return {
        "Dataset": name,
        "Condition_number": condition,
        "Max_abs_pairwise_cosine": float(
            np.max(np.abs(cosine - np.eye(len(PARAMETERS))))
        ),
        "Rank_relative_1e-2": int(np.sum(singular / singular[0] >= 1.0e-2)),
        "Rank_relative_1e-3": int(np.sum(singular / singular[0] >= 1.0e-3)),
        **{"RMS_{}_mV_per_ln".format(p): float(v) for p, v in zip(PARAMETERS, rms)},
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    hppc_targets = load_hppc_targets()
    hppc_times = hppc_target_times()
    hppc_exp = hppc_targets["Voltage_V"].to_numpy(float)

    baseline, stage, anode, cathode, experiment, qcell = comparison.current_p0()
    base_transport = comparison.with_transport_set(baseline, "ORegan2022")
    base_model = apply_parameters(base_transport, INITIAL)
    first_rest_voltage = float(hppc_targets.iloc[0]["Voltage_V"])
    hppc_z0 = priority.monotone_voltage_inverse(
        base_model["v_qocv"], base_model["soc"], first_rest_voltage
    )
    hppc_stage = priority.stage_at_soc(stage, hppc_z0, False)

    print("Baseline dynamic", flush=True)
    baseline_dynamic_runs = run_dynamic(base_model, stage, anode, cathode)
    dynamic_grids = make_dynamic_grids(experiment, baseline_dynamic_runs, qcell)
    dynamic_exp = dynamic_vector(experiment, dynamic_grids, qcell, experiment=True)
    dynamic_base = dynamic_vector(baseline_dynamic_runs, dynamic_grids, qcell)
    print("Baseline HPPC", flush=True)
    hppc_base = run_hppc(base_model, hppc_stage, anode, cathode, hppc_times["Time_s"])

    dynamic_columns = []
    hppc_columns = []
    ranking_rows = []
    for parameter in PARAMETERS:
        outputs = {}
        for sign, label in ((-1.0, "minus"), (1.0, "plus")):
            values = dict(INITIAL)
            values[parameter] *= np.exp(sign * SENSITIVITY_STEP)
            model = apply_parameters(base_transport, values)
            print("Sensitivity {} {} dynamic".format(parameter, label), flush=True)
            dynamic_runs = run_dynamic(model, stage, anode, cathode)
            outputs[(label, "dynamic")] = dynamic_vector(
                dynamic_runs, dynamic_grids, qcell
            )
            print("Sensitivity {} {} HPPC".format(parameter, label), flush=True)
            outputs[(label, "hppc")] = run_hppc(
                model, hppc_stage, anode, cathode, hppc_times["Time_s"]
            )
        dynamic_derivative = 1000.0 * (
            outputs[("plus", "dynamic")] - outputs[("minus", "dynamic")]
        ) / (2.0 * SENSITIVITY_STEP)
        hppc_derivative = 1000.0 * (
            outputs[("plus", "hppc")] - outputs[("minus", "hppc")]
        ) / (2.0 * SENSITIVITY_STEP)
        dynamic_columns.append(dynamic_derivative)
        hppc_columns.append(hppc_derivative)
        ranking_rows.extend(
            [
                {
                    "Dataset": "dynamic",
                    "Parameter": parameter,
                    "Initial": INITIAL[parameter],
                    "RMS_mV_per_ln_parameter": float(
                        np.sqrt(np.mean(dynamic_derivative**2))
                    ),
                    "Max_abs_mV_per_ln_parameter": float(
                        np.max(np.abs(dynamic_derivative))
                    ),
                },
                {
                    "Dataset": "hppc",
                    "Parameter": parameter,
                    "Initial": INITIAL[parameter],
                    "RMS_mV_per_ln_parameter": float(
                        np.sqrt(np.mean(hppc_derivative**2))
                    ),
                    "Max_abs_mV_per_ln_parameter": float(
                        np.max(np.abs(hppc_derivative))
                    ),
                },
            ]
        )

    dynamic_matrix = np.column_stack(dynamic_columns)
    hppc_matrix = np.column_stack(hppc_columns)
    combined_matrix = np.vstack(
        [
            dynamic_matrix / np.sqrt(dynamic_matrix.shape[0]),
            hppc_matrix / np.sqrt(hppc_matrix.shape[0]),
        ]
    )
    np.save(OUT / "dynamic_sensitivity_matrix.npy", dynamic_matrix)
    np.save(OUT / "hppc_sensitivity_matrix.npy", hppc_matrix)
    np.save(OUT / "combined_equal_dataset_sensitivity_matrix.npy", combined_matrix)
    ranking = pd.DataFrame(ranking_rows)
    ranking.to_csv(OUT / "sensitivity_ranking.csv", index=False, encoding="utf-8-sig")
    reports = pd.DataFrame(
        [
            matrix_report("dynamic", dynamic_matrix),
            matrix_report("hppc", hppc_matrix),
            matrix_report("combined_equal_dataset", combined_matrix),
        ]
    )
    reports.to_csv(OUT / "identifiability_summary.csv", index=False, encoding="utf-8-sig")

    baseline_metrics = {
        "Dynamic_RMSE_mV": float(
            np.sqrt(np.mean((1000.0 * (dynamic_base - dynamic_exp)) ** 2))
        ),
        "HPPC_RMSE_mV": float(
            np.sqrt(np.mean((1000.0 * (hppc_base - hppc_exp)) ** 2))
        ),
        "HPPC_initial_SOC_pct": 100.0 * hppc_z0,
        "HPPC_points": len(hppc_exp),
        "Dynamic_points": len(dynamic_exp),
    }
    pd.DataFrame([baseline_metrics]).to_csv(
        OUT / "baseline_objective_metrics.csv", index=False, encoding="utf-8-sig"
    )
    hppc_compare = hppc_targets.copy()
    hppc_compare["Model_V"] = hppc_base
    hppc_compare["Error_mV"] = 1000.0 * (hppc_base - hppc_exp)
    hppc_compare.to_csv(OUT / "hppc_baseline_keypoints.csv", index=False, encoding="utf-8-sig")

    fig, axes = plt.subplots(1, 3, figsize=(16.0, 5.2), constrained_layout=True)
    for ax, dataset in zip(axes[:2], ("dynamic", "hppc")):
        frame = ranking[ranking["Dataset"] == dataset].set_index("Parameter").loc[list(PARAMETERS)]
        ax.bar(PARAMETERS, frame["RMS_mV_per_ln_parameter"], color="#0072B2")
        ax.set(
            title="{} sensitivity".format(dataset.upper()),
            ylabel="RMS sensitivity (mV / ln parameter)",
        )
        ax.grid(axis="y", alpha=0.25)
    cosine = pd.read_csv(OUT / "combined_equal_dataset_cosine.csv", index_col=0).to_numpy(float)
    image = axes[2].imshow(cosine, vmin=-1.0, vmax=1.0, cmap="coolwarm")
    axes[2].set_xticks(range(len(PARAMETERS)), PARAMETERS)
    axes[2].set_yticks(range(len(PARAMETERS)), PARAMETERS)
    axes[2].set_title("Combined normalized-sensitivity cosine")
    for i in range(len(PARAMETERS)):
        for j in range(len(PARAMETERS)):
            axes[2].text(j, i, "{:.2f}".format(cosine[i, j]), ha="center", va="center")
    fig.colorbar(image, ax=axes[2], shrink=0.82)
    fig.savefig(OUT / "sensitivity_identifiability.png", dpi=240)
    plt.close(fig)

    manifest = {
        "electrolyte_set": "ORegan2022",
        "brugg_n": BRUGG_N,
        "parameters": PARAMETERS,
        "initial": INITIAL,
        "sensitivity_perturbation": "multiplicative +/-5%",
        "dynamic_data": "cell 6-6, 0.5C/1C/2C charge and discharge; 80 capacity-grid points per branch",
        "hppc_source": str(HPPC_FILE),
        "hppc_data": "cycles 1-10, five key points per cycle",
        "combined_weighting": "dynamic and HPPC blocks divided by sqrt(number of points), giving equal dataset-level least-squares weight",
        "baseline_metrics": baseline_metrics,
    }
    (OUT / "analysis_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("\nBASELINE\n", pd.DataFrame([baseline_metrics]).to_string(index=False))
    print("\nSENSITIVITY\n", ranking.to_string(index=False))
    print("\nIDENTIFIABILITY\n", reports.to_string(index=False))


if __name__ == "__main__":
    main()
