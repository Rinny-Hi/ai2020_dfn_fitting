"""Area ablation and local identifiability for the five requested DFN parameters.

The July BoL three-cell cohort is kept as the modelling cohort.  Capacity is
never used in the local sensitivity residual.  The residual is terminal
voltage at the experimental time points within 10--70% of the measured C/50
reference capacity, with equal numbers of points for every cell and C-rate.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import c50_selected_ocp_dynamic_validation as selected
import comprehensive_prefit_cross_cohort as cross
import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import kinetics_original_reversion_ablation as kinetics
import paper_geometry_area_porosity_ablation as area_study
import prefit_physical_anchor_decision as p0
import priority_initial_state_protocol_recheck as priority


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260927_p0_five_parameter_sensitivity_area"
RESULTS.mkdir(parents=True, exist_ok=True)

R_GAS = 8.314462618
FARADAY = 96485.33212
T_REF_K = 298.15
CE_REF = 1000.0
RCT_P_ASR_OHM_CM2 = 18.78433
AI_K = 1.0e-11 * FARADAY

PARAMETERS = ("Dsn", "Dsp", "kn", "kp", "brugg_n")
INITIAL_COMMON = {
    "Dsn": 2.30e-14,
    "Dsp": 4.36e-14,
    "kn": AI_K,
    "brugg_n": 2.914,
}
BOUNDS = {
    "Dsn": (1.0e-14, 5.0e-14),
    "Dsp": (2.0e-14, 9.0e-14),
    "kn": (3.0e-7, 3.0e-6),
    "kp": (3.0e-7, 1.2e-6),
    "brugg_n": (1.5, 4.0),
}


def load_july():
    curves: dict[str, dict] = {}
    rows: list[dict] = []
    for cell, path in cross.OLD_FILES.items():
        cell_curves, cell_rows = cross.parse_cell(
            "Old July BoL", cell, path, cross.OLD_BRANCH_STEPS
        )
        curves[cell] = cell_curves
        rows.extend(cell_rows)
    qc = pd.DataFrame(rows)
    return curves, qc, cross.cohort_protocol_means(qc)


def geometric_area_m2(model: dict) -> float:
    params = model["params"]
    return float(params["Electrode height [m]"]) * float(
        params["Electrode width [m]"]
    ) * float(params["Number of electrodes connected in parallel to make a cell"])


def kp_from_area_specific_rct(model: dict) -> dict[str, float]:
    """Convert cathode area-specific Rct to the PyBaMM kinetic prefactor.

    Rct is labelled ohm cm2 in the EIS table.  Therefore the macroscopic
    electrode area cancels, while Rp, Lp, eps_s and c_s,max remain in the
    conversion through b=3 eps_s/Rp and S/A=bL.
    """

    params = model["params"]
    eps = float(params["Positive electrode active material volume fraction"])
    radius = float(params["Positive particle radius [m]"])
    thickness = float(params["Positive electrode thickness [m]"])
    cmax = float(params["Maximum concentration in positive electrode [mol.m-3]"])
    b = 3.0 * eps / radius
    rct_asr = RCT_P_ASR_OHM_CM2 * 1.0e-4
    j0 = R_GAS * T_REF_K / (FARADAY * rct_asr * b * thickness)
    kp = 2.0 * j0 / (cmax * np.sqrt(CE_REF))
    return {
        "eps_s_p": eps,
        "Rp_um": radius * 1.0e6,
        "Lp_um": thickness * 1.0e6,
        "csp_max_mol_m3": cmax,
        "specific_surface_area_m2_m3": b,
        "j0_A_m2": j0,
        "kp": kp,
    }


def apply_parameters(base: dict, values: dict[str, float]) -> dict:
    changed = kinetics.with_prefactors(base, values["kn"], values["kp"])
    changed = dict(changed)
    changed["params"] = changed["params"].copy()
    changed["params"].update(
        {
            "Negative particle diffusivity [m2.s-1]": values["Dsn"],
            "Positive particle diffusivity [m2.s-1]": values["Dsp"],
            "Negative electrode Bruggeman coefficient (electrolyte)": values[
                "brugg_n"
            ],
        },
        check_already_exists=False,
    )
    return changed


def simulate(model, stage, anode, cathode, protocol, charge: bool):
    runs = {}
    direction = "Charge" if charge else "Discharge"
    for rate in ga.RATES:
        row = protocol.loc[("Old July BoL", rate, direction)]
        z0 = priority.monotone_voltage_inverse(
            model["v_qocv"], model["soc"], float(row.Rest_end_V)
        )
        local_stage = priority.stage_at_soc(stage, z0, charge)
        sim_rate = float(row.Measured_current_A) / cross.NOMINAL_CAPACITY_AH
        print(f"{rate:g}C {direction}", flush=True)
        runs[rate] = ehq.run_dfn(
            sim_rate, charge, local_stage, model, anode, cathode, False, False
        )
    return runs


def summarize_runs(name, runs_by_direction, curves, qcell):
    rows = []
    for charge, runs in runs_by_direction.items():
        direction = "Charge" if charge else "Discharge"
        for rate, sim in runs.items():
            for cell, branches in curves.items():
                metric = priority.curve_metrics(
                    branches[(rate, charge)], sim, qcell, charge
                )
                rows.append(
                    {
                        "Area_policy": name,
                        "Cell": cell,
                        "C_rate": rate,
                        "Direction": direction,
                        **metric,
                        "Capacity_error_mAh": 1000.0
                        * (metric["Q_end_model_Ah"] - metric["Q_end_exp_Ah"]),
                    }
                )
    return pd.DataFrame(rows)


def objective_grids(curves, qcell, baseline_runs, n_points=80):
    grids = {}
    # The most resistive +5% candidate reaches the 2C cutoff before 70%.
    # Restrict the identification window to the common mid-range instead of
    # padding/extrapolating a voltage after cutoff.
    lo_q, hi_q = 0.10 * qcell, 0.55 * qcell
    for rate in ga.RATES:
        for cell, branches in curves.items():
            obs = branches[(rate, True)]
            q = np.asarray(obs["Q_Ah"], float)
            t = np.asarray(obs["t_min"], float)
            mask = (q >= lo_q) & (q <= min(hi_q, 0.995 * q[-1]))
            if np.count_nonzero(mask) < 3:
                raise RuntimeError(f"Insufficient objective range: {cell} {rate:g}C")
            # Keep one immutable time domain for every later candidate.  The
            # baseline model may reach cutoff before the experimental 70% time,
            # so use the shared interval rather than silently extrapolating.
            upper_t = min(t[mask][-1], 0.95 * float(baseline_runs[rate]["t_min"][-1]))
            grids[(rate, cell)] = np.linspace(t[mask][0], upper_t, n_points)
    return grids


def voltage_vector(runs, curves, grids, experiment=False):
    blocks = []
    for rate in ga.RATES:
        sim = None if experiment else runs[rate]
        for cell, branches in curves.items():
            grid = grids[(rate, cell)]
            source = branches[(rate, True)] if experiment else sim
            if grid[-1] > float(np.max(source["t_min"])):
                raise RuntimeError(
                    f"Simulation ended before objective grid: {rate:g}C {cell}"
                )
            blocks.append(np.interp(grid, source["t_min"], source["V"]))
    return np.concatenate(blocks)


def plot_area_comparison(area_runs, curves, area_detail):
    colors = {"DFN area": "#0072B2", "Teardown overlap area": "#D55E00"}
    fig, axes = plt.subplots(2, 3, figsize=(16.2, 8.8), constrained_layout=True)
    for row_i, charge in enumerate((True, False)):
        direction = "Charge" if charge else "Discharge"
        for col_i, rate in enumerate(ga.RATES):
            ax = axes[row_i, col_i]
            for idx, (cell, branches) in enumerate(curves.items()):
                obs = branches[(rate, charge)]
                ax.plot(
                    obs["t_min"],
                    obs["V"],
                    color="black",
                    alpha=0.30,
                    lw=1.0,
                    label="Experiment (July BoL, n=3)" if idx == 0 else None,
                )
            for area_name, by_direction in area_runs.items():
                sim = by_direction[charge][rate]
                sel = area_detail[
                    (area_detail.Area_policy == area_name)
                    & np.isclose(area_detail.C_rate, rate)
                    & (area_detail.Direction == direction)
                ]
                rmse = float(sel.Center10_70_RMSE_mV.mean())
                cap = float(sel.Capacity_error_mAh.mean())
                ax.plot(
                    sim["t_min"],
                    sim["V"],
                    color=colors[area_name],
                    lw=2.0,
                    label=f"{area_name}: {rmse:.1f} mV, {cap:+.1f} mAh",
                )
            ax.set(
                title=f"{rate:g}C {direction}",
                xlabel="Time [min]",
                ylabel="Voltage [V]",
            )
            ax.grid(alpha=0.2)
            ax.legend(fontsize=7.5)
    fig.suptitle("DFN area versus teardown overlap area | no dynamic fitting")
    fig.savefig(RESULTS / "area_time_voltage_comparison.png", dpi=220)
    plt.close(fig)


def main():
    curves, qc, protocol = load_july()
    model0, stage, anode, cathode, _, _, qcell = selected.selected_inputs()
    p0_base = p0.build_model(model0, p0.Scenario("P0 anchored reference"))

    # Compare the current DFN area against the smaller teardown overlap.  The
    # latter rescales c_s,max so Qn/Qp and the OCP window are unchanged.
    geometries = {
        "DFN area": p0_base,
        "Teardown overlap area": area_study.apply_paper_overlap(
            p0_base, preserve_capacity=True
        ),
    }
    area_rows = []
    area_runs = {}
    geometry_meta = {}
    for area_name, geometry in geometries.items():
        kp_meta = kp_from_area_specific_rct(geometry)
        initial = {**INITIAL_COMMON, "kp": kp_meta["kp"]}
        model = apply_parameters(geometry, initial)
        runs = {
            True: simulate(model, stage, anode, cathode, protocol, True),
            False: simulate(model, stage, anode, cathode, protocol, False),
        }
        area_runs[area_name] = runs
        area_rows.append(summarize_runs(area_name, runs, curves, qcell))
        geometry_meta[area_name] = {
            "area_cm2": geometric_area_m2(geometry) * 1.0e4,
            "initial_parameters": initial,
            "kp_recalculation": kp_meta,
        }

    area_detail = pd.concat(area_rows, ignore_index=True)
    area_summary = (
        area_detail.groupby(["Area_policy", "Direction"], as_index=False)
        .agg(
            Full_RMSE_mV=("Full_RMSE_mV", "mean"),
            Center10_70_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
            Capacity_RMSE_pct=(
                "Capacity_error_pct",
                lambda x: float(np.sqrt(np.mean(np.asarray(x, float) ** 2))),
            ),
            Mean_capacity_error_mAh=("Capacity_error_mAh", "mean"),
        )
        .sort_values(["Direction", "Center10_70_RMSE_mV"])
    )
    area_detail.to_csv(RESULTS / "area_detail.csv", index=False, encoding="utf-8-sig")
    area_summary.to_csv(
        RESULTS / "area_summary.csv", index=False, encoding="utf-8-sig"
    )
    plot_area_comparison(area_runs, curves, area_detail)

    # Use the physically measured overlap-area bookkeeping for the sensitivity
    # analysis.  This choice is not made from the voltage error ranking.
    fitting_geometry = geometries["Teardown overlap area"]
    kp_meta = kp_from_area_specific_rct(fitting_geometry)
    initial = {**INITIAL_COMMON, "kp": kp_meta["kp"]}
    base_model = apply_parameters(fitting_geometry, initial)
    base_runs = area_runs["Teardown overlap area"][True]
    grids = objective_grids(curves, qcell, base_runs)
    experiment_vector = voltage_vector(None, curves, grids, experiment=True)
    baseline_vector = voltage_vector(base_runs, curves, grids)
    baseline_rmse = float(
        np.sqrt(np.mean((1000.0 * (baseline_vector - experiment_vector)) ** 2))
    )

    delta = np.log(1.05)
    sensitivity_columns = []
    sensitivity_rows = []
    rate_rows = []
    for parameter in PARAMETERS:
        values_minus = dict(initial)
        values_plus = dict(initial)
        values_minus[parameter] *= np.exp(-delta)
        values_plus[parameter] *= np.exp(delta)
        model_minus = apply_parameters(fitting_geometry, values_minus)
        model_plus = apply_parameters(fitting_geometry, values_plus)
        print(f"Sensitivity {parameter} -", flush=True)
        runs_minus = simulate(model_minus, stage, anode, cathode, protocol, True)
        print(f"Sensitivity {parameter} +", flush=True)
        runs_plus = simulate(model_plus, stage, anode, cathode, protocol, True)
        v_minus = voltage_vector(runs_minus, curves, grids)
        v_plus = voltage_vector(runs_plus, curves, grids)
        derivative = 1000.0 * (v_plus - v_minus) / (2.0 * delta)
        sensitivity_columns.append(derivative)
        sensitivity_rows.append(
            {
                "Parameter": parameter,
                "Initial": initial[parameter],
                "Lower_bound": BOUNDS[parameter][0],
                "Upper_bound": BOUNDS[parameter][1],
                "RMS_mV_per_ln_parameter": float(np.sqrt(np.mean(derivative**2))),
                "Max_abs_mV_per_ln_parameter": float(np.max(np.abs(derivative))),
            }
        )
        offset = 0
        n_block = len(next(iter(grids.values())))
        for rate in ga.RATES:
            n = len(curves) * n_block
            block = derivative[offset : offset + n]
            offset += n
            rate_rows.append(
                {
                    "Parameter": parameter,
                    "C_rate": rate,
                    "RMS_mV_per_ln_parameter": float(np.sqrt(np.mean(block**2))),
                }
            )

    matrix = np.column_stack(sensitivity_columns)
    norms = np.linalg.norm(matrix, axis=0)
    normalized = matrix / np.maximum(norms, np.finfo(float).eps)
    correlation = normalized.T @ normalized
    singular = np.linalg.svd(normalized, compute_uv=False)
    condition = float(singular[0] / singular[-1])

    ranking = pd.DataFrame(sensitivity_rows).sort_values(
        "RMS_mV_per_ln_parameter", ascending=False
    )
    by_rate = pd.DataFrame(rate_rows)
    corr = pd.DataFrame(correlation, index=PARAMETERS, columns=PARAMETERS)
    ranking.to_csv(
        RESULTS / "five_parameter_sensitivity_ranking.csv",
        index=False,
        encoding="utf-8-sig",
    )
    by_rate.to_csv(
        RESULTS / "five_parameter_sensitivity_by_rate.csv",
        index=False,
        encoding="utf-8-sig",
    )
    corr.to_csv(RESULTS / "five_parameter_sensitivity_cosine.csv", encoding="utf-8-sig")
    np.save(RESULTS / "five_parameter_sensitivity_matrix.npy", matrix)

    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.6), constrained_layout=True)
    axes[0].bar(
        ranking.Parameter,
        ranking.RMS_mV_per_ln_parameter,
        color="#0072B2",
    )
    axes[0].set(
        ylabel="RMS sensitivity [mV / ln(parameter)]",
        title="Local voltage sensitivity | charge 0.5C/1C/2C",
    )
    axes[0].tick_params(axis="x", rotation=30)
    image = axes[1].imshow(correlation, vmin=-1, vmax=1, cmap="coolwarm")
    axes[1].set_xticks(range(len(PARAMETERS)), PARAMETERS, rotation=30)
    axes[1].set_yticks(range(len(PARAMETERS)), PARAMETERS)
    axes[1].set_title(f"Normalized sensitivity cosine | cond={condition:.1f}")
    for i in range(len(PARAMETERS)):
        for j in range(len(PARAMETERS)):
            axes[1].text(j, i, f"{correlation[i, j]:.2f}", ha="center", va="center", fontsize=8)
    fig.colorbar(image, ax=axes[1], shrink=0.85)
    fig.savefig(RESULTS / "five_parameter_sensitivity_identifiability.png", dpi=220)
    plt.close(fig)

    manifest = {
        "cohort": "July BoL, cells 6-5/6-6/6-8",
        "objective_for_sensitivity": (
            "voltage residual on experimental time axis, 10-55% of q_meas, "
            "80 points per cell and C-rate; equal weighting; capacity excluded"
        ),
        "q_meas_Ah": qcell,
        "baseline_objective_RMSE_mV": baseline_rmse,
        "geometry": geometry_meta,
        "fitting_geometry": "Teardown overlap area",
        "kp_assumption": (
            "18.78433 ohm cm2 is area-specific Rct; total cell area cancels. "
            "kp changes with Rp, Lp, eps_s,p and c_s,p,max."
        ),
        "normalized_sensitivity_singular_values": singular.tolist(),
        "normalized_condition_number": condition,
        "perturbation": "multiplicative +/-5%",
    }
    (RESULTS / "analysis_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("\nAREA SUMMARY\n", area_summary.to_string(index=False))
    print("\nSENSITIVITY\n", ranking.to_string(index=False))
    print("\nCOSINE\n", corr.to_string())
    print(f"\nNormalized condition number: {condition:.3g}")


if __name__ == "__main__":
    main()
