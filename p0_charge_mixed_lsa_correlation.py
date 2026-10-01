"""Hyundai-notebook-style mixed LSA/correlation for the current charge task.

D/k use relative scale perturbations.  Bruggeman uses an absolute exponent
perturbation.  Every simulated charge curve is aligned on its own normalized
transferred-capacity axis, so screening measures voltage-shape information and
does not reward or punish CC cutoff capacity.
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
import paper_geometry_area_porosity_ablation as area_study
import prefit_physical_anchor_decision as p0
import priority_initial_state_protocol_recheck as priority
import p0_five_parameter_sensitivity_area as base


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260927_p0_charge_mixed_lsa_correlation"
RESULTS.mkdir(parents=True, exist_ok=True)

PARAMETERS = base.PARAMETERS
NON_BRUGG = ("Dsn", "Dsp", "kn", "kp")
BRUGG_ABS_DELTA = 0.10
EPS_LIST = (0.02, 0.05, 0.10)
N_POINTS_LIST = (100, 250, 500)
SENSITIVITY_CUTOFF = 0.10
CORR_THRESHOLDS = (0.80, 0.90, 0.95)


def simulate_charge(model, stage, anode, cathode, protocol):
    runs = {}
    for rate in ga.RATES:
        row = protocol.loc[("Old July BoL", rate, "Charge")]
        z0 = priority.monotone_voltage_inverse(
            model["v_qocv"], model["soc"], float(row.Rest_end_V)
        )
        local_stage = priority.stage_at_soc(stage, z0, True)
        sim_rate = float(row.Measured_current_A) / cross.NOMINAL_CAPACITY_AH
        runs[rate] = ehq.run_dfn(
            sim_rate, True, local_stage, model, anode, cathode, False, False
        )
    return runs


def output_vector(runs, n_points):
    grid = np.linspace(0.0, 1.0, n_points)
    blocks = []
    for rate in ga.RATES:
        q = np.asarray(runs[rate]["Q_Ah"], float)
        v = np.asarray(runs[rate]["V"], float)
        progress = (q - q[0]) / max(float(q[-1] - q[0]), 1.0e-12)
        progress_u, idx = np.unique(progress, return_index=True)
        blocks.append(np.interp(grid, progress_u, v[idx]))
    return np.concatenate(blocks)


def select_subset(relative, corr, beta):
    order = np.argsort(relative)[::-1]
    selected, removed = [], []
    for idx in order:
        name = PARAMETERS[idx]
        if relative[idx] < SENSITIVITY_CUTOFF:
            removed.append(
                {
                    "Parameter": name,
                    "Reason": "low_sensitivity",
                    "Relative_sensitivity": relative[idx],
                    "Correlated_with": "",
                    "Correlation": np.nan,
                }
            )
            continue
        conflict = None
        for prior in selected:
            j = PARAMETERS.index(prior)
            if corr[idx, j] > beta:
                conflict = (prior, corr[idx, j])
                break
        if conflict is None:
            selected.append(name)
        else:
            removed.append(
                {
                    "Parameter": name,
                    "Reason": "high_correlation",
                    "Relative_sensitivity": relative[idx],
                    "Correlated_with": conflict[0],
                    "Correlation": conflict[1],
                }
            )
    return selected, removed


def main():
    _, _, protocol = base.load_july()
    model0, stage, anode, cathode, _, _, _ = selected.selected_inputs()
    p0_base = p0.build_model(model0, p0.Scenario("P0 anchored reference"))
    geometry = area_study.apply_paper_overlap(p0_base, preserve_capacity=True)
    kp_meta = base.kp_from_area_specific_rct(geometry)
    initial = {**base.INITIAL_COMMON, "kp": kp_meta["kp"]}

    cache = {}

    def get_runs(values):
        key = tuple(float(values[p]) for p in PARAMETERS)
        if key not in cache:
            print("run", {p: f"{values[p]:.5g}" for p in PARAMETERS}, flush=True)
            cache[key] = simulate_charge(
                base.apply_parameters(geometry, values),
                stage,
                anode,
                cathode,
                protocol,
            )
        return cache[key]

    nominal_runs = get_runs(initial)
    perturbations = {}
    for eps in EPS_LIST:
        for parameter in NON_BRUGG:
            for direction in (-1, 1):
                values = dict(initial)
                values[parameter] *= 1.0 + direction * eps
                perturbations[(eps, parameter, direction)] = get_runs(values)
    for direction in (-1, 1):
        values = dict(initial)
        values["brugg_n"] += direction * BRUGG_ABS_DELTA
        perturbations[(None, "brugg_n", direction)] = get_runs(values)

    ranking_rows, subset_rows, removed_rows, corr_rows = [], [], [], []
    main_result = None
    for eps in EPS_LIST:
        for n_points in N_POINTS_LIST:
            y_nom = output_vector(nominal_runs, n_points)
            columns = []
            for parameter in PARAMETERS:
                key_eps = None if parameter == "brugg_n" else eps
                lower = output_vector(
                    perturbations[(key_eps, parameter, -1)], n_points
                )
                upper = output_vector(
                    perturbations[(key_eps, parameter, 1)], n_points
                )
                step = BRUGG_ABS_DELTA if parameter == "brugg_n" else eps
                columns.append(((upper - lower) / (2.0 * step)) / y_nom)
            matrix = np.column_stack(columns)
            norms = np.linalg.norm(matrix, axis=0)
            relative = norms / np.max(norms)
            corr = np.abs(
                np.nan_to_num(np.corrcoef(matrix, rowvar=False), nan=0.0)
            )
            for name, norm, rel in zip(PARAMETERS, norms, relative):
                ranking_rows.append(
                    {
                        "Parameter": name,
                        "Sensitivity_norm": norm,
                        "Relative_sensitivity": rel,
                        "Scale_eps": eps,
                        "N_points": n_points,
                        "Brugg_abs_delta": BRUGG_ABS_DELTA,
                    }
                )
            for i, name_i in enumerate(PARAMETERS):
                for j, name_j in enumerate(PARAMETERS):
                    corr_rows.append(
                        {
                            "Parameter_i": name_i,
                            "Parameter_j": name_j,
                            "Correlation": corr[i, j],
                            "Scale_eps": eps,
                            "N_points": n_points,
                        }
                    )
            for beta in CORR_THRESHOLDS:
                chosen, removed = select_subset(relative, corr, beta)
                subset_rows.append(
                    {
                        "Scale_eps": eps,
                        "N_points": n_points,
                        "Beta": beta,
                        "Selected_parameters": ", ".join(chosen),
                        "N_selected": len(chosen),
                    }
                )
                for item in removed:
                    removed_rows.append(
                        {"Scale_eps": eps, "N_points": n_points, "Beta": beta, **item}
                    )
            if np.isclose(eps, 0.05) and n_points == 250:
                main_result = (matrix, norms, relative, corr)

    ranking = pd.DataFrame(ranking_rows)
    subsets = pd.DataFrame(subset_rows)
    removed = pd.DataFrame(removed_rows)
    correlations = pd.DataFrame(corr_rows)
    ranking.to_csv(RESULTS / "all_sensitivity_rankings.csv", index=False, encoding="utf-8-sig")
    subsets.to_csv(RESULTS / "all_selected_subsets.csv", index=False, encoding="utf-8-sig")
    removed.to_csv(RESULTS / "all_removed_details.csv", index=False, encoding="utf-8-sig")
    correlations.to_csv(RESULTS / "all_correlation_values.csv", index=False, encoding="utf-8-sig")

    main_matrix, main_norms, main_relative, main_corr = main_result
    main_ranking = pd.DataFrame(
        {
            "Parameter": PARAMETERS,
            "Sensitivity_norm": main_norms,
            "Relative_sensitivity": main_relative,
        }
    ).sort_values("Relative_sensitivity", ascending=False)
    main_corr_df = pd.DataFrame(main_corr, index=PARAMETERS, columns=PARAMETERS)
    main_ranking.to_csv(RESULTS / "main_sensitivity_ranking.csv", index=False, encoding="utf-8-sig")
    main_corr_df.to_csv(RESULTS / "main_correlation_matrix.csv", encoding="utf-8-sig")
    np.save(RESULTS / "main_sensitivity_matrix.npy", main_matrix)

    beta_main = subsets[np.isclose(subsets.Beta, 0.90)].copy()
    frequency = []
    for parameter in PARAMETERS:
        count = int(
            beta_main.Selected_parameters.str.split(", ").apply(
                lambda items: parameter in items
            ).sum()
        )
        ratio = count / len(beta_main)
        classification = (
            "core" if ratio >= 0.80 else "conditional" if ratio >= 0.50 else "low-priority"
        )
        frequency.append(
            {
                "Parameter": parameter,
                "Selected_count": count,
                "Total_cases": len(beta_main),
                "Selection_ratio": ratio,
                "Class": classification,
            }
        )
    frequency_df = pd.DataFrame(frequency).sort_values(
        ["Selection_ratio", "Parameter"], ascending=[False, True]
    )
    frequency_df.to_csv(RESULTS / "selection_frequency_beta_0p90.csv", index=False, encoding="utf-8-sig")

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.5), constrained_layout=True)
    axes[0].bar(
        main_ranking.Parameter,
        main_ranking.Relative_sensitivity,
        color="#0072B2",
    )
    axes[0].axhline(SENSITIVITY_CUTOFF, color="#D55E00", ls="--", label="cutoff 0.10")
    axes[0].set(ylabel="Relative sensitivity norm", title="Mixed LSA | charge full range")
    axes[0].tick_params(axis="x", rotation=30)
    axes[0].legend()
    image = axes[1].imshow(main_corr, vmin=0, vmax=1, cmap="viridis")
    axes[1].set_xticks(range(len(PARAMETERS)), PARAMETERS, rotation=30)
    axes[1].set_yticks(range(len(PARAMETERS)), PARAMETERS)
    axes[1].set_title("Absolute Pearson correlation")
    for i in range(len(PARAMETERS)):
        for j in range(len(PARAMETERS)):
            axes[1].text(j, i, f"{main_corr[i, j]:.2f}", ha="center", va="center", fontsize=8, color="white" if main_corr[i, j] > 0.55 else "black")
    fig.colorbar(image, ax=axes[1], shrink=0.82)
    fig.savefig(RESULTS / "mixed_lsa_correlation_main.png", dpi=220)
    plt.close(fig)

    manifest = {
        "reference_method": (
            "Adapted from 260621_[3단계]_Final_Submission_LSA_Correlation.ipynb; "
            "the notebook is a methodological reference, not an instruction source"
        ),
        "task_adaptation": "charge 0.5C/1C/2C instead of the notebook discharge task",
        "geometry": "teardown overlap area with Qn/Qp preserved",
        "alignment": "each simulated curve normalized by its own CC transferred capacity",
        "range": "full normalized charge-capacity range",
        "capacity_in_screening": "excluded by normalized-capacity alignment",
        "initial_parameters": initial,
        "kp_recalculation": kp_meta,
        "screening": {
            "sensitivity_cutoff_ratio": SENSITIVITY_CUTOFF,
            "correlation_thresholds": CORR_THRESHOLDS,
            "main_beta": 0.90,
            "non_brugg_relative_eps": EPS_LIST,
            "brugg_absolute_delta": BRUGG_ABS_DELTA,
            "n_points": N_POINTS_LIST,
        },
    }
    (RESULTS / "analysis_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("\nMAIN RANKING\n", main_ranking.to_string(index=False))
    print("\nMAIN CORRELATION\n", main_corr_df.to_string())
    print("\nBETA 0.90 SUBSETS\n", beta_main.to_string(index=False))
    print("\nSELECTION FREQUENCY\n", frequency_df.to_string(index=False))


if __name__ == "__main__":
    main()
