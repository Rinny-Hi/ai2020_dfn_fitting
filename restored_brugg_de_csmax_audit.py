"""Restore Ai2020 Bruggeman, rescreen radii, and audit De/c_s,max corrections."""

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
import ocp_grounded_candidate_search as common
import old_notebook_method_comparison as omc
import ocp_rn_sweep_fixed_transport as rs


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "results" / "260920_mass_independent_nominal_anode"
RESULTS = ROOT / "results" / "260921_restored_brugg_de_csmax_audit"
RESULTS.mkdir(parents=True, exist_ok=True)

NAME = "Ai2020 nominal equilibrium + harvested hysteresis"
RN_VALUES_UM = (2.0, 3.0, 4.0, 5.0)
RP_VALUES_UM = (2.0, 3.0, 4.0, 5.0)
DE_FACTORS = (1e-4, 1e-3, 1e-2, 1e-1, 1.0)
CSN_CORRECTED = 29700.0
CSN_AI2020 = 28700.0


def build_base(bundle, de_factor, csn_max):
    original_factory = ga.make_source_corrected_ai2020

    def factory():
        raw = pybamm.ParameterValues("Ai2020")
        values = raw.copy()
        values.update(
            {
                "Electrolyte diffusivity [m2.s-1]": ga.scale_parameter_function(
                    raw["Electrolyte diffusivity [m2.s-1]"], de_factor
                ),
                "Maximum concentration in negative electrode [mol.m-3]": csn_max,
            }
        )
        return values

    ga.make_source_corrected_ai2020 = factory
    try:
        return ga.build_legacy_model_inputs(bundle)
    finally:
        ga.make_source_corrected_ai2020 = original_factory


def set_dynamic_parameters(model, rn_um, rp_um, de_factor):
    changed = dict(model)
    changed["params"] = model["params"].copy()
    raw = pybamm.ParameterValues("Ai2020")
    changed["params"].update(
        {
            "Negative particle diffusivity [m2.s-1]": rs.DSN,
            "Positive particle diffusivity [m2.s-1]": rs.DSP,
            "Negative particle radius [m]": rn_um * 1e-6,
            "Positive particle radius [m]": rp_um * 1e-6,
            "Electrolyte diffusivity [m2.s-1]": ga.scale_parameter_function(
                raw["Electrolyte diffusivity [m2.s-1]"], de_factor
            ),
            # Explicit restoration requested by the user.
            "Negative electrode Bruggeman coefficient (electrolyte)": 2.914,
            "Positive electrode Bruggeman coefficient (electrolyte)": 1.83,
            "Separator Bruggeman coefficient (electrolyte)": 1.5,
        },
        check_already_exists=False,
    )
    return changed


def run_case(case, model, stage, anode, cathode, experiment, simulations, metadata):
    rows = []
    cathode_dynamic = common.scaled_detail(cathode, 0.50)
    anode_dynamic = common.scaled_detail(anode, 0.50)
    for rate in ga.RATES:
        for charge in (True, False):
            direction = "Charge" if charge else "Discharge"
            print(f"{case} | {rate:g}C | {direction}", flush=True)
            sim = ehq.run_dfn(
                rate,
                charge,
                stage,
                model,
                anode_dynamic,
                cathode_dynamic,
                True,
                True,
            )
            simulations[(case, rate, charge)] = sim
            result = omc.old_time_metrics(experiment[(rate, charge)], sim)
            capacity_error = result["end_time_error_min"] * rate * 2.28 / 60.0
            rows.append(
                {
                    "Case": case,
                    **metadata,
                    "C_rate": rate,
                    "Direction": direction,
                    **result,
                    "capacity_error_Ah": capacity_error,
                }
            )
    return rows


def summarize(frame, group_columns):
    return frame.groupby(group_columns, as_index=False).agg(
        Dynamic_MAE_SOC10_70_mV=("MAE_SOC10_70_mV", "mean"),
        Charge_MAE_SOC10_70_mV=("MAE_SOC10_70_mV", lambda x: float(x[frame.loc[x.index, "Direction"] == "Charge"].mean())),
        Discharge_MAE_SOC10_70_mV=("MAE_SOC10_70_mV", lambda x: float(x[frame.loc[x.index, "Direction"] == "Discharge"].mean())),
        Mean_abs_capacity_error_mAh=("capacity_error_Ah", lambda x: float(np.mean(np.abs(x))) * 1000.0),
        Max_abs_capacity_error_mAh=("capacity_error_Ah", lambda x: float(np.max(np.abs(x))) * 1000.0),
    )


def main():
    bundle = ga.load_legacy_bundle(ga.get_paths(ROOT))
    base = build_base(bundle, 1e-4, CSN_CORRECTED)
    measured, nominal, hybrid, cathode = rs.candidate_details(bundle, base)
    experiment = omc.load_old_dynamic_data(base["dynamic"], base["q_meas"])
    metrics = pd.read_csv(SOURCE / "mass_independent_ocp_candidate_metrics.csv")
    metric_row = metrics[metrics.Candidate == NAME].iloc[0]
    stage = {key: float(metric_row[key]) for key in ("x0", "x100", "y100", "y0")}
    effective = mit.model_for_effective_capacities(base, metric_row)
    simulations = {}

    # 1) Radius screen with restored Ai2020 electrolyte Bruggeman values.
    radius_rows = []
    for rn_um in RN_VALUES_UM:
        for rp_um in RP_VALUES_UM:
            case = f"Rn={rn_um:g},Rp={rp_um:g}"
            model = set_dynamic_parameters(effective, rn_um, rp_um, 1e-4)
            radius_rows.extend(
                run_case(
                    case,
                    model,
                    stage,
                    hybrid,
                    cathode,
                    experiment,
                    simulations,
                    {"Rn_um": rn_um, "Rp_um": rp_um},
                )
            )
    radius_validation = pd.DataFrame(radius_rows)
    radius_summary = summarize(radius_validation, ["Case", "Rn_um", "Rp_um"])
    best = radius_summary.sort_values(
        ["Dynamic_MAE_SOC10_70_mV", "Mean_abs_capacity_error_mAh"]
    ).iloc[0]
    best_rn = float(best.Rn_um)
    best_rp = float(best.Rp_um)

    # 2) Audit the legacy De multiplier at the selected geometry.
    de_rows = []
    for factor in DE_FACTORS:
        case = f"De_factor={factor:g}"
        model = set_dynamic_parameters(effective, best_rn, best_rp, factor)
        de_rows.extend(
            run_case(
                case,
                model,
                stage,
                hybrid,
                cathode,
                experiment,
                simulations,
                {"De_factor_vs_Ai2020": factor},
            )
        )
    de_validation = pd.DataFrame(de_rows)
    de_summary = summarize(de_validation, ["Case", "De_factor_vs_Ai2020"])

    # 3) Audit the Stage-0 c_s,n,max correction through the complete current
    # mass-independent pipeline.  Qn is fixed by the fitted window, so the
    # effective c_s,max should cancel algebraically; only the harvested
    # hysteresis-to-stoichiometry mapping can retain a small dependency.
    cs_rows = []
    cs_meta = []
    hybrid_sources = {}
    for source_csn in (CSN_AI2020, CSN_CORRECTED):
        alt_base = build_base(bundle, 1e-4, source_csn)
        _, _, alt_hybrid, alt_cathode = rs.candidate_details(bundle, alt_base)
        alt_effective = mit.model_for_effective_capacities(alt_base, metric_row)
        alt_model = set_dynamic_parameters(alt_effective, best_rn, best_rp, 1e-4)
        case = f"Stage0_csn={source_csn:.0f}"
        hybrid_sources[source_csn] = alt_hybrid
        cs_rows.extend(
            run_case(
                case,
                alt_model,
                stage,
                alt_hybrid,
                alt_cathode,
                experiment,
                simulations,
                {"Stage0_csn_max_mol_m3": source_csn},
            )
        )
        cs_meta.append(
            {
                "Stage0_csn_max_mol_m3": source_csn,
                "Stage0_Qn_Ah": alt_base["qn_cell"],
                "Effective_csn_after_Qn_fit_mol_m3": alt_effective["csn_max"],
                "Fitted_Qn_Ah": float(metric_row.Q_n_Ah),
            }
        )
    cs_validation = pd.DataFrame(cs_rows)
    cs_summary = summarize(cs_validation, ["Case", "Stage0_csn_max_mol_m3"])
    cs_meta = pd.DataFrame(cs_meta)
    cs_summary = cs_summary.merge(cs_meta, on="Stage0_csn_max_mol_m3", how="left")

    grid = np.linspace(0.02, 0.80, 2000)
    h_raw = hybrid_sources[CSN_AI2020]
    h_corrected = hybrid_sources[CSN_CORRECTED]
    branch_differences = {}
    for branch in ("lithiation", "delithiation"):
        raw_curve = np.interp(grid, h_raw["grid"], h_raw[branch])
        corrected_curve = np.interp(grid, h_corrected["grid"], h_corrected[branch])
        diff = (corrected_curve - raw_curve) * 1000.0
        branch_differences[branch] = {
            "MAE_mV": float(np.mean(np.abs(diff))),
            "RMSE_mV": float(np.sqrt(np.mean(diff**2))),
            "max_abs_mV": float(np.max(np.abs(diff))),
        }

    radius_validation.to_csv(RESULTS / "restored_brugg_radius_validation.csv", index=False, encoding="utf-8-sig")
    radius_summary.to_csv(RESULTS / "restored_brugg_radius_summary.csv", index=False, encoding="utf-8-sig")
    de_validation.to_csv(RESULTS / "de_multiplier_validation.csv", index=False, encoding="utf-8-sig")
    de_summary.to_csv(RESULTS / "de_multiplier_summary.csv", index=False, encoding="utf-8-sig")
    cs_validation.to_csv(RESULTS / "csmax_validation.csv", index=False, encoding="utf-8-sig")
    cs_summary.to_csv(RESULTS / "csmax_summary.csv", index=False, encoding="utf-8-sig")

    pivot = radius_summary.pivot(index="Rn_um", columns="Rp_um", values="Dynamic_MAE_SOC10_70_mV")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), constrained_layout=True)
    im = axes[0].imshow(pivot.to_numpy(), origin="lower", cmap="viridis")
    axes[0].set_xticks(range(len(pivot.columns)), [f"{v:g}" for v in pivot.columns])
    axes[0].set_yticks(range(len(pivot.index)), [f"{v:g}" for v in pivot.index])
    axes[0].set_xlabel("Rp (um)")
    axes[0].set_ylabel("Rn (um)")
    axes[0].set_title("Dynamic MAE with restored Bruggeman")
    for i in range(len(pivot.index)):
        for j in range(len(pivot.columns)):
            axes[0].text(j, i, f"{pivot.iloc[i, j]:.1f}", ha="center", va="center", color="white" if pivot.iloc[i, j] > pivot.to_numpy().mean() else "black")
    fig.colorbar(im, ax=axes[0], label="mV")
    axes[1].semilogx(de_summary.De_factor_vs_Ai2020, de_summary.Dynamic_MAE_SOC10_70_mV, "o-", lw=2, label="Dynamic MAE")
    axes[1].set_xlabel("Electrolyte diffusivity multiplier vs Ai2020")
    axes[1].set_ylabel("Dynamic MAE (mV)")
    axes[1].set_title("Legacy De multiplier audit")
    axes[1].grid(alpha=0.25, which="both")
    fig.suptitle("Restored b_n=2.914, b_p=1.83, b_s=1.5 | Dsn/Dsp user values")
    fig.savefig(RESULTS / "restored_brugg_radius_and_de_audit.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2), constrained_layout=True)
    axes[0].bar(cs_summary.Stage0_csn_max_mol_m3.astype(str), cs_summary.Dynamic_MAE_SOC10_70_mV, color="#4c78a8")
    axes[0].set_xlabel("Stage-0 negative c_s,max (mol/m3)")
    axes[0].set_ylabel("Dynamic MAE (mV)")
    axes[0].set_title("Full current-pipeline sensitivity")
    axes[1].bar(cs_meta.Stage0_csn_max_mol_m3.astype(str), cs_meta.Effective_csn_after_Qn_fit_mol_m3, color="#f58518")
    axes[1].set_xlabel("Stage-0 negative c_s,max (mol/m3)")
    axes[1].set_ylabel("Effective c_s,max after Qn fit (mol/m3)")
    axes[1].set_title("Capacity-consistent cancellation")
    for ax in axes:
        ax.grid(axis="y", alpha=0.2)
    fig.savefig(RESULTS / "csmax_pipeline_audit.png", dpi=220)
    plt.close(fig)

    report = {
        "restored_bruggeman": {"negative": 2.914, "positive": 1.83, "separator": 1.5},
        "fixed_solid_diffusivities": {"Dsn_m2_s": rs.DSN, "Dsp_m2_s": rs.DSP},
        "best_radius": best.to_dict(),
        "de_note": (
            "PyBaMM 26.8 Ai2020 returns about 3.22e-6 at 1000 mol/m3 and "
            "298.15 K despite the parameter label m2/s. The 1e-4 multiplier "
            "converts the source correlation from cm2/s to about 3.22e-10 m2/s; "
            "using a factor of 1 would be physically inconsistent."
        ),
        "de_summary": de_summary.to_dict(orient="records"),
        "csmax_summary": cs_summary.to_dict(orient="records"),
        "hybrid_branch_difference_due_to_stage0_csn": branch_differences,
    }
    (RESULTS / "restored_brugg_de_csmax_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print("\nBEST RADIUS")
    print(best.to_string())
    print("\nDE SUMMARY")
    print(de_summary.to_string(index=False))
    print("\nCSMAX SUMMARY")
    print(cs_summary.to_string(index=False))
    print("\nBRANCH DIFFERENCES", branch_differences)
    print("\nSaved to", RESULTS)


if __name__ == "__main__":
    main()
