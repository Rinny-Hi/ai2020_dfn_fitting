"""Apply the small remeasured coating-mass change to the three Lp/Ln cases."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import c50_selected_ocp_dynamic_validation as selected
import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import kp_eis_recalculation_ablation as kpstudy
import priority_initial_state_protocol_recheck as priority
import three_thickness_measurement_comparison as thickness3


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260927_remeasured_electrode_mass_sensitivity"
RESULTS.mkdir(parents=True, exist_ok=True)

OLD_MASS_MG = {"Negative": 10.40, "Positive": 17.25}
NEW_MASS_MG = {"Negative": 10.10, "Positive": 17.35}
MASS_RATIO = {key: NEW_MASS_MG[key] / OLD_MASS_MG[key] for key in OLD_MASS_MG}


def apply_mass_ratio(model0: dict) -> dict:
    """Scale solid inventory by measured mass ratio; close volume with porosity."""
    changed = thickness3.copy_model(model0)
    p = changed["params"]
    for electrode in ("Negative", "Positive"):
        eps_s_key = f"{electrode} electrode active material volume fraction"
        porosity_key = f"{electrode} electrode porosity"
        old_eps_s = float(p[eps_s_key])
        old_porosity = float(p[porosity_key])
        inactive = 1.0 - old_eps_s - old_porosity
        new_eps_s = old_eps_s * MASS_RATIO[electrode]
        new_porosity = 1.0 - inactive - new_eps_s
        if new_porosity <= 0:
            raise ValueError(f"Non-positive {electrode} porosity: {new_porosity}")
        p.update(
            {eps_s_key: new_eps_s, porosity_key: new_porosity},
            check_already_exists=False,
        )
    return changed


def main() -> None:
    model0, stage, anode, cathode, experiment, _, qcell = selected.selected_inputs()
    model0 = kpstudy.with_kp(model0, kpstudy.KP_EIS)
    protocol, _, _ = priority.load_protocol(model0)
    pidx = protocol.set_index(["C_rate", "Direction"])
    models = {}
    for case, values in thickness3.CASES_UM.items():
        same = thickness3.same_loading_geometry(model0, values["Ln"], values["Lp"])
        models[case] = apply_mass_ratio(same)

    rows = []
    simulations = {}
    for case, model in models.items():
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                rest_v = float(pidx.loc[(rate, direction), "Rest_end_V"])
                initial_soc = priority.monotone_voltage_inverse(model["v_qocv"], model["soc"], rest_v)
                local_stage = priority.stage_at_soc(stage, initial_soc, charge)
                measured_current = float(np.nanmedian(np.abs(experiment[(rate, charge)]["I_A"])))
                sim_rate = measured_current / thickness3.NOMINAL_CAPACITY_AH
                print(f"{case} mass-adjusted | {rate:g}C {direction}", flush=True)
                sim = ehq.run_dfn(
                    sim_rate,
                    charge,
                    local_stage,
                    model,
                    anode,
                    cathode,
                    False,
                    True,
                    positive_initial_h=(-1.0 if charge else 1.0),
                )
                simulations[(case, rate, charge)] = sim
                metric = priority.curve_metrics(experiment[(rate, charge)], sim, qcell, charge)
                rows.append(
                    {
                        "Case": case,
                        "C_rate": rate,
                        "Direction": direction,
                        "Full_RMSE_mV": metric["Full_RMSE_mV"],
                        "Center10_70_RMSE_mV": metric["Center10_70_RMSE_mV"],
                        "Capacity_error_pct": metric["Capacity_error_pct"],
                        "Capacity_error_mAh": 1000.0 * (metric["Q_end_model_Ah"] - metric["Q_end_exp_Ah"]),
                        "eps_s_n": float(model["params"]["Negative electrode active material volume fraction"]),
                        "eps_e_n": float(model["params"]["Negative electrode porosity"]),
                        "eps_s_p": float(model["params"]["Positive electrode active material volume fraction"]),
                        "eps_e_p": float(model["params"]["Positive electrode porosity"]),
                    }
                )

    adjusted = pd.DataFrame(rows)
    adjusted.to_csv(RESULTS / "mass_adjusted_detail.csv", index=False, encoding="utf-8-sig")
    base_path = thickness3.RESULTS / "three_thickness_detail.csv"
    base = pd.read_csv(base_path)
    keys = ["Case", "C_rate", "Direction"]
    metrics = ["Full_RMSE_mV", "Center10_70_RMSE_mV", "Capacity_error_pct", "Capacity_error_mAh"]
    comparison = base[keys + metrics].merge(
        adjusted[keys + metrics], on=keys, suffixes=("_same_mass", "_remeasured_mass")
    )
    for metric in metrics:
        comparison[f"Delta_{metric}"] = (
            comparison[f"{metric}_remeasured_mass"] - comparison[f"{metric}_same_mass"]
        )
    comparison.to_csv(RESULTS / "mass_effect_delta.csv", index=False, encoding="utf-8-sig")

    summary_rows = []
    for label, frame in (("Same mass", base), ("Remeasured mass", adjusted)):
        grouped = frame.groupby(["Case", "Direction"], sort=False)
        for (case, direction), group in grouped:
            summary_rows.append(
                {
                    "Mass_case": label,
                    "Case": case,
                    "Direction": direction,
                    "Mean_full_RMSE_mV": float(group.Full_RMSE_mV.mean()),
                    "Mean_center_RMSE_mV": float(group.Center10_70_RMSE_mV.mean()),
                    "Capacity_RMSE_pct": float(np.sqrt(np.mean(group.Capacity_error_pct**2))),
                    "Mean_abs_capacity_error_mAh": float(np.mean(np.abs(group.Capacity_error_mAh))),
                }
            )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(RESULTS / "mass_effect_summary.csv", index=False, encoding="utf-8-sig")

    colors = {"Ai2020 paper": "#777777", "Thickness gauge": "#0072B2", "SEM cross-section": "#D55E00"}
    fig, axes = plt.subplots(2, 3, figsize=(16.8, 8.8), constrained_layout=True)
    for col, rate in enumerate(ga.RATES):
        for row_idx, charge in enumerate((True, False)):
            ax = axes[row_idx, col]
            direction = "Charge" if charge else "Discharge"
            obs = experiment[(rate, charge)]
            q_exp = priority.transferred_capacity(obs, qcell, charge)
            ax.plot(q_exp, obs["V"], color="black", lw=2.6, label="Experiment")
            for case in thickness3.CASES_UM:
                sim = simulations[(case, rate, charge)]
                result = adjusted[(adjusted.Case == case) & np.isclose(adjusted.C_rate, rate) & (adjusted.Direction == direction)].iloc[0]
                q_sim = priority.transferred_capacity(sim, qcell, charge)
                ax.plot(q_sim, sim["V"], color=colors[case], lw=1.8, label=f"{case}: {result.Capacity_error_mAh:+.1f} mAh")
            ax.set(title=f"{rate:g}C {direction}", xlabel="Transferred capacity [Ah]", ylabel="Voltage [V]")
            ax.grid(alpha=0.22)
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("Three thickness cases with remeasured electrode-mass ratios")
    fig.savefig(RESULTS / "mass_adjusted_capacity_voltage_curves.png", dpi=220)
    plt.close(fig)

    report = {
        "mass_mg": {"old": OLD_MASS_MG, "remeasured": NEW_MASS_MG},
        "mass_ratio": MASS_RATIO,
        "interpretation": "small measured coating-mass change applied as solid-inventory ratio; c_s,max and OCP endpoints held fixed for isolated sensitivity",
        "warning": "Qn/Qp change is not used to re-fit the OCP window in this ablation",
        "summary": summary.to_dict(orient="records"),
        "delta": comparison.to_dict(orient="records"),
    }
    (RESULTS / "mass_sensitivity_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print("\nSUMMARY")
    print(summary.to_string(index=False))
    print("\nDELTA")
    print(comparison[keys + [f"Delta_{m}" for m in metrics]].to_string(index=False))
    print("\nSaved to", RESULTS)


if __name__ == "__main__":
    main()
