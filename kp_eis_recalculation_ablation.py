"""Compare the current kp with the EIS-recalculated cathode kp."""

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
import nominal_radius_current_configuration as current
import priority_initial_state_protocol_recheck as priority


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260924_kp_eis_recalculation"
RESULTS.mkdir(parents=True, exist_ok=True)
NOMINAL_CAPACITY_AH = 2.28
KP_EIS = 4.18e-7
KP_EIS_MODEL_CONSISTENT = 4.32e-7
CASES = {
    "Current kp 3.12e-7": current.KP_PREF,
    "EIS recalculated kp 4.18e-7": KP_EIS,
    "EIS model-consistent kp 4.32e-7": KP_EIS_MODEL_CONSISTENT,
}


def with_kp(model0: dict, kp: float) -> dict:
    changed = dict(model0)
    changed["params"] = model0["params"].copy()
    key = "Positive electrode exchange-current density [A.m-2]"
    changed["params"].update(
        {key: ga.scale_parameter_function(model0["params"][key], kp / current.KP_PREF)},
        check_already_exists=False,
    )
    return changed


def main() -> None:
    model0, stage, anode, cathode, experiment, _, qcell = selected.selected_inputs()
    protocol, _, _ = priority.load_protocol(model0)
    pidx = protocol.set_index(["C_rate", "Direction"])
    rows = []
    simulations = {}

    for case, kp in CASES.items():
        model = with_kp(model0, kp)
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                rest_v = float(pidx.loc[(rate, direction), "Rest_end_V"])
                initial_soc = priority.monotone_voltage_inverse(
                    model["v_qocv"], model["soc"], rest_v
                )
                local_stage = priority.stage_at_soc(stage, initial_soc, charge)
                measured_current = float(
                    np.nanmedian(np.abs(experiment[(rate, charge)]["I_A"]))
                )
                sim_rate = measured_current / NOMINAL_CAPACITY_AH
                print(f"{case} | {rate:g}C {direction}", flush=True)
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
                metric = priority.curve_metrics(
                    experiment[(rate, charge)], sim, qcell, charge
                )
                rows.append(
                    {
                        "Case": case,
                        "kp": kp,
                        "C_rate": rate,
                        "Direction": direction,
                        "measured_current_A": measured_current,
                        "initial_SOC": initial_soc,
                        **metric,
                        "Capacity_error_mAh": 1000
                        * (metric["Q_end_model_Ah"] - metric["Q_end_exp_Ah"]),
                    }
                )

    detail = pd.DataFrame(rows)
    detail.to_csv(RESULTS / "kp_recalculation_detail.csv", index=False, encoding="utf-8-sig")
    summary = detail.groupby(["Case", "Direction"], sort=False, as_index=False).agg(
        Mean_full_RMSE_mV=("Full_RMSE_mV", "mean"),
        Mean_center_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
        Capacity_RMSE_pct=("Capacity_error_pct", lambda x: float(np.sqrt(np.mean(x**2)))),
        Mean_capacity_error_mAh=("Capacity_error_mAh", "mean"),
        Mean_abs_capacity_error_mAh=("Capacity_error_mAh", lambda x: float(np.mean(np.abs(x)))),
    )
    summary.to_csv(RESULTS / "kp_recalculation_summary.csv", index=False, encoding="utf-8-sig")

    colors = {
        "Current kp 3.12e-7": "#777777",
        "EIS recalculated kp 4.18e-7": "#D55E00",
        "EIS model-consistent kp 4.32e-7": "#0072B2",
    }
    fig, axes = plt.subplots(2, 3, figsize=(16.5, 8.8), constrained_layout=True)
    for col, rate in enumerate(ga.RATES):
        for row_index, charge in enumerate((True, False)):
            ax = axes[row_index, col]
            direction = "Charge" if charge else "Discharge"
            obs = experiment[(rate, charge)]
            q_exp = priority.transferred_capacity(obs, qcell, charge)
            ax.plot(q_exp, obs["V"], color="black", lw=2.5, label="Experiment")
            for case in CASES:
                sim = simulations[(case, rate, charge)]
                metric = detail[
                    (detail.Case == case)
                    & np.isclose(detail.C_rate, rate)
                    & (detail.Direction == direction)
                ].iloc[0]
                q_sim = priority.transferred_capacity(sim, qcell, charge)
                ax.plot(
                    q_sim,
                    sim["V"],
                    color=colors[case],
                    lw=1.8,
                    label=(
                        f"{case}: {metric.Full_RMSE_mV:.1f} mV, "
                        f"{metric.Capacity_error_mAh:+.1f} mAh"
                    ),
                )
            ax.set(
                title=f"{rate:g}C {direction}",
                xlabel="Transferred capacity [Ah]",
                ylabel="Voltage [V]",
            )
            ax.grid(alpha=0.22)
    axes[0, 0].legend(fontsize=6.7)
    fig.suptitle("Effect of EIS-recalculated cathode kp")
    fig.savefig(RESULTS / "kp_recalculation_dynamic_curves.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(15.8, 4.7), constrained_layout=True)
    width = 0.25
    x = np.arange(len(ga.RATES))
    for case_index, case in enumerate(CASES):
        offset = (case_index - (len(CASES) - 1) / 2) * width
        charge = detail[(detail.Case == case) & (detail.Direction == "Charge")].set_index("C_rate").loc[list(ga.RATES)]
        discharge = detail[(detail.Case == case) & (detail.Direction == "Discharge")].set_index("C_rate").loc[list(ga.RATES)]
        axes[0].bar(x + offset, charge.Capacity_error_mAh, width, color=colors[case], label=case)
        axes[1].bar(x + offset, discharge.Capacity_error_mAh, width, color=colors[case], label=case)
        axes[2].plot(ga.RATES, charge.Full_RMSE_mV, marker="o", color=colors[case], label=f"{case}, charge")
        axes[2].plot(ga.RATES, discharge.Full_RMSE_mV, marker="s", ls="--", color=colors[case], label=f"{case}, discharge")
    axes[0].set(title="Charge capacity error", ylabel="Model - experiment [mAh]")
    axes[1].set(title="Discharge capacity error", ylabel="Model - experiment [mAh]")
    for ax in axes[:2]:
        ax.set_xticks(x, [f"{r:g}C" for r in ga.RATES])
        ax.axhline(0, color="black", lw=0.8)
        ax.grid(axis="y", alpha=0.22)
    axes[2].set(title="Full-range voltage RMSE", xlabel="C-rate", ylabel="RMSE [mV]")
    axes[2].grid(alpha=0.22)
    axes[0].legend(fontsize=7)
    axes[2].legend(fontsize=6.5)
    fig.suptitle("Capacity and voltage changes from kp update")
    fig.savefig(RESULTS / "kp_recalculation_metric_changes.png", dpi=220)
    plt.close(fig)

    base = detail[detail.Case == "Current kp 3.12e-7"].set_index(["C_rate", "Direction"])
    new = detail[detail.Case == "EIS recalculated kp 4.18e-7"].set_index(["C_rate", "Direction"])
    delta = new[["Full_RMSE_mV", "Center10_70_RMSE_mV", "Capacity_error_mAh"]] - base[["Full_RMSE_mV", "Center10_70_RMSE_mV", "Capacity_error_mAh"]]
    delta = delta.reset_index().rename(
        columns={
            "Full_RMSE_mV": "Delta_full_RMSE_mV",
            "Center10_70_RMSE_mV": "Delta_center_RMSE_mV",
            "Capacity_error_mAh": "Delta_capacity_error_mAh",
        }
    )
    delta.to_csv(RESULTS / "kp_recalculation_delta.csv", index=False, encoding="utf-8-sig")

    report = {
        "current_kp": current.KP_PREF,
        "new_kp": KP_EIS,
        "new_kp_reexpressed_for_current_model_c_s_max": KP_EIS_MODEL_CONSISTENT,
        "new_over_current_ratio": KP_EIS / current.KP_PREF,
        "fixed_parameters": {
            "kn": current.KN_PREF,
            "Rn_um": current.RN_UM,
            "Rp_um": current.RP_UM,
            "Ln_um": current.LN_UM,
            "Lp_um": current.LP_UM,
            "positive_active_material_volume_fraction_in_model": float(
                model0["params"]["Positive electrode active material volume fraction"]
            ),
            "positive_max_concentration_mol_m3": float(
                model0["params"]["Maximum concentration in positive electrode [mol.m-3]"]
            ),
        },
        "protocol": "current selected OCP; measured current; preceding-rest initial SOC; negative hysteresis off; positive hysteresis on",
        "summary": summary.to_dict(orient="records"),
        "delta": delta.to_dict(orient="records"),
    }
    (RESULTS / "kp_recalculation_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("\nDETAIL")
    print(detail.to_string(index=False))
    print("\nSUMMARY")
    print(summary.to_string(index=False))
    print("\nDELTA")
    print(delta.to_string(index=False))


if __name__ == "__main__":
    main()
