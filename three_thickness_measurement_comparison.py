"""Compare three measured electrode-thickness interpretations.

The re-weighed electrodes were not materially different from the previous
measurement.  Therefore the primary comparison holds solid inventory per unit
area constant instead of changing c_s,max.  Active-material fraction scales as
1/L, the baseline inactive fraction is held constant, and porosity closes the
volume balance.  This represents the same coating loading at different measured
or sectioned thicknesses.
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
import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import kp_eis_recalculation_ablation as kpstudy
import priority_initial_state_protocol_recheck as priority


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260927_three_thickness_measurement_comparison"
RESULTS.mkdir(parents=True, exist_ok=True)

NOMINAL_CAPACITY_AH = 2.28
KP_APPLIED = kpstudy.KP_EIS

CASES_UM = {
    "Ai2020 paper": {"Ln": 76.5, "Lp": 68.0},
    "Thickness gauge": {"Ln": (155.0 - 9.0) / 2.0, "Lp": (138.0 - 13.0) / 2.0},
    "SEM cross-section": {"Ln": (172.04 - 9.0) / 2.0, "Lp": (138.93 - 13.0) / 2.0},
}


def copy_model(model: dict) -> dict:
    changed = dict(model)
    changed["params"] = model["params"].copy()
    return changed


def same_loading_geometry(model0: dict, ln_um: float, lp_um: float) -> dict:
    """Change L while preserving L*eps_s and baseline inactive fraction."""
    changed = copy_model(model0)
    p = changed["params"]
    updates = {}
    metadata = {}
    for electrode, new_um in (("Negative", ln_um), ("Positive", lp_um)):
        old_l = float(p[f"{electrode} electrode thickness [m]"])
        new_l = float(new_um) * 1e-6
        old_eps_s = float(p[f"{electrode} electrode active material volume fraction"])
        old_eps_e = float(p[f"{electrode} electrode porosity"])
        inactive = 1.0 - old_eps_s - old_eps_e
        new_eps_s = old_eps_s * old_l / new_l
        new_eps_e = 1.0 - inactive - new_eps_s
        if new_eps_e <= 0:
            raise ValueError(f"Non-positive {electrode} porosity: {new_eps_e}")
        updates.update(
            {
                f"{electrode} electrode thickness [m]": new_l,
                f"{electrode} electrode active material volume fraction": new_eps_s,
                f"{electrode} electrode porosity": new_eps_e,
            }
        )
        metadata[electrode.lower()] = {
            "old_L_um": old_l * 1e6,
            "new_L_um": new_um,
            "old_eps_s": old_eps_s,
            "new_eps_s": new_eps_s,
            "old_porosity": old_eps_e,
            "new_porosity": new_eps_e,
            "inactive_fraction_fixed": inactive,
            "inventory_ratio": new_l * new_eps_s / (old_l * old_eps_s),
        }
    p.update(updates, check_already_exists=False)
    changed["geometry_metadata"] = metadata
    return changed


def main() -> None:
    model0, stage, anode, cathode, experiment, _, qcell = selected.selected_inputs()
    model0 = kpstudy.with_kp(model0, KP_APPLIED)
    protocol, _, _ = priority.load_protocol(model0)
    pidx = protocol.set_index(["C_rate", "Direction"])

    models = {
        name: same_loading_geometry(model0, values["Ln"], values["Lp"])
        for name, values in CASES_UM.items()
    }

    rows = []
    simulations = {}
    for case, model in models.items():
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
                        "Ln_um": CASES_UM[case]["Ln"],
                        "Lp_um": CASES_UM[case]["Lp"],
                        "eps_s_n": float(model["params"]["Negative electrode active material volume fraction"]),
                        "eps_e_n": float(model["params"]["Negative electrode porosity"]),
                        "eps_s_p": float(model["params"]["Positive electrode active material volume fraction"]),
                        "eps_e_p": float(model["params"]["Positive electrode porosity"]),
                        "C_rate": rate,
                        "Direction": direction,
                        "measured_current_A": measured_current,
                        "initial_SOC": initial_soc,
                        **metric,
                        "Capacity_error_mAh": 1000.0
                        * (metric["Q_end_model_Ah"] - metric["Q_end_exp_Ah"]),
                    }
                )

    detail = pd.DataFrame(rows)
    detail.to_csv(RESULTS / "three_thickness_detail.csv", index=False, encoding="utf-8-sig")
    summary = detail.groupby(["Case", "Direction"], sort=False, as_index=False).agg(
        Mean_full_RMSE_mV=("Full_RMSE_mV", "mean"),
        Mean_center_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
        Capacity_RMSE_pct=("Capacity_error_pct", lambda x: float(np.sqrt(np.mean(x**2)))),
        Mean_capacity_error_mAh=("Capacity_error_mAh", "mean"),
        Mean_abs_capacity_error_mAh=("Capacity_error_mAh", lambda x: float(np.mean(np.abs(x)))),
    )
    summary.to_csv(RESULTS / "three_thickness_summary.csv", index=False, encoding="utf-8-sig")

    colors = {
        "Ai2020 paper": "#777777",
        "Thickness gauge": "#0072B2",
        "SEM cross-section": "#D55E00",
    }
    fig, axes = plt.subplots(2, 3, figsize=(16.8, 9.0), constrained_layout=True)
    for col, rate in enumerate(ga.RATES):
        for row_index, charge in enumerate((True, False)):
            ax = axes[row_index, col]
            direction = "Charge" if charge else "Discharge"
            obs = experiment[(rate, charge)]
            q_exp = priority.transferred_capacity(obs, qcell, charge)
            ax.plot(q_exp, obs["V"], color="black", lw=2.6, label="Experiment")
            for case in CASES_UM:
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
    axes[0, 0].legend(fontsize=6.8)
    fig.suptitle("Three electrode-thickness cases | same coating inventory")
    fig.savefig(RESULTS / "three_thickness_capacity_voltage_curves.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.8), constrained_layout=True)
    width = 0.25
    x = np.arange(len(ga.RATES))
    for case_index, case in enumerate(CASES_UM):
        offset = (case_index - 1) * width
        ch = detail[(detail.Case == case) & (detail.Direction == "Charge")].set_index("C_rate").loc[list(ga.RATES)]
        dis = detail[(detail.Case == case) & (detail.Direction == "Discharge")].set_index("C_rate").loc[list(ga.RATES)]
        axes[0].bar(x + offset, ch.Capacity_error_mAh, width, color=colors[case], label=case)
        axes[1].bar(x + offset, dis.Capacity_error_mAh, width, color=colors[case], label=case)
        axes[2].plot(ga.RATES, ch.Full_RMSE_mV, marker="o", color=colors[case], label=f"{case}, charge")
        axes[2].plot(ga.RATES, dis.Full_RMSE_mV, marker="s", ls="--", color=colors[case], label=f"{case}, discharge")
    for ax, title in zip(axes[:2], ("Charge capacity error", "Discharge capacity error")):
        ax.set_title(title)
        ax.set_ylabel("Model - experiment [mAh]")
        ax.set_xticks(x, [f"{r:g}C" for r in ga.RATES])
        ax.axhline(0, color="black", lw=0.8)
        ax.grid(axis="y", alpha=0.22)
    axes[2].set(title="Full-range voltage RMSE", xlabel="C-rate", ylabel="RMSE [mV]")
    axes[2].grid(alpha=0.22)
    axes[0].legend(fontsize=7)
    axes[2].legend(fontsize=6.2)
    fig.suptitle("Effect of Lp/Ln measurement method")
    fig.savefig(RESULTS / "three_thickness_metric_comparison.png", dpi=220)
    plt.close(fig)

    geometry = []
    for case, model in models.items():
        meta = model["geometry_metadata"]
        geometry.append(
            {
                "Case": case,
                "Ln_um": CASES_UM[case]["Ln"],
                "Lp_um": CASES_UM[case]["Lp"],
                "eps_s_n": meta["negative"]["new_eps_s"],
                "eps_e_n": meta["negative"]["new_porosity"],
                "eps_s_p": meta["positive"]["new_eps_s"],
                "eps_e_p": meta["positive"]["new_porosity"],
                "inactive_n": meta["negative"]["inactive_fraction_fixed"],
                "inactive_p": meta["positive"]["inactive_fraction_fixed"],
            }
        )
    pd.DataFrame(geometry).to_csv(RESULTS / "geometry_inputs.csv", index=False, encoding="utf-8-sig")

    report = {
        "method": (
            "Same measured coating inventory: c_s,max and L*epsilon_s fixed; "
            "baseline inactive fraction fixed; porosity closes volume balance."
        ),
        "collectors_um": {"Cu": 9.0, "Al": 13.0},
        "electrode_mass_recheck": "No material change from the previous measurement (user report)",
        "fixed": {
            "OCP_window": "current C/50-selected OCP/window",
            "initial_SOC": "preceding-rest voltage",
            "kp": KP_APPLIED,
            "kn": "legacy 7.40e-7 simulation placeholder; anode Rct anchor rejected",
            "Rn_um": 5.0,
            "Rp_um": 3.0,
            "Dsn_m2_s": 2.1e-14,
            "Dsp_m2_s": 4.4e-14,
            "negative_hysteresis": False,
            "positive_hysteresis": True,
        },
        "geometry": geometry,
        "summary": summary.to_dict(orient="records"),
        "detail": detail.to_dict(orient="records"),
    }
    (RESULTS / "three_thickness_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("\nGEOMETRY")
    print(pd.DataFrame(geometry).to_string(index=False))
    print("\nSUMMARY")
    print(summary.to_string(index=False))
    print("\nDETAIL")
    print(detail[["Case", "C_rate", "Direction", "Full_RMSE_mV", "Center10_70_RMSE_mV", "Capacity_error_mAh", "Capacity_error_pct"]].to_string(index=False))
    print("\nSaved to", RESULTS)


if __name__ == "__main__":
    main()
