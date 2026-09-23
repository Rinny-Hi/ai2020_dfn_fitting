"""Ablate reverting experimental kinetic prefactors to Ai2020 originals."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybamm

import c50_selected_ocp_dynamic_validation as selected
import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import nominal_radius_current_configuration as current
import priority_initial_state_protocol_recheck as priority


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260923_kinetics_original_reversion"
RESULTS.mkdir(parents=True, exist_ok=True)
AI2020_PREF = 1.0e-11 * 96485.33212
NOMINAL_CAPACITY_AH = 2.28

CASES = {
    "Current EIS kn/kp": (current.KN_PREF, current.KP_PREF),
    "Ai2020 kn + current kp": (AI2020_PREF, current.KP_PREF),
    "Current kn + Ai2020 kp": (current.KN_PREF, AI2020_PREF),
    "Ai2020 kn/kp both": (AI2020_PREF, AI2020_PREF),
}


def with_prefactors(model: dict, kn: float, kp: float) -> dict:
    changed = dict(model)
    changed["params"] = model["params"].copy()
    raw = pybamm.ParameterValues("Ai2020")
    changed["params"].update(
        {
            "Negative electrode exchange-current density [A.m-2]": ga.scale_parameter_function(
                raw["Negative electrode exchange-current density [A.m-2]"], kn / AI2020_PREF
            ),
            "Positive electrode exchange-current density [A.m-2]": ga.scale_parameter_function(
                raw["Positive electrode exchange-current density [A.m-2]"], kp / AI2020_PREF
            ),
        },
        check_already_exists=False,
    )
    return changed


def main() -> None:
    model0, stage, anode, cathode, experiment, _, qcell = selected.selected_inputs()
    protocol, _, _ = priority.load_protocol(model0)
    protocol_index = protocol.set_index(["C_rate", "Direction"])
    rows = []

    for case, (kn, kp) in CASES.items():
        model = with_prefactors(model0, kn, kp)
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                rest_v = float(protocol_index.loc[(rate, direction), "Rest_end_V"])
                initial_soc = priority.monotone_voltage_inverse(
                    model["v_qocv"], model["soc"], rest_v
                )
                local_stage = priority.stage_at_soc(stage, initial_soc, charge)
                measured_current = float(
                    np.nanmedian(np.abs(experiment[(rate, charge)]["I_A"]))
                )
                sim_rate = measured_current / NOMINAL_CAPACITY_AH
                print(f"{case} | {rate:g}C {direction}", flush=True)
                simulation = ehq.run_dfn(
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
                metric = priority.curve_metrics(
                    experiment[(rate, charge)], simulation, qcell, charge
                )
                rows.append(
                    {
                        "Case": case,
                        "kn": kn,
                        "kp": kp,
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
    detail.to_csv(RESULTS / "kinetics_reversion_detail.csv", index=False, encoding="utf-8-sig")
    summary = detail.groupby(["Case", "Direction"], sort=False, as_index=False).agg(
        mean_full_RMSE_mV=("Full_RMSE_mV", "mean"),
        mean_center_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
        capacity_RMSE_pct=("Capacity_error_pct", lambda x: float(np.sqrt(np.mean(x**2)))),
        mean_capacity_error_mAh=("Capacity_error_mAh", "mean"),
        mean_abs_capacity_error_mAh=("Capacity_error_mAh", lambda x: float(np.mean(np.abs(x)))),
    )
    summary.to_csv(RESULTS / "kinetics_reversion_summary.csv", index=False, encoding="utf-8-sig")

    cases = list(CASES)
    x = np.arange(len(cases))
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)
    width = 0.36
    for ax, metric, ylabel in (
        (axes[0], "mean_full_RMSE_mV", "Full-range voltage RMSE [mV]"),
        (axes[1], "mean_center_RMSE_mV", "10–70% voltage RMSE [mV]"),
        (axes[2], "mean_abs_capacity_error_mAh", "Mean absolute capacity error [mAh]"),
    ):
        charge = summary[summary.Direction == "Charge"].set_index("Case").loc[cases]
        discharge = summary[summary.Direction == "Discharge"].set_index("Case").loc[cases]
        ax.bar(x - width / 2, charge[metric], width, label="Charge", color="#D55E00")
        ax.bar(x + width / 2, discharge[metric], width, label="Discharge", color="#0072B2")
        ax.set_xticks(x, ["Current", "Ai n", "Ai p", "Ai n+p"], rotation=15)
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", alpha=0.22)
    axes[0].legend()
    fig.suptitle("Reverting kinetic prefactors to Ai2020 originals")
    fig.savefig(RESULTS / "kinetics_reversion_comparison.png", dpi=220)
    plt.close(fig)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
