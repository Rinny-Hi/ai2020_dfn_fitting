"""Validate the C/50-selected OCP on 0.5C/1C/2C charge and discharge data."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import c50_ocp_source_combination_study as c50
import current_thickness_estimate_application as thickness
import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import mass_independent_nominal_anode_trial as mit
import nominal_radius_current_configuration as current
import old_notebook_method_comparison as omc
import priority_initial_state_protocol_recheck as priority
import restored_brugg_de_csmax_audit as audit
import stage2_seven_parameter_sensitivity as stage2


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "results" / "260921_c50_ocp_source_combinations"
RESULTS = ROOT / "results" / "260921_c50_selected_ocp_dynamic_validation"
RESULTS.mkdir(parents=True, exist_ok=True)


def selected_inputs():
    bundle = ga.load_legacy_bundle(ga.get_paths(ROOT))
    base = audit.build_base(bundle, 1e-4, audit.CSN_CORRECTED)
    summary = pd.read_csv(SOURCE / "ocp_source_combination_summary.csv")
    selected = summary[(summary.anode_source == "nominal") & (summary.cathode_source == "old")].iloc[0]
    decision = json.loads((SOURCE / "decision.json").read_text(encoding="utf-8"))
    qcell = float(decision["target"]["qcell_Ah"])

    row = pd.Series(
        {
            "Q_n_Ah": float(selected.Qn_Ah),
            "Q_p_Ah": float(selected.Qp_Ah),
            "delta_x": float(selected.delta_x),
            "delta_y": float(selected.delta_y),
        }
    )
    model = mit.model_for_effective_capacities(base, row)
    model = audit.set_dynamic_parameters(model, current.RN_UM, current.RP_UM, 1e-4)
    model = thickness.set_thickness_preserve_capacity(
        model, current.LN_UM * 1e-6, current.LP_UM * 1e-6
    )
    model = current.apply_current_kinetics_and_bruggeman(model)
    model["q_meas"] = qcell

    anode = mit.nominal_anode_detail()
    old_cathode = ehq.build_cathode_detail(bundle, base)
    # The C/50 branch calibration selected the raw old-cathode gap (scale ~= 1).
    cathode = {
        key: value.copy() if isinstance(value, np.ndarray) else value
        for key, value in old_cathode.items()
    }
    stage = {
        "x0": float(selected.x0),
        "x100": float(selected.x100),
        "y100": float(selected.y100),
        "y0": float(selected.y0),
    }
    soc = np.linspace(0.0, 1.0, 1001)
    model["soc"] = soc
    model["v_qocv"] = c50.predict(
        np.array([stage["x0"], stage["x100"], stage["y100"], stage["y0"]]),
        soc,
        anode,
        cathode,
        "equilibrium",
    )
    experiment = omc.load_old_dynamic_data(base["dynamic"], qcell)
    return model, stage, anode, cathode, experiment, selected, qcell


def current_baseline(qcell: float):
    _, model, stage, anode, cathode, _, _ = stage2.build_current_inputs()
    model = dict(model)
    model["q_meas"] = qcell
    bundle = ga.load_legacy_bundle(ga.get_paths(ROOT))
    experiment = omc.load_old_dynamic_data(bundle["dynamic"], qcell)
    return model, stage, anode, cathode, experiment


def run_all(
    name,
    model,
    stage,
    anode,
    cathode,
    negative_hysteresis,
    positive_hysteresis,
    protocol=None,
    history=False,
):
    runs = {}
    protocol_index = None if protocol is None else protocol.set_index(["C_rate", "Direction"])
    for rate in ga.RATES:
        for charge in (True, False):
            direction = "Charge" if charge else "Discharge"
            local_stage = stage
            neg_h = None
            pos_h = None
            if protocol_index is not None:
                z = float(protocol_index.loc[(rate, direction), "qOCV_equivalent_initial_SOC_pct"]) / 100.0
                local_stage = priority.stage_at_soc(stage, z, charge)
            if history:
                neg_h = 1.0 if charge else -1.0
                pos_h = -1.0 if charge else 1.0
            print(f"{name} | {rate:g}C {direction}", flush=True)
            try:
                runs[(rate, charge)] = ehq.run_dfn(
                    rate,
                    charge,
                    local_stage,
                    model,
                    anode,
                    cathode,
                    negative_hysteresis,
                    positive_hysteresis,
                    negative_initial_h=neg_h,
                    positive_initial_h=pos_h,
                )
            except Exception as exc:
                print(f"FAILED: {exc}", flush=True)
                runs[(rate, charge)] = None
    return runs


def score(scenarios, experiment, qcell):
    rows = []
    for name, runs in scenarios.items():
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                sim = runs[(rate, charge)]
                if sim is None:
                    rows.append(
                        {
                            "Scenario": name,
                            "C_rate": rate,
                            "Direction": direction,
                            "Status": "failed",
                        }
                    )
                    continue
                metrics = priority.curve_metrics(
                    experiment[(rate, charge)], sim, qcell, charge
                )
                rows.append(
                    {
                        "Scenario": name,
                        "C_rate": rate,
                        "Direction": direction,
                        "Status": "ok",
                        **metrics,
                        "Capacity_error_mAh": 1000.0
                        * (metrics["Q_end_model_Ah"] - metrics["Q_end_exp_Ah"]),
                    }
                )
    return pd.DataFrame(rows)


def main():
    model, stage, anode, cathode, experiment, selected, qcell = selected_inputs()
    baseline_model, baseline_stage, baseline_anode, baseline_cathode, _ = current_baseline(qcell)

    protocol, _, _ = priority.load_protocol(model)
    scenarios = {
        "Previous current OCP": run_all(
            "Previous current OCP",
            baseline_model,
            baseline_stage,
            baseline_anode,
            baseline_cathode,
            True,
            True,
        ),
        "C50 OCP endpoint start": run_all(
            "C50 OCP endpoint start",
            model,
            stage,
            anode,
            cathode,
            False,
            True,
        ),
        "C50 OCP rest/history start": run_all(
            "C50 OCP rest/history start",
            model,
            stage,
            anode,
            cathode,
            False,
            True,
            protocol=protocol,
            history=True,
        ),
    }
    detail = score(scenarios, experiment, qcell)
    detail.to_csv(RESULTS / "dynamic_condition_metrics.csv", index=False, encoding="utf-8-sig")
    valid = detail[detail.Status == "ok"].copy()
    summary = valid.groupby(["Scenario", "Direction"], as_index=False).agg(
        Mean_full_RMSE_mV=("Full_RMSE_mV", "mean"),
        Mean_center10_70_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
        Mean_center10_70_MAE_mV=("Center10_70_MAE_mV", "mean"),
        Capacity_RMSE_pct=("Capacity_error_pct", lambda x: float(np.sqrt(np.mean(x**2)))),
        Mean_abs_capacity_error_mAh=("Capacity_error_mAh", lambda x: float(np.mean(np.abs(x)))),
    )
    summary.to_csv(RESULTS / "dynamic_summary.csv", index=False, encoding="utf-8-sig")

    colors = {
        "Previous current OCP": "#777777",
        "C50 OCP endpoint start": "#0072B2",
        "C50 OCP rest/history start": "#D55E00",
    }
    fig, axes = plt.subplots(2, 3, figsize=(16.4, 8.8), constrained_layout=True)
    for col, rate in enumerate(ga.RATES):
        for row, charge in enumerate((True, False)):
            ax = axes[row, col]
            direction = "Charge" if charge else "Discharge"
            exp = experiment[(rate, charge)]
            q_exp = priority.transferred_capacity(exp, qcell, charge)
            ax.plot(q_exp, exp["V"], color="black", lw=2.5, label="Experiment")
            for name, runs in scenarios.items():
                sim = runs[(rate, charge)]
                if sim is None:
                    continue
                metric = valid[
                    (valid.Scenario == name)
                    & np.isclose(valid.C_rate, rate)
                    & (valid.Direction == direction)
                ].iloc[0]
                q_sim = priority.transferred_capacity(sim, qcell, charge)
                ax.plot(
                    q_sim,
                    sim["V"],
                    color=colors[name],
                    lw=1.8,
                    label=f"{name}: {metric.Center10_70_RMSE_mV:.1f} mV",
                )
            ax.set_title(f"{rate:g}C {direction}")
            ax.set_xlabel("Transferred capacity [Ah]")
            ax.set_ylabel("Voltage [V]")
            ax.grid(alpha=0.22)
    axes[0, 0].legend(fontsize=7.2)
    fig.suptitle("0.5C / 1C / 2C validation of C/50-selected OCP")
    fig.savefig(RESULTS / "dynamic_curve_comparison.png", dpi=220)
    plt.close(fig)

    result = {
        "selected_OCP": "Ai2020 nominal anode + old cathode GITT",
        "stage": stage,
        "Qn_Ah": float(selected.Qn_Ah),
        "Qp_Ah": float(selected.Qp_Ah),
        "qcell_C50_Ah": qcell,
        "fixed_dynamic_parameters": {
            "Rn_um": current.RN_UM,
            "Rp_um": current.RP_UM,
            "Dsn_m2_s": 2.1e-14,
            "Dsp_m2_s": 4.4e-14,
            "kn_prefactor": current.KN_PREF,
            "kp_prefactor": current.KP_PREF,
            "Bruggeman_n_p_s": [current.BRUGG_N, current.BRUGG_P, current.BRUGG_S],
            "Ln_um": current.LN_UM,
            "Lp_um": current.LP_UM,
            "Cu_collector_um": 10.0,
            "Al_collector_um": 15.0,
            "collector_source": "SPB655060 teardown paper; user-measured 9/13 um retained only as sensitivity evidence",
            "effective_csn_max_mol_m3": float(
                model["params"]["Maximum concentration in negative electrode [mol.m-3]"]
            ),
            "effective_csp_max_mol_m3": float(
                model["params"]["Maximum concentration in positive electrode [mol.m-3]"]
            ),
            "negative_hysteresis": False,
            "positive_hysteresis": True,
        },
        "summary": summary.to_dict(orient="records"),
    }
    (RESULTS / "dynamic_validation.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("\nDETAIL")
    print(detail.to_string(index=False))
    print("\nSUMMARY")
    print(summary.to_string(index=False))
    print("Saved to", RESULTS)


if __name__ == "__main__":
    main()
