"""Compare OCP candidates while sweeping negative-particle radius.

User-requested fixed transport parameters:
  Dsp = 4.4e-14 m2/s
  Dsn = 2.1e-14 m2/s
  electrolyte-phase Bruggeman coefficients b_n, b_p, b_s = 1.5
  Rp fixed at the Ai2020 value (3 um)
  Rn = 2, 3, 4 um

The stoichiometry window and effective electrode capacities are retained from
the prior mass-independent OCP comparison.  Only the held-out 0.5C/1C/2C
dynamic simulations are rerun, so qOCV metrics are invariant with Rn.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import mass_independent_nominal_anode_trial as mit
import ocp_grounded_candidate_search as common
import old_notebook_method_comparison as omc


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "results" / "260920_mass_independent_nominal_anode"
RESULTS = ROOT / "results" / "260920_ocp_rn_sweep_fixed_transport"
RESULTS.mkdir(parents=True, exist_ok=True)

DSN = 2.1e-14
DSP = 4.4e-14
RN_VALUES_UM = (2.0, 3.0, 4.0)
BRUGGEMAN = 1.5
BRUGGEMAN_KEYS = (
    "Negative electrode Bruggeman coefficient (electrolyte)",
    "Positive electrode Bruggeman coefficient (electrolyte)",
    "Separator Bruggeman coefficient (electrolyte)",
)


def candidate_details(bundle, base):
    measured = ehq.build_anode_detail(bundle, base)
    measured["source"] = "Harvested Enertech anode v1/v2"
    nominal = mit.nominal_anode_detail()
    cathode = ehq.build_cathode_detail(bundle, base)

    grid = nominal["grid"]
    measured_grid = measured["grid"]
    hybrid = dict(nominal)
    hybrid["lithiation"] = nominal["equilibrium"] + np.interp(
        grid, measured_grid, measured["lithiation"] - measured["equilibrium"]
    )
    hybrid["delithiation"] = nominal["equilibrium"] + np.interp(
        grid, measured_grid, measured["delithiation"] - measured["equilibrium"]
    )
    hybrid["source"] = "Ai2020 nominal equilibrium + harvested branch offsets"
    return measured, nominal, hybrid, cathode


def update_transport(model, rn_um: float, rp_m: float):
    changed = dict(model)
    changed["params"] = model["params"].copy()
    update = {
        "Negative particle diffusivity [m2.s-1]": DSN,
        "Positive particle diffusivity [m2.s-1]": DSP,
        "Negative particle radius [m]": rn_um * 1e-6,
        "Positive particle radius [m]": rp_m,
    }
    update.update({key: BRUGGEMAN for key in BRUGGEMAN_KEYS})
    changed["params"].update(update, check_already_exists=False)
    return changed


def main():
    bundle = ga.load_legacy_bundle(ga.get_paths(ROOT))
    base = ga.build_legacy_model_inputs(bundle)
    measured, nominal, hybrid, cathode = candidate_details(bundle, base)
    experiment = omc.load_old_dynamic_data(base["dynamic"], base["q_meas"])
    metrics = pd.read_csv(SOURCE / "mass_independent_ocp_candidate_metrics.csv")

    names = [
        "Previous measured-anode ±5% control",
        "Measured anode, mass-independent fit",
        "Ai2020 nominal anode, mass-independent fit",
        "Ai2020 nominal equilibrium + harvested hysteresis",
        "Ai2020 nominal anode, published stoichiometry window",
    ]
    rp_m = float(base["params"]["Positive particle radius [m]"])
    rows = []
    simulations = {}

    for name in names:
        metric_row = metrics[metrics.Candidate == name].iloc[0]
        effective_model = mit.model_for_effective_capacities(base, metric_row)
        stage = {key: float(metric_row[key]) for key in ("x0", "x100", "y100", "y0")}
        if name == "Ai2020 nominal equilibrium + harvested hysteresis":
            anode = common.scaled_detail(hybrid, 0.50)
            negative_hysteresis = True
        elif name.startswith("Ai2020"):
            anode = nominal
            negative_hysteresis = False
        else:
            anode = common.scaled_detail(measured, 0.50)
            negative_hysteresis = True
        cathode_dynamic = common.scaled_detail(cathode, 0.50)

        for rn_um in RN_VALUES_UM:
            model = update_transport(effective_model, rn_um, rp_m)
            for rate in ga.RATES:
                for charge in (True, False):
                    direction = "Charge" if charge else "Discharge"
                    print(f"{name} | Rn={rn_um:g} um | {rate:g}C | {direction}", flush=True)
                    try:
                        sim = ehq.run_dfn(
                            rate,
                            charge,
                            stage,
                            model,
                            anode,
                            cathode_dynamic,
                            negative_hysteresis,
                            True,
                        )
                        simulations[(name, rn_um, rate, charge)] = sim
                        result = omc.old_time_metrics(experiment[(rate, charge)], sim)
                        capacity_error = result["end_time_error_min"] * rate * 2.28 / 60.0
                        status = "ok"
                    except Exception as exc:
                        print(f"FAILED: {exc}", flush=True)
                        result = {
                            "MAE_all_mV": np.nan,
                            "MAE_SOC10_70_mV": np.nan,
                            "N_all": 0,
                            "N_SOC10_70": 0,
                            "end_time_error_min": np.nan,
                        }
                        capacity_error = np.nan
                        status = f"failed: {type(exc).__name__}: {exc}"
                    rows.append(
                        {
                            "Candidate": name,
                            "Rn_um": rn_um,
                            "Rp_um": rp_m * 1e6,
                            "Dsn_m2_s": DSN,
                            "Dsp_m2_s": DSP,
                            "Bruggeman": BRUGGEMAN,
                            "C_rate": rate,
                            "Direction": direction,
                            **result,
                            "capacity_error_Ah": capacity_error,
                            "status": status,
                        }
                    )

    validation = pd.DataFrame(rows)
    summary = validation.groupby(["Candidate", "Rn_um"], as_index=False).agg(
        Dynamic_MAE_SOC10_70_mV=("MAE_SOC10_70_mV", "mean"),
        Charge_MAE_SOC10_70_mV=("MAE_SOC10_70_mV", lambda x: float(x[validation.loc[x.index, "Direction"] == "Charge"].mean())),
        Discharge_MAE_SOC10_70_mV=("MAE_SOC10_70_mV", lambda x: float(x[validation.loc[x.index, "Direction"] == "Discharge"].mean())),
        Mean_abs_capacity_error_mAh=("capacity_error_Ah", lambda x: float(np.nanmean(np.abs(x))) * 1000.0),
        Max_abs_capacity_error_mAh=("capacity_error_Ah", lambda x: float(np.nanmax(np.abs(x))) * 1000.0),
        Failed_runs=("status", lambda x: int(np.sum(x != "ok"))),
    )
    qcols = [
        "Candidate",
        "qOCV_MAE_2_98_mV",
        "qOCV_RMSE_2_98_mV",
        "qOCV_max_endpoint_abs_mV",
    ]
    summary = summary.merge(metrics[qcols], on="Candidate", how="left")
    summary["tau_n_s"] = (summary.Rn_um * 1e-6) ** 2 / DSN
    summary["relative_negative_surface_area_vs_Rn3"] = 3.0 / summary.Rn_um
    validation.to_csv(RESULTS / "ocp_rn_dynamic_validation.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(RESULTS / "ocp_rn_summary.csv", index=False, encoding="utf-8-sig")

    colors = {
        names[0]: "#7f7f7f",
        names[1]: "#1f77b4",
        names[2]: "#d62728",
        names[3]: "#ff7f0e",
        names[4]: "#2ca02c",
    }
    labels = {
        names[0]: "Measured ±5%",
        names[1]: "Measured free",
        names[2]: "Ai nominal",
        names[3]: "Ai nominal + hyst",
        names[4]: "Ai published window",
    }
    fig, axes = plt.subplots(2, 2, figsize=(14.5, 10.0), constrained_layout=True)
    for name in names:
        frame = summary[summary.Candidate == name].sort_values("Rn_um")
        axes[0, 0].plot(frame.Rn_um, frame.Dynamic_MAE_SOC10_70_mV, "o-", lw=2, color=colors[name], label=labels[name])
        axes[0, 1].plot(frame.Rn_um, frame.Mean_abs_capacity_error_mAh, "o-", lw=2, color=colors[name], label=labels[name])
        axes[1, 0].plot(frame.Rn_um, frame.Charge_MAE_SOC10_70_mV, "o-", lw=2, color=colors[name], label=labels[name])
        axes[1, 1].plot(frame.Rn_um, frame.Discharge_MAE_SOC10_70_mV, "o-", lw=2, color=colors[name], label=labels[name])
    axes[0, 0].set(title="Held-out dynamic voltage error", ylabel="Mean MAE, SOC 10-70% (mV)")
    axes[0, 1].set(title="Termination-capacity error", ylabel="Mean absolute error (mAh)")
    axes[1, 0].set(title="Charge voltage error", ylabel="Mean MAE, SOC 10-70% (mV)")
    axes[1, 1].set(title="Discharge voltage error", ylabel="Mean MAE, SOC 10-70% (mV)")
    for ax in axes.flat:
        ax.set_xlabel("Negative particle radius Rn (um)")
        ax.set_xticks(RN_VALUES_UM)
        ax.grid(alpha=0.25)
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("OCP x Rn sweep | Dsn=2.1e-14, Dsp=4.4e-14 m2/s | Rp=3 um | electrolyte Bruggeman=1.5", fontsize=13)
    fig.savefig(RESULTS / "ocp_rn_sweep_summary.png", dpi=210)
    plt.close(fig)

    admissible = summary[(summary.qOCV_max_endpoint_abs_mV <= 10.0) & (summary.Failed_runs == 0)].copy()
    # Keep the selection transparent: among acceptable qOCV endpoints, rank
    # voltage MAE first and report capacity error independently.
    selected = admissible.sort_values(
        ["Dynamic_MAE_SOC10_70_mV", "Mean_abs_capacity_error_mAh"]
    ).iloc[0]
    selected_name = str(selected.Candidate)

    fig, axes = plt.subplots(2, 3, figsize=(16.0, 8.5), constrained_layout=True)
    for col, rate in enumerate(ga.RATES):
        for row_index, charge in enumerate((True, False)):
            ax = axes[row_index, col]
            exp = experiment[(rate, charge)]
            ax.plot(exp["t_min"], exp["V"], color="black", lw=2.3, label="Experiment")
            for rn_um, color in zip(RN_VALUES_UM, ("#1f77b4", "#ff7f0e", "#2ca02c")):
                sim = simulations.get((selected_name, rn_um, rate, charge))
                if sim is not None:
                    ax.plot(sim["t_min"], sim["V"], lw=1.7, color=color, label=f"Rn={rn_um:g} um")
            ax.set_title(f"{rate:g}C {'Charge' if charge else 'Discharge'}")
            ax.set_xlabel("Time (min)")
            ax.set_ylabel("Voltage (V)")
            ax.grid(alpha=0.2)
    axes[0, 0].legend(fontsize=8)
    fig.suptitle(f"Best voltage candidate across Rn sweep: {selected_name}", fontsize=13)
    fig.savefig(RESULTS / "best_ocp_rn_dynamic_curves.png", dpi=210)
    plt.close(fig)

    best_by_ocp = summary.loc[summary.groupby("Candidate")["Dynamic_MAE_SOC10_70_mV"].idxmin()].copy()
    best_by_ocp = best_by_ocp.sort_values("Dynamic_MAE_SOC10_70_mV")
    best_by_ocp.to_csv(RESULTS / "best_rn_by_ocp.csv", index=False, encoding="utf-8-sig")

    report = {
        "fixed_parameters": {
            "Dsn_m2_s": DSN,
            "Dsp_m2_s": DSP,
            "Rp_um": rp_m * 1e6,
            "Bruggeman_electrolyte_n_p_s": BRUGGEMAN,
            "Bruggeman_solid_phase": "retained Ai2020 value (0.0)",
            "Rn_um": list(RN_VALUES_UM),
        },
        "selection_rule": "qOCV endpoint max <= 10 mV and no failed run; then minimum held-out dynamic MAE. Capacity errors are reported separately, not hidden in an arbitrary weight.",
        "selected_voltage_candidate": selected.to_dict(),
        "best_Rn_by_OCP": best_by_ocp.to_dict(orient="records"),
        "notes": [
            "qOCV is invariant to Dsn, Dsp, Bruggeman and particle radius in this equilibrium comparison.",
            "Changing Rn with fixed Dsn changes both Rn^2/Dsn and the active interfacial area 3*eps_s/Rn.",
            "The direct Ai2020 published window is retained as a negative control but excluded by its endpoint error.",
            "The 260918 new GITT data are excluded.",
        ],
    }
    (RESULTS / "ocp_rn_sweep_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("\nBEST BY OCP")
    print(best_by_ocp[["Candidate", "Rn_um", "Dynamic_MAE_SOC10_70_mV", "Mean_abs_capacity_error_mAh", "Max_abs_capacity_error_mAh"]].to_string(index=False))
    print("\nSELECTED")
    print(selected.to_string())
    print(f"\nSaved to {RESULTS}")


if __name__ == "__main__":
    main()
