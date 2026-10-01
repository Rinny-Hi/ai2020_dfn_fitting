"""Plot pre-fit time-voltage curves and ablate nominal/current k and Ds pairs."""

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
import priority_initial_state_protocol_recheck as priority


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260921_prefit_kinetics_diffusion_ablation"
RESULTS.mkdir(parents=True, exist_ok=True)

K_N = "Negative electrode exchange-current density [A.m-2]"
K_P = "Positive electrode exchange-current density [A.m-2]"
D_N = "Negative particle diffusivity [m2.s-1]"
D_P = "Positive particle diffusivity [m2.s-1]"


def changed_model(model: dict, updates: dict) -> dict:
    changed = dict(model)
    changed["params"] = model["params"].copy()
    changed["params"].update(updates, check_already_exists=False)
    return changed


def main() -> None:
    base, stage, anode, cathode, experiment, _, qcell = selected.selected_inputs()
    raw = pybamm.ParameterValues("Ai2020")

    current_kn = base["params"][K_N]
    current_kp = base["params"][K_P]
    current_dn = base["params"][D_N]
    current_dp = base["params"][D_P]
    nominal_kn = raw[K_N]
    nominal_kp = raw[K_P]
    nominal_dn = raw[D_N]
    nominal_dp = raw[D_P]

    scenarios = {
        "Current baseline": base,
        "k: nominal kn + nominal kp": changed_model(base, {K_N: nominal_kn, K_P: nominal_kp}),
        "k: nominal kn + current kp": changed_model(base, {K_N: nominal_kn, K_P: current_kp}),
        "k: current kn + nominal kp": changed_model(base, {K_N: current_kn, K_P: nominal_kp}),
        "Ds: nominal Dsn + nominal Dsp": changed_model(base, {D_N: nominal_dn, D_P: nominal_dp}),
        "Ds: nominal Dsn + current Dsp": changed_model(base, {D_N: nominal_dn, D_P: current_dp}),
        "Ds: current Dsn + nominal Dsp": changed_model(base, {D_N: current_dn, D_P: nominal_dp}),
    }

    rows = []
    simulations = {}
    for name, model in scenarios.items():
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                print(f"{name} | {rate:g}C {direction}", flush=True)
                sim = ehq.run_dfn(
                    rate, charge, stage, model, anode, cathode, False, True
                )
                simulations[(name, rate, charge)] = sim
                metric = priority.curve_metrics(
                    experiment[(rate, charge)], sim, qcell, charge
                )
                rows.append(
                    {
                        "Scenario": name,
                        "C_rate": rate,
                        "Direction": direction,
                        **metric,
                        "Capacity_error_mAh": 1000
                        * (metric["Q_end_model_Ah"] - metric["Q_end_exp_Ah"]),
                    }
                )
    detail = pd.DataFrame(rows)
    summary = detail.groupby(["Scenario", "Direction"], as_index=False, sort=False).agg(
        mean_center_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
        mean_capacity_error_mAh=("Capacity_error_mAh", "mean"),
        mean_abs_capacity_error_mAh=("Capacity_error_mAh", lambda x: float(np.mean(np.abs(x)))),
        max_abs_capacity_error_mAh=("Capacity_error_mAh", lambda x: float(np.max(np.abs(x)))),
    )
    detail.to_csv(RESULTS / "ablation_condition_metrics.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(RESULTS / "ablation_summary.csv", index=False, encoding="utf-8-sig")

    # Current pre-fit time-voltage curves.
    fig, axes = plt.subplots(2, 3, figsize=(16.0, 8.5), constrained_layout=True)
    for col, rate in enumerate(ga.RATES):
        for row, charge in enumerate((True, False)):
            ax = axes[row, col]
            direction = "Charge" if charge else "Discharge"
            exp = experiment[(rate, charge)]
            sim = simulations[("Current baseline", rate, charge)]
            metric = detail[
                (detail.Scenario == "Current baseline")
                & np.isclose(detail.C_rate, rate)
                & (detail.Direction == direction)
            ].iloc[0]
            ax.plot(exp["t_min"], exp["V"], color="black", lw=2.4, label="Experiment")
            ax.plot(sim["t_min"], sim["V"], color="#D55E00", lw=2.0, label="Pre-fit model")
            ax.axvline(exp["t_min"][-1], color="black", ls=":", alpha=0.55)
            ax.axvline(sim["t_min"][-1], color="#D55E00", ls=":", alpha=0.65)
            ax.set_title(
                f"{rate:g}C {direction}\n"
                f"RMSE(10–70%)={metric.Center10_70_RMSE_mV:.1f} mV, "
                f"ΔQ={metric.Capacity_error_mAh:+.0f} mAh"
            )
            ax.set_xlabel("Time [min]")
            ax.set_ylabel("Voltage [V]")
            ax.grid(alpha=0.22)
    axes[0, 0].legend()
    fig.suptitle("Pre-fitting time–voltage curves | current baseline")
    fig.savefig(RESULTS / "prefit_time_voltage_curves.png", dpi=220)
    plt.close(fig)

    # Capacity error by condition for k and Ds combinations.
    kinetic_names = [
        "Current baseline",
        "k: nominal kn + nominal kp",
        "k: nominal kn + current kp",
        "k: current kn + nominal kp",
    ]
    diffusion_names = [
        "Current baseline",
        "Ds: nominal Dsn + nominal Dsp",
        "Ds: nominal Dsn + current Dsp",
        "Ds: current Dsn + nominal Dsp",
    ]
    short = {
        "Current baseline": "Current",
        "k: nominal kn + nominal kp": "kn nom / kp nom",
        "k: nominal kn + current kp": "kn nom / kp current",
        "k: current kn + nominal kp": "kn current / kp nom",
        "Ds: nominal Dsn + nominal Dsp": "Dsn nom / Dsp nom",
        "Ds: nominal Dsn + current Dsp": "Dsn nom / Dsp current",
        "Ds: current Dsn + nominal Dsp": "Dsn current / Dsp nom",
    }
    marker = ["o", "s", "^", "D"]
    palette = ["#222222", "#0072B2", "#D55E00", "#009E73"]
    fig, axes = plt.subplots(2, 2, figsize=(14.0, 8.4), constrained_layout=True)
    for row, (names, title) in enumerate(
        [(kinetic_names, "Kinetic prefactor combinations"), (diffusion_names, "Solid diffusivity combinations")]
    ):
        for col, direction in enumerate(("Charge", "Discharge")):
            ax = axes[row, col]
            for idx, name in enumerate(names):
                part = detail[(detail.Scenario == name) & (detail.Direction == direction)].sort_values("C_rate")
                ax.plot(
                    part.C_rate,
                    part.Capacity_error_mAh,
                    marker=marker[idx],
                    color=palette[idx],
                    lw=1.8,
                    ms=6,
                    label=short[name],
                )
            ax.axhline(0, color="black", lw=0.8)
            ax.set_xticks(ga.RATES, [f"{r:g}C" for r in ga.RATES])
            ax.set_ylabel("Model − experiment [mAh]")
            ax.set_title(f"{title} | {direction}")
            ax.grid(alpha=0.22)
            ax.legend(fontsize=8)
    fig.suptitle("Pre-fitting cutoff-capacity error ablation")
    fig.savefig(RESULTS / "capacity_error_ablation.png", dpi=220)
    plt.close(fig)

    print("\nSUMMARY")
    print(summary.to_string(index=False))
    print("\nSaved", RESULTS)


if __name__ == "__main__":
    main()
