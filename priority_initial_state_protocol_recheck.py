"""Priority recheck of initial state, hysteresis history, and R^2/Ds.

This is a diagnostic, not a new physical parameter recommendation.  It audits
the raw BoL protocol, infers each branch's initial SOC from its pre-step rest
voltage, and then separates initial-SOC, hysteresis-history, and solid-
diffusion-time effects before a small kn re-fit.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybamm

import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import nominal_radius_current_configuration as current
import stage2_seven_parameter_sensitivity as stage2


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260921_priority_initial_state_protocol_recheck"
RESULTS.mkdir(parents=True, exist_ok=True)

DYNAMIC_RAW = Path(
    r"C:\Users\user\OneDrive\01. 적응형 충전 프로토콜\4. Enertech 파우치셀 실험데이터"
    r"\260808 Enertech셀 열화데이터\260808 Eneterch셀 열화데이터 v1 6-6.xlsx"
)
LOWRATE_RAW = Path(
    r"C:\Users\user\OneDrive\01. 적응형 충전 프로토콜\4. Enertech 파우치셀 실험데이터"
    r"\260808 Eneterch셀 초기 저율 데이터\260808 Enertech셀 초기 저율 데이터 v2.xlsx"
)

BRANCH_STEPS = {
    (0.5, True): 7,
    (0.5, False): 10,
    (1.0, True): 12,
    (1.0, False): 15,
    (2.0, True): 17,
    (2.0, False): 20,
}

ATTACHED_RN_UM = 2.853
ATTACHED_RP_UM = 1.035
ATTACHED_DSN = 4.124e-13
ATTACHED_DSP = 1.194e-12


def monotone_voltage_inverse(voltage_grid, soc_grid, voltage):
    order = np.argsort(voltage_grid, kind="stable")
    v_unique, index = np.unique(np.asarray(voltage_grid)[order], return_index=True)
    z_unique = np.asarray(soc_grid)[order][index]
    return float(np.interp(voltage, v_unique, z_unique))


def load_protocol(model_inputs):
    if not DYNAMIC_RAW.exists() or not LOWRATE_RAW.exists():
        raise FileNotFoundError("Required raw BoL workbook is missing")

    step = pd.read_excel(DYNAMIC_RAW, sheet_name="step")
    record = pd.read_excel(DYNAMIC_RAW, sheet_name="record")
    for column in ("Oneset Date", "End Date"):
        step[column] = pd.to_datetime(step[column])
    record["Date"] = pd.to_datetime(record["Date"])

    rows = []
    rest_curves = {}
    for (rate, charge), step_index in BRANCH_STEPS.items():
        branch = step[(step["Cycle Index"] == 1) & (step["Step Index"] == step_index)].iloc[0]
        previous = step[
            (step["Cycle Index"] == 1)
            & (step["Step Number"] == int(branch["Step Number"]) - 1)
        ].iloc[0]
        segment = record[
            (record["Date"] >= previous["Oneset Date"])
            & (record["Date"] < previous["End Date"])
            & (record["Step Type"].astype(str) == "Rest")
        ].copy()
        segment["rest_time_min"] = (
            segment["Date"] - segment["Date"].iloc[0]
        ).dt.total_seconds() / 60.0
        tail = segment[
            segment["rest_time_min"] >= segment["rest_time_min"].max() - 3.0
        ]
        slope = float(
            np.polyfit(tail["rest_time_min"], tail["Voltage(V)"], 1)[0] * 1000.0
        )
        v_start = float(segment["Voltage(V)"].iloc[0])
        v_end = float(segment["Voltage(V)"].iloc[-1])
        initial_soc = monotone_voltage_inverse(
            model_inputs["v_qocv"], model_inputs["soc"], v_end
        )
        key = (rate, charge)
        rest_curves[key] = segment[["rest_time_min", "Voltage(V)"]].copy()
        rows.append(
            {
                "C_rate": rate,
                "Direction": "Charge" if charge else "Discharge",
                "Branch_step": int(step_index),
                "Previous_step": int(previous["Step Index"]),
                "Previous_step_type": str(previous["Step Type"]),
                "Rest_duration_min": float(
                    (previous["End Date"] - previous["Oneset Date"]).total_seconds() / 60.0
                ),
                "Rest_start_V": v_start,
                "Rest_end_V": v_end,
                "Rest_relaxation_mV": (v_end - v_start) * 1000.0,
                "Last_3min_slope_mV_per_min": slope,
                "qOCV_equivalent_initial_SOC_pct": initial_soc * 100.0,
                "Branch_onset_V": float(branch["Oneset Volt.(V)"]),
                "Branch_end_V": float(branch["End Voltage(V)"]),
                "Branch_capacity_Ah": float(branch["Capacity(Ah)"]),
            }
        )

    low_step = pd.read_excel(LOWRATE_RAW, sheet_name="step")
    qocv_rows = []
    # The integrated primary qOCV uses Step 10 charge and Step 9 discharge.
    for step_number, role in ((9, "qOCV discharge primary"), (10, "qOCV charge primary")):
        row = low_step[low_step["Step Number"] == step_number].iloc[0]
        previous = low_step[low_step["Step Number"] == step_number - 1].iloc[0]
        gap_s = (
            pd.to_datetime(row["Oneset Date"]) - pd.to_datetime(previous["End Date"])
        ).total_seconds()
        qocv_rows.append(
            {
                "Role": role,
                "Step_number": step_number,
                "Step_type": str(row["Step Type"]),
                "Duration_h": pd.to_timedelta(str(row["Step Time"])).total_seconds() / 3600.0,
                "Capacity_Ah": float(row["Capacity(Ah)"]),
                "Onset_V": float(row["Oneset Volt.(V)"]),
                "End_V": float(row["End Voltage(V)"]),
                "Previous_step_type": str(previous["Step Type"]),
                "Rest_before_step_min": gap_s / 60.0,
            }
        )
    return pd.DataFrame(rows), rest_curves, pd.DataFrame(qocv_rows)


def stage_at_soc(stage, z, charge):
    changed = dict(stage)
    x = stage["x0"] + z * (stage["x100"] - stage["x0"])
    y = stage["y0"] - z * (stage["y0"] - stage["y100"])
    if charge:
        changed["x0"], changed["y0"] = float(x), float(y)
    else:
        changed["x100"], changed["y100"] = float(x), float(y)
    return changed


def model_with_transport(model_inputs, rn_um, rp_um, dsn, dsp):
    changed = dict(model_inputs)
    changed["params"] = model_inputs["params"].copy()
    changed["params"].update(
        {
            "Negative particle radius [m]": rn_um * 1e-6,
            "Positive particle radius [m]": rp_um * 1e-6,
            "Negative particle diffusivity [m2.s-1]": dsn,
            "Positive particle diffusivity [m2.s-1]": dsp,
        },
        check_already_exists=False,
    )
    return changed


def model_with_kn(model_inputs, kn_value):
    changed = dict(model_inputs)
    changed["params"] = model_inputs["params"].copy()
    nominal = model_inputs["params"][stage2.PARAMETER_KEYS["kn"]]
    changed["params"].update(
        {
            stage2.PARAMETER_KEYS["kn"]: ga.scale_parameter_function(
                nominal, float(kn_value) / current.KN_PREF
            )
        },
        check_already_exists=False,
    )
    return changed


def run_scenario(
    name,
    model_inputs,
    stage,
    anode,
    cathode,
    protocol,
    use_rest_soc,
    negative_hysteresis,
    positive_hysteresis,
    history_consistent,
):
    runs = {}
    protocol_index = protocol.set_index(["C_rate", "Direction"])
    for rate in ga.RATES:
        for charge in (True, False):
            direction = "Charge" if charge else "Discharge"
            z = (
                float(protocol_index.loc[(rate, direction), "qOCV_equivalent_initial_SOC_pct"])
                / 100.0
                if use_rest_soc
                else (0.0 if charge else 1.0)
            )
            local_stage = stage_at_soc(stage, z, charge)
            neg_h = None
            pos_h = None
            if history_consistent:
                # The preceding CC/CV branch has the opposite direction.
                neg_h = 1.0 if charge else -1.0
                pos_h = -1.0 if charge else 1.0
            print(f"{name} | {rate:g}C {direction} | z0={z:.5f}", flush=True)
            runs[(rate, charge)] = ehq.run_dfn(
                rate,
                charge,
                local_stage,
                model_inputs,
                anode,
                cathode,
                negative_hysteresis,
                positive_hysteresis,
                negative_initial_h=neg_h,
                positive_initial_h=pos_h,
            )
    return runs


def transferred_capacity(data, q_meas, charge):
    soc = np.asarray(data["SOC"], dtype=float)
    return soc * q_meas if charge else (1.0 - soc) * q_meas


def curve_metrics(experiment, simulation, q_meas, charge):
    qe = transferred_capacity(experiment, q_meas, charge)
    qm = transferred_capacity(simulation, q_meas, charge)
    ve = np.asarray(experiment["V"], dtype=float)
    vm = np.asarray(simulation["V"], dtype=float)
    oe = np.argsort(qe, kind="stable")
    om = np.argsort(qm, kind="stable")
    qe_u, ie = np.unique(qe[oe], return_index=True)
    qm_u, im = np.unique(qm[om], return_index=True)
    q_end_exp = float(qe_u[-1])
    q_end_model = float(qm_u[-1])
    upper = 0.995 * min(q_end_exp, q_end_model)
    full_grid = np.linspace(0.0, upper, 400)
    center_lo = 0.10 * q_meas
    center_hi = min(0.70 * q_meas, upper)
    center_grid = np.linspace(center_lo, center_hi, 160)

    def errors(grid):
        exp_v = np.interp(grid, qe_u, ve[oe][ie])
        mod_v = np.interp(grid, qm_u, vm[om][im])
        return (mod_v - exp_v) * 1000.0

    full = errors(full_grid)
    center = errors(center_grid)
    return {
        "Full_MAE_mV": float(np.mean(np.abs(full))),
        "Full_RMSE_mV": float(np.sqrt(np.mean(full**2))),
        "Full_Bias_mV": float(np.mean(full)),
        "Center10_70_MAE_mV": float(np.mean(np.abs(center))),
        "Center10_70_RMSE_mV": float(np.sqrt(np.mean(center**2))),
        "Initial_voltage_error_mV": float((vm[0] - ve[0]) * 1000.0),
        "Q_end_exp_Ah": q_end_exp,
        "Q_end_model_Ah": q_end_model,
        "Capacity_error_pct": (q_end_model - q_end_exp) / q_end_exp * 100.0,
    }


def score_scenarios(scenarios, experiment, q_meas):
    rows = []
    for name, runs in scenarios.items():
        for rate in ga.RATES:
            for charge in (True, False):
                rows.append(
                    {
                        "Scenario": name,
                        "C_rate": rate,
                        "Direction": "Charge" if charge else "Discharge",
                        **curve_metrics(
                            experiment[(rate, charge)],
                            runs[(rate, charge)],
                            q_meas,
                            charge,
                        ),
                    }
                )
    detail = pd.DataFrame(rows)
    summary = detail.groupby(["Scenario", "Direction"], as_index=False).agg(
        Mean_full_MAE_mV=("Full_MAE_mV", "mean"),
        Mean_full_RMSE_mV=("Full_RMSE_mV", "mean"),
        Mean_center_MAE_mV=("Center10_70_MAE_mV", "mean"),
        Mean_center_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
        Mean_abs_initial_error_mV=("Initial_voltage_error_mV", lambda x: float(np.mean(np.abs(x)))),
        Capacity_RMSE_pct=("Capacity_error_pct", lambda x: float(np.sqrt(np.mean(x**2)))),
    )
    return detail, summary


def scan_kn(model_inputs, stage, anode, cathode, protocol, experiment):
    values = np.unique(
        np.r_[np.linspace(0.50e-6, 1.80e-6, 14), current.KN_PREF]
    )
    protocol_index = protocol.set_index(["C_rate", "Direction"])
    rows = []
    cache = {}
    q_meas = float(model_inputs["q_meas"])
    for kn in values:
        changed = model_with_kn(model_inputs, kn)
        runs = {}
        for rate in ga.RATES:
            z = float(
                protocol_index.loc[(rate, "Charge"), "qOCV_equivalent_initial_SOC_pct"]
            ) / 100.0
            local_stage = stage_at_soc(stage, z, True)
            print(f"kn scan {kn:.4e} | {rate:g}C Charge", flush=True)
            runs[(rate, True)] = ehq.run_dfn(
                rate,
                True,
                local_stage,
                changed,
                anode,
                cathode,
                True,
                True,
                negative_initial_h=1.0,
                positive_initial_h=-1.0,
            )
        cache[float(kn)] = runs
        metrics = [
            curve_metrics(experiment[(rate, True)], runs[(rate, True)], q_meas, True)
            for rate in ga.RATES
        ]
        rows.append(
            {
                "kn": kn,
                "Charge_center_voltage_MSE_mV2": float(
                    np.mean([m["Center10_70_RMSE_mV"] ** 2 for m in metrics])
                ),
                "Charge_center_RMSE_mV": float(
                    np.sqrt(np.mean([m["Center10_70_RMSE_mV"] ** 2 for m in metrics]))
                ),
                "Charge_full_RMSE_mV": float(np.mean([m["Full_RMSE_mV"] for m in metrics])),
                "Charge_capacity_RMSE_pct": float(
                    np.sqrt(np.mean([m["Capacity_error_pct"] ** 2 for m in metrics]))
                ),
            }
        )
    scan = pd.DataFrame(rows)
    voltage_only = scan.loc[scan["Charge_center_voltage_MSE_mV2"].idxmin()].copy()
    voltage_only["Selection"] = "Voltage-only"
    feasible = scan[scan["Charge_capacity_RMSE_pct"] <= 3.0]
    if feasible.empty:
        capacity_constrained = scan.loc[scan["Charge_capacity_RMSE_pct"].idxmin()].copy()
        capacity_constrained["Selection"] = "Minimum-capacity-error (3% infeasible)"
    else:
        capacity_constrained = feasible.loc[
            feasible["Charge_center_voltage_MSE_mV2"].idxmin()
        ].copy()
        capacity_constrained["Selection"] = "Capacity-constrained 3%"
    selections = pd.DataFrame([voltage_only, capacity_constrained])

    validation_sets = {}
    for _, selected in selections.iterrows():
        label = str(selected["Selection"])
        kn = float(selected["kn"])
        changed = model_with_kn(model_inputs, kn)
        validation_runs = dict(cache[kn])
        for rate in ga.RATES:
            z = float(
                protocol_index.loc[(rate, "Discharge"), "qOCV_equivalent_initial_SOC_pct"]
            ) / 100.0
            local_stage = stage_at_soc(stage, z, False)
            print(
                f"{label} kn {kn:.4e} | {rate:g}C Discharge validation",
                flush=True,
            )
            validation_runs[(rate, False)] = ehq.run_dfn(
                rate,
                False,
                local_stage,
                changed,
                anode,
                cathode,
                True,
                True,
                negative_initial_h=-1.0,
                positive_initial_h=1.0,
            )
        validation_sets[label] = validation_runs
    return scan, selections, validation_sets


def plot_rest(protocol, rest_curves):
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.2), constrained_layout=True)
    for col, rate in enumerate(ga.RATES):
        for row, charge in enumerate((True, False)):
            ax = axes[row, col]
            curve = rest_curves[(rate, charge)]
            p = protocol[
                np.isclose(protocol.C_rate, rate)
                & (protocol.Direction == ("Charge" if charge else "Discharge"))
            ].iloc[0]
            ax.plot(curve.rest_time_min, curve["Voltage(V)"], color="#0072B2", lw=2)
            ax.set_title(
                f"Before {rate:g}C {p.Direction}\n"
                f"dV={p.Rest_relaxation_mV:+.1f} mV, final slope={p.Last_3min_slope_mV_per_min:+.2f} mV/min"
            )
            ax.set_xlabel("Rest time (min)")
            ax.set_ylabel("Voltage (V)")
            ax.grid(alpha=0.25)
    fig.suptitle("Actual 10-minute pre-branch rest in Enertech BoL RPT")
    fig.savefig(RESULTS / "bol_prebranch_rest_diagnostics.png", dpi=220)
    plt.close(fig)


def plot_scenarios(experiment, scenarios, detail, q_meas):
    selected_names = [
        "A Endpoint + default hysteresis",
        "B Rest-SOC + default hysteresis",
        "C Rest-SOC + history hysteresis",
        "D Rest-SOC + mean OCP",
        "E Rest-SOC + attached R2/Ds diagnostic",
    ]
    colors = ["#777777", "#0072B2", "#D55E00", "#009E73", "#CC79A7"]
    fig, axes = plt.subplots(2, 3, figsize=(16.2, 8.8), constrained_layout=True)
    for col, rate in enumerate(ga.RATES):
        for row, charge in enumerate((True, False)):
            ax = axes[row, col]
            exp = experiment[(rate, charge)]
            ax.plot(
                transferred_capacity(exp, q_meas, charge),
                exp["V"],
                color="black",
                lw=2.5,
                label="Experiment",
            )
            for name, color in zip(selected_names, colors):
                sim = scenarios[name][(rate, charge)]
                metric = detail[
                    (detail.Scenario == name)
                    & np.isclose(detail.C_rate, rate)
                    & (detail.Direction == ("Charge" if charge else "Discharge"))
                ].iloc[0]
                ax.plot(
                    transferred_capacity(sim, q_meas, charge),
                    sim["V"],
                    color=color,
                    lw=1.35,
                    label=f"{name[0]}: {metric.Full_RMSE_mV:.1f} mV",
                )
            ax.set_title(f"{rate:g}C {'Charge' if charge else 'Discharge'}")
            ax.set_xlabel("Transferred capacity (Ah)")
            ax.set_ylabel("Voltage (V)")
            ax.grid(alpha=0.2)
    axes[0, 0].legend(fontsize=7.5)
    fig.suptitle("Priority ablation under a common full-overlap metric")
    fig.savefig(RESULTS / "priority_scenario_curves.png", dpi=220)
    plt.close(fig)


def plot_summary(summary, scan):
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.6), constrained_layout=True)
    pivot = summary.pivot(index="Scenario", columns="Direction", values="Mean_full_RMSE_mV")
    x = np.arange(len(pivot))
    width = 0.36
    axes[0].bar(x - width / 2, pivot["Charge"], width, label="Charge")
    axes[0].bar(x + width / 2, pivot["Discharge"], width, label="Discharge")
    axes[0].set_xticks(x, [name[0] for name in pivot.index])
    axes[0].set_ylabel("Mean full-overlap RMSE (mV)")
    axes[0].set_xlabel("Scenario A-E")
    axes[0].set_title("Initial-state / hysteresis / diffusion ablation")
    axes[0].grid(axis="y", alpha=0.25)
    axes[0].legend()

    axes[1].semilogx(scan.kn, scan.Charge_center_RMSE_mV, "o-", label="10-70% voltage RMSE")
    ax2 = axes[1].twinx()
    ax2.semilogx(
        scan.kn,
        scan.Charge_capacity_RMSE_pct,
        "s--",
        color="#D55E00",
        label="capacity RMSE",
    )
    axes[1].set_xlabel("kn")
    axes[1].set_ylabel("Charge voltage RMSE (mV)")
    ax2.set_ylabel("Charge capacity RMSE (%)", color="#D55E00")
    axes[1].set_title("kn re-fit after state correction")
    axes[1].grid(alpha=0.25)
    fig.savefig(RESULTS / "priority_ablation_and_kn_refit.png", dpi=220)
    plt.close(fig)


def main():
    warnings.filterwarnings(
        "ignore", message="The definition of the hysteresis decay rate parameter has changed"
    )
    _, model_inputs, stage, anode, cathode, experiment, ocp_row = stage2.build_current_inputs()
    q_meas = float(model_inputs["q_meas"])
    protocol, rest_curves, qocv_protocol = load_protocol(model_inputs)
    protocol.to_csv(RESULTS / "bol_dynamic_protocol_audit.csv", index=False, encoding="utf-8-sig")
    qocv_protocol.to_csv(RESULTS / "bol_qocv_protocol_audit.csv", index=False, encoding="utf-8-sig")
    plot_rest(protocol, rest_curves)

    attached_transport = model_with_transport(
        model_inputs,
        ATTACHED_RN_UM,
        ATTACHED_RP_UM,
        ATTACHED_DSN,
        ATTACHED_DSP,
    )
    scenarios = {
        "A Endpoint + default hysteresis": run_scenario(
            "A Endpoint + default hysteresis",
            model_inputs,
            stage,
            anode,
            cathode,
            protocol,
            False,
            True,
            True,
            False,
        ),
        "B Rest-SOC + default hysteresis": run_scenario(
            "B Rest-SOC + default hysteresis",
            model_inputs,
            stage,
            anode,
            cathode,
            protocol,
            True,
            True,
            True,
            False,
        ),
        "C Rest-SOC + history hysteresis": run_scenario(
            "C Rest-SOC + history hysteresis",
            model_inputs,
            stage,
            anode,
            cathode,
            protocol,
            True,
            True,
            True,
            True,
        ),
        "D Rest-SOC + mean OCP": run_scenario(
            "D Rest-SOC + mean OCP",
            model_inputs,
            stage,
            anode,
            cathode,
            protocol,
            True,
            False,
            False,
            False,
        ),
        "E Rest-SOC + attached R2/Ds diagnostic": run_scenario(
            "E Rest-SOC + attached R2/Ds diagnostic",
            attached_transport,
            stage,
            anode,
            cathode,
            protocol,
            True,
            True,
            True,
            True,
        ),
    }
    detail, summary = score_scenarios(scenarios, experiment, q_meas)
    detail.to_csv(RESULTS / "priority_scenario_condition_metrics.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(RESULTS / "priority_scenario_summary.csv", index=False, encoding="utf-8-sig")
    plot_scenarios(experiment, scenarios, detail, q_meas)

    scan, selections, selected_runs = scan_kn(
        model_inputs, stage, anode, cathode, protocol, experiment
    )
    scan.to_csv(RESULTS / "rest_history_kn_scan.csv", index=False, encoding="utf-8-sig")
    selections.to_csv(
        RESULTS / "rest_history_kn_selections.csv", index=False, encoding="utf-8-sig"
    )
    selected_detail, selected_summary = score_scenarios(
        {f"Selected kn: {name}": runs for name, runs in selected_runs.items()},
        experiment,
        q_meas,
    )
    selected_detail.to_csv(
        RESULTS / "selected_kn_condition_metrics.csv", index=False, encoding="utf-8-sig"
    )
    selected_summary.to_csv(
        RESULTS / "selected_kn_summary.csv", index=False, encoding="utf-8-sig"
    )
    plot_summary(summary, scan)

    # Low-rate charge/discharge separation is reconstructed from the source bundle.
    bundle = ga.load_legacy_bundle(ga.get_paths(ROOT))
    low = bundle["full_c20"]
    qc, vc = ga._extract_qv_branch(low, "source_version", "v2", "charge")
    qd, vd = ga._extract_qv_branch(low, "source_version", "v2", "discharge")
    z = np.linspace(0.0, 1.0, 1001)
    v_charge = ga._interp_curve(ga._clean_curve(qc / qc.max(), vc), z)
    v_discharge = ga._interp_curve(ga._clean_curve(1.0 - qd / qd.max(), vd), z)
    mask = (z >= 0.10) & (z <= 0.90)
    branch_gap = (v_charge - v_discharge) * 1000.0

    tau_current = {
        "negative_s": (current.RN_UM * 1e-6) ** 2 / 2.1e-14,
        "positive_s": (current.RP_UM * 1e-6) ** 2 / 4.4e-14,
    }
    tau_attached = {
        "negative_s": (ATTACHED_RN_UM * 1e-6) ** 2 / ATTACHED_DSN,
        "positive_s": (ATTACHED_RP_UM * 1e-6) ** 2 / ATTACHED_DSP,
    }
    report = {
        "qOCV_MAE_mV": float(ocp_row.qOCV_MAE_2_98_mV),
        "raw_protocol": {
            "dynamic_source": str(DYNAMIC_RAW),
            "qOCV_source": str(LOWRATE_RAW),
            "dynamic_prebranch_rest_min": 10.0,
            "qOCV_primary_steps": [9, 10],
            "qOCV_rest_before_primary_steps_min": qocv_protocol[
                "Rest_before_step_min"
            ].tolist(),
            "qOCV_charge_discharge_gap_10_90_MAE_mV": float(
                np.mean(np.abs(branch_gap[mask]))
            ),
            "qOCV_charge_discharge_gap_10_90_max_mV": float(
                np.max(np.abs(branch_gap[mask]))
            ),
        },
        "diffusion_time_R2_over_D": {
            "current": tau_current,
            "attached_notebook": tau_attached,
            "ratio_current_over_attached": {
                key: tau_current[key] / tau_attached[key] for key in tau_current
            },
        },
        "scenario_summary": summary.to_dict(orient="records"),
        "selected_kn_after_state_correction": selections.to_dict(orient="records"),
        "selected_kn_validation": selected_summary.to_dict(orient="records"),
        "interpretation_guardrail": (
            "Scenario E copies the attached notebook's R and Ds only as a causal "
            "R^2/D diagnostic; it is not a parameter recommendation for the Enertech cell."
        ),
    }
    (RESULTS / "priority_recheck_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print("\nPROTOCOL")
    print(protocol.to_string(index=False))
    print("\nSCENARIO SUMMARY")
    print(summary.to_string(index=False))
    print("\nSELECTED KN")
    print(selections.to_string(index=False))
    print("\nSELECTED KN VALIDATION")
    print(selected_summary.to_string(index=False))


if __name__ == "__main__":
    main()
