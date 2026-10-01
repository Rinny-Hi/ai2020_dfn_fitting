"""Compare 3 electrode-thickness cases x 2 particle-radius cases."""

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
RESULTS = ROOT / "results" / "260927_thickness_radius_six_case_comparison"
RESULTS.mkdir(parents=True, exist_ok=True)

RADIUS_CASES_UM = {
    "Nominal R": {"Rn": 5.0, "Rp": 3.0},
    "SEM ImageJ R": {"Rn": 3.7542, "Rp": 4.2387},
}


def with_radius(model0: dict, rn_um: float, rp_um: float) -> dict:
    changed = thickness3.copy_model(model0)
    changed["params"].update(
        {
            "Negative particle radius [m]": rn_um * 1e-6,
            "Positive particle radius [m]": rp_um * 1e-6,
        },
        check_already_exists=False,
    )
    return changed


def main() -> None:
    model0, stage, anode, cathode, experiment, _, qcell = selected.selected_inputs()
    model0 = kpstudy.with_kp(model0, kpstudy.KP_EIS)
    protocol, _, _ = priority.load_protocol(model0)
    pidx = protocol.set_index(["C_rate", "Direction"])

    models = {}
    for thickness_case, lv in thickness3.CASES_UM.items():
        geometry = thickness3.same_loading_geometry(model0, lv["Ln"], lv["Lp"])
        for radius_case, rv in RADIUS_CASES_UM.items():
            models[(thickness_case, radius_case)] = with_radius(
                geometry, rv["Rn"], rv["Rp"]
            )

    rows = []
    simulations = {}
    for (thickness_case, radius_case), model in models.items():
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
                sim_rate = measured_current / thickness3.NOMINAL_CAPACITY_AH
                print(
                    f"{thickness_case} | {radius_case} | {rate:g}C {direction}",
                    flush=True,
                )
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
                simulations[(thickness_case, radius_case, rate, charge)] = sim
                metric = priority.curve_metrics(
                    experiment[(rate, charge)], sim, qcell, charge
                )
                rows.append(
                    {
                        "Thickness_case": thickness_case,
                        "Radius_case": radius_case,
                        "Ln_um": thickness3.CASES_UM[thickness_case]["Ln"],
                        "Lp_um": thickness3.CASES_UM[thickness_case]["Lp"],
                        "Rn_um": RADIUS_CASES_UM[radius_case]["Rn"],
                        "Rp_um": RADIUS_CASES_UM[radius_case]["Rp"],
                        "C_rate": rate,
                        "Direction": direction,
                        "initial_SOC": initial_soc,
                        **metric,
                        "Capacity_error_mAh": 1000.0
                        * (metric["Q_end_model_Ah"] - metric["Q_end_exp_Ah"]),
                    }
                )

    detail = pd.DataFrame(rows)
    detail.to_csv(RESULTS / "six_case_detail.csv", index=False, encoding="utf-8-sig")
    summary = detail.groupby(
        ["Thickness_case", "Radius_case", "Direction"], sort=False, as_index=False
    ).agg(
        Mean_full_RMSE_mV=("Full_RMSE_mV", "mean"),
        Mean_center_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
        Capacity_RMSE_pct=(
            "Capacity_error_pct",
            lambda x: float(np.sqrt(np.mean(x**2))),
        ),
        Mean_capacity_error_mAh=("Capacity_error_mAh", "mean"),
        Mean_abs_capacity_error_mAh=(
            "Capacity_error_mAh", lambda x: float(np.mean(np.abs(x)))
        ),
        Max_abs_capacity_error_mAh=(
            "Capacity_error_mAh", lambda x: float(np.max(np.abs(x)))
        ),
    )
    summary.to_csv(RESULTS / "six_case_summary.csv", index=False, encoding="utf-8-sig")

    overall = detail.groupby(
        ["Thickness_case", "Radius_case"], sort=False, as_index=False
    ).agg(
        Mean_full_RMSE_mV=("Full_RMSE_mV", "mean"),
        Mean_center_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
        Capacity_RMSE_pct=(
            "Capacity_error_pct",
            lambda x: float(np.sqrt(np.mean(x**2))),
        ),
        Mean_abs_capacity_error_mAh=(
            "Capacity_error_mAh", lambda x: float(np.mean(np.abs(x)))
        ),
        Max_abs_capacity_error_mAh=(
            "Capacity_error_mAh", lambda x: float(np.max(np.abs(x)))
        ),
    )
    overall.to_csv(RESULTS / "six_case_overall.csv", index=False, encoding="utf-8-sig")

    colors = {
        "Ai2020 paper": "#777777",
        "Thickness gauge": "#0072B2",
        "SEM cross-section": "#D55E00",
    }
    for radius_case in RADIUS_CASES_UM:
        fig, axes = plt.subplots(2, 3, figsize=(16.8, 9.0), constrained_layout=True)
        for col, rate in enumerate(ga.RATES):
            for row_idx, charge in enumerate((True, False)):
                ax = axes[row_idx, col]
                direction = "Charge" if charge else "Discharge"
                obs = experiment[(rate, charge)]
                q_exp = priority.transferred_capacity(obs, qcell, charge)
                ax.plot(q_exp, obs["V"], color="black", lw=2.6, label="Experiment")
                for thickness_case in thickness3.CASES_UM:
                    sim = simulations[(thickness_case, radius_case, rate, charge)]
                    result = detail[
                        (detail.Thickness_case == thickness_case)
                        & (detail.Radius_case == radius_case)
                        & np.isclose(detail.C_rate, rate)
                        & (detail.Direction == direction)
                    ].iloc[0]
                    q_sim = priority.transferred_capacity(sim, qcell, charge)
                    ax.plot(
                        q_sim,
                        sim["V"],
                        color=colors[thickness_case],
                        lw=1.8,
                        label=(
                            f"{thickness_case}: {result.Full_RMSE_mV:.1f} mV, "
                            f"{result.Capacity_error_mAh:+.1f} mAh"
                        ),
                    )
                ax.set(
                    title=f"{rate:g}C {direction}",
                    xlabel="Transferred capacity [Ah]",
                    ylabel="Voltage [V]",
                )
                ax.grid(alpha=0.22)
        axes[0, 0].legend(fontsize=6.8)
        fig.suptitle(f"Thickness comparison | {radius_case}")
        slug = "nominal" if radius_case == "Nominal R" else "sem_imagej"
        fig.savefig(RESULTS / f"capacity_voltage_{slug}_radius.png", dpi=220)
        plt.close(fig)

    # Heatmaps make the 3 x 2 factorial comparison explicit.
    metric_specs = [
        ("Mean_full_RMSE_mV", "Mean full RMSE [mV]"),
        ("Mean_center_RMSE_mV", "Mean 10-70% RMSE [mV]"),
        ("Capacity_RMSE_pct", "Capacity RMSE [%]"),
        ("Mean_abs_capacity_error_mAh", "Mean |capacity error| [mAh]"),
    ]
    fig, axes = plt.subplots(2, 4, figsize=(18.0, 8.3), constrained_layout=True)
    for row_idx, direction in enumerate(("Charge", "Discharge")):
        local = summary[summary.Direction == direction]
        for col_idx, (metric, title) in enumerate(metric_specs):
            ax = axes[row_idx, col_idx]
            matrix = local.pivot(
                index="Thickness_case", columns="Radius_case", values=metric
            ).loc[list(thickness3.CASES_UM), list(RADIUS_CASES_UM)]
            image = ax.imshow(matrix.values, cmap="viridis_r", aspect="auto")
            for i in range(matrix.shape[0]):
                for j in range(matrix.shape[1]):
                    ax.text(j, i, f"{matrix.iloc[i, j]:.2f}", ha="center", va="center", color="white" if matrix.iloc[i, j] > np.nanmean(matrix.values) else "black", fontsize=9)
            ax.set_xticks(range(matrix.shape[1]), matrix.columns, rotation=15, ha="right")
            ax.set_yticks(range(matrix.shape[0]), matrix.index)
            ax.set_title(f"{direction} | {title}")
            fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle("3 thickness x 2 radius cases")
    fig.savefig(RESULTS / "six_case_metric_heatmaps.png", dpi=220)
    plt.close(fig)

    nominal = detail[detail.Radius_case == "Nominal R"].set_index(
        ["Thickness_case", "C_rate", "Direction"]
    )
    sem = detail[detail.Radius_case == "SEM ImageJ R"].set_index(
        ["Thickness_case", "C_rate", "Direction"]
    )
    delta = (
        sem[["Full_RMSE_mV", "Center10_70_RMSE_mV", "Capacity_error_mAh"]]
        - nominal[["Full_RMSE_mV", "Center10_70_RMSE_mV", "Capacity_error_mAh"]]
    ).reset_index()
    delta = delta.rename(
        columns={
            "Full_RMSE_mV": "Delta_full_RMSE_mV_SEM_minus_nominal",
            "Center10_70_RMSE_mV": "Delta_center_RMSE_mV_SEM_minus_nominal",
            "Capacity_error_mAh": "Delta_capacity_error_mAh_SEM_minus_nominal",
        }
    )
    delta.to_csv(RESULTS / "sem_radius_minus_nominal_delta.csv", index=False, encoding="utf-8-sig")

    report = {
        "thickness_cases_um": thickness3.CASES_UM,
        "radius_cases_um": RADIUS_CASES_UM,
        "radius_source": (
            "SEM ImageJ arithmetic mean of per-particle 2D equivalent-circle radii; "
            "20 cathode and 44 anode ROIs"
        ),
        "fixed": {
            "Dsn_m2_s": 2.1e-14,
            "Dsp_m2_s": 4.4e-14,
            "kp": kpstudy.KP_EIS,
            "kn": "legacy 7.40e-7 comparison placeholder; fit later",
            "OCP_window": "current C/50-selected window",
            "initial_SOC": "preceding-rest voltage",
            "negative_hysteresis": False,
            "positive_hysteresis": True,
            "thickness_inventory": "L*epsilon_s fixed; inactive fraction fixed; porosity closes balance",
        },
        "summary": summary.to_dict(orient="records"),
        "overall": overall.to_dict(orient="records"),
        "sem_minus_nominal_delta": delta.to_dict(orient="records"),
    }
    (RESULTS / "six_case_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("\nSUMMARY")
    print(summary.to_string(index=False))
    print("\nOVERALL")
    print(overall.to_string(index=False))
    print("\nSEM RADIUS MINUS NOMINAL")
    print(delta.to_string(index=False))
    print("\nSaved to", RESULTS)


if __name__ == "__main__":
    main()
