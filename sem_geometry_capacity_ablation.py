"""Ablate the 2026-09-22 SEM geometry before dry-loading is available.

The raw cases hold porosity and maximum concentration fixed, so changing
electrode thickness changes the implied active-material inventory.  The
Q-preserved cases rescale c_s,max inversely with thickness and therefore
isolate the transport/geometric effect from the unknown coating loading.
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
import current_thickness_estimate_application as thickness
import gitt_ocp_analysis as ga
import priority_initial_state_protocol_recheck as priority


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260922_sem_geometry_capacity_ablation"
RESULTS.mkdir(parents=True, exist_ok=True)

# Cross-section SEM totals with independently measured collectors removed.
LN_SEM_UM = (172.04 - 9.0) / 2.0
LP_SEM_UM = (138.93 - 13.0) / 2.0
RN_SEM_UM = 6.0
RP_SEM_UM = 5.0
CU_UM = 9.0
AL_UM = 13.0


def copy_model(model):
    changed = dict(model)
    changed["params"] = model["params"].copy()
    return changed


def set_raw_sem_thickness(model):
    """Change geometry while holding porosity and c_s,max fixed."""
    changed = copy_model(model)
    changed["params"].update(
        {
            "Negative electrode thickness [m]": LN_SEM_UM * 1e-6,
            "Positive electrode thickness [m]": LP_SEM_UM * 1e-6,
            "Negative current collector thickness [m]": CU_UM * 1e-6,
            "Positive current collector thickness [m]": AL_UM * 1e-6,
        },
        check_already_exists=False,
    )
    return changed


def set_q_preserved_sem_thickness(model):
    """Change thickness but preserve fitted Qn/Qp via inverse c_s,max scaling."""
    changed = thickness.set_thickness_preserve_capacity(
        model, LN_SEM_UM * 1e-6, LP_SEM_UM * 1e-6
    )
    changed["params"].update(
        {
            "Negative current collector thickness [m]": CU_UM * 1e-6,
            "Positive current collector thickness [m]": AL_UM * 1e-6,
        },
        check_already_exists=False,
    )
    return changed


def set_sem_radii(model, preserve_diffusion_time=False):
    changed = copy_model(model)
    old_rn = float(changed["params"]["Negative particle radius [m]"])
    old_rp = float(changed["params"]["Positive particle radius [m]"])
    new_rn = RN_SEM_UM * 1e-6
    new_rp = RP_SEM_UM * 1e-6
    updates = {
        "Negative particle radius [m]": new_rn,
        "Positive particle radius [m]": new_rp,
    }
    if preserve_diffusion_time:
        updates.update(
            {
                "Negative particle diffusivity [m2.s-1]": float(
                    changed["params"]["Negative particle diffusivity [m2.s-1]"]
                )
                * (new_rn / old_rn) ** 2,
                "Positive particle diffusivity [m2.s-1]": float(
                    changed["params"]["Positive particle diffusivity [m2.s-1]"]
                )
                * (new_rp / old_rp) ** 2,
            }
        )
    changed["params"].update(updates, check_already_exists=False)
    return changed


def inventory_ratio(model, reference, electrode):
    if electrode == "negative":
        prefix = "Negative"
    else:
        prefix = "Positive"
    length_key = f"{prefix} electrode thickness [m]"
    eps_key = f"{prefix} electrode active material volume fraction"
    cs_key = f"Maximum concentration in {electrode} electrode [mol.m-3]"
    numerator = (
        float(model["params"][length_key])
        * float(model["params"][eps_key])
        * float(model["params"][cs_key])
    )
    denominator = (
        float(reference["params"][length_key])
        * float(reference["params"][eps_key])
        * float(reference["params"][cs_key])
    )
    return numerator / denominator


def main():
    base, stage, anode, cathode, experiment, _, qcell = selected.selected_inputs()
    protocol, _, _ = priority.load_protocol(base)

    raw_thickness = set_raw_sem_thickness(base)
    q_preserved = set_q_preserved_sem_thickness(base)
    scenarios = {
        "Current baseline": base,
        "Nominal thickness + SEM Rn only": set_sem_radii(base),
        "Nominal thickness + SEM Rp only": set_sem_radii(base),
        "Nominal thickness + SEM radii": set_sem_radii(base),
        "Nominal thickness + SEM radii tau-preserved": set_sem_radii(
            base, preserve_diffusion_time=True
        ),
        "SEM thickness raw": raw_thickness,
        "SEM thickness + radii raw": set_sem_radii(raw_thickness),
        "SEM thickness + radii raw tau-preserved": set_sem_radii(
            raw_thickness, preserve_diffusion_time=True
        ),
        "SEM thickness Q-preserved": q_preserved,
        "SEM thickness + radii Q-preserved": set_sem_radii(q_preserved),
        "SEM thickness + radii Q-preserved tau-preserved": set_sem_radii(
            q_preserved, preserve_diffusion_time=True
        ),
    }
    scenarios["Nominal thickness + SEM Rn only"]["params"].update(
        {"Positive particle radius [m]": float(base["params"]["Positive particle radius [m]"])},
        check_already_exists=False,
    )
    scenarios["Nominal thickness + SEM Rp only"]["params"].update(
        {"Negative particle radius [m]": float(base["params"]["Negative particle radius [m]"])},
        check_already_exists=False,
    )

    runs = {}
    for name, model in scenarios.items():
        runs[name] = selected.run_all(
            name,
            model,
            stage,
            anode,
            cathode,
            False,
            True,
            protocol=protocol,
            history=True,
        )

    detail = selected.score(runs, experiment, qcell)
    detail.to_csv(RESULTS / "sem_geometry_branch_metrics.csv", index=False, encoding="utf-8-sig")
    valid = detail[detail.Status == "ok"].copy()
    summary = valid.groupby("Scenario", as_index=False, sort=False).agg(
        mean_center_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
        capacity_RMSE_pct=("Capacity_error_pct", lambda x: float(np.sqrt(np.mean(x**2)))),
        mean_capacity_error_mAh=("Capacity_error_mAh", "mean"),
        mean_abs_capacity_error_mAh=("Capacity_error_mAh", lambda x: float(np.mean(np.abs(x)))),
        max_abs_capacity_error_mAh=("Capacity_error_mAh", lambda x: float(np.max(np.abs(x)))),
    )
    directional = valid.groupby(["Scenario", "Direction"], as_index=False, sort=False).agg(
        mean_center_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
        capacity_RMSE_pct=("Capacity_error_pct", lambda x: float(np.sqrt(np.mean(x**2)))),
        mean_capacity_error_mAh=("Capacity_error_mAh", "mean"),
        mean_abs_capacity_error_mAh=("Capacity_error_mAh", lambda x: float(np.mean(np.abs(x)))),
    )

    inventory = []
    for name, model in scenarios.items():
        inventory.append(
            {
                "Scenario": name,
                "Ln_um": float(model["params"]["Negative electrode thickness [m]"]) * 1e6,
                "Lp_um": float(model["params"]["Positive electrode thickness [m]"]) * 1e6,
                "Rn_um": float(model["params"]["Negative particle radius [m]"]) * 1e6,
                "Rp_um": float(model["params"]["Positive particle radius [m]"]) * 1e6,
                "negative_inventory_ratio": inventory_ratio(model, base, "negative"),
                "positive_inventory_ratio": inventory_ratio(model, base, "positive"),
                "Dsn_m2_s": float(model["params"]["Negative particle diffusivity [m2.s-1]"]),
                "Dsp_m2_s": float(model["params"]["Positive particle diffusivity [m2.s-1]"]),
            }
        )
    inventory = pd.DataFrame(inventory)
    summary = summary.merge(inventory, on="Scenario", how="left")
    summary.to_csv(RESULTS / "sem_geometry_summary.csv", index=False, encoding="utf-8-sig")
    directional.to_csv(RESULTS / "sem_geometry_direction_summary.csv", index=False, encoding="utf-8-sig")

    plot_names = ["Current baseline", "SEM thickness raw", "SEM thickness Q-preserved"]
    colors = {
        "Current baseline": "#666666",
        "SEM thickness raw": "#D55E00",
        "SEM thickness Q-preserved": "#0072B2",
    }
    fig, axes = plt.subplots(2, 3, figsize=(16.5, 8.8), constrained_layout=True)
    for col, rate in enumerate(ga.RATES):
        for row, charge in enumerate((True, False)):
            direction = "Charge" if charge else "Discharge"
            ax = axes[row, col]
            exp = experiment[(rate, charge)]
            q_exp = priority.transferred_capacity(exp, qcell, charge)
            ax.plot(q_exp, exp["V"], color="black", lw=2.4, label="Experiment")
            for name in plot_names:
                sim = runs[name][(rate, charge)]
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
                    lw=1.7,
                    label=f"{name}: cap {metric.Capacity_error_mAh:+.0f} mAh",
                )
            ax.set_title(f"{rate:g}C {direction}")
            ax.set_xlabel("Transferred capacity [Ah]")
            ax.set_ylabel("Voltage [V]")
            ax.grid(alpha=0.22)
    axes[0, 0].legend(fontsize=7.0)
    fig.suptitle("SEM thickness ablation before coating-loading measurement")
    fig.savefig(RESULTS / "sem_thickness_dynamic_curves.png", dpi=220)
    plt.close(fig)

    x = np.arange(len(summary))
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)
    axes[0].bar(x, summary.mean_abs_capacity_error_mAh, color="#E69F00")
    axes[0].set_ylabel("Mean |capacity error| [mAh]")
    axes[1].bar(x, summary.capacity_RMSE_pct, color="#CC79A7")
    axes[1].set_ylabel("Capacity RMSE [%]")
    axes[2].bar(x, summary.mean_center_RMSE_mV, color="#56B4E9")
    axes[2].set_ylabel("Mean voltage RMSE, 10-70% [mV]")
    for ax in axes:
        ax.set_xticks(x, [s.replace("SEM thickness ", "SEM ") for s in summary.Scenario], rotation=24, ha="right")
        ax.grid(axis="y", alpha=0.22)
    fig.savefig(RESULTS / "sem_geometry_metric_comparison.png", dpi=220)
    plt.close(fig)

    report = {
        "sem_inputs": {
            "negative_double_sided_total_um": 172.04,
            "Cu_um": CU_UM,
            "Ln_um": LN_SEM_UM,
            "positive_double_sided_total_um": 138.93,
            "Al_um": AL_UM,
            "Lp_um": LP_SEM_UM,
            "Rn_um_provisional": RN_SEM_UM,
            "Rp_um_provisional": RP_SEM_UM,
        },
        "fixed": {
            "OCP": "C/50-selected Ai2020 nominal anode + old cathode",
            "initial_SOC": "per-test preceding-rest equivalent SOC",
            "negative_hysteresis": False,
            "positive_hysteresis": True,
            "porosity": "unchanged",
            "kinetics_diffusion_bruggeman": "unchanged",
        },
        "interpretation": {
            "raw": "c_s,max fixed; implied Qn/Qp change with thickness and can absorb loading error",
            "Q_preserved": "c_s,max inversely scaled; Qn/Qp fixed and only geometry/transport changes",
        },
        "summary": summary.to_dict(orient="records"),
        "direction_summary": directional.to_dict(orient="records"),
    }
    (RESULTS / "sem_geometry_capacity_ablation.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print("\nSUMMARY")
    print(summary.to_string(index=False))
    print("\nDIRECTION")
    print(directional.to_string(index=False))
    print("\nDETAIL")
    print(valid[["Scenario", "C_rate", "Direction", "Center10_70_RMSE_mV", "Capacity_error_mAh", "Capacity_error_pct"]].to_string(index=False))
    print("\nSaved to", RESULTS)


if __name__ == "__main__":
    main()
