from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybamm


PROJECT = Path(__file__).resolve().parents[1]
REPO = Path(
    r"C:\Users\user\Documents\Codex\2026-09-18"
    r"\https-github-com-rinny-hi-ai2020\work\ai2020_dfn_fitting"
)
OUT = PROJECT / "outputs" / "260928_july_nominal_reversion_audit"
sys.path.insert(0, str(REPO))

import c50_ocp_source_combination_study as c50
import charge_physical_time_common as common
import comprehensive_prefit_cross_cohort as cross
import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import p0_five_parameter_sensitivity_area as legacy
import priority_initial_state_protocol_recheck as priority


FINAL = {
    "Dsn": 7.466200473885764e-14,
    "Dsp": 6.888175323807748e-14,
    "kn": 1.4502357402082455e-6,
    "brugg_n": 2.914,
}
AI_K = 1.0e-11 * 96485.33212
BEST_CELL = "6-8"


def copy_model(model: dict) -> dict:
    changed = dict(model)
    changed["params"] = model["params"].copy()
    return changed


def set_params(model: dict, updates: dict) -> dict:
    changed = copy_model(model)
    changed["params"].update(updates, check_already_exists=False)
    return changed


def make_scenarios(context: common.Context):
    raw = pybamm.ParameterValues("Ai2020")
    baseline = legacy.apply_parameters(
        context.model_geometry,
        {
            "Dsn": FINAL["Dsn"],
            "Dsp": FINAL["Dsp"],
            "kn": FINAL["kn"],
            "kp": context.kp,
            "brugg_n": FINAL["brugg_n"],
        },
    )
    nominal_cathode = c50.nominal_detail(raw["Positive electrode OCP [V]"])

    scenarios = {"Current July final": (baseline, context.cathode, "selected baseline")}
    scenarios["De raw PyBaMM"] = (
        set_params(
            baseline,
            {"Electrolyte diffusivity [m2.s-1]": raw["Electrolyte diffusivity [m2.s-1]"]},
        ),
        context.cathode,
        "removes the current x1e-4 unit-conversion factor",
    )
    scenarios["Dsn nominal"] = (
        set_params(
            baseline,
            {"Negative particle diffusivity [m2.s-1]": raw["Negative particle diffusivity [m2.s-1]"]},
        ),
        context.cathode,
        "Ai2020 graphite Dsn function",
    )
    scenarios["Dsp nominal"] = (
        set_params(
            baseline,
            {"Positive particle diffusivity [m2.s-1]": raw["Positive particle diffusivity [m2.s-1]"]},
        ),
        context.cathode,
        "Ai2020 LiCoO2 Dsp function",
    )
    scenarios["kn nominal"] = (
        set_params(
            baseline,
            {
                "Negative electrode exchange-current density [A.m-2]": raw[
                    "Negative electrode exchange-current density [A.m-2]"
                ]
            },
        ),
        context.cathode,
        "Ai2020 negative exchange-current function",
    )
    scenarios["kp nominal"] = (
        set_params(
            baseline,
            {
                "Positive electrode exchange-current density [A.m-2]": raw[
                    "Positive electrode exchange-current density [A.m-2]"
                ]
            },
        ),
        context.cathode,
        "Ai2020 positive exchange-current function",
    )
    scenarios["Rn nominal"] = (
        set_params(baseline, {"Negative particle radius [m]": raw["Negative particle radius [m]"]}),
        context.cathode,
        "Ai2020 negative particle radius",
    )
    scenarios["Rp nominal"] = (
        set_params(baseline, {"Positive particle radius [m]": raw["Positive particle radius [m]"]}),
        context.cathode,
        "Ai2020 positive particle radius",
    )
    scenarios["Electrode width nominal"] = (
        set_params(baseline, {"Electrode width [m]": raw["Electrode width [m]"]}),
        context.cathode,
        "Ai2020 width; changes common reaction area without capacity preservation",
    )
    csn = copy_model(baseline)
    csn["csn_max"] = float(raw["Maximum concentration in negative electrode [mol.m-3]"])
    csn["params"].update(
        {"Maximum concentration in negative electrode [mol.m-3]": csn["csn_max"]},
        check_already_exists=False,
    )
    scenarios["csn_max nominal"] = (csn, context.cathode, "Ai2020 negative maximum concentration")
    csp = copy_model(baseline)
    csp["csp_max"] = float(raw["Maximum concentration in positive electrode [mol.m-3]"])
    csp["params"].update(
        {"Maximum concentration in positive electrode [mol.m-3]": csp["csp_max"]},
        check_already_exists=False,
    )
    scenarios["csp_max nominal"] = (csp, context.cathode, "Ai2020 positive maximum concentration")

    ocp_model = copy_model(baseline)
    stage_vector = np.array(
        [context.stage["x0"], context.stage["x100"], context.stage["y100"], context.stage["y0"]]
    )
    ocp_model["v_qocv"] = c50.predict(
        stage_vector,
        ocp_model["soc"],
        context.anode,
        nominal_cathode,
        "equilibrium",
    )
    scenarios["Positive OCP nominal"] = (
        ocp_model,
        nominal_cathode,
        "Ai2020 cathode OCP on the currently selected stoichiometric window",
    )

    all_dynamic = set_params(
        baseline,
        {
            "Electrolyte diffusivity [m2.s-1]": raw["Electrolyte diffusivity [m2.s-1]"],
            "Negative particle diffusivity [m2.s-1]": raw["Negative particle diffusivity [m2.s-1]"],
            "Positive particle diffusivity [m2.s-1]": raw["Positive particle diffusivity [m2.s-1]"],
            "Negative electrode exchange-current density [A.m-2]": raw[
                "Negative electrode exchange-current density [A.m-2]"
            ],
            "Positive electrode exchange-current density [A.m-2]": raw[
                "Positive electrode exchange-current density [A.m-2]"
            ],
            "Negative particle radius [m]": raw["Negative particle radius [m]"],
            "Positive particle radius [m]": raw["Positive particle radius [m]"],
        },
    )
    scenarios["All dynamic nominal"] = (
        all_dynamic,
        context.cathode,
        "De, Dsn, Dsp, kn, kp, Rn, Rp reverted together; OCP/window/inventory retained",
    )
    return scenarios, raw


def simulate_model(context, model, cathode):
    runs = {}
    for rate in ga.RATES:
        for charge in (True, False):
            direction = "Charge" if charge else "Discharge"
            row = context.protocol.loc[("Old July BoL", rate, direction)]
            z0 = priority.monotone_voltage_inverse(
                model["v_qocv"], model["soc"], float(row.Rest_end_V)
            )
            local_stage = priority.stage_at_soc(context.stage, z0, charge)
            sim_rate = float(row.Measured_current_A) / cross.NOMINAL_CAPACITY_AH
            runs[(rate, charge)] = ehq.run_dfn(
                sim_rate,
                charge,
                local_stage,
                model,
                context.anode,
                cathode,
                False,
                False,
            )
    return runs


def branch_metric(experiment: dict, simulation: dict, qcell: float, charge: bool) -> dict:
    view = common._cc_view(experiment)
    grid = np.linspace(0.0, float(view["t_min"][-1]), common.N_POINTS)
    sim_t = np.asarray(simulation["t_min"], float)
    sim_v = np.asarray(simulation["V"], float)
    margin = float(sim_t[-1] - grid[-1])
    overlap_end = min(float(grid[-1]), float(sim_t[-1]))
    overlap_grid = np.linspace(0.0, overlap_end, common.N_POINTS)
    overlap_exp = np.interp(overlap_grid, view["t_min"], view["V"])
    overlap_mod = np.interp(overlap_grid, sim_t, sim_v)
    overlap_error = 1000.0 * (overlap_mod - overlap_exp)
    feasible = bool(margin >= -1e-9)
    time_rmse = np.nan
    time_bias = np.nan
    if feasible:
        exp_v = np.interp(grid, view["t_min"], view["V"])
        mod_v = np.interp(grid, sim_t, sim_v)
        error = 1000.0 * (mod_v - exp_v)
        time_rmse = float(np.sqrt(np.mean(error**2)))
        time_bias = float(np.mean(error))
    cap = priority.curve_metrics(experiment, simulation, qcell, charge)
    return {
        "Full_time_feasible": feasible,
        "Time_margin_min": margin,
        "Physical_time_RMSE_mV": time_rmse,
        "Physical_time_bias_mV": time_bias,
        "Overlap_time_RMSE_mV": float(np.sqrt(np.mean(overlap_error**2))),
        "Overlap_time_bias_mV": float(np.mean(overlap_error)),
        "Capacity_aligned_full_RMSE_mV": cap["Full_RMSE_mV"],
        "Capacity_error_pct": cap["Capacity_error_pct"],
        "Capacity_error_mAh": 1000.0 * (cap["Q_end_model_Ah"] - cap["Q_end_exp_Ah"]),
    }


def summarize(detail: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (scenario, direction), part in detail.groupby(["Scenario", "Direction"], sort=False):
        feasible = bool(part.Full_time_feasible.all())
        rows.append(
            {
                "Scenario": scenario,
                "Direction": direction,
                "All_full_time_feasible": feasible,
                "Minimum_time_margin_min": float(part.Time_margin_min.min()),
                "Physical_time_RMSE_mV": float(
                    np.sqrt(np.mean(part.Physical_time_RMSE_mV.to_numpy(float) ** 2))
                )
                if feasible
                else np.nan,
                "Overlap_time_RMSE_mV": float(
                    np.sqrt(np.mean(part.Overlap_time_RMSE_mV.to_numpy(float) ** 2))
                ),
                "Capacity_aligned_full_RMSE_mV": float(
                    np.sqrt(np.mean(part.Capacity_aligned_full_RMSE_mV.to_numpy(float) ** 2))
                ),
                "Capacity_RMSE_pct": float(
                    np.sqrt(np.mean(part.Capacity_error_pct.to_numpy(float) ** 2))
                ),
                "Capacity_MAE_mAh": float(np.mean(np.abs(part.Capacity_error_mAh))),
            }
        )
    return pd.DataFrame(rows)


def parameter_audit(context, raw) -> pd.DataFrame:
    p = context.model_geometry["params"]
    return pd.DataFrame(
        [
            {"Parameter": "De @1000 mol/m3, 298K", "Current": raw["Electrolyte diffusivity [m2.s-1]"](1000, 298.15) * 1e-4, "Nominal": raw["Electrolyte diffusivity [m2.s-1]"](1000, 298.15), "Status": "different; current includes unit conversion"},
            {"Parameter": "Dsn @x=0.5, 298K", "Current": FINAL["Dsn"], "Nominal": raw["Negative particle diffusivity [m2.s-1]"](0.5, 298.15), "Status": "different"},
            {"Parameter": "Dsp @y=0.7, 298K", "Current": FINAL["Dsp"], "Nominal": raw["Positive particle diffusivity [m2.s-1]"](0.7, 298.15), "Status": "different"},
            {"Parameter": "kn prefactor", "Current": FINAL["kn"], "Nominal": AI_K, "Status": "different"},
            {"Parameter": "kp prefactor", "Current": context.kp, "Nominal": AI_K, "Status": "different"},
            {"Parameter": "Rn [m]", "Current": p["Negative particle radius [m]"], "Nominal": raw["Negative particle radius [m]"], "Status": "different"},
            {"Parameter": "Rp [m]", "Current": p["Positive particle radius [m]"], "Nominal": raw["Positive particle radius [m]"], "Status": "different"},
            {"Parameter": "Ln [m]", "Current": p["Negative electrode thickness [m]"], "Nominal": raw["Negative electrode thickness [m]"], "Status": "already nominal"},
            {"Parameter": "Lp [m]", "Current": p["Positive electrode thickness [m]"], "Nominal": raw["Positive electrode thickness [m]"], "Status": "already nominal"},
            {"Parameter": "brugg_n", "Current": p["Negative electrode Bruggeman coefficient (electrolyte)"], "Nominal": raw["Negative electrode Bruggeman coefficient (electrolyte)"], "Status": "already nominal"},
            {"Parameter": "brugg_p", "Current": p["Positive electrode Bruggeman coefficient (electrolyte)"], "Nominal": raw["Positive electrode Bruggeman coefficient (electrolyte)"], "Status": "already nominal"},
            {"Parameter": "brugg_s", "Current": p["Separator Bruggeman coefficient (electrolyte)"], "Nominal": raw["Separator Bruggeman coefficient (electrolyte)"], "Status": "already nominal"},
            {"Parameter": "Electrode width [m]", "Current": p["Electrode width [m]"], "Nominal": raw["Electrode width [m]"], "Status": "different"},
            {"Parameter": "csn_max [mol/m3]", "Current": context.model_geometry["csn_max"], "Nominal": raw["Maximum concentration in negative electrode [mol.m-3]"], "Status": "different; current preserves selected Qn"},
            {"Parameter": "csp_max [mol/m3]", "Current": context.model_geometry["csp_max"], "Nominal": raw["Maximum concentration in positive electrode [mol.m-3]"], "Status": "different; current preserves selected Qp"},
        ]
    )


def plot_summary(summary: pd.DataFrame) -> None:
    scenarios = [name for name in summary.Scenario.unique() if name != "Current July final"]
    baseline = summary[summary.Scenario == "Current July final"].set_index("Direction")
    fig, axes = plt.subplots(2, 1, figsize=(14.5, 10.2), constrained_layout=True)
    colors = {"Charge": "#0072B2", "Discharge": "#D55E00"}
    width = 0.36
    x = np.arange(len(scenarios))
    for offset, direction in ((-width / 2, "Charge"), (width / 2, "Discharge")):
        part = summary[summary.Direction == direction].set_index("Scenario").loc[scenarios]
        delta = part.Overlap_time_RMSE_mV.to_numpy(float) - float(
            baseline.loc[direction, "Overlap_time_RMSE_mV"]
        )
        axes[0].bar(x + offset, delta, width, color=colors[direction], label=direction)
        cap_delta = part.Capacity_RMSE_pct.to_numpy(float) - float(
            baseline.loc[direction, "Capacity_RMSE_pct"]
        )
        axes[1].bar(x + offset, cap_delta, width, color=colors[direction], label=direction)
        for index, (_, row) in enumerate(part.iterrows()):
            if not bool(row.All_full_time_feasible):
                axes[0].text(x[index] + offset, delta[index], "x", ha="center", va="bottom", fontsize=10, fontweight="bold")
    short = [
        name.replace(" nominal", " nom").replace("Electrode width", "Width").replace("Positive OCP", "Cathode OCP").replace("All dynamic", "All dyn")
        for name in scenarios
    ]
    for ax in axes:
        ax.axhline(0.0, color="black", lw=1)
        ax.set_xticks(x, short, rotation=28, ha="right")
        ax.grid(axis="y", alpha=0.24)
    axes[0].set(ylabel="Change in overlap-time RMSE [mV]", title="Voltage error change vs current July final (x = full-time infeasible)")
    axes[1].set(ylabel="Change in capacity RMSE [percentage point]", title="Post-hoc capacity error change vs current July final")
    axes[0].legend()
    fig.suptitle("One-at-a-time rollback to Ai2020 nominal values", fontsize=15)
    fig.savefig(OUT / "nominal_reversion_error_attribution.png", dpi=220)
    plt.close(fig)


def plot_de(context, run_cache):
    fig, axes = plt.subplots(2, 3, figsize=(16.2, 9.5), constrained_layout=True)
    for row, charge in enumerate((True, False)):
        direction = "Charge" if charge else "Discharge"
        for col, rate in enumerate(ga.RATES):
            ax = axes[row, col]
            exp = common._cc_view(context.curves[BEST_CELL][(rate, charge)])
            base = run_cache[("Current July final", rate, charge)]
            raw = run_cache[("De raw PyBaMM", rate, charge)]
            ax.plot(exp["t_min"], exp["V"], color="black", lw=2.3, label=f"Experiment ({BEST_CELL})")
            ax.plot(base["t_min"], base["V"], color="#0072B2", lw=2.2, label="Current De x1e-4")
            ax.plot(raw["t_min"], raw["V"], color="#D55E00", lw=2.1, ls="--", label="Raw PyBaMM De")
            ax.set(title=f"{rate:g}C {direction}", xlabel="Time from 95% current [min]", ylabel="Voltage [V]")
            ax.grid(alpha=0.23)
            if row == 0 and col == 0:
                ax.legend(fontsize=8)
    fig.suptitle("Electrolyte diffusivity rollback: current unit-corrected De vs raw PyBaMM function", fontsize=14)
    fig.savefig(OUT / "de_raw_pybamm_curve_comparison.png", dpi=220)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    context = common.build_context()
    scenarios, raw = make_scenarios(context)
    audit = parameter_audit(context, raw)
    audit.to_csv(OUT / "current_vs_nominal_parameters.csv", index=False, encoding="utf-8-sig")
    rows = []
    failures = []
    run_cache = {}
    for scenario, (model, cathode, note) in scenarios.items():
        print(scenario, flush=True)
        try:
            runs = simulate_model(context, model, cathode)
        except Exception as exc:
            failures.append({"Scenario": scenario, "Error": f"{type(exc).__name__}: {exc}"})
            continue
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                run_cache[(scenario, rate, charge)] = runs[(rate, charge)]
                for cell, branches in context.curves.items():
                    rows.append(
                        {
                            "Scenario": scenario,
                            "Note": note,
                            "Cell": cell,
                            "C_rate": rate,
                            "Direction": direction,
                            **branch_metric(
                                branches[(rate, charge)],
                                runs[(rate, charge)],
                                context.qcell,
                                charge,
                            ),
                        }
                    )
    detail = pd.DataFrame(rows)
    summary = summarize(detail)
    detail.to_csv(OUT / "nominal_reversion_detail.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(OUT / "nominal_reversion_summary.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(failures).to_csv(OUT / "solver_failures.csv", index=False, encoding="utf-8-sig")
    plot_summary(summary)
    if ("Current July final", 0.5, True) in run_cache and ("De raw PyBaMM", 0.5, True) in run_cache:
        plot_de(context, run_cache)
    manifest = {
        "baseline": FINAL,
        "kp": context.kp,
        "nominal_kinetic_prefactor": AI_K,
        "De_interpretation": "current model uses raw Ai2020 function x1e-4 as cm2/s to m2/s unit conversion; raw-function rollback is diagnostic, not recommended physics",
        "already_nominal_and_not_resimulated": ["Ln", "Lp", "porosity_n", "porosity_p", "active_fraction_n", "active_fraction_p", "brugg_n", "brugg_p", "brugg_s"],
        "capacity_in_objective": False,
        "scenarios": list(scenarios),
        "failures": failures,
    }
    (OUT / "analysis_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(summary.to_string(index=False), flush=True)
    if failures:
        print(pd.DataFrame(failures).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
