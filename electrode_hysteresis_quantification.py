"""Quantify anode/cathode OCP hysteresis and test both in the DFN.

The 260918 new GITT data is deliberately excluded.  This analysis uses the
legacy anode low-rate v1/v2 data, legacy cathode GITT 7-7/7-8 rest-end data,
and the endpoint-constrained stoichiometry window.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybamm
from scipy.integrate import cumulative_trapezoid

import gitt_ocp_analysis as ga
import old_notebook_method_comparison as omc


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260920_electrode_hysteresis_quantification"
RESULTS.mkdir(parents=True, exist_ok=True)


def build_anode_detail(
    bundle: dict[str, Any], model: dict[str, Any]
) -> dict[str, Any]:
    diameter = float(bundle["prior"]("disk_diameter_m"))
    area = np.pi * (diameter / 2.0) ** 2
    params = model["params"]
    qn_disk = (
        ga.FARADAY_CONSTANT
        * model["csn_max"]
        * float(params["Negative electrode active material volume fraction"])
        * float(params["Negative electrode thickness [m]"])
        * area
        / 3600.0
    )
    frame = bundle["anode_lowrate"]
    lith = {}
    delith = {}
    for source in ("v1", "v2"):
        q_chg, v_chg = ga._extract_qv_branch(frame, "source_version", source, "charge")
        q_dch, v_dch = ga._extract_qv_branch(frame, "source_version", source, "discharge")
        lith[source] = ga._clean_curve(q_dch / qn_disk, v_dch)
        span = q_chg.max() / qn_disk
        delith[source] = ga._clean_curve(span - q_chg / qn_disk, v_chg)
    return combine_directional_curves(lith, delith)


def build_cathode_detail(
    bundle: dict[str, Any], model: dict[str, Any]
) -> dict[str, Any]:
    diameter = float(bundle["prior"]("disk_diameter_m"))
    area = np.pi * (diameter / 2.0) ** 2
    params = model["params"]
    qp_disk = (
        ga.FARADAY_CONSTANT
        * model["csp_max"]
        * float(params["Positive electrode active material volume fraction"])
        * float(params["Positive electrode thickness [m]"])
        * area
        / 3600.0
    )
    frame = bundle["cathode_gitt"]
    lith = {}
    delith = {}
    for source in ("7-7", "7-8"):
        branch = ga._extract_gitt_branch(frame, "source_channel", source, "discharge_pulse")
        span = branch[:, 0].max() / qp_disk
        lith[source] = ga._clean_curve((1.0 - span) + branch[:, 0] / qp_disk, branch[:, 1])
        branch = ga._extract_gitt_branch(frame, "source_channel", source, "charge_pulse")
        delith[source] = ga._clean_curve(1.0 - branch[:, 0] / qp_disk, branch[:, 1])
    return combine_directional_curves(lith, delith)


def combine_directional_curves(
    lith_replicates: dict[str, tuple[np.ndarray, np.ndarray]],
    delith_replicates: dict[str, tuple[np.ndarray, np.ndarray]],
) -> dict[str, Any]:
    x_lith, u_lith_0 = ga._mean_overlap(list(lith_replicates.values()), n=3000)
    x_delith, u_delith_0 = ga._mean_overlap(list(delith_replicates.values()), n=3000)
    lo = max(x_lith.min(), x_delith.min())
    hi = min(x_lith.max(), x_delith.max())
    grid = np.linspace(lo, hi, 3000)
    u_lith = np.interp(grid, x_lith, u_lith_0)
    u_delith = np.interp(grid, x_delith, u_delith_0)
    return {
        "grid": grid,
        "lithiation": u_lith,
        "delithiation": u_delith,
        "equilibrium": 0.5 * (u_lith + u_delith),
        "lithiation_replicates": lith_replicates,
        "delithiation_replicates": delith_replicates,
    }


def difference_metrics(values_mV: np.ndarray) -> dict[str, float]:
    absolute = np.abs(values_mV)
    return {
        "signed_mean_mV": float(np.mean(values_mV)),
        "MAE_mV": float(np.mean(absolute)),
        "RMSE_mV": float(np.sqrt(np.mean(values_mV**2))),
        "P95_abs_mV": float(np.quantile(absolute, 0.95)),
        "maximum_abs_mV": float(np.max(absolute)),
    }


def branch_metrics(
    electrode: str,
    detail: dict[str, Any],
    operating_lo: float,
    operating_hi: float,
) -> list[dict[str, Any]]:
    grid = detail["grid"]
    gap_mV = (detail["lithiation"] - detail["delithiation"]) * 1000.0
    rows = []
    for region, mask in (
        ("full_common_coverage", np.ones_like(grid, dtype=bool)),
        ("endpoint_operating_window", (grid >= operating_lo) & (grid <= operating_hi)),
        (
            "operating_window_central_10_90",
            (grid >= operating_lo + 0.10 * (operating_hi - operating_lo))
            & (grid <= operating_lo + 0.90 * (operating_hi - operating_lo)),
        ),
    ):
        rows.append(
            {
                "Electrode": electrode,
                "Comparison": "lithiation_minus_delithiation",
                "Region": region,
                "stoichiometry_min": float(grid[mask].min()),
                "stoichiometry_max": float(grid[mask].max()),
                **difference_metrics(gap_mV[mask]),
            }
        )

    for direction in ("lithiation", "delithiation"):
        replicates = detail[f"{direction}_replicates"]
        names = list(replicates)
        lo = max(replicates[name][0].min() for name in names)
        hi = min(replicates[name][0].max() for name in names)
        compare_grid = np.linspace(lo, hi, 3000)
        first = np.interp(compare_grid, *replicates[names[0]])
        second = np.interp(compare_grid, *replicates[names[1]])
        rows.append(
            {
                "Electrode": electrode,
                "Comparison": f"replicate_difference_{direction}_{names[0]}_minus_{names[1]}",
                "Region": "replicate_common_coverage",
                "stoichiometry_min": float(lo),
                "stoichiometry_max": float(hi),
                **difference_metrics((first - second) * 1000.0),
            }
        )
    return rows


def make_parameter_values(
    model_inputs: dict[str, Any],
    stage: dict[str, Any],
    anode: dict[str, Any],
    cathode: dict[str, Any],
    charge: bool,
    negative_hysteresis: bool,
    positive_hysteresis: bool,
    negative_initial_h: float | None = None,
    positive_initial_h: float | None = None,
    negative_decay_rate: float = 3.2,
    positive_decay_rate: float = 3.2,
):
    values = model_inputs["params"].copy()
    x_init = stage["x0"] if charge else stage["x100"]
    y_init = stage["y0"] if charge else stage["y100"]
    update = {
        "Initial concentration in negative electrode [mol.m-3]": x_init
        * model_inputs["csn_max"],
        "Initial concentration in positive electrode [mol.m-3]": y_init
        * model_inputs["csp_max"],
        "Negative electrode OCP [V]": ga.pybamm_interp(
            anode["grid"], anode["equilibrium"], "Experimental negative equilibrium OCP"
        ),
        "Positive electrode OCP [V]": ga.pybamm_interp(
            cathode["grid"], cathode["equilibrium"], "Experimental positive equilibrium OCP"
        ),
        "Initial temperature [K]": model_inputs["temperature"],
        "Ambient temperature [K]": model_inputs["temperature"],
    }
    if negative_hysteresis:
        update.update(
            {
                "Negative electrode lithiation OCP [V]": ga.pybamm_interp(
                    anode["grid"], anode["lithiation"], "Negative lithiation OCP"
                ),
                "Negative electrode delithiation OCP [V]": ga.pybamm_interp(
                    anode["grid"], anode["delithiation"], "Negative delithiation OCP"
                ),
                "Negative particle lithiation hysteresis decay rate": lambda sto, temp: negative_decay_rate
                + 0 * sto,
                "Negative particle delithiation hysteresis decay rate": lambda sto, temp: negative_decay_rate
                + 0 * sto,
                "Initial hysteresis state in negative electrode": (
                    (-1.0 if charge else 1.0)
                    if negative_initial_h is None
                    else float(negative_initial_h)
                ),
            }
        )
    if positive_hysteresis:
        update.update(
            {
                "Positive electrode lithiation OCP [V]": ga.pybamm_interp(
                    cathode["grid"], cathode["lithiation"], "Positive lithiation OCP"
                ),
                "Positive electrode delithiation OCP [V]": ga.pybamm_interp(
                    cathode["grid"], cathode["delithiation"], "Positive delithiation OCP"
                ),
                "Positive particle lithiation hysteresis decay rate": lambda sto, temp: positive_decay_rate
                + 0 * sto,
                "Positive particle delithiation hysteresis decay rate": lambda sto, temp: positive_decay_rate
                + 0 * sto,
                "Initial hysteresis state in positive electrode": (
                    (1.0 if charge else -1.0)
                    if positive_initial_h is None
                    else float(positive_initial_h)
                ),
            }
        )
    values.update(update, check_already_exists=False)
    return values


def run_dfn(
    rate: float,
    charge: bool,
    stage: dict[str, Any],
    model_inputs: dict[str, Any],
    anode: dict[str, Any],
    cathode: dict[str, Any],
    negative_hysteresis: bool,
    positive_hysteresis: bool,
    negative_initial_h: float | None = None,
    positive_initial_h: float | None = None,
    negative_decay_rate: float = 3.2,
    positive_decay_rate: float = 3.2,
) -> dict[str, np.ndarray]:
    negative_option = "one-state hysteresis" if negative_hysteresis else "single"
    positive_option = "one-state hysteresis" if positive_hysteresis else "single"
    model = pybamm.lithium_ion.DFN(
        {
            "open-circuit potential": (negative_option, positive_option),
            "thermal": "isothermal",
        }
    )
    step = f"{'Charge' if charge else 'Discharge'} at {rate}C until {4.2 if charge else 3.0} V"
    simulation = pybamm.Simulation(
        model,
        parameter_values=make_parameter_values(
            model_inputs,
            stage,
            anode,
            cathode,
            charge,
            negative_hysteresis,
            positive_hysteresis,
            negative_initial_h,
            positive_initial_h,
            negative_decay_rate,
            positive_decay_rate,
        ),
        experiment=pybamm.Experiment([step], period="10 seconds"),
        var_pts={"x_n": 15, "x_s": 15, "x_p": 15, "r_n": 15, "r_p": 15},
        solver=pybamm.IDAKLUSolver(rtol=1e-6, atol=1e-8),
    )
    solution = simulation.solve()
    time_s = np.asarray(solution["Time [s]"].entries, dtype=float)
    voltage = np.asarray(solution["Terminal voltage [V]"].entries, dtype=float)
    current = np.asarray(solution["Current [A]"].entries, dtype=float)
    capacity = cumulative_trapezoid(np.abs(current), time_s, initial=0.0) / 3600.0
    soc = capacity / model_inputs["q_meas"] if charge else 1.0 - capacity / model_inputs["q_meas"]
    negative_surface = np.asarray(
        solution["Negative particle surface concentration [mol.m-3]"].entries,
        dtype=float,
    ) / model_inputs["csn_max"]
    positive_surface = np.asarray(
        solution["Positive particle surface concentration [mol.m-3]"].entries,
        dtype=float,
    ) / model_inputs["csp_max"]
    return {
        "t_min": time_s / 60.0,
        "V": voltage,
        "SOC": np.clip(soc, 0.0, 1.0),
        "negative_surface_sto_min": float(np.nanmin(negative_surface)),
        "negative_surface_sto_max": float(np.nanmax(negative_surface)),
        "positive_surface_sto_min": float(np.nanmin(positive_surface)),
        "positive_surface_sto_max": float(np.nanmax(positive_surface)),
    }


def main() -> None:
    paths = ga.get_paths(ROOT)
    bundle = ga.load_legacy_bundle(paths)
    model = ga.build_legacy_model_inputs(bundle)
    anode = build_anode_detail(bundle, model)
    cathode = build_cathode_detail(bundle, model)

    comparison_model = dict(model)
    comparison_model["y_exp"] = cathode["grid"]
    comparison_model["up_exp"] = cathode["equilibrium"]
    anode_for_window = {
        "x": anode["grid"],
        "equilibrium": anode["equilibrium"],
    }
    stage, _ = omc.endpoint_window(comparison_model, anode_for_window)

    rows = []
    rows.extend(branch_metrics("Anode", anode, stage["x0"], stage["x100"]))
    rows.extend(branch_metrics("Cathode", cathode, stage["y100"], stage["y0"]))
    branch_table = pd.DataFrame(rows)

    for electrode in ("Anode", "Cathode"):
        h = branch_table[
            (branch_table["Electrode"] == electrode)
            & (branch_table["Comparison"] == "lithiation_minus_delithiation")
            & (branch_table["Region"] == "operating_window_central_10_90")
        ].iloc[0]
        repeat = branch_table[
            (branch_table["Electrode"] == electrode)
            & branch_table["Comparison"].str.startswith("replicate_difference")
        ]["MAE_mV"].mean()
        branch_table.loc[
            (branch_table["Electrode"] == electrode)
            & (branch_table["Comparison"] == "lithiation_minus_delithiation"),
            "mean_replicate_MAE_mV",
        ] = repeat
        branch_table.loc[
            (branch_table["Electrode"] == electrode)
            & (branch_table["Comparison"] == "lithiation_minus_delithiation"),
            "hysteresis_to_repeatability_ratio",
        ] = float(h["MAE_mV"] / repeat)

    branch_table.to_csv(
        RESULTS / "electrode_branch_difference_metrics.csv", index=False, encoding="utf-8-sig"
    )

    exp = omc.load_old_dynamic_data(model["dynamic"], model["q_meas"])
    cases = {
        "No hysteresis": (False, False),
        "Anode only": (True, False),
        "Cathode only": (False, True),
        "Both electrodes": (True, True),
    }
    colors = {
        "No hysteresis": "#7f7f7f",
        "Anode only": "#1f77b4",
        "Cathode only": "#ff7f0e",
        "Both electrodes": "#2ca02c",
    }
    simulations = {}
    dynamic_rows = []
    for name, (negative_hysteresis, positive_hysteresis) in cases.items():
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                print(f"{name} | {rate:g}C | {direction}", flush=True)
                sim = run_dfn(
                    rate,
                    charge,
                    stage,
                    comparison_model,
                    anode,
                    cathode,
                    negative_hysteresis,
                    positive_hysteresis,
                )
                simulations[(name, rate, charge)] = sim
                dynamic_rows.append(
                    {
                        "Model": name,
                        "C_rate": rate,
                        "Direction": direction,
                        **omc.old_time_metrics(exp[(rate, charge)], sim),
                    }
                )
    dynamic = pd.DataFrame(dynamic_rows)
    dynamic.to_csv(
        RESULTS / "dynamic_hysteresis_model_comparison.csv", index=False, encoding="utf-8-sig"
    )
    summary = (
        dynamic.groupby("Model", as_index=False)
        .agg(
            Mean_MAE_all_mV=("MAE_all_mV", "mean"),
            Mean_MAE_SOC10_70_mV=("MAE_SOC10_70_mV", "mean"),
            Mean_abs_end_time_error_min=("end_time_error_min", lambda x: np.mean(np.abs(x))),
        )
        .sort_values("Mean_MAE_SOC10_70_mV")
    )
    summary.to_csv(RESULTS / "dynamic_hysteresis_summary.csv", index=False, encoding="utf-8-sig")

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.5))
    for col, (electrode, detail, lo, hi) in enumerate(
        (
            ("Anode", anode, stage["x0"], stage["x100"]),
            ("Cathode", cathode, stage["y100"], stage["y0"]),
        )
    ):
        axes[0, col].plot(detail["grid"], detail["lithiation"], lw=2, label="Lithiation")
        axes[0, col].plot(detail["grid"], detail["delithiation"], lw=2, label="Delithiation")
        axes[0, col].axvspan(lo, hi, color="#999999", alpha=0.12, label="Endpoint operating window")
        axes[0, col].set(title=f"{electrode} directional OCP", xlabel="Stoichiometry", ylabel="OCP [V]")
        axes[0, col].grid(alpha=0.25)
        axes[0, col].legend(fontsize=8)
        gap = (detail["lithiation"] - detail["delithiation"]) * 1000.0
        axes[1, col].plot(detail["grid"], gap, color="#9467bd", lw=2)
        axes[1, col].axhline(0, color="black", lw=1)
        axes[1, col].axvspan(lo, hi, color="#999999", alpha=0.12)
        axes[1, col].set(title=f"{electrode} lithiation − delithiation", xlabel="Stoichiometry", ylabel="OCP difference [mV]")
        axes[1, col].grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(RESULTS / "electrode_branch_hysteresis.png", dpi=190)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.2))
    for col, rate in enumerate(ga.RATES):
        for row, charge in enumerate((True, False)):
            axis = axes[row, col]
            axis.plot(exp[(rate, charge)]["t_min"], exp[(rate, charge)]["V"], color="black", lw=2.3, label="Experiment")
            for name in cases:
                sim = simulations[(name, rate, charge)]
                metric = dynamic[
                    (dynamic["Model"] == name)
                    & np.isclose(dynamic["C_rate"], rate)
                    & (dynamic["Direction"] == ("Charge" if charge else "Discharge"))
                ].iloc[0]
                axis.plot(sim["t_min"], sim["V"], color=colors[name], lw=1.5, label=f"{name} ({metric['MAE_SOC10_70_mV']:.1f} mV)")
            axis.set(title=f"{rate:g}C {'Charge' if charge else 'Discharge'}", xlabel="Time [min]", ylabel="Terminal voltage [V]")
            axis.grid(alpha=0.25)
            if row == 0 and col == 0:
                axis.legend(fontsize=7.5)
    fig.suptitle("Endpoint window: electrode hysteresis model comparison", fontsize=15)
    fig.tight_layout()
    fig.savefig(RESULTS / "dynamic_electrode_hysteresis_comparison.png", dpi=190)
    plt.close(fig)

    pivot = dynamic.pivot_table(
        index=["C_rate", "Direction"], columns="Model", values="MAE_SOC10_70_mV"
    ).reset_index()
    fig, ax = plt.subplots(figsize=(11.3, 5.6))
    positions = np.arange(len(pivot))
    width = 0.20
    for offset, name in zip((-1.5 * width, -0.5 * width, 0.5 * width, 1.5 * width), cases):
        bars = ax.bar(positions + offset, pivot[name], width, color=colors[name], label=name)
        ax.bar_label(bars, fmt="%.1f", fontsize=8, padding=2)
    ax.set_xticks(positions, [f"{r:g}C\n{d}" for r, d in zip(pivot["C_rate"], pivot["Direction"])])
    ax.set_ylabel("MAE in experimental SOC 10–70% [mV]")
    ax.set_title("Does cathode hysteresis improve the DFN?")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(RESULTS / "dynamic_hysteresis_mae.png", dpi=190)
    plt.close(fig)

    manifest = {
        "new_260918_GITT_used": False,
        "window": "endpoint constrained; delta x and delta y fixed by measured full-cell capacity",
        "anode_sources": "legacy low-rate v1/v2",
        "cathode_sources": "legacy GITT 7-7/7-8 rest-end",
        "hysteresis_decay_rate_test_value": 3.2,
        "kinetics": "Ai2020 nominal",
        "pybamm_version": pybamm.__version__,
    }
    (RESULTS / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print("\nEndpoint stage")
    print(pd.Series(stage).to_string())
    print("\nBranch metrics")
    print(branch_table.to_string(index=False))
    print("\nDynamic summary")
    print(summary.to_string(index=False))
    print("\nSaved to", RESULTS)


if __name__ == "__main__":
    main()
