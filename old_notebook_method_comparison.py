"""Compare the 260818 notebook method with the current legacy-OCP recheck.

The new 260918 GITT measurements are intentionally excluded.  The comparison
separates the effect of the stoichiometry-window objective from the effect of
the one-state anode hysteresis model used in the old notebook.
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
from scipy.optimize import differential_evolution

import gitt_ocp_analysis as ga


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260920_old_notebook_method_comparison"
RESULTS.mkdir(parents=True, exist_ok=True)


def build_old_anode_ocp(
    bundle: dict[str, Any], model: dict[str, Any]
) -> dict[str, np.ndarray]:
    """Reproduce cell 10 of the 260818 notebook."""
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

    anode = bundle["anode_lowrate"]
    lith_curves = []
    delith_curves = []
    for source in ("v1", "v2"):
        q_chg, v_chg = ga._extract_qv_branch(
            anode, "source_version", source, "charge"
        )
        q_dch, v_dch = ga._extract_qv_branch(
            anode, "source_version", source, "discharge"
        )
        lith_curves.append(ga._clean_curve(q_dch / qn_disk, v_dch))
        span = q_chg.max() / qn_disk
        delith_curves.append(ga._clean_curve(span - q_chg / qn_disk, v_chg))

    x_lith, un_lith_0 = ga._mean_overlap(lith_curves, n=3000)
    x_delith, un_delith_0 = ga._mean_overlap(delith_curves, n=3000)
    x = np.linspace(max(x_lith.min(), x_delith.min()), min(x_lith.max(), x_delith.max()), 3000)
    un_lith = np.interp(x, x_lith, un_lith_0)
    un_delith = np.interp(x, x_delith, un_delith_0)
    return {
        "x": x,
        "equilibrium": 0.5 * (un_lith + un_delith),
        "lithiation": un_lith,
        "delithiation": un_delith,
    }


def window_voltage(
    x0: float,
    y100: float,
    model: dict[str, Any],
    anode: dict[str, np.ndarray],
) -> tuple[np.ndarray, tuple[float, float, float, float]]:
    x100 = x0 + model["delta_x"]
    y0 = y100 + model["delta_y"]
    x = x0 + model["soc"] * model["delta_x"]
    y = y0 - model["soc"] * model["delta_y"]
    voltage = np.interp(y, model["y_exp"], model["up_exp"]) - np.interp(
        x, anode["x"], anode["equilibrium"]
    )
    return voltage, (x0, x100, y100, y0)


def solve_old_window(
    model: dict[str, Any], anode: dict[str, np.ndarray]
) -> tuple[dict[str, Any], np.ndarray]:
    mask = (model["soc"] >= 0.02) & (model["soc"] <= 0.98)

    def objective(z: np.ndarray) -> float:
        voltage, (x0, x100, y100, y0) = window_voltage(z[0], z[1], model, anode)
        if not (0 < x0 < x100 < 1 and 0 < y100 < y0 < 1):
            return 100.0
        if x0 < anode["x"].min() or x100 > anode["x"].max():
            return 50.0
        if y100 < model["y_exp"].min() or y0 > model["y_exp"].max():
            return 50.0
        error = voltage[mask] - model["v_qocv"][mask]
        return float(np.sqrt(np.mean(error**2)))

    result = differential_evolution(
        objective,
        bounds=[(1e-4, 0.08), (0.35, 0.50)],
        seed=42,
        maxiter=600,
        tol=1e-10,
        polish=True,
    )
    voltage, (x0, x100, y100, y0) = window_voltage(
        float(result.x[0]), float(result.x[1]), model, anode
    )
    error_mV = (voltage - model["v_qocv"]) * 1000.0
    central = error_mV[mask]
    row = {
        "Method": "old_notebook_window",
        "x0": x0,
        "x100": x100,
        "y100": y100,
        "y0": y0,
        "qOCV_RMSE_2_98_mV": float(np.sqrt(np.mean(central**2))),
        "qOCV_MAE_2_98_mV": float(np.mean(np.abs(central))),
        "qOCV_bias_2_98_mV": float(np.mean(central)),
        "SOC_0_endpoint_error_mV": float(error_mV[0]),
        "SOC_100_endpoint_error_mV": float(error_mV[-1]),
        "optimizer_success": bool(result.success),
    }
    return row, voltage


def endpoint_window(
    model: dict[str, Any], anode: dict[str, np.ndarray]
) -> tuple[dict[str, Any], np.ndarray]:
    source, _ = ga.solve_stage1_case(
        "endpoint", anode["x"], anode["equilibrium"], model
    )
    voltage, (x0, x100, y100, y0) = window_voltage(
        float(source["x0"]), float(source["y100"]), model, anode
    )
    error_mV = (voltage - model["v_qocv"]) * 1000.0
    mask = (model["soc"] >= 0.02) & (model["soc"] <= 0.98)
    central = error_mV[mask]
    return {
        "Method": "current_endpoint_window",
        "x0": x0,
        "x100": x100,
        "y100": y100,
        "y0": y0,
        "qOCV_RMSE_2_98_mV": float(np.sqrt(np.mean(central**2))),
        "qOCV_MAE_2_98_mV": float(np.mean(np.abs(central))),
        "qOCV_bias_2_98_mV": float(np.mean(central)),
        "SOC_0_endpoint_error_mV": float(error_mV[0]),
        "SOC_100_endpoint_error_mV": float(error_mV[-1]),
        "optimizer_success": True,
    }, voltage


def hysteresis_parameter_values(
    model_inputs: dict[str, Any],
    stage: dict[str, Any],
    anode: dict[str, np.ndarray],
    charge: bool,
):
    values = model_inputs["params"].copy()
    x_init = stage["x0"] if charge else stage["x100"]
    y_init = stage["y0"] if charge else stage["y100"]
    values.update(
        {
            "Initial concentration in negative electrode [mol.m-3]": x_init
            * model_inputs["csn_max"],
            "Initial concentration in positive electrode [mol.m-3]": y_init
            * model_inputs["csp_max"],
            "Negative electrode OCP [V]": ga.pybamm_interp(
                anode["x"], anode["equilibrium"], "Experimental negative OCP"
            ),
            "Negative electrode lithiation OCP [V]": ga.pybamm_interp(
                anode["x"], anode["lithiation"], "Experimental negative lithiation OCP"
            ),
            "Negative electrode delithiation OCP [V]": ga.pybamm_interp(
                anode["x"], anode["delithiation"], "Experimental negative delithiation OCP"
            ),
            "Negative particle lithiation hysteresis decay rate": lambda sto, temp: 3.2
            + 0 * sto,
            "Negative particle delithiation hysteresis decay rate": lambda sto, temp: 3.2
            + 0 * sto,
            "Initial hysteresis state in negative electrode": -1.0 if charge else 1.0,
            "Positive electrode OCP [V]": ga.pybamm_interp(
                model_inputs["y_exp"], model_inputs["up_exp"], "Experimental positive OCP"
            ),
            "Initial temperature [K]": model_inputs["temperature"],
            "Ambient temperature [K]": model_inputs["temperature"],
        },
        check_already_exists=False,
    )
    return values


def run_old_hysteresis_dfn(
    rate: float,
    charge: bool,
    stage: dict[str, Any],
    model_inputs: dict[str, Any],
    anode: dict[str, np.ndarray],
) -> dict[str, np.ndarray]:
    model = pybamm.lithium_ion.DFN(
        {
            "open-circuit potential": ("one-state hysteresis", "single"),
            "thermal": "isothermal",
        }
    )
    step = (
        f"Charge at {rate}C until 4.2 V"
        if charge
        else f"Discharge at {rate}C until 3.0 V"
    )
    simulation = pybamm.Simulation(
        model,
        parameter_values=hysteresis_parameter_values(
            model_inputs, stage, anode, charge
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
    return {"t_min": time_s / 60.0, "V": voltage, "SOC": np.clip(soc, 0.0, 1.0)}


def old_time_metrics(
    experiment: dict[str, np.ndarray], simulation: dict[str, np.ndarray]
) -> dict[str, float]:
    te = np.asarray(experiment["t_min"], dtype=float)
    ve = np.asarray(experiment["V"], dtype=float)
    soc = np.asarray(experiment["SOC"], dtype=float)
    ts = np.asarray(simulation["t_min"], dtype=float)
    vs = np.asarray(simulation["V"], dtype=float)
    overlap = (
        np.isfinite(te)
        & np.isfinite(ve)
        & np.isfinite(soc)
        & (te >= max(np.nanmin(te), np.nanmin(ts)))
        & (te <= min(np.nanmax(te), np.nanmax(ts)))
    )
    order = np.argsort(ts)
    ts_u, index = np.unique(ts[order], return_index=True)
    vm = np.interp(te[overlap], ts_u, vs[order][index])
    error_mV = np.abs(vm - ve[overlap]) * 1000.0
    soc_use = soc[overlap]
    mid = (soc_use >= 0.10) & (soc_use <= 0.70)
    return {
        "MAE_all_mV": float(np.mean(error_mV)),
        "MAE_SOC10_70_mV": float(np.mean(error_mV[mid])),
        "N_all": int(overlap.sum()),
        "N_SOC10_70": int(mid.sum()),
        "end_time_error_min": float(ts[-1] - te[-1]),
    }


def load_old_dynamic_data(
    frame: pd.DataFrame, q_meas: float
) -> dict[tuple[float, bool], dict[str, np.ndarray]]:
    """Use the old notebook's Q/Q_MEAS SOC definition, not branch normalization."""
    output = {}
    for rate in ga.RATES:
        for charge in (True, False):
            direction = "charge" if charge else "discharge"
            selected = frame[
                np.isclose(pd.to_numeric(frame["rate_C"], errors="coerce"), rate)
                & (frame["direction"].astype(str).str.lower() == direction)
            ].copy()
            capacity = pd.to_numeric(selected["Q_Ah"], errors="coerce").to_numpy(dtype=float)
            capacity = capacity - capacity[0]
            soc = capacity / q_meas if charge else 1.0 - capacity / q_meas
            output[(rate, charge)] = {
                "t_min": pd.to_numeric(selected["t_min"], errors="coerce").to_numpy(dtype=float),
                "V": pd.to_numeric(selected["V_V"], errors="coerce").to_numpy(dtype=float),
                "I_A": pd.to_numeric(selected["I_A"], errors="coerce").to_numpy(dtype=float),
                # Keep the measured branch capacity itself.  Reconstructing
                # it later from clipped SOC can silently cap Q at q_meas.
                "Q_Ah": capacity,
                "SOC": np.clip(soc, 0.0, 1.0),
            }
    return output


def main() -> None:
    paths = ga.get_paths(ROOT)
    bundle = ga.load_legacy_bundle(paths)
    model = ga.build_legacy_model_inputs(bundle)
    anode = build_old_anode_ocp(bundle, model)

    endpoint, endpoint_voltage = endpoint_window(model, anode)
    old_window, old_voltage = solve_old_window(model, anode)
    stage_table = pd.DataFrame([endpoint, old_window])
    stage_table.to_csv(
        RESULTS / "stage1_method_comparison.csv", index=False, encoding="utf-8-sig"
    )

    exp = load_old_dynamic_data(model["dynamic"], model["q_meas"])
    cases = {
        "Current: endpoint + mean OCP": (endpoint, False),
        "Recommended test: endpoint + hysteresis": (endpoint, True),
        "Old window only: RMSE + mean OCP": (old_window, False),
        "Old notebook exact: RMSE + hysteresis": (old_window, True),
    }
    simulations = {}
    dynamic_rows = []
    for name, (stage, use_hysteresis) in cases.items():
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                print(f"{name} | {rate:g}C | {direction}", flush=True)
                if use_hysteresis:
                    sim = run_old_hysteresis_dfn(
                        rate, charge, stage, model, anode
                    )
                else:
                    sim = ga.run_cc_dfn(
                        rate,
                        charge,
                        stage,
                        anode["x"],
                        anode["equilibrium"],
                        model,
                    )
                simulations[(name, rate, charge)] = sim
                dynamic_rows.append(
                    {
                        "Method": name,
                        "C_rate": rate,
                        "Direction": direction,
                        **old_time_metrics(exp[(rate, charge)], sim),
                    }
                )
    dynamic = pd.DataFrame(dynamic_rows)
    dynamic.to_csv(
        RESULTS / "dynamic_method_comparison.csv", index=False, encoding="utf-8-sig"
    )

    colors = {
        "Current: endpoint + mean OCP": "#1f77b4",
        "Recommended test: endpoint + hysteresis": "#9467bd",
        "Old window only: RMSE + mean OCP": "#ff7f0e",
        "Old notebook exact: RMSE + hysteresis": "#2ca02c",
    }
    fig, ax = plt.subplots(figsize=(8.8, 5.2))
    ax.plot(model["soc"], model["v_qocv"], color="black", lw=2.5, label="Measured C/20 qOCV")
    ax.plot(model["soc"], endpoint_voltage, color=colors["Current: endpoint + mean OCP"], lw=2, label="Current endpoint window")
    ax.plot(model["soc"], old_voltage, color=colors["Old window only: RMSE + mean OCP"], lw=2, label="Old RMSE window")
    ax.set(xlabel="SOC", ylabel="Cell OCV [V]", title="Stoichiometry-window method: old vs current")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(RESULTS / "qocv_old_vs_current.png", dpi=190)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.2))
    for col, rate in enumerate(ga.RATES):
        for row, charge in enumerate((True, False)):
            axis = axes[row, col]
            axis.plot(exp[(rate, charge)]["t_min"], exp[(rate, charge)]["V"], color="black", lw=2.2, label="Experiment")
            for name in cases:
                sim = simulations[(name, rate, charge)]
                metric = dynamic[
                    (dynamic["Method"] == name)
                    & np.isclose(dynamic["C_rate"], rate)
                    & (dynamic["Direction"] == ("Charge" if charge else "Discharge"))
                ].iloc[0]
                axis.plot(sim["t_min"], sim["V"], color=colors[name], lw=1.6, label=f"{name} ({metric['MAE_SOC10_70_mV']:.1f} mV)")
            axis.set_title(f"{rate:g}C {'Charge' if charge else 'Discharge'}")
            axis.set_xlabel("Time [min]")
            axis.set_ylabel("Terminal voltage [V]")
            axis.grid(alpha=0.25)
            if row == 0 and col == 0:
                axis.legend(fontsize=7.5)
    fig.suptitle("Dynamic validation: separating window and hysteresis effects", fontsize=15)
    fig.tight_layout()
    fig.savefig(RESULTS / "dynamic_old_vs_current.png", dpi=190)
    plt.close(fig)

    pivot = dynamic.pivot_table(
        index=["C_rate", "Direction"], columns="Method", values="MAE_SOC10_70_mV"
    ).reset_index()
    fig, ax = plt.subplots(figsize=(10.5, 5.5))
    positions = np.arange(len(pivot))
    width = 0.20
    for offset, name in zip((-1.5 * width, -0.5 * width, 0.5 * width, 1.5 * width), cases):
        bars = ax.bar(positions + offset, pivot[name], width, color=colors[name], label=name)
        ax.bar_label(bars, fmt="%.1f", fontsize=8.5, padding=2)
    ax.set_xticks(positions, [f"{r:g}C\n{d}" for r, d in zip(pivot["C_rate"], pivot["Direction"])])
    ax.set_ylabel("MAE in experimental SOC 10–70% [mV]")
    ax.set_title("Why the old notebook looked better")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(RESULTS / "mae_old_vs_current.png", dpi=190)
    plt.close(fig)

    summary = (
        dynamic.groupby("Method", as_index=False)
        .agg(
            Mean_MAE_all_mV=("MAE_all_mV", "mean"),
            Mean_MAE_SOC10_70_mV=("MAE_SOC10_70_mV", "mean"),
        )
    )
    summary.to_csv(RESULTS / "dynamic_summary.csv", index=False, encoding="utf-8-sig")

    manifest = {
        "new_260918_GITT_used": False,
        "old_notebook": r"C:\Users\user\OneDrive\02. 현대자동차 미팅자료\★ 추후 제출 코드\260818_0단계~2단계_진행_결과.ipynb",
        "pybamm_old_notebook": "26.7.1.0",
        "pybamm_recheck": pybamm.__version__,
        "old_metric": "model interpolated on experimental time points; experimental SOC 10-70%",
        "old_dynamic_model": "DFN one-state hysteresis; nominal Ai2020 k",
    }
    (RESULTS / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print("\nStage 1")
    print(stage_table.to_string(index=False))
    print("\nDynamic")
    print(dynamic.to_string(index=False))
    print("\nSummary")
    print(summary.to_string(index=False))
    print("\nSaved to", RESULTS)


if __name__ == "__main__":
    main()
