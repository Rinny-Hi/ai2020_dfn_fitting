"""Recheck legacy experimental OCP with two stoichiometry-window methods.

The new 260918 GITT data is deliberately excluded.  This script compares the
current endpoint-constrained Stage 1 window with a full-qOCV least-squares
window while keeping the source workbook, OCP curves, DFN, parameters, mesh,
solver, and dynamic data fixed.
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
from scipy.optimize import least_squares

import gitt_ocp_analysis as ga


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260920_legacy_ocp_window_recheck"
RESULTS.mkdir(parents=True, exist_ok=True)


def stage1_metrics(
    method: str,
    x0: float,
    y100: float,
    model: dict[str, Any],
) -> tuple[dict[str, Any], np.ndarray]:
    x100 = x0 + model["delta_x"]
    y0 = y100 + model["delta_y"]
    return explicit_window_metrics(method, x0, x100, y100, y0, model)


def explicit_window_metrics(
    method: str,
    x0: float,
    x100: float,
    y100: float,
    y0: float,
    model: dict[str, Any],
) -> tuple[dict[str, Any], np.ndarray]:
    soc = model["soc"]
    delta_x = x100 - x0
    delta_y = y0 - y100
    un = np.interp(x0 + soc * delta_x, model["x_legacy"], model["un_legacy"])
    up = np.interp(y0 - soc * delta_y, model["y_exp"], model["up_exp"])
    voltage = up - un
    error_mV = (voltage - model["v_qocv"]) * 1000.0
    central = (soc >= 0.02) & (soc <= 0.98)
    endpoint_error = np.array([error_mV[0], error_mV[-1]])
    x_lower = max(1e-6, float(model["x_legacy"].min()))
    x_upper = min(1.0 - model["delta_x"] - 1e-6, float(model["x_legacy"].max()) - model["delta_x"])
    y_lower = max(1e-6, float(model["y_exp"].min()))
    y_upper = min(1.0 - model["delta_y"] - 1e-6, float(model["y_exp"].max()) - model["delta_y"])
    tolerance = 1e-5
    bound_hits = []
    if x0 - x_lower <= tolerance:
        bound_hits.append("x0_lower")
    if x_upper - x0 <= tolerance:
        bound_hits.append("x0_upper")
    if y100 - y_lower <= tolerance:
        bound_hits.append("y100_lower")
    if y_upper - y100 <= tolerance:
        bound_hits.append("y100_upper")
    row = {
        "Method": method,
        "x0": x0,
        "x100": x100,
        "y100": y100,
        "y0": y0,
        "Delta_x": delta_x,
        "Delta_y": delta_y,
        "SOC_0_endpoint_error_mV": endpoint_error[0],
        "SOC_100_endpoint_error_mV": endpoint_error[1],
        "qOCV_RMSE_full_mV": float(np.sqrt(np.mean(error_mV**2))),
        "qOCV_MAE_full_mV": float(np.mean(np.abs(error_mV))),
        "qOCV_RMSE_2_98_mV": float(np.sqrt(np.mean(error_mV[central] ** 2))),
        "qOCV_MAE_2_98_mV": float(np.mean(np.abs(error_mV[central]))),
        "optimizer_boundary_hit": bool(bound_hits),
        "optimizer_boundary": ";".join(bound_hits) if bound_hits else "none",
    }
    return row, voltage


def solve_full_qocv_window(model: dict[str, Any]) -> tuple[dict[str, Any], np.ndarray]:
    x_lower = max(1e-6, float(model["x_legacy"].min()))
    x_upper = min(1.0 - model["delta_x"] - 1e-6, float(model["x_legacy"].max()) - model["delta_x"])
    y_lower = max(1e-6, float(model["y_exp"].min()))
    y_upper = min(1.0 - model["delta_y"] - 1e-6, float(model["y_exp"].max()) - model["delta_y"])

    def residual(offsets: np.ndarray) -> np.ndarray:
        x0, y100 = offsets
        _, voltage = stage1_metrics("full_qocv_least_squares", x0, y100, model)
        return voltage - model["v_qocv"]

    endpoint_row, _ = ga.solve_stage1_case(
        "Legacy", model["x_legacy"], model["un_legacy"], model
    )
    initial = np.array([endpoint_row["x0"], endpoint_row["y100"]])
    result = least_squares(
        residual,
        x0=initial,
        bounds=([x_lower, y_lower], [x_upper, y_upper]),
        xtol=1e-13,
        ftol=1e-13,
        gtol=1e-13,
        max_nfev=500,
    )
    row, voltage = stage1_metrics(
        "full_qocv_least_squares", float(result.x[0]), float(result.x[1]), model
    )
    row["Solver_success"] = bool(result.success)
    row["Solver_message"] = str(result.message)
    return row, voltage


def interp_time(time: np.ndarray, values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    order = np.argsort(time)
    time = np.asarray(time, dtype=float)[order]
    values = np.asarray(values, dtype=float)[order]
    time, index = np.unique(time, return_index=True)
    return np.interp(grid, time, values[index])


def time_metrics(
    experiment: dict[str, np.ndarray],
    simulation: dict[str, np.ndarray],
) -> dict[str, float]:
    t0 = max(float(experiment["t_min"].min()), float(simulation["t_min"].min()))
    t1 = min(float(experiment["t_min"].max()), float(simulation["t_min"].max()))
    grid = np.linspace(t0, t1, 400)
    exp_v = interp_time(experiment["t_min"], experiment["V"], grid)
    sim_v = interp_time(simulation["t_min"], simulation["V"], grid)
    error = (sim_v - exp_v) * 1000.0

    mask = (experiment["SOC"] >= 0.10) & (experiment["SOC"] <= 0.70)
    mid0 = max(float(experiment["t_min"][mask].min()), float(simulation["t_min"].min()))
    mid1 = min(float(experiment["t_min"][mask].max()), float(simulation["t_min"].max()))
    mid_grid = np.linspace(mid0, mid1, 300)
    exp_mid = interp_time(experiment["t_min"], experiment["V"], mid_grid)
    sim_mid = interp_time(simulation["t_min"], simulation["V"], mid_grid)
    mid_error = (sim_mid - exp_mid) * 1000.0
    return {
        "RMSE_time_full_mV": float(np.sqrt(np.mean(error**2))),
        "MAE_time_full_mV": float(np.mean(np.abs(error))),
        "RMSE_time_SOC10_70_mV": float(np.sqrt(np.mean(mid_error**2))),
        "MAE_time_SOC10_70_mV": float(np.mean(np.abs(mid_error))),
    }


def save_figures(
    model: dict[str, Any],
    stage_voltages: dict[str, np.ndarray],
    simulations: dict[tuple[str, float, bool], dict[str, np.ndarray]],
    experiment: dict[tuple[float, bool], dict[str, np.ndarray]],
    dynamic: pd.DataFrame,
) -> None:
    colors = {
        "endpoint_constrained": "#1f77b4",
        "full_qocv_least_squares": "#ff7f0e",
        "slide_reported_window": "#2ca02c",
    }
    labels = {
        "endpoint_constrained": "Endpoint-constrained",
        "full_qocv_least_squares": "Full-qOCV least squares",
        "slide_reported_window": "Slide-reported window",
    }

    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    ax.plot(model["soc"], model["v_qocv"], color="black", lw=2.6, label="Measured C/20 qOCV")
    for method, voltage in stage_voltages.items():
        ax.plot(model["soc"], voltage, color=colors[method], lw=2.0, label=labels[method])
    ax.set(xlabel="SOC", ylabel="Cell OCV [V]", title="Legacy OCP: stoichiometry-window comparison")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(RESULTS / "legacy_qocv_window_comparison.png", dpi=190)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(15, 8.2), sharex=False)
    for col, rate in enumerate(ga.RATES):
        for row_index, charge in enumerate((True, False)):
            axis = axes[row_index, col]
            exp = experiment[(rate, charge)]
            axis.plot(exp["t_min"], exp["V"], color="black", lw=2.0, label="Experiment")
            for method in labels:
                sim = simulations[(method, rate, charge)]
                axis.plot(sim["t_min"], sim["V"], color=colors[method], lw=1.7, label=labels[method])
                metric = dynamic[
                    (dynamic["Method"] == method)
                    & np.isclose(dynamic["C_rate"], rate)
                    & (dynamic["direction"] == ("Charge" if charge else "Discharge"))
                ].iloc[0]
                text_offsets = {
                    "endpoint_constrained": 0.04,
                    "full_qocv_least_squares": 0.12,
                    "slide_reported_window": 0.20,
                }
                axis.text(
                    0.02,
                    text_offsets[method],
                    f"{labels[method]} MAE10-70={metric['MAE_time_SOC10_70_mV']:.1f} mV",
                    transform=axis.transAxes,
                    color=colors[method],
                    fontsize=8.5,
                )
            axis.set_title(f"{rate:g}C {'Charge' if charge else 'Discharge'}")
            axis.set_xlabel("Time [min]")
            axis.set_ylabel("Terminal voltage [V]")
            axis.grid(alpha=0.25)
            if row_index == 0 and col == 0:
                axis.legend(fontsize=8)
    fig.suptitle("Legacy OCP dynamic baseline by stoichiometry-window method", fontsize=15)
    fig.tight_layout()
    fig.savefig(RESULTS / "legacy_dynamic_window_comparison.png", dpi=190)
    plt.close(fig)

    pivot = dynamic.pivot_table(
        index=["C_rate", "direction"],
        columns="Method",
        values="MAE_time_SOC10_70_mV",
    ).reset_index()
    fig, ax = plt.subplots(figsize=(9.2, 5.2))
    positions = np.arange(len(pivot))
    width = 0.25
    for offset, method in (
        (-width, "endpoint_constrained"),
        (0.0, "full_qocv_least_squares"),
        (width, "slide_reported_window"),
    ):
        bars = ax.bar(
            positions + offset,
            pivot[method],
            width=width,
            color=colors[method],
            label=labels[method],
        )
        ax.bar_label(bars, fmt="%.1f", padding=2, fontsize=9)
    ax.set_xticks(
        positions,
        [f"{rate:g}C\n{direction}" for rate, direction in zip(pivot["C_rate"], pivot["direction"])],
    )
    ax.set_ylabel("Time-domain MAE, SOC 10–70% [mV]")
    ax.set_title("Legacy OCP dynamic error")
    ax.grid(axis="y", alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(RESULTS / "legacy_dynamic_mae_10_70.png", dpi=190)
    plt.close(fig)


def main() -> None:
    paths = ga.get_paths(ROOT)
    bundle = ga.load_legacy_bundle(paths)
    model = ga.build_legacy_model_inputs(bundle)

    endpoint_row, endpoint_voltage = ga.solve_stage1_case(
        "Legacy", model["x_legacy"], model["un_legacy"], model
    )
    endpoint_row = {
        "Method": "endpoint_constrained",
        **{key: value for key, value in endpoint_row.items() if key != "Case"},
    }
    full_row, full_voltage = solve_full_qocv_window(model)
    slide_row, slide_voltage = explicit_window_metrics(
        "slide_reported_window",
        0.0001,
        0.7723,
        0.46062,
        0.9665,
        model,
    )
    slide_row["Solver_success"] = True
    slide_row["Solver_message"] = "Values transcribed from the supplied prior-result slide"

    stage1 = pd.DataFrame([endpoint_row, full_row, slide_row])
    stage1.to_csv(RESULTS / "legacy_stage1_window_comparison.csv", index=False, encoding="utf-8-sig")

    methods = {
        "endpoint_constrained": endpoint_row,
        "full_qocv_least_squares": full_row,
        "slide_reported_window": slide_row,
    }
    experiment = ga.load_dynamic_data(model["dynamic"])
    simulations: dict[tuple[str, float, bool], dict[str, np.ndarray]] = {}
    rows = []
    for method, stage in methods.items():
        for rate in ga.RATES:
            for charge in (True, False):
                direction = "Charge" if charge else "Discharge"
                print(f"{method}: {rate:g}C {direction}", flush=True)
                sim = ga.run_cc_dfn(
                    rate,
                    charge,
                    stage,
                    model["x_legacy"],
                    model["un_legacy"],
                    model,
                )
                simulations[(method, rate, charge)] = sim
                metrics = ga._dynamic_error_metrics(experiment[(rate, charge)], sim)
                metrics.update(time_metrics(experiment[(rate, charge)], sim))
                rows.append(
                    {
                        "Method": method,
                        "C_rate": rate,
                        "direction": direction,
                        **metrics,
                        "simulation_success": True,
                    }
                )
    dynamic = pd.DataFrame(rows)
    dynamic.to_csv(RESULTS / "legacy_dynamic_window_comparison.csv", index=False, encoding="utf-8-sig")

    save_figures(
        model,
        {
            "endpoint_constrained": endpoint_voltage,
            "full_qocv_least_squares": full_voltage,
            "slide_reported_window": slide_voltage,
        },
        simulations,
        experiment,
        dynamic,
    )

    manifest = {
        "new_gitt_used": False,
        "source_workbook": str(paths.legacy_book),
        "negative_ocp": "legacy low-rate v1/v2 mean",
        "positive_ocp": "legacy cathode GITT 7-7/7-8 rest-end mean",
        "dynamic_parameters": "Ai2020 source-corrected nominal values; no Stage 2-4 fitting",
        "methods": list(methods),
    }
    (RESULTS / "analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print("\nStage 1")
    print(stage1.to_string(index=False))
    print("\nDynamic")
    print(dynamic.to_string(index=False))
    print("\nSaved to", RESULTS)


if __name__ == "__main__":
    main()
