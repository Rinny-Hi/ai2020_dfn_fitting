"""Compare old/new/nominal electrode OCP source combinations against C/50 full-cell data."""

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
from scipy.optimize import differential_evolution, minimize, minimize_scalar

import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga


warnings.filterwarnings("ignore", message="Workbook contains no default style")

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260921_c50_ocp_source_combinations"
RESULTS.mkdir(parents=True, exist_ok=True)

JUNE_FILES = [
    Path(r"C:\Users\user\OneDrive\01. 적응형 충전 프로토콜\4. Enertech 파우치셀 실험데이터\260610 파우치셀 데이터\260610 저율데이터 1.xlsx"),
    Path(r"C:\Users\user\OneDrive\01. 적응형 충전 프로토콜\4. Enertech 파우치셀 실험데이터\260610 파우치셀 데이터\260610 저율데이터 2.xlsx"),
]
NEW_ANODE_FILES = [
    Path(r"C:\Users\user\OneDrive\01. 적응형 충전 프로토콜\4. Enertech 파우치셀 실험데이터\260918 GITT 실험데이터\260918 Anode GITT v1 1-5.xlsx"),
    Path(r"C:\Users\user\OneDrive\01. 적응형 충전 프로토콜\4. Enertech 파우치셀 실험데이터\260918 GITT 실험데이터\260918 Anode GITT v2 1-6.xlsx"),
    Path(r"C:\Users\user\OneDrive\01. 적응형 충전 프로토콜\4. Enertech 파우치셀 실험데이터\260918 GITT 실험데이터\260918 Anode GITT v3 1-7.xlsx"),
]
NEW_CATHODE_FILES = [
    Path(r"C:\Users\user\OneDrive\01. 적응형 충전 프로토콜\4. Enertech 파우치셀 실험데이터\260918 GITT 실험데이터\260918 Cathode GITT v1 2-6.xlsx"),
    Path(r"C:\Users\user\OneDrive\01. 적응형 충전 프로토콜\4. Enertech 파우치셀 실험데이터\260918 GITT 실험데이터\260918 Cathode GITT v2 2-7.xlsx"),
    Path(r"C:\Users\user\OneDrive\01. 적응형 충전 프로토콜\4. Enertech 파우치셀 실험데이터\260918 GITT 실험데이터\260918 Cathode GITT v3 2-8.xlsx"),
]


def clean_curve(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return ga._clean_curve(np.asarray(x, dtype=float), np.asarray(y, dtype=float))


def interp_curve(curve: tuple[np.ndarray, np.ndarray], grid: np.ndarray) -> np.ndarray:
    return np.interp(grid, curve[0], curve[1])


def fullcell_c50_target(soc: np.ndarray) -> dict[str, object]:
    charge_curves: list[np.ndarray] = []
    discharge_curves: list[np.ndarray] = []
    charge_caps: list[float] = []
    discharge_caps: list[float] = []
    source_rows: list[dict[str, object]] = []

    for file_name in JUNE_FILES:
        frame = pd.read_excel(file_name, sheet_name="record")
        step_type = frame["Step Type"].astype(str)
        group = step_type.ne(step_type.shift()).cumsum()
        branches: dict[str, list[tuple[float, np.ndarray]]] = {"CC Chg": [], "CC DChg": []}
        for _, part in frame.groupby(group):
            kind = str(part["Step Type"].iloc[0])
            if kind not in branches:
                continue
            q = pd.to_numeric(part["Capacity(Ah)"], errors="coerce").to_numpy(dtype=float)
            v = pd.to_numeric(part["Voltage(V)"], errors="coerce").to_numpy(dtype=float)
            current = float(pd.to_numeric(part["Current(A)"], errors="coerce").median())
            valid = np.isfinite(q) & np.isfinite(v)
            q, v = q[valid], v[valid]
            if len(q) < 50 or float(np.max(q)) < 2.2 or abs(current) > 0.06:
                continue
            qmax = float(np.max(q))
            z = q / qmax if kind == "CC Chg" else 1.0 - q / qmax
            curve = clean_curve(z, v)
            branches[kind].append((qmax, interp_curve(curve, soc)))

        # First discharge follows the partial initial charge. Match complete cycles 2-5.
        selected_charge = branches["CC Chg"]
        selected_discharge = branches["CC DChg"][1:]
        for kind, selected in (("charge", selected_charge), ("discharge", selected_discharge)):
            for index, (capacity, values) in enumerate(selected, start=1):
                source_rows.append(
                    {
                        "file": file_name.name,
                        "direction": kind,
                        "selected_cycle_index": index,
                        "capacity_Ah": capacity,
                    }
                )
                if kind == "charge":
                    charge_curves.append(values)
                    charge_caps.append(capacity)
                else:
                    discharge_curves.append(values)
                    discharge_caps.append(capacity)

    charge = np.mean(charge_curves, axis=0)
    discharge = np.mean(discharge_curves, axis=0)
    qocv = 0.5 * (charge + discharge)
    qcell = 0.5 * (float(np.mean(charge_caps)) + float(np.mean(discharge_caps)))
    pd.DataFrame(source_rows).to_csv(
        RESULTS / "c50_selected_branches.csv", index=False, encoding="utf-8-sig"
    )
    return {
        "charge": charge,
        "discharge": discharge,
        "qocv": qocv,
        "qcell": qcell,
        "charge_capacity_mean": float(np.mean(charge_caps)),
        "discharge_capacity_mean": float(np.mean(discharge_caps)),
        "n_charge": len(charge_caps),
        "n_discharge": len(discharge_caps),
    }


def new_halfcell_capacity_qc(files: list[Path], electrode: str) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for file_name in files:
        frame = pd.read_excel(file_name, sheet_name="step")
        seconds = pd.to_timedelta(frame["Step Time"].astype(str), errors="coerce").dt.total_seconds()
        frame = frame.assign(seconds=seconds)
        pulse = frame[
            frame["Step Type"].isin(["CC Chg", "CC DChg"])
            & frame.seconds.between(500, 700)
            & (frame["Capacity(Ah)"] > 0)
        ]
        sums = pulse.groupby("Step Type")["Capacity(Ah)"].sum()
        rows.append(
            {
                "electrode": electrode,
                "source": file_name.stem,
                "charge_pulse_count": int((pulse["Step Type"] == "CC Chg").sum()),
                "discharge_pulse_count": int((pulse["Step Type"] == "CC DChg").sum()),
                "charge_capacity_mAh": float(sums.get("CC Chg", np.nan) * 1000),
                "discharge_capacity_mAh": float(sums.get("CC DChg", np.nan) * 1000),
                "last_voltage_V": float(frame["End Voltage(V)"].iloc[-1]),
                "last_step_type": str(frame["Step Type"].iloc[-1]),
                "last_rest_minutes": float(seconds.iloc[-1] / 60),
            }
        )
    return pd.DataFrame(rows)


def extract_new_cathode_delith() -> tuple[dict[str, object], pd.DataFrame]:
    replicates: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    qc_rows: list[dict[str, object]] = []
    for file_name in NEW_CATHODE_FILES:
        frame = pd.read_excel(file_name, sheet_name="step")
        seconds = pd.to_timedelta(frame["Step Time"].astype(str), errors="coerce").dt.total_seconds()
        frame = frame.assign(seconds=seconds)
        cycle = frame[frame["Cycle Index"] == 2].reset_index(drop=True)
        first_rest_candidates = cycle[
            (cycle["Step Type"] == "Rest") & cycle.seconds.between(7100, 7300)
        ]
        first_rest_index = int(first_rest_candidates.index[0])
        gitt = cycle.loc[first_rest_index:].reset_index(drop=True)

        cumulative = {"CC Chg": 0.0, "CC DChg": 0.0}
        points = {"CC Chg": [(0.0, float(gitt.iloc[0]["End Voltage(V)"]))], "CC DChg": []}
        last_branch: str | None = None
        for _, row in gitt.iterrows():
            kind = str(row["Step Type"])
            if kind in cumulative and 500 <= float(row["seconds"]) <= 700:
                cumulative[kind] += float(row["Capacity(Ah)"])
                last_branch = kind
            elif kind == "Rest" and last_branch is not None and float(row["seconds"]) >= 6000:
                points[last_branch].append((cumulative[last_branch], float(row["End Voltage(V)"])))

        delith = np.asarray(points["CC Chg"], dtype=float)
        qmax = float(delith[:, 0].max())
        y = 1.0 - delith[:, 0] / qmax
        replicates[file_name.stem] = clean_curve(y, delith[:, 1])
        qc_rows.append(
            {
                "source": file_name.stem,
                "delithiation_capacity_mAh": 1000 * cumulative["CC Chg"],
                "lithiation_capacity_mAh": 1000 * cumulative["CC DChg"],
                "reverse_branch_completion_fraction": cumulative["CC DChg"] / cumulative["CC Chg"],
                "last_relaxed_voltage_V": float(frame["End Voltage(V)"].iloc[-1]),
            }
        )

    grid, mean = ga._mean_overlap(list(replicates.values()), n=3001)
    return {
        "grid": grid,
        "equilibrium": mean,
        "lithiation": mean.copy(),  # reverse branch is incomplete; do not claim it is measured
        "delithiation": mean.copy(),
        "replicates": replicates,
        "directional_complete": False,
        "absolute_anchor": False,
    }, pd.DataFrame(qc_rows)


def nominal_detail(function) -> dict[str, object]:
    grid = np.linspace(0.001, 0.999, 3001)
    values = np.array(
        [float(function(pybamm.Scalar(float(value))).evaluate()) for value in grid],
        dtype=float,
    )
    return {
        "grid": grid,
        "equilibrium": values,
        "lithiation": values.copy(),
        "delithiation": values.copy(),
        "directional_complete": False,
        "absolute_anchor": True,
    }


def new_anode_detail() -> dict[str, object]:
    paths = ga.get_paths(ROOT)
    gitt = ga.load_new_gitt(paths)
    built = ga.build_new_ocp(gitt)
    curves = built["curves"]

    def mean_branch(branch: str) -> tuple[np.ndarray, np.ndarray]:
        return ga._mean_overlap([curves["v2"][branch], curves["v3"][branch]], n=3001)

    lith = mean_branch("lithiation")
    delith = mean_branch("delithiation")
    lo = max(lith[0].min(), delith[0].min())
    hi = min(lith[0].max(), delith[0].max())
    grid = np.linspace(lo, hi, 3001)
    lith_v = interp_curve(lith, grid)
    delith_v = interp_curve(delith, grid)
    return {
        "grid": grid,
        "equilibrium": 0.5 * (lith_v + delith_v),
        "lithiation": lith_v,
        "delithiation": delith_v,
        "directional_complete": True,
        "absolute_anchor": False,
    }


def old_detail(detail: dict[str, object]) -> dict[str, object]:
    return {
        "grid": np.asarray(detail["grid"], dtype=float),
        "equilibrium": np.asarray(detail["equilibrium"], dtype=float),
        "lithiation": np.asarray(detail["lithiation"], dtype=float),
        "delithiation": np.asarray(detail["delithiation"], dtype=float),
        "directional_complete": True,
        "absolute_anchor": True,
    }


def predict(
    z: np.ndarray,
    soc: np.ndarray,
    anode: dict[str, object],
    cathode: dict[str, object],
    branch: str = "equilibrium",
) -> np.ndarray:
    x0, x100, y100, y0 = map(float, z)
    x = x0 + soc * (x100 - x0)
    y = y0 - soc * (y0 - y100)
    if branch == "charge":
        un_key, up_key = "lithiation", "delithiation"
    elif branch == "discharge":
        un_key, up_key = "delithiation", "lithiation"
    else:
        un_key = up_key = "equilibrium"
    un = np.interp(x, anode["grid"], anode[un_key])
    up = np.interp(y, cathode["grid"], cathode[up_key])
    return up - un


def fit_window(
    target: np.ndarray,
    soc: np.ndarray,
    qcell: float,
    anode: dict[str, object],
    cathode: dict[str, object],
    seed: int,
) -> tuple[np.ndarray, bool, str]:
    xmin, xmax = float(np.min(anode["grid"])), float(np.max(anode["grid"]))
    ymin, ymax = float(np.min(cathode["grid"])), float(np.max(cathode["grid"]))
    min_width = qcell / 6.0
    max_width = qcell / 2.4
    mask = (soc >= 0.02) & (soc <= 0.98)

    def voltage_mse(z: np.ndarray) -> float:
        err = predict(z, soc, anode, cathode) - target
        return float(np.mean(err[mask] ** 2))

    def violation(z: np.ndarray) -> float:
        x0, x100, y100, y0 = z
        dx, dy = x100 - x0, y0 - y100
        endpoint = predict(z, soc[[0, -1]], anode, cathode) - target[[0, -1]]
        terms = [
            max(0.0, min_width - dx), max(0.0, dx - max_width),
            max(0.0, min_width - dy), max(0.0, dy - max_width),
            max(0.0, abs(float(endpoint[0])) - 0.010),
            max(0.0, abs(float(endpoint[1])) - 0.010),
        ]
        return float(np.sum(np.square(terms)))

    def penalized(z: np.ndarray) -> float:
        return voltage_mse(z) + 1e4 * violation(z)

    bounds = [(xmin, xmax), (xmin, xmax), (ymin, ymax), (ymin, ymax)]
    de = differential_evolution(
        penalized,
        bounds=bounds,
        seed=seed,
        popsize=16,
        maxiter=350,
        tol=1e-10,
        polish=False,
        workers=1,
    )

    constraints = [
        {"type": "ineq", "fun": lambda z: (z[1] - z[0]) - min_width},
        {"type": "ineq", "fun": lambda z: max_width - (z[1] - z[0])},
        {"type": "ineq", "fun": lambda z: (z[3] - z[2]) - min_width},
        {"type": "ineq", "fun": lambda z: max_width - (z[3] - z[2])},
        {"type": "ineq", "fun": lambda z: 0.010 - abs(float(predict(z, soc[[0]], anode, cathode)[0] - target[0]))},
        {"type": "ineq", "fun": lambda z: 0.010 - abs(float(predict(z, soc[[-1]], anode, cathode)[0] - target[-1]))},
    ]
    local = minimize(
        voltage_mse,
        x0=de.x,
        method="SLSQP",
        bounds=bounds,
        constraints=constraints,
        options={"maxiter": 3000, "ftol": 1e-14, "disp": False},
    )
    z = local.x if violation(local.x) <= violation(de.x) + 1e-12 else de.x
    feasible = violation(z) < 1e-10
    return np.asarray(z, dtype=float), bool(feasible), str(local.message)


def main() -> None:
    paths = ga.get_paths(ROOT)
    bundle = ga.load_legacy_bundle(paths)
    model = ga.build_legacy_model_inputs(bundle)
    target_soc = np.linspace(0.0, 1.0, 1001)
    target = fullcell_c50_target(target_soc)

    old_anode = old_detail(ehq.build_anode_detail(bundle, model))
    old_cathode = old_detail(ehq.build_cathode_detail(bundle, model))
    new_anode = new_anode_detail()
    new_cathode, new_cathode_qc = extract_new_cathode_delith()
    params = pybamm.ParameterValues("Ai2020")
    nominal_anode = nominal_detail(params["Negative electrode OCP [V]"])
    nominal_cathode = nominal_detail(params["Positive electrode OCP [V]"])

    anodes = {"old": old_anode, "new": new_anode, "nominal": nominal_anode}
    cathodes = {"old": old_cathode, "new": new_cathode, "nominal": nominal_cathode}

    qc = pd.concat(
        [
            new_halfcell_capacity_qc(NEW_ANODE_FILES, "negative"),
            new_halfcell_capacity_qc(NEW_CATHODE_FILES, "positive"),
        ],
        ignore_index=True,
    )
    qc.to_csv(RESULTS / "new_gitt_capacity_qc.csv", index=False, encoding="utf-8-sig")
    new_cathode_qc.to_csv(
        RESULTS / "new_cathode_completion_qc.csv", index=False, encoding="utf-8-sig"
    )

    rows: list[dict[str, object]] = []
    curves: dict[tuple[str, str], dict[str, np.ndarray]] = {}
    for ai, (anode_name, anode) in enumerate(anodes.items()):
        for ci, (cathode_name, cathode) in enumerate(cathodes.items()):
            print(f"Fitting {anode_name}/{cathode_name}", flush=True)
            z, feasible, message = fit_window(
                np.asarray(target["qocv"]),
                target_soc,
                float(target["qcell"]),
                anode,
                cathode,
                seed=260921 + 10 * ai + ci,
            )
            eq = predict(z, target_soc, anode, cathode, "equilibrium")
            raw_charge = predict(z, target_soc, anode, cathode, "charge")
            raw_discharge = predict(z, target_soc, anode, cathode, "discharge")
            mask_10_90 = (target_soc >= 0.10) & (target_soc <= 0.90)

            def directional_objective(scale: float) -> float:
                scaled_charge = eq + scale * (raw_charge - eq)
                scaled_discharge = eq + scale * (raw_discharge - eq)
                return float(
                    np.mean((scaled_charge[mask_10_90] - target["charge"][mask_10_90]) ** 2)
                    + np.mean((scaled_discharge[mask_10_90] - target["discharge"][mask_10_90]) ** 2)
                )

            scale_result = minimize_scalar(
                directional_objective, bounds=(0.0, 1.0), method="bounded"
            )
            hysteresis_scale = float(scale_result.x)
            charge = eq + hysteresis_scale * (raw_charge - eq)
            discharge = eq + hysteresis_scale * (raw_discharge - eq)
            curves[(anode_name, cathode_name)] = {
                "equilibrium": eq,
                "charge": charge,
                "discharge": discharge,
            }
            mask_2_98 = (target_soc >= 0.02) & (target_soc <= 0.98)
            eq_error = (eq - target["qocv"]) * 1000
            charge_error = (charge - target["charge"]) * 1000
            discharge_error = (discharge - target["discharge"]) * 1000
            raw_charge_error = (raw_charge - target["charge"]) * 1000
            raw_discharge_error = (raw_discharge - target["discharge"]) * 1000
            dx, dy = z[1] - z[0], z[3] - z[2]
            primary_eligible = feasible and cathode_name != "new"
            rows.append(
                {
                    "anode_source": anode_name,
                    "cathode_source": cathode_name,
                    "x0": z[0], "x100": z[1], "y100": z[2], "y0": z[3],
                    "delta_x": dx, "delta_y": dy,
                    "Qn_Ah": float(target["qcell"]) / dx,
                    "Qp_Ah": float(target["qcell"]) / dy,
                    "qOCV_MAE_2_98_mV": float(np.mean(np.abs(eq_error[mask_2_98]))),
                    "qOCV_RMSE_2_98_mV": float(np.sqrt(np.mean(eq_error[mask_2_98] ** 2))),
                    "qOCV_MAE_10_90_mV": float(np.mean(np.abs(eq_error[mask_10_90]))),
                    "SOC0_error_mV": float(eq_error[0]),
                    "SOC100_error_mV": float(eq_error[-1]),
                    "max_endpoint_abs_mV": float(max(abs(eq_error[0]), abs(eq_error[-1]))),
                    "optimal_hysteresis_scale": hysteresis_scale,
                    "C50_raw_branch_mean_MAE_10_90_mV": float(
                        0.5 * (np.mean(np.abs(raw_charge_error[mask_10_90])) + np.mean(np.abs(raw_discharge_error[mask_10_90])))
                    ),
                    "C50_charge_MAE_10_90_mV": float(np.mean(np.abs(charge_error[mask_10_90]))),
                    "C50_discharge_MAE_10_90_mV": float(np.mean(np.abs(discharge_error[mask_10_90]))),
                    "C50_branch_mean_MAE_10_90_mV": float(
                        0.5 * (np.mean(np.abs(charge_error[mask_10_90])) + np.mean(np.abs(discharge_error[mask_10_90])))
                    ),
                    "solver_feasible": feasible,
                    "primary_eligible": primary_eligible,
                    "anode_absolute_anchor": bool(anode["absolute_anchor"]),
                    "cathode_absolute_anchor": bool(cathode["absolute_anchor"]),
                    "anode_directional_complete": bool(anode["directional_complete"]),
                    "cathode_directional_complete": bool(cathode["directional_complete"]),
                    "solver_message": message,
                }
            )

    summary = pd.DataFrame(rows).sort_values(
        ["primary_eligible", "qOCV_RMSE_2_98_mV", "C50_branch_mean_MAE_10_90_mV"],
        ascending=[False, True, True],
    )
    summary.to_csv(RESULTS / "ocp_source_combination_summary.csv", index=False, encoding="utf-8-sig")

    eligible = summary[summary["primary_eligible"]].copy()
    selected = eligible.sort_values(
        ["qOCV_RMSE_2_98_mV", "C50_branch_mean_MAE_10_90_mV"]
    ).iloc[0]
    selected_key = (str(selected["anode_source"]), str(selected["cathode_source"]))

    # Heatmaps.
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    for ax, value, title in (
        (axes[0], "qOCV_RMSE_2_98_mV", "C/50 qOCV RMSE [mV]"),
        (axes[1], "C50_branch_mean_MAE_10_90_mV", "C/50 branch mean MAE [mV]"),
    ):
        matrix = summary.pivot(index="anode_source", columns="cathode_source", values=value).reindex(
            index=["old", "new", "nominal"], columns=["old", "new", "nominal"]
        )
        image = ax.imshow(matrix.to_numpy(), cmap="viridis_r")
        for i in range(3):
            for j in range(3):
                label = f"{matrix.iloc[i, j]:.1f}"
                if matrix.index[i] == selected_key[0] and matrix.columns[j] == selected_key[1]:
                    label += "*"
                ax.text(j, i, label, ha="center", va="center", color="white" if matrix.iloc[i, j] > matrix.to_numpy().mean() else "black")
        ax.set_xticks(range(3), matrix.columns)
        ax.set_yticks(range(3), matrix.index)
        ax.set_xlabel("Cathode OCP source")
        ax.set_ylabel("Anode OCP source")
        ax.set_title(title)
        fig.colorbar(image, ax=ax, shrink=0.8)
    fig.suptitle("Old / new / nominal OCP source combinations (* selected eligible)")
    fig.tight_layout()
    fig.savefig(RESULTS / "ocp_source_combination_heatmaps.png", dpi=220)
    plt.close(fig)

    # Top eligible curves and the diagnostic new-cathode best case.
    eligible_keys = [
        (str(row.anode_source), str(row.cathode_source))
        for row in eligible.sort_values("qOCV_RMSE_2_98_mV").head(4).itertuples()
    ]
    diagnostic_new = summary[summary["cathode_source"] == "new"].sort_values("qOCV_RMSE_2_98_mV").iloc[0]
    diagnostic_key = (str(diagnostic_new["anode_source"]), "new")
    plot_keys = list(dict.fromkeys(eligible_keys + [diagnostic_key]))

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))
    observed = [target["qocv"], target["charge"], target["discharge"]]
    branch_names = ["equilibrium", "charge", "discharge"]
    titles = ["C/50 mean qOCV", "C/50 charge", "C/50 discharge"]
    for ax, obs, branch, title in zip(axes, observed, branch_names, titles):
        ax.plot(target_soc, obs, color="black", lw=2.7, label="Measured C/50")
        for key in plot_keys:
            is_selected = key == selected_key
            is_diagnostic = key == diagnostic_key and key not in eligible_keys
            ax.plot(
                target_soc,
                curves[key][branch],
                lw=2.5 if is_selected else 1.4,
                ls="--" if is_diagnostic else "-",
                alpha=1.0 if is_selected else 0.72,
                label=f"{key[0]}/{key[1]}" + (" selected" if is_selected else "") + (" incomplete cathode" if is_diagnostic else ""),
            )
        ax.set(title=title, xlabel="Normalized full-cell SOC", ylabel="Voltage [V]")
        ax.grid(alpha=0.25)
    axes[0].legend(fontsize=7)
    fig.suptitle("C/50 OCP source comparison")
    fig.tight_layout()
    fig.savefig(RESULTS / "top_ocp_source_curves.png", dpi=220)
    plt.close(fig)

    decision = {
        "target": {
            "source": "260610 C/50 cells 1 and 2, complete cycles 2-5",
            "qcell_Ah": float(target["qcell"]),
            "charge_capacity_mean_Ah": float(target["charge_capacity_mean"]),
            "discharge_capacity_mean_Ah": float(target["discharge_capacity_mean"]),
            "n_charge_branches": int(target["n_charge"]),
            "n_discharge_branches": int(target["n_discharge"]),
        },
        "selected": selected.to_dict(),
        "selection_rule": "minimum qOCV RMSE among endpoint-feasible candidates, excluding incomplete new-cathode reverse GITT",
        "new_cathode_reverse_branch_completion_mean": float(new_cathode_qc["reverse_branch_completion_fraction"].mean()),
        "new_cathode_last_relaxed_voltage_mean_V": float(new_cathode_qc["last_relaxed_voltage_V"].mean()),
    }
    (RESULTS / "decision.json").write_text(
        json.dumps(decision, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("\nCombination summary")
    print(summary.to_string(index=False))
    print("\nSelected", selected_key)
    print("Saved to", RESULTS)


if __name__ == "__main__":
    main()
