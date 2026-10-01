"""Fit passed charge subsets and validate them on untouched discharge data.

Capacity is never part of the optimizer objective.  It is calculated only in
the post-fit reporting section.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution, dual_annealing, least_squares
from scipy.stats import qmc

import charge_physical_time_common as common
import gitt_ocp_analysis as ga


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "260927_charge_subset_fit_validation"
LSA = ROOT / "results" / "260927_charge_physical_time_lsa" / "subset_identifiability.csv"
CURRENT_10PCT_DSN = 2.381074e-14


def make_starts(parameters: tuple[str, ...]) -> list[tuple[str, np.ndarray]]:
    starts = [
        ("all-range GITT", common.values_to_normalized(common.INITIAL, parameters)),
        (
            "10pct-floor GITT",
            common.values_to_normalized(
                {**common.INITIAL, "Dsn": CURRENT_10PCT_DSN}, parameters
            ),
        ),
        ("bound center", np.full(len(parameters), 0.5)),
    ]
    sampler = qmc.LatinHypercube(d=len(parameters), seed=260927)
    for index, vector in enumerate(sampler.random(2), start=1):
        starts.append((f"log-space LHS {index}", vector))
    return starts


class ModelFitter:
    def __init__(self, context: common.Context, model_name: str):
        self.context = context
        self.model_name = model_name
        self.parameters = common.MODEL_PARAMETERS[model_name]
        self.result_dir = OUT / model_name
        self.result_dir.mkdir(parents=True, exist_ok=True)
        self.cache: dict[tuple[float, ...], tuple[np.ndarray, dict | None, dict, dict]] = {}
        self.history: list[dict] = []
        self.final_rows: list[dict] = []
        self.final_vectors: dict[str, np.ndarray] = {}
        self.final_results: dict[str, object] = {}

    def evaluate(self, coordinates) -> tuple[np.ndarray, dict | None, dict, dict]:
        vector = np.clip(np.asarray(coordinates, float), 0.0, 1.0)
        key = tuple(np.round(vector, 12))
        if key in self.cache:
            return self.cache[key]
        values = common.normalized_to_values(vector, self.parameters)
        started = time.perf_counter()
        status = "ok"
        runs = None
        try:
            runs = common.simulate_direction(self.context, values, True)
            exact, feasibility = common.charge_residual(self.context, runs)
            if exact is None:
                residual, _ = common.charge_residual(
                    self.context, runs, invalid_penalty_mV=1000.0
                )
                status = "invalid_early_cutoff"
            else:
                residual = exact
        except Exception as exc:
            residual = np.full(3 * len(self.context.curves) * common.N_POINTS, 2000.0)
            feasibility = {
                "feasible": False,
                "minimum_time_margin_min": float("nan"),
                "time_margins_min": {},
            }
            status = f"failed:{type(exc).__name__}"
        elapsed = time.perf_counter() - started
        rmse = float(np.sqrt(np.mean(residual**2)))
        row = {
            "Evaluation": len(self.history) + 1,
            "Objective_RMSE_mV": rmse,
            "Status": status,
            "Feasible": bool(feasibility["feasible"]),
            "Minimum_time_margin_min": feasibility["minimum_time_margin_min"],
            "Elapsed_s": elapsed,
            **values,
        }
        self.history.append(row)
        if len(self.history) == 1 or len(self.history) % 10 == 0:
            print(
                f"{self.model_name} evaluation {len(self.history)}: "
                f"{rmse:.3f} mV ({status})",
                flush=True,
            )
            self.save_history()
        result = (residual, runs, feasibility, values)
        self.cache[key] = result
        return result

    def save_history(self):
        pd.DataFrame(self.history).to_csv(
            self.result_dir / "optimizer_evaluation_history.csv",
            index=False,
            encoding="utf-8-sig",
        )

    def residual(self, coordinates):
        return self.evaluate(coordinates)[0]

    def scalar(self, coordinates):
        residual = self.residual(coordinates)
        return float(np.mean(residual**2))

    def record(self, label: str, vector, result, algorithm: str, started: float):
        coordinates = np.clip(np.asarray(vector, float), 0.0, 1.0)
        residual, _, feasibility, values = self.evaluate(coordinates)
        self.final_vectors[label] = coordinates
        self.final_results[label] = result
        self.final_rows.append(
            {
                "Candidate": label,
                "Algorithm": algorithm,
                "Success": bool(getattr(result, "success", False)),
                "Message": str(getattr(result, "message", "")),
                "NFEV_reported": int(getattr(result, "nfev", -1)),
                "Wall_s": time.perf_counter() - started,
                "Charge_objective_RMSE_mV": float(np.sqrt(np.mean(residual**2))),
                "Feasible": bool(feasibility["feasible"]),
                "Minimum_time_margin_min": feasibility["minimum_time_margin_min"],
                **values,
            }
        )
        pd.DataFrame(self.final_rows).to_csv(
            self.result_dir / "optimizer_candidates.csv",
            index=False,
            encoding="utf-8-sig",
        )

    def local_fit(self, start, max_nfev):
        return least_squares(
            self.residual,
            np.asarray(start, float),
            bounds=(np.zeros(len(self.parameters)), np.ones(len(self.parameters))),
            method="trf",
            x_scale="jac",
            diff_step=1.0e-3,
            max_nfev=max_nfev,
            ftol=1.0e-4,
            xtol=1.0e-4,
            gtol=1.0e-4,
        )

    def run(self):
        max_local = 80 if self.model_name == "main" else 60
        for index, (start_name, start) in enumerate(make_starts(self.parameters), start=1):
            print(f"{self.model_name}: TRF start {index}/5 ({start_name})", flush=True)
            begun = time.perf_counter()
            result = self.local_fit(start, max_local)
            self.record(
                f"TRF {index} - {start_name}", result.x, result, "Multi-start TRF", begun
            )

        bounds = [(0.0, 1.0)] * len(self.parameters)
        print(f"{self.model_name}: differential evolution", flush=True)
        begun = time.perf_counter()
        global_result = differential_evolution(
            self.scalar,
            bounds,
            popsize=6,
            maxiter=8,
            seed=260927,
            polish=False,
            workers=1,
            updating="immediate",
            tol=0.01,
        )
        polished = self.local_fit(global_result.x, max_local)
        self.record("DE -> TRF", polished.x, polished, "DE + TRF", begun)

        print(f"{self.model_name}: dual annealing", flush=True)
        begun = time.perf_counter()
        global_result = dual_annealing(
            self.scalar,
            bounds,
            maxfun=180 if self.model_name == "main" else 120,
            seed=260927,
            no_local_search=True,
            x0=np.full(len(self.parameters), 0.5),
        )
        polished = self.local_fit(global_result.x, max_local)
        self.record("DA -> TRF", polished.x, polished, "DA + TRF", begun)
        self.save_history()
        return self.postprocess()

    def postprocess(self):
        candidates = pd.DataFrame(self.final_rows)
        feasible = candidates[candidates.Feasible.astype(bool)].sort_values(
            "Charge_objective_RMSE_mV"
        )
        if feasible.empty:
            raise RuntimeError(f"No feasible optimizer result for {self.model_name}")
        best_name = str(feasible.iloc[0].Candidate)
        best_values = common.normalized_to_values(
            self.final_vectors[best_name], self.parameters
        )
        charge_runs = common.simulate_direction(self.context, best_values, True)
        discharge_runs = common.simulate_direction(self.context, best_values, False)
        detail = common.postfit_detail(
            self.context, best_name, charge_runs, discharge_runs
        )
        detail.to_csv(
            self.result_dir / "best_postfit_detail_by_cell.csv",
            index=False,
            encoding="utf-8-sig",
        )
        summary = (
            detail.groupby(["Direction", "C_rate"], as_index=False)
            .agg(
                Full_RMSE_mV=("Full_RMSE_mV", "mean"),
                Center10_70_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
                Capacity_error_mAh=("Capacity_error_mAh", "mean"),
                Capacity_MAE_mAh=("Capacity_error_mAh", lambda x: float(np.mean(np.abs(x)))),
                Capacity_RMSE_pct=(
                    "Capacity_error_pct",
                    lambda x: float(np.sqrt(np.mean(np.asarray(x, float) ** 2))),
                ),
            )
        )
        summary.to_csv(
            self.result_dir / "best_postfit_summary.csv",
            index=False,
            encoding="utf-8-sig",
        )

        fig, axes = plt.subplots(2, 3, figsize=(16, 9), constrained_layout=True)
        for row, (is_charge, runs) in enumerate(((True, charge_runs), (False, discharge_runs))):
            direction = "Charge" if is_charge else "Discharge"
            for column, rate in enumerate(ga.RATES):
                ax = axes[row, column]
                for idx, (cell, branches) in enumerate(self.context.curves.items()):
                    observed = branches[(rate, is_charge)]
                    ax.plot(
                        observed["t_min"],
                        observed["V"],
                        color="black",
                        alpha=0.30,
                        lw=1.0,
                        label="Experiment, n=3" if idx == 0 else None,
                    )
                ax.plot(runs[rate]["t_min"], runs[rate]["V"], color="#0072B2", lw=2.2, label="DFN")
                ax.set(title=f"{rate:g}C {direction}", xlabel="Time [min]", ylabel="Voltage [V]")
                ax.grid(alpha=0.2)
                ax.legend(fontsize=8)
        fig.suptitle(f"{self.model_name}: charge fit and untouched discharge validation")
        fig.savefig(self.result_dir / "best_charge_fit_discharge_validation.png", dpi=220)
        plt.close(fig)

        # Recompute a stable physical-time Jacobian at the selected optimum.
        center = self.final_vectors[best_name]
        center_residual, _, center_feasibility, _ = self.evaluate(center)
        if not center_feasibility["feasible"]:
            raise RuntimeError("Selected optimum is not full-grid feasible")
        columns = []
        fd_step = 0.01
        for index, parameter in enumerate(self.parameters):
            low = center.copy()
            high = center.copy()
            low[index] = max(0.0, center[index] - fd_step)
            high[index] = min(1.0, center[index] + fd_step)
            low_residual, _, low_feasibility, _ = self.evaluate(low)
            high_residual, _, high_feasibility, _ = self.evaluate(high)
            if low_feasibility["feasible"] and high_feasibility["feasible"]:
                derivative = (high_residual - low_residual) / (high[index] - low[index])
            elif high_feasibility["feasible"]:
                derivative = (high_residual - center_residual) / (high[index] - center[index])
            elif low_feasibility["feasible"]:
                derivative = (center_residual - low_residual) / (center[index] - low[index])
            else:
                raise RuntimeError(f"No feasible finite-difference side for {parameter}")
            columns.append(derivative)
        jacobian = np.column_stack(columns)
        unit = jacobian / np.maximum(np.linalg.norm(jacobian, axis=0), np.finfo(float).eps)
        cosine = unit.T @ unit
        singular = np.linalg.svd(unit, compute_uv=False)
        condition = float(singular[0] / singular[-1])
        pd.DataFrame(cosine, index=self.parameters, columns=self.parameters).to_csv(
            self.result_dir / "best_postfit_sensitivity_cosine.csv", encoding="utf-8-sig"
        )
        np.save(self.result_dir / "best_postfit_jacobian.npy", jacobian)

        result = {
            "Model": self.model_name,
            "Best_candidate": best_name,
            "Charge_objective_RMSE_mV": float(feasible.iloc[0].Charge_objective_RMSE_mV),
            "Minimum_time_margin_min": float(feasible.iloc[0].Minimum_time_margin_min),
            **best_values,
            "Postfit_normalized_condition_number": condition,
            "Postfit_max_abs_pairwise_cosine": float(
                np.max(np.abs(cosine - np.eye(len(self.parameters))))
            ),
            "Charge_full_RMSE_mV": float(summary[summary.Direction == "Charge"].Full_RMSE_mV.mean()),
            "Discharge_full_RMSE_mV": float(summary[summary.Direction == "Discharge"].Full_RMSE_mV.mean()),
            "Charge_capacity_RMSE_pct": float(
                np.sqrt(
                    np.mean(
                        detail[detail.Direction == "Charge"].Capacity_error_pct.to_numpy(float)
                        ** 2
                    )
                )
            ),
            "Discharge_capacity_RMSE_pct": float(
                np.sqrt(
                    np.mean(
                        detail[detail.Direction == "Discharge"].Capacity_error_pct.to_numpy(float)
                        ** 2
                    )
                )
            ),
        }
        (self.result_dir / "best_result.json").write_text(
            json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
        )
        return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--model", choices=("main", "conservative", "all", "summarize"), default="all"
    )
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    lsa = pd.read_csv(LSA)
    passed = (
        lsa.groupby("Model").Pass.apply(lambda x: bool(np.all(x.astype(bool)))).to_dict()
    )
    requested = (
        list(common.MODEL_PARAMETERS)
        if args.model in {"all", "summarize"}
        else [args.model]
    )
    context = None if args.model == "summarize" else common.build_context()
    comparison = []
    for model_name in requested:
        existing = OUT / model_name / "best_result.json"
        if args.model == "summarize" and existing.exists():
            comparison.append(json.loads(existing.read_text(encoding="utf-8")))
            continue
        if not passed.get(model_name, False):
            print(f"skip {model_name}: physical-time identifiability did not pass", flush=True)
            continue
        comparison.append(ModelFitter(context, model_name).run())
        pd.DataFrame(comparison).to_csv(
            OUT / "model_comparison.csv", index=False, encoding="utf-8-sig"
        )
    pd.DataFrame(comparison).to_csv(
        OUT / "model_comparison.csv", index=False, encoding="utf-8-sig"
    )
    manifest = {
        "objective": "charge voltage ordinary L2/RMSE only, 100 physical-time points per cell and C-rate",
        "capacity_in_objective": False,
        "invalid_candidate": "simulation cutoff before any fixed experimental grid endpoint; penalized only as a feasibility violation",
        "validation": "July BoL discharge 0.5C/1C/2C never used by optimizers",
        "current_mass_corrected_GITT": {
            "Dsn": common.CURRENT_MASS_DSN,
            "Dsp": common.CURRENT_MASS_DSP,
        },
        "bounds": common.BOUNDS,
        "algorithms": {
            "multi_start_TRF": "5 starts; main max_nfev=80, conservative=60",
            "differential_evolution": "popsize=6, maxiter=8, seed=260927, then TRF",
            "dual_annealing": "maxfun main=180/conservative=120, seed=260927, then TRF",
        },
    }
    (OUT / "analysis_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print("\nMODEL COMPARISON\n", pd.DataFrame(comparison).to_string(index=False))


if __name__ == "__main__":
    main()
