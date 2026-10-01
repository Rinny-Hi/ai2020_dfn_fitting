"""Finalize the completed September profile run with feasible-side Jacobian."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
WORK = PROJECT / "work"
OUT = PROJECT / "outputs" / "260927_september_recommended_dsn_kn_fit"
sys.path.insert(0, str(WORK))
import september_kn_only_after_identifiability as fit
import september_recommended_dsn_kn_fit as base


def main() -> None:
    context = base.build_context()
    dsn_audit = pd.read_csv(OUT / "dsn_kn_prefit_identifiability.csv")
    candidates = pd.read_csv(OUT / "optimizer_candidates.csv")
    profile = pd.read_csv(OUT / "dsn_fixed_profile_summary.csv")
    best = candidates[candidates.Feasible.astype(bool)].sort_values("Objective_RMSE_mV").iloc[0]
    values = {"Dsn": fit.GITT_DSN, "kn": float(best.kn)}

    evaluator = fit.KnEvaluator(context, fit.GITT_DSN)
    center = fit.kn_to_z(values["kn"])
    ident_result, ident = fit.kn_identifiability(evaluator, center, "post-fit kn-only")
    ident.to_csv(OUT / "kn_only_identifiability.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame([ident_result]).to_csv(
        OUT / "kn_only_identifiability_summary.csv", index=False, encoding="utf-8-sig"
    )
    if not ident_result["Pass"]:
        raise RuntimeError(f"Feasible-side kn identifiability failed: {ident_result}")

    charge_runs = base.simulate_direction(context, values, "Charge")
    discharge_runs = base.simulate_direction(context, values, "Discharge")
    metrics = pd.concat(
        [
            base.branch_metrics(context, charge_runs, "Charge"),
            base.branch_metrics(context, discharge_runs, "Discharge"),
        ],
        ignore_index=True,
    )
    metrics.to_csv(OUT / "branch_metrics.csv", index=False, encoding="utf-8-sig")
    by_rate = (
        metrics.groupby(["Direction", "C_rate"], as_index=False)
        .agg(
            N_cells=("Cell", "count"),
            All_physical_time_feasible=("Physical_time_feasible", "all"),
            Mean_physical_time_RMSE_mV=("Physical_time_RMSE_mV", "mean"),
            Mean_physical_time_bias_mV=("Physical_time_bias_mV", "mean"),
            Capacity_RMSE_pct=(
                "Capacity_error_pct",
                lambda x: float((x.astype(float).pow(2).mean()) ** 0.5),
            ),
            Capacity_MAE_mAh=("Capacity_error_mAh", lambda x: float(x.astype(float).abs().mean())),
        )
    )
    by_rate.to_csv(OUT / "summary_by_rate.csv", index=False, encoding="utf-8-sig")
    base.plot_fit(context, charge_runs, discharge_runs)
    summary = fit.validation_summary(metrics)
    profile.loc[
        profile.Fixed_Dsn.sub(fit.GITT_DSN).abs().idxmin(),
        list(summary.keys()),
    ] = list(summary.values())
    profile.loc[
        ~profile.Discharge_all_time_feasible.astype(bool),
        "Discharge_physical_time_RMSE_mV",
    ] = float("nan")
    profile.to_csv(OUT / "dsn_fixed_profile_summary.csv", index=False, encoding="utf-8-sig")
    fit.profile_plot(profile)

    manifest = {
        "decision": "Dsn+kn rejected before fitting because Dsn relative sensitivity stayed below 0.10; fit kn only with Dsn fixed to independent mass-corrected GITT",
        "selected_ocp": "Ai2020 nominal graphite equilibrium + old measured cathode equilibrium",
        "window_policy": "retain previously identified window; update each September branch initial SOC from its own 2-hour rest endpoint voltage",
        "stage": {key: float(context.stage[key]) for key in ("x0", "x100", "y100", "y0")},
        "qcell_Ah": context.qcell,
        "fixed": {**base.FIXED, "Dsn": fit.GITT_DSN, "kp": context.kp},
        "fitted": {"kn": values["kn"]},
        "objective": "charge CC ordinary voltage RMSE; 2 cells x 3 rates x 100 fixed physical-time points",
        "capacity_in_objective": False,
        "cv_in_objective": False,
        "validation": summary,
        "kn_identifiability": ident_result,
    }
    (OUT / "analysis_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    fit.write_report(
        context,
        dsn_audit,
        ident,
        candidates,
        profile,
        values,
        summary,
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2), flush=True)
    print(profile.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
