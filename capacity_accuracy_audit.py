"""Reproducible audit of capacity accounting and capacity/voltage trade-offs.

This script does not identify new material properties. It checks the current
capacity bookkeeping and summarises already-computed feasibility/refinement
results using the measured CC current of every branch.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import c50_selected_ocp_dynamic_validation as selected
import gitt_ocp_analysis as ga


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260922_capacity_accuracy_audit"
RESULTS.mkdir(parents=True, exist_ok=True)
F = 96485.33212
NOMINAL_CAPACITY_AH = 2.28


def parameter_scalar(params, name: str) -> float:
    return float(params.evaluate(params[name]))


def main() -> None:
    model, stage, _, _, experiment, selected_row, qcell = selected.selected_inputs()
    params = model["params"]
    one_face_area = parameter_scalar(params, "Electrode height [m]") * parameter_scalar(
        params, "Electrode width [m]"
    )
    n_parallel = parameter_scalar(
        params, "Number of electrodes connected in parallel to make a cell"
    )
    area = one_face_area * n_parallel
    ln = parameter_scalar(params, "Negative electrode thickness [m]")
    lp = parameter_scalar(params, "Positive electrode thickness [m]")
    eps_n = parameter_scalar(params, "Negative electrode active material volume fraction")
    eps_p = parameter_scalar(params, "Positive electrode active material volume fraction")
    csn = parameter_scalar(params, "Maximum concentration in negative electrode [mol.m-3]")
    csp = parameter_scalar(params, "Maximum concentration in positive electrode [mol.m-3]")
    qn_formula = F * area * ln * eps_n * csn / 3600.0
    qp_formula = F * area * lp * eps_p * csp / 3600.0

    accounting = pd.DataFrame(
        [
            {
                "quantity": "Q_cell_C50_target_Ah",
                "direct_calculation": qcell,
                "model_value": qcell,
                "difference": 0.0,
            },
            {
                "quantity": "Q_n_Ah",
                "direct_calculation": qn_formula,
                "model_value": float(selected_row.Qn_Ah),
                "difference": qn_formula - float(selected_row.Qn_Ah),
            },
            {
                "quantity": "Q_p_Ah",
                "direct_calculation": qp_formula,
                "model_value": float(selected_row.Qp_Ah),
                "difference": qp_formula - float(selected_row.Qp_Ah),
            },
        ]
    )
    accounting.to_csv(RESULTS / "capacity_formula_crosscheck.csv", index=False, encoding="utf-8-sig")

    current_rows = []
    for rate in ga.RATES:
        for charge in (True, False):
            branch = experiment[(rate, charge)]
            current_a = float(np.nanmedian(np.abs(branch["I_A"])))
            setpoint_a = float(rate * NOMINAL_CAPACITY_AH)
            current_rows.append(
                {
                    "C_rate_label": rate,
                    "Direction": "Charge" if charge else "Discharge",
                    "measured_CC_current_A": current_a,
                    "nominal_setpoint_A": setpoint_a,
                    "relative_error_pct": 100.0 * (current_a / setpoint_a - 1.0),
                }
            )
    current_audit = pd.DataFrame(current_rows)
    current_audit.to_csv(RESULTS / "measured_current_crosscheck.csv", index=False, encoding="utf-8-sig")

    refinement_dir = ROOT / "results" / "260921_capacity_consistent_qn_qp_kn_refinement"
    detail = pd.read_csv(refinement_dir / "capacity_consistent_dynamic_detail.csv")
    candidates = pd.read_csv(refinement_dir / "capacity_candidate_selection.csv")
    kn = pd.read_csv(refinement_dir / "kn_candidate_selection.csv")
    chosen_candidates = ["Qn1.0000_Qp1.0000", "Qn1.0400_Qp1.0400"]
    detail_out = detail[detail.Candidate.isin(chosen_candidates)].copy()
    detail_out.to_csv(RESULTS / "baseline_vs_capacity_feasibility_detail.csv", index=False, encoding="utf-8-sig")
    candidates[candidates.Candidate.isin(chosen_candidates)].to_csv(
        RESULTS / "baseline_vs_capacity_feasibility_summary.csv", index=False, encoding="utf-8-sig"
    )
    kn.to_csv(RESULTS / "kn_capacity_voltage_tradeoff.csv", index=False, encoding="utf-8-sig")

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.8), constrained_layout=True)
    rates = np.asarray(ga.RATES, dtype=float)
    width = 0.18
    colors = {"Charge": "#D55E00", "Discharge": "#0072B2"}
    for j, candidate in enumerate(chosen_candidates):
        label = "Baseline" if j == 0 else "Qn,Qp +4% (feasibility)"
        for i, direction in enumerate(("Charge", "Discharge")):
            rows = detail_out[(detail_out.Candidate == candidate) & (detail_out.Direction == direction)]
            axes[0].bar(
                rates + (j * 2 + i - 1.5) * width,
                rows.Capacity_error_mAh.to_numpy(),
                width,
                color=colors[direction],
                alpha=1.0 if j == 0 else 0.48,
                hatch=None if j == 0 else "//",
                label=f"{label}, {direction}",
            )
    axes[0].axhline(0.0, color="black", lw=0.8)
    axes[0].set(
        xticks=rates,
        xticklabels=[f"{x:g}C" for x in rates],
        ylabel="Model - experiment capacity [mAh]",
        title="Endpoint capacity residual",
    )
    axes[0].legend(fontsize=7, ncol=2)

    selected_candidates = candidates[candidates.Candidate.isin(chosen_candidates)].copy()
    labels = ["Baseline", "Qn,Qp +4%"]
    x = np.arange(2)
    axes[1].bar(x - 0.18, selected_candidates.charge_capacity_RMSE_pct, 0.36, label="Charge capacity", color="#D55E00")
    axes[1].bar(x + 0.18, selected_candidates.discharge_capacity_RMSE_pct, 0.36, label="Discharge capacity", color="#0072B2")
    axes[1].set(xticks=x, xticklabels=labels, ylabel="Capacity RMSE [%]", title="Capacity inventory feasibility")
    ax1b = axes[1].twinx()
    ax1b.plot(x, selected_candidates.qOCV_RMSE_2_98_mV, "ko--", label="qOCV RMSE")
    ax1b.set_ylabel("qOCV RMSE [mV]")
    axes[1].legend(loc="upper right", fontsize=8)
    ax1b.legend(loc="center right", fontsize=8)

    axes[2].plot(kn.kn * 1e7, kn.charge_capacity_RMSE_pct, "o-", color="#D55E00", label="Charge capacity")
    ax2b = axes[2].twinx()
    ax2b.plot(kn.kn * 1e7, kn.charge_mean_center_RMSE_mV, "s--", color="#5E3C99", label="Charge voltage")
    axes[2].set(xlabel=r"$k_n$ [$10^{-7}$ scale]", ylabel="Charge capacity RMSE [%]", title=r"Why $k_n$ must not be fit to capacity")
    ax2b.set_ylabel("10–70% voltage RMSE [mV]")
    axes[2].legend(loc="upper right", fontsize=8)
    ax2b.legend(loc="center right", fontsize=8)

    for ax in axes:
        ax.grid(alpha=0.22)
    fig.savefig(RESULTS / "capacity_error_diagnosis.png", dpi=220)
    plt.close(fig)

    print(accounting.to_string(index=False))
    print(current_audit.to_string(index=False))
    print("Saved to", RESULTS)


if __name__ == "__main__":
    main()
