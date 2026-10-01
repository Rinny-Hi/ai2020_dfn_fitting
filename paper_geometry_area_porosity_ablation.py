"""Ablate paper collector, overlap-area and source-porosity choices.

All dynamic cases use the selected capacity-consistent OCP, measured branch
currents and rest-voltage initial SOC.  The raw-area case is diagnostic only:
it changes geometric inventory without re-identifying the OCP window.  The
capacity-preserved case rescales c_s,max inversely with area so Qn/Qp and the
OCP window remain consistent.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

import capacity_consistent_qn_qp_kn_refinement as refinement
import capacity_consistent_refined_configuration as refined
import paper_teardown_geometry as paper
import priority_initial_state_protocol_recheck as priority


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260921_paper_geometry_area_porosity_ablation"
RESULTS.mkdir(parents=True, exist_ok=True)


def copy_model(model):
    changed = dict(model)
    changed["params"] = model["params"].copy()
    return changed


def apply_paper_collectors(model):
    changed = copy_model(model)
    changed["params"].update(
        {
            "Negative current collector thickness [m]": paper.NEGATIVE_CURRENT_COLLECTOR_M,
            "Positive current collector thickness [m]": paper.POSITIVE_CURRENT_COLLECTOR_M,
        },
        check_already_exists=False,
    )
    return changed


def area(model):
    p = model["params"]
    return (
        float(p["Electrode height [m]"])
        * float(p["Electrode width [m]"])
        * float(p["Number of electrodes connected in parallel to make a cell"])
    )


def apply_paper_overlap(model, preserve_capacity):
    old_area = area(model)
    changed = paper.apply_overlap_area(model)
    new_area = area(changed)
    if preserve_capacity:
        scale = old_area / new_area
        changed["csn_max"] = float(model["csn_max"]) * scale
        changed["csp_max"] = float(model["csp_max"]) * scale
        changed["params"].update(
            {
                "Maximum concentration in negative electrode [mol.m-3]": changed["csn_max"],
                "Maximum concentration in positive electrode [mol.m-3]": changed["csp_max"],
            },
            check_already_exists=False,
        )
    return changed


def apply_source_porosity(model):
    changed = copy_model(model)
    changed["params"].update(
        {
            "Negative electrode porosity": paper.NEGATIVE_POROSITY,
            "Positive electrode porosity": paper.POSITIVE_POROSITY,
        },
        check_already_exists=False,
    )
    return changed


def main():
    model, stage, anode, cathode, experiment, report, qcell = refined.refined_inputs()
    protocol, _, _ = priority.load_protocol(model)
    baseline = apply_paper_collectors(model)
    raw_area = apply_paper_overlap(baseline, preserve_capacity=False)
    preserved = apply_paper_overlap(baseline, preserve_capacity=True)
    source_porosity = apply_source_porosity(preserved)

    models = {
        "Current area + paper collectors": baseline,
        "Paper overlap raw (diagnostic)": raw_area,
        "Paper overlap + Q preserved": preserved,
        "Paper overlap + Q preserved + source porosity": source_porosity,
    }
    detail_frames = []
    for name, local_model in models.items():
        detail, _ = refinement.run_candidate(
            name, local_model, stage, anode, cathode, experiment, protocol, qcell
        )
        detail_frames.append(detail)
    detail = pd.concat(detail_frames, ignore_index=True)
    summary = refinement.summarize(detail)
    detail.to_csv(RESULTS / "dynamic_detail.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(RESULTS / "dynamic_summary.csv", index=False, encoding="utf-8-sig")

    old_area = area(baseline)
    new_area = area(preserved)
    metadata = {
        "old_common_area_cm2": old_area * 1e4,
        "paper_negative_total_area_cm2": paper.NEGATIVE_TOTAL_AREA_M2 * 1e4,
        "paper_positive_total_area_cm2": paper.POSITIVE_TOTAL_AREA_M2 * 1e4,
        "paper_dfn_overlap_area_cm2": new_area * 1e4,
        "area_ratio_new_over_old": new_area / old_area,
        "current_density_ratio_new_over_old": old_area / new_area,
        "raw_geometric_capacity_change_pct": 100.0 * (new_area / old_area - 1.0),
        "capacity_preserving_csmax_scale": old_area / new_area,
        "paper_collectors_um": {"Cu": 10.0, "Al": 15.0},
        "collector_only_effect": (
            "none in the isothermal electrochemical DFN because Ai2020 already uses 10/15 um"
        ),
        "source_porosity": {"negative": 0.32, "positive": 0.33},
        "active_fraction_status": (
            "kept at n=0.61,p=0.62; source values are estimates from measured porosity "
            "and assumed inactive fractions"
        ),
    }
    (RESULTS / "report.json").write_text(
        json.dumps(
            {"metadata": metadata, "summary": summary.to_dict(orient="records")},
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )

    metrics = ["mean_center_RMSE_mV", "capacity_RMSE_pct"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    for ax, metric, ylabel in zip(
        axes, metrics, ["Mean 10-70% RMSE [mV]", "Capacity RMSE [%]"]
    ):
        pivot = summary.pivot(index="Candidate", columns="Direction", values=metric)
        pivot.plot(kind="bar", ax=ax, color=["#3b82f6", "#f97316"])
        ax.set_ylabel(ylabel)
        ax.set_xlabel("")
        ax.grid(axis="y", alpha=0.25)
        ax.tick_params(axis="x", rotation=20)
    fig.suptitle("SPB655060 paper geometry ablation")
    fig.savefig(RESULTS / "paper_geometry_ablation.png", dpi=220)
    plt.close(fig)
    print(summary.to_string(index=False))
    print(json.dumps(metadata, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
