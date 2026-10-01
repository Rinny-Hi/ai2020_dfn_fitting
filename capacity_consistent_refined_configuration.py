"""Reusable provisional configuration from the capacity-consistent refinement."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import c50_ocp_source_combination_study as c50
import c50_selected_ocp_dynamic_validation as selected
import capacity_consistent_qn_qp_kn_refinement as refinement


ROOT = Path(__file__).resolve().parent
REPORT = (
    ROOT
    / "results"
    / "260921_capacity_consistent_qn_qp_kn_refinement"
    / "capacity_consistent_refinement_report.json"
)


def refined_inputs():
    """Return the +2% Qn/Qp provisional model while preserving capacity-window equality."""
    model, _, anode, cathode, experiment, selected_row, qcell = selected.selected_inputs()
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    chosen = report["selected_capacity_candidate"]
    stage = {key: float(report["selected_stage"][key]) for key in ("x0", "x100", "y100", "y0")}
    model = refinement.model_with_capacities(
        model,
        float(chosen["Qn_Ah"]),
        float(chosen["Qp_Ah"]),
        float(selected_row.Qn_Ah),
        float(selected_row.Qp_Ah),
    )
    soc = np.asarray(model["soc"], dtype=float)
    model["v_qocv"] = c50.predict(
        np.asarray([stage["x0"], stage["x100"], stage["y100"], stage["y0"]]),
        soc,
        anode,
        cathode,
        "equilibrium",
    )
    model["delta_x"] = stage["x100"] - stage["x0"]
    model["delta_y"] = stage["y0"] - stage["y100"]
    return model, stage, anode, cathode, experiment, report, qcell

