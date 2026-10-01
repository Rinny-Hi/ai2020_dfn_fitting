"""Shared, final-design charge-fitting setup on fixed physical-time grids."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np

import c50_selected_ocp_dynamic_validation as selected
import comprehensive_prefit_cross_cohort as cross
import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import paper_geometry_area_porosity_ablation as area_study
import p0_five_parameter_sensitivity_area as legacy
import prefit_physical_anchor_decision as p0
import priority_initial_state_protocol_recheck as priority


N_POINTS = 100
CURRENT_MASS_DSN = 4.170774145968247e-14
CURRENT_MASS_DSP = 6.888175323807748e-14

INITIAL = {
    "Dsn": CURRENT_MASS_DSN,
    "kn": 9.648533212e-7,
    "brugg_n": 2.914,
}
BOUNDS = {
    "Dsn": (1.0e-14, 8.0e-14),
    "kn": (3.0e-7, 3.0e-6),
    "brugg_n": (1.5, 3.5),
}
FIXED = {"Dsp": CURRENT_MASS_DSP, "brugg_p": 1.83, "brugg_s": 1.5}

MODEL_PARAMETERS = {
    "main": ("Dsn", "kn", "brugg_n"),
    "conservative": ("Dsn", "kn"),
}


@dataclass
class Context:
    curves: dict
    qc: object
    protocol: object
    model_geometry: dict
    stage: dict
    anode: object
    cathode: object
    qcell: float
    kp: float
    grids: dict
    experimental_blocks: dict


def _cc_view(branch: dict) -> dict[str, np.ndarray]:
    """Start at the first point at or above 95% of the median CC current."""
    time = np.asarray(branch["t_min"], float)
    voltage = np.asarray(branch["V"], float)
    current = np.abs(np.asarray(branch["I_A"], float))
    capacity = np.asarray(branch["Q_Ah"], float)
    command = float(np.nanmedian(current))
    matches = np.flatnonzero(current >= 0.95 * command)
    if not len(matches):
        raise RuntimeError("No sample reaches 95% of the CC command")
    first = int(matches[0])
    return {
        "t_min": time[first:] - time[first],
        "V": voltage[first:],
        "I_A": current[first:],
        "Q_Ah": capacity[first:] - capacity[first],
    }


def build_context(n_points: int = N_POINTS) -> Context:
    curves, qc, protocol = legacy.load_july()
    model0, stage, anode, cathode, _, _, qcell = selected.selected_inputs()
    anchored = p0.build_model(model0, p0.Scenario("P0 anchored reference"))
    geometry = area_study.apply_paper_overlap(anchored, preserve_capacity=True)
    kp = float(legacy.kp_from_area_specific_rct(geometry)["kp"])
    grids = {}
    experimental = {}
    for rate in ga.RATES:
        for cell, branches in curves.items():
            view = _cc_view(branches[(rate, True)])
            grid = np.linspace(0.0, float(view["t_min"][-1]), n_points)
            grids[(rate, cell)] = grid
            experimental[(rate, cell)] = np.interp(grid, view["t_min"], view["V"])
    return Context(
        curves=curves,
        qc=qc,
        protocol=protocol,
        model_geometry=geometry,
        stage=stage,
        anode=anode,
        cathode=cathode,
        qcell=float(qcell),
        kp=kp,
        grids=grids,
        experimental_blocks=experimental,
    )


def complete_values(values: dict[str, float], context: Context) -> dict[str, float]:
    return {
        "Dsn": float(values.get("Dsn", INITIAL["Dsn"])),
        "Dsp": FIXED["Dsp"],
        "kn": float(values.get("kn", INITIAL["kn"])),
        "kp": context.kp,
        "brugg_n": float(values.get("brugg_n", INITIAL["brugg_n"])),
    }


def simulate_direction(context: Context, values: dict[str, float], charge: bool) -> dict:
    model = legacy.apply_parameters(context.model_geometry, complete_values(values, context))
    runs = {}
    direction = "Charge" if charge else "Discharge"
    for rate in ga.RATES:
        row = context.protocol.loc[("Old July BoL", rate, direction)]
        z0 = priority.monotone_voltage_inverse(
            model["v_qocv"], model["soc"], float(row.Rest_end_V)
        )
        local_stage = priority.stage_at_soc(context.stage, z0, charge)
        sim_rate = float(row.Measured_current_A) / cross.NOMINAL_CAPACITY_AH
        runs[rate] = ehq.run_dfn(
            sim_rate,
            charge,
            local_stage,
            model,
            context.anode,
            context.cathode,
            False,
            False,
        )
    return runs


def charge_residual(
    context: Context, runs: dict, invalid_penalty_mV: float | None = None
) -> tuple[np.ndarray | None, dict]:
    blocks = []
    margins = {}
    missing_any = False
    for rate in ga.RATES:
        sim_t = np.asarray(runs[rate]["t_min"], float)
        sim_v = np.asarray(runs[rate]["V"], float)
        for cell in context.curves:
            grid = context.grids[(rate, cell)]
            margin = float(sim_t[-1] - grid[-1])
            margins[(rate, cell)] = margin
            if margin < -1.0e-9:
                missing_any = True
            prediction = np.interp(np.minimum(grid, sim_t[-1]), sim_t, sim_v)
            residual = 1000.0 * (prediction - context.experimental_blocks[(rate, cell)])
            if margin < -1.0e-9 and invalid_penalty_mV is not None:
                missing = np.maximum(grid - sim_t[-1], 0.0) / max(grid[-1], 1.0e-12)
                residual = residual + invalid_penalty_mV * (1.0 + missing)
            blocks.append(residual)
    vector = np.concatenate(blocks)
    if missing_any and invalid_penalty_mV is None:
        vector = None
    return vector, {
        "feasible": not missing_any,
        "minimum_time_margin_min": min(margins.values()),
        "time_margins_min": {
            f"{rate:g}C_{cell}": value for (rate, cell), value in margins.items()
        },
    }


def normalized_to_values(
    coordinates: Iterable[float], parameters: tuple[str, ...]
) -> dict[str, float]:
    result = dict(INITIAL)
    for name, coordinate in zip(parameters, coordinates):
        lo, hi = BOUNDS[name]
        z = float(coordinate)
        if name in {"Dsn", "kn"}:
            result[name] = float(np.exp(np.log(lo) + z * (np.log(hi) - np.log(lo))))
        else:
            result[name] = float(lo + z * (hi - lo))
    return result


def values_to_normalized(values: dict[str, float], parameters: tuple[str, ...]) -> np.ndarray:
    coordinates = []
    for name in parameters:
        lo, hi = BOUNDS[name]
        value = float(values[name])
        if name in {"Dsn", "kn"}:
            coordinates.append((np.log(value) - np.log(lo)) / (np.log(hi) - np.log(lo)))
        else:
            coordinates.append((value - lo) / (hi - lo))
    return np.asarray(coordinates, float)


def postfit_detail(context: Context, name: str, charge_runs: dict, discharge_runs: dict):
    return legacy.summarize_runs(
        name,
        {True: charge_runs, False: discharge_runs},
        context.curves,
        context.qcell,
    ).rename(columns={"Area_policy": "Candidate"})
