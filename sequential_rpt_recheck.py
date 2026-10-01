"""Simulate the first BoL RPT as one continuous experiment.

This avoids estimating SOC from a still-relaxing 10-minute rest voltage and
lets concentration gradients and hysteresis carry naturally from one branch
to the next.
"""

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
from scipy.integrate import cumulative_trapezoid

import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import priority_initial_state_protocol_recheck as priority
import stage2_seven_parameter_sensitivity as stage2


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260921_priority_initial_state_protocol_recheck"


def run_sequence(model_inputs, stage, anode, cathode):
    model = pybamm.lithium_ion.DFN(
        {
            "open-circuit potential": (
                "one-state hysteresis",
                "one-state hysteresis",
            ),
            "thermal": "isothermal",
        }
    )
    steps = (
        "Discharge at 0.05C until 3.0 V",
        "Rest for 10 minutes",
        "Charge at 0.5C until 4.2 V",
        "Hold at 4.2 V until C/20",
        "Rest for 10 minutes",
        "Discharge at 0.5C until 3.0 V",
        "Rest for 10 minutes",
        "Charge at 1C until 4.2 V",
        "Hold at 4.2 V until C/20",
        "Rest for 10 minutes",
        "Discharge at 1C until 3.0 V",
        "Rest for 10 minutes",
        "Charge at 2C until 4.2 V",
        "Hold at 4.2 V until C/20",
        "Rest for 10 minutes",
        "Discharge at 2C until 3.0 V",
    )
    values = ehq.make_parameter_values(
        model_inputs,
        stage,
        anode,
        cathode,
        charge=False,
        negative_hysteresis=True,
        positive_hysteresis=True,
        # The sequence starts from a discharge branch at the charged endpoint.
        negative_initial_h=1.0,
        positive_initial_h=-1.0,
    )
    simulation = pybamm.Simulation(
        model,
        parameter_values=values,
        experiment=pybamm.Experiment([steps], period="10 seconds"),
        var_pts={"x_n": 15, "x_s": 15, "x_p": 15, "r_n": 15, "r_p": 15},
        solver=pybamm.IDAKLUSolver(rtol=1e-6, atol=1e-8),
    )
    solution = simulation.solve()
    cycle = solution.cycles[0]
    if len(cycle.steps) != len(steps):
        raise RuntimeError(
            f"Expected {len(steps)} sequential steps, obtained {len(cycle.steps)}"
        )
    branch_positions = {
        (0.5, True): 2,
        (0.5, False): 5,
        (1.0, True): 7,
        (1.0, False): 10,
        (2.0, True): 12,
        (2.0, False): 15,
    }
    runs = {}
    q_meas = float(model_inputs["q_meas"])
    for key, position in branch_positions.items():
        rate, charge = key
        step_solution = cycle.steps[position]
        time_s = np.asarray(step_solution["Time [s]"].entries, dtype=float)
        time_s = time_s - time_s[0]
        voltage = np.asarray(step_solution["Terminal voltage [V]"].entries, dtype=float)
        current = np.asarray(step_solution["Current [A]"].entries, dtype=float)
        capacity = cumulative_trapezoid(np.abs(current), time_s, initial=0.0) / 3600.0
        soc = capacity / q_meas if charge else 1.0 - capacity / q_meas
        runs[key] = {
            "t_min": time_s / 60.0,
            "V": voltage,
            "SOC": np.clip(soc, 0.0, 1.0),
        }
    return runs


def main():
    warnings.filterwarnings(
        "ignore", message="The definition of the hysteresis decay rate parameter has changed"
    )
    _, model_inputs, stage, anode, cathode, experiment, _ = stage2.build_current_inputs()
    print("Running the complete first RPT as one continuous DFN experiment", flush=True)
    runs = run_sequence(model_inputs, stage, anode, cathode)
    detail, summary = priority.score_scenarios(
        {"F Continuous original RPT sequence": runs},
        experiment,
        float(model_inputs["q_meas"]),
    )
    detail.to_csv(
        RESULTS / "sequential_rpt_condition_metrics.csv", index=False, encoding="utf-8-sig"
    )
    summary.to_csv(
        RESULTS / "sequential_rpt_summary.csv", index=False, encoding="utf-8-sig"
    )

    fig, axes = plt.subplots(2, 3, figsize=(16, 8.5), constrained_layout=True)
    for col, rate in enumerate(ga.RATES):
        for row, charge in enumerate((True, False)):
            ax = axes[row, col]
            exp = experiment[(rate, charge)]
            sim = runs[(rate, charge)]
            metric = detail[
                np.isclose(detail.C_rate, rate)
                & (detail.Direction == ("Charge" if charge else "Discharge"))
            ].iloc[0]
            ax.plot(
                priority.transferred_capacity(exp, model_inputs["q_meas"], charge),
                exp["V"],
                color="black",
                lw=2.5,
                label="Experiment",
            )
            ax.plot(
                priority.transferred_capacity(sim, model_inputs["q_meas"], charge),
                sim["V"],
                color="#0072B2",
                lw=2.0,
                label=f"Continuous DFN: {metric.Full_RMSE_mV:.1f} mV",
            )
            ax.set_title(f"{rate:g}C {'Charge' if charge else 'Discharge'}")
            ax.set_xlabel("Transferred capacity (Ah)")
            ax.set_ylabel("Voltage (V)")
            ax.grid(alpha=0.22)
    axes[0, 0].legend()
    fig.suptitle("Sequential reproduction of the original Enertech BoL RPT")
    fig.savefig(RESULTS / "sequential_rpt_curves.png", dpi=220)
    plt.close(fig)

    report = {
        "method": "one continuous PyBaMM experiment with CC, CV, and 10-minute rests in original order",
        "initial_state": "charged endpoint; discharge hysteresis state",
        "kinetics_refitted": False,
        "summary": summary.to_dict(orient="records"),
        "conditions": detail.to_dict(orient="records"),
    }
    (RESULTS / "sequential_rpt_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(summary.to_string(index=False))
    print(detail.to_string(index=False))


if __name__ == "__main__":
    main()
