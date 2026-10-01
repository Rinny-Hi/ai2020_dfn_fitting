from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pybamm
from scipy.integrate import cumulative_trapezoid


REPO = Path(
    r"C:\Users\user\Documents\Codex\2026-09-18\https-github-com-rinny-hi-ai2020"
    r"\work\ai2020_dfn_fitting"
)
OUT = Path(
    r"C:\Users\user\Documents\Codex\2026-09-27"
    r"\github-260921-ocp-eis-current-progress\outputs\260927_ocp_initial_state_deep_dive"
)
sys.path.insert(0, str(REPO))

import charge_physical_time_common as common
import comprehensive_prefit_cross_cohort as cross
import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import p0_five_parameter_sensitivity_area as legacy


FINAL = {"Dsn": 7.466200473885764e-14, "kn": 1.4502357402082455e-6, "brugg_n": 2.914}


def entries(step, name):
    return np.asarray(step[name].entries, dtype=float)


def main():
    context = common.build_context()
    model_inputs = legacy.apply_parameters(
        context.model_geometry, common.complete_values(FINAL, context)
    )
    rates = {
        float(rate): float(
            context.protocol.loc[("Old July BoL", rate, "Charge")].Measured_current_A
        )
        / cross.NOMINAL_CAPACITY_AH
        for rate in ga.RATES
    }
    steps = ["Discharge at C/50 until 3.0 V", "Rest for 10 minutes"]
    charge_step_indices = {}
    for rate in ga.RATES:
        charge_step_indices[float(rate)] = len(steps)
        steps.extend(
            [
                f"Charge at {rates[float(rate)]:.8f}C until 4.2 V",
                "Hold at 4.2 V until C/20",
                "Rest for 10 minutes",
                f"Discharge at {rates[float(rate)]:.8f}C until 3.0 V",
                "Rest for 10 minutes",
            ]
        )
    experiment = pybamm.Experiment([tuple(steps)], period="10 seconds")
    parameters = ehq.make_parameter_values(
        model_inputs,
        context.stage,
        context.anode,
        context.cathode,
        False,
        False,
        False,
    )
    simulation = pybamm.Simulation(
        pybamm.lithium_ion.DFN({"open-circuit potential": ("single", "single"), "thermal": "isothermal"}),
        parameter_values=parameters,
        experiment=experiment,
        var_pts={"x_n": 15, "x_s": 15, "x_p": 15, "r_n": 15, "r_p": 15},
        solver=pybamm.IDAKLUSolver(rtol=1e-6, atol=1e-8),
    )
    solution = simulation.solve()
    cycle = solution.cycles[0]
    step_rows = []
    for index, step_solution in enumerate(cycle.steps):
        step_time = entries(step_solution, "Time [s]")
        step_voltage = entries(step_solution, "Terminal voltage [V]")
        step_current = np.abs(entries(step_solution, "Current [A]"))
        step_capacity = cumulative_trapezoid(
            step_current, step_time - step_time[0], initial=0.0
        ) / 3600.0
        step_rows.append(
            {
                "step_index": index,
                "protocol_step": steps[index],
                "duration_min": float((step_time[-1] - step_time[0]) / 60.0),
                "start_V": float(step_voltage[0]),
                "end_V": float(step_voltage[-1]),
                "absolute_throughput_Ah": float(step_capacity[-1]),
            }
        )
    rows = []
    curve_rows = []
    component_rows = []
    component_names = [
        "Terminal voltage [V]",
        "Battery open-circuit voltage [V]",
        "Battery particle concentration overpotential [V]",
        "X-averaged battery reaction overpotential [V]",
        "X-averaged battery negative reaction overpotential [V]",
        "X-averaged battery positive reaction overpotential [V]",
        "X-averaged battery solid phase ohmic losses [V]",
        "X-averaged battery electrolyte ohmic losses [V]",
        "X-averaged battery concentration overpotential [V]",
    ]
    for rate, step_index in charge_step_indices.items():
        step = cycle.steps[step_index]
        previous_rest = cycle.steps[step_index - 1]
        t = entries(step, "Time [s]")
        t = (t - t[0]) / 60.0
        v = entries(step, "Terminal voltage [V]")
        current = np.abs(entries(step, "Current [A]"))
        q = cumulative_trapezoid(current, t * 60.0, initial=0.0) / 3600.0
        component_record = {"C_rate": rate}
        for name in component_names:
            label = (
                name.replace("X-averaged battery ", "")
                .replace("Battery ", "")
                .replace(" [V]", "")
                .replace(" ", "_")
            )
            rest_values = entries(previous_rest, name)
            charge_values = entries(step, name)
            component_record[f"rest_end_{label}_mV"] = 1000.0 * float(rest_values[-1])
            component_record[f"charge_start_{label}_mV"] = 1000.0 * float(charge_values[0])
            component_record[f"step_change_{label}_mV"] = 1000.0 * float(
                charge_values[0] - rest_values[-1]
            )
        component_rows.append(component_record)
        for cell in sorted(context.curves):
            grid = context.grids[(rate, cell)]
            exp = context.experimental_blocks[(rate, cell)]
            pred = np.interp(np.minimum(grid, t[-1]), t, v)
            err = 1000.0 * (pred - exp)
            rows.append(
                {
                    "scenario": "continuous_history",
                    "cell": cell,
                    "C_rate": rate,
                    "model_charge_start_V": float(v[0]),
                    "experiment_charge_start_V": float(exp[0]),
                    "initial_error_mV": float(err[0]),
                    "full_RMSE_mV": float(np.sqrt(np.mean(err**2))),
                    "full_bias_mV": float(np.mean(err)),
                    "early_0_10_RMSE_mV": float(np.sqrt(np.mean(err[:10] ** 2))),
                    "middle_10_70_RMSE_mV": float(np.sqrt(np.mean(err[10:70] ** 2))),
                    "late_70_100_RMSE_mV": float(np.sqrt(np.mean(err[70:] ** 2))),
                    "model_time_margin_min": float(t[-1] - grid[-1]),
                    "model_CC_capacity_Ah": float(q[-1]),
                }
            )
        for ti, vi, qi in zip(t, v, q):
            curve_rows.append({"C_rate": rate, "t_min": ti, "V": vi, "Q_Ah": qi})
    detail = pd.DataFrame(rows)
    summary = (
        detail.groupby("C_rate", as_index=False)
        .agg(
            initial_error_mV=("initial_error_mV", "mean"),
            full_RMSE_mV=("full_RMSE_mV", lambda x: float(np.sqrt(np.mean(np.asarray(x) ** 2)))),
            full_bias_mV=("full_bias_mV", "mean"),
            early_0_10_RMSE_mV=("early_0_10_RMSE_mV", lambda x: float(np.sqrt(np.mean(np.asarray(x) ** 2)))),
            middle_10_70_RMSE_mV=("middle_10_70_RMSE_mV", lambda x: float(np.sqrt(np.mean(np.asarray(x) ** 2)))),
            late_70_100_RMSE_mV=("late_70_100_RMSE_mV", lambda x: float(np.sqrt(np.mean(np.asarray(x) ** 2)))),
            min_time_margin_min=("model_time_margin_min", "min"),
            model_CC_capacity_Ah=("model_CC_capacity_Ah", "mean"),
        )
    )
    detail.to_csv(OUT / "continuous_protocol_detail.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(OUT / "continuous_protocol_by_rate.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(curve_rows).to_csv(OUT / "continuous_protocol_charge_curves.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(step_rows).to_csv(OUT / "continuous_protocol_steps.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(component_rows).to_csv(OUT / "continuous_protocol_onset_components.csv", index=False, encoding="utf-8-sig")
    print(detail.to_string(index=False))
    print("\nSUMMARY")
    print(summary.to_string(index=False))
    print("\nSTEPS")
    print(pd.DataFrame(step_rows).to_string(index=False))
    print("\nONSET COMPONENTS")
    print(pd.DataFrame(component_rows).to_string(index=False))


if __name__ == "__main__":
    main()
