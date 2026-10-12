"""Evaluate anode potential for ORegan2022 with b_n = 2.5."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybamm
from scipy.integrate import cumulative_trapezoid

import comprehensive_prefit_cross_cohort as cross
import electrode_hysteresis_quantification as ehq
import gitt_ocp_analysis as ga
import priority_initial_state_protocol_recheck as priority

import electrolyte_transport_set_comparison as comparison


def at_separator(entries: np.ndarray, n_time: int) -> np.ndarray:
    """Return a time series from the separator-side edge of the negative electrode."""
    values = np.squeeze(np.asarray(entries, dtype=float))
    if values.ndim == 1:
        return values
    if values.shape[-1] == n_time:
        return values[-1, :]
    if values.shape[0] == n_time:
        return values[:, -1]
    raise RuntimeError("Unexpected variable shape: {}".format(values.shape))


def time_fraction(time_s: np.ndarray, condition: np.ndarray) -> float:
    if len(time_s) < 2 or time_s[-1] <= time_s[0]:
        return float("nan")
    duration = np.diff(time_s)
    occupancy = 0.5 * (
        condition[:-1].astype(float) + condition[1:].astype(float)
    )
    return float(np.sum(duration * occupancy) / (time_s[-1] - time_s[0]))


def run_branch(
    rate: float,
    charge: bool,
    stage: dict,
    model_inputs: dict,
    anode: dict,
    cathode: dict,
) -> dict:
    model = pybamm.lithium_ion.DFN(
        {"open-circuit potential": ("single", "single"), "thermal": "isothermal"}
    )
    direction = "Charge" if charge else "Discharge"
    step = "{} at {}C until {} V".format(
        direction, rate, 4.2 if charge else 3.0
    )
    simulation = pybamm.Simulation(
        model,
        parameter_values=ehq.make_parameter_values(
            model_inputs, stage, anode, cathode, charge, False, False
        ),
        experiment=pybamm.Experiment([step], period="10 seconds"),
        var_pts={"x_n": 15, "x_s": 15, "x_p": 15, "r_n": 15, "r_p": 15},
        solver=pybamm.IDAKLUSolver(rtol=1e-6, atol=1e-8),
    )
    solution = simulation.solve()
    time_s = np.asarray(solution["Time [s]"].entries, dtype=float)
    terminal_voltage = np.asarray(
        solution["Terminal voltage [V]"].entries, dtype=float
    )
    current = np.asarray(solution["Current [A]"].entries, dtype=float)
    capacity = cumulative_trapezoid(np.abs(current), time_s, initial=0.0) / 3600.0
    anode_potential = np.squeeze(
        np.asarray(
            solution[
                "Negative electrode surface potential difference at separator interface [V]"
            ].entries,
            dtype=float,
        )
    )
    anode_ocp = at_separator(
        solution["Negative electrode open-circuit potential [V]"].entries,
        len(time_s),
    )
    surface_sto = at_separator(
        solution["Negative particle surface stoichiometry"].entries,
        len(time_s),
    )
    return {
        "time_s": time_s,
        "capacity_Ah": capacity,
        "Q_Ah": capacity,
        "V": terminal_voltage,
        "potential_V": anode_potential,
        "ocp_V": anode_ocp,
        "surface_sto": surface_sto,
    }


def main() -> None:
    comparison.OUT.mkdir(parents=True, exist_ok=True)
    baseline, stage, anode, cathode, _, qcell = comparison.current_p0()
    model = comparison.with_transport_set(baseline, "ORegan2022")
    model["params"].update(
        {"Negative electrode Bruggeman coefficient (electrolyte)": 2.5}
    )

    fig, axes = plt.subplots(
        2,
        3,
        figsize=(15.6, 8.5),
        sharex=True,
        sharey="row",
        constrained_layout=True,
    )
    rows = []
    curve_rows = []
    for column, rate in enumerate(ga.RATES):
        for row, charge in enumerate((True, False)):
            direction = "Charge" if charge else "Discharge"
            rest_v, measured_current = comparison.PROTOCOL[(rate, charge)]
            z0 = priority.monotone_voltage_inverse(
                model["v_qocv"], model["soc"], rest_v
            )
            local_stage = priority.stage_at_soc(stage, z0, charge)
            output = run_branch(
                measured_current / cross.NOMINAL_CAPACITY_AH,
                charge,
                local_stage,
                model,
                anode,
                cathode,
            )
            potential = output["potential_V"]
            ocp = output["ocp_V"]
            capacity = output["capacity_Ah"]
            minimum_index = int(np.nanargmin(potential))
            below_zero_fraction = time_fraction(
                output["time_s"], potential < 0.0
            )
            rows.append(
                {
                    "C_rate": rate,
                    "Direction": direction,
                    "Initial_SOC_pct": 100.0 * z0,
                    "Minimum_anode_potential_V": float(potential[minimum_index]),
                    "Capacity_at_minimum_Ah": float(capacity[minimum_index]),
                    "Surface_stoichiometry_at_minimum": float(
                        output["surface_sto"][minimum_index]
                    ),
                    "End_anode_potential_V": float(potential[-1]),
                    "Minimum_anode_OCP_V": float(np.nanmin(ocp)),
                    "Below_0V_time_fraction_pct": 100.0 * below_zero_fraction,
                }
            )

            ax = axes[row, column]
            ax.plot(
                capacity,
                potential,
                color="#0072B2",
                linewidth=2.2,
                label="DFN interfacial anode potential",
            )
            ax.plot(
                capacity,
                ocp,
                color="#555555",
                linewidth=1.7,
                linestyle="--",
                label="Local equilibrium OCP",
            )
            ax.axhline(0.0, color="#D55E00", linewidth=1.2, linestyle=":")
            ax.scatter(
                [capacity[minimum_index]],
                [potential[minimum_index]],
                color="#0072B2",
                s=28,
                zorder=4,
            )
            ax.set_title("{}C {}".format(rate, direction), fontsize=12, weight="bold")
            ax.set_xlim(0.0, 2.42)
            ax.set_ylim((-0.12, 0.35) if charge else (-0.03, 0.85))
            ax.grid(True, color="#B0B0B0", alpha=0.25, linewidth=0.7)
            ax.text(
                0.03,
                0.05,
                "min {:+.1f} mV\n< 0 V: {:.1f}% of time".format(
                    1000.0 * potential[minimum_index],
                    100.0 * below_zero_fraction,
                ),
                transform=ax.transAxes,
                fontsize=9.5,
                va="bottom",
                bbox={
                    "boxstyle": "round,pad=0.3",
                    "facecolor": "white",
                    "edgecolor": "#B0B0B0",
                    "alpha": 0.86,
                },
            )
            if row == 1:
                ax.set_xlabel("Transferred capacity (Ah)")
            if column == 0:
                ax.set_ylabel("Anode potential (V vs. Li/Li+)")

            for index in range(len(capacity)):
                curve_rows.append(
                    {
                        "C_rate": rate,
                        "Direction": direction,
                        "Transferred_capacity_Ah": float(capacity[index]),
                        "Anode_potential_V": float(potential[index]),
                        "Local_equilibrium_OCP_V": float(ocp[index]),
                        "Negative_surface_stoichiometry": float(
                            output["surface_sto"][index]
                        ),
                    }
                )

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=2,
        frameon=False,
        bbox_to_anchor=(0.5, 1.02),
    )
    fig.suptitle(
        "Cell 6-6 anode-potential prediction\n"
        "ORegan2022 electrolyte transport + negative-electrode $b_n=2.5$",
        fontsize=15,
        weight="bold",
        y=1.08,
    )
    figure_path = comparison.OUT / "oregan2022_brugg_n_2p5_anode_potential.png"
    summary_path = comparison.OUT / "oregan2022_brugg_n_2p5_anode_potential_summary.csv"
    curve_path = comparison.OUT / "oregan2022_brugg_n_2p5_anode_potential_curves.csv"
    fig.savefig(figure_path, dpi=240, bbox_inches="tight")
    plt.close(fig)
    summary = pd.DataFrame(rows)
    summary.to_csv(summary_path, index=False, encoding="utf-8-sig")
    pd.DataFrame(curve_rows).to_csv(curve_path, index=False, encoding="utf-8-sig")
    print(figure_path)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
