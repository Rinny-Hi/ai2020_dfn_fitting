"""Compare GITT and Ai2020 Dsn at the ORegan b_n=2.914-start fit point."""

from pathlib import Path

import pandas as pd

import cumulative_six_stage_model_comparison as cumulative
import electrolyte_transport_set_comparison as transport
import s7_pso_voltage_fitting as s7


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "261020_final_dsn_ai2020_vs_gitt"
COMMON = {
    "kn": 8.726541044267879e-7,
    "brugg_n": 2.533583406033751,
}
SCENARIOS = {
    "GITT Dsn": 4.1707741459682473e-14,
    "Ai2020 Dsn": 3.9e-14,
}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    _, _, _, _, experiment, qcell = transport.current_p0()
    entry = cumulative.build_stages(qcell)["S6"]
    model = entry["model"]
    stage = entry["stage"]
    anode = entry["anode"]
    cathode = entry["cathode"]

    details = []
    potentials = []
    for label, dsn in SCENARIOS.items():
        values = {"Dsn": dsn, **COMMON}
        changed = s7.apply_values(model, values)
        runs = s7.run_model(changed, stage, anode, cathode)
        result = {**values, "runs": runs}
        branch_rows, anode_rows = s7.scenario_metrics(
            label, result, experiment, qcell
        )
        details.extend(branch_rows)
        potentials.extend(anode_rows)

    detail = pd.DataFrame(details)
    summary = s7.aggregate_metrics(detail)
    anode_summary = pd.DataFrame(potentials)
    detail.to_csv(OUT / "branch_metrics.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(OUT / "summary.csv", index=False, encoding="utf-8-sig")
    anode_summary.to_csv(
        OUT / "anode_potential.csv", index=False, encoding="utf-8-sig"
    )
    print(summary.to_string(index=False))
    print(anode_summary.to_string(index=False))


if __name__ == "__main__":
    main()
