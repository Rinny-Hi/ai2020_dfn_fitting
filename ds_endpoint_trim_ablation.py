"""Ablate no-trim GITT representatives against the current neutral Ds values."""

from pathlib import Path

import pandas as pd

import c50_selected_ocp_dynamic_validation as selected
import paper_geometry_area_porosity_ablation as area_study
import prefit_physical_anchor_decision as p0
import p0_five_parameter_sensitivity_area as sens


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260927_ds_endpoint_trim_ablation"
RESULTS.mkdir(parents=True, exist_ok=True)


def main():
    curves, _, protocol = sens.load_july()
    model0, stage, anode, cathode, _, _, qcell = selected.selected_inputs()
    geometry = area_study.apply_paper_overlap(
        p0.build_model(model0, p0.Scenario("P0 anchored reference")),
        preserve_capacity=True,
    )
    kp = sens.kp_from_area_specific_rct(geometry)["kp"]
    common = {**sens.INITIAL_COMMON, "kp": kp}
    cases = {
        "10% neutral reference": (2.30e-14, 4.36e-14),
        "0% trim Dsn only": (4.422222641191311e-14, 4.36e-14),
        "0% trim Dsp only": (2.30e-14, 6.809001552344235e-14),
        "0% trim Dsn+Dsp": (4.422222641191311e-14, 6.809001552344235e-14),
    }
    frames = []
    for name, (dsn, dsp) in cases.items():
        model = sens.apply_parameters(geometry, {**common, "Dsn": dsn, "Dsp": dsp})
        charge = sens.simulate(model, stage, anode, cathode, protocol, True)
        discharge = sens.simulate(model, stage, anode, cathode, protocol, False)
        frame = sens.summarize_runs(
            name, {True: charge, False: discharge}, curves, qcell
        ).rename(columns={"Area_policy": "Case"})
        frame["Dsn_m2_s"] = dsn
        frame["Dsp_m2_s"] = dsp
        frames.append(frame)
    detail = pd.concat(frames, ignore_index=True)
    summary = (
        detail.groupby(["Case", "Direction", "Dsn_m2_s", "Dsp_m2_s"], as_index=False)
        .agg(
            Full_RMSE_mV=("Full_RMSE_mV", "mean"),
            Center10_70_RMSE_mV=("Center10_70_RMSE_mV", "mean"),
            Capacity_RMSE_pct=(
                "Capacity_error_pct",
                lambda x: float((sum(float(v) ** 2 for v in x) / len(x)) ** 0.5),
            ),
            Mean_capacity_error_mAh=("Capacity_error_mAh", "mean"),
        )
        .sort_values(["Direction", "Center10_70_RMSE_mV"])
    )
    detail.to_csv(RESULTS / "detail.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(RESULTS / "summary.csv", index=False, encoding="utf-8-sig")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
