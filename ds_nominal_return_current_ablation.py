"""Current-baseline ablation of GITT-neutral versus Ai2020 nominal Ds."""

from pathlib import Path

import pandas as pd

import c50_selected_ocp_dynamic_validation as selected
import paper_geometry_area_porosity_ablation as area_study
import prefit_physical_anchor_decision as p0
import p0_five_parameter_sensitivity_area as sens


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "260927_ds_nominal_return_current_ablation"
RESULTS.mkdir(parents=True, exist_ok=True)


def main():
    curves, _, protocol = sens.load_july()
    model0, stage, anode, cathode, _, _, qcell = selected.selected_inputs()
    p0_base = p0.build_model(model0, p0.Scenario("P0 anchored reference"))
    geometry = area_study.apply_paper_overlap(p0_base, preserve_capacity=True)
    kp = sens.kp_from_area_specific_rct(geometry)["kp"]
    common = {**sens.INITIAL_COMMON, "kp": kp}
    cases = {
        "GITT-neutral Dsn/Dsp": (2.30e-14, 4.36e-14),
        "Ai2020 Dsn only": (3.90e-14, 4.36e-14),
        "Ai2020 Dsp only": (2.30e-14, 5.387e-15),
        "Ai2020 Dsn+Dsp": (3.90e-14, 5.387e-15),
    }
    details = []
    for name, (dsn, dsp) in cases.items():
        values = {**common, "Dsn": dsn, "Dsp": dsp}
        model = sens.apply_parameters(geometry, values)
        print(name, "charge", flush=True)
        charge = sens.simulate(model, stage, anode, cathode, protocol, True)
        print(name, "discharge", flush=True)
        discharge = sens.simulate(model, stage, anode, cathode, protocol, False)
        frame = sens.summarize_runs(
            name, {True: charge, False: discharge}, curves, qcell
        ).rename(columns={"Area_policy": "Case"})
        frame["Dsn_m2_s"] = dsn
        frame["Dsp_m2_s"] = dsp
        details.append(frame)
    detail = pd.concat(details, ignore_index=True)
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
