"""Check Ds alternatives at fixed measured kinetic constants, one factor at a time."""

import json
from pathlib import Path

import numpy as np
import pandas as pd

from run_stage2_two_ocp import cases, kinetic_constants, scope, stage_with_values


OUT = Path(__file__).resolve().parents[1] / "outputs"
raw = scope["raw_parameters"]
def scalar(value):
    return float(value.evaluate()) if hasattr(value, "evaluate") else float(value)


nominal_n = scalar(raw["Negative particle diffusivity [m2.s-1]"](0.5, 298.15))
nominal_p = scalar(raw["Positive particle diffusivity [m2.s-1]"](0.7, 298.15))
conditions = tuple((rate, charge) for rate in (0.5, 1.0, 2.0) for charge in (True, False))
rows = []
for case_name, case in cases().items():
    values = kinetic_constants(case)
    variants = {
        "nominal_both": (nominal_n, nominal_p),
        "GITT_n_only": (values["Dsn"], nominal_p),
        "GITT_p_only": (nominal_n, values["Dsp"]),
        "GITT_both": (values["Dsn"], values["Dsp"]),
    }
    for label, (dsn, dsp) in variants.items():
        adjusted = {**values, "Dsn": dsn, "Dsp": dsp}
        stage = stage_with_values(case, adjusted)
        metrics = []
        for rate, charge in conditions:
            run = scope["run_dfn_branch"](stage, rate, charge)
            metrics.append(scope["compare_with_experiment"](scope["experimental_curves"][(rate, charge)], run))
        rows.append({
            "OCP_case": case_name,
            "Ds_variant": label,
            "Dsn_m2_s": dsn,
            "Dsp_m2_s": dsp,
            "Full_RMSE_mV": float(np.sqrt(np.mean([m["full_RMSE_mV"] ** 2 for m in metrics]))),
            "Center_10_70_RMSE_mV": float(np.sqrt(np.mean([m["center_10_70_RMSE_mV"] ** 2 for m in metrics]))),
        })
result = pd.DataFrame(rows)
result.to_csv(OUT / "Stage2_Ds_nominal_vs_GITT_fixed_k.csv", index=False, encoding="utf-8-sig")
print(result.to_string(index=False))
