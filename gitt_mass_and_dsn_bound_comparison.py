"""Compare the legacy and remeasured coating masses in the GITT Ds result.

The Neware/Weppner--Huggins expression used by the project is quadratic in
the coating mass.  All voltage, area, molar-volume and pulse selections are
therefore unchanged by a mass-only correction, so the existing audited
pulse/aggregate table can be rescaled exactly without rereading raw xlsx
files.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "results" / "260927_gitt_ds_revalidation" / "trim_aggregation_audit.csv"
OUT = ROOT / "results" / "260927_gitt_mass_and_dsn_bound_comparison"

OLD_MASS_MG = {"Anode": 10.40, "Cathode": 17.25}
CURRENT_MASS_MG = {"Anode": 10.10, "Cathode": 17.35}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    source = pd.read_csv(SOURCE)
    rows: list[dict] = []
    for _, row in source.iterrows():
        electrode = str(row["electrode"])
        old_mass = OLD_MASS_MG[electrode]
        current_mass = CURRENT_MASS_MG[electrode]
        factor = (current_mass / old_mass) ** 2
        old_d = float(row["D_m2_s"])
        rows.append(
            {
                **row.to_dict(),
                "old_mass_mg": old_mass,
                "current_mass_mg": current_mass,
                "current_over_old_D_factor": factor,
                "current_mass_D_m2_s": old_d * factor,
                "D_change_pct": 100.0 * (factor - 1.0),
            }
        )
    comparison = pd.DataFrame(rows)
    comparison.to_csv(OUT / "mass_rescaled_trim_aggregation.csv", index=False, encoding="utf-8-sig")

    selected = comparison[
        (comparison["rounding"] == "floor")
        & comparison["trim_fraction"].isin([0.0, 0.05, 0.10, 0.15, 0.20])
    ].copy()
    selected.to_csv(OUT / "selected_mass_comparison.csv", index=False, encoding="utf-8-sig")

    manifest = {
        "equation_scaling": "D_current = D_old * (m_current / m_old)^2",
        "unchanged_inputs": [
            "pulse voltages",
            "pulse duration",
            "molar mass",
            "molar volume",
            "10 mm geometric punch area",
            "endpoint aggregation",
        ],
        "mass_definition": "single-face composite coating mass = (punched electrode - current collector) / 2",
        "old_mass_mg": OLD_MASS_MG,
        "current_mass_mg": CURRENT_MASS_MG,
        "recommended_fitting_reference": {
            "Dsn_m2_s": float(
                selected[
                    (selected.electrode == "Anode")
                    & (selected.trim_fraction == 0.0)
                    & (selected.direction == "Combined")
                ].iloc[0].current_mass_D_m2_s
            ),
            "Dsp_m2_s": float(
                selected[
                    (selected.electrode == "Cathode")
                    & (selected.trim_fraction == 0.0)
                    & (selected.direction == "Combined")
                ].iloc[0].current_mass_D_m2_s
            ),
        },
        "interpretation": (
            "Use current masses for the mass-corrected apparent-GITT fitting reference, "
            "but retain the old-mass result as a measurement uncertainty case because "
            "replicate punch statistics and active-material mass fractions are unavailable."
        ),
    }
    (OUT / "analysis_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(selected.to_string(index=False))
    print(json.dumps(manifest["recommended_fitting_reference"], indent=2))


if __name__ == "__main__":
    main()
