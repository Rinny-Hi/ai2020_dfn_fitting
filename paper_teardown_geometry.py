"""Geometry bookkeeping from the SPB655060 teardown paper.

The teardown reports different positive and negative sheet footprints.  A 1D
DFN, however, has one through-plane current-carrying area.  We therefore use
the cathode footprint (the smaller overlap) as the DFN reaction area and retain
the full negative footprint separately for mass/inventory bookkeeping.
"""

from __future__ import annotations

from typing import Any


N_COATED_FACES = 34.0

NEGATIVE_WIDTH_M = 46.4e-3
NEGATIVE_LENGTH_M = 52.0e-3
POSITIVE_WIDTH_M = 45.4e-3
POSITIVE_LENGTH_M = 51.0e-3

NEGATIVE_CURRENT_COLLECTOR_M = 10.0e-6
POSITIVE_CURRENT_COLLECTOR_M = 15.0e-6

# Values reported by the teardown source.  Porosity was measured by mercury
# porosimetry; active-material fractions were then estimated from porosity and
# assumed inactive-material volume fractions (5% positive, 7% negative).
NEGATIVE_POROSITY = 0.32
POSITIVE_POROSITY = 0.33
NEGATIVE_ACTIVE_MATERIAL_FRACTION = 0.61
POSITIVE_ACTIVE_MATERIAL_FRACTION = 0.62

NEGATIVE_TOTAL_AREA_M2 = NEGATIVE_WIDTH_M * NEGATIVE_LENGTH_M * N_COATED_FACES
POSITIVE_TOTAL_AREA_M2 = POSITIVE_WIDTH_M * POSITIVE_LENGTH_M * N_COATED_FACES
DFN_OVERLAP_AREA_M2 = POSITIVE_TOTAL_AREA_M2


def geometry_metadata() -> dict[str, float | str]:
    return {
        "negative_total_area_m2": NEGATIVE_TOTAL_AREA_M2,
        "positive_total_area_m2": POSITIVE_TOTAL_AREA_M2,
        "dfn_overlap_area_m2": DFN_OVERLAP_AREA_M2,
        "negative_total_area_cm2": NEGATIVE_TOTAL_AREA_M2 * 1e4,
        "positive_total_area_cm2": POSITIVE_TOTAL_AREA_M2 * 1e4,
        "dfn_overlap_area_cm2": DFN_OVERLAP_AREA_M2 * 1e4,
        "number_of_coated_faces": N_COATED_FACES,
        "negative_current_collector_um": NEGATIVE_CURRENT_COLLECTOR_M * 1e6,
        "positive_current_collector_um": POSITIVE_CURRENT_COLLECTOR_M * 1e6,
        "source_negative_porosity_measured": NEGATIVE_POROSITY,
        "source_positive_porosity_measured": POSITIVE_POROSITY,
        "source_negative_active_fraction_estimated": NEGATIVE_ACTIVE_MATERIAL_FRACTION,
        "source_positive_active_fraction_estimated": POSITIVE_ACTIVE_MATERIAL_FRACTION,
        "area_policy": (
            "Use distinct teardown footprints for dry-mass/inventory bookkeeping; "
            "use the smaller cathode footprint as the common 1D DFN current area."
        ),
    }


def apply_overlap_area(model: dict[str, Any]) -> dict[str, Any]:
    """Apply the paper cathode footprint as the common 1D DFN area."""
    changed = dict(model)
    changed["params"] = model["params"].copy()
    changed["params"].update(
        {
            "Electrode height [m]": POSITIVE_LENGTH_M,
            "Electrode width [m]": POSITIVE_WIDTH_M,
            "Number of electrodes connected in parallel to make a cell": N_COATED_FACES,
            "Negative current collector thickness [m]": NEGATIVE_CURRENT_COLLECTOR_M,
            "Positive current collector thickness [m]": POSITIVE_CURRENT_COLLECTOR_M,
        },
        check_already_exists=False,
    )
    changed["paper_geometry"] = geometry_metadata()
    return changed
