"""Build the QC-only GITT preparation notebook."""

from pathlib import Path

import nbformat as nbf


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "260918 GITT OCP QC and Stage1 Preparation.ipynb"


def markdown(text: str):
    return nbf.v4.new_markdown_cell(text.strip())


def code(text: str):
    return nbf.v4.new_code_cell(text.strip())


notebook = nbf.v4.new_notebook()
notebook["metadata"] = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": "3.12"},
}
notebook["cells"] = [
    markdown(
        """
# 260918 GITT OCP QC and Stage 1 Preparation

This notebook performs parser and QC preparation only. It keeps v2 and v3 independent until normalized comparisons are calculated.

Guardrails:

- `Q_NORM_NEW` is a normalized position, not absolute stoichiometry.
- The new candidate is not injected into Stage 1.
- No parameter sensitivity or Stage 2–4 fitting is run.
- The raw Neware workbooks are read-only inputs.
"""
    ),
    code(
        """
from pathlib import Path
from IPython.display import Image, display
import pandas as pd

from gitt_ocp_qc import ANODE_OCP_CASE, run_qc, select_anode_ocp_source

ROOT = Path.cwd().resolve()
QC = run_qc(ROOT)
RESULTS = QC["paths"].results
print("Results:", RESULTS)
"""
    ),
    markdown("## Section A — Raw data parser"),
    code("display(QC['parser_status'])"),
    markdown(
        "The parser locates GITT steps from step type, cycle, duration, and sequence. Record rows are assigned by timestamp boundaries rather than hard-coded row numbers."
    ),
    markdown("## Section B — GITT pulse / rest identification"),
    code(
        """
display(QC["pulse_rest"].head(12))
display(QC["pulse_rest"].groupby(["cell", "direction"]).size().rename("pulse_count").reset_index())
"""
    ),
    markdown("## Section C — 1 h vs 2 h relaxation comparison"),
    code(
        """
display(QC["relaxation_metrics"])
display(Image(filename=str(RESULTS / "anode_rest_1h_vs_2h.png")))
"""
    ),
    markdown("## Section D — Rest equilibrium quality"),
    code(
        """
print("Method: linear-regression slope over the final 10 minutes of each 2 h rest.")
display(QC["equilibrium_metrics"])
display(Image(filename=str(RESULTS / "anode_rest_equilibrium.png")))
"""
    ),
    markdown("## Section E — Capacity QC"),
    code(
        """
display(QC["capacity_qc"])
print("The 2.20–2.25 mAh reference is a full-cell operating-window sanity check only. No correction is applied.")
"""
    ),
    markdown("## Section F — v2 / v3 normalized reproducibility"),
    code(
        """
display(QC["reproducibility"])
display(Image(filename=str(RESULTS / "anode_v2_v3_normalized_reproducibility.png")))
"""
    ),
    markdown("## Section G — Representative new anode OCP candidate"),
    code(
        """
Q_NORM_NEW = QC["candidate"]["Q_NORM_NEW"].to_numpy()
UN_NEW_NORM_CANDIDATE = QC["candidate"]["UN_NEW_NORM_CANDIDATE_V"].to_numpy()
display(QC["candidate"].head())
print("Axis status: normalized position only; absolute stoichiometry is not assigned.")
"""
    ),
    markdown("## Section H — Legacy vs new normalized OCP shape"),
    code(
        """
display(QC["legacy_metrics"])
display(Image(filename=str(RESULTS / "legacy_vs_new_anode_ocp_normalized.png")))
"""
    ),
    markdown("## Section I — OCP source switch safeguard"),
    code(
        """
ANODE_OCP_CASE = "legacy"
display(pd.DataFrame([select_anode_ocp_source(ANODE_OCP_CASE)]))
display(pd.DataFrame([select_anode_ocp_source("new_2h_candidate")]))
"""
    ),
    markdown(
        "`new_2h_candidate` remains `NOT_READY` unless an explicit absolute stoichiometry mapping is supplied."
    ),
    markdown("## Section J — Stage 1 comparison skeleton"),
    code("display(QC['stage1_skeleton'])"),
    markdown(
        "Case B–D are intentionally not executed. The existing endpoint-constrained Stage 1 method is not changed."
    ),
    markdown("## Cathode 2 h GITT placeholder"),
    code("display(QC['cathode_placeholder'])"),
    markdown("## Saved outputs"),
    code(
        """
for path in sorted(RESULTS.iterdir()):
    print(path.name)
"""
    ),
]

nbf.write(notebook, OUTPUT)
print(OUTPUT)
