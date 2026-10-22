import os
from pathlib import Path

import nbformat
from nbclient import NotebookClient


ROOT = Path(__file__).resolve().parent.parent
NOTEBOOK = ROOT / "outputs" / "261020_Ai2020_1단계_OCP_Stoichiometry_Measured_Anode_비교.ipynb"
RUN_DIR = ROOT / "work" / "run_inputs"
VENV = ROOT / "work" / "venv-pybamm2671"

os.environ["JUPYTER_PATH"] = str(VENV / "share" / "jupyter")
notebook = nbformat.read(NOTEBOOK, as_version=4)
notebook.cells[2].source = notebook.cells[2].source.replace(
    "%pip install -q pybamm==26.7.1.0 openpyxl\n", ""
)
notebook.metadata.kernelspec = {
    "display_name": "PyBaMM 26.7.1",
    "language": "python",
    "name": "pybamm2671",
}
client = NotebookClient(
    notebook,
    timeout=3600,
    kernel_name="pybamm2671",
    resources={"metadata": {"path": str(RUN_DIR)}},
    allow_errors=False,
)

try:
    client.execute()
except Exception:
    partial = ROOT / "work" / "measured_anode_partial.ipynb"
    nbformat.write(notebook, partial)
    print(f"Partial notebook: {partial}", flush=True)
    raise
else:
    nbformat.write(notebook, NOTEBOOK)
    print(f"Executed notebook: {NOTEBOOK}", flush=True)
