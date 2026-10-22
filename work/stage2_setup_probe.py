from pathlib import Path
import json
import os
import time

ROOT = Path(__file__).resolve().parents[1]
os.chdir(ROOT / "work" / "stage2_inputs")
notebook = json.loads(next((ROOT / "outputs").glob("*Measured_Anode*ipynb")).read_text(encoding="utf-8"))
scope = {"__name__": "__main__"}
for index in (2, 3, 5, 7, 8, 9, 12, 13):
    code = "".join(notebook["cells"][index]["source"])
    if index == 2:
        code = "\n".join(line for line in code.splitlines() if not line.startswith("%pip"))
        code = code.replace('REQUIRED_PYBAMM = "26.7.1.0"', 'REQUIRED_PYBAMM = "26.8.0.0"')
    if index == 8:
        code = '''
def build_cathode_gitt_ocp_compact(book_name, parameters):
    frame = pd.read_excel(book_name, sheet_name="CATHODE_GITT", header=8)
    disk_area_m2 = np.pi * (0.01 / 2.0) ** 2
    disk_capacity_ah = (
        FARADAY_CONSTANT
        * float(parameters["Maximum concentration in positive electrode [mol.m-3]"])
        * float(parameters["Positive electrode active material volume fraction"])
        * float(parameters["Positive electrode thickness [m]"])
        * disk_area_m2 / 3600.0
    )
    lithiation_curves, delithiation_curves = [], []
    for channel in ("7-7", "7-8"):
        block = frame.loc[frame["source_channel"].astype(str).eq(channel)]
        charge = block.loc[block["pulse_direction"].eq("charge_pulse"), ["Qcum_Ah", "Vrest_V"]].to_numpy(float)
        discharge = block.loc[block["pulse_direction"].eq("discharge_pulse"), ["Qcum_Ah", "Vrest_V"]].to_numpy(float)
        delithiation_curves.append(clean_curve(1.0 - charge[:, 0] / disk_capacity_ah, charge[:, 1]))
        discharge_span = discharge[:, 0].max() / disk_capacity_ah
        lithiation_curves.append(clean_curve((1.0 - discharge_span) + discharge[:, 0] / disk_capacity_ah, discharge[:, 1]))
    return mean_overlap([mean_overlap(lithiation_curves), mean_overlap(delithiation_curves)])

CATHODE_GITT_OCP = build_cathode_gitt_ocp_compact(MEASURED_ANODE_BOOK, raw_parameters)
''' + code[code.index('def build_measured_anode_lowrate_ocp('):]
    start = time.monotonic()
    print("EXEC CELL", index, flush=True)
    exec(compile(code, f"notebook_cell_{index}", "exec"), scope)
    print("CELL DONE", index, round(time.monotonic() - start, 2), flush=True)

stage = {
    "parameters": scope["improved_parameters"],
    "window": scope["IMPROVED_WINDOW"],
    "c_s_max_n": scope["improved_csmax_n"],
    "c_s_max_p": scope["improved_csmax_p"],
    "anode_ocp": scope["AI2020_OCP"]["negative"],
    "cathode_ocp": scope["CATHODE_GITT_OCP"],
}
if __name__ == "__main__":
    start = time.monotonic()
    result = scope["run_dfn_branch"](stage, 1.0, True)
    print("SIM DONE", round(time.monotonic() - start, 2), len(result["time_min"]), result["voltage_v"][0], result["voltage_v"][-1], flush=True)
