"""Six-branch no-voltage-fit baseline with current experimental k_n and b_n=2.5."""
from pathlib import Path
import pandas as pd

from oregan_cross_direction_fit import evaluate_both
from run_stage2_two_ocp import cases, kinetic_constants

OUT = Path(__file__).resolve().parents[1] / "outputs"
CASE_NAME = "Ai2020_anode__measured_cathode"
CASE = cases()[CASE_NAME]
kn = kinetic_constants(CASE)["kn"]
fit = {"fit_type": "EIS_kn_bn2p5_before_fit", "fit_direction": "none", "kn": kn, "brugg_n": 2.5}
metrics, curves = evaluate_both(CASE_NAME, CASE, fit)
for row in metrics:
    row["Data_role"] = "Before voltage fit"
pd.DataFrame(metrics).to_csv(OUT / "ORegan2022_Ai2020_bn2.5_피팅전_조건별결과.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(curves).to_csv(OUT / "ORegan2022_Ai2020_bn2.5_피팅전_곡선.csv", index=False, encoding="utf-8-sig")
