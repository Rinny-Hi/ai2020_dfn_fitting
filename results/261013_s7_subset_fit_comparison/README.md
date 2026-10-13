# S7 subset fitting comparison

All fits minimize the same fixed-physical-time terminal-voltage RMSE over the
six Cell 6-6 branches. Capacity and anode potential are validation metrics,
not fitting targets. `Dsn`, `Dsp`, and `brugg_s` remain fixed.

The standard starting point is hybrid: `kn = 9.648533212e-7` from Ai2020
(numerically almost identical to the negative-electrode EIS v2 estimate) and
`kp = 5.8995199067e-7` from positive-electrode EIS.

## Main comparison

| Scenario | Fitted variables | Objective RMSE | Full RMSE | Center 10-70% RMSE | Capacity RMSE |
|---|---|---:|---:|---:|---:|
| S7a | `brugg_n + kn` | 37.437 mV | 37.146 mV | 20.967 mV | 1.981% |
| S7c | `brugg_n + kn + brugg_p` | 37.216 mV | 37.017 mV | 20.585 mV | 2.953% |
| S7d | `brugg_n + kp + kn + brugg_p` | 37.139 mV | 36.766 mV | 15.427 mV | 2.696% |
| S7e | same four variables, `kn = kp = Ai2020` start | 37.166 mV | 36.928 mV | 20.063 mV | 2.887% |

S7c gives only a marginal voltage improvement and worsens the endpoint-capacity
metric. S7d gives the lowest raw voltage error, but independent starts produce
very different `kn` and `kp` combinations at similar objective values. This is
consistent with the 0.991 `kn`-`kp` sensitivity cosine and means the individual
kinetic parameters are not identifiable from terminal voltage alone. S7a is
therefore retained as the preferred practical model.

Scripts:

- `s7_fit_subset_comparison.py`: S7a/S7c/S7d comparison
- `s7d_second_seed_check.py`: independent S7d seed
- `s7d_ai2020_start_fit.py`: S7e with both kinetics initialized at Ai2020

