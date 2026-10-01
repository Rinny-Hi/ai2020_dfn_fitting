# Stage 2 five-parameter sensitivity and identifiability

Current OCP/stoichiometry/geometry were fixed. The candidate set was `kn, kp, Dsn, Dsp, brugg_n`.
Positive parameters were perturbed by ±5% in log-effect coordinates. `brugg_n` was perturbed to create the same ±5% change in `eps_n**b`, making its RMS effect comparable to the other candidates.

## Selected subsets

- C-rate combined: `brugg_n, kn`
- C-rate charge: `kn, brugg_n, Dsn`
- C-rate discharge: `brugg_n, kp`
- HPPC: `kn, brugg_n, Dsn`

## Ranking

### C-rate combined

| Parameter | RMS (mV/log-effect) | Relative | Selected |
|---|---:|---:|---|
| brugg_n | 71.125 | 1.000 | yes |
| kn | 43.305 | 0.609 | yes |
| kp | 43.155 | 0.607 | no |
| Dsn | 6.804 | 0.096 | no |
| Dsp | 1.234 | 0.017 | no |

### HPPC

| Parameter | RMS (mV/log-effect) | Relative | Selected |
|---|---:|---:|---|
| kn | 28.119 | 1.000 | yes |
| kp | 27.014 | 0.961 | no |
| brugg_n | 8.410 | 0.299 | yes |
| Dsn | 5.926 | 0.211 | yes |
| Dsp | 0.931 | 0.033 | no |

## Interpretation guardrails

- Selection is local to the current parameter point and the 5% effect scale.
- HPPC raw pulse records contain only pulse start/end points; diffusion identifiability must therefore be interpreted conservatively.
- A parameter excluded from one dataset can still be retained in the other dataset's independent fitting path.
- Stage 3 should fit only the subset selected from its own training dataset to prevent validation leakage.
