# Stage 2 seven-parameter sensitivity and identifiability

Current OCP/stoichiometry/geometry were fixed. The candidate set was `kn, kp, Dsn, Dsp, brugg_n, brugg_p, brugg_s`.
Positive parameters were perturbed by ±5% in log-effect coordinates. Each Bruggeman coefficient was perturbed to create the same ±5% change in its local `eps_i**b_i`, making RMS effects comparable across all candidates.

## Selected subsets

- C-rate combined: `brugg_n, kn, brugg_p`
- C-rate charge: `kn, brugg_n, Dsn`
- C-rate discharge: `brugg_n, kp`
- HPPC: `kn, brugg_n, Dsn, brugg_p`

## Ranking

### C-rate combined

| Parameter | RMS (mV/log-effect) | Relative | Selected |
|---|---:|---:|---|
| brugg_n | 71.125 | 1.000 | yes |
| kn | 43.305 | 0.609 | yes |
| kp | 43.155 | 0.607 | no |
| brugg_p | 15.698 | 0.221 | yes |
| brugg_s | 7.691 | 0.108 | no |
| Dsn | 6.804 | 0.096 | no |
| Dsp | 1.234 | 0.017 | no |

### HPPC

| Parameter | RMS (mV/log-effect) | Relative | Selected |
|---|---:|---:|---|
| kn | 28.119 | 1.000 | yes |
| kp | 27.014 | 0.961 | no |
| brugg_n | 8.410 | 0.299 | yes |
| Dsn | 5.926 | 0.211 | yes |
| brugg_p | 3.139 | 0.112 | yes |
| brugg_s | 1.654 | 0.059 | no |
| Dsp | 0.931 | 0.033 | no |

## Interpretation guardrails

- Selection is local to the current parameter point and the 5% effect scale.
- HPPC raw pulse records contain only pulse start/end points; diffusion identifiability must therefore be interpreted conservatively.
- A parameter excluded from one dataset can still be retained in the other dataset's independent fitting path.
- Stage 3 should fit only the subset selected from its own training dataset to prevent validation leakage.
