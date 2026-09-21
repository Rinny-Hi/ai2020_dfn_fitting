# Experimental kinetics fixed + Bruggeman screen

## Fixed values

- `kn = 7.40e-7` (N4-only EIS)
- `kp = 3.12e-7` (positive-electrode EIS)
- `Rn = 5 um`, `Rp = 3 um`
- `Dsn = 2.1e-14 m2/s`, `Dsp = 4.4e-14 m2/s`
- OCP and stoichiometry window unchanged (`qOCV MAE = 4.64 mV`)

## Screen result

With `bs=1.5`, the Pareto electrode-Bruggeman candidates were:

| bn | bp | Mean voltage MAE (mV) | Mean absolute capacity error (mAh) | Max absolute capacity error (mAh) |
|---:|---:|---:|---:|---:|
| 2.5 | 1.50 | 15.40 | 40.98 | 118.02 |
| 2.5 | 1.83 | 13.67 | 43.00 | 123.06 |
| 2.5 | 2.00 | 13.41 | 49.40 | 126.13 |

The lowest grid voltage MAE was obtained at `bn=2.5, bp=2.0, bs=1.5`, but its capacity error was larger. Separator refinement produced a marginally lower voltage MAE near `bn=2.5, bp=2.0, bs=1.25`, but `bs=1.25` is less defensible than the conventional 1.5 without independent tortuosity evidence.

## Numerical compromise found in this screen

`bn=2.5, bp=1.83, bs=1.5` was the grounded numerical compromise within this screen. It keeps the measured EIS kinetics and the established positive/separator Bruggeman values, while making only a modest change to the negative electrode coefficient.

- Mean MAE: `13.67 mV`
- Charge mean MAE: `10.74 mV`
- Discharge mean MAE: `16.60 mV`
- Mean absolute capacity error: `43.00 mAh`
- Maximum absolute capacity error: `123.06 mAh` (0.5C charge)

The remaining limitation is the charge endpoint/capacity mismatch, especially at 0.5C. It should not be removed by unconstrained Bruggeman tuning because Bruggeman primarily controls electrolyte transport, whereas the persistent endpoint error also depends on usable capacity, stoichiometric endpoints, and voltage cutoffs.

> Current-status note: after this screen, `bn` was restored to the Ai2020 nominal value `2.914` by user instruction. This document preserves the screen result and is not the final current configuration; see `260921_CURRENT_PROGRESS.md`.
