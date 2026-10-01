# 260918 GITT OCP Stage 1 comparison

## 1. Data QC

- v2: full experiment 254.07 h; GITT section 196.81 h; 90 pulses and 91 two-hour rests. Last-10-min |dV/dt| median/mean/P90/max = 0.0088/0.0174/0.0234/0.3330 mV/min.
- v3: full experiment 260.02 h; GITT section 203.36 h; 93 pulses and 94 two-hour rests. Last-10-min |dV/dt| median/mean/P90/max = 0.0086/0.0119/0.0169/0.3060 mV/min.
- GITT usable capacity: v2 lithiation/delithiation = 1.790/1.760 mAh; v3 = 1.860/1.830 mAh. Replicate differences are 3.84% and 3.90%.
- Replicate mean-OCP difference: full-range MAE/RMSE = 3.19/8.61 mV; 10–90% MAE/RMSE = 1.32/2.06 mV. The full-range maximum is 64.4 mV and is endpoint-driven.
- The legacy bundle contains only aggregated 1 h GITT rest-end points, not the underlying final-10-min record data. A like-for-like 1 h vs 2 h |dV/dt| comparison is therefore not available.

## 2. OCP comparison

- v2 lithiation-delithiation hysteresis over 10–90%: MAE 11.86 mV, RMSE 12.56 mV, maximum 21.38 mV.
- v3 lithiation-delithiation hysteresis over 10–90%: MAE 12.53 mV, RMSE 13.30 mV, maximum 22.54 mV.
- v2 and v3 were normalized by their own usable branch capacity before interpolation and averaging. Raw GITT capacity was not converted directly to absolute stoichiometry.
- For this first-pass candidate, the normalized usable range was mapped to candidate x=0..1. This improves nominal measured coverage but does not establish the absolute stoichiometry anchor.

## 3. Stage 1 comparison

- Legacy: x0/x100 = 0.001417/0.773635; y100/y0 = 0.445100/0.953391; qOCV RMSE full/2–98% = 18.64/18.04 mV.
- New GITT: x0/x100 = 0.000001/0.772220; y100/y0 = 0.414786/0.923078; qOCV RMSE full/2–98% = 32.86/28.73 mV.
- New GITT selected the x0 lower bound and left endpoint errors of 9.40 mV (SOC 0) and -0.81 mV (SOC 100).

## 4. Dynamic baseline

- Legacy completed 6/6 nominal DFN cases.
- New GITT completed 3/3 discharge cases and 0/3 charge cases. All New GITT charge cases failed at initialization with x0 near zero.
- 0.5C discharge full-overlap RMSE: Legacy 61.82 mV; New GITT 45.42 mV.
- 1C discharge full-overlap RMSE: Legacy 68.48 mV; New GITT 54.48 mV.
- 2C discharge full-overlap RMSE: Legacy 54.82 mV; New GITT 44.29 mV.

## 5. Anode potential

- Legacy 0.5C charge: minimum 0.0184 V at SOC 0.910 and 111.99 min; time fraction below 0 V = 0.000.
- Legacy 1C charge: minimum -0.0220 V at SOC 0.823 and 50.61 min; time fraction below 0 V = 0.169.
- Legacy 2C charge: minimum -0.0805 V at SOC 0.651 and 20.03 min; time fraction below 0 V = 0.680.
- New GITT charge anode-potential metrics are unavailable because all charge simulations failed at initialization.
- A modeled anode potential below 0 V is treated only as a DFN electrochemical-state/risk indicator, not proof of lithium plating.

## 6. Recommendation

**Do not adopt yet**

The 2 h GITT experiment has good central-range replicate reproducibility and the New GITT OCP improves the three discharge baseline RMSE values. However, it worsens Stage 1 qOCV error, drives x0 to the lower bound, fails to satisfy the SOC 0 endpoint closely, and makes all nominal charge DFN simulations fail at initialization. The absolute stoichiometry anchor of the normalized capacity axis also remains unverified.

Before Stage 2–4 is considered:

1. establish an independent absolute stoichiometry/capacity anchor for the new half-cells;
2. investigate the low-x endpoint mismatch and endpoint replicate divergence without tuning dynamic parameters;
3. resolve the x0 boundary solution and reproduce all three charge baselines;
4. then reassess surface-stoichiometry coverage and anode potential for New GITT.

Stage 2–4 fitting was not run.