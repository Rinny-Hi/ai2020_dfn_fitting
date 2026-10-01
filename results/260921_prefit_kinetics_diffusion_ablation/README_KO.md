# Fitting 전 시간-전압 및 k/Ds ablation

## 기준

- OCP: Ai2020 nominal 음극 + old cathode GITT
- Stoichiometry window와 Qn/Qp 고정
- Ln/Lp/Rn/Rp 및 Bruggeman: Ai2020 nominal
- 초기상태: endpoint start
- 음극 hysteresis OFF, 양극 hysteresis ON

## 현재 시간-전압 결과

| 조건 | 10–70% RMSE | 용량오차 |
|---|---:|---:|
| 0.5C Charge | 16.27 mV | -89.26 mAh |
| 0.5C Discharge | 7.46 mV | +27.43 mAh |
| 1C Charge | 12.10 mV | -86.70 mAh |
| 1C Discharge | 12.28 mV | +21.56 mAh |
| 2C Charge | 18.14 mV | -102.58 mAh |
| 2C Discharge | 14.54 mV | -18.21 mAh |

## k 조합 평균

| 조합 | Charge RMSE | Charge 평균오차 | Charge 평균절대오차 | Discharge RMSE | Discharge 평균절대오차 |
|---|---:|---:|---:|---:|---:|
| current kn/current kp | 15.50 mV | -92.85 mAh | 92.85 mAh | 11.43 mV | 22.40 mAh |
| nominal kn/nominal kp | 47.18 mV | +28.36 mAh | 42.88 mAh | 42.28 mV | 20.95 mAh |
| nominal kn/current kp | 16.18 mV | -59.71 mAh | 59.71 mAh | 7.27 mV | 22.25 mAh |
| current kn/nominal kp | 36.08 mV | -0.65 mAh | 32.37 mAh | 31.98 mV | 21.07 mAh |

두 k를 모두 nominal로 되돌리면 용량오차 절대값은 감소하지만 전압 RMSE가 약 3–4배 증가한다. nominal kp가 충전 cutoff 용량을 크게 늘리는 주된 항이지만 전압 형상을 훼손한다. nominal kn/current kp 조합은 용량오차를 줄이면서 전압 RMSE를 거의 유지하므로 fitting 초기 후보로 사용할 수 있다.

## Ds 조합 평균

여기서 nominal Ds는 Ai2020의 농도 의존 함수이며 current Ds는 Dsn=2.1e-14, Dsp=4.4e-14 m2/s 상수이다.

| 조합 | Charge RMSE | Charge 평균오차 | Charge 평균절대오차 | Discharge RMSE | Discharge 평균절대오차 |
|---|---:|---:|---:|---:|---:|
| current Dsn/current Dsp | 15.50 mV | -92.85 mAh | 92.85 mAh | 11.43 mV | 22.40 mAh |
| nominal Dsn/nominal Dsp | 15.96 mV | -134.37 mAh | 134.37 mAh | 24.66 mV | 25.50 mAh |
| nominal Dsn/current Dsp | 14.17 mV | -86.11 mAh | 86.11 mAh | 10.47 mV | 33.69 mAh |
| current Dsn/nominal Dsp | 17.57 mV | -141.07 mAh | 141.07 mAh | 25.60 mV | 23.84 mAh |

Nominal Dsp 함수가 포함된 두 조합은 충전 용량오차를 크게 악화한다. 현재 Dsp=4.4e-14 m2/s를 유지하는 것이 낫다. Dsn은 nominal 함수로 되돌릴 때 전압 RMSE와 충전 용량오차가 소폭 개선되지만 방전 용량오차가 증가하므로 아직 확정하지 않는다.

## 판단

1. kn/kp를 둘 다 nominal로 복귀하지 않는다.
2. nominal kn/current kp를 fitting 초기 후보로 추가한다.
3. Dsp는 현재 4.4e-14 m2/s를 유지한다.
4. Dsn은 current 상수와 nominal 함수 모두 초기 후보로 남긴다.
5. 용량오차만 최소화하지 않고 전압 RMSE와 cutoff 용량오차를 함께 평가한다.
