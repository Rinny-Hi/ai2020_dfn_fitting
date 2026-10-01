# 동적 fitting 전 baseline 정리

## 범위

- OCP 및 stoichiometry window는 C/50 데이터로 선정된 상태이다.
- 아래 0.5C/1C/2C 결과에는 동적 데이터에 대한 parameter optimization을 아직 적용하지 않았다.
- 두께, 반경, Bruggeman은 Ai2020 nominal을 사용한다.
- Dsn/Dsp와 kn/kp는 지정값/EIS 값을 고정했다.

## OCP 및 용량

| 항목 | 값 |
|---|---:|
| OCP | Ai2020 nominal anode + old cathode GITT |
| Qcell | 2.35653 Ah |
| Qn | 2.47947 Ah |
| Qp | 4.37120 Ah |
| x0 / x100 | 0.003806 / 0.954221 |
| y100 / y0 | 0.431642 / 0.970746 |
| qOCV MAE / RMSE | 4.69 / 6.24 mV |
| C/50 charge / discharge MAE | 9.12 / 4.48 mV |

## 동적 사전 검증

| C_rate | Direction | Center10_70_MAE_mV | Center10_70_RMSE_mV | Full_RMSE_mV | Initial_voltage_error_mV | Q_end_exp_Ah | Q_end_model_Ah | Capacity_error_mAh |
|---|---|---|---|---|---|---|---|---|
| 0.500 | Charge | 14.363 | 16.269 | 39.181 | 151.250 | 2.146 | 2.057 | -89.258 |
| 0.500 | Discharge | 5.103 | 7.461 | 35.436 | -77.059 | 2.290 | 2.318 | 27.432 |
| 1.000 | Charge | 9.456 | 12.097 | 50.393 | 51.549 | 1.909 | 1.822 | -86.705 |
| 1.000 | Discharge | 10.137 | 12.279 | 32.774 | -115.221 | 2.261 | 2.282 | 21.555 |
| 2.000 | Charge | 13.987 | 18.135 | 56.252 | -45.002 | 1.462 | 1.359 | -102.582 |
| 2.000 | Discharge | 14.053 | 14.541 | 25.421 | -142.054 | 2.214 | 2.196 | -18.208 |


## 해석

- 중앙 SOC 전압 형상은 6조건 모두 약 7–18 mV RMSE 수준이다.
- 충전 CC cutoff 용량은 모든 C-rate에서 모델이 작으며, 동적 fitting의 핵심 잔차이다.
- 방전 용량오차는 충전보다 작으므로 Qn/Qp 전체를 동시에 확대하는 방식은 우선하지 않는다.
- 다음 fitting에서는 OCP/Qn/Qp를 우선 고정하고 actual current, cathode-side capacity/kinetics, 초기 이력을 분리한다.
