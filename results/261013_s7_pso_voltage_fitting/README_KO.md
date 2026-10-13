# S6 marker plot and S7 PSO voltage fitting

## 정의

- S6: 현재 no-fit 기준모델.
- S7a: S6 국소 민감도/식별도 선별 후 `kn+brugg_n`만 PSO fitting.
- S7b: 선별 없이 `Dsn+kn+brugg_n`을 모두 PSO fitting.
- 목적함수: Cell 6-6의 0.5C/1C/2C 충·방전 6곡선 고정 시간격자 terminal-voltage RMSE.
- 용량, cutoff 시간, anode potential은 fitting에 사용하지 않고 사후평가했다.

## 민감도/식별도

| Parameter | Relative sensitivity | 판단 |
|---|---:|---|
| Dsn | 0.069 | 고정 |
| kn | 0.450 | fit 후보 |
| brugg_n | 1.000 | fit 후보 |

`Dsn`의 상대 민감도 0.069가 기준 0.10보다 낮아 GITT anchor로 고정했다. `kn+brugg_n` 부분집합은 full rank, condition number 2.29, 최대 sensitivity cosine 0.679로 기준을 통과했다.

## fitted parameter

| Scenario | Dsn [m2/s] | kn prefactor | brugg_n |
|---|---:|---:|---:|
| S6 | 4.170774e-14 | 9.648533e-07 | 2.500000 |
| S7a | 4.170774e-14 | 7.746510e-07 | 2.494938 |
| S7b | 4.097036e-14 | 7.482399e-07 | 2.510917 |

## 사후 성능

| Scenario | 전체 RMSE | 10–70% RMSE | 용량 RMSE | 최대 용량오차 |
|---|---:|---:|---:|---:|
| S6 | 38.59 mV | 26.30 mV | 1.02% | 29.8 mAh |
| S7a | 37.15 mV | 20.97 mV | 1.98% | 55.1 mAh |
| S7b | 37.14 mV | 19.94 mV | 2.36% | 68.4 mAh |

## 결론

- S7a와 S7b의 전체 전압 RMSE 차이는 약 0.01 mV로 사실상 동일하다.
- S7b에서 Dsn과 brugg_n은 S6 대비 거의 움직이지 않아 세 파라미터 동시 fitting의 추가 실익이 작다.
- fitting 모델 중에는 실험 Dsn을 보존하고 더 적은 파라미터를 사용하는 S7a가 우선이다.
- 그러나 S6 대비 전체 RMSE 개선은 약 1.44 mV에 그치고, 용량 RMSE는 1.02%에서 1.98%로 증가한다. 따라서 현재 기준모델은 S6 유지가 더 보수적이다.
- PSO는 7 particles × 7 iterations, 단일 seed의 탐색 비교이며 전역수렴 또는 optimizer 불확실성을 증명하지 않는다.
