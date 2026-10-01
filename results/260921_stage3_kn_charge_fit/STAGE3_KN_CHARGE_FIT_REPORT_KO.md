# 1차 fitting: charge 기반 kn 추정과 discharge validation

## 목적과 고정 조건

- fitting 파라미터: `kn` 하나
- 초기 EIS 값: `kn = 7.40e-7`
- 고정값: `kp=3.12e-7`, `Dsn=2.10e-14 m2/s`, `Dsp=4.40e-14 m2/s`, `brugg_n/p/s=2.914/1.83/1.5`
- OCP와 stoichiometry window 고정, qOCV MAE 4.64 mV 유지
- fitting 데이터: 0.5C, 1C, 2C charge
- validation 데이터: 0.5C, 1C, 2C discharge
- prior 없음. `kn` 범위 안에서 직접 탐색했다.

## 표준 전압 목적함수

세 charge 곡선을 각각 같은 비중으로 둔 ordinary least squares를 사용했다.

`J_V = mean_over_rates(mean_over_Q((V_model - V_exp)^2))`

각 조건은 고정된 실제 통과용량 grid에서 비교했으며, 점 개수가 많은 저율 곡선이 목적함수를 지배하지 않도록 조건별 MSE를 먼저 계산한 뒤 평균했다. 최적화는 MSE로 수행하고 결과는 RMSE와 MAE로 보고했다.

## 용량오차 처리 방법

### 1. Voltage-only LS

용량오차를 목적함수에 넣지 않고 별도 지표로 보고한다. DFN 전압 fitting에서 가장 일반적이고 해석이 단순한 기준이다.

### 2. 다목적 Pareto

`J_V`와 charge 종료용량의 상대 RMSE인 `J_Q`를 별도 목적함수로 둔다. 전압과 용량 어느 한쪽도 동시에 개선할 수 없는 비지배해를 Pareto front로 보고했다. 표시한 knee point는 정규화된 두 목적의 utopia distance가 최소인 설명용 절충점이며 물리적으로 유일한 정답은 아니다.

### 3. 용량 제약

전압 MSE를 최소화하되 charge 용량 RMSE가 0.5%, 1%, 2%, 3%, 5% 이하인 경우만 허용했다. 이 tolerance는 현재 반복실험 산포로부터 정한 값이 아니라 feasibility를 확인하기 위한 시나리오다. 최종 tolerance는 반복 capacity 측정 또는 coulomb-counting 불확도로 정해야 한다.

## Charge fitting 결과

| 방식 | kn | EIS 대비 | Charge 전압 RMSE | Charge 용량 RMSE | 용량 MAE |
|---|---:|---:|---:|---:|---:|
| EIS baseline | 7.400e-7 | 1.00x | 20.73 mV | 7.04% | 126.74 mAh |
| Voltage-only LS | 1.076e-6 | 1.45x | **13.94 mV** | 4.25% | 78.99 mAh |
| Pareto knee | 1.352e-6 | 1.83x | 16.86 mV | 3.01% | 52.65 mAh |
| Capacity constraint 3% | 1.391e-6 | 1.88x | 17.50 mV | **2.91%** | 49.64 mAh |
| Capacity-only reference | 1.706e-6 | 2.31x | 22.78 mV | **2.55%** | 46.12 mAh |

0.5%, 1%, 2% capacity-RMSE 제약은 `kn` 단독 변경으로는 feasible solution이 없었다. 용량만 최소화해도 2.55% 아래로 내려가지 않았고 전압 RMSE는 baseline보다 나빠졌다. 따라서 charge 용량 불일치 전체를 `kn`으로 보정하면 안 된다.

## Discharge blind validation

| 방식 | Discharge 전압 RMSE | Discharge 용량 RMSE | 용량 MAE |
|---|---:|---:|---:|
| EIS baseline | 26.32 mV | 0.97% | 14.60 mAh |
| Voltage-only LS | 12.80 mV | 0.85% | 14.39 mAh |
| Pareto knee | 7.77 mV | 0.80% | 14.30 mAh |
| Capacity constraint 3% | **7.53 mV** | 0.79% | 14.29 mAh |
| Capacity-only reference | 8.91 mV | **0.75%** | 14.21 mAh |

3% 제약 해의 조건별 discharge RMSE는 0.5C 8.71 mV, 1C 9.09 mV, 2C 3.43 mV이다. baseline의 23.75, 27.69, 27.34 mV보다 모두 개선됐다.

## 해석과 1차 선택

- 일반적인 단일 목적 fitting 결과는 `kn = 1.076e-6`이다. 이것이 표준 voltage least-squares 해다.
- 용량을 함께 고려하면 Pareto 구간은 대략 `1.08e-6`에서 `1.71e-6` 사이로 나타난다.
- 현재 시나리오 중 charge 용량 RMSE 3% 제약 해인 `kn = 1.391e-6`이 discharge 전압 validation에서 가장 좋았다.
- 하지만 3% tolerance가 반복실험에서 정해진 값이 아니므로 `1.391e-6`을 최종 물리 `kn`으로 확정할 수는 없다.
- fitted `kn`이 EIS 값의 1.45–1.88배라는 점은 `kn`이 다른 kinetics 또는 전해질 수송 오차를 일부 흡수하고 있을 가능성을 뜻한다.
- 여러 방식을 discharge 결과로 비교했기 때문에 discharge는 이제 모델 선택용 validation이다. 최종 선택 후 HPPC를 추가 외부 validation으로 사용해야 편향되지 않은 최종 검증이 된다.

현재 권고는 **표준 결과로 voltage-only `kn=1.076e-6`을 유지하면서, validation 우수 후보로 3% 제약 `kn=1.391e-6`을 병렬 보존하고 HPPC에서 둘을 비교하는 것**이다.

## 선행연구와의 대응

- Forman et al.은 실험과 DFN 전압 궤적 사이의 L2 제곱오차를 최소화했다: https://doi.org/10.1016/j.jpowsour.2012.03.009
- Reddy et al.도 측정·계산 전압곡선의 least-squares 차이를 최소화했다: https://doi.org/10.1108/COMPEL-12-2018-0533
- Jin et al.은 thermodynamic/kinetic 파라미터를 단계적으로 나누고 pulse 데이터와 별도 실험 데이터로 검증했다: https://doi.org/10.1002/er.4022
- DFN grouping 연구는 전압 fitting 정확도가 높아도 추정 파라미터가 물리적으로 의미 있지 않을 수 있으며, teardown과 추정을 결합해야 한다고 지적한다: https://doi.org/10.1016/j.jpowsour.2021.229901

