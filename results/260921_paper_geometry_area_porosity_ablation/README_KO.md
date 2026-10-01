# SPB655060 원 논문 geometry 적용 영향

## 비교 조건

모든 후보는 현재 선택된 capacity-consistent OCP, 실제 평균전류, 직전 rest 전압 기반 초기 SOC, 음극 hysteresis off/양극 hysteresis on 조건에서 비교했다.

1. 현재 공통면적 `814.98 cm2` + 논문 집전체 `Cu/Al=10/15 um`
2. 논문 overlap 면적 `787.236 cm2`를 그대로 적용한 진단 후보
3. overlap 면적을 적용하되 `c_s,max`를 면적에 반비례시켜 `Qn/Qp` 보존
4. 3번에 원 teardown porosity `n=0.32,p=0.33` 추가

## Geometry 변화량

| 항목 | 값 |
|---|---:|
| 기존 DFN 공통면적 | 814.980 cm2 |
| 원 논문 cathode/overlap 면적 | 787.236 cm2 |
| 면적 변화 | -3.404% |
| 동일 전류에서 전류밀도 변화 | +3.524% |
| Q 보존 시 필요한 `c_s,max` scale | 1.035242 |

집전체 `10/15 um`는 Ai2020 기본값과 이미 동일하다. 현재 isothermal electrochemical DFN에서 collector 두께만 변경한 효과는 0이다.

## 동적 검증 결과

| 후보 | Charge 중앙 RMSE | Discharge 중앙 RMSE | Charge 용량 RMSE | Discharge 용량 RMSE |
|---|---:|---:|---:|---:|
| 현재 면적 + 논문 집전체 | 20.86 mV | 10.60 mV | 3.05% | 0.92% |
| 논문 overlap raw | 17.88 mV | 10.81 mV | 7.47% | 4.19% |
| 논문 overlap + Q 보존 | 19.88 mV | 9.18 mV | 3.71% | 1.01% |
| 논문 overlap + Q 보존 + 원 porosity | 19.11 mV | 10.31 mV | 5.04% | 1.31% |

### 논문 overlap + Q 보존의 C-rate별 결과

| 조건 | 중앙 RMSE | 용량오차 |
|---|---:|---:|
| 0.5C charge | 28.69 mV | -66.69 mAh |
| 1C charge | 17.64 mV | -67.08 mAh |
| 2C charge | 13.32 mV | -64.07 mAh |
| 0.5C discharge | 11.91 mV | +1.65 mAh |
| 1C discharge | 11.44 mV | -4.55 mAh |
| 2C discharge | 4.19 mV | -38.61 mAh |

현재 면적 대비 중앙 전압 RMSE는 charge `-0.98 mV`, discharge `-1.42 mV` 개선되지만, charge capacity RMSE는 `+0.66%p`, discharge는 `+0.09%p` 악화된다.

## 해석

- collector 두께 변경은 결과 변화 원인이 아니다.
- overlap 면적을 Q 재식별 없이 바꾸면 active inventory가 3.4% 감소하므로 용량오차가 크게 증가한다. 이 후보는 물리적으로 불완전한 진단용이다.
- Q를 보존하면 전류밀도 증가에 의해 전압 형상은 소폭 개선되지만 cutoff가 빨라져 충전 용량오차가 커진다.
- 원 논문 porosity를 추가하면 charge 중앙 RMSE는 더 낮아지지만 charge capacity RMSE가 5.04%로 악화된다. 전압 RMSE만 보고 채택하면 안 된다.
- active fraction `n=0.61,p=0.62`는 그대로 유지했다. 원 논문에서도 직접 측정값이 아니라 measured porosity와 inactive fraction 가정에서 계산한 값이므로 loading 측정 전에는 다시 조정하지 않았다.

## 결정

원 teardown 면적은 physical bookkeeping에 채택하되, 현재 DFN 기본 결과를 즉시 paper-overlap 결과로 교체하지 않는다. 내일 dry loading으로 `Qn/Qp`를 독립 계산한 뒤 paper overlap area와 함께 OCP window를 재식별해야 한다.

현재 성능만 기준으로는 기존 공통면적 구성이 용량오차가 가장 작다. paper-overlap 후보는 전압 형상은 조금 낫지만 OCP까지 맞춘 종합 결과에서는 아직 기존 구성을 이기지 못한다.
