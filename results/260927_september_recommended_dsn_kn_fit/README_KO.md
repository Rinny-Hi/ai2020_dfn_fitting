# September BoL 권장 설정 실행 결과

## 최종 결론

- 9월 charge 데이터에서 사전 `Dsn+kn` 민감도/식별도를 재검증했으나, Dsn 상대민감도가 기준 0.10에 미달했다. 따라서 두 파라미터 동시 fitting은 수행하지 않았다.
- 통과한 축소 모델은 **Dsn 고정 + kn 단일 fitting**이다. 독립 mass-corrected GITT 값 `Dsn=4.170774e-14 m2/s`를 고정하고 `kn=5.348413e-07`를 얻었다.
- charge physical-time RMSE는 57.80 mV, untouched discharge RMSE는 75.80 mV다.
- capacity는 목적함수에 넣지 않았다. 사후 capacity RMSE는 charge 15.82%, discharge 2.09%다.

## 권장 설정 적용

- OCP: Ai2020 nominal graphite equilibrium + old measured cathode equilibrium
- hysteresis: OFF
- `brugg_n/p/s=2.914/1.83/1.5`
- 초기 SOC: charge 1.66--6.68%, discharge 97.33--97.55% (각 cell·rate의 2시간 rest 종점전압으로 개별 산정)
- 유지 OCP window: `x0=0.0038062`, `x100=0.9542213`, `y100=0.4316424`, `y0=0.9707458`

9월에는 독립 full-qOCV가 없어 QLi와 4개 window endpoint를 새로 동시에 식별할 수 없다. 근거 없는 window 재최적화 대신, 검증된 window를 유지하고 branch별 초기 상태만 9월 rest로 갱신했다.

## Dsn+kn 사전 식별도

| Anchor | step | 통과 | 최소 상대민감도 | max abs(cosine) | condition |
|---|---:|---:|---:|---:|---:|
| GITT anchor | 0.010 | False | 0.064 | 0.631 | 2.10 |
| GITT anchor | 0.025 | False | 0.064 | 0.632 | 2.10 |
| GITT anchor | 0.050 | False | 0.064 | 0.632 | 2.11 |
| July-selected anchor | 0.025 | False | 0.045 | 0.613 | 2.04 |
| Bound center | 0.025 | False | 0.084 | 0.670 | 2.25 |
| High-kn anchor | 0.025 | False | 0.082 | 0.655 | 2.19 |

full rank와 양호한 condition에도 Dsn 신호가 일관되게 약하므로, 수치적 2변수 optimum을 물성값으로 채택하지 않았다.

## kn 단일 모델 식별도와 최적화

- 1/2.5/5% normalized perturbation에서 kn RMS sensitivity는 104.84--105.65 mV이며 모두 통과했다.
- 최선 후보는 `Bounded scalar -> TRF`이고 full-time feasibility를 만족했다.

| 후보 | charge RMSE | kn | feasible |
|---|---:|---:|---:|
| Bounded scalar -> TRF | 57.801 mV | 5.348413e-07 | True |
| DE -> TRF | 57.801 mV | 5.348451e-07 | True |
| TRF July-kn start | 57.802 mV | 5.348542e-07 | True |

## Dsn 고정값 profile

| Fixed Dsn | fitted kn | charge RMSE | discharge RMSE | discharge time feasible | charge cap RMSE | discharge cap RMSE |
|---:|---:|---:|---:|---:|---:|---:|
| 1.0000e-14 | 5.5950e-07 | 56.48 mV | invalid (early cutoff) | False | 14.97% | 0.52% |
| 4.1708e-14 | 5.3484e-07 | 57.80 mV | 75.80 mV | True | 15.82% | 2.09% |
| 7.4662e-14 | 5.3174e-07 | 58.02 mV | 77.14 mV | True | 15.90% | 2.35% |
| 8.0000e-14 | 5.3147e-07 | 58.03 mV | 77.25 mV | True | 15.91% | 2.37% |

profile은 Dsn을 데이터에서 추정한 것이 아니라 고정 가정 변화에 대한 구조 민감도다. 최종값은 independent GITT anchor를 사용한다.

## 목적함수

- charge CC 0.5C/1C/2C, 2셀, 각 100개 physical-time voltage point의 ordinary RMSE
- capacity, CV, voltage slope/shape, prior는 목적함수에서 제외
- candidate가 실험 cutoff time 전에 종료되면 fixed-time residual을 정의할 수 없어 infeasible 처리
- discharge와 capacity는 fitting 이후에만 평가

## 해석 제한

- 이번 9월 charge만으로 Dsn을 재추정하지 못했다. 이는 실패가 아니라 파라미터 수를 줄여야 한다는 식별도 결과다.
- 선택 kn의 최소 charge time margin은 0.0021초로 feasibility 경계에 매우 가깝다. 따라서 kn도 독립 kinetic 물성값보다 cutoff를 포함한 effective parameter로 해석한다.
- 2C charge capacity 오차는 두 셀에서 +25.0/+29.3%이고, discharge voltage bias는 전 branch에서 +54--66 mV다. 단일 kn 조정만으로 고율 charge polarization과 charge/discharge 비대칭을 동시에 설명하지 못한다.
- 9월 전용 QLi/window를 얻으려면 같은 셀의 저율 full qOCV 또는 전극 half-cell OCP가 필요하다.
- 2C charge는 CV 비중이 약 52.5%이므로 CC cutoff capacity를 inventory fitting 신호로 사용하지 않았다.
