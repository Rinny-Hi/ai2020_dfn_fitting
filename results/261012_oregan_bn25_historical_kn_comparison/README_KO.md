# 과거 k_n 후보 동일 조건 비교

## 비교 조건

- ORegan2022 electrolyte transport
- `b_n=2.5`
- 현재 geometry, OCP, stoichiometry window, `Dsn`, `Dsp`, `kp` 고정
- cell 6-6의 0.5C/1C/2C 충·방전 및 HPPC 사용
- September 전용 fitting 값은 제외

## 후보의 성격

| 후보 | PyBaMM prefactor | 표기용 k_n | 근거와 현재 상태 |
|---|---:|---:|---|
| EIS N3+N4 | 2.48e-7 | 2.570e-12 | 과거 음극 half-cell Rct/EIS 해석값, 재현성 문제로 폐기 |
| Past 2.59 check | 2.59e-7 | 2.684e-12 | 과거 대조검증 값, 실험 출처는 현재 기록에서 확인되지 않음 |
| EIS N4-only | 7.40e-7 | 7.670e-12 | 과거 음극 half-cell Rct/EIS 해석값, 재현성 문제로 폐기 |
| Ai2020 anchor | 9.6485e-7 | 1.000e-11 | 현재 기준값 |
| Current combined fit | 1.0606e-6 | 1.099e-11 | 현재 Dynamic+HPPC 진단 fitting, 미채택 |
| July prior fit | 1.4502e-6 | 1.503e-11 | 과거 모델 구성의 effective fitting 값, 실험값 아님 |

## 동일 조건 결과

| 후보 | Dynamic fit-window RMSE | Dynamic full RMSE | 용량 RMSE | HPPC pulse RMSE | Dynamic+HPPC 동일가중 RMSE |
|---|---:|---:|---:|---:|---:|
| EIS N3+N4 | 42.10 mV | 58.79 mV | 9.41% | 87.76 mV | 68.83 mV |
| Past 2.59 check | 40.25 mV | 57.12 mV | 9.06% | 85.65 mV | 66.92 mV |
| EIS N4-only | **19.91 mV** | **36.25 mV** | 2.28% | 38.79 mV | 30.83 mV |
| Ai2020 anchor | 26.28 mV | 37.60 mV | 1.02% | 29.39 mV | 27.88 mV |
| Current combined fit | 28.87 mV | 38.63 mV | **0.71%** | 26.45 mV | **27.69 mV** |
| July prior fit | 37.50 mV | 43.58 mV | 1.25% | **18.79 mV** | 29.66 mV |

## 해석

`7.40e-7`은 Dynamic 전압만 보면 가장 좋지만, Ai2020 기준보다 용량 RMSE가 1.02%에서 2.28%로 증가하고 HPPC pulse RMSE도 29.39 mV에서 38.79 mV로 악화된다. 충전 음극 최저전위도 1C에서 -6.01 mV에서 -14.60 mV, 2C에서 -53.93 mV에서 -61.97 mV로 더 낮아진다.

`2.48e-7`과 `2.59e-7`은 Dynamic, 용량, HPPC 및 음극전위를 모두 크게 악화시켜 현재 조건에서는 배제할 수 있다. 특히 `2.48e-7`의 용량 RMSE는 9.41%이며 2C 충전 최저 음극전위는 -96.06 mV다.

따라서 과거 실험 해석값 중에서는 `7.40e-7`이 상대적으로 낫지만, Dynamic 한 데이터셋의 보상값에 가깝고 HPPC 및 음극전위까지 포함하면 채택 근거가 부족하다. 음극 half-cell Rct 자체도 재현성 문제로 이미 폐기됐으므로 현재 최종값은 Ai2020 anchor `1.0e-11`을 유지하는 것이 가장 방어 가능하다.
