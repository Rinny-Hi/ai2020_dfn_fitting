# ORegan2022, b_n=2.5, k_n 단독 fitting 결과

## 고정 조건

- 전해질 transport: ORegan2022
- 음극 Bruggeman 계수: `b_n = 2.5`
- 고정값: `Dsn = 4.170774e-14 m2/s`, `Dsp = 6.888175e-14 m2/s`
- 고정값: `kp = 5.899520e-7` (PyBaMM prefactor)
- fitting 변수: 음극 반응속도 prefactor `kn` 한 개
- 기준 `kn = 9.648533e-7` (표기용 반응속도계수 `1.000000e-11`)

## 목적함수

- Dynamic: cell 6-6의 0.5C, 1C, 2C 충·방전 곡선에서 공통 10–70% 용량 구간 전압오차
- HPPC: 각 SOC의 휴지전압을 기준으로 한 30초 방전 pulse drop 및 10초 충전 pulse rise, 총 20점
- Combined: Dynamic MSE와 HPPC pulse MSE를 데이터셋별 동일 가중치로 결합
- HPPC 절대전압 50점은 검증용으로만 사용했다. OCV 및 초기 SOC 편향이 포함되므로 `kn` fitting 목적함수로 사용하지 않았다.

## 결과

| 구분 | kn 배율 | 반응속도계수 | Dynamic RMSE | HPPC pulse RMSE | HPPC 절대전압 RMSE |
|---|---:|---:|---:|---:|---:|
| Pre-fit | 1.000 | 1.000e-11 | 26.28 mV | 29.39 mV | 55.98 mV |
| Dynamic-only | 0.604 | 6.045e-12 | 17.38 mV | 48.36 mV | 62.01 mV |
| HPPC-pulse-only | 2.256 | 2.256e-11 | 47.59 mV | 15.11 mV | 52.48 mV |
| Combined | 1.099 | 1.099e-11 | 28.87 mV | 26.45 mV | 55.22 mV |

Dynamic-only 값은 중앙구간 전압을 개선하지만 전체곡선 RMSE는 37.60 → 37.23 mV로 거의 변하지 않고, 용량 RMSE는 1.02 → 3.59%로 악화된다. HPPC-only 값은 pulse RMSE를 크게 줄이지만 Dynamic 오차를 크게 악화시킨다. Combined 값도 목적함수는 27.88 → 27.69 mV로 0.19 mV만 개선되며, 전체곡선 RMSE는 37.60 → 38.63 mV로 악화된다.

## 판단

Dynamic과 HPPC가 서로 반대 방향의 `kn` 보정을 요구하므로 두 데이터에 공통인 `kn` 최적값은 현재 상태에서 강하게 식별되지 않는다. `kn`이 다른 파라미터 또는 OCV/초기화 오차를 대신 보정하는 현상으로 해석하는 것이 타당하다.

따라서 최종 권고는 **`kn` fitting 값을 채택하지 않고 EIS anchor인 `1.0e-11`을 유지**하는 것이다. Combined 값 `1.099e-11`은 민감도 비교용 결과로만 보관한다.

Combined 값의 음극전위 변화도 작다. 1C 충전 최저전위는 -6.01 → -3.09 mV, 2C 충전은 -53.93 → -51.11 mV로 소폭 개선되지만 0 V 미만 시간 비율은 각각 3.01 → 1.45%, 47.68 → 43.86%에 머문다.
