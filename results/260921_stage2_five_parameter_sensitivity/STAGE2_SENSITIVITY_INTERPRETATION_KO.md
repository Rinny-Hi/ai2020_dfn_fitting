# 2단계 민감도·식별성 분석 결과

## 분석 조건

- 고정 항목: 현재 OCP, stoichiometry window, 전극 두께/용량 보정, `Rn=5 um`, `Rp=3 um`, `De x 1e-4`
- 기준값: `kn=7.40e-7`, `kp=3.12e-7`, `Dsn=2.10e-14 m2/s`, `Dsp=4.40e-14 m2/s`, `brugg_n=2.914`
- 대상 파라미터: `kn`, `kp`, `Dsn`, `Dsp`, `brugg_n`
- 각 양의 파라미터에 대해 로그 좌표 중앙차분으로 ±5%를 적용했다.
- `brugg_n`은 지수 자체를 ±5% 변경하지 않고 실제 음극 유효 전해질 수송항 `eps_n**brugg_n`이 ±5% 변하도록 했다. 음극 공극률 0.33에서 사용된 지수는 2.870과 2.958이다.
- C-rate 전압 형상은 실험 SOC 10–70% 중 모델 cutoff 전 공통 구간에서 비교했다. cutoff 시각과 용량 변화는 별도의 endpoint 민감도로 보존했다.
- HPPC는 10개 SOC 수준의 `rest_pre_end`, `dchg_end`, `rest_30s`, `rest_end`, `chg_end` 총 50점을 사용했다.

## ±5% 변화가 만드는 실제 RMS 전압 변화

| 데이터 | brugg_n | kn | kp | Dsn | Dsp |
|---|---:|---:|---:|---:|---:|
| C-rate 충·방전 통합 | 3.47 mV | 2.11 mV | 2.11 mV | 0.33 mV | 0.06 mV |
| HPPC | 0.41 mV | 1.37 mV | 1.32 mV | 0.29 mV | 0.04 mV |

따라서 C-rate에서는 `brugg_n`의 영향이 가장 크고, HPPC pulse 끝점에서는 `kn`, `kp`가 가장 크다. `Dsp`는 두 데이터 모두에서 현재 실험 분해능으로 식별하기 어렵다.

## 조건별 해석

- `brugg_n`은 C-rate 전체에서 1순위이며 특히 2C 방전 전압 형상에서 가장 강하다. ±5% 유효수송 변화에 대한 RMS 영향은 약 7.32 mV이다.
- `kn`과 `kp`는 charge/discharge 양쪽에서 모두 크지만, 전체 C-rate 민감도 벡터의 cosine correlation이 0.998이고 HPPC에서도 0.998이다. 즉 두 값을 동시에 최적화하면 각 전극 반응속도를 독립적으로 구분하기보다 같은 전압 오차를 서로 보상할 가능성이 매우 크다.
- `Dsn`은 C-rate 통합 상대 민감도 0.096으로 경계값이지만, HPPC에서는 0.211이며 방전 pulse 후 relaxation과 2C 방전에 정보가 집중된다.
- `Dsp` 상대 민감도는 C-rate 0.017, HPPC 0.033뿐이다. 이번 데이터로 최적화하면 물리값보다 잡음·다른 파라미터 오차를 흡수할 위험이 크다.
- ±5% 중앙차분의 비선형성 비율은 모든 조합에서 0.047 이하이므로 이번 국소 선형 순위는 안정적이다.

## Endpoint 영향

2C에서 ±5% 변화가 만드는 cutoff 용량 변화의 국소 크기는 다음과 같다.

| 파라미터 | 2C Charge | 2C Discharge |
|---|---:|---:|
| brugg_n | 18.52 mAh | 3.65 mAh |
| kn | 9.53 mAh | 0.67 mAh |
| kp | 9.04 mAh | 0.69 mAh |
| Dsn | 1.53 mAh | 4.33 mAh |
| Dsp | 0.52 mAh | 0.06 mAh |

이 결과는 전압 형상만 볼 때보다 2C charge cutoff가 `brugg_n`, `kn`, `kp`에 더 민감하다는 뜻이다. 따라서 3단계 목적함수에는 공통 시간 구간 전압 RMSE와 endpoint/capacity 항을 분리해서 함께 넣어야 한다.

## 3단계에 넘길 권고 subset

1. EIS 실험값을 신뢰하는 기본안: `kn`, `kp`를 각각 현재 실험값에 고정하고 `brugg_n`, `Dsn`만 최적화한다. `Dsp`는 현재값에 고정한다.
2. EIS 값의 보정을 허용하는 대안: `kn`, `kp`를 독립 변수로 동시에 풀지 않고 공통 kinetic scale 1개와 `brugg_n`, `Dsn`을 최적화한다.
3. C-rate 데이터만 학습할 경우 자동 선별 subset은 `brugg_n + kn`이며 정규화 민감도 행렬 조건수는 2.84이다.
4. HPPC 데이터만 학습할 경우 자동 선별 subset은 `kn + brugg_n + Dsn`이며 조건수는 2.71이다. 다만 실제 pulse 내부 표본이 시작/끝 중심이라 `Dsn` 추정은 보수적으로 해석해야 한다.

최종 권고는 **우선 `kn`, `kp`를 EIS 값으로 고정하고 `brugg_n + Dsn`을 C-rate 또는 HPPC 한쪽에서 학습한 뒤 다른 쪽으로 검증**하는 것이다. 이 방식이 현재 데이터에서 가장 식별 가능하고, EIS 측정값을 DFN 전압 오차에 과도하게 재흡수시키지 않는다.

