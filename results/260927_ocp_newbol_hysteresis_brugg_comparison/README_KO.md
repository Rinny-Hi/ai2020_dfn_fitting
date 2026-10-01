# OCP 선택, 2시간 rest, hysteresis + Bruggeman 비교

## 결론

1. OCP는 **Ai2020 nominal graphite 음극 + 기존 측정 cathode 평형 OCP**를 유지한다. 유효 후보 중 qOCV RMSE가 6.24 mV이고 C/50 branch MAE가 6.80 mV다. 기존/기존 조합은 8.22/9.08 mV로 더 나쁘다.
2. nominal/new cathode 조합의 qOCV RMSE는 5.03 mV로 더 작지만, new cathode는 reverse branch가 73.6%만 완성되어 absolute/directional anchor가 없다. 그래서 수치만 보고 선택하지 않았다.
3. September BoL의 120분 rest는 July 10분 rest보다 초기 SOC/OCP 정의를 **확실히 개선**한다. 다만 0.5C charge 직전 저전압 rest에는 약 0.30 mV/min의 잔류 relaxation이 있어 완전 평형으로 간주하면 안 된다.
4. `brugg_n=1.5`는 negative electrode electrolyte Bruggeman exponent만 2.914에서 1.5로 바꾼 경우다. positive=1.83, separator=1.5, fitted Dsn/kn은 고정했다. nominal anode에는 방향별 hysteresis 데이터가 없으므로 hysteresis 비교는 측정 방향성이 있는 cathode에만 적용했다.

## 1. OCP 선택 근거

- 음극: Ai2020 nominal은 절대 stoichiometry anchor가 있고 선택된 full-cell window와 직접 결합할 수 있다. 새 음극 측정치는 양방향 형상은 있으나 절대 anchor가 없어 주 OCP로 쓰지 않았다.
- 양극: 기존 측정 cathode는 절대 anchor와 lithiation/delithiation 양방향이 모두 있다. 새 cathode는 역방향 미완성이라 diagnostic 후보로만 두었다.
- 최종 선택은 단일 qOCV 오차가 아니라 C/50 charge/discharge branch 재현성과 데이터 완결성을 함께 본 결과다.

## 2. 2시간 rest가 초기 상태 문제를 얼마나 줄였나

| Rate | Direction | July slope | September slope | 감소 배수 | September inferred SOC |
|---:|---|---:|---:|---:|---:|
| 0.5C | Charge | 2.467 mV/min | 0.298 mV/min | 8.3x | 1.67% |
| 0.5C | Discharge | 0.322 mV/min | 0.009 mV/min | 34.4x | 97.38% |
| 1C | Charge | 5.691 mV/min | 0.197 mV/min | 28.9x | 2.93% |
| 1C | Discharge | 0.345 mV/min | 0.012 mV/min | 29.2x | 97.48% |
| 2C | Charge | 6.381 mV/min | 0.071 mV/min | 90.0x | 6.30% |
| 2C | Discharge | 0.309 mV/min | 0.016 mV/min | 18.7x | 97.50% |

- Charge 직전 rest 기준 기울기 감소는 0.5/1/2C에서 각각 약 8/29/90배다. Discharge 직전 고전압 rest도 약 19--34배로 평탄해졌다.
- 따라서 July의 큰 initial voltage error를 OCP fitting으로 흡수할 위험은 크게 줄었다. 하지만 긴 rest가 cell inventory나 protocol history 차이를 없애지는 않는다. September는 rate별 rest 전압과 inferred SOC가 다르며, July보다 capacity도 작으므로 cohort-specific initial SOC/QLi/window가 여전히 필요하다.

## 3. Hysteresis와 brugg_n=1.5 비교

July는 원본 full voltage curve에 대해 RMSE를 계산했고, September는 현재 원본 파일이 연결되지 않아 저장된 read-only audit의 rest/current/cutoff 정보로 endpoint capacity만 비교했다. Capacity는 fitting 목적함수가 아니라 여기서도 사후 지표다.

### July full-curve 결과

| Scenario | Direction | Full RMSE | Center RMSE | Initial | Capacity RMSE | Physical-time feasible |
|---|---|---:|---:|---:|---:|---|
| Baseline | Charge | 53.69 mV | 55.64 mV | 133.68 mV | 3.11% | yes |
| Baseline | Discharge | 49.51 mV | 37.61 mV | 66.21 mV | 1.09% | yes |
| +Cathode hysteresis | Charge | 55.95 mV | 59.68 mV | 110.44 mV | 3.25% | no |
| +Cathode hysteresis | Discharge | 51.52 mV | 41.65 mV | 60.94 mV | 1.08% | yes |
| b_n = 1.5 | Charge | 82.51 mV | 89.35 mV | 122.60 mV | 14.69% | yes |
| b_n = 1.5 | Discharge | 101.18 mV | 81.51 mV | 55.20 mV | 1.88% | yes |
| Hysteresis + b_n = 1.5 | Charge | 85.45 mV | 93.33 mV | 99.35 mV | 14.71% | yes |
| Hysteresis + b_n = 1.5 | Discharge | 103.26 mV | 85.47 mV | 49.92 mV | 1.87% | yes |

### September endpoint-only 결과

| Scenario | Direction | Capacity RMSE | Mean abs capacity error | Duration RMSE |
|---|---|---:|---:|---:|
| Baseline | Charge | 24.14% | 209.1 mAh | 4.43 min |
| Baseline | Discharge | 2.83% | 58.4 mAh | 1.59 min |
| +Cathode hysteresis | Charge | 24.42% | 210.5 mAh | 4.44 min |
| +Cathode hysteresis | Discharge | 2.82% | 58.2 mAh | 1.58 min |
| b_n = 1.5 | Charge | 40.34% | 361.4 mAh | 7.74 min |
| b_n = 1.5 | Discharge | 3.70% | 72.3 mAh | 1.77 min |
| Hysteresis + b_n = 1.5 | Charge | 40.43% | 361.4 mAh | 7.72 min |
| Hysteresis + b_n = 1.5 | Discharge | 3.69% | 72.1 mAh | 1.76 min |

### 판단

- 결합 변경은 July charge full RMSE를 53.69→85.45 mV, discharge를 49.51→103.26 mV로 바꾼다.
- hysteresis는 rest 직후 offset을 일부 이동시키지만 전체 형상/중간 SOC가 같이 좋아지는지는 별개다. 양극 old-GITT hysteresis를 1:1로 적용하는 것은 transfer test이지 확정 calibration이 아니다.
- brugg_n=1.5는 electrolyte effective transport를 크게 높이는 구조적 변경이다. 기존 fitted Dsn/kn과 함께 바꾸면 파라미터 보상이 깨지므로, 성능이 좋아 보여도 새 bound 고정값으로 바로 채택하지 말고 재식별해야 한다.

## 제한

- September 원본 `E:\예린\260927 6-1/6-2 ...xlsx`가 현재 세션에서 마운트되지 않았다. 2시간 rest 진단은 프로젝트에 이미 저장된 record-level QC 결과를 사용했다.
- September의 scenario 결과는 full voltage RMSE가 아니라 cutoff endpoint 비교다. 원본이 다시 연결되면 동일한 physical-time grid로 full rerun해야 한다.
- 모든 scenario는 current selected OCP와 final conservative Dsn/kn을 고정한 국소 비교다.
