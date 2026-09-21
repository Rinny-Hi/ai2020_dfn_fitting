# C/50 선정 OCP의 0.5C/1C/2C 검증

## 결론

`Ai2020 nominal 음극 + old 양극 GITT` OCP는 기존 현재 OCP보다 동적 검증 결과가 좋다.

| 방향 | 기존 평균 10–70% RMSE | C/50 OCP 평균 10–70% RMSE | 변화 |
|---|---:|---:|---:|
| Charge | 19.95 mV | 14.49 mV | -27.4% |
| Discharge | 26.13 mV | 8.46 mV | -67.6% |

충전 평균 절대 용량오차도 126.7 mAh에서 75.6 mAh로 감소했다. 방전 평균 절대 용량오차는 14.6 mAh에서 20.7 mAh로 증가했지만, 세 조건 모두 절대값 28 mAh 이내이고 capacity RMSE는 약 1%이다.

현재 자료에서는 **C/50 OCP + endpoint start**를 우선 적용한다. 기존 10분 rest 전압을 평형 SOC로 역산한 `rest/history start`는 충전 용량오차를 더 줄이지만 중앙 전압 RMSE를 악화시키므로 채택하지 않는다.

## 조건별 결과

| C-rate | 방향 | 기존 RMSE 10–70% | C/50 OCP RMSE 10–70% | C/50 OCP MAE 10–70% | 용량오차 |
|---:|---|---:|---:|---:|---:|
| 0.5C | Charge | 14.92 mV | 17.25 mV | 15.16 mV | -84.0 mAh |
| 1C | Charge | 19.49 mV | 11.74 mV | 9.72 mV | -75.3 mAh |
| 2C | Charge | 25.44 mV | 14.47 mV | 9.82 mV | -67.4 mAh |
| 0.5C | Discharge | 23.60 mV | 7.15 mV | 4.88 mV | +27.6 mAh |
| 1C | Discharge | 27.52 mV | 10.87 mV | 8.58 mV | +22.3 mAh |
| 2C | Discharge | 27.28 mV | 7.35 mV | 6.62 mV | -12.2 mAh |

0.5C charge 중앙 RMSE는 기존보다 2.33 mV 증가했지만 1C/2C charge 및 세 discharge 조건이 크게 개선됐다. 따라서 전체 C-rate 일반화 성능은 C/50 선정 OCP가 우수하다.

## Full-range 결과

| 방향 | 기존 평균 full RMSE | C/50 OCP 평균 full RMSE |
|---|---:|---:|
| Charge | 52.02 mV | 47.44 mV |
| Discharge | 32.30 mV | 30.24 mV |

full-range 오차가 중앙 구간만큼 감소하지 않는 이유는 시작전압 오차가 여전히 크기 때문이다. 기존 BoL branch 직전 rest가 10분에 불과하여 평형 초기상태를 제공하지 못한다. 새 BoL 시험의 60분 rest 결과가 확보되면 초기상태와 OCP 오차를 다시 분리해야 한다.

## 고정한 동적 파라미터

- Rn/Rp: 5/3 µm
- Dsn/Dsp: 2.1e-14 / 4.4e-14 m²/s
- kn/kp prefactor: 7.40e-7 / 3.12e-7
- Bruggeman n/p/s: 2.914 / 1.83 / 1.5
- Ln/Lp: 72.5/61.5 µm
- OCP window: x0/x100 = 0.003806/0.954221, y100/y0 = 0.431642/0.970746

동적 파라미터를 재최적화하지 않고 OCP만 교체했으므로, 개선은 OCP/window 선택의 효과로 해석할 수 있다.

## 판단

1. 현재 OCP는 `nominal 음극 + old 양극 GITT`로 변경한다.
2. 초기상태는 기존 10분 rest 전압 역산보다 endpoint start를 우선 사용한다.
3. 새 BoL 결과가 들어오면 60분 rest를 포함한 연속 프로토콜로 다시 검증한다.
4. 그 다음 charge 용량 부족을 kn 또는 usable capacity로 분리해 보정한다.
