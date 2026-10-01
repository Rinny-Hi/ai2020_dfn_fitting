# 현재 P0의 0.5C/1C/2C 곡선과 fitting 대상 판단

## 현재 성능

July BoL 3셀의 CC 구간만 사용했다. 검은 선은 셀별 실험, 색 선은 fitting 전
P0 모델이다.

| C-rate | 방향 | 10-70% RMSE | 용량오차 |
|---:|---|---:|---:|
| 0.5C | Charge | 37.79 mV | -38.88 mAh (-1.81%) |
| 0.5C | Discharge | 15.98 mV | +4.24 mAh (+0.19%) |
| 1C | Charge | 29.64 mV | -28.32 mAh (-1.48%) |
| 1C | Discharge | 17.84 mV | +7.73 mAh (+0.34%) |
| 2C | Charge | 20.78 mV | -22.64 mAh (-1.51%) |
| 2C | Discharge | 9.47 mV | -10.44 mAh (-0.47%) |

## 추가 R0 진단

직전 rest 전압과 CC 첫 전압의 차이를 실험과 DFN에서 비교했다.

| 조건 | 실험 점프 | DFN 점프 | 실험-DFN | apparent extra R0 |
|---|---:|---:|---:|---:|
| 0.5C charge | 44.87 mV | 193.32 mV | -148.45 mV | -130.22 mOhm |
| 0.5C discharge | 28.50 mV | 92.59 mV | -64.09 mV | -56.22 mOhm |
| 1C charge | 75.80 mV | 251.94 mV | -176.14 mV | -77.26 mOhm |
| 1C discharge | 57.93 mV | 156.64 mV | -98.71 mV | -43.29 mOhm |
| 2C charge | 137.30 mV | 322.90 mV | -185.60 mV | -40.70 mOhm |
| 2C discharge | 113.17 mV | 238.82 mV | -125.65 mV | -27.56 mOhm |

모델이 이미 실험보다 큰 초기 분극을 계산하므로 양의 추가 직렬저항을 넣으면 오차가
더 커진다. 필요한 R0가 전 조건에서 음수라는 것은 R0가 누락된 것이 아니라 현재의
반응/수송 분극 또는 초기상태 표현이 과도하다는 뜻이다. 따라서 `Contact resistance`
는 0으로 유지하고 fitting 대상에서 제외한다.

## 다음 fitting

1. OCP/window, geometry, measured R, measured Ds, cathode-EIS kp와 Bruggeman은 고정한다.
2. 신뢰 가능한 anode Rct가 없고 Ai2020 임시값을 사용 중인 `kn` 하나만 0.5C/1C/2C
   charge에서 fitting한다.
3. 목적함수는 C-rate별 동일 가중의 voltage SSE/RMSE이며, CC 용량오차는 각 branch
   2% 이내 제약으로 사용한다.
4. 같은 kn으로 discharge를 validation한다.
5. 잔차가 고율/SOC 의존적으로 남을 때만 `tau_n=Rn^2/Dsn` 또는 `b_n` 중 하나를
   추가한다.

## 파일

- `p0_time_voltage.png`
- `p0_capacity_voltage.png`
- `p0_capacity_error.png`
- `p0_curve_summary.csv`
- `r0_initial_jump_summary.csv`

