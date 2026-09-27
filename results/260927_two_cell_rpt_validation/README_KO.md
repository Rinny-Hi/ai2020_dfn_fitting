# 260927 two-cell 0.5C/1C/2C 데이터 검증

## 결론

- 6-1과 6-2의 반복성은 좋다. CC 용량 차이는 모든 조건에서 평균 대비 2.2% 이내이며, 방전은 0.2~0.7% 수준이다.
- 고율 충전에서 CC 용량이 작아지는 것은 셀 총용량 손실이 아니라 CV 구간으로 용량이 이동하기 때문이다. CC+CV 충전용량은 같은 C-rate의 방전용량과 약 1~4 mAh 이내로 일치한다.
- 현재 고정 확산계수와 SEM 입자반경을 사용할 때 `R^2/D`는 음극 11.19분, 양극 6.81분이다. C-rate 응답에서 농도구배가 생길 수 있는 물리적으로 가능한 크기다.
- 그러나 10분 GITT pulse의 Fourier number가 음극 0.89, 양극 1.47이므로, 10분 전체를 반무한 확산식으로 처리해 얻은 `D_s`라면 그 가정은 성립하지 않는다. 초기 약 40~70초의 sqrt(t) 선형부 또는 유한구/DFN pulse fitting으로 재검증해야 한다.
- 120분 rest는 SEM 반경 기준 `R^2/D`의 약 11~18배다. 그럼에도 방전 직후 rest 말단에 양의 전압 기울기가 남고 방향 비대칭이 크므로, 장시간 relaxation을 단일 상수 `D_s`만으로 설명하면 안 된다.

## 시험 데이터 QC

| C-rate | 방향 | 6-1 CC (Ah) | 6-2 CC (Ah) | 셀간 범위/평균 |
|---:|---|---:|---:|---:|
| 0.5C | Charge | 2.0272 | 2.0444 | 0.85% |
| 0.5C | Discharge | 2.2359 | 2.2517 | 0.71% |
| 1C | Charge | 1.7426 | 1.7637 | 1.20% |
| 1C | Discharge | 2.2079 | 2.2213 | 0.61% |
| 2C | Charge | 1.0064 | 1.0281 | 2.13% |
| 2C | Discharge | 2.1345 | 2.1387 | 0.20% |

평균 충전용량의 CC/CV 분해는 다음과 같다.

| C-rate | CC (Ah) | CV (Ah) | CC+CV (Ah) | CV 비율 | 평균 방전 (Ah) |
|---:|---:|---:|---:|---:|---:|
| 0.5C | 2.0358 | 0.2114 | 2.2472 | 9.4% | 2.2438 |
| 1C | 1.7531 | 0.4624 | 2.2155 | 20.9% | 2.2146 |
| 2C | 1.0173 | 1.1234 | 2.1406 | 52.5% | 2.1366 |

따라서 이번 데이터에서 모델 검증용 charge capacity는 CC만 비교할 때와 CC+CV를 비교할 때를 반드시 구분해야 한다. 현재 DFN 비교는 CC branch의 cutoff까지를 같은 정의로 맞췄다.

## 확산 시간척도

사용한 정의는 `tau_D = R^2 / D_s`이며, 구형입자의 가장 느린 고유모드 시간은 대략 `R^2 / (pi^2 D_s)`이다.

| 반경 | 전극 | R (um) | Ds (m2/s) | R2/D (min) | R2/(pi2 D) (min) | 10분 pulse Fo |
|---|---|---:|---:|---:|---:|---:|
| Nominal | Negative | 5.0000 | 2.1e-14 | 19.84 | 2.01 | 0.50 |
| Nominal | Positive | 3.0000 | 4.4e-14 | 3.41 | 0.35 | 2.93 |
| SEM ImageJ | Negative | 3.7542 | 2.1e-14 | 11.19 | 1.13 | 0.89 |
| SEM ImageJ | Positive | 4.2387 | 4.4e-14 | 6.81 | 0.69 | 1.47 |

SEM 반경 기준 평균 CC 지속시간과 `R^2/D`의 비는 다음과 같다.

- 0.5C charge: 음극 9.58 tau, 양극 15.74 tau
- 1C charge: 음극 4.13 tau, 양극 6.78 tau
- 2C charge: 음극 1.20 tau, 양극 1.97 tau
- 0.5C discharge: 음극 10.56 tau, 양극 17.35 tau
- 1C discharge: 음극 5.21 tau, 양극 8.56 tau
- 2C discharge: 음극 2.51 tau, 양극 4.13 tau

이는 2C 충전에서 고체 확산 제한이 눈에 띌 수 있음을 뜻한다. 다만 `tau`보다 시간이 길다는 사실이 정상상태 또는 균일 농도를 뜻하지는 않는다. 전류가 계속 흐르는 동안에는 농도구배가 지속적으로 형성된다.

## 모델 비교

두 셀, 충·방전 6개 branch 전체를 동일 조건으로 비교한 결과다.

| L case | R case | Full RMSE (mV) | 10~70% RMSE (mV) | Capacity RMSE (%) | 최대 용량오차 (mAh) |
|---|---|---:|---:|---:|---:|
| Thickness gauge | SEM ImageJ | 50.83 | 45.33 | 6.63 | 178.51 |
| Thickness gauge | Nominal | 50.84 | 46.78 | 6.66 | 178.70 |
| Ai2020 | Nominal | 66.94 | 63.01 | 12.03 | 315.96 |
| Ai2020 | SEM ImageJ | 67.43 | 61.51 | 11.88 | 311.26 |
| SEM cross-section | Nominal | 74.89 | 69.76 | 14.73 | 382.06 |
| SEM cross-section | SEM ImageJ | 75.60 | 68.24 | 14.53 | 376.05 |

수치상 두께측정기+SEM 반경이 가장 좋지만, 전체 RMSE 약 51 mV는 여전히 크다. 실험 대비 모델은 충전에서 대체로 낮고 방전에서 높아서, 전극 두께 하나로 해결되는 오차가 아니다. 두께 case가 낮은 porosity를 통해 미반영 분극을 대신 흡수했을 가능성이 있으므로 이를 최종 geometry 선정 근거로 쓰면 안 된다.

최선 case에서도 가장 큰 실패는 2C charge이며 셀별 CC 용량오차는 +178.5 mAh, +141.8 mAh다. 반면 2C discharge는 -59.1, -61.3 mAh로 재현성이 높다. 이는 데이터 노이즈보다는 방향 의존 분극, OCP/hysteresis, 계면·접촉저항 또는 kinetics 쪽의 구조적 오차를 우선 의심하게 한다.

## 다음 검증 우선순위

1. `D_s` 산출에 사용한 GITT fitting window를 확인하고 초기 sqrt(t) 구간 또는 finite-sphere/DFN pulse fitting으로 재산출한다.
2. full-rest 전압에서 전류 인가 직후의 `Delta V / I`를 이용해 series/contact resistance를 분리한다. 이번 데이터는 대략 discharge 31~33 mOhm, charge 34~41 mOhm 수준의 apparent resistance를 보인다.
3. 동일 초기 SOC에서 평균 OCP/no-hysteresis와 현재 OCP/hysteresis를 다시 비교해 방향별 중간영역 편차를 확인한다.
4. `k_n`, contact resistance, electrolyte transport를 각각 ablation하되 geometry와 `D_s`를 동시에 움직이지 않는다.
5. `D_s(SOC)` 또는 particle-size distribution이 필요한지는 위 항목 이후 잔차의 SOC 의존성으로 판단한다.

## 산출물

- `new_rpt_branch_qc.csv`: branch별 용량, 시간, rest 전압/기울기
- `two_cell_repeatability.csv`: 6-1/6-2 반복성
- `diffusion_timescales.csv`: 반경·확산계수 조합별 시간척도
- `branch_duration_vs_diffusion_time.csv`: branch 시간 대비 tau
- `six_case_new_data_detail.csv`: 72개 simulation 상세 결과
- `six_case_new_data_summary.csv`, `six_case_new_data_overall.csv`: 집계 결과
- `new_rpt_two_cell_repeatability.png`: 두 셀 실험곡선 비교
- `new_rpt_pareto_candidate_curves.png`: 대표 geometry 후보와 실험 비교
