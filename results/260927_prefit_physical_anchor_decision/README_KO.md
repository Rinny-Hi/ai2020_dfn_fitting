# July BoL fitting 전 물리 파라미터 선정

## 결론

동적 파라미터 최적화 전 공식 baseline은 `P0 anchored reference`를 권장한다.

- OCP/window: C/50에서 선정한 nominal anode + old cathode
- `Ln/Lp`: Ai2020 76.5/68 um
- `Rn/Rp`: SEM 3.7542/4.2387 um
- `Dsn/Dsp`: 소재 기반 2.1e-14/4.4e-14 m2/s
- `kn`: Ai2020 9.6485e-7
- `kp`: cathode EIS 재계산 4.18e-7
- Bruggeman n/p/s: Ai2020 2.914/1.83/1.5
- negative/positive hysteresis: 모두 off
- initial state: branch 직전 rest 전압으로 OCP 역산

이 설정은 July BoL 3셀, 0.5C/1C/2C 충방전 18개 곡선에서 fitting 없이 다음 성능을 보였다.

| 지표 | P0 |
|---|---:|
| full MAE / RMSE | 26.77 / 37.95 mV |
| 10-70% MAE / RMSE | 19.72 / 21.92 mV |
| capacity RMSE | 1.32% |
| 평균 절대 용량오차 | 18.71 mAh |
| 최대 절대 용량오차 | 53.26 mAh |

과거 notebook-like 후보는 10-70% RMSE가 15.53 mV로 더 낮지만 capacity RMSE는 2.90%이며, 폐기한 anode-Rct `kn`, legacy `kp`, endpoint initialization을 사용한다. 따라서 과거 결과의 낮은 전압오차는 독립적으로 신뢰할 수 있는 물성 개선이 아니라 오차 상쇄를 포함한다.

## 과거 코드 해석

`260818_0단계~2단계_진행_결과.ipynb`의 Stage 1은 experimental OCP와 stoichiometry window를 먼저 재식별한 결과이다. 완전한 nominal pre-fit이 아니다.

- qOCV RMSE/MAE: 10.99/9.37 mV
- 고율 평균 10-70% MAE: 22.59 mV
- charge 평균: 12.18 mV
- discharge 평균: 33.00 mV
- 전체 시간 MAE: 31.29 mV

당시 코드는 실험 시간축에서 simulation을 interpolation하고 공통 overlap만 평가했다. cutoff에 일찍 도달해 사라진 시간 및 용량 차이는 전압오차에 직접 포함되지 않았고, capacity error 표도 함께 제시하지 않았다. 따라서 charge plot은 매우 좋아 보이지만 discharge와 cutoff capacity를 합친 전체 baseline은 현재보다 우수하다고 할 수 없다.

현재 P0는 더 엄격한 절대 Ah 평가와 3셀 평균에서 center MAE 19.72 mV, capacity RMSE 1.32%이다. 같은 의미의 평균값으로 보면 과거 Stage 1의 22.59 mV보다 작다.

## 직접 측정값 적용 판단

| 파라미터 | 비교 결과 | 판단 |
|---|---|---|
| `Rn/Rp` | SEM 사용 시 center RMSE 21.92 mV, capacity RMSE 1.32%; Ai2020 5/3 um은 24.00 mV, 1.48% | SEM 적용 |
| `Dsp` | 측정 4.4e-14에서 21.92 mV, 1.32%; Ai2020 함수로 복귀하면 24.42 mV, 6.25% | 측정값 적용 |
| `Dsn` | 측정값 21.92 mV, 1.32%; Ai2020 함수 22.67 mV, 1.27% | 차이가 작으므로 측정값 유지, nominal은 uncertainty 후보 |
| `kp` | EIS 4.18e-7에서 21.92 mV, 1.32%; legacy 3.12e-7은 15.77 mV, 2.55% | EIS 적용. legacy의 낮은 전압오차는 오차 상쇄 |
| `kn` | 신뢰 가능한 anode Rct가 없으므로 Ai2020 사용 | 이후 fitting 대상 |
| `Ln/Lp` | Ai2020 21.92 mV, 1.32%; gauge 19.66 mV, 5.39%; SEM 28.37 mV, 1.54%; consensus 19.45 mV, 1.86% | 공식 baseline은 Ai2020, L uncertainty 별도 보고 |
| positive hysteresis | on에서 25.15 mV, 1.26%; off에서 21.92 mV, 1.32% | 근거가 추가되기 전 off |
| 초기화 | rest 역산 21.92 mV, 1.32%; endpoint 25.82 mV, 1.56% | rest 역산 적용 |

`Ln/Lp`를 변경할 때는 fitted `Qn/Qp`를 유지하도록 `c_s,max`를 함께 변환했다. 따라서 두께 비교는 열역학적 용량을 임의로 바꾼 비교가 아니라 transport/kinetic geometry 효과 비교이다.

## P0 C-rate별 결과

| 조건 | center RMSE | 용량오차 |
|---|---:|---:|
| 0.5C charge | 37.79 mV | -38.88 mAh (-1.81%) |
| 0.5C discharge | 15.98 mV | +4.24 mAh (+0.19%) |
| 1C charge | 29.64 mV | -28.32 mAh (-1.48%) |
| 1C discharge | 17.84 mV | +7.73 mAh (+0.34%) |
| 2C charge | 20.78 mV | -22.64 mAh (-1.51%) |
| 2C discharge | 9.47 mV | -10.44 mAh (-0.47%) |

남은 구조적 오차는 0.5C/1C charge 중앙 형상에 집중되어 있다. 전체 용량을 바꾸기보다 다음 단계에서 contact/series resistance를 먼저 분리하고, 그 후 `kn` 하나만 charge fitting하여 discharge를 validation하는 순서가 적절하다.

## 최적화 전 제안 후보

### 공식 baseline: P0

물리적 근거, 전압 형상, 용량오차의 균형이 가장 좋다. 발표자료와 fitting의 시작점으로 사용한다.

### 두께 uncertainty 보조 후보: L2

`Ln/Lp=76.5/62.7325 um`에서 center RMSE 19.45 mV, capacity RMSE 1.86%이다. 전압 형상은 P0보다 2.46 mV 좋아지지만 용량오차가 증가하므로 최종값이 아니라 두께 불확실성 범위로만 보고한다.

### 수치상 하한: K0/OLD

legacy `kp` 또는 폐기된 `kn`을 쓰면 center RMSE 15-16 mV까지 낮아진다. 그러나 capacity RMSE가 2.55-2.90%로 커지고 독립 물성 근거가 약하므로 최종 파라미터로 채택하지 않는다. 이 값은 향후 fitting으로 얻을 수 있는 전압오차의 대략적 수치 하한으로만 사용한다.

## 산출물

- `prefit_anchor_summary.csv`: 모든 고정 후보 종합 지표
- `prefit_anchor_detail.csv`: 셀/C-rate/방향별 상세 지표
- `selected_branch_summary.csv`: P0, endpoint, old-notebook-like 비교
- `prefit_anchor_tradeoff.png`: 전압-용량 Pareto 비교
- `selected_prefit_capacity_voltage.png`: 절대 Ah 축 곡선 비교
- `selected_prefit_time_voltage.png`: 시간-전압 곡선 비교
- `analysis_manifest.json`: 고정값과 평가 조건
