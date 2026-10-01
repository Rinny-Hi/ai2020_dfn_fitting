# Enertech old/new RPT 및 pre-fit 후보 종합 비교

## 최종 판단

1. September RPT는 120분 rest와 두 셀 반복시험을 사용하므로 초기상태 정의와 시험 재현성 측면에서 더 좋은 동특성 데이터다.
2. July BoL과 September RPT는 동일한 parameterization target으로 바로 합치면 안 된다. September 셀은 July BoL 대비 0.5C 기준 약 2.5%, 2C 기준 약 3.6~5.4% 낮은 usable capacity와 더 큰 충·방전 전압 분리를 보인다.
3. 현재 GitHub 최신 평가는 절대 전달용량 축을 사용하므로 동특성 정규화 오류는 없다. 각 branch를 0~1로 정규화하는 방식은 형상 진단에만 사용해야 하고 최종 fitting 목적함수로 사용하면 안 된다.
4. 16개 고정 후보 중 두 데이터군에서 전압오차와 용량오차가 모두 작은 단일 후보는 없었다. 낮은 `k`와 gauge thickness를 사용하면 September 곡선은 좋아지지만 July BoL 용량오차가 악화한다. 이는 `k`가 cohort capacity/contact resistance 차이를 흡수하고 있다는 신호다.
5. 최종 fitting은 shared intrinsic parameter와 cohort-specific state parameter를 분리해야 한다. OCP shape, measured R, D는 공유하고 `Q_Li/window`와 contact resistance는 cohort별로 둔다.

## 데이터 품질 비교

| 항목 | July BoL | September RPT | 판단 |
|---|---:|---:|---|
| 셀 수 | 3 | 2 | July가 반복수 우세 |
| branch 전 rest | 10분 | 120분 | September 우세 |
| 방전 후 rest 말단 기울기 | +2.47 / +5.69 / +6.38 mV/min | +0.30 / +0.20 / +0.07 mV/min | September가 초기 평형에 훨씬 가까움 |
| 충전 후 rest 말단 기울기 | 약 -0.31~-0.34 mV/min | 약 -0.01~-0.02 mV/min | September 우세 |
| 방전 용량 셀간 범위 | 0.65 / 2.12 / 4.62 mAh | 15.83 / 13.48 / 4.22 mAh | 둘 다 양호, July 저율 우세 |
| 충전 CC+CV/방전 일치 | 양호 | 1~4 mAh 수준 | September 내부 일관성 매우 양호 |

September 데이터의 2시간 rest에도 0.5C 방전 후 저전압 rest tail이 완전히 평평하지 않으므로 엄밀한 equilibrium OCV로 보지는 않는다. 다만 기존 10분 rest보다 초기 SOC 정의에는 훨씬 유리하다.

## old-new capacity 차이

| C-rate | 방향 | July | September | September-July |
|---:|---|---:|---:|---:|
| 0.5C | CC charge | 2.1459 Ah | 2.0358 Ah | -110.1 mAh (-5.13%) |
| 0.5C | CC+CV charge | 2.3145 Ah | 2.2472 Ah | -67.3 mAh (-2.91%) |
| 0.5C | discharge | 2.2899 Ah | 2.2438 Ah | -46.1 mAh (-2.02%) |
| 1C | CC+CV charge | 2.2886 Ah | 2.2155 Ah | -73.0 mAh (-3.19%) |
| 1C | discharge | 2.2613 Ah | 2.2146 Ah | -46.7 mAh (-2.06%) |
| 2C | CC+CV charge | 2.2637 Ah | 2.1406 Ah | -123.0 mAh (-5.43%) |
| 2C | discharge | 2.2160 Ah | 2.1366 Ah | -79.4 mAh (-3.58%) |

0.5C에서 charge total과 discharge를 평균해 cohort capacity ratio를 구하면 September/July=0.97535다. 기존 C/50 `Qcell=2.35653 Ah`에 이 비율을 적용한 September의 사전 eSOH 용량 추정값은 약 2.298 Ah이다. 이 값은 새로운 C/50 실측이 아니라 cohort-specific prior로만 사용해야 한다.

## 시간-전압 형상 차이

각 branch를 자체 진행률 0~1로 정규화해도 September는 중앙부에서 July보다 다음과 같이 이동했다.

| C-rate | charge bias | discharge bias | 해석 |
|---:|---:|---:|---|
| 0.5C | +21.2 mV | -41.3 mV | September 전압 분리 증가 |
| 1C | +35.8 mV | -51.3 mV | 분극 증가 |
| 2C | +52.5 mV | -63.4 mV | 분극 증가 |

따라서 old-new 차이는 branch 용량 정규화만으로 사라지지 않는다. 충전은 위로, 방전은 아래로 이동하므로 추가 series/contact resistance, hysteresis 또는 cohort 상태 차이가 존재한다. 중간 SOC 차이를 전류로 나눈 apparent added resistance는 대략 12~30 mOhm 범위다.

## 정규화 검토

- OCP/qOCV: C/50 branch별 용량 정규화와 절대 Ah 재구성 차이는 qOCV RMSE 0.021 mV였다. C/50 용량 반복성이 높아 현재 OCP 선택에는 문제가 없다.
- dynamic voltage: 최신 `priority.curve_metrics`는 실험과 모델의 원래 `Q_Ah`를 사용하고 공통 절대용량 구간에서 비교한다. 현재 방식이 맞다.
- normalized progress: 서로 다른 용량의 전압 형상을 진단하는 보조 그림으로만 사용한다.
- capacity: CC endpoint와 CC+CV total을 분리한다. CC-only DFN과 CC+CV 실험 총량을 직접 비교하지 않는다.
- branch weighting: 시간 표본 수를 그대로 합치면 0.5C가 과대 가중된다. 각 branch를 동일한 수의 Ah grid로 재표본하고 branch별 metric을 같은 가중치로 평균한다.

## 현대차 셀과 비교

현대차 노트북의 좋은 점은 0.05C pseudo-OCV로 window를 먼저 정하고 0.5C/1C/2C를 따로 검증하며, 직전 rest 전압과 실제 전류를 사용하는 구조다. 이 구조는 유지할 가치가 있다.

다만 현대차 결과는 branch별 0~1 정규화와 overlap voltage error를 사용한다. 저장 결과의 CC-only 평균 전압 MAE는 약 14.63 mV, 용량오차는 2.63%이며 CC+CV 총량오차 0.48%는 CC-only DFN 검증과 직접 비교할 수 없다. 현재 Enertech old cohort는 선택 후보에 따라 중앙 MAE 약 16~23 mV와 capacity RMSE 약 1.3~3%까지 가능해 현대차와 크게 동떨어지지 않는다. September cohort는 가장 좋은 수치 후보에서도 중앙 MAE 약 30.9 mV, capacity RMSE 약 5.0%로 명확히 더 어렵다.

즉 현대차 코드가 특별히 더 정확해서라기보다 데이터군의 OCV-capacity 일관성과 분극 크기가 더 유리했고, metric도 일부 더 낙관적이다.

## 16개 pre-fit 후보 결과

주요 후보만 요약한다. 모든 후보는 current OCP, `Dsn=2.1e-14`, `Dsp=4.4e-14`, nominal Bruggeman, negative hysteresis off를 사용했다.

| 후보 | 핵심 설정 | July center RMSE | Sep center RMSE | July cap RMSE | Sep cap RMSE | 판단 |
|---|---|---:|---:|---:|---:|---|
| A | Ai L, SEM R, Ai kn, kp=4.18e-7 | 25.15 mV | 71.29 mV | 1.26% | 13.49% | BoL 물리 baseline |
| C | gauge L, SEM R, Ai kn, kp=4.18e-7 | 20.18 | 55.03 | 5.19% | 8.17% | current-state 물리 baseline |
| D | gauge L, SEM R, kn=7.4e-7, kp=4.18e-7 | 20.81 | 45.35 | 6.40% | 6.63% | measured kp 유지 numerical start |
| O | gauge L, SEM R, Ai kn, kp=3.12e-7 | 21.21 | 42.56 | 6.52% | 6.52% | legacy kp가 결손분극 흡수 가능 |
| P | gauge L, SEM R, kn=7.4e-7, kp=3.12e-7 | 26.07 | 33.11 | 7.83% | 5.00% | 전압상 최선이나 물성 anchor 약함 |

P가 전체 전압에서 가장 좋아 보여도 최종 물성값으로 선택하면 안 된다. 음극 Rct 기반 `kn=7.4e-7`은 이미 폐기되었고, `kp=3.12e-7`도 재계산값 4.18e-7보다 낮다. 두 kinetic 값을 낮춰 September의 추가 저항과 cohort capacity 차이를 대신 표현한 결과다.

## 선택 가능한 세 방향

### 1. BoL parameterization

- fitting: July BoL 3셀 평균
- validation: cell leave-one-out + September external/current-state validation
- 시작값: 후보 A
- 고정: current OCP, SEM R, Dsn/Dsp, kp=4.18e-7, nominal Bruggeman
- fitting 포함: contact resistance, kn
- 이후 조건부: brugg_n 또는 tau_n 중 하나

이 방향은 논문/Ai2020과 일관된 BoL 모델이 목적일 때 가장 적절하다.

### 2. 현재 September 셀 직접 parameterization

- 시작 geometry: gauge `Ln=73 um`, `Lp=62.5 um`, SEM `Rn=3.7542 um`, `Rp=4.2387 um`
- OCP shape: 현재 nominal anode + old cathode 유지
- cohort eSOH prior: `Qcell approximately 2.298 Ah`
- 첫 단계: `Q_Li/window`와 contact resistance만 식별
- 두 번째: kn fitting, kp=4.18e-7 고정
- validation: charge fitting -> discharge validation, 6-1/6-2 leave-one-cell-out

이 방향이 이번 새 곡선을 가장 잘 설명할 가능성이 높다. 단순히 후보 P를 채택하는 것보다 물리적이다.

### 3. old-new 공동 모델

- shared: OCP shape, Ln/Lp uncertainty case, SEM R, Dsn/Dsp, kp
- cohort-specific: QLi/stoichiometry window, contact resistance, 필요 시 hysteresis initial state
- objective: cohort/branch별 voltage RMSE와 capacity error를 별도 Pareto metric으로 유지

가장 일반화 가능하지만 parameter 수가 늘어 identifiability 검사가 필수다.

## 권장 fitting 순서

1. OCP/window stage
   - current OCP source를 고정한다.
   - July `Qcell=2.35653 Ah`, September prior 약 2.298 Ah를 분리한다.
   - qOCV RMSE와 endpoint 오차를 hard constraint로 둔다.
2. ohmic stage
   - initial voltage jump와 HPPC로 cell-level contact/series resistance를 식별한다.
   - September에 필요한 추가 저항은 우선 0~25 mOhm 범위에서 진단한다.
3. kinetic stage
   - kp=4.18e-7 고정 후 kn만 fitting한다.
   - kn 진단범위는 약 0.5e-6~1.3e-6으로 제한하고 bound 도달 여부를 확인한다.
4. electrolyte/solid transport stage
   - C-rate 잔차가 남으면 brugg_n 하나를 우선한다.
   - Rn과 Dsn을 따로 fitting하지 말고 `tau_n=Rn^2/Dsn`으로 profile한다.
   - Dsp는 기존 민감도가 매우 낮아 초기 fitting에서 제외한다.
5. hysteresis structure
   - positive hysteresis on/off를 discrete model comparison으로 먼저 결정한다.
   - negative hysteresis는 계속 off로 둔다.

`kn`, `kp`, `brugg_n`, `brugg_p`, `brugg_s`, `Dsn`, `Dsp`를 한 번에 fitting하면 안 된다. 기존 민감도 분석에서 `kn-kp` correlation은 약 0.998, `brugg_p-brugg_s`는 약 0.99였고 charge sensitivity matrix condition number는 약 300이었다.

## 최종 목적함수와 validation

가중합 하나로 합치지 않고 아래를 따로 보고한다.

- full-range voltage RMSE/MAE
- 10~70% 절대 Ah 구간 RMSE/MAE
- C-rate/방향별 CC capacity error
- CC+CV total capacity reconciliation
- initial jump error
- cutoff 인접 tail error

validation은 다음 두 축을 동시에 사용한다.

1. direction split: charge fitting -> discharge validation
2. cell split: 6-1 fitting -> 6-2 validation, 이후 반대 방향 교차검증

마지막에는 각 branch 독립 초기화 결과와 전체 RPT sequential simulation 결과를 모두 보고한다.

## 데이터별 역할

| 데이터 | 권장 역할 |
|---|---|
| 260610 C/50/저율 | OCP, qOCV, eSOH/window |
| 260703 BoL 6-5/6-6/6-8 | BoL C-rate fitting/validation |
| 260927 6-1/6-2 | current-state fitting 또는 external validation |
| HPPC | R0/contact와 kinetics 분리 |
| GITT | OCP/Ds prior, 10분 pulse finite-diffusion 재검증 필요 |
| cathode EIS | kp prior 4.18e-7 |
| anode EIS | washing/reproducibility 문제로 kn anchor에서 제외 |
| 260210/260610 장기 cycle 파일 | degradation validation, 단일 BoL fitting target에서 제외 |

## 산출물

- `old_vs_new_time_voltage.png`
- `old_vs_new_capacity_voltage.png`
- `cross_cohort_candidate_tradeoff.png`
- `branch_qc_all_cells.csv`
- `cohort_branch_summary.csv`
- `old_vs_new_capacity_shift.csv`
- `old_vs_new_normalized_shape_difference.csv`
- `candidate_detail.csv`
- `candidate_by_cohort.csv`
- `candidate_robust_comparison.csv`
- `analysis_manifest.json`
