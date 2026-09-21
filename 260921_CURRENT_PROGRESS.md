# 260921 DFN/OCP 재검토 진행 현황

## 1. 범위와 현재 기준

신규 GITT는 재현성과 용량 차이 검토가 끝날 때까지 현재 모델 선정에서 제외했다. 기존 저율 full-cell 데이터, 기존 half-cell OCP, 0.5C/1C/2C 충·방전 데이터를 사용해 OCP, stoichiometry window, transport/kinetic parameter의 영향을 분리했다.

현재 코드는 OCP 적합도만 최소화하지 않고 다음 두 지표를 함께 기록한다.

- 10–70% SOC 전압 MAE
- 종료시간으로 환산한 용량오차

## 2. OCP와 stoichiometry window

현재 OCP 조합은 **Ai2020 nominal equilibrium + harvested-electrode hysteresis offset**이다.

- 충전: 음극 lithiation / 양극 delithiation
- 방전: 음극 delithiation / 양극 lithiation
- `Qn = 2.5651 Ah`, `Qp = 4.3010 Ah`
- `x0 = 0.003948`, `x100 = 0.915316`
- `y100 = 0.429778`, `y0 = 0.973308`
- qOCV MAE: **4.64 mV**
- qOCV RMSE: **5.59 mV**
- endpoint 최대 절대오차: **2.01 mV**

이 window는 punched-electrode 질량에 직접 의존하지 않는 effective-capacity fit이다. 따라서 이후 teardown 형상을 반영할 때는 이미 맞춘 `Qn/Qp`를 보존하도록 `c_s,max`를 함께 환산했다.

## 3. 현재 반영 파라미터

| 구분 | 현재값 | 근거/상태 |
|---|---:|---|
| `Ln` | 72.5 µm | 음극 양면 총두께 155 µm, Cu 10 µm 가정의 잠정값 |
| `Lp` | 61.5 µm | 양극 양면 총두께 138 µm, Al 15 µm 가정의 잠정값 |
| `Rn` | 5 µm | Ai2020 nominal |
| `Rp` | 3 µm | Ai2020 nominal |
| `Dsn` | 2.1e-14 m²/s | 지정값 |
| `Dsp` | 4.4e-14 m²/s | 지정값 |
| `kn` prefactor | 7.40e-7 | EIS N4-only 값 |
| `kp` prefactor | 3.12e-7 | EIS 1차값 |
| `bn` | 2.914 | 사용자 지시에 따라 Ai2020 nominal로 복원 |
| `bp` | 1.83 | Ai2020 nominal |
| `bs` | 1.5 | Ai2020 nominal |
| electrolyte `De` | raw Ai2020 함수 ×1e-4 | 단위 환산 유지 |
| effective `csn,max` | 26,553.8 mol/m³ | `Ln` 변경 후 `Qn` 보존 환산값 |
| effective `csp,max` | 51,641.4 mol/m³ | `Lp` 변경 후 `Qp` 보존 환산값 |

`c_s,max` 환산값은 소재 고유 최대농도가 새로 측정되었다는 뜻이 아니라, OCP에서 식별한 active-material inventory를 새로운 두께 표현에서도 유지하기 위한 effective parameter이다.

## 4. EIS kinetic parameter 검토

PyBaMM Ai2020/Dualfoil 교환전류 식은 `m_ref = 1e-11 F = 9.6485e-7`을 사용한다. 첨부 EIS 값을 같은 prefactor 단위로 해석해 적용했다.

| 조합 | 평균 전압 MAE | 평균 절대 용량오차 |
|---|---:|---:|
| Ai2020 `kn`, `kp` | 31.24 mV | 30.92 mAh |
| `kn=2.48e-7`만 | 31.94 mV | 106.40 mAh |
| `kp=3.12e-7`만 | 15.47 mV | 62.15 mAh |
| `kn=2.48e-7`, `kp=3.12e-7` | 71.24 mV | 179.89 mAh |
| `kn=7.40e-7`만 | 21.05 mV | 32.56 mAh |
| `kn=7.40e-7`, `kp=3.12e-7` | 22.17 mV | 79.69 mAh |

현재는 사용자의 요청에 따라 `kn=7.40e-7`, `kp=3.12e-7`을 모두 실험값으로 유지한다.

## 5. Bruggeman 탐색과 최종 선택

두 EIS kinetic 값을 고정하고 `bn/bp`를 격자 탐색했다. 수치상 Pareto 후보는 `bn=2.5` 부근이었다.

| `bn` | `bp` | `bs` | 평균 전압 MAE | 평균 절대 용량오차 |
|---:|---:|---:|---:|---:|
| 2.5 | 1.50 | 1.5 | 15.40 mV | 40.98 mAh |
| 2.5 | 1.83 | 1.5 | 13.67 mV | 43.00 mAh |
| 2.5 | 2.00 | 1.5 | 13.41 mV | 49.40 mAh |

하지만 Bruggeman을 단순 fitting parameter로 확정하지 않고, **현재 최종 구성은 사용자 지시에 따라 nominal `bn=2.914`, `bp=1.83`, `bs=1.5`**를 사용한다.

## 6. 잠정 전극 두께 반영

단면 코팅 두께는 다음과 같이 계산했다.

```text
Ln = (155 - assumed Cu 10) / 2 = 72.5 µm
Lp = (138 - assumed Al 15) / 2 = 61.5 µm
```

집전체 실측 전까지 잠정값이다. 두께만 변경해 active-material inventory가 바뀌지 않도록 `c_s,max`를 두께에 반비례시켰다.

`bn=2.5` 탐색 상태에서 비교했을 때 nominal 두께 대비 전압 MAE는 13.67→14.89 mV로 소폭 증가했고 평균 절대 용량오차는 43.00→39.95 mAh로 소폭 감소했다. 즉 두께 반영의 동적 영향은 작았으며, 형상 근거를 우선해 잠정값을 현재 구성에 유지했다.

## 7. 현재 최종 동적 검증 결과

현재 설정은 실험 `kn/kp`, 잠정 `Ln/Lp`, nominal Bruggeman을 동시에 사용한다.

| C-rate | 방향 | 10–70% SOC MAE | 용량오차 |
|---:|---|---:|---:|
| 0.5C | Charge | 11.70 mV | -138.94 mAh |
| 0.5C | Discharge | 16.98 mV | +6.89 mAh |
| 1C | Charge | 15.36 mV | -123.00 mAh |
| 1C | Discharge | 18.07 mV | +0.62 mAh |
| 2C | Charge | 21.67 mV | -117.90 mAh |
| 2C | Discharge | 26.20 mV | -36.31 mAh |

- 평균 전압 MAE: **18.33 mV**
- 충전 평균 MAE: **16.24 mV**
- 방전 평균 MAE: **20.42 mV**
- 평균 절대 용량오차: **70.61 mAh**
- 최대 절대 용량오차: **138.94 mAh**

중간 SOC 전압 형상은 개선됐지만 충전 종료 용량오차가 남아 있다. 이 오차를 Bruggeman만으로 줄이면 kinetic/transport와 capacity/endpoint 영향이 섞일 수 있으므로, 다음 단계에서는 usable capacity, stoichiometric endpoint, cutoff 부근 OCP를 분리해야 한다.

## 8. 주요 실행 파일과 결과

- 현재 설정 실행: `nominal_radius_current_configuration.py`
- EIS kinetic ablation: `eis_first_kinetics_application.py`
- Bruggeman grid: `eis_bruggeman_screen.py`
- separator Bruggeman refinement: `eis_bruggeman_separator_refine.py`
- 잠정 두께 A/B: `current_thickness_estimate_application.py`
- transport ablation: `transport_parameter_ablation.py`
- `Rn/Rp` 조합: `ocp_rn_rp_combination_screen.py`
- OCP 후보/stoichiometry 분석: `mass_independent_nominal_anode_trial.py`, `finalize_ocp_recommendation.py`

주요 산출물:

- `results/260921_nominal_radius_current_configuration/`
- `results/260921_eis_first_kinetics/`
- `results/260921_eis_bruggeman_screen/`
- `results/260921_current_thickness_estimate/`
- `results/260921_transport_parameter_ablation/`
- `results/260920_mass_independent_nominal_anode/`

## 9. 재현

```powershell
.\.venv\Scripts\python.exe nominal_radius_current_configuration.py
.\.venv\Scripts\python.exe eis_first_kinetics_application.py
.\.venv\Scripts\python.exe eis_bruggeman_screen.py
.\.venv\Scripts\python.exe current_thickness_estimate_application.py
```

집전체 두께가 측정되면 `Ln=(155-LCu)/2`, `Lp=(138-LAl)/2`로 갱신하고, 동일한 방식으로 `Qn/Qp` 보존 환산 후 재검증해야 한다.
