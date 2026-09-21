# 260921 DFN/OCP 재검토 진행 현황

## 1. 범위와 현재 기준

신규 GITT는 재현성과 용량 차이 검토가 끝날 때까지 현재 모델 선정에서 제외했다. 기존 저율 full-cell 데이터, 기존 half-cell OCP, 0.5C/1C/2C 충·방전 데이터를 사용해 OCP, stoichiometry window, transport/kinetic parameter의 영향을 분리했다.

현재 코드는 OCP 적합도만 최소화하지 않고 다음 두 지표를 함께 기록한다.

- 10–70% SOC 전압 MAE
- 종료시간으로 환산한 용량오차

## 2. OCP와 stoichiometry window

현재 OCP 조합은 **Ai2020 nominal 음극 + 기존(old) 양극 GITT**이다.

- 음극: nominal equilibrium, hysteresis 미사용
- 양극 충전: delithiation / 방전: lithiation
- `Qn = 2.47947 Ah`, `Qp = 4.37120 Ah`
- `x0 = 0.003806`, `x100 = 0.954221`
- `y100 = 0.431642`, `y0 = 0.970746`
- C/50 qOCV MAE: **4.69 mV**
- C/50 qOCV RMSE: **6.24 mV**
- endpoint 최대 절대오차: **10.00 mV**

이 window는 punched-electrode 질량에 직접 의존하지 않는 effective-capacity fit이다. 따라서 이후 teardown 형상을 반영할 때는 이미 맞춘 `Qn/Qp`를 보존하도록 `c_s,max`를 함께 환산했다.

## 3. 현재 반영 파라미터

| 구분 | 현재값 | 근거/상태 |
|---|---:|---|
| `Ln` | 76.5 µm | Ai2020 nominal 기준으로 복귀 |
| `Lp` | 68.0 µm | Ai2020 nominal 기준으로 복귀 |
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
| effective `csn,max` | 24,325.5 mol/m³ | 현재 `Qn`과 `Ln`을 보존한 환산값 |
| effective `csp,max` | 47,467.3 mol/m³ | 현재 `Qp`와 `Lp`를 보존한 환산값 |

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

## 6. 전극 두께 반영

현재 단면 코팅 두께는 측정법 불확실성이 정리될 때까지 Ai2020 nominal 값을 적용한다.

```text
Ln = 76.5 µm
Lp = 68.0 µm
```

기존 OCP fitting에서 식별한 active-material inventory가 바뀌지 않도록 `c_s,max`를 두께에 반비례시켰다. 따라서 이 변경은 기하학적 두께만 반영하며 `Qn/Qp` 자체를 증가시키지 않는다.

기존 가정 두께로 수행한 A/B에서는 동적 영향이 작았다. 현재 기준 구성은 Ai2020 nominal `76.5/68.0 µm`를 유지하고, 집전체 실측 기반 두께는 sensitivity 결과로만 보관한다.

## 7. 현재 최종 동적 검증 결과

현재 설정은 실험 `kn/kp`, Ai2020 nominal `Ln/Lp`, nominal Bruggeman을 동시에 사용한다.

| C-rate | 방향 | 10–70% SOC MAE | 용량오차 |
|---:|---|---:|---:|
| 0.5C | Charge | 14.36 mV | -89.26 mAh |
| 0.5C | Discharge | 5.10 mV | +27.43 mAh |
| 1C | Charge | 9.46 mV | -86.70 mAh |
| 1C | Discharge | 10.14 mV | +21.56 mAh |
| 2C | Charge | 13.99 mV | -102.58 mAh |
| 2C | Discharge | 14.05 mV | -18.21 mAh |

- 평균 전압 MAE: **11.18 mV**
- 충전 평균 MAE: **12.60 mV**
- 방전 평균 MAE: **9.76 mV**
- 평균 절대 용량오차: **57.62 mAh**
- 최대 절대 용량오차: **102.58 mAh**

음극 hysteresis scale을 0/0.25/0.5/0.75/1.0으로 ablation한 결과, C/50 branch 평균 MAE는 6.80/9.82/17.79/28.11/39.15 mV로 단조 악화했다. 동적 검증에서도 시험한 모든 두께 조건에서 같은 방향이므로 현재 구성에서는 음극 hysteresis를 사용하지 않는다.

중간 SOC 전압 형상은 개선됐지만 충전 종료 용량오차가 남아 있다. 이 오차를 Bruggeman만으로 줄이면 kinetic/transport와 capacity/endpoint 영향이 섞일 수 있으므로, 다음 단계에서는 usable capacity, stoichiometric endpoint, cutoff 부근 OCP를 분리해야 한다.

원시험 로그의 용량을 전류-시간으로 재적분한 결과 cycler 기록과의 차이는 주요 CC/CV step 모두 0.1% 미만이었다. 실제 평균전류를 적용하면 충전 평균 용량오차는 -92.85→-83.03 mAh로 감소하지만 잔차가 남는다. 양·음극 용량을 함께 4% 늘리면 충전은 맞지만 방전이 +78.66 mAh로 악화한다. 양극 유효용량만 4% 늘리면 충전 평균 오차는 -38.95 mAh로 줄고 방전 평균 절대오차는 23.31 mAh로 유지되어, cathode effective capacity/OCP mapping 재검토가 우선이다.

두께측정기와 SEM 차이에 대한 sensitivity로 `Ln/Lp=74/66 µm`, Ai2020 `76.5/68 µm`, 총두께-집전체 `73/62.5 µm`, 현차셀 SEM/두께측정기 비율을 이식한 `89.36/72.54 µm`를 비교했다. 6조건 평균 중앙 RMSE는 각각 12.35/13.46/11.69/19.85 mV였다. 다른 셀의 SEM 보정비를 Enertech에 이식할 근거는 없으므로 현재는 Ai2020 nominal `76.5/68 µm`를 유지한다.

## 8. Capacity-consistent `Qn/Qp-c_s,max-window` 재식별

현차셀 코드의 장점인 전극용량-화학양론 폭-`c_s,max` 연결을 반영하되, 다음 equality를 모든 후보에서 정확히 유지했다.

```text
delta_x = Qcell / Qn
delta_y = Qcell / Qp
Qn,p = F * A * L_n,p * eps_s,n,p * c_s,max,n,p / 3600
```

따라서 `c_s,max`만 독립적으로 fitting하지 않고 `Qn/Qp`, `delta_x/delta_y`, endpoint를 함께 갱신했다. 실제 C/50 qOCV를 target으로 사용하고, 고율 검증은 실제 branch 전류와 직전 rest 종점전압 기반 초기 SOC를 적용했다.

가중합 대신 다음 제약을 적용했다: qOCV RMSE 악화 ≤0.5 mV, 충전/방전 중앙 RMSE 악화 ≤1/2 mV, 충전/방전 용량 RMSE ≤3.1/1.0%. 그 안에서 `Qn/Qp` 변경량이 가장 작은 후보를 선택했다.

| 항목 | 기존 | 보수적 개선 후보 |
|---|---:|---:|
| `Qn` | 2.47947 Ah | 2.52906 Ah (+2%) |
| `Qp` | 4.37120 Ah | 4.45862 Ah (+2%) |
| `csn,max` | 24,325.5 | 24,812.0 mol/m³ |
| `csp,max` | 47,467.3 | 48,416.6 mol/m³ |
| `delta_x` | 0.950415 | 0.931779 |
| `delta_y` | 0.539103 | 0.528533 |
| qOCV RMSE | 6.302 mV | 6.494 mV |
| 충전 중앙 RMSE | 20.28 mV | 20.86 mV |
| 방전 중앙 RMSE | 8.74 mV | 10.60 mV |
| 충전 용량 RMSE | 3.75% | 3.05% |
| 방전 용량 RMSE | 1.24% | 0.92% |
| 충전 평균 절대 용량오차 | 68.50 mAh | 56.01 mAh |
| 방전 평균 절대 용량오차 | 22.19 mAh | 13.67 mAh |

`kn` 재탐색에서는 현재 EIS 값 `7.40e-7`이 3.1% 용량 제약을 만족하면서 충전 전압 RMSE가 가장 낮았다. 따라서 `kn`은 변경하지 않는다. 이 후보는 질량/조성 독립 측정 전까지 **보수적 provisional refinement**이며, 기존 기준 결과를 삭제하거나 물성 확정값으로 간주하지 않는다.

이 equality는 C/50 qOCV/저율 reversible capacity 단계에만 적용한다. 0.5C/1C/2C cutoff 용량에는 강제하지 않는다. 고율 종점은 아직 확정되지 않은 `Rct/k`, `Ds`, electrolyte/ohmic polarization의 영향을 받기 때문이다. `csmax` 변경은 교환전류식의 `j0`에도 영향을 주므로, 향후 SOC별 EIS에서 charge-transfer `Rct`가 분리되면 intrinsic `k` 고정과 EIS-`j0` 보존(`k ∝ 1/csmax`) 두 해석을 비교해야 한다.

## 9. SPB655060 원 논문 geometry 적용 원칙

- 집전체 두께는 원 논문값 `Cu=10 um`, `Al=15 um`를 기준값으로 사용한다. 사용자 측정 `9/13 um`는 sensitivity 기록으로만 남긴다.
- teardown 면적은 `A_n=820.352 cm2`, `A_p=787.236 cm2`로 기록한다.
- 표준 1D DFN의 current/reaction area는 작은 cathode footprint인 `A_overlap=A_p`를 사용한다. 서로 다른 `A_n/A_p`를 양·음극 전류식에 직접 넣지 않는다.
- 원 논문의 porosity는 mercury porosimetry 측정값 `eps_e,p=0.33`, `eps_e,n=0.32`다.
- active-material fraction `eps_s,p=0.62`, `eps_s,n=0.61`은 직접 측정이 아니라 inactive solid `0.05/0.07` 가정에서 계산한 값이다.

집전체 `10/15 um`는 Ai2020 기본값과 같아 현재 isothermal DFN 결과를 바꾸지 않는다. DFN overlap 면적을 `814.98 -> 787.236 cm2`로 바꾸면 전류밀도가 3.52% 증가한다. 면적만 바꾼 raw 후보는 charge/discharge capacity RMSE가 `7.47/4.19%`로 악화했다. `c_s,max`를 면적에 반비례시켜 Q를 보존하면 중앙 전압 RMSE는 `20.86/10.60 -> 19.88/9.18 mV`로 개선되지만 capacity RMSE는 `3.05/0.92 -> 3.71/1.01%`로 악화한다. 원 porosity까지 적용하면 capacity RMSE가 `5.04/1.31%`로 더 커진다. 따라서 paper geometry는 bookkeeping에 채택하되, 내일 dry loading으로 Q를 재계산하기 전에는 기존 동적 baseline을 교체하지 않는다.

## 10. 주요 실행 파일과 결과

- 현재 설정 실행: `nominal_radius_current_configuration.py`
- EIS kinetic ablation: `eis_first_kinetics_application.py`
- Bruggeman grid: `eis_bruggeman_screen.py`
- separator Bruggeman refinement: `eis_bruggeman_separator_refine.py`
- 전극 두께 A/B: `current_thickness_estimate_application.py`
- 음극 hysteresis/실측 집전체 재검토: `anode_hysteresis_thickness_recheck.py`
- 용량오차 원인 ablation: `capacity_error_cause_audit.py`
- 두께 측정법 sensitivity: `thickness_measurement_method_sensitivity.py`
- fitting 전 시간-전압 및 k/Ds 조합: `prefit_kinetics_diffusion_ablation.py`
- transport ablation: `transport_parameter_ablation.py`
- `Rn/Rp` 조합: `ocp_rn_rp_combination_screen.py`
- OCP 후보/stoichiometry 분석: `mass_independent_nominal_anode_trial.py`, `finalize_ocp_recommendation.py`
- capacity-consistent 개선: `capacity_consistent_qn_qp_kn_refinement.py`
- 개선 후보 재사용 설정: `capacity_consistent_refined_configuration.py`
- 원 논문 geometry 상수/정책: `paper_teardown_geometry.py`
- paper area/porosity ablation: `paper_geometry_area_porosity_ablation.py`

주요 산출물:

- `results/260921_nominal_radius_current_configuration/`
- `results/260921_eis_first_kinetics/`
- `results/260921_eis_bruggeman_screen/`
- `results/260921_current_thickness_estimate/`
- `results/260921_transport_parameter_ablation/`
- `results/260920_mass_independent_nominal_anode/`
- `results/260921_capacity_consistent_qn_qp_kn_refinement/`
- `results/260921_ai2020_geometry_source_audit/`
- `results/260921_coating_loading_precheck/`
- `results/260921_paper_geometry_area_porosity_ablation/`

## 11. 재현

```powershell
.\.venv\Scripts\python.exe nominal_radius_current_configuration.py
.\.venv\Scripts\python.exe eis_first_kinetics_application.py
.\.venv\Scripts\python.exe eis_bruggeman_screen.py
.\.venv\Scripts\python.exe current_thickness_estimate_application.py
.\.venv\Scripts\python.exe capacity_consistent_qn_qp_kn_refinement.py
.\.venv\Scripts\python.exe paper_geometry_area_porosity_ablation.py
```

현재 모델에는 Ai2020 nominal 단면 코팅 두께 `Ln=76.5 µm`, `Lp=68 µm`를 적용했다.
