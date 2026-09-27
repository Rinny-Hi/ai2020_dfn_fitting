# 260927 Enertech DFN fitting 인수인계

> 이 문서는 2026-09-27 작업 종료 시점의 최신 결정이다. 이전 `260921_CURRENT_PROGRESS.md`, `260925_HANDOFF_CURRENT_STATE_KO.md` 및 과거 result 문서와 충돌하면 이 문서와 `results/260927_charge_fitting_design/README_KO.md`를 우선한다.

## 1. 현재 상태

- GitHub: <https://github.com/Rinny-Hi/ai2020_dfn_fitting>
- 작업 브랜치: `260921-ocp-eis-current-progress`
- 다음 단계는 fitting 실행이 아니라, 먼저 **최종 physical-time 목적함수와 동일한 grid에서 선택 subset의 민감도/식별도를 재확인**하는 것이다.
- 5개 파라미터 동시 fitting은 사용자 요청으로 중단했다. `results/260927_p0_charge_fit_discharge_validation/optimizer_evaluation_history.csv`는 최종값이 아닌 중간 이력이다.

## 2. 최신 물리 baseline

| 항목 | 적용값 | 상태/근거 |
|---|---:|---|
| overlap area | `Ap=787.236 cm2` | 양극 45.4 mm x 51 mm x 34면, 1D DFN 공통 반응면적 |
| negative physical area | `An=820.352 cm2` | 음극 46.4 mm x 52 mm x 34면, inventory/N/P bookkeeping용 |
| `Ln/Lp` | `76.5/68 um` | Ai2020 nominal 유지 |
| `Rn/Rp` | `3.7542/4.2387 um` | SEM ImageJ 대표값 |
| `eps_s,n/eps_s,p` | `0.61/0.62` | Ai2020 |
| `brugg_n/p/s` | `2.914/1.83/1.5` | 현재 prefit baseline |
| `c_s,n,max` | `25182.812 mol/m3` | teardown overlap area에서 Qn 보존한 effective 값 |
| `c_s,p,max` | `49140.112 mol/m3` | teardown overlap area에서 Qp 보존한 effective 값 |
| `kp` | 약 `5.90e-7` | cathode P4 `Rct_ASR=18.78433 ohm cm2` 환산 apparent EIS 값 |
| `kn` | `9.6485e-7` 초기값 | 신뢰 가능한 anode Rct가 없어 Ai2020 nominal 환산값 사용, fitting 대상 |

`c_s,max`는 독립적으로 측정된 intrinsic 최대농도가 아니라, 현재 OCP window와 electrode inventory를 공통 overlap area에 보존한 effective parameter다. `An` 전체 면적과 현재 `c_s,n,max`를 동시에 용량식에 다시 넣으면 음극 inventory가 이중 반영된다.

## 3. OCP / stoichiometry baseline

- OCP 조합: nominal graphite anode + old measured cathode
- `Qcell=2.356527656 Ah`
- `x0=0.0038062`, `x100=0.9542213`
- `y100=0.4316424`, `y0=0.9707458`
- `delta_x=0.9504150`, `delta_y=0.5391034`
- `Qn=2.4794722 Ah`, `Qp=4.3711981 Ah`
- 각 0.5C/1C/2C branch의 초기 SOC는 직전 rest 종점전압으로 재설정한다.

## 4. 최신 GITT 결정

고정 10% endpoint 제거에는 보편적 문헌 규칙이 없으므로, 재현 baseline은 전체 pulse를 사용한 **all-range apparent GITT reference**로 변경했다.

| 전극 | 0% charge | 0% discharge | 0% 결합값 | 단위 |
|---|---:|---:|---:|---|
| anode | 6.792 | 2.879 | `4.422` | `1e-14 m2/s` |
| cathode | 4.558 | 10.171 | `6.809` | `1e-14 m2/s` |

- 최신 fitting 초기값: `Dsn=4.422222641e-14`
- 1차 charge fitting 고정값: `Dsp=6.809001552e-14`
- 10% 결과는 폐기하지 않고 endpoint 민감도/강건성 비교로 유지한다.
- 이 값들은 intrinsic crystal diffusivity 확정값이 아니다. 600 s pulse, composite coating mass 100 wt% active 가정, 전압 resolution, finite-particle 효과를 포함한 apparent 값이다.
- 상세: `results/260927_gitt_ds_revalidation/README_KO.md`

## 5. 민감도/식별도 결과

charge 0.5C/1C/2C mixed LSA 결과:

| parameter | relative sensitivity | 해석 |
|---|---:|---|
| `kn` | 1.000 | 유지 |
| `kp` | 0.983 | `kn`과 correlation 약 0.996이므로 동시 fitting 제외 |
| `brugg_n` | 0.771 | 주 모델에 유지 |
| `Dsn` | 0.284 | 유지 |
| `Dsp` | 0.023 | charge 1차 fitting에서 제외/고정 |

선택 subset:

1. 주 모델: `Dsn + kn + brugg_n`
2. 보수 모델: `Dsn + kn`

주의: 이 mixed LSA는 normalized transferred-capacity 축을 사용했다. 실제 fitting은 physical-time voltage RMSE이므로, fitting 전에 동일 time grid에서 두 subset의 Jacobian/correlation을 재계산해야 한다.

## 6. fitting 목적함수

charge July BoL 3셀의 0.5C/1C/2C를 사용한다.

\[
J_V(\theta)=
\sqrt{\frac{1}{3N_{cell}N_t}
\sum_{r\in\{0.5,1,2\}}\sum_j\sum_i
\left[1000\left(V_{DFN,r}(t_i;\theta)-V_{exp,r,j}(t_i)\right)\right]^2}
\]

- 각 C-rate와 셀을 동일 가중한다.
- 각 branch를 physical-time 기준 100 point로 재표본화한다.
- voltage ordinary L2/RMSE만 최소화한다.
- capacity, prior, voltage slope/shape 항은 목적함수에 넣지 않는다.
- 실험 time grid 전에 model cutoff가 발생하면 유효하지 않은 candidate로 처리한다.
- capacity error, full-range RMSE, 10--70% RMSE는 fitting 이후 별도 평가한다.
- charge fitting 결과는 untouched discharge 0.5C/1C/2C로 validation한다.

## 7. 최종 1차 bound

| parameter | initial | 1차 bound | 근거 |
|---|---:|---:|---|
| `Dsn [m2/s]` | `4.422e-14` | `1.0e-14--8.0e-14` | combined endpoint 최소 1.689e-14와 0% charge 6.792e-14를 모두 포함 |
| `kn` | `9.6485e-7` | `3.0e-7--3.0e-6` | Ai2020 `kref=1e-11 m/s`를 `F*kref`로 환산 후 약 +/-0.5 log10 decade |
| `brugg_n` | `2.914` | `1.5--3.5` | 이론 하한 1.5, graphite tortuosity 환산 약 2.73, Ai2020 2.914 포함 |

중요 수정:

- 이전 `Dsn` 상한 `5.0e-14`는 charge 0% GITT 값 `6.792e-14`를 배제하므로 `8.0e-14`로 수정했다.
- 이전 `brugg_n` 상한 `4.0`은 1차 bound의 직접 근거가 약해 `3.5`로 줄였다.
- boundary test: `Dsn -> 1.2e-13`, `kn -> 2.0e-7--5.0e-6`, `brugg_n -> 4.0`.
- `Dsn`, `kn`은 log space, `brugg_n`은 linear space에서 최적화한다.

## 8. 최적화 알고리즘 설계

두 모델에 동일한 목적함수와 algorithm ladder를 적용한다.

1. Multi-start TRF least-squares
   - 0% GITT, 10% GITT, bound 중앙 및 log-space 시작점 총 5개
   - 주 모델 `max_nfev=80`, 보수 모델 `max_nfev=60`
2. Differential Evolution -> TRF polish
   - `popsize=6`, `maxiter=8`, `seed=260927`
3. Dual Annealing -> TRF polish
   - 주 모델 `maxfun=180`, 보수 모델 `maxfun=120`, `seed=260927`

선택 기준은 charge training RMSE 단독 최저가 아니다. 반복/알고리즘 간 수렴성, bound 비접촉, discharge validation, 사후 capacity error를 함께 본다. 주 모델이 validation에서 이득이 없거나 `brugg_n`이 경계에 붙으면 보수 모델을 선택한다.

## 9. 다음 실행 순서

1. 코드의 fitting 전용 config를 최신 0% GITT와 위 bound로 분리한다. 과거 `p0_five_parameter_sensitivity_area.py`의 old initial/bound를 그대로 재사용하지 않는다.
2. physical-time 목적함수로 `Dsn+kn+brugg_n` 및 `Dsn+kn`의 민감도/correlation/condition number를 재검증한다.
3. 식별도 통과 시 두 모델에 Multi-start TRF, DE->TRF, DA->TRF를 실행한다.
4. charge fitting 결과를 discharge validation하고, capacity는 사후 지표로만 비교한다.
5. optimizer별 parameter, objective history, Experiment/Fit 곡선, validation 곡선, capacity error를 저장한다.

## 10. 핵심 파일

- fitting 설계: `results/260927_charge_fitting_design/README_KO.md`
- GITT 재검증: `results/260927_gitt_ds_revalidation/README_KO.md`
- endpoint ablation: `results/260927_ds_endpoint_trim_ablation/summary.csv`
- 민감도/식별도: `results/260927_p0_charge_mixed_lsa_correlation/README_KO.md`
- area 및 5 parameter LSA: `results/260927_p0_five_parameter_sensitivity_area/`
- prefit 물리 anchor: `results/260927_prefit_physical_anchor_decision/README_KO.md`
- July/new RPT 비교: `results/260927_comprehensive_prefit_cross_cohort/README_KO.md`
- 발표 그림: `results/260927_presentation_figures/`
- 중단된 5 parameter fitting 상태: `results/260927_p0_charge_fit_discharge_validation/STATUS_KO.md`

## 11. 새 채팅 시작 문구

> GitHub의 `260921-ocp-eis-current-progress` 브랜치에서 `260927_HANDOFF_CURRENT_STATE_KO.md`를 먼저 읽고 이어서 진행해줘. 이전 문서와 충돌하면 이 문서와 `results/260927_charge_fitting_design/README_KO.md`를 우선해. 현재 baseline은 all-range apparent GITT `Dsn=4.422e-14`, `Dsp=6.809e-14`이고, 먼저 physical-time 목적함수로 주 모델 `Dsn+kn+brugg_n`과 보수 모델 `Dsn+kn`의 민감도/식별도를 재검증해. 통과한 모델만 charge 0.5C/1C/2C fitting 후 discharge validation해. capacity는 목적함수에 넣지 말고 사후 지표로 평가해.

