# 260927 Enertech DFN fitting 인수인계

> 이 문서는 2026-09-27 작업 종료 시점의 최신 결정이다. 이전 `260921_CURRENT_PROGRESS.md`, `260925_HANDOFF_CURRENT_STATE_KO.md` 및 과거 result 문서와 충돌하면 이 문서와 `results/260927_charge_fitting_design/README_KO.md`를 우선한다.

## 0. 최신 후속 실행 결과 (이 절이 아래 계획보다 우선)

- GITT 전극 질량을 이전 `10.40/17.25 mg`에서 최신 `10.10/17.35 mg`로
  보정했다. 최신 all-range apparent reference는
  `Dsn=4.170774146e-14`, `Dsp=6.888175324e-14 m2/s`이다.
- 질량 보정은 GITT 식에만 적용한다. 단일 펀칭 시편 질량으로 full-cell `Qn/Qp`를
  변경하지 않는다.
- graphite 문헌 및 선행 fitting bound를 비교한 뒤 1차 Dsn bound
  `1e-14--8e-14`를 유지했다. `1.2e-13` 상한 확장 test에서도 RMSE 개선은
  0.005 mV 미만이었다.
- final physical-time grid에서 주 모델 `Dsn+kn+brugg_n`과 보수 모델 `Dsn+kn`은
  모두 local 식별도 검사를 통과했다.
- 두 모델의 charge fitting 및 untouched discharge validation을 완료했다.
  주 모델은 charge objective `53.676 mV`, discharge full RMSE `62.08 mV`;
  보수 모델은 `53.997/49.51 mV`였다.
- 최종 선택은 **보수 모델**이다: `Dsn=7.4662e-14`, `kn=1.45024e-6`,
  `brugg_n=2.914` 고정. capacity는 목적함수에 사용하지 않았으며 사후 charge/
  discharge RMSE는 `3.11/1.09%`였다.
- fitted Dsn은 optimizer별 `5.74--7.47e-14`로 분산되어 intrinsic 물성 확정값으로
  해석하지 않는다.
- 전체 근거와 산출물: `results/260927_charge_physical_time_final/README_KO.md`

### 0.1 OCP/초기상태와 7월·9월 후속 검증 (2026-09-28)

- OCP는 **Ai2020 nominal graphite equilibrium + 기존 측정 cathode equilibrium**을
  유지한다. 새 cathode는 qOCV 수치만 일부 개선됐지만 reverse branch와 absolute
  anchor가 불완전해 진단용으로만 남겼다.
- July의 10분 rest는 평형이 아니며, 연속 protocol에서도 1C/2C rest 종점은 실험과
  4 mV 이내지만 charge onset은 약 140--145 mV 높다. 시작부 오차의 주원인은 OCP
  선택보다 저 SOC reaction overpotential 과대 예측이다.
- September의 2시간 rest는 July보다 잔류 slope를 charge 기준 약 8--90배 줄여
  초기 SOC/OCP 불확실성을 크게 낮춘다. 그래도 cohort별 inventory/window 차이는
  남는다.
- cathode hysteresis를 켜면 July 시작점은 일부 개선되지만 전체 RMSE가 악화하며,
  `brugg_n=1.5` 변경은 charge/discharge RMSE를 `85.45/103.26 mV`까지 악화한다.
- September charge에서는 `Dsn+kn` 중 Dsn 민감도가 기준을 통과하지 못해,
  `Dsn=4.170774e-14`를 독립 GITT 값으로 고정하고 `kn=5.348413e-7`만 fitting했다.
  결과는 charge/discharge `57.80/75.80 mV`, 사후 capacity RMSE
  `15.82/2.09%`다. 7월 모델보다 종합 성능이 낮다.
- 최종 7월 모델에서 `kn=2.59e-7`만 적용하면 charge/discharge 모두 조기 cutoff되어
  채택할 수 없다.
- Ai2020 nominal one-at-a-time 복원에서 raw PyBaMM `De`가 가장 큰 악화를 만들었다.
  현재 `De=raw Ai2020 function x 1e-4` 단위변환을 유지한다. raw De 적용 시
  charge/discharge RMSE는 `90.38/112.76 mV`다.
- 최신 상세 문서:
  - `results/260927_ocp_initial_state_deep_dive/README_KO.md`
  - `results/260927_ocp_newbol_hysteresis_brugg_comparison/README_KO.md`
  - `results/260927_september_recommended_dsn_kn_fit/README_KO.md`
  - `results/260928_july_final_result_graph/README_KO.md`
  - `results/260928_july_kn_2p59e7_check/README_KO.md`
  - `results/260928_july_nominal_reversion_audit/README_KO.md`

## 1. 현재 상태

- GitHub: <https://github.com/Rinny-Hi/ai2020_dfn_fitting>
- 작업 브랜치: `260921-ocp-eis-current-progress`
- 최종 physical-time 민감도/식별도, charge fitting, discharge validation까지 완료했다. 상세 결과는 `results/260927_charge_physical_time_final/README_KO.md`에 있다.
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
| anode | 6.406 | 2.715 | `4.171` | `1e-14 m2/s` |
| cathode | 4.611 | 10.290 | `6.888` | `1e-14 m2/s` |

- 최신 fitting 초기값: `Dsn=4.170774146e-14`
- 1차 charge fitting 고정값: `Dsp=6.888175324e-14`
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
| `Dsn [m2/s]` | `4.171e-14` | `1.0e-14--8.0e-14` | 최신 질량 기준 combined endpoint 최소 1.593e-14와 0% charge 6.406e-14를 모두 포함 |
| `kn` | `9.6485e-7` | `3.0e-7--3.0e-6` | Ai2020 `kref=1e-11 m/s`를 `F*kref`로 환산 후 약 +/-0.5 log10 decade |
| `brugg_n` | `2.914` | `1.5--3.5` | 이론 하한 1.5, graphite tortuosity 환산 약 2.73, Ai2020 2.914 포함 |

중요 수정:

- 이전 `Dsn` 상한 `5.0e-14`는 최신 질량 기준 charge 0% GITT 값 `6.406e-14`를 배제하므로 `8.0e-14`로 유지했다.
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

1. [완료] fitting 전용 config를 최신 질량 보정 GITT와 bound로 분리했다.
2. [완료] physical-time 목적함수로 두 subset의 민감도/correlation/condition number를 재검증했다.
3. [완료] 두 모델에 Multi-start TRF, DE->TRF, DA->TRF를 실행했다.
4. [완료] charge fitting 결과를 discharge validation하고 capacity를 사후 비교했다.
5. [완료] optimizer history, 곡선, validation 및 capacity 지표를 저장했다.
6. 다음 연구 단계는 constant D를 추가 최적화하기보다 SOC-dependent D/finite-sphere pulse fitting 또는 EIS·dynamic-current 기반 독립 식별이다.

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
