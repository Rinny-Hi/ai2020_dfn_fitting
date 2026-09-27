# 0.5C/1C/2C Charge fitting 설계

## 1. 비교할 모델

| 모델 | fitting parameter | 고정 parameter | 목적 |
|---|---|---|---|
| 주 모델 | `Dsn`, `kn`, `brugg_n` | `Dsp`, `kp`, `brugg_p`, `brugg_s` | 고체 확산·반응속도·음극 전해질 수송을 함께 식별 |
| 보수 모델 | `Dsn`, `kn` | `brugg_n`, `Dsp`, `kp`, `brugg_p`, `brugg_s` | 상관성과 과적합 위험을 낮춘 최소 subset |

고정값은 `Dsp=6.809001552e-14 m2/s`, `kp≈5.90e-7`, `brugg_n=2.914`(보수 모델만), `brugg_p=1.83`, `brugg_s=1.5`를 사용한다. `Dsp`는 all-range apparent GITT reference이고 `kp`는 cathode P4를 area-specific Rct로 해석해 환산한 apparent EIS 초기값이다. 두 값 모두 확정 intrinsic 물성이라는 의미는 아니다.

## 2. 목적함수

각 C-rate 및 각 셀이 동일한 기여를 갖는 ordinary voltage least-squares를 사용한다.

\[
J_V(\theta)=
\sqrt{\frac{1}{N_RN_CN_t}
\sum_{r\in\{0.5,1,2\}}
\sum_{c=1}^{N_C}
\sum_{i=1}^{N_t}
\left[1000\left(V_{\mathrm{DFN},r}(t_i;\theta)-V_{\mathrm{exp},r,c}(t_i)\right)\right]^2}
\]

- 단위: mV
- `N_R=3`, 현재 July BoL cohort는 `N_C=3`
- 각 셀·C-rate를 `N_t=100`개의 동일 간격 physical-time point로 재표본화한다.
- 각 시험의 CC 시작 시점을 `t=0`으로 두고, 직전 rest 전압으로 초기 SOC를 정한다.
- 전류가 command의 95%에 도달한 첫 record부터 실험 CC cutoff 직전까지를 사용한다.
- C-rate 가중치, 형상 미분항, capacity 항, prior/regularization 항을 추가하지 않는다.
- candidate가 고정 실험 time grid를 덮기 전에 cutoff에 도달하거나 solver가 실패하면 유효하지 않은 candidate로 처리한다. 이는 별도의 capacity 목적함수가 아니라 time-voltage 잔차를 정의할 수 없다는 feasibility 처리이다.
- fitting 후에는 full-range RMSE, 10--70% RMSE, cutoff capacity error를 별도로 산출한다.

현재 mixed LSA는 normalized transferred-capacity 축으로 수행되었으므로, 실제 fitting 직전에 위 physical-time objective와 동일한 grid로 선택된 subset의 Jacobian/correlation을 한 번 재계산한다. subset이 유지될 때만 최적화를 시작한다.

## 3. 초기값과 bound

양수 파라미터 `Dsn`, `kn`은 log space에서, `brugg_n`은 linear space에서 최적화한다.

| Parameter | Initial | Lower | Upper | 근거/처리 |
|---|---:|---:|---:|---|
| `Dsn [m2/s]` | `4.422222641e-14` | `1.0e-14` | `8.0e-14` | endpoint 처리 0--20%의 combined 최소 `1.689e-14`와 0% charge 대표값 `6.792e-14`를 모두 포함; charge-only fitting에서 방향별 실험값을 bound 밖으로 배제하지 않음 |
| `kn` | `9.648533212e-7` | `3.0e-7` | `3.0e-6` | Ai2020 `kref=1e-11 m/s`를 현재 Butler--Volmer 구현의 `F*kref=9.6485e-7`로 환산하고, 이를 중심으로 약 `±0.5 log10 decade`(약 ÷3.16, ×3.16)를 탐색 |
| `brugg_n` | `2.914` | `1.5` | `3.5` | 이상 구형 입자 Bruggeman 하한 `1.5`, graphite tortuosity 실험 환산 약 `2.73`, Ai2020의 `1+alpha_B=2.914`를 모두 포함 |

`Dsn` 최적값이 상한의 log-distance 기준 5% 이내에 붙으면 결과를 채택하지 않고 상한을 `1.2e-13`으로 한 번 확장해 boundary test를 수행한다. `kn`이 경계에 근접하면 `2.0e-7--5.0e-6`, `brugg_n`이 3.5에 근접하면 4.0까지 확장한다. 이 확장값들은 1차 fitting bound가 아니라 경계 절단 여부를 확인하기 위한 진단 범위다.

### Bound 근거의 성격

- `Dsn`: 현재 실험 결과로 직접 정한 **measurement-envelope bound**다. 10% 제외값만이 아니라 0/5/10/15/20% endpoint 민감도와 방향별 결과를 사용했다. 특히 charge fitting 대상인데 0% charge 값 `6.792e-14`를 제외하는 기존 `5.0e-14` 상한은 논리적으로 맞지 않아 `8.0e-14`로 수정했다.
- `kn`: 동일한 Butler--Volmer 정의를 쓰는 Ai2020 nominal을 중심으로 정한 **nominal-centered uncertainty bound**다. 문헌의 `k`는 교환전류 정의와 단위가 달라 직접 최소·최댓값을 섞지 않았다. 재현성이 낮아 폐기한 anode Rct 값 `2.48e-7/7.40e-7`은 bound의 anchor가 아니라 사후 비교점으로만 남긴다.
- `brugg_n`: 고전 Bruggeman 값과 실제 graphite tortuosity 측정, Ai2020 값을 함께 포함한 **physics/literature bound**다. Ai2020 식 `Psi_eff=epsilon^(1+alpha_B)Psi_0`에서 음극 `alpha_B=1.914`, 즉 PyBaMM식 exponent는 `2.914`다. Landesfeind 등의 graphite 결과 `epsilon≈0.43`, `tau≈4.3`을 `tau=epsilon^(1-b)`로 환산하면 `b≈2.73`이다.

따라서 이 bound들은 통계적 95% 신뢰구간이 아니다. optimizer가 물리·실험 anchor 주변을 탐색하기 위한 범위이며, 최적값이 경계에 붙으면 그 값을 물성으로 채택하지 않고 bound 또는 식별도 문제로 판단한다.

## 4. 최적화 알고리즘

두 모델에 동일한 알고리즘 구조와 random seed를 적용한다.

### A. Multi-start TRF least-squares

- ordinary linear loss
- 시작점: 0% GITT, 10% GITT, bound 중심 및 log-space Latin-hypercube 시작점으로 총 5개
- 주 모델 `max_nfev=80`, 보수 모델 `max_nfev=60`
- `ftol=xtol=gtol=1e-4`
- 잔차 vector를 직접 사용하므로 최종 Jacobian, correlation 및 condition number를 계산할 수 있다.

### B. Differential Evolution + TRF polish

- 모든 변수를 동일한 normalized search space로 변환
- `popsize=6`, `maxiter=8`, `seed=260927`, optimizer 내부 polish는 끔
- 예상 global evaluation: 주 모델 약 162회, 보수 모델 약 108회
- DE 최적점을 TRF로 local polish

### C. Dual Annealing + TRF polish

- 주 모델 `maxfun=180`, 보수 모델 `maxfun=120`, `seed=260927`
- global optimum 재현성 비교용이며 최종 해는 TRF로 polish

`L-BFGS-B`는 event/cutoff에 따른 목적함수 비평활성 때문에 주 알고리즘으로 사용하지 않고 필요할 때만 scalar-objective 진단 비교에 사용한다.

## 5. 모델 및 해 선택 기준

1. charge fitting objective RMSE
2. 서로 다른 시작점/알고리즘에서의 parameter 수렴성
3. 최적값의 bound 접촉 여부
4. fitting에 사용하지 않은 discharge 0.5C/1C/2C time-voltage RMSE
5. fitting 후 cutoff capacity error
6. 셀별 결과 및 leave-one-cell-out 방향 일관성

주 모델의 charge 오차가 더 낮아도 `brugg_n`이 bound에 붙거나 discharge validation 개선이 반복되지 않으면 보수 모델을 선택한다. 반대로 주 모델이 여러 시작점에서 유사한 parameter로 수렴하고 세 C-rate 및 discharge에서 일관되게 개선되면 주 모델을 채택한다. 단순히 training RMSE가 가장 낮은 한 번의 run은 선택 근거로 사용하지 않는다.

## 6. 출력물

- optimizer별 최적 parameter와 목적함수 history
- charge 0.5C/1C/2C Experiment/Fit time-voltage overlay
- discharge 0.5C/1C/2C Experiment/Validation overlay
- C-rate별 full-range 및 10--70% RMSE
- cutoff capacity signed error, MAE, RMSE(%)
- Jacobian sensitivity, parameter correlation, singular value/condition number
- bound 접촉 및 solver failure log
