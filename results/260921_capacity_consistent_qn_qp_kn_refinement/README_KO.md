# Capacity-consistent 개선 결과

## 방법

기존처럼 `c_s,max`만 바꾸지 않고 아래 관계를 모든 후보에서 강제했다.

```text
Qn * delta_x = Qp * delta_y = Qcell
c_s,max,n / c_s,max,n,0 = Qn / Qn,0
c_s,max,p / c_s,max,p,0 = Qp / Qp,0
```

각 `Qn/Qp` 후보마다 window 폭을 먼저 고정하고, C/50 실험 qOCV에 대해 `x0`, `y100` 위치만 다시 최적화했다. 0.5C/1C/2C 검증에는 실제 평균전류와 직전 rest 종점전압 기반 초기 SOC를 사용했다. 음극 hysteresis는 끄고 양극 hysteresis는 직전 branch history와 함께 유지했다.

선택은 가중합이 아니라 다음 제약을 사용했다.

- qOCV RMSE: 기준 대비 +0.5 mV 이내
- 충전/방전 중앙 RMSE: 기준 대비 +1/+2 mV 이내
- 충전/방전 용량 RMSE: 3.1/1.0% 이하
- 위 조건에서 `Qn/Qp` 변경량 최소

## 선택 결과

| 항목 | 기준 | 개선 후보 | 변화 |
|---|---:|---:|---:|
| `Qn` | 2.47947 Ah | 2.52906 Ah | +2.0% |
| `Qp` | 4.37120 Ah | 4.45862 Ah | +2.0% |
| `csn,max` | 24,325.5 | 24,812.0 mol/m³ | +2.0% |
| `csp,max` | 47,467.3 | 48,416.6 mol/m³ | +2.0% |
| `delta_x` | 0.950415 | 0.931779 | -2.0% |
| `delta_y` | 0.539103 | 0.528533 | -2.0% |
| `x0 → x100` | 0.00381 → 0.95422 | 0.00325 → 0.93503 | 재최적화 |
| `y100 → y0` | 0.43164 → 0.97075 | 0.43297 → 0.96150 | 재최적화 |
| qOCV MAE / RMSE | 4.78 / 6.30 mV | 5.07 / 6.49 mV | +0.29 / +0.19 mV |

`kn`은 7.4e-7~1.35e-6 범위에서 다시 평가했다. 현재 EIS 값 `7.40e-7`이 충전 용량 RMSE 3.1% 제약을 만족하는 후보 중 전압 RMSE가 가장 낮으므로 유지했다.

## 동적 검증

| C-rate | 방향 | 중앙 RMSE | 용량오차 |
|---:|---|---:|---:|
| 0.5C | Charge | 29.47 mV | -63.50 mAh |
| 1C | Charge | 18.81 mV | -60.27 mAh |
| 2C | Charge | 14.30 mV | -44.25 mAh |
| 0.5C | Discharge | 12.66 mV | +1.76 mAh |
| 1C | Discharge | 12.56 mV | -4.16 mAh |
| 2C | Discharge | 6.59 mV | -35.11 mAh |

동일한 실제전류/rest 초기화 기준에서 평균 결과는 다음과 같다.

| 지표 | 기준 | 개선 후보 |
|---|---:|---:|
| 충전 중앙 RMSE | 20.28 mV | 20.86 mV |
| 방전 중앙 RMSE | 8.74 mV | 10.60 mV |
| 충전 용량 RMSE | 3.75% | 3.05% |
| 방전 용량 RMSE | 1.24% | 0.92% |
| 충전 평균 절대 용량오차 | 68.50 mAh | 56.01 mAh |
| 방전 평균 절대 용량오차 | 22.19 mAh | 13.67 mAh |

용량오차는 양 방향 모두 감소하지만 전압 RMSE가 소폭 증가한다. 따라서 측정 질량/조성 없이 더 큰 `c_s,max` 변경은 채택하지 않고 +2% 후보를 provisional refinement로 둔다.

## 현차셀 방식과의 차이

현차셀도 최종적으로 `Qcell`, window 폭, 전극 형상을 이용해 `c_s,max`를 역산한다는 점은 같다. 차이는 다음과 같다.

- 현차셀: tear-down `c_s,max`로 span을 계산하지만 endpoint를 span 내부에서 독립 최적화한 뒤 `c_s,max`를 다시 역산한다. 초기 tear-down 값은 실질적으로 bound/anchor다.
- 현재 개선: 후보 `Qn/Qp`가 정해지면 `delta_x=Qcell/Qn`, `delta_y=Qcell/Qp`를 equality로 고정한다. 따라서 capacity와 window가 수치적으로 분리되지 않는다.
- 현차셀은 실측 질량·조성·면적·층수 anchor가 있지만, Enertech는 아직 이 독립 anchor가 없다. 그러므로 현재 `c_s,max`는 소재 고유값 확정치가 아니라 effective capacity parameter다.

최종 확정에는 양·음극 양면 코팅 질량에서 집전체 질량을 뺀 값, 활물질 질량분율, 실제 전극 면적/적층 수가 필요하다.

## Equality 강제 범위와 Rct 주의사항

`Qn*delta_x = Qp*delta_y = Qcell`은 저율에서 동일한 reversible lithium throughput을 표현하는 물질수지다. 다음 조건에서는 논리적으로 타당하다.

- formation 이후 Coulombic efficiency가 거의 100%인 안정 사이클
- 같은 온도와 같은 상·하한 전압 사이의 저율 CC 용량
- side reaction, lithium plating, active-material isolation이 무시 가능한 경우
- `Qcell`이 CV tail이나 고율 polarization으로 정의되지 않은 경우

반대로 0.5C/1C/2C cutoff 용량에 이 equality를 직접 강제하면 안 된다. 고율에서는 `Rct/k`, `Ds`, electrolyte transport와 ohmic resistance 때문에 표면농도가 먼저 cutoff에 도달하며 electrode inventory 전체가 사용되지 않는다. 현재 개선 코드는 equality를 C/50 qOCV/window 단계에만 적용하고 고율 용량은 validation으로만 사용했다.

또한 `j0 = k*sqrt(ce)*sqrt(cs)*sqrt(csmax-cs)` 형태에서는 고정 stoichiometry에서 `j0`가 대략 `k*csmax`에 비례한다. 따라서 `csmax`를 바꾸면 두 가지 해석이 가능하다.

1. intrinsic `k`가 독립 물성이라고 보면 `k`를 유지하고 `Rct` 변화는 모델 결과로 둔다.
2. EIS에서 얻은 `j0/Rct`를 고정하고 `csmax`만 재표현한 것이라면 `k_new = k_EIS*csmax_EIS/csmax_new`로 보정한다.

현재 Enertech의 `Rct`는 확정값이 아니므로 1번을 사용했고, +2% `csmax` 변경과 kinetics를 완전히 분리했다고 주장하지 않는다. EIS 등가회로에서 charge-transfer arc가 분리되면 2번 보정을 별도 검증해야 한다.

현차셀 코드는 2번을 사용한다. `KAPPA_N_RAW/P_RAW`를 EIS 값으로 두고 최종 역산된 `CSN/CSP`에 대해 `KAPPA = KAPPA_RAW*CSMAX_AT_EIS/CSMAX_FINAL`로 재조정한다. 즉 현차 결과의 좋은 일치는 `csmax` 변경 시 EIS `j0`를 보존한다는 추가 가정도 포함한다.

## Enertech에서 추가로 확보할 값

| 우선순위 | 측정/확인값 | 사용 목적 |
|---:|---|---|
| 1 | 같은 면적의 양·음극 양면 총질량과 집전체 질량, 세척·건조 후 반복값 | coating loading 및 `Qn/Qp` 독립 계산 |
| 2 | 양·음극 활물질 질량분율, binder/carbon 조성 | coating 질량에서 active-material 질량 분리 |
| 3 | 단면 coating 두께, 집전체 두께, 전극 기공률 | `eps_s*L`과 effective capacity 계산 |
| 4 | 실제 단면 전극 면적과 병렬 electrode-pair 수 | 총 active area 계산; 양면/층수 이중계산 방지 |
| 5 | 동일 셀·온도의 C/50 또는 충분히 낮은 C-rate CC 충·방전과 rest | `Qcell`, pseudo-OCV, endpoint 기준 |
| 6 | 양·음극 half-cell OCP/GITT 반복 측정 | OCP source와 hysteresis 검증 |
| 7 | SOC별 EIS에서 `R_ohm`, `R_SEI/contact`, `Rct` 분리 및 반복 | `kn/kp` 및 `csmax`-kinetics coupling 결정 |

회수 음극에 전해질/SEI가 남아 있으면 coating mass가 과대평가될 수 있다. 가능하면 pristine sibling electrode를 사용한다. 회수 전극만 가능하면 동일 세척·건조 protocol로 여러 시편을 반복 측정하고, 같은 면적의 박리된 Cu/Al 집전체 질량을 빼며 그 분산을 uncertainty로 보고한다. 양면 코팅 시편은 `(총질량-집전체질량)/2`로 한 면 loading을 구하되, 모델의 병렬 면적/층수 정의에서 다시 2배 하지 않도록 한다.
