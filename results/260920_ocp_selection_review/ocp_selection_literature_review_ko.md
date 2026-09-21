# Teardown 이후 OCP 선택 및 파라미터 식별 기준

## 결론

OCP의 목적은 0.5C–2C 전압 MAE를 직접 최소화하는 것이 아니라, 전극의 열역학적 전위 함수와 full-cell의 lithium inventory 및 electrode balancing을 결정하는 것이다. 동적 파라미터를 식별하기 전에 이 기준선이 고정되어야 한다.

현재 ±5% 해는 수치적으로 가장 좋은 qOCV 후보지만 최종 물리 파라미터가 아니다. teardown에서 전극 두께와 활성물질 관련 값이 확보되면 Qn, Qp를 다시 계산하고 window를 재식별해야 한다. ±5% 경계는 그때 제거하고 측정 불확도를 prior로 사용한다.

## 파라미터별 OCP 영향

| 파라미터 | 평형 OCP 곡선 형상 | stoichiometry window/용량 | 동적 전압 | 처리 순서 |
|---|---|---|---|---|
| Ln, Lp | 직접 변화 없음 | Qn, Qp에 선형 영향 | 전해질 수송거리와 반응분포 영향 | teardown 후 OCP window 재계산 |
| Rn, Rp | 변화 없음 | 직접 영향 없음 | 비표면적 및 확산시간 영향 | OCP 고정 후 반영 |
| Dsn, Dsp | 변화 없음 | 영향 없음 | 입자 내부 농도구배와 relaxation 영향 | R과 함께 식별 |
| kn, kp | 변화 없음 | 영향 없음 | Butler–Volmer 과전압 영향 | OCP 고정 후 식별 |

전극 용량은 다음 관계를 따른다.

Qn = F A Ln eps_n c_s,max,n / 3600

Qp = F A Lp eps_p c_s,max,p / 3600

따라서 Ln과 Lp가 바뀌면 delta_x=Qcell/Qn 및 delta_y=Qcell/Qp가 바뀐다. 반면 R, Ds, k는 전류가 0이고 충분히 평형에 도달한 상태의 U(theta)를 바꾸지 않는다.

Rn과 Dsn, Rp와 Dsp는 각각 diffusion time constant R^2/D의 형태로 강하게 결합된다. 입자 반경이 teardown으로 변경되면 기존 GITT D 값을 그대로 사용할 수 없고 다시 계산해야 한다. 또한 a_s=3 eps_s/R이므로 R 변경은 k 또는 exchange-current 해석에도 영향을 준다.

## 선행연구의 OCP 선택 기준

1. Chen et al. 2020은 pseudo-OCV보다 충분히 이완된 GITT-OCV를 선택했다. 특히 낮은 SOC에서 pseudo-OCV가 relaxation OCV를 정확히 나타내지 못한다고 보고했다. 전극 window는 half-cell과 three-electrode full-cell의 특징점을 least-squares로 맞췄다. 동적 방전 MAE를 기준으로 OCP를 선택하지 않았다.
2. Ecker et al. 2015는 동일 전극에서 만든 half-cell/full-cell OCV의 변곡점과 dV/dQ 특징을 대응시켜 balancing을 계산했다. balancing은 셀마다 고유하며 문헌값을 가져올 수 없다고 명시했다.
3. Birkl et al. 2015는 half-cell과 full-cell OCV를 동시에 사용하고, charge/discharge와 온도별 hysteresis를 별도로 고려했다. 단일 평균곡선의 RMSE만으로 선택하지 않았다.
4. Hamed et al. 2024는 전압차만 최소화한 뒤 미분 특징이 맞지 않는 문제를 확인해 OCV 오차와 dU/dQ 오차를 함께 사용했다. 선택 결과가 charge/discharge branch 조합에 민감하다고 명시했다.

## OCP 후보 선택 우선순위

### 1순위: 데이터와 물리 근거의 적합성

- 실제 LCO/graphite 전극에서 측정했는가
- 충분한 relaxation 조건을 만족했는가
- charge/discharge hysteresis branch가 사용 목적과 일치하는가
- 동일 전극 lot와 세척/재조립 조건인가
- 필요한 stoichiometry 범위를 포함하는가

이 조건을 통과하지 못한 후보는 MAE가 낮아도 선택하지 않는다.

### 2순위: teardown 기반 물리 제약

- coating mass, thickness, porosity, active fraction, electrode area로 Qn/Qp 및 c_s,max 계산
- Qcell=Qn delta_x=Qp delta_y 강제
- 0 < x0 < x100 < 1, 0 < y100 < y0 < 1 강제
- N/P ratio와 lithium inventory가 실험 범위 안에 있는지 확인

### 3순위: 특징점 및 미분곡선

- OCV plateau 및 변곡점 위치
- dV/dQ 또는 dQ/dV peak 위치와 상대 크기
- 상하한 전압 endpoint

### 4순위: 전압 오차

- 중앙 SOC의 qOCV MAE/RMSE
- endpoint 또는 가파른 구간이 전체 점수에 과도하게 영향을 주지 않도록 구간별 확인
- 가중치는 임의값보다 측정 반복성 및 relaxation 불확도의 역분산으로 설정

### 5순위: 독립 검증

- 다른 half-cell replicate 및 다른 full-cell 저율 시험에서 window가 유지되는지 확인
- 0.5C–2C 곡선은 OCP 선택 기준이 아니라 OCP를 고정한 후 transport/kinetic 파라미터 검증에 사용

## 현재 결과의 해석

| 후보 | qOCV MAE | 동적 평균 MAE | 최대 용량오차 | 해석 |
|---|---:|---:|---:|---|
| 기존 고정 | 16.21 mV | 37.87 mV | 0.0833 Ah | Ai2020 물성 의존성이 큼 |
| ±5% window | 7.13 mV | 36.12 mV | 0.0568 Ah | qOCV와 용량 균형이 가장 좋지만 prior 경계에 걸림 |
| Chen mass/capacity | 12.66 mV | 33.40 mV | 0.0926 Ah | 충전 동적전압은 개선되지만 qOCV와 종료용량 악화 |

Chen 해의 동적 평균 MAE 감소는 주로 충전에서 발생했다. 충전 MAE는 19.8, 16.2, 15.7 mV였지만 방전은 43.4, 58.3, 46.9 mV였다. 따라서 이 개선을 평형 OCP가 더 정확해졌다는 증거로 볼 수 없다. 충방전 비대칭은 hysteresis branch, charge/discharge kinetic asymmetry, resistance 또는 concentration-dependent transport 오차의 특성이 더 강하다.

qOCV residual에서 SOC 약 0.5–0.7의 넓은 양의 hump가 후보 전체에 남는다. rate와 무관한 동일 SOC 위치의 형상오차이므로 전극 OCP 형상 또는 plateau/feature alignment 문제일 가능성이 높다. 반면 rate에 따라 증가하거나 charge/discharge에서 부호가 반전되는 오차는 R, Ds, k 또는 hysteresis 문제로 분류해야 한다.

저 SOC의 급격한 residual은 음극 OCP의 가파른 구간, cutoff 정의 및 충분하지 않은 relaxation에 민감하다. 고 SOC endpoint 오차는 양극 window 상단과 full-cell cutoff alignment에 민감하다.

## 권장 식별 순서

1. teardown으로 Ln, Lp, electrode area/layer count, porosity, coating mass와 active fraction을 확정한다.
2. 이 값으로 c_s,max, Qn, Qp와 측정 불확도를 계산한다.
3. chemistry-specific half-cell GITT OCP를 유지하고, full-cell 저율 데이터와 dV/dQ 특징으로 x0, x100, y100, y0 및 QLi를 식별한다.
4. charge/discharge branch별 OCP 또는 equilibrium OCP+hysteresis 모델을 비교한다.
5. OCP와 window를 고정한다.
6. teardown Rn/Rp를 반영하고 GITT로 Dsn/Dsp를 다시 계산한다.
7. EIS/pulse 및 0.5C–2C 데이터로 kn/kp와 저항/전해질 수송 파라미터를 식별한다.
8. 마지막에만 모든 파라미터를 측정 prior 안에서 제한적으로 joint refinement한다.

## 현재 선택

현재 기준 모델은 ±5% 해를 유지하되 이름을 '최종 물리 OCP'가 아니라 'teardown 전 provisional OCP window'로 변경한다. Chen 해는 물리 prior 후보로 유지한다. teardown 이후에는 ±5% 제한을 제거하고 측정된 Ln/Lp와 용량 불확도로 window를 재식별한다.

최종 OCP 선택은 qOCV MAE 최소값이 아니라 다음 순서로 결정한다.

1. protocol/chemistry validity
2. teardown capacity consistency
3. dV/dQ feature alignment
4. endpoint and hysteresis consistency
5. qOCV MAE/RMSE
6. held-out low-rate reproducibility

0.5C–2C 오차는 그 다음 단계의 R, Ds, k 식별 성능으로 평가한다.
