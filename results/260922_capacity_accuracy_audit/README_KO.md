# 용량 계산 정확도 및 선행연구 기반 개선방향 재검토

## 1. 최종 판단

현재 용량 계산식, 총 활성면적, 병렬 전극수, 실제 시험전류, 적산용량 처리에는 큰 계산 오류가 없다. 남은 오차는 하나가 아니라 다음 두 성분으로 분리된다.

- **정적 성분:** OCP reference/window, Qn/Qp, cyclable lithium `QLi`, formation loss, 초기 SOC
- **동적 성분:** Rct/kinetics, solid/electrolyte transport, ohmic/contact 저항 때문에 cutoff에 일찍 도달하는 양

현재 결과에서 충전 용량오차는 0.5C/1C/2C 모두 약 `-77/-77/-79 mAh`로 거의 일정하다. 이는 단순한 C-rate 의존 kinetics보다 먼저 **branch 공통의 window/QLi/OCP offset**을 의심해야 한다는 패턴이다. 방전은 `-8/-15/-52 mAh`로 C-rate에 따라 커지므로 **동적 polarization** 성분이 뚜렷하다.

따라서 `Dsn`, `Dsp`, `kn`, `kp`를 용량에 맞춰 움직이지 않는다. Rct 기반 `kn/kp`가 들어오면 고정하고, loading 기반 Qn/Qp와 저율 전압 기반 QLi/window를 별도로 식별하는 것이 가장 근거가 확실한 방향이다.

## 2. 계산 교차검산

PyBaMM의 전극 최대 용량과 동일한 식을 직접 계산했다.

`Q_k = F * A_face * N_parallel * L_k * epsilon_s,k * c_s,max,k / 3600`

| 항목 | 직접 계산 | 모델 값 | 차이 |
|---|---:|---:|---:|
| C/50 target Qcell | 2.3565277 Ah | 2.3565277 Ah | 0 |
| Qn | 2.4794722 Ah | 2.4794722 Ah | 0 |
| Qp | 4.3711981 Ah | 4.3711981 Ah | 0 |

여기서 `A_face=51×47 mm²`, `N_parallel=34`이다. 한 면 면적만 쓰고 병렬 전극수를 누락하면 Q가 1/34로 계산되므로 둘을 반드시 함께 관리한다.

실제 시험 전류도 2.28 Ah nominal C-rate setpoint와 최대 0.005% 이내로 일치했다.

| 조건 | 실제 CC 전류 | 2.28 Ah 기준 setpoint |
|---|---:|---:|
| 0.5C charge/discharge | 1.139948 / 1.139974 A | 1.140 A |
| 1C charge/discharge | 2.279965 / 2.279931 A | 2.280 A |
| 2C charge/discharge | 4.559894 / 4.560041 A | 4.560 A |

`Nominal cell capacity`는 C-rate를 전류로 바꾸는 값일 뿐 실제 전극 inventory가 아니다. 실제 용량은 위 Qn/Qp, QLi, OCP window와 cutoff dynamics가 결정한다.

적산용량도 이제 clipped SOC를 역산하지 않고 simulation/experiment의 원래 `Q_Ah`를 직접 사용한다. 이번 데이터는 clipping 범위를 넘지 않았기 때문에 수치 결과는 바뀌지 않았지만, 이후 후보에서 용량이 qcell을 넘을 때 생길 수 있는 숨은 cap을 제거했다.

## 3. 실제 endpoint 용량오차

| C-rate | charge | discharge |
|---:|---:|---:|
| 0.5C | -76.92 mAh | -8.23 mAh |
| 1C | -76.61 mAh | -14.82 mAh |
| 2C | -79.33 mAh | -52.12 mAh |

모두 `model - experiment`이다. 충전의 일정한 offset과 방전의 rate dependence를 같은 파라미터 하나로 해결하려 하면 전압 형상을 훼손한다.

## 4. 초기 rest 불확실성의 크기

모든 동적 branch 직전 rest는 10분이었다. charge 시작 전 rest의 마지막 3분 기울기는 약 `2.88/6.01/6.56 mV/min`이라 완전 평형으로 보기 어렵다. qOCV-equivalent SOC는 charge `0.074/0.472/1.094%`, discharge `98.87/98.84/98.91%`였다.

이를 극단적으로 정확한 0/100%로 보정해도 가능한 용량 이동량은 대략 charge `1.8–25.8 mAh`, discharge `25.8–27.3 mAh` 수준이다. 따라서 초기 SOC 불확실성은 일부 오차를 설명하지만 0.5C와 1C charge의 약 77 mAh 전체를 설명하지 못한다.

원본 RPT 순서를 rest까지 포함해 연속 DFN으로 모사한 기존 검토에서도 charge/discharge 용량 RMSE가 `6.17/3.18%`였고 충전 부족이 해소되지 않았다. 즉 branch 독립 초기화와 이력 누락은 오차 원인 중 하나지만 주원인 전체는 아니다.

## 5. capacity inventory 조정 feasibility

실제 CC 전류로 Qn/Qp grid를 다시 계산했다.

| 지표 | baseline | Qn,Qp +4% feasibility |
|---|---:|---:|
| qOCV RMSE | 6.302 mV | 6.886 mV |
| charge 용량 RMSE | 4.412% | 3.140% |
| discharge 용량 RMSE | 1.426% | 0.933% |
| charge 평균 절대오차 | 77.62 mAh | 55.27 mAh |
| discharge 평균 절대오차 | 25.05 mAh | 14.63 mAh |

용량은 개선되지만 qOCV와 charge 전압은 악화되고, 사전 제약을 모두 통과한 후보는 없었다. 이 결과는 loading 측정 전 `c_s,max`를 +4% 하라는 뜻이 아니라, **inventory 오차가 실제 원인의 일부일 가능성**만 보여준다.

`kn` scan에서는 `kn=1.075e-6`일 때 charge 용량 RMSE가 0.80%까지 내려갔지만 charge 전압 RMSE는 20.49→32.31 mV로 악화됐다. 용량을 맞추기 위해 kinetics를 조정하면 안 된다는 정량적 반례다.

## 6. 선행연구와 PyBaMM/PyBOP가 제시하는 방향

- PyBaMM eSOH는 `Qn`, `Qp`, `QLi`, 전압 상·하한으로 `x0,x100,y0,y100`과 usable cell capacity를 함께 푼다. 즉 window 폭을 Qcell에만 맞추는 것보다 `QLi`를 독립 상태량으로 두는 것이 정석에 가깝다.
- PyBOP의 OCP 방법은 charge-capacity/stoichiometry domain에서 RMSE를 기본 cost로 사용하고, OCP branch 정렬에 shift와 stretch를 사용한다. 용량 endpoint를 임의 가중치로 voltage SSE에 더하는 것이 필수 표준은 아니다.
- Weng 등의 DVA 방법은 저율 full-cell 전압과 half-cell reference로 electrode capacities와 cyclable lithium/formation loss를 분리한다. 또한 제한된 half-cell OCP 구간 때문에 fitted capacity가 true design capacity와 다를 수 있는 `inaccessible lithium` 문제를 명시한다.
- Chen2020은 OCP/stoichiometry를 GITT와 three-electrode 자료로, kinetics를 EIS로, geometry/조성을 tear-down으로 각각 측정한 뒤 P2D를 검증한다. 서로 다른 오차를 한 종류의 시험으로 모두 흡수시키지 않는 구조다.

## 7. 권장 fitting/validation 순서

1. **데이터 정리:** 각 branch의 실제 전류를 사용하고, `Q_Ah=∫|I|dt/3600`을 원본 그대로 유지한다.
2. **loading prior:** 양·음극 dry coating loading, 활물질 분율, 실제 면적/병렬 face 수로 Qn/Qp의 측정 prior와 불확도를 만든다.
3. **저율 eSOH/DVA:** C/50 charge와 discharge를 capacity domain에서 각각 확인하고, Qn/Qp prior 안에서 `Qn,Qp,QLi` 또는 동등한 window shift/stretch를 fitting한다. 전압 RMSE와 dV/dQ feature mismatch를 별도 보고한다.
4. **OCP admissibility:** qOCV RMSE, endpoint 오차, stoichiometry bounds, practical NPR을 제약으로 사용한다. 임의의 0.7/0.3 가중합으로 후보를 한 줄 세우지 않는다.
5. **동역학:** Rct 기반 kn/kp를 고정하고 0.5C/1C/2C 전압 residual의 rate/direction pattern으로 Ds, electrolyte transport, contact/ohmic 성분을 검증한다.
6. **검증 분리:** charge fitting을 한다면 discharge는 validation으로 남기고, rate별 full-range/10–70% voltage RMSE와 endpoint capacity error를 각각 표시한다.
7. **불확도 보고:** loading/OCP/초기 SOC 허용범위에서 predicted capacity band를 만들고, 측정 곡선이 band 밖인지 확인한다.

## 8. 지금 바꾸지 않을 값과 다음 입력

- `Dsn`, `Dsp`: 소재 기반값 유지
- `kn`, `kp`: 임의 scan값을 채택하지 않고 사용자가 제공할 Rct 기반값 적용
- `Rn`, `Rp`, `Ln`, `Lp`: 이번 SEM 결과는 feasibility 최종본으로 공유하되 baseline 강제 변경 근거로 쓰지 않음
- `c_s,max`: loading prior 전에는 현 값 유지
- 다음 핵심 입력: dry coating loading/활물질 분율과 SOC별 Rct

## 9. 근거 자료

- PyBaMM electrode SOH: https://docs.pybamm.org/en/latest/source/examples/notebooks/models/electrode-state-of-health.html
- PyBaMM capacity discussion: https://github.com/pybamm-team/PyBaMM/discussions/4211
- PyBaMM capacity/geometry discussion: https://github.com/pybamm-team/PyBaMM/discussions/1635
- PyBOP OCP methods: https://pybop-docs.readthedocs.io/en/latest/_modules/pybop/applications/ocp_methods.html
- Weng et al., DVA for manufacturing: https://doi.org/10.3389/fenrg.2023.1087269
- Chen et al. 2020 parameterisation: https://doi.org/10.1149/1945-7111/ab9050

재현 파일:

- `capacity_accuracy_audit.py`
- `capacity_consistent_qn_qp_kn_refinement.py`
- `capacity_error_diagnosis.png`
- `capacity_formula_crosscheck.csv`
- `measured_current_crosscheck.csv`
