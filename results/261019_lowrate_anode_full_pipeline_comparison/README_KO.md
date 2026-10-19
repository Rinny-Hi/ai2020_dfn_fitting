# Ai2020 vs measured low-rate anode OCP through S3-S7a

양극 OCP, geometry, Rn/Rp, kp/Dsn/Dsp, ORegan2022 electrolyte, b_n=2.5, 실험 데이터와 오차 정의를 동일하게 유지했다. 각 음극 OCP에는 C/50에서 독립 선정한 window와 Qn/Qp, capacity-consistent cs,max를 사용했다.

## S3-S6 동일 변경 적용 결과

| Scenario | Overall_full_time_RMSE_mV | Overall_full_time_MAE_mV | Overall_center_time_RMSE_mV | Overall_center_time_MAE_mV | Capacity_RMSE_pct | Max_abs_capacity_error_mAh | Minimum_charge_anode_potential_mV | Stage | Anode_OCP | qOCV_RMSE_2_98_mV |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| S3 Ai2020 | 49.0325 | 41.1389 | 42.2536 | 38.5998 | 0.9624 | 35.9338 | -82.3918 | S3 | Ai2020 | 6.2427 |
| S3 Low-rate | 58.0552 | 51.2688 | 52.851 | 50.7682 | 3.2977 | 100.6336 | -74.1777 | S3 | Low-rate | 8.2196 |
| S4 Ai2020 | 46.2686 | 38.2546 | 38.8694 | 34.2723 | 3.2228 | 99.9781 | -55.7578 | S4 | Ai2020 | 6.2427 |
| S4 Low-rate | 51.817 | 43.9632 | 46.5876 | 42.7925 | 1.8749 | 94.7642 | -50.9098 | S4 | Low-rate | 8.2196 |
| S5 Ai2020 | 43.9037 | 35.6134 | 35.7445 | 32.8838 | 1.0752 | 30.9454 | -70.2805 | S5 | Ai2020 | 6.2427 |
| S5 Low-rate | 52.2044 | 45.5938 | 45.1596 | 43.5022 | 4.336 | 130.5264 | -61.9974 | S5 | Low-rate | 8.2196 |
| S6 Ai2020 | 38.5878 | 27.6823 | 26.3025 | 21.8436 | 1.0155 | 29.8283 | -53.9321 | S6 | Ai2020 | 6.2427 |
| S6 Low-rate | 44.119 | 35.6026 | 34.2431 | 32.0647 | 2.5038 | 66.8626 | -49.6803 | S6 | Low-rate | 8.2196 |

## S7a 동일 민감도/식별도 및 PSO fitting

- 저율 OCP 기준으로 선별된 fitting subset: `kn+brugg_n`
- Ai2020 S6→S7a 전체 RMSE: 38.59→37.15 mV
- 저율 S6→S7a 전체 RMSE: 44.12→37.69 mV

| Scenario | Overall_full_time_RMSE_mV | Overall_full_time_MAE_mV | Overall_center_time_RMSE_mV | Overall_center_time_MAE_mV | Capacity_RMSE_pct | Max_abs_capacity_error_mAh | Minimum_charge_anode_potential_mV |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Ai2020 S6 | 38.5878 | 27.6823 | 26.3025 | 21.8436 | 1.0155 | 29.8283 | -53.9321 |
| Ai2020 S7a | 37.1455 | 24.3829 | 20.9666 | 16.1645 | 1.9807 | 55.0905 | -60.172 |
| Low-rate S6 | 44.119 | 35.6026 | 34.2431 | 32.0647 | 2.5038 | 66.8626 | -49.6803 |
| Low-rate S7a | 37.6856 | 23.8188 | 18.4585 | 14.6704 | 1.8288 | 59.6468 | -70.3282 |

| Scenario | Fitted_parameters | Dsn_m2_s | kn_prefactor | brugg_n | Objective_RMSE_mV |
| --- | --- | --- | --- | --- | --- |
| Ai2020 S6 | none | 4.170774e-14 | 9.648533e-07 | 2.500000 | 38.933 |
| Ai2020 S7a | kn+brugg_n | 4.170774e-14 | 7.746510e-07 | 2.494938 | 37.437 |
| Low-rate S6 | none | 4.170774e-14 | 9.648533e-07 | 2.500000 | 44.665 |
| Low-rate S7a | kn+brugg_n | 4.170774e-14 | 5.955163e-07 | 2.607941 | 38.236 |

## 해석

- 저율 음극 OCP에도 이후 변경과 fitting을 동일하게 적용할 수 있다.
- 다만 terminal-voltage fitting이 OCP 차이를 kinetic/transport parameter로 보상할 수 있으므로, fitted RMSE만으로 음극 OCP의 물리적 우수성을 판단하지 않는다.
- 저율 OCP의 absolute x 환산은 Ai2020 c_s,max, active fraction, thickness에 의존하므로 독립적인 실험 anchor가 아니다.
- 일부 고율 분기에서 저율 OCP 보간 상한의 미세 외삽 경고가 발생했다.
- 음극 전위 절대값은 3전극 검증 전까지 후보 간 경향으로만 해석한다.

## 실험 저율 OCP의 x 범위가 좁은 이유

실험 저율 OCP에서 직접 측정한 값은 half-cell 전압과 용량이다. stoichiometry는 다음 식으로 환산했다.

`delta_x = Q_half-cell / Q_disk,assumed`

`Q_disk,assumed = F A_disk L_n epsilon_s,n c_s,max,n / 3600`

현재 구현은 직경 10 mm disk 면적과 Ai2020의 음극 두께, 활물질 부피분율, 최대 고체상 농도를 사용한다. 따라서 실험 저율 OCP의 좁은 x 범위는 전압 측정 자체의 오류라기보다 다음 항목이 합쳐진 **absolute stoichiometry calibration 불확실성**으로 해석한다.

- 실제 disk dry loading과 Ai2020 inventory의 차이
- 활물질 중량/부피분율, 코팅 두께, 유효면적의 불확실성
- half-cell cutoff, formation/SEI 및 조립 상태에 따른 가역용량 감소
- 채취 위치의 coating 편차와 teardown 영향
- 양극 OCP와 함께 full-cell qOCV를 맞추는 과정에서 발생한 window 보상

| Anode OCP | x0 | x100 | delta_x | Qn (Ah) | qOCV RMSE (mV) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Ai2020 + 현재 qOCV fitting | 0.003806 | 0.954221 | 0.950415 | 2.4795 | 6.2427 |
| 실험 저율 + 현재 qOCV fitting | 0.002074 | 0.714487 | 0.712413 | 3.3078 | 8.2196 |

좁은 `delta_x`가 실제로 더 큰 음극 여유용량을 의미할 가능성은 있으나, disk별 dry loading과 활물질 분율로 Qn을 독립 검증하기 전에는 실제 음극 사용범위로 확정하지 않는다.

## Ai2020 OCP 출처의 해석

Ai et al.은 LCO와 graphite OCV 곡선을 자체 GITT로 신규 측정했다고 기술하지 않고, Rieger et al.의 선행 문헌에서 가져왔다고 명시한다. Ai2020의 실험 검증은 Enertech SPB655060 full-cell의 0.5C, 1C, 2C 방전 전압·온도·두께 변화 비교였다.

Rieger et al.의 연구도 상용 pouch cell의 방전 전압, 온도, 변위와 pulse excitation 기반 열팽창 검증이 중심이다. 따라서 PyBaMM의 `graphite_ocp_Enertech_Ai2020.csv`는 **문헌 기반 equilibrium-potential lookup table**로 표현하고, Enertech 해체 전극을 Ai2020에서 GITT로 직접 측정한 결과라고 표현하지 않는다.

- Ai et al. (2020): https://doi.org/10.1149/2.0122001JES
- Rieger et al. (2016): https://doi.org/10.1149/2.0441608jes
- PyBaMM Ai2020 implementation: https://github.com/pybamm-team/PyBaMM/blob/main/packages/pybamm/src/pybamm/input/parameters/lithium_ion/Ai2020.py

## 현재 선정

현재 기준모델에는 **실험 GITT 양극 + Ai2020 graphite 음극 OCP**를 유지한다. 이는 Ai2020 OCP가 현차셀의 직접 측정값이기 때문이 아니라, 현재 확보된 후보 중 qOCV와 S3–S6 동적 전압 재현성이 더 좋고 절대 x축의 추가적인 조건부 환산이 필요하지 않기 때문이다. 실험 저율 OCP는 dry loading 기반 Qn과 절대 stoichiometry anchor가 확보되면 다시 평가한다.
