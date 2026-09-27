# 초기 SOC/OCP 심층 진단

## 결론

현재 적용 중인 **Ai2020 nominal graphite 음극 + old measured cathode GITT 양극의 equilibrium OCP**를 유지하는 것이 가장 타당하다.

다만 충전 시작부가 맞지 않는 주원인은 OCP 곡선의 선택 자체가 아니다. 1C와 2C에서는 현재 OCP와 연속 protocol simulation이 충전 직전 10분 rest 종점전압을 각각 약 `+3.8/-3.3 mV`로 재현했지만, 전류 인가 직후 전압 상승을 `214/286 mV`로 예측했다. 실험 상승량은 `75.8/137.3 mV`이다. 시작 오차의 대부분은 **저 SOC에서 과대 예측된 reaction overpotential**이다.

따라서 OCP를 다른 후보로 바꾸기보다 다음 단계에서 저 SOC exchange-current/kinetics와 연속 초기상태를 먼저 재검토해야 한다.

## 사용 데이터

- Full-cell OCP 선택: 260610 C/50 셀 1·2의 완전 cycle 2–5, charge 8개와 discharge 8개
- 동특성 진단: 260703 July BoL 셀 `6-5`, `6-6`, `6-8`의 0.5C/1C/2C charge
- 기존 음극/양극: `260902 Enertech셀 데이터 모음.xlsx`의 `ANODE_LOWRATE`, `CATHODE_GITT`
- 신규 음극/양극: 260918 GITT 각 3개 replicate
- 현재 동특성 파라미터: 선택된 보수 모델 `Dsn=7.4662e-14`, `kn=1.45024e-6`, `Dsp=6.88818e-14`, `brugg_n=2.914`

## 1. 충전 직전 rest는 평형 OCV가 아니다

모든 충전 branch 직전 rest는 10분이다. 마지막 구간에서도 전압이 계속 상승한다.

| 다음 charge | rest 종점 평균 | 마지막 rest slope 범위 | 실험 onset-rest 상승 |
|---:|---:|---:|---:|
| 0.5C | 3.0574 V | 2.33–2.72 mV/min | 44.9 mV |
| 1C | 3.2068 V | 5.52–5.95 mV/min | 75.8 mV |
| 2C | 3.3235 V | 6.20–6.69 mV/min | 137.3 mV |

평형 판정값으로 보기에는 slope가 매우 크다. 따라서 `rest_end_V -> qOCV inverse -> uniform electrode concentration` 방식은 물리적으로 완전한 초기상태가 아니다. 이 전압에는 직전 discharge에서 남은 입자 농도구배와 전해질 농도구배가 포함된다.

현재 OCP로 rest 전압을 역산한 평균 SOC는 `0.089/0.543/1.155%`이다. 셀별 충·방전량 bookkeeping으로 얻은 상대 SOC는 `0/1.043/2.200%`이다. SOC 차이는 1.05%p 이하이지만 저 SOC qOCV가 매우 가파르므로 1C와 2C에서 약 `100/124 mV`의 전압 차이로 증폭된다.

이 결과는 10분 rest 전압만으로 bulk SOC를 고유하게 정할 수 없다는 뜻이다. rest 전압 역산과 coulomb counting 어느 하나도 독립 branch의 균일 농도 초기조건으로 그대로 사용하면 안 된다.

## 2. 초기화 방식만 바꾼 대조실험

최종 `Dsn+kn` 파라미터를 고정하고 초기화만 변경했다.

| 초기화 | 시작점 RMSE | 전체 RMSE | 초기 0–10% RMSE | 중간 10–70% RMSE | 판단 |
|---|---:|---:|---:|---:|---|
| C-rate 평균 rest 역산, 현재 방식 | 134.7 mV | 54.00 mV | 74.13 mV | 55.71 mV | 기준 |
| 셀별 rest 역산 | 134.6 mV | 53.99 mV | 74.09 mV | 55.71 mV | 차이 없음 |
| SOC=0 endpoint | 85.6 mV | 52.19 mV | 39.02 mV | 57.21 mV | 수치상 일부 개선, 물리 불일치 |
| cycle coulomb bookkeeping | 201.7 mV | 57.34 mV | 103.88 mV | 54.56 mV | 균일 농도 적용 부적절 |
| 평균 rest + 양극 hysteresis | 111.6 mV | 56.03 mV | 65.03 mV | 60.03 mV | 시작만 개선, 전체 악화 |

셀별 rest 전압 사용은 평균 rest 사용과 `0.01 mV` 수준 차이뿐이다. 셀 편차가 원인이 아니다.

SOC=0 endpoint는 전체 RMSE를 약 `1.8 mV` 줄이지만, qOCV가 실제 rest보다 0.5C/1C/2C에서 평균 `47.6/197.0/313.6 mV` 낮아진다. 즉 시작 전압을 수치적으로 상쇄한 것이지 물리적인 초기상태 개선이 아니다.

## 3. 연속 protocol simulation

C/50 discharge부터 0.5C/1C/2C charge–CV–discharge–rest를 연속으로 simulation하여 입자·전해질 농도 이력을 보존했다.

| 다음 charge | 모델 rest 종점 | 실험 rest 종점 | 모델 charge onset | 실험 onset | onset 오차 |
|---:|---:|---:|---:|---:|---:|
| 0.5C | 3.0147 V | 3.0574 V | 3.1772 V | 3.1023 V | +74.9 mV |
| 1C | 3.2106 V | 3.2068 V | 3.4247 V | 3.2826 V | +142.1 mV |
| 2C | 3.3202 V | 3.3235 V | 3.6060 V | 3.4608 V | +145.3 mV |

1C와 2C rest 종점은 이미 4 mV 이내로 맞는다. 그럼에도 charge onset이 약 140–145 mV 높다. 따라서 이 두 조건의 시작 오차를 OCP 또는 SOC 초기화 오차로 설명하기 어렵다.

## 4. 전류 인가 순간 과전압 분해

연속 simulation에서 전류 인가 직후의 전압 상승을 분해했다.

| Charge | 실험 onset-rest | 모델 onset-rest | 음극 reaction 기여 절대값 | 양극 reaction 기여 | 전해질 ohmic 기여 |
|---:|---:|---:|---:|---:|---:|
| 0.5C | 44.9 mV | 162.5 mV | 89.1 mV | 63.4 mV | 10.0 mV |
| 1C | 75.8 mV | 214.1 mV | 101.3 mV | 93.3 mV | 19.4 mV |
| 2C | 137.3 mV | 285.9 mV | 123.3 mV | 125.5 mV | 37.0 mV |

실험 apparent step resistance는 약 `39.4/33.2/30.1 mOhm`, 모델은 `142.5/93.9/62.7 mOhm`이다. 모델이 약 `3.6/2.8/2.1배` 크다.

모델 reaction overpotential만으로도 실험 전체 전압상승보다 크다. 별도의 양의 series/contact resistance를 추가하면 시작 오차는 더 악화한다. 현재 우선순위는 저 SOC에서의 `kn`, `kp`, exchange-current 농도의존식, 그리고 OCP window가 만드는 endpoint stoichiometry가 서로 일관되는지 확인하는 것이다.

음극과 양극 reaction 기여가 비슷한 크기이므로 `kn` 하나만으로 시작부를 교정하기 어렵다. 다만 `kn`과 `kp`는 voltage-only fitting에서 강하게 상관되므로 둘을 자유롭게 동시에 fitting하기보다 EIS/pulse 근거 또는 하나의 공통 kinetic scale 진단이 필요하다.

## 5. 음극 OCP 비교

양극을 유효한 old cathode GITT로 고정한 결과다.

| 음극 | C/50 qOCV RMSE | C/50 charge/discharge 평균 MAE | 데이터 상태 | 판단 |
|---|---:|---:|---|---|
| Ai2020 nominal graphite | **6.24 mV** | **6.80 mV** | absolute stoichiometry 지원, directional branch 없음 | 적용 유지 |
| old measured anode | 8.22 mV | 9.08 mV | absolute 및 양방향 완전 | 보조 control |
| new 260918 anode | 16.42 mV | 13.07 mV | 양방향 형상은 있으나 absolute anchor 없음 | 적용 보류 |

old measured anode의 중앙 operating window에서 lithiation-delithiation gap은 약 `88.6 mV`인데 replicate 차이는 약 `2.5 mV`다. 그러나 full-cell C/50가 요구한 old/old hysteresis scale은 `0.260`뿐이다. 측정된 half-cell gap 전체를 full-cell graphite hysteresis로 적용하면 July 동특성이 악화했다.

따라서 현재 음극은 **Ai2020 nominal equilibrium OCP, hysteresis off**가 가장 안정적이다. old measured anode는 형상 불확실성 control로 남기고, new anode는 absolute stoichiometry/areal-capacity anchor가 확보될 때 재검토한다.

## 6. 양극 OCP 비교

음극을 Ai2020 nominal graphite로 고정한 결과다.

| 양극 | C/50 qOCV RMSE | C/50 charge/discharge 평균 MAE | 데이터 상태 | 판단 |
|---|---:|---:|---|---|
| old measured cathode GITT | 6.24 mV | **6.80 mV** | absolute 및 양방향 완전 | 적용 유지 |
| Ai2020 nominal cathode | 16.35 mV | 17.08 mV | 현재 셀 chemistry 형상 불일치 | 제외 |
| new 260918 cathode | **5.03 mV** | 14.80 mV | reverse branch 평균 73.6%만 완료, 마지막 rest 약 3.912 V | 진단용만 사용 |

new cathode는 equilibrium qOCV 숫자만 보면 가장 낮지만, reverse lithiation이 3.0 V까지 완료되지 않았고 absolute anchor도 없다. 또한 양방향 C/50 오차는 old cathode보다 두 배 이상 크다. 현재 데이터로는 채택할 수 없다.

old cathode의 equilibrium 형상은 유지한다. positive hysteresis를 켜면 시작점 RMSE는 `134.7 -> 111.6 mV`로 줄지만 전체 RMSE는 `54.0 -> 56.0 mV`, 중간 RMSE는 `55.7 -> 60.0 mV`로 악화한다. 따라서 현재 fitting baseline에서는 **양극 equilibrium OCP 사용, hysteresis off**가 타당하다.

## 최종 권고

1. OCP source는 변경하지 않는다.
   - 음극: `Ai2020 nominal graphite equilibrium OCP`
   - 양극: `old measured cathode GITT equilibrium OCP`
   - hysteresis: 음극/양극 모두 우선 off
2. branch마다 rest 전압을 균일 농도 SOC로 역산하는 방식을 최종 물리 모델로 해석하지 않는다. fitting은 가능하면 직전 discharge와 10분 rest를 포함한 연속 protocol 상태에서 수행한다.
3. 현재 50 mV대 fitting 한계의 첫 번째 교정 대상은 OCP 교체가 아니라 저 SOC reaction overpotential이다.
4. `Dsn+kn` fitting을 다시 하기 전에 시작 1–10초의 전극별 reaction overpotential, exchange-current density 및 EIS/pulse 일관성을 검증한다.
5. OCP를 새 측정값으로 교체하려면 새 양극 reverse GITT를 3.0 V까지 완료하고, 새 음극에는 absolute stoichiometry 또는 areal-capacity anchor를 추가한다.

## 파일

- `ocp_candidate_summary.csv`: 9개 음극/양극 조합과 데이터 완전성
- `ocp_candidate_initial_soc_detail.csv`: 후보별 rest 역산 SOC와 coulomb bookkeeping 비교
- `initialization_dynamic_overall.csv`: 초기화 방식 전체 비교
- `initialization_dynamic_by_rate.csv`: C-rate별 비교
- `continuous_protocol_steps.csv`: 연속 protocol simulation 단계별 결과
- `continuous_protocol_onset_components.csv`: 전류 인가 순간 과전압 분해
- `rest_non_equilibrium_diagnostics.png`: rest slope와 실험 onset jump
- `ocp_candidate_qocv_comparison.png`: 유효 OCP 조합 비교
- `initialization_scenario_comparison.png`: 초기화 방식별 시작/전체 RMSE

이 분석은 원본 Excel과 저장소 baseline을 변경하지 않은 read-only 진단이다.
