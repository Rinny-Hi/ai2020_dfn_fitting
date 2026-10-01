# C/50 기반 old/new/nominal OCP 조합 비교

## 결론

현재 채택 후보는 **Ai2020 nominal 음극 + 기존(old) 양극 GITT**이다.

- C/50 qOCV MAE, SOC 2–98%: **4.69 mV**
- C/50 qOCV RMSE, SOC 2–98%: **6.24 mV**
- C/50 charge MAE, SOC 10–90%: **9.12 mV**
- C/50 discharge MAE, SOC 10–90%: **4.48 mV**
- 충·방전 평균 MAE: **6.80 mV**
- 최대 endpoint 오차: **10.00 mV**

새 양극을 사용한 `nominal/new` 조합은 qOCV RMSE가 5.03 mV로 더 낮지만, 새 양극 역방향 GITT가 평균 73.6%만 진행되었고 마지막 이완전압이 약 3.912 V이다. 따라서 낮은 RMSE만으로 최종 채택하지 않았다.

## 비교 기준

- Full-cell 목표: 260610 C/50 셀 1·2의 완전 cycle 2–5
- 사용 branch: charge 8개, discharge 8개
- 평균 full-cell 용량: 2.35653 Ah
- OCP window fitting: SOC 2–98% qOCV MSE 최소화
- 제약조건: SOC 0/100 endpoint 오차 각각 10 mV 이내
- 전극용량 허용범위: 2.4–6.0 Ah
- 충·방전 분기 scale은 equilibrium window를 고정한 후 0–1 범위에서 별도 fitting

## 주요 후보

| 음극 | 양극 | qOCV MAE (mV) | qOCV RMSE (mV) | C/50 branch 평균 MAE (mV) | Hysteresis scale | 판단 |
|---|---|---:|---:|---:|---:|---|
| nominal | old | 4.69 | 6.24 | 6.80 | 1.000 | 채택 후보 |
| old | old | 6.40 | 8.22 | 9.08 | 0.260 | 완전 측정 control |
| old | nominal | 10.75 | 12.50 | 11.54 | 0.311 | qOCV 열세 |
| nominal | nominal | 13.71 | 16.35 | 17.08 | 1.000 | 열세 |
| new | old | 13.29 | 16.42 | 13.07 | 0.903 | 음극 절대 anchor 없음 |
| nominal | new | 3.66 | 5.03 | 14.80 | 1.000 | 새 양극 미완료로 제외 |

## 선정 window

| 항목 | 값 |
|---|---:|
| x0 | 0.003806 |
| x100 | 0.954221 |
| y100 | 0.431642 |
| y0 | 0.970746 |
| delta x | 0.950415 |
| delta y | 0.539103 |
| Qn | 2.47947 Ah |
| Qp | 4.37120 Ah |

SOC 0 endpoint가 10 mV 제약 경계에 있으므로 최종 확정 전 새 BoL 데이터에서 endpoint를 다시 확인한다.

## 새 GITT QC

### 음극

- v1 charge/discharge pulse capacity: 약 1.68/1.68 mAh
- v2: 약 1.76/1.76 mAh
- v3: 약 1.80/1.84 mAh
- v2/v3 중앙 형상은 사용할 수 있으나 절대 stoichiometry anchor와 셀 간 용량 차이가 남아 있다.

### 양극

- delithiation capacity: 2.80–2.84 mAh
- reverse lithiation capacity: 2.04–2.12 mAh
- reverse branch completion: 71.8–75.7%
- 마지막 이완전압: 3.911–3.913 V

따라서 새 양극은 완전한 양방향 OCP 자료가 아니다.

## 추가 실험 판단

현재 OCP 모델링과 BoL fitting을 진행하기 위한 추가 저율 full-cell 시험은 필수가 아니다. 기존 C/50 두 셀의 재현성이 충분하고, nominal 음극 + old 양극 조합이 qOCV와 양방향 C/50 곡선을 모두 잘 재현한다.

다만 다음 목적에는 추가 실험이 필요하다.

1. 새 양극 GITT를 최종 OCP로 채택하려면 reverse lithiation을 3.0 V까지 완료하고 각 pulse 뒤 2 h rest를 확보한다.
2. 새 음극을 절대 stoichiometry 기준으로 채택하려면 세척·건조 조건을 통일하고 dry coating mass 또는 areal loading을 확보한다.
3. 선정 OCP의 DFN 적용 확정에는 동일 초기상태로 측정한 새 0.5C/1C/2C BoL 데이터가 필요하다.

## 즉시 실행할 BoL 프로토콜

| 순서 | 단계 | 조건 | 종료조건 | Rest | 목적 |
|---:|---|---|---|---:|---|
| 0 | 온도 안정화 | 25 ± 1 °C, 무부하 | 2 h | - | 열·전압 안정화 |
| 1 | 초기 상태 정렬 | 0.5C discharge | 3.0 V | 60 min | 첫 charge 초기 SOC 확정 |
| 2 | 0.5C charge | CC 0.5C → CV 4.2 V | CV current ≤ C/20 | 60 min | charge fitting |
| 3 | 0.5C discharge | CC 0.5C | 3.0 V | 60 min | discharge validation |
| 4 | 1C charge | CC 1C → CV 4.2 V | CV current ≤ C/20 | 60 min | charge fitting |
| 5 | 1C discharge | CC 1C | 3.0 V | 60 min | discharge validation |
| 6 | 2C charge | CC 2C → CV 4.2 V | CV current ≤ C/20 | 60 min | charge fitting |
| 7 | 2C discharge | CC 2C | 3.0 V | 60 min | discharge validation |

- 동일 셀을 사용하고 0.5C → 1C → 2C 순서로 진행한다.
- CC/CV 기록 간격은 1 s, rest는 5–10 s로 한다.
- 각 rest의 시작·종료 전압을 보존한다.
- 예상 총시간은 약 17–20 h이다.

## 결과 파일

- `ocp_source_combination_summary.csv`: 9개 조합과 window/오차
- `new_gitt_capacity_qc.csv`: 새 양·음극 pulse capacity QC
- `new_cathode_completion_qc.csv`: 새 양극 역방향 completion QC
- `ocp_source_combination_heatmaps.png`: 조합별 qOCV/branch 오차
- `top_ocp_source_curves.png`: C/50 곡선 비교
