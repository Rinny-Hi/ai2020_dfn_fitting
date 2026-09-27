# September BoL 원본 full-curve 복구 및 protocol QC

## 결론

- 6-1/6-2 원본 workbook의 `step`과 `record`를 읽어 0.5C/1C/2C charge·discharge CC 곡선을 모두 복구했다.
- 선택된 12개 branch에서 Neware `Total Time`은 모두 단조 증가하며, step 누락이나 잘못된 CC/CV 순서는 발견되지 않았다.
- 6-2 0.5C discharge의 마지막 1초 cutoff record는 `Date`와 step `Time`이 직전 record와 같게 표시되지만 `Total Time`에는 1초 차이가 있다. full physical-time 축은 이 record를 보존하는 `Total Time`을 사용한다.
- 각 branch 직전 rest는 120분이다. 고전압 rest 말단은 거의 평형이나, 저전압 charge 직전 rest는 0.5C와 1C에서 여전히 0.30/0.20 mV/min의 relaxation이 남는다.
- 2C charge는 CC보다 CV에서 들어가는 용량이 더 크다. 따라서 2C charge의 낮은 CC cutoff capacity를 cell inventory 감소로만 해석하면 안 된다.

## Branch 요약

| Rate | Direction | Rest tail slope | Load jump | CC time | CC capacity | CV capacity | CV fraction | Cell range |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 0.5C | Charge | +0.298 mV/min | +46.1 mV | 107.15 min | 2.0358 Ah | 0.2114 Ah | 9.4% | 17.2 mAh |
| 0.5C | Discharge | -0.009 mV/min | -36.9 mV | 118.09 min | 2.2438 Ah | - | - | 15.8 mAh |
| 1C | Charge | +0.197 mV/min | +87.0 mV | 46.13 min | 1.7531 Ah | 0.4624 Ah | 20.9% | 21.1 mAh |
| 1C | Discharge | -0.012 mV/min | -70.8 mV | 58.28 min | 2.2146 Ah | - | - | 13.5 mAh |
| 2C | Charge | +0.071 mV/min | +162.8 mV | 13.38 min | 1.0173 Ah | 1.1234 Ah | 52.5% | 21.7 mAh |
| 2C | Discharge | -0.016 mV/min | -141.0 mV | 28.11 min | 2.1366 Ah | - | - | 4.2 mAh |

## 해석

1. **초기 상태 품질**: discharge 직전 고전압 rest slope는 절대값 0.009--0.016 mV/min으로 안정적이다. charge 직전 저전압 rest는 2시간 후에도 0.5C 0.298, 1C 0.197, 2C 0.071 mV/min이다. July 10분 rest보다 훨씬 낫지만 0.5C/1C 저전압 상태를 완전 평형 OCP로 간주하지 않는다.
2. **전류 인가 응답**: 평균 rest-to-load jump는 charge에서 +46.1/+87.0/+162.8 mV, discharge에서 -36.9/-70.8/-141.0 mV다. rate에 따라 거의 증가하므로 OCP offset보다 kinetics·ohmic·transport 성분이 크다.
3. **CC/CV 분리**: charge CV 용량 비중은 0.5C 약 9.4%, 1C 약 20.9%, 2C 약 52.5%다. 2C CC 구간은 polarization 때문에 4.2 V에 일찍 도달하며, 이후 CV에서 약 1.12 Ah가 추가된다.
4. **CC→CV 전환 품질**: 두 셀·세 rate 모두 record gap은 0초이고 전압 불연속은 0.0--0.2 mV다. protocol 전환 문제로 볼 증거는 없다.
5. **셀 반복성**: 두 셀의 CC capacity 범위는 branch별 4--22 mAh 수준이다. cohort 내 반복성보다 July--September cohort 차이가 더 크다.
6. **fitting 입력**: physical-time charge fitting에는 CC record만 사용하고, CV capacity는 목적함수에 넣지 않는다. CV 및 total capacity는 post-hoc validation으로 남긴다.

## 이전 추출값 교차검증

이전 저장 QC와 raw rerun의 최대 차이는 `Rest_end_V=0.000e+00 V`, `rest slope=9.374e-14 mV/min`, `CC duration=3.333e-02 min`, `CC capacity=0.000e+00 Ah`이다. 기존 September 요약 추출은 재현된다.

## 산출물

- `september_full_cc_records.csv`: fitting용 CC full physical-time records
- `september_charge_cccv_records.csv`: charge CC+CV 연속 record
- `september_prebranch_rest_records.csv`: 각 branch 직전 120분 rest record
- `september_branch_qc.csv`, `september_protocol_summary.csv`: branch 및 cohort QC
- PNG 3종: full CC, rest-to-load, CC-CV 전환
