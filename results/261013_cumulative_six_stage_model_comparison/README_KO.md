# Ai2020에서 현재 기준모델까지 누적 6단계 비교

## 비교 정의

- **S1 Ai2020 raw**: Raw PyBaMM Ai2020 parameters, nominal OCPs, and published window
- **S2 S1 + corrected De**: Only Ai2020 electrolyte diffusivity is multiplied by 1e-4
- **S3 S2 + selected OCP/window**: Ai2020 graphite + existing LCO GITT; C/50 window and capacity-consistent inventory
- **S4 S3 + measured Rn/Rp**: Rn=3.7542 um and Rp=4.2387 um
- **S5 S4 + kp/Dsn/Dsp**: Measured-anchor kp=5.8995e-7, Dsn=4.1708e-14, Dsp=6.8882e-14
- **S6 S5 + ORegan2022/bn=2.5**: ORegan2022 electrolyte transport and negative-electrode Bruggeman 2.5

S3에서 OCP/window와 함께 effective cs,max 및 공통 overlap 면적을 적용했다. 이는 선택된 window의 Qn/Qp를 실제 DFN inventory로 구현하기 위한 연동 bookkeeping이며 동적 fitting은 아니다.

모든 단계는 Cell 6-6 July-BOL 0.5C/1C/2C 충·방전 실험의 각 rest-end 전압에서 독립 초기화했으며, 별도 fitting 없이 4.2/3.0 V cutoff까지 계산했다.

## 종합 결과

| 단계 | 성공/6 | 전체 RMSE/MAE | 10–70% RMSE/MAE | 용량 RMSE | 최대 용량 오차 |
|---|---:|---:|---:|---:|---:|
| S1 Ai2020 raw | 6/6 | 118.20/105.14 mV | 112.31/101.76 mV | 14.40% | 444.4 mAh |
| S2 S1 + corrected De | 6/6 | 76.28/63.85 mV | 72.72/60.61 mV | 7.92% | 214.1 mAh |
| S3 S2 + selected OCP/window | 6/6 | 49.03/41.14 mV | 42.25/38.60 mV | 0.96% | 35.9 mAh |
| S4 S3 + measured Rn/Rp | 6/6 | 46.27/38.25 mV | 38.87/34.27 mV | 3.22% | 100.0 mAh |
| S5 S4 + kp/Dsn/Dsp | 6/6 | 43.90/35.61 mV | 35.74/32.88 mV | 1.08% | 30.9 mAh |
| S6 S5 + ORegan2022/bn=2.5 | 6/6 | 38.59/27.68 mV | 26.30/21.84 mV | 1.02% | 29.8 mAh |

## C-rate별 전체 시간 전압 RMSE

아래 값은 `충전 / 방전` 순서다.

| 단계 | 0.5C | 1C | 2C |
|---|---:|---:|---:|
| S1 Ai2020 raw | 107.0 / 72.6 mV | 126.6 / 95.0 mV | 157.0 / 131.9 mV |
| S2 S1 + corrected De | 89.9 / 60.4 mV | 93.7 / 62.3 mV | 95.3 / 37.9 mV |
| S3 S2 + selected OCP/window | 47.2 / 47.6 mV | 53.7 / 54.5 mV | 54.5 / 33.3 mV |
| S4 S3 + measured Rn/Rp | 47.0 / 46.5 mV | 50.3 / 48.6 mV | 47.7 / 36.0 mV |
| S5 S4 + kp/Dsn/Dsp | 42.3 / 43.6 mV | 48.4 / 48.4 mV | 49.9 / 26.4 mV |
| S6 S5 + ORegan2022/bn=2.5 | 38.8 / 39.5 mV | 43.5 / 40.3 mV | 44.7 / 18.9 mV |

## 해석

- S1→S2에서 De 단위 보정만으로 전체 전압 RMSE가 118.20→76.28 mV, 용량 RMSE가 14.40→7.92%로 감소했다.
- S2→S3의 선택 OCP/window 및 연동 inventory 반영이 가장 큰 추가 개선으로, 전체 전압 RMSE는 49.03 mV, 용량 RMSE는 0.96%가 됐다.
- S3→S4의 측정 입자반경 반영은 전압 RMSE를 소폭 줄였지만 용량 RMSE를 3.22%로 일시 증가시켰다. 입자반경이 분극과 cutoff 도달시점에 함께 영향을 준 결과다.
- S4→S5의 kp/Dsn/Dsp 반영은 전압 및 용량 오차를 다시 감소시켰다.
- S6가 전체 전압 RMSE 38.59 mV, 중앙 RMSE 26.30 mV로 가장 우수하며 용량 RMSE도 1.02%로 유지된다. 기존 현재 기준모델 결과와 수치상 일치한다.

## 주의

- 전체 시간 RMSE는 실험과 모델이 모두 존재하는 공통 경과시간 구간에서 계산했다. 종료용량 차이는 용량 RMSE와 최대 용량 오차로 별도 평가했다.
- 단계 간 차이는 누적 변경의 효과이며 개별 파라미터의 독립 인과효과로 해석하지 않는다.
- S1의 raw Ai2020 De와 S2의 corrected De 차이는 Ai2020 구현의 단위 해석 검증을 위한 비교다.
- S4 1C 방전 계산 중 선택된 양극 OCP 보간 상한의 미세 외삽 경고가 1회 발생했다. 계산은 정상 종료됐으며 최종 S6에서는 발생하지 않았다.
