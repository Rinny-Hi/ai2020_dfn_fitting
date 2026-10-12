# OCP 조합별 동적 교차검증

## 고정 조건

- Cell 6-6 July-BOL 0.5C/1C/2C 충·방전 6개 CC branch를 사용했다.
- ORegan2022 전해질 수송식과 음극 Bruggeman 계수 `b_n=2.5`를 사용했다.
- geometry, Dsn/Dsp, kn/kp는 현재 no-fit 기준값으로 고정했다.
- 각 OCP 조합에는 해당 조합의 C/50 fitted stoichiometry window, Qn/Qp, cs,max를 적용했다.
- 동적 파라미터는 OCP 조합별로 다시 fitting하지 않았다.
- 전압 RMSE는 공통 경과시간 구간에서 계산하고 용량 오차는 cutoff 종료용량으로 별도 평가했다.

## 종합 결과

| OCP 조합 | 신뢰 구분 | 성공/6 | 전체 시간 RMSE | 중앙 시간 RMSE | 용량 RMSE | 최대 용량 오차 | 충전 최저 음극전위 |
|---|---|---:|---:|---:|---:|---:|---:|
| New GITT anode + existing GITT cathode | sensitivity-only | 6/6 | 36.71 mV | 27.65 mV | 1.46% | 34.5 mAh | -54.9 mV |
| Ai2020 anode + existing GITT cathode | decision-grade comparator | 6/6 | 38.59 mV | 26.30 mV | 1.02% | 29.8 mAh | -53.9 mV |
| Ai2020 anode + new GITT cathode | sensitivity-only | 6/6 | 39.39 mV | 20.22 mV | 2.97% | 76.5 mAh | -44.7 mV |
| Existing low-rate anode + existing GITT cathode | decision-grade comparator | 6/6 | 44.12 mV | 34.24 mV | 2.50% | 66.9 mAh | -49.7 mV |
| Ai2020 anode + Ai2020 cathode | decision-grade comparator | 6/6 | 46.63 mV | 39.15 mV | 1.51% | 41.4 mAh | -60.6 mV |
| Existing low-rate anode + Ai2020 cathode | decision-grade comparator | 6/6 | 48.36 mV | 39.47 mV | 1.51% | 52.9 mAh | -78.7 mV |

## 판단

- 절대 stoichiometry가 확보된 조합 중 `Ai2020 anode + existing GITT cathode`가 전체 시간 RMSE, 중앙 시간 RMSE, 용량 RMSE, 최대 용량 오차를 종합했을 때 가장 우수했다.
- `existing low-rate anode + existing GITT cathode`는 C/50 qOCV는 경쟁력이 있었지만 동적 검증에서 전체 시간 RMSE 44.12 mV, 용량 RMSE 2.50%, 최대 용량 오차 66.9 mAh로 현재 조합보다 열세였다.
- `new GITT anode + existing GITT cathode`는 전체 시간 RMSE가 36.71 mV로 가장 낮았지만 절대 stoichiometry anchor가 없고 용량 RMSE도 현재 조합보다 커서 최종 선정 근거로 사용하지 않는다.
- 따라서 현재 OCP인 `Ai2020 nominal graphite + existing measured LCO GITT`를 유지한다.
- 실행 중 기존 저율 음극 OCP 조합 일부에서 측정 OCP 보간범위 상한 extrapolation 경고가 발생했다. 이는 고율 동적 상태에서 표면 stoichiometry가 측정 support를 벗어날 수 있음을 뜻하므로 추가적인 물리적 주의사항이다.

## 해석 제한

- `new GITT anode`는 절대 stoichiometry anchor가 확립되지 않아 sensitivity-only이다.
- `new GITT cathode`는 reverse lithiation branch가 미완결이므로 sensitivity-only이다.
- 동적 전압 RMSE가 낮더라도 측정 완결성과 절대 stoichiometry가 부족한 조합을 최종 OCP로 선정하지 않는다.
- OCP 변경에 따른 음극전위 차이는 모델 예측 민감도이며 3전극 측정에 의한 절대 검증값이 아니다.

## C/50 입력 window

| OCP 조합 | x0 | x100 | y0 | y100 | Qn | Qp | qOCV RMSE | primary eligible |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| Ai2020 anode + existing GITT cathode | 0.003806 | 0.954221 | 0.970746 | 0.431642 | 2.4795 Ah | 4.3712 Ah | 6.24 mV | True |
| Existing low-rate anode + existing GITT cathode | 0.002074 | 0.714487 | 0.976235 | 0.448271 | 3.3078 Ah | 4.4634 Ah | 8.22 mV | True |
| Existing low-rate anode + Ai2020 cathode | 0.001675 | 0.840542 | 0.920101 | 0.452709 | 2.8092 Ah | 5.0419 Ah | 12.50 mV | True |
| Ai2020 anode + Ai2020 cathode | 0.003180 | 0.985066 | 0.867768 | 0.441824 | 2.4000 Ah | 5.5325 Ah | 16.35 mV | True |
| New GITT anode + existing GITT cathode | 0.004843 | 0.986729 | 0.975692 | 0.429007 | 2.4000 Ah | 4.3106 Ah | 16.42 mV | True |
| Ai2020 anode + new GITT cathode | 0.002955 | 0.941869 | 0.935516 | 0.070386 | 2.5098 Ah | 2.7239 Ah | 5.03 mV | False |
