# Capacity-consistent Qn/Qp feasibility 재검토

## 결론

실제 각 branch의 CC 전류 중앙값으로 다시 계산했다. 이전 버전의 하드코딩 평균전류는 실제 CC 전류보다 특히 2C에서 약 2–3% 작았으므로 제거했다.

현재 grid에서 용량오차 합이 가장 작은 후보는 `Qn,Qp +4%`였지만, 사전에 정한 qOCV/동적 전압/용량 제약을 모두 통과한 후보는 없었다. 따라서 이 후보는 최종 물성값이 아니라 **용량 inventory를 늘렸을 때 얻을 수 있는 개선 상한을 보는 feasibility**로만 해석한다.

## 계산 조건

- `Qn*delta_x = Qp*delta_y = Qcell`을 C/50 window 단계에서만 강제
- 각 Qn/Qp 후보마다 `x0`, `y100`을 C/50 qOCV에 재최적화
- 동적 시험 전류: 각 0.5C/1C/2C branch의 실제 CC 전류 중앙값
- 초기 SOC: 각 시험 직전 10분 rest 종점전압의 qOCV inverse
- 음극 hysteresis off, 양극 history branch 사용
- `Dsn`, `Dsp`, `kn`, `kp`, Bruggeman, geometry는 현 baseline 유지

## baseline 대 용량 후보

| 지표 | baseline | Qn,Qp +4% feasibility |
|---|---:|---:|
| Qn / Qp | 2.4795 / 4.3712 Ah | 2.5787 / 4.5460 Ah |
| qOCV RMSE (2–98%) | 6.302 mV | 6.886 mV |
| 충전 10–70% 전압 RMSE | 19.99 mV | 20.49 mV |
| 방전 10–70% 전압 RMSE | 10.32 mV | 10.89 mV |
| 충전 용량 RMSE | 4.412% | 3.140% |
| 방전 용량 RMSE | 1.426% | 0.933% |
| 충전 평균 절대 용량오차 | 77.62 mAh | 55.27 mAh |
| 방전 평균 절대 용량오차 | 25.05 mAh | 14.63 mAh |

용량 inventory를 4% 키워도 충전 오차 약 55 mAh가 남고 qOCV와 충전 전압은 나빠졌다. 따라서 현재 용량오차를 `c_s,max` 또는 Qn/Qp 하나로 전부 흡수시키면 안 된다.

## k_n scan의 의미

Qn,Qp +4% 후보에서 `kn`을 키우면 충전 용량 RMSE는 감소했지만 전압 형상은 반대로 악화됐다.

| kn | 충전 용량 RMSE | 충전 10–70% 전압 RMSE |
|---:|---:|---:|
| 7.40e-7 | 3.140% | 20.49 mV |
| 8.50e-7 | 2.117% | 24.33 mV |
| 9.65e-7 | 1.282% | 28.53 mV |
| 1.075e-6 | 0.803% | 32.31 mV |
| 1.20e-6 | 0.911% | 36.23 mV |
| 1.35e-6 | 1.497% | 40.38 mV |

따라서 endpoint 용량만으로 `kn`을 선택하면 잘못된 보정이 된다. `kn/kp`는 사용자가 제공할 SOC별 Rct 기반 값으로 고정하고, 용량 inventory와 분리해 검증한다.

## 최종 판단

1. baseline Qn/Qp와 `c_s,max`를 임의 변경하지 않는다.
2. `Qn,Qp +4%`는 dry loading이 허용하는지 확인할 민감도/feasibility 결과로만 남긴다.
3. loading 측정 뒤 Qn, Qp prior를 만들고 C/50 qOCV/DVA에서 QLi와 window를 먼저 식별한다.
4. 0.5C/1C/2C endpoint 용량은 저율 inventory fitting 목표가 아니라 동역학 validation 지표로 둔다.
5. charge와 discharge의 전압 RMSE 및 endpoint 용량오차를 각각 보고하고, 임의 가중합으로 하나의 점수로 합치지 않는다.

재현 스크립트: `capacity_consistent_qn_qp_kn_refinement.py`
