# 7월 최종 모델의 Ai2020 nominal one-at-a-time 복원 점검

## 결론

- 현재 7월 최종 모델에서 파라미터를 하나씩 Ai2020 nominal로 복원하면, 가장 큰 전압 오차 증가는 raw PyBaMM `De`에서 발생한다.
- 현재 `De`는 raw Ai2020 함수에 `1e-4`를 곱한다. 1000 mol/m3, 298.15 K에서 raw 함수는 `3.223e-6`, 현재값은 `3.223e-10 m2/s`다. 이 배율은 기존 검토에서 `cm2/s -> m2/s` 단위변환으로 판단됐으므로 raw 함수 직접 사용은 물리적 권장안이 아니라 진단 비교다.
- `Ln/Lp`, porosity, active fraction, `brugg_n/p/s`는 이미 Ai2020 nominal과 같아서 별도 simulation을 하지 않았다.
- `Dsp`, `kn`, `Rn` nominal 복원은 겹치는 전압구간의 RMSE를 낮추지만 charge가 실험시간 전에 cutoff되어 채택할 수 없다.
- 현재 OCP/window와 inventory를 유지한 채 dynamic parameter를 모두 nominal로 복원하면 charge/discharge RMSE가 `86.23/97.21 mV`로 악화한다.

## 기준

- 기준 모델: `Dsn=7.4662e-14`, `Dsp=6.8882e-14`, `kn=1.4502e-6`, `kp=5.8995e-7`, `brugg_n=2.914`
- 기준 physical-time RMSE: charge `54.00 mV`, discharge `52.30 mV`
- 기준 사후 capacity RMSE: charge `3.11%`, discharge `1.09%`
- 각 scenario는 한 항목만 nominal로 복원했다. capacity는 목적함수에 넣지 않았다.

## 주요 결과

| 복원 항목 | Charge RMSE | Discharge RMSE | Charge capacity RMSE | Discharge capacity RMSE | 판정 |
|---|---:|---:|---:|---:|---|
| Current July final | 54.00 mV | 52.30 mV | 3.11% | 1.09% | 기준 |
| De raw PyBaMM | 90.38 mV | 112.76 mV | 15.62% | 1.95% | 크게 악화 |
| Dsn nominal | invalid | 50.17 mV | 3.02% | 0.80% | charge 조기 cutoff |
| Dsp nominal | invalid | invalid | 4.19% | 2.70% | 양방향 조기 cutoff |
| kn nominal | invalid | 43.24 mV | 1.27% | 0.94% | charge 조기 cutoff |
| kp nominal | 66.49 mV | 67.98 mV | 5.32% | 1.29% | 양방향 악화 |
| Rn nominal | invalid | 44.10 mV | 1.64% | 0.76% | charge 조기 cutoff |
| Rp nominal | 63.28 mV | 64.35 mV | 4.99% | 1.27% | 양방향 악화 |
| Electrode width nominal | 63.10 mV | 84.31 mV | 7.45% | 4.81% | inventory/area 불일치로 악화 |
| csn_max nominal | 63.52 mV | 98.03 mV | 7.35% | 7.24% | discharge 크게 악화 |
| csp_max nominal | 56.04 mV | 55.73 mV | 4.08% | 1.27% | 소폭 악화 |
| Positive OCP nominal | 76.87 mV | invalid | 4.65% | 0.73% | charge 악화, discharge 조기 cutoff |
| All dynamic nominal | 86.23 mV | 97.21 mV | 11.28% | 0.38% | 전압 크게 악화 |

`invalid`는 적어도 한 셀·rate에서 모델이 실험 physical-time 끝보다 먼저 cutoff되어 정식 RMSE를 정의할 수 없다는 뜻이다. 그래프의 overlap RMSE 감소만 보고 해당 항목을 개선으로 판단하면 안 된다.

## 오차 증가의 중심

1. `De raw PyBaMM`: charge/discharge 전압 형상을 가장 크게 악화한다. 현재 `x1e-4` 처리를 제거하면 안 된다.
2. 음극 inventory/면적: `csn_max`, electrode width를 nominal로 단독 복원하면 특히 discharge와 capacity가 크게 악화한다. 현재 OCP window에 맞춘 Qn 보존과 묶여 있기 때문이다.
3. 양극 동특성: `kp`, `Rp`, positive OCP nominal 복원은 charge·discharge 전압 오차를 증가시킨다.
4. `Dsp`, `kn`, `Rn`: 일부 전압구간은 좋아지지만 charge cutoff capacity가 부족해진다. 이 항목들은 전압과 cutoff feasibility의 trade-off를 만든다.

## 산출물

- `nominal_reversion_error_attribution.png`: one-at-a-time 오차 변화
- `de_raw_pybamm_curve_comparison.png`: raw PyBaMM De와 현재 De 곡선 비교
- `nominal_reversion_summary.csv`, `nominal_reversion_detail.csv`: 전체 수치
- `current_vs_nominal_parameters.csv`: 현재값과 nominal 값 대조
