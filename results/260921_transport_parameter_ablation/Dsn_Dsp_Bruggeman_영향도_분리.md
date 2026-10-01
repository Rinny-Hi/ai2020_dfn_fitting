# Dsn, Dsp, Bruggeman 영향도 분리

## 고정 조건

- OCP: Ai2020 nominal equilibrium + 측정 hysteresis
- Stoichiometry window 및 Qn/Qp: 이전 질량 비의존 fitting 결과 고정
- Rn = 5 µm, Rp = 3 µm 고정
- 0.5C, 1C, 2C charge/discharge held-out 검증
- Dsn, Dsp, electrolyte Bruggeman의 2^3=8개 모든 조합 계산

## 단독 변경 결과

| 조건 | 동적 MAE (mV) | 기존 대비 변화 (mV) | Charge MAE | Discharge MAE | 평균 절대 용량오차 (mAh) |
|---|---:|---:|---:|---:|---:|
| 기존 Ai2020 transport | 27.09 | 0.00 | 27.99 | 26.19 | 37.77 |
| Dsn만 2.1e-14로 변경 | 25.12 | -1.97 | 26.21 | 24.03 | 39.70 |
| Dsp만 4.4e-14로 변경 | 33.84 | +6.75 | 33.60 | 34.08 | 35.76 |
| Bruggeman만 1.5로 변경 | 73.09 | +46.00 | 65.47 | 80.71 | 62.94 |
| Dsn+Dsp+Bruggeman 변경 | 77.73 | +50.64 | 69.88 | 85.57 | 78.00 |

## 상호작용을 포함한 정확한 기여도

모든 변경 순서를 평균한 exact Shapley allocation 결과:

| 항목 | 동적 MAE 변화 기여 |
|---|---:|
| Dsn | -2.31 mV |
| Dsp | +6.81 mV |
| Bruggeman | +46.14 mV |
| 합계 | +50.64 mV |

양수는 오차 악화, 음수는 오차 개선을 의미한다. Bruggeman 변경이 전체 악화의 약 91%를 설명한다.

## Bruggeman이 큰 이유

Ai2020의 electrolyte Bruggeman coefficient는 negative 2.914, positive 1.83, separator 1.5다. 모두 1.5로 통일하면 실제로 변경되는 것은 negative와 positive이며 separator는 변하지 않는다.

Effective electrolyte transport는 대략 epsilon^b에 비례한다. Ai2020 porosity negative=0.33, positive=0.32, separator=0.50를 적용하면:

- negative: b 2.914→1.5에서 effective transport 약 4.80배 증가
- positive: b 1.83→1.5에서 effective transport 약 1.46배 증가
- separator: b=1.5로 동일하여 변화 없음

따라서 electrolyte concentration gradient와 ohmic polarization이 크게 줄어든다. 현재 실험은 baseline 모델보다 더 큰 분극을 요구하므로, Bruggeman을 1.5로 낮추면 charge 전압은 너무 낮고 discharge 전압은 너무 높아져 오차가 증가한다.

## 판단

전압 오차가 크게 증가한 주원인은 Bruggeman 통일이다. Dsp 증가는 두 번째로 작지만 분명한 악화 요인이고, Dsn 감소는 오히려 약 2 mV 개선 방향이었다. 따라서 Dsn/Dsp 값 자체를 먼저 폐기할 필요는 없으며, Bruggeman을 Ai2020 값 또는 실험적으로 근거 있는 값으로 되돌려 검증하는 것이 우선이다.
