# Ln/Lp 3종 x Rn/Rp 2종 비교

## 비교 조건

### 전극 두께

| 구분 | Ln | Lp |
|---|---:|---:|
| Ai2020 원논문 | 76.50 um | 68.00 um |
| 두께측정기 | 73.00 um | 62.50 um |
| SEM 단면 | 81.52 um | 62.965 um |

두께 변화 시 coating inventory `L*epsilon_s`와 inactive fraction을 유지하고
porosity로 체적수지를 닫았다.

### 입자 반경

| 구분 | Rn | Rp |
|---|---:|---:|
| Nominal | 5.000 um | 3.000 um |
| SEM ImageJ 평균 | 3.7542 um | 4.2387 um |

SEM 값은 음극 44개, 양극 20개 ROI의 2D 등가원 반경 산술평균이다.

### 고정값

- `Dsn=2.1e-14`, `Dsp=4.4e-14 m2/s`
- `kp=4.18e-7`
- `kn=7.40e-7`: 최종 물성값이 아니라 비교용 placeholder
- 현재 C/50 OCP/window
- 각 시험의 직전 rest 전압 기반 초기 SOC
- 음극 hysteresis 미사용, 양극 hysteresis 사용
- 전극 재측정 질량비는 셀 전체 inventory에 강제하지 않음

## 6개 조합 종합 결과

충전/방전 6개 branch 전체를 합친 값이다.

| Ln/Lp | Rn/Rp | full RMSE | 10-70% RMSE | 용량 RMSE | 평균 절대 용량오차 |
|---|---|---:|---:|---:|---:|
| Ai2020 | Nominal | 37.12 mV | 19.51 mV | 2.10% | 35.97 mAh |
| **Ai2020** | **SEM ImageJ** | **35.87 mV** | **18.31 mV** | 1.97% | 28.90 mAh |
| 두께측정기 | Nominal | 43.10 mV | 19.57 mV | 6.17% | 85.93 mAh |
| 두께측정기 | SEM ImageJ | 42.18 mV | 19.92 mV | 6.06% | 80.54 mAh |
| SEM 단면 | Nominal | 39.86 mV | 25.65 mV | 1.56% | 26.67 mAh |
| **SEM 단면** | **SEM ImageJ** | 39.20 mV | 24.32 mV | **1.26%** | **17.99 mAh** |

## 반경 변경 효과

SEM ImageJ 반경은 세 두께 조합 모두에서 종합 full RMSE와 용량오차를 개선했다.

- Ai2020 두께: 평균 절대 용량오차 `35.97 -> 28.90 mAh`
- 두께측정기: `85.93 -> 80.54 mAh`
- SEM 두께: `26.67 -> 17.99 mAh`

특히 방전에서 개선이 크다.

| Ln/Lp | R | 방전 용량 RMSE | 방전 평균 절대 용량오차 |
|---|---|---:|---:|
| Ai2020 | Nominal | 1.31% | 22.66 mAh |
| Ai2020 | SEM ImageJ | **0.36%** | **7.05 mAh** |
| SEM 단면 | Nominal | 1.09% | 19.49 mAh |
| SEM 단면 | SEM ImageJ | **0.20%** | **3.79 mAh** |

반면 방전 full RMSE는 Ai2020 두께에서 `26.03 -> 26.82 mV`, SEM 두께에서
`31.10 -> 32.99 mV`로 소폭 악화했다. 즉 cutoff 용량과 중앙부 형상은 좋아지지만
일부 초기/말단 전압 잔차가 증가한다.

## 우선 후보별 C-rate 결과

### 균형형: Ai2020 Ln/Lp + SEM ImageJ Rn/Rp

| C-rate | 충전 full RMSE | 충전 용량오차 | 방전 full RMSE | 방전 용량오차 |
|---:|---:|---:|---:|---:|
| 0.5C | 36.79 mV | -58.4 mAh | 28.69 mV | +2.5 mAh |
| 1C | 46.16 mV | -53.3 mAh | 30.95 mV | +6.4 mAh |
| 2C | 51.85 mV | -40.6 mAh | 20.81 mV | -12.3 mAh |

### 용량형: SEM Ln/Lp + SEM ImageJ Rn/Rp

| C-rate | 충전 full RMSE | 충전 용량오차 | 방전 full RMSE | 방전 용량오차 |
|---:|---:|---:|---:|---:|
| 0.5C | 37.88 mV | -44.5 mAh | 30.85 mV | +2.7 mAh |
| 1C | 46.70 mV | -25.3 mAh | 36.48 mV | +7.2 mAh |
| 2C | 51.64 mV | +26.8 mAh | 31.63 mV | -1.5 mAh |

## 판단

비가중 Pareto 기준으로 최종적으로 남는 후보는 두 개다.

1. **전압 형상까지 포함한 균형형 후보:** Ai2020 `Ln/Lp=76.5/68 um` +
   SEM ImageJ `Rn/Rp=3.754/4.239 um`
2. **cutoff 용량 우선 후보:** SEM `Ln/Lp=81.52/62.965 um` +
   SEM ImageJ `Rn/Rp=3.754/4.239 um`

두께측정기 조합은 SEM 반경을 적용해도 2C 용량오차가 크므로 현재 후보에서 제외한다.

SEM 반경은 nominal 반경보다 현재 데이터에서 일관되게 좋아졌으므로 **예비 적용값으로
승격할 근거가 있다.** 다만 `D_s`와 `k`를 고정한 결과이며, ImageJ 2D 등가반경이
실제 DFN 확산반경과 동일하다는 의미는 아니다. 새 충방전 데이터에서 두 Pareto 후보를
독립 validation한 뒤 geometry를 선택하고, 마지막에 `kn`을 fitting한다.

