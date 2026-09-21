# 질량 비의존 OCP 재검토와 Ai2020 nominal 음극 비교

## 결론

**OCP/동적 전압의 1순위 후보는 Ai2020 nominal equilibrium + measured hysteresis**이다. qOCV MAE는 기존 control보다 2.49 mV, held-out 동적 MAE는 9.03 mV 개선된다. 다만 평균 절대 용량오차는 25.25→37.77 mAh, 최대 절대 용량오차는 56.81→94.59 mAh로 증가하므로 아직 최종 확정 모델은 아니다. 현재 보수적 기준 모델은 **Previous measured-anode ±5% control**로 유지하고, hybrid를 다음 보정의 1순위 후보로 둔다.

| 후보 | qOCV MAE (mV) | qOCV RMSE (mV) | 최대 endpoint 오차 (mV) | 동적 평균 MAE (mV) | 평균 절대 용량오차 (mAh) |
|---|---:|---:|---:|---:|---:|
| 기존 측정 음극 ±5% | 7.13 | 9.36 | 3.80 | 36.12 | 25.25 |
| 측정 음극, 질량 비의존 자유 fitting | 6.15 | 8.29 | 2.98 | 37.80 | 27.28 |
| Ai2020 nominal 음극, 자유 fitting | 4.64 | 5.59 | 2.01 | 48.71 | 24.97 |
| Ai2020 nominal equilibrium + 측정 hysteresis | 4.64 | 5.59 | 2.01 | 27.09 | 37.77 |
| Ai2020 nominal published window | 6.00 | 7.42 | 151.63 | 53.37 | 39.01 |

## 질량을 사용하지 않은 방식

펀칭 질량을 목적함수에서 완전히 제외했다. Full-cell 측정 용량을 hard constraint로 두고 x0, x100, y100, y0를 qOCV, dV/dQ 및 endpoint 전압에 맞췄다. Qn과 Qp는 각각 Qcell/(x100-x0), Qcell/(y0-y100)으로 계산했다. 0.5C, 1C, 2C 데이터는 fitting에 사용하지 않고 검증에만 사용했다.

## 음극 데이터 재현성과 nominal OCP

- 기록된 v1/v2 charge/discharge 네 branch의 평균 용량은 2.4747 mAh이다.
- branch 전체 CV는 0.80%, 범위는 평균 대비 1.77%이다.
- x=0.02~0.80에서 회수 음극 equilibrium OCP는 Ai2020 nominal보다 평균 32.0 mV 낮다.

용량 적분값 자체의 v1/v2 재현성은 현재 파일 안에서는 양호하다. 더 큰 불확실성은 젖은 전극의 질량, 초기 lithiation 상태, branch-to-stoichiometry 변환과 OCP 형상에 있다.

## Formation, SEI와 잔류 전해질 해석

- 잔류 전해질은 펀칭 질량을 증가시키므로 질량 기반 c_s,max와 specific capacity를 왜곡할 수 있다.
- Formation 후 SEI는 irreversible lithium inventory와 첫 cycle 효율, 계면 저항에 영향을 줄 수 있다.
- 그러나 충분히 낮은 전류와 이완 조건에서 graphite staging 전위 자체가 SEI 때문에 수십 mV 이동한다고 단정할 근거는 부족하다. 따라서 현재 약 32 mV의 systematic OCP 차이를 SEI 하나로 설명하지 않는다.
- 우선 확인할 실험 요인은 DMC washing 및 건조 조건, 해체 SOC, 펀칭 위치/면, coin-cell 압력과 wetting time, Li counter-electrode 상태, 첫 cycle 제외 여부다.

## 문헌 위치

- Ai et al., 2020, DOI: https://doi.org/10.1149/2.0122001JES
- Lu et al., 2021, half-cell OCP data processing, DOI: https://doi.org/10.1149/1945-7111/ac11a4
- PyBaMM Ai2020 parameter-set documentation: https://docs.pybamm.org/en/v25.6.0/source/examples/notebooks/models/Validating_mechanical_models_Enertech_DFN.html
