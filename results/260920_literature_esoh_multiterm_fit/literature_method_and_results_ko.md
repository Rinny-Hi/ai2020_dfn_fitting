# Stoichiometry window / OCP 다목적 최적화 검토

## 결론

제시한 형태의 다목적 함수는 선행연구의 방향과 부합한다. 다만 문헌에서는 네 항을 모두 임의 가중치로 더하기보다, OCV 오차와 dV/dQ 오차를 목적함수로 두고 전지 용량 및 전압 끝점은 물리식 또는 경계조건으로 강제하는 경우가 더 명확하다. 따라서 이번 계산에서는 다음과 같이 구성했다.

\[
J=w_1J_{\mathrm{OCV}}+w_2J_{dV/dQ}+w_3J_{\mathrm{endpoint}},
\qquad
J_{\mathrm{capacity}}=0\;\text{(hard constraint)}
\]

전지 실측용량 \(Q_{cell}=2.3377272\,\mathrm{Ah}\)에 대해

\[
Q_{cell}=Q_n(x_{100}-x_0)=Q_p(y_0-y_{100})
\]

를 항상 만족하도록 했으므로 용량 항을 별도 패널티로 중복 계산하지 않았다.

## 선행연구와의 관계

- 64 Ah 파우치셀 OCV 진단 연구는 전극 OCV의 수평 이동과 운용창을 맞추면서 전압 오차와 dU/dQ 오차를 함께 사용하는 비용함수를 제시한다. 즉 \(J_{OCV}+J_{dV/dQ}\) 구성은 문헌에 직접 근거가 있다. [논문 PDF](https://documentserver.uhasselt.be/bitstream/1942/42480/1/Experimental%20Investigation%20of%20a%2064%20Ah%20Lithium-Ion%20Pouch%20Cell.pdf)
- Birkl 등의 parametric OCV 연구는 먼저 전극 유효용량 범위를 추정한 뒤 half-cell OCV와 full-cell OCV를 함께 맞춘다. 용량 범위와 stoichiometry 창을 분리하지 않고 물리적으로 연결하는 접근이다. [논문 PDF](https://www.researchgate.net/profile/Christoph-Birkl/publication/281550775_A_Parametric_Open_Circuit_Voltage_Model_for_Lithium_Ion_Batteries/links/55ed649f08ae3e12184802a6/A-Parametric-Open-Circuit-Voltage-Model-for-Lithium-Ion-Batteries.pdf)
- PyBaMM electrode-SOH 식도 \(Q_n,Q_p,Q_{Li}\), 전압 끝점, 전지 용량을 연립해 \(x_0,x_{100},y_0,y_{100}\)을 결정한다. [PyBaMM electrode-SOH](https://docs.pybamm.org/en/stable/_modules/pybamm/models/full_battery_models/lithium_ion/electrode_soh.html)
- Chen 2020의 DFN 파라미터화는 \(c_{s,max}\)를 전극 loading, 두께, 활물질 분율 등으로 계산하고, OCV 특징을 이용한 stoichiometry mapping과 구분한다. 따라서 OCP만 잘 맞추기 위해 \(c_{s,max}\)를 자유롭게 조정하는 것보다 \(Q_n,Q_p\) 또는 유효 활물질 분율을 먼저 조정하는 편이 해석상 안전하다. [Chen 2020](https://doi.org/10.1149/1945-7111/AB9050)

## 이번 최적화 정의

- 신규 260918 GITT는 사용하지 않았다.
- 기존 half-cell OCP와 full-cell qOCV만 사용했다.
- 최적화 변수: \(Q_n,Q_p,x_0,y_{100}\)
- 종속 변수: \(x_{100}=x_0+Q_{cell}/Q_n\), \(y_0=y_{100}+Q_{cell}/Q_p\)
- 균형 목적함수: 정규화된 OCV, dV/dQ, endpoint 항에 동일 가중치 \((1,1,1)\)
- OCP 후보 선택 후에만 0.5C/1C/2C 충·방전 곡선으로 외부 검증했다.
- 동적 비교에서는 모든 후보에 동일한 음극 hysteresis scale 0.50을 적용해 window/capacity 차이만 비교했다. qOCV 최적화 자체에는 이 scale이 들어가지 않는다.
- \(c_{s,max,n}=29,700\,\mathrm{mol\,m^{-3}}\), \(c_{s,max,p}=49,943\,\mathrm{mol\,m^{-3}}\)는 고정했다.

## 핵심 결과

| 항목 | 기존 고정 \(Q_n,Q_p\) | 선정: 균형목적 ±5% | 진단용: 균형목적 ±15% |
|---|---:|---:|---:|
| \(Q_n\) (Ah) | 3.0273 | 3.1787 (+5.00%) | 3.3360 (+10.20%) |
| \(Q_p\) (Ah) | 4.5992 | 4.3692 (-5.00%) | 4.3784 (-4.80%) |
| \(Q_{Li}\) (Ah) | 4.3891 | 4.2950 | 4.2933 |
| \(x_0 \rightarrow x_{100}\) | 0.0014→0.7736 | 0.0023→0.7378 | 0.0022→0.7029 |
| \(y_{100} \rightarrow y_0\) | 0.4451→0.9534 | 0.4463→0.9813 | 0.4450→0.9789 |
| qOCV MAE, 2–98% (mV) | 16.21 | **7.13** | 6.15 |
| qOCV RMSE, 2–98% (mV) | 18.04 | **9.36** | 8.29 |
| 최대 끝점 오차 (mV) | 0.00 | **3.80** | 2.98 |
| 동적 평균 MAE, SOC 10–70% (mV) | 37.87 | **36.12** | 37.80 |
| 평균 절대 용량오차 (Ah) | 0.0415 | **0.0252** | 0.0273 |
| 최대 절대 용량오차 (Ah) | 0.0833 | **0.0568** | 0.0583 |

선정안은 qOCV MAE를 16.21 mV에서 7.13 mV로 56% 줄였고, 동적 평균 MAE도 37.87 mV에서 36.12 mV로 소폭 개선했다. ±15% 해는 qOCV만 보면 더 좋지만 동적 오차가 다시 커져 OCP 과적합 가능성이 있다. 따라서 현재 자료에서는 ±5% 제한 해를 권장한다.

## 0.5C/1C/2C 선정안 검증

| 조건 | SOC 10–70% MAE (mV) | 종료시간 오차 (min) | 용량오차 (Ah) |
|---|---:|---:|---:|
| 0.5C Charge | 25.11 | -2.990 | -0.0568 |
| 0.5C Discharge | 40.15 | +0.879 | +0.0167 |
| 1C Charge | 26.43 | -0.327 | -0.0124 |
| 1C Discharge | 55.10 | +0.399 | +0.0151 |
| 2C Charge | 27.65 | +0.355 | +0.0270 |
| 2C Discharge | 42.26 | -0.308 | -0.0234 |

1C 방전의 전압 오차가 가장 크므로, 다음 단계에서는 stoichiometry/OCP를 다시 움직이기보다 전달계수·확산계수·저항·열 파라미터를 식별해야 한다.

## 파라미터 반영 권고

현재 해를 DFN에 반영할 때는 \(c_{s,max}\)를 바꾸지 않고 유효 전극용량에 해당하는 활물질 분율을 임시로 다음과 같이 매핑했다.

- 음극 유효 활물질 분율: 0.6100 → 0.6405
- 양극 유효 활물질 분율: 0.6200 → 0.5890

단, 두 값이 물리적인 활물질 분율 변화라는 뜻은 아니다. 면적, 두께, loading, porosity 오차가 모두 동일한 \(Q_n,Q_p\) 변화로 나타날 수 있기 때문이다. 또한 최적해가 ±5% 경계에 닿았으므로, 실제 전극 loading·두께·활물질 함량 자료가 확보되기 전에는 보수적 후보값으로 사용해야 한다.

## 최종 판단

1. 제시한 다목적 접근은 선행연구와 일치한다.
2. 다만 capacity는 soft penalty보다 물리적 hard constraint로 두는 편이 타당하다.
3. 현재 권장안은 균형 목적함수 \((1,1,1)\)와 \(Q_n,Q_p\) ±5% 제한 해이다.
4. \(c_{s,max}\)는 현 단계에서 변경하지 않는다.
5. 다음 식별은 OCP 창을 고정한 뒤 동적 파라미터로 진행한다.
