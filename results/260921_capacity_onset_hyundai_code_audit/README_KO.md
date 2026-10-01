# 용량오차 발생 원인 및 현차셀 코드 신뢰성 검토

## 결론

현재 충전 용량 부족의 가장 큰 **변경 파라미터 계열**은 `k_p`와 OCP fitting에서 함께 결정된 양극 유효용량(`Q_p`, `csp,max`, `Δy`)이다. `Dsn`, `Dsp`, 반경, 두께, Bruggeman, `De × 10^-4`가 주원인이라는 근거는 없다.

- kinetics 중에는 `k_p=3.12e-7`가 가장 큰 용량 lever다. nominal `k_p`로만 복원하면 충전 평균 절대 용량오차가 92.85 → 32.37 mAh로 감소하지만 중앙 전압 RMSE는 15.50 → 36.08 mV로 악화된다. 따라서 nominal 복원은 채택할 수 없다.
- `k_n` nominal 복원은 용량오차를 92.85 → 59.71 mAh로 줄이면서 전압 RMSE가 15.50 → 16.18 mV에 머문다. 현재 우선 fitting 대상으로 `k_n`을 택한 판단은 유지할 수 있다.
- 양극 `csp,max`만 Ai2020 nominal로 복원하면 충전 용량오차가 23.24 mAh로 감소하고 방전은 거의 유지된다. 그러나 이는 OCP fitting에서 얻은 `Q_p=4.371 Ah`, `Δy=0.5391` 관계를 깨뜨린다. 직접 교체할 값이 아니라 **양극 활물질량/활물질분율/두께/면적/층수와 stoichiometry window를 함께 다시 식별해야 한다는 진단 신호**다.
- `Dsp` nominal 복원은 충전 오차를 141.07 mAh로 악화한다. 현재 `Dsp=4.4e-14`는 용량오차 증가 원인이 아니다.
- branch별 실제 전류 사용은 92.85 → 83.03 mAh, rest/history 초기조건은 92.85 → 77.54 mAh로 개선하지만 단독 해결책은 아니다.
- 양극 hysteresis를 끄면 66.66 mAh까지 개선되지만 충전 중앙 RMSE가 18.43 mV로 증가한다. hysteresis 선택도 목적함수에 포함해야 한다.

`De ×10^-4`는 오류가 아니다. PyBaMM 26.8의 Ai2020 함수는 1000 mol/m3, 298.15 K에서 약 `3.22e-6`을 반환한다. `×10^-4` 적용 후 약 `3.22e-10 m2/s`가 되므로 원 상관식의 `cm2/s → m2/s` 단위 환산이다. `×1`은 물리적으로 부적절하다.

## 현재 오차가 단일 mass/capacity 오차가 아닌 이유

현재 모델은 충전에서는 세 조건 모두 약 87–103 mAh 짧지만, 방전은 -18–+27 mAh로 비교적 가깝다. 전극 총용량만 작게 산정됐다면 충전과 방전 종점이 같은 방향으로 움직여야 한다. 실제로 양 전극 용량을 4% 늘리면 충전은 -12.42 mAh로 좋아지지만 방전은 +78.66 mAh로 악화된다.

따라서 남은 오차는 다음의 조합으로 해석하는 것이 타당하다.

1. 양극 반응속도와 양극 유효용량/window의 trade-off
2. C/50 OCP 셀(2.3565 Ah)과 동특성 시험 셀/사이클의 용량 차이
3. 독립 branch endpoint 초기화와 실제 rest/history 차이
4. 양극 hysteresis 모델 및 PyBaMM 25.10 이후 decay-rate 정의 변경 영향
5. nominal C-rate 전류와 실제 전류의 소규모 차이

## 현차셀 노트북 신뢰성

검토 대상은 `현차셀_Stoichiometry (8).ipynb`의 코드와 저장된 출력이다. 원본 Hyundai 데이터 파일은 노트북의 Colab 경로에만 있고 이 PC에서 발견되지 않아, 완전 재실행 검증은 불가능했다.

### 신뢰할 수 있는 부분

- 0.05C fitting 후 0.5C/1C/2C를 별도 평가하므로 완전한 same-data fitting은 아니다.
- 각 branch의 실제 전류를 사용하고, 직전 rest 종점전압으로 초기 SOC를 정한다.
- CC-only와 CC-CV를 분리 보고한다.
- CC-only 결과는 평균 전압 MAE 14.63 mV, 용량오차 2.63%로 내부적으로 일관된다.
- CC-CV 결과도 CC 용량오차 1.53%, CC-CV 총량오차 0.48%, 전압 MAE 18.26 mV로 저장된 출력끼리 일관된다.

### 최종 근거로 사용하기 전에 고쳐야 할 부분

1. 첫 설치 셀이 `pip install pybamm`으로 버전을 고정하지 않는다. 아래 주석에는 26.8.0.0이 있으나 실제 실행 보장은 없다.
2. 모든 셀의 `execution_count`가 비어 있고 유사한 통합 코드가 여러 셀에 중복되어 있다. 저장 출력이 현재 소스와 같은 실행에서 나온 것인지 notebook 자체로 증명되지 않는다.
3. 파일 선택이 `glob(...)[0]`이다. 폴더에 파일이 추가되면 다른 CSV/OCP가 조용히 선택될 수 있다.
4. 주석은 `cs,max` 분자량 오류 수정값을 33106/47489라고 쓰지만 실제 상수는 계속 28032/45488이다. 코드와 설명이 충돌한다.
5. `Δθ`를 tear-down `cs,max`로 고정한다고 설명하지만 실제 최적화는 네 endpoint를 각 span 내부에서 독립적으로 움직인 뒤 `cs,max`를 다시 역산한다. 즉 tear-down span은 hard equality가 아니라 bound 역할만 한다.
6. 0.05C 충·방 곡선을 각자의 최대용량으로 정규화한 뒤 평균한다. OCV 형상에는 유용하지만 충·방 용량차를 축에서 제거하므로 capacity consistency를 검증하지 않는다.
7. 0.05C 평균은 pseudo-OCV이지 평형 OCV가 아니다. 평균화가 저항 과전압 일부를 상쇄하지만 확산·비대칭 hysteresis까지 제거한다고 보장할 수 없다.
8. `z0_of`는 같은 pseudo-OCV로 rest 전압을 SOC로 역변환한다. 초기 SOC 개선에는 유용하지만 OCP fitting과 초기화가 부분적으로 순환적이다.
9. CC-only 모델은 전역 cutoff를 실험 2.5/4.4 V가 아니라 2.45/4.45 V로 넓혀 둔다. 보고된 CC-only 종점 용량은 엄밀한 동일 프로토콜 비교가 아니다.
10. CC-CV 모델은 window 밖 OCP를 음극 -0.5 V, 양극 6 V까지 인위적으로 기울인다. CV 폭주 방지용 수치 guard지만, endpoint capacity에 영향을 줄 수 있어 guard 민감도 보고가 필요하다.
11. 전압 MAE는 실험·모델의 겹치는 용량 구간 99.5%에서만 계산한다. 곡선 형상 비교로는 타당하지만 unmatched tail은 용량오차에서만 보인다.
12. 각 branch를 독립적으로 초기화해 이전 C-rate 시험의 농도 구배와 hysteresis history를 제거한다. 완전한 protocol validation은 아니다.

종합 평가는 **탐색용 baseline으로는 신뢰 가능하지만, 최종 parameter identification 근거로는 보완이 필요**하다. 특히 “CC-CV 용량오차 0.48%”만 인용하면 과도하게 낙관적이며, 비교 가능한 대표값은 CC-only 2.63%와 branch별 CC 오차다.

## 현차셀 코드와 현재 코드 비교

| 항목 | 현차셀 코드 | 현재 Enertech 코드 | 판단 |
|---|---|---|---|
| base set | Chen2020 | Ai2020 | 화학계가 달라 단순 우열 없음 |
| OCV target | 0.05C 충·방 평균 pOCV | C/50 qOCV + endpoint/OCP source QC | 현재 코드가 OCP 선택 근거는 더 엄격 |
| half-cell OCP | 양·음극 모두 충·방 평균 | nominal 음극 + old cathode GITT, 양극 hysteresis | 현재 방식이 데이터 품질 차이를 반영하지만 hysteresis 검증 필요 |
| window/capacity | tear-down span을 bound로 쓰고 csmax 재역산 | Qn/Qp와 endpoint/qOCV 동시 식별 | 둘 다 csmax와 window 결합성이 남음 |
| 초기 SOC | 직전 rest 전압 역산 | endpoint와 rest/history 두 시나리오 | 현차 방식의 branch 초기화는 참고 가치 있음 |
| current | branch 실제 A | nominal C-rate | 현재 코드에 실제 전류 반영 필요 |
| protocol | CC-only와 CC-CV 둘 다 | 현재 핵심 비교는 CC cutoff | 목적별로 둘 다 보고해야 함 |
| error metric | overlap MAE + 별도 capacity | full/10–70% RMSE·MAE + capacity | 현재 코드가 더 보수적이고 상세 |
| validation independence | 고율 branch held-out이나 독립 초기화 | 고율 held-out이나 독립 초기화 | 둘 다 sequential validation 추가 필요 |
| reproducibility | glob·unversioned install·중복 셀 | pinned requirements·CSV 결과 저장 | 현재 코드 우세, 단 config 중복은 개선 필요 |

현재 모델의 고율 surface stoichiometry는 음극 0.0038–0.9807, 양극 0.4316–0.9735이고 OCP 표 범위는 음극 0–1, 양극 0.4067–0.9914다. 따라서 현차 코드처럼 극단적 OCP extrapolation guard를 추가할 필요는 현재 결과에서 확인되지 않는다.

## 적용 우선순위

1. 실제 branch 전류와 직전 rest 기반 초기 SOC를 기본 validation에 반영한다.
2. `Qp/csp,max/Δy`를 독립 변경하지 말고 한 묶음으로 재식별한다. qOCV RMSE, endpoint, C/50 용량을 hard constraint 또는 별도 목적함수로 유지하고 `Qp`를 ±6% 범위에서 진단한다.
3. 그 상태에서 `k_n`을 먼저 fitting한다. 현재 결과상 전압 형상을 거의 훼손하지 않고 용량오차를 약 33 mAh 줄이는 방향이다.
4. 남은 충전 종점 오차에 대해서만 `k_p`를 좁은 범위에서 추가 탐색한다. nominal까지 자유롭게 올리면 전압 형상이 크게 무너진다.
5. positive hysteresis on/off를 동일 목적함수로 비교하고, PyBaMM 26.8 decay-rate 정의에 맞춘 단위/식 검증을 추가한다.
6. 최종 평가는 CC-only를 주 지표로 하고 CC-CV는 별도 보조 지표로 제시한다. 마지막에는 전체 시험 sequence를 연속 simulation하여 독립 branch 초기화 이점을 제거한다.

