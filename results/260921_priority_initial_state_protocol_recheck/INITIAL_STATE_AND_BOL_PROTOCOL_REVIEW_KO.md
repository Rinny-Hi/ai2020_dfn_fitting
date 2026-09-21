# 초기상태 및 BoL 프로토콜 재검토

## 결론

현재 전압오차는 단순한 초기 SOC 숫자보다 시험 이력과 hysteresis 초기상태의 영향이 크다. Enertech C-rate 데이터에는 branch 직전 10분 휴지가 있으나, 충전 전 휴지는 평형에 도달하지 않았다. 현재 qOCV 데이터도 선택된 C/20 charge/discharge branch 직전에 휴지가 없다.

따라서 기존 데이터는 원본 RPT 순서를 연속 모사하는 용도로 계속 사용할 수 있지만, 평형 OCP와 stoichiometry window를 높은 신뢰도로 확정하려면 별도의 저율/GITT 데이터를 다시 측정하는 것이 좋다.

## 실제 BoL C-rate 프로토콜

모델 검증에 사용한 데이터는 `260808 Eneterch셀 열화데이터 v1 6-6.xlsx`의 첫 RPT Step 7/10/12/15/17/20이다. 모든 branch 직전에 10분 Rest가 있다.

| 다음 branch | 10분 동안 전압 변화 | 마지막 3분 기울기 | qOCV 역산 초기 SOC |
|---|---:|---:|---:|
| 0.5C Charge | +54.0 mV | +2.88 mV/min | 0.074% |
| 1C Charge | +164.3 mV | +6.01 mV/min | 0.472% |
| 2C Charge | +249.4 mV | +6.56 mV/min | 1.094% |
| 0.5C Discharge | -15.7 mV | -0.36 mV/min | 98.868% |
| 1C Discharge | -16.0 mV | -0.39 mV/min | 98.840% |
| 2C Discharge | -15.5 mV | -0.36 mV/min | 98.905% |

충전 전 10분 Rest는 종료 시점에도 전압이 빠르게 상승한다. 이 전압을 평형 qOCV에 직접 역매핑하면 잔류 분극을 SOC로 잘못 해석한다.

## 초기상태 ablation

모든 값은 모델과 실험의 공통 용량구간 99.5%에서 계산한 C-rate 평균 RMSE이다.

| Scenario | Charge RMSE | Discharge RMSE | 해석 |
|---|---:|---:|---|
| A. endpoint + 새 방향 hysteresis | 52.0 mV | 32.3 mV | 기존 방식 |
| B. 10분 rest 전압으로 SOC만 보정 | 59.0 mV | 38.7 mV | 비평형 전압 역산으로 악화 |
| C. SOC 보정 + 직전 방향 hysteresis | 45.4 mV | 17.3 mV | 시험 이력 반영 효과가 큼 |
| D. SOC 보정 + 평균 OCP | 49.7 mV | 22.7 mV | hysteresis 제거도 방전 개선 |
| E. 첨부 노트북의 R²/Ds 적용 | 72.1 mV | 76.1 mV | 해당 물성을 복사하면 악화 |

Scenario E에서 시작전압은 개선되지만 곡선 전체가 크게 틀어진다. 첨부 노트북의 빠른 확산조건은 해당 현차셀에는 잘 맞지만 Enertech셀의 파라미터로 직접 사용할 근거는 없다. 현재 느린 확산조건은 다른 모델오차를 보상하고 있을 가능성이 있다.

## 원본 RPT 연속 모사

`C/20 방전 → Rest → 0.5C 충전/CV → Rest → 0.5C 방전 → ... → 2C`를 하나의 DFN Experiment로 연속 모사했다.

| 방향 | Full RMSE | 10–70% RMSE | 용량 RMSE |
|---|---:|---:|---:|
| Charge | 57.7 mV | 19.1 mV | 6.17% |
| Discharge | 28.6 mV | 11.2 mV | 3.18% |

연속 모사는 방전 중앙 형상을 크게 개선했지만, 충전 시작전압과 충전용량 부족은 해결하지 못했다. 따라서 시험 이력 누락은 실제 원인 중 하나지만 유일한 원인은 아니다.

## kn 재검토

초기 SOC와 직전 hysteresis 방향을 반영한 뒤 charge를 재평가했다.

| 선택 | kn | Charge 10–70% RMSE | Charge full RMSE | Charge 용량 RMSE | Discharge full RMSE |
|---|---:|---:|---:|---:|---:|
| 전압만 최소화 | 5.0e-7 | 17.4 mV | 50.5 mV | 6.86% | 22.8 mV |
| 용량 RMSE 3% 제약 | 9.0e-7 | 28.7 mV | 45.3 mV | 2.82% | 19.8 mV |

전압만 최소화한 5.0e-7은 용량오차를 키우므로 채택하지 않는다. 현재 정보에서는 9.0e-7이 더 균형적이지만, 초기 휴지구간을 모델이 충분히 재현하지 못하므로 최종 물성값으로 확정하지 않는다.

## qOCV 데이터 취득 방식

현재 Stage 1 qOCV는 `260808 Enertech셀 초기 저율 데이터 v2.xlsx`의 Step 9 C/20 discharge와 Step 10 C/20 charge를 SOC 기준으로 평균해 만들었다. 두 branch 모두 직전 단계 후 Rest가 0분이다.

- SOC 10–90% charge-discharge 간격 평균: 41.4 mV
- 최대 간격: 58.3 mV
- 따라서 현재 곡선은 pseudo-OCV로 사용할 수 있지만 equilibrium OCV로 해석하면 안 된다.

## 데이터 재측정 판단

### 재측정 없이 가능한 작업

- 기존 C-rate 데이터를 원본 순서대로 연속 모사
- 각 Rest 구간까지 fitting/validation에 포함
- 전압 형상과 종료용량을 별도 지표로 평가
- HPPC는 원래의 10분 Rest–pulse–40초 Rest 순서를 유지해 사용

### 재측정을 권하는 작업

- OCP와 stoichiometry window 최종 확정
- Ds 또는 장시간 확산/완화 파라미터 식별
- C-rate 간 공정한 독립 validation

재측정 시에는 모든 C-rate 조건에서 같은 preconditioning을 사용하고, charge와 discharge 시작 전에 충분한 Rest를 넣어야 한다. 현차셀과 직접 비교하려면 우선 동일한 2시간 Rest를 사용하는 것이 가장 명확하다. 가능하면 고정 시간만 쓰지 말고 Rest 말단의 전압 변화율도 함께 저장해 평형 도달 여부를 판단한다. qOCV는 C/50 또는 C/20 양방향 저율곡선과 endpoint Rest를 확보하고, 더 엄밀한 OCP가 필요하면 SOC step별 GITT/relaxation 데이터를 별도로 취득한다.
