# 260922 geometry·capacity feasibility 최종 요약

## 최종 채택 상태

현재 baseline을 유지한다. 이번 SEM 기반 `Rn/Rp/Ln/Lp`와 Qn/Qp scan 결과는 모두 feasibility로 공유하며 물성 확정값으로 강제 적용하지 않는다.

| 구분 | 현재 판단 |
|---|---|
| `Dsn`, `Dsp` | 소재 기반값 유지 |
| `kn`, `kp` | 현재 실험값 유지, SOC별 Rct 산출 후 갱신 |
| `c_s,max,n/p` | 현 값 유지, dry loading prior 확보 후 재계산 |
| `Rn`, `Rp` | SEM 임시 6/5 µm는 feasibility만 기록 |
| `Ln`, `Lp` | SEM 81.52/62.965 µm는 geometry 후보로 기록, baseline 강제 변경 안 함 |
| OCP/window | C/50 기반 현 baseline 유지, loading 후 Qn/Qp–QLi eSOH/DVA 재식별 |

## SEM geometry feasibility

| 조건 | 10–70% 전압 RMSE | 용량 RMSE | 평균 절대 용량오차 |
|---|---:|---:|---:|
| 현재 baseline | 15.16 mV | 3.28% | 51.24 mAh |
| nominal L + SEM Rn만 | 17.70 mV | 4.84% | 78.41 mAh |
| nominal L + SEM Rp만 | 26.97 mV | 6.53% | 90.09 mAh |
| nominal L + SEM Rn/Rp | 34.91 mV | 8.20% | 118.52 mAh |
| SEM Ln/Lp raw | 18.00 mV | 6.37% | 94.34 mAh |
| SEM Ln/Lp, Qn/Qp 보존 진단 | 15.77 mV | 3.73% | 56.58 mAh |

해석은 “큰 입자가 나쁘다”가 아니라, 기존 `Ds/k/a_s`와 결합된 baseline에 반지름만 교체하면 유효 동역학이 동시에 바뀐다는 뜻이다. 사용자가 소재 기반으로 확보한 `Ds/k`를 반지름에 맞춰 임의 rescale하지 않는다.

## 용량 계산 검증

- 직접식과 PyBaMM 내부 Qn/Qp: `2.479472/4.371198 Ah`로 일치
- C/50 target Qcell: `2.356528 Ah`
- 0.5C/1C/2C 실제 CC 전류와 2.28 Ah nominal setpoint: 최대 0.005% 이내 일치
- simulation/experiment의 원래 적산 `Q_Ah`를 사용하도록 수정; clipped SOC 역산 제거

현재 baseline endpoint 오차(`model-experiment`)는 다음과 같다.

| C-rate | charge | discharge |
|---:|---:|---:|
| 0.5C | -76.92 mAh | -8.23 mAh |
| 1C | -76.61 mAh | -14.82 mAh |
| 2C | -79.33 mAh | -52.12 mAh |

충전의 약 77–79 mAh 공통 offset은 window/QLi/OCP 쪽 정적 성분을, 방전의 rate-dependent 증가는 동적 polarization 성분을 우선 의심하게 한다.

## Qn/Qp 및 kn 진단 결과

`Qn,Qp +4%`는 현재 grid에서 용량오차를 가장 줄였지만 모든 제약을 통과하지 못했다.

| 지표 | baseline | Qn,Qp +4% feasibility |
|---|---:|---:|
| qOCV RMSE | 6.302 mV | 6.886 mV |
| charge 용량 RMSE | 4.412% | 3.140% |
| discharge 용량 RMSE | 1.426% | 0.933% |
| charge 평균 절대 용량오차 | 77.62 mAh | 55.27 mAh |

`kn`을 키우면 endpoint 용량은 줄일 수 있었지만 charge 전압 RMSE가 20.49→32.31 mV 이상으로 악화됐다. 따라서 `kn`은 용량 보정 노브로 사용하지 않는다.

## 다음 갱신 순서

1. dry coating loading·활물질 분율로 Qn/Qp prior와 불확도를 계산한다.
2. C/50 charge/discharge를 capacity domain에서 fitting해 Qn/Qp prior 안의 QLi/window를 식별한다.
3. 사용자가 산출한 SOC별 Rct 기반 `kn/kp`를 적용한다.
4. 0.5C/1C/2C를 charge fitting–discharge validation으로 다시 평가한다.
5. 전압 RMSE, endpoint 용량오차, qOCV/DVA feature를 각각 보고하고 임의 가중합 하나로 합치지 않는다.

세부 보고서:

- `../260922_sem_geometry_capacity_ablation/README_KO.md`
- `../260922_capacity_accuracy_audit/README_KO.md`
- `../260921_capacity_consistent_qn_qp_kn_refinement/README_KO.md`
