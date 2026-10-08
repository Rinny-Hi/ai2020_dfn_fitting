# 2026-10-08 집 Codex 인수인계

이 문서는 집에서 현재 작업을 그대로 이어가기 위한 최신 기준이다. 이번 작업 범위는 **September 데이터 제외**, **stoichiometry 및 물리 baseline 확정**, **parameter fitting 직전까지**다. 과거 post-fit 문서와 목적이 다르므로 이번 작업에서는 이 문서와 `261008_Enertech_Ai2020_피팅전_P0_2-98.ipynb`를 우선한다.

## 바로 확인할 파일

- 실행·검토본: `261008_Enertech_Ai2020_피팅전_P0_2-98.ipynb`
- stoichiometry 계산 코드: `c50_ocp_source_combination_study.py`
- 선택 결과: `results/260921_c50_ocp_source_combinations/decision.json`
- 공통 면적 근거 코드: `paper_teardown_geometry.py`
- Python 환경: `requirements-gitt-ocp.txt`

`261007_Enertech_Ai2020_피팅전_P0_fullSOC.ipynb`와 방식비교 노트북은 현재 최종 실행본이 아니다. 최종 검토 기준은 SOC 2–98% fitting이다.

## 확정한 stoichiometry 기준

- C/50 full-cell 데이터: 260610 cell-1, cell-2의 complete cycle 2–5
- 평균 대상: charge 8개 곡선과 discharge 8개 곡선
- `Qcell = (charge 평균 + discharge 평균) / 2 = 2.3565276563125 Ah`
- OCP 조합: Ai2020 nominal graphite anode + 기존 측정 cathode
- 목적함수 구간: full-cell SOC 2–98%
- 별도 endpoint 제약: SOC 0%, 100% 전압 오차 절댓값 각각 10 mV 이하
- 전극 용량 bounds: `2.4–6.0 Ah`

| 값 | 결과 |
|---|---:|
| `x0` | 0.0038062237 |
| `x100` | 0.9542212582 |
| `y100` | 0.4316424031 |
| `y0` | 0.9707457852 |
| `delta_x` | 0.9504150345 |
| `delta_y` | 0.5391033821 |
| `Qn` | 2.4794721998 Ah |
| `Qp` | 4.3711980569 Ah |
| qOCV RMSE, SOC 2–98% | 6.242657 mV |
| SOC 0% endpoint 오차 | +10.000 mV |
| SOC 100% endpoint 오차 | -0.875 mV |

2–98%는 full-cell SOC endpoint의 급격한 오차·미분 민감도를 목적함수에서 제외하기 위한 robustness 선택이다. 전극 stoichiometry 범위를 2–98%로 잘라 정의한 것이 아니다. `x0/x100/y0/y100`은 SOC 0/100%까지 선형 사상하며 endpoint 전압 제약으로 양 끝단도 통제한다.

## Geometry 및 effective concentration 기준

1D DFN은 양·음극에 공통인 through-plane 반응면적 하나를 사용한다. 양극 footprint가 더 작으므로 실제 겹치는 면적을 공통 면적으로 둔다.

| 항목 | 값 | 근거/의미 |
|---|---:|---|
| 음극 물리 총면적 | 820.352 cm2 | 46.4 mm × 52.0 mm × 34면; 별도 inventory 참고 |
| 공통 DFN overlap 면적 | 787.236 cm2 | 45.4 mm × 51.0 mm × 34면 |
| `Ln / Lp` | 76.5 / 68.0 um | Ai2020 |
| `Rn / Rp` | 3.7542 / 4.2387 um | SEM ImageJ |
| `eps_s,n / eps_s,p` | 0.61 / 0.62 | teardown 문헌 기반 추정값; 독립 실측 확정값 아님 |
| effective `cs,max,n` | 25,182.8125 mol/m3 | `Qn` 보존 역산값 |
| effective `cs,max,p` | 49,140.1119 mol/m3 | `Qp` 보존 역산값 |

계산식은 다음과 같다.

`Qk = F A_overlap Lk eps_s,k cs,max,k / 3600`

따라서 현재 데이터와 가정으로 직접 식별되는 것은 `eps_s * cs,max`의 곱이다. `eps_s`를 문헌 기반 가정으로 고정한 뒤 `Qn/Qp`를 보존하도록 `cs,max`를 역산했다. 이 값을 teardown으로 독립 측정한 intrinsic 최대농도라고 표현하면 안 된다.

## fitting 전 고정 입력과 재실행 결과

노트북에는 재스케일링 식이 아니라 계산 완료값을 직접 넣었다.

| 항목 | 값 | 상태 |
|---|---:|---|
| `Dsn` | 4.17077414597e-14 m2/s | all-range apparent GITT |
| `Dsp` | 6.88817532381e-14 m2/s | all-range apparent GITT |
| `kn` | 9.648533212e-7 | Ai2020 anchor, 향후 fitting 후보 |
| `kp` | 5.89951990671e-7 | cathode P4 EIS apparent |
| `brugg_n/p/s` | 2.914 / 1.83 / 1.5 | 고정 baseline |

공통 overlap 면적과 effective `cs,max`, 위 GITT/EIS 값을 적용해 July BoL 3셀의 0.5C/1C/2C charge/discharge를 다시 실행한 무피팅 결과:

| 지표 | 결과 |
|---|---:|
| Full MAE / RMSE | 34.68 / 42.56 mV |
| SOC 10–70% MAE / RMSE | 31.53 / 33.45 mV |
| Capacity RMSE | 1.00% |
| Mean / max absolute capacity error | 16.37 / 31.92 mAh |

이 결과가 현재 parameter fitting 전 기준이다. parameter fitting은 아직 실행하지 않는다.

## 집에서 Codex에 붙여 넣을 문장

> GitHub `https://github.com/Rinny-Hi/ai2020_dfn_fitting`의 `codex/home-handoff-261008` 브랜치를 열고 `AGENTS.md`와 `261008_HOME_CODEX_HANDOFF_KO.md`를 먼저 모두 읽어줘. 현재 기준은 September 제외, C/50 SOC 2–98% stoichiometry fitting + SOC 0/100% endpoint ±10 mV 제약, 공통 overlap 면적 787.236 cm2, effective cs,max,n/p = 25182.8125/49140.1119 mol/m3이고 parameter fitting은 아직 하지 않는다. `261008_Enertech_Ai2020_피팅전_P0_2-98.ipynb`의 저장 출력과 입력 파일 경로부터 검증해서 현재 상태를 짧게 보고하고, 다음 작업은 내 지시를 기다려줘. 원본 파일이 없으면 임의 대체하지 말고 정확한 파일명만 알려줘.

## 입력 파일 처리

통합 Excel은 저장소에 포함돼 있다. July BoL 원본 3개는 저장소에 없으므로 같은 OneDrive가 동기화돼 있으면 Codex가 검색해 환경변수로 연결하게 한다.

- `260703 BoL 데이터 6-5.xlsx`
- `260703 BoL 데이터 6-6.xlsx`
- `260703 BoL 데이터 6-8.xlsx`

찾지 못한 경우에만 해당 3개 파일을 집 PC로 옮기면 된다. 회사에서 지금 별도로 할 작업은 없다.
