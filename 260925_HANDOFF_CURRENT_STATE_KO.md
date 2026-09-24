# 260925 Enertech DFN fitting 인수인계

이 문서는 다른 PC와 새 Codex 대화에서 현재 작업을 그대로 이어가기 위한 기준 문서다. 과거 문서의 값과 충돌할 경우 이 문서를 우선한다.

## 1. 저장소와 브랜치

- GitHub: <https://github.com/Rinny-Hi/ai2020_dfn_fitting>
- 작업 브랜치: `260921-ocp-eis-current-progress`
- 2026-09-24까지의 주요 결과 커밋:
  - `00eed3a`: 절대 Ah OCP/eSOH 및 `kn` trade-off
  - `7d3700d`: cathode EIS 기반 `kp` 재계산 비교

새 PC에서 시작할 때:

```powershell
git clone https://github.com/Rinny-Hi/ai2020_dfn_fitting.git
cd ai2020_dfn_fitting
git checkout 260921-ocp-eis-current-progress
git pull origin 260921-ocp-eis-current-progress
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements-gitt-ocp.txt
```

일부 분석 스크립트는 `C:\Users\user\OneDrive\...`의 원본 Excel 파일을 절대경로로 참조한다. 새 PC의 OneDrive 경로가 다르면 소스 경로를 먼저 수정해야 한다.

## 2. 새 대화에 보낼 시작 문장

아래 문장을 GitHub 링크와 함께 새 대화에 보내면 된다.

> `260921-ocp-eis-current-progress` 브랜치의 `260925_HANDOFF_CURRENT_STATE_KO.md`를 먼저 읽고 이어서 진행해줘. 최신 결정은 음극 half-cell Rct/kn 값을 폐기하고 kn은 나중에 fitting하는 것이다. 현재 geometry 적용값과 새로 측정한 Ln, Lp, Rn, Rp를 비교한 뒤, 확정 geometry를 반영하고 cathode kp 및 kn fitting/validation을 순서대로 진행해줘. 과거 md와 충돌하면 260925 handoff 문서를 우선해줘.

## 3. 현재 geometry 정리

### 3.1 코드에 현재 적용된 값

| 파라미터 | 현재 코드값 | 상태 | 다음 조치 |
|---|---:|---|---|
| `Ln` | 76.5 um | Ai2020/SPB655060 기준 임시 적용 | 새 단면 측정값과 비교 후 결정 |
| `Lp` | 68.0 um | Ai2020/SPB655060 기준 임시 적용 | 새 단면 측정값과 비교 후 결정 |
| `Rn` | 5.0 um | nominal/feasibility 값 | 면적등가반경 분포 확보 후 결정 |
| `Rp` | 3.0 um | nominal/feasibility 값 | 면적등가반경 분포 확보 후 결정 |
| Cu 집전체 | 사용자 측정 9 um | geometry 기록값 | 전기화학 coating 두께와 분리 관리 |
| Al 집전체 | 사용자 측정 13 um | geometry 기록값 | 전기화학 coating 두께와 분리 관리 |
| `A_n` | 820.352 cm2 | 논문 stack geometry | dry loading 계산 시 사용 |
| `A_p` | 787.236 cm2 | 논문 stack geometry | DFN overlap 기준 후보 |

코드 상수는 `nominal_radius_current_configuration.py`에 있다.

```text
RN_UM = 5.0
RP_UM = 3.0
LN_UM = 76.5
LP_UM = 68.0
```

### 3.2 SEM에서 얻은 초기 참고값

표면 SEM의 표시 길이를 단순 직경으로 간주한 초기 산술평균은 다음과 같다.

- 음극 표시 길이 평균: 약 12.36 um → 단순 반경 약 6.18 um
- 양극 표시 길이 평균: 약 9.73 um → 단순 반경 약 4.87 um

이 값은 입자 면적 segmentation 결과가 아니며 응집체와 1차 입자를 구분하지 못한다. 따라서 현재 `Rn/Rp`를 교체하지 않았다. 최종값은 입자별 면적등가반경

```text
R_eq = sqrt(A_particle / pi)
```

의 D10/D50/D90과 분석 대상이 1차 입자인지 2차 입자인지를 함께 기록해 결정한다.

### 3.3 두께 적용 원칙

양면 전극 단면에서는 다음 식으로 단면 coating 두께를 계산한다.

```text
L_one_side = (total electrode thickness - current collector thickness) / 2
```

상·하부 coating이 비대칭이면 각각 기록하고 평균과 표준편차를 함께 보관한다. 새 측정값을 바로 확정값으로 넣지 말고 `76.5/68.0 um` baseline과 A/B 검증한다.

## 4. 현재 OCP와 stoichiometry 기준

신규 260918 GITT는 재현성 문제가 정리될 때까지 최종 OCP 선정에서 제외했다.

- 음극 OCP: Ai2020 nominal graphite equilibrium
- 음극 hysteresis: 사용하지 않음
- 양극 OCP: 기존 old cathode GITT
- 양극 hysteresis: charge delithiation / discharge lithiation 사용
- C/50 complete cycle 2–5의 charge/discharge 평균 사용

| 항목 | 현재값 |
|---|---:|
| `Qcell` | 2.356528 Ah |
| `Qn` | 2.479472 Ah |
| `Qp` | 4.371198 Ah |
| `x0` | 0.003806 |
| `x100` | 0.954221 |
| `y100` | 0.431642 |
| `y0` | 0.970746 |
| qOCV RMSE | 약 6.24–6.31 mV |

개별 C/50 곡선을 SOC로 정규화한 방식과 절대 Ah 방식의 qOCV 차이는 RMSE 0.021 mV에 불과했다. 약 77 mAh의 충전 용량오차는 C/50 정규화 때문에 생긴 것이 아니다. `Qn/Qp`를 ±5% 범위에서 재적합해도 qOCV 개선은 약 0.021 mV뿐이었고 충전 용량오차가 악화했으므로 기존 window를 유지한다.

주요 결과:

- `results/260923_absolute_capacity_esoh_recheck/README_KO.md`
- `absolute_capacity_esoh_recheck.py`

## 5. transport와 effective capacity 파라미터

| 파라미터 | 현재값 | 상태 |
|---|---:|---|
| `Dsn` | 2.1e-14 m2/s | 소재 기반 지정값 유지 |
| `Dsp` | 4.4e-14 m2/s | 소재 기반 지정값 유지 |
| `b_n` | 2.914 | Ai2020 nominal 유지 |
| `b_p` | 1.83 | Ai2020 nominal 유지 |
| `b_s` | 1.5 | Ai2020 nominal 유지 |
| effective `csn,max` | 약 24,325.5 mol/m3 | fitted Qn과 geometry를 표현하는 effective 값 |
| effective `csp,max` | 약 47,467.3 mol/m3 | fitted Qp와 geometry를 표현하는 effective 값 |

`c_s,max`는 현재 소재 고유값이 아니라 loading, 활물질분율, 두께 불확실성을 일부 흡수한 effective parameter다. 독립 loading 없이 `L`, `epsilon_s`, `c_s,max`를 동시에 fitting하지 않는다.

## 6. kinetic parameter의 최신 결정

### 6.1 음극: Rct와 kn 폐기

음극 half-cell은 재현성이 매우 낮고, 세척 시간을 늘려도 전극 표면에 자국이 남는다. 잔류 전해액, SEI, 회수/세척 이력에 의한 오염 가능성이 있으므로 음극 Rct에서 계산한 `kn`은 물성 anchor로 사용하지 않는다.

- 과거 음극 Rct 기반 `kn=2.48e-7`, `7.40e-7` 해석: **폐기**
- 코드의 `KN_PREF=7.40e-7`: 과거 결과 재현을 위한 **임시 baseline만 유지**
- 보고서나 발표에서 `kn=7.40e-7`을 실험 확정값으로 표현하지 않는다.
- 최종 `kn`: geometry, OCP, `kp`, loading 관련 값이 정리된 뒤 fitting

권장 `kn` fitting/validation 절차:

1. 0.5C/1C/2C charge를 fitting set으로 사용
2. discharge는 독립 validation으로 유지
3. full-range RMSE, 10–70% RMSE, C-rate별 용량오차를 별도로 보고
4. 임의 가중합 하나로 선택하지 말고 Pareto 후보를 비교
5. bound 도달 여부와 C-rate별 일관성을 확인

과거 `kn` 스캔은 trade-off 확인용으로만 남긴다. Rct가 폐기되었으므로 그 결과를 물리적 최종값 선정으로 해석하지 않는다.

### 6.2 양극: cathode Rct 기반 kp

cathode half-cell SOC 50%의 면적정규화 Rct:

```text
Rct,p = 18.78433 ohm cm2
```

Chen2020 방식:

```text
b_p = 3 epsilon_s,p / R_p
Rct,ASR = R T / (j0 F b_p L_p)
j0 = kp sqrt(c_e (c_s,max-c_s) c_s)
SOC 50%: c_s = c_s,max / 2
```

현재 조건에서:

| 표현 | epsilon_s,p | c_s,max | kp |
|---|---:|---:|---:|
| Chen/물성 표현 | 0.620 | 49,943 mol/m3 | 4.11e-7 |
| fitted Qp에서 역산한 표현 | 약 0.610 | 49,943 mol/m3 | 4.18e-7 |
| 현재 DFN effective 표현 | 0.620 | 47,467 mol/m3 | 4.32e-7 |

세 값은 서로 다른 `epsilon_s-c_s,max` 표현에서 같은 Rct를 나타내는 조건부 값이다.

- 현재 코드 기본 `KP_PREF`: 3.12e-7, legacy baseline
- 다음 우선 후보: 4.18e-7
- 현재 effective DFN에 Rct를 엄밀히 재표현한 sensitivity: 4.32e-7
- loading과 intrinsic `c_s,max`가 확정되기 전까지 4.18/4.32를 최종 상수로 고정하지 않는다.

`kp=4.18e-7` 결과:

| 방향 | 기존 3.12e-7 | 신규 4.18e-7 |
|---|---:|---:|
| 충전 평균 절대 용량오차 | 77.54 mAh | 49.27 mAh |
| 충전 용량 RMSE | 4.41% | 2.66% |
| 충전 full RMSE | 50.07 mV | 48.21 mV |
| 충전 10–70% RMSE | 19.99 mV | 25.88 mV |
| 방전 평균 절대 용량오차 | 24.94 mAh | 22.66 mAh |
| 방전 용량 RMSE | 1.42% | 1.31% |
| 방전 full RMSE | 25.48 mV | 26.03 mV |

용량과 고율 full-range 전압은 개선하지만 저율/중앙부 형상은 악화한다. 자세한 내용:

- `results/260924_kp_eis_recalculation/README_KO.md`
- `kp_eis_recalculation_ablation.py`

## 7. porosity와 active-material fraction

porosity `epsilon_e`와 활물질 부피분율 `epsilon_s`는 같지 않다.

```text
epsilon_s + epsilon_e + epsilon_binder/carbon = 1
```

질량분율이 기밀이라면 전극 질량과 두께만으로 coating bulk density는 구할 수 있지만 `epsilon_s`와 porosity를 각각 유일하게 분리할 수는 없다. 가능한 방법:

1. scraped coating의 He pycnometer skeletal density로 porosity 계산
2. 불가능하면 nominal `epsilon_s`를 유지하고 범위 sensitivity 수행
3. fitted `Qp`, 면적, 두께, intrinsic `c_s,max`에서 effective `epsilon_s,p` 역산

현재 `Qp`, `A_p`, `Lp=68 um`, `c_s,max=49,943`으로 역산한 `epsilon_s,p`는 약 0.610이며 reference 0.620과 가깝다.

## 8. 사용자가 다음에 준비할 값

가장 먼저 아래 표만 채우면 된다.

| 항목 | 새 측정값 | 단위 | 측정/해석 메모 |
|---|---:|---|---|
| `Ln` |  | um | 단면 coating 한 면 두께 |
| `Lp` |  | um | 단면 coating 한 면 두께 |
| `Rn D50` |  | um | 면적등가반경, 1차/2차 입자 구분 |
| `Rn D10/D90` |  | um | 분포 범위 |
| `Rp D50` |  | um | 면적등가반경, 1차/2차 입자 구분 |
| `Rp D10/D90` |  | um | 분포 범위 |
| 양극 dry coating loading |  | g/m2/side | 집전체 질량 차감 |
| 음극 dry coating loading |  | g/m2/side | 세척 자국으로 신뢰도 낮을 수 있음 |
| cathode Rct 온도 |  | degC | 18.78433 ohm cm2의 측정 온도 |
| cathode Rct 정규화 확인 |  |  | raw ohm인지 ohm cm2인지 확인 |

음극은 washing 후 자국이 계속 남으므로, formed/harvested 음극 질량과 Rct를 독립 물성값으로 강제하지 않는다. 가능하면 미사용 dry electrode를 우선한다.

## 9. 다음 계산 순서

1. 새 `Ln/Lp/Rn/Rp`를 현재 baseline과 A/B 비교
2. dry loading이 있으면 `epsilon_s`, `Qn/Qp`, effective `c_s,max`의 일관성 확인
3. cathode `kp=4.18e-7`과 current-model 표현 `4.32e-7` 비교 후 provisional 선택
4. 위 값을 고정한 뒤 `kn` charge fitting
5. discharge 및 0.5C/1C/2C 전체 validation
6. 최종 파라미터를 `measured / literature / fitted / effective`로 구분해 표 작성

## 10. 실행 파일과 결과 위치

최신 핵심 실행 파일:

- `c50_selected_ocp_dynamic_validation.py`: 현재 OCP 동적 검증
- `absolute_capacity_esoh_recheck.py`: 절대 Ah OCP/eSOH 검증
- `kp_eis_recalculation_ablation.py`: 3.12/4.18/4.32e-7 `kp` 비교
- `kn_pareto_scan_current_ocp.py`: `kn` trade-off 진단, 최종 물성 선정용 아님
- `kinetics_original_reversion_ablation.py`: Ai2020 원 kinetic 비교
- `paper_geometry_area_porosity_ablation.py`: 논문 geometry/porosity 영향
- `sem_geometry_capacity_ablation.py`: SEM geometry와 capacity 보존 비교

최신 결과:

- `results/260923_absolute_capacity_esoh_recheck/`
- `results/260923_kn_pareto_current_ocp/`
- `results/260923_kinetics_original_reversion/`
- `results/260924_kp_eis_recalculation/`

## 11. 재현성 주의사항

- 2026-09-24 `kp` 비교는 로컬 `.venv`의 SciPy DLL 오류 때문에 Anaconda Python의 PyBaMM 25.12.2로 완료했다. `requirements-gitt-ocp.txt`는 PyBaMM 26.8.0.0을 명시한다. 최종 채택 전 새 PC의 고정 환경에서 한 번 재실행한다.
- PyBaMM 25.10 이후 hysteresis decay-rate 정의 변경 경고가 발생한다. 같은 버전으로 비교해야 한다.
- 원본 Excel 경로와 데이터 선택 cycle을 바꾸면 결과가 달라지므로 스크립트의 절대경로와 cycle filtering을 확인한다.
- `outputs/`는 현재 Git에 포함하지 않은 로컬 산출물 폴더다.

