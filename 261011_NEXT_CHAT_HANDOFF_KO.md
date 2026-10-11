# 2026-10-11 다음 대화 인수인계: OCP·window·hysteresis 비교

## 먼저 읽을 파일과 현재 상태

1. [`AGENTS.md`](AGENTS.md): 입력·결과 검증 및 기존 결과 보존 원칙.
2. [10월 11일 July BoL 단계별 결과](results/261011_july_bol_stage_reproduction/README_KO.md)와 [재현 노트북](261011%20Ai2020%20July%20BoL%20DFN%20단계별%20재현.ipynb).
3. [사용자가 제공한 OCP·window·hysteresis 비교 실행계획](docs/261011_OCP_WINDOW_HYSTERESIS_PLAN.md). **계획이며 아직 새 후보 비교 결과가 아니다.**
4. [9월 27일 기존 최종 physical-time 피팅](results/260927_charge_physical_time_final/README_KO.md). 이 피팅과 10월 11일 단계별 진단은 설정과 지표가 다르다.

현재 기준 OCP는 Ai2020 graphite + 기존 측정 양극 GITT이고, window는 `x0=0.003806`, `x100=0.954221`, `y0=0.970746`, `y100=0.431642`이다. 6월 C/50 qOCV SoC 2–98% RMSE는 6.243 mV다. July 3셀 순차 DFN의 S3(Geometry+OCP+window)는 0.5C/1C/2C 충전 RMSE 46.51/54.83/59.41 mV, 방전 42.90/49.33/47.78 mV다. 이는 겹침 CC 구간 점수이며 전구간 또는 SoC 10–80% 점수가 아니다. 음극 절대 전위 정확성은 검증되지 않았다.

8월 실험 음극+실험 양극 OCP의 저장된 C/20 SoC 2–98% RMSE 10.99 mV는 그 데이터에 대한 in-sample 근거다. 원본 8월 노트북은 `xlsx_files` 미정의로 처음부터 실행되지 않고, Stage 2 동적 출력이 없다. 실험 음극 OCP는 다음 비교의 후보로 남긴다. 저장된 과거 수치를 새 검증 결과로 옮겨 쓰지 않는다.

## 다음 대화에서 준비할 파일

### 필수: 저장소에 없는 5개 원본

이 파일들은 같은 로컬 PC의 `01. 적응형 충전 프로토콜/4. Enertech 파우치셀 실험데이터/` 아래에 존재함을 2026-10-11 확인했다. 새 Codex **로컬** 대화가 같은 PC에 접근하면 파일을 다시 업로드할 필요 없이 아래 이름과 하위 폴더를 지정하면 된다. 다른 호스트나 클라우드 대화라면 다섯 파일을 첨부해야 한다.

| 하위 폴더 | 파일 | 용도 |
| --- | --- | --- |
| `260610 파우치셀 데이터` | `260610 저율데이터 1.xlsx` | C/50 셀 1 qOCV·window |
| `260610 파우치셀 데이터` | `260610 저율데이터 2.xlsx` | C/50 셀 2 qOCV·독립 확인 |
| `260703 실험데이터` | `260703 BoL 데이터 6-5.xlsx` | July BoL 순차 0.5C/1C/2C |
| `260703 실험데이터` | `260703 BoL 데이터 6-6.xlsx` | July BoL 순차 0.5C/1C/2C |
| `260703 실험데이터` | `260703 BoL 데이터 6-8.xlsx` | July BoL 순차 0.5C/1C/2C |

저장소 루트의 `260902 Enertech셀 데이터 모음.xlsx`에는 과거 음·양극 OCP/GITT, C/20, 동적 데이터의 통합 시트가 이미 있으므로 첫 실행에 중복 업로드할 필요가 없다. 신규 July 노트북은 이 통합 Excel과 위 5개를 동일 데이터 폴더에 둔 뒤 `/content` 또는 `AI2020_DATA_DIR`에서 찾는다. 다음 구현은 별도 설정 파일이나 입력 경로 인자로 분리하는 것이 좋다.

### 원본 수준의 GITT/legacy 감사에 추가 권장

통합 Excel만으로는 각 GITT 휴지 구간의 시간별 전압 기울기와 원본 추출 품질을 모두 확인할 수 없다. 비교 실행계획의 rest-end QC, 절대 축, legacy 재현까지 하려면 아래도 첨부하거나 로컬 경로로 연결한다.

| 하위 폴더 | 파일 |
| --- | --- |
| `260808 GITT 실험데이터` | `260816 Anode GITT + OCP v1.xlsx` |
| `260808 GITT 실험데이터` | `260816 Anode GITT + OCP v2.xlsx` |
| `260808 GITT 실험데이터` | `260816 Positive GITT 7-7.xlsx` |
| `260808 GITT 실험데이터` | `260816 Positive GITT 7-8.xlsx` |
| `260808 Eneterch셀 초기 저율 데이터` | `260808 Enertech셀 초기 저율 데이터 v2.xlsx` |

과거 노트북의 셀별 구현을 직접 대조하려면 사용자 보유 `260818_0단계~2단계_진행_결과_RctN2반영.ipynb`도 다음 대화에 첨부한다. 해당 노트북은 최종 결과가 아닌 historical reference다.

September cohort까지 독립 검증하려면 `260927 BoL 실험데이터/260927 6-1 0.5C 1C 2C.xlsx`, `260927 6-2 0.5C 1C 2C.xlsx`가 추가로 필요하다. July와 September의 `6-*`는 같은 종류의 **서로 다른 실물 셀**이며 번호는 위치 번호다. 같은 inventory/window를 자동 공유하지 않는다.

원본 Excel을 이 공개 GitHub 저장소에 무조건 커밋할 필요는 없다. 먼저 저장소 통합 시트와 기존 CSV가 필요한 raw 필드·출처를 보존하는지 확인하고, 부족한 원본만 다음 작업 환경에 제공한다. 원본을 공개 추가할지의 결정은 별도로 한다.

## 평가 기준에서 먼저 확인할 차이

사용자의 이전 직접 요청은 **SoC 10–80%의 시간–전압 오차가 20 mV 이하면 충분**하다는 것이었다. 첨부 비교 실행계획은 더 엄격하게 **각 셀·rate·방향의 full physical-time RMSE ≤20 mV**를 목표로 적었다. 다음 작업에서는 두 지표를 모두 계산하고, 어느 것이 최종 합격 기준인지 명시한다. 사용자가 기준을 바꾸지 않았다면 직접 요청한 SoC 10–80%를 우선 판정하고 full branch는 별도 강건성 지표로 보고한다. 조기 cutoff가 나면 미계산 구간을 버린 RMSE로 달성을 주장하지 않는다.

## 다음 대화에 그대로 붙여 넣을 요청문

```text
https://github.com/Rinny-Hi/ai2020_dfn_fitting 의 최신 main을 확인하고
AGENTS.md, 261011_NEXT_CHAT_HANDOFF_KO.md,
docs/261011_OCP_WINDOW_HYSTERESIS_PLAN.md,
results/261011_july_bol_stage_reproduction/README_KO.md를 먼저 읽어줘.

첨부 계획의 OCP·stoichiometry window·hysteresis 공정 비교를 구현하고
가능한 기존 데이터로 실제 실행해줘. 7월 BoL 3셀을 우선 대상으로 하고,
현재 Ai2020 음극+기존 실험 양극 OCP와 8월 실험 음극+실험 양극 OCP를
같은 geometry, 같은 nominal kinetics/transport, 같은 실제 전류·이력에서 비교해줘.
Kinetics/transport fitting 결과를 이번 OCP 비교 baseline에 섞지 마.

내 기존 목표는 SoC 10–80% 시간–전압 오차 20 mV 이하이고,
첨부 계획에는 full physical-time RMSE 20 mV 이하가 적혀 있어.
두 지표를 모두 계산하고 평가 구간·coverage·조기 cutoff를 명확히 구분해줘.
음극 전위는 local 최소값과 위치·시간, OCP coverage, 초기 상태·mesh 민감도를
보고하되 기준전극 실측 검증을 했다고 표현하지 마.

원본 파일은 이 대화에 첨부하거나 같은 PC의
01. 적응형 충전 프로토콜/4. Enertech 파우치셀 실험데이터 하위 폴더에서 찾아줘.
필수: 260610 저율데이터 1.xlsx, 2.xlsx 및
260703 BoL 데이터 6-5.xlsx, 6-6.xlsx, 6-8.xlsx.
저장소 통합 Excel에 없는 원본 정보가 필요하면 이 문서의 추가 권장 목록을 확인해줘.
기존 결과를 덮어쓰지 말고 새 결과 폴더에 입력 manifest, 비교 설정,
조건별 전압·용량·음극 전위 지표, 선정 근거를 남겨줘.
결과가 20 mV 목표를 못 맞추면 가장 타당한 후보와 남은 오차를 그대로 보고해줘.
```

## 완료 전 확인

- 원본이 있는지, 시트와 방향명이 기대한 것과 맞는지 확인한다.
- 실험 음극 GITT 전체와 저율 full-range OCP를 혼동하지 않는다.
- OCP source를 바꿔도 공통 geometry, 전류, nominal kinetics/transport fingerprint가 같은지 확인한다.
- source별 window를 같은 용량·coverage 규칙에서 식별하고, hysteresis ON/OFF에서는 window를 고정한다.
- first-rest 전압 역산만으로 내부 농도장·h가 알려진다고 주장하지 않는다.
- 후보 선택에 쓴 cell/branch와 검증 cell/branch를 구분한다.
- `TARGET_MET`, `BEST_VALID_PREFIT`, `INFEASIBLE_COVERAGE` 등 계획의 판정을 실제 결과에만 붙인다.
