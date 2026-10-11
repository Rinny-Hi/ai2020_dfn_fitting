> 사용자 제공 비교 실행계획을 보존한 문서입니다. 이 계획의 새 비교와 20 mV 목표는 아직 실행·검증되지 않았습니다. 현재 저장 결과는 [10월 11일 단계별 결과](../results/261011_july_bol_stage_reproduction/README_KO.md)와 9월 최종 피팅 문서를 구분해 확인하세요.

# DFN OCP·stoichiometry window·hysteresis 비교 실행계획

작성일: 2026-10-11<br>
대상 저장소: https://github.com/Rinny-Hi/ai2020_dfn_fitting<br>
검토한 저장소 HEAD: `aa003e9bd2d058de3672ad514043f79bf1ade8c7`<br>
용도: Codex에서 비교 코드를 구현하고 실행할 때 사용하는 작업 명세

## 1 목표와 이번 작업 범위

추가 실험 없이 기존 측정 데이터와 코드를 활용해, 선택 과정을 설명할 수 있으면서 전압 예측이 가장 좋은 OCP·window·hysteresis 조합을 선정한다

최신 geometry를 모든 주 비교 모델에 동일하게 적용한다<br>
kinetics와 transport는 아직 실험값 적용 또는 fitting을 하지 않은 공통 기준으로 유지한다<br>
이번 문서는 실행계획이며, 아래 과거 수치를 제외한 새 모델 성능이나 목표 달성을 주장하지 않는다

| 항목 | 기준 |
|---|---|
| 전압 최종 목표 | 0.5C·1C·2C 충전과 방전에서 full physical-time RMSE ≤ 20 mV |
| 20 mV 판정 | 원칙적으로 각 cell·rate·direction 단위에서 확인하고, 평균만 통과한 경우 따로 표시 |
| 음극 전위 목표 | 기준전극 정의, 공간적 최소값, 방향별 경향, 농도 범위와 수치 수렴을 확인한 물리적으로 설명 가능한 예측 |
| 함께 확인할 지표 | 종료용량 오차, 조기 cutoff, 시작 과전압, OCP coverage, 초기 상태 의존성 |
| 추가 실험 | 요구하지 않음 |
| 이번 최적화 허용 | window offset, 사전 정의한 hysteresis 관련 소수 변수 |
| 이번 최적화 제외 | kn·kp·Dsn·Dsp·Bruggeman·전해질 물성 함수·geometry |

모든 후보가 20 mV 목표를 못 맞추면, 가장 타당한 prefit 후보와 목표까지 남은 오차를 보고한다<br>
OCP 축이나 물리 용량을 근거 없이 바꿔 목표를 달성시키지 않는다<br>
kinetic/transport fitting은 후속 단계로 분리한다

## 2 기존 코드에서 확인한 핵심 차이

### 2.1 hysteresis는 과거부터 같은 방식으로 들어간 것이 아니다

| 근거 파일 | 실제 구조 | 비교에서 주의할 점 |
|---|---|---|
| 첨부 `260818_0단계~2단계_진행_결과_RctN2반영(1).ipynb`, 셀 10·21·22 | 음극 lithiation/delithiation 평균을 중심 OCP로 사용하고 음극만 one-state hysteresis, gamma=3.2 | 양극은 평균 OCP인 single 모델 |
| `old_notebook_method_comparison.py` | 과거 RMSE window와 endpoint window, 음극 hysteresis 유무를 분리 | 과거 결과 재현용 참고 |
| `electrode_hysteresis_quantification.py` | OFF·음극만·양극만·양쪽 one-state hysteresis, gamma 기본 3.2 | 전극별 옵션 순서는 음극, 양극 |
| `nominal_radius_current_configuration.py` | Ai2020 음극 중심 OCP에 실험 방향차를 결합한 후보, 양·음극 branch gap scale=0.50 | kinetics와 물성 보정까지 포함한 builder를 그대로 사용하면 이번 비교 조건이 달라짐 |
| `c50_ocp_source_combination_study.py` | 정적 방향별 곡선에 gap scale을 최적화 | one-state 상태방정식을 푼 동적 비교가 아님 |
| `ocp_newbol_hysteresis_brugg_analysis.py` | 기존 fitted Dsn·kn을 고정한 상태에서 양극 hysteresis 비교 | 이번 nominal kinetics/transport 비교와 구분 |

이 문서의 중심 OCP는 방향별 평균을 사용하는 운영상 equilibrium 근사다<br>
저율 음극 곡선 및 1 h GITT rest-end 양극 곡선이 완전한 열역학적 평형이라는 의미는 아니다

### 2.2 과거 좋은 결과를 동일한 목표 지표로 다시 확인해야 한다

`results/260920_electrode_hysteresis_quantification/dynamic_hysteresis_summary.csv`의 저장 수치

| hysteresis | 전 구간 MAE의 조건 평균 [mV] | SOC 10–70% MAE의 조건 평균 [mV] |
|---|---:|---:|
| OFF | 56.37 | 62.57 |
| 음극만 | 29.63 | 25.50 |
| 양극만 | 48.12 | 52.43 |
| 양·음극 | 27.70 | 20.42 |

양·음극 hysteresis가 당시 조건에서는 크게 개선했다는 근거가 있다<br>
그러나 20.42 mV는 일부 SOC 구간의 평균 MAE이며, 전 구간 RMSE 20 mV 이하 달성 근거로 사용할 수 없다

`results/260921_c50_ocp_source_combinations/ocp_source_combination_summary.csv`의 저장 수치

| 음극/양극 OCP | qOCV RMSE, SOC 2–98% [mV] | C/50 branch 평균 MAE, SOC 10–90% [mV] | gap scale | Qn [Ah] | Qp [Ah] |
|---|---:|---:|---:|---:|---:|
| Ai2020/기존 실험 | 6.243 | 6.800 | 약 1.000 | 2.479 | 4.371 |
| 기존 실험/기존 실험 | 8.220 | 9.084 | 0.260 | 3.308 | 4.463 |
| Ai2020/Ai2020 | 16.354 | 17.085 | 약 1.000 | 2.400 | 5.532 |

여기서 branch 평균 MAE는 gap scale 최적화 후의 값이다<br>
실험/실험 조합의 9.084 mV를 원시 half-cell branch를 그대로 적용한 성능으로 설명하면 안 된다<br>
Ai2020/Ai2020은 방향차가 0이므로 그 경우 gap scale 자체에 식별 가능한 의미가 없다

해당 코드의 `fit_window`는 네 endpoint를 fitting하고 `Qn=Qcell/Δx`, `Qp=Qcell/Δy`로 용량을 역산한다<br>
따라서 이 용량은 OCP 조합에서 추정된 값이며 독립적인 실제 전극 용량 측정값이 아니다<br>
용량 차이는 OCP 축·inventory·window 식별이 서로 얽혀 있다는 신호로 다뤄야 한다

### 2.3 현재 값도 무조건적인 물리 prior로 쓰면 안 된다

현재 인수인계에는 `Qn=2.4794722 Ah`, `Qp=4.3711981 Ah`, 음극 `Δx=0.9504150`이 기록돼 있다<br>
그러나 `c_s,max`에는 기존 식별 용량을 geometry 변경 후 보존한 effective 값이 포함돼 있다

첨부 과거 노트북에서는 음극 공통 coverage가 약 `0–0.840542`, 양극이 약 `0.406734–0.991381`이다<br>
현재 Δx를 이 과거 축에 그대로 요구하면 음극 coverage를 초과한다<br>
이 사실은 실험 음극이 틀렸다는 증명도, Ai2020가 독립적으로 검증됐다는 증명도 아니다

**실행 원칙: OCP source별로 물리적 inventory 가정을 동일하게 유지하되, 그 공통 가정 자체의 근거와 불확실성을 먼저 감사한다**

## 3 활용할 데이터와 코드

### 3.1 OCP 데이터 계보

| 대상 | 원래 측정 파일 | 저장소에서 우선 확인할 입력 | 추출 규칙 |
|---|---|---|---|
| 실험 음극 | `260816 Anode GITT + OCP v1.xlsx`, `v2.xlsx` | `260902 Enertech셀 데이터 모음.xlsx`, `ANODE_LOWRATE` | `source_version` v1/v2, charge/discharge 이름을 실제 Li 출입 방향으로 확인 |
| 실험 양극 | `260816 Positive GITT 7-7.xlsx`, `7-8.xlsx` | 같은 묶음 파일의 `CATHODE_GITT` | `source_channel` 7-7/7-8, charge/discharge pulse의 rest-end |
| 과거 full-cell C/20 | `260808 Enertech셀 초기 저율 데이터 v2.xlsx` | 같은 묶음 파일의 `FULLCELL_LOWRATE` | charge/discharge branch 및 용량 |
| 과거 동적 데이터 | 묶음 데이터 원본 기록을 추적 | 같은 묶음 파일의 `DYNAMIC` | 0.5C·1C·2C, 충·방전, 실제 전류·시간 |
| C/50 full-cell | `260610 저율데이터 1.xlsx`, `260610 저율데이터 2.xlsx` | `c50_ocp_source_combination_study.py`의 입력 경로와 저장된 target CSV | 기존 코드가 사용한 complete cycles 2–5와 cell 구분 확인 |
| Ai2020 OCP | 버전 고정된 PyBaMM Ai2020 | installed parameter function 및 OCP CSV | `graphite_ocp_Enertech_Ai2020.csv`, `lico2_ocp_Ai2020.csv`와 일치 여부 |

실험 음극 GITT 전체를 그대로 OCP로 쓰는 방식과, 같은 workbook의 저율 full-range branch를 쓰는 방식은 다르다<br>
첫 주 비교는 과거에 사용한 음극 저율 v1/v2와 양극 GITT 7-7/7-8을 기준으로 한다

실행 전에 실제 workbook의 sheet·column·단위·source별 행 수를 확인한다<br>
부족한 원본을 가상 데이터로 채우지 않는다<br>
저장 CSV가 필요한 정보를 충분히 보존하면 활용할 수 있으나, raw 기반 분석인지 aggregate 기반 분석인지 manifest에 기록한다

### 3.2 첨부 노트북 역할

| 첨부 노트북 | 활용할 부분 |
|---|---|
| `260817 [1단계] Ai2020 vs 실험 OCP 선정 과정.ipynb` | OCP source 선택, 저율 음극과 GITT 양극 구성, capacity-constrained MAE window |
| `260818_0단계~2단계_진행_결과_RctN2반영(1).ipynb` | RMSE window, 음극 one-state hysteresis, 실제 초기 h 설정 |
| 다른 260818 결과·RctN2 변형 | source와 parameter builder 차이를 확인한 뒤 과거 재현 모델 구성 |
| `260903 0~4단계 HPPC fitting C-rate validation (3).ipynb` | 두 endpoint만 맞추는 window 방식, 기존 nominal baseline 분리 |
| `260908 Charge Fitting Discharge Validation Feasibility.ipynb` | 충전 selection/fitting과 방전 검증 분리, feasibility 처리 참고 |
| `260610_[2단계]_stoichiometry_최적화_reparameterization.ipynb` 및 최종 window 노트북 | Qn/Qp를 고정한 2변수 window 구성 |
| 260621 LSA 노트북 | 후속 kinetic/transport 단계 참고, 이번 후보에 fitted 값 주입 금지 |
| 저율 실험-시뮬레이션 오차 비교 노트북 | 비교 그림과 오차 정의 참고 |

`(최종 검증 완료) ... stoichiometry window`의 두 첨부는 SHA-256이 같고, 최종 LSA의 두 첨부도 동일하다<br>
중복 파일을 독립 실험 또는 독립적인 재현 검증으로 세지 않는다<br>
파일 이름의 ‘최종 검증 완료’만으로 현재 조건에서 재현됐다고 주장하지 않는다

### 3.3 저장소 코드 재사용 시 위험

`c50_selected_ocp_dynamic_validation.selected_inputs()`에는 current kinetics 적용이 들어 있다<br>
`restored_brugg_de_csmax_audit.set_dynamic_parameters()`에는 Dsn/Dsp 변경이 들어 있다<br>
이를 그대로 부르면 OCP-only 비교가 아니다

후속 Codex는 읽기·정리 함수와 결과 스코어 함수를 재사용하고, 이번 비교 전용의 순수 baseline builder를 새로 만든다<br>
import만으로 디렉터리 생성 등 부작용이 발생하는 기존 코드도 확인한다

## 4 Stage 0 — geometry와 공통 파라미터 잠금

### 4.1 최신 geometry의 의미부터 명시한다

저장소 9월 27일 인수인계의 참고값은 overlap area 787.236 cm², Ln/Lp 76.5/68 μm, SEM Rn/Rp 3.7542/4.2387 μm다<br>
이 값들을 오늘의 확정 최신 측정값이라고 자동 채택하지 않는다<br>
더 최근의 파일·설정이 있으면 source와 단위를 확인하고 이번 run의 공통 geometry로 고정한다

| 항목 | 필수 기록 |
|---|---|
| Ln/Lp | coating 단면 두께인지, 집전체와 양면이 어떻게 제거됐는지 |
| Rn/Rp | 반경/직경 구분, 1차/2차 입자, 평균/D50 등 집계 규칙 |
| 반응 면적 | 1D DFN overlap 면적과 전체 음극 inventory 면적 구분 |
| 병렬 전극 수 | 면적에 이미 포함했는지 확인해 이중 계산 방지 |
| εs·εe | 측정·문헌·유도값 여부, 체적분율 합 검증 |
| cs,max | 문헌 intrinsic 가정인지 용량 보존용 effective 표현인지 |

### 4.2 물리 용량 계산

공통 반응 면적을 Aeff라 하면

```text
Qn = F · csn,max · εs,n · Ln · Aeff / 3600
Qp = F · csp,max · εs,p · Lp · Aeff / 3600
Δx = Qref / Qn
Δy = Qref / Qp
```

독립 active loading이 있으면 loading × active fraction × specific capacity 계산도 교차 확인한다<br>
증거가 없는 단일 펀칭 질량을 full-cell 전체 질량으로 확대하거나, GITT 질량 보정으로 full-cell Qn/Qp를 자동 변경하지 않는다

Ai2020 OCP 함수의 사용과 Ai2020 geometry·용량의 사용은 각각 별도의 선택이다<br>
Ai2020 OCP를 후보 함수로 사용하는 것은 가능하다<br>
그 선택이 실제 셀의 용량·absolute stoichiometry·음극 전위를 검증해주지는 않는다

### 4.3 두 종류의 공통 capacity anchor

| anchor | 용도 | 적용 |
|---|---|---|
| G — geometry 및 독립 근거 기반 Qn/Qp | 우선 주 비교 | 독립 근거가 실제로 있는 범위에서 공통값을 고정 |
| E — 기존 fitted inventory 보존 Qn/Qp | 현재 모델과의 연속성 및 민감도 진단 | source별로 다시 조정하지 않고 동일한 기존값을 적용 |

G도 cs,max·active fraction을 문헌에 의존한다면 ‘독립 실측 용량’ 대신 ‘공통 물리 가정 기반 용량’으로 표시한다<br>
독립 근거가 불충분하면 E 결과를 폐기할 필요는 없지만 effective inventory 조건부 모델로 설명한다

기존 fitted Qn/Qp를 고정한 채 B가 불가능하다는 이유만으로 A의 물리적 우월성을 선언하지 않는다<br>
공통 uncertainty 범위가 기존 측정에 있으면 그 범위의 작은 격자를 동일하게 A/B에 적용해 feasible 영역을 제시한다<br>
독립 근거 없는 ±10% 또는 2.4–6 Ah 범위를 물리 prior로 만들지 않는다

### 4.4 kinetic/transport baseline 잠금

| 항목 | 이번 처리 |
|---|---|
| 음·양극 exchange-current 함수 | 버전 고정 Ai2020 nominal 함수 |
| 고체 확산 함수 | 같은 nominal 함수, GITT/fitted 값 주입 금지 |
| Bruggeman | nominal 기준값 고정, 기존 기록 2.914/1.83/1.5와 installed 값 대조 |
| 전해질 κ·De·t+ | 같은 source와 같은 단위로 고정 |
| 접촉저항·film 저항 등 | 기존 보정 여부 공개, 이번 후보별 조정 금지 |
| 온도·열 모델·solver·mesh | 모든 후보에 동일 |
| C-rate 기준용량 | 공통 Qrated 및 실제 전류를 확인, 후보별 Qn에 따라 인가전류 변경 금지 |

저장소에는 `De × 1e-4`를 source 단위변환으로 유지한 이력이 있다<br>
이 문서만으로 그 변환의 정확성을 새로 보증하지 않는다<br>
설치 버전의 원함수·원단위·기존 correction을 추적하고, 확인한 단위 보정만 공통 baseline에 넣는다<br>
`source correction`과 `성능 개선용 multiplier`를 구별하고 불명확하면 두 baseline의 별도 진단 결과로 남긴다

각 후보의 processed parameter fingerprint를 비교해 허용된 OCP·window·hysteresis 외 차이가 없는지 검증한다

## 5 Stage 1 — 실험 OCP 품질과 축 검증

### 5.1 추가 측정 없이 확인할 수 있는 재현성

기존 `electrode_branch_difference_metrics.csv`의 replicate-common-coverage 수치

| 전극과 branch | 반복쌍 MAE [mV] | 반복쌍 RMSE [mV] |
|---|---:|---:|
| 음극 lithiation v1/v2 | 2.947 | 5.038 |
| 음극 delithiation v1/v2 | 2.082 | 3.128 |
| 양극 lithiation 7-7/7-8 | 2.007 | 3.878 |
| 양극 delithiation 7-7/7-8 | 7.734 | 9.437 |

이는 기존 두 곡선 사이의 일치도를 보여준다<br>
독립 시편 수, 동일 시편 반복 여부, 세척 잔류물, 열화, 장비 offset을 모두 배제한 재현성 보증은 아니다<br>
특히 음극 약 89 mV 수준의 중앙 구간 방향차가 반복쌍 차이보다 크다고 해서 순수 intrinsic hysteresis라고 확정하지 않는다

이번 실행에서는 최신 축에서 사용 window와 공통 coverage를 기준으로 위 지표를 다시 계산한다<br>
평균 OCP뿐 아니라 v1-only/v2-only, 7-7-only/7-8-only를 최종 후보에 적용해 성능과 음극 전위 민감도를 확인한다<br>
평균곡선만 유독 좋아지거나 특정 replicate에서 문제가 커지면 재현성 제한을 명시한다

### 5.2 capacity 축과 절대 anchor

기존 음극 축은 `x=q/Qn,disk`, 역방향은 `x=qmax/Qn,disk−q/Qn,disk`다<br>
기존 양극 축은 lithiation `y=1−qmax/Qp,disk+q/Qp,disk`, delithiation `y=1−q/Qp,disk`다

이 규칙은 방향별 endpoint에 각각 anchor를 두므로 두 branch의 다른 용량이 상대 정렬에 영향을 줄 수 있다<br>
모든 측정 endpoint가 정확히 x=0 또는 y=1인 것도 자동으로 성립하지 않는다

기존 축 유지 버전과, 기존 측정으로 독립 근거가 있는 경우에 한해 재계산한 축 버전을 분리한다<br>
full-cell 전압만 맞추기 위한 자유로운 axis stretching은 주 비교에서 제외한다<br>
필요하면 affine 축 보정은 진단으로 실행하고 scale·offset·inventory 자유도 증가를 공개한다

재보정에는 각 half-cell 시편 자체의 disk geometry/loading을 쓴다<br>
새 full-cell geometry로 과거 half-cell 축을 무조건 다시 만드는 방식은 사용하지 않는다

### 5.3 데이터 처리 규칙

중복 stoichiometry, 비단조 capacity, rest-end 누락, 이상값을 확인하고 제외 규칙과 행 수를 저장한다<br>
replicate끼리 먼저 같은 branch에서 정렬·평균하고, 그 다음 lith/delith 중심과 gap을 구성한다<br>
공통 coverage 밖 보간·외삽을 금지한다<br>
단순 `np.interp`는 범위 밖 값을 endpoint로 고정할 수 있으므로 별도 bounds 검사가 필요하다

선형 보간을 공통 기준으로 사용하고, 강한 smoothing이나 고차 보간의 overshoot를 성능 개선으로 채택하지 않는다<br>
full-cell 충·방전을 각각 자기 cutoff 용량으로 0–1 정규화하면 서로 다른 상태가 같은 SOC로 묶일 수 있다<br>
기존 정규화 방식은 legacy 재현에 사용하고, 주 분석은 공통 absolute Ah 및 이력에 맞는 alignment를 우선한다<br>
alignment 자체의 민감도도 저장한다

## 6 Stage 2 — stoichiometry window 비교

### 6.1 공통 capacity가 고정된 2변수 문제

충전 기준 SOC z에 대해

```text
x(z) = x0 + z · Δx
y(z) = y100 + (1−z) · Δy
x100 = x0 + Δx
y0 = y100 + Δy
Vcenter(z) = Up,center(y(z)) − Un,center(x(z))
```

주 비교에서는 x0와 y100만 식별한다<br>
source를 바꿀 때 center OCP가 달라지므로 offset은 source별로 재식별하되 Qn/Qp와 Qref는 동일하다<br>
그 후 hysteresis OFF/ON 비교에서는 같은 source·window의 endpoint를 고정한다

각 전극 coverage가 `[xmin,xmax]`, `[ymin,ymax]`라면

```text
x0 ∈ [max(δ,xmin), min(1−δ−Δx,xmax−Δx)]
y100 ∈ [max(δ,ymin), min(1−δ−Δy,ymax−Δy)]
```

δ는 수치 singularity 방지용 작은 margin이며 물리 fitting 값이 아니다<br>
구간이 비면 `INFEASIBLE_COVERAGE`로 기록하고 optimizer를 돌리지 않는다<br>
정적 평균 stoichiometry뿐 아니라 DFN의 local surface stoichiometry도 coverage를 만족해야 한다

### 6.2 window 방법

| 코드 | 방식 | 역할 |
|---|---|---|
| W1 | full qOCV RMSE 최소화, 용량·coverage 제약 유지 | 권장 주 비교 |
| W2 | 두 endpoint 잔차 최소화, 용량·coverage 제약 유지 | 과거 endpoint 방식과 비교 |
| W3 | MAE 최소화, 같은 제약 유지 | 과거 MAE 방식 진단 |
| W4 | 네 endpoint 자유 + 용량 prior | 독립 prior가 있을 때만 별도 진단 |

W1은 원칙적으로 full curve를 사용한다<br>
측정 또는 보간 불확실성이 큰 양 끝점은 그 근거를 기록한 뒤 제외/가중할 수 있으며, 기존 SOC 2–98% 방식은 따로 재현한다<br>
full-range·2–98% 지표 모두 공개해 끝단을 숨기지 않는다

W2의 두 endpoint가 정확히 맞아도 내부 형상·inventory가 맞다는 뜻은 아니다<br>
endpoint qOCV에 저율 분극과 충·방전 alignment 오차가 포함될 수 있다<br>
정확한 두 점 통과보다 validation 성능과 음극 전위 안정성으로 선택한다

W4는 기존 9월 21일 자유 window 결과를 설명하는 데 유용하다<br>
OCP source마다 전극 용량을 자유롭게 바꾼 결과를 동일 physical model의 source-only 비교로 부르지 않는다<br>
E anchor에서 불가능한 B를 살리기 위해 B에만 W4를 허용하지 않는다

### 6.3 식별 안정성

동일 목적함수에 최소 3개 시작점 또는 seed를 사용하고, 최적점만 아니라 근접한 해의 window·전압·음극 전위 범위를 저장한다<br>
거의 같은 full-cell 오차에서 음극 전위가 크게 달라지면 비식별성을 보고한다<br>
현재 endpoint를 초기값으로 쓸 수 있으나 좁은 임의 bounds로 결과를 고정하지 않는다

## 7 Stage 3 — OCP × hysteresis 비교

### 7.1 최소 주 비교 16개

아래 8개 조합을 W1과 W2 각각에 적용한다<br>
feasible 조합만 동적 simulation으로 진행한다

| 모델 | 음극 center | 양극 center | 음극 hysteresis | 양극 hysteresis |
|---|---|---|---|---|
| A0 | Ai2020 | 기존 실험 | OFF | OFF |
| A1 | Ai2020 | 기존 실험 | OFF | ON |
| A2 | Ai2020 | 기존 실험 | ON, 실험 gap 결합 | OFF |
| A3 | Ai2020 | 기존 실험 | ON, 실험 gap 결합 | ON |
| B0 | 기존 실험 | 기존 실험 | OFF | OFF |
| B1 | 기존 실험 | 기존 실험 | OFF | ON |
| B2 | 기존 실험 | 기존 실험 | ON | OFF |
| B3 | 기존 실험 | 기존 실험 | ON | ON |

Ai2020/Ai2020 + OFF는 별도 기준 모델로 둔다<br>
G/E anchor 비교가 필요하면 동일 matrix를 anchor별로 반복하며 결과를 섞지 않는다

### 7.2 중심을 유지하는 branch 구성

각 전극에서 실험 branch의 반차를 정의한다

```text
Uexp,center = (Uexp,lith + Uexp,delith) / 2
Hexp = (Uexp,delith − Uexp,lith) / 2
Ulith,new = Ucenter − α · Hexp
Udelith,new = Ucenter + α · Hexp
```

B는 실험 center를 쓴다<br>
A의 음극은 Ai2020 center를 쓰고 측정 Hexp를 결합한다<br>
A2/A3는 실험 기반의 방향차를 이전한 가설 모델이며 ‘실험에서 측정한 Ai2020 branch’라고 부르지 않는다

raw 실험 branch를 A의 lith/delith 함수에 그대로 넣으면 평균 center가 Ai2020에서 실험으로 바뀔 수 있다<br>
반드시 `(Ulith,new+Udelith,new)/2 == Ucenter`를 수치 확인한다<br>
gap은 signed 차이로 계산하고 임의 절댓값 처리로 방향을 바꾸지 않는다

A2/A3는 실험 gap의 축 anchor와 Ai2020 축을 대응시킬 근거가 있어야 한다<br>
전체 운전 surface stoichiometry가 gap coverage를 벗어나면 주 비교에서는 infeasible로 표시한다<br>
gap 외삽이나 normalized-axis 전이는 별도 진단 모델이며, 사용하면 mapping 가정과 영향 범위를 공개한다

### 7.3 amplitude와 state dynamics를 분리

| 변수 | 1차 기준 | 후속 민감도 |
|---|---|---|
| gap scale αn·αp | ON 전극에 1, OFF에 0 | 우수 후보만 0.25·0.5·1 비교 |
| decay γn·γp | 기존 3.2 | 우수 후보만 1·3.2·10, 연구용 탐색값으로 표시 |
| 초기 h | 실제 직전 이력으로 결정 | 이력 불완전 시 h=−1·0·+1 범위 진단 |

3.2는 과거 코드의 재현 기준이며 intrinsic 물성으로 확정된 값이 아니다<br>
정적 C/50 branch gap calibration과 동적 gamma fitting은 서로 다른 작업이다<br>
특히 B의 과거 C/50 α≈0.26은 낮은 α 후보를 확인할 이유이며, 새 조건의 고정 정답이 아니다

α·γ를 추정할 때는 selection 데이터만 쓴다<br>
전극별 αn/αp 또는 γn/γp를 동시에 크게 자유화하기 전에 상관성과 식별도를 확인한다<br>
가장 단순한 공동 α 또는 고정 γ 모델을 우선하며, 개선이 명확할 때만 자유도를 늘린다

### 7.4 PyBaMM 구현과 초기 상태

모델 옵션 예시

```python
options = {
    "thermal": "isothermal",
    "open-circuit potential": (
        "one-state hysteresis" if hyst_n else "single",
        "one-state hysteresis" if hyst_p else "single",
    ),
}
```

설치 버전의 실제 상태방정식과 branch 혼합식을 확인한다<br>
설명용 `Ucenter+hH`의 h 부호를 라이브러리에 그대로 가정하지 않는다<br>
작은 charge/discharge simulation으로 음극·양극이 올바른 lith/delith branch로 수렴하는지 검증한다

과거 builder의 초기 h는 charge에서 음극 −1, 양극 +1, discharge에서 음극 +1, 양극 −1이다<br>
후속 코드에는 직전 반대 방향 이력을 반영하려고 이를 뒤집은 경우가 있다<br>
어느 쪽이 작은 오차를 만드는지 branch별로 선택하는 방식은 금지한다

가장 좋은 처리 방식은 기록된 선행 charge/discharge/rest를 연속 simulation해 h와 농도장을 전달하는 것이다<br>
현재 branch 방향만으로 직전 h를 초기화하지 않는다<br>
휴지시간이 길다고 h=0 또는 균일 농도가 보장되는지도 가정하지 않는다

선행 이력이 부족하면 공통의 대표 초기화 규칙을 선택하고 가능한 h 범위에서 결과를 보고한다<br>
branch마다 h0를 voltage 전체에 fitting해 얻은 성능은 주 비교가 아닌 진단이다

## 8 Stage 4 — 공정한 동적 검증

### 8.1 데이터 분할

| 데이터 | 기본 용도 |
|---|---|
| C/50 저율 full-cell | OCP·window 식별, α의 정적 진단 |
| C/20 저율 full-cell | 기존 결과 재현 및 rate transfer 확인 |
| July 동적 charge | 후보 selection, 필요한 경우 γ의 제한적 calibration |
| 같은 cohort의 discharge | calibration에서 제외한 반대 방향 검증 |
| 독립 cell·cycle·September cohort | 가능한 범위에서 최종 transfer 검증 |

파일 존재와 데이터 양을 확인한 뒤 실행 전 split manifest를 고정한다<br>
같은 곡선의 인접 점을 무작위로 나누지 않고 cell 또는 cycle 단위로 분리한다<br>
validation을 보며 다시 후보를 고르면 그 데이터는 selection에 사용된 것으로 표시한다<br>
독립 test가 없으면 내부 검증이라고 보고하고 일반화 성능을 확정하지 않는다

다른 cohort의 inventory/window는 같다고 가정하지 않는다<br>
필요한 cohort별 저율 window 재식별은 저율 데이터만 사용하고, 재식별 유무를 별도 기록한다

### 8.2 실제 시간·전류·초기 상태

인가전류는 실제 기록을 우선하고, 고정 C-rate 명령을 쓰면 그 기준용량을 공개한다<br>
모델과 실험은 current onset을 같은 t=0으로 정렬한다<br>
모델별로 종료시간이나 cutoff 용량을 0–1로 정규화해 전압 오차를 평가하지 않는다

SOC0는 coulomb counting과 알려진 이력을 우선한다<br>
직전 rest 전압만 사용할 경우 해당 모델의 center/h 상태와 함께 SOC0를 추정하고 ambiguity를 기록한다<br>
rest 전압 역산의 SOC0를 후보마다 공개하고 과도하게 변하는 경우 판정에 반영한다<br>
전체 동적 곡선을 맞추도록 SOC0를 자유 fitting하지 않는다

공통 구간 점수는 진단으로 제공할 수 있다<br>
그러나 model cutoff가 실험 평가 종료 전에 발생하면 full-branch 성능은 infeasible로 표시한다<br>
잔여 실험 구간을 버리고 작은 RMSE만 보고하지 않는다

### 8.3 지표 정의

```text
ei = 1000 · [Vmodel(ti) − Vexp(ti)]
RMSE = sqrt(sum(wi · ei²) / sum(wi))
MAE = sum(wi · abs(ei)) / sum(wi)
Bias = sum(wi · ei) / sum(wi)
Qerror[%] = 100 · (Qmodel,end − Qexp,end) / Qexp,end
```

시간가중 적분 또는 같은 physical-time grid를 사용해 raw sampling rate의 영향을 통제한다<br>
rate·direction·cell을 동등 가중한 macro RMSE를 함께 보고한다<br>
capacity는 실제 전류 적분의 unbounded transferred Ah를 사용하고 plot용 clipped SOC로 복원하지 않는다

| 필수 보고 | 목적 |
|---|---|
| full physical-time RMSE·MAE·Bias | 최종 20 mV 목표 판정 |
| 조건별 RMSE와 worst-cell/worst-branch | 평균에 가려지는 실패 확인 |
| SOC 2–98%, 10–70% 지표 | 과거 비교 연결, 최종 기준 대체 금지 |
| 첫 10 s·60 s 오차 | 초기 h·농도·kinetic mismatch 진단 |
| 종료용량·시간 오차 | 형상만 잘 맞는 후보 배제 |
| 평가 가능한 시간 비율 | 조기 cutoff로 인한 유리한 점수 방지 |

용량오차 3%는 연구용 우선 기준으로 둘 수 있으나 보편적 합격 규칙은 아니다<br>
기존 데이터의 cell 간 편차가 확인되면 그에 맞춰 실행 전에 기준을 정하고, 결과를 본 뒤 문턱을 완화하지 않는다

## 9 음극 전위 검증 — full-cell 전압과 별도로 수행

### 9.1 무엇을 저장할 것인가

| 변수 | 의미 |
|---|---|
| Un,center 및 Un,hyst | equilibrium 근사와 hysteresis 적용 OCP |
| φs,n−φe,n | local negative electrode interfacial potential difference |
| ηn=(φs,n−φe,n)−Un,hyst | reaction overpotential, 모델 정의와 대조 |
| separator interface 전위 | 고율 충전에서 주의할 위치의 추세 |
| 전 음극 공간 최소값 | 공간평균에 숨겨지는 낮은 전위 확인 |
| surface stoichiometry | 0–1과 OCP·gap coverage 확인 |
| electrolyte concentration | 양수 및 수치적으로 정상인지 확인 |

PyBaMM 변수 이름은 실행 버전에서 검색한다<br>
기존 `gitt_ocp_analysis.py`는 `Negative electrode surface potential difference at separator interface [V]`를 사용한다<br>
새 runner에서는 평균만 아니라 시간×공간 배열을 저장한다<br>
기존 hysteresis runner는 surface stoichiometry 최소·최대만 반환하므로 전위 저장을 추가해야 한다

### 9.2 판정 규칙

0 V vs Li/Li+ 기준의 해석이 일관되는지 먼저 확인한다<br>
충전에서 음극 lithiation이 진행될 때 OCP 변화와 cathodic η가 설명 가능한지, 방전에서는 그 반대 경향이 설명 가능한지 확인한다

local φs−φe가 0 V 아래로 내려가는 것은 plating 위험의 진단 신호다<br>
그 값이 무조건 수학적으로 불가능하거나 모델 오류라는 뜻은 아니다<br>
반대로 전위가 높거나 모두 양수라는 사실만으로 예측이 정확하다고 판단하지 않는다

사용자의 ‘정상적 음극 전위’ 목표는, 안전한 충전 조건을 설명하려는 경우 전 공간 최소값이 양수인지도 확인하는 것으로 운영한다<br>
0 V 통과가 있으면 시간·위치·rate·지속시간과 초기화/mesh 민감도를 보고하고 근거 없이 채택하지 않는다<br>
전압을 잘 맞추기 위해 OCP offset을 조정해 전위를 강제로 양수화하지 않는다

같은 full-cell 전압을 맞추는 여러 window 해가 다른 음극 전위를 예측하면 모두 저장하고 범위를 보여준다<br>
reference electrode 측정이 없으므로 최종 표현은 ‘물리적 일관성과 강건성을 확인한 음극 전위 예측’이다<br>
‘음극 전위 실측 검증 완료’ 또는 ‘정확한 절대 전위 보장’이라고 쓰지 않는다

## 10 최종 후보 선택 순서

| 순서 | 조건 |
|---|---|
| 1 | 데이터 source·축·단위·geometry와 파라미터 fingerprint가 설명 가능한가 |
| 2 | coverage·농도·h 범위·solver 수렴·charge conservation을 통과하는가 |
| 3 | 초기 h·replicate·근접 window 해에서 음극 전위와 voltage가 안정적인가 |
| 4 | 선택에 쓰지 않은 branch/cell에서 full physical-time 전압 성능이 좋은가 |
| 5 | 종료용량과 조기 cutoff가 함께 허용 가능한가 |
| 6 | 20 mV 목표를 만족하는 후보 중 자유도와 가정이 가장 적은가 |

먼저 validity와 제한을 판정하고, 통과한 후보 중 validation 성능으로 순위를 정한다<br>
현재 모델을 반드시 유지하거나 실험/실험 모델을 반드시 선택한다는 결론을 미리 두지 않는다

가능한 최종 판정

| 판정 | 의미 |
|---|---|
| `TARGET_MET` | 사전에 정한 모든 평가 조건에서 20 mV·물리 검증·용량 기준 만족 |
| `BEST_VALID_PREFIT` | 유효한 후보 중 최선이나 20 mV 목표 미달 |
| `CONDITIONAL_EFFECTIVE_MODEL` | fitted inventory 또는 axis/gap transfer 가정에 조건부인 후보 |
| `INFEASIBLE_COVERAGE` | OCP/gap 범위와 필요한 운전 범위 불일치 |
| `EARLY_CUTOFF` | 실험 평가시간을 충분히 커버하지 못함 |
| `NONIDENTIFIABLE_ANODE_POTENTIAL` | 비슷한 voltage 해가 서로 다른 음극 전위를 예측 |
| `NOT_RUN_MISSING_INPUT` | 필요한 기존 데이터 부재로 미실행 |

α=0에 가까운 dual-hysteresis 모델은 사실상 OFF 모델과의 차이를 확인한다<br>
차이가 cell 간 편차나 수치 오차보다 작으면 단순 모델을 선택한다

## 11 구현 순서와 산출물

새 작업 폴더 제안: `results/261011_ocp_window_hysteresis_comparison/`<br>
기존 결과는 덮어쓰지 않는다

| 순서 | 작업 | 산출물 |
|---|---|---|
| 1 | 실제 branch·HEAD·AGENTS·입력 inventory 확인 | `input_manifest.json` |
| 2 | historical builder 차이와 baseline 잠금 | `historical_method_audit.csv`, `locked_baseline.json` |
| 3 | raw branch 정리·replicate·axis·coverage QC | `ocp_source_qc.csv`, `ocp_axis_audit.csv`, `ocp_branches.csv` |
| 4 | capacity anchor와 feasible window 확인 | `capacity_anchor_audit.csv`, `coverage_feasibility.csv` |
| 5 | A/B × W1/W2 window 식별 | `window_comparison.csv`, `window_near_optimal_solutions.csv` |
| 6 | A0–A3/B0–B3 동일 baseline 동적 비교 | `voltage_metrics_by_condition.csv`, `capacity_metrics.csv` |
| 7 | 우수 후보 α·γ·h·replicate 민감도 | `hysteresis_sensitivity.csv`, `replicate_transfer_metrics.csv` |
| 8 | local 음극 전위와 수치 검증 | `anode_potential_metrics.csv`, `numerical_convergence.csv` |
| 9 | 선정과 근거 정리 | `decision_table.csv`, `comparison_report.md`, `selected_model.json` |

위 파일들은 새 구현의 예상 산출물이며 현재 생성됐다는 의미가 아니다

예상 runner: `compare_ocp_window_hysteresis.py`<br>
설정: `configs/ocp_window_hysteresis_comparison.json`<br>
코드는 저장소 내부 구조와 AGENTS 지침을 확인한 뒤 그에 맞춰 배치한다

```text
입력 감사 → baseline 잠금 → OCP QC → capacity/coverage 검사
→ window 식별 → frozen-window hysteresis 비교
→ 우수 후보의 상태·분기차 민감도 → 음극 전위 검증 → 선정
```

legacy 재현은 별도 테이블에 둔다<br>
과거 geometry·kinetics·window·h를 그대로 쓰는 재현 모델과 최신 geometry의 주 비교 모델을 같은 순위표에 합치지 않는다

### 최소 수치 검증

center 보존과 branch 방향 수렴을 먼저 검증한다<br>
최종 1–2개 후보에서 mesh를 한 단계 늘리고 tolerance를 조여 terminal voltage와 local minimum negative potential 변화를 확인한다<br>
약 1 mV보다 큰 수치 변화 또는 0 V 통과 여부 변화가 있으면 해결 전 결론을 확정하지 않는다<br>
1 mV는 이번 20 mV 목표를 위한 작업용 수치 기준이며 보편 규칙이 아니다

## 12 최종 결과 표 양식

| 모델 | anchor | window | αn/αp | γn/γp | 초기 이력 | full validation RMSE | worst branch RMSE | 최대 용량오차 | 최소 음극 전위 | coverage | 판정 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A0 | 미실행 | 미실행 | 0/0 | 해당 없음 | 미실행 | 미실행 | 미실행 | 미실행 | 미실행 | 미실행 | 미실행 |
| A3 | 미실행 | 미실행 | 미실행 | 미실행 | 미실행 | 미실행 | 미실행 | 미실행 | 미실행 | 미실행 | 미실행 |
| B0 | 미실행 | 미실행 | 0/0 | 해당 없음 | 미실행 | 미실행 | 미실행 | 미실행 | 미실행 | 미실행 | 미실행 |
| B3 | 미실행 | 미실행 | 미실행 | 미실행 | 미실행 | 미실행 | 미실행 | 미실행 | 미실행 | 미실행 | 미실행 |

음극 전위는 mV vs Li/Li+로 표시하고 해당 시간·공간 위치를 함께 제공한다<br>
모든 A/B 중간 모델도 실제 결과 표에는 포함한다

선정 이유는 다음 연결이 드러나도록 작성한다

```text
OCP source를 고른 데이터 근거
→ capacity/axis 가정과 불확실성
→ window 방법을 고른 검증 근거
→ hysteresis 추가가 필요한 방향별 오차
→ 초기 상태와 음극 전위의 강건성
→ 독립 또는 내부 검증에서 실제 목표 달성 여부
```

## 13 Codex에 전달할 실행 요청문

```text
이 문서를 기준으로 ai2020_dfn_fitting에서 OCP·stoichiometry window·hysteresis
비교 코드를 구현하고, 이용 가능한 기존 데이터로 실행해줘

추가 실험은 하지 않아
목표는 전압 full physical-time RMSE 20 mV 이하와 설명 가능한 음극 전위 예측이야
최신 geometry는 모든 주 비교에 동일하게 반영하고, kinetic/transport는
실험값 또는 fitted 값을 넣지 않은 공통 nominal baseline으로 잠가줘

먼저 실제 HEAD·AGENTS·입력 파일·parameter builder를 확인해줘
현재 Qn/Qp와 cs,max가 독립 측정값인지 fitted inventory 표현인지 구별해줘
Ai2020 음극/실험 양극과 실험 음극/실험 양극을 비교하되,
각 source의 중심 OCP를 유지하고 W1 RMSE와 W2 endpoint window를 비교해줘
window는 source·방법별로 식별하고 같은 window의 hysteresis OFF/ON에서는 고정해줘

A0–A3, B0–B3를 각각 W1/W2에서 비교해줘
Ai2020 중심에 실험 음극 gap을 결합하는 A2/A3는 가설 모델로 표시하고
gap coverage와 축 정렬 근거가 부족하면 주 비교에서 infeasible로 남겨줘
우수 후보만 amplitude, gamma, 초기 h와 replicate 민감도를 추가로 확인해줘
실제 직전 이력으로 h와 농도 상태를 전달하는 방식을 우선해줘

cell·rate·direction별 실제 시간축 전압 RMSE, cutoff 용량오차,
local 최소 음극 전위, surface stoichiometry, OCP coverage를 저장해줘
조기 cutoff 구간을 제거하거나 각 curve를 cutoff 기준으로 다시 정규화해
성능을 좋게 보이게 하지 마

현재 결과를 덮어쓰지 말고 새 날짜 폴더에 manifest·config·곡선·지표·선정 근거를 남겨줘
선정에 사용한 데이터와 validation 데이터를 구분해줘
원본이 없으면 저장 CSV의 정보 충분성을 확인하고, 불가능한 항목은 미실행으로 남겨줘
모두 20 mV를 넘으면 가장 타당한 prefit 후보와 남은 오차를 보고해줘
추가 kinetic/transport fitting은 이번 작업과 분리해줘
음극 전위를 실측 검증했다고 표현하지 말고 예측의 일관성과 불확실성을 보고해줘
```

## 14 근거와 확인 범위

아래 저장소 파일·결과 CSV와 첨부 노트북의 source 및 저장 출력을 읽어 이 명세를 작성했다<br>
이번 작성에서는 새 DFN simulation·optimizer·replicate 재분석을 실행하지 않았다<br>
외부 경로를 참조하는 C/50·July·September 원본이 모두 현재 환경에 있다는 의미도 아니다

| 근거 | 링크 또는 위치 |
|---|---|
| C/50 비교 코드 | [c50_ocp_source_combination_study.py](https://github.com/Rinny-Hi/ai2020_dfn_fitting/blob/aa003e9bd2d058de3672ad514043f79bf1ade8c7/c50_ocp_source_combination_study.py) |
| C/50 저장 결과 | [ocp_source_combination_summary.csv](https://github.com/Rinny-Hi/ai2020_dfn_fitting/blob/aa003e9bd2d058de3672ad514043f79bf1ade8c7/results/260921_c50_ocp_source_combinations/ocp_source_combination_summary.csv) |
| hysteresis 비교 코드 | [electrode_hysteresis_quantification.py](https://github.com/Rinny-Hi/ai2020_dfn_fitting/blob/aa003e9bd2d058de3672ad514043f79bf1ade8c7/electrode_hysteresis_quantification.py) |
| 방향차·반복쌍 저장 결과 | [electrode_branch_difference_metrics.csv](https://github.com/Rinny-Hi/ai2020_dfn_fitting/blob/aa003e9bd2d058de3672ad514043f79bf1ade8c7/results/260920_electrode_hysteresis_quantification/electrode_branch_difference_metrics.csv) |
| hysteresis 동적 저장 결과 | [dynamic_hysteresis_summary.csv](https://github.com/Rinny-Hi/ai2020_dfn_fitting/blob/aa003e9bd2d058de3672ad514043f79bf1ade8c7/results/260920_electrode_hysteresis_quantification/dynamic_hysteresis_summary.csv) |
| 과거 방식 분리 비교 | [old_notebook_method_comparison.py](https://github.com/Rinny-Hi/ai2020_dfn_fitting/blob/aa003e9bd2d058de3672ad514043f79bf1ade8c7/old_notebook_method_comparison.py) |
| 현재 모델과 effective inventory 설명 | [260927_HANDOFF_CURRENT_STATE_KO.md](https://github.com/Rinny-Hi/ai2020_dfn_fitting/blob/aa003e9bd2d058de3672ad514043f79bf1ade8c7/260927_HANDOFF_CURRENT_STATE_KO.md) |
| 후속 fitted 조건의 hysteresis 결과 | [README_KO.md](https://github.com/Rinny-Hi/ai2020_dfn_fitting/blob/aa003e9bd2d058de3672ad514043f79bf1ade8c7/results/260927_ocp_newbol_hysteresis_brugg_comparison/README_KO.md) |
| PyBaMM one-state 모델 정의 | [공식 문서 v26.4.0](https://docs.pybamm.org/en/v26.4.0/source/api/models/submodels/interface/open_circuit_potential/one_state_hysteresis_ocp.html) |

공식 문서는 one-state hysteresis가 local surface stoichiometry와 volumetric interfacial current에 의존하는 상태방정식을 사용하는 모델이라고 설명한다<br>
실제 실행에서는 저장소 requirements와 설치된 PyBaMM 버전을 확인하고, 그 버전의 구현을 기준으로 옵션·파라미터·변수 이름을 검증한다
