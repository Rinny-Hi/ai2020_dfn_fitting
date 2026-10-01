# 파일 지도

이 문서는 저장소 **현재 배치**를 기준으로 작성했습니다. 연구 순서는 파일명의 날짜와 항상 같지 않습니다. 최신 결정은 루트의 [`260927_HANDOFF_CURRENT_STATE_KO.md`](../260927_HANDOFF_CURRENT_STATE_KO.md)와 [`260927_charge_physical_time_final/README_KO.md`](../results/260927_charge_physical_time_final/README_KO.md)를 우선합니다.

## 디렉터리와 입력

| 경로 | 역할 | 재실행 조건 |
| --- | --- | --- |
| [`260902 Enertech셀 데이터 모음.xlsx`](../260902%20Enertech셀%20데이터%20모음.xlsx) | 통합 half-cell OCP/GITT, full-cell C/20 qOCV, 동적 데이터, prior | 저장소에 포함. 여러 초기 노트북과 GITT 스크립트가 루트에서 읽음 |
| [`data/raw/260918_Anode_GITT/`](../data/raw/260918_Anode_GITT/) | 음극 GITT v2, v3, 결합 원본 | 저장소에 포함. `gitt_ocp_qc.py`, `gitt_ocp_analysis.py` 입력 |
| `260808 Eneterch셀 열화데이터 v1 6-6(5).xlsx` | 초기 HPPC 노트북의 BoL HPPC 원본 | 저장소에 없음. 초기 노트북 실행 위치에 별도 준비 |
| July BoL 3셀, September BoL 2셀 원본 Excel | 후속 charge/discharge 피팅·cohort 비교 입력 | 저장소에 없음. 재실행 담당 에이전트가 파일을 찾아 `OLD_FILES`·`NEW_FILES` 등과 연결해야 함 |
| [`results/`](../results/) | 실행별 CSV, JSON, PNG, Markdown 보고서 | 과거 실행 산출물 포함. 재실행 시 같은 폴더를 덮어쓸 수 있음 |
| [`outputs/`](../outputs/) | 공유용 결과 압축 파일 | 분석의 주 산출물 위치는 `results/` |

통합 Excel의 `SOURCE_INDEX` 시트에서 묶인 원본을 확인할 수 있습니다. `data/raw/`의 세 원본은 새 음극 GITT QC용입니다. 저장소에 없는 외부 원본을 통합 Excel의 값으로 임의 대체하지 마세요.

스크립트가 읽는 CSV는 대부분 `results/`의 이전 단계 산출물입니다. 로컬 CSV 경로를 사용자에게 묻기 전에 저장소에서 파일명과 생산 스크립트를 검색하세요. 경로 추적과 수정 원칙은 루트 [`AGENTS.md`](../AGENTS.md)에 적었습니다.

## 노트북 4개

| 파일 | 성격 | 필요한 입력 |
| --- | --- | --- |
| [`260903 0~4단계 HPPC fitting C-rate validation.ipynb`](../260903%200~4단계%20HPPC%20fitting%20C-rate%20validation.ipynb) | 초기 Stage 0–4: OCP, HPPC 민감도·피팅, C-rate 검증 | 루트 통합 Excel과 저장소에 없는 HPPC 원본. 첫 셀의 PyBaMM 버전은 당시 환경 |
| [`260908 Charge Fitting Discharge Validation.ipynb`](../260908%20Charge%20Fitting%20Discharge%20Validation.ipynb) | 초기 charge-only 피팅, discharge 검증 | 루트 통합 Excel. 기본값은 `SMOKE_TEST=True` |
| [`260918 GITT OCP QC and Stage1 Preparation.ipynb`](../260918%20GITT%20OCP%20QC%20and%20Stage1%20Preparation.ipynb) | 새 음극 GITT의 파싱·QC 결과 확인 | `gitt_ocp_qc.py` 및 포함된 raw 원본 |
| [`260918 GITT OCP Stage1 Comparison.ipynb`](../260918%20GITT%20OCP%20Stage1%20Comparison.ipynb) | legacy/new OCP 및 Stage 1 비교 결과 열람 | `gitt_ocp_analysis.py`가 만든 `results/260918_gitt_ocp_comparison/` |

이 노트북들은 당시 실험 기록입니다. `260927` 최종 physical-time 피팅의 실행 파일은 아래 Python 스크립트입니다.

## 최신 July 피팅 흐름

```text
통합 Excel + 외부 July BoL 원본
  └─ comprehensive_prefit_cross_cohort.py / 기존 OCP·geometry 분석
      └─ charge_physical_time_common.py (공통 입력·시간 grid·목적함수)
          ├─ p0_charge_physical_time_lsa.py (사전 식별도)
          └─ p0_charge_subset_fit_validation.py (charge 피팅 / discharge 검증)
              └─ p0_charge_postfit_identifiability_audit.py (사후 식별도)
```

| 파일 | 주요 산출물·의존성 |
| --- | --- |
| [`charge_physical_time_common.py`](../charge_physical_time_common.py) | 100개 physical-time point, 파라미터 범위, 시뮬레이션·잔차 공통 함수. July 원본을 `comprehensive_prefit_cross_cohort.py`를 통해 읽음 |
| [`p0_charge_physical_time_lsa.py`](../p0_charge_physical_time_lsa.py) | `results/260927_charge_physical_time_lsa/subset_identifiability.csv` 등 |
| [`p0_charge_subset_fit_validation.py`](../p0_charge_subset_fit_validation.py) | 사전 식별도 CSV를 읽고 주/보수 모델을 피팅. `results/260927_charge_subset_fit_validation/`에 optimizer·검증 결과 저장 |
| [`p0_charge_postfit_identifiability_audit.py`](../p0_charge_postfit_identifiability_audit.py) | 선택 후보의 경계 근처 Jacobian 재검사 |
| [`p0_charge_dsn_boundary_test.py`](../p0_charge_dsn_boundary_test.py) | `Dsn` 상한 확장 비교 |
| [`results/260927_charge_physical_time_final/README_KO.md`](../results/260927_charge_physical_time_final/README_KO.md) | 위 단계의 최종 해석과 선택 근거 |

[`p0_charge_fit_discharge_validation.py`](../p0_charge_fit_discharge_validation.py)는 **중단된 5파라미터 동시 피팅**입니다. 해당 결과의 중간 optimizer 이력은 최종 파라미터가 아닙니다.

## 검증 및 후속 작업

| 작업 | 실행 파일 | 저장된 판단 근거 |
| --- | --- | --- |
| 음극 GITT 원본 파싱·QC | [`gitt_ocp_qc.py`](../gitt_ocp_qc.py) | [`260918_gitt_ocp_qc/qc_manifest.json`](../results/260918_gitt_ocp_qc/qc_manifest.json) |
| legacy/new GITT OCP 비교 | [`gitt_ocp_analysis.py`](../gitt_ocp_analysis.py) | [`260918_gitt_ocp_comparison/comparison_report.md`](../results/260918_gitt_ocp_comparison/comparison_report.md) |
| 9월 raw 곡선 추출·QC | [`september_raw_fullcurve_audit.py`](../september_raw_fullcurve_audit.py) | [`260927_september_raw_fullcurve_audit/README_KO.md`](../results/260927_september_raw_fullcurve_audit/README_KO.md) |
| 9월 Dsn 식별도와 kn 피팅 | [`september_recommended_dsn_kn_fit.py`](../september_recommended_dsn_kn_fit.py) | [`260927_september_recommended_dsn_kn_fit/README_KO.md`](../results/260927_september_recommended_dsn_kn_fit/README_KO.md) |
| 새 2셀 RPT geometry 검증 | [`validate_260927_two_cell_rpt.py`](../validate_260927_two_cell_rpt.py) | [`260927_two_cell_rpt_validation/README_KO.md`](../results/260927_two_cell_rpt_validation/README_KO.md) |
| July 최종값의 Ai2020 nominal 복원 비교 | [`july_one_at_a_time_nominal_reversion.py`](../july_one_at_a_time_nominal_reversion.py) | [`260928_july_nominal_reversion_audit/README_KO.md`](../results/260928_july_nominal_reversion_audit/README_KO.md) |

여기서 “검증”은 측정 데이터와 모델 결과를 비교하는 연구 절차입니다. 별도의 `pytest`/`unittest` 회귀 테스트 모음은 없습니다.

## 나머지 루트 스크립트 찾기

루트 `.py` 파일은 서로 import하는 연구 스크립트이므로 물리적으로 폴더를 나누지 않았습니다. 이름으로 찾을 때 아래 묶음을 사용하세요.

- **OCP·GITT·hysteresis:** `c50_*`, `gitt_*`, `ocp_*`, `legacy_ocp_window_recheck.py`, `mass_independent_nominal_anode_trial.py`, `finalize_ocp_recommendation.py`, `electrode_hysteresis_quantification.py`, `hysteresis_initial_state_sensitivity.py`, `final_endpoint_dual_hysteresis_baseline.py`, `absolute_capacity_esoh_recheck.py`, `literature_esoh_multiterm_fit.py`
- **Geometry·EIS·transport:** `paper_*`, `eis_*`, `*thickness*`, `*radius*`, `*mass*`, `*brugg*`, `*kinetics*`, `*diffusion*`, `*ds_*`, `transport_parameter_ablation.py`, `sem_geometry_capacity_ablation.py`
- **용량·cohort·초기 상태:** `capacity_*`, `comprehensive_prefit_cross_cohort.py`, `priority_initial_state_protocol_recheck.py`, `sequential_rpt_recheck.py`, `continuous_protocol_probe.py`, `current_p0_curve_and_r0_diagnostic.py`
- **민감도·피팅:** `p0_*`, `stage2_seven_parameter_sensitivity.py`, `stage3_kn_charge_fit_discharge_validation.py`, `september_*`
- **요약·그림:** `build_prefit_summary.py`, `build_gitt_qc_notebook.py`, `create_gitt_and_lsa_slide_figures.py`, `plot_july_*.py`, `check_july_kn_2p59e7.py`

결과 폴더는 `YYMMDD_주제` 형식입니다. 해당 폴더의 `README_KO.md` 또는 `analysis_manifest.json`·`*_report.json`에서 입력, 조건, 해석을 확인하세요. 앞선 결과 문서와 충돌할 때는 최신 인수인계 및 최종 결과 문서를 우선합니다.
