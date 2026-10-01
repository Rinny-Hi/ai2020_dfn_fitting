# Enertech 셀 DFN 파라미터 피팅

Enertech 파우치셀의 OCP, GITT, EIS, HPPC, 충·방전 데이터를 이용해 PyBaMM의 Ai2020 DFN 기준 모델을 검토하고 파라미터를 식별한 연구 저장소입니다. 초기 노트북, 후속 분석 스크립트, 원본 데이터 일부, 날짜별 결과를 함께 보존합니다.

**현재 결과를 읽는 출발점은 [2026-09-27 인수인계](260927_HANDOFF_CURRENT_STATE_KO.md)와 [최종 physical-time 피팅 결과](results/260927_charge_physical_time_final/README_KO.md)입니다.** 이전 노트북과 결과 폴더는 실험 당시의 가정을 기록하므로 최신 권장 설정과 다를 수 있습니다.

## 현재 권장 설정

2026-09-28까지 저장된 결과 기준, July BoL 3셀의 0.5C·1C·2C **charge** 전압으로 `Dsn`과 `kn`을 피팅하고 **discharge**를 독립 검증한 보수 모델이 선택됐습니다. `Dsn=7.4662e-14 m²/s`, `kn=1.45024e-6`, `brugg_n=2.914`(고정)입니다. Charge objective RMSE는 53.997 mV, discharge full RMSE는 49.51 mV입니다. 용량은 목적함수에 넣지 않고 사후에 평가했습니다. 값의 의미와 한계는 [최종 결과 문서](results/260927_charge_physical_time_final/README_KO.md)를 확인하세요.

September BoL 2셀은 [별도 검증](results/260927_september_recommended_dsn_kn_fit/README_KO.md)입니다. 해당 데이터에서는 `Dsn+kn` 동시 피팅의 사전 식별도 기준을 통과하지 못해 GITT 기반 `Dsn`을 고정하고 `kn`만 피팅했습니다. July 결과를 September 셀의 확정 물성값으로 옮겨 쓰지 않습니다.

## 파일 구성

```text
ai2020_dfn_fitting/
├── README.md                         프로젝트 안내와 실행 경로
├── AGENTS.md                         후속 Codex/에이전트의 검증·수정 지침
├── docs/FILE_MAP.md                  노트북·스크립트·검증 파일 상세 지도
├── requirements-gitt-ocp.txt        후속 Python 분석 환경의 버전 목록
├── 260902 Enertech셀 데이터 모음.xlsx  통합 OCP/GITT/qOCV/동적 데이터
├── 2609*.ipynb                      초기 피팅 및 GITT 검토 노트북 4개
├── *.py                              분석, 피팅, 검증 및 그림 생성 스크립트
├── 2609*.md                          날짜별 진행·인수인계 기록
├── data/raw/260918_Anode_GITT/      저장소에 포함된 음극 GITT 원본 3개
├── results/                          분석별 CSV·JSON·그림·보고서
└── outputs/                          별도 공유용 결과 압축 파일
```

루트의 Python 파일은 서로 같은 폴더에서 import하고 일부는 `__file__` 또는 현재 작업 디렉터리를 기준으로 입력·결과 경로를 만듭니다. 파일을 다른 폴더로 옮기려면 import와 경로를 함께 수정해야 합니다. 역할별 핵심 파일, 입력 조건, 산출물은 [파일 지도](docs/FILE_MAP.md)에 정리했습니다.

## 설치와 입력 확인

Python 3.12를 권장합니다. 저장소 루트에서 가상 환경을 만들고 의존성을 설치하세요. `requirements-gitt-ocp.txt`는 후속 분석에서 사용한 버전(`pybamm==26.8.0.0` 포함)을 기록합니다.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-gitt-ocp.txt
```

Windows PowerShell에서는 `python -m venv .venv` 후 `.venv\Scripts\Activate.ps1`을 실행합니다. 초기 `260903` 노트북의 첫 셀에는 당시 사용한 `pybamm==26.7.1.0` 설치 명령이 별도로 남아 있습니다. 후속 분석 버전에서 그 노트북을 재실행하려면 환경 차이를 먼저 확인하세요.

통합 Excel의 `SOURCE_INDEX`는 원본 출처를, `ANODE_LOWRATE`, `ANODE_GITT`, `CATHODE_GITT`, `CATHODE_LOWRATE`, `FULLCELL_LOWRATE`, `PARAMETER_PRIOR`, `DYNAMIC`은 분석 입력을 담습니다. 통합 Excel은 **저장소 루트**에 있습니다. `260908 Charge Fitting Discharge Validation.ipynb`도 이 위치를 읽습니다.

## 실행 경로

### 저장소에 포함된 입력으로 실행

아래 명령은 저장소 루트에서 실행합니다. GITT QC에는 `data/raw/260918_Anode_GITT/`의 v2·v3 원본과 루트의 통합 Excel이 필요합니다.

```bash
python gitt_ocp_qc.py
```

QC 결과는 `results/260918_gitt_ocp_qc/`에 기록됩니다. 같은 내용을 단계별로 보려면 `260918 GITT OCP QC and Stage1 Preparation.ipynb`를 사용하세요. `new_2h_candidate`는 절대 stoichiometry 기준이 확정되지 않아 Stage 1 투입 준비 완료 상태가 아닙니다. [QC manifest](results/260918_gitt_ocp_qc/qc_manifest.json)에 범위가 명시돼 있습니다.

`260908 Charge Fitting Discharge Validation.ipynb`는 통합 Excel만 사용하며 기본 `SMOKE_TEST=True`입니다. 이 설정의 출력은 실행 경로 확인용이며 최종 피팅 결과가 아닙니다. 전체 실행을 하려면 `SMOKE_TEST=False`로 바꾸고, 결과를 현재 권장 모델과 별도로 비교하세요.

### 추가 원본이 필요한 분석

`260903 0~4단계 HPPC fitting C-rate validation.ipynb`의 Stage 2–3은 저장소에 없는 `260808 Eneterch셀 열화데이터 v1 6-6(5).xlsx`가 필요합니다. 노트북 설명대로 통합 Excel과 함께 작업 위치에 두어야 합니다.

최신 July/September 피팅을 다시 계산하는 스크립트는 저장소에 없는 July BoL 3셀 및 September BoL 2셀의 원본 Excel을 참조합니다. 일부 스크립트에는 작성자 PC의 Windows 절대경로가 남아 있습니다. `comprehensive_prefit_cross_cohort.py`의 `OLD_FILES`/`NEW_FILES`와 September 분석 스크립트의 `REPO`/`PROJECT`/`SOURCE`가 예입니다. CSV 입력 상당수는 이미 저장소의 `results/`에 있으므로, 재실행할 때 먼저 해당 파일과 이를 만든 스크립트를 찾고 경로를 연결해야 합니다. **클론만으로 전체 결과를 재생성할 수는 없습니다.** 후속 Codex가 이 경로 조사와 수정을 맡도록 [에이전트 지침](AGENTS.md)에 명시했습니다.

원본을 준비한 뒤 July의 핵심 분석 순서는 다음과 같습니다. 각 단계의 산출물은 이름이 대응하는 `results/260927_*/` 폴더에 저장됩니다.

```bash
python p0_charge_physical_time_lsa.py
python p0_charge_subset_fit_validation.py
python p0_charge_postfit_identifiability_audit.py
```

첫 명령은 고정 physical-time grid의 민감도·식별도를 검사합니다. 둘째는 통과한 subset을 charge로 피팅하고 discharge로 검증하며 첫 명령의 `subset_identifiability.csv`를 읽습니다. 마지막은 피팅 후 Jacobian을 재검사합니다. 계산 비용이 큰 DFN 최적화이므로 결과를 확인하려는 경우에는 먼저 [저장된 최종 보고서](results/260927_charge_physical_time_final/README_KO.md)와 [피팅 설계](results/260927_charge_fitting_design/README_KO.md)를 읽는 편이 좋습니다.

## 검증 파일과 결과 해석

이 저장소에는 독립된 `tests/` 또는 `pytest` 테스트 모음이 없습니다. `*_validation.py`, `*_audit.py`, `*_recheck.py`, `*_qc.py` 및 해당 노트북은 실험 데이터에 대한 분석·검증 실행 파일입니다. 각 결과 폴더의 `README_KO.md`, CSV, JSON이 판단 근거입니다. 대표적으로:

| 검토 내용 | 실행 파일 | 결과 |
| --- | --- | --- |
| 음극 GITT 파싱·재현성·평형 QC | [`gitt_ocp_qc.py`](gitt_ocp_qc.py) | [`260918_gitt_ocp_qc`](results/260918_gitt_ocp_qc/) |
| July charge 피팅과 discharge 검증 | [`p0_charge_subset_fit_validation.py`](p0_charge_subset_fit_validation.py) | [`260927_charge_subset_fit_validation`](results/260927_charge_subset_fit_validation/) |
| July 피팅 후 식별도 | [`p0_charge_postfit_identifiability_audit.py`](p0_charge_postfit_identifiability_audit.py) | [`260927_charge_postfit_identifiability_audit`](results/260927_charge_postfit_identifiability_audit/) |
| September 권장 설정 검증 | [`september_recommended_dsn_kn_fit.py`](september_recommended_dsn_kn_fit.py) | [`260927_september_recommended_dsn_kn_fit`](results/260927_september_recommended_dsn_kn_fit/) |
| 7월 최종값과 Ai2020 nominal 비교 | [`july_one_at_a_time_nominal_reversion.py`](july_one_at_a_time_nominal_reversion.py) | [`260928_july_nominal_reversion_audit`](results/260928_july_nominal_reversion_audit/) |

`results/260927_p0_charge_fit_discharge_validation/`의 5개 파라미터 동시 피팅은 중단됐습니다. 해당 폴더의 optimizer history를 최종 결과로 사용하지 마세요. [상태 기록](results/260927_p0_charge_fit_discharge_validation/STATUS_KO.md)을 확인하세요.

## 참고 문서

- [현재 인수인계와 결정](260927_HANDOFF_CURRENT_STATE_KO.md)
- [최종 physical-time 피팅 결과](results/260927_charge_physical_time_final/README_KO.md)
- [September BoL 별도 검증](results/260927_september_recommended_dsn_kn_fit/README_KO.md)
- [날짜별 파일·실행 지도](docs/FILE_MAP.md)

측정 데이터와 연구 결과를 변경할 때는 기존 `results/`를 덮어쓰기 전에 사용 입력과 설정을 기록하고 새 결과 폴더에서 비교하세요.
