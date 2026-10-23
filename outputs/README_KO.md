# 2026-10-21까지의 DFN 결과 스냅샷

> **최신 결정:** 이 문서는 과거 결과 기록이다. 2026-10-23 이후 후속 충전 프로토콜 최적화에는 [최종 A모델](261023_final_a_model/README_KO.md)을 사용한다. 아래 ORegan2022 조합은 최종 기준이 아니다.

이 폴더는 2026-10-21까지 생성한 간결본, 현대차 미팅 PPT, 그림, 수치 CSV 및 검토 보고서를 보존한다. 숫자를 읽을 때는 **최종 권고 조합**과 이후 발표 자료용 **Ai2020 전해질 중간 단계**를 구분해야 한다.

## 최종 권고 조합

- Ai2020 음극 OCP + 측정 양극 GITT OCP + ORegan2022 전해질.
- GITT `D_sn=4.170774×10⁻¹⁴`, `D_sp=6.888175×10⁻¹⁴ m²/s`; 최신 양극 `k_p=4.768215585×10⁻⁷` 고정.
- `b_n=2.5` 고정 후 0.5/1/2C **충전**으로 `k_n=9.7028×10⁻⁷` 피팅. 방전 3곡선은 독립 검증.
- 6곡선 전구간 RMSE 37.14 mV, 방전 검증 RMSE 31.74 mV, 2C 방전 RMSE 21.41 mV. 1C 충전 계산 최저 음극 전위는 −3.33 mV로, SoC 약 79% 이후 짧은 구간에서 0 mV를 하회한다.
- 판단 근거: [최종 선정 보고서](ORegan2022_Ai2020_최종선정_및_SoC음극전위.md), [최신 방법·오차 검토](ORegan2022_Ai2020_최종결과_방법검토.md), [조건별 오차 CSV](ORegan2022_Ai2020_최종_3단계_조건별_RMSE_MAE.csv), [충방전 그래프](ORegan2022_Ai2020_최종_충방전6곡선.png), [음극 전위 그래프](ORegan2022_Ai2020_최종_SoC별음극전위.png).
- 열람용 간결본: [양방향 검증 간결본](261021_Ai2020_DFN_ORegan_양방향검증_간결본.ipynb), [전체 단계 간결본](261021_Ai2020_DFN_최종_간결본.ipynb). 발표용: [최신 검증 PPT](261021_미팅자료_Ai2020_ORegan_최신검증본.pptx).

## 후속 발표 자료용 Ai2020 전해질 중간 단계

- Stage 1의 De/입자반경 → OCP/window → 측정 양극 OCP 변경은 [영문 6곡선 그래프](Stage1_dynamic_voltage_0p5C_1C_2C_English.png), [전체 오차표](Stage1_dynamic_voltage_full_error_table_English.png), [중앙 오차표](Stage1_dynamic_voltage_central_error_table_English.png) 및 [수치 CSV](Stage1_De_RpRn_OCP_window_양극OCP_3단계_동적전압_RMSE_MAE.csv)에 있다. 이 Stage 1의 “중앙”은 **전달 용량 10–70%** 기준이다.
- Stage 2에서 사용자가 지정한 EIS 초기 `k_n=5.54×10⁻⁷`, `b_n=2.914` 고정, 최신 `k_p=4.768215585×10⁻⁷`으로 시작하여 **`k_n`만** 충전 피팅했다. 피팅값은 `6.977111568×10⁻⁷`이다. `b_n`은 피팅하지 않았다.
- 해당 중간 단계의 충전 전구간 RMSE는 43.81→42.48 mV, 독립 방전은 30.41→31.03 mV, 6곡선 중앙 RMSE는 16.03→19.81 mV이다. 여기서 중앙은 전달 용량 10–70%의 고정 시간 구간이다. 충전 점수만으로 개선을 주장하지 않으며, 이 결과는 위 최종 권고 조합을 대체하지 않는다.
- [단계별 영문 표](Stage2_intermediate_kn554_to_charge_fit_table_English.png), [수치 CSV](Stage2_intermediate_kn554_to_charge_fit_summary_English.csv), [초기/피팅 음극 전위 그래프](Stage2_kn554_initial_vs_kn_only_fit_anode_potential_SoC.png), [전위 요약 CSV](Stage2_kn554_initial_vs_kn_only_fit_anode_potential_summary.csv).

## 재현 범위

그림/표 생성 및 DFN 실행에 사용한 스크립트는 [`../work/`](../work/)에 있다. 저장된 CSV 기반 영문 그림·표는 `work/make_stage1_dynamic_english.py`, `work/make_stage2_user_kn554_fit_tables.py` 등으로 재생성할 수 있다. DFN 전 구간을 처음부터 다시 계산하려면 일부 외부 원본 Excel이 필요하다. 특히 `work/stage2_setup_probe.py`는 `work/stage2_inputs/`의 저율 데이터, Cell 6-6 동적 데이터, 양극 GITT 7-7/7-8 원본을 읽으며 이 큰 원본들은 이 Git 스냅샷에 포함하지 않았다. 결과 CSV/그림은 당시 실행 결과이지, 이 클론만으로 재계산을 완료했다는 주장이 아니다.

최신 미팅 PPT는 기존 슬라이드를 기반으로 만든 한 개 버전만 올렸다. PPT 그림/표와 계산 근거를 확인할 때는 같은 이름의 CSV와 위 보고서를 함께 참조한다.
