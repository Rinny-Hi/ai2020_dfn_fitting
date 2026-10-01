# 7월 최종 모델 시각화

## 목적

선택된 7월 보수 모델 `Dsn+kn`의 fitting 전후 charge 및 untouched discharge
validation을, 셀 반복 곡선을 중복 표시하지 않고 각 rate에서 전압 RMSE가 가장 작은
실험 셀 하나와 비교했다.

## 적용 모델

- OCP: Ai2020 nominal graphite equilibrium + 기존 측정 cathode equilibrium
- `Dsn=7.466200474e-14 m2/s`
- `Dsp=6.888175324e-14 m2/s`
- `kn=1.450235740e-6`
- `kp=5.899519907e-7`
- `brugg_n/p/s=2.914/1.83/1.5`
- capacity는 fitting 목적함수에 포함하지 않고 사후 지표로만 평가

## 파일

- `july_prefit_postfit_best_cell.png`: rate별 최저오차 실험 셀과 fitting 전후 비교
- `july_final_charge_fit_discharge_validation.png`: 최종 charge/discharge 검증 곡선
- `july_final_rate_metrics.png`: rate별 전압 및 capacity 지표
- `july_experiment_cell_selection.csv`: 대표 실험 셀 선택 근거
- `july_prefit_postfit_best_cell_metrics.csv`: fitting 전후 수치
