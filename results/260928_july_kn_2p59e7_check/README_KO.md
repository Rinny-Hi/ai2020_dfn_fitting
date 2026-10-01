# 7월 모델 `kn=2.59e-7` 대조 검증

## 결론

7월 최종 모델에서 `kn`만 `2.59e-7`로 변경하면 charge와 discharge 모두 실험
physical-time 끝 전에 cutoff되어 공통 모델로 채택할 수 없다.

| Scenario | Direction | Full-time feasible | RMSE | Overlap RMSE | Capacity RMSE |
|---|---|---|---:|---:|---:|
| July final | Charge | yes | 54.00 mV | 54.00 mV | 3.11% |
| July final | Discharge | yes | 52.30 mV | 52.30 mV | 1.09% |
| `kn=2.59e-7` | Charge | no | invalid | 56.17 mV | 9.49% |
| `kn=2.59e-7` | Discharge | no | invalid | 45.28 mV | 0.62% |

방전 overlap 오차와 capacity 오차는 낮아지지만, charge는 최소 `7.38 min`,
discharge는 최소 `0.19 min` 일찍 종료된다. 따라서 겹치는 구간의 낮은 RMSE만으로
개선이라고 판단하면 안 된다.

## 파일

- `july_kn_2p59e7_comparison.png`: rate별 실험/기준/변경 모델 비교
- `kn_2p59e7_summary.csv`: 방향별 요약
- `kn_2p59e7_detail.csv`: 셀·rate별 상세
- `analysis_manifest.json`: 실행 조건과 파라미터
