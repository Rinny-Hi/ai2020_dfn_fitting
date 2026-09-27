# 재측정 coating loading 반영 후 Lp/Ln 3종 비교

## 입력 해석

- 양극 `17.35 mg`, 음극 `10.10 mg`은 각각
  `(펀칭 전극 전체 질량 - 집전체 질량) / 2`로 얻은 **한 면 coating 질량**이다.
- 10 mm 펀칭 면적 `0.785398 cm2` 기준 loading은 양극 `22.09 mg/cm2`,
  음극 `12.86 mg/cm2`이다.
- 이전 측정 `17.25/10.40 mg` 대비 양극 inventory는 `+0.58%`, 음극은
  `-2.88%`로 반영했다.
- 활물질/도전재/바인더 질량분율은 없으므로 절대 porosity를 새로 단정하지 않고,
  기존 effective solid inventory에 재측정 질량비를 적용했다.

## 물리적으로 일관된 반영 순서

1. 재측정 질량비로 `Qn/Qp`와 solid inventory를 변경했다.
2. `delta_x=Qcell/Qn`, `delta_y=Qcell/Qp`로 window 폭을 고정했다.
3. C/50 qOCV에 대해 `x0`, `y100` 위치만 재정렬했으며 양 끝점 오차를
   각각 10 mV 이내로 제한했다.
4. 세 두께에서 같은 coating inventory를 사용하고 두께에 따라 porosity가
   체적수지를 닫도록 했다.
5. 현재 P0의 SEM `Rn/Rp`, 측정 `Dsn/Dsp`, EIS `kp`, Ai2020 `kn`,
   Bruggeman 계수와 hysteresis OFF를 공통 적용했다.

재식별 결과:

| 항목 | 기존 P0 | loading 반영 |
|---|---:|---:|
| Qn | 2.47947 Ah | 2.40795 Ah |
| Qp | 4.37120 Ah | 4.39654 Ah |
| delta x | 0.95042 | 0.97865 |
| delta y | 0.53910 | 0.53600 |
| x0 / x100 | 0.00381 / 0.95422 | 0.00374 / 0.98239 |
| y100 / y0 | 0.43164 / 0.97075 | 0.43703 / 0.97302 |
| qOCV RMSE | 6.24 mV | 7.27 mV |

## 두께별 결과

July BoL 3셀, 0.5C/1C/2C CC 충·방전 18개 곡선 평균이다.

| 두께 | Ln/Lp | full RMSE | 10-70% RMSE | 용량 RMSE | 평균 절대 용량오차 | 최대 용량오차 |
|---|---:|---:|---:|---:|---:|---:|
| Ai2020 | 76.5/68.0 um | 42.64 mV | 26.91 mV | **1.12%** | **13.46 mAh** | **42.00 mAh** |
| 두께측정기 | 73.0/62.5 um | **40.16 mV** | **18.95 mV** | 2.96% | 37.86 mAh | 121.37 mAh |
| SEM | 81.52/62.965 um | 46.31 mV | 31.23 mV | 2.20% | 20.32 mAh | 89.08 mAh |

### Ai2020 두께의 C-rate별 용량오차

| 조건 | 용량오차 |
|---|---:|
| 0.5C charge | -26.63 mAh (-1.24%) |
| 0.5C discharge | +4.61 mAh (+0.20%) |
| 1C charge | -7.24 mAh (-0.38%) |
| 1C discharge | +9.26 mAh (+0.41%) |
| 2C charge | +26.44 mAh (+1.81%) |
| 2C discharge | +0.25 mAh (+0.01%) |

## 기존 P0 대비 판단

| 모델 | 10-70% RMSE | 용량 RMSE | 평균 절대 용량오차 |
|---|---:|---:|---:|
| 기존 P0, Ai2020 두께 | **21.92 mV** | 1.32% | 18.71 mAh |
| loading+window 재식별, Ai2020 두께 | 26.91 mV | **1.12%** | **13.46 mAh** |

loading을 반영하면 용량 RMSE는 `1.32 -> 1.12%`, 평균 절대 용량오차는
`18.71 -> 13.46 mAh`로 개선되지만 중앙 전압 RMSE는 `21.92 -> 26.91 mV`로
약 5 mV 악화한다. 따라서 새 안은 용량 측면에서는 개선되지만 기존 P0를 전 지표에서
지배하지는 않는다.

## 선택

세 두께 중에서는 **Ai2020 `Ln/Lp=76.5/68 um`를 선택**한다.

- 두께측정기 안은 중앙 전압은 가장 좋지만 2C charge에서 `-90.75 mAh`, 전체
  최대 `121.37 mAh` 오차가 발생해 geometry 후보로 채택하기 어렵다.
- SEM 안은 Ai2020 대비 전압과 용량 모두 열세다.
- Ai2020 안은 18개 곡선 용량 RMSE가 가장 낮고 모든 C-rate에서 용량오차가
  약 +/-1.8% 이내로 유지된다.

최종적으로는 두 모델을 다음처럼 분리한다.

- **물리 loading 우선 후보:** loading+window 재식별, Ai2020 두께
- **동특성 fitting 시작점:** 기존 P0와 위 후보를 모두 남기고, 동일한 contact
  resistance/`kn` fitting 후 방전 validation 결과로 최종 선택

loading 기반 후보의 qOCV endpoint가 허용 한계 `+/-10 mV`에 닿고 있으므로,
활물질 질량분율 또는 porosity 없이 이를 유일한 확정안으로 강제하지 않는다.

## 산출물

- `capacity_voltage.png`: 세 두께의 C-rate별 용량-전압 곡선
- `qocv_window_error.png`: 기존/재식별 window의 qOCV 오차
- `summary.csv`: 전체 지표
- `branch_summary.csv`: C-rate/방향별 지표
- `detail.csv`: 셀별 상세 지표
- `analysis.json`: window와 계산 조건

