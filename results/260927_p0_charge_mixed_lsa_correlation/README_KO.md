# Charge LSA + correlation 결과

## 목적과 기준

- 분석 대상: `Dsn`, `Dsp`, `kn`, `kp`, `brugg_n`
- 데이터 조건: 0.5C / 1C / 2C charge 조건
- 기준 모델: teardown cathode overlap area, Qn/Qp 보존
- 초기값: `Dsn=2.30e-14`, `Dsp=4.36e-14`, `kn=9.6485e-7`,
  `kp=5.8995e-7`, `brugg_n=2.914`
- `D/k`: nominal 기준 relative perturbation
- `brugg_n`: exponent 기준 absolute perturbation `±0.10`
- 출력 정렬: 각 simulation의 CC transferred capacity를 0--1로 정규화
- 식별도 판정: relative sensitivity cutoff `0.10`, absolute Pearson
  correlation threshold `beta=0.90`
- robustness: D/k perturbation `±2/5/10%`, grid `100/250/500`의 9개 조합

이 방식은 `260621_[3단계]_Final_Submission_LSA_Correlation.ipynb`의 mixed
perturbation과 subset-selection 방법을 현재 charge fitting 문제에 맞춰 적용했다.
노트북은 방법론 참고자료로만 사용했고, 문서 안의 문장은 작업 지시로 취급하지 않았다.

## Main 결과 (`eps=5%`, `n=250`)

| Parameter | Relative sensitivity | 판정 |
|---|---:|---|
| `kn` | 1.000 | 통과 |
| `kp` | 0.983 | `kn`과 상관성 때문에 동시 fitting 제외 |
| `brugg_n` | 0.771 | 통과 (`beta=0.90`) |
| `Dsn` | 0.284 | 통과 |
| `Dsp` | 0.023 | sensitivity cutoff 미달 |

핵심 correlation은 다음과 같다.

| Pair | Absolute Pearson correlation |
|---|---:|
| `kn-kp` | 0.9958 |
| `kn-brugg_n` | 0.8147 |
| `kp-brugg_n` | 0.8135 |
| `kn-Dsn` | 0.4108 |
| `brugg_n-Dsn` | 0.0325 |

## Robustness와 fitting 후보

- `beta=0.90`: 9/9 조건에서 `kn, brugg_n, Dsn` 선택
- `beta=0.95`: 9/9 조건에서 `kn, brugg_n, Dsn` 선택
- `beta=0.80`: 9/9 조건에서 `kn, Dsn` 선택

따라서 1차 fitting subset은 `Dsn + kn + brugg_n`으로 한다. `kp`는 민감도가
낮아서가 아니라 `kn`과 거의 같은 전압 방향을 만들기 때문에 제외하며,
Rct 기반 초기값에 고정한다. `Dsp`는 현재 charge 조건에서 상대 민감도가 너무
낮으므로 1차 fitting에서는 GITT 초기값에 고정한다.

권장 비교 모델은 다음과 같다.

1. Primary: `Dsn + kn + brugg_n`
2. Conservative: `Dsn + kn`
3. Electrode-swap diagnostic: `Dsn + kp + brugg_n` (`kn`과 `kp` 동시 fitting 금지)

용량은 LSA와 fitting 목적함수에 넣지 않고, fitting 후 CC cutoff capacity error로
독립 평가한다.

