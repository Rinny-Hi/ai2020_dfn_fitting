# 전극 질량 보정, Dsn bound 검토, physical-time fitting 최종 결과

## 결론

1. GITT 식에 최신 한 면 coating 질량 `negative=10.10 mg`, `positive=17.35 mg`을
   적용한다. 다만 한 개 펀칭 시편의 질량을 full-cell `Qn/Qp` 변경으로 전파하지
   않고, **GITT apparent D 계산에만** 적용한다.
2. mass-corrected all-range reference는
   `Dsn=4.170774146e-14`, `Dsp=6.888175324e-14 m2/s`이다.
3. graphite 문헌과 현재 GITT endpoint envelope를 함께 보면 1차
   `Dsn=1e-14--8e-14 m2/s`는 유지할 수 있다. 이는 신뢰구간이 아니라 상수
   effective-D 탐색범위다.
4. physical-time 민감도/식별도에서는 주 모델과 보수 모델이 모두 통과했다.
5. charge fitting과 untouched discharge validation 후에는 **보수 모델
   `Dsn+kn`을 선택**한다. 주 모델의 charge training 이득은 0.32 mV에 불과한데
   discharge full RMSE와 사후 capacity error가 크게 악화했다.
6. 선택 모델의 fitted `Dsn`은 optimizer별 분산이 남으므로 intrinsic graphite
   diffusivity 확정값으로 해석하지 않는다.

## 1. 이전 질량과 최신 질량 비교

프로젝트의 Neware/Weppner--Huggins 식에서 다른 항을 고정하면

\[
D_{s,new}=D_{s,old}\left(\frac{m_{new}}{m_{old}}\right)^2
\]

이다.

| 전극 | 이전 질량 | 최신 질량 | D 배율 | all-range 이전 D | all-range 최신 D |
|---|---:|---:|---:|---:|---:|
| Anode | 10.40 mg | 10.10 mg | 0.943140 | 4.422e-14 | **4.171e-14** |
| Cathode | 17.25 mg | 17.35 mg | 1.011628 | 6.809e-14 | **6.888e-14** |

방향별 최신 값은 anode charge/discharge `6.406/2.715e-14`, cathode
charge/discharge `4.611/10.290e-14 m2/s`이다. 음극 20% endpoint trim combined
최솟값도 `1.689e-14 -> 1.593e-14`로 내려간다.

두 값 모두 `(전극 원판-집전체)/2`의 한 면 **composite coating 질량**이며,
활물질 질량분율은 미확정이다. 실제 GITT 시편의 최신 계량값이라는 전제에서 최신
질량을 nominal로 쓰되, 이전 질량 결과를 약 `+5.7%/-1.2%`의 mass-choice
uncertainty로 보존한다. endpoint 선택과 active fraction 불확실성이 이 차이보다
더 크므로 질량 변경만으로 intrinsic D 정확도가 확보되는 것은 아니다.

## 2. graphite Dsn 문헌 및 선행연구 bound 비교

| 출처/방법 | graphite D 또는 탐색범위 | 이번 작업에서의 의미 |
|---|---:|---|
| Pham et al., commercial graphite GITT | `1e-15--4e-14 m2/s` | 직접 GITT 범위. 최신 combined 값 4.17e-14는 상단과 비슷함 |
| Chen et al. 2020 / PyBaMM Chen2020 | `3.3e-14 m2/s` | 상수 DFN reference |
| Ai2020/Marquis2019 PyBaMM implementation | `3.9e-14 m2/s` | 상수 graphite reference |
| Ecker2015 PyBaMM 함수, 이번 x-window/298 K 평가 | `8.98e-15--8.86e-13`, median `1.31e-14 m2/s` | SOC 함수는 국부 peak 때문에 단일 상수보다 훨씬 넓음 |
| Röder et al., MCMB graphite 정리 | 약 `1e-15--1e-13 m2/s` | graphite 종류·입자구조에 따른 넓은 물성범위 |
| Le Roux/Khalik grouped DFN bound | `Dhat_sn=1.3e-4--1.6e-3 s-1`; 현재 Rn으로 환산 시 `1.83e-15--2.26e-14 m2/s` | 다른 셀·grouped 정의이므로 직접 bound 복사는 부적절 |
| Zhao & Jossen time-domain identification | anode D posterior가 상한으로 접근하여 상한만으로는 고유값 식별 불가 | voltage fit에서 큰 D가 rate-limiting이 아니면 upper plateau가 생김 |

문헌 링크:

- Pham et al., *Investigation of Lithium Ion Diffusion of Graphite Anode by GITT*:
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC8397968/>
- Röder et al., *Simulating the Impact of Particle Size Distribution on Graphite Electrodes*:
  <https://doi.org/10.1002/ente.201600232>
- Chen2020 graphite diffusivity의 PyBaMM 구현과 원 논문 인용:
  <https://github.com/pybamm-team/PyBaMM/blob/main/packages/pybamm/src/pybamm/input/parameters/lithium_ion/OKane2022.py>
- Le Roux et al.의 grouped parameter range와 log/linear normalization:
  <https://doi.org/10.1016/j.ifacol.2023.10.714>
- Khalik et al., teardown과 voltage fitting 결합 필요성:
  <https://doi.org/10.1016/j.jpowsour.2021.229901>
- Zhao & Jossen, time/frequency-domain identifiability 비교:
  <https://doi.org/10.3390/batteries8110222>
- Forman et al., L2 voltage objective와 Fisher identifiability:
  <https://doi.org/10.1016/j.jpowsour.2012.03.009>

### Bound 결정

- 1차 `1e-14--8e-14`는 최신 mass-corrected endpoint envelope
  `1.593e-14--6.406e-14`를 포함한다.
- 문헌의 더 낮은 값은 SOC-resolved/local diffusion에서 가능하지만, 이번 fitting은
  0.5C/1C/2C 전 구간의 **상수 effective D**를 찾으므로 측정 envelope를 우선한다.
- 상한 `8e-14`는 direct GITT와 상수 DFN reference보다 넓고, graphite/MCMB의
  `~1e-13` 범위와도 인접한다.
- 설계대로 `1.2e-13` 확장 test를 수행했다. 주/보수 모델 RMSE 개선은 각각
  `0.0045/0.0032 mV`뿐이고 Dsn은 `7.682/7.466e-14`에 머물렀다. 따라서 1차
  상한이 local optimum을 절단했다는 증거는 없다.

## 3. physical-time 민감도/식별도 재검증

July BoL 3셀의 charge 0.5C/1C/2C 각 branch를 100개 physical-time point로
재표본화했다. CC command의 95%에 처음 도달한 record를 `t=0`, 실험 CC cutoff를
끝점으로 사용했으며 capacity는 Jacobian에 넣지 않았다.

최신 GITT baseline은 일부 실험 cutoff보다 최대 1.214 min 먼저 종료하여 full-grid
Jacobian을 정의할 수 없었다. 따라서 `Dsn=4.171e-14`, `brugg_n=2.914`를 유지하고
`kn=1.90e-6`으로 올린 가장 가까운 단순 feasible anchor에서 1/2.5/5% normalized
finite difference를 계산했다.

| 모델 | relative sensitivity | max |cosine| | condition | 판정 |
|---|---|---:|---:|---|
| `Dsn+kn+brugg_n` | `0.126/0.670/1.000` | 0.891 (`kn-brugg_n`) | 5.12 | 통과 |
| `Dsn+kn` | `0.187/1.000` | 0.507 | 1.75 | 통과 |

세 perturbation 크기에서 같은 판정이었다. 통과 규칙은 full rank, 최소 상대민감도
0.10, 최대 pairwise cosine 0.95, normalized condition number 20 이하이다.

## 4. charge fitting 및 discharge validation

두 모델에 5-start TRF, DE->TRF, DA->TRF를 동일하게 적용했다. 조기 cutoff는
capacity objective가 아니라 fixed time-voltage residual을 정의할 수 없는 feasibility
위반으로 처리했다.

| 지표 | 주 모델 | 보수 모델 |
|---|---:|---:|
| best optimizer | DA -> TRF | TRF, all-range start |
| Dsn | 7.682e-14 | **7.466e-14** |
| kn | 8.926e-7 | **1.450e-6** |
| brugg_n | 2.524 | **2.914 고정** |
| charge objective RMSE | **53.676 mV** | 53.997 mV |
| charge full RMSE | **52.89 mV** | 53.69 mV |
| discharge full RMSE | 62.08 mV | **49.51 mV** |
| charge capacity RMSE, post hoc | 6.45% | **3.11%** |
| discharge capacity RMSE, post hoc | 1.51% | **1.09%** |

주 모델의 charge objective 이득은 `0.321 mV`뿐인 반면 discharge full RMSE는
`+12.58 mV`, charge capacity RMSE는 `+3.34%p` 악화했다. 따라서 설계의 선택
규칙에 따라 보수 모델을 선택한다.

최적점은 full-time feasibility 경계에 매우 가깝다. cutoff penalty를 미분한 잘못된
Jacobian을 피하기 위해 feasible 방향의 one-sided voltage derivative로 다시 검사했다.
주 모델 condition `6.15`, max cosine `0.937`; 보수 모델 condition `1.62`, max
cosine `0.451`로 두 모델 모두 local test는 통과했다. 그러나 보수 모델이 훨씬
안정적이다.

## 5. 최종 해석과 제한

- 선택 predictive model: `Dsn+kn`, `Dsn=7.466e-14`, `kn=1.450e-6`,
  `brugg_n=2.914` 고정.
- optimizer별 보수 모델 Dsn은 `5.74--7.47e-14`로 퍼졌고 objective 차이는 작다.
  따라서 fitted Dsn은 unique intrinsic property가 아니라 full-time voltage와 cutoff
  feasibility가 정한 effective parameter다.
- 최신 질량을 GITT 식에 쓰는 것과 full-cell solid inventory를 바꾸는 것은 별개다.
  3--5개 이상 펀칭 시편의 평균/표준편차와 active fraction이 확보되기 전에는 최신
  한 개 질량으로 `Qn/Qp`를 재설정하지 않는다.
- 다음 개선은 constant D bound를 더 넓히는 것보다 SOC-dependent D 또는
  finite-sphere pulse fitting, 그리고 current/EIS가 포함된 별도 식별 실험이 우선이다.

## 산출물

- mass rescaling: `../260927_gitt_mass_and_dsn_bound_comparison/`
- physical-time LSA: `../260927_charge_physical_time_lsa/`
- optimizer/validation: `../260927_charge_subset_fit_validation/`
- Dsn upper-bound test: `../260927_charge_dsn_boundary_test/`
- post-fit feasible-side Jacobian: `../260927_charge_postfit_identifiability_audit/`
