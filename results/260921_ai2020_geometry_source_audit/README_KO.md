# Ai2020 geometry 출처 및 용량 영향 검토

## 1. 세 논문의 관계

1. `Multi-scale investigation of thickness changes in a commercial pouch type lithium-ion battery`가 Enertech SPB655060을 실제 해체해 geometry를 얻은 원 실험 논문이다.
2. `A New Method to Model the Thickness Change of a Commercial Pouch Cell during Discharge`는 위 논문을 Ref. 15로 인용해 geometry를 P2D 모델에 적용했다.
3. `Electrochemical Thermal-Mechanical Modelling of Stress Inhomogeneity in Lithium-Ion Pouch Cells`는 Ai et al. (2020) 논문이며, 앞선 Rieger 논문의 geometry와 전기화학 파라미터를 다시 사용했다. PyBaMM `Ai2020` parameter set의 직접적인 이름 출처다.

즉 Ai2020 geometry는 fitting 결과가 아니라 동일 모델 셀의 teardown 측정값을 전승한 것이다.

## 2. 원 실험 논문에서 geometry를 얻은 방법

셀을 해체하여 내부 구조를 직접 확인했다.

| 항목 | 원 논문 값 | 방법/성격 |
|---|---:|---|
| anode sheet | 양면 17장 | 해체 후 직접 계수 |
| cathode sheet | 양면 16장 + 단면 2장 | 해체 후 직접 계수 |
| electrochemically active coated faces | 양·음극 각각 34면 | 위 sheet 수에서 계산 |
| anode coating thickness | `77 +/- 0.5 um` | measured, 구체적인 두께 측정 장비는 논문에 미기재 |
| cathode coating thickness | `68 +/- 0.5 um` | measured, 구체적인 두께 측정 장비는 논문에 미기재 |
| separator thickness | `25 +/- 0.5 um` | measured |
| Cu current collector | `10 +/- 0.5 um` | measured |
| Al current collector | `15 +/- 0.5 um` | measured |
| anode footprint | `46.4 x 52 mm` | measured |
| cathode footprint | `45.4 x 51 mm` | measured |
| cathode porosity | `33 +/- 0.5%` | mercury porosimetry |
| anode porosity | `32 +/- 1%` | mercury porosimetry |
| cathode active fraction | `62 +/- 0.5%` | estimated, inactive material 5 vol% 가정 |
| anode active fraction | `61 +/- 1%` | estimated, inactive material 7 vol% 가정 |

원 논문의 active-material fraction은 독립적인 영상분할/밀도 측정값이 아니다.

```text
positive eps_s = 1 - porosity(0.33) - inactive(0.05) = 0.62
negative eps_s = 1 - porosity(0.32) - inactive(0.07) = 0.61
```

따라서 `eps_s=0.62/0.61`의 신뢰도는 mercury porosimetry뿐 아니라 비활물질 체적분율 가정에도 의존한다.

## 3. 논문 사이의 porosity 불일치

원 teardown 논문:

```text
positive porosity = 0.33
negative porosity = 0.32
```

후속 Rieger P2D 논문과 Ai2020 parameter set:

```text
positive porosity = 0.32
negative porosity = 0.33
```

양·음극 값이 서로 뒤바뀌어 있다. 반면 active fractions `0.62/0.61`은 원 논문 값을 유지했다. 후속 논문은 각 상의 합을 양·음극 모두 `0.94`로 만들고 inactive fraction을 각각 `0.06`처럼 취급하는 형태가 된다.

이 0.01 차이는 `eps_s`를 고정하면 이론용량에 직접 들어가지 않지만, `D_e,eff=D_e*eps_e^b`와 `kappa_eff=kappa*eps_e^b`를 통해 고율 수송에는 영향을 준다. 현재처럼 `b_n=2.914`, `b_p=1.83`이면 원 논문 porosity로 되돌릴 때 유효수송 변화는 대략:

```text
negative: (0.32/0.33)^2.914 = 0.914  -> 약 -8.6%
positive: (0.33/0.32)^1.83  = 1.058  -> 약 +5.8%
```

따라서 원 논문 porosity 복원은 수 퍼센트 이상의 고율 전압/용량 변화를 만들 수 있으나, `eps_s`를 같이 바꾸지 않으면 저율 이론용량은 그대로다.

## 4. 전극 면적과 layer 수

원 teardown 치수를 그대로 사용한 총 단면 active area:

```text
A_n = 46.4 mm * 52 mm * 34 = 820.352 cm2
A_p = 45.4 mm * 51 mm * 34 = 787.236 cm2
```

PyBaMM Ai2020은 양·음극 공통 면적으로 다음 값을 사용한다.

```text
A_common = 47 mm * 51 mm * 34 = 814.980 cm2
```

따라서 공통 면적은 원 geometry 대비:

- anode 면적보다 `0.65%` 작음
- cathode 면적보다 `3.52%` 큼

P2D가 동일한 through-plane 면적을 요구하므로 단순화한 것으로 보이지만, 전극 inventory 계산에는 cathode 기준 약 3.5% 차이를 만든다.

## 5. Ai2020 값으로 계산한 Q

```text
Q_k = F*A_k*L_k*eps_s,k*c_s,max,k/3600
Q_cell = Q_n*delta_x = Q_p*delta_y
```

Ai2020 nominal 값:

- `c_s,max,n=28,700 mol/m3`
- `c_s,max,p=49,943 mol/m3`
- `delta_x=0.84-0.0065=0.8335`
- `delta_y=0.9651-0.435=0.5301`

| 면적 정의 | `Qn` full inventory | `Qn*delta_x` | `Qp` full inventory | `Qp*delta_y` |
|---|---:|---:|---:|---:|
| Ai2020 공통 면적 814.98 cm2 | 2.9254 Ah | 2.4383 Ah | 4.5992 Ah | 2.4380 Ah |
| 원 anode 면적 820.35 cm2 | 2.9446 Ah | 2.4544 Ah | - | - |
| 원 cathode 면적 787.24 cm2 | - | - | 4.4426 Ah | 2.3550 Ah |

공통 면적과 endpoint를 쓰면 양·음극 accessible capacity가 약 `2.438 Ah`로 거의 정확하게 맞춰져 있다. 반면 실제 전극별 footprint를 사용하면 cathode accessible capacity가 `2.355 Ah`로 제한된다.

현재 Enertech C/50 기준 용량 `2.35653 Ah`는 원 cathode geometry로 계산한 `2.35503 Ah`와 약 `1.50 mAh (0.064%)` 차이다. 이는 다음 조합이 현재 저율용량과 매우 일관됨을 의미한다.

```text
measured cathode area 45.4*51 mm
34 coated faces
Lp=68 um
eps_s,p=0.62
c_s,max,p=49,943 mol/m3
Ai2020 cathode window delta_y=0.5301
```

완전히 독립적인 검증은 아니지만, 적어도 양극 geometry와 nominal `c_s,max,p`를 크게 변경할 근거는 약하다. 반대로 anode는 같은 계산에서 `2.454 Ah`로 약 4.1%의 여유를 가지며 commercial cell의 unused graphite capacity 설명과 일치한다.

## 6. 논문에서 electrode balancing을 처리한 방법

원 논문은 half-cell graphite와 full-cell differential-potential peak를 맞추기 위해 graphite 데이터를 면적용량 축에서 `-0.32 mAh/cm2` 이동했다. 상용 셀에서 Li plating을 피하기 위해 graphite가 완전히 lithiation되지 않는 unused capacity를 제거하기 위한 조치다.

즉 논문도 geometry만으로 full-cell capacity를 결정하지 않았다.

1. teardown geometry와 porosity로 physical scale을 설정
2. half-cell 전극 곡선 측정
3. full-cell differential potential과 graphite peak를 정렬
4. 사용되지 않는 graphite capacity를 window에서 제외

현재의 `Qn/Qp-window` 동시 재식별과 기본 논리는 같다.

## 7. 두께 상태와 측정법의 영향

원 논문은 셀을 15 cycle pre-cycling한 뒤 사용했고, dilatometry용 전극은 2.8 V 셀에서 회수해 DMC 10분 세척, vacuum 60분 건조 후 측정했다. Dilatometer에서는 working electrode에 1 N 하중을 적용했다.

또한 operating window 내 graphite coating의 reversible thickness change를 약 `5.2%`, 약 `4 um`로 보고했다. 따라서 teardown 시점 SOC, 세척·건조, 하중 유무에 따라 수 um 차이는 가능하다. SEM이 무하중·국부 단면을 측정하고 두께측정기가 넓은 면적을 압축한다면 SEM 값이 더 클 수도 있다.

그러나 이 효과는 현재 음극의 `73-76.5 um`과 질량 역산 `~115 um` 사이의 40 um 전후 차이를 단독으로 설명하기에는 부족하다.

## 8. 현재 Enertech 모델에 대한 결론

1. `Lp=68 um`은 동일 셀 teardown measured value이며, 원 cathode footprint까지 적용하면 현재 C/50 용량과 거의 정확히 일치한다. 우선 유지할 근거가 강하다.
2. `Ln=76.5 um`도 동일 셀 measured value `77 +/- 0.5 um`에서 온 값이다. SEM이 더 크게 나올 수 있지만 현재 음극 질량 역산값 `115 um`를 바로 채택할 근거는 없다.
3. `eps_s,p=0.62`, `eps_s,n=0.61`은 직접 측정이 아니라 porosity와 inactive-volume 가정의 잔차다. coating loading과 조성이 확보되면 재계산해야 한다.
4. porosity는 원 논문과 Ai2020에서 양·음극이 뒤바뀌어 있다. 이후 ablation에서 `source porosity (p=0.33,n=0.32)`와 `Ai2020 porosity (p=0.32,n=0.33)`를 비교할 가치가 있다.
5. 총 active area는 양·음극 공통 면적을 쓰는 현재 구현보다 원 전극별 footprint가 용량 해석에 더 적합하다. 특히 cathode area `787.24 cm2`를 적용하면 nominal 양극 inventory가 실제 저율용량과 일치한다.
6. 내일 loading을 받으면 `A`, `L`, `eps_s`, `c_s,max`를 각각 독립 fitting하지 않고 `m_AM,areal=rho_AM*L*eps_s` 질량수지로 묶어야 한다.

집전체 두께는 이후 기준 구성에서 원 논문값 `Cu=10 um`, `Al=15 um`를 사용한다. 사용자 측정 `9/13 um`는 측정법 차이를 확인하기 위한 sensitivity 기록으로만 보존한다.

부피분율의 측정/추정 관계는 다음과 같다.

- electrolyte volume fraction `eps_e` = porosity: mercury porosimetry로 **측정** (`positive=0.33`, `negative=0.32`)
- inactive solid fraction: 직접 정량값이 아니라 **가정** (`positive=0.05`, `negative=0.07`)
- active-material volume fraction `eps_s`: `1-eps_e-eps_inactive`로 **추정/계산** (`positive=0.62`, `negative=0.61`)

따라서 `eps_s`가 측정값이고 `eps_e`가 추정값인 것은 아니다. 원 논문에서는 반대로 `eps_e`가 측정값이고 `eps_s`가 가정에 의존한 계산값이다.

## 9. 이후 적용할 면적 정책

Teardown 원 치수는 다음과 같이 그대로 기록한다.

```text
negative total coated area = 46.4 mm * 52 mm * 34 = 820.352 cm2
positive total coated area = 45.4 mm * 51 mm * 34 = 787.236 cm2
```

다만 표준 1D DFN은 through-plane 전류가 통과하는 공통 단면적 하나만 가질 수 있다. 음극 overhang까지 양극과 다른 전류면적으로 넣으면 양·음극 전류보존이 깨진다. 따라서 앞으로 다음처럼 구분한다.

- dry mass 및 전체 전극 inventory bookkeeping: `A_n=820.352 cm2`, `A_p=787.236 cm2`
- 1D DFN current/reaction area: 작은 cathode footprint인 `A_overlap=A_p=787.236 cm2`
- 음극 overhang이 고율에서도 모두 즉시 사용 가능하다고 가정하지 않음

이 정책은 `paper_teardown_geometry.py`에 재사용 가능한 형태로 저장했다. 기존 결과와의 비교 가능성을 위해 현재 baseline을 즉시 덮어쓰지 않고, 다음 OCP 재식별부터 paper-area candidate에 적용한다.
