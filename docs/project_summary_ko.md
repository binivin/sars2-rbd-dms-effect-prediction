# RBD DMS 효과 예측 및 일반화 한계 분석

## 프로젝트 개요

이 프로젝트는 RBD deep mutational scanning(DMS) 데이터를 이용해 아미노산 치환이 receptor binding 및 expression/folding score에 미치는 영향을 예측하는 생물정보학 분석이다.

핵심 목표는 단순히 높은 예측 성능을 얻는 것이 아니라, 모델이 실제로 어디까지 일반화되는지를 검증하는 것이다. 이를 위해 mutation-level feature, 구조 feature, reference ESM-2 embedding, mutant-sequence delta ESM feature를 단계적으로 추가하고 여러 validation setting을 비교하였다.

## 연구 질문

이 프로젝트의 핵심 질문은 다음과 같다.

```text
RBD 아미노산 치환의 DMS score를 sequence, structure, protein language model feature로 예측할 수 있는가?

그리고 모델은 이미 관측된 residue position 안에서만 잘 보간하는가,
아니면 완전히 새로운 residue position까지 일반화할 수 있는가?
```

## 데이터셋

분석용 데이터셋은 public RBD DMS 자료를 통합해 구축하였다.

```text
22,572 amino-acid substitutions
5 measured backgrounds: Wuhan-Hu-1, Delta, Beta, E484K, N501Y
2 prediction targets: receptor-binding score, expression/folding score
```

Raw score는 source/background마다 scale이 다르기 때문에, 각 source/background group 안에서 z-score normalization을 수행하였다.

## Feature 구성

사용한 feature는 네 단계로 확장되었다.

1. Mutation-level physicochemical feature
   - residue position
   - wild-type / mutant amino-acid identity
   - hydropathy, volume, charge 변화
   - proline/cysteine indicator

2. Structure-derived feature
   - 6M0J 구조 기반 receptor-binding-motif annotation
   - receptor chain과의 최소 거리
   - 5Å interface 여부
   - 접촉 residue 수
   - 구조상 해당 site 존재 여부

3. Reference ESM-2 residue embedding
   - `facebook/esm2_t6_8M_UR50D` 모델을 이용한 residue-level sequence-context embedding

4. Mutant-sequence delta ESM feature
   - mutant sequence와 background wild-type sequence의 ESM representation 차이

```text
delta_embedding = mutant_sequence_embedding - background_wildtype_sequence_embedding
```

## Validation 전략

이 프로젝트의 가장 중요한 부분은 validation setting을 세분화한 것이다.

### Random split

가장 쉬운 검증이다. 같은 site 또는 유사한 mutation pattern이 train/test에 함께 포함될 수 있으므로, 높은 성능이 실제 일반화를 의미한다고 보기 어렵다.

### Background-held-out

하나의 measured background를 통째로 test set으로 제외한다. 모델이 다른 background로 mutation effect pattern을 이전할 수 있는지 확인한다.

### Mutation-held-out

특정 substitution identity를 모든 background에서 제외한다. 같은 residue position의 다른 substitution은 train에 남아 있을 수 있으므로, known site 안에서 새로운 mutation을 예측하는 능력을 평가한다.

### Site-held-out

특정 residue position 전체를 train에서 제외한다. 완전히 처음 보는 site의 mutation effect를 예측해야 하므로 가장 어려운 검증이다.

## 주요 결과

### Background-held-out 결과

Background-held-out에서는 ESM feature를 포함한 RandomForest 모델이 높은 성능을 보였다.

```text
Binding R2_mean    about 0.91
Expression R2_mean about 0.94
```

이는 측정된 background 사이에서 mutation effect pattern이 어느 정도 이전될 수 있음을 의미한다.

### Mutation-held-out 결과

Mutation-held-out에서도 reference ESM + structure + background feature를 포함한 RandomForest 모델이 높은 성능을 보였다.

```text
Binding R2_mean    about 0.81
Expression R2_mean about 0.81
```

이는 모델이 이미 관측된 residue position 안에서 새로운 substitution effect를 잘 예측할 수 있음을 의미한다.

### Site-held-out 결과

Site-held-out은 가장 어려운 검증이었다. Delta ESM + structure + background feature가 가장 좋은 결과를 냈지만, 성능은 여전히 낮았다.

```text
Binding R2_mean    about 0.03
Expression R2_mean about 0.12
```

즉, delta ESM은 완전히 새로운 site 예측에 약간 도움을 주지만, unseen-site prediction 문제를 해결하기에는 부족했다.

## Feature importance 해석

Grouped permutation importance 결과, mutation-held-out validation에서는 reference ESM feature가 가장 중요한 feature group이었다.

Binding prediction에서는 reference ESM을 섞었을 때 R2 drop이 가장 컸고, expression/folding prediction에서도 같은 경향이 나타났다.

이는 모델이 단순한 amino-acid property보다 ESM-2가 담고 있는 residue-level sequence context를 강하게 활용했음을 의미한다.

반면 site-held-out에서는 reference ESM과 delta ESM이 일부 정보를 제공했지만, 전체 성능 자체가 낮았기 때문에 완전한 일반화에는 한계가 있었다.

## Difficult-site 분석

Site-held-out error는 RBD 전체에 균일하게 분포하지 않았다.

Combined error 기준으로 어려운 site는 다음과 같았다.

```text
454, 442, 355, 398, 379, 467, 490, 350, 423, 461
```

Binding prediction에서는 receptor-binding-motif 관련 위치에서 큰 error가 많이 나타났다.

```text
490, 504, 442, 454, 487, 483, 499
```

Expression/folding prediction에서는 receptor-binding-motif뿐 아니라 non-receptor-binding-motif 위치에서도 어려운 site가 나타났다.

```text
379, 355, 461, 398, 467, 454, 380, 430, 492, 417
```

이는 binding effect와 expression/folding effect가 서로 다른 생물학적 민감 영역에 의해 영향을 받을 수 있음을 시사한다.

## 최종 해석

이 프로젝트의 핵심 결론은 다음과 같다.

```text
모델은 이미 관측된 residue position 안에서의 새로운 substitution은 비교적 잘 예측하지만,
완전히 처음 보는 residue position의 mutation effect는 안정적으로 예측하지 못한다.
```

따라서 이 모델은 완전한 일반 예측 모델이라기보다, DMS effect reproduction 및 interpolation model로 해석하는 것이 적절하다.

이 프로젝트의 강점은 단순히 높은 성능을 제시한 것이 아니라, validation setting을 세분화하여 모델이 실제로 어떤 상황에서 일반화되고 어떤 상황에서 실패하는지를 분석했다는 점이다.

## 한계

1. Structure feature는 6M0J 단일 정적 구조와 단순 거리 기준에 의존하였다.
2. ESM-2는 작은 모델을 사용했으므로 더 큰 protein language model을 사용할 경우 성능이 달라질 수 있다.
3. Reference ESM은 known-site interpolation에는 강했지만 unseen-site prediction을 해결하지 못했다.
4. Delta ESM은 site-held-out 성능을 소폭 개선했지만 충분하지 않았다.
5. 이 모델을 완전히 새로운 residue position의 효과를 안정적으로 예측하는 모델로 해석해서는 안 된다.

## 후속 연구

후속 연구에서는 다음 방향을 고려할 수 있다.

- 더 큰 protein language model 적용
- 여러 구조 conformation 사용
- solvent accessibility 추가
- residue conservation score 추가
- structural dynamics feature 추가
- high-error site 중심 생물학적 해석 강화
- 최종 보고서 및 발표 자료 제작
