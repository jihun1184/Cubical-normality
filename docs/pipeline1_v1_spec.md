# Pipeline 1 v1 — `check_embedded_cubical_normality`

## 1. 목적과 비목적

**목적**: doubled-coordinate 3D maximal-voxel 입력이 주어졌을 때, Cubical
Theorem 3의 네 조건(purity / facet degree∈{1,2} / facet-adjacency
connectivity / codim≥2 link pseudomanifold성)을 검사하고, 구조화된
pass/fail 증명서를 낸다.

**비목적**:
- abstract regular cubical complex 전체 Q5 범주를 인증하지 않는다 — 이
  backend는 Q5의 embedded-grid 특수 케이스만 다룬다.
- malformed abstract incidence poset을 복구하거나 normalize하지 않는다.
- `F(K)`가 n-PCM/discrete surface인지는 검사하지 않는다 — Pipeline 1은
  `K` 자체(와 그 link들)가 cubical (normal) pseudomanifold인지만 검사한다.
  `F(K)`의 PCM/surface 판정은 별도 pipeline의 몫이다.

## 2. 입력 계약

- nonempty finite set of 3D integer voxel coordinates `(i, j, k) ∈ Z^3`.
- **정수 좌표 계약 (2026-패키징 확정)**: 각 좌표는 `numbers.Integral`의
  non-Boolean 인스턴스여야 한다. NumPy 정수 스칼라를 포함한 accepted
  integral scalar는 dedup 및 기하 처리 전에 built-in `int`로 정규화된다.
  `bool`은 Python에서 `int`의 subclass이지만 명시적으로 reject한다.
  실수형은 값이 정수처럼 보여도 (`1.0`, `np.float64(1.0)` 등) reject한다.
  (`np.bool_`는 별도 예외처리가 필요 없음 — `numbers.Integral`의
  인스턴스가 아니므로 이미 위 규칙으로 걸러짐, 직접 검증 완료.)
- 입력 원소는 **maximal 3-cell**(top cell)로 해석된다 — face poset 전체를
  입력받지 않는다.
- 중복 정책: **deduplicate** (voxel 좌표 집합이므로 `set()`으로 자동
  중복 제거; 결과 보고서에 `duplicates_removed` 개수를 기록한다. reject는
  하지 않는다 — 중복 좌표 자체는 공리 위반이 아니라 표현상의 중복이므로).
- 구조 위반(길이가 3이 아닌 tuple, 정수가 아닌 좌표, 빈 입력)은 **reject**
  — 조용히 고치지 않는다.
- `closure()`로 전체 face-closed complex를 내부에서 생성한다.
- 임의의 구조적 repair/normalization은 하지 않는다.
- (이 representation에서는 서로 다른 voxel 좌표가 자동으로 서로 다른
  top cell을 주므로, top cell 간 interior overlap은 애초에 표현 불가능
  — 별도 검사 대상이 아니다.)

## 3. 수학적 표현

- **Cell ID**: doubled-coordinate 정수 tuple 그 자체
  (`voxel_to_cell(v) = tuple(2x+1 for x in v)`).
- **`dim(cell)`**: 홀수 좌표 개수. Rank는 dim과 그대로 identify (shift 없음).
- **`is_face(sub, sup)`**: `sup`의 짝수좌표는 `sub`와 일치해야 하고,
  홀수좌표는 `sub`가 ±1 이내여야 한다.
- **closure**: cell 자신을 포함한 모든 face의 집합.
- **facet**: top cell(dim=3)의 dim=2인 face.
- **`facet_parents`**: `facet -> [top cells having this facet]`, degree와
  connectivity 검사가 공유하는 단일 인덱스.
- **simplicial link (L2)**: `A_h(c) := {a ∈ K : h ⋖ a ≤ c}`,
  `Lk_K(h) := {A_h(c) : c ≥ h}`. `c=h`일 때 `A_h(h)=∅`(phantom empty
  face)이므로 **`F^×` convention**(빈 face 제거)을 쓴다. 효율을 위해
  maximal coface(h를 포함하는 top cell)에서 얻은 link facet들의
  하향폐포(downward closure)로 전체 link face를 생성한다 — L2의
  order-preservation(`c1≤c2 ⟹ A_h(c1)⊆A_h(c2)`)에 의해 이것으로 충분함이
  보장된다.

## 4. 네 검사 조건

| # | Mathematical condition | Implementation rule | Pass criterion | Failure certificate |
|---|---|---|---|---|
| 1 | Purity of `K` | 모든 face는 어떤 top cell의 `closure()`에서 나옴 (입력이 §2를 통과하고 모든 voxel이 동일 ambient dim(=3)이면 자동) | 항상 PASS (입력 validation 통과 시) | 해당 없음 — 위반은 입력 validation 단계에서 이미 reject됨 |
| 2 | Facet degree ∈ {1,2} | `facet_parents`의 각 facet에 대해 `1 ≤ len(parents) ≤ 2` | 모든 facet이 만족 | `(facet, parent_cells)` — 위반 facet과 그 parent 목록 (**v1에서는 도달 불가**, 아래 참조) |
| 3 | Facet-adjacency connectivity | `facet_parents`에서 degree=2인 facet만으로 top-cell 그래프의 edge를 만들고 connected 여부 검사 (`related()`/`connected_components_via_theta()` **미사용**) | 컴포넌트 1개 | `components` — top cell을 컴포넌트별로 묶은 리스트 |
| 4 | codim≥2 link pseudomanifold성 | 각 `h`(dim≤1)에 대해 `Lk_K(h)`를 L2 공식으로 생성한 뒤, **그 link 자체**의 purity(자동)/ridge-degree∈{1,2}/ridge-adjacency connectivity **세 조건만** 직접 검사 (normality 재귀 없음) | 모든 `h`에서 세 조건 만족 | `(h, failed_link_condition, witness)` |

**Facet degree의 실질적 지위 (v1 backend 전용 사실, 2025-검증)**: doubled-coordinate
unit-voxel 표현에서, dim=2인 facet은 짝수좌표 축 하나를 ±1 옮긴 **정확히 2개**의
후보 top cell만 가질 수 있다(다른 두 좌표는 이미 홀수로 고정돼 있어 나머지 후보가
없음). 따라서 `facet_parents`의 parent 수는 입력이 무엇이든 **구조적으로 2를
초과할 수 없다** — 구조적 논증과 200회 랜덤 입력 실증 검사(둘 다
`test_check_embedded_cubical_normality.py` 참조)로 확인됨. 즉 네 조건의 실질
검증력은 다음과 같이 나뉜다:

- purity: `PASS_BY_CONSTRUCTION` (입력 계약 + closure)
- **facet degree: v1에서 사실상 `PASS_BY_CONSTRUCTION`** — 코드의 실제 카운팅
  로직 자체는 틀리지 않았고 그대로 유지하지만(방어적 검사), abstract Q5
  backend에서는 이게 비자명한 조건이라는 점과 v1에서는 표현 자체가 이를
  구조적으로 보장한다는 점의 차이를 여기 명시해둔다.
- facet connectivity: **실제 비자명 검사**
- codim≥2 link pseudomanifold성: **실제 비자명 검사**

## 5. 명시적 비재사용 항목

> ⚠️ **Do not use `connected_components_via_theta()` on the top-cell set.**
> 서로 다른 두 top cell은 (같은 차원이므로) 절대 `is_face` 관계가 아니라서,
> 이 함수를 top-cell 집합에 직접 적용하면 공통 facet을 공유해도 항상
> 비연결로 잘못 판정된다.
>
> ⚠️ **Do not recursively test normality inside a link.**
> D10-candidate(cubically normal)는 각 `Lk_K(σ)`가 `[A25]` **Definition 6**
> (pseudomanifold, 단일 레벨)을 만족하라는 것이지, **Definition 10**
> (normal pseudomanifold, link의 link까지 재귀)을 만족하라는 게 아니다.
>
> ⚠️ **Do not describe this backend as covering all abstract Q5 complexes.**
> 함수명은 `check_embedded_cubical_normality`로 두고, `check_cubical_
> normality`처럼 Q5 전체를 인증하는 것으로 오해될 수 있는 이름은 쓰지 않는다.

## 6. 출력 계약

```json
{
    "valid_input": true,
    "duplicates_removed": 0,
    "purity": {"status": "PASS_BY_CONSTRUCTION"},
    "facet_degree": {"status": "PASS", "failures": []},
    "facet_connectivity": {"status": "PASS", "components": [[...]]},
    "links": {"status": "PASS", "checked_faces": 0, "failures": []},
    "overall": "PASS"
}
```

## 6-A. 패키징 레이어 (2026-패키징 확정): `build_certificate`

§6의 출력 예시는 `check_embedded_cubical_normality()`의 raw 반환값이며,
이 함수의 계약은 회귀테스트 78개가 직접 검증하므로 변경하지 않는다.
논문/부록용 산출물은 별도 함수 `build_certificate(voxels)`가 담당하며,
다음을 추가한다:

- `schema_version`, `backend` 필드.
- `overall`을 `PASS` / `FAIL` / `INVALID_INPUT` 세 값으로 확장
  (`INVALID_INPUT`은 complex 자체를 구성 못한 경우, `FAIL`은 유효한
  complex가 조건을 위반한 경우 — 서로 다른 의미이므로 분리). 네 번째 값
  `INTERNAL_ERROR`는 `build_certificate()` 자체가 아니라
  `build_certificate_cli()` wrapper가 예기치 못한 내부 예외를 잡을 때만
  추가되는 값이다 (아래 참조) — 이 구분을 명시하지 않은 이전 문구는
  정정한다.
- 구조 위반 시 첫 번째 오류만이 아니라 발견된 모든 오류를 `validation.
  failures`에 code/item_index/coordinate_index/value_repr/message로
  구조화해 기록.
- `facet_degree`/`facet_connectivity`/`links`의 실패도 code가 붙은
  구조화된 witness로 기록 (`FACET_DEGREE_OUT_OF_RANGE`,
  `TOP_CELL_FACET_DISCONNECTED`류 대신 실제 구현은
  `LINK_RIDGE_DEGREE_OUT_OF_RANGE`/`LINK_RIDGE_DISCONNECTED` 코드 사용 —
  현재 v1 백엔드에서 실제로 발생 가능한 실패 종류는 facet-connectivity와
  ridge-degree/ridge-connectivity 뿐이며, facet-degree 실패는 §4의
  구조적 논증대로 이 backend에서 도달 불가능하지만 필드 자체는 abstract
  backend 호환을 위해 유지).
- 모든 cell/component/witness를 tuple/set/frozenset이 아닌 정렬된
  list로 직렬화 — 입력 순서에 무관하게 동일한 JSON이 나옴을 회귀테스트로
  확인(`test_certificate_deterministic_under_input_order_permutation`).
- `build_certificate_cli(voxels)`: 예기치 못한 내부 예외만
  `INTERNAL_ERROR`로 포장하는 CLI 경계용 wrapper. 라이브러리 사용자는
  `build_certificate()`를 직접 호출해 실제 버그는 그대로 raise되게 둔다.

## 7. 구현 순서

1. 입력 validation
2. `facet_parents` 구축
3. degree와 connectivity 검사
4. simplicial link 생성 (L2 공식, maximal-coface + downward closure)
5. simplicial pseudomanifold checker (link 전용, 비재귀)
6. 통합 결과와 witness 출력
7. 작은 수작업 사례 unit test (단일 voxel / facet-공유 두 voxel / edge만
   공유하는 두 voxel — pinch 사례)
