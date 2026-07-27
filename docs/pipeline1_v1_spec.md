# Pipeline 1 v1 — `check_embedded_cubical_normality`

## 1. Purpose and Non-Goals

**Purpose**: Given a doubled-coordinate 3D maximal-voxel input, check the four conditions of Cubical Theorem 3 (purity / facet degree ∈ {1,2} / facet-adjacency connectivity / codim≥2 link pseudomanifold-ness) and emit a structured pass/fail certificate.

**Non-goals**:
- This backend does not certify the entire abstract regular cubical complex category Q5 — it only handles the embedded-grid special case of Q5.
- It does not repair or normalize a malformed abstract incidence poset.
- It does not check whether `F(K)` is an n-PCM/discrete surface — Pipeline 1 only checks whether `K` itself (and its links) is a cubical (normal) pseudomanifold. Determining whether `F(K)` is a PCM/surface is the job of a separate pipeline.

## 2. Input Contract

- A nonempty finite set of 3D integer voxel coordinates `(i, j, k) ∈ Z^3`.
- **Integer coordinate contract (finalized 2026 packaging)**: each coordinate must be a non-Boolean instance of `numbers.Integral`. Accepted integral scalars, including NumPy integer scalars, are normalized to a built-in `int` before dedup and geometric processing. `bool` is a subclass of `int` in Python but is explicitly rejected. Floating-point types are rejected even when the value looks integral (`1.0`, `np.float64(1.0)`, etc.). (`np.bool_` needs no separate exception handling — it is not an instance of `numbers.Integral`, so it is already filtered out by the above rule; this has been directly verified.)
- Input elements are interpreted as **maximal 3-cells** (top cells) — the entire face poset is not taken as input.
- Duplicate policy: **deduplicate** (since the input is a set of voxel coordinates, duplicates are automatically removed via `set()`; the number of `duplicates_removed` is recorded in the resulting report. Duplicates are not rejected — a duplicate coordinate is not itself an axiom violation, merely a redundancy of representation).
- Structural violations (tuples whose length is not 3, non-integer coordinates, empty input) are **rejected** — they are not silently fixed.
- `closure()` internally generates the full face-closed complex.
- No arbitrary structural repair/normalization is performed.
- (In this representation, distinct voxel coordinates automatically give distinct top cells, so interior overlap between top cells is inherently unrepresentable — it is not a separate check target.)

## 3. Mathematical Representation

- **Cell ID**: the doubled-coordinate integer tuple itself (`voxel_to_cell(v) = tuple(2x+1 for x in v)`).
- **`dim(cell)`**: the number of odd coordinates. Rank is identified directly with dim (no shift).
- **`is_face(sub, sup)`**: the even coordinates of `sup` must match `sub`, and the odd coordinates must be within ±1 of `sub`.
- **closure**: the set of all faces including the cell itself.
- **facet**: a dim=2 face of a top cell (dim=3).
- **`facet_parents`**: `facet -> [top cells having this facet]`, a single index shared by both the degree check and the connectivity check.
- **simplicial link (L2)**: `A_h(c) := {a ∈ K : h ⋖ a ≤ c}`, `Lk_K(h) := {A_h(c) : c ≥ h}`. Since `A_h(h) = ∅` (a phantom empty face) when `c=h`, the **`F^×` convention** (empty face removed) is used. For efficiency, the full set of link faces is generated as the downward closure of the link facets obtained from the maximal cofaces (top cells containing h) — this suffices due to L2's order-preservation property (`c1 ≤ c2 ⟹ A_h(c1) ⊆ A_h(c2)`).

## 4. The Four Check Conditions

| # | Mathematical condition | Implementation rule | Pass criterion | Failure certificate |
|---|---|---|---|---|
| 1 | Purity of `K` | Every face arises from the `closure()` of some top cell (automatic once the input passes §2 and all voxels share the same ambient dimension, =3) | Always PASS (once input validation passes) | N/A — a violation would already have been rejected at the input validation stage |
| 2 | Facet degree ∈ {1,2} | For each facet in `facet_parents`, `1 ≤ len(parents) ≤ 2` | All facets satisfy this | `(facet, parent_cells)` — the violating facet and its list of parents (**unreachable in v1**, see below) |
| 3 | Facet-adjacency connectivity | Using only degree=2 facets from `facet_parents`, build edges of the top-cell graph and check connectedness (`related()` / `connected_components_via_theta()` are **not used**) | Exactly 1 component | `components` — a list of top cells grouped by component |
| 4 | codim≥2 link pseudomanifold-ness | For each `h` (dim≤1), generate `Lk_K(h)` via the L2 formula, then directly check only **three** conditions on **that link itself**: purity (automatic) / ridge-degree ∈ {1,2} / ridge-adjacency connectivity (no recursion into normality) | All three conditions hold for every `h` | `(h, failed_link_condition, witness)` |

**Effective status of facet degree (v1-backend-only fact, verified 2025)**: in the doubled-coordinate unit-voxel representation, a dim=2 facet can have **exactly 2** candidate top cells at most — obtained by shifting one even-coordinate axis by ±1 (the other two coordinates are already fixed as odd, leaving no further candidates). Therefore the number of parents in `facet_parents` **cannot structurally exceed 2**, regardless of input — confirmed both by structural argument and by empirical testing over 200 random inputs (both documented in `test_check_embedded_cubical_normality.py`). In other words, the actual verification power of the four conditions breaks down as follows:

- purity: `PASS_BY_CONSTRUCTION` (input contract + closure)
- **facet degree: effectively `PASS_BY_CONSTRUCTION` in v1** — the actual counting logic in the code is not incorrect and is kept as-is (as a defensive check), but it is noted here that while this is a non-trivial condition in the abstract Q5 backend, in v1 the representation itself structurally guarantees it.
- facet connectivity: **a genuinely non-trivial check**
- codim≥2 link pseudomanifold-ness: **a genuinely non-trivial check**

## 5. Explicit Non-Reuse Items

> ⚠️ **Do not use `connected_components_via_theta()` on the top-cell set.**
> Two distinct top cells (being of the same dimension) are never in an `is_face` relation, so applying this function directly to the set of top cells will always incorrectly judge them as disconnected even when they share a common facet.
>
> ⚠️ **Do not recursively test normality inside a link.**
> Being D10-candidate (cubically normal) means each `Lk_K(σ)` satisfies `[A25]` **Definition 6** (pseudomanifold, single level), not that it satisfies **Definition 10** (normal pseudomanifold, recursing into the link's link).
>
> ⚠️ **Do not describe this backend as covering all abstract Q5 complexes.**
> The function is kept named `check_embedded_cubical_normality`; a name like `check_cubical_normality`, which could be misread as certifying the entirety of Q5, is not used.

## 6. Output Contract

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

## 6-A. Packaging Layer (finalized 2026 packaging): `build_certificate`

The output example in §6 is the raw return value of `check_embedded_cubical_normality()`, and since this function's contract is directly verified by 78 regression tests, it is not changed. The paper/appendix-facing output is handled by a separate function `build_certificate(voxels)`, which adds the following:

- `schema_version` and `backend` fields.
- Extends `overall` to three values: `PASS` / `FAIL` / `INVALID_INPUT` (`INVALID_INPUT` is for when the complex itself could not be constructed; `FAIL` is for when a valid complex violates a condition — these are distinct meanings, hence separated). A fourth value, `INTERNAL_ERROR`, is added not by `build_certificate()` itself but only by the `build_certificate_cli()` wrapper when it catches an unexpected internal exception (see below) — earlier wording that did not make this distinction explicit is hereby corrected.
- On structural violations, records not just the first error but **all** discovered errors in `validation.failures`, structured with code/item_index/coordinate_index/value_repr/message.
- Failures of `facet_degree`/`facet_connectivity`/`links` are also recorded as structured witnesses with codes attached (instead of things like `FACET_DEGREE_OUT_OF_RANGE`/`TOP_CELL_FACET_DISCONNECTED`, the actual implementation uses the codes `LINK_RIDGE_DEGREE_OUT_OF_RANGE` / `LINK_RIDGE_DISCONNECTED` — the only failure types actually possible in the current v1 backend are facet-connectivity and ridge-degree/ridge-connectivity; facet-degree failure is structurally unreachable in this backend per §4's argument, but the field itself is kept for compatibility with the abstract backend).
- All cells/components/witnesses are serialized as sorted lists rather than tuples/sets/frozensets — confirmed by regression test (`test_certificate_deterministic_under_input_order_permutation`) that the same JSON results regardless of input order.
- `build_certificate_cli(voxels)`: a CLI-boundary wrapper that wraps only unexpected internal exceptions as `INTERNAL_ERROR`. Library users are expected to call `build_certificate()` directly, letting real bugs raise normally.

## 7. Implementation Order

1. Input validation
2. Build `facet_parents`
3. Check degree and connectivity
4. Generate simplicial link (L2 formula, maximal-coface + downward closure)
5. Simplicial pseudomanifold checker (link only, non-recursive)
6. Output combined results and witnesses
7. Small hand-crafted unit tests (single voxel / two voxels sharing a facet / two voxels sharing only an edge — pinch case)
