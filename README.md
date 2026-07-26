# Cubical Normality Checker

Reference implementation for the paper's computational checker of the embedded-grid special case of cubical normal pseudomanifold conditions.

## Scope

The input is a nonempty finite collection of maximal 3D voxels `(i, j, k) ∈ Z³`. Internally, voxels are mapped to doubled coordinates and expanded to the face-closed cubical complex. The checker evaluates:

1. purity (guaranteed by construction for validated maximal-voxel input),
2. facet degree in `{1, 2}`,
3. connectivity of top cells through shared facets, and
4. single-level pseudomanifold conditions for links of codimension-at-least-two faces.

This release **does not** certify arbitrary abstract regular cubical complexes, recursively test normality inside links, or decide whether the face poset is a PCM/discrete surface. See `docs/pipeline1_v1_spec.md` for the precise contract.

## Installation

```bash
python -m pip install -e ".[test]"
```

## Python API

```python
from cubical_normality import build_certificate

voxels = [(0, 0, 0), (1, 0, 0)]
certificate = build_certificate(voxels)
print(certificate["overall"])  # PASS
```

`build_certificate()` returns a deterministic, JSON-serializable, schema-versioned certificate. For library use, unexpected implementation errors are intentionally propagated. `build_certificate_cli()` converts such errors into an `INTERNAL_ERROR` certificate at the command-line boundary.

## Command line

```bash
cubical-normality examples/pass_two_face_adjacent_voxels.json
cubical-normality examples/fail_diagonal_edge_pinch.json
cat examples/pass_two_face_adjacent_voxels.json | cubical-normality -
```

The command exits with code `0` only for an overall `PASS`; valid mathematical failures and invalid inputs exit with code `1`.

## Tests

```bash
pytest
```

The included suite covers input validation, deterministic certificates, translation and axis-permutation equivariance, disconnected configurations, vertex-only contacts, and the edge-quadrant path/cycle/pinch family.

## Repository contents

- `src/cubical_normality/checker.py`: canonical mathematical implementation
- `src/cubical_normality/cli.py`: thin JSON command-line interface
- `tests/test_checker.py`: canonical regression tests
- `docs/pipeline1_v1_spec.md`: exact mathematical and software contract
- `examples/`: one passing and one failing minimal input

## Citation and license

This project is licensed under the BSD 3-Clause License. See [`LICENSE`](LICENSE) for details. Citation metadata for release `v1.0.0` is provided in [`CITATION.cff`](CITATION.cff); update its `date-released` field to match the date of the GitHub release/tag it accompanies.
