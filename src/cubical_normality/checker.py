"""
Pipeline 1 v1 — check_embedded_cubical_normality

Implements pipeline1_v1_spec.md exactly. Scope: doubled-coordinate 3D
maximal-voxel cubical complexes only (the embedded-grid special case of
Q5), NOT abstract regular cubical complexes in general.

Does NOT reuse connected_components_via_theta()/related() for facet
connectivity (see spec §5 — those test face-inclusion adjacency, which
never holds between two distinct top cells). Does NOT recursively test
normality inside a link (spec §5 — only single-level pseudomanifold-ness
of each link is checked).
"""

from itertools import product, combinations
from numbers import Integral


# ---------------------------------------------------------------------
# Section 3: mathematical representation
# ---------------------------------------------------------------------

def voxel_to_cell(v):
    """Integer voxel coordinate (i, j, k) -> doubled coordinate (2i+1, 2j+1, 2k+1)."""
    return tuple(2 * x + 1 for x in v)


def dim(cell):
    """Dimension of a cell = number of odd (free) coordinates."""
    return sum(1 for c in cell if c % 2 == 1)


def is_face(sub, sup):
    """True iff sub is a face of sup (sub == sup counts as a face)."""
    for s, S in zip(sub, sup):
        if S % 2 == 0:
            if s != S:
                return False
        else:
            if abs(s - S) > 1:
                return False
    return True


def closure(cell):
    """All faces of `cell`, including itself."""
    options = []
    for c in cell:
        if c % 2 == 1:
            options.append([c - 1, c, c + 1])
        else:
            options.append([c])
    return set(product(*options))


# ---------------------------------------------------------------------
# Section 2: input contract
# ---------------------------------------------------------------------

class InputValidationError(ValueError):
    pass


def _is_valid_coordinate(x):
    """
    Coordinate contract (packaging session 37): a coordinate is valid iff
    it is a non-Boolean instance of numbers.Integral. This accepts
    built-in int and integral scalar types from other libraries (e.g.
    NumPy integer scalars), and rejects bool even though bool is a
    subclass of int in Python, and rejects any float even when
    numerically integral (float is not a numbers.Integral instance).

    Note: numpy.bool_ does not need special-casing here — it is not an
    instance of numbers.Integral in the first place (verified
    empirically), so `isinstance(x, Integral)` alone already excludes
    it without a separate numpy-specific check.
    """
    return isinstance(x, Integral) and not isinstance(x, bool)


def validate_and_normalize_input(voxels):
    """
    Enforce spec §2 (as extended in the packaging session). Returns
    (deduped_voxel_list, duplicates_removed). Raises InputValidationError
    on structural violations. Never repairs a structural violation
    silently.

    Accepted coordinates are normalized to builtin `int` before
    deduplication and all downstream geometric processing, so that
    e.g. a NumPy integer scalar and the equal-valued builtin int are
    always treated as the same coordinate.
    """
    if voxels is None:
        raise InputValidationError("input is None")

    voxel_list = list(voxels)
    if len(voxel_list) == 0:
        raise InputValidationError("input voxel set is empty")

    normalized = []
    for v in voxel_list:
        if not isinstance(v, tuple) or len(v) != 3:
            raise InputValidationError(
                f"voxel {v!r} is not a length-3 tuple")
        if not all(_is_valid_coordinate(x) for x in v):
            raise InputValidationError(
                f"voxel {v!r} has non-integer coordinates "
                f"(Boolean and floating-point coordinates are rejected; "
                f"only non-Boolean numbers.Integral values are accepted)")
        normalized.append(tuple(int(x) for x in v))

    deduped = sorted(set(normalized))
    duplicates_removed = len(normalized) - len(deduped)
    return deduped, duplicates_removed


# ---------------------------------------------------------------------
# facet_parents index (shared by conditions 2 and 3)
# ---------------------------------------------------------------------

def build_facet_parents(top_cells):
    """
    top_cells: list of doubled-coordinate top cells (dim 3).
    Returns dict: facet (dim-2 cell) -> list of top cells having it as a facet.
    """
    facet_parents = {}
    for c in top_cells:
        for f in closure(c):
            if dim(f) == 2:
                facet_parents.setdefault(f, []).append(c)
    return facet_parents


# ---------------------------------------------------------------------
# Condition 2: facet degree
# ---------------------------------------------------------------------

def check_facet_degree(facet_parents):
    failures = []
    for facet, parents in facet_parents.items():
        if not (1 <= len(parents) <= 2):
            failures.append((facet, list(parents)))
    status = "PASS" if not failures else "FAIL"
    return {"status": status, "failures": failures}


# ---------------------------------------------------------------------
# Condition 3: facet-adjacency connectivity
# NOTE: deliberately does NOT use related()/connected_components_via_theta().
# ---------------------------------------------------------------------

def check_facet_connectivity(top_cells, facet_parents):
    parent_of = {c: c for c in top_cells}

    def find(x):
        while parent_of[x] != x:
            parent_of[x] = parent_of[parent_of[x]]
            x = parent_of[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent_of[ra] = rb

    for facet, parents in facet_parents.items():
        if len(parents) == 2:
            union(parents[0], parents[1])

    components = {}
    for c in top_cells:
        components.setdefault(find(c), []).append(c)
    components_list = list(components.values())

    status = "PASS" if len(components_list) == 1 else "FAIL"
    return {"status": status, "components": components_list}


# ---------------------------------------------------------------------
# Condition 1: purity — PASS_BY_CONSTRUCTION given validated input
# ---------------------------------------------------------------------

def check_purity():
    return {
        "status": "PASS_BY_CONSTRUCTION",
        "justification": (
            "every generated face is a member of closure(voxel_to_cell(v)) "
            "for some validated input voxel v, hence a face of a top cell "
            "by construction; no separate structural check is needed"
        ),
    }


# ---------------------------------------------------------------------
# Condition 4: codim>=2 link pseudomanifold-ness (single level, no recursion)
# ---------------------------------------------------------------------

def atoms_above(h, all_cells):
    """Covers of h (dim(h)+1 elements a with h a face of a) within all_cells."""
    return [a for a in all_cells if dim(a) == dim(h) + 1 and is_face(h, a)]


def build_link_faces(h, all_cells, top_cells):
    """
    Lk_K(h) via L2: A_h(c) = {a in atoms(h) : a <= c}, for c a maximal
    coface of h (top cell containing h). Downward closure of these facets
    gives the full link (F-times convention: no empty face).
    Returns (link_facets: list[frozenset], all_link_faces: set[frozenset]).
    """
    atoms = atoms_above(h, all_cells)
    cofaces = [c for c in top_cells if is_face(h, c)]

    link_facets = []
    for c in cofaces:
        facet = frozenset(a for a in atoms if is_face(a, c))
        if facet:
            link_facets.append(facet)
    link_facets = list(set(link_facets))

    all_faces = set()
    for facet in link_facets:
        k = len(facet)
        for size in range(1, k + 1):
            for sub in combinations(sorted(facet), size):
                all_faces.add(frozenset(sub))
    return link_facets, all_faces


def check_link_pseudomanifold(h, link_facets):
    """
    Single-level pseudomanifold check on the simplicial link generated at h:
    purity (by construction), ridge-degree in {1,2}, ridge-adjacency
    connectivity. Does NOT recurse into normality of the link.
    """
    if not link_facets:
        return {"status": "PASS", "reason": "empty link (h is itself a top cell boundary edge case)"}

    facet_size = len(link_facets[0])
    if any(len(f) != facet_size for f in link_facets):
        return {"status": "FAIL", "reason": "link facets of inconsistent size (purity violation)"}

    ridge_parents = {}
    for facet in link_facets:
        if facet_size == 1:
            # rank-0 link: nothing further to check (single point per facet)
            continue
        for ridge in combinations(sorted(facet), facet_size - 1):
            ridge_parents.setdefault(frozenset(ridge), []).append(facet)

    if facet_size == 1:
        return {"status": "PASS", "reason": "rank-0 link, trivially a pseudomanifold"}

    degree_failures = [(r, ps) for r, ps in ridge_parents.items() if not (1 <= len(ps) <= 2)]
    if degree_failures:
        return {"status": "FAIL", "reason": "ridge_degree", "witness": degree_failures}

    parent_of = {f: f for f in link_facets}

    def find(x):
        while parent_of[x] != x:
            parent_of[x] = parent_of[parent_of[x]]
            x = parent_of[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent_of[ra] = rb

    for ridge, ps in ridge_parents.items():
        if len(ps) == 2:
            union(ps[0], ps[1])

    n_components = len({find(f) for f in link_facets})
    if n_components != 1:
        return {"status": "FAIL", "reason": "ridge_connectivity", "witness": n_components}

    return {"status": "PASS"}


def check_links(all_cells, top_cells):
    codim_ge_2_faces = [h for h in all_cells if dim(h) <= 1]
    failures = []
    for h in codim_ge_2_faces:
        link_facets, _ = build_link_faces(h, all_cells, top_cells)
        result = check_link_pseudomanifold(h, link_facets)
        if result["status"] == "FAIL":
            failures.append((h, result.get("reason"), result.get("witness")))
    status = "PASS" if not failures else "FAIL"
    return {"status": status, "checked_faces": len(codim_ge_2_faces), "failures": failures}


# ---------------------------------------------------------------------
# Top level
# ---------------------------------------------------------------------

def check_embedded_cubical_normality(voxels):
    try:
        deduped, duplicates_removed = validate_and_normalize_input(voxels)
    except InputValidationError as e:
        return {"valid_input": False, "error": str(e), "overall": "FAIL"}

    top_cells = [voxel_to_cell(v) for v in deduped]
    all_cells = set()
    for c in top_cells:
        all_cells |= closure(c)

    facet_parents = build_facet_parents(top_cells)

    purity_result = check_purity()
    degree_result = check_facet_degree(facet_parents)
    connectivity_result = check_facet_connectivity(top_cells, facet_parents)
    links_result = check_links(all_cells, top_cells)

    overall = "PASS" if all(
        r["status"] in ("PASS", "PASS_BY_CONSTRUCTION")
        for r in (purity_result, degree_result, connectivity_result, links_result)
    ) else "FAIL"

    return {
        "valid_input": True,
        "duplicates_removed": duplicates_removed,
        "purity": purity_result,
        "facet_degree": degree_result,
        "facet_connectivity": connectivity_result,
        "links": links_result,
        "overall": overall,
    }


# ---------------------------------------------------------------------
# Packaging layer (session 37, §5): JSON-serializable certificate.
#
# This is additive on top of check_embedded_cubical_normality() above --
# it does not change that function's return contract (the 78-test
# regression suite tests that contract directly and must keep passing
# unmodified). build_certificate() re-runs the same core checks and
# repackages them into a schema-versioned, deterministic, JSON-safe
# format suitable for a paper's computational appendix.
# ---------------------------------------------------------------------

CERTIFICATE_SCHEMA_VERSION = "1.0"
CERTIFICATE_BACKEND = "embedded_doubled_coordinate_3d"

_LINK_FAILURE_CODES = {
    "ridge_degree": "LINK_RIDGE_DEGREE_OUT_OF_RANGE",
    "ridge_connectivity": "LINK_RIDGE_DISCONNECTED",
    "purity": "LINK_PURITY_VIOLATION",
}


def _empty_certificate_shell(received_count, failures):
    return {
        "schema_version": CERTIFICATE_SCHEMA_VERSION,
        "backend": CERTIFICATE_BACKEND,
        "valid_input": False,
        "input": {
            "received_count": received_count,
            "normalized_count": 0,
            "duplicates_removed": 0,
        },
        "validation": {"status": "FAIL", "failures": failures},
        "purity": {"status": "NOT_RUN"},
        "facet_degree": {"status": "NOT_RUN", "failures": []},
        "facet_connectivity": {"status": "NOT_RUN", "component_count": 0, "components": []},
        "links": {"status": "NOT_RUN", "checked_faces": 0, "failures": []},
        "overall": "INVALID_INPUT",
    }


def build_certificate(voxels):
    """
    Build the schema-versioned JSON-safe certificate described in
    pipeline1_v1_spec.md's packaging notes.

    Differences from check_embedded_cubical_normality()'s raw dict:
      - "overall" distinguishes INVALID_INPUT (complex could not be
        built) from FAIL (a valid complex violated a condition).
      - validation collects ALL structural failures, not just the
        first one encountered.
      - failures are structured records with a `code` field, not bare
        tuples/strings.
      - all coordinates/cells are serialized as sorted, JSON-safe lists
        (no tuple/set/frozenset ever reaches the returned dict), so the
        certificate is deterministic under input reordering.
    """
    try:
        voxel_list = list(voxels) if voxels is not None else None
    except TypeError:
        voxel_list = None
    received_count = len(voxel_list) if voxel_list is not None else 0

    failures = []
    if voxel_list is None:
        failures.append({
            "code": "NULL_INPUT", "item_index": None, "coordinate_index": None,
            "value_repr": repr(voxels), "message": "input is None",
        })
    elif len(voxel_list) == 0:
        failures.append({
            "code": "EMPTY_INPUT", "item_index": None, "coordinate_index": None,
            "value_repr": "[]", "message": "input voxel set is empty",
        })
    else:
        for i, v in enumerate(voxel_list):
            if not isinstance(v, tuple) or len(v) != 3:
                failures.append({
                    "code": "NOT_LENGTH_3_TUPLE", "item_index": i, "coordinate_index": None,
                    "value_repr": repr(v), "message": f"voxel {v!r} is not a length-3 tuple",
                })
                continue
            for j, x in enumerate(v):
                if not _is_valid_coordinate(x):
                    code = "BOOLEAN_COORDINATE" if isinstance(x, bool) else "NON_INTEGRAL_COORDINATE"
                    msg = ("Boolean coordinates are not accepted as integer voxel coordinates."
                           if code == "BOOLEAN_COORDINATE" else
                           f"coordinate {x!r} is not a non-Boolean numbers.Integral value")
                    failures.append({
                        "code": code, "item_index": i, "coordinate_index": j,
                        "value_repr": repr(x), "message": msg,
                    })

    if failures:
        return _empty_certificate_shell(received_count, failures)

    # --- valid input: normalize, dedupe, run the core checks ---
    normalized = [tuple(int(x) for x in v) for v in voxel_list]
    deduped = sorted(set(normalized))
    duplicates_removed = len(normalized) - len(deduped)

    top_cells = [voxel_to_cell(v) for v in deduped]
    all_cells = set()
    for c in top_cells:
        all_cells |= closure(c)

    facet_parents = build_facet_parents(top_cells)
    purity_result = check_purity()
    degree_result = check_facet_degree(facet_parents)
    connectivity_result = check_facet_connectivity(top_cells, facet_parents)
    links_result = check_links(all_cells, top_cells)

    degree_failures = sorted(
        (
            {
                "code": "FACET_DEGREE_OUT_OF_RANGE",
                "facet": list(facet),
                "degree": len(parents),
                "parent_cells": sorted(list(p) for p in parents),
            }
            for facet, parents in degree_result["failures"]
        ),
        key=lambda f: (f["facet"], f["code"]),
    )

    components_serialized = sorted(
        (sorted(list(c) for c in comp) for comp in connectivity_result["components"]),
        key=lambda comp: comp[0] if comp else [],
    )

    link_failures = []
    for h, reason, witness in links_result["failures"]:
        entry = {
            "code": _LINK_FAILURE_CODES.get(reason, "LINK_FAIL_UNKNOWN"),
            "base_face": list(h),
            "base_face_dimension": dim(h),
            "link_condition": reason,
            "witness": None,
        }
        if reason == "ridge_degree" and witness:
            entry["witness"] = sorted(
                (
                    {
                        "ridge": sorted(list(r)),
                        "degree": len(parent_facets),
                        "parent_facets": sorted(sorted(list(f)) for f in parent_facets),
                    }
                    for r, parent_facets in witness
                ),
                key=lambda w: w["ridge"],
            )
        elif reason == "ridge_connectivity":
            entry["witness"] = {"component_count": witness}
        link_failures.append(entry)

    link_failures.sort(key=lambda e: (e["base_face"], e["code"]))

    overall = "PASS" if all(
        r["status"] in ("PASS", "PASS_BY_CONSTRUCTION")
        for r in (purity_result, degree_result, connectivity_result, links_result)
    ) else "FAIL"

    return {
        "schema_version": CERTIFICATE_SCHEMA_VERSION,
        "backend": CERTIFICATE_BACKEND,
        "valid_input": True,
        "input": {
            "received_count": received_count,
            "normalized_count": len(deduped),
            "duplicates_removed": duplicates_removed,
        },
        "validation": {"status": "PASS", "failures": []},
        "purity": purity_result,
        "facet_degree": {
            "status": degree_result["status"],
            "checked_facets": len(facet_parents),
            "failures": degree_failures,
        },
        "facet_connectivity": {
            "status": connectivity_result["status"],
            "component_count": len(connectivity_result["components"]),
            "components": components_serialized,
        },
        "links": {
            "status": links_result["status"],
            "checked_faces": links_result["checked_faces"],
            "failures": link_failures,
        },
        "overall": overall,
    }


def build_certificate_cli(voxels):
    """
    CLI-mode wrapper around build_certificate(): converts any
    *unexpected* internal exception into an INTERNAL_ERROR certificate
    instead of propagating a raw traceback. Structural input problems
    are NOT exceptions here (build_certificate already reports those as
    "overall": "INVALID_INPUT"); this only catches genuine implementation
    bugs. Library callers should call build_certificate() directly and
    let such bugs raise, per packaging §6.
    """
    try:
        return build_certificate(voxels)
    except Exception as e:  # intentional catch-all: this is the CLI boundary
        return {
            "schema_version": CERTIFICATE_SCHEMA_VERSION,
            "backend": CERTIFICATE_BACKEND,
            "valid_input": None,
            "input": {"received_count": None, "normalized_count": 0, "duplicates_removed": 0},
            "validation": {"status": "NOT_RUN", "failures": []},
            "purity": {"status": "NOT_RUN"},
            "facet_degree": {"status": "NOT_RUN", "failures": []},
            "facet_connectivity": {"status": "NOT_RUN", "component_count": 0, "components": []},
            "links": {"status": "NOT_RUN", "checked_faces": 0, "failures": []},
            "overall": "INTERNAL_ERROR",
            "internal_error": {"type": type(e).__name__, "message": str(e)},
        }


if __name__ == "__main__":
    print("Regression tests moved to test_check_embedded_cubical_normality.py "
          "(run: pytest test_check_embedded_cubical_normality.py -v)")
