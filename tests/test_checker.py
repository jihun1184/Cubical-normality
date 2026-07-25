"""
Regression tests for check_embedded_cubical_normality (Pipeline 1 v1).

Priority per this session's plan: harden the decision boundary of the
current 3D checker with small adversarial/positive pairs, before widening
mathematical scope (n in {0,1} base cases) or moving to a large
hand-verified complex.

Sections:
  1. Basic pass/reject sanity (single voxel, malformed input)
  2. Facet-degree structural bound (v1 backend only, see spec)
  3. The edge-quadrant family: 2-adjacent / 2-opposite (pinch) /
     vertex-only / far-apart / 3-quadrant (path) / 4-quadrant (cycle)
  4. Invariance/equivariance: input-order permutation, translation,
     axis permutation
"""
import itertools
import random

import pytest

import json

from cubical_normality.checker import (
    check_embedded_cubical_normality,
    build_facet_parents,
    voxel_to_cell,
    build_certificate,
    build_certificate_cli,
)


# ---------------------------------------------------------------------
# 1. Basic sanity
# ---------------------------------------------------------------------

def test_single_voxel_passes_everything():
    r = check_embedded_cubical_normality([(0, 0, 0)])
    assert r["overall"] == "PASS"
    assert r["purity"]["status"] == "PASS_BY_CONSTRUCTION"
    assert r["facet_degree"]["status"] == "PASS"
    assert r["facet_connectivity"]["status"] == "PASS"
    assert len(r["facet_connectivity"]["components"]) == 1
    assert r["links"]["status"] == "PASS"
    # 8 vertices + 12 edges of a single cube
    assert r["links"]["checked_faces"] == 20


def test_malformed_input_is_rejected_not_repaired():
    r = check_embedded_cubical_normality([(0, 0, 0), (1, 0)])
    assert r["valid_input"] is False
    assert r["overall"] == "FAIL"
    assert "length-3" in r["error"]


def test_duplicate_voxels_are_deduplicated_not_rejected():
    r = check_embedded_cubical_normality([(0, 0, 0), (0, 0, 0), (1, 0, 0)])
    assert r["valid_input"] is True
    assert r["duplicates_removed"] == 1
    assert r["overall"] == "PASS"


def test_empty_input_rejected():
    r = check_embedded_cubical_normality([])
    assert r["valid_input"] is False


def test_boolean_coordinate_is_rejected():
    # bool is a subclass of int in Python, so a naive isinstance(x, int)
    # check would silently accept this. The packaging-session contract
    # (numbers.Integral, non-bool) must reject it explicitly.
    r = check_embedded_cubical_normality([(True, 0, 0)])
    assert r["valid_input"] is False
    assert "non-integer" in r["error"] or "Boolean" in r["error"]


def test_numpy_integer_coordinate_is_accepted_and_normalized():
    np = pytest.importorskip("numpy")
    r_numpy = check_embedded_cubical_normality(
        [(np.int64(0), np.int64(0), np.int64(0))])
    r_plain = check_embedded_cubical_normality([(0, 0, 0)])
    assert r_numpy["valid_input"] is True
    # Results must be identical to the equal-valued builtin-int input --
    # i.e. numpy integer scalars are genuinely normalized, not merely
    # tolerated as a different code path.
    assert r_numpy == r_plain


def test_numpy_and_builtin_int_of_equal_value_are_deduplicated_together():
    np = pytest.importorskip("numpy")
    r = check_embedded_cubical_normality(
        [(0, 0, 0), (np.int64(0), np.int64(0), np.int64(0)), (1, 0, 0)])
    assert r["valid_input"] is True
    # The numpy-int voxel and the equal-valued builtin-int voxel must be
    # recognized as the same coordinate and deduplicated, not kept as
    # two distinct top cells.
    assert r["duplicates_removed"] == 1


def test_numpy_float_coordinate_is_rejected_even_if_integral_valued():
    np = pytest.importorskip("numpy")
    r = check_embedded_cubical_normality(
        [(np.float64(0.0), np.float64(0.0), np.float64(0.0))])
    assert r["valid_input"] is False


# ---------------------------------------------------------------------
# 2. Facet-degree structural bound (v1 backend claim)
# ---------------------------------------------------------------------

def test_facet_parent_count_never_exceeds_two_structurally():
    """
    A dim=2 facet has exactly one even coordinate; the only two cells
    that can contain it as a facet are obtained by moving that axis to
    +/-1 while the other two (already odd) coordinates match exactly.
    So no input can ever produce a facet with 3+ parents in this
    representation -- checked here across many random configurations,
    not just asserted.
    """
    random.seed(0)
    for _ in range(200):
        n = random.randint(1, 15)
        voxels = list({
            (random.randint(0, 3), random.randint(0, 3), random.randint(0, 3))
            for _ in range(n)
        })
        top_cells = [voxel_to_cell(v) for v in voxels]
        fp = build_facet_parents(top_cells)
        for parents in fp.values():
            assert len(parents) <= 2


def test_facet_degree_status_is_always_pass_in_v1():
    random.seed(1)
    for _ in range(50):
        n = random.randint(1, 10)
        voxels = list({
            (random.randint(0, 3), random.randint(0, 3), random.randint(0, 3))
            for _ in range(n)
        })
        r = check_embedded_cubical_normality(voxels)
        assert r["facet_degree"]["status"] == "PASS"
        assert r["facet_degree"]["failures"] == []


# ---------------------------------------------------------------------
# 3. Edge-quadrant family
#
# Four voxels (0,0,k), (1,0,k), (0,1,k), (1,1,k) all touch the shared
# vertical edge at doubled coordinate (2, 2, 2k+1). Subsets of these
# four voxels reconstruct, at that shared edge, either: two disjoint
# simplicial edges (opposite quadrants only), a path (3 quadrants), or
# a cycle (4 quadrants) in the link -- exactly Definition 6's
# pseudomanifold connectivity/degree clause, one dimension down.
# ---------------------------------------------------------------------

Q00, Q10, Q01, Q11 = (0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0)
SHARED_EDGE = (2, 2, 1)
SHARED_EDGE_ENDPOINTS = {(2, 2, 0), (2, 2, 2)}


def test_far_apart_voxels_fail_connectivity_but_links_pass():
    r = check_embedded_cubical_normality([(0, 0, 0), (10, 10, 10)])
    assert r["facet_connectivity"]["status"] == "FAIL"
    assert len(r["facet_connectivity"]["components"]) == 2
    # no shared cells at all -> each local link is independently fine
    assert r["links"]["status"] == "PASS"
    assert r["overall"] == "FAIL"


def test_vertex_only_sharing_fails_connectivity_and_the_shared_vertex_link():
    # (0,0,0) and (1,1,1) touch only at the single corner vertex (2,2,2)
    r = check_embedded_cubical_normality([(0, 0, 0), (1, 1, 1)])
    assert r["facet_connectivity"]["status"] == "FAIL"
    assert len(r["facet_connectivity"]["components"]) == 2
    assert r["links"]["status"] == "FAIL"
    failed_faces = {h for h, _reason, _witness in r["links"]["failures"]}
    assert (2, 2, 2) in failed_faces  # the shared corner vertex


def test_edge_two_adjacent_quadrants_pass():
    # Q00, Q10 share a full facet (already share the x=2 plane) -- this
    # is the "adjacent quadrants" case, not touching only along the edge.
    r = check_embedded_cubical_normality([Q00, Q10])
    assert r["facet_connectivity"]["status"] == "PASS"
    assert r["links"]["status"] == "PASS"
    assert r["overall"] == "PASS"


def test_edge_two_opposite_quadrants_is_a_pinch():
    # Q00, Q11 touch only along the shared edge (2,2,z) -- diagonal pinch.
    r = check_embedded_cubical_normality([Q00, Q11])
    assert r["facet_connectivity"]["status"] == "FAIL"
    assert len(r["facet_connectivity"]["components"]) == 2
    assert r["links"]["status"] == "FAIL"
    failed_faces = {h for h, _reason, _witness in r["links"]["failures"]}
    assert SHARED_EDGE in failed_faces
    assert SHARED_EDGE_ENDPOINTS <= failed_faces
    # the failure mode at the shared edge/endpoints must specifically be
    # a connectivity split (two disjoint simplicial edges in the link),
    # not a degree violation
    reasons = {reason for h, reason, _w in r["links"]["failures"]
               if h == SHARED_EDGE}
    assert reasons == {"ridge_connectivity"}


def test_edge_three_quadrants_is_a_path_and_passes():
    r = check_embedded_cubical_normality([Q00, Q10, Q01])
    assert r["facet_connectivity"]["status"] == "PASS"
    assert r["links"]["status"] == "PASS"
    assert r["overall"] == "PASS"


def test_edge_four_quadrants_is_a_cycle_and_passes():
    r = check_embedded_cubical_normality([Q00, Q10, Q01, Q11])
    assert r["facet_connectivity"]["status"] == "PASS"
    assert r["links"]["status"] == "PASS"
    assert r["overall"] == "PASS"


# ---------------------------------------------------------------------
# 5. Larger hand-verified cases
#
# These re-confirm the small edge-quadrant results simultaneously inside
# one real 3D chunk (2x2x2), and probe a case where the *correct*
# mathematical expectation is PASS despite an apparent "hole" (3x3x3
# with the center voxel removed) -- checking that Pipeline 1 enforces
# only the local pseudomanifold conditions and does not secretly require
# global properties like a connected border or simple-connectivity.
# ---------------------------------------------------------------------

def _independent_codim_ge2_count(voxels):
    """
    Recount codim>=2 faces (dim <= 1) by direct enumeration over the
    union of closures, independently of check_embedded_cubical_normality's
    own internal bookkeeping -- so a bug shared between the checker and
    its own counting logic cannot silently agree with itself.
    """
    from cubical_normality.checker import closure as _closure, dim as _dim
    all_cells = set()
    for v in voxels:
        all_cells |= _closure(voxel_to_cell(v))
    return sum(1 for c in all_cells if _dim(c) <= 1)


def test_2x2x2_solid_block_passes_with_all_local_link_types():
    voxels = [(i, j, k) for i in range(2) for j in range(2) for k in range(2)]
    result = check_embedded_cubical_normality(voxels)

    assert result["facet_degree"]["status"] == "PASS"
    assert result["facet_connectivity"]["status"] == "PASS"
    assert result["links"]["status"] == "PASS"
    assert result["overall"] == "PASS"

    assert len(result["facet_connectivity"]["components"]) == 1
    assert result["links"]["failures"] == []
    # 3^3 vertices + 3*(2*3*3) edges = 27 + 54 = 81, pinned as a regression
    # baseline and cross-checked against an independent recount.
    assert result["links"]["checked_faces"] == 81
    assert result["links"]["checked_faces"] == _independent_codim_ge2_count(voxels)


def test_3x3x3_center_removed_cavity_shell_passes():
    """
    A 1-voxel cavity strictly interior to a 3x3x3 block: the border
    Delta F(K) has two disjoint components (outer surface + inner
    cavity surface) after this removal, but Definition 6/10 pseudo-
    manifold-ness is a purely local condition and never requires the
    border to be connected (Boutry's own "smooth n-PCM" explicitly
    allows a *separated union* of surfaces as a border) -- so the
    mathematically correct expectation here is PASS, not FAIL.
    """
    voxels = [(i, j, k) for i in range(3) for j in range(3) for k in range(3)
              if (i, j, k) != (1, 1, 1)]
    result = check_embedded_cubical_normality(voxels)

    assert result["facet_degree"]["status"] == "PASS"
    assert result["facet_connectivity"]["status"] == "PASS"
    assert result["links"]["status"] == "PASS"
    assert result["overall"] == "PASS"

    assert len(result["facet_connectivity"]["components"]) == 1
    assert result["links"]["failures"] == []
    assert result["links"]["checked_faces"] == 208
    assert result["links"]["checked_faces"] == _independent_codim_ge2_count(voxels)


# ---------------------------------------------------------------------
# 6. Invariance / equivariance
# ---------------------------------------------------------------------

BASE_CASES = [
    [(0, 0, 0)],
    [Q00, Q10],
    [Q00, Q11],
    [Q00, Q10, Q01],
    [Q00, Q10, Q01, Q11],
    [(0, 0, 0), (1, 1, 1)],
]


def _summarize(r):
    """Order/labeling-independent summary of a result, for comparison."""
    return (
        r["overall"],
        r["facet_degree"]["status"],
        r["facet_connectivity"]["status"],
        len(r["facet_connectivity"]["components"]),
        r["links"]["status"],
        r["links"]["checked_faces"],
        len(r["links"]["failures"]),
    )


@pytest.mark.parametrize("voxels", BASE_CASES)
def test_result_invariant_under_input_order_permutation(voxels):
    baseline = _summarize(check_embedded_cubical_normality(voxels))
    random.seed(42)
    for _ in range(5):
        shuffled = list(voxels)
        random.shuffle(shuffled)
        assert _summarize(check_embedded_cubical_normality(shuffled)) == baseline


@pytest.mark.parametrize("voxels", BASE_CASES)
@pytest.mark.parametrize("shift", [(1, 0, 0), (0, 5, 0), (3, -2, 7)])
def test_result_invariant_under_translation(voxels, shift):
    baseline = _summarize(check_embedded_cubical_normality(voxels))
    shifted = [(v[0] + shift[0], v[1] + shift[1], v[2] + shift[2]) for v in voxels]
    assert _summarize(check_embedded_cubical_normality(shifted)) == baseline


@pytest.mark.parametrize("voxels", BASE_CASES)
@pytest.mark.parametrize("perm", list(itertools.permutations(range(3))))
def test_result_invariant_under_axis_permutation(voxels, perm):
    baseline = _summarize(check_embedded_cubical_normality(voxels))
    permuted = [(v[perm[0]], v[perm[1]], v[perm[2]]) for v in voxels]
    assert _summarize(check_embedded_cubical_normality(permuted)) == baseline


# ---------------------------------------------------------------------
# 5. Certificate layer (packaging session 37, §5): schema, determinism,
#    JSON round-trip, invalid-input unification.
# ---------------------------------------------------------------------

def test_certificate_valid_input_schema_and_counts():
    cert = build_certificate([(0, 0, 0), (1, 0, 0)])
    assert cert["schema_version"] == "1.0"
    assert cert["backend"] == "embedded_doubled_coordinate_3d"
    assert cert["valid_input"] is True
    assert cert["overall"] == "PASS"
    assert cert["input"] == {
        "received_count": 2, "normalized_count": 2, "duplicates_removed": 0,
    }
    assert cert["facet_degree"]["checked_facets"] > 0
    assert cert["facet_connectivity"]["component_count"] == 1


def test_certificate_checked_faces_matches_hand_verified_cases():
    block = list(itertools.product(range(2), range(2), range(2)))
    cube = list(itertools.product(range(3), range(3), range(3)))
    cavity = [v for v in cube if v != (1, 1, 1)]
    assert build_certificate(block)["links"]["checked_faces"] == 81
    assert build_certificate(cavity)["links"]["checked_faces"] == 208


def test_certificate_invalid_input_uses_unified_schema_not_raw_error_string():
    cert = build_certificate([(0, 0, 0), (1, 0)])
    assert cert["valid_input"] is False
    assert cert["overall"] == "INVALID_INPUT"  # distinct from FAIL
    assert cert["validation"]["status"] == "FAIL"
    assert cert["validation"]["failures"][0]["code"] == "NOT_LENGTH_3_TUPLE"
    assert cert["purity"]["status"] == "NOT_RUN"
    assert cert["links"]["status"] == "NOT_RUN"


def test_certificate_invalid_input_collects_all_failures_not_just_first():
    cert = build_certificate([(0, 0), (True, 0, 0), (1.5, 0, 0)])
    codes = [f["code"] for f in cert["validation"]["failures"]]
    assert "NOT_LENGTH_3_TUPLE" in codes
    assert "BOOLEAN_COORDINATE" in codes
    assert "NON_INTEGRAL_COORDINATE" in codes
    assert len(codes) == 3  # all three collected, not just the first


def test_certificate_is_fully_json_serializable_no_tuples_or_sets():
    cert = build_certificate([(0, 0, 0), (1, 0, 0), (2, 0, 0)])
    roundtripped = json.loads(json.dumps(cert))
    assert roundtripped == cert


def test_certificate_deterministic_under_input_order_permutation():
    voxels = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)]
    baseline = build_certificate(voxels)
    random.seed(7)
    for _ in range(5):
        shuffled = list(voxels)
        random.shuffle(shuffled)
        assert build_certificate(shuffled) == baseline


def test_certificate_link_failures_canonically_sorted_across_multiple_witnesses():
    # Vertex-only-sharing case produces several independent link
    # failures (several disconnected base faces). This exercises the
    # outer-list sort of `link_failures` (by base_face, code) and not
    # just single-witness cases.
    voxels = [(0, 0, 0), (1, 1, 0)]
    baseline = build_certificate(voxels)
    assert len(baseline["links"]["failures"]) > 1
    base_faces = [f["base_face"] for f in baseline["links"]["failures"]]
    assert base_faces == sorted(base_faces)  # outer list is canonically ordered

    random.seed(11)
    for _ in range(5):
        shuffled = list(voxels)
        random.shuffle(shuffled)
        assert build_certificate(shuffled) == baseline


def test_certificate_cli_wrapper_passes_through_on_success():
    voxels = [(0, 0, 0)]
    assert build_certificate_cli(voxels) == build_certificate(voxels)


def test_certificate_cli_wrapper_reports_internal_error_not_traceback():
    # Force an unexpected internal exception (not a structural input
    # problem) and confirm the CLI wrapper reports it structurally
    # instead of propagating a raw traceback.
    class ExplodesOnLen:
        def __iter__(self):
            raise RuntimeError("simulated internal failure")

    cert = build_certificate_cli(ExplodesOnLen())
    assert cert["overall"] == "INTERNAL_ERROR"
    assert cert["internal_error"]["type"] == "RuntimeError"


def test_regression_original_contract_unaffected_by_certificate_layer():
    # Sanity: the raw checker's contract (tested by the other tests in
    # this file) is untouched by adding build_certificate on top.
    r = check_embedded_cubical_normality([(0, 0, 0)])
    assert set(r.keys()) == {
        "valid_input", "duplicates_removed", "purity", "facet_degree",
        "facet_connectivity", "links", "overall",
    }
