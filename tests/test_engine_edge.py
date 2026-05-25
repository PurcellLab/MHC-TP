"""Edge-case tests for the numba search engine.

Covers kernels.compute_all_correlations, cache.build_reference_array, and
search.search. The first numba-compiled call is slow (JIT compile); that is
expected.
"""

import numpy as np
import pandas as pd
import pytest

from mhc_tp.constants import N_AMINO_ACIDS
from mhc_tp.engine.cache import build_reference_array
from mhc_tp.engine.kernels import compute_all_correlations
from mhc_tp.engine.search import search


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _ramp_matrix(p=1):
    """(p, 20) ramp matrix; cell 0 is 0.0 so masked out, leaving 20*p-1 cells."""
    return np.arange(p * N_AMINO_ACIDS, dtype=np.float32).reshape(p, N_AMINO_ACIDS)


def _ref_df(formatted, matrices, n_positions=None):
    """Build a reference DataFrame from a list of (p,20) matrices."""
    if n_positions is None:
        n_positions = [m.shape[0] for m in matrices]
    return pd.DataFrame(
        {
            "allotype": [f"A*{i:02d}:01" for i in range(len(formatted))],
            "formatted": list(formatted),
            "mhc_class": ["I"] * len(formatted),
            "locus": ["A"] * len(formatted),
            "n_positions": list(n_positions),
            "matrix": [
                np.asarray(m, dtype=np.float32).ravel().tolist() for m in matrices
            ],
            "source": ["test"] * len(formatted),
        }
    )


# ===========================================================================
# kernels.compute_all_correlations
# ===========================================================================
def test_kernel_identical_corr_near_one():
    g = _ramp_matrix(1).reshape(1, 1, N_AMINO_ACIDS)
    refs = _ramp_matrix(1).reshape(1, 1, N_AMINO_ACIDS)
    mask = np.ones(1, dtype=np.bool_)
    corr, invalid = compute_all_correlations(g, refs, mask, 0.5)
    assert invalid[0] == 0
    assert abs(corr[0, 0] - 1.0) < 1e-5


def test_kernel_inverted_below_threshold_stays_sentinel():
    # inverted ramp -> corr ~ -1, which is below any positive threshold -> -1.0 sentinel
    g = _ramp_matrix(1).reshape(1, 1, N_AMINO_ACIDS)
    inv = _ramp_matrix(1)[:, ::-1].copy().reshape(1, 1, N_AMINO_ACIDS)
    mask = np.ones(1, dtype=np.bool_)
    corr, invalid = compute_all_correlations(g, inv, mask, 0.5)
    assert invalid[0] == 0
    assert corr[0, 0] == np.float32(-1.0)


def test_kernel_inverted_with_low_threshold_stores_negative_corr():
    # With threshold very low (-2.0), a ~-1 correlation IS >= threshold and is stored.
    g = _ramp_matrix(1).reshape(1, 1, N_AMINO_ACIDS)
    inv = _ramp_matrix(1)[:, ::-1].copy().reshape(1, 1, N_AMINO_ACIDS)
    mask = np.ones(1, dtype=np.bool_)
    corr, invalid = compute_all_correlations(g, inv, mask, -2.0)
    assert invalid[0] == 0
    assert abs(corr[0, 0] - (-1.0)) < 1e-5


def test_kernel_too_few_valid_cells_invalid():
    g = np.zeros((1, 1, N_AMINO_ACIDS), dtype=np.float32)
    g[0, 0, :9] = np.arange(1, 10, dtype=np.float32)  # 9 non-zero cells < 10
    refs = np.ones((1, 1, N_AMINO_ACIDS), dtype=np.float32)
    mask = np.ones(1, dtype=np.bool_)
    corr, invalid = compute_all_correlations(g, refs, mask, 0.5)
    assert invalid[0] == 1
    assert corr[0, 0] == np.float32(-1.0)


def test_kernel_exactly_ten_valid_cells_not_invalid():
    # 10 non-zero cells == boundary: needs >= 10, so this should NOT be flagged.
    g = np.zeros((1, 1, N_AMINO_ACIDS), dtype=np.float32)
    g[0, 0, :10] = np.arange(1, 11, dtype=np.float32)
    refs = np.zeros((1, 1, N_AMINO_ACIDS), dtype=np.float32)
    refs[0, 0, :10] = np.arange(1, 11, dtype=np.float32)
    mask = np.ones(1, dtype=np.bool_)
    corr, invalid = compute_all_correlations(g, refs, mask, 0.5)
    assert invalid[0] == 0
    assert abs(corr[0, 0] - 1.0) < 1e-5


def test_kernel_constant_gibbs_zero_std_invalid():
    # gibbs all 7.0 -> >=10 valid cells but zero std -> invalid
    g = np.full((1, 1, N_AMINO_ACIDS), 7.0, dtype=np.float32)
    refs = _ramp_matrix(1).reshape(1, 1, N_AMINO_ACIDS)
    mask = np.ones(1, dtype=np.bool_)
    corr, invalid = compute_all_correlations(g, refs, mask, -2.0)
    assert invalid[0] == 1


def test_kernel_constant_ref_zero_std_skipped_not_invalid():
    # constant ref -> r_std == 0 -> column skipped (stays -1) but gibbs NOT invalid
    g = _ramp_matrix(1).reshape(1, 1, N_AMINO_ACIDS)
    refs = np.full((1, 1, N_AMINO_ACIDS), 3.0, dtype=np.float32)
    mask = np.ones(1, dtype=np.bool_)
    corr, invalid = compute_all_correlations(g, refs, mask, -2.0)
    assert invalid[0] == 0
    assert corr[0, 0] == np.float32(-1.0)


def test_kernel_hla_mask_false_column_stays_sentinel():
    g = _ramp_matrix(1).reshape(1, 1, N_AMINO_ACIDS)
    refs = np.stack([_ramp_matrix(1), _ramp_matrix(1)]).astype(  # both identical to g
        np.float32
    )
    mask = np.array([False, True], dtype=np.bool_)
    corr, invalid = compute_all_correlations(g, refs, mask, 0.5)
    assert invalid[0] == 0
    assert corr[0, 0] == np.float32(-1.0)  # masked out
    assert abs(corr[0, 1] - 1.0) < 1e-5  # included


def test_kernel_nan_cells_handled():
    # NaNs in gibbs are excluded from the valid mask. Place NaNs and verify the
    # correlation matches the non-NaN/non-zero subset (identical to ref there).
    g = _ramp_matrix(1).reshape(1, 1, N_AMINO_ACIDS).copy()
    g[0, 0, 5] = np.nan
    g[0, 0, 11] = np.nan
    refs = _ramp_matrix(1).reshape(1, 1, N_AMINO_ACIDS)
    mask = np.ones(1, dtype=np.bool_)
    corr, invalid = compute_all_correlations(g, refs, mask, 0.5)
    assert invalid[0] == 0
    assert not np.isnan(corr[0, 0])
    assert abs(corr[0, 0] - 1.0) < 1e-5


def test_kernel_threshold_gating_just_below_and_above():
    # Build a gibbs/ref pair with a known intermediate correlation, then gate.
    rng = np.random.default_rng(0)
    base = rng.random(N_AMINO_ACIDS).astype(np.float32) + 0.1  # all non-zero
    noise = rng.random(N_AMINO_ACIDS).astype(np.float32) + 0.1
    g = base.reshape(1, 1, N_AMINO_ACIDS)
    ref = (0.5 * base + 0.5 * noise).reshape(1, 1, N_AMINO_ACIDS).astype(np.float32)
    mask = np.ones(1, dtype=np.bool_)
    # measure the true correlation with a permissive threshold
    c, _ = compute_all_correlations(g, ref, mask, -2.0)
    true_corr = float(c[0, 0])
    assert -1.0 < true_corr < 1.0  # genuinely intermediate
    # threshold just below -> stored
    c_below, _ = compute_all_correlations(g, ref, mask, true_corr - 0.01)
    assert abs(c_below[0, 0] - true_corr) < 1e-5
    # threshold just above -> sentinel
    c_above, _ = compute_all_correlations(g, ref, mask, true_corr + 0.01)
    assert c_above[0, 0] == np.float32(-1.0)


# ===========================================================================
# cache.build_reference_array
# ===========================================================================
def test_build_ref_different_positions_padded_to_max():
    m1 = _ramp_matrix(2)  # 2 positions
    m3 = _ramp_matrix(3)  # 3 positions
    df = _ref_df(["A", "B"], [m1, m3])
    arr, max_positions = build_reference_array(df)
    assert max_positions == 3
    assert arr.shape == (2, 3, N_AMINO_ACIDS)
    # row 0 (2 positions): position index 2 must be zero-padded
    assert np.all(arr[0, 2, :] == 0.0)
    # values preserved for the real rows
    assert np.allclose(arr[0, :2, :], m1)
    assert np.allclose(arr[1, :3, :], m3)


def test_build_ref_empty_df():
    df = _ref_df([], [])
    arr, max_positions = build_reference_array(df)
    assert max_positions == 0
    assert arr.shape == (0, 0, N_AMINO_ACIDS)


def test_build_ref_single_row_values_preserved():
    m = _ramp_matrix(2)
    df = _ref_df(["A"], [m])
    arr, max_positions = build_reference_array(df)
    assert max_positions == 2
    assert arr.shape == (1, 2, N_AMINO_ACIDS)
    assert np.allclose(arr[0], m)


def test_build_ref_reshape_correct_layout():
    # Flat values 0..(2*20-1) should reshape row-major into (2, 20).
    m = np.arange(2 * N_AMINO_ACIDS, dtype=np.float32).reshape(2, N_AMINO_ACIDS)
    df = _ref_df(["A"], [m])
    arr, _ = build_reference_array(df)
    assert arr[0, 0, 0] == 0.0
    assert arr[0, 0, N_AMINO_ACIDS - 1] == np.float32(N_AMINO_ACIDS - 1)
    assert arr[0, 1, 0] == np.float32(N_AMINO_ACIDS)


# ===========================================================================
# search.search
# ===========================================================================
def test_search_empty_gibbs_dict():
    df = _ref_df(["A0201"], [_ramp_matrix(1)])
    assert search(df, {}, threshold=0.5, top_n=3) == {}


def test_search_no_matches_above_threshold():
    df = _ref_df(["A0201"], [_ramp_matrix(1)])
    # gibbs is the inverted ramp -> corr ~ -1, well below 0.9
    gibbs = {"g.mat": _ramp_matrix(1)[:, ::-1].copy()}
    assert search(df, gibbs, threshold=0.9, top_n=3) == {}


def test_search_finds_best_highest_corr():
    ramp = _ramp_matrix(1)
    inv = ramp[:, ::-1].copy()
    df = _ref_df(["A0201", "A0101"], [ramp, inv])
    gibbs = {"1of1.mat": ramp.copy()}  # matches A0201 exactly
    res = search(df, gibbs, threshold=0.5, top_n=2)
    assert res  # non-empty
    best = max(res.items(), key=lambda kv: kv[1])
    assert best[0] == ("1of1.mat", "A0201")
    assert abs(best[1] - 1.0) < 1e-5


def test_search_hla_filter_restricts_results():
    ramp = _ramp_matrix(1)
    df = _ref_df(["A0201", "A0101"], [ramp.copy(), ramp.copy()])  # both match
    gibbs = {"g.mat": ramp.copy()}
    res = search(df, gibbs, threshold=0.5, top_n=5, hla_filter=["A0101"])
    assert set(k[1] for k in res) == {"A0101"}


def test_search_top_n_caps_hits():
    ramp = _ramp_matrix(1)
    # four refs all identical to the gibbs -> all corr ~1; top_n=2 should cap.
    df = _ref_df(["R0", "R1", "R2", "R3"], [ramp.copy() for _ in range(4)])
    gibbs = {"g.mat": ramp.copy()}
    res = search(df, gibbs, threshold=0.5, top_n=2)
    hits_for_g = [k for k in res if k[0] == "g.mat"]
    assert len(hits_for_g) == 2


def test_search_gibbs_shorter_than_reference_padding():
    # reference has 2 positions, gibbs has 1 -> gibbs padded with zeros (masked out).
    ref2 = _ramp_matrix(2)
    df = _ref_df(["A0201"], [ref2])
    gibbs = {"g.mat": _ramp_matrix(1)}  # 1 position
    res = search(df, gibbs, threshold=-2.0, top_n=3)
    # gibbs has 19 valid cells (>=10); correlation computed against ref's first row
    assert ("g.mat", "A0201") in res


def test_search_gibbs_longer_than_reference_raises():
    # gibbs has MORE positions than any reference -> padded array sized to ref
    # max_positions, so the assignment overflows. Documents real behavior.
    ref1 = _ramp_matrix(1)
    df = _ref_df(["A0201"], [ref1])
    gibbs = {"g.mat": _ramp_matrix(3)}  # 3 positions > ref max 1
    with pytest.raises(ValueError):
        search(df, gibbs, threshold=0.5, top_n=3)
