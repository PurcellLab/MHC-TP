import numpy as np
from mhc_tp.engine.kernels import compute_all_correlations


def test_correlation_identical_and_inverted():
    # one gibbs matrix shape (1, 1, 20): a ramp 0..19 (index 0 == 0.0 is excluded
    # by the non-zero mask, leaving 19 valid cells >= 10).
    g = np.arange(20, dtype=np.float32).reshape(1, 1, 20)
    refs = np.stack(
        [
            np.arange(20, dtype=np.float32).reshape(1, 20),  # identical -> corr ~1
            np.arange(20, dtype=np.float32)[::-1].reshape(
                1, 20
            ),  # inverted  -> corr ~-1 (below thr)
        ]
    ).astype(np.float32)
    mask = np.ones(2, dtype=np.bool_)
    corr, invalid = compute_all_correlations(g, refs, mask, 0.5)
    assert invalid[0] == 0
    assert abs(corr[0, 0] - 1.0) < 1e-5
    assert corr[0, 1] == np.float32(-1.0)  # below-threshold sentinel preserved


def test_too_few_valid_cells_flagged_invalid():
    g = np.zeros((1, 1, 20), dtype=np.float32)
    g[0, 0, 0] = 5.0  # only 1 non-zero cell -> < 10 -> invalid
    refs = np.ones((1, 1, 20), dtype=np.float32)
    mask = np.ones(1, dtype=np.bool_)
    corr, invalid = compute_all_correlations(g, refs, mask, 0.5)
    assert invalid[0] == 1
