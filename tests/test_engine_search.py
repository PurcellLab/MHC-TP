import numpy as np
import pandas as pd
from mhc_tp.engine.search import search


def _ref_df():
    ramp = np.arange(20, dtype=np.float32)  # one position, 20 AA values
    inv = ramp[::-1].copy()
    return pd.DataFrame(
        {
            "allotype": ["A*02:01", "A*01:01"],
            "formatted": ["A0201", "A0101"],
            "mhc_class": ["I", "I"],
            "locus": ["A", "A"],
            "n_positions": [1, 1],
            "matrix": [ramp.tolist(), inv.tolist()],
            "source": ["test", "test"],
        }
    )


def test_search_finds_best():
    ref = _ref_df()
    gibbs = {"1of1.mat": np.arange(20, dtype=np.float32).reshape(1, 20)}  # == A0201
    cd = search(ref, gibbs, threshold=0.5, top_n=2)
    best = max(cd.items(), key=lambda kv: kv[1])
    assert best[0][0] == "1of1.mat"
    assert best[0][1] == "A0201"
    assert abs(best[1] - 1.0) < 1e-5


def test_search_threshold_gated_drops_below():
    # A0101 is the anti-correlated motif (PCC ~ -1) → excluded by the threshold.
    ref = _ref_df()
    gibbs = {"1of1.mat": np.arange(20, dtype=np.float32).reshape(1, 20)}
    cd = search(ref, gibbs, threshold=0.5, top_n=2)
    assert {k[1] for k in cd} == {"A0201"}


def _ref_df_with_intermediate():
    """A0201 (PCC≈1) plus A8888, a 'tent' motif with a weak (sub-0.5) PCC."""
    ramp = np.arange(20, dtype=np.float32)
    tent = np.concatenate([np.arange(10), np.arange(10)[::-1]]).astype(np.float32)
    return pd.DataFrame(
        {
            "allotype": ["A*02:01", "A*88:88"],
            "formatted": ["A0201", "A8888"],
            "mhc_class": ["I", "I"],
            "locus": ["A", "A"],
            "n_positions": [1, 1],
            "matrix": [ramp.tolist(), tent.tolist()],
            "source": ["test", "test"],
        }
    )


def test_search_always_top_n_keeps_below_threshold():
    # always_top_n returns the top-N per cluster regardless of threshold, so a
    # genuine but weak match (A8888, PCC < 0.5) is still returned, ranked last.
    ref = _ref_df_with_intermediate()
    gibbs = {"1of1.mat": np.arange(20, dtype=np.float32).reshape(1, 20)}

    measured = search(ref, gibbs, threshold=-1.0, top_n=2, always_top_n=True)
    weak = measured[("1of1.mat", "A8888")]
    assert -1.0 < weak < 0.5  # genuinely intermediate, not the -1.0 sentinel

    gated = search(ref, gibbs, threshold=0.5, top_n=2)
    assert {k[1] for k in gated} == {"A0201"}  # weak match dropped

    full = search(ref, gibbs, threshold=0.5, top_n=2, always_top_n=True)
    assert {k[1] for k in full} == {"A0201", "A8888"}  # weak match kept
    assert full[("1of1.mat", "A0201")] > full[("1of1.mat", "A8888")]


def test_search_always_top_n_caps_at_top_n():
    # With more refs than top_n, always_top_n still returns exactly top_n best.
    ref = _ref_df_with_intermediate()
    gibbs = {"1of1.mat": np.arange(20, dtype=np.float32).reshape(1, 20)}
    one = search(ref, gibbs, threshold=0.9, top_n=1, always_top_n=True)
    assert len(one) == 1
    assert next(iter(one))[1] == "A0201"  # the single best


def test_search_always_top_n_backfills_to_top_n():
    # Only A0201 clears a high threshold (default → 1 row); always_top_n
    # backfills the weak A8888 to reach the requested top_n=2.
    ref = _ref_df_with_intermediate()
    gibbs = {"1of1.mat": np.arange(20, dtype=np.float32).reshape(1, 20)}
    gated = search(ref, gibbs, threshold=0.9, top_n=2)
    assert {k[1] for k in gated} == {"A0201"}
    full = search(ref, gibbs, threshold=0.9, top_n=2, always_top_n=True)
    assert {k[1] for k in full} == {"A0201", "A8888"}


def test_search_always_top_n_respects_hla_filter():
    ref = _ref_df_with_intermediate()
    gibbs = {"1of1.mat": np.arange(20, dtype=np.float32).reshape(1, 20)}
    full = search(
        ref, gibbs, threshold=0.5, top_n=2, always_top_n=True, hla_filter=["A8888"]
    )
    assert {k[1] for k in full} == {"A8888"}  # only the filtered allotype considered
