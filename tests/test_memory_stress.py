"""Memory-leak and high-memory-stress tests for the search pipeline.

Run with:  pixi run pytest tests/test_memory_stress.py -v

Measurement uses STDLIB only:
  * ``tracemalloc`` for Python-object allocation growth (cross-platform).
  * ``resource.getrusage(RUSAGE_SELF).ru_maxrss`` for peak RSS.

``resource`` is Unix-only, so it is guarded.  On Linux ``ru_maxrss`` is in
kilobytes; on macOS it is in bytes -- we normalise to bytes below.

Synthetic data only (no reading of data/ref_data/*.parquet).  Sizes are kept
modest so the whole module stays well under ~60s once numba has JIT-compiled.
"""

from __future__ import annotations

import gc
import sys
import tracemalloc

import numpy as np
import pandas as pd
import pytest

from hla_pepclust.constants import N_AMINO_ACIDS
from hla_pepclust.engine.cache import build_reference_array
from hla_pepclust.engine.search import search

# resource is Unix-only; skip the whole module on Windows.
if sys.platform.startswith("win"):
    pytest.skip("resource module unavailable on Windows", allow_module_level=True)

resource = pytest.importorskip("resource")


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _peak_rss_bytes() -> int:
    """Peak resident-set size of this process, normalised to bytes.

    Linux reports ``ru_maxrss`` in kilobytes, macOS/BSD in bytes.
    """
    maxrss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        return int(maxrss)
    return int(maxrss) * 1024


def make_ref(n_alleles: int, n_positions: int = 9, seed: int = 0) -> pd.DataFrame:
    """Build a synthetic reference DataFrame with the real schema.

    Columns: allotype, formatted, mhc_class, locus, n_positions, matrix, source.
    ``matrix`` holds a flattened (n_positions * 20) list of random float32.
    """
    rng = np.random.default_rng(seed)
    matrices = [
        rng.standard_normal(n_positions * N_AMINO_ACIDS).astype(np.float32).tolist()
        for _ in range(n_alleles)
    ]
    return pd.DataFrame(
        {
            "allotype": [f"A*{i:04d}" for i in range(n_alleles)],
            "formatted": [f"A{i:06d}" for i in range(n_alleles)],
            "mhc_class": ["I"] * n_alleles,
            "locus": ["A"] * n_alleles,
            "n_positions": [n_positions] * n_alleles,
            "matrix": matrices,
            "source": ["synthetic"] * n_alleles,
        }
    )


def make_gibbs(
    n_matrices: int, n_positions: int = 9, seed: int = 1
) -> dict[str, np.ndarray]:
    """Build a dict of synthetic Gibbs matrices shaped (n_positions, 20)."""
    rng = np.random.default_rng(seed)
    return {
        f"{k + 1}of{n_matrices}.mat": rng.standard_normal(
            (n_positions, N_AMINO_ACIDS)
        ).astype(np.float32)
        for k in range(n_matrices)
    }


# --------------------------------------------------------------------------- #
# 1. Leak test -- repeated search must not grow Python-allocated memory.
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_repeated_search_no_leak() -> None:
    ref = make_ref(500, n_positions=9, seed=42)
    gibbs = make_gibbs(4, n_positions=9, seed=43)

    warmup = 3
    total = 20
    growth_bound_mb = 50.0

    tracemalloc.start()
    snapshot_before = None
    peak_after = 0

    for it in range(total):
        result = search(ref, gibbs, threshold=0.5, top_n=3)
        assert isinstance(result, dict)
        del result
        gc.collect()

        if it == warmup - 1:
            # Baseline taken after warmup (lets caches / numba settle).
            tracemalloc.reset_peak()
            snapshot_before = tracemalloc.take_snapshot()
        elif it >= warmup:
            _cur, peak = tracemalloc.get_traced_memory()
            peak_after = max(peak_after, peak)

    snapshot_after = tracemalloc.take_snapshot()
    assert snapshot_before is not None

    stats = snapshot_after.compare_to(snapshot_before, "filename")
    extra_bytes = sum(s.size_diff for s in stats if s.size_diff > 0)
    extra_mb = extra_bytes / 1e6
    peak_mb = peak_after / 1e6

    tracemalloc.stop()

    print(
        f"\n[leak] post-warmup peak={peak_mb:.1f} MB  "
        f"net positive growth={extra_mb:.1f} MB across {total - warmup} iters"
    )

    leak_detected = extra_mb >= growth_bound_mb
    assert not leak_detected, (
        f"Possible per-call leak: Python allocations grew {extra_mb:.1f} MB "
        f"(>= {growth_bound_mb} MB) across {total - warmup} post-warmup searches"
    )

    del ref, gibbs
    gc.collect()


# --------------------------------------------------------------------------- #
# 2. High-memory stress -- large reference, bounded peak RSS.
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_large_reference_search_bounded_rss() -> None:
    # Human scale: 20k alleles x 9 positions x 20 AA.  Drop to 10k if too slow.
    n_alleles = 20_000
    n_positions = 9
    rss_bound_bytes = 4 * 1024**3  # 4 GB

    ref = make_ref(n_alleles, n_positions=n_positions, seed=7)
    gibbs = make_gibbs(5, n_positions=n_positions, seed=8)

    result = search(ref, gibbs, threshold=0.7, top_n=3)
    assert isinstance(result, dict)
    # Keys are (gibbs_name, ref_formatted); values float correlations.
    for (gname, fmt), corr in result.items():
        assert isinstance(gname, str)
        assert isinstance(fmt, str)
        assert isinstance(corr, float)

    peak = _peak_rss_bytes()
    peak_gb = peak / 1024**3
    print(
        f"\n[stress] {n_alleles} alleles x {n_positions} pos -> "
        f"{len(result)} hits, peak RSS={peak_gb:.2f} GB"
    )

    blowup_detected = peak >= rss_bound_bytes
    assert not blowup_detected, (
        f"Peak RSS {peak_gb:.2f} GB exceeded {rss_bound_bytes / 1024**3:.0f} GB "
        f"bound -- possible O(n^2) memory blowup"
    )

    # Cleanup so the stress test does not bloat the rest of the suite.
    del ref, gibbs, result
    gc.collect()


# --------------------------------------------------------------------------- #
# 3. build_reference_array -- correct dtype/shape and expected memory order.
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_build_reference_array_memory() -> None:
    n_alleles = 20_000
    n_positions = 9

    ref = make_ref(n_alleles, n_positions=n_positions, seed=11)

    arr, max_positions = build_reference_array(ref)

    assert arr.dtype == np.float32
    assert max_positions == n_positions
    assert arr.shape == (n_alleles, n_positions, N_AMINO_ACIDS)

    expected_bytes = n_alleles * n_positions * N_AMINO_ACIDS * 4
    assert arr.nbytes == expected_bytes
    # Order-of-magnitude sanity: ~14.4 MB for 20k x 9 x 20 x 4 bytes.
    assert arr.nbytes < 1024**3  # comfortably under 1 GB

    print(
        f"\n[build] array shape={arr.shape} dtype={arr.dtype} "
        f"nbytes={arr.nbytes / 1e6:.1f} MB"
    )

    del ref, arr
    gc.collect()
