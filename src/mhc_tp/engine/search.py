"""Search Gibbs cluster matrices against the reference array.

Method (how a cluster is matched to an allotype)
------------------------------------------------
1. Each GibbsCluster motif and each reference allotype is represented as a
   position-specific scoring matrix (PSSM): ``n_positions x 20`` amino-acid
   weights. Gibbs matrices are zero-padded to a common position count so they
   can be batched.
2. For every (cluster, allotype) pair the two PSSMs are compared by **Pearson
   correlation** of their flattened weights, computed only over the cells
   ``V = {k : g_k != 0 and g_k not NaN}`` that are informative in the *cluster*
   matrix (so padding and empty positions do not dilute the score)::

       PCC(g, r) = Σ_{k∈V}(g_k - ḡ)(r_k - r̄) / ( |V| · σ_g · σ_r )   ∈ [-1, 1]

   with means/std taken over ``V``. It is scale- and offset-invariant, so it
   scores motif *shape*, not absolute magnitudes. Full derivation and numerical
   guards: :func:`mhc_tp.engine.kernels.compute_all_correlations`.
3. Per cluster the allotypes are ranked by correlation (PCC, ``-1..1``;
   ``1.0`` = identical motif). Selection is then either threshold-gated
   (default) or pure top-N (``always_top_n``); see :func:`search`.

The correlation is a *motif-shape* similarity: it rewards matching the relative
preference pattern across positions, not the absolute weight magnitudes.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from mhc_tp.constants import N_AMINO_ACIDS
from mhc_tp.engine.cache import build_reference_array
from mhc_tp.engine.kernels import compute_all_correlations


def search(
    reference: pd.DataFrame,
    gibbs_matrices: dict[str, np.ndarray],
    threshold: float = 0.70,
    top_n: int = 3,
    hla_filter: list[str] | None = None,
    always_top_n: bool = False,
) -> dict[tuple[str, str], float]:
    """Return ``{(gibbs_name, ref_formatted): correlation}`` for top-N hits.

    By default a hit must score ``>= threshold`` to be returned, so a cluster
    may yield fewer than ``top_n`` rows (or none). When ``always_top_n`` is set,
    every cluster returns its ``top_n`` best matches regardless of threshold —
    the threshold then only annotates confidence downstream, it never drops a
    row.
    """
    ref_arr, max_positions = build_reference_array(reference)
    names = list(gibbs_matrices.keys())

    padded = np.zeros((len(names), max_positions, N_AMINO_ACIDS), dtype=np.float32)
    for i, name in enumerate(names):
        m = gibbs_matrices[name]
        padded[i, : m.shape[0], : m.shape[1]] = m

    mask = np.ones(len(reference), dtype=np.bool_)
    if hla_filter:
        mask = reference["formatted"].isin(hla_filter).to_numpy()

    # In always-top-N mode, store every valid correlation (kernel keeps a -1.0
    # sentinel for cells below its threshold), then rank in Python.
    kernel_threshold = -2.0 if always_top_n else threshold
    corr, _invalid = compute_all_correlations(
        padded, ref_arr.astype(np.float32), mask, kernel_threshold
    )

    formatted = reference["formatted"].to_numpy()
    out: dict[tuple[str, str], float] = {}
    for i, name in enumerate(names):
        row = corr[i, :]
        order = np.argsort(row)[::-1]
        if always_top_n:
            # Top-N among computed (non-sentinel, unmasked) cells, any score.
            hits = [j for j in order if mask[j] and row[j] > -1.0][:top_n]
        else:
            hits = [j for j in order if row[j] >= threshold][:top_n]
        for j in hits:
            out[(name, str(formatted[j]))] = float(row[j])
    return out
