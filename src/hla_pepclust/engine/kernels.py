"""Numba JIT correlation kernel (ported from the proven NumbaSearch engine)."""

import numpy as np
from numba import jit, prange


@jit(nopython=True, parallel=True, fastmath=True)
def compute_all_correlations(gibbs_matrices, ref_matrices, hla_mask, threshold):
    """All-pairs flattened Pearson correlation, parallel over Gibbs matrices.

    Correlates over the cells that are non-zero and non-NaN in the Gibbs
    matrix. Stores a correlation only when it is >= threshold; otherwise the
    cell keeps the -1.0 sentinel. Returns (correlations, invalid_flags).
    """
    n_gibbs = gibbs_matrices.shape[0]
    n_refs = ref_matrices.shape[0]
    correlations = np.full((n_gibbs, n_refs), -1.0, dtype=np.float32)
    invalid_flags = np.zeros(n_gibbs, dtype=np.int32)

    for i in prange(n_gibbs):
        gibbs_flat = gibbs_matrices[i].flatten()
        valid = ~(np.isnan(gibbs_flat) | (gibbs_flat == 0.0))
        gibbs_clean = gibbs_flat[valid]
        if len(gibbs_clean) < 10:
            invalid_flags[i] = 1
            continue
        g_mean = np.mean(gibbs_clean)
        g_std = np.std(gibbs_clean)
        if g_std == 0.0:
            invalid_flags[i] = 1
            continue
        for j in range(n_refs):
            if not hla_mask[j]:
                continue
            ref_clean = ref_matrices[j].flatten()[valid]
            r_mean = np.mean(ref_clean)
            r_std = np.std(ref_clean)
            if r_std == 0.0:
                continue
            num = np.mean((gibbs_clean - g_mean) * (ref_clean - r_mean))
            corr = num / (g_std * r_std)
            if corr >= threshold:
                correlations[i, j] = corr
    return correlations, invalid_flags
