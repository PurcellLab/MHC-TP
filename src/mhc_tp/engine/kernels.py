"""Numba JIT correlation kernel (ported from the proven NumbaSearch engine)."""

import numpy as np
from numba import jit, prange


# NOTE: no fastmath — it would let LLVM assume no-NaN and silently disable the
# np.isnan() masking below (NaN cells must be excluded, per the docstring).
@jit(nopython=True, parallel=True, cache=True)
def compute_all_correlations(gibbs_matrices, ref_matrices, hla_mask, threshold):
    r"""All-pairs flattened Pearson correlation, parallel over Gibbs matrices.

    Each PSSM is flattened to a vector. Only the cells that are informative in
    the Gibbs matrix are scored: the valid set is

        V = { k : g_k != 0 and g_k is not NaN }

    Restricted to those cells, the score for a (Gibbs g, reference r) pair is
    the Pearson correlation coefficient

                  (1/|V|) * Σ_{k in V} (g_k - ḡ)(r_k - r̄)
        PCC(g,r) = ----------------------------------------
                                σ_g · σ_r

    where ḡ, r̄ are the means and σ_g, σ_r the population standard deviations
    taken over V:

        ḡ   = (1/|V|) Σ g_k ,                σ_g = sqrt( (1/|V|) Σ (g_k - ḡ)^2 )
        r̄   = (1/|V|) Σ r_k ,                σ_r = sqrt( (1/|V|) Σ (r_k - r̄)^2 )

    PCC lies in [-1, 1] (1 = identical motif shape) and is scale/offset
    invariant, so it measures the *pattern* of position preferences rather than
    absolute weight magnitudes.

    Guards: a Gibbs matrix with |V| < 10 or σ_g = 0 is flagged invalid (its row
    is skipped); a reference with σ_r = 0 is skipped for that pair. A score is
    stored only when PCC >= ``threshold``; otherwise the cell keeps the -1.0
    sentinel. Returns (correlations, invalid_flags).
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
