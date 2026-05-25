"""Search Gibbs cluster matrices against the reference array."""

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
) -> dict[tuple[str, str], float]:
    """Return ``{(gibbs_name, ref_formatted): correlation}`` for top-N hits."""
    ref_arr, max_positions = build_reference_array(reference)
    names = list(gibbs_matrices.keys())

    padded = np.zeros((len(names), max_positions, N_AMINO_ACIDS), dtype=np.float32)
    for i, name in enumerate(names):
        m = gibbs_matrices[name]
        padded[i, : m.shape[0], : m.shape[1]] = m

    mask = np.ones(len(reference), dtype=np.bool_)
    if hla_filter:
        mask = reference["formatted"].isin(hla_filter).to_numpy()

    corr, _invalid = compute_all_correlations(
        padded, ref_arr.astype(np.float32), mask, threshold
    )

    formatted = reference["formatted"].to_numpy()
    out: dict[tuple[str, str], float] = {}
    for i, name in enumerate(names):
        row = corr[i, :]
        order = np.argsort(row)[::-1]
        hits = [j for j in order if row[j] >= threshold][:top_n]
        for j in hits:
            out[(name, str(formatted[j]))] = float(row[j])
    return out
