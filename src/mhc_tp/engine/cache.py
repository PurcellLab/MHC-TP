"""Build the in-memory padded 3D reference array from a reference DataFrame."""

from __future__ import annotations

import numpy as np
import pandas as pd

from mhc_tp.constants import N_AMINO_ACIDS


def build_reference_array(df: pd.DataFrame) -> tuple[np.ndarray, int]:
    """Return ((n_allotypes, max_positions, 20) float32, max_positions)."""
    n = len(df)
    positions = df["n_positions"].to_numpy()
    max_positions = int(positions.max()) if n else 0
    arr = np.zeros((n, max_positions, N_AMINO_ACIDS), dtype=np.float32)
    for i, (flat, npos) in enumerate(zip(df["matrix"], positions)):
        m = np.asarray(flat, dtype=np.float32).reshape(int(npos), N_AMINO_ACIDS)
        arr[i, : int(npos), :] = m
    return arr, max_positions
