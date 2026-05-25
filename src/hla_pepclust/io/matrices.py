"""Parse PSSM matrix text files (NetMHCpan / GibbsCluster layout)."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

import numpy as np

from hla_pepclust.constants import N_AMINO_ACIDS


def parse_matrix(file_path: str | os.PathLike) -> Optional[np.ndarray]:
    """Read a PSSM file into an (n_positions, 20) float32 array.

    Skips blank lines, comment lines (``#``) and the amino-acid header.
    For each data row, the last 20 whitespace-separated fields are the
    amino-acid scores in canonical order. Returns ``None`` if the file is
    missing or has no parseable rows.
    """
    path = Path(file_path)
    if not path.exists():
        return None

    rows: list[list[float]] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "A R N D" in line:
            continue
        parts = line.split()
        if len(parts) < N_AMINO_ACIDS:
            continue
        try:
            rows.append(
                [float(parts[-(N_AMINO_ACIDS - i)]) for i in range(N_AMINO_ACIDS)]
            )
        except (ValueError, IndexError):
            continue

    return np.array(rows, dtype=np.float32) if rows else None
