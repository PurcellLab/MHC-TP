"""DEV-ONLY: build a per-species reference Parquet from raw matrices.

Not part of the user runtime. Run via the ``clust-search build-db`` subcommand.
"""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd

from hla_pepclust.io.matrices import parse_matrix
from hla_pepclust.io.naming import format_allotype
from hla_pepclust.refdata.parquet_io import write_reference

# Defaults reflect the real <species>.db columns (verified in Step 0).
DEFAULT_ALLOTYPE_COL = "allotypes"
DEFAULT_MATRIX_PATH_COL = "matrices_path"


def build_species_parquet(
    db_csv: str | os.PathLike,
    matrix_root: str | os.PathLike,
    species: str,
    out_parquet: str | os.PathLike,
    source: str,
    allotype_col: str = DEFAULT_ALLOTYPE_COL,
    matrix_path_col: str = DEFAULT_MATRIX_PATH_COL,
) -> None:
    """Read the legacy <species>.db CSV + matrix files -> write <species>.parquet."""
    db = pd.read_csv(db_csv)
    matrix_root = Path(matrix_root)
    rows = []
    for _, r in db.iterrows():
        mat = parse_matrix(matrix_root / str(r[matrix_path_col]))
        if mat is None:
            continue
        info = format_allotype(str(r[allotype_col]), species=species)
        rows.append({
            "allotype": info.raw,
            "formatted": info.formatted,
            "mhc_class": info.mhc_class,
            "locus": info.locus,
            "n_positions": int(mat.shape[0]),
            "matrix": mat.reshape(-1).tolist(),
            "source": source,
        })
    write_reference(pd.DataFrame(rows), out_parquet)
