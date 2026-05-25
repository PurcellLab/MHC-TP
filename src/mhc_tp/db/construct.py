"""DEV-ONLY: build a per-species reference Parquet from raw matrices.

Not part of the user runtime. Run via the ``mhc-tp build-db`` subcommand.
"""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd

from mhc_tp.io.matrices import parse_matrix
from mhc_tp.io.naming import format_allotype
from mhc_tp.refdata.parquet_io import write_reference

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
    with_logos: bool = False,
    seq2logo_path: str | None = None,
    seq2logo_python: str | None = None,
) -> None:
    """Read the legacy <species>.db CSV + matrix files -> write <species>.parquet.

    When ``with_logos`` is set, generate a Seq2Logo reference logo per allele
    (via the external Seq2Logo install) and embed it in a ``logo`` column.
    """
    db = pd.read_csv(db_csv)
    matrix_root = Path(matrix_root)
    rows = []
    for _, r in db.iterrows():
        matrix_file = matrix_root / str(r[matrix_path_col])
        mat = parse_matrix(matrix_file)
        if mat is None:
            continue
        info = format_allotype(str(r[allotype_col]), species=species)
        row = {
            "allotype": info.raw,
            "formatted": info.formatted,
            "mhc_class": info.mhc_class,
            "locus": info.locus,
            "n_positions": int(mat.shape[0]),
            "matrix": mat.reshape(-1).tolist(),
            "source": source,
        }
        if with_logos:
            from mhc_tp.db.logos import reference_logo_bytes

            row["logo"] = reference_logo_bytes(
                matrix_file,
                seq2logo_path=seq2logo_path,
                python_exe=seq2logo_python,
                title=info.raw,
            )
        rows.append(row)
    write_reference(pd.DataFrame(rows), out_parquet)
