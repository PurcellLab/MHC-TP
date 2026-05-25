"""Read/write the per-species reference Parquet (zstd-compressed)."""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd

from mhc_tp.refdata.schema import COLUMNS


def write_reference(df: pd.DataFrame, path: str | os.PathLike) -> None:
    """Write the reference DataFrame to a zstd Parquet with the canonical columns.

    The optional ``logo`` column (Seq2Logo PNG bytes) is written when present.
    """
    missing = set(COLUMNS) - set(df.columns)
    if missing:
        raise ValueError(f"reference df missing columns: {sorted(missing)}")
    cols = list(COLUMNS) + (["logo"] if "logo" in df.columns else [])
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df[cols].to_parquet(path, engine="pyarrow", compression="zstd", index=False)


def read_reference(
    path: str | os.PathLike, columns: list[str] | None = None
) -> pd.DataFrame:
    """Load a reference Parquet. Pass ``columns`` to skip heavy ones (e.g. logo).

    For a search, pass ``columns=COLUMNS`` to avoid loading the large ``logo``
    blob column into memory.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)
    return pd.read_parquet(p, engine="pyarrow", columns=columns)


def load_logos(path: str | os.PathLike, formatted_ids) -> dict[str, bytes]:
    """Return ``{formatted: logo_bytes}`` for the given ids, reading ONLY the
    logo column of the matching rows (predicate pushdown) — not the whole file.

    Empty dict if the parquet has no ``logo`` column or no ids match.
    """
    import pyarrow.parquet as pq

    p = Path(path)
    ids = list(dict.fromkeys(formatted_ids))
    if not p.exists() or not ids:
        return {}
    if "logo" not in pq.ParquetFile(p).schema_arrow.names:
        return {}
    table = pq.read_table(
        p, columns=["formatted", "logo"], filters=[("formatted", "in", ids)]
    )
    cols = table.to_pydict()
    return {fmt: blob for fmt, blob in zip(cols["formatted"], cols["logo"]) if blob}
