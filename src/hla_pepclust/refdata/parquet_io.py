"""Read/write the per-species reference Parquet (zstd-compressed)."""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd

from hla_pepclust.refdata.schema import COLUMNS


def write_reference(df: pd.DataFrame, path: str | os.PathLike) -> None:
    """Write the reference DataFrame to a zstd Parquet with the canonical columns."""
    missing = set(COLUMNS) - set(df.columns)
    if missing:
        raise ValueError(f"reference df missing columns: {sorted(missing)}")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df[COLUMNS].to_parquet(path, engine="pyarrow", compression="zstd", index=False)


def read_reference(path: str | os.PathLike) -> pd.DataFrame:
    """Load a reference Parquet into a DataFrame."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)
    return pd.read_parquet(p, engine="pyarrow")
