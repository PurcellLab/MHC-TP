"""Export embedded Seq2Logo reference logos from a parquet to PNG files.

Opt-in (the logos normally live only inside the parquet ``logo`` column and are
inlined into the HTML report). Useful for downstream / publication figure use.
"""

from __future__ import annotations

import os
from pathlib import Path

from hla_pepclust.refdata.parquet_io import read_reference


def export_logos(parquet_path: str | os.PathLike, out_dir: str | os.PathLike) -> int:
    """Write each row's embedded logo to ``<out_dir>/<formatted>.png``.

    Returns the number of logos written. Returns 0 if the parquet has no
    ``logo`` column (i.e. was built without ``--with-logos``).
    """
    df = read_reference(parquet_path)
    if "logo" not in df.columns:
        return 0
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    written = 0
    for row in df.itertuples():
        data = getattr(row, "logo", None)
        if data:
            (out / f"{row.formatted}.png").write_bytes(data)
            written += 1
    return written
