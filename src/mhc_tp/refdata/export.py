"""Export embedded Seq2Logo reference logos from a parquet to PNG files.

Opt-in (the logos normally live only inside the parquet ``logo`` column and are
inlined into the HTML report). Useful for downstream / publication figure use.
"""

from __future__ import annotations

import os
import re
from pathlib import Path

from mhc_tp.refdata.parquet_io import read_reference

_PREFIX_RE = re.compile(r"^(hla|h2)")


def _norm(name) -> str:
    """Loose match key: lowercase, alphanumerics only (drops ``-``/``*``/``:``…)."""
    return re.sub(r"[^a-z0-9]", "", str(name).lower())


def _row_keys(row) -> set[str]:
    """Match keys for a reference row: the raw allotype and the formatted key,
    each also with the leading ``HLA``/``H2`` species prefix stripped — so a
    user can pass ``HLA-B*39:124``, ``HLAB39124`` or ``B39124`` interchangeably.
    """
    norm_allo, norm_fmt = _norm(row.allotype), _norm(row.formatted)
    return {
        norm_allo,
        norm_fmt,
        _PREFIX_RE.sub("", norm_allo),
        _PREFIX_RE.sub("", norm_fmt),
    }


def export_logos(
    parquet_path: str | os.PathLike,
    out_dir: str | os.PathLike,
    allotypes: str | list[str] | None = None,
) -> int:
    """Write embedded logos to ``<out_dir>/<formatted>.png``.

    ``allotypes`` selects which logos to export:

    * ``None`` (default) — export every allotype that has an embedded logo.
    * a name or list of names — export only those. Matching is forgiving:
      ``HLA-B*39:124``, ``HLAB39124`` and ``B39124`` all resolve to the same row.

    The parquet is expected to carry ``allotype``, ``formatted`` and ``logo``
    columns. Returns the number of logos written; returns 0 if the parquet has
    no ``logo`` column (i.e. was built without ``--with-logos``). Raises
    ``ValueError`` if a requested allotype is not present.
    """
    df = read_reference(parquet_path)
    if "logo" not in df.columns:
        return 0
    missing_cols = {"allotype", "formatted", "logo"} - set(df.columns)
    if missing_cols:
        raise ValueError(
            f"Reference parquet is missing required columns: {sorted(missing_cols)}"
        )

    if allotypes is None:
        wanted = list(df.itertuples())
    else:
        if isinstance(allotypes, str):
            allotypes = [allotypes]
        by_key: dict[str, object] = {}
        for row in df.itertuples():
            for key in _row_keys(row):
                by_key.setdefault(key, row)
        wanted, missing, seen = [], [], set()
        for a in allotypes:
            row = by_key.get(_norm(a)) or by_key.get(_PREFIX_RE.sub("", _norm(a)))
            if row is None:
                missing.append(a)
            elif id(row) not in seen:
                seen.add(id(row))
                wanted.append(row)
        if missing:
            available = sorted(str(v) for v in df["allotype"].dropna().unique())
            preview = ", ".join(available[:10])
            more = f" (and {len(available) - 10} more)" if len(available) > 10 else ""
            raise ValueError(
                f"Requested allotype(s) not found: {missing}; "
                f"available e.g.: {preview}{more}"
            )

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    written = 0
    for row in wanted:
        data = getattr(row, "logo", None)
        if data:
            (out / f"{row.formatted}.png").write_bytes(data)
            written += 1
    return written
