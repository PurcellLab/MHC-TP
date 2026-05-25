"""Parse the GibbsCluster ``gibbs.KLDvsClusters.tab`` file."""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd


def read_kld(file_path: str | os.PathLike) -> pd.DataFrame:
    """Return a DataFrame with columns ``cluster``, ``group1..groupN``, ``total``.

    Raises FileNotFoundError if the file is absent.
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(path)

    data: dict[str, list] = {"cluster": []}
    lines = path.read_text().splitlines()
    for line in lines[1:]:  # skip header
        parts = line.strip().split("\t")
        if not parts or parts == [""]:
            continue
        cluster = int(parts[0])
        values = [max(float(v), 0.0) for v in parts[1:]]
        for i in range(1, len(values) + 1):
            data.setdefault(f"group{i}", [])
        data["cluster"].append(cluster)
        for i, v in enumerate(values, start=1):
            data[f"group{i}"].append(v)

    df = pd.DataFrame(data).fillna(0.0)
    group_cols = [c for c in df.columns if c.startswith("group")]
    df["total"] = df[group_cols].sum(axis=1)
    return df.reset_index(drop=True)
