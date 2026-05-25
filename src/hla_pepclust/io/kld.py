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
        data["cluster"].append(cluster)
        n = len(data["cluster"])
        for i, v in enumerate(values, start=1):
            col = f"group{i}"
            if col not in data:  # new group column: back-fill earlier rows with 0
                data[col] = [0.0] * (n - 1)
            data[col].append(v)
        # ragged rows: pad any group column this row lacked so all stay aligned
        for col, col_values in data.items():
            if col != "cluster" and len(col_values) < n:
                col_values.append(0.0)

    df = pd.DataFrame(data).fillna(0.0)
    group_cols = [c for c in df.columns if c.startswith("group")]
    df["total"] = df[group_cols].sum(axis=1)
    return df.reset_index(drop=True)
