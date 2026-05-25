"""Shape search results into report-ready structures."""

from __future__ import annotations

import re

import pandas as pd

_CLUSTER_RE = re.compile(r"(\d+)of(\d+)")


def parse_cluster_id(gibbs_name: str) -> tuple[str, int, int]:
    """``gibbs.1of3.mat`` -> ("1of3", group=1, n_clusters=3)."""
    m = _CLUSTER_RE.search(gibbs_name)
    if not m:
        return (gibbs_name, 0, 0)
    group, nclust = int(m.group(1)), int(m.group(2))
    return (f"{group}of{nclust}", group, nclust)


def pcc_records(correlation_dict: dict[tuple[str, str], float]) -> list[dict]:
    """Records for the D3 heatmap: {Cluster, HLA, Correlation}."""
    out = []
    for (gibbs_name, hla), corr in correlation_dict.items():
        cid, _, _ = parse_cluster_id(gibbs_name)
        out.append({"Cluster": cid, "HLA": hla, "Correlation": round(float(corr), 4)})
    return out


def datatable_rows(correlation_dict, kld_df: pd.DataFrame | None) -> list[dict]:
    """Rows for the results DataTable, sorted by correlation desc, with KLD."""
    rows = []
    for (gibbs_name, hla), corr in correlation_dict.items():
        cid, group, nclust = parse_cluster_id(gibbs_name)
        kld = _lookup_kld(kld_df, group, nclust)
        rows.append({"cluster": cid, "hla": hla,
                     "correlation": round(float(corr), 4), "kld": kld})
    rows.sort(key=lambda r: r["correlation"], reverse=True)
    return rows


def _lookup_kld(kld_df, group: int, nclust: int):
    if kld_df is None or nclust == 0:
        return None
    col = f"group{group}"
    sub = kld_df[kld_df["cluster"] == nclust]
    if sub.empty or col not in kld_df.columns:
        return None
    val = sub[col].iloc[0]
    return round(float(val), 4)
