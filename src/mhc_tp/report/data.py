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


def pcc_records(
    correlation_dict: dict[tuple[str, str], float],
    name_map: dict[str, str] | None = None,
) -> list[dict]:
    """Records for the D3 heatmap: {Cluster, HLA, Correlation}.

    ``name_map`` ({formatted: display}) supplies pretty allele labels.
    """
    name_map = name_map or {}
    out = []
    for (gibbs_name, hla), corr in correlation_dict.items():
        cid, _, _ = parse_cluster_id(gibbs_name)
        out.append(
            {
                "Cluster": cid,
                "HLA": name_map.get(hla, hla),
                "Correlation": round(float(corr), 4),
            }
        )
    return out


def datatable_rows(
    correlation_dict,
    kld_df: pd.DataFrame | None,
    name_map: dict[str, str] | None = None,
    threshold: float = 0.70,
) -> list[dict]:
    """Rows for the results DataTable, sorted by correlation desc, with KLD.

    ``name_map`` ({formatted: display}) supplies pretty allele labels.
    ``below`` flags rows whose correlation is under ``threshold`` (only possible
    when the search ran in always-top-N mode).
    """
    name_map = name_map or {}
    rows = []
    for (gibbs_name, hla), corr in correlation_dict.items():
        cid, group, nclust = parse_cluster_id(gibbs_name)
        kld = _lookup_kld(kld_df, group, nclust)
        rows.append(
            {
                "cluster": cid,
                "hla": name_map.get(hla, hla),
                "correlation": round(float(corr), 4),
                "below": float(corr) < threshold,
                "kld": kld,
            }
        )
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
