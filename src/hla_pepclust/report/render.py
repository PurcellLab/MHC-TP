"""Assemble the standalone HTML report."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape

from hla_pepclust.constants import N_AMINO_ACIDS
from hla_pepclust.report.assets import find_cluster_logo, png_bytes_to_data_uri
from hla_pepclust.report.data import datatable_rows, parse_cluster_id, pcc_records
from hla_pepclust.report.logos import render_logo

_TEMPLATES = Path(__file__).parent / "templates"


def render_report(
    correlation_dict,
    reference_df,
    gibbs_matrices,
    output_dir,
    kld_df: pd.DataFrame | None = None,
    version: str = "",
    gibbs_dir: str | None = None,
) -> str:
    """Write <output_dir>/clust_result/clust-search-result.html and return its path."""
    ref_by_fmt = {r.formatted: r for r in reference_df.itertuples()}

    table_rows = datatable_rows(correlation_dict, kld_df)
    pcc_json = json.dumps(pcc_records(correlation_dict))

    # Best HLA per (cluster, class), with both logos rendered from the matrices.
    groups: dict[str, list] = {}
    seen: set[tuple[str, str]] = set()
    for (gibbs_name, hla), corr in sorted(
        correlation_dict.items(), key=lambda kv: -kv[1]
    ):
        cid, group, nclust = parse_cluster_id(gibbs_name)
        ref = ref_by_fmt.get(hla)
        if ref is None:
            continue
        key = (cid, ref.mhc_class)
        if key in seen:
            continue
        seen.add(key)
        # Reference logo: embedded Seq2Logo PNG from the parquet if present,
        # else the logomaker fallback rendered from the matrix.
        ref_logo_bytes = getattr(ref, "logo", None)
        if ref_logo_bytes:
            ref_logo = png_bytes_to_data_uri(ref_logo_bytes)
        else:
            ref_mat = np.asarray(ref.matrix, dtype=np.float32).reshape(
                int(ref.n_positions), N_AMINO_ACIDS
            )
            ref_logo = render_logo(ref_mat, title=hla)

        # Cluster logo: GibbsCluster's own Seq2Logo output if available, else fallback.
        cluster_logo = find_cluster_logo(gibbs_dir, cid) if gibbs_dir else None
        if not cluster_logo:
            gibbs_mat = gibbs_matrices.get(gibbs_name)
            cluster_logo = (
                render_logo(gibbs_mat, title=cid) if gibbs_mat is not None else ""
            )

        groups.setdefault(ref.mhc_class, []).append(
            {
                "cluster_id": cid,
                "hla": hla,
                "correlation": round(float(corr), 3),
                "kld": _kld(kld_df, group, nclust),
                "ref_logo": ref_logo,
                "cluster_logo": cluster_logo,
            }
        )

    class_groups = [{"mhc_class": c, "clusters": groups[c]} for c in sorted(groups)]

    env = Environment(
        loader=FileSystemLoader(str(_TEMPLATES)),
        autoescape=select_autoescape(["html", "j2"]),
    )
    template = env.get_template("report.html.j2")
    html = template.render(
        version=version,
        pcc_json=pcc_json,
        table_rows=table_rows,
        class_groups=class_groups,
    )

    out_dir = Path(output_dir) / "clust_result"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "clust-search-result.html"
    out_path.write_text(html)
    return str(out_path)


def _kld(kld_df, group, nclust):
    if kld_df is None or nclust == 0:
        return None
    col = f"group{group}"
    sub = kld_df[kld_df["cluster"] == nclust]
    if sub.empty or col not in kld_df.columns:
        return None
    return round(float(sub[col].iloc[0]), 4)
