"""Assemble the standalone HTML report."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape

from mhc_tp.constants import N_AMINO_ACIDS
from mhc_tp.report.assets import find_cluster_logo, png_bytes_to_data_uri
from mhc_tp.report.data import datatable_rows, parse_cluster_id, pcc_records

_TEMPLATES = Path(__file__).parent / "templates"


def _render_logo(*args, **kwargs):
    """Lazy logomaker fallback: import matplotlib/logomaker (~1.6s) only when a
    logo actually has to be drawn (i.e. no embedded Seq2Logo / gibbs logo)."""
    from mhc_tp.report.logos import render_logo

    return render_logo(*args, **kwargs)


def render_report(
    correlation_dict,
    reference_df,
    gibbs_matrices,
    output_dir,
    kld_df: pd.DataFrame | None = None,
    version: str = "",
    gibbs_dir: str | None = None,
    logo_map: dict | None = None,
    name_map: dict | None = None,
    top_n: int = 3,
    threshold: float = 0.70,
    always_top_n: bool = False,
) -> str:
    """Write <output_dir>/clust_result/mhc-tp-result.html and return its path.

    ``logo_map`` ({formatted: png_bytes}) supplies reference logos when the
    reference DataFrame was loaded without the heavy ``logo`` column.
    ``name_map`` ({formatted: display}) supplies pretty allele labels.
    """
    logo_map = logo_map or {}
    name_map = name_map or {}
    ref_by_fmt = {r.formatted: r for r in reference_df.itertuples()}

    table_rows = datatable_rows(correlation_dict, kld_df, name_map, threshold)
    pcc_json = json.dumps(pcc_records(correlation_dict, name_map))

    # Best HLA per cluster id, grouped by the number of clusters N, with both
    # logos rendered from the matrices.
    sections: dict[int, list] = {}
    seen: set[str] = set()
    for (gibbs_name, hla), corr in sorted(
        correlation_dict.items(), key=lambda kv: -kv[1]
    ):
        cid, group, nclust = parse_cluster_id(gibbs_name)
        ref = ref_by_fmt.get(hla)
        if ref is None:
            continue
        if cid in seen:
            continue
        seen.add(cid)
        hla_display = name_map.get(hla, hla)
        # Reference logo: embedded Seq2Logo PNG from the parquet if present,
        # else the logomaker fallback rendered from the matrix.
        ref_logo_bytes = logo_map.get(hla) or getattr(ref, "logo", None)
        if ref_logo_bytes:
            ref_logo = png_bytes_to_data_uri(ref_logo_bytes)
        else:
            ref_mat = np.asarray(ref.matrix, dtype=np.float32).reshape(
                int(ref.n_positions), N_AMINO_ACIDS
            )
            ref_logo = _render_logo(ref_mat, title=hla_display)

        # Cluster logo: GibbsCluster's own Seq2Logo output if available, else fallback.
        cluster_logo = find_cluster_logo(gibbs_dir, cid) if gibbs_dir else None
        if not cluster_logo:
            gibbs_mat = gibbs_matrices.get(gibbs_name)
            cluster_logo = (
                _render_logo(gibbs_mat, title=cid) if gibbs_mat is not None else ""
            )

        sections.setdefault(nclust, []).append(
            {
                "cid": cid,
                "group": group,
                "hla": hla_display,
                "correlation": round(float(corr), 3),
                "below": float(corr) < threshold,
                "kld": _kld(kld_df, group, nclust),
                "ref_logo": ref_logo,
                "cluster_logo": cluster_logo,
            }
        )

    cluster_sections = [
        {
            "n_clusters": n,
            "groups": sorted(sections[n], key=lambda c: c["group"]),
        }
        for n in sorted(sections)
    ]

    env = Environment(
        loader=FileSystemLoader(str(_TEMPLATES)),
        autoescape=select_autoescape(["html", "j2"]),
    )
    template = env.get_template("report.html.j2")
    html = template.render(
        version=version,
        pcc_json=pcc_json,
        table_rows=table_rows,
        cluster_sections=cluster_sections,
        top_n=top_n,
        threshold=threshold,
        always_top_n=always_top_n,
    )

    out_dir = Path(output_dir) / "clust_result"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "mhc-tp-result.html"
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
