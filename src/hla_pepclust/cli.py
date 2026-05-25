"""Command-line interface for HLA-PepClust."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import pandas as pd

from hla_pepclust import __version__
from hla_pepclust.io.matrices import parse_matrix
from hla_pepclust.refdata.parquet_io import read_reference


def _load_gibbs_matrices(gibbs_dir: str, n_clusters: str = "all") -> dict:
    matrices_dir = Path(gibbs_dir) / "matrices"
    out = {}
    for fn in os.listdir(matrices_dir):
        if not fn.endswith(".mat"):
            continue
        if n_clusters != "all" and n_clusters.isdigit() and not fn.endswith(f"of{n_clusters}.mat"):
            continue
        m = parse_matrix(matrices_dir / fn)
        if m is not None:
            out[fn] = m
    return out


def run_search(gibbs_dir, reference, species, output, threshold=0.70, top_n=3,
               hla_filter=None, make_html=True):
    """Run the search; write correlations.csv and (default) the HTML report."""
    # Lazy import: pulls numba (~1.4s) only when a search actually runs, so
    # `clust-search --version` / `build-db` stay instant.
    from hla_pepclust.engine.search import search

    ref = read_reference(reference)
    gibbs = _load_gibbs_matrices(gibbs_dir)
    cd = search(ref, gibbs, threshold=threshold, top_n=top_n, hla_filter=hla_filter)

    rows = [
        {"cluster": name.replace(".mat", ""), "hla": hla, "correlation": round(corr, 4)}
        for (name, hla), corr in cd.items()
    ]
    df = pd.DataFrame(rows).sort_values("correlation", ascending=False) if rows else pd.DataFrame(columns=["cluster", "hla", "correlation"])
    out_dir = Path(output) / "clust_result"
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / "correlations.csv", index=False)

    if make_html:
        from hla_pepclust.io.kld import read_kld
        from hla_pepclust.report.render import render_report

        kld = None
        try:
            kld = read_kld(Path(gibbs_dir) / "images" / "gibbs.KLDvsClusters.tab")
        except FileNotFoundError:
            kld = None
        render_report(cd, ref, gibbs, output, kld_df=kld, version=__version__)

    return df


def main(argv=None):
    parser = argparse.ArgumentParser(prog="clust-search", description="HLA-PepClust")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command")

    s = sub.add_parser("search", help="run the cluster search (default)")
    s.add_argument("gibbs_folder")
    s.add_argument("-r", "--reference", required=True, help="path to <species>.parquet")
    s.add_argument("-s", "--species", default="human")
    s.add_argument("-t", "--threshold", type=float, default=0.70)
    s.add_argument("--topNHits", type=int, default=3)
    s.add_argument("-o", "--output", default="output")
    s.add_argument("--no-html", action="store_true", help="skip the HTML report (CSV only)")

    b = sub.add_parser("build-db", help="DEV: build a reference parquet")
    b.add_argument("db_csv")
    b.add_argument("matrix_root")
    b.add_argument("species")
    b.add_argument("out_parquet")
    b.add_argument("--source", default="NetMHCpan-4.2")

    args = parser.parse_args(argv)
    if args.command == "build-db":
        from hla_pepclust.db.construct import build_species_parquet
        build_species_parquet(args.db_csv, args.matrix_root, args.species, args.out_parquet, args.source)
        return
    if args.command == "search":
        run_search(args.gibbs_folder, args.reference, args.species, args.output,
                   threshold=args.threshold, top_n=args.topNHits, make_html=not args.no_html)
        return
    parser.print_help()


if __name__ == "__main__":
    main()
