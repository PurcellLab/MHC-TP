"""Command-line interface for HLA-PepClust."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from hla_pepclust.runtime import apply_thread_env

# Bound BLAS/OpenMP/numba thread pools BEFORE numpy/pandas are imported, so the
# tool stays a good citizen on shared servers (e.g. the Immunolyser backend).
apply_thread_env()

import pandas as pd  # noqa: E402
from hla_pepclust import __version__  # noqa: E402
from hla_pepclust.io.matrices import parse_matrix  # noqa: E402
from hla_pepclust.refdata.parquet_io import read_reference  # noqa: E402
from hla_pepclust.refdata.schema import COLUMNS  # noqa: E402
from hla_pepclust.tui import (  # noqa: E402
    banner,
    configure_logging,
    results_table,
    save_console_log,
)


def _load_gibbs_matrices(gibbs_dir: str, n_clusters: str = "all") -> dict:
    matrices_dir = Path(gibbs_dir) / "matrices"
    out = {}
    for fn in os.listdir(matrices_dir):
        if not fn.endswith(".mat"):
            continue
        if (
            n_clusters != "all"
            and n_clusters.isdigit()
            and not fn.endswith(f"of{n_clusters}.mat")
        ):
            continue
        m = parse_matrix(matrices_dir / fn)
        if m is not None:
            out[fn] = m
    return out


def run_search(
    gibbs_dir,
    reference,
    species,
    output,
    threshold=0.70,
    top_n=3,
    hla_filter=None,
    make_html=True,
    log_level="info",
    log_to_file=False,
    threads=None,
):
    """Run the search; write correlations.csv and (default) the HTML report."""
    # Lazy import: pulls numba (~1.4s) only when a search actually runs, so
    # `clust-search --version` / `build-db` stay instant.
    from hla_pepclust.engine.search import search
    from hla_pepclust.runtime import apply_numba_threads

    log = configure_logging(log_level, log_to_file)
    n_threads = apply_numba_threads(threads)
    log.debug("thread budget: %d", n_threads)
    out_dir = Path(output) / "clust_result"
    out_dir.mkdir(parents=True, exist_ok=True)

    log.info("[bold]Stage 1/4[/bold] loading reference …")
    # Load matrices + metadata only (skip the heavy logo blob column) for the search.
    ref = read_reference(reference, columns=COLUMNS)
    log.info(
        "Reference: [bold]%d[/bold] %s allotypes (class I+II) from %s",
        len(ref),
        species,
        reference,
    )

    log.info("[bold]Stage 2/4[/bold] parsing Gibbs cluster matrices …")
    gibbs = _load_gibbs_matrices(gibbs_dir)
    log.info("Loaded [bold]%d[/bold] Gibbs matrices from %s", len(gibbs), gibbs_dir)

    log.info(
        "[bold]Stage 3/4[/bold] numba correlation search (threshold %.2f) …", threshold
    )
    cd = search(ref, gibbs, threshold=threshold, top_n=top_n, hla_filter=hla_filter)
    log.info("[green]Found %d matches above threshold[/green]", len(cd))

    # Map the raw `formatted` key to the Immunolyser display name; keep both so
    # the CSV stays joinable on `formatted` while showing pretty `hla`.
    from hla_pepclust.naming import pretty_allele

    name_map = {r.formatted: pretty_allele(r.allotype) for r in ref.itertuples()}
    rows = [
        {
            "cluster": name.replace(".mat", ""),
            "hla": name_map.get(hla, hla),
            "formatted": hla,
            "correlation": round(corr, 4),
        }
        for (name, hla), corr in cd.items()
    ]
    df = (
        pd.DataFrame(rows).sort_values("correlation", ascending=False)
        if rows
        else pd.DataFrame(columns=["cluster", "hla", "formatted", "correlation"])
    )
    df.to_csv(out_dir / "correlations.csv", index=False)
    if len(df):
        results_table(df, threshold=threshold)

    log.info("[bold]Stage 4/4[/bold] rendering outputs …")
    if make_html:
        from hla_pepclust.io.kld import read_kld
        from hla_pepclust.report.render import render_report

        kld = None
        try:
            kld = read_kld(Path(gibbs_dir) / "images" / "gibbs.KLDvsClusters.tab")
        except FileNotFoundError:
            log.debug("no KLD file found; skipping KLD column")
        # Load embedded Seq2Logo logos ONLY for the matched alleles (targeted read).
        from hla_pepclust.refdata.parquet_io import load_logos

        logo_map = load_logos(reference, {hla for (_, hla) in cd})
        render_report(
            cd,
            ref,
            gibbs,
            output,
            kld_df=kld,
            version=__version__,
            gibbs_dir=gibbs_dir,
            logo_map=logo_map,
            name_map=name_map,
        )
        log.info("HTML report → %s", out_dir / "clust-search-result.html")
    log.info("[green]Done. Results in %s[/green]", out_dir)

    if log_to_file:
        save_console_log(str(out_dir / "search_cluster.log"))

    return df


def main(argv=None):
    from rich_argparse import RichHelpFormatter

    fmt = {"formatter_class": RichHelpFormatter}
    parser = argparse.ArgumentParser(
        prog="clust-search",
        description="HLA-PepClust: cluster immunopeptidomics peptides by HLA/MHC motif",
        **fmt,
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    sub = parser.add_subparsers(dest="command")

    s = sub.add_parser("search", help="run the cluster search (default)", **fmt)
    s.add_argument("gibbs_folder")
    s.add_argument(
        "-r",
        "--reference",
        default=None,
        help="path to <species>.parquet (default: fetched data dir; see `fetch`)",
    )
    s.add_argument("-s", "--species", default="human")
    s.add_argument("-t", "--threshold", type=float, default=0.70)
    s.add_argument("--topNHits", type=int, default=3)
    s.add_argument("-o", "--output", default="output")
    s.add_argument(
        "--no-html", action="store_true", help="skip the HTML report (CSV only)"
    )
    s.add_argument(
        "-l",
        "--log",
        action="store_true",
        help="save the coloured session log to the output dir",
    )
    s.add_argument(
        "--log-level",
        default="info",
        choices=["debug", "info", "warning", "error", "critical"],
        help="logging verbosity (default: info)",
    )
    s.add_argument(
        "--threads",
        type=int,
        default=None,
        help="max CPU threads for the search (default: $HLA_PEPCLUST_THREADS or 4)",
    )

    b = sub.add_parser("build-db", help="DEV: build a reference parquet", **fmt)
    b.add_argument("db_csv")
    b.add_argument("matrix_root")
    b.add_argument("species")
    b.add_argument("out_parquet")
    b.add_argument("--source", default="NetMHCpan-4.2")
    b.add_argument(
        "--with-logos",
        action="store_true",
        help="embed Seq2Logo reference logos (needs SEQ2LOGO_PATH)",
    )
    b.add_argument(
        "--seq2logo-path",
        default=os.environ.get("SEQ2LOGO_PATH"),
        help="seq2logo-2.1 dir (default: $SEQ2LOGO_PATH)",
    )
    b.add_argument(
        "--seq2logo-python",
        default=os.environ.get("SEQ2LOGO_PYTHON"),
        help="python2.7 interpreter for Seq2Logo (default: $SEQ2LOGO_PYTHON)",
    )

    br = sub.add_parser(
        "build-ref",
        help="DEV: build a combined <species>.parquet from class I + II packs",
        **fmt,
    )
    br.add_argument("species")
    br.add_argument("class_i_pack", help="extracted NetMHCpan class I pack dir")
    br.add_argument("class_ii_pack", help="extracted NetMHCIIpan class II pack dir")
    br.add_argument("out_parquet")
    br.add_argument(
        "--with-logos",
        action="store_true",
        help="embed Seq2Logo reference logos (needs SEQ2LOGO_PATH; slow)",
    )
    br.add_argument("--seq2logo-path", default=os.environ.get("SEQ2LOGO_PATH"))
    br.add_argument("--seq2logo-python", default=os.environ.get("SEQ2LOGO_PYTHON"))
    br.add_argument(
        "--workers",
        type=int,
        default=1,
        help="parallel Seq2Logo render workers (default 1)",
    )

    fp = sub.add_parser(
        "fetch", help="download prebuilt reference parquets to the data dir", **fmt
    )
    fp.add_argument("-s", "--species", default="all", choices=["human", "mouse", "all"])
    fp.add_argument("-d", "--dest", default=None, help="override the data dir")

    ep = sub.add_parser(
        "export-logos",
        help="export embedded reference logos from a parquet to PNG files",
        **fmt,
    )
    ep.add_argument(
        "-r", "--reference", required=True, help="path to <species>.parquet"
    )
    ep.add_argument(
        "-o", "--output", required=True, help="output directory for PNG files"
    )
    ep.add_argument(
        "-a",
        "--allotypes",
        nargs="+",
        default=None,
        metavar="ALLOTYPE",
        help="allotypes to export, space- and/or comma-separated, e.g. "
        "'HLA-B*39:124 HLAB3942' or 'HLA-B39:124,A0201' (default: all). "
        "Matching ignores prefixes/separators.",
    )

    args = parser.parse_args(argv)
    if args.command == "fetch":
        from hla_pepclust.refdata.fetch import data_dir, fetch

        paths = fetch(args.species, args.dest)
        print(f"fetched {len(paths)} file(s) to {args.dest or data_dir()}:")
        for p in paths:
            print(f"  {p}")
        return
    if args.command == "export-logos":
        from hla_pepclust.refdata.export import export_logos

        # Accept both space-separated (nargs) and comma-separated tokens.
        allotypes = None
        if args.allotypes:
            allotypes = [
                a.strip()
                for chunk in args.allotypes
                for a in chunk.split(",")
                if a.strip()
            ] or None
        n = export_logos(args.reference, args.output, allotypes)
        print(f"exported {n} reference logos to {args.output}")
        return
    if args.command == "build-db":
        from hla_pepclust.db.construct import build_species_parquet

        build_species_parquet(
            args.db_csv,
            args.matrix_root,
            args.species,
            args.out_parquet,
            args.source,
            with_logos=args.with_logos,
            seq2logo_path=args.seq2logo_path,
            seq2logo_python=args.seq2logo_python,
        )
        return
    if args.command == "build-ref":
        from hla_pepclust.db.pack import build_species_reference

        n_i, n_ii = build_species_reference(
            args.species,
            args.class_i_pack,
            args.class_ii_pack,
            args.out_parquet,
            with_logos=args.with_logos,
            seq2logo_path=args.seq2logo_path,
            seq2logo_python=args.seq2logo_python,
            workers=args.workers,
        )
        print(f"{args.species}: class I={n_i}, class II={n_ii} -> {args.out_parquet}")
        return
    if args.command == "search":
        from hla_pepclust.refdata.fetch import resolve_reference

        banner()
        reference = str(resolve_reference(args.species, args.reference))
        run_search(
            args.gibbs_folder,
            reference,
            args.species,
            args.output,
            threshold=args.threshold,
            top_n=args.topNHits,
            make_html=not args.no_html,
            log_level=args.log_level,
            log_to_file=args.log,
            threads=args.threads,
        )
        return
    parser.print_help()


if __name__ == "__main__":
    main()
