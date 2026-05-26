# MHC-TP

[![PyPI](https://img.shields.io/pypi/v/mhc-tp.svg)](https://pypi.org/project/mhc-tp/)
[![Python](https://img.shields.io/pypi/pyversions/mhc-tp.svg)](https://pypi.org/project/mhc-tp/)
[![Docs](https://img.shields.io/badge/docs-mkdocs--material-526CFE.svg)](https://purcelllab.github.io/MHC-TP/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Cluster immunopeptidomics peptides by their **HLA/MHC binding motif** and get a
ranked table plus a standalone, interactive HTML report.

`mhc-tp` takes a **GibbsCluster** output folder, correlates each cluster's
position-specific scoring matrix against a reference of HLA/MHC **class I + II**
binding motifs (human & mouse), and reports the best allele match per cluster.

> [!TIP]
> 📖 **[Documentation](https://purcelllab.github.io/MHC-TP/)** &nbsp;·&nbsp;
> 🔬 **[Live example report](https://purcelllab.github.io/MHC-TP/example-report.html)** &nbsp;·&nbsp;
> 📦 **[PyPI](https://pypi.org/project/mhc-tp/)**

---

## Quick start

```bash
pip install mhc-tp
mhc-tp fetch -s human                                  # download reference motifs (once)
mhc-tp search <gibbscluster_output_dir> -s human -o results/
```

Open `results/clust_result/mhc-tp-result.html` in any browser — see what it looks like in the
**[live example report](https://purcelllab.github.io/MHC-TP/example-report.html)**.

> [!NOTE]
> **Requirements:** Python 3.9–3.11. A virtual environment is recommended
> (`python -m venv .venv && source .venv/bin/activate`).

---

## Install

From PyPI (recommended):

```bash
pip install mhc-tp
```

<details>
<summary>Or install editable from source</summary>

So that `git pull` updates the tool:

```bash
git clone https://github.com/PurcellLab/MHC-TP.git
cd MHC-TP
pip install -e .
```

One-liner without cloning: `pip install "git+https://github.com/PurcellLab/MHC-TP.git"`

</details>

## Download the reference data (once)

The reference motifs are fetched from the GitHub release, not bundled:

```bash
mhc-tp fetch -s human     # or:  mouse  |  all
```

## Run a search

```bash
mhc-tp search <gibbscluster_output_dir> -s human -o results/
```

`<gibbscluster_output_dir>` is a GibbsCluster run folder (it must contain a `matrices/` subdirectory).

**Outputs** land in `results/clust_result/`:

| file | what it is |
|------|------------|
| `correlations.csv` | every cluster→allele match (`hla` = display name, `formatted` = raw key, `correlation` = PCC) |
| `mhc-tp-result.html` | standalone interactive report — open it in any browser |

### Options

| flag | meaning | default |
|------|---------|---------|
| `-s, --species` | `human` or `mouse` | `human` |
| `-c, --class` | restrict the reference to MHC class `I`, `II`, or `all` | `all` |
| `-r, --reference` | path to a `<species>.parquet` (otherwise the fetched one is used) | auto |
| `-t, --threshold` | minimum Pearson correlation (PCC) to report | `0.70` |
| `--topNHits` | allele matches to keep per cluster | `3` |
| `--always-top-n` | keep each cluster's top-N even below threshold (flagged in the report) | off |
| `-o, --output` | output directory | `output` |
| `--threads` | max CPU threads (also `$MHC_TP_THREADS`) | `4` |
| `--no-html` | write only the CSV | off |
| `-l, --log` | also save the coloured session log | off |

Run `mhc-tp search --help` for the full list.

### Examples

```bash
# Class I only, keep the top 5 matches per cluster
mhc-tp search runs/sampleA -s human -c I --topNHits 5 -o results/

# Guarantee a top-3 for every cluster (weak matches tagged "below cutoff")
mhc-tp search runs/sampleA -s human --always-top-n -o results/
```

> [!IMPORTANT]
> By default a match must score `≥ --threshold`, so a cluster can return fewer than
> `--topNHits` rows (or none). `--always-top-n` returns the best N regardless — the
> threshold then only **annotates** confidence and nothing is dropped.

---

## For contributors / developers

<details>
<summary>Dev environment, tests, and docs (click to expand)</summary>

The project uses [pixi](https://pixi.sh) for a reproducible dev environment (Python 3.11)
and a `src/` layout packaged with hatchling.

```bash
git clone https://github.com/PurcellLab/MHC-TP.git
cd MHC-TP
pixi install            # create the dev env from pixi.lock
pixi run dev-install    # editable-install the package (run once)

pixi run test           # pytest
pixi run lint           # ruff
pixi run fmt            # black
```

> [!WARNING]
> Always run via `pixi run …` — a bare `python` may pick up a different interpreter
> without the pinned dependencies. CI enforces `black --check`, so run `pixi run fmt` before pushing.

### Preview the docs site

```bash
pip install -e ".[docs]"
mkdocs serve            # live preview at http://127.0.0.1:8000
mkdocs build            # static site in ./site
```

### Rebuilding the reference data (dev only)

End users never do this. The per-species parquets are built once from the
NetMHCpan / NetMHCIIpan packs and uploaded to the release. Embedding the
Seq2Logo reference logos (`--with-logos`) needs a separate Python 2.7 env and is slow:

```bash
mhc-tp build-ref <species> <classI_pack> <classII_pack> <out.parquet> --with-logos --workers 16
```

### Layout

```text
src/mhc_tp/
  cli.py            entry point (mhc-tp)
  engine/           numba correlation search
  refdata/          reference parquet read/write, fetch, schema
  report/           HTML report rendering (data, logos, templates)
  db/               DEV-ONLY reference-pack ingestion
  tui/              Rich console banner, logging, results table
tests/              pytest suite
docs/               MkDocs site
```

</details>

---

## How it works

For each GibbsCluster motif, every reference allotype motif is scored by the **Pearson
correlation** of their flattened position-weight matrices, computed only over the
informative cells of the cluster motif. Per cluster the allotypes are ranked by PCC
(`1.0` = identical motif shape). Full method and formula:
the **[API reference](https://purcelllab.github.io/MHC-TP/api/)**.

---

## Citation

If you use MHC-TP in your work, please cite:

> Munday PR, Krishna SSG, Fehring J, Croft NP, Purcell AW, Li C, Braun A.
> *Immunolyser 2.0: An advanced computational pipeline for comprehensive analysis of
> immunopeptidomic data.* Comput Struct Biotechnol J. 2025;29:296–304.
> doi:[10.1016/j.csbj.2025.10.007](https://doi.org/10.1016/j.csbj.2025.10.007).
> PMID: [41209766](https://pubmed.ncbi.nlm.nih.gov/41209766/); PMCID: PMC12590289.

<details>
<summary>BibTeX</summary>

```bibtex
@article{Munday2025Immunolyser2,
  title   = {Immunolyser 2.0: An advanced computational pipeline for comprehensive analysis of immunopeptidomic data},
  author  = {Munday, Prithvi Raj and Krishna, Sanjay S. G. and Fehring, Joshua and Croft, Nathan P. and Purcell, Anthony W. and Li, Chen and Braun, Asolina},
  journal = {Computational and Structural Biotechnology Journal},
  volume  = {29},
  pages   = {296--304},
  year    = {2025},
  doi     = {10.1016/j.csbj.2025.10.007},
  pmid    = {41209766},
  pmcid   = {PMC12590289}
}
```

</details>
