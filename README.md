# MHC-TP

Cluster immunopeptidomics peptides by their HLA/MHC binding motif and get a
ranked table plus a standalone interactive HTML report.

`mhc-tp` takes a **GibbsCluster** output folder, correlates each cluster's
position-specific scoring matrix against a reference of HLA/MHC **class I + II**
binding motifs (human & mouse), and writes the best allele match per cluster.

---

## For users

**Requirements:** Python 3.9–3.11.

### 1. Install

Clone the repo and install it (editable, so `git pull` updates the tool):

```bash
git clone https://github.com/PurcellLab/MHC-TP.git
cd MHC-TP
pip install -e .
```

> Prefer a one-liner without cloning? `pip install git+https://github.com/PurcellLab/MHC-TP.git`
> A virtual environment (`python -m venv .venv && source .venv/bin/activate`) is recommended.

### 2. Download the reference data (once)

The reference motifs are fetched from the GitHub release, not bundled:

```bash
mhc-tp fetch -s human     # or:  mouse  |  all
```

### 3. Run a search

```bash
mhc-tp search <gibbscluster_output_dir> -s human -o results/
```

`<gibbscluster_output_dir>` is a GibbsCluster run folder (it must contain a
`matrices/` subdirectory).

**Outputs** land in `results/clust_result/`:

| file | what it is |
|------|------------|
| `correlations.csv` | every cluster→allele match above the threshold (`hla` = display name, `formatted` = raw key, `correlation` = PCC) |
| `mhc-tp-result.html` | standalone report — open it in any browser |

### Common options

| flag | meaning | default |
|------|---------|---------|
| `-s, --species` | `human` or `mouse` | `human` |
| `-r, --reference` | path to a `<species>.parquet` (otherwise the fetched one is used) | auto |
| `-t, --threshold` | minimum Pearson correlation to report | `0.70` |
| `-o, --output` | output directory | `output` |
| `--threads` | max CPU threads (also `$MHC_TP_THREADS`) | `4` |
| `--no-html` | write only the CSV | off |
| `-l, --log` | also save the coloured session log | off |

Run `mhc-tp search --help` for the full list.

---

## For contributors / developers

The project uses [pixi](https://pixi.sh) for a reproducible dev environment
(Python 3.11) and a `src/` layout packaged with hatchling.

```bash
git clone https://github.com/PurcellLab/MHC-TP.git
cd MHC-TP
pixi install            # create the dev env from pixi.lock
pixi run dev-install    # editable-install the package into the env (run once)

pixi run test           # pytest
pixi run lint           # ruff
pixi run fmt            # black
```

Always run via `pixi run …` — a bare `python` may pick up a different
interpreter without the pinned dependencies.

### Rebuilding the reference data (dev only)

End users never do this. The per-species parquets are built once from the
NetMHCpan / NetMHCIIpan packs and uploaded to the release. Embedding the
Seq2Logo reference logos (`--with-logos`) needs a separate Python 2.7 env and
is slow — run it on a cluster:

```bash
mhc-tp build-ref <species> <classI_pack> <classII_pack> <out.parquet> \
    --with-logos --workers 16
# Seq2Logo itself runs in its own env:  pixi run -e seq2logo ...
```

### Layout

```
src/mhc_tp/
  cli.py            entry point (mhc-tp)
  engine/           numba correlation search
  refdata/          reference parquet read/write, fetch, schema
  report/           HTML report rendering (data, logos, templates)
  db/               DEV-ONLY reference-pack ingestion
  tui/              Rich console banner, logging, results table
tests/              pytest suite
```

---

## Citation

If you use MHC-TP in your work, please cite:

> Munday PR, Krishna SSG, Fehring J, Croft NP, Purcell AW, Li C, Braun A.
> Immunolyser 2.0: An advanced computational pipeline for comprehensive analysis
> of immunopeptidomic data. *Comput Struct Biotechnol J.* 2025;29:296–304.
> doi:[10.1016/j.csbj.2025.10.007](https://doi.org/10.1016/j.csbj.2025.10.007).
> PMID: [41209766](https://pubmed.ncbi.nlm.nih.gov/41209766/); PMCID: PMC12590289.

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

