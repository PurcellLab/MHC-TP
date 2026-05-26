# LLM guide — complete operating instructions

> This page is a single, self-contained reference written for an LLM assistant to
> help an end user run **every** MHC-TP use case correctly. It restates the
> install, inputs, commands, flags, recipes, outputs, and the matching method so
> no other page is required.
>
> **Raw text:** a machine-readable copy of this page is served at
> [`/llm.txt`](../llm.txt). Use the **Copy page as text** button (top right) to
> copy the whole guide to your clipboard.

## What MHC-TP does

`mhc-tp` is a command-line tool. Given a **GibbsCluster** output folder, it
correlates each peptide cluster's position-specific scoring matrix (PSSM) against
a reference library of HLA/MHC **class I + II** binding motifs (human and mouse)
and reports, per cluster, the best-matching allotype(s). It writes a ranked CSV
and a standalone interactive HTML report.

- Package / import name: `mhc-tp` / `mhc_tp`. Console command: `mhc-tp`.
- Python 3.9–3.11. No internet needed at search time (only once, to fetch the reference).

## Install

```bash
pip install mhc-tp
```

Editable from source (for development):

```bash
git clone https://github.com/PurcellLab/MHC-TP.git
cd MHC-TP
pip install -e .
```

## Reference data (required, fetched once)

Reference motif parquets are downloaded from the GitHub release, not bundled.

```bash
mhc-tp fetch -s human     # human  |  mouse  |  all
```

`fetch` options: `-s/--species {human,mouse,all}` (default `all`), `-d/--dest DIR`
(override the data dir; otherwise a per-user data dir, also settable via
`$MHC_TP_DATA_DIR`). After fetching, `mhc-tp search` finds the parquet
automatically — you only pass `-r/--reference` to use a custom file.

## Input: the GibbsCluster folder

The positional `gibbs_folder` argument is a GibbsCluster run directory. It **must
contain a `matrices/` subdirectory** with the per-cluster matrix files
(`gibbs.<g>of<N>.mat`). If a `images/gibbs.KLDvsClusters.tab` file is present, the
report also shows each cluster's KLD (information content). GibbsCluster's own
logo PNGs, if present, are reused in the report.

## Command: `mhc-tp search`

```text
mhc-tp search [-h] [-r REFERENCE] [-s SPECIES] [-c {I,II,all}]
              [-t THRESHOLD] [--topNHits TOPNHITS] [--always-top-n]
              [-o OUTPUT] [--no-html] [-l]
              [--log-level {debug,info,warning,error,critical}]
              [--threads THREADS]
              gibbs_folder
```

| flag | meaning | default |
|------|---------|---------|
| `gibbs_folder` | GibbsCluster run dir (must have `matrices/`) | required |
| `-s, --species` | `human` or `mouse` | `human` |
| `-c, --class` | restrict reference to MHC class `I`, `II`, or `all` | `all` |
| `-r, --reference` | path to a `<species>.parquet` (else the fetched one) | auto |
| `-t, --threshold` | minimum Pearson correlation (PCC) to report | `0.70` |
| `--topNHits` | allotype matches to keep per cluster | `3` |
| `--always-top-n` | keep each cluster's top-N even below threshold (flagged "below cutoff") | off |
| `-o, --output` | output directory | `output` |
| `--no-html` | write only the CSV, skip the HTML report | off |
| `-l, --log` | also save the coloured session log to the output dir | off |
| `--log-level` | `debug`/`info`/`warning`/`error`/`critical` | `info` |
| `--threads` | max CPU threads (also `$MHC_TP_THREADS`) | `4` |

## Use-case recipes (copy-paste)

```bash
# 1. Basic human search (class I + II reference)
mhc-tp search runs/sampleA -s human -o results/

# 2. Mouse sample
mhc-tp fetch -s mouse
mhc-tp search runs/mouseA -s mouse -o results_mouse/

# 3. Restrict to MHC class I only (faster, class-I immunopeptidome)
mhc-tp search runs/sampleA -s human -c I -o results/

# 4. Class II only
mhc-tp search runs/sampleA -s human -c II -o results/

# 5. Keep the top 5 matches per cluster instead of 3
mhc-tp search runs/sampleA -s human --topNHits 5 -o results/

# 6. Guarantee a top-3 for EVERY cluster even if below threshold
#    (weak matches are tagged "below cutoff" in the report)
mhc-tp search runs/sampleA -s human --always-top-n -o results/

# 7. Stricter / looser confidence cutoff
mhc-tp search runs/sampleA -s human -t 0.80 -o results/

# 8. CSV only (no HTML report)
mhc-tp search runs/sampleA -s human --no-html -o results/

# 9. Use a custom / local reference parquet
mhc-tp search runs/sampleA -s human -r path/to/human.parquet -o results/

# 10. Limit CPU threads and save a log
mhc-tp search runs/sampleA -s human --threads 8 -l -o results/
```

> [!IMPORTANT]
> **Threshold vs top-N selection order:** per cluster, all allotypes are ranked by
> PCC, then the threshold is applied, then the top-N are taken. So by default a
> cluster can return fewer than `--topNHits` rows (or none). With `--always-top-n`,
> the top-N are taken regardless of threshold and the threshold only annotates
> confidence — every cluster keeps its best N.

## Outputs

Written to `<output>/clust_result/`:

| file | description |
|------|-------------|
| `correlations.csv` | columns: `cluster` (e.g. `gibbs.2of5`), `hla` (canonical display name, e.g. `HLA-B39:124`), `formatted` (raw join key, e.g. `HLAB39124`), `correlation` (PCC) |
| `mhc-tp-result.html` | standalone interactive report (open in any browser) |

The report contains: a motif-comparison carousel (per cluster solution, choose via
dropdown), a paginated + searchable Top-N table, and a correlation analysis section
(heatmap + force-directed network, with a PCC threshold slider that filters the
view only). Rows/motifs below the search threshold (only present with
`--always-top-n`) are tagged **below cutoff**.

## How matching works (method)

Each motif is a PSSM (`n_positions × 20`). For a cluster motif `g` and reference
`r`, the score is the **Pearson correlation** over the cells `V` that are
informative in the cluster motif (`g_k ≠ 0` and not NaN):

```text
PCC(g, r) = Σ_{k∈V}(g_k − ḡ)(r_k − r̄) / ( |V| · σ_g · σ_r )    ∈ [−1, 1]
```

`1.0` = identical motif shape. It is scale/offset invariant (scores the pattern of
position preferences, not absolute magnitudes). Allele display names are returned
verbatim from the reference (matching the embedded logo titles), e.g. `HLA-A25:08`.

## Other commands

```bash
# Export embedded reference logos to PNGs
mhc-tp export-logos -r human.parquet -o logos/ -a "HLA-B39:124,A0201"   # -a optional; default all
```

Developer-only (end users never need these): `mhc-tp build-db` and
`mhc-tp build-ref` rebuild the reference parquets from NetMHCpan / NetMHCIIpan
packs; embedding Seq2Logo logos (`--with-logos`) needs a separate Python 2.7 env.

## Troubleshooting

- **"no class II allotypes" / 0 matches with `-c II`**: the sample is likely class I
  (short 8–11mers); use `-c I` or `-c all`.
- **No matches**: lower `-t` (e.g. `-t 0.5`) or use `--always-top-n` to force the
  best N per cluster.
- **Reference not found**: run `mhc-tp fetch -s <species>` first, or pass
  `-r path/to/<species>.parquet`.
- **`gibbs_folder` rejected**: ensure it contains a `matrices/` subdirectory.

## Citation

Munday PR, Krishna SSG, Fehring J, Croft NP, Purcell AW, Li C, Braun A.
*Immunolyser 2.0: An advanced computational pipeline for comprehensive analysis of
immunopeptidomic data.* Comput Struct Biotechnol J. 2025;29:296–304.
doi:10.1016/j.csbj.2025.10.007. PMID: 41209766.
