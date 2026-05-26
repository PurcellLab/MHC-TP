# Usage

## Install

**Requirements:** Python 3.9–3.11.

```bash
git clone https://github.com/PurcellLab/MHC-TP.git
cd MHC-TP
pip install -e .
```

> One-liner without cloning: `pip install git+https://github.com/PurcellLab/MHC-TP.git`
> A virtual environment is recommended.

## Download the reference data (once)

The reference motifs are fetched from the GitHub release, not bundled:

```bash
mhc-tp fetch -s human     # or:  mouse  |  all
```

## Run a search

```bash
mhc-tp search <gibbscluster_output_dir> -s human -o results/
```

`<gibbscluster_output_dir>` is a GibbsCluster run folder (it must contain a
`matrices/` subdirectory).

**Outputs** land in `results/clust_result/`:

| file | what it is |
|------|------------|
| `correlations.csv` | every cluster→allele match (`hla` = display name, `formatted` = raw key, `correlation` = PCC) |
| `mhc-tp-result.html` | standalone interactive report — open it in any browser |

## Options

| flag | meaning | default |
|------|---------|---------|
| `-s, --species` | `human` or `mouse` | `human` |
| `-c, --class` | restrict the reference to MHC class `I`, `II`, or `all` | `all` |
| `-r, --reference` | path to a `<species>.parquet` (otherwise the fetched one is used) | auto |
| `-t, --threshold` | minimum Pearson correlation to report | `0.70` |
| `--topNHits` | number of allotype matches to keep per cluster | `3` |
| `--always-top-n` | keep each cluster's top-N even if below `--threshold` (flagged in the report) | off |
| `-o, --output` | output directory | `output` |
| `--threads` | max CPU threads (also `$MHC_TP_THREADS`) | `4` |
| `--no-html` | write only the CSV | off |
| `-l, --log` | also save the coloured session log | off |

Run `mhc-tp search --help` for the full list.

### Examples

Restrict to MHC class I and keep the top 5 per cluster:

```bash
mhc-tp search runs/sampleA -s human -c I --topNHits 5 -o results/
```

Guarantee a top-3 for **every** cluster, even weak ones (sub-threshold matches
are tagged *below cutoff* in the report):

```bash
mhc-tp search runs/sampleA -s human --always-top-n -o results/
```

## How matches are selected

Per cluster, all reference allotypes are ranked by PCC. By default a hit must
score `>= --threshold`, so a cluster can yield fewer than `--topNHits` rows (or
none). With `--always-top-n`, every cluster returns its best `--topNHits`
regardless of threshold; the threshold then only annotates confidence. See the
[API reference](api.md#mhc_tp.engine.search) for the maths.
