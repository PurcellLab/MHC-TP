# MHC-TP

Cluster immunopeptidomics peptides by their **HLA/MHC binding motif** and get a
ranked table plus a standalone interactive HTML report.

`mhc-tp` takes a **GibbsCluster** output folder, correlates each cluster's
position-specific scoring matrix against a reference of HLA/MHC **class I + II**
binding motifs (human & mouse), and writes the best allele match per cluster.

[Get started](usage.md){ .md-button .md-button--primary }
[View an example report](example-report.html){ .md-button }

## Install

```bash
git clone https://github.com/PurcellLab/MHC-TP.git
cd MHC-TP
pip install -e .
```

Then fetch the reference motifs (once) and run a search:

```bash
mhc-tp fetch -s human
mhc-tp search <gibbscluster_output_dir> -s human -o results/
```

See the [Usage](usage.md) guide for all options and the
[API reference](api.md) for the search engine internals (including the
correlation maths).

## How it works

For each GibbsCluster motif, every reference allotype motif is scored by the
**Pearson correlation** of their flattened position-weight matrices, computed
only over the informative cells of the cluster motif. Per cluster the allotypes
are ranked by this correlation (PCC); see
[`mhc_tp.engine.search`](api.md#mhc_tp.engine.search) for the full method.

## Citation

If you use MHC-TP, please cite:

> Munday PR, Krishna SSG, Fehring J, Croft NP, Purcell AW, Li C, Braun A.
> *Immunolyser 2.0: An advanced computational pipeline for comprehensive
> analysis of immunopeptidomic data.* Comput Struct Biotechnol J. 2025;29:296–304.
> doi:[10.1016/j.csbj.2025.10.007](https://doi.org/10.1016/j.csbj.2025.10.007).
