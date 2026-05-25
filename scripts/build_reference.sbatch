#!/bin/bash
# Build the class I+II reference parquets WITH embedded Seq2Logo logos.
# This is the slow (~hours) one-time release build: ~21k Seq2Logo renders.
#
# Submit with your scheduler account/partition (SLURM directives cannot read
# env vars, so pass them on the command line):
#   sbatch --account=<HPC_ACCOUNT> --partition=<HPC_PARTITION> scripts/build_reference.sbatch
#
# Prerequisites:
#   - export SEQ2LOGO_PATH=/path/to/seq2logo-2.1   (your own licensed install)
#   - pixi env 'seq2logo' installed:  pixi install -e seq2logo
#   - packs extracted:  bash scripts/fetch_sources.sh
#
#SBATCH --job-name=hlapep-refbuild
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=slurm-refbuild-%j.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

: "${SEQ2LOGO_PATH:?export SEQ2LOGO_PATH to your seq2logo-2.1 directory}"
PY27="$PWD/.pixi/envs/seq2logo/bin/python"
C1="data/ref_sources/extracted/netmhcpan42_all_logos"
C2="data/ref_sources/extracted/netmhciipan43_all_logos"
OUT="${OUT_DIR:-data/ref_data}"
mkdir -p "$OUT"

for SP in human mouse; do
    echo ">>> building ${SP}.parquet (--with-logos)"
    pixi run clust-search build-ref "$SP" "$C1" "$C2" "$OUT/${SP}.parquet" \
        --with-logos --seq2logo-path "$SEQ2LOGO_PATH" --seq2logo-python "$PY27"
done

echo ">>> done. recompute these into src/hla_pepclust/refdata/reference_manifest.tsv:"
sha256sum "$OUT"/human.parquet "$OUT"/mouse.parquet
