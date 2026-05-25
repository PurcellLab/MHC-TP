#!/usr/bin/env bash
# Fetch the dev-only reference source packs (NetMHCpan-4.2 / NetMHCIIpan-4.3).
#
# Manifest + checksums live in data/ref_sources/sources.tsv (committed).
# Raw downloads land in data/ref_sources/raw/ and extracted archives in
# data/ref_sources/extracted/ — both gitignored (large third-party files).
# A file whose on-disk SHA-256 matches the manifest is skipped. A mismatch is
# fatal. A manifest sha256 of "-" means unpinned: download and PRINT the
# computed checksum (seed the manifest, then commit it).
#
# Usage:
#   scripts/fetch_sources.sh             # fetch + verify + extract all
#   scripts/fetch_sources.sh --no-extract
#   HLA_PEPCLUST_RECORD_CHECKSUMS=1 scripts/fetch_sources.sh   # print all checksums

set -Eeuo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC_DIR="${REPO_ROOT}/data/ref_sources"
RAW_DIR="${SRC_DIR}/raw"
EXTRACT_DIR="${SRC_DIR}/extracted"
MANIFEST="${SRC_DIR}/sources.tsv"

DO_EXTRACT=1
[[ "${1:-}" == "--no-extract" ]] && DO_EXTRACT=0
RECORD="${HLA_PEPCLUST_RECORD_CHECKSUMS:-0}"

[[ -f "${MANIFEST}" ]] || { echo "ERROR: manifest not found at ${MANIFEST}" >&2; exit 1; }
mkdir -p "${RAW_DIR}" "${EXTRACT_DIR}"
sha256_of() { sha256sum "$1" | cut -d' ' -f1; }

while IFS=$'\t' read -r name filename bytes sha256 extract url; do
    [[ -z "${name}" || "${name}" == \#* || "${name}" == "name" ]] && continue
    dest="${RAW_DIR}/${filename}"

    if [[ -f "${dest}" && "${sha256}" != "-" ]]; then
        if [[ "$(sha256_of "${dest}")" == "${sha256}" ]]; then
            echo "==> ${name}: up to date"; continue
        fi
        echo "==> ${name}: checksum drift, re-fetching"; rm -f "${dest}"
    fi

    if [[ ! -f "${dest}" ]]; then
        echo "==> ${name}: downloading"
        curl -sS --fail --location --max-time 3600 -o "${dest}" "${url}" \
            || { echo "    FAIL — download error for ${url}" >&2; exit 1; }
    fi

    have="$(sha256_of "${dest}")"
    if [[ "${sha256}" == "-" || "${RECORD}" == "1" ]]; then
        printf '    %s: bytes=%s sha256=%s  (seed the manifest)\n' "${name}" "$(stat -c%s "${dest}")" "${have}"
    elif [[ "${have}" != "${sha256}" ]]; then
        echo "    FAIL — SHA-256 mismatch for ${filename}: expected ${sha256}, got ${have}" >&2
        exit 1
    fi

    if [[ "${DO_EXTRACT}" == "1" && "${extract}" == "tar.gz" ]]; then
        target="${EXTRACT_DIR}/${name}"
        if [[ -d "${target}" ]]; then echo "    extracted present (${target})";
        else echo "    extracting -> ${target}"; mkdir -p "${target}"; tar -xzf "${dest}" -C "${target}"; fi
    fi
done < "${MANIFEST}"

echo; echo "done. raw=${RAW_DIR} (gitignored), extracted=${EXTRACT_DIR} (gitignored)"
