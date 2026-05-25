"""Resolve and fetch prebuilt reference parquets to a per-user data dir.

End users run ``clust-search fetch`` to download the prebuilt class I+II
reference parquets (with embedded Seq2Logo reference logos) instead of building
them. The download source + checksums live in the packaged
``reference_manifest.tsv``; the maintainer fills them in on each release.
"""

from __future__ import annotations

import hashlib
import os
import urllib.request
from importlib.resources import files
from pathlib import Path

import platformdirs

_APP = "hla_pepclust"


def data_dir() -> Path:
    """User data dir for reference files. Overridable via ``HLA_PEPCLUST_DATA_DIR``."""
    override = os.environ.get("HLA_PEPCLUST_DATA_DIR")
    return Path(override) if override else Path(platformdirs.user_data_dir(_APP))


def reference_path(species: str) -> Path:
    """Expected path of a species reference parquet in the data dir."""
    return data_dir() / f"{species.lower()}.parquet"


def resolve_reference(species: str, override: str | None = None) -> Path:
    """Return the reference parquet path, raising a helpful error if absent."""
    if override:
        p = Path(override)
        if not p.exists():
            raise FileNotFoundError(f"reference not found: {p}")
        return p
    p = reference_path(species)
    if not p.exists():
        raise FileNotFoundError(
            f"No {species} reference at {p}. Run: clust-search fetch --species {species} "
            f"(or pass --reference <path>)."
        )
    return p


def load_manifest() -> list[dict]:
    """Parse the packaged reference manifest (species/filename/sha256/url)."""
    text = files("hla_pepclust.refdata").joinpath("reference_manifest.tsv").read_text()
    rows = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("species\t"):
            continue
        species, filename, sha256, url = line.split("\t")
        rows.append(
            {"species": species, "filename": filename, "sha256": sha256, "url": url}
        )
    return rows


def fetch(species: str = "all", dest: str | None = None) -> list[Path]:
    """Download the reference parquet(s) into the data dir; verify checksums."""
    out_dir = Path(dest) if dest else data_dir()
    out_dir.mkdir(parents=True, exist_ok=True)
    fetched = []
    for row in load_manifest():
        if species != "all" and row["species"] != species:
            continue
        if not row["url"] or row["url"] == "-":
            raise RuntimeError(
                f"No download URL configured for the {row['species']} reference yet. "
                f"Build it locally with `clust-search build-ref`, or point --reference "
                f"at a parquet."
            )
        target = out_dir / row["filename"]
        urllib.request.urlretrieve(
            row["url"], target
        )  # noqa: S310 (trusted manifest URL)
        expected = row["sha256"]
        if expected and expected != "-":
            actual = hashlib.sha256(target.read_bytes()).hexdigest()
            if actual != expected:
                target.unlink(missing_ok=True)
                raise RuntimeError(f"checksum mismatch for {row['filename']}")
        fetched.append(target)
    return fetched
