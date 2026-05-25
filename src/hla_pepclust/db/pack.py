"""DEV-ONLY: ingest NetMHCpan reference packs into per-species parquets."""
from __future__ import annotations
import re
from dataclasses import dataclass
from pathlib import Path
import pandas as pd
from hla_pepclust.io.matrices import parse_matrix
from hla_pepclust.refdata.parquet_io import write_reference


@dataclass(frozen=True)
class PackAllotype:
    allotype: str
    formatted: str
    mhc_class: str   # "I" | "II"
    locus: str
    species: str     # "human" | "mouse" | "other"


def classify_allele(name: str) -> PackAllotype | None:
    """Classify a pack allele name; return None if not human/mouse."""
    n = name.strip()
    up = n.upper()
    # mouse: H-2-/H2- prefix
    if up.startswith(("H-2-", "H2-", "H-2", "H2_")):
        core = re.sub(r"^H-?2[-_]?", "", n)              # e.g. Db, Kb, IAb, IEk
        if core[:2].upper() in ("IA", "IE"):
            return PackAllotype(n, "H2" + re.sub(r"[^A-Za-z0-9]", "", core), "II", core[:2].upper(), "mouse")
        if core[:1].upper() in ("K", "D", "L", "Q"):
            return PackAllotype(n, "H2" + re.sub(r"[^A-Za-z0-9]", "", core), "I", core[:1].upper(), "mouse")
        return None
    # human HLA- prefixed
    if up.startswith("HLA-"):
        body = n[4:]
        if body[:1] in ("A", "B", "C", "E", "G") and (len(body) < 2 or not body[1].isalpha()):
            return PackAllotype(n, "HLA" + re.sub(r"[^A-Za-z0-9]", "", body), "I", body[0], "human")
        loc = body[:2].upper()
        if loc in ("DR", "DQ", "DP"):
            return PackAllotype(n, "HLA" + re.sub(r"[^A-Za-z0-9]", "", body), "II", loc, "human")
        return None
    # bare class II human (e.g. DRB1_0101)
    if up[:2] in ("DR", "DQ", "DP"):
        return PackAllotype(n, re.sub(r"[^A-Za-z0-9]", "", n), "II", up[:2], "human")
    return None


def _pairs_class_i(base: Path):
    for nf in (base / "mhc_names").glob("*_names.txt"):
        for line in nf.read_text().splitlines():
            parts = line.split()
            if len(parts) >= 2:
                yield parts[0], parts[1]   # (allele, pseudo)


def _pairs_class_ii(base: Path):
    for line in (base / "pseudo_mhc_list").read_text().splitlines():
        parts = line.split()
        if len(parts) >= 2:
            yield parts[1], parts[0]       # (allele, pseudo)


def build_pack_parquet(pack_dir, mhc_class, species, out_parquet, source):
    """Build a parquet from a pack for one species. mhc_class in {"I","II"}.

    Uses the LOG-ODDS (halfbits) score matrices so the reference is in the same
    representation as GibbsCluster output (the search query) — class I
    ``score_mat_el/<pseudo>.txt``, class II ``log_odds/<pseudo>.txt``. The
    frequency matrices are for logos, not correlation.
    """
    base = Path(pack_dir) / "all_logos"
    if mhc_class == "I":
        mat_dir = base / "score_mat_el"
        pairs = _pairs_class_i(base)
    else:
        mat_dir = base / "log_odds"
        pairs = _pairs_class_ii(base)
    rows, seen = [], set()
    for allele, pseudo in pairs:
        info = classify_allele(allele)
        if info is None or info.species != species or info.mhc_class != mhc_class:
            continue
        if info.formatted in seen:
            continue
        mat = parse_matrix(mat_dir / (pseudo + ".txt"))
        if mat is None:
            continue
        seen.add(info.formatted)
        rows.append({
            "allotype": info.allotype, "formatted": info.formatted,
            "mhc_class": info.mhc_class, "locus": info.locus,
            "n_positions": int(mat.shape[0]), "matrix": mat.reshape(-1).tolist(),
            "source": source,
        })
    write_reference(pd.DataFrame(rows), out_parquet)
    return len(rows)
