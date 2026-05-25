"""DEV-ONLY: ingest NetMHCpan reference packs into per-species parquets."""

from __future__ import annotations
import re
from dataclasses import dataclass
from pathlib import Path
import pandas as pd
from hla_pepclust.io.matrices import parse_matrix
from hla_pepclust.refdata.parquet_io import read_reference, write_reference


@dataclass(frozen=True)
class PackAllotype:
    allotype: str
    formatted: str
    mhc_class: str  # "I" | "II"
    locus: str
    species: str  # "human" | "mouse" | "other"


def classify_allele(name: str) -> PackAllotype | None:
    """Classify a pack allele name; return None if not human/mouse."""
    n = name.strip()
    up = n.upper()
    # mouse: H-2-/H2- prefix
    if up.startswith(("H-2-", "H2-", "H-2", "H2_")):
        core = re.sub(r"^H-?2[-_]?", "", n)  # e.g. Db, Kb, IAb, IEk
        if core[:2].upper() in ("IA", "IE"):
            return PackAllotype(
                n,
                "H2" + re.sub(r"[^A-Za-z0-9]", "", core),
                "II",
                core[:2].upper(),
                "mouse",
            )
        if core[:1].upper() in ("K", "D", "L", "Q"):
            return PackAllotype(
                n,
                "H2" + re.sub(r"[^A-Za-z0-9]", "", core),
                "I",
                core[:1].upper(),
                "mouse",
            )
        return None
    # human HLA- prefixed
    if up.startswith("HLA-"):
        body = n[4:]
        if body[:1] in ("A", "B", "C", "E", "G") and (
            len(body) < 2 or not body[1].isalpha()
        ):
            return PackAllotype(
                n, "HLA" + re.sub(r"[^A-Za-z0-9]", "", body), "I", body[0], "human"
            )
        loc = body[:2].upper()
        if loc in ("DR", "DQ", "DP"):
            return PackAllotype(
                n, "HLA" + re.sub(r"[^A-Za-z0-9]", "", body), "II", loc, "human"
            )
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
                yield parts[0], parts[1]  # (allele, pseudo)


def _pairs_class_ii(base: Path):
    for line in (base / "pseudo_mhc_list").read_text().splitlines():
        parts = line.split()
        if len(parts) >= 2:
            yield parts[1], parts[0]  # (allele, pseudo)


def build_pack_parquet(
    pack_dir,
    mhc_class,
    species,
    out_parquet,
    source,
    with_logos=False,
    seq2logo_path=None,
    seq2logo_python=None,
):
    """Build a parquet from a pack for one species. mhc_class in {"I","II"}.

    Correlation ``matrix`` comes from the LOG-ODDS (halfbits) score matrices so
    the reference matches GibbsCluster output — class I ``score_mat_el/<pseudo>.txt``,
    class II ``log_odds/<pseudo>.txt``. When ``with_logos`` is set, a Seq2Logo
    reference logo is rendered from the FREQUENCY matrix (class I
    ``freq_mat_el/<pseudo>_freq.mat``, class II ``freq_mat/<pseudo>_freq.mat``)
    and embedded in a ``logo`` column.
    """
    base = Path(pack_dir) / "all_logos"
    if mhc_class == "I":
        mat_dir, freq_dir = base / "score_mat_el", base / "freq_mat_el"
        pairs = _pairs_class_i(base)
    else:
        mat_dir, freq_dir = base / "log_odds", base / "freq_mat"
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
        row = {
            "allotype": info.allotype,
            "formatted": info.formatted,
            "mhc_class": info.mhc_class,
            "locus": info.locus,
            "n_positions": int(mat.shape[0]),
            "matrix": mat.reshape(-1).tolist(),
            "source": source,
        }
        if with_logos:
            from hla_pepclust.db.logos import reference_logo_bytes

            row["logo"] = reference_logo_bytes(
                freq_dir / (pseudo + "_freq.mat"),
                seq2logo_path=seq2logo_path,
                python_exe=seq2logo_python,
                title=info.allotype,
            )
        rows.append(row)
    write_reference(pd.DataFrame(rows), out_parquet)
    return len(rows)


def build_species_reference(
    species,
    class_i_pack,
    class_ii_pack,
    out_parquet,
    source_i="NetMHCpan-4.2",
    source_ii="NetMHCIIpan-4.3",
    with_logos=False,
    seq2logo_path=None,
    seq2logo_python=None,
):
    """Build one ``<species>.parquet`` combining class I + class II from the packs.

    Returns (n_class_i, n_class_ii). Pass ``with_logos`` (+ Seq2Logo paths) to
    embed reference logos.
    """
    import tempfile

    tmp = Path(tempfile.mkdtemp())
    pi, pii = tmp / "i.parquet", tmp / "ii.parquet"
    logo_kw = dict(
        with_logos=with_logos,
        seq2logo_path=seq2logo_path,
        seq2logo_python=seq2logo_python,
    )
    n_i = build_pack_parquet(class_i_pack, "I", species, pi, source_i, **logo_kw)
    n_ii = build_pack_parquet(class_ii_pack, "II", species, pii, source_ii, **logo_kw)
    frames = []
    if n_i:
        frames.append(read_reference(pi))
    if n_ii:
        frames.append(read_reference(pii))
    combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    write_reference(combined, out_parquet)
    return n_i, n_ii
