"""Edge-case tests for the refdata + db.pack layers.

Covers:
  * parquet_io  -> write_reference / read_reference / load_logos
  * fetch       -> data_dir / reference_path / resolve_reference / load_manifest
  * export      -> export_logos
  * db.pack     -> classify_allele / build_pack_parquet / build_species_reference

All fixtures use tmp_path + synthetic data only; the real
``data/ref_data/*.parquet`` files are never read (a build job writes there).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from mhc_tp.constants import N_AMINO_ACIDS
from mhc_tp.db.pack import (
    build_pack_parquet,
    build_species_reference,
    classify_allele,
)
from mhc_tp.refdata.export import export_logos
from mhc_tp.refdata.fetch import (
    data_dir,
    load_manifest,
    reference_path,
    resolve_reference,
)
from mhc_tp.refdata.parquet_io import (
    load_logos,
    read_reference,
    write_reference,
)
from mhc_tp.refdata.schema import COLUMNS

# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

_MAT_HEADER = (
    "#cmd\n"
    "Last position-specific scoring matrix computed, values are in halfbits\n"
    "A R N D C Q E G H I L K M F P S T W Y V\n"
)


def _score_rows(n_pos: int, val: str = "0.5") -> str:
    """n_pos data rows in the NetMHCpan/Gibbs matrix layout."""
    return "".join(
        f"{i + 1} A " + " ".join([val] * N_AMINO_ACIDS) + "\n" for i in range(n_pos)
    )


def _matrix_file(n_pos: int = 2, val: str = "0.5") -> str:
    return _MAT_HEADER + _score_rows(n_pos, val)


def _ref_df(n_rows: int = 1, *, with_logo: bool = False) -> pd.DataFrame:
    """A minimal valid reference DataFrame (all required COLUMNS present)."""
    rows = []
    for i in range(n_rows):
        n_pos = 2
        matrix = list(np.arange(n_pos * N_AMINO_ACIDS, dtype=np.float32))
        rows.append(
            {
                "allotype": f"A*02:{i:02d}",
                "formatted": f"A02{i:02d}",
                "mhc_class": "I",
                "locus": "A",
                "n_positions": n_pos,
                "matrix": matrix,
                "source": "synthetic",
            }
        )
    df = pd.DataFrame(rows)
    if with_logo:
        df["logo"] = [f"PNGBYTES-{r.formatted}".encode() for r in df.itertuples()]
    return df


# --------------------------------------------------------------------------- #
# parquet_io: write / read roundtrip
# --------------------------------------------------------------------------- #


def test_write_read_roundtrip_preserves_matrix_and_columns(tmp_path):
    df = _ref_df(2)
    path = tmp_path / "ref.parquet"
    write_reference(df, path)

    back = read_reference(path)
    assert list(back.columns) == list(COLUMNS)  # no logo column
    assert len(back) == 2
    # matrix list<float32> survives the roundtrip element-for-element
    for orig, got in zip(df["matrix"], back["matrix"]):
        assert list(orig) == list(got)
    assert list(back["formatted"]) == list(df["formatted"])


def test_write_creates_parent_dirs(tmp_path):
    df = _ref_df(1)
    path = tmp_path / "nested" / "deep" / "ref.parquet"
    write_reference(df, path)
    assert path.exists()


@pytest.mark.parametrize("missing_col", COLUMNS)
def test_write_rejects_missing_required_column(tmp_path, missing_col):
    df = _ref_df(1).drop(columns=[missing_col])
    with pytest.raises(ValueError, match="missing columns"):
        write_reference(df, tmp_path / "bad.parquet")


def test_write_extra_columns_are_dropped(tmp_path):
    df = _ref_df(1)
    df["extra_junk"] = ["nope"]
    path = tmp_path / "ref.parquet"
    write_reference(df, path)
    back = read_reference(path)
    assert "extra_junk" not in back.columns


# --------------------------------------------------------------------------- #
# parquet_io: logo column handling
# --------------------------------------------------------------------------- #


def test_write_read_roundtrip_with_logo_column(tmp_path):
    df = _ref_df(2, with_logo=True)
    path = tmp_path / "ref.parquet"
    write_reference(df, path)

    back = read_reference(path)
    assert "logo" in back.columns
    assert list(back.columns) == list(COLUMNS) + ["logo"]
    assert list(back["logo"]) == list(df["logo"])


def test_read_with_columns_subset_skips_logo(tmp_path):
    df = _ref_df(2, with_logo=True)
    path = tmp_path / "ref.parquet"
    write_reference(df, path)

    back = read_reference(path, columns=list(COLUMNS))
    assert "logo" not in back.columns
    assert list(back.columns) == list(COLUMNS)


def test_read_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_reference(tmp_path / "does_not_exist.parquet")


# --------------------------------------------------------------------------- #
# parquet_io: load_logos predicate pushdown
# --------------------------------------------------------------------------- #


def test_load_logos_returns_only_requested_and_present(tmp_path):
    df = _ref_df(3, with_logo=True)  # formatted: A0200, A0201, A0202
    path = tmp_path / "ref.parquet"
    write_reference(df, path)

    got = load_logos(path, ["A0200", "A0202", "NOTHERE"])
    assert set(got) == {"A0200", "A0202"}
    assert got["A0200"] == b"PNGBYTES-A0200"
    assert got["A0202"] == b"PNGBYTES-A0202"


def test_load_logos_no_logo_column_returns_empty(tmp_path):
    df = _ref_df(2)  # no logo column
    path = tmp_path / "ref.parquet"
    write_reference(df, path)
    assert load_logos(path, ["A0200", "A0201"]) == {}


def test_load_logos_empty_ids_returns_empty(tmp_path):
    df = _ref_df(2, with_logo=True)
    path = tmp_path / "ref.parquet"
    write_reference(df, path)
    assert load_logos(path, []) == {}


def test_load_logos_ids_not_present_returns_empty(tmp_path):
    df = _ref_df(2, with_logo=True)
    path = tmp_path / "ref.parquet"
    write_reference(df, path)
    assert load_logos(path, ["ZZZ999", "QQQ000"]) == {}


def test_load_logos_missing_file_returns_empty(tmp_path):
    assert load_logos(tmp_path / "nope.parquet", ["A0200"]) == {}


def test_load_logos_dedups_requested_ids(tmp_path):
    df = _ref_df(2, with_logo=True)
    path = tmp_path / "ref.parquet"
    write_reference(df, path)
    got = load_logos(path, ["A0200", "A0200", "A0200"])
    assert got == {"A0200": b"PNGBYTES-A0200"}


def test_load_logos_skips_empty_blobs(tmp_path):
    # A row whose logo blob is empty bytes must be filtered out.
    df = _ref_df(2, with_logo=True)
    df.loc[0, "logo"] = b""
    path = tmp_path / "ref.parquet"
    write_reference(df, path)
    got = load_logos(path, ["A0200", "A0201"])
    assert "A0200" not in got
    assert got["A0201"] == b"PNGBYTES-A0201"


# --------------------------------------------------------------------------- #
# fetch
# --------------------------------------------------------------------------- #


def test_data_dir_honours_env_override(tmp_path, monkeypatch):
    monkeypatch.setenv("MHC_TP_DATA_DIR", str(tmp_path))
    assert data_dir() == tmp_path


def test_data_dir_default_without_env(monkeypatch):
    monkeypatch.delenv("MHC_TP_DATA_DIR", raising=False)
    d = data_dir()
    assert "mhc_tp" in str(d)


def test_reference_path_lowercases_species(tmp_path, monkeypatch):
    monkeypatch.setenv("MHC_TP_DATA_DIR", str(tmp_path))
    assert reference_path("HUMAN") == tmp_path / "human.parquet"
    assert reference_path("Mouse") == tmp_path / "mouse.parquet"


def test_resolve_reference_override_existing_returned(tmp_path):
    p = tmp_path / "custom.parquet"
    p.write_bytes(b"x")
    assert resolve_reference("human", str(p)) == p


def test_resolve_reference_override_missing_raises(tmp_path):
    missing = tmp_path / "missing.parquet"
    with pytest.raises(FileNotFoundError, match="reference not found"):
        resolve_reference("human", str(missing))


def test_resolve_reference_present_in_data_dir(tmp_path, monkeypatch):
    monkeypatch.setenv("MHC_TP_DATA_DIR", str(tmp_path))
    (tmp_path / "human.parquet").write_bytes(b"x")
    assert resolve_reference("human") == tmp_path / "human.parquet"


def test_resolve_reference_missing_mentions_fetch(tmp_path, monkeypatch):
    monkeypatch.setenv("MHC_TP_DATA_DIR", str(tmp_path))
    with pytest.raises(FileNotFoundError, match="fetch"):
        resolve_reference("human")


def test_load_manifest_rows_have_expected_keys():
    rows = load_manifest()
    assert rows  # at least one row
    for r in rows:
        assert set(r) == {"species", "filename", "sha256", "url"}
        assert r["species"]
        assert r["filename"]
        assert r["sha256"]
        assert r["url"]


# --------------------------------------------------------------------------- #
# export_logos
# --------------------------------------------------------------------------- #


def test_export_logos_writes_png_per_row(tmp_path):
    df = _ref_df(3, with_logo=True)
    parquet = tmp_path / "ref.parquet"
    write_reference(df, parquet)

    out_dir = tmp_path / "out"
    n = export_logos(parquet, out_dir)
    assert n == 3
    for fmt, blob in zip(df["formatted"], df["logo"]):
        png = out_dir / f"{fmt}.png"
        assert png.exists()
        assert png.read_bytes() == blob


def test_export_logos_no_logo_column_returns_zero(tmp_path):
    df = _ref_df(2)  # no logo column
    parquet = tmp_path / "ref.parquet"
    write_reference(df, parquet)

    out_dir = tmp_path / "out"
    n = export_logos(parquet, out_dir)
    assert n == 0
    # Nothing written (dir not even created, and certainly no pngs).
    assert not out_dir.exists() or not list(out_dir.glob("*.png"))


def test_export_logos_skips_empty_blob_rows(tmp_path):
    df = _ref_df(2, with_logo=True)
    df.loc[0, "logo"] = b""
    parquet = tmp_path / "ref.parquet"
    write_reference(df, parquet)

    out_dir = tmp_path / "out"
    n = export_logos(parquet, out_dir)
    assert n == 1
    assert not (out_dir / "A0200.png").exists()
    assert (out_dir / "A0201.png").exists()


# --------------------------------------------------------------------------- #
# classify_allele
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "name,locus,formatted",
    [
        ("HLA-A02:01", "A", "HLAA0201"),
        ("HLA-B07:02", "B", "HLAB0702"),
        ("HLA-C07:01", "C", "HLAC0701"),
        ("HLA-E01:01", "E", "HLAE0101"),
        ("HLA-G01:01", "G", "HLAG0101"),
    ],
)
def test_classify_human_class_i(name, locus, formatted):
    a = classify_allele(name)
    assert a is not None
    assert a.mhc_class == "I"
    assert a.locus == locus
    assert a.species == "human"
    assert a.formatted == formatted


@pytest.mark.parametrize(
    "name,locus",
    [
        ("DRB1_0101", "DR"),
        ("DQA1_0501", "DQ"),
        ("DPB1_0401", "DP"),
        ("HLA-DRB1*01:01", "DR"),
        ("HLA-DQB1*06:02", "DQ"),
        ("HLA-DPA1*01:03", "DP"),
    ],
)
def test_classify_human_class_ii(name, locus):
    a = classify_allele(name)
    assert a is not None
    assert a.mhc_class == "II"
    assert a.locus == locus
    assert a.species == "human"


@pytest.mark.parametrize(
    "name,locus",
    [
        ("H-2-Db", "D"),
        ("H-2-Kb", "K"),
        ("H-2-Ld", "L"),
        ("H-2-Qa1", "Q"),
        ("H2-Kb", "K"),
    ],
)
def test_classify_mouse_class_i(name, locus):
    a = classify_allele(name)
    assert a is not None
    assert a.mhc_class == "I"
    assert a.locus == locus
    assert a.species == "mouse"
    assert a.formatted.startswith("H2")


@pytest.mark.parametrize(
    "name,locus",
    [
        ("H-2-IAb", "IA"),
        ("H-2-IEk", "IE"),
        ("H2-IAd", "IA"),
    ],
)
def test_classify_mouse_class_ii(name, locus):
    a = classify_allele(name)
    assert a is not None
    assert a.mhc_class == "II"
    assert a.locus == locus
    assert a.species == "mouse"


@pytest.mark.parametrize(
    "name",
    [
        "Mamu-A01:01",
        "BoLA-DRB3_00101",
        "Patr-A01:01",
        "SLA-1*01:01",
        "",
        "   ",
        "garbage",
        "XYZ-123",
        "H-2-Zz",  # mouse prefix but unknown locus letter
        "HLA-X01:01",  # HLA but unknown class I locus
    ],
)
def test_classify_returns_none(name):
    assert classify_allele(name) is None


def test_classify_strips_whitespace():
    a = classify_allele("  HLA-A02:01  ")
    assert a is not None and a.locus == "A"


# --------------------------------------------------------------------------- #
# build_pack_parquet
# --------------------------------------------------------------------------- #


def _make_class_i_pack(base_dir, names_text, mat_files):
    """mat_files: {pseudo: matrix_text}."""
    base = base_dir / "all_logos"
    (base / "mhc_names").mkdir(parents=True)
    (base / "score_mat_el").mkdir()
    (base / "mhc_names" / "PS_names.txt").write_text(names_text)
    for pseudo, text in mat_files.items():
        (base / "score_mat_el" / f"{pseudo}.txt").write_text(text)
    return base_dir


def _make_class_ii_pack(base_dir, pseudo_list_text, mat_files):
    base = base_dir / "all_logos"
    (base / "log_odds").mkdir(parents=True)
    (base / "pseudo_mhc_list").write_text(pseudo_list_text)
    for pseudo, text in mat_files.items():
        (base / "log_odds" / f"{pseudo}.txt").write_text(text)
    return base_dir


def test_build_pack_class_i_basic(tmp_path):
    pack = _make_class_i_pack(
        tmp_path / "p",
        "HLA-A02:01 PS\n",
        {"PS": _matrix_file(2)},
    )
    out = tmp_path / "i.parquet"
    n = build_pack_parquet(pack, "I", "human", out, "NetMHCpan-4.2")
    assert n == 1
    d = read_reference(out)
    assert d["formatted"].iloc[0] == "HLAA0201"
    assert d["mhc_class"].iloc[0] == "I"
    assert d["n_positions"].iloc[0] == 2
    assert len(d["matrix"].iloc[0]) == 2 * N_AMINO_ACIDS
    assert "logo" not in d.columns  # with_logos=False -> no logo column


def test_build_pack_class_ii_basic(tmp_path):
    # pseudo_mhc_list lines are "<pseudo> <allele>".
    pack = _make_class_ii_pack(
        tmp_path / "p",
        "Q DRB1_0101\n",
        {"Q": _matrix_file(3)},
    )
    out = tmp_path / "ii.parquet"
    n = build_pack_parquet(pack, "II", "human", out, "NetMHCIIpan-4.3")
    assert n == 1
    d = read_reference(out)
    assert d["formatted"].iloc[0] == "DRB10101"
    assert d["mhc_class"].iloc[0] == "II"
    assert d["locus"].iloc[0] == "DR"
    assert d["n_positions"].iloc[0] == 3


def test_build_pack_species_filtering(tmp_path):
    # Mixed names: human + non-human; only the human row should be kept.
    pack = _make_class_i_pack(
        tmp_path / "p",
        "HLA-A02:01 PS1\nMamu-A01:01 PS2\nBoLA-x PS3\n",
        {"PS1": _matrix_file(2), "PS2": _matrix_file(2), "PS3": _matrix_file(2)},
    )
    out = tmp_path / "i.parquet"
    n = build_pack_parquet(pack, "I", "human", out, "src")
    assert n == 1
    assert read_reference(out)["formatted"].iloc[0] == "HLAA0201"


def test_build_pack_class_filtering(tmp_path):
    # A class II allele appearing in a class I pack listing must be skipped
    # when building class I.
    pack = _make_class_i_pack(
        tmp_path / "p",
        "HLA-A02:01 PS1\nDRB1_0101 PS2\n",
        {"PS1": _matrix_file(2), "PS2": _matrix_file(2)},
    )
    out = tmp_path / "i.parquet"
    n = build_pack_parquet(pack, "I", "human", out, "src")
    assert n == 1
    assert read_reference(out)["formatted"].iloc[0] == "HLAA0201"


def test_build_pack_dedups_identical_formatted(tmp_path):
    # Two listings collapse to the same formatted key -> only one row.
    pack = _make_class_i_pack(
        tmp_path / "p",
        "HLA-A02:01 PS1\nHLA-A02:01 PS2\n",
        {"PS1": _matrix_file(2), "PS2": _matrix_file(2)},
    )
    out = tmp_path / "i.parquet"
    n = build_pack_parquet(pack, "I", "human", out, "src")
    assert n == 1


def test_build_pack_missing_matrix_file_skips_row(tmp_path):
    # PS2 is referenced but its matrix file does not exist -> parse_matrix None
    # -> row skipped.
    pack = _make_class_i_pack(
        tmp_path / "p",
        "HLA-A02:01 PS1\nHLA-B07:02 PS2\n",
        {"PS1": _matrix_file(2)},  # PS2.txt deliberately absent
    )
    out = tmp_path / "i.parquet"
    n = build_pack_parquet(pack, "I", "human", out, "src")
    assert n == 1
    assert read_reference(out)["formatted"].iloc[0] == "HLAA0201"


def test_build_pack_no_matches_writes_empty(tmp_path):
    # All listings are non-human -> zero rows, but parquet still written.
    pack = _make_class_i_pack(
        tmp_path / "p",
        "Mamu-A01:01 PS1\n",
        {"PS1": _matrix_file(2)},
    )
    out = tmp_path / "i.parquet"
    n = build_pack_parquet(pack, "I", "human", out, "src")
    assert n == 0
    assert out.exists()


# --------------------------------------------------------------------------- #
# build_species_reference
# --------------------------------------------------------------------------- #


def test_build_species_reference_combines_class_i_and_ii(tmp_path):
    c1 = _make_class_i_pack(
        tmp_path / "c1",
        "HLA-A02:01 P\n",
        {"P": _matrix_file(2)},
    )
    c2 = _make_class_ii_pack(
        tmp_path / "c2",
        "Q DRB1_0101\n",
        {"Q": _matrix_file(3)},
    )
    out = tmp_path / "human.parquet"
    n_i, n_ii = build_species_reference("human", c1, c2, out)
    assert (n_i, n_ii) == (1, 1)
    d = read_reference(out)
    assert len(d) == 2
    assert set(d["mhc_class"]) == {"I", "II"}


def test_build_species_reference_class_i_only(tmp_path):
    c1 = _make_class_i_pack(
        tmp_path / "c1",
        "HLA-A02:01 P\n",
        {"P": _matrix_file(2)},
    )
    # class II pack: valid structure but no matching (human class II) rows.
    c2 = _make_class_ii_pack(
        tmp_path / "c2",
        "Q Mamu-A01:01\n",
        {"Q": _matrix_file(3)},
    )
    out = tmp_path / "human.parquet"
    n_i, n_ii = build_species_reference("human", c1, c2, out)
    assert (n_i, n_ii) == (1, 0)
    d = read_reference(out)
    assert len(d) == 1
    assert set(d["mhc_class"]) == {"I"}
