"""Edge-case tests for the io layer.

Covers:
* ``hla_pepclust.io.matrices.parse_matrix``
* ``hla_pepclust.io.kld.read_kld``
* ``hla_pepclust.io.naming.format_allotype``

These complement the happy-path tests in ``test_io_matrices.py``,
``test_io_kld.py`` and ``test_io_naming.py`` by hammering the boundary
conditions: empty/garbage input, ragged rows, sign preservation, etc.
"""

import textwrap

import numpy as np
import pytest

from hla_pepclust.io.kld import read_kld
from hla_pepclust.io.matrices import parse_matrix
from hla_pepclust.io.naming import format_allotype

HEADER = "Pos A R N D C Q E G H I L K M F P S T W Y V"


def _zeros20():
    return ["0.0"] * 20


# --------------------------------------------------------------------------- #
# parse_matrix
# --------------------------------------------------------------------------- #


def test_parse_matrix_missing_file_returns_none(tmp_path):
    assert parse_matrix(tmp_path / "does_not_exist.txt") is None


def test_parse_matrix_empty_file_returns_none(tmp_path):
    p = tmp_path / "empty.txt"
    p.write_text("")
    assert parse_matrix(p) is None


def test_parse_matrix_only_comments_and_header_returns_none(tmp_path):
    p = tmp_path / "headeronly.txt"
    p.write_text(textwrap.dedent(f"""\
            # a comment
            # another comment
            {HEADER}
            """))
    assert parse_matrix(p) is None


def test_parse_matrix_only_blank_lines_returns_none(tmp_path):
    p = tmp_path / "blanks.txt"
    p.write_text("\n\n   \n\t\n")
    assert parse_matrix(p) is None


def test_parse_matrix_rows_with_too_few_fields_skipped(tmp_path):
    # First "row" only has 19 numeric fields -> skipped; second is valid.
    short = "1 P " + " ".join(["0.0"] * 19)
    good_fields = _zeros20()
    good_fields[0] = "0.3"
    good = "2 P " + " ".join(good_fields)
    p = tmp_path / "short.txt"
    p.write_text(f"# c\n{HEADER}\n{short}\n{good}\n")
    m = parse_matrix(p)
    assert m is not None
    assert m.shape == (1, 20)
    assert m[0, 0] == np.float32(0.3)


def test_parse_matrix_valid_two_positions_shape_dtype_values(tmp_path):
    # Use distinct values in the trailing 20 columns to verify ordering.
    row1 = _zeros20()
    row1[0] = "0.1"  # A
    row1[19] = "0.7"  # V
    row2 = _zeros20()
    row2[9] = "0.5"  # I
    row2[10] = "0.5"  # L
    p = tmp_path / "valid.txt"
    p.write_text(
        f"# comment\n{HEADER}\n" f"1 P {' '.join(row1)}\n" f"2 P {' '.join(row2)}\n"
    )
    m = parse_matrix(p)
    assert m is not None
    assert m.shape == (2, 20)
    assert m.dtype == np.float32
    assert m[0, 0] == np.float32(0.1)
    assert m[0, 19] == np.float32(0.7)
    assert m[1, 9] == np.float32(0.5)
    assert m[1, 10] == np.float32(0.5)


def test_parse_matrix_uses_last_20_columns(tmp_path):
    # Many leading metadata columns; only the LAST 20 fields are scores.
    scores = [str(float(i)) for i in range(20)]  # 0.0 .. 19.0
    leading = "Pos PepLen SomeJunk MoreJunk"
    p = tmp_path / "leading.txt"
    p.write_text(f"{HEADER}\n{leading} {' '.join(scores)}\n")
    m = parse_matrix(p)
    assert m is not None
    assert m.shape == (1, 20)
    np.testing.assert_array_equal(
        m[0], np.array([float(i) for i in range(20)], dtype=np.float32)
    )


def test_parse_matrix_nonnumeric_trailing_field_skips_row(tmp_path):
    # A row whose trailing-20 window contains junk should be skipped.
    bad = _zeros20()
    bad[5] = "NOTANUMBER"
    good = _zeros20()
    good[3] = "1.5"
    p = tmp_path / "junk.txt"
    p.write_text(f"{HEADER}\n" f"1 P {' '.join(bad)}\n" f"2 P {' '.join(good)}\n")
    m = parse_matrix(p)
    assert m is not None
    assert m.shape == (1, 20)
    assert m[0, 3] == np.float32(1.5)


def test_parse_matrix_all_rows_junk_returns_none(tmp_path):
    bad = _zeros20()
    bad[0] = "junk"
    p = tmp_path / "alljunk.txt"
    p.write_text(f"{HEADER}\n1 P {' '.join(bad)}\n")
    assert parse_matrix(p) is None


def test_parse_matrix_negative_and_zero_values_preserved(tmp_path):
    row = _zeros20()
    row[0] = "-3.5"
    row[1] = "0.0"
    row[2] = "-0.0"
    row[19] = "2.25"
    p = tmp_path / "neg.txt"
    p.write_text(f"{HEADER}\n1 P {' '.join(row)}\n")
    m = parse_matrix(p)
    assert m is not None
    assert m[0, 0] == np.float32(-3.5)
    assert m[0, 1] == np.float32(0.0)
    assert m[0, 2] == np.float32(0.0)  # -0.0 == 0.0
    assert m[0, 19] == np.float32(2.25)


def test_parse_matrix_tolerates_blank_lines_between_rows(tmp_path):
    r1 = _zeros20()
    r1[0] = "1.0"
    r2 = _zeros20()
    r2[0] = "2.0"
    p = tmp_path / "blanky.txt"
    p.write_text(
        f"{HEADER}\n\n" f"1 P {' '.join(r1)}\n" f"\n  \n" f"2 P {' '.join(r2)}\n\n"
    )
    m = parse_matrix(p)
    assert m is not None
    assert m.shape == (2, 20)
    assert m[0, 0] == np.float32(1.0)
    assert m[1, 0] == np.float32(2.0)


# --------------------------------------------------------------------------- #
# read_kld
# --------------------------------------------------------------------------- #


def test_read_kld_normal_file(tmp_path):
    f = tmp_path / "gibbs.KLDvsClusters.tab"
    f.write_text("Number\tg1\tg2\n1\t2.0\t0.5\n2\t1.0\t1.0\n")
    df = read_kld(f)
    assert list(df["cluster"]) == [1, 2]
    assert {"group1", "group2", "total"}.issubset(df.columns)
    assert df.loc[df["cluster"] == 1, "group1"].iloc[0] == 2.0
    assert df.loc[df["cluster"] == 1, "group2"].iloc[0] == 0.5
    assert df.loc[df["cluster"] == 1, "total"].iloc[0] == 2.5
    assert df.loc[df["cluster"] == 2, "total"].iloc[0] == 2.0


def test_read_kld_negative_values_clamped_to_zero(tmp_path):
    f = tmp_path / "neg.tab"
    f.write_text("Number\tg1\tg2\n1\t-3.0\t2.0\n")
    df = read_kld(f)
    assert df.loc[0, "group1"] == 0.0  # clamped from -3.0
    assert df.loc[0, "group2"] == 2.0
    assert df.loc[0, "total"] == 2.0


def test_read_kld_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_kld(tmp_path / "nope.tab")


def test_read_kld_trailing_blank_lines_tolerated(tmp_path):
    f = tmp_path / "trailing.tab"
    f.write_text("Number\tg1\tg2\n1\t2.0\t1.0\n2\t1.0\t0.5\n\n\n   \n")
    df = read_kld(f)
    assert list(df["cluster"]) == [1, 2]
    assert len(df) == 2


def test_read_kld_uniform_rows_total_is_row_sum(tmp_path):
    f = tmp_path / "uniform.tab"
    f.write_text("Number\tg1\tg2\tg3\n" "1\t1.0\t2.0\t3.0\n" "2\t0.0\t0.0\t0.0\n")
    df = read_kld(f)
    assert df.loc[df["cluster"] == 1, "total"].iloc[0] == 6.0
    assert df.loc[df["cluster"] == 2, "total"].iloc[0] == 0.0
    assert {"group1", "group2", "group3"}.issubset(df.columns)


def test_read_kld_ragged_rows_filled_with_zero(tmp_path):
    # Row 1 has 2 groups, row 2 has only 1 group.
    f = tmp_path / "ragged.tab"
    f.write_text("Number\tg1\tg2\n1\t2.0\t1.0\n2\t1.5\n")
    df = read_kld(f)
    # Expected (per spec): the missing group2 cell is filled with 0.0.
    assert list(df["cluster"]) == [1, 2]
    assert df.loc[df["cluster"] == 2, "group2"].iloc[0] == 0.0
    assert df.loc[df["cluster"] == 2, "total"].iloc[0] == 1.5


# --------------------------------------------------------------------------- #
# format_allotype
# --------------------------------------------------------------------------- #


def test_format_allotype_human_class_i():
    info = format_allotype("HLA-A*02:01", species="human")
    assert info.raw == "HLA-A*02:01"
    assert info.formatted == "A0201"
    assert info.mhc_class == "I"
    assert info.locus == "A"


def test_format_allotype_human_class_ii():
    info = format_allotype("HLA-DRB1*15:01", species="human")
    assert info.mhc_class == "II"
    assert info.locus == "DR"
    assert info.formatted == "DRB11501"


def test_format_allotype_mouse_h2dash():
    info = format_allotype("H2-Kb", species="mouse")
    assert info.mhc_class == "I"
    assert info.locus == "K"
    assert info.formatted == "H2Kb"


def test_format_allotype_mouse_hdash2dash():
    # "H-2-Db": "H-2" -> "H2" then "H2-" -> "H2" leaves "H2Db".
    info = format_allotype("H-2-Db", species="mouse")
    assert info.mhc_class == "I"
    assert info.locus == "D"
    assert info.formatted == "H2Db"


def test_format_allotype_already_bare_name():
    # No prefixes/separators to strip.
    info = format_allotype("B0702", species="human")
    assert info.formatted == "B0702"
    assert info.mhc_class == "I"
    assert info.locus == "B"


def test_format_allotype_lowercase_input_not_normalized():
    # The formatter only strips the uppercase "HLA-" token, so lowercase
    # input is passed through unchanged (documents current behaviour, which
    # leaves the prefix in and lowercases the locus).
    info = format_allotype("hla-a*02:01", species="human")
    assert info.formatted == "hlaa0201"
    assert info.locus == "h"
    assert info.mhc_class == "I"


def test_format_allotype_garbage_does_not_crash():
    info = format_allotype("garbage123", species="human")
    assert info.raw == "garbage123"
    assert info.formatted == "garbage123"
    assert info.locus == "g"
    assert info.mhc_class == "I"


def test_format_allotype_empty_string_does_not_crash():
    info = format_allotype("", species="human")
    assert info.raw == ""
    assert info.formatted == ""
    assert info.locus == ""
    assert info.mhc_class == "I"


def test_format_allotype_mouse_class_ii():
    assert format_allotype("H2-IAb", species="mouse").locus == "IA"
    assert format_allotype("H2-IAb", species="mouse").mhc_class == "II"
    assert format_allotype("H2-IEd", species="mouse").locus == "IE"
