import numpy as np
import pandas as pd
import pytest

from hla_pepclust.refdata.export import export_logos
from hla_pepclust.refdata.parquet_io import write_reference


def _ref_with_logos(tmp_path, with_logo=True):
    row = {
        "allotype": "HLA-A*02:01",
        "formatted": "A0201",
        "mhc_class": "I",
        "locus": "A",
        "n_positions": 1,
        "matrix": np.zeros(20, np.float32).tolist(),
        "source": "test",
    }
    if with_logo:
        row["logo"] = b"\x89PNG\r\n\x1a\nfake-png-bytes"
    p = tmp_path / "ref.parquet"
    write_reference(pd.DataFrame([row]), p)
    return p


def _multi_ref(tmp_path):
    """A few rows mimicking the real parquet keys (formatted carries the prefix)."""
    rows = [
        ("HLA-B39:124", "HLAB39124", "I", "B"),
        ("HLA-B39:42", "HLAB3942", "I", "B"),
        ("H-2-Kb", "H2Kb", "I", "K"),
    ]
    df = pd.DataFrame(
        [
            {
                "allotype": allo,
                "formatted": fmt,
                "mhc_class": cls,
                "locus": loc,
                "n_positions": 1,
                "matrix": np.zeros(20, np.float32).tolist(),
                "source": "test",
                "logo": b"\x89PNG\r\n\x1a\n" + fmt.encode(),
            }
            for allo, fmt, cls, loc in rows
        ]
    )
    p = tmp_path / "multi.parquet"
    write_reference(df, p)
    return p


def test_export_logos_writes_png(tmp_path):
    ref = _ref_with_logos(tmp_path, with_logo=True)
    out = tmp_path / "logos"
    n = export_logos(ref, out)
    assert n == 1
    assert (out / "A0201.png").read_bytes().startswith(b"\x89PNG")


def test_export_logos_none_when_no_logo_column(tmp_path):
    ref = _ref_with_logos(tmp_path, with_logo=False)
    assert export_logos(ref, tmp_path / "logos") == 0


def test_export_all_by_default(tmp_path):
    """allotypes=None (the default) must export every logo, not crash."""
    ref = _multi_ref(tmp_path)
    out = tmp_path / "all"
    assert export_logos(ref, out) == 3
    assert {p.name for p in out.glob("*.png")} == {
        "HLAB39124.png",
        "HLAB3942.png",
        "H2Kb.png",
    }


@pytest.mark.parametrize(
    "requested,expected_file",
    [
        ("HLA-B39:124", "HLAB39124.png"),  # raw allotype form
        ("HLAB39124", "HLAB39124.png"),  # formatted key
        ("B39124", "HLAB39124.png"),  # prefix-stripped short form
        ("HLA-B*39:124", "HLAB39124.png"),  # IMGT star form
        ("H-2-Kb", "H2Kb.png"),  # mouse raw
        ("Kb", "H2Kb.png"),  # mouse prefix-stripped
    ],
)
def test_export_filter_forms_all_match(tmp_path, requested, expected_file):
    ref = _multi_ref(tmp_path)
    out = tmp_path / "one"
    assert export_logos(ref, out, requested) == 1
    assert (out / expected_file).exists()


def test_export_filter_list_multiple(tmp_path):
    ref = _multi_ref(tmp_path)
    out = tmp_path / "two"
    n = export_logos(ref, out, ["HLAB39124", "H-2-Kb"])
    assert n == 2
    assert {p.name for p in out.glob("*.png")} == {"HLAB39124.png", "H2Kb.png"}


def test_export_duplicate_request_dedups(tmp_path):
    ref = _multi_ref(tmp_path)
    out = tmp_path / "dup"
    # Two spellings of the same allotype -> written once.
    assert export_logos(ref, out, ["HLAB39124", "B39124"]) == 1


def test_export_unknown_allotype_raises_cleanly(tmp_path):
    """The unknown-allotype error path must build its message without crashing
    (previously sliced a set, raising TypeError instead of the ValueError)."""
    ref = _multi_ref(tmp_path)
    with pytest.raises(ValueError, match="not found"):
        export_logos(ref, tmp_path / "x", "HLA-Z99:99")
