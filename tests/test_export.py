import numpy as np
import pandas as pd

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


def test_export_logos_writes_png(tmp_path):
    ref = _ref_with_logos(tmp_path, with_logo=True)
    out = tmp_path / "logos"
    n = export_logos(ref, out)
    assert n == 1
    assert (out / "A0201.png").read_bytes().startswith(b"\x89PNG")


def test_export_logos_none_when_no_logo_column(tmp_path):
    ref = _ref_with_logos(tmp_path, with_logo=False)
    assert export_logos(ref, tmp_path / "logos") == 0
