import numpy as np
import pandas as pd
from mhc_tp.refdata.parquet_io import write_reference, read_reference


def _sample_df():
    return pd.DataFrame(
        {
            "allotype": ["A*02:01", "A*01:01"],
            "formatted": ["A0201", "A0101"],
            "mhc_class": ["I", "I"],
            "locus": ["A", "A"],
            "n_positions": [2, 2],
            "matrix": [
                np.zeros(40, np.float32).tolist(),
                np.ones(40, np.float32).tolist(),
            ],
            "source": ["NetMHCpan-4.2", "NetMHCpan-4.2"],
        }
    )


def test_roundtrip(tmp_path):
    df = _sample_df()
    out = tmp_path / "human.parquet"
    write_reference(df, out)
    back = read_reference(out)
    assert list(back["formatted"]) == ["A0201", "A0101"]
    assert back["matrix"].iloc[1][0] == 1.0
    assert out.exists()


def test_write_rejects_missing_columns(tmp_path):
    import pytest

    bad = pd.DataFrame({"allotype": ["X"]})
    with pytest.raises(ValueError):
        write_reference(bad, tmp_path / "bad.parquet")
