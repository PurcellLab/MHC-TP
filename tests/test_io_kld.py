import textwrap
import pytest
from hla_pepclust.io.kld import read_kld


def test_read_kld(tmp_path):
    f = tmp_path / "gibbs.KLDvsClusters.tab"
    f.write_text(textwrap.dedent(
        """\
        Number\tg1\tg2
        1\t2.0\t0
        2\t1.5\t1.0
        """
    ))
    df = read_kld(f)
    assert list(df["cluster"]) == [1, 2]
    assert df.loc[df["cluster"] == 2, "total"].iloc[0] == 2.5
    assert "group1" in df.columns and "group2" in df.columns


def test_read_kld_missing_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_kld(tmp_path / "nope.tab")
