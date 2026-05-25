import numpy as np
import pandas as pd
from hla_pepclust.cli import run_search
from hla_pepclust.refdata.parquet_io import write_reference


def test_run_search_writes_csv(tmp_path):
    ramp = np.arange(20, dtype=np.float32)
    ref = pd.DataFrame(
        {
            "allotype": ["A*02:01"],
            "formatted": ["A0201"],
            "mhc_class": ["I"],
            "locus": ["A"],
            "n_positions": [1],
            "matrix": [ramp.tolist()],
            "source": ["test"],
        }
    )
    ref_path = tmp_path / "human.parquet"
    write_reference(ref, ref_path)

    mdir = tmp_path / "gibbs" / "matrices"
    mdir.mkdir(parents=True)
    (mdir / "1of1.mat").write_text(
        "Pos A R N D C Q E G H I L K M F P S T W Y V\n"
        "1 P " + " ".join(str(float(x)) for x in ramp) + "\n"
    )

    out = tmp_path / "out"
    run_search(
        str(tmp_path / "gibbs"),
        str(ref_path),
        "human",
        str(out),
        threshold=0.5,
        top_n=3,
        make_html=False,
    )
    res = pd.read_csv(out / "clust_result" / "correlations.csv")
    assert {"cluster", "hla", "correlation"} <= set(res.columns)
    assert (res["hla"] == "A0201").any()


def test_run_search_writes_html(tmp_path):
    ramp = np.arange(20, dtype=np.float32)
    ref = pd.DataFrame(
        {
            "allotype": ["A*02:01"],
            "formatted": ["A0201"],
            "mhc_class": ["I"],
            "locus": ["A"],
            "n_positions": [1],
            "matrix": [ramp.tolist()],
            "source": ["test"],
        }
    )
    write_reference(ref, tmp_path / "human.parquet")
    mdir = tmp_path / "gibbs" / "matrices"
    mdir.mkdir(parents=True)
    (mdir / "gibbs.1of1.mat").write_text(
        "Pos A R N D C Q E G H I L K M F P S T W Y V\n"
        "1 P " + " ".join(str(float(x)) for x in ramp) + "\n"
    )
    out = tmp_path / "out"
    run_search(
        str(tmp_path / "gibbs"),
        str(tmp_path / "human.parquet"),
        "human",
        str(out),
        threshold=0.5,
        top_n=3,
    )
    assert (out / "clust_result" / "clust-search-result.html").exists()
