"""End-to-end CLI test on committed fixtures.

Exercises the full pipeline (parse real-format matrices -> numba search ->
CSV + HTML report) using only stdlib path handling, so it runs identically on
Linux, macOS, Windows, and WSL2 (covered by the CI OS matrix).
"""

from pathlib import Path

import pandas as pd

from mhc_tp.cli import run_search
from mhc_tp.io.matrices import parse_matrix
from mhc_tp.refdata.parquet_io import write_reference

DATA = Path(__file__).parent / "data"


def _build_reference(tmp_path):
    rows = []
    for name in ("A0201", "B0702"):
        m = parse_matrix(DATA / "ref_matrices" / f"{name}.txt")
        assert m is not None
        rows.append(
            {
                "allotype": f"HLA-{name}",
                "formatted": name,
                "mhc_class": "I",
                "locus": name[0],
                "n_positions": int(m.shape[0]),
                "matrix": m.reshape(-1).tolist(),
                "source": "test",
            }
        )
    ref = tmp_path / "ref.parquet"
    write_reference(pd.DataFrame(rows), ref)
    return ref


def test_end_to_end_search(tmp_path):
    ref = _build_reference(tmp_path)
    out = tmp_path / "out"
    run_search(
        str(DATA / "gibbs_sample"),
        str(ref),
        "human",
        str(out),
        threshold=0.3,
        top_n=2,
    )
    csv = out / "clust_result" / "correlations.csv"
    html = out / "clust_result" / "mhc-tp-result.html"
    assert csv.exists() and html.exists()

    res = pd.read_csv(csv)
    assert len(res) >= 1
    # gibbs.1of1 was built from the A0201 motif + small noise -> A0201 should win.
    top = res.sort_values("correlation", ascending=False).iloc[0]
    assert top["hla"] == "HLA-A0201"  # canonical display = raw allotype
    assert top["formatted"] == "A0201"  # raw key preserved
    assert top["correlation"] > 0.8
