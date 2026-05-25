import numpy as np
import pandas as pd
from hla_pepclust.engine.search import search


def _ref_df():
    ramp = np.arange(20, dtype=np.float32)  # one position, 20 AA values
    inv = ramp[::-1].copy()
    return pd.DataFrame(
        {
            "allotype": ["A*02:01", "A*01:01"],
            "formatted": ["A0201", "A0101"],
            "mhc_class": ["I", "I"],
            "locus": ["A", "A"],
            "n_positions": [1, 1],
            "matrix": [ramp.tolist(), inv.tolist()],
            "source": ["test", "test"],
        }
    )


def test_search_finds_best():
    ref = _ref_df()
    gibbs = {"1of1.mat": np.arange(20, dtype=np.float32).reshape(1, 20)}  # == A0201
    cd = search(ref, gibbs, threshold=0.5, top_n=2)
    best = max(cd.items(), key=lambda kv: kv[1])
    assert best[0][0] == "1of1.mat"
    assert best[0][1] == "A0201"
    assert abs(best[1] - 1.0) < 1e-5
