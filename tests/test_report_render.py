import numpy as np
import pandas as pd

from mhc_tp.report.render import render_report


def _ref_df():
    return pd.DataFrame(
        {
            "allotype": ["H2-Kb"],
            "formatted": ["H2Kb"],
            "mhc_class": ["I"],
            "locus": ["K"],
            "n_positions": [9],
            "matrix": [
                np.random.default_rng(1).random(9 * 20).astype(np.float32).tolist()
            ],
            "source": ["NetMHCpan-4.2"],
        }
    )


def test_render_report_writes_standalone_html(tmp_path):
    ref = _ref_df()
    cd = {("gibbs.1of1.mat", "H2Kb"): 0.81}
    gibbs = {
        "gibbs.1of1.mat": np.random.default_rng(2).random((9, 20)).astype(np.float32)
    }
    out = tmp_path / "out"
    path = render_report(cd, ref, gibbs, str(out), kld_df=None, version="2.0.0-dev")
    html = open(path).read()
    assert "correlation_table" in html
    assert "data:image/png;base64," in html  # logos inlined
    assert "H2Kb" in html
    assert "const PCC =" in html  # pcc data inlined
    assert "clusters" in html  # per-cluster-count section heading
    assert path.endswith("mhc-tp-result.html")
