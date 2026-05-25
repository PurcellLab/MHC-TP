import numpy as np
from hla_pepclust.report.logos import render_logo


def test_render_logo_returns_data_uri():
    rng = np.random.default_rng(0)
    m = rng.random((9, 20)).astype(np.float32)
    uri = render_logo(m)
    assert isinstance(uri, str)
    assert uri.startswith("data:image/png;base64,")
    assert len(uri) > 100


def test_render_logo_handles_all_zero_row():
    m = np.zeros((9, 20), dtype=np.float32)
    m[0, 0] = 1.0
    uri = render_logo(m)  # must not raise
    assert uri.startswith("data:image/png;base64,")
