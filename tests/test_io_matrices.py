import numpy as np
from mhc_tp.io.matrices import parse_matrix


def test_parse_matrix_shape_and_values(pssm_text):
    m = parse_matrix(pssm_text)
    assert m.shape == (2, 20)
    assert m.dtype == np.float32
    assert m[0, 0] == np.float32(0.1)
    assert m[0, 19] == np.float32(0.7)
    assert m[1, 9] == np.float32(0.5)
    assert m[1, 10] == np.float32(0.5)


def test_parse_matrix_missing_file_returns_none(tmp_path):
    assert parse_matrix(tmp_path / "nope.txt") is None
