import base64

from hla_pepclust.report.assets import (
    find_cluster_logo,
    png_bytes_to_data_uri,
    png_file_to_data_uri,
)


def test_png_bytes_to_data_uri():
    uri = png_bytes_to_data_uri(b"hello")
    assert uri == "data:image/png;base64," + base64.b64encode(b"hello").decode()


def test_png_file_to_data_uri(tmp_path):
    p = tmp_path / "x.png"
    p.write_bytes(b"\x89PNG\r\n")
    assert png_file_to_data_uri(p).startswith("data:image/png;base64,")
    assert png_file_to_data_uri(tmp_path / "missing.png") is None


def test_find_cluster_logo_png(tmp_path):
    logos = tmp_path / "logos"
    logos.mkdir()
    (logos / "gibbs_logos_1of3-001.png").write_bytes(b"\x89PNG\r\nfake")
    uri = find_cluster_logo(tmp_path, "1of3")
    assert uri and uri.startswith("data:image/png;base64,")


def test_find_cluster_logo_missing(tmp_path):
    (tmp_path / "logos").mkdir()
    assert find_cluster_logo(tmp_path, "9of9") is None
