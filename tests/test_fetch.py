import pytest

from hla_pepclust.refdata.fetch import (
    data_dir,
    load_manifest,
    reference_path,
    resolve_reference,
)


def test_data_dir_env_override(tmp_path, monkeypatch):
    monkeypatch.setenv("HLA_PEPCLUST_DATA_DIR", str(tmp_path))
    assert data_dir() == tmp_path
    assert reference_path("Human") == tmp_path / "human.parquet"


def test_resolve_override(tmp_path):
    p = tmp_path / "r.parquet"
    p.write_bytes(b"x")
    assert resolve_reference("human", str(p)) == p


def test_resolve_missing_raises_with_guidance(tmp_path, monkeypatch):
    monkeypatch.setenv("HLA_PEPCLUST_DATA_DIR", str(tmp_path))
    with pytest.raises(FileNotFoundError, match="clust-search fetch"):
        resolve_reference("human")


def test_manifest_loads():
    rows = load_manifest()
    species = {r["species"] for r in rows}
    assert {"human", "mouse"} <= species
    assert all(r["url"].startswith("https://") for r in rows)
    assert all(len(r["sha256"]) == 64 for r in rows)
