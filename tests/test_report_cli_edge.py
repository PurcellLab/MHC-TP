"""Edge-case tests for the report, CLI, and runtime layers.

Run: ``pixi run pytest tests/test_report_cli_edge.py -v``

These exercise boundary behaviour of the thread-budget resolver, the report
data-shaping helpers, the image/asset helpers, the HTML renderer, and the CLI
glue. Everything uses tmp_path + synthetic data; nothing reads the packaged
parquet reference data.
"""

from __future__ import annotations

import base64

import numpy as np
import pandas as pd
import pytest

from hla_pepclust import runtime
from hla_pepclust.report import assets, data
from hla_pepclust.report.render import render_report

# ---------------------------------------------------------------------------
# runtime
# ---------------------------------------------------------------------------


def test_resolve_threads_explicit_wins_over_env(monkeypatch):
    monkeypatch.setenv("HLA_PEPCLUST_THREADS", "9")
    # An explicit argument beats the env var.
    assert runtime.resolve_threads(2) == 2


def test_resolve_threads_env_honoured(monkeypatch):
    monkeypatch.setenv("HLA_PEPCLUST_THREADS", "7")
    assert runtime.resolve_threads(None) == 7


def test_resolve_threads_default_when_no_env(monkeypatch):
    monkeypatch.delenv("HLA_PEPCLUST_THREADS", raising=False)
    assert runtime.resolve_threads(None) == 4


def test_resolve_threads_non_numeric_env_falls_back_to_default(monkeypatch):
    # isdigit() guard: a non-numeric env value must not raise; falls back to 4.
    monkeypatch.setenv("HLA_PEPCLUST_THREADS", "lots")
    assert runtime.resolve_threads(None) == 4


def test_resolve_threads_empty_env_falls_back_to_default(monkeypatch):
    monkeypatch.setenv("HLA_PEPCLUST_THREADS", "")
    assert runtime.resolve_threads(None) == 4


def test_resolve_threads_floors_at_one(monkeypatch):
    monkeypatch.delenv("HLA_PEPCLUST_THREADS", raising=False)
    assert runtime.resolve_threads(0) == 1
    assert runtime.resolve_threads(-5) == 1


def test_apply_thread_env_sets_all_vars(monkeypatch):
    monkeypatch.delenv("HLA_PEPCLUST_THREADS", raising=False)
    for var in runtime._THREAD_ENV_VARS:
        monkeypatch.delenv(var, raising=False)
    n = runtime.apply_thread_env()
    assert n == 4
    for var in runtime._THREAD_ENV_VARS:
        assert runtime.os.environ[var] == "4"


def test_apply_thread_env_does_not_override_existing(monkeypatch):
    # setdefault: an operator-set value must survive untouched.
    monkeypatch.delenv("HLA_PEPCLUST_THREADS", raising=False)
    for var in runtime._THREAD_ENV_VARS:
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setenv("OMP_NUM_THREADS", "13")
    runtime.apply_thread_env()
    assert runtime.os.environ["OMP_NUM_THREADS"] == "13"
    # The other (unset) vars still get the default.
    assert runtime.os.environ["MKL_NUM_THREADS"] == "4"


def test_apply_numba_threads_caps_at_compiled_max(monkeypatch):
    pytest.importorskip("numba")
    import numba

    # Ask for a huge number; numba must cap at its compiled maximum.
    requested = numba.config.NUMBA_NUM_THREADS + 100
    returned = runtime.apply_numba_threads(requested)
    assert returned == runtime.resolve_threads(requested)
    assert numba.get_num_threads() <= numba.config.NUMBA_NUM_THREADS


def test_apply_numba_threads_swallows_failure(monkeypatch):
    # If numba is missing/raises, the function must not propagate.
    import builtins

    real_import = builtins.__import__

    def boom(name, *a, **k):
        if name == "numba":
            raise ImportError("no numba")
        return real_import(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", boom)
    assert runtime.apply_numba_threads(3) == 3


# ---------------------------------------------------------------------------
# report.data
# ---------------------------------------------------------------------------


def test_parse_cluster_id_basic():
    assert data.parse_cluster_id("gibbs.1of3.mat") == ("1of3", 1, 3)


def test_parse_cluster_id_other():
    assert data.parse_cluster_id("2of5.mat") == ("2of5", 2, 5)


def test_parse_cluster_id_no_pattern():
    assert data.parse_cluster_id("somefile.mat") == ("somefile.mat", 0, 0)


def test_parse_cluster_id_first_match_wins():
    # The regex search picks the first NofM occurrence.
    assert data.parse_cluster_id("gibbs.10of12.foo") == ("10of12", 10, 12)


def test_pcc_records_keys_and_values():
    cd = {("gibbs.1of3.mat", "A0201"): 0.812345, ("gibbs.2of3.mat", "B0702"): 0.5}
    recs = data.pcc_records(cd)
    assert {r["Cluster"] for r in recs} == {"1of3", "2of3"}
    for r in recs:
        assert set(r.keys()) == {"Cluster", "HLA", "Correlation"}
    # rounded to 4 dp
    rec1 = next(r for r in recs if r["Cluster"] == "1of3")
    assert rec1["HLA"] == "A0201"
    assert rec1["Correlation"] == 0.8123


def test_datatable_rows_sorted_desc_kld_none():
    cd = {
        ("gibbs.1of3.mat", "A0201"): 0.4,
        ("gibbs.2of3.mat", "B0702"): 0.9,
        ("gibbs.3of3.mat", "C0701"): 0.6,
    }
    rows = data.datatable_rows(cd, kld_df=None)
    corrs = [r["correlation"] for r in rows]
    assert corrs == sorted(corrs, reverse=True)
    assert all(r["kld"] is None for r in rows)
    assert set(rows[0].keys()) == {"cluster", "hla", "correlation", "kld"}


def test_datatable_rows_kld_looked_up():
    cd = {("gibbs.1of3.mat", "A0201"): 0.81, ("gibbs.2of3.mat", "B0702"): 0.7}
    kld = pd.DataFrame(
        {"cluster": [3], "group1": [1.5], "group2": [2.5], "total": [4.0]}
    )
    rows = data.datatable_rows(cd, kld_df=kld)
    by_cluster = {r["cluster"]: r["kld"] for r in rows}
    assert by_cluster["1of3"] == 1.5
    assert by_cluster["2of3"] == 2.5


def test_datatable_rows_kld_missing_group_column():
    cd = {("gibbs.1of3.mat", "A0201"): 0.81}
    # kld_df present but only has group2; group1 lookup => None handled gracefully.
    kld = pd.DataFrame({"cluster": [3], "group2": [2.5]})
    rows = data.datatable_rows(cd, kld_df=kld)
    assert rows[0]["kld"] is None


def test_datatable_rows_kld_no_matching_cluster():
    cd = {("gibbs.1of3.mat", "A0201"): 0.81}
    kld = pd.DataFrame({"cluster": [5], "group1": [1.0]})
    rows = data.datatable_rows(cd, kld_df=kld)
    assert rows[0]["kld"] is None


def test_datatable_rows_no_pattern_name_kld_none():
    # nclust==0 -> _lookup_kld short-circuits to None even with a kld_df.
    cd = {("plain_name", "A0201"): 0.5}
    kld = pd.DataFrame({"cluster": [0], "group0": [1.0]})
    rows = data.datatable_rows(cd, kld_df=kld)
    assert rows[0]["kld"] is None
    assert rows[0]["cluster"] == "plain_name"


# ---------------------------------------------------------------------------
# report.assets
# ---------------------------------------------------------------------------

# A 1x1 transparent PNG.
_PNG_1X1 = base64.b64decode(
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNk+M8AAAMBAQDJ/pLv"
    "AAAAAElFTkSuQmCC"
)


def test_png_bytes_to_data_uri_prefix():
    uri = assets.png_bytes_to_data_uri(_PNG_1X1)
    assert uri.startswith("data:image/png;base64,")
    payload = uri.split(",", 1)[1]
    assert base64.b64decode(payload) == _PNG_1X1


def test_png_file_to_data_uri_roundtrip(tmp_path):
    p = tmp_path / "logo.png"
    p.write_bytes(_PNG_1X1)
    uri = assets.png_file_to_data_uri(p)
    assert uri is not None and uri.startswith("data:image/png;base64,")


def test_png_file_to_data_uri_missing_returns_none(tmp_path):
    assert assets.png_file_to_data_uri(tmp_path / "nope.png") is None


def test_find_cluster_logo_finds_png(tmp_path):
    logos = tmp_path / "logos"
    logos.mkdir()
    (logos / "gibbs_logos_1of3-001.png").write_bytes(_PNG_1X1)
    uri = assets.find_cluster_logo(tmp_path, "1of3")
    assert uri is not None and uri.startswith("data:image/png;base64,")


def test_find_cluster_logo_alt_png_name(tmp_path):
    logos = tmp_path / "logos"
    logos.mkdir()
    # the second candidate name (without -001)
    (logos / "gibbs_logos_2of2.png").write_bytes(_PNG_1X1)
    uri = assets.find_cluster_logo(tmp_path, "2of2")
    assert uri is not None and uri.startswith("data:image/png;base64,")


def test_find_cluster_logo_missing_returns_none(tmp_path):
    (tmp_path / "logos").mkdir()
    assert assets.find_cluster_logo(tmp_path, "1of3") is None


def test_find_cluster_logo_no_logos_dir_returns_none(tmp_path):
    # No logos/ subdir at all -> None (png_file_to_data_uri handles missing path).
    assert assets.find_cluster_logo(tmp_path, "1of3") is None


# ---------------------------------------------------------------------------
# report.render
# ---------------------------------------------------------------------------


def _ref_df(formatted="A0201", n_positions=9, with_logo=False):
    cols = {
        "allotype": ["A*02:01"],
        "formatted": [formatted],
        "mhc_class": ["I"],
        "locus": ["A"],
        "n_positions": [n_positions],
        "matrix": [
            np.random.default_rng(0)
            .random(n_positions * 20)
            .astype(np.float32)
            .tolist()
        ],
        "source": ["NetMHCpan-4.2"],
    }
    if with_logo:
        cols["logo"] = [_PNG_1X1]
    return pd.DataFrame(cols)


def test_render_report_uses_logo_map(tmp_path):
    ref = _ref_df(with_logo=False)  # no logo column
    cd = {("gibbs.1of1.mat", "A0201"): 0.9}
    gibbs = {
        "gibbs.1of1.mat": np.random.default_rng(1).random((9, 20)).astype(np.float32)
    }
    logo_map = {"A0201": _PNG_1X1}
    path = render_report(cd, ref, gibbs, str(tmp_path), logo_map=logo_map)
    html = open(path).read()
    # The supplied PNG was embedded as a base64 data URI.
    assert "data:image/png;base64," in html
    assert base64.b64encode(_PNG_1X1).decode("ascii") in html
    assert "const PCC =" in html
    assert path.endswith("clust-search-result.html")


def test_render_report_empty_correlation_dict(tmp_path):
    ref = _ref_df()
    path = render_report({}, ref, {}, str(tmp_path))
    html = open(path).read()
    # Valid HTML produced, no crash.
    assert path.endswith("clust-search-result.html")
    assert "const PCC =" in html
    assert "<html" in html.lower()


def test_render_report_logomaker_fallback_no_logo(tmp_path):
    # reference_df has no 'logo' column and no logo_map -> logomaker fallback.
    ref = _ref_df(with_logo=False)
    cd = {("gibbs.1of1.mat", "A0201"): 0.77}
    gibbs = {
        "gibbs.1of1.mat": np.random.default_rng(2).random((9, 20)).astype(np.float32)
    }
    path = render_report(cd, ref, gibbs, str(tmp_path))  # no logo_map
    html = open(path).read()
    assert "data:image/png;base64," in html  # fallback logo rendered
    assert "A0201" in html


def test_render_report_sections_grouped_by_n_clusters(tmp_path):
    # Two cluster-counts (N=2 and N=3); two refs so each cluster matches one.
    ref = pd.DataFrame(
        {
            "allotype": ["A*02:01", "B*07:02"],
            "formatted": ["A0201", "B0702"],
            "mhc_class": ["I", "I"],
            "locus": ["A", "B"],
            "n_positions": [9, 9],
            "matrix": [
                np.random.default_rng(3).random(9 * 20).astype(np.float32).tolist(),
                np.random.default_rng(4).random(9 * 20).astype(np.float32).tolist(),
            ],
            "source": ["NetMHCpan-4.2", "NetMHCpan-4.2"],
        }
    )
    cd = {
        ("gibbs.1of2.mat", "A0201"): 0.95,
        ("gibbs.2of2.mat", "B0702"): 0.90,
        ("gibbs.1of3.mat", "A0201"): 0.85,
    }
    gibbs = {
        "gibbs.1of2.mat": np.random.default_rng(5).random((9, 20)).astype(np.float32),
        "gibbs.2of2.mat": np.random.default_rng(6).random((9, 20)).astype(np.float32),
        "gibbs.1of3.mat": np.random.default_rng(7).random((9, 20)).astype(np.float32),
    }
    path = render_report(cd, ref, gibbs, str(tmp_path))
    html = open(path).read()
    # Both N=2 and N=3 cluster-count groupings appear.
    assert "2 clusters" in html
    assert "3 clusters" in html


def test_render_report_kld_shown(tmp_path):
    ref = _ref_df()
    cd = {("gibbs.1of3.mat", "A0201"): 0.81}
    gibbs = {
        "gibbs.1of3.mat": np.random.default_rng(8).random((9, 20)).astype(np.float32)
    }
    kld = pd.DataFrame(
        {"cluster": [3], "group1": [1.5], "group2": [0.0], "group3": [0.0]}
    )
    path = render_report(cd, ref, gibbs, str(tmp_path), kld_df=kld)
    html = open(path).read()
    # Template renders the table KLD column with %.3f and the section with %.2f.
    assert "1.500" in html  # KLD surfaced in the results table
    assert "KLD 1.50" in html  # KLD surfaced in the per-cluster section


def test_render_report_skips_unmatched_hla(tmp_path):
    # correlation references an HLA not in the reference_df -> skipped, no crash.
    ref = _ref_df(formatted="A0201")
    cd = {("gibbs.1of1.mat", "GHOST"): 0.99}
    gibbs = {
        "gibbs.1of1.mat": np.random.default_rng(9).random((9, 20)).astype(np.float32)
    }
    path = render_report(cd, ref, gibbs, str(tmp_path))
    assert path.endswith("clust-search-result.html")


# ---------------------------------------------------------------------------
# cli
# ---------------------------------------------------------------------------


def test_main_version_exits_zero():
    from hla_pepclust import cli

    with pytest.raises(SystemExit) as exc:
        cli.main(["--version"])
    assert exc.value.code == 0


def test_main_search_missing_reference_raises():
    from hla_pepclust import cli

    # resolve_reference raises FileNotFoundError for a non-existent override.
    with pytest.raises(FileNotFoundError):
        cli.main(["search", "nonexistent_gibbs", "-r", "nope.parquet"])


def _make_mat_file(path, arr):
    """Write a PSSM .mat file parseable by parse_matrix."""
    header = "Pos A R N D C Q E G H I L K M F P S T W Y V"
    lines = [header]
    for i, row in enumerate(arr, start=1):
        scores = " ".join(f"{v:.4f}" for v in row)
        lines.append(f"{i} P {scores}")
    path.write_text("\n".join(lines) + "\n")


def test_load_gibbs_matrices_filters_mat_and_n_clusters(tmp_path):
    from hla_pepclust import cli

    matrices = tmp_path / "matrices"
    matrices.mkdir()
    arr = np.random.default_rng(0).random((9, 20))
    _make_mat_file(matrices / "gibbs.1of3.mat", arr)
    _make_mat_file(matrices / "gibbs.2of3.mat", arr)
    _make_mat_file(matrices / "gibbs.1of2.mat", arr)
    # non-.mat files must be ignored
    (matrices / "notes.txt").write_text("ignore me")
    (matrices / "gibbs.1of3.png").write_bytes(_PNG_1X1)

    # n_clusters="all" -> all three .mat files
    all_m = cli._load_gibbs_matrices(str(tmp_path), n_clusters="all")
    assert set(all_m) == {"gibbs.1of3.mat", "gibbs.2of3.mat", "gibbs.1of2.mat"}

    # n_clusters="3" -> only the *of3.mat files
    only3 = cli._load_gibbs_matrices(str(tmp_path), n_clusters="3")
    assert set(only3) == {"gibbs.1of3.mat", "gibbs.2of3.mat"}


def test_run_search_no_html_writes_correlations_csv(tmp_path):
    from hla_pepclust import cli
    from hla_pepclust.refdata.parquet_io import write_reference

    # Build a tiny synthetic reference parquet.
    n = 9
    mat = np.random.default_rng(0).random(n * 20).astype(np.float32)
    ref = pd.DataFrame(
        {
            "allotype": ["A*02:01"],
            "formatted": ["A0201"],
            "mhc_class": ["I"],
            "locus": ["A"],
            "n_positions": [n],
            "matrix": [mat.tolist()],
            "source": ["NetMHCpan-4.2"],
        }
    )
    ref_path = tmp_path / "ref.parquet"
    write_reference(ref, ref_path)

    # Gibbs matrix that is an exact copy of the reference -> correlation ~1.0.
    gibbs_dir = tmp_path / "gibbs"
    matrices = gibbs_dir / "matrices"
    matrices.mkdir(parents=True)
    _make_mat_file(matrices / "gibbs.1of1.mat", mat.reshape(n, 20))

    out = tmp_path / "out"
    df = cli.run_search(
        str(gibbs_dir),
        str(ref_path),
        "human",
        str(out),
        threshold=0.1,
        make_html=False,
    )
    csv_path = out / "clust_result" / "correlations.csv"
    assert csv_path.exists()
    written = pd.read_csv(csv_path)
    assert list(written.columns) == ["cluster", "hla", "formatted", "correlation"]
    assert (written["hla"] == "HLA_A0201").any()
    assert (written["formatted"] == "A0201").any()
    # No HTML produced when make_html=False.
    assert not (out / "clust_result" / "clust-search-result.html").exists()
    assert isinstance(df, pd.DataFrame)


def test_run_search_no_matches_writes_empty_csv(tmp_path):
    from hla_pepclust import cli
    from hla_pepclust.refdata.parquet_io import write_reference

    n = 9
    rng = np.random.default_rng(11)
    ref = pd.DataFrame(
        {
            "allotype": ["A*02:01"],
            "formatted": ["A0201"],
            "mhc_class": ["I"],
            "locus": ["A"],
            "n_positions": [n],
            "matrix": [rng.random(n * 20).astype(np.float32).tolist()],
            "source": ["NetMHCpan-4.2"],
        }
    )
    ref_path = tmp_path / "ref.parquet"
    write_reference(ref, ref_path)

    gibbs_dir = tmp_path / "gibbs"
    matrices = gibbs_dir / "matrices"
    matrices.mkdir(parents=True)
    # A different random matrix; threshold ~1.0 ensures no hits.
    _make_mat_file(matrices / "gibbs.1of1.mat", rng.random((n, 20)))

    out = tmp_path / "out"
    df = cli.run_search(
        str(gibbs_dir),
        str(ref_path),
        "human",
        str(out),
        threshold=0.999999,
        make_html=False,
    )
    csv_path = out / "clust_result" / "correlations.csv"
    assert csv_path.exists()
    written = pd.read_csv(csv_path)
    assert list(written.columns) == ["cluster", "hla", "formatted", "correlation"]
    assert len(written) == 0
    assert len(df) == 0
