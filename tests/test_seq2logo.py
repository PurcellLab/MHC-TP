import pytest

from mhc_tp.report.seq2logo import Seq2LogoNotConfigured, Seq2LogoRenderer


def test_not_configured(monkeypatch):
    monkeypatch.delenv("SEQ2LOGO_PATH", raising=False)
    with pytest.raises(Seq2LogoNotConfigured):
        Seq2LogoRenderer()


def test_missing_script(tmp_path, monkeypatch):
    monkeypatch.delenv("SEQ2LOGO_PATH", raising=False)
    with pytest.raises(Seq2LogoNotConfigured):
        Seq2LogoRenderer(seq2logo_path=str(tmp_path))  # no Seq2Logo.py inside


def test_build_command(tmp_path):
    (tmp_path / "Seq2Logo.py").write_text("# stub")
    r = Seq2LogoRenderer(seq2logo_path=str(tmp_path), python_exe="python2")
    cmd = r.build_command("in.mat", tmp_path / "out", title="HLA-A*02:01")
    assert cmd[0] == "python2"
    assert cmd[1].endswith("Seq2Logo.py")
    assert "-f" in cmd and "in.mat" in cmd
    assert "-I" in cmd and "2" in cmd  # KL logo
    assert "--format" in cmd and "PNG" in cmd
    assert "-t" in cmd and "HLA-A*02:01" in cmd


def test_env_var_path(tmp_path, monkeypatch):
    (tmp_path / "Seq2Logo.py").write_text("# stub")
    monkeypatch.setenv("SEQ2LOGO_PATH", str(tmp_path))
    r = Seq2LogoRenderer()
    assert r.script.name == "Seq2Logo.py"
