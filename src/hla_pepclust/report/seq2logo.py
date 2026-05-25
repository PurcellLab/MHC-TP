"""Seq2Logo 2.1 plugin — render Kullback-Leibler sequence logos from PSSMs.

Wraps the external Seq2Logo 2.1 tool (python 2.7). The tool location is taken
from the ``SEQ2LOGO_PATH`` environment variable (or passed explicitly) — never
hardcoded — so any user can point it at their own install. Run inside the
project's python2.7 env, e.g. ``pixi run -e seq2logo``.

Modelled on the project's original Seq2LogoRunner. Used at DEV time to
pre-generate reference logos that get embedded in the reference Parquet; not a
runtime dependency of the end-user package.
"""

from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path

# Kullback-Leibler logo, bits, beta=50, no clustering, PNG (matches the
# NetMHCpan / reference-pack logo style).
_KL_PARAMS = {
    "-I": "2",
    "-u": "Bits",
    "-b": "50",
    "-C": "0",
    "-l": "1",
    "-S": "1",
    "-i": "1",
    "-p": "5333x4000",
    "-s": "40",
    "--format": "PNG",
}


class Seq2LogoNotConfigured(RuntimeError):
    """Raised when the Seq2Logo tool path is not configured or not found."""


class Seq2LogoRenderer:
    """Render sequence-logo PNGs from PSSM matrix files via Seq2Logo 2.1.

    Args:
        seq2logo_path: directory containing ``Seq2Logo.py``. Defaults to the
            ``SEQ2LOGO_PATH`` environment variable. No path is hardcoded.
        python_exe: python 2.7 interpreter to run the tool. Defaults to the
            ``SEQ2LOGO_PYTHON`` env var, else ``"python"`` (assumes the call is
            made inside the python2.7 env, e.g. ``pixi run -e seq2logo``).
    """

    def __init__(self, seq2logo_path: str | None = None, python_exe: str | None = None):
        path = seq2logo_path or os.environ.get("SEQ2LOGO_PATH")
        if not path:
            raise Seq2LogoNotConfigured(
                "Seq2Logo path not set. Pass seq2logo_path=... or set SEQ2LOGO_PATH "
                "to your seq2logo-2.1 directory."
            )
        self.seq2logo_path = Path(path)
        self.script = self.seq2logo_path / "Seq2Logo.py"
        if not self.script.exists():
            raise Seq2LogoNotConfigured(f"Seq2Logo.py not found under {self.seq2logo_path}")
        self.python_exe = python_exe or os.environ.get("SEQ2LOGO_PYTHON") or "python"

    def build_command(self, matrix_file: str | os.PathLike, output_path: str | os.PathLike,
                      title: str = "") -> list[str]:
        """Construct the Seq2Logo command (KL logo, PNG). Pure — no side effects."""
        params = dict(_KL_PARAMS)
        params["-f"] = str(matrix_file)
        params["-o"] = str(output_path)
        if title:
            params["-t"] = title
        cmd = [self.python_exe, str(self.script)]
        for key, value in params.items():
            cmd.extend([key, str(value)])
        return cmd

    def render(self, matrix_file: str | os.PathLike, out_dir: str | os.PathLike,
               name: str | None = None, title: str = "", timeout: int = 90) -> Path | None:
        """Render a logo PNG; return its path (Seq2Logo appends ``-001``) or None."""
        matrix_file = Path(matrix_file)
        if not matrix_file.exists():
            return None
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        name = name or matrix_file.stem
        output_path = out_dir / name
        cmd = self.build_command(matrix_file, output_path, title)

        # Seq2Logo imports Seq2Logo_module relative to its own directory.
        result = subprocess.run(
            cmd, cwd=str(self.seq2logo_path),
            capture_output=True, text=True, timeout=timeout,
        )
        if result.returncode != 0:
            return None
        for cand in (out_dir / f"{name}-001.png", out_dir / f"{name}.png"):
            if cand.exists():
                return cand
        return None

    def render_png_bytes(self, matrix_file: str | os.PathLike, title: str = "") -> bytes | None:
        """Render and return the PNG bytes (for embedding), or None on failure."""
        with tempfile.TemporaryDirectory() as tmp:
            png = self.render(matrix_file, tmp, name="logo", title=title)
            if png is None:
                return None
            return Path(png).read_bytes()
