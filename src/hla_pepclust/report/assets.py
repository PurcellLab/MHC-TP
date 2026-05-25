"""Image helpers for the report: base64 data URIs + locating gibbs cluster logos."""

from __future__ import annotations

import base64
import subprocess
import tempfile
from pathlib import Path


def png_bytes_to_data_uri(data: bytes) -> str:
    """Encode PNG bytes as a ``data:image/png;base64,...`` URI."""
    return "data:image/png;base64," + base64.b64encode(data).decode("ascii")


def png_file_to_data_uri(path: str | Path) -> str | None:
    """Read a PNG file and return a data URI, or None if absent."""
    p = Path(path)
    if not p.exists():
        return None
    return png_bytes_to_data_uri(p.read_bytes())


def _eps_to_png_bytes(eps_path: Path) -> bytes | None:
    """Convert an EPS to PNG bytes via the system ghostscript (`gs`)."""
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "logo.png"
        cmd = [
            "gs",
            "-q",
            "-dSAFER",
            "-dBATCH",
            "-dNOPAUSE",
            "-sDEVICE=pngalpha",
            "-r150",
            "-dEPSCrop",
            f"-sOutputFile={out}",
            str(eps_path),
        ]
        try:
            r = subprocess.run(cmd, capture_output=True, timeout=60)
        except (OSError, subprocess.TimeoutExpired):
            return None
        if r.returncode != 0 or not out.exists():
            return None
        return out.read_bytes()


def find_cluster_logo(gibbs_dir: str | Path, cluster_id: str) -> str | None:
    """Return a data URI for the GibbsCluster logo of ``cluster_id`` (e.g. ``1of3``).

    Prefers the PNG GibbsCluster emits; falls back to converting its EPS via gs.
    """
    logos = Path(gibbs_dir) / "logos"
    for name in (f"gibbs_logos_{cluster_id}-001.png", f"gibbs_logos_{cluster_id}.png"):
        uri = png_file_to_data_uri(logos / name)
        if uri:
            return uri
    for name in (f"gibbs_logos_{cluster_id}-001.eps", f"gibbs_logos_{cluster_id}.eps"):
        eps = logos / name
        if eps.exists():
            data = _eps_to_png_bytes(eps)
            if data:
                return png_bytes_to_data_uri(data)
    return None
