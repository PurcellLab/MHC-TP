"""DEV-ONLY: generate reference sequence logos via Seq2Logo for the parquet.

Invokes the user's own (external) Seq2Logo install through the plugin. Never
ships Seq2Logo; used only at DB-construction time.
"""

from __future__ import annotations

import os


def reference_logo_bytes(
    matrix_file: str | os.PathLike,
    seq2logo_path: str | None = None,
    python_exe: str | None = None,
    title: str = "",
) -> bytes | None:
    """Return Seq2Logo PNG bytes for ``matrix_file``, or None if unavailable.

    Returns None (rather than raising) when Seq2Logo is not configured, so the
    build degrades gracefully to a logo-less parquet.
    """
    from mhc_tp.report.seq2logo import Seq2LogoNotConfigured, Seq2LogoRenderer

    try:
        renderer = Seq2LogoRenderer(seq2logo_path=seq2logo_path, python_exe=python_exe)
    except Seq2LogoNotConfigured:
        return None
    return renderer.render_png_bytes(matrix_file, title=title)
