"""Render sequence logos from PSSM matrices via logomaker (base64 PNG)."""

from __future__ import annotations

import base64
import io

import matplotlib

matplotlib.use("Agg")  # headless

import logomaker  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from hla_pepclust.constants import AMINO_ACIDS


def render_logo(matrix: np.ndarray, title: str | None = None) -> str:
    """Render an information-content logo for an (n_positions, 20) PSSM.

    Returns a ``data:image/png;base64,...`` URI suitable for inlining in HTML.
    """
    arr = np.asarray(matrix, dtype=float)
    df = pd.DataFrame(arr, columns=AMINO_ACIDS)
    prob = df.clip(lower=0)
    row_sums = prob.sum(axis=1).replace(0, np.nan)
    prob = prob.div(row_sums, axis=0).fillna(0.0)
    try:
        info = logomaker.transform_matrix(
            prob, from_type="probability", to_type="information"
        )
    except Exception:
        info = prob  # fallback: plot the probability matrix directly

    fig, ax = plt.subplots(figsize=(max(2.0, 0.5 * len(info)), 2.0))
    logomaker.Logo(info, ax=ax, color_scheme="skylign_protein")
    ax.set_xticks(range(len(info)))
    ax.set_xticklabels([str(i + 1) for i in range(len(info))])
    if title:
        ax.set_title(title, fontsize=9)
    ax.set_ylabel("bits")
    fig.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=100)
    plt.close(fig)
    encoded = base64.b64encode(buf.getvalue()).decode("ascii")
    return f"data:image/png;base64,{encoded}"
