"""Display formatting for allele names.

The reference parquet stores each allotype in the very string that was burned
into its Seq2Logo motif title (e.g. ``HLA-A25:08``, ``HLA-A3301``,
``DRB1_0101``, ``H-2-IAb``). To keep the report's text labels consistent with
those embedded reference logo titles, the canonical display name is the source
allotype verbatim.
"""

from __future__ import annotations


def pretty_allele(allotype: str | None) -> str:
    """Canonical display name for an allotype, matching its reference logo title.

    The source allotype is already the name shown on the embedded Seq2Logo
    motif (``HLA-A25:08``, ``HLA-A3301``, ``DRB1_0101`` ...), so it is returned
    verbatim apart from surrounding whitespace.

        ``HLA-A25:08`` -> ``HLA-A25:08``
        ``HLA-A3301``  -> ``HLA-A3301``
        ``DRB1_0101``  -> ``DRB1_0101``
        ``H-2-IAb``    -> ``H-2-IAb``
    """
    return (allotype or "").strip()
