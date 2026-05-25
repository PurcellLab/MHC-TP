"""Display formatting for allele names (Immunolyser report style)."""

from __future__ import annotations

import re

_MOUSE_PREFIX = re.compile(r"^H-?2[-_]?", re.IGNORECASE)
_HUMAN_PREFIX = re.compile(r"^HLA[-_]?", re.IGNORECASE)
_NON_ALNUM = re.compile(r"[^A-Za-z0-9]")


def pretty_allele(allotype: str | None) -> str:
    """Render a source allotype in the Immunolyser display style.

    The species prefix and field separators are normalised while the locus and
    allele-field digits are preserved verbatim (no field-split guessing):

        ``HLA-B39:124`` -> ``HLA_B39124``
        ``HLA-B3942``   -> ``HLA_B3942``
        ``DRB1_0101``   -> ``HLA_DRB10101``
        ``H-2-Kb``      -> ``H2_Kb``
        ``H-2-IAb``     -> ``H2_IAb``

    Already-formatted display names are stable (idempotent).
    """
    n = (allotype or "").strip()
    if not n:
        return ""
    if re.match(r"^H-?2", n, re.IGNORECASE):
        core = _MOUSE_PREFIX.sub("", n)
        return "H2_" + _NON_ALNUM.sub("", core)
    core = _HUMAN_PREFIX.sub("", n)
    return "HLA_" + _NON_ALNUM.sub("", core)
