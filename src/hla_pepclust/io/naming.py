"""Normalize HLA/MHC allotype names for class I and II, human and mouse.

Ported from the formatting logic in ``cli/cluster_search.py``
(``formate_HLA_DB`` / ``formate_HLA_user_in``) and the
``formatted_allotypes`` column produced for the reference databases
(``data/ref_data/{human,mouse}.db``).

Reference examples reproduced here:

* human  ``HLA-C*07:485`` -> ``C07485``
* mouse  ``MHC*H2_Db``    -> ``H2Db`` (the ``H2`` prefix is retained)
"""

from __future__ import annotations

import re
from dataclasses import dataclass

_CLASS_II_HUMAN = ("DR", "DQ", "DP")
_CLASS_II_MOUSE = ("IA", "IE")


@dataclass(frozen=True)
class AllotypeInfo:
    raw: str
    formatted: str
    mhc_class: str
    locus: str


def format_allotype(name: str, species: str) -> AllotypeInfo:
    """Normalize an allotype name into an :class:`AllotypeInfo`.

    ``formatted`` is the match key used against the reference database
    (no ``*``/``:``/spaces/dots/underscores).
    """
    raw = name
    species = species.lower()

    # Strip the species/database prefixes the legacy formatters removed.
    core = name.replace("MHC*", "").replace("MHC_", "")
    core = (
        core.replace("HLA-", "")
        .replace("HLA_", "")
        .replace("H-2", "H2")
        .replace("H2-", "H2")
        .replace("H2_", "H2")
    )

    if species == "human":
        locus = next(
            (loc for loc in _CLASS_II_HUMAN if core.startswith(loc)),
            core[:1],
        )
        mhc_class = "II" if locus in _CLASS_II_HUMAN else "I"
    else:
        # Mouse names keep the ``H2`` prefix; the locus follows it.
        rest = core[2:] if core.startswith("H2") else core
        locus = next(
            (loc for loc in _CLASS_II_MOUSE if rest.startswith(loc)),
            rest[:1],
        )
        mhc_class = "II" if locus in _CLASS_II_MOUSE else "I"

    # Match key: drop separators while keeping the alphanumeric body.
    formatted = re.sub(r"[*:\s._-]", "", core)

    return AllotypeInfo(
        raw=raw,
        formatted=formatted,
        mhc_class=mhc_class,
        locus=locus,
    )
