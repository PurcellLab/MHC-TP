"""Bound CPU/thread usage so MHC-TP is a good citizen on shared servers.

Run as a backend (e.g. the Immunolyser server), MHC-TP must not grab every
core (numba ``parallel=True`` + BLAS default to all cores) and starve
co-located programs. The thread budget resolves from, in order: an explicit
argument, the ``HLA_PEPCLUST_THREADS`` env var, else a conservative default.

``apply_thread_env()`` is called at import time of the CLI *before* numpy/numba
are imported, so the BLAS/OpenMP caps take effect.
"""

from __future__ import annotations

import os

_DEFAULT_THREADS = 4
_THREAD_ENV_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "NUMBA_NUM_THREADS",
)


def resolve_threads(threads: int | None = None) -> int:
    """Thread budget: explicit arg > ``HLA_PEPCLUST_THREADS`` env > default (4)."""
    if threads is None:
        env = os.environ.get("HLA_PEPCLUST_THREADS")
        threads = int(env) if env and env.isdigit() else _DEFAULT_THREADS
    return max(1, int(threads))


def apply_thread_env() -> int:
    """Cap BLAS/OpenMP/numba thread pools via env vars. Call before importing numpy.

    Only sets a var if not already set (respects an explicit operator override).
    Returns the resolved thread budget.
    """
    n = resolve_threads()
    for var in _THREAD_ENV_VARS:
        os.environ.setdefault(var, str(n))
    return n


def apply_numba_threads(threads: int | None = None) -> int:
    """Cap numba's active thread count at runtime (≤ its compiled maximum)."""
    n = resolve_threads(threads)
    try:
        import numba

        numba.set_num_threads(min(n, numba.config.NUMBA_NUM_THREADS))
    except Exception:
        pass
    return n
