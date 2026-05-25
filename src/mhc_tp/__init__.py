"""MHC-TP: cluster immunopeptidomics peptides by HLA binding motif."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("mhc-tp")
except PackageNotFoundError:  # not installed (e.g. running from a source checkout)
    __version__ = "0.0.0+dev"

__all__ = ["__version__"]
