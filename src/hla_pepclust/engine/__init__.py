"""Numba search engine.

Import the public entry point lazily to avoid paying the numba import cost
(~1.4s) unless the engine is actually used:

    from hla_pepclust.engine.search import search
"""
