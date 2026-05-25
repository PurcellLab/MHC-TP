"""Reference Parquet schema (one row per allotype)."""

COLUMNS = [
    "allotype",    # str display name, e.g. A*02:01
    "formatted",   # str match key, e.g. A0201
    "mhc_class",   # "I" | "II"
    "locus",       # A/B/C | DR/DQ/DP | K/D/L | IA/IE
    "n_positions", # int motif length
    "matrix",      # list<float32> flattened n_positions x 20, canonical AA order
    "source",      # provenance, e.g. NetMHCpan-4.2 / NetMHCIIpan-4.3
]
