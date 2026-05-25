"""Rich terminal UI: banner, colored logging, and the results table."""

from __future__ import annotations

import pandas as pd
from rich.console import Console
from rich.table import Table
from rich.text import Text

from mhc_tp import __version__

CONSOLE = Console(record=True)  # record=True so logs can be exported to a file

_BANNER = r"""
 __  __ _   _  ____    _____ ____  
|  \/  | | | |/ ___|  |_   _|  _ \ 
| |\/| | |_| | |   _____| | | |_) |
| |  | |  _  | |__|_____| | |  __/ 
|_|  |_|_| |_|\____|    |_| |_|                            
"""


def banner() -> None:
    """Print the MHC-TP banner + credits."""
    CONSOLE.print(_BANNER, style="bold cyan")
    text = Text()
    text.append("MHC-TP ", style="bold")
    text.append(f"v{__version__}", style="cyan")
    text.append("  cluster immunopeptidomics peptides by HLA/MHC motif\n", style="dim")
    text.append("Li Lab / Purcell Lab, Monash University", style="dim")
    CONSOLE.print(text)


def results_table(df: pd.DataFrame, threshold: float = 0.70, top: int = 25) -> None:
    """Print the top cluster→HLA matches as a colored Rich table."""
    table = Table(title="Top cluster → HLA/MHC matches", title_style="bold")
    table.add_column("#", justify="right", style="dim")
    table.add_column("Cluster", style="bold")
    table.add_column("Best HLA/MHC")
    table.add_column("PCC", justify="right")
    if "kld" in df.columns:
        table.add_column("KLD", justify="right")
    for i, (_, r) in enumerate(df.head(top).iterrows(), 1):
        corr = float(r["correlation"])
        style = "green" if corr >= threshold else "yellow" if corr >= 0.5 else "red"
        cells = [
            str(i),
            str(r["cluster"]),
            str(r["hla"]),
            Text(f"{corr:.3f}", style=style),
        ]
        if "kld" in df.columns:
            kld = r["kld"]
            cells.append("NA" if pd.isna(kld) else f"{float(kld):.3f}")
        table.add_row(*cells)
    CONSOLE.print(table)
