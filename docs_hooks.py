"""MkDocs build hooks.

Emit a raw plain-text copy of the LLM guide at the site root (``/llm.txt``) so
machine/LLM consumers can fetch the instructions without HTML. The source of
truth is ``docs/llm.md``; this copies it verbatim on every build.
"""

from __future__ import annotations

import pathlib


def on_post_build(config, **kwargs) -> None:
    site = pathlib.Path(config["site_dir"])
    src = pathlib.Path(config["docs_dir"]) / "llm.md"
    if src.exists():
        (site / "llm.txt").write_text(src.read_text(encoding="utf-8"), encoding="utf-8")
