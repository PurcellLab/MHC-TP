"""Standard error messages for the CLI (modernized from cli/error_logs.py)."""

from __future__ import annotations


def file_not_found(path: str) -> str:
    return f"File not found: {path}"


def png_not_found(path: str) -> str:
    return f"PNG file not found: {path}"


def reference_not_found(path: str) -> str:
    return f"HLA/MHC reference not found: {path}"
