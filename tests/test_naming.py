"""Tests for canonical allele display names (match reference logo titles)."""

import pytest

from mhc_tp.naming import pretty_allele


@pytest.mark.parametrize(
    "allotype,expected",
    [
        # The source allotype is the reference logo title, returned verbatim.
        ("HLA-A25:08", "HLA-A25:08"),
        ("HLA-A3301", "HLA-A3301"),
        ("HLA-B39:124", "HLA-B39:124"),
        ("HLA-C14:90", "HLA-C14:90"),
        # class II + mouse keep their source spelling
        ("DRB1_0101", "DRB1_0101"),
        ("H-2-IAb", "H-2-IAb"),
        ("H-2-Kb", "H-2-Kb"),
        # surrounding whitespace is trimmed
        ("  HLA-A25:08  ", "HLA-A25:08"),
    ],
)
def test_pretty_allele(allotype, expected):
    assert pretty_allele(allotype) == expected


def test_pretty_allele_empty():
    assert pretty_allele("") == ""
    assert pretty_allele(None) == ""
