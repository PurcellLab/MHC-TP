"""Tests for Immunolyser-style allele display names."""

import pytest

from hla_pepclust.naming import pretty_allele


@pytest.mark.parametrize(
    "allotype,expected",
    [
        # human class I — colon and colon-less pack forms collapse the same way
        ("HLA-B39:124", "HLA_B39124"),
        ("HLA-B3942", "HLA_B3942"),
        ("HLA-B5309", "HLA_B5309"),
        ("HLA-A68117", "HLA_A68117"),
        ("HLA-A*02:01", "HLA_A0201"),
        # human class II — bare DRB/DQ/DP pack names gain the HLA_ prefix
        ("DRB1_0101", "HLA_DRB10101"),
        ("HLA-DRB1*01:01", "HLA_DRB10101"),
        # mouse class I
        ("H-2-Kb", "H2_Kb"),
        ("H-2-Db", "H2_Db"),
        ("H2-Kb", "H2_Kb"),
        # mouse class II
        ("H-2-IAb", "H2_IAb"),
        ("H-2-IEk", "H2_IEk"),
    ],
)
def test_pretty_allele(allotype, expected):
    assert pretty_allele(allotype) == expected


def test_pretty_allele_empty():
    assert pretty_allele("") == ""
    assert pretty_allele(None) == ""


def test_pretty_allele_idempotent_on_display_form():
    # Re-formatting an already-pretty human name is stable.
    assert pretty_allele("HLA_B39124") == "HLA_B39124"
    assert pretty_allele("H2_Kb") == "H2_Kb"
