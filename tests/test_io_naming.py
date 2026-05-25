from mhc_tp.io.naming import format_allotype


def test_class_i_human():
    info = format_allotype("HLA-A*02:01", species="human")
    assert info.formatted == "A0201"
    assert info.mhc_class == "I"
    assert info.locus == "A"


def test_class_i_human_other_loci():
    assert format_allotype("HLA-B*07:02", species="human").locus == "B"
    assert format_allotype("HLA-C*07:485", species="human").formatted == "C07485"
    assert format_allotype("HLA-C*07:485", species="human").locus == "C"


def test_class_ii_human():
    info = format_allotype("HLA-DRB1*15:01", species="human")
    assert info.mhc_class == "II"
    assert info.locus == "DR"


def test_class_ii_human_other_loci():
    assert format_allotype("HLA-DQB1*06:02", species="human").locus == "DQ"
    assert format_allotype("HLA-DPB1*04:01", species="human").locus == "DP"


def test_class_i_mouse():
    info = format_allotype("H2-Kb", species="mouse")
    assert info.mhc_class == "I"
    assert info.locus == "K"


def test_class_i_mouse_formatted_matches_db():
    # Reference DB stores MHC*H2_Db -> H2Db; the H2 prefix is retained.
    info = format_allotype("MHC*H2_Db", species="mouse")
    assert info.formatted == "H2Db"
    assert info.mhc_class == "I"
    assert info.locus == "D"


def test_class_i_mouse_other_loci():
    assert format_allotype("H2-Ld", species="mouse").locus == "L"
    assert format_allotype("H2-Qa1", species="mouse").locus == "Q"


def test_class_ii_mouse():
    info = format_allotype("H2-IAb", species="mouse")
    assert info.mhc_class == "II"
    assert info.locus == "IA"
    assert format_allotype("H2-IEd", species="mouse").locus == "IE"


def test_raw_preserved():
    info = format_allotype("HLA-A*02:01", species="human")
    assert info.raw == "HLA-A*02:01"
