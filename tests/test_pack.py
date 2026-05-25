from hla_pepclust.db.pack import (
    build_pack_parquet,
    build_species_reference,
    classify_allele,
)
from hla_pepclust.refdata.parquet_io import read_reference

_HDR = ("#cmd\nLast position-specific scoring matrix computed, values are in halfbits\n"
        "A R N D C Q E G H I L K M F P S T W Y V\n")


def _score_row(v="0.5"):
    return "1 A " + " ".join([v] * 20) + "\n"


def test_classify():
    a = classify_allele("HLA-A02:352"); assert a.mhc_class=="I" and a.locus=="A" and a.species=="human" and a.formatted=="HLAA02352"
    assert classify_allele("H-2-Db").locus=="D" and classify_allele("H-2-Db").mhc_class=="I"
    assert classify_allele("H-2-IAb").locus=="IA" and classify_allele("H-2-IAb").mhc_class=="II"
    assert classify_allele("DRB1_0101").mhc_class=="II" and classify_allele("DRB1_0101").locus=="DR" and classify_allele("DRB1_0101").species=="human"
    assert classify_allele("HLA-E01:01").locus=="E"
    assert classify_allele("BoLA-DRB3_00101") is None
    assert classify_allele("Mamu-A01:01") is None


def test_build_pack_parquet(tmp_path):
    # Class I correlation matrices come from score_mat_el/<pseudo>.txt (halfbits).
    base = tmp_path / "all_logos"; (base/"mhc_names").mkdir(parents=True); (base/"score_mat_el").mkdir()
    (base/"mhc_names"/"PS_names.txt").write_text("HLA-A02:01 PS\nBoLA-x PS\n")
    (base/"score_mat_el"/"PS.txt").write_text(
        "#cmd\nLast position-specific scoring matrix computed, values are in halfbits\n"
        "A R N D C Q E G H I L K M F P S T W Y V\n"
        "1 A " + " ".join(["0.5"]*20) + "\n")
    n = build_pack_parquet(tmp_path, "I", "human", tmp_path/"h.parquet", "NetMHCpan-4.2")
    assert n == 1
    d = read_reference(tmp_path/"h.parquet")
    assert d["formatted"].iloc[0]=="HLAA0201" and d["mhc_class"].iloc[0]=="I"


def test_build_species_reference(tmp_path):
    c1 = tmp_path / "c1" / "all_logos"; (c1 / "mhc_names").mkdir(parents=True); (c1 / "score_mat_el").mkdir()
    (c1 / "mhc_names" / "P_names.txt").write_text("HLA-A02:01 P\n")
    (c1 / "score_mat_el" / "P.txt").write_text(_HDR + _score_row("0.5"))
    c2 = tmp_path / "c2" / "all_logos"; (c2 / "log_odds").mkdir(parents=True)
    (c2 / "pseudo_mhc_list").write_text("Q DRB1_0101\n")
    (c2 / "log_odds" / "Q.txt").write_text(_HDR + _score_row("0.3"))

    n_i, n_ii = build_species_reference("human", tmp_path / "c1", tmp_path / "c2", tmp_path / "human.parquet")
    assert n_i == 1 and n_ii == 1
    d = read_reference(tmp_path / "human.parquet")
    assert set(d["mhc_class"]) == {"I", "II"} and len(d) == 2
