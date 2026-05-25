import textwrap
import pandas as pd
from hla_pepclust.db.construct import build_species_parquet
from hla_pepclust.refdata.parquet_io import read_reference


def test_build_species_parquet(tmp_path):
    mat = tmp_path / "A_02_01.txt"
    mat.write_text(textwrap.dedent("""\
        Pos A R N D C Q E G H I L K M F P S T W Y V
        1 P 0.1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.9
        """))
    db = pd.DataFrame({"allotypes": ["HLA-A*02:01"], "matrices_path": ["A_02_01.txt"]})
    db_csv = tmp_path / "human.db"
    db.to_csv(db_csv, index=False)

    out = tmp_path / "human.parquet"
    build_species_parquet(db_csv, tmp_path, "human", out, source="NetMHCpan-4.2")
    ref = read_reference(out)
    assert ref["formatted"].iloc[0] == "A0201"
    assert ref["mhc_class"].iloc[0] == "I"
    assert ref["n_positions"].iloc[0] == 1
    assert len(ref["matrix"].iloc[0]) == 20
