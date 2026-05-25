from mhc_tp.report.data import parse_cluster_id, pcc_records, datatable_rows


def test_parse_cluster_id():
    assert parse_cluster_id("gibbs.1of3.mat") == ("1of3", 1, 3)
    assert parse_cluster_id("2of5.mat") == ("2of5", 2, 5)


def test_pcc_records_and_table():
    cd = {("gibbs.1of1.mat", "H2Kb"): 0.81, ("gibbs.1of2.mat", "H2Db"): 0.77}
    recs = pcc_records(cd)
    assert {"Cluster", "HLA", "Correlation"} <= set(recs[0])
    assert any(r["Cluster"] == "1of1" and r["HLA"] == "H2Kb" for r in recs)

    rows = datatable_rows(cd, kld_df=None)
    assert rows[0]["correlation"] >= rows[-1]["correlation"]  # sorted desc
    assert "kld" in rows[0]
