from boltzscan.pwmmap import align


def test_parse_blast_rejects_low_coverage_and_cross_pfam():
    rows = [
        # full-length, same pfam, 50/55 identical -> kept, pct = 50/55
        ["qA__PF00010__0", "rB__PF00010__3", "90.9", "55", "50", "60", "58"],
        # short HSP (len 20 of qlen 60) -> coverage <0.8 -> dropped
        ["qA__PF00010__0", "rC__PF00010__1", "100.0", "20", "20", "60", "58"],
        # cross-pfam -> dropped
        ["qA__PF00010__0", "rD__PF00046__2", "95.0", "55", "52", "60", "58"],
    ]
    hits = align.parse_blast_rows(rows, min_cov=0.8)
    assert len(hits) == 1
    h = hits[0]
    assert h.query_tf == "qA" and h.ref_dbd_id == "rB__PF00010__3"
    assert h.query_pfam == "PF00010" and round(h.pct_id, 4) == round(50/55, 4)


import pytest


@pytest.mark.network  # not network, but external binary; excluded from default run
def test_real_blast_roundtrip(tmp_path):
    ref = tmp_path/"ref.fasta"; ref.write_text(">rX__PF00010__0\nMKRAHHNALERRRRDHIKDSFSSLRDSVPSLQ\n")
    qry = tmp_path/"q.fasta"; qry.write_text(">qX__PF00010__0\nMKRAHHNALERRRRDHIKDSFSSLRDSVPSLQ\n")
    bp, mk = align.resolve_blast_bins()
    db = align.make_blast_db(ref, tmp_path/"db", mk)
    hits = align.blast_dbd_pct_id(qry, db, bp, min_cov=0.8)
    assert hits and hits[0].pct_id == 1.0
