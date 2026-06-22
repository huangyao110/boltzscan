import pytest

@pytest.fixture
def tiny_domtbl(tmp_path):
    # domtblout columns: target(0) ... query_name(3) query_acc(4) ... dom_score(13)
    # ... env_from(19) env_to(20). PF00010 is a DBD (bHLH); PF99999 is not.
    rows = [
        "geneA - 120 HLH PF00010.30 54 1e-20 70 0 1 1 1e-22 1e-18 65 0 1 54 11 64 10 65 0.95 -",
        "geneB - 90 SomethingElse PF99999.1 40 1 5 0 1 1 1 1 4 0 1 40 5 44 4 45 0.9 -",
    ]
    p = tmp_path / "pfam.domtbl"
    p.write_text("# comment\n" + "\n".join(" ".join(r.split()) for r in rows) + "\n")
    return p

@pytest.fixture
def tiny_proteins(tmp_path):
    p = tmp_path / "prot.fasta"
    p.write_text(">geneA\n" + "M"*9 + "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRS" + "G"*60 + "\n"
                 ">geneB\nMKKKKKKKKK\n")
    return p
