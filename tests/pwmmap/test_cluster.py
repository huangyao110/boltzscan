from boltzscan.pwmmap import cluster


def test_parse_tomtom_edges_symmetric_excludes_self_and_comments():
    tsv = (
        "Query_ID\tTarget_ID\toff\tp\tE\tq\n"
        "A\tA\t0\t1\t1\t1\n"      # self -> excluded
        "A\tB\t0\t1\t1\t0.01\n"   # edge
        "# Computing q-values\n"   # comment -> skipped
        "C\tD\t0\t1\t1\t0.02\n"
    )
    adj = cluster.parse_tomtom_edges(tsv)
    assert adj == {"A": {"B"}, "B": {"A"}, "C": {"D"}, "D": {"C"}}


def test_greedy_cluster_assigns_to_highest_support_rep_no_chaining():
    motifs = ["A", "B", "C", "D"]
    support = {"A": 10, "B": 5, "C": 8, "D": 1}
    # A~B and A~D similar; C is isolated; B and D are NOT directly similar
    adj = {"A": {"B", "D"}, "B": {"A"}, "D": {"A"}, "C": set()}
    rep = cluster.greedy_cluster(motifs, support, adj)
    # order by support: A(10), C(8), B(5), D(1)
    assert rep == {"A": "A", "B": "A", "C": "C", "D": "A"}


def test_build_motif_meta_family_and_support(tmp_path):
    idx = tmp_path / "ref_index.tsv"
    idx.write_text(
        "ref_id\tsource\tspecies\tfamily\tpfam_acc\tdbd_seq_id\tmotif_ids\n"
        "r1\tcisbp\tsp\tAP2\tPF00847\tr1__PF00847__0\tM1;M2\n"
        "r2\tcisbp\tsp\tAP2\tPF00847\tr2__PF00847__0\tM1\n"
        "r3\tjaspar\tsp\tbHLH\tPF00010\tr3__PF00010__0\tM3\n"
    )
    fam, support = cluster.build_motif_meta(idx)
    assert fam == {"M1": "AP2", "M2": "AP2", "M3": "bHLH"}
    assert support == {"M1": 2, "M2": 1, "M3": 1}


def test_load_clusters_missing_returns_empty(tmp_path):
    assert cluster.load_clusters(tmp_path) == {}
