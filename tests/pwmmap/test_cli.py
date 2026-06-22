from boltzscan.cli import _build_parser


def test_map_pwm_parses():
    p = _build_parser()
    ns = p.parse_args(["map-pwm", "-f", "cm.fasta", "-o", "data/pwms/cm_pwms",
                       "--domtbl", "x/pfam.domtbl", "--threshold-mode", "family"])
    assert ns.command == "map-pwm" and ns.proteins == "cm.fasta"


def test_build_pwm_refs_parses():
    p = _build_parser()
    ns = p.parse_args(["build-pwm-refs", "--refs", "data/pwms/_refs", "-c", "16"])
    assert ns.command == "build-pwm-refs" and ns.cpu == 16
