"""Pfam DBD accession -> TF family, and family -> DBD %ID transfer cutoff.

Cutoffs follow the cisBP / Weirauch et al. 2014 family-specific DBD identity
thresholds for motif inference (fraction identity over the aligned domain).
Families without a calibrated value use DEFAULT_CUTOFF.
"""

# Canonical Pfam-acc -> family. Seeded from boltzscan.utils.find_tf.DBD_ACCS so
# both modules agree on which Pfam is which family.
PFAM_FAMILY = {
    "PF00046": "Homeobox", "PF05920": "Homeobox", "PF00249": "MYB",
    "PF13837": "MYB", "PF00010": "bHLH", "PF00847": "AP2", "PF02365": "NAC",
    "PF02362": "B3", "PF00096": "C2H2", "PF13894": "C2H2", "PF13912": "C2H2",
    "PF13465": "C2H2", "PF03106": "WRKY", "PF03514": "GRAS", "PF03101": "FAR1",
    "PF00642": "C3H", "PF00170": "bZIP", "PF07716": "bZIP", "PF00319": "MADS",
    "PF03195": "LBD", "PF02701": "Dof", "PF03634": "TCP", "PF00447": "HSF",
    "PF03110": "SBP", "PF00320": "GATA", "PF02319": "E2F", "PF00808": "NF-Y",
    "PF02045": "NF-Y", "PF08879": "GRF", "PF04770": "ZF-HD", "PF02042": "RWP-RK",
    "PF00643": "BBX", "PF06203": "CCT", "PF03859": "CAMTA", "PF04690": "YABBY",
    "PF04504": "GeBP", "PF03638": "CPP", "PF04873": "EIL", "PF06943": "LSD",
    "PF08536": "Whirly", "PF01422": "NF-X1", "PF04689": "S1Fa", "PF01698": "LFY",
    "PF17538": "LFY", "PF05142": "SRS", "PF05687": "BES1", "PF06217": "BBR-BPC",
    "PF28235": "VOZ",
}

DEFAULT_CUTOFF = 0.70

# Family-specific DBD %ID cutoffs (fraction). Values from cisBP/Weirauch 2014;
# families not listed inherit DEFAULT_CUTOFF.
FAMILY_CUTOFF = {
    "Homeobox": 0.70, "bHLH": 0.69, "bZIP": 0.70, "MYB": 0.62, "AP2": 0.59,
    "NAC": 0.62, "B3": 0.60, "C2H2": 0.78, "WRKY": 0.66, "Dof": 0.70,
    "GATA": 0.70, "MADS": 0.70, "HSF": 0.67, "NF-Y": 0.70, "E2F": 0.70,
    "SBP": 0.70, "TCP": 0.66, "GRAS": 0.60, "ZF-HD": 0.70, "LBD": 0.62,
}

def family_for_pfam(pfam_acc):
    return PFAM_FAMILY.get(pfam_acc)

def cutoff_for(pfam_acc, mode="family", global_thr=0.70):
    if mode == "global":
        return global_thr
    fam = family_for_pfam(pfam_acc)
    return FAMILY_CUTOFF.get(fam, DEFAULT_CUTOFF)
