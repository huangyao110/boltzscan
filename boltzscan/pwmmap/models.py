from dataclasses import dataclass, field


@dataclass
class RefTf:
    ref_id: str
    source: str            # "cisbp" | "jaspar"
    species: str
    family: str
    dbid: str
    motif_ids: list = field(default_factory=list)
    uniprot_ids: list = field(default_factory=list)
