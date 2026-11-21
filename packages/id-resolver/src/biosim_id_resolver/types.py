"""
Canonical identifier types for the biosim platform.

Every biological entity in the platform has exactly one canonical ID.
All external ontology IDs (KEGG, Reactome, UniProt, ChEBI, HMDB) are
cross-referenced to this canonical ID. No module except this one is
permitted to hold raw external IDs — all inter-module communication
uses CanonicalID instances.

Canonical ID format:
    BIOSIM:<entity_type>:<zero-padded-integer>
    e.g. BIOSIM:PROTEIN:000142
         BIOSIM:METABOLITE:004891
         BIOSIM:PATHWAY:000023
         BIOSIM:GENE:001107
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional


class EntityType(str, Enum):
    PROTEIN = "PROTEIN"
    METABOLITE = "METABOLITE"
    PATHWAY = "PATHWAY"
    GENE = "GENE"
    REACTION = "REACTION"
    COMPARTMENT = "COMPARTMENT"


class OntologySource(str, Enum):
    """
    Supported external ontology sources, in priority order.

    Priority determines which source wins when resolving ambiguous mappings:
        UNIPROT > CHEBI > HMDB > REACTOME > KEGG

    KEGG is cross-reference only — never a primary resolution source.
    """
    UNIPROT = "uniprot"       # Primary: proteins
    CHEBI = "chebi"           # Primary: small molecules, metabolites
    HMDB = "hmdb"             # Primary: human metabolome, edge-case metabolites
    REACTOME = "reactome"     # Primary: pathways
    KEGG = "kegg"             # Cross-reference only — never primary


# Resolution priority order (index 0 = highest priority)
RESOLUTION_PRIORITY: list[OntologySource] = [
    OntologySource.UNIPROT,
    OntologySource.CHEBI,
    OntologySource.HMDB,
    OntologySource.REACTOME,
    OntologySource.KEGG,
]


@dataclass(frozen=True)
class ExternalID:
    """
    An identifier from an external ontology system.

    Examples:
        ExternalID(source=OntologySource.UNIPROT, accession="P06213")   # Insulin receptor
        ExternalID(source=OntologySource.CHEBI, accession="CHEBI:17234") # Glucose
        ExternalID(source=OntologySource.HMDB, accession="HMDB0000122")  # Glucose (HMDB)
        ExternalID(source=OntologySource.KEGG, accession="C00031")       # Glucose (KEGG)
    """
    source: OntologySource
    accession: str

    def __post_init__(self) -> None:
        if not self.accession or not self.accession.strip():
            raise ValueError(f"Accession cannot be empty for source {self.source}")
        # Normalize accession to uppercase stripped string
        object.__setattr__(self, "accession", self.accession.strip().upper())

    def __str__(self) -> str:
        return f"{self.source.value}:{self.accession}"


@dataclass(frozen=True)
class CanonicalID:
    """
    The platform's internal canonical identifier for a biological entity.

    This is the only ID format that crosses module boundaries.
    All modules receive and return CanonicalID instances.
    Raw external accessions never leave this package.
    """
    entity_type: EntityType
    numeric_id: int

    def __post_init__(self) -> None:
        if self.numeric_id < 0:
            raise ValueError(f"CanonicalID numeric_id must be non-negative, got {self.numeric_id}")

    def __str__(self) -> str:
        return f"BIOSIM:{self.entity_type.value}:{self.numeric_id:06d}"

    @classmethod
    def from_string(cls, s: str) -> "CanonicalID":
        """
        Parse a canonical ID string back into a CanonicalID instance.

        Args:
            s: A string of the form "BIOSIM:<TYPE>:<NUMBER>"

        Raises:
            ValueError: If the string is not a valid canonical ID.
        """
        parts = s.strip().split(":")
        if len(parts) != 3 or parts[0] != "BIOSIM":
            raise ValueError(f"Invalid canonical ID format: '{s}'. Expected 'BIOSIM:<TYPE>:<NUMBER>'")
        try:
            entity_type = EntityType(parts[1])
        except ValueError:
            raise ValueError(f"Unknown entity type '{parts[1]}' in canonical ID '{s}'")
        try:
            numeric_id = int(parts[2])
        except ValueError:
            raise ValueError(f"Non-integer numeric ID '{parts[2]}' in canonical ID '{s}'")
        return cls(entity_type=entity_type, numeric_id=numeric_id)


@dataclass
class ResolutionRecord:
    """
    The full resolution record for a canonical entity.

    Stores the canonical ID, all known external cross-references,
    human-readable name, and provenance metadata.
    """
    canonical_id: CanonicalID
    name: str
    description: str
    cross_references: dict[OntologySource, list[str]] = field(default_factory=dict)
    deprecated_canonical_ids: list[CanonicalID] = field(default_factory=list)

    def get_external_id(self, source: OntologySource) -> Optional[str]:
        """Return the primary accession for a given ontology source, or None."""
        refs = self.cross_references.get(source, [])
        return refs[0] if refs else None

    def all_external_accessions(self) -> list[ExternalID]:
        """Return all external IDs across all sources as a flat list."""
        result = []
        for source, accessions in self.cross_references.items():
            for acc in accessions:
                result.append(ExternalID(source=source, accession=acc))
        return result
