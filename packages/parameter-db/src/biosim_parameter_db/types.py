"""
Parameter types for the biosim parameter database.

Every kinetic parameter in the platform is represented as a Parameter instance.
Parameters are immutable once created — updates create new versions, preserving
full history. No parameter is ever deleted, only superseded.

Confidence tiers:
    GOLD   — expert-validated against primary source
    SILVER — auto-extracted, passes plausibility checks, not yet expert-reviewed
    BRONZE — inferred, estimated, or from indirect evidence

Unit standard:
    Concentrations : micromolar (μM)
    Time           : seconds (s)
    Rate constants : per second (s⁻¹) or μM⁻¹s⁻¹ for bimolecular
    Dimensionless  : Hill coefficients, stoichiometric ratios
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional
from datetime import datetime, timezone


class ParameterKind(str, Enum):
    KM          = "Km"           # Michaelis constant (μM)
    VMAX        = "Vmax"         # Maximum reaction rate (μM/s)
    KD          = "Kd"           # Dissociation constant (μM)
    KI          = "Ki"           # Inhibition constant (μM)
    EC50        = "EC50"         # Half-maximal effective concentration (μM)
    IC50        = "IC50"         # Half-maximal inhibitory concentration (μM)
    KCAT        = "kcat"         # Catalytic rate constant (s⁻¹)
    HILL        = "hill"         # Hill coefficient (dimensionless)
    RATE        = "rate"         # Generic rate constant
    FOLD_CHANGE = "fold_change"  # Expression or activity fold change (dimensionless)
    HALF_LIFE   = "half_life"    # Molecular half-life (s)
    ABUNDANCE   = "abundance"    # Basal concentration (μM)
    CUSTOM      = "custom"       # Explicitly named custom parameter


class ConfidenceTier(str, Enum):
    GOLD   = "gold"    # Expert-validated
    SILVER = "silver"  # Auto-extracted, plausibility-checked
    BRONZE = "bronze"  # Inferred or estimated


class ExperimentalSystem(str, Enum):
    IN_VIVO       = "in_vivo"
    IN_VITRO      = "in_vitro"
    IN_SILICO     = "in_silico"
    EX_VIVO       = "ex_vivo"
    CLINICAL      = "clinical"
    UNKNOWN       = "unknown"


class Organism(str, Enum):
    HUMAN         = "homo_sapiens"
    MOUSE         = "mus_musculus"
    RAT           = "rattus_norvegicus"
    YEAST         = "saccharomyces_cerevisiae"
    ECOLI         = "escherichia_coli"
    OTHER         = "other"
    UNKNOWN       = "unknown"


@dataclass(frozen=True)
class ExperimentalContext:
    """
    The experimental conditions under which a parameter was measured.

    This is critical — a Km measured in HEK293 cells is not the same
    as one measured in primary hepatocytes. Context is stored alongside
    every parameter value and propagated to uncertainty calculations.
    """
    organism: Organism = Organism.UNKNOWN
    system: ExperimentalSystem = ExperimentalSystem.UNKNOWN
    cell_type: Optional[str] = None        # e.g. "HEK293", "primary hepatocyte"
    tissue: Optional[str] = None           # e.g. "liver", "skeletal muscle"
    temperature_celsius: Optional[float] = None
    ph: Optional[float] = None
    notes: str = ""


class VerificationStatus(str, Enum):
    """
    Whether a stored value was checked against its cited source's text.

    VERIFIED: the fetched text states the value (quote recorded).
    DERIVED: the value follows from numbers the fetched text states.
    PAPER_ON_TOPIC_VALUE_NOT_FOUND: right paper, value not in available text.
    WRONG_PAPER: the cited paper is about something else.
    CONTRADICTED: the paper states a different value for the same quantity.
    UNVERIFIED: no real source is cited (pubmed_id="UNVERIFIED").
    UNCHECKED: nobody has checked yet (the default).
    """
    VERIFIED = "VERIFIED"
    DERIVED = "DERIVED"
    PAPER_ON_TOPIC_VALUE_NOT_FOUND = "PAPER-ON-TOPIC-VALUE-NOT-FOUND"
    WRONG_PAPER = "WRONG-PAPER"
    CONTRADICTED = "CONTRADICTED"
    UNVERIFIED = "UNVERIFIED"
    UNCHECKED = "UNCHECKED"


@dataclass(frozen=True)
class Provenance:
    """
    Full provenance record for a parameter value.

    Every parameter must have a source. No parameter enters the database
    without a traceable origin — this is what makes the model reproducible.
    """
    pubmed_id: Optional[str] = None        # PubMed ID (preferred)
    doi: Optional[str] = None              # DOI as fallback
    arxiv_id: Optional[str] = None         # arXiv ID for preprints
    source_description: str = ""           # Human-readable source description
    extracted_by: str = ""                 # "auto" or curator name/ID
    extracted_at: datetime = field(
        default_factory=lambda: datetime.now(timezone.utc)
    )
    figure_or_table: Optional[str] = None  # e.g. "Table 2", "Figure 3B"
    is_primary_source: bool = True         # False if from a review paper
    # VerificationStatus value, plus the exact source sentence when the
    # status is VERIFIED/DERIVED/CONTRADICTED.
    verification_status: str = VerificationStatus.UNCHECKED.value
    verification_quote: str = ""

    def __post_init__(self) -> None:
        valid = {s.value for s in VerificationStatus}
        if self.verification_status not in valid:
            raise ValueError(
                f"verification_status must be one of {sorted(valid)}, "
                f"got {self.verification_status!r}"
            )
        if not any([self.pubmed_id, self.doi, self.arxiv_id, self.source_description]):
            raise ValueError(
                "Provenance requires at least one of: pubmed_id, doi, arxiv_id, "
                "or source_description. Parameters without traceable sources "
                "cannot enter the database."
            )


@dataclass(frozen=True)
class ParameterValue:
    """
    A single measured or estimated parameter value with full uncertainty.

    Values are always stored in the platform's canonical units.
    Unit conversion happens at ingestion time, never at query time.
    """
    value: float
    unit: str                              # Canonical unit string e.g. "uM", "s-1", "uM-1s-1"
    uncertainty_lower: Optional[float] = None   # Lower bound of confidence interval
    uncertainty_upper: Optional[float] = None   # Upper bound of confidence interval
    uncertainty_type: str = "none"         # "95ci", "sd", "sem", "range", "none"

    def __post_init__(self) -> None:
        if self.value < 0 and self.unit not in ("dimensionless", "fold_change"):
            # Negative kinetic parameters are physically meaningless
            # Exception: fold changes and dimensionless ratios can be negative
            raise ValueError(
                f"Parameter value {self.value} is negative. "
                f"Kinetic parameters must be non-negative. "
                f"If this is intentional (e.g. a log-transformed value), "
                f"set unit='dimensionless' explicitly."
            )

    @property
    def has_uncertainty(self) -> bool:
        return self.uncertainty_lower is not None and self.uncertainty_upper is not None

    @property
    def uncertainty_range(self) -> Optional[float]:
        if self.has_uncertainty:
            return self.uncertainty_upper - self.uncertainty_lower
        return None


@dataclass
class Parameter:
    """
    A kinetic parameter linking biological entities to a measured value.

    Every parameter connects:
        subject_canonical_id  — the entity being characterized (e.g. insulin receptor)
        object_canonical_id   — the entity it acts on or binds (e.g. insulin)
        reaction_canonical_id — the reaction context (optional)

    Parameters are versioned. When a parameter is updated (e.g. a better
    measurement is published), a new Parameter is created with
    supersedes_parameter_id pointing to the old one. The old one is never
    deleted — only marked as superseded.
    """
    parameter_id: str                      # UUID assigned at creation
    kind: ParameterKind
    value: ParameterValue
    subject_canonical_id: str             # CanonicalID string of the subject entity
    provenance: Provenance
    experimental_context: ExperimentalContext
    confidence_tier: ConfidenceTier = ConfidenceTier.SILVER

    object_canonical_id: Optional[str] = None    # CanonicalID string, if applicable
    reaction_canonical_id: Optional[str] = None  # CanonicalID string, if applicable
    custom_kind_name: Optional[str] = None        # Required if kind == CUSTOM
    supersedes_parameter_id: Optional[str] = None
    curator_notes: str = ""
    created_at: datetime = field(
        default_factory=lambda: datetime.now(timezone.utc)
    )
    is_superseded: bool = False

    def __post_init__(self) -> None:
        if self.kind == ParameterKind.CUSTOM and not self.custom_kind_name:
            raise ValueError(
                "Parameters with kind=CUSTOM must provide a custom_kind_name. "
                "This name must be unique within the parameter database and "
                "documented in docs/unit-standards.md."
            )

    def supersede(self, new_parameter_id: str) -> "Parameter":
        """Mark this parameter as superseded by a newer measurement."""
        self.is_superseded = True
        return self
