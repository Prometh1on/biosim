"""
Fluent builder for Parameter instances.

Using the builder avoids constructing deeply nested dataclasses directly,
which is error-prone and hard to read. The builder validates as it goes
and gives clear error messages when required fields are missing.

Usage:
    param = (
        ParameterBuilder()
        .for_subject("BIOSIM:PROTEIN:000001")
        .kind(ParameterKind.KM)
        .value(45.0, unit="uM")
        .with_uncertainty(lower=30.0, upper=65.0, uncertainty_type="95ci")
        .from_pubmed("12345678")
        .in_human_liver()
        .in_vitro()
        .confidence(ConfidenceTier.GOLD)
        .build()
    )
"""

from __future__ import annotations

import uuid
from typing import Optional

from .types import (
    ConfidenceTier,
    ExperimentalContext,
    ExperimentalSystem,
    Organism,
    Parameter,
    ParameterKind,
    ParameterValue,
    Provenance,
)


class ParameterBuilder:
    def __init__(self) -> None:
        self._subject_canonical_id: Optional[str] = None
        self._object_canonical_id: Optional[str] = None
        self._reaction_canonical_id: Optional[str] = None
        self._kind: Optional[ParameterKind] = None
        self._custom_kind_name: Optional[str] = None
        self._value: Optional[float] = None
        self._unit: Optional[str] = None
        self._uncertainty_lower: Optional[float] = None
        self._uncertainty_upper: Optional[float] = None
        self._uncertainty_type: str = "none"
        self._pubmed_id: Optional[str] = None
        self._doi: Optional[str] = None
        self._arxiv_id: Optional[str] = None
        self._source_description: str = ""
        self._extracted_by: str = "auto"
        self._figure_or_table: Optional[str] = None
        self._is_primary_source: bool = True
        self._verification_status: str = "UNCHECKED"
        self._verification_quote: str = ""
        self._organism: Organism = Organism.UNKNOWN
        self._system: ExperimentalSystem = ExperimentalSystem.UNKNOWN
        self._cell_type: Optional[str] = None
        self._tissue: Optional[str] = None
        self._temperature_celsius: Optional[float] = None
        self._ph: Optional[float] = None
        self._context_notes: str = ""
        self._confidence: ConfidenceTier = ConfidenceTier.SILVER
        self._curator_notes: str = ""
        self._supersedes: Optional[str] = None

    def for_subject(self, canonical_id: str) -> "ParameterBuilder":
        self._subject_canonical_id = canonical_id
        return self

    def for_object(self, canonical_id: str) -> "ParameterBuilder":
        self._object_canonical_id = canonical_id
        return self

    def for_reaction(self, canonical_id: str) -> "ParameterBuilder":
        self._reaction_canonical_id = canonical_id
        return self

    def kind(self, kind: ParameterKind, custom_name: Optional[str] = None) -> "ParameterBuilder":
        self._kind = kind
        self._custom_kind_name = custom_name
        return self

    def value(self, v: float, unit: str) -> "ParameterBuilder":
        self._value = v
        self._unit = unit
        return self

    def with_uncertainty(
        self,
        lower: float,
        upper: float,
        uncertainty_type: str = "95ci",
    ) -> "ParameterBuilder":
        self._uncertainty_lower = lower
        self._uncertainty_upper = upper
        self._uncertainty_type = uncertainty_type
        return self

    def from_pubmed(self, pubmed_id: str, figure_or_table: Optional[str] = None) -> "ParameterBuilder":
        self._pubmed_id = pubmed_id
        self._figure_or_table = figure_or_table
        return self

    def verification(self, status, quote: str = "") -> "ParameterBuilder":
        """Record whether the value was checked against the source text."""
        self._verification_status = getattr(status, "value", status)
        self._verification_quote = quote
        return self

    def from_doi(self, doi: str) -> "ParameterBuilder":
        self._doi = doi
        return self

    def from_description(self, description: str) -> "ParameterBuilder":
        self._source_description = description
        return self

    def extracted_by(self, name: str) -> "ParameterBuilder":
        self._extracted_by = name
        return self

    def secondary_source(self) -> "ParameterBuilder":
        self._is_primary_source = False
        return self

    def organism(self, organism: Organism) -> "ParameterBuilder":
        self._organism = organism
        return self

    def in_human(self) -> "ParameterBuilder":
        self._organism = Organism.HUMAN
        return self

    def in_human_liver(self) -> "ParameterBuilder":
        self._organism = Organism.HUMAN
        self._tissue = "liver"
        return self

    def in_human_muscle(self) -> "ParameterBuilder":
        self._organism = Organism.HUMAN
        self._tissue = "skeletal muscle"
        return self

    def in_human_adipose(self) -> "ParameterBuilder":
        self._organism = Organism.HUMAN
        self._tissue = "adipose"
        return self

    def cell_type(self, cell_type: str) -> "ParameterBuilder":
        self._cell_type = cell_type
        return self

    def experimental_system(self, system: ExperimentalSystem) -> "ParameterBuilder":
        self._system = system
        return self

    def in_vitro(self) -> "ParameterBuilder":
        self._system = ExperimentalSystem.IN_VITRO
        return self

    def in_vivo(self) -> "ParameterBuilder":
        self._system = ExperimentalSystem.IN_VIVO
        return self

    def clinical(self) -> "ParameterBuilder":
        self._system = ExperimentalSystem.CLINICAL
        return self

    def temperature(self, celsius: float) -> "ParameterBuilder":
        self._temperature_celsius = celsius
        return self

    def ph(self, ph: float) -> "ParameterBuilder":
        self._ph = ph
        return self

    def confidence(self, tier: ConfidenceTier) -> "ParameterBuilder":
        self._confidence = tier
        return self

    def gold(self) -> "ParameterBuilder":
        return self.confidence(ConfidenceTier.GOLD)

    def silver(self) -> "ParameterBuilder":
        return self.confidence(ConfidenceTier.SILVER)

    def bronze(self) -> "ParameterBuilder":
        return self.confidence(ConfidenceTier.BRONZE)

    def notes(self, notes: str) -> "ParameterBuilder":
        self._curator_notes = notes
        return self

    def experimental_notes(self, notes: str) -> "ParameterBuilder":
        """
        Set ExperimentalContext.notes (documented there as "critical").

        Regression: `_context_notes` was initialized in __init__ and read
        by build() into ExperimentalContext(notes=self._context_notes,
        ...), but no builder method ever set it -- every other
        ExperimentalContext field (organism, system, cell_type, tissue,
        temperature, ph) had one. Every Parameter built through this
        builder had its experimental-conditions notes permanently frozen
        at "" -- confirmed via grep, no assignment site existed anywhere
        in the package before this method.
        """
        self._context_notes = notes
        return self

    def supersedes(self, parameter_id: str) -> "ParameterBuilder":
        self._supersedes = parameter_id
        return self

    def build(self) -> Parameter:
        errors = []
        if not self._subject_canonical_id:
            errors.append("subject canonical ID is required (call .for_subject(...))")
        if not self._kind:
            errors.append("parameter kind is required (call .kind(...))")
        if self._value is None:
            errors.append("parameter value is required (call .value(...))")
        if not self._unit:
            errors.append("unit is required (call .value(v, unit=...))")
        if not any([self._pubmed_id, self._doi, self._arxiv_id, self._source_description]):
            errors.append("provenance is required (call .from_pubmed(), .from_doi(), or .from_description(...))")
        if errors:
            raise ValueError(
                f"Cannot build Parameter — {len(errors)} required field(s) missing:\n"
                + "\n".join(f"  - {e}" for e in errors)
            )

        return Parameter(
            parameter_id=str(uuid.uuid4()),
            kind=self._kind,
            custom_kind_name=self._custom_kind_name,
            value=ParameterValue(
                value=self._value,
                unit=self._unit,
                uncertainty_lower=self._uncertainty_lower,
                uncertainty_upper=self._uncertainty_upper,
                uncertainty_type=self._uncertainty_type,
            ),
            subject_canonical_id=self._subject_canonical_id,
            object_canonical_id=self._object_canonical_id,
            reaction_canonical_id=self._reaction_canonical_id,
            provenance=Provenance(
                pubmed_id=self._pubmed_id,
                doi=self._doi,
                arxiv_id=self._arxiv_id,
                source_description=self._source_description,
                extracted_by=self._extracted_by,
                figure_or_table=self._figure_or_table,
                is_primary_source=self._is_primary_source,
                verification_status=self._verification_status,
                verification_quote=self._verification_quote,
            ),
            experimental_context=ExperimentalContext(
                organism=self._organism,
                system=self._system,
                cell_type=self._cell_type,
                tissue=self._tissue,
                temperature_celsius=self._temperature_celsius,
                ph=self._ph,
                notes=self._context_notes,
            ),
            confidence_tier=self._confidence,
            curator_notes=self._curator_notes,
            supersedes_parameter_id=self._supersedes,
        )
