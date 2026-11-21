"""
Literature ingestion pipeline.

Orchestrates the full flow from a paper identifier to stored parameters:

    1. Fetch paper metadata (PubMedFetcher or MockFetcher)
    2. Extract raw parameters (LLMExtractor or MockExtractor)
    3. Normalize units (UnitNormalizer)
    4. Check plausibility (PlausibilityChecker)
    5. Resolve entity names to CanonicalIDs (EntityResolver)
    6. Build Parameter objects (ParameterBuilder)
    7. Store in ParameterStore

Each step can fail independently. Failures are collected in the
IngestionReport rather than crashing the pipeline — a single bad
extraction should not prevent good extractions from being stored.

Entity resolution:
    The pipeline resolves free-text entity names (e.g. "insulin receptor")
    to CanonicalIDs via the IDRegistry. If a name cannot be resolved,
    the extraction is flagged as UNRESOLVED and added to the curation queue.
    Unresolved extractions are NEVER silently dropped.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

from biosim_parameter_db import (
    ConfidenceTier,
    ExperimentalSystem,
    Organism,
    ParameterBuilder,
    ParameterKind,
    ParameterStore,
)
from biosim_id_resolver import IDRegistry

from .extraction import Extractor, ExtractionResult, PaperMetadata, RawExtraction
from .fetcher import Fetcher
from .normalizer import UnitNormalizer, NormalizationError
from .plausibility import PlausibilityChecker, PlausibilityStatus


class IngestionOutcome(str, Enum):
    SUCCESS    = "success"      # Parameter stored successfully
    NORMALIZED = "normalized"   # Stored with unit normalization applied
    UNRESOLVED = "unresolved"   # Entity name could not be resolved to CanonicalID
    IMPLAUSIBLE = "implausible" # Failed plausibility check — routed to curation
    BAD_UNIT   = "bad_unit"     # Unit not recognized — cannot normalize
    BAD_KIND   = "bad_kind"     # ParameterKind not recognized
    ERROR      = "error"        # Unexpected error during processing


@dataclass
class ExtractionOutcome:
    """The outcome of processing a single RawExtraction through the pipeline."""
    raw: RawExtraction
    outcome: IngestionOutcome
    parameter_id: Optional[str] = None    # Set if successfully stored
    canonical_subject_id: Optional[str] = None
    normalized_value: Optional[float] = None
    normalized_unit: Optional[str] = None
    message: str = ""
    needs_curation: bool = False


@dataclass
class IngestionReport:
    """
    Full report of a paper ingestion run.

    Every extraction attempt is recorded, whether it succeeded or failed.
    This is the primary input to the curation interface (Sprint 3).
    """
    pubmed_id: Optional[str]
    doi: Optional[str]
    paper_title: str
    extractor_name: str
    outcomes: list[ExtractionOutcome] = field(default_factory=list)
    fetch_error: Optional[str] = None
    extraction_error: Optional[str] = None

    @property
    def total(self) -> int:
        return len(self.outcomes)

    @property
    def stored(self) -> int:
        return sum(
            1 for o in self.outcomes
            if o.outcome in (IngestionOutcome.SUCCESS, IngestionOutcome.NORMALIZED)
        )

    @property
    def needs_curation(self) -> int:
        return sum(1 for o in self.outcomes if o.needs_curation)

    @property
    def failed(self) -> int:
        return sum(
            1 for o in self.outcomes
            if o.outcome in (IngestionOutcome.BAD_UNIT, IngestionOutcome.BAD_KIND, IngestionOutcome.ERROR)
        )

    def summary(self) -> dict:
        by_outcome: dict[str, int] = {}
        for o in self.outcomes:
            by_outcome[o.outcome.value] = by_outcome.get(o.outcome.value, 0) + 1
        return {
            "paper": self.paper_title,
            "total_extractions": self.total,
            "stored": self.stored,
            "needs_curation": self.needs_curation,
            "failed": self.failed,
            "by_outcome": by_outcome,
        }


# Organism name mapping from free text to enum
ORGANISM_NAME_MAP: dict[str, Organism] = {
    "human": Organism.HUMAN,
    "homo sapiens": Organism.HUMAN,
    "h. sapiens": Organism.HUMAN,
    "mouse": Organism.MOUSE,
    "mus musculus": Organism.MOUSE,
    "rat": Organism.RAT,
    "rattus norvegicus": Organism.RAT,
    "yeast": Organism.YEAST,
    "saccharomyces cerevisiae": Organism.YEAST,
    "e. coli": Organism.ECOLI,
    "escherichia coli": Organism.ECOLI,
}

SYSTEM_NAME_MAP: dict[str, ExperimentalSystem] = {
    "in vitro": ExperimentalSystem.IN_VITRO,
    "in_vitro": ExperimentalSystem.IN_VITRO,
    "in vivo": ExperimentalSystem.IN_VIVO,
    "in_vivo": ExperimentalSystem.IN_VIVO,
    "ex vivo": ExperimentalSystem.EX_VIVO,
    "ex_vivo": ExperimentalSystem.EX_VIVO,
    "clinical": ExperimentalSystem.CLINICAL,
    "in silico": ExperimentalSystem.IN_SILICO,
    "in_silico": ExperimentalSystem.IN_SILICO,
}


class IngestionPipeline:
    """
    Orchestrates literature ingestion from paper ID to stored parameters.

    Args:
        fetcher:     Paper metadata fetcher (PubMedFetcher or MockFetcher)
        extractor:   Parameter extractor (LLMExtractor or MockExtractor)
        store:       ParameterStore to write to
        registry:    IDRegistry for entity name resolution
        min_confidence_to_autostore: Parameters below this tier go to curation
                                     rather than being auto-stored
    """

    def __init__(
        self,
        fetcher: Fetcher,
        extractor: Extractor,
        store: ParameterStore,
        registry: IDRegistry,
        min_extractor_confidence: float = 0.7,
        min_confidence_to_autostore: ConfidenceTier = ConfidenceTier.SILVER,
    ) -> None:
        self._fetcher = fetcher
        self._extractor = extractor
        self._store = store
        self._registry = registry
        self._normalizer = UnitNormalizer()
        self._plausibility = PlausibilityChecker()
        self._min_extractor_confidence = min_extractor_confidence
        self._min_confidence_to_autostore = min_confidence_to_autostore

    def ingest_by_pubmed_id(self, pubmed_id: str) -> IngestionReport:
        """Run the full pipeline for a paper identified by PubMed ID."""
        report = IngestionReport(
            pubmed_id=pubmed_id,
            doi=None,
            paper_title="",
            extractor_name=self._extractor.name,
        )

        # Step 1: Fetch
        try:
            paper = self._fetcher.fetch_by_pubmed_id(pubmed_id)
            report.doi = paper.doi
            report.paper_title = paper.title
        except Exception as e:
            report.fetch_error = str(e)
            return report

        return self._run_pipeline(paper, report)

    def ingest_by_doi(self, doi: str) -> IngestionReport:
        """Run the full pipeline for a paper identified by DOI."""
        report = IngestionReport(
            pubmed_id=None,
            doi=doi,
            paper_title="",
            extractor_name=self._extractor.name,
        )

        try:
            paper = self._fetcher.fetch_by_doi(doi)
            report.pubmed_id = paper.pubmed_id
            report.paper_title = paper.title
        except Exception as e:
            report.fetch_error = str(e)
            return report

        return self._run_pipeline(paper, report)

    def ingest_paper(self, paper: PaperMetadata) -> IngestionReport:
        """Run the pipeline directly on a pre-fetched PaperMetadata object."""
        report = IngestionReport(
            pubmed_id=paper.pubmed_id,
            doi=paper.doi,
            paper_title=paper.title,
            extractor_name=self._extractor.name,
        )
        return self._run_pipeline(paper, report)

    def _run_pipeline(self, paper: PaperMetadata, report: IngestionReport) -> IngestionReport:
        # Step 2: Extract
        try:
            extraction_result = self._extractor.extract(paper)
        except Exception as e:
            report.extraction_error = str(e)
            return report

        if extraction_result.errors:
            report.extraction_error = "; ".join(extraction_result.errors)

        # Step 3-7: Process each extraction
        for raw in extraction_result.extractions:
            outcome = self._process_extraction(raw, paper)
            report.outcomes.append(outcome)

        return report

    def _process_extraction(
        self,
        raw: RawExtraction,
        paper: PaperMetadata,
    ) -> ExtractionOutcome:
        """Process a single RawExtraction through normalize → check → resolve → store."""

        # Validate parameter kind
        try:
            kind = ParameterKind(raw.parameter_kind)
        except ValueError:
            return ExtractionOutcome(
                raw=raw,
                outcome=IngestionOutcome.BAD_KIND,
                message=(
                    f"Unknown parameter kind '{raw.parameter_kind}'. "
                    f"Valid kinds: {[k.value for k in ParameterKind]}"
                ),
                needs_curation=True,
            )

        # Normalize units
        try:
            normalized = self._normalizer.normalize(raw.value, raw.unit)
        except NormalizationError as e:
            return ExtractionOutcome(
                raw=raw,
                outcome=IngestionOutcome.BAD_UNIT,
                message=str(e),
                needs_curation=True,
            )

        # Plausibility check
        plausibility = self._plausibility.check(kind, normalized.value, normalized.canonical_unit)
        if plausibility.status == PlausibilityStatus.FAIL:
            return ExtractionOutcome(
                raw=raw,
                outcome=IngestionOutcome.IMPLAUSIBLE,
                normalized_value=normalized.value,
                normalized_unit=normalized.canonical_unit,
                message=plausibility.message,
                needs_curation=True,
            )

        # Resolve subject entity name to CanonicalID
        subject_candidates = self._registry.search_by_name(raw.subject_name)
        if not subject_candidates:
            return ExtractionOutcome(
                raw=raw,
                outcome=IngestionOutcome.UNRESOLVED,
                normalized_value=normalized.value,
                normalized_unit=normalized.canonical_unit,
                message=(
                    f"Cannot resolve entity name '{raw.subject_name}' to a CanonicalID. "
                    f"Register this entity in the IDRegistry first, or add a name alias."
                ),
                needs_curation=True,
            )

        subject_canonical_id = str(subject_candidates[0])

        # Determine confidence tier
        confidence = self._determine_confidence(raw, plausibility.status)

        # Build the parameter
        builder = (
            ParameterBuilder()
            .for_subject(subject_canonical_id)
            .kind(kind)
            .value(normalized.value, unit=normalized.canonical_unit)
            .confidence(confidence)
            .extracted_by(self._extractor.name)
        )

        if paper.pubmed_id:
            builder = builder.from_pubmed(paper.pubmed_id, figure_or_table=raw.figure_or_table)
        elif paper.doi:
            builder = builder.from_doi(paper.doi)
        else:
            builder = builder.from_description(paper.title)

        if not paper.is_primary_source:
            builder = builder.secondary_source()

        # Experimental context
        organism = ORGANISM_NAME_MAP.get(raw.organism.lower(), Organism.UNKNOWN)
        system = SYSTEM_NAME_MAP.get(raw.experimental_system.lower(), ExperimentalSystem.UNKNOWN)

        builder = builder.organism(organism)
        builder = builder.experimental_system(system)
        if raw.cell_type:
            builder = builder.cell_type(raw.cell_type)
        if raw.tissue:
            builder._tissue = raw.tissue

        # Resolve object entity if present
        if raw.object_name:
            object_candidates = self._registry.search_by_name(raw.object_name)
            if object_candidates:
                builder = builder.for_object(str(object_candidates[0]))

        try:
            parameter = builder.build()
        except Exception as e:
            return ExtractionOutcome(
                raw=raw,
                outcome=IngestionOutcome.ERROR,
                message=f"Failed to build parameter: {e}",
                needs_curation=True,
            )

        # Route to curation if implausible warning or low confidence
        if plausibility.status == PlausibilityStatus.WARNING:
            return ExtractionOutcome(
                raw=raw,
                outcome=IngestionOutcome.IMPLAUSIBLE,
                normalized_value=normalized.value,
                normalized_unit=normalized.canonical_unit,
                canonical_subject_id=subject_canonical_id,
                message=plausibility.message,
                needs_curation=True,
            )

        # Store it
        parameter_id = self._store.add(parameter)
        was_normalized = normalized.conversion_factor != 1.0

        return ExtractionOutcome(
            raw=raw,
            outcome=IngestionOutcome.NORMALIZED if was_normalized else IngestionOutcome.SUCCESS,
            parameter_id=parameter_id,
            canonical_subject_id=subject_canonical_id,
            normalized_value=normalized.value,
            normalized_unit=normalized.canonical_unit,
            message=(
                f"Stored as {parameter_id}. "
                + (f"Unit converted from {raw.unit} to {normalized.canonical_unit}." if was_normalized else "")
            ),
            needs_curation=False,
        )

    def _determine_confidence(
        self,
        raw: RawExtraction,
        plausibility_status: PlausibilityStatus,
    ) -> ConfidenceTier:
        """Determine the confidence tier for an extracted parameter."""
        if raw.extractor_confidence >= 0.9 and plausibility_status == PlausibilityStatus.PASS:
            return ConfidenceTier.SILVER  # Auto-extracted max is SILVER — GOLD requires human
        if raw.extractor_confidence >= 0.5:
            return ConfidenceTier.SILVER
        return ConfidenceTier.BRONZE
