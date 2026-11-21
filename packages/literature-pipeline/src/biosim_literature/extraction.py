"""
Extraction types and the Extractor interface.

A RawExtraction is what comes out of the extraction step — either from an LLM
or from a human curator. It uses free-text entity names and raw units, both of
which are resolved and normalized downstream.

The Extractor interface is abstract. Two implementations exist:
    MockExtractor  — returns pre-defined results; used in tests
    LLMExtractor   — calls the Anthropic API; used in production
                     (requires API key; not used in tests)

This separation means the entire pipeline can be tested without any API calls.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class PaperMetadata:
    """
    Metadata for a scientific paper fetched from PubMed or similar.
    """
    pubmed_id: Optional[str]
    doi: Optional[str]
    title: str
    abstract: str
    authors: list[str]
    year: Optional[int]
    journal: Optional[str]
    full_text: Optional[str] = None   # Available for open-access papers
    is_primary_source: bool = True

    @property
    def best_id(self) -> str:
        """Return the best available identifier for provenance tracking."""
        return self.pubmed_id or self.doi or f"title:{self.title[:50]}"


@dataclass
class RawExtraction:
    """
    A single parameter extracted from a paper before normalization.

    Entity names are free text (e.g. "insulin receptor", "glucose") and
    must be resolved to CanonicalIDs by the pipeline's entity resolver.

    Units are raw (e.g. "nM", "min-1") and must be normalized by the
    UnitNormalizer before entering the ParameterStore.
    """
    parameter_kind: str              # String matching ParameterKind enum value
    subject_name: str                # Free-text entity name (to be resolved)
    value: float                     # Raw numeric value
    unit: str                        # Raw unit string (to be normalized)

    object_name: Optional[str] = None      # Free-text entity name (if applicable)
    organism: str = "unknown"
    cell_type: Optional[str] = None
    tissue: Optional[str] = None
    experimental_system: str = "unknown"
    temperature_celsius: Optional[float] = None
    figure_or_table: Optional[str] = None
    raw_text: str = ""                     # The source sentence(s) from the paper
    extractor_confidence: float = 1.0      # 0.0 to 1.0

    def __post_init__(self) -> None:
        if not (0.0 <= self.extractor_confidence <= 1.0):
            raise ValueError(
                f"extractor_confidence must be between 0.0 and 1.0, "
                f"got {self.extractor_confidence}"
            )


@dataclass
class ExtractionResult:
    """
    The result of running an extractor on a paper.
    """
    paper: PaperMetadata
    extractions: list[RawExtraction] = field(default_factory=list)
    extractor_name: str = ""
    errors: list[str] = field(default_factory=list)

    @property
    def succeeded(self) -> bool:
        return len(self.extractions) > 0 and len(self.errors) == 0

    @property
    def extraction_count(self) -> int:
        return len(self.extractions)


class Extractor(ABC):
    """
    Abstract base class for parameter extractors.

    All extractors take a PaperMetadata and return an ExtractionResult.
    The extraction step is deliberately isolated behind this interface so
    the rest of the pipeline is independent of extraction strategy.
    """

    @abstractmethod
    def extract(self, paper: PaperMetadata) -> ExtractionResult:
        """
        Extract raw parameters from a paper.

        Args:
            paper: Fetched paper metadata and optional full text

        Returns:
            ExtractionResult containing all found parameters
        """
        ...

    @property
    @abstractmethod
    def name(self) -> str:
        """Human-readable name of this extractor (for provenance records)."""
        ...


class MockExtractor(Extractor):
    """
    Mock extractor for use in tests and development.

    Returns pre-defined extractions regardless of paper content.
    Used to test the full pipeline without API calls.
    """

    def __init__(self, extractions: list[RawExtraction], errors: Optional[list[str]] = None):
        self._extractions = extractions
        self._errors = errors or []

    def extract(self, paper: PaperMetadata) -> ExtractionResult:
        return ExtractionResult(
            paper=paper,
            extractions=list(self._extractions),
            extractor_name=self.name,
            errors=list(self._errors),
        )

    @property
    def name(self) -> str:
        return "mock-extractor-v1"
