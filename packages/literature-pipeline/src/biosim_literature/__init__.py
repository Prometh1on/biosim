"""
biosim_literature — literature ingestion pipeline.

Public API:

    from biosim_literature import (
        IngestionPipeline,
        IngestionReport,
        IngestionOutcome,
        MockExtractor,
        MockFetcher,
        PubMedFetcher,
        RawExtraction,
        PaperMetadata,
        UnitNormalizer,
        NormalizationError,
        PlausibilityChecker,
        PlausibilityStatus,
    )
"""

from .curation_bridge import BridgeResult, enqueue_needs_curation_outcomes
from .extraction import (
    Extractor,
    ExtractionResult,
    MockExtractor,
    PaperMetadata,
    RawExtraction,
)
from .fetcher import FetchError, Fetcher, MockFetcher, PubMedFetcher
from .llm_extractor import LLMExtractor
from .normalizer import NormalizationError, NormalizedValue, UnitNormalizer
from .pipeline import (
    IngestionOutcome,
    IngestionPipeline,
    IngestionReport,
)
from .plausibility import PlausibilityChecker, PlausibilityResult, PlausibilityStatus
from .verification import CitationVerification, ClaimCheck, verify_citation

__all__ = [
    "BridgeResult",
    "CitationVerification",
    "ClaimCheck",
    "Extractor",
    "ExtractionResult",
    "FetchError",
    "Fetcher",
    "IngestionOutcome",
    "IngestionPipeline",
    "IngestionReport",
    "LLMExtractor",
    "MockExtractor",
    "MockFetcher",
    "NormalizationError",
    "NormalizedValue",
    "PaperMetadata",
    "PlausibilityChecker",
    "PlausibilityResult",
    "PlausibilityStatus",
    "PubMedFetcher",
    "RawExtraction",
    "UnitNormalizer",
    "enqueue_needs_curation_outcomes",
    "verify_citation",
]

__version__ = "0.1.0"
