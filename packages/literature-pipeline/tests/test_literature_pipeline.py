"""
Test suite for biosim_literature.

Coverage:
    - Unit normalization: all unit families, conversion accuracy, unknown units
    - Plausibility: pass/warning/fail for all parameter kinds, boundary conditions
    - MockExtractor and MockFetcher behaviour
    - Full pipeline: success, unresolved entity, bad unit, implausible value
    - IngestionReport statistics
    - Pipeline with real IDRegistry and ParameterStore integration
"""

import pytest

from biosim_literature import (
    IngestionOutcome,
    IngestionPipeline,
    MockExtractor,
    MockFetcher,
    NormalizationError,
    PaperMetadata,
    PlausibilityChecker,
    PlausibilityStatus,
    RawExtraction,
    UnitNormalizer,
)
from biosim_parameter_db import ExperimentalSystem, ParameterKind, ParameterStore
from biosim_id_resolver import EntityType, IDRegistry, OntologySource


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def normalizer() -> UnitNormalizer:
    return UnitNormalizer()


@pytest.fixture
def checker() -> PlausibilityChecker:
    return PlausibilityChecker()


@pytest.fixture
def registry_with_insulin_axis() -> IDRegistry:
    """IDRegistry pre-seeded with core insulin axis entities."""
    r = IDRegistry()
    r.register(
        name="insulin receptor",
        entity_type=EntityType.PROTEIN,
        cross_references={OntologySource.UNIPROT: ["P06213"]},
        force_numeric_id=1,
    )
    r.register(
        name="insulin",
        entity_type=EntityType.PROTEIN,
        cross_references={OntologySource.UNIPROT: ["P01308"]},
        force_numeric_id=2,
    )
    r.register(
        name="D-Glucose",
        entity_type=EntityType.METABOLITE,
        cross_references={OntologySource.CHEBI: ["CHEBI:17234"]},
        force_numeric_id=1,
    )
    r.register(
        name="IGF-1",
        entity_type=EntityType.PROTEIN,
        cross_references={OntologySource.UNIPROT: ["P05019"]},
        force_numeric_id=3,
    )
    return r


@pytest.fixture
def empty_store() -> ParameterStore:
    return ParameterStore()


@pytest.fixture
def pipeline(registry_with_insulin_axis, empty_store) -> IngestionPipeline:
    return IngestionPipeline(
        fetcher=MockFetcher(),
        extractor=MockExtractor([]),   # No extractions by default — override per test
        store=empty_store,
        registry=registry_with_insulin_axis,
    )


def make_pipeline(extractions, registry, store):
    return IngestionPipeline(
        fetcher=MockFetcher(),
        extractor=MockExtractor(extractions),
        store=store,
        registry=registry,
    )


def make_extraction(**kwargs) -> RawExtraction:
    defaults = dict(
        parameter_kind="Km",
        subject_name="insulin receptor",
        value=0.45,
        unit="nM",
        organism="human",
        experimental_system="in_vitro",
        extractor_confidence=0.9,
    )
    defaults.update(kwargs)
    return RawExtraction(**defaults)


# ---------------------------------------------------------------------------
# UnitNormalizer — concentration
# ---------------------------------------------------------------------------

class TestNormalizerConcentration:
    def test_um_passthrough(self, normalizer):
        result = normalizer.normalize_concentration(45.0, "uM")
        assert result.value == pytest.approx(45.0)
        assert result.canonical_unit == "uM"
        assert result.conversion_factor == 1.0

    def test_nm_to_um(self, normalizer):
        result = normalizer.normalize_concentration(450.0, "nM")
        assert result.value == pytest.approx(0.45)
        assert result.canonical_unit == "uM"

    def test_mm_to_um(self, normalizer):
        result = normalizer.normalize_concentration(1.0, "mM")
        assert result.value == pytest.approx(1000.0)
        assert result.canonical_unit == "uM"

    def test_m_to_um(self, normalizer):
        result = normalizer.normalize_concentration(1e-6, "M")
        assert result.value == pytest.approx(1.0)
        assert result.canonical_unit == "uM"

    def test_pm_to_um(self, normalizer):
        result = normalizer.normalize_concentration(1000.0, "pM")
        assert result.value == pytest.approx(1e-3)
        assert result.canonical_unit == "uM"

    def test_original_value_preserved(self, normalizer):
        result = normalizer.normalize_concentration(450.0, "nM")
        assert result.original_value == 450.0
        assert result.original_unit == "nM"


# ---------------------------------------------------------------------------
# UnitNormalizer — time
# ---------------------------------------------------------------------------

class TestNormalizerTime:
    def test_s_passthrough(self, normalizer):
        result = normalizer.normalize_time(60.0, "s")
        assert result.value == pytest.approx(60.0)
        assert result.canonical_unit == "s"

    def test_min_to_s(self, normalizer):
        result = normalizer.normalize_time(1.0, "min")
        assert result.value == pytest.approx(60.0)

    def test_h_to_s(self, normalizer):
        result = normalizer.normalize_time(1.0, "h")
        assert result.value == pytest.approx(3600.0)

    def test_ms_to_s(self, normalizer):
        result = normalizer.normalize_time(500.0, "ms")
        assert result.value == pytest.approx(0.5)

    def test_day_to_s(self, normalizer):
        result = normalizer.normalize_time(1.0, "d")
        assert result.value == pytest.approx(86400.0)


# ---------------------------------------------------------------------------
# UnitNormalizer — rate constants
# ---------------------------------------------------------------------------

class TestNormalizerRates:
    def test_s_minus_1_passthrough(self, normalizer):
        result = normalizer.normalize_first_order_rate(0.023, "s-1")
        assert result.value == pytest.approx(0.023)
        assert result.canonical_unit == "s-1"

    def test_min_minus_1_to_s_minus_1(self, normalizer):
        result = normalizer.normalize_first_order_rate(1.0, "min-1")
        assert result.value == pytest.approx(1/60)

    def test_h_minus_1_to_s_minus_1(self, normalizer):
        result = normalizer.normalize_first_order_rate(3600.0, "h-1")
        assert result.value == pytest.approx(1.0)

    def test_M_minus_1_s_minus_1_to_uM(self, normalizer):
        result = normalizer.normalize_second_order_rate(1e6, "M-1s-1")
        assert result.value == pytest.approx(1.0)
        assert result.canonical_unit == "uM-1s-1"

    def test_reaction_rate_uM_per_min(self, normalizer):
        result = normalizer.normalize_reaction_rate(60.0, "uM/min")
        assert result.value == pytest.approx(1.0)
        assert result.canonical_unit == "uM/s"

    def test_reaction_rate_mM_per_s(self, normalizer):
        result = normalizer.normalize_reaction_rate(1.0, "mM/s")
        assert result.value == pytest.approx(1000.0)


# ---------------------------------------------------------------------------
# UnitNormalizer — dimensionless and auto-detect
# ---------------------------------------------------------------------------

class TestNormalizerDimensionless:
    def test_dimensionless_passthrough(self, normalizer):
        result = normalizer.normalize_dimensionless(2.5, "dimensionless")
        assert result.value == pytest.approx(2.5)
        assert result.canonical_unit == "dimensionless"

    def test_fold_change_passthrough(self, normalizer):
        result = normalizer.normalize_dimensionless(3.2, "fold_change")
        assert result.value == pytest.approx(3.2)

    def test_unknown_dimensionless_raises(self, normalizer):
        with pytest.raises(NormalizationError):
            normalizer.normalize_dimensionless(1.0, "arbitrary_unit")


class TestNormalizerAutoDetect:
    def test_auto_detects_nM(self, normalizer):
        result = normalizer.normalize(450.0, "nM")
        assert result.canonical_unit == "uM"
        assert result.value == pytest.approx(0.45)

    def test_auto_detects_s_minus_1(self, normalizer):
        result = normalizer.normalize(0.023, "s-1")
        assert result.canonical_unit == "s-1"

    def test_auto_detects_dimensionless(self, normalizer):
        result = normalizer.normalize(2.5, "dimensionless")
        assert result.canonical_unit == "dimensionless"

    def test_unknown_unit_raises(self, normalizer):
        with pytest.raises(NormalizationError, match="Unknown unit"):
            normalizer.normalize(1.0, "parsecs")

    def test_is_known_unit_true(self, normalizer):
        assert normalizer.is_known_unit("nM") is True
        assert normalizer.is_known_unit("uM/s") is True
        assert normalizer.is_known_unit("dimensionless") is True

    def test_is_known_unit_false(self, normalizer):
        assert normalizer.is_known_unit("parsecs") is False

    def test_canonical_unit_for(self, normalizer):
        assert normalizer.canonical_unit_for("nM") == "uM"
        assert normalizer.canonical_unit_for("min-1") == "s-1"
        assert normalizer.canonical_unit_for("unknown") is None


# ---------------------------------------------------------------------------
# PlausibilityChecker
# ---------------------------------------------------------------------------

class TestPlausibilityChecker:
    def test_km_in_range_passes(self, checker):
        result = checker.check(ParameterKind.KM, 45.0, "uM")
        assert result.status == PlausibilityStatus.PASS

    def test_km_zero_fails(self, checker):
        result = checker.check(ParameterKind.KM, 0.0, "uM")
        assert result.status == PlausibilityStatus.FAIL

    def test_km_negative_fails(self, checker):
        result = checker.check(ParameterKind.KM, -1.0, "uM")
        assert result.status == PlausibilityStatus.FAIL

    def test_km_extremely_large_warns(self, checker):
        result = checker.check(ParameterKind.KM, 1e10, "uM")
        assert result.status == PlausibilityStatus.WARNING

    def test_hill_in_range_passes(self, checker):
        result = checker.check(ParameterKind.HILL, 2.0, "dimensionless")
        assert result.status == PlausibilityStatus.PASS

    def test_hill_zero_fails(self, checker):
        result = checker.check(ParameterKind.HILL, 0.0, "dimensionless")
        assert result.status == PlausibilityStatus.FAIL

    def test_hill_very_large_warns(self, checker):
        result = checker.check(ParameterKind.HILL, 50.0, "dimensionless")
        assert result.status == PlausibilityStatus.WARNING

    def test_vmax_in_range_passes(self, checker):
        result = checker.check(ParameterKind.VMAX, 0.5, "uM/s")
        assert result.status == PlausibilityStatus.PASS

    def test_kcat_in_range_passes(self, checker):
        result = checker.check(ParameterKind.KCAT, 100.0, "s-1")
        assert result.status == PlausibilityStatus.PASS

    def test_fold_change_can_be_negative(self, checker):
        result = checker.check(ParameterKind.FOLD_CHANGE, -5.0, "fold_change")
        assert result.status == PlausibilityStatus.PASS

    def test_passed_property(self, checker):
        result = checker.check(ParameterKind.KM, 45.0, "uM")
        assert result.passed is True

    def test_needs_review_on_warning(self, checker):
        result = checker.check(ParameterKind.KM, 1e10, "uM")
        assert result.needs_review is True

    def test_needs_review_on_fail(self, checker):
        result = checker.check(ParameterKind.KM, -1.0, "uM")
        assert result.needs_review is True

    def test_batch_check(self, checker):
        items = [
            (ParameterKind.KM, 45.0, "uM"),
            (ParameterKind.KM, -1.0, "uM"),
            (ParameterKind.HILL, 2.0, "dimensionless"),
        ]
        results = checker.check_batch(items)
        assert len(results) == 3
        assert results[0].status == PlausibilityStatus.PASS
        assert results[1].status == PlausibilityStatus.FAIL
        assert results[2].status == PlausibilityStatus.PASS

    def test_summary(self, checker):
        items = [
            (ParameterKind.KM, 45.0, "uM"),
            (ParameterKind.KM, -1.0, "uM"),
            (ParameterKind.KM, 1e10, "uM"),
        ]
        results = checker.check_batch(items)
        summary = checker.summary(results)
        assert summary["passed"] == 1
        assert summary["failed"] == 1
        assert summary["warnings"] == 1


# ---------------------------------------------------------------------------
# Full pipeline integration
# ---------------------------------------------------------------------------

class TestPipelineSuccess:
    def test_successful_ingestion(self, registry_with_insulin_axis, empty_store):
        extraction = make_extraction(
            subject_name="insulin receptor",
            value=0.45,
            unit="nM",    # Will be normalized to 0.00045 uM
        )
        pipeline = make_pipeline([extraction], registry_with_insulin_axis, empty_store)
        report = pipeline.ingest_by_pubmed_id("12345678")

        assert report.stored == 1
        assert report.total == 1
        assert report.needs_curation == 0

    def test_unit_normalization_applied(self, registry_with_insulin_axis, empty_store):
        extraction = make_extraction(value=450.0, unit="nM")
        pipeline = make_pipeline([extraction], registry_with_insulin_axis, empty_store)
        report = pipeline.ingest_by_pubmed_id("12345678")

        outcome = report.outcomes[0]
        assert outcome.outcome == IngestionOutcome.NORMALIZED
        assert outcome.normalized_value == pytest.approx(0.45)
        assert outcome.normalized_unit == "uM"

    def test_no_normalization_when_already_canonical(self, registry_with_insulin_axis, empty_store):
        extraction = make_extraction(value=45.0, unit="uM")
        pipeline = make_pipeline([extraction], registry_with_insulin_axis, empty_store)
        report = pipeline.ingest_by_pubmed_id("12345678")

        outcome = report.outcomes[0]
        assert outcome.outcome == IngestionOutcome.SUCCESS

    def test_multiple_extractions_from_one_paper(self, registry_with_insulin_axis, empty_store):
        extractions = [
            make_extraction(subject_name="insulin receptor", parameter_kind="Km", value=0.45, unit="nM"),
            make_extraction(subject_name="insulin receptor", parameter_kind="kcat", value=0.023, unit="s-1"),
        ]
        pipeline = make_pipeline(extractions, registry_with_insulin_axis, empty_store)
        report = pipeline.ingest_by_pubmed_id("12345678")

        assert report.stored == 2

    def test_parameter_stored_in_store(self, registry_with_insulin_axis, empty_store):
        extraction = make_extraction(value=45.0, unit="uM")
        pipeline = make_pipeline([extraction], registry_with_insulin_axis, empty_store)
        pipeline.ingest_by_pubmed_id("12345678")

        stats = empty_store.stats()
        assert stats["active_parameters"] == 1

    def test_experimental_system_is_stored_not_discarded(self, registry_with_insulin_axis, empty_store):
        # Regression: _process_extraction() computed `system` from
        # raw.experimental_system via SYSTEM_NAME_MAP but never called
        # builder.experimental_system(system) -- every parameter was
        # silently stored as ExperimentalSystem.UNKNOWN regardless of what
        # was actually extracted (found via ruff F841, the same
        # "computed but never wired in" bug class as the gi_meals kabs
        # and exercise-personalisation bugs -- see KNOWN BUGS FIXED #16
        # and the discarded-computation sweep in CLAUDE.md, which this
        # instance slipped past). make_extraction()'s own default
        # (experimental_system="in_vitro") means every existing pipeline
        # test had been silently exercising this bug for a long time
        # without ever asserting on the stored value.
        extraction = make_extraction(experimental_system="clinical", value=45.0, unit="uM")
        pipeline = make_pipeline([extraction], registry_with_insulin_axis, empty_store)
        pipeline.ingest_by_pubmed_id("12345678")

        subject_id = str(registry_with_insulin_axis.search_by_name("insulin receptor")[0])
        stored = empty_store.query(subject_canonical_id=subject_id, kind=ParameterKind.KM)
        assert len(stored) == 1
        assert stored[0].experimental_context.system == ExperimentalSystem.CLINICAL


class TestPipelineFailures:
    def test_unresolvable_entity_flagged(self, registry_with_insulin_axis, empty_store):
        extraction = make_extraction(subject_name="completely unknown protein XYZ")
        pipeline = make_pipeline([extraction], registry_with_insulin_axis, empty_store)
        report = pipeline.ingest_by_pubmed_id("12345678")

        assert report.stored == 0
        assert report.needs_curation == 1
        outcome = report.outcomes[0]
        assert outcome.outcome == IngestionOutcome.UNRESOLVED

    def test_bad_unit_flagged(self, registry_with_insulin_axis, empty_store):
        extraction = make_extraction(unit="parsecs")
        pipeline = make_pipeline([extraction], registry_with_insulin_axis, empty_store)
        report = pipeline.ingest_by_pubmed_id("12345678")

        assert report.stored == 0
        outcome = report.outcomes[0]
        assert outcome.outcome == IngestionOutcome.BAD_UNIT

    def test_unknown_parameter_kind_flagged(self, registry_with_insulin_axis, empty_store):
        extraction = make_extraction(parameter_kind="made_up_kind")
        pipeline = make_pipeline([extraction], registry_with_insulin_axis, empty_store)
        report = pipeline.ingest_by_pubmed_id("12345678")

        outcome = report.outcomes[0]
        assert outcome.outcome == IngestionOutcome.BAD_KIND

    def test_implausible_value_flagged(self, registry_with_insulin_axis, empty_store):
        extraction = make_extraction(value=-999.0, unit="uM")
        pipeline = make_pipeline([extraction], registry_with_insulin_axis, empty_store)
        report = pipeline.ingest_by_pubmed_id("12345678")

        assert report.stored == 0
        outcome = report.outcomes[0]
        assert outcome.outcome == IngestionOutcome.IMPLAUSIBLE

    def test_fetch_error_recorded(self, registry_with_insulin_axis, empty_store):
        from biosim_literature.fetcher import Fetcher, FetchError

        class FailingFetcher(Fetcher):
            def fetch_by_pubmed_id(self, pubmed_id):
                raise FetchError("Network unreachable")
            def fetch_by_doi(self, doi):
                raise FetchError("Network unreachable")

        pipeline = IngestionPipeline(
            fetcher=FailingFetcher(),
            extractor=MockExtractor([]),
            store=empty_store,
            registry=registry_with_insulin_axis,
        )
        report = pipeline.ingest_by_pubmed_id("99999999")
        assert report.fetch_error is not None
        assert report.stored == 0

    def test_mixed_success_and_failure(self, registry_with_insulin_axis, empty_store):
        extractions = [
            make_extraction(subject_name="insulin receptor", value=45.0, unit="uM"),  # success
            make_extraction(subject_name="unknown entity xyz", value=45.0, unit="uM"),  # unresolved
            make_extraction(subject_name="insulin receptor", value=45.0, unit="parsecs"),  # bad unit
        ]
        pipeline = make_pipeline(extractions, registry_with_insulin_axis, empty_store)
        report = pipeline.ingest_by_pubmed_id("12345678")

        assert report.stored == 1
        assert report.needs_curation == 2
        assert report.failed == 1


class TestIngestionReport:
    def test_summary_structure(self, registry_with_insulin_axis, empty_store):
        extraction = make_extraction(value=45.0, unit="uM")
        pipeline = make_pipeline([extraction], registry_with_insulin_axis, empty_store)
        report = pipeline.ingest_by_pubmed_id("12345678")
        summary = report.summary()

        assert "total_extractions" in summary
        assert "stored" in summary
        assert "needs_curation" in summary
        assert "by_outcome" in summary

    def test_ingest_by_doi(self, registry_with_insulin_axis, empty_store):
        extraction = make_extraction(value=45.0, unit="uM")
        pipeline = make_pipeline([extraction], registry_with_insulin_axis, empty_store)
        report = pipeline.ingest_by_doi("10.0000/test.001")
        assert report.stored == 1

    def test_ingest_paper_directly(self, registry_with_insulin_axis, empty_store):
        extraction = make_extraction(value=45.0, unit="uM")
        pipeline = make_pipeline([extraction], registry_with_insulin_axis, empty_store)
        paper = PaperMetadata(
            pubmed_id="11111111",
            doi=None,
            title="Test paper",
            abstract="An abstract.",
            authors=["Test Author"],
            year=2024,
            journal="Test Journal",
        )
        report = pipeline.ingest_paper(paper)
        assert report.paper_title == "Test paper"
        assert report.stored == 1
