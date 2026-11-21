# SAVE TO: packages\literature-pipeline\tests\test_literature_curation_bridge.py
"""
Tests for biosim_literature.curation_bridge (MASTER_PLAN 5.2c).

Named test_literature_curation_bridge.py (not test_curation_bridge.py)
to avoid a pytest "import file mismatch" collision with
packages/parameter-db/tests/test_curation_bridge.py -- this repo's
flat test layout has no __init__.py per test directory to disambiguate
same-named modules, so basenames must be unique repo-wide.

Closes the third gap found scoping CurationQueue's own workflow (5.2):
IngestionPipeline already flags problematic extractions with
needs_curation=True, but nothing surfaced them to the curation queue.
These tests confirm the bridge actually does that, end to end through
a real IngestionPipeline (MockFetcher/MockExtractor -- no API key
needed) and a real CurationQueue -- not mocked.
"""
import pytest

from biosim_curation import CurationItemType, CurationPriority, CurationQueue
from biosim_id_resolver import EntityType, IDRegistry
from biosim_literature import (
    IngestionPipeline,
    MockExtractor,
    MockFetcher,
    RawExtraction,
    enqueue_needs_curation_outcomes,
)
from biosim_parameter_db import ParameterStore


@pytest.fixture
def registry_with_glucose():
    r = IDRegistry()
    r.register(
        name="glucose",
        entity_type=EntityType.METABOLITE,
        cross_references={},
        force_numeric_id=1,
    )
    return r


def make_pipeline(extractions, registry):
    return IngestionPipeline(
        fetcher=MockFetcher(),
        extractor=MockExtractor(extractions=extractions),
        store=ParameterStore(),
        registry=registry,
    )


class TestEnqueueNeedsCurationOutcomes:
    def test_flagged_extraction_is_enqueued(self, registry_with_glucose):
        extractions = [
            RawExtraction(
                parameter_kind="Km",
                subject_name="an unregistered protein",
                value=5.0,
                unit="uM",
            ),
        ]
        pipeline = make_pipeline(extractions, registry_with_glucose)
        report = pipeline.ingest_by_pubmed_id("11111111")
        queue = CurationQueue()

        result = enqueue_needs_curation_outcomes(report, queue)

        assert result.n_enqueued == 1
        assert len(queue.pending()) == 1

    def test_successful_extraction_is_not_enqueued(self, registry_with_glucose):
        extractions = [
            RawExtraction(
                # 5.0 uM, well within ABUNDANCE's plausible range
                # (1e-9 to 1e3 uM per plausibility.py's own
                # PLAUSIBLE_RANGES) -- a value like 5.0 mM (5000 uM)
                # would itself be flagged IMPLAUSIBLE, which is the
                # opposite of what this test needs to check.
                parameter_kind="abundance",
                subject_name="glucose",
                value=5.0,
                unit="uM",
            ),
        ]
        pipeline = make_pipeline(extractions, registry_with_glucose)
        report = pipeline.ingest_by_pubmed_id("11111111")
        queue = CurationQueue()

        result = enqueue_needs_curation_outcomes(report, queue)

        assert report.stored == 1
        assert result.n_enqueued == 0
        assert queue.pending() == []

    def test_enqueued_item_carries_the_real_extraction_details(self, registry_with_glucose):
        extractions = [
            RawExtraction(
                parameter_kind="not-a-real-kind",
                subject_name="glucose",
                value=42.0,
                unit="mM",
            ),
        ]
        pipeline = make_pipeline(extractions, registry_with_glucose)
        report = pipeline.ingest_by_pubmed_id("11111111")
        queue = CurationQueue()

        enqueue_needs_curation_outcomes(report, queue)

        item = queue.pending()[0]
        assert item.item_type == CurationItemType.EXTRACTION
        assert item.subject_name == "glucose"
        assert item.parameter_kind == "not-a-real-kind"
        assert item.value == 42.0
        assert item.unit == "mM"
        assert "bad_kind" in item.reason

    def test_unresolved_is_high_priority(self, registry_with_glucose):
        extractions = [
            RawExtraction(
                parameter_kind="Km",
                subject_name="an unregistered protein",
                value=5.0,
                unit="uM",
            ),
        ]
        pipeline = make_pipeline(extractions, registry_with_glucose)
        report = pipeline.ingest_by_pubmed_id("11111111")
        queue = CurationQueue()

        enqueue_needs_curation_outcomes(report, queue)

        item = queue.pending()[0]
        assert item.priority == CurationPriority.HIGH

    def test_bad_unit_is_enqueued_at_medium_priority_with_raw_value(self, registry_with_glucose):
        # Found by mandatory code review: the original test suite never
        # exercised BAD_UNIT/IMPLAUSIBLE at all, only UNRESOLVED/BAD_KIND
        # -- yet the value/unit fallback logic (normalized-if-present,
        # else raw) genuinely differs across these branches. BAD_UNIT
        # never reaches normalization at all (NormalizationError raised
        # before any normalized_value/normalized_unit exists), so the
        # enqueued item must carry the RAW value/unit, not None/None.
        extractions = [
            RawExtraction(
                parameter_kind="Km",
                subject_name="glucose",
                value=5.0,
                unit="not-a-real-unit",
            ),
        ]
        pipeline = make_pipeline(extractions, registry_with_glucose)
        report = pipeline.ingest_by_pubmed_id("11111111")
        queue = CurationQueue()

        enqueue_needs_curation_outcomes(report, queue)

        item = queue.pending()[0]
        assert item.priority == CurationPriority.MEDIUM
        assert item.value == 5.0
        assert item.unit == "not-a-real-unit"
        assert "bad_unit" in item.reason

    def test_implausible_fail_is_enqueued_at_high_priority_with_normalized_value(self, registry_with_glucose):
        # The other real code-review gap: IMPLAUSIBLE has two distinct
        # sub-cases in plausibility.py (a PHYSICAL_LIMITS FAIL -- e.g. a
        # negative Km, "almost certainly a data entry error" -- and a
        # PLAUSIBLE_RANGES WARNING -- a value outside the expected
        # range but not physically impossible). Both map to the same
        # IngestionOutcome.IMPLAUSIBLE and both DO reach normalization
        # first (unlike BAD_UNIT/BAD_KIND), so both correctly carry a
        # normalized_value/normalized_unit -- verified here for the FAIL
        # sub-case specifically (a negative Km, which fails the
        # PHYSICAL_LIMITS "must be > 0" check).
        extractions = [
            RawExtraction(
                parameter_kind="Km",
                subject_name="glucose",
                value=-5.0,
                unit="uM",
            ),
        ]
        pipeline = make_pipeline(extractions, registry_with_glucose)
        report = pipeline.ingest_by_pubmed_id("11111111")
        queue = CurationQueue()

        enqueue_needs_curation_outcomes(report, queue)

        item = queue.pending()[0]
        assert item.priority == CurationPriority.HIGH
        assert item.value == -5.0
        assert item.unit == "uM"
        assert "implausible" in item.reason

    def test_multiple_flagged_extractions_all_enqueued(self, registry_with_glucose):
        extractions = [
            RawExtraction(parameter_kind="Km", subject_name="unregistered A", value=1.0, unit="uM"),
            RawExtraction(parameter_kind="not-a-kind", subject_name="glucose", value=2.0, unit="mM"),
            RawExtraction(parameter_kind="abundance", subject_name="glucose", value=5.0, unit="uM"),  # succeeds
        ]
        pipeline = make_pipeline(extractions, registry_with_glucose)
        report = pipeline.ingest_by_pubmed_id("11111111")
        queue = CurationQueue()

        result = enqueue_needs_curation_outcomes(report, queue)

        assert result.n_enqueued == 2
        assert len(queue.pending()) == 2

    def test_idempotent_on_repeat_calls(self, registry_with_glucose):
        # Re-running the SAME report through the bridge twice must not
        # duplicate the same flagged extraction as two separate items --
        # item_id is built from the paper's own stable identifier plus
        # (subject_name, parameter_kind), not a list index.
        extractions = [
            RawExtraction(
                parameter_kind="Km",
                subject_name="an unregistered protein",
                value=5.0,
                unit="uM",
            ),
        ]
        pipeline = make_pipeline(extractions, registry_with_glucose)
        report = pipeline.ingest_by_pubmed_id("11111111")
        queue = CurationQueue()

        enqueue_needs_curation_outcomes(report, queue)
        enqueue_needs_curation_outcomes(report, queue)

        assert len(queue.pending()) == 1

    def test_empty_report_enqueues_nothing(self, registry_with_glucose):
        pipeline = make_pipeline([], registry_with_glucose)
        report = pipeline.ingest_by_pubmed_id("11111111")
        queue = CurationQueue()

        result = enqueue_needs_curation_outcomes(report, queue)

        assert result.n_enqueued == 0
        assert queue.pending() == []

    def test_summary_reports_n_enqueued(self, registry_with_glucose):
        # Matches the sibling bridges' own summary() shape (found by
        # mandatory code review as a consistency gap, added here).
        extractions = [
            RawExtraction(parameter_kind="Km", subject_name="unregistered A", value=1.0, unit="uM"),
        ]
        pipeline = make_pipeline(extractions, registry_with_glucose)
        report = pipeline.ingest_by_pubmed_id("11111111")
        queue = CurationQueue()

        result = enqueue_needs_curation_outcomes(report, queue)

        assert result.summary() == {"n_enqueued": 1}
