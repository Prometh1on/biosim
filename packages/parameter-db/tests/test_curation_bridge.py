# SAVE TO: packages\parameter-db\tests\test_curation_bridge.py
"""
Tests for biosim_parameter_db.curation_bridge (MASTER_PLAN 5.2a).

Closes the gap found scoping CurationQueue's own workflow (5.2):
ParameterStore already detects real conflicts, but nothing surfaced
them to the curation queue. These tests confirm the bridge actually
does that, end to end through real ParameterStore/CurationQueue/
IDRegistry instances -- not mocked.
"""
import pytest

from biosim_curation import CurationQueue
from biosim_id_resolver import EntityType, IDRegistry
from biosim_parameter_db import (
    ConfidenceTier,
    Organism,
    ParameterBuilder,
    ParameterKind,
    ParameterStore,
    enqueue_unresolved_conflicts,
)

INSULIN_RECEPTOR_ID = "BIOSIM:PROTEIN:000001"


def make_km(value: float, pubmed_id: str = "12345678"):
    return (
        ParameterBuilder()
        .for_subject(INSULIN_RECEPTOR_ID)
        .kind(ParameterKind.KM)
        .value(value, unit="uM")
        .from_pubmed(pubmed_id)
        .organism(Organism.HUMAN)
        .confidence(ConfidenceTier.SILVER)
        .build()
    )


@pytest.fixture
def registry_with_insulin_receptor():
    r = IDRegistry()
    r.register(
        name="insulin receptor",
        entity_type=EntityType.PROTEIN,
        cross_references={},
        force_numeric_id=1,
    )
    return r


class TestEnqueueUnresolvedConflicts:
    def test_enqueues_a_real_conflict(self, registry_with_insulin_receptor):
        store = ParameterStore()
        store.add(make_km(10.0, pubmed_id="11111111"))
        store.add(make_km(50.0, pubmed_id="22222222"))  # 5x -- above threshold
        queue = CurationQueue()

        result = enqueue_unresolved_conflicts(store, queue, registry_with_insulin_receptor)

        assert result.n_conflicts_enqueued == 1
        pending = queue.pending()
        assert len(pending) == 1

    def test_subject_name_resolved_via_registry_not_raw_id(self, registry_with_insulin_receptor):
        store = ParameterStore()
        store.add(make_km(10.0, pubmed_id="11111111"))
        store.add(make_km(50.0, pubmed_id="22222222"))
        queue = CurationQueue()

        enqueue_unresolved_conflicts(store, queue, registry_with_insulin_receptor)

        item = queue.pending()[0]
        assert item.subject_name == "insulin receptor"
        assert item.subject_name != INSULIN_RECEPTOR_ID

    def test_unit_mismatch_skipped_not_silently_mislabeled(self, registry_with_insulin_receptor):
        # Mandatory code-reviewer finding (MEDIUM): ParameterStore's own
        # conflict detection never checks that the two compared
        # parameters share a unit -- constructing that exact scenario
        # directly (bypassing the literature-pipeline's own unit
        # normalization, which ParameterBuilder itself doesn't enforce).
        store = ParameterStore()
        store.add(
            ParameterBuilder()
            .for_subject(INSULIN_RECEPTOR_ID)
            .kind(ParameterKind.KM)
            .value(10.0, unit="uM")
            .from_pubmed("11111111")
            .organism(Organism.HUMAN)
            .confidence(ConfidenceTier.SILVER)
            .build()
        )
        store.add(
            ParameterBuilder()
            .for_subject(INSULIN_RECEPTOR_ID)
            .kind(ParameterKind.KM)
            .value(50.0, unit="nM")  # different unit, same subject/kind
            .from_pubmed("22222222")
            .organism(Organism.HUMAN)
            .confidence(ConfidenceTier.SILVER)
            .build()
        )
        queue = CurationQueue()

        result = enqueue_unresolved_conflicts(store, queue, registry_with_insulin_receptor)

        assert result.n_conflicts_enqueued == 0
        assert len(result.conflict_ids_skipped_unit_mismatch) == 1
        assert queue.pending() == []

    def test_falls_back_to_raw_canonical_id_when_malformed(self, registry_with_insulin_receptor):
        # Exercises the ValueError branch specifically (CanonicalID.
        # from_string() rejects a malformed string) -- distinct from
        # test_falls_back_to_raw_canonical_id_when_unregistered's
        # ResolutionError branch (a well-formed but unregistered ID).
        # ParameterBuilder.for_subject() performs no format validation,
        # so a malformed canonical ID is genuinely reachable here.
        store = ParameterStore()
        store.add(
            ParameterBuilder()
            .for_subject("not-a-valid-canonical-id")
            .kind(ParameterKind.KM)
            .value(10.0, unit="uM")
            .from_pubmed("11111111")
            .organism(Organism.HUMAN)
            .confidence(ConfidenceTier.SILVER)
            .build()
        )
        store.add(
            ParameterBuilder()
            .for_subject("not-a-valid-canonical-id")
            .kind(ParameterKind.KM)
            .value(50.0, unit="uM")
            .from_pubmed("22222222")
            .organism(Organism.HUMAN)
            .confidence(ConfidenceTier.SILVER)
            .build()
        )
        queue = CurationQueue()

        result = enqueue_unresolved_conflicts(store, queue, registry_with_insulin_receptor)

        assert result.n_conflicts_enqueued == 1
        assert queue.pending()[0].subject_name == "not-a-valid-canonical-id"

    def test_falls_back_to_raw_canonical_id_when_unregistered(self):
        # An empty registry -- the conflict must still be enqueued, not
        # silently dropped just because this one cosmetic lookup fails.
        store = ParameterStore()
        store.add(make_km(10.0, pubmed_id="11111111"))
        store.add(make_km(50.0, pubmed_id="22222222"))
        queue = CurationQueue()
        empty_registry = IDRegistry()

        result = enqueue_unresolved_conflicts(store, queue, empty_registry)

        assert result.n_conflicts_enqueued == 1
        item = queue.pending()[0]
        assert item.subject_name == INSULIN_RECEPTOR_ID

    def test_unit_and_value_and_fold_difference_correct(self, registry_with_insulin_receptor):
        store = ParameterStore()
        store.add(make_km(10.0, pubmed_id="11111111"))
        store.add(make_km(50.0, pubmed_id="22222222"))
        queue = CurationQueue()

        enqueue_unresolved_conflicts(store, queue, registry_with_insulin_receptor)

        item = queue.pending()[0]
        assert item.unit == "uM"
        assert item.parameter_kind == "Km"
        assert "5.0" in item.reason or "5.00" in item.reason

    def test_source_description_uses_pubmed_ids(self, registry_with_insulin_receptor):
        store = ParameterStore()
        store.add(make_km(10.0, pubmed_id="11111111"))
        store.add(make_km(50.0, pubmed_id="22222222"))
        queue = CurationQueue()

        enqueue_unresolved_conflicts(store, queue, registry_with_insulin_receptor)

        item = queue.pending()[0]
        assert "11111111" in item.source_description
        assert "22222222" in item.source_description

    def test_no_conflicts_enqueues_nothing(self, registry_with_insulin_receptor):
        store = ParameterStore()
        store.add(make_km(45.0, pubmed_id="11111111"))
        store.add(make_km(50.0, pubmed_id="22222222"))  # 1.11x -- below threshold
        queue = CurationQueue()

        result = enqueue_unresolved_conflicts(store, queue, registry_with_insulin_receptor)

        assert result.n_conflicts_enqueued == 0
        assert queue.pending() == []

    def test_resolved_conflicts_not_reenqueued(self, registry_with_insulin_receptor):
        store = ParameterStore()
        store.add(make_km(10.0, pubmed_id="11111111"))
        store.add(make_km(50.0, pubmed_id="22222222"))
        conflict_id = store.unresolved_conflicts()[0].conflict_id
        store.resolve_conflict(conflict_id, notes="Agreed on 10.0 -- better assay conditions")
        queue = CurationQueue()

        result = enqueue_unresolved_conflicts(store, queue, registry_with_insulin_receptor)

        assert result.n_conflicts_enqueued == 0
        assert queue.pending() == []

    def test_idempotent_on_repeat_calls(self, registry_with_insulin_receptor):
        store = ParameterStore()
        store.add(make_km(10.0, pubmed_id="11111111"))
        store.add(make_km(50.0, pubmed_id="22222222"))
        queue = CurationQueue()

        enqueue_unresolved_conflicts(store, queue, registry_with_insulin_receptor)
        enqueue_unresolved_conflicts(store, queue, registry_with_insulin_receptor)

        assert len(queue.pending()) == 1

    def test_multiple_conflicts_all_enqueued(self, registry_with_insulin_receptor):
        store = ParameterStore()
        store.add(make_km(10.0, pubmed_id="11111111"))
        store.add(make_km(50.0, pubmed_id="22222222"))  # conflict 1: Km
        store.add(
            ParameterBuilder()
            .for_subject(INSULIN_RECEPTOR_ID)
            .kind(ParameterKind.VMAX)
            .value(5.0, unit="uM/s")
            .from_pubmed("33333333")
            .organism(Organism.HUMAN)
            .confidence(ConfidenceTier.SILVER)
            .build()
        )
        store.add(
            ParameterBuilder()
            .for_subject(INSULIN_RECEPTOR_ID)
            .kind(ParameterKind.VMAX)
            .value(30.0, unit="uM/s")  # conflict 2: Vmax, 6x
            .from_pubmed("44444444")
            .organism(Organism.HUMAN)
            .confidence(ConfidenceTier.SILVER)
            .build()
        )
        queue = CurationQueue()

        result = enqueue_unresolved_conflicts(store, queue, registry_with_insulin_receptor)

        assert result.n_conflicts_enqueued == 2
        assert len(queue.pending()) == 2
