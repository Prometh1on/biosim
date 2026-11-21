"""
Test suite for biosim_parameter_db.

Coverage:
    - Parameter construction via builder and direct instantiation
    - Provenance validation — no orphan parameters
    - Value validation — negative kinetic parameters rejected
    - Store add, query, best, supersede operations
    - Conflict detection at correct fold-difference threshold
    - Priority ordering: GOLD > SILVER > BRONZE, HUMAN > others
    - Full version history preservation
    - Persistence: save and load round-trip
    - Thread safety under concurrent reads and writes
    - Builder validation — missing required fields give clear errors
    - Stats reporting
"""

import json
import threading
import tempfile
from pathlib import Path

import pytest

from biosim_parameter_db import (
    ConfidenceTier,
    ConflictRecord,
    ExperimentalContext,
    ExperimentalSystem,
    Organism,
    Parameter,
    ParameterBuilder,
    ParameterKind,
    ParameterStore,
    ParameterValue,
    Provenance,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

GLUCOSE_ID = "BIOSIM:METABOLITE:000001"
INSULIN_RECEPTOR_ID = "BIOSIM:PROTEIN:000001"
IGF1_ID = "BIOSIM:PROTEIN:000002"
GLYCOLYSIS_ID = "BIOSIM:PATHWAY:000001"


def make_km(
    subject: str = INSULIN_RECEPTOR_ID,
    value: float = 45.0,
    tier: ConfidenceTier = ConfidenceTier.SILVER,
    organism: Organism = Organism.HUMAN,
    system: ExperimentalSystem = ExperimentalSystem.IN_VITRO,
    pubmed_id: str = "12345678",
) -> Parameter:
    return (
        ParameterBuilder()
        .for_subject(subject)
        .kind(ParameterKind.KM)
        .value(value, unit="uM")
        .from_pubmed(pubmed_id)
        .organism(organism)
        .confidence(tier)
        .build()
    )


def make_kd(
    subject: str = INSULIN_RECEPTOR_ID,
    object_id: str = "BIOSIM:METABOLITE:000099",
    value: float = 1.0,
    tier: ConfidenceTier = ConfidenceTier.SILVER,
) -> Parameter:
    return (
        ParameterBuilder()
        .for_subject(subject)
        .for_object(object_id)
        .kind(ParameterKind.KD)
        .value(value, unit="nM")
        .from_pubmed("11112222")
        .in_human()
        .confidence(tier)
        .build()
    )


def make_vmax(
    subject: str = INSULIN_RECEPTOR_ID,
    value: float = 10.0,
    tier: ConfidenceTier = ConfidenceTier.SILVER,
) -> Parameter:
    return (
        ParameterBuilder()
        .for_subject(subject)
        .kind(ParameterKind.VMAX)
        .value(value, unit="uM/s")
        .from_pubmed("99999999")
        .in_human()
        .confidence(tier)
        .build()
    )


# ---------------------------------------------------------------------------
# ParameterValue tests
# ---------------------------------------------------------------------------

class TestParameterValue:
    def test_basic_construction(self):
        v = ParameterValue(value=45.0, unit="uM")
        assert v.value == 45.0
        assert v.unit == "uM"

    def test_has_uncertainty_false_by_default(self):
        v = ParameterValue(value=45.0, unit="uM")
        assert not v.has_uncertainty

    def test_has_uncertainty_true_when_set(self):
        v = ParameterValue(
            value=45.0, unit="uM",
            uncertainty_lower=30.0, uncertainty_upper=60.0,
            uncertainty_type="95ci",
        )
        assert v.has_uncertainty
        assert v.uncertainty_range == pytest.approx(30.0)

    def test_negative_value_rejected_for_kinetic(self):
        with pytest.raises(ValueError, match="negative"):
            ParameterValue(value=-5.0, unit="uM")

    def test_negative_value_allowed_for_fold_change(self):
        v = ParameterValue(value=-1.5, unit="fold_change")
        assert v.value == -1.5

    def test_negative_value_allowed_for_dimensionless(self):
        v = ParameterValue(value=-0.5, unit="dimensionless")
        assert v.value == -0.5


# ---------------------------------------------------------------------------
# Provenance tests
# ---------------------------------------------------------------------------

class TestProvenance:
    def test_pubmed_id_sufficient(self):
        p = Provenance(pubmed_id="12345678")
        assert p.pubmed_id == "12345678"

    def test_doi_sufficient(self):
        p = Provenance(doi="10.1000/xyz123")
        assert p.doi == "10.1000/xyz123"

    def test_description_sufficient(self):
        p = Provenance(source_description="Manual curation from review")
        assert p.source_description

    def test_no_source_raises(self):
        with pytest.raises(ValueError, match="Provenance requires"):
            Provenance()

    def test_arxiv_id_sufficient(self):
        p = Provenance(arxiv_id="2301.00001")
        assert p.arxiv_id == "2301.00001"


# ---------------------------------------------------------------------------
# ParameterBuilder tests
# ---------------------------------------------------------------------------

class TestParameterBuilder:
    def test_builds_valid_parameter(self):
        p = make_km()
        assert p.kind == ParameterKind.KM
        assert p.value.value == 45.0
        assert p.subject_canonical_id == INSULIN_RECEPTOR_ID

    def test_missing_subject_raises(self):
        with pytest.raises(ValueError, match="subject canonical ID"):
            (
                ParameterBuilder()
                .kind(ParameterKind.KM)
                .value(45.0, unit="uM")
                .from_pubmed("12345678")
                .build()
            )

    def test_missing_kind_raises(self):
        with pytest.raises(ValueError, match="parameter kind"):
            (
                ParameterBuilder()
                .for_subject(INSULIN_RECEPTOR_ID)
                .value(45.0, unit="uM")
                .from_pubmed("12345678")
                .build()
            )

    def test_missing_value_raises(self):
        with pytest.raises(ValueError, match="parameter value"):
            (
                ParameterBuilder()
                .for_subject(INSULIN_RECEPTOR_ID)
                .kind(ParameterKind.KM)
                .from_pubmed("12345678")
                .build()
            )

    def test_missing_provenance_raises(self):
        with pytest.raises(ValueError, match="provenance"):
            (
                ParameterBuilder()
                .for_subject(INSULIN_RECEPTOR_ID)
                .kind(ParameterKind.KM)
                .value(45.0, unit="uM")
                .build()
            )

    def test_multiple_missing_fields_all_listed(self):
        with pytest.raises(ValueError) as exc:
            ParameterBuilder().build()
        assert "subject canonical ID" in str(exc.value)
        assert "parameter kind" in str(exc.value)
        assert "provenance" in str(exc.value)

    def test_uncertainty_round_trip(self):
        p = (
            ParameterBuilder()
            .for_subject(INSULIN_RECEPTOR_ID)
            .kind(ParameterKind.KM)
            .value(45.0, unit="uM")
            .with_uncertainty(30.0, 60.0, "95ci")
            .from_pubmed("12345678")
            .build()
        )
        assert p.value.uncertainty_lower == 30.0
        assert p.value.uncertainty_upper == 60.0
        assert p.value.uncertainty_type == "95ci"

    def test_custom_kind_requires_name(self):
        with pytest.raises(ValueError, match="CUSTOM"):
            (
                ParameterBuilder()
                .for_subject(INSULIN_RECEPTOR_ID)
                .kind(ParameterKind.CUSTOM)
                .value(1.0, unit="dimensionless")
                .from_pubmed("12345678")
                .build()
            )

    def test_custom_kind_with_name_succeeds(self):
        p = (
            ParameterBuilder()
            .for_subject(INSULIN_RECEPTOR_ID)
            .kind(ParameterKind.CUSTOM, custom_name="insulin_clearance_rate")
            .value(0.05, unit="s-1")
            .from_pubmed("12345678")
            .build()
        )
        assert p.custom_kind_name == "insulin_clearance_rate"

    def test_convenience_tissue_methods(self):
        p = (
            ParameterBuilder()
            .for_subject(INSULIN_RECEPTOR_ID)
            .kind(ParameterKind.KM)
            .value(45.0, unit="uM")
            .from_pubmed("12345678")
            .in_human_liver()
            .build()
        )
        assert p.experimental_context.tissue == "liver"
        assert p.experimental_context.organism == Organism.HUMAN

    def test_parameter_id_is_unique(self):
        p1 = make_km()
        p2 = make_km()
        assert p1.parameter_id != p2.parameter_id

    def test_gold_silver_bronze_convenience(self):
        p_gold = make_km(tier=ConfidenceTier.GOLD)
        p_silver = make_km(tier=ConfidenceTier.SILVER)
        p_bronze = make_km(tier=ConfidenceTier.BRONZE)
        assert p_gold.confidence_tier == ConfidenceTier.GOLD
        assert p_silver.confidence_tier == ConfidenceTier.SILVER
        assert p_bronze.confidence_tier == ConfidenceTier.BRONZE

    def test_experimental_notes_reaches_experimental_context(self):
        # Regression: _context_notes was initialized and read by build()
        # into ExperimentalContext(notes=...), but no builder method ever
        # set it -- every Parameter built through this builder had its
        # experimental-conditions notes permanently frozen at "".
        p = (
            ParameterBuilder()
            .for_subject(INSULIN_RECEPTOR_ID)
            .kind(ParameterKind.KM)
            .value(45.0, unit="uM")
            .from_pubmed("12345678")
            .experimental_notes("Measured at 37C in serum-free medium.")
            .build()
        )
        assert p.experimental_context.notes == "Measured at 37C in serum-free medium."

    def test_experimental_notes_defaults_to_empty_string(self):
        p = make_km()
        assert p.experimental_context.notes == ""

    def test_experimental_notes_distinct_from_curator_notes(self):
        # notes() (curator_notes) and experimental_notes() (ExperimentalContext
        # .notes) are two different fields -- confirm setting one doesn't
        # silently also set the other.
        p = (
            ParameterBuilder()
            .for_subject(INSULIN_RECEPTOR_ID)
            .kind(ParameterKind.KM)
            .value(45.0, unit="uM")
            .from_pubmed("12345678")
            .notes("curator review pending")
            .experimental_notes("fasted overnight")
            .build()
        )
        assert p.curator_notes == "curator review pending"
        assert p.experimental_context.notes == "fasted overnight"


# ---------------------------------------------------------------------------
# ParameterStore — add and retrieve
# ---------------------------------------------------------------------------

class TestParameterStoreBasics:
    def test_add_and_get(self):
        store = ParameterStore()
        p = make_km()
        pid = store.add(p)
        retrieved = store.get(pid)
        assert retrieved.parameter_id == pid
        assert retrieved.value.value == 45.0

    def test_get_unknown_raises(self):
        store = ParameterStore()
        with pytest.raises(KeyError):
            store.get("nonexistent-id")

    def test_query_returns_matching(self):
        store = ParameterStore()
        store.add(make_km(value=45.0))
        results = store.query(INSULIN_RECEPTOR_ID, ParameterKind.KM)
        assert len(results) == 1
        assert results[0].value.value == 45.0

    def test_query_returns_empty_for_unknown(self):
        store = ParameterStore()
        results = store.query("BIOSIM:PROTEIN:999999", ParameterKind.KM)
        assert results == []

    def test_best_returns_top_priority(self):
        store = ParameterStore()
        store.add(make_km(value=45.0, tier=ConfidenceTier.SILVER))
        store.add(make_km(value=50.0, tier=ConfidenceTier.GOLD))
        best = store.best(INSULIN_RECEPTOR_ID, ParameterKind.KM)
        assert best.value.value == 50.0
        assert best.confidence_tier == ConfidenceTier.GOLD

    def test_best_returns_none_when_empty(self):
        store = ParameterStore()
        result = store.best("BIOSIM:PROTEIN:999999", ParameterKind.KM)
        assert result is None

    def test_best_respects_min_confidence(self):
        store = ParameterStore()
        store.add(make_km(value=45.0, tier=ConfidenceTier.BRONZE))
        result = store.best(
            INSULIN_RECEPTOR_ID,
            ParameterKind.KM,
            min_confidence=ConfidenceTier.SILVER,
        )
        assert result is None

    def test_get_all_for_subject(self):
        store = ParameterStore()
        store.add(make_km())
        store.add(make_vmax())
        all_params = store.get_all_for_subject(INSULIN_RECEPTOR_ID)
        assert ParameterKind.KM.value in all_params
        assert ParameterKind.VMAX.value in all_params


# ---------------------------------------------------------------------------
# Priority ordering
# ---------------------------------------------------------------------------

class TestPriorityOrdering:
    def test_gold_beats_silver(self):
        store = ParameterStore()
        store.add(make_km(value=45.0, tier=ConfidenceTier.SILVER))
        store.add(make_km(value=50.0, tier=ConfidenceTier.GOLD))
        results = store.query(INSULIN_RECEPTOR_ID, ParameterKind.KM)
        assert results[0].confidence_tier == ConfidenceTier.GOLD

    def test_silver_beats_bronze(self):
        store = ParameterStore()
        store.add(make_km(value=45.0, tier=ConfidenceTier.BRONZE))
        store.add(make_km(value=50.0, tier=ConfidenceTier.SILVER))
        results = store.query(INSULIN_RECEPTOR_ID, ParameterKind.KM)
        assert results[0].confidence_tier == ConfidenceTier.SILVER

    def test_human_beats_mouse_at_same_confidence(self):
        store = ParameterStore()
        store.add(make_km(value=45.0, tier=ConfidenceTier.SILVER, organism=Organism.MOUSE))
        store.add(make_km(value=50.0, tier=ConfidenceTier.SILVER, organism=Organism.HUMAN))
        results = store.query(INSULIN_RECEPTOR_ID, ParameterKind.KM)
        assert results[0].experimental_context.organism == Organism.HUMAN

    def test_organism_filter(self):
        store = ParameterStore()
        store.add(make_km(value=45.0, organism=Organism.MOUSE))
        store.add(make_km(value=50.0, organism=Organism.HUMAN))
        results = store.query(
            INSULIN_RECEPTOR_ID,
            ParameterKind.KM,
            organism=Organism.MOUSE,
        )
        assert len(results) == 1
        assert results[0].experimental_context.organism == Organism.MOUSE


# ---------------------------------------------------------------------------
# Conflict detection
# ---------------------------------------------------------------------------

class TestConflictDetection:
    def test_no_conflict_within_threshold(self):
        store = ParameterStore()
        store.add(make_km(value=45.0))
        store.add(make_km(value=60.0))  # 60/45 = 1.33 — below 2x threshold
        assert len(store.unresolved_conflicts()) == 0

    def test_conflict_detected_above_threshold(self):
        store = ParameterStore()
        store.add(make_km(value=10.0))
        store.add(make_km(value=50.0))  # 50/10 = 5x — above threshold
        conflicts = store.unresolved_conflicts()
        assert len(conflicts) == 1
        assert conflicts[0].fold_difference == pytest.approx(5.0)

    def test_conflict_sorted_by_fold_difference(self):
        store = ParameterStore()
        # First conflict: 10x
        store.add(make_km(value=10.0, subject=INSULIN_RECEPTOR_ID))
        store.add(make_km(value=100.0, subject=INSULIN_RECEPTOR_ID))
        # Second conflict: 3x
        store.add(make_km(value=10.0, subject=IGF1_ID))
        store.add(make_km(value=30.0, subject=IGF1_ID))
        conflicts = store.unresolved_conflicts()
        assert conflicts[0].fold_difference > conflicts[1].fold_difference

    def test_resolve_conflict(self):
        store = ParameterStore()
        store.add(make_km(value=10.0))
        store.add(make_km(value=100.0))
        conflicts = store.unresolved_conflicts()
        assert len(conflicts) == 1
        store.resolve_conflict(conflicts[0].conflict_id, notes="Agreed on 10.0 — better assay conditions")
        assert len(store.unresolved_conflicts()) == 0

    def test_both_conflicting_parameters_retained(self):
        store = ParameterStore()
        p1 = make_km(value=10.0)
        p2 = make_km(value=100.0)
        id1 = store.add(p1)
        id2 = store.add(p2)
        # Both should be retrievable
        assert store.get(id1).value.value == 10.0
        assert store.get(id2).value.value == 100.0


class TestCustomKindScopedParameters:
    """
    Regression (MASTER_PLAN 5.1): conflict detection keyed on
    (subject, kind, object, reaction) only, silently ignoring
    custom_kind_name for kind=CUSTOM parameters. ParameterKind.CUSTOM's
    whole reason to exist is to hold many DIFFERENT named quantities
    under one bare enum value (custom_kind_name is the real
    disambiguator) -- two real, semantically unrelated CUSTOM
    parameters for the same subject (e.g. an absorption rate constant
    and a per-Unit scaling factor, both scoped to the same drug/entity)
    would have been falsely flagged as rival measurements of "the
    same" quantity. Found while curating exogenous-insulin PK constants
    into the store: 4 real CUSTOM parameters for subject="insulin"
    produced 5 nonsensical conflicts (fold differences up to ~61,500x)
    before this fix. Exactly the same bug class TestObjectScopedParameters
    already covers for object_canonical_id, just for custom_kind_name.
    """

    def _custom_param(self, custom_name: str, value: float) -> Parameter:
        return (
            ParameterBuilder()
            .for_subject(INSULIN_RECEPTOR_ID)
            .kind(ParameterKind.CUSTOM, custom_name=custom_name)
            .value(value, unit="dimensionless")
            .from_pubmed("12345678")
            .in_human()
            .confidence(ConfidenceTier.BRONZE)
            .build()
        )

    def test_different_custom_kind_names_do_not_conflict(self):
        store = ParameterStore()
        store.add(self._custom_param("rate_constant", 0.01))
        store.add(self._custom_param("scaling_factor", 700.0))  # ~70,000x apart
        assert store.unresolved_conflicts() == []

    def test_same_custom_kind_name_still_conflicts(self):
        # The fix must not disable conflict detection entirely for
        # CUSTOM-kind parameters -- only scope it correctly.
        store = ParameterStore()
        store.add(self._custom_param("rate_constant", 1.0))
        store.add(self._custom_param("rate_constant", 10.0))
        conflicts = store.unresolved_conflicts()
        assert len(conflicts) == 1
        assert conflicts[0].fold_difference == pytest.approx(10.0)

    def test_query_without_custom_kind_name_still_returns_all_matching_kind(self):
        # query()/best()'s own lookup semantics are DELIBERATELY
        # untouched by this fix (ParameterLoader.load_pathway_params()'s
        # real p3 mapping calls store.best() with no custom_kind_name
        # filter at all) -- confirms query(kind=CUSTOM) with no
        # custom_kind_name argument still returns every CUSTOM parameter
        # for that subject, same coarse behaviour as before this fix.
        store = ParameterStore()
        store.add(self._custom_param("rate_constant", 0.01))
        store.add(self._custom_param("scaling_factor", 700.0))
        results = store.query(INSULIN_RECEPTOR_ID, ParameterKind.CUSTOM)
        assert len(results) == 2


class TestObjectScopedParameters:
    """
    Regression: the store's index/conflict-detection/query all keyed on
    (subject, kind) only, silently ignoring object_canonical_id --
    Parameter's own docstring defines it as a first-class disambiguator
    (e.g. a receptor's Kd for one ligand vs. another), and
    ParameterBuilder.for_object() is already used in production seed data
    (the insulin receptor's Kd for insulin). Two real, independent Kd
    measurements for the SAME receptor but DIFFERENT ligands -- e.g. the
    insulin receptor's well-documented cross-reactivity with IGF-1, both
    ligands this platform already tracks as species -- would have been
    falsely flagged as a 2x-fold conflict between rival measurements of
    "the same" quantity, and best()/query() had no way to even ask for
    one specifically: whichever ranked highest would silently win for
    every caller, with no way to request the other.
    """
    def test_different_objects_do_not_conflict(self):
        store = ParameterStore()
        store.add(make_kd(object_id="BIOSIM:PROTEIN:000002", value=1.0))    # Kd for IGF-1
        store.add(make_kd(object_id="BIOSIM:METABOLITE:000099", value=50.0))  # Kd for insulin -- 50x apart
        assert store.unresolved_conflicts() == []

    def test_same_object_still_conflicts(self):
        # Scoping by object must not silently disable conflict detection
        # when the object genuinely matches.
        store = ParameterStore()
        store.add(make_kd(object_id="BIOSIM:PROTEIN:000002", value=1.0))
        store.add(make_kd(object_id="BIOSIM:PROTEIN:000002", value=50.0))
        conflicts = store.unresolved_conflicts()
        assert len(conflicts) == 1
        assert conflicts[0].fold_difference == pytest.approx(50.0)

    def test_query_without_object_does_not_return_object_scoped_parameter(self):
        store = ParameterStore()
        store.add(make_kd(object_id="BIOSIM:PROTEIN:000002"))
        assert store.query(INSULIN_RECEPTOR_ID, ParameterKind.KD) == []
        assert store.best(INSULIN_RECEPTOR_ID, ParameterKind.KD) is None

    def test_query_with_matching_object_id_retrieves_it(self):
        store = ParameterStore()
        store.add(make_kd(object_id="BIOSIM:PROTEIN:000002", value=1.0))
        results = store.query(
            INSULIN_RECEPTOR_ID, ParameterKind.KD, object_canonical_id="BIOSIM:PROTEIN:000002",
        )
        assert len(results) == 1
        assert results[0].value.value == pytest.approx(1.0)

    def test_query_disambiguates_two_objects_for_the_same_subject(self):
        store = ParameterStore()
        store.add(make_kd(object_id="BIOSIM:PROTEIN:000002", value=1.0))     # IGF-1
        store.add(make_kd(object_id="BIOSIM:METABOLITE:000099", value=50.0))  # insulin
        igf1_kd = store.best(
            INSULIN_RECEPTOR_ID, ParameterKind.KD, object_canonical_id="BIOSIM:PROTEIN:000002",
        )
        insulin_kd = store.best(
            INSULIN_RECEPTOR_ID, ParameterKind.KD, object_canonical_id="BIOSIM:METABOLITE:000099",
        )
        assert igf1_kd.value.value == pytest.approx(1.0)
        assert insulin_kd.value.value == pytest.approx(50.0)

    def test_save_load_roundtrip_preserves_object_scoped_index(self, tmp_path):
        store = ParameterStore()
        store.add(make_kd(object_id="BIOSIM:PROTEIN:000002", value=1.0))
        store.add(make_kd(object_id="BIOSIM:METABOLITE:000099", value=50.0))
        path = tmp_path / "store.json"
        store.save(path)
        loaded = ParameterStore.load(path)
        assert loaded.unresolved_conflicts() == []
        assert loaded.best(
            INSULIN_RECEPTOR_ID, ParameterKind.KD, object_canonical_id="BIOSIM:PROTEIN:000002",
        ).value.value == pytest.approx(1.0)

    def test_get_all_for_subject_returns_every_scope_not_just_the_last(self):
        # Regression: get_all_for_subject() indexed its result dict by the
        # bare kind string and did `result[kind_str] = sorted(...)` inside a
        # loop over ALL index keys for the subject -- a plain assignment,
        # not an accumulation. A subject+kind spanning multiple
        # object/reaction scopes (a real, supported case per
        # ParameterBuilder.for_object()) silently kept only whichever
        # scope's bucket the loop visited last, dropping the others with
        # no warning. Reproduced directly: 3 stored KD params for one
        # subject (one subject-only, two object-scoped) -> 1 returned.
        store = ParameterStore()
        store.add(make_kd(object_id=None, value=1.0))
        store.add(make_kd(object_id="BIOSIM:PROTEIN:000002", value=2.0))
        store.add(make_kd(object_id="BIOSIM:METABOLITE:000099", value=3.0))
        all_params = store.get_all_for_subject(INSULIN_RECEPTOR_ID)
        assert ParameterKind.KD.value in all_params
        values = sorted(p.value.value for p in all_params[ParameterKind.KD.value])
        assert values == [1.0, 2.0, 3.0]

    def test_get_all_for_subject_ties_break_by_global_insertion_order(self):
        # A naive fix (accumulate, then `sorted(key=self._priority_score,
        # reverse=True)`, with no explicit insertion-sequence tiebreak)
        # still fails this -- NOT because `reverse=True` reverses ties
        # (Python's sort is stable either way: verified directly,
        # `sorted([...], key=..., reverse=True)` preserves tied elements'
        # original relative order, same as negating the key and sorting
        # ascending). The real divergence is upstream of the sort
        # entirely: the accumulation loop iterates
        # `self._subject_kind_index.items()` -- grouped by when each
        # distinct (subject, kind, object, reaction) scope key was FIRST
        # created, with each key's own id-list in its own append order --
        # not a flat merge in true global add() order. INTERLEAVED
        # insertion into different scopes (add scope A, then scope B,
        # then scope A again) exposes this: global insertion order is
        # [a1, b, a2], but grouped-by-key accumulation visits key A's
        # full list [a1, a2] before key B's [b], giving [a1, a2, b] --
        # wrong. The real fix's explicit `self._insertion_seq` tiebreak
        # (not reliant on any sort's tie-stability) restores true global
        # order regardless of scope-key grouping.
        store = ParameterStore()
        id_a1 = store.add(make_kd(object_id="BIOSIM:PROTEIN:000002", value=1.0))
        id_b = store.add(make_kd(object_id="BIOSIM:METABOLITE:000099", value=2.0))
        id_a2 = store.add(make_kd(object_id="BIOSIM:PROTEIN:000002", value=3.0))
        all_params = store.get_all_for_subject(INSULIN_RECEPTOR_ID)
        ids = [p.parameter_id for p in all_params[ParameterKind.KD.value]]
        assert ids == [id_a1, id_b, id_a2]

    def test_get_all_for_subject_orders_by_confidence_then_organism_then_system(self):
        # Non-tie ranking must still hold once ties break on insertion
        # order -- adding the seq tiebreak must not disturb rank ordering.
        store = ParameterStore()
        store.add(make_kd(object_id="BIOSIM:PROTEIN:000002", value=1.0, tier=ConfidenceTier.BRONZE))
        store.add(make_kd(object_id="BIOSIM:METABOLITE:000099", value=2.0, tier=ConfidenceTier.GOLD))
        all_params = store.get_all_for_subject(INSULIN_RECEPTOR_ID)
        results = all_params[ParameterKind.KD.value]
        assert results[0].confidence_tier == ConfidenceTier.GOLD
        assert results[1].confidence_tier == ConfidenceTier.BRONZE

    def test_get_all_for_subject_include_superseded_filters_across_scopes(self):
        store = ParameterStore()
        old_id = store.add(make_kd(object_id="BIOSIM:PROTEIN:000002", value=1.0))
        store.supersede(old_id, make_kd(object_id="BIOSIM:METABOLITE:000099", value=2.0))
        default_view = store.get_all_for_subject(INSULIN_RECEPTOR_ID)
        assert len(default_view[ParameterKind.KD.value]) == 1
        full_view = store.get_all_for_subject(INSULIN_RECEPTOR_ID, include_superseded=True)
        assert len(full_view[ParameterKind.KD.value]) == 2

    def test_get_all_for_subject_groups_by_kind_across_scopes(self):
        store = ParameterStore()
        store.add(make_km(value=45.0))
        store.add(make_kd(object_id="BIOSIM:PROTEIN:000002", value=1.0))
        store.add(make_kd(object_id="BIOSIM:METABOLITE:000099", value=2.0))
        all_params = store.get_all_for_subject(INSULIN_RECEPTOR_ID)
        assert len(all_params[ParameterKind.KM.value]) == 1
        assert len(all_params[ParameterKind.KD.value]) == 2


# ---------------------------------------------------------------------------
# Versioning / supersession
# ---------------------------------------------------------------------------

class TestVersioning:
    def test_supersede_marks_old_as_superseded(self):
        store = ParameterStore()
        old = make_km(value=45.0)
        old_id = store.add(old)
        new = make_km(value=48.0, tier=ConfidenceTier.GOLD)
        store.supersede(old_id, new)
        assert store.get(old_id).is_superseded

    def test_superseded_excluded_from_query_by_default(self):
        store = ParameterStore()
        old = make_km(value=45.0)
        old_id = store.add(old)
        new = make_km(value=48.0, tier=ConfidenceTier.GOLD)
        store.supersede(old_id, new)
        results = store.query(INSULIN_RECEPTOR_ID, ParameterKind.KM)
        values = [r.value.value for r in results]
        assert 45.0 not in values
        assert 48.0 in values

    def test_superseded_included_when_requested(self):
        store = ParameterStore()
        old = make_km(value=45.0)
        old_id = store.add(old)
        new = make_km(value=48.0)
        store.supersede(old_id, new)
        results = store.query(
            INSULIN_RECEPTOR_ID, ParameterKind.KM, include_superseded=True
        )
        values = [r.value.value for r in results]
        assert 45.0 in values
        assert 48.0 in values

    def test_supersede_unknown_raises(self):
        store = ParameterStore()
        with pytest.raises(KeyError):
            store.supersede("nonexistent-id", make_km())


# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------

class TestPersistence:
    def test_save_and_load_round_trip(self, tmp_path):
        store = ParameterStore()
        store.add(make_km(value=45.0, tier=ConfidenceTier.GOLD))
        store.add(make_vmax(value=10.0))
        path = tmp_path / "params.json"
        store.save(path)

        loaded = ParameterStore.load(path)
        best_km = loaded.best(INSULIN_RECEPTOR_ID, ParameterKind.KM)
        assert best_km is not None
        assert best_km.value.value == 45.0
        assert best_km.confidence_tier == ConfidenceTier.GOLD

    def test_conflicts_persist(self, tmp_path):
        store = ParameterStore()
        store.add(make_km(value=10.0))
        store.add(make_km(value=100.0))
        path = tmp_path / "params.json"
        store.save(path)
        loaded = ParameterStore.load(path)
        assert len(loaded.unresolved_conflicts()) == 1

    def test_superseded_status_persists(self, tmp_path):
        store = ParameterStore()
        old = make_km(value=45.0)
        old_id = store.add(old)
        store.supersede(old_id, make_km(value=48.0))
        path = tmp_path / "params.json"
        store.save(path)
        loaded = ParameterStore.load(path)
        assert loaded.get(old_id).is_superseded

    def test_created_at_survives_save_load_round_trip(self, tmp_path):
        # Regression: _serialize_parameter() writes created_at, but
        # _deserialize_parameter() never passed it to the Parameter
        # constructor, so load() silently regenerated it via
        # default_factory=datetime.now(timezone.utc) -- destroying real
        # creation-time provenance every time a store is saved and
        # reloaded. Real consequence for Track 3: parameter.created_at is
        # timestamptz NOT NULL in the append-only Postgres schema, so a
        # migration backfill reading via load() would have permanently
        # stamped every row with the migration's wall clock instead of
        # the parameter's real creation time, in a table where UPDATE is
        # revoked and trigger-blocked -- unfixable after the fact.
        store = ParameterStore()
        p_id = store.add(make_km(value=45.0))
        original_created_at = store.get(p_id).created_at
        path = tmp_path / "params.json"
        store.save(path)
        loaded = ParameterStore.load(path)
        assert loaded.get(p_id).created_at == original_created_at

    def test_insertion_seq_persists_and_survives_reload(self, tmp_path):
        # Found during security review of the created_at fix above: seq
        # was re-derived purely from the parameters array's JSON order on
        # load(), not persisted -- correct as long as the array is never
        # reordered, but silently wrong (with no error) the moment it
        # is (a hand-edited snapshot, a merge of two partial exports).
        # Persist seq explicitly (same .get()-fallback pattern as
        # created_at/extracted_at) so a reloaded store's tie order is
        # robust rather than merely "correct by convention." Verified via
        # the same interleaved-scope scenario that discriminates the
        # get_all_for_subject() fix itself, round-tripped through save/load.
        store = ParameterStore()
        id_a1 = store.add(make_kd(object_id="BIOSIM:PROTEIN:000002", value=1.0))
        id_b = store.add(make_kd(object_id="BIOSIM:METABOLITE:000099", value=2.0))
        id_a2 = store.add(make_kd(object_id="BIOSIM:PROTEIN:000002", value=3.0))
        path = tmp_path / "params.json"
        store.save(path)
        raw = json.loads(path.read_text())
        assert all("seq" in p for p in raw["parameters"])
        loaded = ParameterStore.load(path)
        all_params = loaded.get_all_for_subject(INSULIN_RECEPTOR_ID)
        ids = [p.parameter_id for p in all_params[ParameterKind.KD.value]]
        assert ids == [id_a1, id_b, id_a2]

    def test_query_and_best_use_seq_not_array_order_after_a_reordered_snapshot(self, tmp_path):
        # Regression (MASTER_PLAN 3.1): query()/best() are the primary
        # lookup path the ODE engine actually uses (best()'s own
        # docstring), but only get_all_for_subject() had this tiebreak
        # fixed -- query()/best() relied on the unstated assumption that
        # a single index key's own id-list is "already a subsequence of
        # insertion order," true for the in-memory add()-only path but
        # false the moment a snapshot's "parameters" array isn't in true
        # seq order (a hand-edited snapshot, a merge of two partial
        # exports -- load()'s own docstring already anticipates this).
        # Reproduced directly here: two tied-priority parameters added
        # A-then-B, saved, the JSON array reversed while each
        # parameter's own persisted "seq" is left untouched, reloaded --
        # query()/best() must still return A first (seq order), matching
        # get_all_for_subject(), not B (array order).
        store = ParameterStore()
        id_a = store.add(make_km(value=1.0))
        id_b = store.add(make_km(value=9.0))
        path = tmp_path / "params.json"
        store.save(path)
        data = json.loads(path.read_text())
        data["parameters"].reverse()
        path.write_text(json.dumps(data))

        loaded = ParameterStore.load(path)
        queried_ids = [p.parameter_id for p in loaded.query(INSULIN_RECEPTOR_ID, ParameterKind.KM)]
        assert queried_ids == [id_a, id_b]
        assert loaded.best(INSULIN_RECEPTOR_ID, ParameterKind.KM).parameter_id == id_a
        # Confirms this now genuinely agrees with the already-fixed
        # get_all_for_subject(), not just independently "happens to work".
        all_params = loaded.get_all_for_subject(INSULIN_RECEPTOR_ID)
        assert [p.parameter_id for p in all_params[ParameterKind.KM.value]] == [id_a, id_b]

    def test_conflict_detected_at_survives_save_load_round_trip(self, tmp_path):
        # Same default_factory=datetime.now(timezone.utc)-regenerated-on
        # -load bug as Parameter.created_at/Provenance.extracted_at, for
        # ConflictRecord.detected_at -- load()'s conflict reconstruction
        # only overwrote conflict_id/resolved/resolution_notes from the
        # JSON, never detected_at.
        store = ParameterStore()
        store.add(make_km(value=10.0))
        store.add(make_km(value=100.0))
        conflict = store.unresolved_conflicts()[0]
        original_detected_at = conflict.detected_at
        path = tmp_path / "params.json"
        store.save(path)
        loaded = ParameterStore.load(path)
        loaded_conflict = loaded.unresolved_conflicts()[0]
        assert loaded_conflict.detected_at == original_detected_at

    def test_extracted_at_survives_save_load_round_trip(self, tmp_path):
        # Same bug, same fix, a second field: Provenance.extracted_at has
        # the identical default_factory=datetime.now(timezone.utc) gap --
        # _deserialize_parameter() builds a fresh Provenance() without
        # passing extracted_at through either. Found while fixing
        # created_at above; not previously flagged.
        store = ParameterStore()
        p_id = store.add(make_km(value=45.0))
        original_extracted_at = store.get(p_id).provenance.extracted_at
        path = tmp_path / "params.json"
        store.save(path)
        loaded = ParameterStore.load(path)
        assert loaded.get(p_id).provenance.extracted_at == original_extracted_at


# ---------------------------------------------------------------------------
# Thread safety
# ---------------------------------------------------------------------------

class TestThreadSafety:
    def test_concurrent_reads(self):
        store = ParameterStore()
        store.add(make_km(value=45.0))
        errors = []

        def read():
            try:
                for _ in range(100):
                    store.best(INSULIN_RECEPTOR_ID, ParameterKind.KM)
            except Exception as e:
                errors.append(e)

        threads = [threading.Thread(target=read) for _ in range(20)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        assert errors == []

    def test_concurrent_writes(self):
        store = ParameterStore()
        errors = []

        def write(i):
            try:
                for j in range(20):
                    store.add(make_km(value=float(i * 100 + j), subject=f"BIOSIM:PROTEIN:{i:06d}"))
            except Exception as e:
                errors.append(e)

        threads = [threading.Thread(target=write, args=(i,)) for i in range(10)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        assert errors == []
        assert store.stats()["total_parameters"] == 200


# ---------------------------------------------------------------------------
# Stats
# ---------------------------------------------------------------------------

class TestStats:
    def test_stats_counts_parameters(self):
        store = ParameterStore()
        store.add(make_km(tier=ConfidenceTier.GOLD))
        store.add(make_vmax(tier=ConfidenceTier.SILVER))
        stats = store.stats()
        assert stats["total_parameters"] == 2
        assert stats["active_parameters"] == 2

    def test_stats_counts_superseded(self):
        store = ParameterStore()
        old = make_km(value=45.0)
        old_id = store.add(old)
        store.supersede(old_id, make_km(value=48.0))
        stats = store.stats()
        assert stats["superseded_parameters"] == 1
        assert stats["active_parameters"] == 1

    def test_stats_counts_by_tier(self):
        store = ParameterStore()
        store.add(make_km(tier=ConfidenceTier.GOLD))
        store.add(make_vmax(tier=ConfidenceTier.SILVER))
        stats = store.stats()
        assert stats["by_confidence_tier"]["gold"] == 1
        assert stats["by_confidence_tier"]["silver"] == 1


class TestVerificationStatus:
    """Provenance records whether the stored value was checked against the
    cited text, and how (2026-11-15 provenance audit)."""

    def test_default_status_is_unchecked(self):
        assert make_km().provenance.verification_status == "UNCHECKED"
        assert make_km().provenance.verification_quote == ""

    def test_builder_records_status_and_quote(self):
        from biosim_parameter_db import VerificationStatus
        p = (
            ParameterBuilder()
            .for_subject(GLUCOSE_ID)
            .kind(ParameterKind.KM)
            .value(5.0, unit="mM")
            .from_pubmed("12345678")
            .verification(VerificationStatus.VERIFIED, "Km was 5.0 mM")
            .build()
        )
        assert p.provenance.verification_status == "VERIFIED"
        assert p.provenance.verification_quote == "Km was 5.0 mM"

    def test_unknown_status_is_rejected(self):
        with pytest.raises(ValueError):
            Provenance(pubmed_id="1", verification_status="LOOKS_FINE")

    def test_status_survives_save_and_load(self, tmp_path):
        from biosim_parameter_db import VerificationStatus
        store = ParameterStore()
        p = (
            ParameterBuilder()
            .for_subject(GLUCOSE_ID)
            .kind(ParameterKind.KM)
            .value(5.0, unit="mM")
            .from_pubmed("12345678")
            .verification(VerificationStatus.DERIVED, "10.13 mM")
            .build()
        )
        pid = store.add(p)
        path = tmp_path / "store.json"
        store.save(path)
        loaded = ParameterStore.load(path).get(pid).provenance
        assert loaded.verification_status == "DERIVED"
        assert loaded.verification_quote == "10.13 mM"
