# SAVE TO: packages\context-registry\tests\test_context_registry.py
"""
Test suite for biosim_context_registry.

New package (2026-08-18, a project architectural decision): a
deliberately separate, deliberately simpler registry for non-molecular
CONTEXT entities (disease states, physiological axes) than
biosim_id_resolver.IDRegistry -- see types.py's own module docstring
for the full reasoning.
"""
import dataclasses

import pytest

from biosim_context_registry import (
    ContextEntityType,
    ContextID,
    ContextRecord,
    ContextRegistry,
    ContextResolutionError,
    normalize_name,
)


class TestNormalizeName:
    def test_lowercases(self):
        assert normalize_name("Chronic Kidney Disease") == "chronic_kidney_disease"

    def test_strips_whitespace(self):
        assert normalize_name("  sleep restriction  ") == "sleep_restriction"

    def test_already_normalized_is_idempotent(self):
        assert normalize_name("chronic_kidney_disease") == "chronic_kidney_disease"

    def test_collapses_repeated_whitespace(self):
        # "chronic  kidney disease" (double space) must normalize to the
        # SAME string as "chronic kidney disease" (single space) -- before
        # this fix, only single spaces were replaced, so the two would
        # normalize to different strings and evade register()'s duplicate
        # check entirely.
        assert normalize_name("chronic  kidney disease") == normalize_name("chronic kidney disease")

    def test_rejects_colon(self):
        # ':' is ContextID's own string-format delimiter
        # (CONTEXT:<type>:<name>) -- a name containing one would produce
        # an ID that from_string() could never parse back.
        with pytest.raises(ValueError, match="a-z0-9_"):
            normalize_name("foo:bar")

    def test_rejects_name_exceeding_max_length(self):
        with pytest.raises(ValueError, match="exceeds"):
            normalize_name("a" * 201)

    def test_accepts_name_at_max_length(self):
        assert normalize_name("a" * 200) == "a" * 200


class TestContextID:
    def test_str_format(self):
        cid = ContextID(entity_type=ContextEntityType.DISEASE_STATE, name="chronic_kidney_disease")
        assert str(cid) == "CONTEXT:DISEASE_STATE:chronic_kidney_disease"

    def test_empty_name_raises(self):
        with pytest.raises(ValueError):
            ContextID(entity_type=ContextEntityType.DISEASE_STATE, name="")

    def test_unnormalized_name_raises(self):
        # Direct construction requires an already-normalized name --
        # callers should go through ContextRegistry.register()/
        # resolve_by_name(), which normalize automatically.
        with pytest.raises(ValueError):
            ContextID(entity_type=ContextEntityType.DISEASE_STATE, name="Chronic Kidney Disease")

    def test_from_string_round_trips(self):
        cid = ContextID(entity_type=ContextEntityType.PHYSIOLOGICAL_AXIS, name="sleep_restriction")
        assert ContextID.from_string(str(cid)) == cid

    def test_from_string_rejects_bad_prefix(self):
        with pytest.raises(ValueError):
            ContextID.from_string("BIOSIM:DISEASE_STATE:chronic_kidney_disease")

    def test_from_string_rejects_unknown_entity_type(self):
        with pytest.raises(ValueError):
            ContextID.from_string("CONTEXT:NOT_A_REAL_TYPE:foo")

    def test_from_string_rejects_wrong_part_count(self):
        with pytest.raises(ValueError):
            ContextID.from_string("CONTEXT:DISEASE_STATE")

    def test_construction_rejects_colon_in_name(self):
        # Direct construction goes through normalize_name()'s own
        # round-trip check in __post_init__, so a colon-containing name
        # is rejected here too, not just via ContextRegistry.register().
        with pytest.raises(ValueError):
            ContextID(entity_type=ContextEntityType.DISEASE_STATE, name="foo:bar")


class TestContextRecord:
    def test_is_frozen(self):
        record = ContextRecord(
            context_id=ContextID(ContextEntityType.DISEASE_STATE, "x"),
            display_name="x",
        )
        with pytest.raises(dataclasses.FrozenInstanceError):
            record.description = "mutated"


class TestContextRegistry:
    def test_register_returns_context_id(self):
        r = ContextRegistry()
        cid = r.register("chronic_kidney_disease", ContextEntityType.DISEASE_STATE)
        assert cid == ContextID(ContextEntityType.DISEASE_STATE, "chronic_kidney_disease")

    def test_register_normalizes_the_name(self):
        r = ContextRegistry()
        cid = r.register("Chronic Kidney Disease", ContextEntityType.DISEASE_STATE)
        assert cid.name == "chronic_kidney_disease"

    def test_duplicate_registration_raises(self):
        r = ContextRegistry()
        r.register("chronic_kidney_disease", ContextEntityType.DISEASE_STATE)
        with pytest.raises(ValueError, match="already registered"):
            r.register("chronic kidney disease", ContextEntityType.DISEASE_STATE)

    def test_duplicate_registration_error_does_not_promise_nonexistent_rename_method(self):
        # This registry has no deregister()/rename() method -- the error
        # message must not tell the caller to "remove the old
        # registration" (a capability that doesn't exist), only point at
        # constructing a fresh registry or editing the seed file.
        r = ContextRegistry()
        r.register("chronic_kidney_disease", ContextEntityType.DISEASE_STATE)
        with pytest.raises(ValueError) as exc_info:
            r.register("chronic_kidney_disease", ContextEntityType.DISEASE_STATE)
        assert "remove the old registration" not in str(exc_info.value)

    def test_register_rejects_name_with_colon(self):
        r = ContextRegistry()
        with pytest.raises(ValueError):
            r.register("foo:bar", ContextEntityType.DISEASE_STATE)

    def test_same_name_different_entity_type_does_not_collide(self):
        # entity_type is part of the identity -- a disease state and a
        # physiological axis coincidentally sharing a name (not a real
        # case today, but a real design property to lock in) must not
        # collide.
        r = ContextRegistry()
        r.register("x", ContextEntityType.DISEASE_STATE)
        r.register("x", ContextEntityType.PHYSIOLOGICAL_AXIS)   # must not raise

    def test_resolve_by_name(self):
        r = ContextRegistry()
        cid = r.register("sleep_restriction", ContextEntityType.PHYSIOLOGICAL_AXIS)
        resolved = r.resolve_by_name("sleep_restriction", ContextEntityType.PHYSIOLOGICAL_AXIS)
        assert resolved == cid

    def test_resolve_by_name_is_case_and_whitespace_insensitive(self):
        r = ContextRegistry()
        r.register("sleep_restriction", ContextEntityType.PHYSIOLOGICAL_AXIS)
        resolved = r.resolve_by_name(" Sleep Restriction ", ContextEntityType.PHYSIOLOGICAL_AXIS)
        assert resolved.name == "sleep_restriction"

    def test_resolve_unregistered_name_raises(self):
        r = ContextRegistry()
        with pytest.raises(ContextResolutionError):
            r.resolve_by_name("not_a_real_thing", ContextEntityType.DISEASE_STATE)

    def test_get_record(self):
        r = ContextRegistry()
        cid = r.register(
            "chronic_kidney_disease", ContextEntityType.DISEASE_STATE,
            description="Moderate-severe non-diabetic CKD.",
        )
        record = r.get_record(cid)
        assert record.display_name == "chronic_kidney_disease"
        assert record.description == "Moderate-severe non-diabetic CKD."

    def test_get_record_unregistered_raises(self):
        r = ContextRegistry()
        cid = ContextID(entity_type=ContextEntityType.DISEASE_STATE, name="nonexistent")
        with pytest.raises(ContextResolutionError):
            r.get_record(cid)

    def test_is_registered(self):
        r = ContextRegistry()
        assert not r.is_registered("chronic_kidney_disease", ContextEntityType.DISEASE_STATE)
        r.register("chronic_kidney_disease", ContextEntityType.DISEASE_STATE)
        assert r.is_registered("chronic_kidney_disease", ContextEntityType.DISEASE_STATE)

    def test_stats(self):
        r = ContextRegistry()
        r.register("chronic_kidney_disease", ContextEntityType.DISEASE_STATE)
        r.register("sleep_restriction", ContextEntityType.PHYSIOLOGICAL_AXIS)
        stats = r.stats()
        assert stats["total_entities"] == 2
        assert stats["by_entity_type"] == {"DISEASE_STATE": 1, "PHYSIOLOGICAL_AXIS": 1}
