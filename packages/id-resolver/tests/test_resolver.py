"""
Test suite for biosim_id_resolver.

Coverage targets:
    - Happy path resolution for all ontology sources
    - Resolution priority ordering
    - Ambiguous resolution detection (data integrity violation)
    - Deprecated canonical ID handling
    - Cross-reference translation
    - Name-based search
    - Serialization round-trip (to_dict / from_dict)
    - File-based loading
    - Thread safety under concurrent reads and writes
    - Input validation and edge cases
    - Registry statistics
"""

import json
import threading
import tempfile
from pathlib import Path

import pytest

from biosim_id_resolver import (
    AmbiguousResolutionError,
    CanonicalID,
    EntityType,
    ExternalID,
    IDRegistry,
    OntologySource,
    ResolutionError,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def empty_registry() -> IDRegistry:
    return IDRegistry()


@pytest.fixture
def seeded_registry() -> IDRegistry:
    """
    Registry pre-seeded with a small set of real biological entities
    from the glucose/insulin axis — the v1 pathway cluster.
    """
    r = IDRegistry()

    # D-Glucose — cross-referenced across ChEBI, HMDB, KEGG
    r.register(
        name="D-Glucose",
        entity_type=EntityType.METABOLITE,
        description="Primary blood sugar; substrate for glycolysis",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:17234"],
            OntologySource.HMDB: ["HMDB0000122"],
            OntologySource.KEGG: ["C00031"],
        },
        force_numeric_id=1,
    )

    # Insulin receptor — UniProt primary
    r.register(
        name="Insulin receptor",
        entity_type=EntityType.PROTEIN,
        description="Receptor tyrosine kinase; mediates cellular insulin response",
        cross_references={
            OntologySource.UNIPROT: ["P06213"],
            OntologySource.KEGG: ["K04527"],
        },
        force_numeric_id=1,
    )

    # IGF-1 — UniProt primary
    r.register(
        name="Insulin-like growth factor 1",
        entity_type=EntityType.PROTEIN,
        description="Mediates growth hormone signalling; structural homology to insulin",
        cross_references={
            OntologySource.UNIPROT: ["P05019"],
            OntologySource.HMDB: ["HMDB0005000"],
        },
        force_numeric_id=2,
    )

    # Glycolysis pathway — Reactome primary
    r.register(
        name="Glycolysis",
        entity_type=EntityType.PATHWAY,
        description="Cytoplasmic conversion of glucose to pyruvate",
        cross_references={
            OntologySource.REACTOME: ["R-HSA-70171"],
            OntologySource.KEGG: ["hsa00010"],
        },
        force_numeric_id=1,
    )

    return r


# ---------------------------------------------------------------------------
# CanonicalID and ExternalID type tests
# ---------------------------------------------------------------------------

class TestCanonicalID:
    def test_string_representation(self):
        cid = CanonicalID(EntityType.PROTEIN, 42)
        assert str(cid) == "BIOSIM:PROTEIN:000042"

    def test_zero_padded(self):
        cid = CanonicalID(EntityType.METABOLITE, 1)
        assert str(cid) == "BIOSIM:METABOLITE:000001"

    def test_large_id(self):
        cid = CanonicalID(EntityType.GENE, 999999)
        assert str(cid) == "BIOSIM:GENE:999999"

    def test_from_string_roundtrip(self):
        original = CanonicalID(EntityType.PATHWAY, 123)
        parsed = CanonicalID.from_string(str(original))
        assert parsed == original

    def test_from_string_invalid_prefix(self):
        with pytest.raises(ValueError, match="Invalid canonical ID format"):
            CanonicalID.from_string("NOTBIOSIM:PROTEIN:000001")

    def test_from_string_unknown_type(self):
        with pytest.raises(ValueError, match="Unknown entity type"):
            CanonicalID.from_string("BIOSIM:WIDGET:000001")

    def test_from_string_non_integer(self):
        with pytest.raises(ValueError, match="Non-integer numeric ID"):
            CanonicalID.from_string("BIOSIM:PROTEIN:abc")

    def test_negative_id_rejected(self):
        with pytest.raises(ValueError):
            CanonicalID(EntityType.PROTEIN, -1)

    def test_frozen(self):
        cid = CanonicalID(EntityType.PROTEIN, 1)
        with pytest.raises(Exception):
            cid.numeric_id = 2  # type: ignore


class TestExternalID:
    def test_string_representation(self):
        eid = ExternalID(OntologySource.UNIPROT, "P06213")
        assert str(eid) == "uniprot:P06213"

    def test_accession_normalized_uppercase(self):
        eid = ExternalID(OntologySource.CHEBI, "chebi:17234")
        assert eid.accession == "CHEBI:17234"

    def test_accession_whitespace_stripped(self):
        eid = ExternalID(OntologySource.HMDB, "  HMDB0000122  ")
        assert eid.accession == "HMDB0000122"

    def test_empty_accession_rejected(self):
        with pytest.raises(ValueError):
            ExternalID(OntologySource.UNIPROT, "")

    def test_whitespace_only_accession_rejected(self):
        with pytest.raises(ValueError):
            ExternalID(OntologySource.UNIPROT, "   ")

    def test_frozen(self):
        eid = ExternalID(OntologySource.UNIPROT, "P06213")
        with pytest.raises(Exception):
            eid.accession = "OTHER"  # type: ignore


# ---------------------------------------------------------------------------
# Registration tests
# ---------------------------------------------------------------------------

class TestRegistration:
    def test_register_returns_canonical_id(self, empty_registry):
        cid = empty_registry.register(
            name="D-Glucose",
            entity_type=EntityType.METABOLITE,
        )
        assert isinstance(cid, CanonicalID)
        assert cid.entity_type == EntityType.METABOLITE

    def test_auto_increment_ids(self, empty_registry):
        cid1 = empty_registry.register(name="Entity A", entity_type=EntityType.PROTEIN)
        cid2 = empty_registry.register(name="Entity B", entity_type=EntityType.PROTEIN)
        assert cid2.numeric_id == cid1.numeric_id + 1

    def test_id_spaces_independent_per_type(self, empty_registry):
        p1 = empty_registry.register(name="Protein 1", entity_type=EntityType.PROTEIN)
        m1 = empty_registry.register(name="Metabolite 1", entity_type=EntityType.METABOLITE)
        assert p1.numeric_id == m1.numeric_id  # both start at 1

    def test_cross_reference_collision_raises(self, empty_registry):
        empty_registry.register(
            name="D-Glucose",
            entity_type=EntityType.METABOLITE,
            cross_references={OntologySource.CHEBI: ["CHEBI:17234"]},
        )
        with pytest.raises(ValueError, match="Cross-reference collision"):
            empty_registry.register(
                name="Beta-D-Glucose",
                entity_type=EntityType.METABOLITE,
                cross_references={OntologySource.CHEBI: ["CHEBI:17234"]},
            )

    def test_force_numeric_id(self, empty_registry):
        cid = empty_registry.register(
            name="Test entity",
            entity_type=EntityType.GENE,
            force_numeric_id=9999,
        )
        assert cid.numeric_id == 9999


# ---------------------------------------------------------------------------
# Resolution tests
# ---------------------------------------------------------------------------

class TestResolution:
    def test_resolve_by_chebi(self, seeded_registry):
        eid = ExternalID(OntologySource.CHEBI, "CHEBI:17234")
        cid = seeded_registry.resolve(eid)
        assert cid.entity_type == EntityType.METABOLITE
        assert str(cid) == "BIOSIM:METABOLITE:000001"

    def test_resolve_by_hmdb(self, seeded_registry):
        eid = ExternalID(OntologySource.HMDB, "HMDB0000122")
        cid = seeded_registry.resolve(eid)
        assert str(cid) == "BIOSIM:METABOLITE:000001"

    def test_resolve_by_kegg(self, seeded_registry):
        eid = ExternalID(OntologySource.KEGG, "C00031")
        cid = seeded_registry.resolve(eid)
        assert str(cid) == "BIOSIM:METABOLITE:000001"

    def test_resolve_by_uniprot(self, seeded_registry):
        eid = ExternalID(OntologySource.UNIPROT, "P06213")
        cid = seeded_registry.resolve(eid)
        assert cid.entity_type == EntityType.PROTEIN

    def test_resolve_by_reactome(self, seeded_registry):
        eid = ExternalID(OntologySource.REACTOME, "R-HSA-70171")
        cid = seeded_registry.resolve(eid)
        assert cid.entity_type == EntityType.PATHWAY

    def test_resolve_unknown_raises(self, seeded_registry):
        eid = ExternalID(OntologySource.UNIPROT, "UNKNOWN999")
        with pytest.raises(ResolutionError, match="Cannot resolve"):
            seeded_registry.resolve(eid)

    def test_resolve_case_insensitive(self, seeded_registry):
        eid_upper = ExternalID(OntologySource.UNIPROT, "P06213")
        eid_lower = ExternalID(OntologySource.UNIPROT, "p06213")
        assert seeded_registry.resolve(eid_upper) == seeded_registry.resolve(eid_lower)

    def test_resolve_any_finds_across_sources(self, seeded_registry):
        results = seeded_registry.resolve_any("HMDB0000122")
        assert len(results) == 1
        assert results[0].entity_type == EntityType.METABOLITE

    def test_resolve_any_returns_empty_for_unknown(self, seeded_registry):
        results = seeded_registry.resolve_any("COMPLETELYFAKEID")
        assert results == []

    def test_ambiguous_resolution_raises(self, empty_registry):
        """
        Two entities registered with the same cross-reference should be caught
        at registration time, preventing ambiguous resolution at query time.
        """
        empty_registry.register(
            name="Entity A",
            entity_type=EntityType.METABOLITE,
            cross_references={OntologySource.HMDB: ["HMDB9999999"]},
        )
        with pytest.raises(ValueError, match="Cross-reference collision"):
            empty_registry.register(
                name="Entity B",
                entity_type=EntityType.METABOLITE,
                cross_references={OntologySource.HMDB: ["HMDB9999999"]},
            )


# ---------------------------------------------------------------------------
# Record retrieval tests
# ---------------------------------------------------------------------------

class TestRecordRetrieval:
    def test_get_record_name(self, seeded_registry):
        cid = seeded_registry.resolve(ExternalID(OntologySource.CHEBI, "CHEBI:17234"))
        record = seeded_registry.get_record(cid)
        assert record.name == "D-Glucose"

    def test_get_record_description(self, seeded_registry):
        cid = seeded_registry.resolve(ExternalID(OntologySource.CHEBI, "CHEBI:17234"))
        record = seeded_registry.get_record(cid)
        assert "glycolysis" in record.description.lower()

    def test_get_record_unknown_raises(self, empty_registry):
        fake_cid = CanonicalID(EntityType.PROTEIN, 99999)
        with pytest.raises(ResolutionError):
            empty_registry.get_record(fake_cid)

    def test_cross_reference_translation(self, seeded_registry):
        cid = seeded_registry.resolve(ExternalID(OntologySource.CHEBI, "CHEBI:17234"))
        kegg_id = seeded_registry.cross_reference(cid, OntologySource.KEGG)
        assert kegg_id == "C00031"

    def test_cross_reference_returns_none_if_missing(self, seeded_registry):
        cid = seeded_registry.resolve(ExternalID(OntologySource.CHEBI, "CHEBI:17234"))
        reactome_id = seeded_registry.cross_reference(cid, OntologySource.REACTOME)
        assert reactome_id is None

    def test_search_by_name(self, seeded_registry):
        results = seeded_registry.search_by_name("D-Glucose")
        assert len(results) == 1
        assert results[0].entity_type == EntityType.METABOLITE

    def test_search_by_name_case_insensitive(self, seeded_registry):
        results_lower = seeded_registry.search_by_name("d-glucose")
        results_upper = seeded_registry.search_by_name("D-Glucose")
        assert results_lower == results_upper

    def test_search_by_name_not_found(self, seeded_registry):
        results = seeded_registry.search_by_name("Unobtanium")
        assert results == []


# ---------------------------------------------------------------------------
# Deprecated ID tests
# ---------------------------------------------------------------------------

class TestDeprecation:
    def test_deprecated_id_raises_with_migration_hint(self, empty_registry):
        old_cid = CanonicalID(EntityType.PROTEIN, 999)
        new_cid = empty_registry.register(
            name="Migrated protein",
            entity_type=EntityType.PROTEIN,
            deprecated_canonical_ids=[old_cid],
        )
        with pytest.raises(ResolutionError, match="deprecated"):
            empty_registry.get_record(old_cid)

    def test_deprecated_id_error_contains_current_id(self, empty_registry):
        old_cid = CanonicalID(EntityType.PROTEIN, 888)
        new_cid = empty_registry.register(
            name="Migrated protein",
            entity_type=EntityType.PROTEIN,
            deprecated_canonical_ids=[old_cid],
        )
        try:
            empty_registry.get_record(old_cid)
        except ResolutionError as e:
            assert str(new_cid) in str(e)

    def test_current_id_still_resolves_after_deprecation(self, empty_registry):
        old_cid = CanonicalID(EntityType.PROTEIN, 777)
        new_cid = empty_registry.register(
            name="Migrated protein",
            entity_type=EntityType.PROTEIN,
            deprecated_canonical_ids=[old_cid],
        )
        record = empty_registry.get_record(new_cid)
        assert record.name == "Migrated protein"

    def test_merging_a_real_entity_with_shared_xref_does_not_raise(self, empty_registry):
        # Regression: deprecated_canonical_ids is the documented mechanism
        # for merging two entities that legitimately share a real-world
        # identifier -- but the collision check used to fire
        # unconditionally, with no exemption for accessions belonging to
        # an entity being deprecated by this same call. Since register()
        # is the only way to add data (no separate update()/merge()
        # method exists), this made the merge workflow entirely
        # non-functional. Every OTHER test in this class deprecates an ID
        # that was never actually registered, which is why this was never
        # caught -- this test registers the OLD entity first, for real,
        # with the SAME cross-reference the new merged entity will carry.
        old = empty_registry.register(
            name="D-Glucose (old)",
            entity_type=EntityType.METABOLITE,
            cross_references={OntologySource.HMDB: ["HMDB0000122"]},
        )
        merged = empty_registry.register(
            name="D-Glucose (merged)",
            entity_type=EntityType.METABOLITE,
            cross_references={OntologySource.HMDB: ["HMDB0000122"]},
            deprecated_canonical_ids=[old],
        )
        assert merged != old

    def test_merged_entity_resolves_unambiguously_not_a_second_candidate(self, empty_registry):
        # A second latent bug alongside the collision check: even once the
        # merge is allowed, nothing used to remove the deprecated entity's
        # OLD cross-reference index entry, so resolve() would then return
        # TWO candidates for the shared accession and start raising
        # AmbiguousResolutionError for what should now be an unambiguous,
        # correctly-merged lookup.
        old = empty_registry.register(
            name="D-Glucose (old)",
            entity_type=EntityType.METABOLITE,
            cross_references={OntologySource.HMDB: ["HMDB0000122"]},
        )
        merged = empty_registry.register(
            name="D-Glucose (merged)",
            entity_type=EntityType.METABOLITE,
            cross_references={OntologySource.HMDB: ["HMDB0000122"]},
            deprecated_canonical_ids=[old],
        )
        resolved = empty_registry.resolve(ExternalID(OntologySource.HMDB, "HMDB0000122"))
        assert resolved == merged

    def test_merge_inherits_cross_references_not_restated_in_the_call(self, empty_registry):
        # Regression: a merge used to only re-point cross-references
        # EXPLICITLY RESTATED in the new register() call's
        # cross_references argument -- any accession the old entity owned
        # under a source the caller didn't happen to restate was left
        # silently pointing at the now-stale, deprecated canonical ID
        # forever (resolve() doesn't consult the deprecated index at all,
        # only get_record() does, so this failed with no warning).
        old = empty_registry.register(
            name="D-Glucose (old)",
            entity_type=EntityType.METABOLITE,
            cross_references={
                OntologySource.HMDB: ["HMDB0000122"],
                OntologySource.KEGG: ["C00031"],
            },
        )
        merged = empty_registry.register(
            name="D-Glucose (merged)",
            entity_type=EntityType.METABOLITE,
            cross_references={OntologySource.HMDB: ["HMDB0000122"]},  # KEGG not restated
            deprecated_canonical_ids=[old],
        )
        resolved = empty_registry.resolve(ExternalID(OntologySource.KEGG, "C00031"))
        assert resolved == merged

    def test_merge_producing_two_accessions_for_one_source_keeps_the_new_ones_own_as_primary(self, empty_registry):
        # Regression/documentation lock for MASTER_PLAN Track 3.4 (the
        # Neo4j invariant-4 decision, infra/MIGRATION_PLAN.md question 8):
        # merging an entity into a new one that restates a DIFFERENT
        # accession from the SAME ontology source is a genuine, legitimate
        # scenario (e.g. two curated records for the same protein initially
        # assigned different UniProt accessions, later consolidated) --
        # not a data-integrity violation. The resulting active entity
        # legitimately carries BOTH accessions, and cross_reference()
        # must deterministically return the NEW entity's own restated
        # accession (list position 0) as "primary", not the inherited one,
        # regardless of alphabetical order between the two. This exact
        # contract is what the Neo4j schema's `IDENTIFIES.rank` property
        # and its rewritten invariant 4 depend on matching.
        old = empty_registry.register(
            name="Protein isoform (old)",
            entity_type=EntityType.PROTEIN,
            cross_references={OntologySource.UNIPROT: ["P00001"]},
        )
        merged = empty_registry.register(
            name="Protein isoform (consolidated)",
            entity_type=EntityType.PROTEIN,
            cross_references={OntologySource.UNIPROT: ["P99999"]},  # alphabetically AFTER P00001
            deprecated_canonical_ids=[old],
        )
        record = empty_registry.get_record(merged)
        assert record.cross_references[OntologySource.UNIPROT] == ["P99999", "P00001"]
        assert record.get_external_id(OntologySource.UNIPROT) == "P99999"
        assert empty_registry.cross_reference(merged, OntologySource.UNIPROT) == "P99999"
        # Both accessions must still resolve -- the merge doesn't make the
        # deprecated entity's own accession unresolvable.
        assert empty_registry.resolve(ExternalID(OntologySource.UNIPROT, "P00001")) == merged
        assert empty_registry.resolve(ExternalID(OntologySource.UNIPROT, "P99999")) == merged

    def test_merge_does_not_double_count_the_deprecated_entity_in_stats(self, empty_registry):
        # Regression: nothing removed the deprecated entity's own entry
        # from the internal canonical index, so stats() (used for
        # monitoring/health checks) silently double-counted every
        # merged-away entity as still a separate, live one.
        before = empty_registry.stats()["total_entities"]
        old = empty_registry.register(
            name="D-Glucose (old)",
            entity_type=EntityType.METABOLITE,
            cross_references={OntologySource.HMDB: ["HMDB0000122"]},
        )
        empty_registry.register(
            name="D-Glucose (merged)",
            entity_type=EntityType.METABOLITE,
            cross_references={OntologySource.HMDB: ["HMDB0000122"]},
            deprecated_canonical_ids=[old],
        )
        after = empty_registry.stats()["total_entities"]
        assert after == before + 1  # one net new entity, not two

    def test_genuine_third_party_collision_still_raises_during_a_merge(self, empty_registry):
        # The collision-check exemption must be scoped to the entities
        # actually being deprecated by THIS call -- an unrelated third
        # entity sharing the same accession is a real data-integrity
        # problem and must still raise.
        unrelated = empty_registry.register(
            name="Unrelated entity",
            entity_type=EntityType.METABOLITE,
            cross_references={OntologySource.HMDB: ["HMDB0000122"]},
        )
        old = empty_registry.register(
            name="D-Glucose (old)",
            entity_type=EntityType.METABOLITE,
            cross_references={OntologySource.CHEBI: ["CHEBI:1"]},
        )
        with pytest.raises(ValueError, match="collision"):
            empty_registry.register(
                name="D-Glucose (merged)",
                entity_type=EntityType.METABOLITE,
                cross_references={OntologySource.HMDB: ["HMDB0000122"]},
                deprecated_canonical_ids=[old],
            )

    def test_transitive_merge_chain_reports_the_final_live_id_not_a_dead_intermediate(self, empty_registry):
        # Regression: _deprecated_index[old] is set once, at the moment
        # `old` is deprecated, and never updated if the entity it points
        # at is LATER deprecated too (A->B->C). get_record(A) used to
        # report "current canonical ID is B" as final truth, even though
        # B was itself already deprecated into C by the time A is looked
        # up -- a caller trusting that single hint lands on another dead
        # id. Confirmed directly before fixing.
        a = CanonicalID(EntityType.PROTEIN, 501)
        b = empty_registry.register(
            name="Intermediate", entity_type=EntityType.PROTEIN,
            deprecated_canonical_ids=[a],
        )
        c = empty_registry.register(
            name="Final", entity_type=EntityType.PROTEIN,
            deprecated_canonical_ids=[b],
        )
        with pytest.raises(ResolutionError) as exc_info:
            empty_registry.get_record(a)
        assert str(c) in str(exc_info.value)
        assert str(b) not in str(exc_info.value)

    def test_search_by_old_name_still_resolves_after_merge(self, empty_registry):
        # Regression: _name_index was never touched when an entity was
        # deprecated -- search_by_name() for the OLD entity's own
        # original name still returned its now-dead canonical id, and
        # get_record() on that result raised ResolutionError, even
        # though the entity is still perfectly findable, just under a
        # new id. Confirmed directly before fixing (with the fix
        # reverted, this test's own assertion below fails with
        # ResolutionError instead of finding the merged record).
        old = empty_registry.register(name="Old Name", entity_type=EntityType.PROTEIN)
        merged = empty_registry.register(
            name="Merged Name", entity_type=EntityType.PROTEIN,
            deprecated_canonical_ids=[old],
        )
        found = empty_registry.search_by_name("Old Name")
        assert len(found) == 1
        assert found[0] == merged
        record = empty_registry.get_record(found[0])
        assert record.name == "Merged Name"

    def test_search_by_old_name_does_not_duplicate_if_names_collide(self, empty_registry):
        # If the new (merging) entity happens to share the old entity's
        # exact name, repointing the old name-index entry must not leave
        # two identical entries for the same canonical id.
        old = empty_registry.register(name="Shared Name", entity_type=EntityType.PROTEIN)
        merged = empty_registry.register(
            name="Shared Name", entity_type=EntityType.PROTEIN,
            deprecated_canonical_ids=[old],
        )
        found = empty_registry.search_by_name("Shared Name")
        assert found == [merged]

    def test_transitive_merge_chain_survives_serialization_roundtrip(self, empty_registry):
        # Regression: to_dict() only serialized each SURVIVING record's
        # own deprecated_canonical_ids (a direct, one-hop list captured
        # at that record's own registration) -- in a transitive chain
        # A->B->C, C's own record lists only [B]; A's deprecation (which
        # only ever lived on B's now-removed record) was nowhere in the
        # serialized form at all. Reloading from that dict made A
        # permanently unresolvable ("No record found"), even though it
        # had a real, resolvable migration path before the round-trip.
        a = CanonicalID(EntityType.PROTEIN, 601)
        b = empty_registry.register(
            name="Intermediate", entity_type=EntityType.PROTEIN,
            deprecated_canonical_ids=[a],
        )
        c = empty_registry.register(
            name="Final", entity_type=EntityType.PROTEIN,
            deprecated_canonical_ids=[b],
        )
        reloaded = IDRegistry.from_dict(empty_registry.to_dict())
        with pytest.raises(ResolutionError) as exc_info:
            reloaded.get_record(a)
        assert str(c) in str(exc_info.value)
        assert reloaded.get_record(c).name == "Final"

    def test_force_numeric_id_self_loop_raises(self, empty_registry):
        # Code-reviewer finding: force_numeric_id was never validated,
        # even though both real deserialization entry points
        # (from_dict()/from_mapping_files()) call register() with it
        # unconditionally for every loaded record -- a corrupted/hand-
        # edited snapshot could reach this, not just adversarial test
        # code. A self-loop (deprecating X while also forcing this new
        # registration's own id to equal X) is nonsensical and must be
        # rejected at registration, not silently create a cycle for
        # _resolve_deprecation_chain()'s guard to catch later.
        a = empty_registry.register(name="A", entity_type=EntityType.PROTEIN)
        with pytest.raises(ValueError, match="collision"):
            empty_registry.register(
                name="Self-loop",
                entity_type=EntityType.PROTEIN,
                deprecated_canonical_ids=[a],
                force_numeric_id=a.numeric_id,
            )

    def test_force_numeric_id_colliding_with_existing_deprecated_id_raises(self, empty_registry):
        # A different collision from the self-loop above: force_numeric_id
        # reconstructs a canonical_str that is ALREADY a stale
        # _deprecated_index key from an earlier, unrelated merge --
        # get_record() checks _deprecated_index before _canonical_index,
        # so a naive load would make a live, just-registered entity
        # incorrectly report itself as deprecated.
        old = CanonicalID(EntityType.PROTEIN, 701)
        empty_registry.register(
            name="Already merged", entity_type=EntityType.PROTEIN,
            deprecated_canonical_ids=[old],
        )
        with pytest.raises(ValueError, match="collision"):
            empty_registry.register(
                name="Colliding new entity",
                entity_type=EntityType.PROTEIN,
                force_numeric_id=old.numeric_id,
            )

    def test_normal_from_dict_roundtrip_does_not_false_positive_on_the_new_validation(self, empty_registry):
        # Confirms the new force_numeric_id validation doesn't reject
        # ORDINARY, non-adversarial reconstruction -- from_dict()'s own
        # real, everyday usage (not just the transitive-chain edge case
        # above) must keep working.
        empty_registry.register(
            name="Ordinary entity", entity_type=EntityType.PROTEIN,
            cross_references={OntologySource.UNIPROT: ["P00001"]},
        )
        reloaded = IDRegistry.from_dict(empty_registry.to_dict())
        found = reloaded.search_by_name("Ordinary entity")
        assert len(found) == 1


# ---------------------------------------------------------------------------
# Serialization tests
# ---------------------------------------------------------------------------

class TestSerialization:
    def test_to_dict_from_dict_roundtrip(self, seeded_registry):
        data = seeded_registry.to_dict()
        restored = IDRegistry.from_dict(data)

        # Resolve the same entity in both
        eid = ExternalID(OntologySource.CHEBI, "CHEBI:17234")
        original_cid = seeded_registry.resolve(eid)
        restored_cid = restored.resolve(eid)
        assert original_cid == restored_cid

    def test_serialization_preserves_names(self, seeded_registry):
        data = seeded_registry.to_dict()
        restored = IDRegistry.from_dict(data)
        cid = restored.resolve(ExternalID(OntologySource.UNIPROT, "P06213"))
        record = restored.get_record(cid)
        assert record.name == "Insulin receptor"

    def test_from_mapping_files(self, seeded_registry, tmp_path):
        data = seeded_registry.to_dict()
        mapping_file = tmp_path / "metabolites.json"
        metabolite_records = [
            r for r in data["records"]
            if "METABOLITE" in r["canonical_id"]
        ]
        mapping_file.write_text(json.dumps({"records": metabolite_records}))

        loaded = IDRegistry.from_mapping_files(tmp_path)
        cid = loaded.resolve(ExternalID(OntologySource.CHEBI, "CHEBI:17234"))
        assert cid.entity_type == EntityType.METABOLITE

    def test_transitive_merge_chain_survives_mapping_files_roundtrip(self, empty_registry, tmp_path):
        # Regression (MASTER_PLAN 3.3): mirrors
        # test_transitive_merge_chain_survives_serialization_roundtrip
        # (TestDeprecation, above), but through from_mapping_files() --
        # the class's OTHER documented production deserialization entry
        # point, which only restored ONE-HOP deprecations via register()
        # and never merged the deeper deprecated_index the way
        # from_dict() already did. In a chain A->B->C, C's own record
        # lists only [B]; A's own deprecation only ever lived in the
        # registry's flat _deprecated_index, which this method never
        # read back before this fix.
        a = CanonicalID(EntityType.PROTEIN, 801)
        b = empty_registry.register(
            name="Intermediate", entity_type=EntityType.PROTEIN,
            deprecated_canonical_ids=[a],
        )
        c = empty_registry.register(
            name="Final", entity_type=EntityType.PROTEIN,
            deprecated_canonical_ids=[b],
        )
        mapping_file = tmp_path / "proteins.json"
        mapping_file.write_text(json.dumps(empty_registry.to_dict()))

        loaded = IDRegistry.from_mapping_files(tmp_path)
        with pytest.raises(ResolutionError) as exc_info:
            loaded.get_record(a)
        assert str(c) in str(exc_info.value)
        assert loaded.get_record(c).name == "Final"

    def test_from_mapping_files_without_deprecated_index_key_still_loads(self, seeded_registry, tmp_path):
        # Backward compatibility: a mapping file written before this fix
        # (or by anything emitting only the bare records list, matching
        # from_mapping_files()'s own original docstring) has no
        # "deprecated_index" key at all -- must still load correctly,
        # just without any deeper merge-chain history restored, the same
        # graceful degradation from_dict()'s own data.get(..., {})
        # fallback already established.
        data = seeded_registry.to_dict()
        mapping_file = tmp_path / "metabolites.json"
        metabolite_records = [
            r for r in data["records"] if "METABOLITE" in r["canonical_id"]
        ]
        mapping_file.write_text(json.dumps({"records": metabolite_records}))

        loaded = IDRegistry.from_mapping_files(tmp_path)
        cid = loaded.resolve(ExternalID(OntologySource.CHEBI, "CHEBI:17234"))
        assert cid.entity_type == EntityType.METABOLITE


# ---------------------------------------------------------------------------
# Thread safety tests
# ---------------------------------------------------------------------------

class TestThreadSafety:
    def test_concurrent_reads_are_safe(self, seeded_registry):
        errors = []
        eid = ExternalID(OntologySource.CHEBI, "CHEBI:17234")

        def read_loop():
            try:
                for _ in range(100):
                    seeded_registry.resolve(eid)
            except Exception as e:
                errors.append(e)

        threads = [threading.Thread(target=read_loop) for _ in range(20)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()

        assert errors == [], f"Thread safety errors: {errors}"

    def test_concurrent_registrations_are_safe(self, empty_registry):
        errors = []

        def register_entities(start_idx: int):
            try:
                for i in range(50):
                    empty_registry.register(
                        name=f"Protein {start_idx}_{i}",
                        entity_type=EntityType.PROTEIN,
                    )
            except Exception as e:
                errors.append(e)

        threads = [threading.Thread(target=register_entities, args=(t,)) for t in range(10)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()

        assert errors == [], f"Thread safety errors: {errors}"
        # 10 threads × 50 registrations = 500 entities
        stats = empty_registry.stats()
        assert stats["total_entities"] == 500


# ---------------------------------------------------------------------------
# Diagnostics tests
# ---------------------------------------------------------------------------

class TestDiagnostics:
    def test_stats_counts_entities(self, seeded_registry):
        stats = seeded_registry.stats()
        assert stats["total_entities"] == 4  # glucose, insulin-r, IGF-1, glycolysis

    def test_stats_counts_by_type(self, seeded_registry):
        stats = seeded_registry.stats()
        assert stats["by_entity_type"]["METABOLITE"] == 1
        assert stats["by_entity_type"]["PROTEIN"] == 2
        assert stats["by_entity_type"]["PATHWAY"] == 1

    def test_stats_counts_by_source(self, seeded_registry):
        stats = seeded_registry.stats()
        assert stats["by_ontology_source"]["kegg"] == 3  # glucose, insulin-r, glycolysis
        assert stats["by_ontology_source"]["uniprot"] == 2  # insulin-r, IGF-1

    def test_is_registered_true(self, seeded_registry):
        eid = ExternalID(OntologySource.CHEBI, "CHEBI:17234")
        assert seeded_registry.is_registered(eid) is True

    def test_is_registered_false(self, seeded_registry):
        eid = ExternalID(OntologySource.CHEBI, "CHEBI:99999999")
        assert seeded_registry.is_registered(eid) is False
