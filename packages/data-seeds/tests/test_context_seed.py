# SAVE TO: packages\data-seeds\tests\test_context_seed.py
"""
Test suite for biosim_seeds.context_seed.

MASTER_PLAN 5.1, chunks 11/12 -- the first ContextRegistry-backed
seeds, deliberately separate from registry_seed.py (IDRegistry,
molecular entities only). See biosim_context_registry.types's own
module docstring for the architectural reasoning.
"""
import pytest

from biosim_context_registry import ContextEntityType, ContextID
from biosim_seeds import create_seeded_context_registry


class TestContextSeed:
    @pytest.fixture(scope="class")
    @classmethod
    def registry(cls):
        return create_seeded_context_registry()

    def test_registry_has_entities(self, registry):
        stats = registry.stats()
        assert stats["total_entities"] > 0

    def test_chronic_kidney_disease_registered(self, registry):
        cid = registry.resolve_by_name(
            "chronic_kidney_disease", ContextEntityType.DISEASE_STATE,
        )
        assert cid == ContextID(ContextEntityType.DISEASE_STATE, "chronic_kidney_disease")

    def test_sleep_restriction_registered(self, registry):
        cid = registry.resolve_by_name(
            "sleep_restriction", ContextEntityType.PHYSIOLOGICAL_AXIS,
        )
        assert cid == ContextID(ContextEntityType.PHYSIOLOGICAL_AXIS, "sleep_restriction")

    def test_gck_mody_registered(self, registry):
        cid = registry.resolve_by_name("gck_mody", ContextEntityType.DISEASE_STATE)
        assert cid == ContextID(ContextEntityType.DISEASE_STATE, "gck_mody")

    def test_gestational_diabetes_registered(self, registry):
        cid = registry.resolve_by_name("gestational_diabetes", ContextEntityType.DISEASE_STATE)
        assert cid == ContextID(ContextEntityType.DISEASE_STATE, "gestational_diabetes")

    def test_type_1_diabetes_registered(self, registry):
        cid = registry.resolve_by_name("type_1_diabetes", ContextEntityType.DISEASE_STATE)
        assert cid == ContextID(ContextEntityType.DISEASE_STATE, "type_1_diabetes")

    def test_type_2_diabetes_registered(self, registry):
        cid = registry.resolve_by_name("type_2_diabetes", ContextEntityType.DISEASE_STATE)
        assert cid == ContextID(ContextEntityType.DISEASE_STATE, "type_2_diabetes")

    def test_insulin_resistant_registered(self, registry):
        cid = registry.resolve_by_name("insulin_resistant", ContextEntityType.DISEASE_STATE)
        assert cid == ContextID(ContextEntityType.DISEASE_STATE, "insulin_resistant")

    def test_hypothyroidism_registered(self, registry):
        cid = registry.resolve_by_name("hypothyroidism", ContextEntityType.DISEASE_STATE)
        assert cid == ContextID(ContextEntityType.DISEASE_STATE, "hypothyroidism")

    def test_chronic_kidney_disease_and_sleep_restriction_are_distinct(self, registry):
        ckd = registry.resolve_by_name("chronic_kidney_disease", ContextEntityType.DISEASE_STATE)
        sleep = registry.resolve_by_name("sleep_restriction", ContextEntityType.PHYSIOLOGICAL_AXIS)
        assert ckd != sleep

    def test_gck_mody_and_ckd_are_distinct_despite_same_entity_type(self, registry):
        ckd = registry.resolve_by_name("chronic_kidney_disease", ContextEntityType.DISEASE_STATE)
        gck_mody = registry.resolve_by_name("gck_mody", ContextEntityType.DISEASE_STATE)
        assert ckd != gck_mody

    def test_stats_reflect_both_entity_types(self, registry):
        stats = registry.stats()
        assert stats["by_entity_type"]["DISEASE_STATE"] >= 7
        assert stats["by_entity_type"]["PHYSIOLOGICAL_AXIS"] >= 1

    def test_gck_mody_record_has_a_real_description(self, registry):
        cid = registry.resolve_by_name("gck_mody", ContextEntityType.DISEASE_STATE)
        record = registry.get_record(cid)
        assert "11119017" in record.description   # Guenat et al. 2000

    def test_gestational_diabetes_record_has_a_real_description(self, registry):
        cid = registry.resolve_by_name("gestational_diabetes", ContextEntityType.DISEASE_STATE)
        record = registry.get_record(cid)
        assert "9353615" in record.description   # Kautzky-Willer et al. 1997
        assert "20190296" in record.description  # IADPSG 2010

    def test_type_1_diabetes_record_has_a_real_description(self, registry):
        cid = registry.resolve_by_name("type_1_diabetes", ContextEntityType.DISEASE_STATE)
        record = registry.get_record(cid)
        assert "3100365" in record.description   # Simonson et al. 1987

    def test_type_2_diabetes_record_has_a_real_description(self, registry):
        cid = registry.resolve_by_name("type_2_diabetes", ContextEntityType.DISEASE_STATE)
        record = registry.get_record(cid)
        assert "15823385" in record.description   # Stumvoll et al. 2005

    def test_hypothyroidism_record_discloses_the_unresolved_conflict(self, registry):
        cid = registry.resolve_by_name("hypothyroidism", ContextEntityType.DISEASE_STATE)
        record = registry.get_record(cid)
        assert "12519082" in record.description   # Rochon et al. 2003
        assert "18429948" in record.description   # Handisurya et al. 2008
        assert "1175861" in record.description    # Shah et al. 1975 (contradicts)

    def test_ckd_record_has_a_real_description(self, registry):
        cid = registry.resolve_by_name(
            "chronic_kidney_disease", ContextEntityType.DISEASE_STATE,
        )
        record = registry.get_record(cid)
        assert "26823551" in record.description   # de Boer et al. 2016

    def test_sleep_restriction_record_has_a_real_description(self, registry):
        cid = registry.resolve_by_name(
            "sleep_restriction", ContextEntityType.PHYSIOLOGICAL_AXIS,
        )
        record = registry.get_record(cid)
        assert "20585000" in record.description   # Buxton et al. 2010
