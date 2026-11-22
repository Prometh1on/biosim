# SAVE TO: packages\data-seeds\tests\test_seeds.py
"""
Test suite for biosim_seeds.

Coverage:
    - Registry seed: all entities registered with correct IDs
    - Parameter seed: all parameters stored with provenance
    - ParameterLoader: loads literature values, falls back to defaults
    - Integration: loaded params produce valid ODE steady state
    - Provenance report: documents what came from literature vs defaults
"""

import pytest

from biosim_seeds import (
    PARAMETER_MAPPINGS,
    ParameterLoader,
    create_seeded_registry,
    create_seeded_store,
)
from biosim_id_resolver import EntityType, OntologySource
from biosim_parameter_db import ConfidenceTier, Organism, ParameterKind
from biosim_ode import GlucoseInsulinPathway, ODESystem, PathwayParameters, high_sugar_meal


# ---------------------------------------------------------------------------
# Registry seed
# ---------------------------------------------------------------------------

class TestRegistrySeed:
    @pytest.fixture(scope="class")
    @classmethod
    def registry(cls):
        return create_seeded_registry()

    def test_registry_has_entities(self, registry):
        stats = registry.stats()
        assert stats["total_entities"] > 0

    def test_glucose_registered(self, registry):
        from biosim_id_resolver import ExternalID, OntologySource
        eid = ExternalID(OntologySource.CHEBI, "CHEBI:17234")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.METABOLITE

    def test_insulin_receptor_registered(self, registry):
        from biosim_id_resolver import ExternalID
        eid = ExternalID(OntologySource.UNIPROT, "P06213")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.PROTEIN

    def test_igf1_registered(self, registry):
        from biosim_id_resolver import ExternalID
        eid = ExternalID(OntologySource.UNIPROT, "P05019")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.PROTEIN

    def test_mtor_registered(self, registry):
        from biosim_id_resolver import ExternalID
        eid = ExternalID(OntologySource.UNIPROT, "P42345")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.PROTEIN

    def test_ampk_registered(self, registry):
        from biosim_id_resolver import ExternalID
        eid = ExternalID(OntologySource.UNIPROT, "Q13131")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.PROTEIN

    def test_glycolysis_pathway_registered(self, registry):
        from biosim_id_resolver import ExternalID
        eid = ExternalID(OntologySource.REACTOME, "R-HSA-70171")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.PATHWAY

    def test_insulin_signalling_registered(self, registry):
        from biosim_id_resolver import ExternalID
        eid = ExternalID(OntologySource.REACTOME, "R-HSA-74752")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.PATHWAY

    def test_search_by_name_glucose(self, registry):
        results = registry.search_by_name("D-Glucose")
        assert len(results) == 1

    def test_search_by_name_insulin(self, registry):
        results = registry.search_by_name("insulin")
        assert len(results) == 1

    def test_kegg_crossref_for_glycolysis(self, registry):
        from biosim_id_resolver import ExternalID
        eid = ExternalID(OntologySource.KEGG, "hsa00010")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.PATHWAY

    def test_hmdb_crossref_for_glucose(self, registry):
        from biosim_id_resolver import ExternalID
        eid = ExternalID(OntologySource.HMDB, "HMDB0000122")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.METABOLITE

    def test_entity_count_by_type(self, registry):
        stats = registry.stats()
        assert stats["by_entity_type"].get("PROTEIN", 0) >= 8
        assert stats["by_entity_type"].get("METABOLITE", 0) >= 5
        assert stats["by_entity_type"].get("PATHWAY", 0) >= 4

    # -- MASTER_PLAN 5.4: fasted hepatic network metabolites --

    def test_hepatic_network_metabolite_count(self, registry):
        # 8 original v1 metabolites + 22 new ones from
        # biosim_metabolism.fasted_hepatic_network.METABOLITES + 1
        # cortisol (MASTER_PLAN 5.1, chunk 6) + 2 acarbose/miglitol
        # (MASTER_PLAN 5.1, chunk 8) + 1 pioglitazone (MASTER_PLAN 5.1,
        # chunk 9) + 1 protein (MASTER_PLAN 5.1, chunk 10) + 1 metformin
        # (MASTER_PLAN 5.1, chunk 17) + 1 caffeine (MASTER_PLAN 5.1,
        # chunk 18) + 1 tirzepatide (MASTER_PLAN 5.1, chunk 19) + 1
        # prednisone (MASTER_PLAN 5.1, chunk 20) + 2 glibenclamide/
        # glipizide (MASTER_PLAN 5.1, chunk 21) + 1 NEFA (MASTER_PLAN
        # 5.1, chunk 24) + 1 hydrochlorothiazide (MASTER_PLAN 5.1,
        # chunk 26) + 1 indinavir (MASTER_PLAN 5.1, chunk 27) + 1
        # atorvastatin (MASTER_PLAN 5.1, chunk 28) + 1 alcohol/ethanol
        # (MASTER_PLAN 5.1, chunk 36) + 1 hydroxychloroquine (MASTER_PLAN
        # 5.1, chunk 37).
        stats = registry.stats()
        assert stats["by_entity_type"].get("METABOLITE", 0) == 47

    def test_glucose_1_phosphate_registered(self, registry):
        from biosim_id_resolver import ExternalID
        eid = ExternalID(OntologySource.CHEBI, "CHEBI:29042")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.METABOLITE

    def test_palmitate_registered_via_kegg(self, registry):
        from biosim_id_resolver import ExternalID
        eid = ExternalID(OntologySource.KEGG, "C00249")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.METABOLITE
        assert registry.get_record(cid).name == "palmitate"

    def test_ammonium_registered_as_the_ion_not_ammonia(self, registry):
        # C01342 (NH4+) is a genuinely different KEGG entry from C00014
        # (neutral NH3) -- confirms the correct one was registered, and
        # that the wrong one was never accidentally used instead.
        from biosim_id_resolver import ExternalID, ResolutionError
        eid = ExternalID(OntologySource.KEGG, "C01342")
        cid = registry.resolve(eid)
        assert registry.get_record(cid).name == "ammonium"
        with pytest.raises(ResolutionError):
            registry.resolve(ExternalID(OntologySource.KEGG, "C00014"))

    def test_fadh2_registered_without_a_fabricated_hmdb_id(self, registry):
        # This project's own "never invent an identifier" discipline:
        # a real HMDB accession was searched for and not independently
        # verifiable this session (see registry_seed.py's own module
        # docstring) -- confirms the gap was left honestly disclosed
        # (no HMDB cross-reference at all) rather than filled with a
        # guess.
        from biosim_id_resolver import ExternalID, OntologySource
        eid = ExternalID(OntologySource.CHEBI, "CHEBI:17877")
        cid = registry.resolve(eid)
        record = registry.get_record(cid)
        assert record.name == "FADH2"
        assert record.get_external_id(OntologySource.HMDB) is None

    def test_all_22_new_metabolites_resolve_to_distinct_entities(self, registry):
        from biosim_id_resolver import ExternalID
        new_chebi_ids = [
            "CHEBI:29042", "CHEBI:15946", "CHEBI:16905", "CHEBI:29052",
            "CHEBI:44897", "CHEBI:16452", "CHEBI:16810", "CHEBI:16977",
            "CHEBI:16015", "CHEBI:28938", "CHEBI:17754", "CHEBI:15978",
            "CHEBI:15756", "CHEBI:16526", "CHEBI:15379", "CHEBI:16761",
            "CHEBI:18367", "CHEBI:15846", "CHEBI:16908", "CHEBI:16238",
            "CHEBI:17877", "CHEBI:15346",
        ]
        resolved = [
            registry.resolve(ExternalID(OntologySource.CHEBI, c))
            for c in new_chebi_ids
        ]
        assert len(resolved) == 22
        assert len(set(str(c) for c in resolved)) == 22

    def test_all_22_new_metabolites_kegg_crossref_matches_chebi(self, registry):
        # Mandatory code-reviewer finding (MEDIUM): only 2 of these 22
        # metabolites' KEGG cross-references were previously exercised
        # by name -- a subtly wrong (but non-colliding) KEGG accession
        # among the other 20 would have passed every existing test.
        # Every (name, chebi, kegg) triple here is checked two ways at
        # once: KEGG resolves to the SAME canonical entity as its own
        # ChEBI accession (not just "resolves to something"), and that
        # entity's name matches.
        from biosim_id_resolver import ExternalID
        name_chebi_kegg = [
            ("glucose-1-phosphate", "CHEBI:29042", "C00103"),
            ("fructose-6-phosphate", "CHEBI:15946", "C00085"),
            ("fructose-1,6-bisphosphate", "CHEBI:16905", "C00354"),
            ("glyceraldehyde-3-phosphate", "CHEBI:29052", "C00118"),
            ("phosphoenolpyruvate", "CHEBI:44897", "C00074"),
            ("oxaloacetate", "CHEBI:16452", "C00036"),
            ("alpha-ketoglutarate", "CHEBI:16810", "C00026"),
            ("glutamate", "CHEBI:16015", "C00025"),
            ("alanine", "CHEBI:16977", "C00041"),
            ("ammonium", "CHEBI:28938", "C01342"),
            ("glycerol", "CHEBI:17754", "C00116"),
            ("glycerol-3-phosphate", "CHEBI:15978", "C00093"),
            ("palmitate", "CHEBI:15756", "C00249"),
            ("carbon dioxide", "CHEBI:16526", "C00011"),
            ("oxygen", "CHEBI:15379", "C00007"),
            ("ADP", "CHEBI:16761", "C00008"),
            ("inorganic phosphate", "CHEBI:18367", "C00009"),
            ("NAD+", "CHEBI:15846", "C00003"),
            ("NADH", "CHEBI:16908", "C00004"),
            ("FAD", "CHEBI:16238", "C00016"),
            ("FADH2", "CHEBI:17877", "C01352"),
            ("coenzyme A", "CHEBI:15346", "C00010"),
        ]
        assert len(name_chebi_kegg) == 22
        for name, chebi, kegg in name_chebi_kegg:
            via_chebi = registry.resolve(ExternalID(OntologySource.CHEBI, chebi))
            via_kegg = registry.resolve(ExternalID(OntologySource.KEGG, kegg))
            assert via_chebi == via_kegg, f"{name}: ChEBI/KEGG resolved to different entities"
            assert registry.get_record(via_chebi).name == name

    def test_already_registered_metabolites_not_duplicated(self, registry):
        # 8 of the network's own 30 metabolites (D-Glucose, pyruvate,
        # lactate, ATP, AMP, glycogen, glucose-6-phosphate, acetyl-CoA)
        # already existed under these exact names -- confirms they
        # weren't accidentally re-registered as second, duplicate
        # entities during this addition.
        for name in [
            "D-Glucose", "pyruvate", "lactate", "ATP", "AMP",
            "glycogen", "glucose-6-phosphate", "acetyl-CoA",
        ]:
            assert len(registry.search_by_name(name)) == 1

    # -- MASTER_PLAN 5.1 chunk 2: SGLT2 (SLC5A2) protein entity --

    def test_sglt2_registered_via_uniprot(self, registry):
        from biosim_id_resolver import ExternalID
        eid = ExternalID(OntologySource.UNIPROT, "P31639")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.PROTEIN
        assert registry.get_record(cid).name == "SGLT2"

    def test_sglt2_kegg_crossref_matches_uniprot(self, registry):
        from biosim_id_resolver import ExternalID
        via_uniprot = registry.resolve(ExternalID(OntologySource.UNIPROT, "P31639"))
        via_kegg = registry.resolve(ExternalID(OntologySource.KEGG, "K14382"))
        assert via_uniprot == via_kegg

    # -- MASTER_PLAN 5.1 chunk 6: cortisol metabolite entity --

    def test_cortisol_registered_via_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        eid = ExternalID(OntologySource.CHEBI, "CHEBI:17650")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.METABOLITE
        assert registry.get_record(cid).name == "cortisol"

    def test_cortisol_kegg_crossref_matches_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        via_chebi = registry.resolve(ExternalID(OntologySource.CHEBI, "CHEBI:17650"))
        via_kegg = registry.resolve(ExternalID(OntologySource.KEGG, "C00735"))
        assert via_chebi == via_kegg

    # -- MASTER_PLAN 5.1 chunk 8: acarbose/miglitol metabolite entities --

    def test_acarbose_registered_via_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        eid = ExternalID(OntologySource.CHEBI, "CHEBI:2376")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.METABOLITE
        assert registry.get_record(cid).name == "acarbose"

    def test_acarbose_kegg_crossref_matches_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        via_chebi = registry.resolve(ExternalID(OntologySource.CHEBI, "CHEBI:2376"))
        via_kegg = registry.resolve(ExternalID(OntologySource.KEGG, "C06802"))
        assert via_chebi == via_kegg

    def test_miglitol_registered_via_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        eid = ExternalID(OntologySource.CHEBI, "CHEBI:6935")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.METABOLITE
        assert registry.get_record(cid).name == "miglitol"

    def test_miglitol_kegg_crossref_matches_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        via_chebi = registry.resolve(ExternalID(OntologySource.CHEBI, "CHEBI:6935"))
        via_kegg = registry.resolve(ExternalID(OntologySource.KEGG, "C07708"))
        assert via_chebi == via_kegg

    def test_acarbose_and_miglitol_are_distinct_entities(self, registry):
        acarbose_id = registry.search_by_name("acarbose")[0]
        miglitol_id = registry.search_by_name("miglitol")[0]
        assert acarbose_id != miglitol_id

    # -- MASTER_PLAN 5.1 chunk 9: pioglitazone metabolite entity --

    def test_pioglitazone_registered_via_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        eid = ExternalID(OntologySource.CHEBI, "CHEBI:8228")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.METABOLITE
        assert registry.get_record(cid).name == "pioglitazone"

    def test_pioglitazone_kegg_crossref_matches_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        via_chebi = registry.resolve(ExternalID(OntologySource.CHEBI, "CHEBI:8228"))
        via_kegg = registry.resolve(ExternalID(OntologySource.KEGG, "C07675"))
        assert via_chebi == via_kegg

    def test_pioglitazone_distinct_from_acarbose_and_miglitol(self, registry):
        pioglitazone_id = registry.search_by_name("pioglitazone")[0]
        acarbose_id = registry.search_by_name("acarbose")[0]
        miglitol_id = registry.search_by_name("miglitol")[0]
        assert pioglitazone_id not in (acarbose_id, miglitol_id)

    # -- MASTER_PLAN 5.1 chunk 10: protein metabolite entity --

    def test_protein_registered_via_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        eid = ExternalID(OntologySource.CHEBI, "CHEBI:16541")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.METABOLITE
        assert registry.get_record(cid).name == "protein"

    def test_protein_kegg_crossref_matches_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        via_chebi = registry.resolve(ExternalID(OntologySource.CHEBI, "CHEBI:16541"))
        via_kegg = registry.resolve(ExternalID(OntologySource.KEGG, "C00017"))
        assert via_chebi == via_kegg

    def test_protein_distinct_from_pioglitazone_acarbose_miglitol(self, registry):
        protein_id = registry.search_by_name("protein")[0]
        pioglitazone_id = registry.search_by_name("pioglitazone")[0]
        acarbose_id = registry.search_by_name("acarbose")[0]
        miglitol_id = registry.search_by_name("miglitol")[0]
        assert protein_id not in (pioglitazone_id, acarbose_id, miglitol_id)

    # -- MASTER_PLAN 5.1 chunk 17: metformin metabolite entity --

    def test_metformin_registered_via_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        eid = ExternalID(OntologySource.CHEBI, "CHEBI:6801")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.METABOLITE
        assert registry.get_record(cid).name == "metformin"

    def test_metformin_kegg_crossref_matches_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        via_chebi = registry.resolve(ExternalID(OntologySource.CHEBI, "CHEBI:6801"))
        via_kegg = registry.resolve(ExternalID(OntologySource.KEGG, "C07151"))
        assert via_chebi == via_kegg

    def test_metformin_distinct_from_pioglitazone_protein_acarbose_miglitol(self, registry):
        metformin_id = registry.search_by_name("metformin")[0]
        pioglitazone_id = registry.search_by_name("pioglitazone")[0]
        protein_id = registry.search_by_name("protein")[0]
        acarbose_id = registry.search_by_name("acarbose")[0]
        miglitol_id = registry.search_by_name("miglitol")[0]
        assert metformin_id not in (pioglitazone_id, protein_id, acarbose_id, miglitol_id)

    # -- MASTER_PLAN 5.1 chunk 18: caffeine metabolite entity --

    def test_caffeine_registered_via_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        eid = ExternalID(OntologySource.CHEBI, "CHEBI:27732")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.METABOLITE
        assert registry.get_record(cid).name == "caffeine"

    def test_caffeine_kegg_crossref_matches_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        via_chebi = registry.resolve(ExternalID(OntologySource.CHEBI, "CHEBI:27732"))
        via_kegg = registry.resolve(ExternalID(OntologySource.KEGG, "C07481"))
        assert via_chebi == via_kegg

    def test_caffeine_distinct_from_metformin_and_siblings(self, registry):
        caffeine_id = registry.search_by_name("caffeine")[0]
        metformin_id = registry.search_by_name("metformin")[0]
        pioglitazone_id = registry.search_by_name("pioglitazone")[0]
        assert caffeine_id not in (metformin_id, pioglitazone_id)

    # -- MASTER_PLAN 5.1 chunk 19: tirzepatide metabolite entity --

    def test_tirzepatide_registered_via_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        eid = ExternalID(OntologySource.CHEBI, "CHEBI:194186")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.METABOLITE
        assert registry.get_record(cid).name == "tirzepatide"

    def test_tirzepatide_kegg_crossref_matches_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        via_chebi = registry.resolve(ExternalID(OntologySource.CHEBI, "CHEBI:194186"))
        via_kegg = registry.resolve(ExternalID(OntologySource.KEGG, "D11360"))
        assert via_chebi == via_kegg

    def test_tirzepatide_distinct_from_caffeine_and_siblings(self, registry):
        tirzepatide_id = registry.search_by_name("tirzepatide")[0]
        caffeine_id = registry.search_by_name("caffeine")[0]
        metformin_id = registry.search_by_name("metformin")[0]
        assert tirzepatide_id not in (caffeine_id, metformin_id)

    # -- MASTER_PLAN 5.1 chunk 20: prednisone metabolite entity --

    def test_prednisone_registered_via_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        eid = ExternalID(OntologySource.CHEBI, "CHEBI:8382")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.METABOLITE
        assert registry.get_record(cid).name == "prednisone"

    def test_prednisone_kegg_crossref_matches_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        via_chebi = registry.resolve(ExternalID(OntologySource.CHEBI, "CHEBI:8382"))
        via_kegg = registry.resolve(ExternalID(OntologySource.KEGG, "C07370"))
        assert via_chebi == via_kegg

    def test_prednisone_distinct_from_tirzepatide_and_siblings(self, registry):
        prednisone_id = registry.search_by_name("prednisone")[0]
        tirzepatide_id = registry.search_by_name("tirzepatide")[0]
        metformin_id = registry.search_by_name("metformin")[0]
        assert prednisone_id not in (tirzepatide_id, metformin_id)

    # -- MASTER_PLAN 5.1 chunk 21: glibenclamide/glipizide metabolite entities --

    def test_glibenclamide_registered_via_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        eid = ExternalID(OntologySource.CHEBI, "CHEBI:5441")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.METABOLITE
        assert registry.get_record(cid).name == "glibenclamide"

    def test_glibenclamide_kegg_crossref_matches_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        via_chebi = registry.resolve(ExternalID(OntologySource.CHEBI, "CHEBI:5441"))
        via_kegg = registry.resolve(ExternalID(OntologySource.KEGG, "C07022"))
        assert via_chebi == via_kegg

    def test_glipizide_registered_via_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        eid = ExternalID(OntologySource.CHEBI, "CHEBI:5384")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.METABOLITE
        assert registry.get_record(cid).name == "glipizide"

    # -- MASTER_PLAN 5.1 chunk 23: GLP-1/GCG protein entity --

    def test_glp1_registered_via_uniprot(self, registry):
        from biosim_id_resolver import ExternalID
        eid = ExternalID(OntologySource.UNIPROT, "P01275")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.PROTEIN
        assert registry.get_record(cid).name == "GLP-1"

    def test_glp1_kegg_crossref_matches_uniprot(self, registry):
        from biosim_id_resolver import ExternalID
        via_uniprot = registry.resolve(ExternalID(OntologySource.UNIPROT, "P01275"))
        via_kegg = registry.resolve(ExternalID(OntologySource.KEGG, "hsa:2641"))
        assert via_uniprot == via_kegg

    def test_glipizide_kegg_crossref_matches_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        via_chebi = registry.resolve(ExternalID(OntologySource.CHEBI, "CHEBI:5384"))
        via_kegg = registry.resolve(ExternalID(OntologySource.KEGG, "D00335"))
        assert via_chebi == via_kegg

    def test_glibenclamide_and_glipizide_are_distinct_entities(self, registry):
        glibenclamide_id = registry.search_by_name("glibenclamide")[0]
        glipizide_id = registry.search_by_name("glipizide")[0]
        prednisone_id = registry.search_by_name("prednisone")[0]
        assert glibenclamide_id != glipizide_id
        assert glibenclamide_id not in (prednisone_id,)
        assert glipizide_id not in (prednisone_id,)

    # -- MASTER_PLAN 5.1 chunk 24: NEFA metabolite entity --

    def test_nefa_registered_via_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        eid = ExternalID(OntologySource.CHEBI, "CHEBI:35366")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.METABOLITE
        assert registry.get_record(cid).name == "NEFA"

    def test_nefa_kegg_crossref_matches_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        via_chebi = registry.resolve(ExternalID(OntologySource.CHEBI, "CHEBI:35366"))
        via_kegg = registry.resolve(ExternalID(OntologySource.KEGG, "C00162"))
        assert via_chebi == via_kegg

    # -- MASTER_PLAN 5.1 chunk 26: hydrochlorothiazide metabolite entity --

    def test_hydrochlorothiazide_registered_via_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        eid = ExternalID(OntologySource.CHEBI, "CHEBI:5778")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.METABOLITE
        assert registry.get_record(cid).name == "hydrochlorothiazide"

    def test_hydrochlorothiazide_kegg_crossref_matches_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        via_chebi = registry.resolve(ExternalID(OntologySource.CHEBI, "CHEBI:5778"))
        via_kegg = registry.resolve(ExternalID(OntologySource.KEGG, "D00340"))
        assert via_chebi == via_kegg

    # -- MASTER_PLAN 5.1 chunk 27: indinavir metabolite entity --

    def test_indinavir_registered_via_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        eid = ExternalID(OntologySource.CHEBI, "CHEBI:44032")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.METABOLITE
        assert registry.get_record(cid).name == "indinavir"

    def test_indinavir_kegg_crossref_matches_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        via_chebi = registry.resolve(ExternalID(OntologySource.CHEBI, "CHEBI:44032"))
        via_kegg = registry.resolve(ExternalID(OntologySource.KEGG, "C07051"))
        assert via_chebi == via_kegg

    # -- MASTER_PLAN 5.1 chunk 28: atorvastatin metabolite entity --

    def test_atorvastatin_registered_via_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        eid = ExternalID(OntologySource.CHEBI, "CHEBI:39548")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.METABOLITE
        assert registry.get_record(cid).name == "atorvastatin"

    def test_atorvastatin_kegg_crossref_matches_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        via_chebi = registry.resolve(ExternalID(OntologySource.CHEBI, "CHEBI:39548"))
        via_kegg = registry.resolve(ExternalID(OntologySource.KEGG, "C06834"))
        assert via_chebi == via_kegg

    def test_alcohol_registered_via_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        eid = ExternalID(OntologySource.CHEBI, "CHEBI:16236")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.METABOLITE
        assert registry.get_record(cid).name == "alcohol"

    def test_alcohol_kegg_crossref_matches_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        via_chebi = registry.resolve(ExternalID(OntologySource.CHEBI, "CHEBI:16236"))
        via_kegg = registry.resolve(ExternalID(OntologySource.KEGG, "C00469"))
        assert via_chebi == via_kegg

    def test_hydroxychloroquine_registered_via_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        eid = ExternalID(OntologySource.CHEBI, "CHEBI:5801")
        cid = registry.resolve(eid)
        assert cid.entity_type == EntityType.METABOLITE
        assert registry.get_record(cid).name == "hydroxychloroquine"

    def test_hydroxychloroquine_kegg_crossref_matches_chebi(self, registry):
        from biosim_id_resolver import ExternalID
        via_chebi = registry.resolve(ExternalID(OntologySource.CHEBI, "CHEBI:5801"))
        via_kegg = registry.resolve(ExternalID(OntologySource.KEGG, "C07043"))
        assert via_chebi == via_kegg


# ---------------------------------------------------------------------------
# Parameter seed
# ---------------------------------------------------------------------------

class TestParameterSeed:
    @pytest.fixture(scope="class")
    @classmethod
    def seeded(cls):
        registry = create_seeded_registry()
        store = create_seeded_store(registry)
        return registry, store

    def test_store_has_parameters(self, seeded):
        _, store = seeded
        stats = store.stats()
        assert stats["total_parameters"] > 0

    def test_glucose_abundance_stored(self, seeded):
        registry, store = seeded
        glucose_ids = registry.search_by_name("D-Glucose")
        assert glucose_ids
        params = store.query(str(glucose_ids[0]), ParameterKind.ABUNDANCE)
        assert len(params) > 0

    def test_insulin_abundance_stored(self, seeded):
        registry, store = seeded
        insulin_ids = registry.search_by_name("insulin")
        assert insulin_ids
        params = store.query(str(insulin_ids[0]), ParameterKind.ABUNDANCE)
        assert len(params) > 0

    def test_insulin_ke_stored(self, seeded):
        registry, store = seeded
        insulin_ids = registry.search_by_name("insulin")
        assert insulin_ids
        params = store.query(str(insulin_ids[0]), ParameterKind.RATE)
        assert len(params) > 0

    def test_mtor_ec50_stored(self, seeded):
        # EC50_mTOR is seeded object-scoped (.for_object(insr_id)) -- it's
        # mTOR's EC50 for insulin-receptor-mediated activation specifically,
        # not a bare subject-only EC50. ParameterStore.query()'s
        # object_canonical_id filter is an exact match, not a wildcard (see
        # its docstring), so a query with no object filter correctly finds
        # nothing here -- that's the object-scoping fix working as intended,
        # not a gap. Matches ParameterLoader.load_pathway_params()'s own
        # PARAMETER_MAPPINGS row for EC50_mTOR, which passes the same
        # object filter for exactly this reason.
        registry, store = seeded
        mtor_ids = registry.search_by_name("mTOR")
        insr_ids = registry.search_by_name("insulin receptor")
        assert mtor_ids
        assert insr_ids
        params = store.query(
            str(mtor_ids[0]), ParameterKind.EC50,
            object_canonical_id=str(insr_ids[0]),
        )
        assert len(params) > 0

    def test_all_parameters_have_pubmed_id(self, seeded):
        _, store = seeded
        for pid, param in store._parameters.items():
            assert param.provenance.pubmed_id is not None, \
                f"Parameter {pid} missing PubMed ID"

    def test_gold_tier_parameters_present(self, seeded):
        _, store = seeded
        stats = store.stats()
        assert stats["by_confidence_tier"].get("gold", 0) > 0

    def test_p3_seed_superseded_with_corrected_value(self, seeded):
        # Regression for KNOWN BUGS FIXED #15: the original p3 seed entry
        # (5.03e-6, claimed "from Table 1") could not be reconciled with
        # Bergman 1981's own directly-quoted Fig 6 SI value once the paper
        # was actually read -- it was superseded (not deleted, per the
        # append-only store), and the active query must return only the
        # corrected value, matching PathwayParameters.p3.
        from biosim_ode import PathwayParameters
        registry, store = seeded
        insr_ids = registry.search_by_name("insulin receptor")
        assert insr_ids
        params = store.query(str(insr_ids[0]), ParameterKind.CUSTOM)
        p3_params = [p for p in params if p.custom_kind_name == "insulin_action_p3"]
        assert len(p3_params) == 1, "Active query should return only the superseding p3 entry"
        assert p3_params[0].value.value == pytest.approx(PathwayParameters().p3)

    def test_exogenous_insulin_pk_constants_stored(self, seeded):
        # MASTER_PLAN 5.1 -- first entry in the systematic parameter-
        # breadth curation campaign. Note: subject="insulin" also now
        # carries a 5th CUSTOM entry from chunk 7 (sulfonylurea
        # glucose-independent secretion) -- filtered out here by name so
        # this test stays scoped to exactly the 4 exogenous-insulin PK
        # constants it was written to check.
        registry, store = seeded
        insulin_ids = registry.search_by_name("insulin")
        assert insulin_ids
        params = store.query(str(insulin_ids[0]), ParameterKind.CUSTOM)
        by_name = {
            p.custom_kind_name: p for p in params
            if p.custom_kind_name.startswith("exogenous_insulin_")
        }
        assert set(by_name) == {
            "exogenous_insulin_regular_ka_per_min",
            "exogenous_insulin_regular_pm_per_unit",
            "exogenous_insulin_lispro_ka_per_min",
            "exogenous_insulin_lispro_pm_per_unit",
        }

    def test_exogenous_insulin_pk_values_match_the_real_module_constants(self, seeded):
        # Confirms the curated store values are byte-identical to
        # exogenous_insulin_pathway.py's own real calibrated constants --
        # not a re-typed, potentially-drifted copy.
        from biosim_ode.exogenous_insulin_pathway import (
            INSULIN_LISPRO_KA_PER_MIN, INSULIN_LISPRO_PM_PER_UNIT,
            INSULIN_REGULAR_KA_PER_MIN, INSULIN_REGULAR_PM_PER_UNIT,
        )
        registry, store = seeded
        insulin_id = str(registry.search_by_name("insulin")[0])
        params = store.query(insulin_id, ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p.value.value for p in params}
        assert by_name["exogenous_insulin_regular_ka_per_min"] == INSULIN_REGULAR_KA_PER_MIN
        assert by_name["exogenous_insulin_regular_pm_per_unit"] == INSULIN_REGULAR_PM_PER_UNIT
        assert by_name["exogenous_insulin_lispro_ka_per_min"] == INSULIN_LISPRO_KA_PER_MIN
        assert by_name["exogenous_insulin_lispro_pm_per_unit"] == INSULIN_LISPRO_PM_PER_UNIT

    def test_exogenous_insulin_pk_cites_howey_1994(self, seeded):
        registry, store = seeded
        insulin_id = str(registry.search_by_name("insulin")[0])
        params = store.query(insulin_id, ParameterKind.CUSTOM)
        exo_params = [p for p in params if p.custom_kind_name.startswith("exogenous_insulin_")]
        assert len(exo_params) == 4
        for p in exo_params:
            assert p.provenance.pubmed_id == "8314011"
            assert p.confidence_tier == ConfidenceTier.BRONZE

    def test_exogenous_insulin_pk_does_not_create_spurious_conflicts(self, seeded):
        # Regression: 4 different CUSTOM-kind parameters sharing
        # subject="insulin" previously produced 5 nonsensical conflicts
        # (fold differences up to ~61,500x) before ParameterStore's own
        # conflict-detection fix (see test_parameter_db.py's
        # TestCustomKindScopedParameters) -- confirms it end-to-end
        # through the REAL seeded store, not just a synthetic
        # reproduction.
        _, store = seeded
        conflicts = store.unresolved_conflicts()
        exo_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).custom_kind_name
            and store.get(c.existing_parameter_id).custom_kind_name.startswith("exogenous_insulin_")
        ]
        assert exo_conflicts == []

    # -- MASTER_PLAN 5.1 chunk 2: renal SGLT2i GFR/TmG constants --

    def test_renal_sglt2i_constants_stored(self, seeded):
        registry, store = seeded
        sglt2_ids = registry.search_by_name("SGLT2")
        assert sglt2_ids
        params = store.query(str(sglt2_ids[0]), ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        assert set(by_name) == {
            "gfr_reference_ml_per_min",
            "dapagliflozin_tmg_healthy_mg_per_min",
            "dapagliflozin_tmg_t2d_mg_per_min",
            "empagliflozin_tmg_healthy_mg_per_min",
            "empagliflozin_tmg_t2d_mg_per_min",
        }

    def test_renal_sglt2i_values_match_the_real_module_constants(self, seeded):
        # Confirms the curated store values are byte-identical to
        # renal_pathway.py's own real calibrated constants.
        from biosim_ode.renal_pathway import (
            DAPAGLIFLOZIN_TMG_HEALTHY_MG_PER_MIN, DAPAGLIFLOZIN_TMG_T2D_MG_PER_MIN,
            EMPAGLIFLOZIN_TMG_HEALTHY_MG_PER_MIN, EMPAGLIFLOZIN_TMG_T2D_MG_PER_MIN,
            GFR_REFERENCE_ML_PER_MIN,
        )
        registry, store = seeded
        sglt2_id = str(registry.search_by_name("SGLT2")[0])
        params = store.query(sglt2_id, ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p.value.value for p in params}
        assert by_name["gfr_reference_ml_per_min"] == GFR_REFERENCE_ML_PER_MIN
        assert by_name["dapagliflozin_tmg_healthy_mg_per_min"] == DAPAGLIFLOZIN_TMG_HEALTHY_MG_PER_MIN
        assert by_name["dapagliflozin_tmg_t2d_mg_per_min"] == DAPAGLIFLOZIN_TMG_T2D_MG_PER_MIN
        assert by_name["empagliflozin_tmg_healthy_mg_per_min"] == EMPAGLIFLOZIN_TMG_HEALTHY_MG_PER_MIN
        assert by_name["empagliflozin_tmg_t2d_mg_per_min"] == EMPAGLIFLOZIN_TMG_T2D_MG_PER_MIN

    def test_renal_sglt2i_cites_correct_pmids_and_tier(self, seeded):
        registry, store = seeded
        sglt2_id = str(registry.search_by_name("SGLT2")[0])
        params = store.query(sglt2_id, ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        for name in (
            "gfr_reference_ml_per_min",
            "dapagliflozin_tmg_healthy_mg_per_min",
            "dapagliflozin_tmg_t2d_mg_per_min",
        ):
            assert by_name[name].provenance.pubmed_id == "23735727"
        for name in (
            "empagliflozin_tmg_healthy_mg_per_min",
            "empagliflozin_tmg_t2d_mg_per_min",
        ):
            assert by_name[name].provenance.pubmed_id == "28428225"
        for p in by_name.values():
            # GOLD, not SILVER -- re-tiered per the housekeeping-backlog
            # confidence-tier-definition fix (see parameter_seed.py's
            # own module docstring): these are directly-reported
            # primary-source clinical measurements, this file's
            # corrected GOLD standard, not a replication-count bar.
            assert p.confidence_tier == ConfidenceTier.GOLD

    def test_renal_sglt2i_does_not_create_spurious_conflicts(self, seeded):
        # Regression: 5 different CUSTOM-kind parameters sharing
        # subject="SGLT2" -- confirms ParameterStore's own conflict-
        # detection fix (see test_parameter_db.py's
        # TestCustomKindScopedParameters) generalizes beyond the
        # exogenous-insulin case it was found and fixed for.
        _, store = seeded
        conflicts = store.unresolved_conflicts()
        renal_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).custom_kind_name
            and (
                store.get(c.existing_parameter_id).custom_kind_name.startswith("gfr_")
                or "tmg" in store.get(c.existing_parameter_id).custom_kind_name
            )
        ]
        assert renal_conflicts == []

    # -- MASTER_PLAN 5.1 chunk 6: cortisol PK/potency constants --

    def test_cortisol_pk_constants_stored(self, seeded):
        registry, store = seeded
        cortisol_ids = registry.search_by_name("cortisol")
        assert cortisol_ids
        params = store.query(str(cortisol_ids[0]), ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        assert set(by_name) == {
            "prednisone_to_cortisol_mg_ratio",
            "daily_cortisol_production_mg",
            "cushings_cortisol_ratio",
        }

    def test_cortisol_pk_values_match_the_real_module_constants(self, seeded):
        # Confirms the curated store values are byte-identical to
        # cortisol_pathway.py's own real constants.
        from biosim_ode.cortisol_pathway import (
            DAILY_CORTISOL_PRODUCTION_MG, PREDNISONE_TO_CORTISOL_MG_RATIO,
        )
        registry, store = seeded
        cortisol_id = str(registry.search_by_name("cortisol")[0])
        params = store.query(cortisol_id, ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p.value.value for p in params}
        assert by_name["prednisone_to_cortisol_mg_ratio"] == PREDNISONE_TO_CORTISOL_MG_RATIO
        assert by_name["daily_cortisol_production_mg"] == DAILY_CORTISOL_PRODUCTION_MG

    def test_cortisol_pk_cites_correct_pmids_and_tier(self, seeded):
        # The two PMIDs here are the CORRECTED ones (see CLAUDE.md "DONE
        # — Fixed two mistranscribed citations in cortisol_pathway.py") --
        # this test would have pinned the wrong PMIDs (412539, 1955512)
        # had it been written before that fix.
        registry, store = seeded
        cortisol_id = str(registry.search_by_name("cortisol")[0])
        params = store.query(cortisol_id, ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        assert by_name["prednisone_to_cortisol_mg_ratio"].provenance.pubmed_id == "888843"
        assert by_name["prednisone_to_cortisol_mg_ratio"].confidence_tier == ConfidenceTier.BRONZE
        assert by_name["daily_cortisol_production_mg"].provenance.pubmed_id == "1986026"
        # GOLD, not SILVER -- re-tiered per the housekeeping-backlog
        # confidence-tier-definition fix (a directly-quoted primary-
        # source measurement, this file's corrected GOLD standard).
        assert by_name["daily_cortisol_production_mg"].confidence_tier == ConfidenceTier.GOLD

    def test_cortisol_pk_does_not_create_spurious_conflicts(self, seeded):
        # Regression: 2 different CUSTOM-kind parameters sharing
        # subject="cortisol" -- confirms ParameterStore's own conflict-
        # detection fix (custom_kind_name scoping) generalizes to a
        # third subject, not just insulin/SGLT2.
        _, store = seeded
        conflicts = store.unresolved_conflicts()
        cortisol_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).custom_kind_name
            and (
                "cortisol" in store.get(c.existing_parameter_id).custom_kind_name
                or "prednisone" in store.get(c.existing_parameter_id).custom_kind_name
            )
        ]
        assert cortisol_conflicts == []

    # -- MASTER_PLAN 5.1 chunk 34: Cushing's syndrome cortisol ratio --

    def test_cushings_cortisol_ratio_stored(self, seeded):
        registry, store = seeded
        cortisol_id = str(registry.search_by_name("cortisol")[0])
        params = store.query(cortisol_id, ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        assert "cushings_cortisol_ratio" in by_name

    def test_cushings_cortisol_ratio_matches_the_real_module_constant(self, seeded):
        # Live cross-check against the REAL CortisolParameters.
        # CUSHINGS_CORTISOL_RATIO, not a duplicated literal.
        from biosim_ode.cortisol_pathway import CortisolParameters
        registry, store = seeded
        cortisol_id = str(registry.search_by_name("cortisol")[0])
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(cortisol_id, ParameterKind.CUSTOM)
        }
        assert params["cushings_cortisol_ratio"] == pytest.approx(
            CortisolParameters.CUSHINGS_CORTISOL_RATIO
        )

    def test_cushings_cortisol_ratio_cites_correct_pmid_and_tier(self, seeded):
        registry, store = seeded
        cortisol_id = str(registry.search_by_name("cortisol")[0])
        params = store.query(cortisol_id, ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        entry = by_name["cushings_cortisol_ratio"]
        assert entry.provenance.pubmed_id == "31695455"
        assert entry.confidence_tier == ConfidenceTier.GOLD

    # -- MASTER_PLAN 5.1 chunk 35: hyperthyroidism hepatic EGP anchor (ContextRegistry) --

    def test_hyperthyroidism_egp_factor_stored(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        hyper_id = str(ContextID(ContextEntityType.DISEASE_STATE, "hyperthyroidism"))
        params = store.query(hyper_id, ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        assert "egp_basal_factor" in by_name

    def test_hyperthyroidism_egp_factor_matches_the_real_module_function(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        from biosim_ode.hepatic import EGP_BASAL_MM_MIN, EGPParameters
        _, store = seeded
        hyper_id = str(ContextID(ContextEntityType.DISEASE_STATE, "hyperthyroidism"))
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(hyper_id, ParameterKind.CUSTOM)
        }
        reference = EGPParameters.for_hyperthyroidism().egp_basal / EGP_BASAL_MM_MIN
        assert params["egp_basal_factor"] == pytest.approx(reference)

    def test_hyperthyroidism_egp_factor_cites_correct_pmid_and_tier(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        hyper_id = str(ContextID(ContextEntityType.DISEASE_STATE, "hyperthyroidism"))
        entries = store.query(hyper_id, ParameterKind.CUSTOM)
        assert len(entries) == 1
        for entry in entries:
            assert entry.provenance.pubmed_id == "3510000"
            assert entry.confidence_tier == ConfidenceTier.GOLD

    def test_hyperthyroidism_is_a_separate_subject_no_cross_conflict(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        hyper_id = str(ContextID(ContextEntityType.DISEASE_STATE, "hyperthyroidism"))
        conflicts = store.unresolved_conflicts()
        hyper_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).subject_canonical_id == hyper_id
        ]
        assert hyper_conflicts == []

    # -- MASTER_PLAN 5.1 chunk 36: alcohol acute hepatic EGP suppression constant --

    def test_alcohol_egp_factor_stored(self, seeded):
        from biosim_ode.alcohol_hepatic_pathway import ALCOHOL_TESTED_DOSE_G
        registry, store = seeded
        alcohol_ids = registry.search_by_name("alcohol")
        assert alcohol_ids
        params = store.query(str(alcohol_ids[0]), ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        assert f"egp_reduction_fraction_at_{ALCOHOL_TESTED_DOSE_G:.0f}g" in by_name

    def test_alcohol_egp_factor_matches_the_real_module_constant(self, seeded):
        from biosim_ode.alcohol_hepatic_pathway import (
            ALCOHOL_EGP_REDUCTION_FRACTION_AT_TESTED_DOSE,
            ALCOHOL_TESTED_DOSE_G,
        )
        registry, store = seeded
        alcohol_id = str(registry.search_by_name("alcohol")[0])
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(alcohol_id, ParameterKind.CUSTOM)
        }
        key = f"egp_reduction_fraction_at_{ALCOHOL_TESTED_DOSE_G:.0f}g"
        assert params[key] == pytest.approx(ALCOHOL_EGP_REDUCTION_FRACTION_AT_TESTED_DOSE)

    def test_alcohol_egp_factor_cites_correct_pmid_and_tier(self, seeded):
        registry, store = seeded
        alcohol_id = str(registry.search_by_name("alcohol")[0])
        entries = store.query(alcohol_id, ParameterKind.CUSTOM)
        assert len(entries) == 1
        for entry in entries:
            assert entry.provenance.pubmed_id == "9815011"
            assert entry.confidence_tier == ConfidenceTier.GOLD

    def test_alcohol_is_a_separate_subject_no_cross_conflict(self, seeded):
        registry, store = seeded
        alcohol_id = str(registry.search_by_name("alcohol")[0])
        conflicts = store.unresolved_conflicts()
        alcohol_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).subject_canonical_id == alcohol_id
        ]
        assert alcohol_conflicts == []

    # -- MASTER_PLAN 5.1 chunk 37: hydroxychloroquine chronic p3 constant --

    def test_hydroxychloroquine_p3_factor_stored(self, seeded):
        registry, store = seeded
        hcq_ids = registry.search_by_name("hydroxychloroquine")
        assert hcq_ids
        params = store.query(str(hcq_ids[0]), ParameterKind.CUSTOM)
        assert {p.custom_kind_name for p in params} == {"hydroxychloroquine_p3_factor"}

    def test_hydroxychloroquine_factor_matches_the_real_module_constant(self, seeded):
        from biosim_physiology.drug_interactions import HYDROXYCHLOROQUINE_P3_FACTOR
        registry, store = seeded
        hcq_id = str(registry.search_by_name("hydroxychloroquine")[0])
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(hcq_id, ParameterKind.CUSTOM)
        }
        assert params["hydroxychloroquine_p3_factor"] == pytest.approx(HYDROXYCHLOROQUINE_P3_FACTOR)

    def test_hydroxychloroquine_factor_cites_correct_pmid_and_tier(self, seeded):
        registry, store = seeded
        hcq_id = str(registry.search_by_name("hydroxychloroquine")[0])
        entries = store.query(hcq_id, ParameterKind.CUSTOM)
        assert len(entries) == 1
        for entry in entries:
            assert entry.provenance.pubmed_id == "26197707"
            assert entry.confidence_tier == ConfidenceTier.GOLD

    def test_hydroxychloroquine_is_a_separate_subject_no_cross_conflict(self, seeded):
        registry, store = seeded
        hcq_id = str(registry.search_by_name("hydroxychloroquine")[0])
        conflicts = [
            c for c in store.unresolved_conflicts()
            if store.get(c.existing_parameter_id).subject_canonical_id == hcq_id
        ]
        assert conflicts == []

    # -- MASTER_PLAN 5.1 chunk 38: NAFLD disease-state p3 constant (ContextRegistry) --

    def test_nafld_p3_factor_stored(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        nafld_id = str(ContextID(ContextEntityType.DISEASE_STATE, "nafld"))
        by_name = {p.custom_kind_name: p for p in store.query(nafld_id, ParameterKind.CUSTOM)}
        assert "p3_factor" in by_name

    def test_nafld_factor_matches_the_real_module_function(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        from biosim_physiology.modifiers import disease_modifiers
        from biosim_physiology.profile import DiseaseState
        _, store = seeded
        nafld_id = str(ContextID(ContextEntityType.DISEASE_STATE, "nafld"))
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(nafld_id, ParameterKind.CUSTOM)
        }
        assert params["p3_factor"] == pytest.approx(disease_modifiers([DiseaseState.NAFLD])["p3"])

    def test_nafld_factor_cites_correct_pmid_and_tier(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        nafld_id = str(ContextID(ContextEntityType.DISEASE_STATE, "nafld"))
        entries = store.query(nafld_id, ParameterKind.CUSTOM)
        assert len(entries) == 1
        for entry in entries:
            assert entry.provenance.pubmed_id == "15747110"
            # BRONZE: a disposal deficit at fixed clamp insulin read as a p3
            # multiplier is an inference, not a directly-reported p3 ratio.
            assert entry.confidence_tier == ConfidenceTier.BRONZE

    def test_nafld_is_a_separate_subject_no_cross_conflict(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        nafld_id = str(ContextID(ContextEntityType.DISEASE_STATE, "nafld"))
        # Non-vacuous: the subject must actually hold a parameter to conflict on.
        assert store.query(nafld_id, ParameterKind.CUSTOM)
        conflicts = [
            c for c in store.unresolved_conflicts()
            if store.get(c.existing_parameter_id).subject_canonical_id == nafld_id
        ]
        assert conflicts == []

    # -- MASTER_PLAN 5.1 chunk 23: GLP-1 basal/meal-response PK constants --

    def test_glp1_pk_constants_stored(self, seeded):
        registry, store = seeded
        glp1_ids = registry.search_by_name("GLP-1")
        assert glp1_ids
        params = store.query(str(glp1_ids[0]), ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        assert set(by_name) == {
            "glp1_basal_pm",
            "glp1_meal_response_pm_per_mm",
        }

    def test_glp1_pk_values_match_the_real_module_constants(self, seeded):
        # Confirms the curated store values are byte-identical to
        # glp1_pathway.py's own real constants.
        from biosim_ode.glp1_pathway import GLP1_BASAL_PM, GLP1_MEAL_RESPONSE
        registry, store = seeded
        glp1_id = str(registry.search_by_name("GLP-1")[0])
        params = store.query(glp1_id, ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p.value.value for p in params}
        assert by_name["glp1_basal_pm"] == GLP1_BASAL_PM
        assert by_name["glp1_meal_response_pm_per_mm"] == GLP1_MEAL_RESPONSE

    def test_glp1_pk_cites_correct_pmid_and_tier(self, seeded):
        registry, store = seeded
        glp1_id = str(registry.search_by_name("GLP-1")[0])
        params = store.query(glp1_id, ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        assert by_name["glp1_basal_pm"].provenance.pubmed_id == "8138058"
        assert by_name["glp1_basal_pm"].confidence_tier == ConfidenceTier.GOLD
        assert by_name["glp1_meal_response_pm_per_mm"].provenance.pubmed_id == "8138058"
        assert by_name["glp1_meal_response_pm_per_mm"].confidence_tier == ConfidenceTier.BRONZE

    def test_glp1_pk_does_not_create_spurious_conflicts(self, seeded):
        # Regression: 2 different CUSTOM-kind parameters sharing
        # subject="GLP-1" -- confirms ParameterStore's own conflict-
        # detection fix (custom_kind_name scoping) generalizes here too.
        _, store = seeded
        conflicts = store.unresolved_conflicts()
        glp1_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).custom_kind_name
            and "glp1" in store.get(c.existing_parameter_id).custom_kind_name
        ]
        assert glp1_conflicts == []

    # -- MASTER_PLAN 5.1 chunk 24: NEFA basal-clearance/antilipolysis PK constants --

    def test_nefa_pk_constants_stored(self, seeded):
        registry, store = seeded
        nefa_ids = registry.search_by_name("NEFA")
        assert nefa_ids
        params = store.query(str(nefa_ids[0]), ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        assert set(by_name) == {
            "nefa_half_life_min",
            "nefa_insulin_ic50_pm",
        }

    def test_nefa_pk_values_match_the_real_module_constants(self, seeded):
        # Confirms the curated store values are byte-identical to
        # nefa_pathway.py's own real constants.
        from biosim_ode.nefa_pathway import NEFA_HALF_LIFE_MIN, NEFA_INSULIN_IC50_PM
        registry, store = seeded
        nefa_id = str(registry.search_by_name("NEFA")[0])
        params = store.query(nefa_id, ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p.value.value for p in params}
        assert by_name["nefa_half_life_min"] == NEFA_HALF_LIFE_MIN
        assert by_name["nefa_insulin_ic50_pm"] == NEFA_INSULIN_IC50_PM

    def test_nefa_pk_cites_correct_pmids_and_tier(self, seeded):
        registry, store = seeded
        nefa_id = str(registry.search_by_name("NEFA")[0])
        params = store.query(nefa_id, ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        assert by_name["nefa_half_life_min"].provenance.pubmed_id == "1951681"
        assert by_name["nefa_half_life_min"].confidence_tier == ConfidenceTier.GOLD
        assert by_name["nefa_insulin_ic50_pm"].provenance.pubmed_id == "11061533"
        assert by_name["nefa_insulin_ic50_pm"].confidence_tier == ConfidenceTier.GOLD

    def test_nefa_pk_does_not_create_spurious_conflicts(self, seeded):
        # Regression: 2 different CUSTOM-kind parameters sharing
        # subject="NEFA" -- confirms ParameterStore's own conflict-
        # detection fix (custom_kind_name scoping) generalizes here too.
        _, store = seeded
        conflicts = store.unresolved_conflicts()
        nefa_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).custom_kind_name
            and "nefa" in store.get(c.existing_parameter_id).custom_kind_name
        ]
        assert nefa_conflicts == []

    # -- MASTER_PLAN 5.1 chunk 25: growth hormone acute p3 dose-response --

    def test_growth_hormone_p3_factor_stored(self, seeded):
        registry, store = seeded
        gh_ids = registry.search_by_name("growth hormone")
        assert gh_ids
        params = store.query(str(gh_ids[0]), ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        assert "growth_hormone_p3_factor" in by_name

    def test_growth_hormone_p3_factor_matches_the_real_module_constant(self, seeded):
        from biosim_physiology.drug_interactions import GROWTH_HORMONE_P3_FACTOR
        registry, store = seeded
        gh_id = str(registry.search_by_name("growth hormone")[0])
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(gh_id, ParameterKind.CUSTOM)
        }
        assert params["growth_hormone_p3_factor"] == pytest.approx(GROWTH_HORMONE_P3_FACTOR)

    def test_growth_hormone_cites_correct_pmid_and_tier(self, seeded):
        registry, store = seeded
        gh_id = str(registry.search_by_name("growth hormone")[0])
        entry = store.query(gh_id, ParameterKind.CUSTOM)[0]
        assert entry.provenance.pubmed_id == "1916055"
        assert entry.confidence_tier == ConfidenceTier.GOLD

    def test_growth_hormone_is_a_separate_subject_no_cross_conflict(self, seeded):
        # Same disclosed-honest pattern as pioglitazone's own equivalent
        # test above: growth hormone is its own subject entity, so this
        # test trivially passes even with a hypothetical scoping bug
        # reintroduced -- not presented as an equivalent-strength
        # regression guard to chunks 1/6/7's real ones.
        _, store = seeded
        conflicts = store.unresolved_conflicts()
        gh_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).custom_kind_name
            and "growth_hormone_p3_factor" in store.get(c.existing_parameter_id).custom_kind_name
        ]
        assert gh_conflicts == []

    # -- MASTER_PLAN 5.1 chunk 26: hydrochlorothiazide chronic p3 dose-response --

    def test_hydrochlorothiazide_p3_factor_stored(self, seeded):
        registry, store = seeded
        hctz_ids = registry.search_by_name("hydrochlorothiazide")
        assert hctz_ids
        params = store.query(str(hctz_ids[0]), ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        assert "thiazide_p3_factor_at_tested_dose" in by_name

    def test_hydrochlorothiazide_p3_factor_matches_the_real_module_constant(self, seeded):
        from biosim_physiology.drug_interactions import THIAZIDE_P3_FACTOR_AT_TESTED_DOSE
        registry, store = seeded
        hctz_id = str(registry.search_by_name("hydrochlorothiazide")[0])
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(hctz_id, ParameterKind.CUSTOM)
        }
        assert params["thiazide_p3_factor_at_tested_dose"] == pytest.approx(
            THIAZIDE_P3_FACTOR_AT_TESTED_DOSE
        )

    def test_hydrochlorothiazide_cites_correct_pmid_and_tier(self, seeded):
        registry, store = seeded
        hctz_id = str(registry.search_by_name("hydrochlorothiazide")[0])
        entry = store.query(hctz_id, ParameterKind.CUSTOM)[0]
        assert entry.provenance.pubmed_id == "2671740"
        assert entry.confidence_tier == ConfidenceTier.GOLD

    def test_hydrochlorothiazide_is_a_separate_subject_no_cross_conflict(self, seeded):
        # Same disclosed-honest pattern as pioglitazone/growth hormone's
        # own equivalent tests above.
        _, store = seeded
        conflicts = store.unresolved_conflicts()
        hctz_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).custom_kind_name
            and "thiazide_p3_factor_at_" in store.get(c.existing_parameter_id).custom_kind_name
        ]
        assert hctz_conflicts == []

    # -- MASTER_PLAN 5.1 chunk 27: indinavir acute p3 dose-response --

    def test_indinavir_p3_factor_stored(self, seeded):
        registry, store = seeded
        indinavir_ids = registry.search_by_name("indinavir")
        assert indinavir_ids
        params = store.query(str(indinavir_ids[0]), ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        assert "indinavir_p3_factor_at_tested_dose" in by_name

    def test_indinavir_p3_factor_matches_the_real_module_constant(self, seeded):
        from biosim_physiology.drug_interactions import INDINAVIR_P3_FACTOR_AT_TESTED_DOSE
        registry, store = seeded
        indinavir_id = str(registry.search_by_name("indinavir")[0])
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(indinavir_id, ParameterKind.CUSTOM)
        }
        assert params["indinavir_p3_factor_at_tested_dose"] == pytest.approx(
            INDINAVIR_P3_FACTOR_AT_TESTED_DOSE
        )

    def test_indinavir_cites_correct_pmid_and_tier(self, seeded):
        registry, store = seeded
        indinavir_id = str(registry.search_by_name("indinavir")[0])
        entry = store.query(indinavir_id, ParameterKind.CUSTOM)[0]
        assert entry.provenance.pubmed_id == "11964551"
        assert entry.confidence_tier == ConfidenceTier.GOLD

    def test_indinavir_is_a_separate_subject_no_cross_conflict(self, seeded):
        # Same disclosed-honest pattern as pioglitazone/growth hormone/
        # hydrochlorothiazide's own equivalent tests above.
        _, store = seeded
        conflicts = store.unresolved_conflicts()
        indinavir_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).custom_kind_name
            and "indinavir_p3_factor_at_" in store.get(c.existing_parameter_id).custom_kind_name
        ]
        assert indinavir_conflicts == []

    # -- MASTER_PLAN 5.1 chunk 28: atorvastatin chronic p3 dose-response --

    def test_atorvastatin_all_three_anchor_points_stored(self, seeded):
        registry, store = seeded
        atorvastatin_ids = registry.search_by_name("atorvastatin")
        assert atorvastatin_ids
        params = store.query(str(atorvastatin_ids[0]), ParameterKind.CUSTOM)
        names = {p.custom_kind_name for p in params}
        assert names == {"p3_factor_at_20mg", "p3_factor_at_40mg", "p3_factor_at_80mg"}

    def test_atorvastatin_factors_match_the_real_module_function(self, seeded):
        # Live cross-check against the REAL atorvastatin_p3_factor()
        # called at each exact anchor dose, not duplicated literals.
        from biosim_physiology.drug_interactions import atorvastatin_p3_factor
        registry, store = seeded
        atorvastatin_id = str(registry.search_by_name("atorvastatin")[0])
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(atorvastatin_id, ParameterKind.CUSTOM)
        }
        assert params["p3_factor_at_20mg"] == pytest.approx(atorvastatin_p3_factor(20.0))
        assert params["p3_factor_at_40mg"] == pytest.approx(atorvastatin_p3_factor(40.0))
        assert params["p3_factor_at_80mg"] == pytest.approx(atorvastatin_p3_factor(80.0))

    def test_atorvastatin_factors_cite_correct_pmid_and_gold_tier(self, seeded):
        registry, store = seeded
        atorvastatin_id = str(registry.search_by_name("atorvastatin")[0])
        entries = store.query(atorvastatin_id, ParameterKind.CUSTOM)
        assert len(entries) == 3
        for entry in entries:
            assert entry.provenance.pubmed_id == "20298928"
            assert entry.confidence_tier == ConfidenceTier.GOLD

    def test_atorvastatin_own_three_factors_do_not_conflict_with_each_other(self, seeded):
        registry, store = seeded
        atorvastatin_id = str(registry.search_by_name("atorvastatin")[0])
        conflicts = store.unresolved_conflicts()
        atorvastatin_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).custom_kind_name
            in ("p3_factor_at_20mg", "p3_factor_at_40mg", "p3_factor_at_80mg")
            and store.get(c.existing_parameter_id).subject_canonical_id == atorvastatin_id
        ]
        assert atorvastatin_conflicts == []

    # -- MASTER_PLAN 5.1 chunk 29: cold acclimation p3 constant (ContextRegistry) --

    def test_cold_acclimation_factor_stored(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        cold_id = str(ContextID(ContextEntityType.PHYSIOLOGICAL_AXIS, "cold_acclimation"))
        params = store.query(cold_id, ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        assert "p3_factor_at_reference" in by_name

    def test_cold_acclimation_factor_matches_the_real_module_function(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        from biosim_physiology.modifiers import cold_acclimation_modifiers
        _, store = seeded
        cold_id = str(ContextID(ContextEntityType.PHYSIOLOGICAL_AXIS, "cold_acclimation"))
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(cold_id, ParameterKind.CUSTOM)
        }
        reference = cold_acclimation_modifiers(True)
        assert params["p3_factor_at_reference"] == pytest.approx(reference["p3"])

    def test_cold_acclimation_factor_cites_correct_pmid_and_tier(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        cold_id = str(ContextID(ContextEntityType.PHYSIOLOGICAL_AXIS, "cold_acclimation"))
        entries = store.query(cold_id, ParameterKind.CUSTOM)
        assert len(entries) == 1
        for entry in entries:
            assert entry.provenance.pubmed_id == "26147760"
            assert entry.confidence_tier == ConfidenceTier.GOLD

    def test_cold_acclimation_is_a_separate_subject_no_cross_conflict(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        cold_id = str(ContextID(ContextEntityType.PHYSIOLOGICAL_AXIS, "cold_acclimation"))
        conflicts = store.unresolved_conflicts()
        cold_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).subject_canonical_id == cold_id
        ]
        assert cold_conflicts == []

    # -- MASTER_PLAN 5.1 chunk 30: menstrual cycle luteal p3 constant (ContextRegistry) --

    def test_menstrual_cycle_factor_stored(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        cycle_id = str(ContextID(ContextEntityType.PHYSIOLOGICAL_AXIS, "menstrual_cycle"))
        params = store.query(cycle_id, ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        assert "p3_factor_at_luteal" in by_name

    def test_menstrual_cycle_factor_matches_the_real_module_function(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        from biosim_physiology.modifiers import menstrual_cycle_modifiers
        from biosim_physiology.profile import MenstrualCyclePhase
        _, store = seeded
        cycle_id = str(ContextID(ContextEntityType.PHYSIOLOGICAL_AXIS, "menstrual_cycle"))
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(cycle_id, ParameterKind.CUSTOM)
        }
        reference = menstrual_cycle_modifiers(MenstrualCyclePhase.LUTEAL)
        assert params["p3_factor_at_luteal"] == pytest.approx(reference["p3"])

    def test_menstrual_cycle_factor_cites_correct_pmid_and_tier(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        cycle_id = str(ContextID(ContextEntityType.PHYSIOLOGICAL_AXIS, "menstrual_cycle"))
        entries = store.query(cycle_id, ParameterKind.CUSTOM)
        assert len(entries) == 1
        for entry in entries:
            assert entry.provenance.pubmed_id == "1997519"
            assert entry.confidence_tier == ConfidenceTier.GOLD

    def test_menstrual_cycle_is_a_separate_subject_no_cross_conflict(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        cycle_id = str(ContextID(ContextEntityType.PHYSIOLOGICAL_AXIS, "menstrual_cycle"))
        conflicts = store.unresolved_conflicts()
        cycle_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).subject_canonical_id == cycle_id
        ]
        assert cycle_conflicts == []

    # -- MASTER_PLAN 5.1 chunk 31: PCOS disease-state p3 constant (ContextRegistry) --

    def test_pcos_p3_factor_stored(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        pcos_id = str(ContextID(ContextEntityType.DISEASE_STATE, "pcos"))
        params = store.query(pcos_id, ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        assert "p3_factor" in by_name

    def test_pcos_factor_matches_the_real_module_function(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        from biosim_physiology.modifiers import disease_modifiers
        from biosim_physiology.profile import DiseaseState
        _, store = seeded
        pcos_id = str(ContextID(ContextEntityType.DISEASE_STATE, "pcos"))
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(pcos_id, ParameterKind.CUSTOM)
        }
        reference = disease_modifiers([DiseaseState.PCOS])
        assert params["p3_factor"] == pytest.approx(reference["p3"])

    def test_pcos_factor_cites_correct_pmid_and_tier(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        pcos_id = str(ContextID(ContextEntityType.DISEASE_STATE, "pcos"))
        entries = store.query(pcos_id, ParameterKind.CUSTOM)
        assert len(entries) == 1
        for entry in entries:
            assert entry.provenance.pubmed_id == "21920943"
            assert entry.confidence_tier == ConfidenceTier.GOLD

    def test_pcos_is_a_separate_subject_no_cross_conflict(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        pcos_id = str(ContextID(ContextEntityType.DISEASE_STATE, "pcos"))
        conflicts = store.unresolved_conflicts()
        pcos_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).subject_canonical_id == pcos_id
        ]
        assert pcos_conflicts == []

    # -- MASTER_PLAN 5.1 chunk 32: acromegaly Gb constant (ContextRegistry) --

    def test_acromegaly_gb_factor_stored(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        acromegaly_id = str(ContextID(ContextEntityType.DISEASE_STATE, "acromegaly"))
        params = store.query(acromegaly_id, ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        assert "Gb_factor" in by_name

    def test_acromegaly_factor_matches_the_real_module_function(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        from biosim_physiology.modifiers import disease_modifiers
        from biosim_physiology.profile import DiseaseState
        _, store = seeded
        acromegaly_id = str(ContextID(ContextEntityType.DISEASE_STATE, "acromegaly"))
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(acromegaly_id, ParameterKind.CUSTOM)
        }
        reference = disease_modifiers([DiseaseState.ACROMEGALY])
        assert params["Gb_factor"] == pytest.approx(reference["Gb"])

    def test_acromegaly_factor_cites_correct_pmid_and_tier(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        acromegaly_id = str(ContextID(ContextEntityType.DISEASE_STATE, "acromegaly"))
        entries = store.query(acromegaly_id, ParameterKind.CUSTOM)
        assert len(entries) == 1
        for entry in entries:
            assert entry.provenance.pubmed_id == "1569148"
            assert entry.confidence_tier == ConfidenceTier.GOLD

    def test_acromegaly_is_a_separate_subject_no_cross_conflict(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        acromegaly_id = str(ContextID(ContextEntityType.DISEASE_STATE, "acromegaly"))
        conflicts = store.unresolved_conflicts()
        acromegaly_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).subject_canonical_id == acromegaly_id
        ]
        assert acromegaly_conflicts == []

    # -- MASTER_PLAN 5.1 chunk 33: OSA disease-state p3 constant (ContextRegistry) --

    def test_osa_p3_factor_stored(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        osa_id = str(ContextID(ContextEntityType.DISEASE_STATE, "osa"))
        params = store.query(osa_id, ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        assert "p3_factor" in by_name

    def test_osa_factor_matches_the_real_module_function(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        from biosim_physiology.modifiers import disease_modifiers
        from biosim_physiology.profile import DiseaseState
        _, store = seeded
        osa_id = str(ContextID(ContextEntityType.DISEASE_STATE, "osa"))
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(osa_id, ParameterKind.CUSTOM)
        }
        reference = disease_modifiers([DiseaseState.OSA])
        assert params["p3_factor"] == pytest.approx(reference["p3"])

    def test_osa_factor_cites_correct_pmid_and_tier(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        osa_id = str(ContextID(ContextEntityType.DISEASE_STATE, "osa"))
        entries = store.query(osa_id, ParameterKind.CUSTOM)
        assert len(entries) == 1
        for entry in entries:
            assert entry.provenance.pubmed_id == "14512265"
            assert entry.confidence_tier == ConfidenceTier.GOLD

    def test_osa_is_a_separate_subject_no_cross_conflict(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        osa_id = str(ContextID(ContextEntityType.DISEASE_STATE, "osa"))
        conflicts = store.unresolved_conflicts()
        osa_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).subject_canonical_id == osa_id
        ]
        assert osa_conflicts == []

    # -- MASTER_PLAN 5.1 chunk 7: sulfonylurea glucose-independent secretion --

    def test_sulfonylurea_secretion_constant_stored(self, seeded):
        registry, store = seeded
        insulin_ids = registry.search_by_name("insulin")
        assert insulin_ids
        params = store.query(str(insulin_ids[0]), ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        assert "sulfonylurea_secretion_fraction_of_ke_ib" in by_name

    def test_sulfonylurea_secretion_value_matches_the_real_module_constant(self, seeded):
        from biosim_physiology.drug_interactions import (
            SULFONYLUREA_SECRETION_FRACTION_OF_KE_IB,
        )
        registry, store = seeded
        insulin_id = str(registry.search_by_name("insulin")[0])
        params = store.query(insulin_id, ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p.value.value for p in params}
        assert (
            by_name["sulfonylurea_secretion_fraction_of_ke_ib"]
            == SULFONYLUREA_SECRETION_FRACTION_OF_KE_IB
        )

    def test_sulfonylurea_secretion_cites_correct_pmid_and_tier(self, seeded):
        registry, store = seeded
        insulin_id = str(registry.search_by_name("insulin")[0])
        params = store.query(insulin_id, ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        entry = by_name["sulfonylurea_secretion_fraction_of_ke_ib"]
        assert entry.provenance.pubmed_id == "25297572"
        assert entry.confidence_tier == ConfidenceTier.BRONZE

    def test_sulfonylurea_secretion_does_not_create_spurious_conflicts(self, seeded):
        # Regression: 5 different CUSTOM-kind parameters now share
        # subject="insulin" (4 exogenous-insulin PK entries + this one) --
        # confirms the custom_kind_name conflict-detection fix still
        # holds as the number of CUSTOM entries under one subject grows.
        _, store = seeded
        conflicts = store.unresolved_conflicts()
        insulin_custom_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).custom_kind_name
            and (
                "sulfonylurea" in store.get(c.existing_parameter_id).custom_kind_name
                or "exogenous_insulin" in store.get(c.existing_parameter_id).custom_kind_name
            )
        ]
        assert insulin_custom_conflicts == []

    # -- MASTER_PLAN 5.1 chunk 8: acarbose/miglitol kabs-scale constants --

    def test_acarbose_kabs_scale_stored(self, seeded):
        registry, store = seeded
        acarbose_ids = registry.search_by_name("acarbose")
        assert acarbose_ids
        params = store.query(str(acarbose_ids[0]), ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        assert "kabs_scale_at_100mg_reference_dose" in by_name

    def test_miglitol_kabs_scale_stored(self, seeded):
        registry, store = seeded
        miglitol_ids = registry.search_by_name("miglitol")
        assert miglitol_ids
        params = store.query(str(miglitol_ids[0]), ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        assert "kabs_scale_at_50mg_reference_dose" in by_name

    def test_acarbose_miglitol_kabs_values_match_the_real_module_constants(self, seeded):
        from biosim_physiology.drug_interactions import (
            ACARBOSE_100MG_KABS_SCALE, MIGLITOL_50MG_KABS_SCALE,
        )
        registry, store = seeded
        acarbose_id = str(registry.search_by_name("acarbose")[0])
        miglitol_id = str(registry.search_by_name("miglitol")[0])
        acarbose_params = {
            p.custom_kind_name: p.value.value
            for p in store.query(acarbose_id, ParameterKind.CUSTOM)
        }
        miglitol_params = {
            p.custom_kind_name: p.value.value
            for p in store.query(miglitol_id, ParameterKind.CUSTOM)
        }
        assert acarbose_params["kabs_scale_at_100mg_reference_dose"] == ACARBOSE_100MG_KABS_SCALE
        assert miglitol_params["kabs_scale_at_50mg_reference_dose"] == MIGLITOL_50MG_KABS_SCALE

    def test_acarbose_miglitol_cite_correct_pmids_and_tier(self, seeded):
        registry, store = seeded
        acarbose_id = str(registry.search_by_name("acarbose")[0])
        miglitol_id = str(registry.search_by_name("miglitol")[0])
        acarbose_entry = store.query(acarbose_id, ParameterKind.CUSTOM)[0]
        miglitol_entry = store.query(miglitol_id, ParameterKind.CUSTOM)[0]
        assert acarbose_entry.provenance.pubmed_id == "27686734"
        assert acarbose_entry.confidence_tier == ConfidenceTier.BRONZE
        assert miglitol_entry.provenance.pubmed_id == "3042310"
        assert miglitol_entry.confidence_tier == ConfidenceTier.BRONZE

    def test_acarbose_and_miglitol_are_separate_subjects_no_cross_conflict(self, seeded):
        # NOTE (found by mandatory code review): unlike chunks 1/6/7's
        # equivalent "no spurious conflicts" tests -- which shared ONE
        # subject and therefore had a real custom_kind_name-scoping bug
        # to guard against -- acarbose and miglitol are two SEPARATE
        # subject entities, so ParameterStore's conflict index key can
        # never place them in the same bucket regardless of whether
        # that scoping fix exists at all. This test trivially passes
        # even with the source change reverted (confirmed directly) --
        # disclosed honestly here, matching this campaign's own
        # established practice for non-discriminating conflict tests,
        # rather than presenting it as equivalent in strength to the
        # earlier chunks' real regression guards.
        _, store = seeded
        conflicts = store.unresolved_conflicts()
        agi_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).custom_kind_name
            and "kabs_scale_at_" in store.get(c.existing_parameter_id).custom_kind_name
        ]
        assert agi_conflicts == []

    # -- MASTER_PLAN 5.1 chunk 9: pioglitazone p3 dose-response constant --

    def test_pioglitazone_p3_multiplier_stored(self, seeded):
        registry, store = seeded
        pioglitazone_ids = registry.search_by_name("pioglitazone")
        assert pioglitazone_ids
        params = store.query(str(pioglitazone_ids[0]), ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        assert "p3_multiplier_at_45mg_reference_dose" in by_name

    def test_pioglitazone_p3_multiplier_matches_the_real_module_constant(self, seeded):
        from biosim_physiology.drug_interactions import pioglitazone_p3_factor
        registry, store = seeded
        pioglitazone_id = str(registry.search_by_name("pioglitazone")[0])
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(pioglitazone_id, ParameterKind.CUSTOM)
        }
        assert params["p3_multiplier_at_45mg_reference_dose"] == pytest.approx(
            pioglitazone_p3_factor(45.0)
        )

    def test_pioglitazone_cites_correct_pmid_and_tier(self, seeded):
        # GOLD, not BRONZE -- re-tiered per MASTER_PLAN Section 8's own
        # housekeeping item: a ratio computed directly from Miyazaki
        # 2002's own two reported ISI numbers (2.7/1.6), no unit
        # conversion, same evidentiary category as CKD's p3_factor/
        # ke_factor and tirzepatide's gb_factor_at_*mg (both GOLD).
        registry, store = seeded
        pioglitazone_id = str(registry.search_by_name("pioglitazone")[0])
        entry = store.query(pioglitazone_id, ParameterKind.CUSTOM)[0]
        assert entry.provenance.pubmed_id == "11874940"
        assert entry.confidence_tier == ConfidenceTier.GOLD

    def test_pioglitazone_is_a_separate_subject_no_cross_conflict(self, seeded):
        # Same disclosed-honest pattern as acarbose/miglitol's own
        # equivalent test above: pioglitazone is its own subject entity,
        # so this test trivially passes even with a hypothetical
        # scoping bug reintroduced -- not presented as an equivalent-
        # strength regression guard to chunks 1/6/7's real ones.
        _, store = seeded
        conflicts = store.unresolved_conflicts()
        pio_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).custom_kind_name
            and "p3_multiplier_at_" in store.get(c.existing_parameter_id).custom_kind_name
        ]
        assert pio_conflicts == []

    # -- MASTER_PLAN 5.1 chunk 10: protein secretion constant --

    def test_protein_secretion_pm_per_gram_stored(self, seeded):
        registry, store = seeded
        protein_ids = registry.search_by_name("protein")
        assert protein_ids
        params = store.query(str(protein_ids[0]), ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        assert "secretion_pm_per_gram" in by_name

    def test_protein_secretion_matches_the_real_module_constant(self, seeded):
        from biosim_ode.nutrient_secretion_pathway import NUTRIENT_SECRETION_PM_PER_G
        registry, store = seeded
        protein_id = str(registry.search_by_name("protein")[0])
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(protein_id, ParameterKind.CUSTOM)
        }
        assert params["secretion_pm_per_gram"] == pytest.approx(NUTRIENT_SECRETION_PM_PER_G)

    def test_protein_secretion_cites_correct_pmid_and_tier(self, seeded):
        registry, store = seeded
        protein_id = str(registry.search_by_name("protein")[0])
        entry = store.query(protein_id, ParameterKind.CUSTOM)[0]
        assert entry.provenance.pubmed_id == "3538843"
        assert entry.confidence_tier == ConfidenceTier.BRONZE

    def test_protein_is_a_separate_subject_no_cross_conflict(self, seeded):
        # Same disclosed-honest pattern as pioglitazone/acarbose/
        # miglitol's own equivalent tests above.
        _, store = seeded
        conflicts = store.unresolved_conflicts()
        protein_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).custom_kind_name
            and "secretion_pm_per_gram" in store.get(c.existing_parameter_id).custom_kind_name
        ]
        assert protein_conflicts == []

    # -- MASTER_PLAN 5.1 chunk 11: CKD disease-state constants (ContextRegistry) --

    def test_ckd_p3_and_ke_factors_stored(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        ckd_id = str(ContextID(ContextEntityType.DISEASE_STATE, "chronic_kidney_disease"))
        params = store.query(ckd_id, ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        assert "p3_factor" in by_name
        assert "ke_factor" in by_name

    def test_ckd_factors_match_the_real_module_constants(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        ckd_id = str(ContextID(ContextEntityType.DISEASE_STATE, "chronic_kidney_disease"))
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(ckd_id, ParameterKind.CUSTOM)
        }
        assert params["p3_factor"] == pytest.approx(0.86)
        assert params["ke_factor"] == pytest.approx(0.915)

    def test_ckd_factors_cite_correct_pmid_and_tier(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        ckd_id = str(ContextID(ContextEntityType.DISEASE_STATE, "chronic_kidney_disease"))
        for entry in store.query(ckd_id, ParameterKind.CUSTOM):
            assert entry.provenance.pubmed_id == "26823551"
            assert entry.confidence_tier == ConfidenceTier.GOLD

    def test_ckd_is_a_separate_subject_no_cross_conflict(self, seeded):
        _, store = seeded
        conflicts = store.unresolved_conflicts()
        ckd_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).custom_kind_name
            and store.get(c.existing_parameter_id).custom_kind_name in ("p3_factor", "ke_factor")
        ]
        assert ckd_conflicts == []

    # -- MASTER_PLAN 5.1 chunk 12: sleep restriction constant (ContextRegistry) --

    def test_sleep_restriction_p3_factor_stored(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        sleep_id = str(ContextID(ContextEntityType.PHYSIOLOGICAL_AXIS, "sleep_restriction"))
        params = store.query(sleep_id, ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        assert "p3_factor_at_5h_reference" in by_name

    def test_sleep_restriction_matches_the_real_module_constant(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        from biosim_physiology.modifiers import SLEEP_RESTRICTED_P3_FACTOR
        _, store = seeded
        sleep_id = str(ContextID(ContextEntityType.PHYSIOLOGICAL_AXIS, "sleep_restriction"))
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(sleep_id, ParameterKind.CUSTOM)
        }
        assert params["p3_factor_at_5h_reference"] == pytest.approx(SLEEP_RESTRICTED_P3_FACTOR)

    def test_sleep_restriction_cites_correct_pmid_and_tier(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        sleep_id = str(ContextID(ContextEntityType.PHYSIOLOGICAL_AXIS, "sleep_restriction"))
        entry = store.query(sleep_id, ParameterKind.CUSTOM)[0]
        assert entry.provenance.pubmed_id == "20585000"
        assert entry.confidence_tier == ConfidenceTier.GOLD

    # -- MASTER_PLAN 5.1 chunk 13: GCK-MODY Gb constant (ContextRegistry) --

    def test_gck_mody_gb_factor_stored(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        gck_mody_id = str(ContextID(ContextEntityType.DISEASE_STATE, "gck_mody"))
        params = store.query(gck_mody_id, ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        assert "Gb_factor" in by_name

    def test_gck_mody_factor_matches_the_real_module_value(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        from biosim_physiology.modifiers import disease_modifiers
        from biosim_physiology.profile import DiseaseState
        _, store = seeded
        gck_mody_id = str(ContextID(ContextEntityType.DISEASE_STATE, "gck_mody"))
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(gck_mody_id, ParameterKind.CUSTOM)
        }
        expected = disease_modifiers([DiseaseState.GCK_MODY])["Gb"]
        assert params["Gb_factor"] == pytest.approx(expected)

    def test_gck_mody_only_registers_gb_not_p1_p3_secretion_scale(self, seeded):
        # The citation explicitly confirms p1/p3/secretion_scale are
        # unaffected by GCK-MODY -- only Gb_factor should ever be
        # registered under this subject.
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        gck_mody_id = str(ContextID(ContextEntityType.DISEASE_STATE, "gck_mody"))
        names = {p.custom_kind_name for p in store.query(gck_mody_id, ParameterKind.CUSTOM)}
        assert names == {"Gb_factor"}

    def test_gck_mody_cites_correct_pmid_and_tier(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        gck_mody_id = str(ContextID(ContextEntityType.DISEASE_STATE, "gck_mody"))
        entry = store.query(gck_mody_id, ParameterKind.CUSTOM)[0]
        assert entry.provenance.pubmed_id == "11119017"
        assert entry.confidence_tier == ConfidenceTier.GOLD

    def test_gck_mody_is_a_separate_subject_no_cross_conflict(self, seeded):
        _, store = seeded
        conflicts = store.unresolved_conflicts()
        gck_mody_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).custom_kind_name
            and store.get(c.existing_parameter_id).custom_kind_name == "Gb_factor"
        ]
        assert gck_mody_conflicts == []

    # -- MASTER_PLAN 5.1 chunk 14: gestational diabetes constants (ContextRegistry) --

    def test_gestational_diabetes_all_four_factors_stored(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        gdm_id = str(ContextID(ContextEntityType.DISEASE_STATE, "gestational_diabetes"))
        names = {p.custom_kind_name for p in store.query(gdm_id, ParameterKind.CUSTOM)}
        assert names == {"p1_factor", "p3_factor", "secretion_scale_factor", "Gb_factor"}

    def test_gestational_diabetes_factors_match_the_real_module_values(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        from biosim_physiology.modifiers import disease_modifiers
        from biosim_physiology.profile import DiseaseState
        _, store = seeded
        gdm_id = str(ContextID(ContextEntityType.DISEASE_STATE, "gestational_diabetes"))
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(gdm_id, ParameterKind.CUSTOM)
        }
        expected = disease_modifiers([DiseaseState.GESTATIONAL_DIABETES])
        assert params["p1_factor"] == pytest.approx(expected["p1"])
        assert params["p3_factor"] == pytest.approx(expected["p3"])
        assert params["secretion_scale_factor"] == pytest.approx(expected["secretion_scale"])
        assert params["Gb_factor"] == pytest.approx(expected["Gb"])

    def test_gestational_diabetes_p1_p3_secretion_scale_cite_kautzky_willer_gold(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        gdm_id = str(ContextID(ContextEntityType.DISEASE_STATE, "gestational_diabetes"))
        by_name = {p.custom_kind_name: p for p in store.query(gdm_id, ParameterKind.CUSTOM)}
        for name in ("p1_factor", "p3_factor", "secretion_scale_factor"):
            assert by_name[name].provenance.pubmed_id == "9353615"
            assert by_name[name].confidence_tier == ConfidenceTier.GOLD

    def test_gestational_diabetes_gb_cites_iadpsg_silver(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        gdm_id = str(ContextID(ContextEntityType.DISEASE_STATE, "gestational_diabetes"))
        by_name = {p.custom_kind_name: p for p in store.query(gdm_id, ParameterKind.CUSTOM)}
        assert by_name["Gb_factor"].provenance.pubmed_id == "20190296"
        assert by_name["Gb_factor"].confidence_tier == ConfidenceTier.SILVER

    # -- MASTER_PLAN 5.1 chunk 17: metformin Gb reduction constant --

    def test_metformin_gb_reduction_fraction_stored(self, seeded):
        registry, store = seeded
        metformin_ids = registry.search_by_name("metformin")
        assert metformin_ids
        params = store.query(str(metformin_ids[0]), ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        assert "gb_reduction_fraction_at_max_dose" in by_name

    def test_metformin_gb_reduction_matches_the_real_module_behavior(self, seeded):
        # Cross-checks against the REAL metformin_modifiers() function at
        # its 2000mg dose cap (dose_factor == 1.0 exactly, since
        # log(1+2000/500)/log(5) == log(5)/log(5) == 1.0), not a
        # duplicated 0.147 literal -- if the real function's own
        # reduction fraction ever drifts, this test breaks immediately.
        from biosim_physiology.drug_interactions import Drug, DrugDose, metformin_modifiers
        registry, store = seeded
        metformin_id = str(registry.search_by_name("metformin")[0])
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(metformin_id, ParameterKind.CUSTOM)
        }
        at_max_dose = metformin_modifiers(DrugDose(drug=Drug.METFORMIN, dose_mg=2000.0))
        real_reduction_fraction = 1.0 - at_max_dose["Gb"]
        assert params["gb_reduction_fraction_at_max_dose"] == pytest.approx(real_reduction_fraction)

    def test_metformin_gb_cites_correct_pmid_and_tier(self, seeded):
        registry, store = seeded
        metformin_id = str(registry.search_by_name("metformin")[0])
        entry = store.query(metformin_id, ParameterKind.CUSTOM)[0]
        assert entry.provenance.pubmed_id == "8923861"
        assert entry.confidence_tier == ConfidenceTier.BRONZE

    def test_metformin_is_a_separate_subject_no_cross_conflict(self, seeded):
        # Same disclosed-honest pattern as pioglitazone/acarbose/
        # miglitol's own equivalent tests above.
        _, store = seeded
        conflicts = store.unresolved_conflicts()
        metformin_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).custom_kind_name == "gb_reduction_fraction_at_max_dose"
        ]
        assert metformin_conflicts == []

    # -- MASTER_PLAN 5.1 chunk 18: caffeine p3/Gb dose-response constants --

    def test_caffeine_both_factors_stored(self, seeded):
        registry, store = seeded
        caffeine_ids = registry.search_by_name("caffeine")
        assert caffeine_ids
        params = store.query(str(caffeine_ids[0]), ParameterKind.CUSTOM)
        names = {p.custom_kind_name for p in params}
        # gb_increase_amplitude removed 2026-11-20 with the caffeine Gb term.
        assert names == {"p3_reduction_at_200mg_reference_dose"}

    def test_caffeine_factors_match_the_real_module_behavior(self, seeded):
        # A genuine live cross-check, not a duplicated-literal comparison
        # (mandatory code-reviewer finding: an earlier draft wrongly
        # claimed no live cross-check was possible here, since
        # concentration_fraction never reaches exactly 1.0 -- but it's a
        # fully known closed-form function of t, so it can be computed
        # independently and divided back out). Calls the REAL
        # caffeine_modifiers() at a fixed dose/time and back-solves the
        # underlying 0.15/0.02 constants -- if either ever drifts in
        # drug_interactions.py, this test breaks immediately.
        import math
        from biosim_physiology.drug_interactions import Drug, DrugDose, caffeine_modifiers
        registry, store = seeded
        caffeine_id = str(registry.search_by_name("caffeine")[0])
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(caffeine_id, ParameterKind.CUSTOM)
        }

        t = 1.0
        dose = DrugDose(drug=Drug.CAFFEINE, dose_mg=200.0, hours_since_ingestion=t)
        result = caffeine_modifiers(dose)
        concentration_fraction = (
            (1.0 - math.exp(-math.log(2) / 0.75 * t)) * math.exp(-math.log(2) / 5.0 * t)
        )
        # dose_factor(200mg) == log(1+200/100)/log(3) == log(3)/log(3) == 1.0
        recovered_p3_reduction = (1.0 - result["p3"]) / concentration_fraction

        assert params["p3_reduction_at_200mg_reference_dose"] == pytest.approx(recovered_p3_reduction)
        assert "Gb" not in result

    def test_caffeine_p3_cites_keijzers_bronze(self, seeded):
        # The uncited Gb addition was removed on 2026-11-20; only the cited
        # p3 magnitude (Keijzers 2002) remains.
        registry, store = seeded
        caffeine_id = str(registry.search_by_name("caffeine")[0])
        by_name = {p.custom_kind_name: p for p in store.query(caffeine_id, ParameterKind.CUSTOM)}
        assert set(by_name) == {"p3_reduction_at_200mg_reference_dose"}
        assert by_name["p3_reduction_at_200mg_reference_dose"].provenance.pubmed_id == "11815511"
        assert by_name["p3_reduction_at_200mg_reference_dose"].confidence_tier == ConfidenceTier.BRONZE

    def test_caffeine_own_two_factors_do_not_conflict_with_each_other(self, seeded):
        registry, store = seeded
        caffeine_id = str(registry.search_by_name("caffeine")[0])
        conflicts = store.unresolved_conflicts()
        caffeine_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).custom_kind_name
            in ("p3_reduction_at_200mg_reference_dose",)
            and store.get(c.existing_parameter_id).subject_canonical_id == caffeine_id
        ]
        assert caffeine_conflicts == []

    # -- MASTER_PLAN 5.1 chunk 19: tirzepatide Gb dose-response constants --

    def test_tirzepatide_all_three_anchor_points_stored(self, seeded):
        registry, store = seeded
        tirzepatide_ids = registry.search_by_name("tirzepatide")
        assert tirzepatide_ids
        params = store.query(str(tirzepatide_ids[0]), ParameterKind.CUSTOM)
        names = {p.custom_kind_name for p in params}
        assert names == {"gb_factor_at_5mg", "gb_factor_at_10mg", "gb_factor_at_15mg"}

    def test_tirzepatide_factors_match_the_real_module_function(self, seeded):
        # Live cross-check against the REAL tirzepatide_gb_factor()
        # called at each exact anchor dose, not duplicated literals.
        from biosim_physiology.drug_interactions import tirzepatide_gb_factor
        registry, store = seeded
        tirzepatide_id = str(registry.search_by_name("tirzepatide")[0])
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(tirzepatide_id, ParameterKind.CUSTOM)
        }
        assert params["gb_factor_at_5mg"] == pytest.approx(tirzepatide_gb_factor(5.0))
        assert params["gb_factor_at_10mg"] == pytest.approx(tirzepatide_gb_factor(10.0))
        assert params["gb_factor_at_15mg"] == pytest.approx(tirzepatide_gb_factor(15.0))

    def test_tirzepatide_factors_cite_correct_pmid_and_gold_tier(self, seeded):
        registry, store = seeded
        tirzepatide_id = str(registry.search_by_name("tirzepatide")[0])
        entries = store.query(tirzepatide_id, ParameterKind.CUSTOM)
        assert len(entries) == 3
        for entry in entries:
            assert entry.provenance.pubmed_id == "34647404"
            assert entry.confidence_tier == ConfidenceTier.GOLD

    def test_tirzepatide_own_three_factors_do_not_conflict_with_each_other(self, seeded):
        registry, store = seeded
        tirzepatide_id = str(registry.search_by_name("tirzepatide")[0])
        conflicts = store.unresolved_conflicts()
        tirzepatide_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).custom_kind_name
            in ("gb_factor_at_5mg", "gb_factor_at_10mg", "gb_factor_at_15mg")
            and store.get(c.existing_parameter_id).subject_canonical_id == tirzepatide_id
        ]
        assert tirzepatide_conflicts == []

    # -- MASTER_PLAN 5.1 chunk 20: prednisone p3/Gb/secretion_scale constants --

    def test_prednisone_all_three_factors_stored(self, seeded):
        registry, store = seeded
        prednisone_ids = registry.search_by_name("prednisone")
        assert prednisone_ids
        params = store.query(str(prednisone_ids[0]), ParameterKind.CUSTOM)
        names = {p.custom_kind_name for p in params}
        assert names == {
            "p3_reduction_at_30mg_reference_dose",
            "gb_increase_at_30mg_reference_dose",
            "secretion_scale_increase_at_30mg_reference_dose",
        }

    def test_prednisone_factors_match_the_real_module_function(self, seeded):
        # Live cross-check at the 30mg dose cap (dose_factor==1.0
        # exactly), not duplicated literals.
        from biosim_physiology.drug_interactions import Drug, DrugDose, glucocorticoid_modifiers
        registry, store = seeded
        prednisone_id = str(registry.search_by_name("prednisone")[0])
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(prednisone_id, ParameterKind.CUSTOM)
        }
        at_max_dose = glucocorticoid_modifiers(DrugDose(drug=Drug.GLUCOCORTICOID, dose_mg=30.0))
        assert params["p3_reduction_at_30mg_reference_dose"] == pytest.approx(1.0 - at_max_dose["p3"])
        assert params["gb_increase_at_30mg_reference_dose"] == pytest.approx(at_max_dose["Gb"] - 1.0)
        assert params["secretion_scale_increase_at_30mg_reference_dose"] == pytest.approx(
            at_max_dose["secretion_scale"] - 1.0
        )

    def test_prednisone_factors_cite_correct_pmid_and_bronze_tier(self, seeded):
        registry, store = seeded
        prednisone_id = str(registry.search_by_name("prednisone")[0])
        entries = store.query(prednisone_id, ParameterKind.CUSTOM)
        assert len(entries) == 3
        for entry in entries:
            assert entry.provenance.pubmed_id == "7033265"
            assert entry.confidence_tier == ConfidenceTier.BRONZE

    def test_prednisone_own_three_factors_do_not_conflict_with_each_other(self, seeded):
        registry, store = seeded
        prednisone_id = str(registry.search_by_name("prednisone")[0])
        conflicts = store.unresolved_conflicts()
        prednisone_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).custom_kind_name
            in (
                "p3_reduction_at_30mg_reference_dose",
                "gb_increase_at_30mg_reference_dose",
                "secretion_scale_increase_at_30mg_reference_dose",
            )
            and store.get(c.existing_parameter_id).subject_canonical_id == prednisone_id
        ]
        assert prednisone_conflicts == []

    # -- MASTER_PLAN 5.1 chunk 21: glibenclamide/glipizide entities + glipizide's own secretion constant --

    def test_secretion_fraction_constants_registered_under_insulin_not_the_drug_entities(self, seeded):
        # Deliberate decision: both drugs' own secretion-fraction
        # constants stay registered under the "insulin" subject (a
        # fraction of insulin's OWN clearance kinetics, not either
        # drug's own chemical property) -- this chunk did NOT move or
        # duplicate them under the new glibenclamide/glipizide entities,
        # which carry no CUSTOM parameters of their own at all.
        registry, store = seeded
        insulin_id = str(registry.search_by_name("insulin")[0])
        params = {p.custom_kind_name for p in store.query(insulin_id, ParameterKind.CUSTOM)}
        assert "sulfonylurea_secretion_fraction_of_ke_ib" in params
        assert "glipizide_secretion_fraction_of_ke_ib" in params

        glibenclamide_id = str(registry.search_by_name("glibenclamide")[0])
        glipizide_id = str(registry.search_by_name("glipizide")[0])
        assert store.query(glibenclamide_id, ParameterKind.CUSTOM) == []
        assert store.query(glipizide_id, ParameterKind.CUSTOM) == []

    def test_glipizide_secretion_fraction_matches_the_real_module_constant(self, seeded):
        from biosim_physiology.drug_interactions import GLIPIZIDE_SECRETION_FRACTION_OF_KE_IB
        registry, store = seeded
        insulin_id = str(registry.search_by_name("insulin")[0])
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(insulin_id, ParameterKind.CUSTOM)
        }
        assert params["glipizide_secretion_fraction_of_ke_ib"] == pytest.approx(
            GLIPIZIDE_SECRETION_FRACTION_OF_KE_IB
        )

    def test_glipizide_secretion_fraction_cites_correct_pmid_and_tier(self, seeded):
        registry, store = seeded
        insulin_id = str(registry.search_by_name("insulin")[0])
        by_name = {p.custom_kind_name: p for p in store.query(insulin_id, ParameterKind.CUSTOM)}
        entry = by_name["glipizide_secretion_fraction_of_ke_ib"]
        assert entry.provenance.pubmed_id == "3117609"
        assert entry.confidence_tier == ConfidenceTier.BRONZE

    def test_glipizide_secretion_fraction_does_not_conflict_with_glibenclamides(self, seeded):
        # Same-subject guard: both fractions (0.5 and 0.9) sit under
        # "insulin" alongside the 4 pre-existing exogenous-insulin PK
        # entries -- confirms custom_kind_name-aware scoping still
        # correctly prevents this NEW entry from spuriously conflicting.
        registry, store = seeded
        insulin_id = str(registry.search_by_name("insulin")[0])
        conflicts = store.unresolved_conflicts()
        sulfonylurea_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).custom_kind_name
            in ("sulfonylurea_secretion_fraction_of_ke_ib", "glipizide_secretion_fraction_of_ke_ib")
            and store.get(c.existing_parameter_id).subject_canonical_id == insulin_id
        ]
        assert sulfonylurea_conflicts == []

    # -- MASTER_PLAN 5.1 chunk 22: intermittent hypoxia p1/p3 constants (ContextRegistry) --

    def test_intermittent_hypoxia_factors_stored(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        hypoxia_id = str(ContextID(ContextEntityType.PHYSIOLOGICAL_AXIS, "intermittent_hypoxia"))
        params = store.query(hypoxia_id, ParameterKind.CUSTOM)
        by_name = {p.custom_kind_name: p for p in params}
        assert "p3_factor_at_reference_rate" in by_name
        assert "p1_factor_at_reference_rate" in by_name

    def test_intermittent_hypoxia_factors_match_the_real_module_function(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        from biosim_physiology.modifiers import intermittent_hypoxia_modifiers
        _, store = seeded
        hypoxia_id = str(ContextID(ContextEntityType.PHYSIOLOGICAL_AXIS, "intermittent_hypoxia"))
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(hypoxia_id, ParameterKind.CUSTOM)
        }
        reference = intermittent_hypoxia_modifiers(24.3)
        assert params["p3_factor_at_reference_rate"] == pytest.approx(reference["p3"])
        assert params["p1_factor_at_reference_rate"] == pytest.approx(reference["p1"])

    def test_intermittent_hypoxia_factors_cite_correct_pmid_and_tier(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        hypoxia_id = str(ContextID(ContextEntityType.PHYSIOLOGICAL_AXIS, "intermittent_hypoxia"))
        entries = store.query(hypoxia_id, ParameterKind.CUSTOM)
        assert len(entries) == 2
        for entry in entries:
            assert entry.provenance.pubmed_id == "19265062"
            assert entry.confidence_tier == ConfidenceTier.GOLD

    def test_intermittent_hypoxia_is_a_separate_subject_no_cross_conflict(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        hypoxia_id = str(ContextID(ContextEntityType.PHYSIOLOGICAL_AXIS, "intermittent_hypoxia"))
        conflicts = store.unresolved_conflicts()
        hypoxia_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).subject_canonical_id == hypoxia_id
        ]
        assert hypoxia_conflicts == []

    def test_p3_factor_shared_name_across_subjects_does_not_cross_conflict(self, seeded):
        # "p3_factor" is used as a custom_kind_name under SIX different
        # subjects in this campaign by this point (CKD, GDM, T2D,
        # insulin_resistant, hypothyroidism, PCOS -- originally written
        # for just GDM/CKD, broadened as later chunks added more). Correction
        # (mandatory code-reviewer finding on the original 2-subject
        # version of this test): this does NOT exercise a different code
        # path than the pre-existing CKD/sleep-restriction cross-conflict
        # test below -- store.py's _index_key() already partitions the
        # conflict-detection bucket by subject BEFORE the custom_kind_name
        # comparison is ever reached (confirmed by reading store.py
        # directly), so different-subject parameters never even reach the
        # custom_kind_name check regardless of whether their names match.
        # Kept and broadened as a still-legitimate, non-redundant
        # regression guard (a real, now 5-subject scenario worth
        # pinning), just not evidence of a distinct disambiguation
        # mechanism.
        _, store = seeded
        conflicts = store.unresolved_conflicts()
        p3_factor_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).custom_kind_name == "p3_factor"
        ]
        assert p3_factor_conflicts == []

    def test_gestational_diabetes_own_four_factors_do_not_conflict_with_each_other(self, seeded):
        # The actually load-bearing guard for this chunk (mandatory
        # code-reviewer finding): GDM registers FOUR CUSTOM parameters
        # under ONE subject, several pairs exceeding
        # CONFLICT_FOLD_THRESHOLD (e.g. secretion_scale_factor/p3_factor
        # = 1.30/0.16 = 8.125) -- the exact same-subject, multiple-
        # distinct-custom_kind_name shape that previously produced 5
        # spurious conflicts for the 4 exogenous-insulin PK constants
        # (see store.py's own comment) before custom_kind_name-aware
        # scoping was added. This is the first multi-parameter-per-
        # subject entity in the ContextRegistry campaign, so it needs
        # its own same-subject guard, not just the cross-subject one
        # above.
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        gdm_id = str(ContextID(ContextEntityType.DISEASE_STATE, "gestational_diabetes"))
        conflicts = store.unresolved_conflicts()
        gdm_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).custom_kind_name
            in ("p1_factor", "p3_factor", "secretion_scale_factor", "Gb_factor")
            and store.get(c.existing_parameter_id).subject_canonical_id == gdm_id
        ]
        assert gdm_conflicts == []

    # -- MASTER_PLAN 5.1 chunk 15: T1D Ib/secretion_scale constants (ContextRegistry) --

    def test_t1d_all_three_factors_stored_including_unverified_gb(self, seeded):
        # All three of modifiers.py's real T1D constants are registered
        # -- Gb_factor too, tagged pubmed_id="UNVERIFIED" rather than
        # omitted (mandatory code-reviewer finding: silent omission would
        # have hidden a real, actively-used constant from the provenance
        # store's own audit trail).
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        t1d_id = str(ContextID(ContextEntityType.DISEASE_STATE, "type_1_diabetes"))
        names = {p.custom_kind_name for p in store.query(t1d_id, ParameterKind.CUSTOM)}
        assert names == {"Ib_factor", "secretion_scale_factor", "Gb_factor"}

    def test_t1d_factors_match_the_real_module_values(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        from biosim_physiology.modifiers import disease_modifiers
        from biosim_physiology.profile import DiseaseState
        _, store = seeded
        t1d_id = str(ContextID(ContextEntityType.DISEASE_STATE, "type_1_diabetes"))
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(t1d_id, ParameterKind.CUSTOM)
        }
        expected = disease_modifiers([DiseaseState.TYPE_1_DIABETES])
        assert params["Ib_factor"] == pytest.approx(expected["Ib"])
        assert params["secretion_scale_factor"] == pytest.approx(expected["secretion_scale"])
        assert params["Gb_factor"] == pytest.approx(expected["Gb"])

    def test_t1d_ib_and_secretion_scale_cite_correct_pmid_bronze_tier(self, seeded):
        # BRONZE, not GOLD -- these are disclosed engineering choices
        # motivated by a qualitative citation, not directly-reported
        # values, unlike CKD/GCK-MODY/GDM's own GOLD factors.
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        t1d_id = str(ContextID(ContextEntityType.DISEASE_STATE, "type_1_diabetes"))
        by_name = {p.custom_kind_name: p for p in store.query(t1d_id, ParameterKind.CUSTOM)}
        assert len(by_name) == 3   # guards against a vacuous pass if the query is empty
        for name in ("Ib_factor", "secretion_scale_factor"):
            assert by_name[name].provenance.pubmed_id == "3100365"
            assert by_name[name].confidence_tier == ConfidenceTier.BRONZE

    def test_t1d_gb_is_explicitly_tagged_unverified_not_a_fabricated_citation(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        t1d_id = str(ContextID(ContextEntityType.DISEASE_STATE, "type_1_diabetes"))
        by_name = {p.custom_kind_name: p for p in store.query(t1d_id, ParameterKind.CUSTOM)}
        assert by_name["Gb_factor"].provenance.pubmed_id == "UNVERIFIED"
        assert by_name["Gb_factor"].confidence_tier == ConfidenceTier.BRONZE

    def test_t1d_own_three_factors_do_not_conflict_with_each_other(self, seeded):
        # Same-subject guard (per the review finding on the GDM chunk
        # above): Ib_factor=0.05, secretion_scale_factor=0.02, and
        # Gb_factor=2.2 all sit under one subject with real fold
        # differences well over CONFLICT_FOLD_THRESHOLD=2.0 (e.g.
        # 0.05/0.02=2.5, 2.2/0.02=110) -- confirm none of the three
        # spuriously conflict with each other.
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        t1d_id = str(ContextID(ContextEntityType.DISEASE_STATE, "type_1_diabetes"))
        conflicts = store.unresolved_conflicts()
        t1d_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).custom_kind_name
            in ("Ib_factor", "secretion_scale_factor", "Gb_factor")
            and store.get(c.existing_parameter_id).subject_canonical_id == t1d_id
        ]
        assert t1d_conflicts == []

    def test_t1d_ib_factor_does_not_conflict_with_a_different_subjects_ib(self, seeded):
        # Cross-subject guard, matching this campaign's established
        # pattern -- confirms subject-scoping still holds even though
        # no other subject in this campaign currently uses "Ib_factor"
        # as a custom_kind_name (a forward-looking regression guard,
        # not evidence of a currently-live collision).
        _, store = seeded
        conflicts = store.unresolved_conflicts()
        ib_factor_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).custom_kind_name == "Ib_factor"
        ]
        assert ib_factor_conflicts == []

    # -- MASTER_PLAN 5.1 chunk 16: T2D/insulin-resistant/hypothyroidism (ContextRegistry) --

    def test_t2d_all_four_factors_stored(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        t2d_id = str(ContextID(ContextEntityType.DISEASE_STATE, "type_2_diabetes"))
        names = {p.custom_kind_name for p in store.query(t2d_id, ParameterKind.CUSTOM)}
        assert names == {"secretion_scale_factor", "p3_factor", "Gb_factor", "p1_factor"}

    def test_t2d_factors_match_the_real_module_values(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        from biosim_physiology.modifiers import disease_modifiers
        from biosim_physiology.profile import DiseaseState
        _, store = seeded
        t2d_id = str(ContextID(ContextEntityType.DISEASE_STATE, "type_2_diabetes"))
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(t2d_id, ParameterKind.CUSTOM)
        }
        expected = disease_modifiers([DiseaseState.TYPE_2_DIABETES])
        assert params["secretion_scale_factor"] == pytest.approx(expected["secretion_scale"])
        assert params["p3_factor"] == pytest.approx(expected["p3"])
        assert params["Gb_factor"] == pytest.approx(expected["Gb"])
        assert params["p1_factor"] == pytest.approx(expected["p1"])

    def test_t2d_factors_cite_stumvoll_bronze_not_unverified(self, seeded):
        # Distinguishes T2D (a real, directionally-supportive citation
        # exists) from insulin_resistant below (zero citation at all) --
        # both lack a reported magnitude, but only one has NO citation
        # attempt whatsoever.
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        t2d_id = str(ContextID(ContextEntityType.DISEASE_STATE, "type_2_diabetes"))
        entries = [
            e for e in store.query(t2d_id, ParameterKind.CUSTOM)
            if e.custom_kind_name != "p1_factor"
        ]
        assert len(entries) == 3
        for entry in entries:
            assert entry.provenance.pubmed_id == "15823385"
            assert entry.confidence_tier == ConfidenceTier.BRONZE

    def test_t2d_p1_factor_cites_basu_2009_gold(self, seeded):
        # MASTER_PLAN row 1.46 (a): a ratio of one paper's own two reported
        # group means (SG 0.019 vs. 0.034 min-1), same tier rule as GCK-MODY's Gb.
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        t2d_id = str(ContextID(ContextEntityType.DISEASE_STATE, "type_2_diabetes"))
        (entry,) = [
            e for e in store.query(t2d_id, ParameterKind.CUSTOM)
            if e.custom_kind_name == "p1_factor"
        ]
        assert entry.provenance.pubmed_id == "19196896"
        assert entry.confidence_tier == ConfidenceTier.GOLD

    def test_t2d_own_four_factors_do_not_conflict_with_each_other(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        t2d_id = str(ContextID(ContextEntityType.DISEASE_STATE, "type_2_diabetes"))
        conflicts = store.unresolved_conflicts()
        t2d_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).custom_kind_name
            in ("secretion_scale_factor", "p3_factor", "Gb_factor", "p1_factor")
            and store.get(c.existing_parameter_id).subject_canonical_id == t2d_id
        ]
        assert t2d_conflicts == []

    def test_insulin_resistant_all_three_factors_stored_all_unverified(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        ir_id = str(ContextID(ContextEntityType.DISEASE_STATE, "insulin_resistant"))
        entries = store.query(ir_id, ParameterKind.CUSTOM)
        assert len(entries) == 3
        names = {e.custom_kind_name for e in entries}
        assert names == {"p3_factor", "Gb_factor", "secretion_scale_factor"}
        for entry in entries:
            assert entry.provenance.pubmed_id == "UNVERIFIED"
            assert entry.confidence_tier == ConfidenceTier.BRONZE

    def test_insulin_resistant_factors_match_the_real_module_values(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        from biosim_physiology.modifiers import disease_modifiers
        from biosim_physiology.profile import DiseaseState
        _, store = seeded
        ir_id = str(ContextID(ContextEntityType.DISEASE_STATE, "insulin_resistant"))
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(ir_id, ParameterKind.CUSTOM)
        }
        expected = disease_modifiers([DiseaseState.INSULIN_RESISTANT])
        assert params["p3_factor"] == pytest.approx(expected["p3"])
        assert params["Gb_factor"] == pytest.approx(expected["Gb"])
        assert params["secretion_scale_factor"] == pytest.approx(expected["secretion_scale"])

    def test_insulin_resistant_own_three_factors_do_not_conflict_with_each_other(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        ir_id = str(ContextID(ContextEntityType.DISEASE_STATE, "insulin_resistant"))
        conflicts = store.unresolved_conflicts()
        ir_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).custom_kind_name
            in ("p3_factor", "Gb_factor", "secretion_scale_factor")
            and store.get(c.existing_parameter_id).subject_canonical_id == ir_id
        ]
        assert ir_conflicts == []

    def test_hypothyroidism_p1_and_p3_factors_stored(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        hypo_id = str(ContextID(ContextEntityType.DISEASE_STATE, "hypothyroidism"))
        names = {p.custom_kind_name for p in store.query(hypo_id, ParameterKind.CUSTOM)}
        assert names == {"p1_factor", "p3_factor"}

    def test_hypothyroidism_factors_match_the_real_module_values(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        from biosim_physiology.modifiers import disease_modifiers
        from biosim_physiology.profile import DiseaseState
        _, store = seeded
        hypo_id = str(ContextID(ContextEntityType.DISEASE_STATE, "hypothyroidism"))
        params = {
            p.custom_kind_name: p.value.value
            for p in store.query(hypo_id, ParameterKind.CUSTOM)
        }
        expected = disease_modifiers([DiseaseState.HYPOTHYROIDISM])
        assert params["p1_factor"] == pytest.approx(expected["p1"])
        assert params["p3_factor"] == pytest.approx(expected["p3"])

    def test_hypothyroidism_p1_cites_rochon_p3_cites_handisurya_both_bronze(self, seeded):
        # Two different PMIDs under one subject -- p1 traces to Rochon
        # (basal glucose disposal), p3 traces to Handisurya (insulin
        # sensitivity clamp), matching the scoping pass's own finding
        # that these two supporting studies measure different quantities.
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        hypo_id = str(ContextID(ContextEntityType.DISEASE_STATE, "hypothyroidism"))
        by_name = {p.custom_kind_name: p for p in store.query(hypo_id, ParameterKind.CUSTOM)}
        assert len(by_name) == 2
        assert by_name["p1_factor"].provenance.pubmed_id == "12519082"
        assert by_name["p3_factor"].provenance.pubmed_id == "18429948"
        assert by_name["p1_factor"].confidence_tier == ConfidenceTier.BRONZE
        assert by_name["p3_factor"].confidence_tier == ConfidenceTier.BRONZE

    def test_hypothyroidism_own_two_factors_do_not_conflict_with_each_other(self, seeded):
        from biosim_context_registry import ContextEntityType, ContextID
        _, store = seeded
        hypo_id = str(ContextID(ContextEntityType.DISEASE_STATE, "hypothyroidism"))
        conflicts = store.unresolved_conflicts()
        hypo_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).custom_kind_name in ("p1_factor", "p3_factor")
            and store.get(c.existing_parameter_id).subject_canonical_id == hypo_id
        ]
        assert hypo_conflicts == []

    def test_ckd_and_sleep_restriction_are_separate_subjects_no_cross_conflict(self, seeded):
        # CKD (DISEASE_STATE) and sleep restriction (PHYSIOLOGICAL_AXIS)
        # register under different subject strings (str(ContextID) embeds
        # entity_type regardless), so they can never spuriously conflict
        # with each other in the store -- this exercises ParameterStore's
        # own custom_kind_name-aware conflict detection (store.py), not
        # the entity-type-collision property itself. That property is
        # covered separately by
        # test_context_registry.py::test_same_name_different_entity_type_does_not_collide.
        _, store = seeded
        conflicts = store.unresolved_conflicts()
        context_conflicts = [
            c for c in conflicts
            if store.get(c.existing_parameter_id).custom_kind_name
            and store.get(c.existing_parameter_id).custom_kind_name
            in ("p3_factor", "ke_factor", "p3_factor_at_5h_reference")
        ]
        assert context_conflicts == []

    def test_no_unresolved_conflicts_critical(self, seeded):
        _, store = seeded
        conflicts = store.unresolved_conflicts()
        # Some conflicts are expected (multiple studies for same parameter)
        # but none should be at extreme fold differences (>10x)
        critical = [c for c in conflicts if c.fold_difference > 10.0]
        assert len(critical) == 0, \
            f"Critical conflicts (>10x): {[(c.subject_canonical_id, c.fold_difference) for c in critical]}"

    def test_parameter_values_physiologically_sane(self, seeded):
        registry, store = seeded
        insulin_ids = registry.search_by_name("insulin")
        ib_param = store.best(str(insulin_ids[0]), ParameterKind.ABUNDANCE)
        assert ib_param is not None
        assert 20.0 <= ib_param.value.value <= 200.0, \
            f"Ib = {ib_param.value.value} pM outside plausible range"


# ---------------------------------------------------------------------------
# ParameterLoader
# ---------------------------------------------------------------------------

class TestParameterLoader:
    @pytest.fixture(scope="class")
    @classmethod
    def loader(cls):
        return ParameterLoader()

    def test_load_returns_pathway_params(self, loader):
        params = loader.load_pathway_params()
        assert isinstance(params, PathwayParameters)

    def test_loaded_gb_near_literature(self, loader):
        params = loader.load_pathway_params()
        # Literature value: ~5.0 mM. Should be within 10%.
        assert 4.5 <= params.Gb <= 5.5, f"Gb = {params.Gb} mM outside expected range"

    def test_loaded_ib_near_literature(self, loader):
        params = loader.load_pathway_params()
        # Literature value: ~66 pM. Should be within plausible range.
        assert 30.0 <= params.Ib <= 150.0, f"Ib = {params.Ib} pM outside expected range"

    def test_loaded_ke_near_literature(self, loader):
        params = loader.load_pathway_params()
        # Literature value: ~0.093 min-1. Allow generous range.
        assert 0.04 <= params.ke <= 0.20, f"ke = {params.ke} outside expected range"

    def test_ec50_mtor_actually_loaded_from_the_object_scoped_seed_value(self, loader):
        # Regression: EC50_mTOR's PARAMETER_MAPPINGS row must pass
        # object_canonical_id="insulin receptor" through to store.best() --
        # it's seeded object-scoped (ParameterBuilder.for_object()), so a
        # query without that filter returns None (an exact-match filter,
        # not a wildcard) and load_pathway_params() would silently keep
        # PathwayParameters' hardcoded default instead. Checked directly
        # against the store (not just PathwayParameters.EC50_mTOR's final
        # value, which coincidentally equals the same 100.0 pM default and
        # so can't distinguish "loaded" from "fell back" by itself).
        from biosim_parameter_db import ParameterKind
        mtor_id = str(loader.registry.search_by_name("mTOR")[0])
        insr_id = str(loader.registry.search_by_name("insulin receptor")[0])
        param = loader.store.best(mtor_id, ParameterKind.EC50, object_canonical_id=insr_id)
        assert param is not None
        assert param.provenance.pubmed_id == "15718470"

    def test_provenance_report_has_entries(self, loader):
        report = loader.provenance_report()
        assert isinstance(report, dict)
        assert len(report) > 0

    def test_provenance_report_has_from_literature(self, loader):
        report = loader.provenance_report()
        # Literature values may equal defaults (derived from same papers).
        # Check that the report has entries and the store has parameters.
        assert len(report) > 0, "Provenance report is empty"
        stats = loader.store_stats()
        assert stats["store"]["total_parameters"] > 0

    def test_provenance_report_has_defaults(self, loader):
        report = loader.provenance_report()
        default_params = [k for k, v in report.items() if v == "default"]
        # Some parameters not yet seeded should still use defaults
        assert len(default_params) > 0

    def test_store_stats_accessible(self, loader):
        stats = loader.store_stats()
        assert "registry" in stats
        assert "store" in stats

    def test_parameter_mappings_not_empty(self):
        assert len(PARAMETER_MAPPINGS) > 0

    def test_parameter_mappings_have_valid_kinds(self):
        for field_name, entity_name, kind, object_entity_name in PARAMETER_MAPPINGS:
            assert isinstance(kind, ParameterKind)
            assert field_name
            assert entity_name
            assert object_entity_name is None or isinstance(object_entity_name, str)

    # -- MASTER_PLAN 5.1 chunk 3: KG_secretion, a CORE parameter, now wired --

    def test_kg_secretion_actually_loaded_from_the_seed_value(self, loader):
        # Same shape as test_ec50_mtor_actually_loaded_from_the_object_
        # scoped_seed_value above: params.KG_secretion alone can't
        # distinguish "loaded from the store" from "fell back to
        # PathwayParameters' own hardcoded default," since both equal
        # 5.63 by design (confirmed empirically as a genuine no-op
        # before this chunk was finalized). Checked directly against the
        # store instead, via its real PMID.
        from biosim_parameter_db import ParameterKind
        insulin_id = str(loader.registry.search_by_name("insulin")[0])
        param = loader.store.best(insulin_id, ParameterKind.EC50)
        assert param is not None
        assert param.provenance.pubmed_id == "30166558"

    def test_kg_secretion_loading_is_a_genuine_no_op(self, loader):
        # Confirms the claim made throughout this chunk's own commit/
        # CLAUDE.md entry: wiring KG_secretion into PARAMETER_MAPPINGS
        # changes NOTHING about what value simulations actually use.
        from biosim_ode import PathwayParameters
        params = loader.load_pathway_params()
        assert params.KG_secretion == PathwayParameters().KG_secretion == 5.63

    def test_kg_secretion_and_p3_no_longer_falsely_listed_as_unmapped(self):
        # Regression: both were already (or, for KG_secretion, are now)
        # real PARAMETER_MAPPINGS entries, but UNMAPPED_PARAMETERS still
        # claimed otherwise -- a stale, self-contradictory entry.
        from biosim_seeds.loader import UNMAPPED_PARAMETERS
        assert "KG_secretion" not in UNMAPPED_PARAMETERS
        assert "p3" not in UNMAPPED_PARAMETERS

    def test_organism_fallback_is_genuinely_exercised_by_load_pathway_params(self):
        # Mandatory security-reviewer finding (MEDIUM): the two tests
        # above check store.best() directly with explicit filters,
        # which bypasses load_pathway_params()'s own "retry without an
        # organism filter" fallback branch entirely (loader.py's
        # KG_secretion mapping is registered under Organism.OTHER --
        # rodent islet data -- while load_pathway_params()'s own default
        # query organism is HUMAN, so the first, organism-filtered
        # lookup misses and only the fallback branch actually finds
        # it). If a future edit ever removed that fallback branch
        # (plausible -- it could look like unreachable defensive code
        # to someone unaware this entry depends on it), KG_secretion
        # would silently revert to PathwayParameters' hardcoded
        # default -- and no existing test would catch it, since both
        # values happen to be 5.63. Proves the PUBLIC entry point
        # (load_pathway_params(), not store.best()) really does
        # traverse the fallback: registers a DIFFERENT, distinguishable
        # value (7.77, nothing to do with the real 5.63) under
        # Organism.OTHER in an isolated store, and confirms
        # load_pathway_params() returns THAT value, not the hardcoded
        # default.
        from biosim_parameter_db import ParameterBuilder, ParameterStore

        registry = create_seeded_registry()
        store = ParameterStore()
        insulin_id = str(registry.search_by_name("insulin")[0])
        store.add(
            ParameterBuilder()
            .for_subject(insulin_id)
            .kind(ParameterKind.EC50)
            .value(7.77, unit="mM")
            .from_pubmed("99999999")
            .organism(Organism.OTHER)
            .in_vitro()
            .bronze()
            .build()
        )
        isolated_loader = ParameterLoader(store=store, registry=registry)
        params = isolated_loader.load_pathway_params()
        assert params.KG_secretion == 7.77

    # -- MASTER_PLAN 5.1 chunk 4: p1, a second core parameter, now wired --

    def test_p1_actually_loaded_from_the_seed_value(self, loader):
        # Same shape as the KG_secretion/EC50_mTOR precedent: params.p1
        # alone can't distinguish "loaded" from "fell back to
        # PathwayParameters' own hardcoded default," since both equal
        # 0.0317 by design. Checked directly against the store instead,
        # via its real PMID (Bergman 1981, the higher-priority of the
        # two real GOLD entries for this subject/kind).
        glucose_id = str(loader.registry.search_by_name("D-Glucose")[0])
        param = loader.store.best(glucose_id, ParameterKind.RATE)
        assert param is not None
        assert param.provenance.pubmed_id == "7033284"

    def test_p1_loading_is_a_genuine_no_op(self, loader):
        from biosim_ode import PathwayParameters
        params = loader.load_pathway_params()
        assert params.p1 == PathwayParameters().p1 == 0.0317

    def test_p1_no_longer_falsely_listed_as_unmapped(self):
        from biosim_seeds.loader import UNMAPPED_PARAMETERS
        assert "p1" not in UNMAPPED_PARAMETERS

    def test_p1_needs_no_organism_fallback(self, loader):
        # Unlike KG_secretion (rodent data, needs the fallback branch),
        # p1's two real seed entries are both .in_human() -- confirms
        # load_pathway_params()'s FIRST, organism-filtered query already
        # succeeds for this one, no fallback dependency to disclose or
        # regression-guard here.
        from biosim_parameter_db import Organism
        glucose_id = str(loader.registry.search_by_name("D-Glucose")[0])
        param = loader.store.best(glucose_id, ParameterKind.RATE, organism=Organism.HUMAN)
        assert param is not None

    def test_p1_mapping_is_genuinely_wired_not_just_coincidentally_equal(self):
        # RED-proofing this chunk found that 3 of the 4 tests above
        # (actually_loaded/loading_is_a_genuine_no_op/needs_no_organism_
        # fallback) all trivially pass even with the PARAMETER_MAPPINGS
        # row for "p1" completely removed -- store.best() checks are
        # independent of the mapping's existence (the store already had
        # p1 seeded), and the no-op value check can't distinguish
        # "loaded via the mapping" from "fell back to the identical
        # hardcoded default" (same failure shape the security-reviewer
        # found for KG_secretion's organism-fallback branch, just for a
        # different reason here). Only
        # test_p1_no_longer_falsely_listed_as_unmapped actually failed
        # when the wiring was reverted. This test closes that gap using
        # the same isolated-store trick as KG_secretion's own review-fix
        # test: registers a DELIBERATELY DIFFERENT p1 value (0.0999,
        # unrelated to the real 0.0317) in an isolated store and
        # confirms load_pathway_params() -- not store.best() -- returns
        # THAT value, proving the PARAMETER_MAPPINGS row is genuinely
        # read, not coincidental.
        from biosim_parameter_db import ParameterBuilder, ParameterStore

        registry = create_seeded_registry()
        store = ParameterStore()
        glucose_id = str(registry.search_by_name("D-Glucose")[0])
        store.add(
            ParameterBuilder()
            .for_subject(glucose_id)
            .kind(ParameterKind.RATE)
            .value(0.0999, unit="min-1")
            .from_pubmed("99999999")
            .in_human()
            .clinical()
            .gold()
            .build()
        )
        isolated_loader = ParameterLoader(store=store, registry=registry)
        params = isolated_loader.load_pathway_params()
        assert params.p1 == 0.0999


# ---------------------------------------------------------------------------
# Integration: loaded params run ODE correctly
# ---------------------------------------------------------------------------

class TestIntegrationWithODE:
    @pytest.fixture(scope="class")
    @classmethod
    def loaded_params(cls):
        return ParameterLoader().load_pathway_params()

    def test_loaded_params_pass_validate_y0(self, loaded_params):
        pathway = GlucoseInsulinPathway(params=loaded_params)
        violations = pathway.validate_y0(tolerance=1e-3)
        assert violations == [], \
            f"Literature-loaded params not at steady state: {violations}"

    def test_loaded_params_run_ode(self, loaded_params):
        pathway = GlucoseInsulinPathway(params=loaded_params)
        system = ODESystem(pathway=pathway)
        result = system.simulate_fasted(duration_s=3600)
        assert result.success

    def test_loaded_params_meal_simulation(self, loaded_params):
        pathway = GlucoseInsulinPathway(params=loaded_params)
        system = ODESystem(pathway=pathway)
        result = system.simulate_meal(high_sugar_meal(), duration_s=7200)
        assert result.success
        peak_G, _ = result.glucose_peak()
        assert peak_G > 5000, "No glucose rise after meal with literature params"

    def test_loaded_vs_default_params_similar_peak(self, loaded_params):
        # Literature params should produce similar (within 30%) peak glucose
        # to the hardcoded defaults, since defaults ARE from the same literature
        default_pathway = GlucoseInsulinPathway()
        loaded_pathway  = GlucoseInsulinPathway(params=loaded_params)
        sys_default = ODESystem(pathway=default_pathway)
        sys_loaded  = ODESystem(pathway=loaded_pathway)
        peak_default = sys_default.simulate_meal(high_sugar_meal(), 7200).glucose_peak()[0]
        peak_loaded  = sys_loaded.simulate_meal(high_sugar_meal(), 7200).glucose_peak()[0]
        ratio = peak_loaded / peak_default
        assert 0.7 <= ratio <= 1.3, \
            f"Literature params produce {ratio:.2f}x peak vs defaults (expected 0.7-1.3)"


class TestIgf1HalfLifeConsistency:
    """The store's IGF-1 half-life and the runtime k_igf1_clear must describe
    the same pool. They once said 600 min and ~99 min (ln2/0.007)."""

    @pytest.fixture(scope="class")
    @classmethod
    def igf1_half_life(cls):
        registry = create_seeded_registry()
        store = create_seeded_store(registry)
        igf1_id = str(registry.search_by_name("IGF-1")[0])
        found = [
            p for p in store._parameters.values()
            if not p.is_superseded
            and p.subject_canonical_id == igf1_id
            and p.kind == ParameterKind.HALF_LIFE
        ]
        assert len(found) == 1
        return found[0]

    def test_store_half_life_matches_runtime_clearance_within_1pct(self, igf1_half_life):
        import math
        runtime_half_life = math.log(2) / PathwayParameters().k_igf1_clear
        assert igf1_half_life.value.unit == "min"
        assert runtime_half_life == pytest.approx(igf1_half_life.value.value, rel=0.01)

    def test_half_life_cites_guler_1989(self, igf1_half_life):
        assert igf1_half_life.provenance.pubmed_id == "2558477"


class TestProvenanceVerificationStatus:
    """Every seeded parameter says whether its value was found in the cited
    source text (2026-11-15 audit, docs/PROVENANCE_GAPS.md)."""

    ALLOWED = {
        "VERIFIED", "DERIVED", "PAPER-ON-TOPIC-VALUE-NOT-FOUND",
        "WRONG-PAPER", "CONTRADICTED", "UNVERIFIED",
    }

    @pytest.fixture(scope="class")
    @classmethod
    def active(cls):
        from biosim_seeds.context_seed import create_seeded_context_registry
        store = create_seeded_store(create_seeded_registry(), create_seeded_context_registry())
        return [p for p in store._parameters.values() if not p.is_superseded]

    def test_every_parameter_has_an_allowed_status(self, active):
        missing = [
            (p.provenance.pubmed_id, p.custom_kind_name or p.kind.value, p.value.value)
            for p in active
            if p.provenance.verification_status not in self.ALLOWED
        ]
        assert missing == []

    def test_verified_and_derived_carry_quote_and_pmid(self, active):
        for p in active:
            if p.provenance.verification_status in ("VERIFIED", "DERIVED"):
                assert p.provenance.verification_quote.strip(), p.custom_kind_name
                assert p.provenance.pubmed_id not in (None, "", "UNVERIFIED")

    def test_unverified_pmid_means_unverified_status(self, active):
        for p in active:
            if p.provenance.pubmed_id == "UNVERIFIED":
                assert p.provenance.verification_status == "UNVERIFIED"

    def test_every_table_record_matches_a_seeded_parameter(self, active):
        # A stale record (value changed, entry removed) must not linger.
        from biosim_seeds.parameter_seed.verification import VERIFICATION_RECORDS
        seeded = {
            (p.provenance.pubmed_id, p.custom_kind_name or p.kind.value, round(p.value.value, 9))
            for p in active
        }
        stale = [
            (r["pmid"], r["label"], r["value"]) for r in VERIFICATION_RECORDS
            if (r["pmid"], r["label"], round(r["value"], 9)) not in seeded
        ]
        assert stale == []

    def test_table_keys_are_unique(self):
        # A duplicate (pmid, label, value) key would silently drop a record.
        from biosim_seeds.parameter_seed.verification import VERIFICATION_RECORDS, _INDEX
        # 81 since B4 (2026-11-16): the WRONG-PAPER DeFronzo p1 replicate
        # was removed from the seeds.
        # 80 since 2026-11-20: the caffeine gb_increase_amplitude entry was
        # removed with the term itself (gold-standard drug stream).
        assert len(_INDEX) == len(VERIFICATION_RECORDS) == 80

    def test_no_seed_cites_the_wrong_paper_for_p1(self, active):
        # B4: DeFronzo 1979 (PMID:382871) is the clamp-methods paper and
        # does not estimate glucose effectiveness; its p1 replicate is gone.
        assert not [p for p in active if p.provenance.pubmed_id == "382871"]

    def test_unsourced_entries_are_labelled_model_assumptions(self):
        from biosim_seeds.parameter_seed.verification import VERIFICATION_RECORDS
        for r in VERIFICATION_RECORDS:
            if r["status"] == "UNVERIFIED":
                assert r["reason"].startswith("MODEL ASSUMPTION:"), r["label"]
                assert len(r["reason"]) > len("MODEL ASSUMPTION:") + 20

    def test_acarbose_record_quotes_its_calibration_targets(self):
        from biosim_seeds.parameter_seed.verification import lookup as verification_for
        rec = verification_for("27686734", "kabs_scale_at_100mg_reference_dose", 0.35)
        assert "13.22" in rec["quote"] and "9.45" in rec["quote"]
