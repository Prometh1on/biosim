# SAVE TO: packages\data-seeds\tests\test_session.py
"""
Tests for PlatformSession and extended variant library.
"""
import pytest
from biosim_seeds import PlatformSession
from biosim_genotype.effect_db import extended_t2d_variants, VariantEffectDB


# ---------------------------------------------------------------------------
# PlatformSession
# ---------------------------------------------------------------------------

class TestPlatformSession:
    @pytest.fixture(scope="class")
    def session(self):
        return PlatformSession()

    def test_session_creates(self, session):
        assert session is not None

    def test_base_params_loaded(self, session):
        from biosim_ode import PathwayParameters
        assert isinstance(session.base_params, PathwayParameters)

    def test_store_stats(self, session):
        stats = session.store_stats()
        assert "registry" in stats
        assert "store" in stats
        assert stats["store"]["total_parameters"] > 0

    def test_simulate_fasted(self, session):
        output = session.simulate("fasted")
        assert output.success

    def test_simulate_high_sugar_meal(self, session):
        output = session.simulate("high_sugar_meal")
        assert output.success
        assert output.expert_metrics.glucose_peak_mM > 5.0

    def test_simulate_with_scenario_name(self, session):
        output = session.simulate("light_snack")
        assert output.success

    def test_simulate_actually_uses_base_params(self):
        # Regression: simulate() used to build a pathway from
        # self._base_params and then discard it, running every simulation
        # through a SimulationRunner() that hardcoded its own fresh
        # PathwayParameters() internally -- self._base_params never
        # actually reached the simulated output (only export_sbml() used
        # it). Invisible via the public use_literature_params flag alone
        # today, since the seeded literature values happen to equal the
        # class defaults exactly -- so this reaches into the "private"
        # _base_params attribute directly to prove the wiring, the same
        # thing a genuinely different literature value would exercise.
        from biosim_ode import PathwayParameters
        session = PlatformSession()
        session._base_params = PathwayParameters(Gb=6.5)
        output = session.simulate("fasted")
        assert output.expert_metrics.glucose_baseline_uM == pytest.approx(6500.0)

    def test_compare_scenarios(self, session):
        comparison = session.compare("high_sugar_meal", "light_snack")
        assert comparison.peak_glucose_diff_uM > 0

    def test_report_generation(self, session):
        output = session.simulate("standard_meal")
        report = session.report(output)
        assert report is not None
        assert len(report.glucose_section) > 20

    def test_assess_risk(self, session):
        output = session.simulate("high_sugar_meal")
        risk = session.assess_risk(output)
        assert risk is not None
        assert risk.overall is not None

    def test_export_sbml(self, session):
        sbml = session.export_sbml()
        assert isinstance(sbml, str)
        assert "sbml" in sbml.lower()

    def test_export_csv(self, session):
        output = session.simulate("standard_meal")
        csv_str = session.export_csv(output)
        assert "time_s" in csv_str

    def test_validate(self, session):
        from biosim_validation import ValidationVerdict
        report = session.validate()
        # Known FAIL since 2026-11-16 (B7 + B4): only the circular peak-glucose
        # point fails; anything else failing is a new regression.
        assert [r.name for r in report.results if r.benchmark.is_scored and not r.passed] == ["peak_glucose"], report.as_text()

    def test_bundle(self, session):
        output = session.simulate("fasted")
        bundle = session.bundle(output)
        assert bundle is not None
        assert bundle.platform_version is not None

    def test_audit_privacy(self, session):
        output = session.simulate("high_sugar_meal")
        result = session.audit_privacy(output)
        assert result.passed, f"Privacy issues: {result.issues}"

    def test_provenance(self, session):
        report = session.provenance()
        assert isinstance(report, dict)
        assert len(report) > 0

    def test_cache_stats(self, session):
        session.simulate("fasted")
        session.simulate("fasted")  # should hit cache
        stats = session.cache_stats()
        assert "hit_rate" in stats

    def test_simulate_actually_hits_the_cache(self):
        # Regression: this test's sibling above ("should hit cache") never
        # actually asserted the hit count, so it silently passed both
        # before and after simulate() was wired to bypass the cache
        # entirely (both branches used to construct a fresh
        # SimulationRunner() directly, never touching self._cached_runner
        # -- confirmed empirically: repeated identical simulate() calls
        # left hits=0/misses=0/size=0). A fresh session (not the shared
        # class-scoped fixture, to avoid interference from other tests'
        # simulate() calls) proves the actual count now.
        session = PlatformSession()
        session.simulate("fasted")
        session.simulate("fasted")
        stats = session.cache_stats()
        assert stats["hits"] == 1
        assert stats["misses"] == 1
        assert stats["size"] == 1

    def test_simulate_with_profile(self, session):
        from biosim_physiology import PhysiologicalProfile, DiseaseState
        profile = PhysiologicalProfile(
            age_years=60,
            disease_states=[DiseaseState.TYPE_2_DIABETES],
        )
        output = session.simulate("high_sugar_meal", profile=profile)
        assert output.success

    def test_simulate_with_variants(self, session):
        from biosim_genotype import Genotype, Variant
        variants = [Variant("rs7903146", Genotype.HETEROZYGOUS)]
        output = session.simulate("high_sugar_meal", variants=variants)
        assert output.success
        assert output.genetic_personalization is not None

    def test_t2d_profile_higher_glucose(self, session):
        from biosim_physiology import PhysiologicalProfile, DiseaseState
        healthy = session.simulate("high_sugar_meal")
        t2d = session.simulate(
            "high_sugar_meal",
            profile=PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES]),
        )
        assert t2d.expert_metrics.glucose_peak_uM > healthy.expert_metrics.glucose_peak_uM


# ---------------------------------------------------------------------------
# Curation-from-sensitivity wiring (MASTER_PLAN 5.2b)
# ---------------------------------------------------------------------------

class TestPrioritizeCurationFromSensitivity:
    """
    PlatformSession.prioritize_curation_from_sensitivity() -- the first
    real caller of biosim_uncertainty.curation_bridge.
    prioritise_from_sensitivity(), which existed and was tested but was
    never invoked from any real entry point (MASTER_PLAN 5.2's own
    scoping finding).

    Own fixture (function-scoped, not the class-scoped `session` used
    above) -- each test needs a session with an empty curation queue to
    start from, not one another test in this file already enqueued into.
    """

    @pytest.fixture
    def fresh_session(self):
        return PlatformSession()

    def test_runs_and_enqueues_real_items(self, fresh_session):
        from biosim_uncertainty import BridgeResult

        result = fresh_session.prioritize_curation_from_sensitivity()

        assert isinstance(result, BridgeResult)
        assert result.n_items_enqueued > 0
        assert len(fresh_session.curation_queue.pending()) == result.n_items_enqueued

    def test_curation_queue_property_reflects_the_same_queue(self, fresh_session):
        fresh_session.prioritize_curation_from_sensitivity()
        # The method must enqueue into THIS session's own queue -- not a
        # fresh, disconnected CurationQueue instance the caller can never
        # reach.
        assert len(fresh_session.curation_queue.pending()) > 0

    def test_custom_parameter_names_restricts_analysis(self, fresh_session):
        result = fresh_session.prioritize_curation_from_sensitivity(
            parameter_names=["Gb"],
        )
        assert set(result.parameters_upgraded) <= {"Gb"}

    def test_high_sensitivity_threshold_forwarded(self, fresh_session):
        # An impossibly high threshold (>1.0, no real sensitivity index
        # can exceed it) must suppress every enqueue -- confirms
        # **prioritise_kwargs genuinely reaches prioritise_from_sensitivity(),
        # not silently dropped.
        result = fresh_session.prioritize_curation_from_sensitivity(
            high_sensitivity_threshold=999.0,
        )
        assert result.n_items_enqueued == 0

    def test_repeated_calls_do_not_duplicate_queue_entries(self, fresh_session):
        first = fresh_session.prioritize_curation_from_sensitivity()
        n_after_first = len(fresh_session.curation_queue.pending())
        fresh_session.prioritize_curation_from_sensitivity()
        n_after_second = len(fresh_session.curation_queue.pending())

        assert first.n_items_enqueued > 0
        assert n_after_second == n_after_first

    def test_export_curation_queue_csv_reflects_enqueued_items(self, fresh_session):
        # MASTER_PLAN 5.2d -- the minimal consumer completing the
        # producer -> queue -> consumer loop 5.2's own scoping found
        # missing. Real end-to-end: real sensitivity-driven items,
        # exported through the real CSVExporter, parsed back.
        import csv
        import io

        result = fresh_session.prioritize_curation_from_sensitivity()
        csv_str = fresh_session.export_curation_queue_csv()

        rows = list(csv.DictReader(io.StringIO(csv_str)))
        assert len(rows) == result.n_items_enqueued
        assert "item_id" in rows[0]
        assert "priority" in rows[0]

    def test_export_curation_queue_csv_empty_queue_does_not_crash(self, fresh_session):
        # No prioritize_curation_from_sensitivity() call -- queue starts empty.
        csv_str = fresh_session.export_curation_queue_csv()
        assert isinstance(csv_str, str)
        assert "item_id" in csv_str  # header only


# ---------------------------------------------------------------------------
# Extended T2D variant library
# ---------------------------------------------------------------------------

class TestExtendedVariants:
    @pytest.fixture(scope="class")
    def db(self):
        return extended_t2d_variants()

    def test_extended_db_has_more_than_default(self, db):
        default = VariantEffectDB.default()
        assert db.count() > default.count()

    def test_kcnq1_present(self, db):
        assert db.has("rs2237892")

    def test_cdkn2a_present(self, db):
        assert db.has("rs10811661")

    def test_cdkal1_present(self, db):
        assert db.has("rs7756992")

    def test_igf2bp2_present(self, db):
        assert db.has("rs4402960")

    def test_adcy5_present(self, db):
        assert db.has("rs11708067")

    def test_hhex_present(self, db):
        assert db.has("rs1111875")

    def test_wfs1_protective(self, db):
        records = db.get("rs10010131")
        assert any(r.effect_per_allele > 0 for r in records), \
            "WFS1 A allele should be protective (positive effect)"

    def test_igf2bp2_has_igf1_effect(self, db):
        records = db.get("rs4402960")
        fields = [r.parameter_field for r in records]
        assert "IGF1b" in fields, "IGF2BP2 should have an IGF-1 effect"

    def test_kcnq1_is_eas_population(self, db):
        from biosim_genotype.effect_db import Population
        records = db.get("rs2237892")
        assert any(r.population == Population.EAS for r in records)

    def test_all_variants_have_pubmed_id(self, db):
        for rsid in db.all_rsids():
            for r in db.get(rsid):
                assert r.pubmed_id, f"{rsid} missing PubMed ID"

    def test_extended_db_with_personalizer(self):
        from biosim_genotype import GenotypePersonalizer, Genotype, Variant
        db = extended_t2d_variants()
        personalizer = GenotypePersonalizer(effect_db=db)
        variants = [
            Variant("rs7903146",  Genotype.HETEROZYGOUS),  # TCF7L2
            Variant("rs10811661", Genotype.HETEROZYGOUS),  # CDKN2A
            Variant("rs7756992",  Genotype.HETEROZYGOUS),  # CDKAL1
        ]
        result = personalizer.personalize(variants)
        # Multiple beta-cell variants should compound the secretion reduction
        assert result.personalized_params.secretion_scale < 0.80, \
            "Three secretion-reducing variants should reduce secretion below 80%"

    def test_wfs1_protective_increases_secretion(self):
        from biosim_genotype import GenotypePersonalizer, Genotype, Variant
        db = extended_t2d_variants()
        personalizer = GenotypePersonalizer(effect_db=db)
        variants = [Variant("rs10010131", Genotype.HOMOZYGOUS_ALT)]  # WFS1 protective
        result = personalizer.personalize(variants)
        assert result.personalized_params.secretion_scale > 1.0, \
            "WFS1 protective variant should increase secretion"
