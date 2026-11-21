# SAVE TO: packages\infrastructure\tests\test_infrastructure.py
"""
Test suite for biosim_infra.

Coverage:
    - SimulationCache: get/put, LRU eviction, hit rate, thread safety
    - CachedRunner: second call hits cache, stats
    - PrivacyAuditor: no rsIDs in outputs, genetic summary clean
    - ReproducibilityBundle: to_json/from_json round-trip, citation
    - ModelVersion: citation format
    - create_bundle: produces valid bundle from output
"""

import json
import threading
import pytest

from biosim_infra import (
    CachedRunner, CURRENT_VERSION, MODEL_VERSION, ModelVersion,
    PLATFORM_VERSION, PrivacyAuditResult, PrivacyAuditor,
    ReproducibilityBundle, SimulationCache, create_bundle,
)
from biosim_api import (
    SimulationRunner, high_sugar_meal_scenario,
    fasted_scenario, light_snack_scenario, MealEvent, ScenarioDefinition,
)


# ---------------------------------------------------------------------------
# SimulationCache
# ---------------------------------------------------------------------------

class TestSimulationCache:
    def test_miss_on_empty_cache(self):
        cache = SimulationCache()
        assert cache.get(high_sugar_meal_scenario()) is None

    def test_hit_after_put(self):
        cache = SimulationCache()
        scenario = high_sugar_meal_scenario()
        output = SimulationRunner().run(scenario)
        cache.put(scenario, output)
        result = cache.get(scenario)
        assert result is not None

    def test_same_scenario_same_key(self):
        cache = SimulationCache()
        s1 = high_sugar_meal_scenario()
        s2 = high_sugar_meal_scenario()
        output = SimulationRunner().run(s1)
        cache.put(s1, output)
        assert cache.get(s2) is not None

    def test_different_scenario_different_key(self):
        cache = SimulationCache()
        s_large = high_sugar_meal_scenario()
        s_small = light_snack_scenario()
        output = SimulationRunner().run(s_large)
        cache.put(s_large, output)
        assert cache.get(s_small) is None

    def test_different_drug_timing_different_key(self):
        # Regression: to_dict() omitted DrugEvent.hours_before_simulation,
        # so a caffeine-now vs. caffeine-3h-ago scenario (otherwise
        # identical) hashed to the same cache key and would have silently
        # served one request's cached result for the other's request.
        from biosim_api import DrugEvent
        cache = SimulationCache()
        s_now = ScenarioDefinition(
            name="caffeine",
            drug_events=[DrugEvent(time_s=0.0, drug="caffeine", dose_mg=100.0,
                                    hours_before_simulation=0.0)],
        )
        s_earlier = ScenarioDefinition(
            name="caffeine",
            drug_events=[DrugEvent(time_s=0.0, drug="caffeine", dose_mg=100.0,
                                    hours_before_simulation=3.0)],
        )
        output = SimulationRunner().run(s_now)
        cache.put(s_now, output)
        assert cache.get(s_earlier) is None

    def test_different_glycaemic_index_different_key(self):
        # Regression: to_dict() omitted MealEvent.glycaemic_index/kabs.
        cache = SimulationCache()
        s_gi100 = ScenarioDefinition(
            name="meal", meal_events=[MealEvent(time_s=0.0, carbs_g=75.0, glycaemic_index=100.0)]
        )
        s_gi30 = ScenarioDefinition(
            name="meal", meal_events=[MealEvent(time_s=0.0, carbs_g=75.0, glycaemic_index=30.0)]
        )
        output = SimulationRunner().run(s_gi100)
        cache.put(s_gi100, output)
        assert cache.get(s_gi30) is None

    def test_different_complex_carb_fraction_different_key(self):
        # Regression guard for MASTER_PLAN 1.4b(iii): to_dict() must
        # include MealEvent.complex_carb_fraction, or two scenarios that
        # differ only in it (e.g. a free-sugar drink vs. a starchy meal,
        # same carbs_g) would silently collide in the cache.
        cache = SimulationCache()
        s_starch = ScenarioDefinition(
            name="meal", meal_events=[MealEvent(time_s=0.0, carbs_g=75.0, complex_carb_fraction=1.0)]
        )
        s_free_sugar = ScenarioDefinition(
            name="meal", meal_events=[MealEvent(time_s=0.0, carbs_g=75.0, complex_carb_fraction=0.0)]
        )
        output = SimulationRunner().run(s_starch)
        cache.put(s_starch, output)
        assert cache.get(s_free_sugar) is None

    def test_different_protein_g_different_key(self):
        # Regression guard for MASTER_PLAN 1.28 (dietary protein
        # co-ingestion): to_dict() must include MealEvent.protein_g, or
        # two scenarios differing only in it (same carbs_g, different
        # protein) would silently collide in the cache.
        cache = SimulationCache()
        s_no_protein = ScenarioDefinition(
            name="meal", meal_events=[MealEvent(time_s=0.0, carbs_g=50.0, protein_g=0.0)]
        )
        s_with_protein = ScenarioDefinition(
            name="meal", meal_events=[MealEvent(time_s=0.0, carbs_g=50.0, protein_g=50.0)]
        )
        output = SimulationRunner().run(s_no_protein)
        cache.put(s_no_protein, output)
        assert cache.get(s_with_protein) is None

    def test_different_personalization_different_key(self):
        # Regression: to_dict() didn't include `personalization` at all --
        # a healthy and a Type-2-diabetes-personalized scenario, otherwise
        # identical, hashed to the exact same cache key.
        from biosim_api import PersonalizationSpec
        from biosim_physiology import PhysiologicalProfile, DiseaseState
        cache = SimulationCache()
        s_healthy = ScenarioDefinition(
            name="meal", meal_events=[MealEvent(time_s=0.0, carbs_g=75.0)]
        )
        s_t2d = ScenarioDefinition(
            name="meal", meal_events=[MealEvent(time_s=0.0, carbs_g=75.0)],
            personalization=PersonalizationSpec(
                physiological_profile=PhysiologicalProfile(
                    disease_states=[DiseaseState.TYPE_2_DIABETES]
                )
            ),
        )
        output = SimulationRunner().run(s_healthy)
        cache.put(s_healthy, output)
        assert cache.get(s_t2d) is None

    def test_different_genetic_variant_different_key(self):
        # Regression: to_dict() didn't include genetic_variants at all.
        from biosim_api import PersonalizationSpec
        from biosim_genotype import Variant, Genotype
        cache = SimulationCache()
        s_no_variant = ScenarioDefinition(
            name="meal", meal_events=[MealEvent(time_s=0.0, carbs_g=75.0)]
        )
        s_variant = ScenarioDefinition(
            name="meal", meal_events=[MealEvent(time_s=0.0, carbs_g=75.0)],
            personalization=PersonalizationSpec(
                genetic_variants=[Variant(rsid="rs7903146", genotype=Genotype.HOMOZYGOUS_ALT,
                                           observed_alleles="TT")]
            ),
        )
        output = SimulationRunner().run(s_no_variant)
        cache.put(s_no_variant, output)
        assert cache.get(s_variant) is None

    def test_different_observed_alleles_different_key(self):
        # Regression (found by code review, not caught on the first
        # attempt): the fingerprint originally hashed only
        # rsid+genotype -- a protective-homozygote call ("CC") and a
        # risk-homozygote call ("TT") for the same SNP share the same
        # rsid and the same Genotype.HOMOZYGOUS_ALT, but
        # dosage_for_allele() (which EffectComposer.compose() actually
        # uses) gives opposite results for them (0 vs. 2) -- the same
        # "23andMe homozygous call" bug class already fixed elsewhere in
        # this codebase, see KNOWN BUGS FIXED in CLAUDE.md.
        from biosim_api import PersonalizationSpec
        from biosim_genotype import Variant, Genotype
        cache = SimulationCache()
        s_cc = ScenarioDefinition(
            name="meal", meal_events=[MealEvent(time_s=0.0, carbs_g=75.0)],
            personalization=PersonalizationSpec(
                genetic_variants=[Variant(rsid="rs7903146", genotype=Genotype.HOMOZYGOUS_ALT,
                                           observed_alleles="CC")]
            ),
        )
        s_tt = ScenarioDefinition(
            name="meal", meal_events=[MealEvent(time_s=0.0, carbs_g=75.0)],
            personalization=PersonalizationSpec(
                genetic_variants=[Variant(rsid="rs7903146", genotype=Genotype.HOMOZYGOUS_ALT,
                                           observed_alleles="TT")]
            ),
        )
        output = SimulationRunner().run(s_cc)
        cache.put(s_cc, output)
        assert cache.get(s_tt) is None

    def test_to_dict_never_includes_raw_rsid(self):
        # Privacy: to_dict() is documented as usable for archiving/API
        # transport, and genetic_variants can carry raw rsIDs (rule 7,
        # CLAUDE.md). The cache-key fix must distinguish variant sets
        # without leaking the identifying rsid into that dict.
        from biosim_api import PersonalizationSpec
        from biosim_genotype import Variant, Genotype
        s = ScenarioDefinition(
            name="genetic", meal_events=[MealEvent(time_s=0.0, carbs_g=75.0)],
            personalization=PersonalizationSpec(
                genetic_variants=[Variant(rsid="rs7903146", genotype=Genotype.HETEROZYGOUS)]
            ),
        )
        assert "rs7903146" not in json.dumps(s.to_dict())

    def test_lru_eviction(self):
        cache = SimulationCache(max_size=2)
        runner = SimulationRunner()
        s1 = fasted_scenario(3600)
        s2 = high_sugar_meal_scenario()
        s3 = light_snack_scenario()
        cache.put(s1, runner.run(s1))
        cache.put(s2, runner.run(s2))
        cache.put(s3, runner.run(s3))  # evicts s1
        assert cache.size == 2
        assert cache.get(s1) is None   # evicted
        assert cache.get(s3) is not None

    def test_hit_rate_zero_on_empty(self):
        cache = SimulationCache()
        cache.get(high_sugar_meal_scenario())
        assert cache.hit_rate == pytest.approx(0.0)

    def test_hit_rate_increases_on_hit(self):
        cache = SimulationCache()
        s = high_sugar_meal_scenario()
        output = SimulationRunner().run(s)
        cache.get(s)  # miss - not in cache yet
        cache.put(s, output)
        cache.get(s)  # hit - now in cache
        assert cache.hit_rate == pytest.approx(0.5)  # 1 miss + 1 hit

    def test_clear_resets_cache(self):
        cache = SimulationCache()
        s = high_sugar_meal_scenario()
        cache.put(s, SimulationRunner().run(s))
        cache.clear()
        assert cache.size == 0
        assert cache.get(s) is None

    def test_invalidate(self):
        cache = SimulationCache()
        s = high_sugar_meal_scenario()
        cache.put(s, SimulationRunner().run(s))
        removed = cache.invalidate(s)
        assert removed
        assert cache.get(s) is None

    def test_invalidate_unknown_returns_false(self):
        cache = SimulationCache()
        assert not cache.invalidate(high_sugar_meal_scenario())

    def test_stats_structure(self):
        cache = SimulationCache()
        stats = cache.stats()
        assert "size" in stats
        assert "hits" in stats
        assert "misses" in stats
        assert "hit_rate" in stats

    def test_max_size_one_raises(self):
        cache = SimulationCache(max_size=1)
        assert cache._max_size == 1

    def test_zero_max_size_raises(self):
        with pytest.raises(ValueError):
            SimulationCache(max_size=0)

    def test_thread_safe_concurrent_puts(self):
        cache = SimulationCache(max_size=100)
        errors = []
        scenarios = [
            ScenarioDefinition(name=f"s{i}", meal_events=[MealEvent(0, float(i*5+5))])
            for i in range(20)
        ]
        runner = SimulationRunner()
        outputs = [runner.run(s) for s in scenarios]

        def put_batch(start):
            try:
                for i in range(start, start + 5):
                    cache.put(scenarios[i], outputs[i])
            except Exception as e:
                errors.append(e)

        threads = [threading.Thread(target=put_batch, args=(i,)) for i in range(0, 20, 5)]
        for t in threads: t.start()
        for t in threads: t.join()
        assert errors == []


class TestCachedRunner:
    def test_second_run_hits_cache(self):
        runner = CachedRunner()
        scenario = fasted_scenario(3600)
        runner.run(scenario)  # miss
        runner.run(scenario)  # hit
        assert runner.cache.hit_rate > 0

    def test_cache_stats_accessible(self):
        runner = CachedRunner()
        stats = runner.cache.stats()
        assert "hit_rate" in stats

    def test_base_params_override_changes_output(self):
        # Same fix as SimulationRunner's own base_params (see runner.py
        # and CLAUDE.md) -- CachedRunner needed the identical parameter
        # threaded through before PlatformSession.simulate() could safely
        # route through its own cache without losing self._base_params.
        from biosim_ode import PathwayParameters
        default_out = CachedRunner().run(fasted_scenario(3600))
        custom_out = CachedRunner(
            base_params=PathwayParameters(Gb=6.5)
        ).run(fasted_scenario(3600))
        assert custom_out.expert_metrics.glucose_baseline_uM == pytest.approx(6500.0)
        assert (
            default_out.expert_metrics.glucose_baseline_uM
            != custom_out.expert_metrics.glucose_baseline_uM
        )


# ---------------------------------------------------------------------------
# PrivacyAuditor
# ---------------------------------------------------------------------------

class TestPrivacyAuditor:
    @pytest.fixture(scope="class")
    def clean_output(self):
        return SimulationRunner().run(high_sugar_meal_scenario())

    @pytest.fixture(scope="class")
    def genetic_output(self):
        from biosim_genotype import Genotype, Variant
        from biosim_api import PersonalizationSpec
        variants = [Variant("rs7903146", Genotype.HETEROZYGOUS)]
        scenario = high_sugar_meal_scenario(
            personalization=PersonalizationSpec(genetic_variants=variants)
        )
        return SimulationRunner().run(scenario)

    def test_clean_output_passes(self, clean_output):
        result = PrivacyAuditor().audit(clean_output)
        assert result.passed, f"Privacy issues: {result.issues}"

    def test_genetic_output_passes(self, genetic_output):
        result = PrivacyAuditor().audit(genetic_output)
        assert result.passed, f"Privacy issues: {result.issues}"

    def test_no_rsid_in_genetic_summary(self, genetic_output):
        summary_str = json.dumps(genetic_output.summary())
        import re
        assert not re.search(r'\brs\d{4,}\b', summary_str), \
            "rsID found in simulation summary"

    def test_no_raw_rsid_in_contributing_variants(self, genetic_output):
        # Regression: ComposedEffect.contributing_variants used to store
        # raw (rsid, dosage) tuples directly on PersonalizationResult
        # .applied_effects -- confirmed a real rsid ("rs7903146") was
        # present there even though PrivacyAuditor.audit() reported PASS
        # (its check was a vacuous hasattr(gp, 'variants') no-op --
        # that attribute never existed). contributing_variants must now
        # hold only privacy-safe fingerprints.
        gp = genetic_output.genetic_personalization
        assert gp is not None
        for effect in gp.applied_effects:
            for fingerprint, _dosage in effect.contributing_variants:
                assert not fingerprint.startswith("rs"), (
                    f"raw-looking rsid {fingerprint!r} found in "
                    f"contributing_variants -- must be a fingerprint"
                )

    def test_audit_catches_a_real_raw_rsid_if_reintroduced(self, genetic_output):
        # Regression: the OLD check (hasattr(gp, 'variants')) could never
        # fire regardless of content -- prove the NEW check actually
        # detects a raw rsid, by simulating exactly the regression this
        # fix closes (a future contributing_variants entry reverting to
        # a raw rsid).
        import dataclasses
        gp = genetic_output.genetic_personalization
        assert gp.applied_effects, "test setup needs at least one applied effect"
        regressed_effect = dataclasses.replace(
            gp.applied_effects[0],
            contributing_variants=[("rs7903146", 1)],
        )
        regressed_gp = dataclasses.replace(
            gp, applied_effects=[regressed_effect] + list(gp.applied_effects[1:]),
        )
        result = PrivacyAuditor().audit_personalization_result(regressed_gp)
        assert not result.passed
        assert any("rsID" in issue for issue in result.issues)

    def test_audit_catches_a_raw_rsid_in_population_warnings_if_reintroduced(self, genetic_output):
        # Regression: a mandatory security-review of the fix above found
        # a sibling leak this pass initially missed -- composer.py's
        # population-mismatch warning used to embed the raw rsid and
        # gene name directly, flowing straight into
        # SimulationOutput.warnings (REST-serializable) with neither the
        # old nor the newly-added Check 3 catching it. Fixed at the
        # source (composer.py no longer embeds rsid/gene); this proves
        # the new defense-in-depth backstop actually detects it if ever
        # reintroduced.
        import dataclasses
        gp = genetic_output.genetic_personalization
        assert gp is not None
        regressed = dataclasses.replace(
            gp, population_warnings=["rs7903146 (TCF7L2): effect may not generalise."],
        )
        result = PrivacyAuditor().audit_personalization_result(regressed)
        assert not result.passed
        assert any("rsID" in issue for issue in result.issues)

    def test_audit_result_structure(self, clean_output):
        result = PrivacyAuditor().audit(clean_output)
        assert hasattr(result, 'passed')
        assert hasattr(result, 'issues')
        assert hasattr(result, 'checked_at')

    def test_audit_result_summary(self, clean_output):
        result = PrivacyAuditor().audit(clean_output)
        s = result.summary()
        assert "passed" in s
        assert "n_issues" in s

    def test_contains_rsid_detection(self):
        assert PrivacyAuditor._contains_rsid("has rs7903146 in text")
        assert not PrivacyAuditor._contains_rsid("no variant here")
        assert not PrivacyAuditor._contains_rsid("rs is a short word")  # too short


# ---------------------------------------------------------------------------
# ReproducibilityBundle
# ---------------------------------------------------------------------------

class TestReproducibilityBundle:
    @pytest.fixture(scope="class")
    def bundle(self):
        output = SimulationRunner().run(high_sugar_meal_scenario())
        return create_bundle(output, random_seed=42)

    def test_create_bundle(self, bundle):
        assert isinstance(bundle, ReproducibilityBundle)

    def test_platform_version_present(self, bundle):
        assert bundle.platform_version == PLATFORM_VERSION

    def test_model_version_present(self, bundle):
        assert bundle.model_version == MODEL_VERSION

    def test_random_seed_preserved(self, bundle):
        assert bundle.random_seed == 42

    def test_python_version_present(self, bundle):
        assert bundle.python_version is not None
        assert "." in bundle.python_version

    def test_numpy_version_present(self, bundle):
        assert bundle.numpy_version is not None

    def test_scipy_version_present(self, bundle):
        assert bundle.scipy_version is not None

    def test_to_json_produces_valid_json(self, bundle):
        json_str = bundle.to_json()
        parsed = json.loads(json_str)
        assert "platform_version" in parsed
        assert "random_seed" in parsed

    def test_from_json_roundtrip(self, bundle):
        json_str = bundle.to_json()
        restored = ReproducibilityBundle.from_json(json_str)
        assert restored.platform_version == bundle.platform_version
        assert restored.random_seed == bundle.random_seed
        assert restored.model_version == bundle.model_version

    def test_citation_returns_string(self, bundle):
        cite = bundle.citation()
        assert isinstance(cite, str)
        assert len(cite) > 20
        assert PLATFORM_VERSION in cite

    def test_scenario_dict_actually_captures_meal_events(self, bundle):
        # Regression: scenario_dict used to be hand-built from just
        # scenario_name/scenario_description -- SimulationOutput never
        # carried the actual ScenarioDefinition, so meals/drugs/
        # personalization were silently absent from every bundle despite
        # the class's own "reproduce a simulation result exactly" claim.
        assert "meal_events" in bundle.scenario_dict
        assert len(bundle.scenario_dict["meal_events"]) == 1
        assert bundle.scenario_dict["meal_events"][0]["carbs_g"] == 75.0

    def test_scenario_dict_captures_drug_events_and_personalization(self):
        from biosim_api import (
            DrugEvent, MealEvent, PersonalizationSpec, ScenarioDefinition,
        )
        from biosim_physiology import DiseaseState, PhysiologicalProfile
        scenario = ScenarioDefinition(
            name="bundle_test",
            meal_events=[MealEvent(time_s=0.0, carbs_g=50.0)],
            drug_events=[DrugEvent(time_s=0.0, drug="metformin", dose_mg=1000.0)],
            personalization=PersonalizationSpec(
                physiological_profile=PhysiologicalProfile(
                    disease_states=[DiseaseState.TYPE_2_DIABETES]
                )
            ),
        )
        output = SimulationRunner().run(scenario)
        bundle = create_bundle(output, random_seed=1)
        assert bundle.scenario_dict["drug_events"][0]["drug"] == "metformin"
        assert bundle.scenario_dict["personalization"]["physiological_profile"][
            "disease_states"
        ] == ["type_2_diabetes"]

    def test_parameter_values_present(self, bundle):
        assert isinstance(bundle.parameter_values, dict)
        assert "p1" in bundle.parameter_values
        assert "ke" in bundle.parameter_values


# ---------------------------------------------------------------------------
# ModelVersion
# ---------------------------------------------------------------------------

class TestModelVersion:
    def test_current_version_exists(self):
        assert CURRENT_VERSION is not None

    def test_current_version_has_platform_version(self):
        assert CURRENT_VERSION.platform_version == PLATFORM_VERSION

    def test_cite_returns_string(self):
        cite = CURRENT_VERSION.cite()
        assert isinstance(cite, str)
        assert PLATFORM_VERSION in cite

    def test_custom_version(self):
        v = ModelVersion(
            platform_version="1.0.0",
            model_id="test_model",
            model_version="2.0.0",
            release_date="2025-01-01",
        )
        cite = v.cite()
        assert "1.0.0" in cite
        assert "test_model" in cite

    def test_model_version_matches_the_actual_ode_engine_not_a_separate_copy(self):
        # Regression: MODEL_VERSION used to be its own hardcoded literal
        # here, an independent copy of ODESystem.MODEL_VERSION rather
        # than a reference to it -- the two could silently drift apart
        # if the ODE model were ever versioned forward without
        # remembering this second copy, making every future
        # ReproducibilityBundle report the wrong model version. Assert
        # they're now the SAME string object/value, not just coincidentally
        # equal today.
        from biosim_ode.system import ODESystem
        assert MODEL_VERSION == ODESystem.MODEL_VERSION
