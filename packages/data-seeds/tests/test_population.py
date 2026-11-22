# SAVE TO: packages\data-seeds\tests\test_population.py
"""Tests for chronobiological scenarios and population simulation."""
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import pytest
import numpy as np
from biosim_seeds.chrono_scenarios import (
    morning_meal, evening_meal, meal_at_hour,
    full_day_scenarios, shift_worker_late_meal,
)
from biosim_seeds.population import PopulationSpec, PopulationRunner, CohortSummary
from biosim_api import SimulationRunner


class TestChronoScenarios:
    def test_morning_meal_creates_scenario(self):
        s = morning_meal()
        assert s.personalization.physiological_profile.time_of_day_hours == 8.0

    def test_evening_meal_creates_scenario(self):
        s = evening_meal()
        assert s.personalization.physiological_profile.time_of_day_hours == 20.0

    def test_meal_at_hour_custom(self):
        s = meal_at_hour(50.0, 14.0)
        assert s.personalization.physiological_profile.time_of_day_hours == 14.0

    def test_full_day_returns_three_scenarios(self):
        assert len(full_day_scenarios()) == 3

    def test_full_day_different_hours(self):
        scenarios = full_day_scenarios()
        hours = [s.personalization.physiological_profile.time_of_day_hours for s in scenarios]
        assert len(set(hours)) == 3

    def test_shift_worker_at_3am(self):
        assert shift_worker_late_meal().personalization.physiological_profile.time_of_day_hours == 3.0

    def test_morning_lower_glucose_than_evening(self):
        runner = SimulationRunner()
        out_m = runner.run(morning_meal(75.0))
        out_e = runner.run(evening_meal(75.0))
        assert out_m.expert_metrics.glucose_peak_mM < out_e.expert_metrics.glucose_peak_mM


class TestPopulationSpec:
    def test_default_construction(self):
        spec = PopulationSpec()
        assert spec.n_individuals == 100

    def test_custom_spec(self):
        spec = PopulationSpec(n_individuals=20, t2d_fraction=0.3)
        assert spec.t2d_fraction == 0.3

    def test_default_age_matches_census_anchor(self):
        # US Census Bureau, vintage-2023 single-year-of-age population
        # counts (ages 18-80): mean=46.50, SD=17.26 -- verified live this
        # session. Previous uncited defaults were 40.0/10.0.
        spec = PopulationSpec()
        assert spec.age_mean == pytest.approx(46.50)
        assert spec.age_sd == pytest.approx(17.26)

    def test_default_bmi_matches_nhanes_anchor(self):
        # Fryar et al. 2025 NCHS anthropometric reference data
        # (PMID:41670955), NHANES Aug 2021-Aug 2023, ages 20+, combined
        # M/F population-weighted: mean ~29.7 kg/m^2, SD ~7.0 (estimated
        # from reported percentiles). Previous uncited defaults were
        # 26.0/4.0.
        spec = PopulationSpec()
        assert spec.bmi_mean == pytest.approx(29.7)
        assert spec.bmi_sd == pytest.approx(7.0)

    def test_default_t2d_fraction_matches_prevalence(self):
        # Menke et al. 2015 JAMA (PMID:26348752), NHANES 2011-2012,
        # combined HbA1c-or-FPG-or-2-hour-PG definition: total diabetes
        # 14.3%. Previous default was 0.00 -- literally zero diabetic
        # individuals in the default cohort, not representative of any
        # real adult population.
        spec = PopulationSpec()
        assert spec.t2d_fraction == pytest.approx(0.143)

    def test_default_insulin_resistant_fraction_matches_prediabetes_prevalence(self):
        # Same source/definition as T2D above (both must share one
        # diagnostic-criteria basis for one coherent synthetic cohort):
        # prediabetes 38.0%. Previous uncited default was 0.10 -- real
        # prevalence is roughly 4x higher.
        spec = PopulationSpec()
        assert spec.insulin_resistant_fraction == pytest.approx(0.380)

    def test_male_fraction_unchanged(self):
        # Confirmed still correct (0.4946 from the same Census data) --
        # no change, verified rather than assumed.
        assert PopulationSpec().male_fraction == pytest.approx(0.50)

    def test_zero_bmi_mean_raises_clear_error(self):
        # Regression: PopulationSpec had no validation at all --
        # PopulationSpec(bmi_mean=0.0) used to construct successfully
        # and only crash later, deep inside PopulationRunner._sample_bmi()
        # (cv=bmi_sd/bmi_mean -> ZeroDivisionError), far from the actual
        # mistake. Found during the mandatory security-reviewer pass.
        with pytest.raises(ValueError, match="bmi_mean"):
            PopulationSpec(bmi_mean=0.0)

    def test_negative_bmi_mean_raises_clear_error(self):
        # Regression: used to construct successfully and only crash
        # later inside _sample_bmi()'s math.log(negative) call
        # (ValueError: math domain error) -- a different, equally
        # confusing exception type depending on the sign of one field.
        with pytest.raises(ValueError, match="bmi_mean"):
            PopulationSpec(bmi_mean=-5.0)

    @pytest.mark.parametrize("kwargs,field", [
        ({"bmi_min": 0.0}, "bmi_min"),          # below MIN_PLAUSIBLE_BMI
        ({"bmi_min": 300.0, "bmi_max": 250.0}, "bmi_min"),   # above MAX alone
        ({"bmi_max": 1000.0}, "bmi_max"),       # above MAX_PLAUSIBLE_BMI
        ({"bmi_max": 1.0, "bmi_min": 5.0}, "bmi_max"),       # below MIN alone
        ({"bmi_min": float("nan")}, "bmi_min"),
        ({"bmi_max": float("nan")}, "bmi_max"),
        ({"bmi_min": float("-inf")}, "bmi_min"),
        ({"bmi_max": float("inf")}, "bmi_max"),
        ({"bmi_mean": float("nan")}, "bmi_mean"),
        ({"bmi_mean": float("inf")}, "bmi_mean"),
        ({"bmi_sd": float("nan")}, "bmi_sd"),
        ({"bmi_sd": float("inf")}, "bmi_sd"),
    ])
    def test_non_finite_or_out_of_range_bmi_fields_raise_at_construction(self, kwargs, field):
        # Regression: these used to construct fine and fail mid-cohort,
        # inside PhysiologicalProfile's own BMI check (or, for
        # bmi_mean=inf, silently pin everyone at bmi_max).
        with pytest.raises(ValueError, match=field):
            PopulationSpec(**kwargs)

    def test_bmi_min_greater_than_bmi_max_raises(self):
        with pytest.raises(ValueError, match="bmi_min"):
            PopulationSpec(bmi_min=40.0, bmi_max=20.0)

    def test_bmi_clip_bounds_at_profile_limits_accepted(self):
        from biosim_physiology import PhysiologicalProfile
        spec = PopulationSpec(
            bmi_min=PhysiologicalProfile.MIN_PLAUSIBLE_BMI,
            bmi_max=PhysiologicalProfile.MAX_PLAUSIBLE_BMI,
        )
        assert spec.bmi_min < spec.bmi_max
        # Boundary-accepted spec must not fail mid-cohort either.
        cohort = PopulationRunner(
            PopulationSpec(n_individuals=20, bmi_min=spec.bmi_min, bmi_max=spec.bmi_max)
        ).generate_cohort()
        assert len(cohort) == 20

    def test_negative_bmi_sd_raises(self):
        with pytest.raises(ValueError, match="bmi_sd"):
            PopulationSpec(bmi_sd=-1.0)

    def test_disease_fractions_summing_over_one_raises(self):
        # Regression: t2d_fraction/insulin_resistant_fraction were
        # completely unbounded -- a sum exceeding 1.0 didn't crash, but
        # silently made HEALTHY unreachable and produced a cohort with
        # unintended proportions, not what the caller asked for. This
        # became a materially bigger risk once the defaults themselves
        # moved from 0.00+0.10=0.10 (0.90 of headroom before an override
        # could exceed 1.0) to the real-prevalence-anchored 0.143+0.380
        # =0.523 (only 0.477 of headroom).
        with pytest.raises(ValueError, match="exceeds"):
            PopulationSpec(t2d_fraction=0.7, insulin_resistant_fraction=0.5)

    def test_disease_fraction_out_of_range_raises(self):
        with pytest.raises(ValueError, match="t2d_fraction"):
            PopulationSpec(t2d_fraction=1.5)
        with pytest.raises(ValueError, match="insulin_resistant_fraction"):
            PopulationSpec(insulin_resistant_fraction=-0.1)

    def test_disease_fractions_summing_to_exactly_one_is_fine(self):
        # 1.0 total is a legitimate (if extreme) edge case -- HEALTHY
        # becomes unreachable, which is exactly what the caller asked
        # for, not an error.
        spec = PopulationSpec(t2d_fraction=0.6, insulin_resistant_fraction=0.4)
        assert spec.t2d_fraction + spec.insulin_resistant_fraction == pytest.approx(1.0)


class TestPopulationRunnerCohort:
    @pytest.fixture(scope="class")
    @classmethod
    def small_spec(cls):
        return PopulationSpec(n_individuals=10, random_seed=42)

    def test_generates_correct_count(self, small_spec):
        assert len(PopulationRunner(small_spec).generate_cohort()) == 10

    def test_age_within_bounds(self, small_spec):
        for p in PopulationRunner(small_spec).generate_cohort():
            assert small_spec.age_min <= p.age_years <= small_spec.age_max

    def test_bmi_within_bounds(self, small_spec):
        for p in PopulationRunner(small_spec).generate_cohort():
            assert small_spec.bmi_min <= p.bmi <= small_spec.bmi_max

    def test_same_seed_same_cohort(self):
        spec = PopulationSpec(n_individuals=10, random_seed=99)
        p1 = PopulationRunner(spec).generate_cohort()
        p2 = PopulationRunner(spec).generate_cohort()
        assert [round(p.age_years, 2) for p in p1] == [round(p.age_years, 2) for p in p2]

    def test_bmi_sampled_log_normal_is_right_skewed(self):
        # Regression: BMI used to be sampled from a symmetric Normal
        # (skew ~0), which understates the real right-skewed shape of
        # adult BMI (a long upper tail into obesity) -- confirmed this
        # session from NHANES's own reported percentiles (see
        # PopulationSpec's docstring). Now sampled log-normal, which
        # must produce genuine positive skew, not just a wider spread.
        spec = PopulationSpec(n_individuals=3000, random_seed=42)
        bmis = np.array([p.bmi for p in PopulationRunner(spec).generate_cohort()])
        skew = np.mean((bmis - bmis.mean()) ** 3) / bmis.std() ** 3
        assert skew > 0.2   # a symmetric Normal would give ~0

    def test_bmi_sample_mean_and_sd_match_spec_reasonably(self):
        # Log-normal is parameterised so its arithmetic mean/SD match
        # PopulationSpec's own bmi_mean/bmi_sd (before clipping) -- a
        # large, mostly-unclipped sample should land close to both.
        spec = PopulationSpec(n_individuals=3000, random_seed=42,
                              bmi_min=5.0, bmi_max=250.0)   # effectively unclipped (profile limits)
        bmis = np.array([p.bmi for p in PopulationRunner(spec).generate_cohort()])
        assert bmis.mean() == pytest.approx(spec.bmi_mean, rel=0.05)
        assert bmis.std() == pytest.approx(spec.bmi_sd, rel=0.10)


class TestPopulationRunnerSimulation:
    @pytest.fixture(scope="class")
    @classmethod
    def small_summary(cls):
        from biosim_api import get_scenario
        spec = PopulationSpec(n_individuals=6, random_seed=42,
                              t2d_fraction=0.0, insulin_resistant_fraction=0.0)
        return PopulationRunner(spec).run(get_scenario("standard_meal"))

    def test_run_produces_summary(self, small_summary):
        assert isinstance(small_summary, CohortSummary)

    def test_all_succeeded(self, small_summary):
        assert small_summary.n_succeeded == small_summary.n_total

    def test_mean_glucose_physiological(self, small_summary):
        assert 5.0 < small_summary.glucose_peak_mean_mM < 20.0

    def test_sd_positive(self, small_summary):
        assert small_summary.glucose_peak_sd_mM > 0

    def test_individual_peaks_match_count(self, small_summary):
        assert len(small_summary.individual_peaks_mM) == small_summary.n_succeeded

    def test_summary_dict(self, small_summary):
        d = small_summary.summary_dict()
        assert "glucose_peak_mM" in d
        assert "pct_above_7p8_mM" in d

    def test_t2d_cohort_higher_peak(self):
        from biosim_api import get_scenario
        scenario = get_scenario("standard_meal")
        healthy = PopulationRunner(PopulationSpec(n_individuals=5, t2d_fraction=0.0,
                                                  insulin_resistant_fraction=0.0, random_seed=42))
        t2d     = PopulationRunner(PopulationSpec(n_individuals=5, t2d_fraction=1.0,
                                                  insulin_resistant_fraction=0.0, random_seed=42))
        assert t2d.run(scenario).glucose_peak_mean_mM > healthy.run(scenario).glucose_peak_mean_mM


class TestGeneticVariantSampling:
    """
    Regression: PopulationSpec.sample_genetic_variants (opt-in, default
    False) drives Hardy-Weinberg genotype sampling at each of the 7
    default T2D-risk variants with a verified population_frequency, then
    threads the sampled Variant list into PersonalizationSpec.genetic_variants
    so it actually reaches GenotypePersonalizer -- previously this data
    existed (population_frequency added this session) but nothing sampled
    or wired it, the same "computed but never wired" bug class CLAUDE.md
    documents repeatedly for this codebase.
    """

    def test_dosage_distribution_matches_hardy_weinberg(self):
        # A large number of draws at a known frequency should match
        # P(0)=(1-q)^2, P(1)=2q(1-q), P(2)=q^2 within sampling noise.
        np_rng = np.random.default_rng(42)
        q = 0.30
        n = 20000
        dosages = [PopulationRunner._sample_genotype_dosage(np_rng, q) for _ in range(n)]
        counts = np.bincount(dosages, minlength=3) / n
        expected = np.array([(1 - q) ** 2, 2 * q * (1 - q), q ** 2])
        assert np.allclose(counts, expected, atol=0.02)

    def test_sampled_variants_cover_all_seven_default_frequency_verified_rsids(self):
        from biosim_genotype import VariantEffectDB
        spec = PopulationSpec(n_individuals=1, sample_genetic_variants=True, random_seed=1)
        np_rng = np.random.default_rng(spec.random_seed + 1)
        variants = PopulationRunner(spec)._sample_genetic_variants_for_individual(np_rng)
        expected_rsids = {
            rsid for rsid in VariantEffectDB.default().all_rsids()
            if VariantEffectDB.default().get(rsid)[0].population_frequency is not None
        }
        assert {v.rsid for v in variants} == expected_rsids
        assert len(variants) == 7

    def test_opt_out_by_default_no_genetic_variants_sampled(self, monkeypatch):
        from biosim_api import get_scenario, PersonalizationSpec

        captured = {}
        real_init = PersonalizationSpec.__init__

        def _capture(self, *args, **kwargs):
            real_init(self, *args, **kwargs)
            captured["genetic_variants"] = list(self.genetic_variants)

        monkeypatch.setattr(PersonalizationSpec, "__init__", _capture)
        spec = PopulationSpec(n_individuals=2, random_seed=42,
                              t2d_fraction=0.0, insulin_resistant_fraction=0.0)
        PopulationRunner(spec).run(get_scenario("standard_meal"))
        assert captured["genetic_variants"] == []

    def test_opt_in_threads_sampled_variants_into_personalization(self, monkeypatch):
        from biosim_api import get_scenario, PersonalizationSpec

        captured = []
        real_init = PersonalizationSpec.__init__

        def _capture(self, *args, **kwargs):
            real_init(self, *args, **kwargs)
            captured.append(list(self.genetic_variants))

        monkeypatch.setattr(PersonalizationSpec, "__init__", _capture)
        spec = PopulationSpec(n_individuals=3, random_seed=42, sample_genetic_variants=True,
                              t2d_fraction=0.0, insulin_resistant_fraction=0.0)
        PopulationRunner(spec).run(get_scenario("standard_meal"))
        assert len(captured) == 3
        assert all(len(variants) == 7 for variants in captured)

    def test_genetically_personalized_outputs_pass_privacy_audit(self, monkeypatch):
        # Rule 7 / security-review regression: this is the first caller
        # to thread real, cohort-scale genetic variant data through
        # PersonalizationSpec via biosim_seeds.population -- previously
        # untested against PrivacyAuditor at this scale. Captures the
        # real per-individual SimulationOutput objects (via a
        # side_effect wrapper around the real SimulationRunner.run, so
        # the actual ODE solve still happens and results stay genuine)
        # and asserts every one passes PrivacyAuditor().audit().
        from biosim_api import get_scenario, SimulationRunner
        from biosim_infra import PrivacyAuditor

        captured_outputs = []
        real_run = SimulationRunner.run

        def _capture(self, scenario):
            output = real_run(self, scenario)
            captured_outputs.append(output)
            return output

        monkeypatch.setattr(SimulationRunner, "run", _capture)
        spec = PopulationSpec(n_individuals=3, random_seed=42, sample_genetic_variants=True,
                              t2d_fraction=0.0, insulin_resistant_fraction=0.0)
        PopulationRunner(spec).run(get_scenario("standard_meal"))

        assert len(captured_outputs) == 3
        auditor = PrivacyAuditor()
        for output in captured_outputs:
            result = auditor.audit(output)
            assert result.passed, result.issues

    def test_exported_cohort_csv_and_json_contain_no_raw_rsids(self):
        # Security-review regression: CohortSummary.failure_reasons and
        # population_export.py's CSV/JSON are the cohort-level surfaces
        # a caller would actually read/export -- confirm no rsid-shaped
        # string ever appears in any of them for a genetically-
        # personalized cohort run.
        import re
        from biosim_api import get_scenario
        from biosim_seeds.population_export import summary_csv, individuals_csv, to_json

        spec = PopulationSpec(n_individuals=4, random_seed=42, sample_genetic_variants=True,
                              t2d_fraction=0.0, insulin_resistant_fraction=0.0)
        summary = PopulationRunner(spec).run(get_scenario("standard_meal"))

        rsid_pattern = re.compile(r'\brs\d{4,}\b')
        for exported in (summary_csv(summary), individuals_csv(summary), to_json(summary)):
            assert not rsid_pattern.search(exported), exported

    def test_reproducible_with_same_seed(self, monkeypatch):
        from biosim_api import get_scenario, PersonalizationSpec

        def _sampled_rsids_per_individual(seed):
            captured = []
            real_init = PersonalizationSpec.__init__

            def _capture(self, *args, **kwargs):
                real_init(self, *args, **kwargs)
                captured.append(tuple(sorted(
                    (v.rsid, v.genotype.value) for v in self.genetic_variants
                )))

            monkeypatch.setattr(PersonalizationSpec, "__init__", _capture)
            spec = PopulationSpec(n_individuals=4, random_seed=seed,
                                  sample_genetic_variants=True,
                                  t2d_fraction=0.0, insulin_resistant_fraction=0.0)
            PopulationRunner(spec).run(get_scenario("standard_meal"))
            monkeypatch.setattr(PersonalizationSpec, "__init__", real_init)
            return captured

        first  = _sampled_rsids_per_individual(seed=7)
        second = _sampled_rsids_per_individual(seed=7)
        assert first == second


class TestPopulationRunnerFailureHandling:
    """
    Regression: individual-profile failures during a cohort run used to be
    silently discarded (bare `except Exception: pass`, flagged by bandit
    B110 -- try/except/pass with no error context). A single bad profile
    (or a real bug affecting many profiles) was indistinguishable from a
    healthy cohort with fewer participants -- n_succeeded < n_total was
    visible, but WHY was not. Fixed by recording a reason per failure in
    CohortSummary.failure_reasons instead of discarding it.
    """

    def test_exception_during_simulation_is_recorded_not_swallowed(self, monkeypatch):
        from biosim_api import get_scenario

        def _raise(self, scenario):
            raise ValueError("synthetic ODE solver failure for testing")

        monkeypatch.setattr(SimulationRunner, "run", _raise)
        spec = PopulationSpec(n_individuals=3, random_seed=42)
        summary = PopulationRunner(spec).run(get_scenario("standard_meal"))

        assert summary.n_succeeded == 0
        assert summary.n_failed == 3
        assert len(summary.failure_reasons) == 3
        assert all("ValueError" in reason and "synthetic ODE solver failure" in reason
                   for reason in summary.failure_reasons)

    def test_unsuccessful_output_is_recorded_not_silently_dropped(self, monkeypatch):
        from biosim_api import get_scenario

        class _FakeOutput:
            success = False
            warnings = ["ODE solver failed: synthetic non-convergence"]

        monkeypatch.setattr(SimulationRunner, "run", lambda self, scenario: _FakeOutput())
        spec = PopulationSpec(n_individuals=2, random_seed=42)
        summary = PopulationRunner(spec).run(get_scenario("standard_meal"))

        assert summary.n_succeeded == 0
        assert len(summary.failure_reasons) == 2
        assert all("synthetic non-convergence" in reason for reason in summary.failure_reasons)

    def test_healthy_cohort_has_no_failure_reasons(self):
        from biosim_api import get_scenario
        spec = PopulationSpec(n_individuals=3, random_seed=42,
                              t2d_fraction=0.0, insulin_resistant_fraction=0.0)
        summary = PopulationRunner(spec).run(get_scenario("standard_meal"))
        assert summary.failure_reasons == []


class TestPopulationRunnerMultiprocessing:
    """
    Track 3 follow-on to BatchRunner's ProcessPoolExecutor fix.
    use_multiprocessing defaults to False here (unlike BatchRunner's
    new default) -- benchmarked directly before choosing this: at
    PopulationSpec's own default n_individuals=100, 8-worker
    ProcessPoolExecutor measured SLOWER than sequential (2.326s vs.
    2.222s, small-n process-spawn overhead not yet amortized); at
    n=400 it measured a genuine 1.73x speedup. See PopulationRunner
    .__init__'s own docstring for the full numbers.
    """

    def test_multiprocessing_defaults_to_false(self):
        assert PopulationRunner().use_multiprocessing is False

    def test_multiprocessing_produces_same_aggregate_stats_as_sequential(self):
        # Cohort generation and genetic-variant sampling are both
        # explicitly documented to stay sequential regardless of
        # use_multiprocessing (RNG-stream reproducibility) -- and each
        # individual's ODE simulation is itself a deterministic pure
        # function of its input PathwayParameters, so switching HOW the
        # simulations are executed must not change the result at all,
        # only how long it takes.
        from biosim_api import get_scenario
        scenario = get_scenario("standard_meal")

        seq_spec = PopulationSpec(n_individuals=6, random_seed=42,
                                  t2d_fraction=0.0, insulin_resistant_fraction=0.0)
        par_spec = PopulationSpec(n_individuals=6, random_seed=42,
                                  t2d_fraction=0.0, insulin_resistant_fraction=0.0)

        seq_summary = PopulationRunner(seq_spec, use_multiprocessing=False).run(scenario)
        par_summary = PopulationRunner(par_spec, use_multiprocessing=True).run(scenario)

        assert par_summary.n_succeeded == seq_summary.n_succeeded
        assert par_summary.glucose_peak_mean_mM == pytest.approx(seq_summary.glucose_peak_mean_mM)
        assert par_summary.individual_peaks_mM == pytest.approx(seq_summary.individual_peaks_mM)

    def test_run_all_captures_real_failure_across_process_boundary(self):
        # A genuine, non-monkeypatched failure -- a types.SimpleNamespace
        # stand-in scenario (has no attribute SimulationRunner.run()
        # needs) fails naturally once actually simulated inside the
        # worker, same technique as BatchRunner's equivalent test.
        # Exercises _run_all() directly (the unit that actually
        # branches on use_multiprocessing) rather than the full run()
        # pipeline, since a genuinely-broken *scenario* is easy to
        # construct but a genuinely-broken *PhysiologicalProfile* that
        # still survives cohort generation is not.
        import types
        from biosim_api import get_scenario
        good = get_scenario("standard_meal")
        bad = types.SimpleNamespace(name="genuinely_invalid")

        runner = PopulationRunner(PopulationSpec(n_individuals=1), use_multiprocessing=True)
        results = runner._run_all([good, bad])

        assert results[0][1] is None and results[0][0] is not None
        assert results[1][0] is None and isinstance(results[1][1], AttributeError)

    def test_failure_reason_recorded_with_multiprocessing_and_thread_pool(self, monkeypatch):
        # Mirrors TestPopulationRunnerFailureHandling's monkeypatch-based
        # failure test, but through the use_multiprocessing=True code
        # path -- must pass executor_cls=ThreadPoolExecutor explicitly,
        # since a spawned (ProcessPoolExecutor) worker never sees this
        # process's monkeypatched SimulationRunner.run.
        from biosim_api import get_scenario

        def _raise(self, scenario):
            raise ValueError("synthetic ODE solver failure for testing")

        monkeypatch.setattr(SimulationRunner, "run", _raise)
        spec = PopulationSpec(n_individuals=3, random_seed=42)
        runner = PopulationRunner(spec, use_multiprocessing=True, executor_cls=ThreadPoolExecutor)
        summary = runner.run(get_scenario("standard_meal"))

        assert summary.n_succeeded == 0
        assert summary.n_failed == 3
        assert len(summary.failure_reasons) == 3
        assert all("ValueError" in reason and "synthetic ODE solver failure" in reason
                   for reason in summary.failure_reasons)

    def test_process_pool_executor_actually_runs_a_full_cohort_correctly(self):
        # End-to-end confirmation through the real default executor
        # (ProcessPoolExecutor) that pickling ScenarioDefinition/
        # PersonalizationSpec/PhysiologicalProfile and the resulting
        # SimulationOutput across the process boundary doesn't corrupt
        # or drop data for a genuinely heterogeneous cohort (mixed
        # disease states, not just healthy).
        from biosim_api import get_scenario
        spec = PopulationSpec(n_individuals=8, random_seed=1)
        summary = PopulationRunner(spec, use_multiprocessing=True).run(get_scenario("standard_meal"))
        assert summary.n_succeeded == 8
        assert summary.glucose_peak_mean_mM > 0

    def test_genetic_variant_sampling_and_multiprocessing_together_match_sequential(self):
        # Code-reviewer finding: the two flags interact in a way neither
        # of the flags' own individual test suites covers on its own --
        # sample_genetic_variants=True means each individual's
        # ScenarioDefinition carries a PersonalizationSpec.genetic_variants
        # list of real Variant objects (rsid + Genotype enum), which must
        # ALSO survive pickling across the ProcessPoolExecutor process
        # boundary alongside everything else, not just the plain-
        # demographic case test_process_pool_executor_actually_runs_a_
        # full_cohort_correctly already covers. Confirms sequential and
        # parallel execution of the SAME genetically-sampled cohort
        # (same seed) produce identical aggregate stats -- if pickling
        # silently dropped or corrupted genetic_variants in one path but
        # not the other, the two arms would diverge (genotype affects
        # secretion_scale/p3/Gb, which measurably shifts glucose_peak_mM).
        from biosim_api import get_scenario
        scenario = get_scenario("standard_meal")

        seq_spec = PopulationSpec(n_individuals=6, random_seed=42, sample_genetic_variants=True,
                                  t2d_fraction=0.0, insulin_resistant_fraction=0.0)
        par_spec = PopulationSpec(n_individuals=6, random_seed=42, sample_genetic_variants=True,
                                  t2d_fraction=0.0, insulin_resistant_fraction=0.0)

        seq_summary = PopulationRunner(seq_spec, use_multiprocessing=False).run(scenario)
        par_summary = PopulationRunner(par_spec, use_multiprocessing=True).run(scenario)

        assert par_summary.n_succeeded == seq_summary.n_succeeded == 6
        assert par_summary.individual_peaks_mM == pytest.approx(seq_summary.individual_peaks_mM)
