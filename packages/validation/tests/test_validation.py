# SAVE TO: packages\validation\tests\test_validation.py
"""
Test suite for biosim_validation.

Coverage:
    - BenchmarkPoint: tolerance bands, check(), pass/fail
    - BenchmarkDataset: construction, add
    - ModelValidator: extract values, validate OGTT, verdict logic
    - ValidationReport: summary, as_text, failed_benchmarks
    - Pre-registered datasets: correct structure
"""

import pytest
from biosim_validation import (
    BenchmarkCheckResult, BenchmarkDataset, BenchmarkPoint,
    ModelValidator, Tolerance, ValidationReport, ValidationVerdict,
    ogtt_75g_healthy_adults,
)
from biosim_api import SimulationRunner, high_sugar_meal_scenario, fasted_scenario


# ---------------------------------------------------------------------------
# BenchmarkPoint
# ---------------------------------------------------------------------------

class TestBenchmarkPoint:
    @pytest.fixture
    def point(self):
        return BenchmarkPoint(
            name="test_glucose", description="Test", reference=7.8,
            unit="mM", tolerance=Tolerance.NORMAL, pubmed_id="12345",
        )

    def test_tolerance_fraction_normal(self, point):
        assert point.tolerance_fraction == pytest.approx(0.25)

    def test_tolerance_fraction_tight(self):
        p = BenchmarkPoint("x","x",5.0,"mM",Tolerance.TIGHT,"0")
        assert p.tolerance_fraction == pytest.approx(0.10)

    def test_tolerance_fraction_wide(self):
        p = BenchmarkPoint("x","x",5.0,"mM",Tolerance.WIDE,"0")
        assert p.tolerance_fraction == pytest.approx(0.50)

    def test_lower_bound(self, point):
        assert point.lower_bound == pytest.approx(7.8 * 0.75)

    def test_upper_bound(self, point):
        assert point.upper_bound == pytest.approx(7.8 * 1.25)

    def test_check_pass_within_band(self, point):
        result = point.check(8.0)
        assert result.passed

    def test_check_fail_above_band(self, point):
        result = point.check(12.0)
        assert not result.passed

    def test_check_fail_below_band(self, point):
        result = point.check(1.0)
        assert not result.passed

    def test_check_at_reference(self, point):
        result = point.check(7.8)
        assert result.passed

    def test_deviation_percent_positive(self, point):
        result = point.check(9.0)
        assert result.deviation_percent > 0

    def test_deviation_percent_negative(self, point):
        result = point.check(6.0)
        assert result.deviation_percent < 0

    def test_status_string(self, point):
        assert point.check(7.8).status == "PASS"
        assert point.check(99.0).status == "MISS"


# ---------------------------------------------------------------------------
# BenchmarkDataset
# ---------------------------------------------------------------------------

class TestBenchmarkDataset:
    def test_construction(self):
        ds = BenchmarkDataset(name="test", description="A test dataset")
        assert len(ds) == 0

    def test_add_point(self):
        ds = BenchmarkDataset(name="test", description="A test dataset")
        ds.add(BenchmarkPoint("p1","d",5.0,"mM",Tolerance.NORMAL,"0"))
        assert len(ds) == 1

    def test_ogtt_dataset_has_points(self):
        ds = ogtt_75g_healthy_adults()
        assert len(ds) >= 5

    def test_ogtt_dataset_has_fasting_glucose(self):
        ds = ogtt_75g_healthy_adults()
        names = [p.name for p in ds.points]
        assert "fasting_glucose" in names

    def test_ogtt_dataset_has_peak_glucose(self):
        ds = ogtt_75g_healthy_adults()
        names = [p.name for p in ds.points]
        assert "peak_glucose" in names

    def test_ogtt_dataset_has_glucose_at_120min(self):
        ds = ogtt_75g_healthy_adults()
        names = [p.name for p in ds.points]
        assert "glucose_at_120min" in names

    def test_all_reference_values_positive(self):
        ds = ogtt_75g_healthy_adults()
        for p in ds.points:
            assert p.reference > 0


# ---------------------------------------------------------------------------
# ModelValidator
# ---------------------------------------------------------------------------

class TestModelValidator:
    @pytest.fixture(scope="class")
    def ogtt_output(self):
        from biosim_api import ogtt_75g
        return SimulationRunner().run(ogtt_75g())

    @pytest.fixture(scope="class")
    def ogtt_report(self, ogtt_output):
        ds = ogtt_75g_healthy_adults()
        return ModelValidator().validate(ds, ogtt_output)

    def test_validate_returns_report(self, ogtt_report):
        assert isinstance(ogtt_report, ValidationReport)

    def test_has_results(self, ogtt_report):
        assert len(ogtt_report.results) > 0

    def test_verdict_is_valid(self, ogtt_report):
        assert ogtt_report.verdict in list(ValidationVerdict)

    def test_fasting_glucose_checked(self, ogtt_report):
        names = [r.name for r in ogtt_report.results]
        assert "fasting_glucose" in names

    def test_peak_glucose_checked(self, ogtt_report):
        names = [r.name for r in ogtt_report.results]
        assert "peak_glucose" in names

    def test_fasting_glucose_passes(self, ogtt_report):
        fg = next(r for r in ogtt_report.results if r.name == "fasting_glucose")
        assert fg.passed, f"Fasting glucose {fg.model_value:.2f} mM outside tolerance"

    def test_pass_rate_between_0_and_1(self, ogtt_report):
        assert 0.0 <= ogtt_report.pass_rate <= 1.0

    def test_n_pass_plus_n_fail_equals_total(self, ogtt_report):
        assert ogtt_report.n_pass + ogtt_report.n_fail == len(ogtt_report.results)

    def test_summary_structure(self, ogtt_report):
        s = ogtt_report.summary()
        assert "verdict" in s
        assert "n_pass" in s
        assert "n_fail" in s
        assert "pass_rate" in s

    def test_as_text_contains_verdict(self, ogtt_report):
        text = ogtt_report.as_text()
        assert ogtt_report.verdict.value in text

    def test_as_text_contains_benchmark_names(self, ogtt_report):
        text = ogtt_report.as_text()
        assert "fasting_glucose" in text

    def test_validate_ogtt_convenience(self):
        report = ModelValidator().validate_ogtt()
        assert isinstance(report, ValidationReport)
        # KNOWN FAIL since 2026-11-16: the in-house OGTT's only scored
        # NORMAL point (peak glucose) is +33% after B7 removed the
        # double-counted clearance its absorption/first-phase constants were
        # fitted around (circular benchmark, NOT re-fitted); the 2-hour point
        # became unscored in B4 (wrong citation) and no longer dilutes it.
        assert report.verdict == ValidationVerdict.FAIL, report.as_text()
        failed = [r.name for r in report.results if r.benchmark.is_scored and not r.passed]
        assert failed == ["peak_glucose"], failed

    def test_ogtt_glucose_timing_within_tight_bound(self, ogtt_report):
        # Regression for the too-fast-absorption bug (CLAUDE.md VALIDATION
        # STATUS): before GutAbsorptionModel.ogtt_calibrated_glucose() was
        # wired in, time_to_peak_glucose_min deviated -46.5% (24.1 vs. 45.0
        # min ref) -- right at the edge of even WIDE tolerance. Assert a
        # much tighter bound than the benchmark's own WIDE (50%) tolerance
        # to actually catch a regression back to that too-fast timescale,
        # not just catch it if it drifts all the way back to failing.
        #
        # B7 (2026-11-16): removing the double-counted resting clearance made
        # the peak LATER (+17.6%), the opposite direction from the regression
        # this guards against. The absorption constants were fitted while the
        # double count was in place, so the circular OGTT fit no longer holds;
        # they are NOT re-fitted. Guard the old direction, pin the new value.
        t_peak = next(r for r in ogtt_report.results if r.name == "time_to_peak_glucose_min")
        assert t_peak.deviation_percent == pytest.approx(17.6, abs=1.0)
        assert t_peak.deviation_percent > -15.0, (
            f"time_to_peak_glucose_min deviation {t_peak.deviation_percent:+.1f}% -- "
            f"regression toward the too-fast monoexponential absorption timing?"
        )

    def test_ogtt_peak_insulin_within_bound(self, ogtt_report):
        # Regression for the uncited KG_secretion=2.0 default (CLAUDE.md
        # VALIDATION STATUS): before it was recalibrated to 5.63 (anchored
        # to Rebelato et al. 2018, PMID:30166558, GSIS EC50=10.13 mM in
        # control islets), peak_insulin deviated -44.9%. Assert a tighter
        # bound than the benchmark's own WIDE (50%) tolerance to catch a
        # regression back toward that undershoot.
        #
        # B7 (2026-11-16): peak insulin now OVERSHOOTS (+37.4%), the opposite
        # direction from the regression this guards against: first_phase_gain
        # was fitted to this OGTT while the double count was in place, and is
        # NOT re-fitted. Guard the old direction, pin the new value.
        peak_i = next(r for r in ogtt_report.results if r.name == "peak_insulin")
        assert peak_i.deviation_percent == pytest.approx(37.4, abs=1.0)
        assert peak_i.deviation_percent > -20.0, (
            f"peak_insulin deviation {peak_i.deviation_percent:+.1f}% -- "
            f"regression toward the uncited KG_secretion=2.0 undershoot, or the "
            f"since-corrected p3, or a lost first_phase_gain contribution?"
        )

    def test_ogtt_insulin_timing_within_tight_bound(self, ogtt_report):
        # Regression for KNOWN BUGS FIXED #15/first-phase secretion: before
        # PathwayParameters.first_phase_gain existed, time_to_peak_insulin_min
        # deviated +41.7% even with the corrected p3. Assert a tight bound
        # so a regression (e.g. first_phase_gain reset to 0) is caught
        # immediately rather than only showing up as a WIDE-tolerance edge case.
        t_peak_i = next(r for r in ogtt_report.results if r.name == "time_to_peak_insulin_min")
        assert abs(t_peak_i.deviation_percent) < 15.0, (
            f"time_to_peak_insulin_min deviation {t_peak_i.deviation_percent:+.1f}% -- "
            f"regression toward missing first-phase secretion?"
        )

    def test_mtor_activation_is_marked_unscored(self, ogtt_report):
        # mtor_activation's own notes admit "no direct clinical comparator
        # exists" -- pubmed_id="UNVERIFIED" should make is_scored False so
        # it can't quietly read as one of "n" genuine validation checks.
        mtor = next(r for r in ogtt_report.results if r.name == "mtor_activation")
        assert mtor.benchmark.is_scored is False

    def test_scored_benchmarks_have_a_real_pubmed_id(self, ogtt_report):
        # Every OTHER OGTT point must have a real citation -- if this ever
        # fails, either a new UNVERIFIED point was added silently (update
        # this test deliberately) or a real citation was accidentally
        # replaced with the UNVERIFIED marker.
        # glucose_at_120min joined the unscored points in B4 (2026-11-16):
        # its only citation was the wrong paper (DeFronzo 1979, PMID:382871).
        unscored = {"mtor_activation", "glucose_at_120min"}
        assert all(not r.benchmark.is_scored for r in ogtt_report.results if r.name in unscored)
        scored = [r for r in ogtt_report.results if r.name not in unscored]
        assert scored, "expected at least one scored benchmark besides mtor_activation"
        for r in scored:
            assert r.benchmark.is_scored is True
            assert r.benchmark.pubmed_id != "UNVERIFIED"

    def test_report_notes_unscored_point_explicitly(self, ogtt_report):
        # _unscored_point_notes() must flag mtor_activation in every
        # report, regardless of whether it happens to pass or fail --
        # otherwise "(7/7 pass)" reads as 7 genuine validations.
        assert any("no real literature anchor" in n for n in ogtt_report.notes)
        assert any("mtor_activation" in n for n in ogtt_report.notes)

    def test_as_text_labels_unscored_point_distinctly(self, ogtt_report):
        text = ogtt_report.as_text()
        mtor_line = next(l for l in text.splitlines() if "mtor_activation" in l and l.strip().startswith("["))
        assert "/INFO" in mtor_line


class TestValidateBergmanIV:
    """
    See CLAUDE.md KNOWN BUGS FIXED #15: this benchmark used to fail
    honestly (glucose_return_time_min deviated -57%) until PathwayParameters
    .p3 was corrected to match Bergman 1981's own directly-quoted SI value
    (Fig 6 caption, PMID:7033284) -- it now genuinely passes.
    """

    @pytest.fixture(scope="class")
    def iv_report(self):
        return ModelValidator().validate_bergman_iv()

    def test_returns_report(self, iv_report):
        assert isinstance(iv_report, ValidationReport)

    def test_has_one_result(self, iv_report):
        assert len(iv_report.results) == 1
        assert iv_report.results[0].name == "glucose_return_time_min"

    def test_verdict_is_pass(self, iv_report):
        assert iv_report.verdict == ValidationVerdict.PASS

    def test_glucose_return_time_actually_passes(self, iv_report):
        # Regression for KNOWN BUGS FIXED #15: this must be a genuine pass
        # (the benchmark point itself, not just the WIDE-failure-tolerant
        # verdict) -- a regression back toward the old p3=5.03e-6 would
        # fail this even though the overall verdict might still say PASS.
        result = iv_report.results[0]
        assert result.passed, (
            f"glucose_return_time_min deviation {result.deviation_percent:+.1f}% -- "
            f"regression toward the old, too-high p3?"
        )


class TestCalibrateOGTT:
    """
    Bayesian calibration of p1/p3/KG_secretion against the OGTT benchmark
    (Sampling-Importance-Resampling). n_prior_samples kept small for test
    speed, matching biosim_uncertainty's own TestMonteCarlo convention.
    """

    @pytest.fixture(scope="class")
    def calibration(self):
        from biosim_uncertainty import CalibrationResult
        result = ModelValidator().calibrate_ogtt(n_prior_samples=40)
        assert isinstance(result, CalibrationResult)
        return result

    def test_returns_calibration_result(self, calibration):
        from biosim_uncertainty import CalibrationResult
        assert isinstance(calibration, CalibrationResult)

    def test_calibrates_p1_p3_kg_secretion(self, calibration):
        assert set(calibration.posterior_parameter_samples.keys()) == {
            "p1", "p3", "KG_secretion",
        }

    def test_mtor_activation_excluded_as_unscored(self, calibration):
        # ogtt_75g_healthy_adults()'s mtor_activation point has
        # pubmed_id="UNVERIFIED" (is_scored=False) -- must never
        # contribute to the likelihood.
        assert "mtor_activation" not in calibration.targets_used

    def test_not_degenerate(self, calibration):
        # A genuine miscalibration (bad prior, bad likelihood_sigma_z
        # choice) would show up as a collapsed effective sample size --
        # this is the single most important health check on this
        # specific calibration target, run for real rather than assumed.
        assert calibration.is_degenerate is False

    def test_posterior_parameters_stay_positive(self, calibration):
        for values in calibration.posterior_parameter_samples.values():
            assert (values > 0).all()

    def test_reproducible_with_same_seed(self):
        r1 = ModelValidator().calibrate_ogtt(n_prior_samples=30, random_seed=7)
        r2 = ModelValidator().calibrate_ogtt(n_prior_samples=30, random_seed=7)
        assert (r1.posterior_parameter_samples["p1"] == r2.posterior_parameter_samples["p1"]).all()


class TestPosteriorPredictiveCheckBergmanIV:
    """
    The statistical-honesty safeguard BayesianCalibrator's own module
    docstring calls for: check the OGTT-calibrated posterior against a
    DIFFERENT benchmark (Bergman IV) that contributed zero weight to the
    likelihood, rather than declaring the posterior correct on the
    strength of the benchmark it was fit to.
    """

    @pytest.fixture(scope="class")
    def pp_report(self):
        calibration = ModelValidator().calibrate_ogtt(n_prior_samples=30)
        return ModelValidator().posterior_predictive_check_bergman_iv(calibration)

    def test_returns_validation_report(self, pp_report):
        assert isinstance(pp_report, ValidationReport)

    def test_scenario_name_marks_this_as_a_posterior_predictive_check(self, pp_report):
        # Must read differently from validate_bergman_iv()'s own
        # scenario_name -- this is NOT the same kind of check (every
        # posterior draw re-simulated, not one point-estimate run) and a
        # caller reading the report text should not be misled.
        assert pp_report.scenario_name == "bergman_iv_posterior_predictive_check"

    def test_has_one_result(self, pp_report):
        assert len(pp_report.results) == 1
        assert pp_report.results[0].name == "glucose_return_time_min"

    def test_notes_document_which_draws_were_used(self, pp_report):
        assert any("posterior parameter draws" in n for n in pp_report.notes)

    def test_all_draws_produce_a_result_for_a_healthy_calibration(self, pp_report):
        # Every posterior draw is a perturbation of literature-anchored
        # defaults -- none should fail to solve or fail to return to
        # baseline within the 14400s window.
        assert any("30/30 posterior draws produced a valid return time" in n for n in pp_report.notes)

    def test_one_raising_posterior_draw_does_not_crash_the_whole_check(self, monkeypatch):
        # Regression: the per-draw loop used to have no exception
        # handling at all, unlike every other per-sample Monte Carlo
        # loop in this codebase (MonteCarloEngine.propagate()) -- found
        # during the mandatory code-reviewer/security-reviewer passes,
        # see CLAUDE.md. Re-simulating a posterior draw under a
        # DIFFERENT scenario (IV bolus) than the one it succeeded under
        # during calibration (OGTT) is not guaranteed to succeed just
        # because the OGTT run did -- one raising draw used to crash
        # this entire posterior-predictive check.
        import biosim_ode.system as system_module

        calibration = ModelValidator().calibrate_ogtt(n_prior_samples=20)

        call_count = {"n": 0}
        real_simulate = system_module.ODESystem.simulate

        def _raise_on_third_call(self, duration_s, meals=None, n_output_points=500, y0_override=None):
            call_count["n"] += 1
            if call_count["n"] == 3:
                raise ValueError("synthetic solver failure for testing")
            return real_simulate(self, duration_s, meals, n_output_points, y0_override)

        monkeypatch.setattr(system_module.ODESystem, "simulate", _raise_on_third_call)
        report = ModelValidator().posterior_predictive_check_bergman_iv(calibration)

        assert isinstance(report, ValidationReport)
        assert any("1/" in n and "raised an exception" in n for n in report.notes)


class TestValidateSemaglutide:
    """
    See CLAUDE.md "DONE -- Built the first validation benchmark for a clinically significant drug pathway": first validation benchmark for a
    clinically significant drug pathway (GLP-1 agonists), not just the
    baseline glucose/insulin model. Compares semaglutide vs. placebo,
    SAME meal size, against Hjerpsted et al. 2018 (PMID:28941314).

    The model used to OVERSHOOT the real effect substantially -- but that
    overshoot was itself found (CLAUDE.md "DONE -- Fixed gastric_emptying_
    modifier mass-non-conservation") to be driven almost entirely by a real
    bug: gastric_emptying_modifier was scaling meal Ra(t)'s OUTPUT rather
    than the underlying absorption rate, which deletes carbohydrate mass
    instead of merely delaying it. With that fixed (mass genuinely
    conserved, confirmed numerically), the model now UNDERSHOOTS the real
    effect instead -- documented and regression-tested here exactly like
    the overshoot was, not silently accepted or forced to pass by
    retuning GLP1 constants against this one benchmark.
    """

    @pytest.fixture(scope="class")
    def semaglutide_report(self):
        return ModelValidator().validate_semaglutide()

    def test_returns_report(self, semaglutide_report):
        assert isinstance(semaglutide_report, ValidationReport)

    def test_has_one_result(self, semaglutide_report):
        assert len(semaglutide_report.results) == 1
        assert semaglutide_report.results[0].name == "glucose_iauc_5h_reduction_mmol_h_per_L"

    def test_semaglutide_genuinely_reduces_glucose_iauc(self, semaglutide_report):
        # Regardless of the magnitude debate, the drug must lower iAUC
        # relative to placebo, not raise it or leave it unchanged --
        # this is the directional sanity check underneath the specific
        # deviation-magnitude question.
        result = semaglutide_report.results[0]
        assert result.model_value > 0, (
            "semaglutide should reduce postprandial glucose iAUC vs. placebo "
            "(a positive reduction), not increase or leave it unchanged"
        )

    def test_verdict_is_pass_and_point_within_1_percent(self, semaglutide_report):
        # B7 (2026-11-16): removing the double-counted resting glucose
        # clearance moved the model to 1.349 vs the cited 1.34 mmol.h/L
        # (+0.7%). A side effect of a fix made for the IVGTT, not a
        # calibration to this benchmark. One point, n=30, on top of an
        # OGTT-fitted meal curve. Reviewed 2026-11-21 and labelled VALIDATED
        # with a population caveat (validation_status.py).
        assert semaglutide_report.verdict == ValidationVerdict.PASS
        assert semaglutide_report.results[0].deviation_percent == pytest.approx(0.7, abs=0.5)

    def test_model_no_longer_undershoots_after_b7(self, semaglutide_report):
        # History: an earlier overshoot was an artifact of deleted meal mass;
        # after the mass-conservation fix the model undershot (0.42, then
        # 0.64 and 0.79 mmol.h/L after the B3 meal-curve change). B7
        # (2026-11-16) removed a double-counted resting glucose clearance,
        # which raises every placebo excursion and so the drug's absolute
        # reduction: 1.35 vs 1.34. Pinned so further drift is looked at.
        result = semaglutide_report.results[0]
        assert result.model_value == pytest.approx(1.349, abs=0.01)


class TestValidateSteroidMealT2D:
    """
    Compares prednisone 40mg vs. no steroid, same meal, in a T2D
    population, against Shah et al. 2020 (PMID:32369480). Deliberately
    targets T2D, not the healthy default -- real evidence (two
    independent studies, PMID:27424036 and PMID:32369480) shows the
    acute glucose-raising effect of short-term prednisone is
    concentrated in T2D (reduced beta-cell reserve) and was NOT
    statistically significant in healthy subjects.
    """

    @pytest.fixture(scope="class")
    def steroid_report(self):
        return ModelValidator().validate_steroid_meal_t2d()

    def test_returns_report(self, steroid_report):
        assert isinstance(steroid_report, ValidationReport)

    def test_has_one_result(self, steroid_report):
        assert len(steroid_report.results) == 1
        assert steroid_report.results[0].name == "glucose_auc_percent_increase_t2d"

    def test_steroid_genuinely_raises_glucose_auc(self, steroid_report):
        # Directional sanity check: prednisone must raise AUC vs. no
        # steroid in this T2D population, not lower or leave it unchanged.
        result = steroid_report.results[0]
        assert result.model_value > 0, (
            "prednisone should raise postprandial glucose AUC vs. no "
            "steroid in a T2D population (a positive percent increase)"
        )

    def test_verdict_is_miss(self, steroid_report):
        assert steroid_report.verdict == ValidationVerdict.MISS

    def test_wide_tolerance_deviation_reflects_the_corrected_model(self, steroid_report):
        # Two independent, deliberate changes moved this WIDE-tier check
        # (see CLAUDE.md's glucocorticoid PK/PD milestone entry):
        #   1. Real PD timing landed on CortisolPathway: the dose is now
        #      2h before the meal (DOSE_TO_MEAL_INTERVAL_HOURS, matching
        #      Shah et al. 2020's own protocol -- "Two hours later,
        #      another OGTT" -- instead of the old, confounded
        #      simultaneous dose+meal), and the glycaemic effect ramps up
        #      from near-zero rather than applying its full magnitude
        #      instantly. This alone shrinks the modeled AUC increase.
        #   2. steroid_t2d_postprandial_glucose()'s own reference value was
        #      WRONG and has been corrected: it claimed "~50% vs. ~19%"
        #      but Shah's own Table 1 gives AUC glucose directly (T2DM
        #      616.8 -> 734.1 mmol/L/min, +19.0%) -- the ~50%/19% figures
        #      were a different metric (HOMA-%B/insulinogenic-index
        #      percentile differences) mistakenly attached to this
        #      BenchmarkPoint. Reference is now 19.0%, not 50.0%.
        # Together: model=8.4% vs. the CORRECTED 19.0% reference, a -56.0%
        # deviation.
        #   3. (2026-10-23, MASTER_PLAN row 1.46 (a)) the T2D preset gained a
        #      cited p1 factor of 0.56 (Basu 2009, PMID:19196896). Lower
        #      glucose effectiveness raises the placebo arm's own excursion,
        #      so the relative steroid increment falls: model 8.4% -> 7.5%,
        #      deviation -56.0% -> -60.5%. Not tuned against.
        #   4. (2026-11-15) default meals became two-compartment: model
        #      7.5% -> 7.6%, deviation -60.5% -> -60.0%. Not tuned against.
        #   5. (2026-11-16, B7) double-counted resting glucose clearance
        #      removed: placebo excursion grows more than the steroid one,
        #      model 7.6% -> 6.0%, deviation -60.0% -> -68.2% (WORSE). Not
        #      tuned against.
        # This is the same "a correct fix can honestly move a
        # benchmark further from a real reference, and that gets
        # documented, not hidden or curve-fit away" precedent as the
        # semaglutide-benchmark and Vg fixes elsewhere in this project.
        # The overall verdict reports this as a MISS (see
        # test_verdict_is_miss); this test pins the size of the gap.
        result = steroid_report.results[0]
        assert result.deviation_percent == pytest.approx(-68.2, abs=0.5), (
            f"glucose_auc_percent_increase_t2d deviation "
            f"{result.deviation_percent:+.1f}% -- model={result.model_value:.1f}% "
            f"vs. real reference={result.benchmark.reference:.1f}% "
            f"(if this moves further, re-check whether it's still explained "
            f"by this fix or by a new, unrelated change)"
        )


class TestValidateMetformin:
    """
    Compares chronic metformin (14-day dosing) vs. placebo, same 4h
    mixed meal test, in a T2D population, against Hansen et al. 2024
    (PMID:39049802). This is the metformin acute-meal-challenge
    benchmark CLAUDE.md's VALIDATION STATUS section flagged as the
    explicit next step. The model currently UNDERSHOOTS the real effect
    substantially -- documented and regression-tested honestly, not
    forced to pass by retuning metformin_modifiers() against this one
    benchmark (same discipline as the semaglutide/steroid benchmarks
    above).
    """

    @pytest.fixture(scope="class")
    def metformin_report(self):
        return ModelValidator().validate_metformin()

    def test_returns_report(self, metformin_report):
        assert isinstance(metformin_report, ValidationReport)

    def test_has_one_result(self, metformin_report):
        assert len(metformin_report.results) == 1
        assert metformin_report.results[0].name == "glucose_iauc_4h_reduction_mmol_min_per_L_t2d"

    def test_metformin_genuinely_reduces_glucose_iauc(self, metformin_report):
        # Directional sanity check: metformin must lower iAUC vs.
        # placebo in this T2D population, not raise or leave it unchanged.
        result = metformin_report.results[0]
        assert result.model_value > 0, (
            "metformin should reduce postprandial glucose iAUC vs. placebo "
            "in a T2D population (a positive reduction), not increase or "
            "leave it unchanged"
        )

    def test_verdict_is_pass_inside_the_wide_band(self, metformin_report):
        # Was MISS (-78%). B7 (2026-11-16) removed a double-counted resting
        # glucose clearance, which raises the T2D placebo excursion and so
        # metformin's absolute iAUC reduction: -30.6%, inside the WIDE +/-50%
        # band only. Not tuned against; the module stays EXPERIMENTAL.
        assert metformin_report.verdict == ValidationVerdict.PASS
        assert metformin_report.results[0].deviation_percent == pytest.approx(-30.6, abs=0.5)

    def test_model_currently_undershoots_the_real_effect(self, metformin_report):
        # Honest snapshot of where this benchmark stands today, not a
        # target: the model's chronic-metformin effect on 4h postprandial
        # glucose iAUC is real and correctly-directioned but far smaller
        # than Hansen et al.'s measured 186 mmol.min/L reduction (see
        # CLAUDE.md's metformin benchmark entry for the actual deviation).
        # If this ever moves to overshoot (crossing the reference value),
        # that's a substantive change worth a deliberate look, not a
        # silent pass-through -- update this test deliberately rather than
        # deleting it, same precedent as the semaglutide/steroid tests
        # above.
        result = metformin_report.results[0]
        assert 0 < result.model_value < result.benchmark.reference, (
            f"model reduction {result.model_value:.2f} mmol.min/L is no "
            f"longer strictly between 0 and the real "
            f"{result.benchmark.reference} mmol.min/L reference -- if this "
            f"is a genuine improvement toward the real value, update this "
            f"test deliberately rather than deleting it"
        )


class TestValidateSGLT2i:
    """
    Compares canagliflozin 300mg vs. placebo, same meal, healthy subjects,
    against Polidori et al. 2013 (PMID:23412078). Canagliflozin's real
    effect at this dose is confounded by intestinal SGLT1 co-inhibition
    (confirmed in the source abstract), which this platform's renal-only
    RenalPathway does not model -- the model is EXPECTED to undershoot,
    documented and regression-tested honestly, same discipline as the
    semaglutide/metformin benchmarks above.
    """

    @pytest.fixture(scope="class")
    def sglt2i_report(self):
        return ModelValidator().validate_sglt2i()

    def test_returns_report(self, sglt2i_report):
        assert isinstance(sglt2i_report, ValidationReport)

    def test_has_two_results(self, sglt2i_report):
        assert len(sglt2i_report.results) == 2
        names = {r.name for r in sglt2i_report.results}
        assert names == {"glucose_auc_2h_percent_reduction", "insulin_auc_2h_percent_reduction"}

    def test_sglt2i_genuinely_reduces_both_glucose_and_insulin_auc(self, sglt2i_report):
        # Directional sanity check: dapagliflozin must lower both AUCs vs.
        # placebo, not raise or leave them unchanged.
        for result in sglt2i_report.results:
            assert result.model_value > 0, (
                f"{result.name} should be a positive reduction vs. placebo, "
                f"got {result.model_value}"
            )

    def test_verdict_is_miss(self, sglt2i_report):
        assert sglt2i_report.verdict == ValidationVerdict.MISS

    def test_model_undershoots_the_confounded_real_effect(self, sglt2i_report):
        # Honest snapshot, not a target: canagliflozin's real 35%/43%
        # reductions are NOT purely renal (confirmed SGLT1/intestinal
        # confound -- see canagliflozin_postprandial_glucose_and_insulin()'s
        # docstring), so this platform's renal-only model is expected to
        # undershoot both. If this ever moves to overshoot (crossing the
        # reference value), that's a substantive change worth a
        # deliberate look, not a silent pass-through -- update this test
        # deliberately rather than deleting it, same precedent as the
        # semaglutide/metformin tests above.
        for result in sglt2i_report.results:
            assert 0 < result.model_value < result.benchmark.reference, (
                f"{result.name}: model value {result.model_value:.2f}% is no "
                f"longer strictly between 0 and the real "
                f"{result.benchmark.reference}% reference -- if this is a "
                f"genuine improvement toward the real value, update this "
                f"test deliberately rather than deleting it"
            )


class TestValidateProteinCoingestion:
    """
    50 g glucose alone vs. with 25 g / 50 g pea protein, healthy subjects,
    180-min iAUCs -- Thondre et al. 2021 (PMID:33515092). Independent of the
    nutrient-secretion calibration sources (Krezowski 1986, Nuttall 1984).
    The model only represents protein-stimulated insulin secretion, not
    protein's slowing of gastric emptying or its glucagon effect, so it is
    EXPECTED to undershoot -- most sharply for glucose.
    """

    @pytest.fixture(scope="class")
    def protein_report(self):
        return ModelValidator().validate_protein_coingestion()

    def test_returns_report_with_four_points(self, protein_report):
        assert isinstance(protein_report, ValidationReport)
        assert {r.name for r in protein_report.results} == {
            "glucose_iauc_3h_percent_reduction_25g_protein",
            "glucose_iauc_3h_percent_reduction_50g_protein",
            "insulin_iauc_3h_percent_increase_25g_protein",
            "insulin_iauc_3h_percent_increase_50g_protein",
        }

    def test_all_points_are_wide_and_scored(self, protein_report):
        for r in protein_report.results:
            assert r.benchmark.tolerance == Tolerance.WIDE
            assert r.benchmark.is_scored
            assert r.benchmark.pubmed_id == "33515092"

    def test_effects_point_the_right_way_and_grow_with_dose(self, protein_report):
        v = {r.name: r.model_value for r in protein_report.results}
        assert 0 < v["glucose_iauc_3h_percent_reduction_25g_protein"]             < v["glucose_iauc_3h_percent_reduction_50g_protein"]
        assert 0 < v["insulin_iauc_3h_percent_increase_25g_protein"]             < v["insulin_iauc_3h_percent_increase_50g_protein"]

    def test_model_undershoots_every_point(self, protein_report):
        # Honest snapshot, not a target (same precedent as the SGLT2i test):
        # if a value ever crosses its reference, look at it deliberately.
        for r in protein_report.results:
            assert 0 < r.model_value < r.benchmark.reference, (
                f"{r.name}: {r.model_value:.2f}% vs reference {r.benchmark.reference}%"
            )

    def test_every_arm_ran_and_was_scored(self, protein_report):
        # A solver failure leaves results empty, which _determine_verdict
        # reports as PASS -- so check that all four points were really scored.
        assert protein_report.n_pass + protein_report.n_fail == 4
        assert not any("solver failed" in n for n in protein_report.notes)


class TestValidateT2DMealInsulin:
    """
    T2D vs. matched non-diabetic subjects, same 75 g glucose drink -- Basu et
    al. 2009 (PMID:19196896). Scored as within-model T2D/healthy ratios. The
    T2D preset secretes far too little insulin relative to real early T2D
    (MASTER_PLAN row 1.46), so every point is expected to MISS; the benchmark
    exists to make that gap measurable, not to be tuned against.
    """

    @pytest.fixture(scope="class")
    def t2d_report(self):
        return ModelValidator().validate_t2d_meal_insulin()

    def test_returns_report_with_five_points(self, t2d_report):
        assert isinstance(t2d_report, ValidationReport)
        assert {r.name for r in t2d_report.results} == {
            "insulin_peak_ratio_t2d_vs_nondiabetic",
            "insulin_iauc_2h_ratio_t2d_vs_nondiabetic",
            "insulin_iauc_6h_ratio_t2d_vs_nondiabetic",
            "insulin_peak_delay_min_t2d_minus_nondiabetic",
            "glucose_iauc_6h_ratio_t2d_vs_nondiabetic",
        }

    def test_all_points_are_wide_scored_and_cite_basu(self, t2d_report):
        for r in t2d_report.results:
            assert r.benchmark.tolerance == Tolerance.WIDE
            assert r.benchmark.is_scored
            assert r.benchmark.pubmed_id == "19196896"
            assert r.benchmark.population == "type_2_diabetes"

    def test_reference_values_are_the_quoted_ratios(self, t2d_report):
        # Basu 2009: peak insulin 382 vs 673 pmol/l at 156 vs 62 min;
        # insulin iAUC 21 vs 52 nmol/l (2 h) and 79 vs 80 nmol/l (6 h);
        # glucose iAUC 1,671 vs 557 over 6 h.
        ref = {r.name: r.benchmark.reference for r in t2d_report.results}
        assert ref["insulin_peak_ratio_t2d_vs_nondiabetic"] == pytest.approx(382 / 673)
        assert ref["insulin_iauc_2h_ratio_t2d_vs_nondiabetic"] == pytest.approx(21 / 52)
        assert ref["insulin_iauc_6h_ratio_t2d_vs_nondiabetic"] == pytest.approx(79 / 80)
        assert ref["insulin_peak_delay_min_t2d_minus_nondiabetic"] == pytest.approx(156 - 62)
        assert ref["glucose_iauc_6h_ratio_t2d_vs_nondiabetic"] == pytest.approx(1671 / 557)

    def test_both_arms_ran_and_every_point_was_scored(self, t2d_report):
        # A solver failure leaves results empty, which _determine_verdict
        # reports as PASS -- so check that all five points were really scored.
        assert t2d_report.n_pass + t2d_report.n_fail == 5
        assert not any("solver failed" in n for n in t2d_report.notes)

    def test_model_undershoots_every_point(self, t2d_report):
        # Honest snapshot of the measured gap (row 1.46), not a target: the
        # preset gives T2D less insulin than the source, with too little
        # delay, and too mild a glucose excursion (the glucose ratio is
        # covered too). If any value crosses its reference, look at it
        # deliberately rather than loosening this test.
        v = {r.name: r for r in t2d_report.results}
        for name in (
            "insulin_peak_ratio_t2d_vs_nondiabetic",
            "insulin_iauc_2h_ratio_t2d_vs_nondiabetic",
            "insulin_iauc_6h_ratio_t2d_vs_nondiabetic",
            "insulin_peak_delay_min_t2d_minus_nondiabetic",
            "glucose_iauc_6h_ratio_t2d_vs_nondiabetic",
        ):
            assert 0 < v[name].model_value < v[name].benchmark.reference, name

    def test_misses_are_surfaced_in_notes_not_hidden(self, t2d_report):
        # 4/5 miss since the T2D p1 factor (row 1.46 (a), 2026-10-23) moved
        # the glucose iAUC ratio 1.11 -> 1.51, just inside its 1.5-4.5 band;
        # every insulin point still misses (the structural gap).
        assert t2d_report.verdict == ValidationVerdict.MISS
        assert t2d_report.n_fail == 4
        passed = {r.name for r in t2d_report.results if r.passed}
        assert passed == {"glucose_iauc_6h_ratio_t2d_vs_nondiabetic"}
        assert any("4 cited benchmark point(s) missed" in n for n in t2d_report.notes)
        assert any("1.46" in n for n in t2d_report.notes)


class TestValidateCaffeineOgtt:
    """
    Caffeine 5 mg/kg one hour before a 75 g OGTT -- Graham et al. 2001
    (PMID:11478588, healthy young men: insulin AUC +60%) and Robinson et al.
    2004 (PMID:15465742, T2D men: insulin AUC +25%). Independent of the
    caffeine calibration source (Keijzers 2002, PMID:11815511). The model
    currently gives the WRONG SIGN (an uncited +2% Gb term raises the delta
    secretion baseline), so both points are expected to miss; the benchmark
    makes that measurable and is not a target to tune the Gb term against.
    """

    @pytest.fixture(scope="class")
    def caffeine_report(self):
        return ModelValidator().validate_caffeine_ogtt()

    def test_returns_report_with_two_points(self, caffeine_report):
        assert isinstance(caffeine_report, ValidationReport)
        assert {r.name for r in caffeine_report.results} == {
            "insulin_auc_percent_increase_healthy_caffeine",
            "insulin_auc_percent_increase_t2d_caffeine",
        }

    def test_points_cite_the_right_papers_wide_and_scored(self, caffeine_report):
        by_name = {r.name: r for r in caffeine_report.results}
        healthy = by_name["insulin_auc_percent_increase_healthy_caffeine"].benchmark
        t2d = by_name["insulin_auc_percent_increase_t2d_caffeine"].benchmark
        assert (healthy.pubmed_id, healthy.reference) == ("11478588", 60.0)
        assert (t2d.pubmed_id, t2d.reference) == ("15465742", 25.0)
        for r in caffeine_report.results:
            assert r.benchmark.tolerance == Tolerance.WIDE
            assert r.benchmark.is_scored

    def test_all_arms_ran_and_every_point_was_scored(self, caffeine_report):
        assert caffeine_report.n_pass + caffeine_report.n_fail == 2
        assert not any("solver failed" in n for n in caffeine_report.notes)

    def test_model_misses_healthy_with_the_right_sign_but_far_too_small(self, caffeine_report):
        # Was the wrong sign (~-5%) until 2026-11-20, when an uncited +2% Gb
        # term that lowered secretion was removed (drug stream). Now +1.3% vs
        # the real +60%: right direction, still a large MISS -- the
        # insulin-raising mechanism is not modelled.
        healthy = next(r for r in caffeine_report.results if "healthy" in r.name)
        assert 0 < healthy.model_value < healthy.benchmark.reference
        assert healthy.model_value == pytest.approx(1.34, abs=0.2)
        assert not healthy.passed

    def test_model_t2d_shows_essentially_no_caffeine_effect(self, caffeine_report):
        # The T2D preset is nearly insulin-independent (MASTER_PLAN row 1.46),
        # so the model's T2D effect is ~-1%, not a robust "wrong sign": assert
        # only that it is far from the observed +25%.
        t2d = next(r for r in caffeine_report.results if "t2d" in r.name)
        assert abs(t2d.model_value) < 5.0
        assert not t2d.passed

    def test_misses_are_surfaced_in_notes(self, caffeine_report):
        assert caffeine_report.verdict == ValidationVerdict.MISS
        assert any("2 cited benchmark point(s) missed" in n for n in caffeine_report.notes)
        assert any("not modelled" in n for n in caffeine_report.notes)


class TestValidateAll:
    """
    validate_all() -- convenience entry point added once four separate
    validate_*() methods existed with no single place to run them all.
    """

    @pytest.fixture(scope="class")
    def all_reports(self):
        return ModelValidator().validate_all()

    def test_returns_dict_of_reports(self, all_reports):
        assert isinstance(all_reports, dict)
        assert all(isinstance(r, ValidationReport) for r in all_reports.values())

    def test_includes_all_nine_registered_benchmarks(self, all_reports):
        # If another benchmark is ever added (validate_*), update this
        # test deliberately -- it exists so validate_all() can't silently
        # drift out of sync with the individual methods.
        assert set(all_reports.keys()) == {
            "OGTT_75g_healthy",
            "Bergman_IV_glucose",
            "Semaglutide_1mg_postprandial",
            "Prednisone_40mg_T2D_postprandial",
            "Metformin_T2D_postprandial",
            "Canagliflozin_300mg_postprandial",
            "Pea_protein_coingestion_50g_glucose",
            "T2D_vs_nondiabetic_meal_insulin",
            "Caffeine_5mgkg_OGTT_insulin",
        }

    def test_no_dataset_fails_outright(self, all_reports):
        # A hard FAIL here (as opposed to a WIDE-tolerance miss within an
        # otherwise-PASS report) would mean a TIGHT/NORMAL benchmark
        # broke -- that should never happen silently.
        # KNOWN FAIL since 2026-11-16: the in-house OGTT (see
        # test_validate_ogtt_convenience). Every other dataset must not FAIL.
        for name, report in all_reports.items():
            if name == "OGTT_75g_healthy":
                continue
            assert report.verdict != ValidationVerdict.FAIL, (
                f"{name} hard-failed:\n{report.as_text()}"
            )

    def test_each_dataset_reports_its_real_verdict(self, all_reports):
        # Since B7 (2026-11-16, double-counted resting glucose clearance
        # removed, nothing re-tuned): the in-house OGTT drops PASS -> WARN
        # (peak glucose +33%, outside its NORMAL band -- its absorption and
        # first-phase constants were fitted with the double count in place);
        # semaglutide (+0.7%) and metformin (-31%, WIDE band) PASS; every
        # other drug/disease benchmark MISSes. If one of these changes, look
        # at it deliberately -- don't loosen the expectation.
        verdicts = {name: r.verdict for name, r in all_reports.items()}
        miss = ValidationVerdict.MISS
        assert verdicts == {
            "OGTT_75g_healthy": ValidationVerdict.FAIL,  # see test_validate_ogtt_convenience
            "Bergman_IV_glucose": ValidationVerdict.PASS,
            "Semaglutide_1mg_postprandial": ValidationVerdict.PASS,
            "Prednisone_40mg_T2D_postprandial": miss,
            "Metformin_T2D_postprandial": ValidationVerdict.PASS,
            "Canagliflozin_300mg_postprandial": miss,
            "Pea_protein_coingestion_50g_glucose": miss,
            "T2D_vs_nondiabetic_meal_insulin": miss,
            "Caffeine_5mgkg_OGTT_insulin": miss,
        }

    def test_verdicts_stay_independent_per_dataset(self, all_reports):
        # Regression: validate_all() must not compute any combined/shared
        # verdict across datasets -- each report's verdict must match
        # calling its own validate_*() method directly.
        validator = ModelValidator()
        assert all_reports["OGTT_75g_healthy"].verdict == validator.validate_ogtt().verdict
        assert all_reports["Bergman_IV_glucose"].verdict == validator.validate_bergman_iv().verdict


# ---------------------------------------------------------------------------
# Verdict logic
# ---------------------------------------------------------------------------

class TestVerdictLogic:
    def _make_report(self, results, dataset_name="test"):
        n_pass = sum(1 for r in results if r.passed)
        n_fail = len(results) - n_pass
        validator = ModelValidator()
        verdict = validator._determine_verdict(results)
        return ValidationReport(
            dataset_name=dataset_name, scenario_name="test",
            results=results, verdict=verdict,
            n_pass=n_pass, n_fail=n_fail,
        )

    def _pass_result(self, tolerance=Tolerance.NORMAL) -> BenchmarkCheckResult:
        bp = BenchmarkPoint("p","d",5.0,"mM",tolerance,"0")
        return bp.check(5.0)

    def _fail_result(self, tolerance=Tolerance.NORMAL) -> BenchmarkCheckResult:
        bp = BenchmarkPoint("p","d",5.0,"mM",tolerance,"0")
        return bp.check(99.0)

    def test_all_pass_gives_pass(self):
        results = [self._pass_result(), self._pass_result()]
        report = self._make_report(results)
        assert report.verdict == ValidationVerdict.PASS

    def test_tight_fail_gives_fail(self):
        results = [self._fail_result(Tolerance.TIGHT)]
        report = self._make_report(results)
        assert report.verdict == ValidationVerdict.FAIL

    def test_normal_fail_minority_gives_warn(self):
        results = [
            self._pass_result(Tolerance.NORMAL),
            self._pass_result(Tolerance.NORMAL),
            self._fail_result(Tolerance.NORMAL),
        ]
        report = self._make_report(results)
        assert report.verdict == ValidationVerdict.WARN

    def test_normal_fail_majority_gives_fail(self):
        results = [
            self._fail_result(Tolerance.NORMAL),
            self._fail_result(Tolerance.NORMAL),
            self._pass_result(Tolerance.NORMAL),
        ]
        report = self._make_report(results)
        assert report.verdict == ValidationVerdict.FAIL

    def test_unscored_normal_point_cannot_dilute_a_failure(self):
        # B4 review (2026-11-16): an UNVERIFIED (unscored) NORMAL point that
        # passes must not count toward the NORMAL-failure ratio. One scored
        # NORMAL failure out of one scored NORMAL point is a FAIL.
        scored_fail = self._fail_result(Tolerance.NORMAL)
        unscored_pass = BenchmarkPoint("u", "d", 5.0, "mM", Tolerance.NORMAL, "UNVERIFIED").check(5.0)
        report = self._make_report([scored_fail, unscored_pass])
        assert report.verdict == ValidationVerdict.FAIL

    def test_unscored_tight_failure_cannot_cause_fail(self):
        unscored_fail = BenchmarkPoint("u", "d", 5.0, "mM", Tolerance.TIGHT, "UNVERIFIED").check(99.0)
        report = self._make_report([self._pass_result(), unscored_fail])
        assert report.verdict == ValidationVerdict.PASS

    def test_scored_wide_miss_gives_miss_not_pass(self):
        # A cited point outside even the +/-50% band is a real miss and must
        # say so -- the old policy reported PASS here.
        results = [
            self._pass_result(Tolerance.NORMAL),
            self._fail_result(Tolerance.WIDE),
        ]
        report = self._make_report(results)
        assert report.verdict == ValidationVerdict.MISS

    def test_no_results_gives_no_data_not_pass(self):
        report = self._make_report([])
        assert report.verdict == ValidationVerdict.NO_DATA
        assert report.verdict.value == "NO DATA"

    def test_only_unscored_points_gives_no_data_not_pass(self):
        # A point with no verified citation is informational only, so a
        # report made only of them has nothing to validate against.
        unscored = BenchmarkPoint("p", "d", 5.0, "mM", Tolerance.NORMAL, "UNVERIFIED")
        report = self._make_report([unscored.check(5.0)])
        assert report.verdict == ValidationVerdict.NO_DATA

    def test_unscored_wide_miss_does_not_give_miss(self):
        # A point with no verified citation is informational only.
        unscored = BenchmarkPoint("p", "d", 5.0, "mM", Tolerance.WIDE, "UNVERIFIED")
        results = [self._pass_result(Tolerance.NORMAL), unscored.check(99.0)]
        report = self._make_report(results)
        assert report.verdict == ValidationVerdict.PASS

    def test_normal_failure_still_outranks_miss(self):
        results = [
            self._pass_result(Tolerance.NORMAL),
            self._pass_result(Tolerance.NORMAL),
            self._fail_result(Tolerance.NORMAL),
            self._fail_result(Tolerance.WIDE),
        ]
        report = self._make_report(results)
        assert report.verdict == ValidationVerdict.WARN

    def test_summary_shows_each_points_deviation_and_status(self):
        results = [self._pass_result(Tolerance.NORMAL), self._fail_result(Tolerance.WIDE)]
        summary = self._make_report(results).summary()
        assert summary["points"] == [
            {"name": "p", "tolerance": "normal", "deviation_percent": 0.0, "status": "PASS"},
            {"name": "p", "tolerance": "wide", "deviation_percent": 1880.0, "status": "MISS"},
        ]


class TestDrugBenchmarksLabelled:
    """
    Drug benchmarks carry an "experimental module" (or, for semaglutide, a
    "validated module" + caveat) label so a reader sees the
    misses in context. The label never changes a verdict (pinned above in
    TestValidateAll.test_each_dataset_reports_its_real_verdict).
    """

    DRUG_DATASETS = {
        "Semaglutide_1mg_postprandial": "semaglutide",
        "Prednisone_40mg_T2D_postprandial": "glucocorticoid",
        "Metformin_T2D_postprandial": "metformin",
        "Canagliflozin_300mg_postprandial": "dapagliflozin",
        "Caffeine_5mgkg_OGTT_insulin": "caffeine",
    }

    @pytest.fixture(scope="class")
    def all_reports(self):
        return ModelValidator().validate_all()

    def test_drug_reports_are_labelled(self, all_reports):
        for name, drug in self.DRUG_DATASETS.items():
            label = all_reports[name].module_label
            expected = (
                "validated module" if drug == "semaglutide"
                else "experimental module"
            )
            assert expected in label, name
            assert drug in label, name
            assert expected in all_reports[name].as_text()

    def test_validated_drug_label_carries_its_caveat(self, all_reports):
        # User decision 2026-11-21: semaglutide VALIDATED, caveat kept.
        from biosim_physiology.drug_interactions import Drug
        from biosim_physiology.drug_interactions.validation_status import DRUG_VALIDATION
        label = all_reports["Semaglutide_1mg_postprandial"].module_label
        assert "status VALIDATED" in label
        assert DRUG_VALIDATION[Drug.SEMAGLUTIDE].caveat in label

    def test_non_drug_reports_are_not_labelled(self, all_reports):
        for name, report in all_reports.items():
            if name not in self.DRUG_DATASETS:
                assert report.module_label == "", name

    def test_status_table_deviations_match_the_live_validator(self, all_reports):
        # The per-drug table quotes each benchmark's % error; if the model
        # moves, this fails so the table (and the docs) get updated.
        from biosim_physiology.drug_interactions import Drug
        from biosim_physiology.drug_interactions.validation_status import DRUG_VALIDATION
        for name, drug in self.DRUG_DATASETS.items():
            record = DRUG_VALIDATION[Drug(drug)]
            assert record.benchmark == name
            live = [r.deviation_percent for r in all_reports[name].results]
            assert record.deviations_percent == pytest.approx(live, abs=0.5), name


def test_no_benchmark_point_cites_the_clamp_methods_paper():
    # B4 (2026-11-16): DeFronzo 1979 (PMID:382871) is the glucose-clamp
    # methods paper; it reports no OGTT 2-hour glucose value.
    from biosim_validation.benchmark import ogtt_75g_healthy_adults
    assert all(p.pubmed_id != "382871" for p in ogtt_75g_healthy_adults().points)
