# SAVE TO: packages\uncertainty\tests\test_calibration.py
"""
Test suite for BayesianCalibrator (Sampling-Importance-Resampling
parameter calibration against clinical benchmark reference values).
"""
import math

import pytest
import numpy as np

from biosim_ode import high_sugar_meal
from biosim_uncertainty import (
    BayesianCalibrator, CalibrationResult, CALIBRATION_METRIC_EXTRACTORS,
)
from biosim_uncertainty.distributions import LogNormalDistribution, FixedValue


class FakeTarget:
    """Minimal object satisfying the CalibrationTarget protocol for tests."""
    def __init__(self, name, reference, tolerance_fraction, is_scored=True):
        self.name = name
        self.reference = reference
        self.tolerance_fraction = tolerance_fraction
        self.is_scored = is_scored


# n_prior_samples kept small throughout for test speed, matching this
# package's existing TestMonteCarlo convention (n_samples=20-30).

class TestCalibrationTargetProtocol:
    def test_real_benchmark_point_satisfies_protocol(self):
        # BenchmarkPoint must duck-type against CalibrationTarget without
        # biosim_uncertainty importing biosim_validation at all (checked
        # here, in a test, which may cross-import for verification).
        from biosim_uncertainty.calibration import CalibrationTarget
        from biosim_validation.benchmark import ogtt_75g_healthy_adults
        point = ogtt_75g_healthy_adults().points[0]
        assert isinstance(point, CalibrationTarget)


class TestBayesianCalibratorValidation:
    def test_non_positive_likelihood_sigma_z_raises(self):
        with pytest.raises(ValueError, match="likelihood_sigma_z"):
            BayesianCalibrator(likelihood_sigma_z=0.0)

    def test_invalid_min_ess_fraction_raises(self):
        with pytest.raises(ValueError, match="min_effective_sample_size_fraction"):
            BayesianCalibrator(min_effective_sample_size_fraction=1.5)
        with pytest.raises(ValueError, match="min_effective_sample_size_fraction"):
            BayesianCalibrator(min_effective_sample_size_fraction=0.0)

    def test_no_scored_targets_raises(self):
        calibrator = BayesianCalibrator(random_seed=42)
        unscored = FakeTarget("mtor_activation", 0.55, 0.50, is_scored=False)
        with pytest.raises(ValueError, match="no scored calibration_targets"):
            calibrator.calibrate(
                distributions={"p1": LogNormalDistribution.from_mean_cv(0.0317, 0.30)},
                calibration_targets=[unscored],
                metric_extractors=CALIBRATION_METRIC_EXTRACTORS,
                meals=[high_sugar_meal()], duration_s=3600, n_prior_samples=20,
            )

    def test_target_reference_zero_raises(self):
        calibrator = BayesianCalibrator(random_seed=42)
        target = FakeTarget("peak_glucose", 0.0, 0.25)
        with pytest.raises(ValueError, match="non-positive"):
            calibrator.calibrate(
                distributions={"p1": LogNormalDistribution.from_mean_cv(0.0317, 0.30)},
                calibration_targets=[target],
                metric_extractors=CALIBRATION_METRIC_EXTRACTORS,
                meals=[high_sugar_meal()], duration_s=3600, n_prior_samples=20,
            )

    def test_all_prior_samples_failing_raises(self, monkeypatch):
        import biosim_ode.system as system_module

        def _raise(self, duration_s, meals=None, n_output_points=500, y0_override=None):
            raise ValueError("synthetic ODE solver failure for testing")

        monkeypatch.setattr(system_module.ODESystem, "simulate", _raise)
        calibrator = BayesianCalibrator(random_seed=42)
        target = FakeTarget("peak_glucose", 7.8, 0.25)
        with pytest.raises(ValueError, match="every prior sample failed"):
            calibrator.calibrate(
                distributions={"p1": LogNormalDistribution.from_mean_cv(0.0317, 0.30)},
                calibration_targets=[target],
                metric_extractors=CALIBRATION_METRIC_EXTRACTORS,
                meals=[high_sugar_meal()], duration_s=3600, n_prior_samples=10,
            )


class TestBayesianCalibratorCalibrate:
    @pytest.fixture(scope="class")
    def result(self) -> CalibrationResult:
        calibrator = BayesianCalibrator(random_seed=42)
        distributions = {
            "p1": LogNormalDistribution.from_mean_cv(0.0317, cv=0.30),
            "p3": LogNormalDistribution.from_mean_cv(9.032e-7, cv=0.30),
        }
        targets = [
            FakeTarget("fasting_glucose", 5.0, 0.10),
            FakeTarget("peak_glucose", 7.8, 0.25),
            FakeTarget("glucose_at_120min", 6.5, 0.25),
            FakeTarget("peak_insulin", 350.0, 0.50),
            FakeTarget("mtor_activation", 0.55, 0.50, is_scored=False),  # unscored, must be ignored
        ]
        return calibrator.calibrate(
            distributions=distributions,
            calibration_targets=targets,
            metric_extractors=CALIBRATION_METRIC_EXTRACTORS,
            meals=[high_sugar_meal()], duration_s=14400, n_prior_samples=40,
        )

    def test_returns_calibration_result(self, result):
        assert isinstance(result, CalibrationResult)

    def test_posterior_bands_cover_all_species(self, result):
        for name in ("G", "I", "X", "IGF1", "mTOR", "AMPK"):
            assert name in result.posterior_bands

    def test_posterior_parameter_samples_present_for_each_distribution(self, result):
        assert set(result.posterior_parameter_samples.keys()) == {"p1", "p3"}
        for values in result.posterior_parameter_samples.values():
            assert len(values) == result.n_posterior_samples
            assert np.all(values > 0)

    def test_unscored_target_excluded_from_targets_used(self, result):
        assert "mtor_activation" not in result.targets_used
        assert "mtor_activation" not in result.likelihood_sigma_used

    def test_effective_sample_size_bounded(self, result):
        assert 0.0 < result.effective_sample_size <= result.n_prior_samples
        assert 0.0 < result.effective_sample_size_fraction <= 1.0

    def test_likelihood_sigma_recorded_for_each_used_target(self, result):
        assert set(result.likelihood_sigma_used.keys()) == {
            "fasting_glucose", "peak_glucose", "glucose_at_120min", "peak_insulin",
        }
        # sigma = reference * tolerance_fraction / z (z=1.0 default)
        assert result.likelihood_sigma_used["peak_glucose"] == pytest.approx(7.8 * 0.25)

    def test_notes_document_the_likelihood_modeling_choice(self, result):
        assert any("modeling choice" in n for n in result.notes)

    def test_not_degenerate_with_default_settings(self, result):
        # A reasonably-matched prior/likelihood combination should not
        # collapse to a handful of samples.
        assert result.is_degenerate is False

    def test_to_uq_result_wraps_posterior_correctly(self, result):
        uq = result.to_uq_result()
        assert uq.n_samples == result.n_posterior_samples
        assert uq.n_failed == 0
        assert uq.bands is result.posterior_bands
        assert uq.sampled_parameter_values is result.posterior_parameter_samples

    def test_same_seed_reproducible(self):
        def run():
            calibrator = BayesianCalibrator(random_seed=7)
            return calibrator.calibrate(
                distributions={"p1": LogNormalDistribution.from_mean_cv(0.0317, 0.30)},
                calibration_targets=[FakeTarget("peak_glucose", 7.8, 0.25)],
                metric_extractors=CALIBRATION_METRIC_EXTRACTORS,
                meals=[high_sugar_meal()], duration_s=3600, n_prior_samples=20,
            )
        r1, r2 = run(), run()
        assert np.array_equal(
            r1.posterior_parameter_samples["p1"], r2.posterior_parameter_samples["p1"]
        )

    def test_posterior_shifts_toward_likelihood_when_prior_is_off_center(self):
        # Regression-style sanity check: a genuinely informative
        # likelihood should move the posterior parameter mean away from
        # the (deliberately mis-centered) prior mean, not leave it
        # unchanged -- confirms the reweighting actually does something,
        # not just that the machinery runs without crashing.
        calibrator = BayesianCalibrator(random_seed=42)
        low_p1_prior = LogNormalDistribution.from_mean_cv(0.015, cv=0.50)  # deliberately low
        result = calibrator.calibrate(
            distributions={"p1": low_p1_prior},
            calibration_targets=[FakeTarget("peak_glucose", 6.5, 0.10)],  # favors higher p1 (lower peak)
            metric_extractors=CALIBRATION_METRIC_EXTRACTORS,
            meals=[high_sugar_meal()], duration_s=3600, n_prior_samples=200,
        )
        prior_mean = low_p1_prior.mean()
        posterior_mean = float(np.mean(result.posterior_parameter_samples["p1"]))
        assert posterior_mean > prior_mean

    def test_target_without_matching_extractor_is_skipped_not_error(self):
        calibrator = BayesianCalibrator(random_seed=42)
        targets = [
            FakeTarget("peak_glucose", 7.8, 0.25),
            FakeTarget("no_such_metric", 1.0, 0.25),  # no entry in metric_extractors
        ]
        result = calibrator.calibrate(
            distributions={"p1": LogNormalDistribution.from_mean_cv(0.0317, 0.30)},
            calibration_targets=targets,
            metric_extractors=CALIBRATION_METRIC_EXTRACTORS,
            meals=[high_sugar_meal()], duration_s=3600, n_prior_samples=20,
        )
        assert result.targets_used == ["peak_glucose"]

    def test_metric_extractor_exception_does_not_crash_calibration(self):
        calibrator = BayesianCalibrator(random_seed=42)

        def _flaky_extractor(times_s, species):
            raise RuntimeError("synthetic extractor failure")

        result = calibrator.calibrate(
            distributions={"p1": LogNormalDistribution.from_mean_cv(0.0317, 0.30)},
            calibration_targets=[FakeTarget("flaky", 7.8, 0.25)],
            metric_extractors={"flaky": _flaky_extractor},
            meals=[high_sugar_meal()], duration_s=3600, n_prior_samples=15,
        )
        assert result.n_posterior_samples == 15
        assert any("RuntimeError" in n for n in result.notes)

    def test_nan_model_value_does_not_poison_every_samples_weight(self):
        # Regression: found during self-verification while the mandatory
        # review-agent pass was still running (see CLAUDE.md). np.max/
        # np.argmax/np.interp propagate NaN silently rather than raising
        # -- a single sample producing a NaN model value (e.g. from a
        # numerically unstable but "successful" ODE solve) used to slip
        # past the try/except entirely (no exception was raised), giving
        # that ONE sample a NaN log-likelihood. logsumexp(log_weights) is
        # NaN if even one entry is NaN, and every OTHER sample's weight
        # then also becomes NaN after normalization -- one bad sample
        # silently broke the entire calibration's weight vector, not
        # just its own contribution.
        calibrator = BayesianCalibrator(random_seed=42)

        call_count = {"n": 0}

        def _extractor_with_one_nan(times_s, species):
            call_count["n"] += 1
            if call_count["n"] == 3:   # poison exactly one sample
                return float("nan")
            return float(np.max(species["G"])) / 1000.0

        result = calibrator.calibrate(
            distributions={"p1": LogNormalDistribution.from_mean_cv(0.0317, cv=0.30)},
            calibration_targets=[FakeTarget("peak_glucose", 7.8, 0.25)],
            metric_extractors={"peak_glucose": _extractor_with_one_nan},
            meals=[high_sugar_meal()], duration_s=3600, n_prior_samples=20,
        )
        weights_are_finite = np.all(
            np.isfinite(result.posterior_parameter_samples["p1"])
        )
        assert weights_are_finite
        assert result.effective_sample_size > 0
        assert not math.isnan(result.effective_sample_size)
        assert any("non-finite model value" in n or "ValueError" in n for n in result.notes)

    def test_degenerate_flag_triggers_with_very_narrow_likelihood(self):
        # A likelihood far narrower than the achievable spread of model
        # values across the prior collapses effective sample size --
        # confirmed this z reliably triggers it before writing the test.
        calibrator = BayesianCalibrator(random_seed=42, likelihood_sigma_z=50.0)
        result = calibrator.calibrate(
            distributions={"p1": LogNormalDistribution.from_mean_cv(0.0317, cv=0.30)},
            calibration_targets=[FakeTarget("peak_glucose", 7.8, 0.25)],
            metric_extractors=CALIBRATION_METRIC_EXTRACTORS,
            meals=[high_sugar_meal()], duration_s=3600, n_prior_samples=40,
        )
        assert result.is_degenerate is True
        assert any("provisional" in n for n in result.notes)

    def test_fixed_value_distribution_gives_maximal_effective_sample_size(self):
        # With no parameter uncertainty at all, every sample is identical
        # -- the likelihood can't discriminate between samples, so ESS
        # should be at (or extremely close to) its maximum, n_successful.
        calibrator = BayesianCalibrator(random_seed=42)
        result = calibrator.calibrate(
            distributions={"p1": FixedValue(0.0317)},
            calibration_targets=[FakeTarget("peak_glucose", 7.8, 0.25)],
            metric_extractors=CALIBRATION_METRIC_EXTRACTORS,
            meals=[high_sugar_meal()], duration_s=3600, n_prior_samples=20,
        )
        assert result.effective_sample_size_fraction == pytest.approx(1.0, abs=1e-9)


class TestCalibrationMetricExtractors:
    """Sanity-check the standard extractors directly against synthetic trajectories."""

    def test_fasting_glucose_reads_initial_value(self):
        from biosim_uncertainty.calibration import _fasting_glucose
        species = {"G": np.array([5000.0, 6000.0, 5500.0])}
        assert _fasting_glucose(np.array([0.0, 60.0, 120.0]), species) == pytest.approx(5.0)

    def test_peak_glucose_reads_max(self):
        from biosim_uncertainty.calibration import _peak_glucose
        species = {"G": np.array([5000.0, 9000.0, 5500.0])}
        assert _peak_glucose(np.array([0.0, 60.0, 120.0]), species) == pytest.approx(9.0)

    def test_time_to_peak_glucose_reads_argmax_time(self):
        from biosim_uncertainty.calibration import _time_to_peak_glucose_min
        species = {"G": np.array([5000.0, 9000.0, 5500.0])}
        result = _time_to_peak_glucose_min(np.array([0.0, 1800.0, 3600.0]), species)
        assert result == pytest.approx(30.0)

    def test_glucose_at_120min_interpolates(self):
        from biosim_uncertainty.calibration import _glucose_at_120min
        species = {"G": np.array([5000.0, 7000.0, 6000.0])}
        result = _glucose_at_120min(np.array([0.0, 7200.0, 14400.0]), species)
        assert result == pytest.approx(7.0)

    def test_peak_insulin_converts_units(self):
        from biosim_uncertainty.calibration import _peak_insulin
        species = {"I": np.array([0.000066, 0.00035, 0.0001])}
        result = _peak_insulin(np.array([0.0, 60.0, 120.0]), species)
        assert result == pytest.approx(350.0)

    def test_glucose_return_time_finds_return_to_baseline(self):
        from biosim_uncertainty.calibration import _glucose_return_time_min
        # Baseline 5000, threshold 5500 (10% above); returns at index 3.
        G = np.array([5000.0, 9000.0, 7000.0, 5400.0, 5100.0])
        times_s = np.array([0.0, 600.0, 1800.0, 3600.0, 5400.0])
        result = _glucose_return_time_min(times_s, {"G": G})
        assert result == pytest.approx(60.0)

    def test_glucose_return_time_falls_back_to_final_time_if_never_returns(self):
        from biosim_uncertainty.calibration import _glucose_return_time_min
        G = np.array([5000.0, 9000.0, 8500.0, 8000.0])
        times_s = np.array([0.0, 600.0, 1800.0, 3600.0])
        result = _glucose_return_time_min(times_s, {"G": G})
        assert result == pytest.approx(60.0)   # final time (3600s) / 60
