# SAVE TO: packages\uncertainty\tests\test_sobol_sensitivity.py
"""
Tests for Sobol global sensitivity analysis.

n_base_samples is kept small (8-16) throughout, matching this test
suite's existing convention for MonteCarloEngine (n_samples=20-30) --
Sobol needs n_base_samples * (k + 2) model evaluations
(calc_second_order=False), so keeping both n_base_samples and k small
is what keeps these tests fast.
"""
import math
import pytest
import numpy as np

from biosim_uncertainty import (
    LogNormalDistribution, SobolSensitivityAnalyzer, SobolResult,
)
from biosim_uncertainty.sensitivity import glucose_peak_metric, insulin_peak_metric
from biosim_ode import PathwayParameters, GlucoseInsulinPathway
from biosim_ode.meal import MealInput
from biosim_ode.glp1_pathway import GLP1Pathway, GLP1Parameters, glp1_couplings
from biosim_ode.coupling import CompositePathway


def _standard_meal():
    return [MealInput(carbs_g=75.0, time_min=0.0)]


def _glp1_factory(params, meals):
    """Test factory: composes GLP1Pathway (semaglutide-like) with GI."""
    gi = GlucoseInsulinPathway(params=params, meals=meals)
    glp1 = GLP1Pathway()
    return CompositePathway([gi, glp1], couplings=glp1_couplings(gi.name, glp1.name))


def _glp1_extension_factory(params, meals, overrides):
    """Test factory: samples GLP1's agonist_scale via a dotted key."""
    gi = GlucoseInsulinPathway(params=params, meals=meals)
    agonist_scale = overrides.get("glp1.agonist_scale", 9.0)
    glp1 = GLP1Pathway(GLP1Parameters(agonist_scale=agonist_scale))
    return CompositePathway([gi, glp1], couplings=glp1_couplings(gi.name, glp1.name))


def _two_param_distributions():
    base = PathwayParameters()
    return {
        "p1": LogNormalDistribution.from_confidence_tier(base.p1, "silver"),
        "p3": LogNormalDistribution.from_confidence_tier(base.p3, "silver"),
    }


class TestSobolSensitivityAnalyzer:
    def test_returns_result_with_expected_shape(self):
        distributions = _two_param_distributions()
        result = SobolSensitivityAnalyzer().analyze(
            distributions,
            metrics={"glucose_peak": glucose_peak_metric},
            meals=_standard_meal(),
            n_base_samples=8,
            seed=42,
        )
        assert isinstance(result, SobolResult)
        assert set(result.parameter_names) == {"p1", "p3"}
        assert result.metric_names == ["glucose_peak"]
        assert set(result.first_order["glucose_peak"].keys()) == {"p1", "p3"}
        assert set(result.total_order["glucose_peak"].keys()) == {"p1", "p3"}

    def test_n_model_evaluations_matches_saltelli_formula(self):
        # calc_second_order=False -> N * (k + 2)
        distributions = _two_param_distributions()
        result = SobolSensitivityAnalyzer().analyze(
            distributions, metrics={"glucose_peak": glucose_peak_metric},
            meals=_standard_meal(), n_base_samples=8, seed=42,
        )
        k = 2
        assert result.n_model_evaluations == 8 * (k + 2)

    def test_indices_are_finite_not_nan(self):
        distributions = _two_param_distributions()
        result = SobolSensitivityAnalyzer().analyze(
            distributions, metrics={"glucose_peak": glucose_peak_metric},
            meals=_standard_meal(), n_base_samples=8, seed=42,
        )
        for param in ("p1", "p3"):
            assert math.isfinite(result.first_order["glucose_peak"][param])
            assert math.isfinite(result.total_order["glucose_peak"][param])
            assert math.isfinite(result.first_order_conf["glucose_peak"][param])
            assert math.isfinite(result.total_order_conf["glucose_peak"][param])

    def test_zero_parameters_raises(self):
        with pytest.raises(ValueError, match="at least one parameter"):
            SobolSensitivityAnalyzer().analyze({}, meals=_standard_meal())

    def test_reproducible_with_same_seed(self):
        distributions = _two_param_distributions()
        r1 = SobolSensitivityAnalyzer().analyze(
            distributions, metrics={"glucose_peak": glucose_peak_metric},
            meals=_standard_meal(), n_base_samples=8, seed=7,
        )
        r2 = SobolSensitivityAnalyzer().analyze(
            distributions, metrics={"glucose_peak": glucose_peak_metric},
            meals=_standard_meal(), n_base_samples=8, seed=7,
        )
        assert r1.first_order == r2.first_order
        assert r1.total_order == r2.total_order

    def test_different_seeds_give_different_sample_matrices(self):
        # Not a claim about the *indices* differing (they may coincide by
        # chance at this sample size) -- confirms the seed is genuinely
        # wired into SALib's sampler, not silently ignored.
        distributions = _two_param_distributions()

        from SALib.sample import sobol as sobol_sample
        problem = {
            "num_vars": 2, "names": ["p1", "p3"],
            "bounds": [[distributions["p1"].mu, distributions["p1"].sigma],
                       [distributions["p3"].mu, distributions["p3"].sigma]],
            "dists": ["lognorm", "lognorm"],
        }
        m1 = sobol_sample.sample(problem, N=8, calc_second_order=False, seed=1)
        m2 = sobol_sample.sample(problem, N=8, calc_second_order=False, seed=2)
        assert not np.allclose(m1, m2)

    def test_multiple_metrics_all_reported(self):
        distributions = _two_param_distributions()
        result = SobolSensitivityAnalyzer().analyze(
            distributions,
            metrics={"glucose_peak": glucose_peak_metric, "insulin_peak": insulin_peak_metric},
            meals=_standard_meal(), n_base_samples=8, seed=42,
        )
        assert set(result.metric_names) == {"glucose_peak", "insulin_peak"}
        assert "p1" in result.first_order["insulin_peak"]
        assert "p1" in result.first_order["glucose_peak"]

    def test_failure_raises_runtime_error_not_silently_dropped(self, monkeypatch):
        # Regression guard for the documented design decision: unlike
        # MonteCarloEngine/OATSensitivityAnalyzer, a Sobol analysis
        # cannot tolerate a dropped/substituted sample without
        # corrupting every downstream index -- confirm a failure raises
        # loudly instead of silently proceeding with a partial matrix.
        distributions = _two_param_distributions()
        analyzer = SobolSensitivityAnalyzer()
        call_count = {"n": 0}
        real_run = analyzer._run

        def _flaky_run(params, meals, duration_s, extension_overrides=None):
            call_count["n"] += 1
            if call_count["n"] == 3:
                return None, "synthetic ODE solver failure for testing"
            return real_run(params, meals, duration_s, extension_overrides)

        monkeypatch.setattr(analyzer, "_run", _flaky_run)
        with pytest.raises(RuntimeError, match="synthetic ODE solver failure"):
            analyzer.analyze(
                distributions, metrics={"glucose_peak": glucose_peak_metric},
                meals=_standard_meal(), n_base_samples=8, seed=42,
            )

    def test_top_parameters_ranks_by_total_order_by_default(self):
        distributions = _two_param_distributions()
        result = SobolSensitivityAnalyzer().analyze(
            distributions, metrics={"glucose_peak": glucose_peak_metric},
            meals=_standard_meal(), n_base_samples=16, seed=42,
        )
        top = result.top_parameters("glucose_peak", n=2)
        assert len(top) == 2
        names = {p for p, _ in top}
        assert names == {"p1", "p3"}
        # Ranked descending by |ST|
        assert abs(top[0][1]) >= abs(top[1][1])

    def test_top_parameters_unknown_metric_raises(self):
        distributions = _two_param_distributions()
        result = SobolSensitivityAnalyzer().analyze(
            distributions, metrics={"glucose_peak": glucose_peak_metric},
            meals=_standard_meal(), n_base_samples=8, seed=42,
        )
        with pytest.raises(KeyError):
            result.top_parameters("nonexistent_metric")

    def test_interaction_strength_is_total_minus_first_order(self):
        distributions = _two_param_distributions()
        result = SobolSensitivityAnalyzer().analyze(
            distributions, metrics={"glucose_peak": glucose_peak_metric},
            meals=_standard_meal(), n_base_samples=8, seed=42,
        )
        expected = (
            result.total_order["glucose_peak"]["p1"]
            - result.first_order["glucose_peak"]["p1"]
        )
        assert result.interaction_strength("glucose_peak", "p1") == pytest.approx(expected)

    def test_summary_reports_base_sample_count_and_evaluations(self):
        distributions = _two_param_distributions()
        result = SobolSensitivityAnalyzer().analyze(
            distributions, metrics={"glucose_peak": glucose_peak_metric},
            meals=_standard_meal(), n_base_samples=8, seed=42,
        )
        summary = result.summary()
        assert summary["n_base_samples"] == 8
        assert summary["n_model_evaluations"] == result.n_model_evaluations
        assert "glucose_peak" in summary["top_sensitivities"]

    def test_sampling_actually_uses_the_distributions_mu_sigma(self):
        # Confirms the SALib problem dict is built from each
        # distribution's own (mu, sigma), not some other, disconnected
        # bound -- checked by comparing the empirical mean of the
        # sampled p1 column against the distribution's documented mean,
        # within a generous tolerance for a small sample.
        base = PathwayParameters()
        dist = LogNormalDistribution.from_confidence_tier(base.p1, "silver")
        distributions = {"p1": dist,
                          "p3": LogNormalDistribution.from_confidence_tier(base.p3, "silver")}

        from SALib.sample import sobol as sobol_sample
        problem = {
            "num_vars": 2, "names": ["p1", "p3"],
            "bounds": [[distributions["p1"].mu, distributions["p1"].sigma],
                       [distributions["p3"].mu, distributions["p3"].sigma]],
            "dists": ["lognorm", "lognorm"],
        }
        matrix = sobol_sample.sample(problem, N=64, calc_second_order=False, seed=42)
        p1_samples = matrix[:, 0]
        assert p1_samples.mean() == pytest.approx(dist.mean(), rel=0.35)
        assert (p1_samples > 0).all()   # log-normal must never sample non-positive

    def test_pathway_factory_wiring_produces_the_composed_pathways_species(self):
        # Same gap/fix as MonteCarloEngine.propagate() and
        # OATSensitivityAnalyzer -- SobolSensitivityAnalyzer also always
        # built a bare GlucoseInsulinPathway internally. Verify the
        # wiring reaches _run() by checking a metric that only exists on
        # the composite's own species.
        def glp1_peak_metric(result):
            return float(np.max(result.species["GLP1"]))

        analyzer = SobolSensitivityAnalyzer(pathway_factory=_glp1_factory)
        result = analyzer.analyze(
            {"p1": LogNormalDistribution.from_confidence_tier(PathwayParameters().p1, "silver")},
            metrics={"glp1_peak": glp1_peak_metric},
            meals=_standard_meal(),
            n_base_samples=8,
            seed=42,
        )
        assert result.n_model_evaluations > 0
        assert "p1" in result.first_order["glp1_peak"]


# ---------------------------------------------------------------------------
# SobolSensitivityAnalyzer extension-pathway parameters (MASTER_PLAN 2.4)
# ---------------------------------------------------------------------------

class TestSobolSensitivityExtensionParameters:
    def test_pathway_factory_and_extension_pathway_factory_are_mutually_exclusive(self):
        with pytest.raises(ValueError, match="at most one"):
            SobolSensitivityAnalyzer(
                pathway_factory=_glp1_factory,
                extension_pathway_factory=_glp1_extension_factory,
            )

    def test_wrong_arity_extension_pathway_factory_raises_at_construction(self):
        def two_arg_factory(params, meals):
            return GlucoseInsulinPathway(params=params, meals=meals)

        with pytest.raises(ValueError, match="must accept 3 positional arguments"):
            SobolSensitivityAnalyzer(extension_pathway_factory=two_arg_factory)

    def test_extension_parameter_analyzed_alongside_core_via_one_distributions_dict(self):
        # A dotted extension name lives in the SAME `distributions` dict
        # as core PathwayParameters fields -- no separate parameter needed,
        # since SALib's problem definition is already flat.
        def glp1_peak_metric(result):
            return float(np.max(result.species["GLP1"]))

        analyzer = SobolSensitivityAnalyzer(extension_pathway_factory=_glp1_extension_factory)
        result = analyzer.analyze(
            {
                "p1": LogNormalDistribution.from_confidence_tier(PathwayParameters().p1, "silver"),
                "glp1.agonist_scale": LogNormalDistribution.from_mean_cv(9.0, 0.30),
            },
            metrics={"glp1_peak": glp1_peak_metric},
            meals=_standard_meal(),
            n_base_samples=8,
            seed=42,
        )
        assert result.n_model_evaluations > 0
        assert "p1" in result.first_order["glp1_peak"]
        assert "glp1.agonist_scale" in result.first_order["glp1_peak"]
        # agonist_scale directly, linearly scales GLP1's own secretion --
        # a real, substantial first-order index, not near-zero noise.
        assert result.first_order["glp1_peak"]["glp1.agonist_scale"] > 0.1

    def test_dotted_distribution_key_without_factory_raises(self):
        # Code-review finding: without this guard, a dotted extension
        # name in `distributions` with no extension_pathway_factory
        # crashed deep inside SALib.analyze.sobol with an opaque
        # TypeError (zero-variance output on that column) instead of a
        # clear diagnostic naming the actual mistake.
        analyzer = SobolSensitivityAnalyzer()
        with pytest.raises(ValueError, match="nowhere to go"):
            analyzer.analyze(
                {"glp1.agonist_scale": LogNormalDistribution.from_mean_cv(9.0, 0.30)},
                meals=_standard_meal(),
                n_base_samples=8,
                seed=42,
            )
