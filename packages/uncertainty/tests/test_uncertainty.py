# SAVE TO: packages\uncertainty\tests\test_uncertainty.py
"""
Test suite for biosim_uncertainty.

Coverage:
    - Distribution sampling and properties
    - LogNormal: mean recovery, positive samples, CV parameterisation
    - Monte Carlo propagation: CI ordering, uncertainty scaling,
      reproducibility with same seed, different results with different seed
    - Sensitivity analysis: p1 is sensitive for glucose peak,
      scores are finite, top_parameters returns ranked list
    - UQResult: band ordering, summary dict, at_time interpolation
"""

import math
import pytest
import numpy as np

from biosim_uncertainty import (
    ConfidenceBand,
    FixedValue,
    LogNormalDistribution,
    MonteCarloEngine,
    NormalDistribution,
    OATSensitivityAnalyzer,
    STANDARD_METRICS,
    UniformDistribution,
    UQResult,
)
from biosim_ode import PathwayParameters, high_sugar_meal, light_snack, GlucoseInsulinPathway
from biosim_ode.glp1_pathway import GLP1Pathway, GLP1Parameters, glp1_couplings
from biosim_ode.renal_pathway import RenalPathway, RenalClearanceParameters, renal_couplings
from biosim_ode.coupling import CompositePathway


# ---------------------------------------------------------------------------
# Distributions
# ---------------------------------------------------------------------------

class TestLogNormalDistribution:
    def test_mean_recovery(self):
        d = LogNormalDistribution.from_mean_cv(45.0, 0.30)
        rng = np.random.default_rng(42)
        samples = d.sample_n(10000, rng)
        assert np.mean(samples) == pytest.approx(45.0, rel=0.05)

    def test_all_samples_positive(self):
        d = LogNormalDistribution.from_mean_cv(0.0317, 0.30)
        rng = np.random.default_rng(42)
        samples = d.sample_n(1000, rng)
        assert np.all(samples > 0)

    def test_cv_parameterisation(self):
        d = LogNormalDistribution.from_mean_cv(100.0, 0.50)
        rng = np.random.default_rng(0)
        samples = d.sample_n(10000, rng)
        empirical_cv = np.std(samples) / np.mean(samples)
        assert empirical_cv == pytest.approx(0.50, rel=0.1)

    def test_public_mu_sigma_match_private_underscore_versions(self):
        # sobol_sensitivity.py is the first consumer of these values from
        # outside this module -- added public mu/sigma properties rather
        # than reaching into the underscore-prefixed _mu/_sigma from a
        # different module. Confirms they're genuinely the same value,
        # not a second, independently-computed copy that could drift.
        d = LogNormalDistribution.from_mean_cv(45.0, 0.30)
        assert d.mu == d._mu
        assert d.sigma == d._sigma

    def test_from_confidence_tier_gold(self):
        d = LogNormalDistribution.from_confidence_tier(45.0, "gold")
        assert d.cv == pytest.approx(0.10)

    def test_from_confidence_tier_silver(self):
        d = LogNormalDistribution.from_confidence_tier(45.0, "silver")
        assert d.cv == pytest.approx(0.30)

    def test_from_confidence_tier_bronze(self):
        d = LogNormalDistribution.from_confidence_tier(45.0, "bronze")
        assert d.cv == pytest.approx(0.60)

    def test_from_confidence_tier_unknown_raises(self):
        with pytest.raises(ValueError):
            LogNormalDistribution.from_confidence_tier(45.0, "platinum")

    def test_from_confidence_tier_default_additional_cv_matches_prior_behaviour(self):
        d = LogNormalDistribution.from_confidence_tier(45.0, "silver")
        assert d.cv == pytest.approx(0.30)

    def test_from_confidence_tier_combines_additional_cv_in_quadrature(self):
        d = LogNormalDistribution.from_confidence_tier(45.0, "silver", additional_cv=0.4)
        assert d.cv == pytest.approx(math.sqrt(0.30 ** 2 + 0.4 ** 2))

    def test_from_confidence_tier_additional_cv_widens_the_distribution(self):
        base = LogNormalDistribution.from_confidence_tier(45.0, "gold")
        widened = LogNormalDistribution.from_confidence_tier(45.0, "gold", additional_cv=0.25)
        assert widened.cv > base.cv

    def test_negative_mean_raises(self):
        with pytest.raises(ValueError):
            LogNormalDistribution(_mean=-1.0, _cv=0.3)

    def test_zero_cv_raises(self):
        with pytest.raises(ValueError):
            LogNormalDistribution(_mean=1.0, _cv=0.0)

    def test_std_equals_mean_times_cv(self):
        d = LogNormalDistribution.from_mean_cv(50.0, 0.20)
        assert d.std() == pytest.approx(50.0 * 0.20)

    def test_fold_range_increases_with_cv(self):
        d_tight = LogNormalDistribution.from_mean_cv(100.0, 0.10)
        d_wide  = LogNormalDistribution.from_mean_cv(100.0, 0.60)
        assert d_wide.fold_range_95() > d_tight.fold_range_95()

    def test_is_not_degenerate(self):
        d = LogNormalDistribution.from_mean_cv(1.0, 0.1)
        assert not d.is_degenerate()


class TestNormalDistribution:
    def test_mean_recovery(self):
        d = NormalDistribution.from_mean_std(50.0, 5.0)
        rng = np.random.default_rng(42)
        samples = d.sample_n(5000, rng)
        assert np.mean(samples) == pytest.approx(50.0, rel=0.05)

    def test_std_recovery(self):
        d = NormalDistribution.from_mean_std(50.0, 5.0)
        rng = np.random.default_rng(42)
        samples = d.sample_n(5000, rng)
        assert np.std(samples) == pytest.approx(5.0, rel=0.1)

    def test_negative_std_raises(self):
        with pytest.raises(ValueError):
            NormalDistribution(_mean=1.0, _std=-1.0)


class TestUniformDistribution:
    def test_mean(self):
        d = UniformDistribution(lower=0.0, upper=10.0)
        assert d.mean() == pytest.approx(5.0)

    def test_samples_within_bounds(self):
        d = UniformDistribution(lower=1.0, upper=3.0)
        rng = np.random.default_rng(42)
        samples = d.sample_n(1000, rng)
        assert np.all(samples >= 1.0)
        assert np.all(samples <= 3.0)

    def test_lower_ge_upper_raises(self):
        with pytest.raises(ValueError):
            UniformDistribution(lower=5.0, upper=5.0)


class TestFixedValue:
    def test_always_returns_same_value(self):
        d = FixedValue(value=42.0)
        rng = np.random.default_rng(0)
        samples = d.sample_n(100, rng)
        assert np.all(samples == 42.0)

    def test_std_is_zero(self):
        assert FixedValue(value=1.0).std() == 0.0

    def test_is_degenerate(self):
        assert FixedValue(value=1.0).is_degenerate()


# ---------------------------------------------------------------------------
# ConfidenceBand
# ---------------------------------------------------------------------------

class TestConfidenceBand:
    @pytest.fixture
    def simple_band(self):
        t = np.linspace(0, 100, 50)
        return ConfidenceBand(
            lower=np.ones(50) * 0.8,
            median=np.ones(50) * 1.0,
            upper=np.ones(50) * 1.2,
        )

    def test_width(self, simple_band):
        assert np.allclose(simple_band.width, 0.4)

    def test_relative_width(self, simple_band):
        assert np.allclose(simple_band.relative_width, 0.4)

    def test_at_time_index(self, simple_band):
        lower, median, upper = simple_band.at_time_index(0)
        assert lower == pytest.approx(0.8)
        assert median == pytest.approx(1.0)
        assert upper == pytest.approx(1.2)

    def test_mean_relative_width(self, simple_band):
        assert simple_band.mean_relative_width() == pytest.approx(0.4)

    def test_mean_relative_width_excludes_near_zero_median_points(self):
        # Regression: relative_width's own median!=0 guard only special-
        # cases EXACT zero -- a median that merely passes NEAR zero (e.g.
        # a decaying bookkeeping compartment like GlucoseInsulinPathway's
        # X species) produces an enormous relative_width at those points
        # from an ordinary-sized numerator over a tiny denominator, which
        # then dominates the mean. Reproduced directly on real X-species
        # trajectories before this fix: mean_relative_width() ~0.25 vs
        # ~0.03-0.08 for physiologically real species, and
        # most_uncertain_species() answered "X" for every ordinary
        # scenario -- a bookkeeping artifact, not a genuine signal.
        # Construct a band that decays from 10 down near 0 (mimicking X):
        # most points are well above the floor, a few late points sit at
        # a tiny median with a real (non-tiny) width -- those should be
        # excluded, not allowed to dominate the average.
        n = 20
        median = np.concatenate([np.linspace(10.0, 1.0, n - 3), np.array([0.05, 0.02, 0.001])])
        width = np.ones(n) * 0.5   # constant absolute width throughout
        band = ConfidenceBand(lower=median - width / 2, median=median, upper=median + width / 2)
        # Without the floor, the tail's tiny-median points would dominate
        # (width/median blows up as median -> 0). With a 5% peak floor
        # (peak=10, floor=0.5), only median >= 0.5 counts -- the mean
        # relative width should stay small and bounded, not blow up.
        assert band.mean_relative_width() < 1.0

    def test_mean_relative_width_floor_is_relative_to_this_bands_own_peak(self):
        # Two bands with identical SHAPE but different absolute scale
        # must get the same mean_relative_width() -- the floor is a
        # fraction of each band's own peak, not a fixed absolute cutoff.
        n = 20
        shape = np.concatenate([np.linspace(10.0, 1.0, n - 3), np.array([0.05, 0.02, 0.001])])
        width_frac = 0.05
        band_a = ConfidenceBand(
            lower=shape * (1 - width_frac), median=shape, upper=shape * (1 + width_frac)
        )
        scaled = shape * 1000.0
        band_b = ConfidenceBand(
            lower=scaled * (1 - width_frac), median=scaled, upper=scaled * (1 + width_frac)
        )
        assert band_a.mean_relative_width() == pytest.approx(band_b.mean_relative_width(), rel=1e-6)


# ---------------------------------------------------------------------------
# Monte Carlo propagation
# ---------------------------------------------------------------------------

class TestMonteCarlo:
    """
    These tests use n_samples=30 for speed.
    Production should use ≥200.
    """

    @pytest.fixture(scope="class")
    def uq_result(self):
        engine = MonteCarloEngine(n_samples=30, random_seed=42)
        distributions = {
            "p1": LogNormalDistribution.from_mean_cv(0.0317, 0.30),
            "ke": LogNormalDistribution.from_mean_cv(0.093,  0.30),
        }
        return engine.propagate(
            distributions=distributions,
            meals=[high_sugar_meal()],
            duration_s=7200,
        )

    def test_propagation_succeeds(self, uq_result):
        assert uq_result.n_successful > 0

    def test_all_species_have_bands(self, uq_result):
        for name in ["G", "I", "X", "IGF1", "mTOR", "AMPK"]:
            assert name in uq_result.bands

    def test_lower_le_median_le_upper(self, uq_result):
        for name, band in uq_result.bands.items():
            assert np.all(band.lower <= band.median + 1e-10), f"{name}: lower > median"
            assert np.all(band.median <= band.upper + 1e-10), f"{name}: median > upper"

    def test_glucose_peak_above_baseline(self, uq_result):
        G_band = uq_result.band("G")
        G_baseline_uM = 5.0 * 1000.0
        assert np.max(G_band.median) > G_baseline_uM

    def test_same_seed_gives_same_result(self):
        distributions = {
            "p1": LogNormalDistribution.from_mean_cv(0.0317, 0.30),
        }
        e1 = MonteCarloEngine(n_samples=20, random_seed=99)
        e2 = MonteCarloEngine(n_samples=20, random_seed=99)
        r1 = e1.propagate(distributions, meals=[high_sugar_meal()], duration_s=3600)
        r2 = e2.propagate(distributions, meals=[high_sugar_meal()], duration_s=3600)
        assert np.allclose(r1.band("G").median, r2.band("G").median)

    def test_different_seed_gives_different_result(self):
        distributions = {
            "p1": LogNormalDistribution.from_mean_cv(0.0317, 0.50),
        }
        e1 = MonteCarloEngine(n_samples=20, random_seed=1)
        e2 = MonteCarloEngine(n_samples=20, random_seed=2)
        r1 = e1.propagate(distributions, meals=[high_sugar_meal()], duration_s=3600)
        r2 = e2.propagate(distributions, meals=[high_sugar_meal()], duration_s=3600)
        assert not np.allclose(r1.band("G").median, r2.band("G").median)

    def test_higher_cv_gives_wider_ci(self):
        meals = [high_sugar_meal()]
        tight = MonteCarloEngine(n_samples=20, random_seed=42).propagate(
            {"p1": LogNormalDistribution.from_mean_cv(0.0317, 0.05)},
            meals=meals, duration_s=3600,
        )
        wide = MonteCarloEngine(n_samples=20, random_seed=42).propagate(
            {"p1": LogNormalDistribution.from_mean_cv(0.0317, 0.60)},
            meals=meals, duration_s=3600,
        )
        tight_width = np.mean(tight.band("G").width)
        wide_width = np.mean(wide.band("G").width)
        assert wide_width > tight_width

    def test_fixed_value_gives_zero_band_width(self):
        engine = MonteCarloEngine(n_samples=20, random_seed=42)
        result = engine.propagate(
            {"p1": FixedValue(0.0317), "ke": FixedValue(0.093)},
            meals=[high_sugar_meal()], duration_s=3600,
        )
        # With all params fixed, band width should be near zero
        G_width = np.mean(result.band("G").width)
        assert G_width < 1.0  # less than 1 μM spread from numerical noise

    def test_n_samples_too_small_raises(self):
        with pytest.raises(ValueError, match="n_samples"):
            MonteCarloEngine(n_samples=5)

    def test_store_trajectories_captures_sampled_parameter_values(self):
        # Regression: propagate() sampled a fresh PathwayParameters per
        # draw, used it to run the ODE, then discarded it -- only the
        # trajectory survived (in raw_trajectories). BayesianCalibrator
        # (calibration.py) needs both the trajectory AND the exact
        # parameter draw that produced it, to resample parameter draws
        # by posterior weight (Sampling-Importance-Resampling).
        engine = MonteCarloEngine(n_samples=20, random_seed=42)
        distributions = {
            "p1": LogNormalDistribution.from_mean_cv(0.0317, 0.30),
            "p3": LogNormalDistribution.from_mean_cv(9.032e-7, 0.30),
        }
        result = engine.propagate(
            distributions, meals=[high_sugar_meal()], duration_s=3600,
            store_trajectories=True,
        )
        assert result.sampled_parameter_values is not None
        assert set(result.sampled_parameter_values.keys()) == {"p1", "p3"}
        for name in ("p1", "p3"):
            values = result.sampled_parameter_values[name]
            assert len(values) == result.n_successful
            assert np.all(values > 0)   # log-normal, always positive
            assert len(set(values.tolist())) > 1   # genuinely varies across samples

    def test_sampled_parameter_values_none_without_store_trajectories(self):
        engine = MonteCarloEngine(n_samples=20, random_seed=42)
        result = engine.propagate(
            {"p1": LogNormalDistribution.from_mean_cv(0.0317, 0.30)},
            meals=[high_sugar_meal()], duration_s=3600,
            store_trajectories=False,
        )
        assert result.sampled_parameter_values is None

    def test_failed_samples_record_a_reason(self, monkeypatch):
        # Regression: a failed sample's exception used to be silently
        # discarded (except Exception: n_failed += 1; continue), only a
        # bare count survived -- a real solver bug and a few benign
        # edge-case draws were indistinguishable. Same bug class as
        # CohortSummary.failure_reasons.
        import biosim_ode.system as system_module

        def _raise(self, duration_s, meals=None, n_output_points=500, y0_override=None,
                   check_clamps=True):
            raise ValueError("synthetic ODE solver failure for testing")

        monkeypatch.setattr(system_module.ODESystem, "simulate", _raise)
        engine = MonteCarloEngine(n_samples=10, random_seed=42)
        result = engine.propagate(
            {"p1": FixedValue(0.0317)}, meals=[high_sugar_meal()], duration_s=3600,
        )
        assert result.n_failed == 10
        assert len(result.failure_reasons) == 10
        assert all("synthetic ODE solver failure" in r for r in result.failure_reasons)

    def test_at_time_returns_all_species(self, uq_result):
        state = uq_result.at_time(1800.0)
        for name in ["G", "I", "mTOR", "AMPK"]:
            assert name in state
            lower, median, upper = state[name]
            assert lower <= median + 1e-10
            assert median <= upper + 1e-10

    def test_success_rate(self, uq_result):
        assert uq_result.success_rate > 0.8

    def test_summary_dict(self, uq_result):
        s = uq_result.uncertainty_summary()
        assert "n_samples" in s
        assert "success_rate" in s
        assert "species_uncertainty" in s

    def test_most_uncertain_species_returns_string(self, uq_result):
        name = uq_result.most_uncertain_species()
        assert name in uq_result.bands


# ---------------------------------------------------------------------------
# pathway_factory -- UQ vs. the point estimate's composed pathway
#
# Regression: propagate() (and OATSensitivityAnalyzer/SobolSensitivityAnalyzer)
# always built a bare GlucoseInsulinPathway, with no way to supply the same
# CompositePathway (GLP1, renal, cortisol, hepatic, ...) SimulationRunner
# composes for the point estimate. For a semaglutide/dapagliflozin scenario
# (both have ZERO static PathwayParameters effect), the reported "confidence
# band" was PROVABLY a no-drug band -- the sampled PathwayParameters were
# identical either way. pathway_factory closes this gap, defaulting to None
# (unchanged, bit-for-bit identical prior behaviour) for every existing
# caller.
# ---------------------------------------------------------------------------

def _glp1_factory(params, meals):
    """Test factory: composes GLP1Pathway (semaglutide-like) with GI."""
    gi = GlucoseInsulinPathway(params=params, meals=meals)
    glp1 = GLP1Pathway()
    return CompositePathway([gi, glp1], couplings=glp1_couplings(gi.name, glp1.name))


class TestMonteCarloPathwayFactory:
    def test_default_none_matches_bare_gi_exactly(self):
        # Explicitly passing pathway_factory=None must be byte-for-byte
        # identical to omitting it -- the parameter's whole point is that
        # every existing caller (calibration.py, validator.py's OGTT/
        # Bergman-IV checks, every test in this file above) is completely
        # unaffected.
        distributions = {"p1": LogNormalDistribution.from_mean_cv(0.0317, 0.30)}
        meals = [high_sugar_meal()]
        r1 = MonteCarloEngine(n_samples=20, random_seed=7).propagate(
            distributions, meals=meals, duration_s=3600,
        )
        r2 = MonteCarloEngine(n_samples=20, random_seed=7).propagate(
            distributions, meals=meals, duration_s=3600, pathway_factory=None,
        )
        assert np.array_equal(r1.band("G").median, r2.band("G").median)
        assert set(r1.bands) == set(r2.bands) == {"G", "I", "X", "IGF1", "mTOR", "AMPK"}

    def test_pathway_factory_adds_the_composed_pathways_own_species(self):
        distributions = {"p1": LogNormalDistribution.from_mean_cv(0.0317, 0.10)}
        result = MonteCarloEngine(n_samples=15, random_seed=1).propagate(
            distributions, meals=[high_sugar_meal()], duration_s=3600,
            pathway_factory=_glp1_factory,
        )
        # Structural proof: the composite's own species (GLP1) now has a
        # real band, not silently absent.
        assert "GLP1" in result.bands
        assert {"G", "I", "X", "IGF1", "mTOR", "AMPK"}.issubset(result.bands)

    def test_pathway_factory_exact_convergence_with_fixed_values(self):
        # Mirrors test_fixed_value_gives_zero_band_width's convention,
        # through a composed pathway: with every sampled parameter fixed,
        # the ensemble is degenerate and the band should collapse to the
        # point estimate regardless of which pathway structure is used.
        result = MonteCarloEngine(n_samples=15, random_seed=42).propagate(
            {"p1": FixedValue(0.0317), "ke": FixedValue(0.093)},
            meals=[high_sugar_meal()], duration_s=3600,
            pathway_factory=_glp1_factory,
        )
        G_width = np.mean(result.band("G").width)
        assert G_width < 1.0

    def test_pathway_factory_produces_a_genuinely_different_band_than_bare_gi(self):
        # The decisive proof this fix closes the real gap: before it, a
        # drug-active composite scenario's band was provably identical to
        # a drug-free bare-GI ensemble's whenever the drug has zero static
        # PathwayParameters effect (true for GLP-1/SGLT2i). With the same
        # distributions/meals/seed, routing through the GLP1-composed
        # factory must now produce a measurably different G band than the
        # bare-GI default -- the ensemble genuinely sees the drug's
        # dynamics (incretin amplification, gastric-emptying slowing),
        # not just static params it doesn't have.
        distributions = {"p1": LogNormalDistribution.from_mean_cv(0.0317, 0.10)}
        meals = [high_sugar_meal()]
        bare = MonteCarloEngine(n_samples=15, random_seed=3).propagate(
            distributions, meals=meals, duration_s=3600,
        )
        composed = MonteCarloEngine(n_samples=15, random_seed=3).propagate(
            distributions, meals=meals, duration_s=3600, pathway_factory=_glp1_factory,
        )
        assert not np.allclose(bare.band("G").median, composed.band("G").median)

    def test_pathway_factory_nondeterministic_species_set_raises_clear_error(self):
        # A pathway_factory that varies structure by something OTHER than
        # (params, meals) violates propagate()'s determinism contract --
        # must fail with a clear, named cause, not a bare KeyError deep in
        # the trajectory-storing loop.
        calls = {"n": 0}

        def flaky_factory(params, meals):
            calls["n"] += 1
            gi = GlucoseInsulinPathway(params=params, meals=meals)
            if calls["n"] == 1:
                return gi
            glp1 = GLP1Pathway()
            return CompositePathway([gi, glp1], couplings=glp1_couplings(gi.name, glp1.name))

        with pytest.raises(ValueError, match="pathway_factory must be deterministic"):
            MonteCarloEngine(n_samples=10, random_seed=1).propagate(
                {"p1": LogNormalDistribution.from_mean_cv(0.0317, 0.10)},
                meals=[high_sugar_meal()], duration_s=3600,
                pathway_factory=flaky_factory,
            )


# ---------------------------------------------------------------------------
# UQResult.species_pathway provenance (MASTER_PLAN 2.3)
# ---------------------------------------------------------------------------
#
# Regression: UQResult.bands had no way to say which pathway produced a
# given species -- purely presentational (species names are already
# globally unique), but useful for a UI grouping bands by pathway instead
# of listing every one flat.

class TestUQResultSpeciesPathway:
    def test_bare_gi_maps_every_species_to_its_own_pathway_name(self):
        result = MonteCarloEngine(n_samples=10, random_seed=1).propagate(
            {"p1": LogNormalDistribution.from_mean_cv(0.0317, 0.10)},
            meals=[high_sugar_meal()], duration_s=3600,
        )
        assert set(result.species_pathway) == set(result.bands)
        for species in result.bands:
            assert result.species_pathway[species] == "glucose_insulin_v1"
            assert result.pathway_for(species) == "glucose_insulin_v1"

    def test_composed_pathway_attributes_each_species_to_its_own_sub_pathway(self):
        result = MonteCarloEngine(n_samples=10, random_seed=1).propagate(
            {"p1": LogNormalDistribution.from_mean_cv(0.0317, 0.10)},
            meals=[high_sugar_meal()], duration_s=3600,
            pathway_factory=_glp1_factory,
        )
        assert result.species_pathway["GLP1"] == "glp1_v1"
        assert result.species_pathway["G"] == "glucose_insulin_v1"
        assert result.species_pathway["I"] == "glucose_insulin_v1"

    def test_all_samples_failing_leaves_species_pathway_empty(self):
        def always_raises(params, meals):
            raise RuntimeError("deliberately broken factory")

        result = MonteCarloEngine(n_samples=10, random_seed=1).propagate(
            {"p1": LogNormalDistribution.from_mean_cv(0.0317, 0.10)},
            meals=[high_sugar_meal()], duration_s=3600,
            pathway_factory=always_raises,
        )
        assert result.species_pathway == {}
        assert result.pathway_for("G") is None

    def test_hand_constructed_uqresult_defaults_to_empty_mapping(self):
        result = UQResult(
            times_s=np.array([0.0, 1.0]),
            bands={},
            n_samples=1,
            n_failed=0,
            overall_confidence=1.0,
        )
        assert result.species_pathway == {}
        assert result.pathway_for("G") is None

    def test_uncertainty_summary_includes_species_pathway(self):
        result = MonteCarloEngine(n_samples=10, random_seed=1).propagate(
            {"p1": LogNormalDistribution.from_mean_cv(0.0317, 0.10)},
            meals=[high_sugar_meal()], duration_s=3600,
            pathway_factory=_glp1_factory,
        )
        summary = result.uncertainty_summary()
        assert summary["species_pathway"]["GLP1"] == "glp1_v1"


# ---------------------------------------------------------------------------
# Extension-pathway parameter sampling (MASTER_PLAN 2.2)
# ---------------------------------------------------------------------------
#
# Regression: pathway_factory (above) made the UQ ensemble see the right
# pathway STRUCTURE, but extension-pathway-specific parameters (e.g.
# RenalClearanceParameters.beta) stayed fixed at their point-estimate value
# in every sample -- a drug scenario's confidence band correctly reflected
# uncertainty in p1/p2/p3/ke/EC50_mTOR but silently treated the drug's own
# calibrated effect size as if it were exactly known.

def _glp1_extension_factory(params, meals, overrides):
    """Test factory: samples GLP1's agonist_scale via a dotted key."""
    gi = GlucoseInsulinPathway(params=params, meals=meals)
    agonist_scale = overrides.get("glp1.agonist_scale", 9.0)
    glp1 = GLP1Pathway(GLP1Parameters(agonist_scale=agonist_scale))
    return CompositePathway([gi, glp1], couplings=glp1_couplings(gi.name, glp1.name))


def _renal_extension_factory(params, meals, overrides):
    """Test factory: samples RenalClearanceParameters.beta via a dotted
    key -- beta's own __init__ validation ([0, 1)) is the real-world
    example of an out-of-range extension draw reaching constructor
    validation instead of being silently floored."""
    gi = GlucoseInsulinPathway(params=params, meals=meals)
    beta = overrides.get("renal.beta", 0.40)
    renal = RenalPathway(RenalClearanceParameters(beta=beta))
    return CompositePathway([gi, renal], couplings=renal_couplings(gi.name, renal))


class TestMonteCarloExtensionParameters:
    def test_pathway_factory_and_extension_pathway_factory_are_mutually_exclusive(self):
        with pytest.raises(ValueError, match="at most one"):
            MonteCarloEngine(n_samples=10, random_seed=1).propagate(
                {"p1": LogNormalDistribution.from_mean_cv(0.0317, 0.10)},
                meals=[high_sugar_meal()], duration_s=3600,
                pathway_factory=_glp1_factory,
                extension_pathway_factory=_glp1_extension_factory,
            )

    def test_extension_distributions_without_factory_raises(self):
        with pytest.raises(ValueError, match="nowhere to go"):
            MonteCarloEngine(n_samples=10, random_seed=1).propagate(
                {"p1": LogNormalDistribution.from_mean_cv(0.0317, 0.10)},
                meals=[high_sugar_meal()], duration_s=3600,
                extension_distributions={
                    "glp1.agonist_scale": LogNormalDistribution.from_mean_cv(9.0, 0.10),
                },
            )

    def test_colliding_extension_key_raises(self):
        with pytest.raises(ValueError, match="collide"):
            MonteCarloEngine(n_samples=10, random_seed=1).propagate(
                {"p1": LogNormalDistribution.from_mean_cv(0.0317, 0.10)},
                meals=[high_sugar_meal()], duration_s=3600,
                extension_distributions={"p1": LogNormalDistribution.from_mean_cv(0.0317, 0.10)},
                extension_pathway_factory=_glp1_extension_factory,
            )

    def test_extension_distributions_default_matches_no_extensions_exactly(self):
        # extension_distributions=None/extension_pathway_factory=None (the
        # default) must reproduce bit-identical output to a call that
        # never mentions either -- every existing caller is unaffected.
        distributions = {"p1": LogNormalDistribution.from_mean_cv(0.0317, 0.10)}
        meals = [high_sugar_meal()]
        r1 = MonteCarloEngine(n_samples=15, random_seed=7).propagate(
            distributions, meals=meals, duration_s=3600,
        )
        r2 = MonteCarloEngine(n_samples=15, random_seed=7).propagate(
            distributions, meals=meals, duration_s=3600,
            extension_distributions=None, extension_pathway_factory=None,
        )
        assert np.array_equal(r1.band("G").median, r2.band("G").median)

    def test_sampling_an_extension_parameter_genuinely_widens_its_own_band(self):
        # The decisive proof: with agonist_scale fixed (via pathway_factory,
        # every sample uses the identical value), the GLP1 band should be
        # narrow; with agonist_scale genuinely SAMPLED via
        # extension_pathway_factory, the GLP1 band must be measurably
        # wider, at the same core distributions/seed.
        core = {"p1": LogNormalDistribution.from_mean_cv(0.0317, 0.05)}
        meals = [high_sugar_meal()]

        fixed = MonteCarloEngine(n_samples=20, random_seed=3).propagate(
            core, meals=meals, duration_s=3600, pathway_factory=_glp1_extension_factory_fixed,
        )
        sampled = MonteCarloEngine(n_samples=20, random_seed=3).propagate(
            core, meals=meals, duration_s=3600,
            extension_distributions={
                "glp1.agonist_scale": LogNormalDistribution.from_mean_cv(9.0, 0.30),
            },
            extension_pathway_factory=_glp1_extension_factory,
        )
        assert sampled.band("GLP1").mean_relative_width() > fixed.band("GLP1").mean_relative_width()

    def test_extension_draw_recorded_under_its_dotted_key(self):
        result = MonteCarloEngine(n_samples=15, random_seed=3).propagate(
            {"p1": LogNormalDistribution.from_mean_cv(0.0317, 0.10)},
            meals=[high_sugar_meal()], duration_s=3600,
            store_trajectories=True,
            extension_distributions={
                "glp1.agonist_scale": LogNormalDistribution.from_mean_cv(9.0, 0.30),
            },
            extension_pathway_factory=_glp1_extension_factory,
        )
        assert "glp1.agonist_scale" in result.parameter_names
        assert "glp1.agonist_scale" in result.sampled_parameter_values
        assert len(result.sampled_parameter_values["glp1.agonist_scale"]) == result.n_successful

    def test_wrong_arity_extension_pathway_factory_raises_once_at_call_entry(self):
        # Security-review finding: without an eager arity check, a
        # wrong-signature factory raises the same TypeError on every
        # single sample (n_samples times) instead of failing once, up
        # front, with a clear diagnostic.
        def two_arg_factory(params, meals):
            return GlucoseInsulinPathway(params=params, meals=meals)

        with pytest.raises(ValueError, match="must accept 3 positional arguments"):
            MonteCarloEngine(n_samples=10, random_seed=1).propagate(
                {"p1": LogNormalDistribution.from_mean_cv(0.0317, 0.10)},
                meals=[high_sugar_meal()], duration_s=3600,
                extension_pathway_factory=two_arg_factory,
            )

    def test_out_of_range_extension_draw_fails_with_a_reasoned_failure_not_a_crash(self):
        # RenalClearanceParameters rejects beta outside [0, 1) -- a wide
        # enough distribution on renal.beta should produce some draws
        # that fail construction, counted honestly in failure_reasons,
        # rather than being silently floored the way a core kinetic
        # parameter's positivity floor would.
        result = MonteCarloEngine(n_samples=40, random_seed=11).propagate(
            {"p1": LogNormalDistribution.from_mean_cv(0.0317, 0.05)},
            meals=[high_sugar_meal()], duration_s=3600,
            extension_distributions={
                "renal.beta": NormalDistribution.from_mean_std(0.5, 0.5),
            },
            extension_pathway_factory=_renal_extension_factory,
        )
        assert result.n_failed > 0
        assert any("beta must be in [0, 1)" in reason for reason in result.failure_reasons)
        assert result.n_successful > 0


def _glp1_extension_factory_fixed(params, meals):
    """2-arg comparator for the widening test above: agonist_scale held
    fixed at the same value _glp1_extension_factory() defaults to."""
    return _glp1_extension_factory(params, meals, {})


# ---------------------------------------------------------------------------
# Sensitivity analysis
# ---------------------------------------------------------------------------

class TestOATSensitivity:
    @pytest.fixture(scope="class")
    def sensitivity_result(self):
        analyzer = OATSensitivityAnalyzer(perturbation_fraction=0.10)
        return analyzer.analyze(
            parameter_names=["p1", "ke", "EC50_mTOR", "k_mTOR"],
            metrics={
                "glucose_peak": STANDARD_METRICS["glucose_peak"],
                "mtor_peak": STANDARD_METRICS["mtor_peak"],
            },
            meals=[high_sugar_meal()],
            duration_s=7200,
        )

    def test_all_params_have_scores(self, sensitivity_result):
        for param in ["p1", "ke", "EC50_mTOR", "k_mTOR"]:
            assert param in sensitivity_result.sensitivity_indices["glucose_peak"]

    def test_p1_is_sensitive_for_glucose(self, sensitivity_result):
        s = sensitivity_result.sensitivity_indices["glucose_peak"]["p1"]
        assert abs(s) > 0.05, f"p1 sensitivity for glucose_peak too low: {s}"

    def test_ec50_mtor_sensitive_for_mtor_peak(self, sensitivity_result):
        s = sensitivity_result.sensitivity_indices["mtor_peak"]["EC50_mTOR"]
        assert abs(s) > 0.05

    def test_top_parameters_returns_ranked(self, sensitivity_result):
        top = sensitivity_result.top_parameters("glucose_peak", n=3)
        assert len(top) == 3
        scores = [abs(s) for _, s in top]
        assert scores == sorted(scores, reverse=True)

    def test_top_parameters_unknown_metric_raises(self, sensitivity_result):
        with pytest.raises(KeyError):
            sensitivity_result.top_parameters("nonexistent_metric")

    def test_base_outputs_positive(self, sensitivity_result):
        for name, val in sensitivity_result.base_outputs.items():
            assert val > 0, f"Base output {name} = {val} is not positive"

    def test_failed_perturbations_record_a_reason(self, monkeypatch):
        # Regression: _run() used to return a bare None on failure --
        # SensitivityResult.n_failed couldn't distinguish a real solver
        # bug from a benign edge-case perturbation. Same bug class as
        # CohortSummary.failure_reasons/UQResult.failure_reasons. Let the
        # baseline call succeed (real ODESystem.simulate), then fail
        # every perturbed call so failure_reasons is populated without
        # raising the "baseline failed" RuntimeError.
        import biosim_ode.system as system_module

        real_simulate = system_module.ODESystem.simulate
        call_count = {"n": 0}

        def _fail_after_baseline(self, duration_s, meals=None, n_output_points=500,
                                 y0_override=None, check_clamps=True):
            call_count["n"] += 1
            if call_count["n"] == 1:
                return real_simulate(
                    self, duration_s, meals=meals, n_output_points=n_output_points,
                    y0_override=y0_override, check_clamps=check_clamps,
                )
            raise ValueError("synthetic perturbation failure for testing")

        monkeypatch.setattr(system_module.ODESystem, "simulate", _fail_after_baseline)
        analyzer = OATSensitivityAnalyzer(perturbation_fraction=0.10)
        result = analyzer.analyze(
            parameter_names=["p1", "ke"],
            meals=[high_sugar_meal()],
            duration_s=3600,
        )
        assert result.n_failed == 2
        assert len(result.failure_reasons) == 4  # upper + lower per parameter
        assert all("synthetic perturbation failure" in r for r in result.failure_reasons)

    def test_summary_structure(self, sensitivity_result):
        s = sensitivity_result.summary()
        assert "top_sensitivities" in s
        assert "glucose_peak" in s["top_sensitivities"]

    def test_invalid_perturbation_raises(self):
        with pytest.raises(ValueError):
            OATSensitivityAnalyzer(perturbation_fraction=1.5)

    def test_pathway_factory_wiring_produces_the_composed_pathways_species(self):
        # Same gap/fix as MonteCarloEngine.propagate() (see
        # TestMonteCarloPathwayFactory above) -- OATSensitivityAnalyzer
        # also always built a bare GlucoseInsulinPathway internally, with
        # no way to analyze sensitivity through a composed extension
        # pathway. Verify the wiring reaches _run() by checking a metric
        # that only exists on the composite's own species.
        def glp1_peak_metric(result):
            return float(np.max(result.species["GLP1"]))

        analyzer = OATSensitivityAnalyzer(
            perturbation_fraction=0.10, pathway_factory=_glp1_factory,
        )
        result = analyzer.analyze(
            parameter_names=["p1"],
            metrics={"glp1_peak": glp1_peak_metric},
            meals=[high_sugar_meal()],
            duration_s=3600,
        )
        assert result.n_failed == 0
        assert result.base_outputs["glp1_peak"] > 0


# ---------------------------------------------------------------------------
# OATSensitivityAnalyzer extension-pathway parameters (MASTER_PLAN 2.4)
# ---------------------------------------------------------------------------
#
# Regression: OATSensitivityAnalyzer's own pathway_factory fix (above) made
# it see the right composed pathway STRUCTURE, but extension-pathway-
# specific parameters (e.g. RenalClearanceParameters.beta) had no way to
# be perturbed/analyzed at all -- parameter_names could only ever name
# PathwayParameters fields. Same additive treatment as
# MonteCarloEngine.propagate() (MASTER_PLAN 2.2).

class TestOATSensitivityExtensionParameters:
    def test_pathway_factory_and_extension_pathway_factory_are_mutually_exclusive(self):
        with pytest.raises(ValueError, match="at most one"):
            OATSensitivityAnalyzer(
                pathway_factory=_glp1_factory,
                extension_pathway_factory=_glp1_extension_factory,
            )

    def test_wrong_arity_extension_pathway_factory_raises_at_construction(self):
        def two_arg_factory(params, meals):
            return GlucoseInsulinPathway(params=params, meals=meals)

        with pytest.raises(ValueError, match="must accept 3 positional arguments"):
            OATSensitivityAnalyzer(extension_pathway_factory=two_arg_factory)

    def test_extension_parameter_gets_a_real_sensitivity_score(self):
        def glp1_peak_metric(result):
            return float(np.max(result.species["GLP1"]))

        analyzer = OATSensitivityAnalyzer(
            perturbation_fraction=0.20,
            extension_pathway_factory=_glp1_extension_factory,
            extension_base_values={"glp1.agonist_scale": 9.0},
        )
        result = analyzer.analyze(
            parameter_names=["glp1.agonist_scale"],
            metrics={"glp1_peak": glp1_peak_metric},
            meals=[high_sugar_meal()],
            duration_s=3600,
        )
        assert result.n_failed == 0
        # agonist_scale directly, linearly scales GLP1's own secretion --
        # a real, substantial sensitivity, not a near-zero numerical
        # artifact.
        assert abs(result.sensitivity_indices["glp1_peak"]["glp1.agonist_scale"]) > 0.5

    def test_core_and_extension_parameters_analyzed_together(self):
        # A single analyze() call must handle a mix of core PathwayParameters
        # fields and dotted extension names in the same parameter_names list.
        analyzer = OATSensitivityAnalyzer(
            perturbation_fraction=0.10,
            extension_pathway_factory=_glp1_extension_factory,
            extension_base_values={"glp1.agonist_scale": 9.0},
        )
        result = analyzer.analyze(
            parameter_names=["p1", "glp1.agonist_scale"],
            meals=[high_sugar_meal()],
            duration_s=3600,
        )
        assert result.n_failed == 0
        assert "p1" in result.sensitivity_indices["glucose_peak"]
        assert "glp1.agonist_scale" in result.sensitivity_indices["glucose_peak"]

    def test_extension_parameter_not_in_base_values_is_silently_skipped(self):
        # Mirrors the existing "unknown core field -> skip" convention
        # (base_val is None -> continue) rather than raising -- a caller
        # analyzing a broad parameter list shouldn't have to filter it
        # per-factory first.
        analyzer = OATSensitivityAnalyzer(
            perturbation_fraction=0.10,
            extension_pathway_factory=_glp1_extension_factory,
            extension_base_values={},  # "glp1.agonist_scale" not present
        )
        result = analyzer.analyze(
            parameter_names=["p1", "glp1.agonist_scale"],
            meals=[high_sugar_meal()],
            duration_s=3600,
        )
        assert "glp1.agonist_scale" not in result.sensitivity_indices["glucose_peak"]
        assert "p1" in result.sensitivity_indices["glucose_peak"]

    def test_extension_base_values_without_factory_raises(self):
        # Code-review finding: without this guard, extension_base_values
        # supplied with no extension_pathway_factory silently produced a
        # misleading "zero influence" result -- the perturbed value never
        # reached the pathway at all (upper/lower runs were identical).
        with pytest.raises(ValueError, match="nowhere to go"):
            OATSensitivityAnalyzer(extension_base_values={"glp1.agonist_scale": 9.0})
