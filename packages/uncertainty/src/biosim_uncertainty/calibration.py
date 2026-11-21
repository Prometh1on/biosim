# SAVE TO: packages\uncertainty\src\biosim_uncertainty\calibration.py
"""
Bayesian parameter calibration against clinical benchmarks.

`ModelValidator` (biosim_validation) checks whether a SINGLE point-estimate
simulation run falls within a fixed tolerance band of each clinical
reference value -- PASS/WARN/FAIL, nothing more. It never uses the
clinical data to update beliefs about which parameter values are more or
less plausible. `BayesianCalibrator` does: given prior distributions over
uncertain parameters (the same `ParameterDistribution` objects
`MonteCarloEngine` already uses) and a set of clinical reference values
with known tolerance, it produces a POSTERIOR sample set -- parameter
draws (and the resulting output trajectories) reweighted by how well each
one's simulated output matches the real clinical data.

Method: Sampling-Importance-Resampling (SIR), not a from-scratch MCMC
sampler. Chosen deliberately for this problem's specific shape (~5-10
uncertain parameters, ~5-20ms per ODE solve, no MCMC library installed in
this project's venv -- only scipy+numpy):
    - SIR reweights samples `MonteCarloEngine.propagate()` already draws
      from the prior -- zero extra ODE solves beyond the existing forward
      Monte Carlo run.
    - A hand-rolled Metropolis-Hastings sampler (the realistic scipy-only
      alternative) needs proposal-covariance tuning, burn-in, and
      convergence diagnostics (Gelman-Rubin, effective sample size) with
      no off-the-shelf library -- meaningfully more engineering for a
      benefit (better mixing when the posterior is much narrower than the
      prior in several correlated dimensions) that matters more at higher
      parameter counts than this problem has.
    - Main risk, tracked explicitly below: importance-weight degeneracy.
      If the likelihood is much narrower than the prior, effective sample
      size (ESS) collapses toward a handful of duplicated points. Always
      report Kish's ESS and flag `is_degenerate` rather than silently
      returning a collapsed posterior.

Likelihood -- a stated MODELING CHOICE, not a measured quantity (same
honesty standard as `HPA_FEEDBACK_GAIN` in `cortisol_pathway.py` and the
GOLD/SILVER/BRONZE confidence-tier CVs in `distributions.py`): each
scored calibration target's `tolerance_fraction` (0.10/0.25/0.50 for
TIGHT/NORMAL/WIDE in biosim_validation's policy) is treated as an implied
Gaussian standard deviation, `sigma = reference * tolerance_fraction /
likelihood_sigma_z`. These tolerance bands are a coarse, pre-registered
PASS/FAIL policy across benchmarks of very different real evidentiary
strength, not per-study confidence intervals (this project's own
metformin benchmark is direct proof: its real 95% CI, 64-307 mmol.min/L
around a reference of 186, is WIDER than the WIDE tier's own +/-50%
policy band around that same reference -- see CLAUDE.md). Treating
`tolerance_fraction` as a literal Gaussian std therefore risks
manufacturing false precision. Mitigate by (a) always checking
`effective_sample_size_fraction` before trusting a posterior, and (b)
validating any calibration with a posterior-predictive check against a
DIFFERENT benchmark that contributed zero weight to the likelihood, not
by declaring the calibrated posterior correct on the strength of the
benchmarks it was fit to.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Optional, Protocol, runtime_checkable

import numpy as np
from scipy.special import logsumexp

from biosim_ode.meal import MealInput

from .distributions import ParameterDistribution
from .propagation import MonteCarloEngine, PathwayFactory
from .result import ConfidenceBand, UQResult


@runtime_checkable
class CalibrationTarget(Protocol):
    """
    Duck-typed against `biosim_validation.benchmark.BenchmarkPoint` --
    deliberately NOT importing biosim_validation here (biosim_uncertainty
    stays below biosim_validation in the dependency graph; any
    `BenchmarkPoint` instance already satisfies this Protocol structurally,
    no adapter needed).
    """
    name: str
    reference: float
    tolerance_fraction: float
    is_scored: bool


# (times_s, {species_name: single-sample trajectory}) -> scalar model value,
# the per-sample analogue of biosim_validation.validator.ModelValidator.
# _extract_model_values() for one Monte Carlo draw instead of one
# SimulationOutput.
MetricExtractor = Callable[[np.ndarray, dict[str, np.ndarray]], float]


@dataclass(frozen=True)
class CalibrationResult:
    """Posterior of a Bayesian calibration run. See module docstring for caveats."""
    times_s:                        np.ndarray
    posterior_bands:                dict[str, ConfidenceBand]
    posterior_parameter_samples:    dict[str, np.ndarray]
    effective_sample_size:          float
    effective_sample_size_fraction: float
    is_degenerate:                  bool
    likelihood_sigma_used:          dict[str, float]
    prior_overall_confidence:       float
    n_prior_samples:                int
    n_posterior_samples:            int
    targets_used:                   list[str]
    notes:                          list[str] = field(default_factory=list)

    def to_uq_result(self) -> UQResult:
        """
        Wrap the posterior as a UQResult so existing UQResult consumers
        (lay reports, etc.) need no changes to accept a calibrated
        posterior instead of a raw forward-propagated prior. Note:
        `overall_confidence` here is the PRIOR run's average watchdog
        confidence over all samples, not recomputed specifically over
        the (reweighted) posterior -- MonteCarloEngine.propagate() only
        returns that average, not a per-sample confidence array, so a
        posterior-specific recomputation isn't available without further
        plumbing. Documented rather than silently presented as exact.
        """
        return UQResult(
            times_s=self.times_s,
            bands=self.posterior_bands,
            n_samples=self.n_posterior_samples,
            n_failed=0,
            overall_confidence=self.prior_overall_confidence,
            parameter_names=list(self.posterior_parameter_samples.keys()),
            raw_trajectories=None,
            sampled_parameter_values=self.posterior_parameter_samples,
            failure_reasons=[],
        )


class BayesianCalibrator:
    """
    Calibrates uncertain PathwayParameters against clinical benchmark
    reference values via Sampling-Importance-Resampling. See module
    docstring for the method, likelihood, and honesty caveats.

    Usage:
        from biosim_uncertainty import BayesianCalibrator
        from biosim_uncertainty.distributions import LogNormalDistribution
        from biosim_validation.benchmark import ogtt_75g_healthy_adults

        calibrator = BayesianCalibrator(random_seed=42)
        distributions = {
            "p1": LogNormalDistribution.from_mean_cv(0.0317, cv=0.30),
            "p3": LogNormalDistribution.from_mean_cv(9.032e-7, cv=0.30),
        }
        dataset = ogtt_75g_healthy_adults()
        result = calibrator.calibrate(
            distributions=distributions,
            calibration_targets=dataset.points,
            metric_extractors=CALIBRATION_METRIC_EXTRACTORS,
            meals=[high_sugar_meal()],
        )
        print(result.effective_sample_size_fraction, result.is_degenerate)
    """

    def __init__(
        self,
        base_params=None,
        random_seed: int = 42,
        likelihood_sigma_z: float = 1.0,
        min_effective_sample_size_fraction: float = 0.10,
    ) -> None:
        """
        Args:
            base_params: Base PathwayParameters. Distributions override
                individual fields (same convention as MonteCarloEngine).
            random_seed: For reproducibility of both the prior draw and
                the posterior resampling step.
            likelihood_sigma_z: How many implied "sigma" the tolerance
                band's edge represents (see module docstring) -- z=1.0
                (the default) means the tolerance-band edge sits at one
                standard deviation, a stated, not measured, choice.
            min_effective_sample_size_fraction: Below this fraction of
                the SUCCESSFUL prior samples (`prior.n_successful`, not
                the requested `n_prior_samples` -- they differ whenever
                any prior draws fail to solve), `CalibrationResult.
                is_degenerate` is True -- the posterior should be
                treated as unreliable (too few effectively-distinct
                samples survived reweighting) rather than trusted at
                face value.
        """
        if likelihood_sigma_z <= 0:
            raise ValueError(f"likelihood_sigma_z must be positive, got {likelihood_sigma_z}")
        if not 0.0 < min_effective_sample_size_fraction <= 1.0:
            raise ValueError(
                "min_effective_sample_size_fraction must be in (0, 1], "
                f"got {min_effective_sample_size_fraction}"
            )
        self.base_params = base_params
        self.random_seed = random_seed
        self.likelihood_sigma_z = likelihood_sigma_z
        self.min_effective_sample_size_fraction = min_effective_sample_size_fraction

    def calibrate(
        self,
        distributions: dict[str, ParameterDistribution],
        calibration_targets: list,
        metric_extractors: dict[str, MetricExtractor],
        meals: Optional[list[MealInput]] = None,
        duration_s: float = 14400.0,
        n_output_points: int = 300,
        n_prior_samples: int = 1000,
        n_posterior_samples: Optional[int] = None,
        pathway_factory: Optional[PathwayFactory] = None,
    ) -> CalibrationResult:
        """
        Run the prior Monte Carlo propagation, compute each successful
        sample's likelihood against `calibration_targets`, and resample
        by posterior weight.

        Args:
            distributions: Prior distributions, same shape as
                MonteCarloEngine.propagate()'s own argument. Not
                restricted to LogNormalDistribution here (any
                ParameterDistribution is accepted), but most
                PathwayParameters fields have no positivity check of
                their own (see pathway.py's __post_init__ -- only Gb,
                Ib, ke, p1, p2, EC50_mTOR, EC50_AMPK, and
                secretion_scale are validated). Per this project's own
                rule 4 ("all kinetic parameters are always-positive,
                use LogNormalDistribution"), a caller supplying e.g. a
                NormalDistribution for an unvalidated field could
                silently sample a physically nonsensical (negative)
                value that still constructs successfully -- use
                LogNormalDistribution for any field without its own
                positivity check, the same convention the rest of this
                codebase follows.
            calibration_targets: Objects satisfying CalibrationTarget
                (e.g. a BenchmarkDataset's .points list). Targets with
                is_scored=False are skipped entirely (same treatment
                ModelValidator gives them), and a target with no matching
                entry in metric_extractors is also skipped (not an
                error -- lets a caller pass a full benchmark dataset and
                only supply extractors for the subset they want to
                calibrate against).
            metric_extractors: Maps target.name -> a function computing
                that target's model value from one sample's (times_s,
                species trajectories).
            n_prior_samples: How many prior draws to propagate. Every
                sample this expensive step succeeds at is reused for
                free by SIR -- no extra ODE solves in the resampling step.
            n_posterior_samples: How many posterior draws to resample.
                Defaults to n_prior_samples.
            pathway_factory: Optional (params, meals) -> BasePathway
                callable, forwarded to MonteCarloEngine.propagate()
                unchanged. Defaults to None (unchanged prior, bare-GI
                behaviour) -- `calibrate_ogtt()` deliberately keeps
                passing nothing, since a healthy 75g OGTT activates zero
                extension pathways; this exists for a future comparative
                -benchmark calibration (e.g. semaglutide) that would need
                the composed pathway UQ propagation is otherwise blind to.

        Returns:
            CalibrationResult with posterior bands, resampled parameter
            draws, and the honesty/degeneracy diagnostics described in
            the module docstring.
        """
        n_posterior_samples = n_posterior_samples or n_prior_samples

        engine = MonteCarloEngine(
            n_samples=n_prior_samples,
            random_seed=self.random_seed,
            base_params=self.base_params,
        )
        prior = engine.propagate(
            distributions=distributions,
            meals=meals,
            duration_s=duration_s,
            n_output_points=n_output_points,
            store_trajectories=True,
            pathway_factory=pathway_factory,
        )

        if prior.n_successful == 0:
            raise ValueError(
                "BayesianCalibrator: every prior sample failed to solve -- "
                f"cannot calibrate. Failure reasons: {prior.failure_reasons[:5]}"
            )

        scored_targets = [
            t for t in calibration_targets
            if t.is_scored and t.name in metric_extractors
        ]
        if not scored_targets:
            raise ValueError(
                "BayesianCalibrator: no scored calibration_targets have a "
                "matching entry in metric_extractors -- nothing to "
                "calibrate against."
            )

        n_successful = prior.n_successful
        log_weights = np.zeros(n_successful)
        likelihood_sigma_used: dict[str, float] = {}
        notes: list[str] = [
            "Likelihood treats each target's tolerance_fraction as an "
            f"implied Gaussian sigma (likelihood_sigma_z={self.likelihood_sigma_z}) -- "
            "a stated modeling choice representing this project's policy "
            "tolerance bands, not a measured per-study confidence interval. "
            "See BayesianCalibrator's module docstring."
        ]

        for target in scored_targets:
            sigma = target.reference * target.tolerance_fraction / self.likelihood_sigma_z
            if sigma <= 0:
                raise ValueError(
                    f"Calibration target '{target.name}' has non-positive "
                    f"implied sigma ({sigma}) -- reference and "
                    "tolerance_fraction must both be positive."
                )
            likelihood_sigma_used[target.name] = sigma
            extractor = metric_extractors[target.name]
            log_const = -math.log(sigma) - 0.5 * math.log(2.0 * math.pi)

            for i in range(n_successful):
                sample_species = {
                    species: traj[i] for species, traj in prior.raw_trajectories.items()
                }
                try:
                    model_value = extractor(prior.times_s, sample_species)
                    if not math.isfinite(model_value):
                        raise ValueError(
                            f"non-finite model value ({model_value!r})"
                        )
                    log_lik = -0.5 * ((model_value - target.reference) / sigma) ** 2 + log_const
                except Exception as e:
                    # A single sample's metric extraction failing (e.g. a
                    # pathological trajectory shape) or returning a
                    # non-finite value (e.g. NaN from a numerically
                    # unstable but "successful" ODE solve -- np.max/
                    # np.argmax/np.interp propagate NaN silently, they
                    # don't raise, so this needed an explicit check, not
                    # just a try/except) must not crash the whole
                    # calibration, and critically must not poison every
                    # OTHER sample's weight too: logsumexp(log_weights)
                    # is NaN if even ONE entry is NaN, and subtracting a
                    # NaN normalizer turns the ENTIRE weights array into
                    # NaN (confirmed directly: a single NaN log-weight
                    # among otherwise-finite ones propagates to all of
                    # them) -- found during self-verification while the
                    # mandatory review-agent pass was still running, see
                    # CLAUDE.md. Exclude just that sample's contribution
                    # for this target (treat as "no information from
                    # this target" for that one sample, not "zero
                    # weight" for the whole run).
                    notes.append(
                        f"metric_extractors['{target.name}'] raised "
                        f"{type(e).__name__} on sample {i}; excluded from "
                        "that target's likelihood contribution only."
                    )
                    log_lik = 0.0
                log_weights[i] += log_lik

        log_weights -= logsumexp(log_weights)
        weights = np.exp(log_weights)

        effective_sample_size = float(1.0 / np.sum(weights ** 2))
        ess_fraction = effective_sample_size / n_successful
        is_degenerate = ess_fraction < self.min_effective_sample_size_fraction
        if is_degenerate:
            notes.append(
                f"Effective sample size fraction ({ess_fraction:.3f}) is below "
                f"min_effective_sample_size_fraction "
                f"({self.min_effective_sample_size_fraction}) -- the "
                "likelihood is much narrower than the prior over these "
                "targets, and this posterior is a handful of duplicated "
                "points, not a reliable distribution. Treat as provisional."
            )

        # Separate RNG from the prior's own (which is seeded and consumed
        # entirely inside MonteCarloEngine.propagate()) -- an offset seed
        # keeps the resampling step reproducible without any appearance
        # of correlation with the prior draw's own random sequence.
        resample_rng = np.random.default_rng(self.random_seed + 1)
        idx = resample_rng.choice(n_successful, size=n_posterior_samples, replace=True, p=weights)

        posterior_bands: dict[str, ConfidenceBand] = {}
        for species, traj in prior.raw_trajectories.items():
            resampled = traj[idx, :]
            posterior_bands[species] = ConfidenceBand(
                lower=np.percentile(resampled, 2.5, axis=0),
                median=np.percentile(resampled, 50.0, axis=0),
                upper=np.percentile(resampled, 97.5, axis=0),
            )

        posterior_parameter_samples = {
            name: values[idx] for name, values in prior.sampled_parameter_values.items()
        }

        return CalibrationResult(
            times_s=prior.times_s,
            posterior_bands=posterior_bands,
            posterior_parameter_samples=posterior_parameter_samples,
            effective_sample_size=effective_sample_size,
            effective_sample_size_fraction=ess_fraction,
            is_degenerate=is_degenerate,
            likelihood_sigma_used=likelihood_sigma_used,
            prior_overall_confidence=prior.overall_confidence,
            n_prior_samples=n_prior_samples,
            n_posterior_samples=n_posterior_samples,
            targets_used=[t.name for t in scored_targets],
            notes=notes,
        )


# ---------------------------------------------------------------------------
# Standard metric extractors -- generic ODE-trajectory metrics matching the
# most common biosim_validation BenchmarkPoint names (species trajectories
# are in canonical units: glucose/insulin/etc. in uM, time in seconds; see
# ModelValidator._extract_model_values() for the equivalent computation
# against a full SimulationOutput rather than one Monte Carlo draw).
# ---------------------------------------------------------------------------

def _fasting_glucose(times_s: np.ndarray, species: dict[str, np.ndarray]) -> float:
    return float(species["G"][0]) / 1000.0   # uM -> mM


def _peak_glucose(times_s: np.ndarray, species: dict[str, np.ndarray]) -> float:
    return float(np.max(species["G"])) / 1000.0


def _time_to_peak_glucose_min(times_s: np.ndarray, species: dict[str, np.ndarray]) -> float:
    return float(times_s[int(np.argmax(species["G"]))]) / 60.0


def _glucose_at_120min(times_s: np.ndarray, species: dict[str, np.ndarray]) -> float:
    return float(np.interp(7200.0, times_s, species["G"])) / 1000.0


def _peak_insulin(times_s: np.ndarray, species: dict[str, np.ndarray]) -> float:
    return float(np.max(species["I"])) * 1e6   # uM -> pM


def _time_to_peak_insulin_min(times_s: np.ndarray, species: dict[str, np.ndarray]) -> float:
    return float(times_s[int(np.argmax(species["I"]))]) / 60.0


def _glucose_return_time_min(times_s: np.ndarray, species: dict[str, np.ndarray]) -> float:
    """
    Time for glucose to return within 10% of its own baseline (first
    value) after its peak -- matches
    biosim_api.output.SimulationOutput.time_to_return_to_baseline()'s
    "G", threshold_fraction=0.10 convention. Returns the trajectory's own
    final time if it never returns (a large but finite value, not None,
    since MetricExtractor must return a plain float for the Gaussian
    likelihood -- a target using this extractor should have a WIDE
    tolerance precisely because "never returned within the window" is a
    real, if crude, signal that the sample is a poor match, not an error).
    """
    G = species["G"]
    baseline = G[0]
    threshold = baseline * 1.10
    peak_idx = int(np.argmax(G))
    for i in range(peak_idx, len(G)):
        if G[i] <= threshold:
            return float(times_s[i]) / 60.0
    return float(times_s[-1]) / 60.0


CALIBRATION_METRIC_EXTRACTORS: dict[str, MetricExtractor] = {
    "fasting_glucose": _fasting_glucose,
    "peak_glucose": _peak_glucose,
    "time_to_peak_glucose_min": _time_to_peak_glucose_min,
    "glucose_at_120min": _glucose_at_120min,
    "peak_insulin": _peak_insulin,
    "time_to_peak_insulin_min": _time_to_peak_insulin_min,
    "glucose_return_time_min": _glucose_return_time_min,
}
