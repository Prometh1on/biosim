# SAVE TO: packages\uncertainty\src\biosim_uncertainty\propagation.py
"""
Monte Carlo uncertainty propagation engine.

Runs N simulations with parameters sampled from their distributions,
then computes confidence bands from the resulting ensemble.

Why Monte Carlo (not analytical propagation):
    - The ODE system is nonlinear — analytical error propagation
      (linearisation) is only valid for small uncertainties
    - Monte Carlo makes no assumptions about the shape of the
      output distribution and handles nonlinear parameter interactions
    - N=200 is sufficient for 95% CI estimation; N=1000 for publishing

Performance:
    Each ODE solve takes ~5-20 ms. N=200 → ~2-4 seconds.
    Simulations are independent so parallelisation is trivial (Sprint 14).
    For now, runs sequentially with a progress callback option.
"""

from __future__ import annotations

import dataclasses
import inspect
from typing import Callable, Mapping, Optional
import numpy as np

from biosim_ode import ODESystem, GlucoseInsulinPathway, PathwayParameters
from biosim_ode.base_pathway import BasePathway
from biosim_ode.meal import MealInput

from .distributions import ParameterDistribution, FixedValue
from .result import ConfidenceBand, UQResult


# A pathway factory lets a caller (e.g. SimulationRunner) supply the SAME
# composed CompositePathway used for the point estimate -- without one,
# every UQ/sensitivity entry point below builds a bare GlucoseInsulinPathway
# regardless of what extension pathways (GLP1, renal/SGLT2i, cortisol,
# hepatic, exercise, meal-absorption) are actually active for the scenario,
# so the reported confidence band can silently describe a structurally
# different ODE system than the one that produced the point estimate.
# Verified directly: for a semaglutide/dapagliflozin scenario (both have
# ZERO static PathwayParameters effect -- _glp1_agonist_modifiers()/
# _sglt2_inhibitor_modifiers() both return {}), the pre-fix band was
# bit-for-bit identical to a drug-free ensemble's band, since the sampled
# PathwayParameters were identical either way.
#
# (params, meals) is the minimum sufficient signature: every extension
# pathway in SimulationRunner._build_pathway() is constructed from the
# scenario alone, never from PathwayParameters, so only the GI pathway
# needs the sampled params. meals must be a factory ARGUMENT, not closed
# over by the caller, because MealAbsorptionPathway derives its own y0
# from meals[i].D_eff_mM and implements no with_meals() -- a factory
# closing over stale meals would silently produce wrong absorption
# initial conditions if meals ever differ from what's passed in here
# (they don't today, but the contract should not rely on that).
PathwayFactory = Callable[[PathwayParameters, list[MealInput]], BasePathway]

# MASTER_PLAN 2.2: pathway_factory (above) closed the "wrong pathway
# STRUCTURE" gap -- the UQ ensemble now sees the same composed extension
# pathways (GLP1, renal/SGLT2i, cortisol, hepatic, ...) the point estimate
# does. A narrower gap remained: extension-pathway-specific parameters
# (e.g. RenalClearanceParameters.beta) were still FIXED at their point-
# estimate value inside every sample, even once the STRUCTURE was right --
# so a drug scenario's confidence band correctly reflected uncertainty in
# p1/p2/p3/ke/EC50_mTOR but silently treated the drug's own calibrated
# effect size as if it were exactly known.
#
# This is deliberately a SEPARATE, additive parameter pair
# (extension_distributions/extension_pathway_factory), not a change to
# PathwayFactory's own 2-argument contract -- extension-specific
# parameters (agonist_scale, beta, ...) do not live on PathwayParameters
# at all (it's a frozen dataclass with a fixed field set); they live on
# separate parameter objects (GLP1Parameters, RenalClearanceParameters)
# that a pathway constructs from at __init__ time. Every existing
# 2-argument pathway_factory (SimulationRunner's real one, calibration.py,
# every pre-existing test) is completely unaffected -- this parameter
# pair defaults to None/None and is validated to not silently no-op (see
# propagate()'s own checks) if partially supplied.
#
# Extension parameter keys are dotted ("renal.beta", "glp1.agonist_scale")
# -- mirrors the "<pathway_name>.<species>" convention _species_pathway_map()
# already uses, which makes collision with a core PathwayParameters field
# name structurally impossible (no core field name contains ".").
ExtensionPathwayFactory = Callable[
    [PathwayParameters, list[MealInput], Mapping[str, float]], BasePathway
]


def _species_pathway_map(pathway: BasePathway) -> dict[str, str]:
    """
    Map each of a pathway's own bare species names to the name of the
    sub-pathway that produced it -- mirrors ODESystem.simulate()'s own
    "<pathway_name>.<species>" prefix-stripping exactly (system.py), so
    UQResult.species_pathway's keys stay consistent with
    SimulationResult.species's own bare keys.

    For a bare (non-composite) pathway, species_names carries no "."
    prefix at all -- every species is attributed to the pathway's own
    .name in that case (matching CompositePathway.species_names' own
    "<pathway.name>.<species>" convention for the single-pathway case).

    Computed once per propagate() call, from the first successful
    sample's own pathway -- safe for every pathway_factory in current
    production use (SimulationRunner's factory depends only on
    `scenario`, fixed for the whole call, so pathway structure cannot
    vary across samples). Not re-validated against later samples: a
    custom pathway_factory that swaps which extension composes based on
    the SAMPLED PARAMETERS (still "deterministic given (params, meals)"
    per propagate()'s own contract) could, in principle, leave this
    mapping attributing a species to a stale pathway if two candidate
    extensions shared a bare species name -- against this platform's own
    global species-name-uniqueness invariant, but not mechanically
    enforced here. species_pathway is presentational only (no
    calculation reads it), so this is a documented limitation, not a
    correctness risk for anything built on this result today.
    """
    mapping: dict[str, str] = {}
    for full_name in pathway.species_names:
        if "." in full_name:
            pathway_name, bare_name = full_name.rsplit(".", 1)
        else:
            pathway_name, bare_name = pathway.name, full_name
        mapping[bare_name] = pathway_name
    return mapping


class MonteCarloEngine:
    """
    Propagates parameter uncertainty through the ODE system via Monte Carlo.

    Usage:
        from biosim_uncertainty import MonteCarloEngine
        from biosim_uncertainty.distributions import LogNormalDistribution

        engine = MonteCarloEngine(n_samples=200, random_seed=42)

        distributions = {
            "p1": LogNormalDistribution.from_mean_cv(0.0317, cv=0.30),
            "ke": LogNormalDistribution.from_mean_cv(0.093,  cv=0.30),
            "EC50_mTOR": LogNormalDistribution.from_mean_cv(100.0, cv=0.50),
        }

        result = engine.propagate(
            distributions=distributions,
            meals=[high_sugar_meal()],
            duration_s=14400,
        )

        glucose_band = result.band("G")
        print(glucose_band.median)   # 50th percentile glucose over time
    """

    def __init__(
        self,
        n_samples: int = 200,
        random_seed: int = 42,
        base_params: Optional[PathwayParameters] = None,
        progress_callback: Optional[Callable[[int, int], None]] = None,
    ) -> None:
        """
        Args:
            n_samples: Number of Monte Carlo samples (200 for exploration, 1000 for publication)
            random_seed: For reproducibility. Every result with the same seed is identical.
            base_params: Base pathway parameters. Distributions override individual fields.
            progress_callback: Called with (current, total) at each sample. For progress bars.
        """
        if n_samples < 10:
            raise ValueError(f"n_samples must be ≥ 10 for meaningful CI, got {n_samples}")
        self.n_samples = n_samples
        self.random_seed = random_seed
        self.base_params = base_params or PathwayParameters()
        self.progress_callback = progress_callback

    def propagate(
        self,
        distributions: dict[str, ParameterDistribution],
        meals: Optional[list[MealInput]] = None,
        duration_s: float = 14400.0,
        n_output_points: int = 300,
        store_trajectories: bool = False,
        pathway_factory: Optional[PathwayFactory] = None,
        extension_distributions: Optional[dict[str, ParameterDistribution]] = None,
        extension_pathway_factory: Optional[ExtensionPathwayFactory] = None,
    ) -> UQResult:
        """
        Run Monte Carlo propagation and return confidence bands.

        Args:
            distributions: Dict mapping PathwayParameters field names to distributions.
                           Fields not listed use their base_params value (FixedValue).
            meals: Meal inputs for the simulation
            duration_s: Simulation duration in seconds
            n_output_points: Number of time points in output
            store_trajectories: If True, store all N trajectories in the result
                               (uses more memory, useful for detailed analysis)
            pathway_factory: Optional (params, meals) -> BasePathway callable,
                used in place of a bare GlucoseInsulinPathway for every
                sample. Defaults to None, which reproduces the EXACT prior
                expression (bit-for-bit identical output for the same seed)
                -- every existing caller is unaffected. Pass a factory that
                reuses SimulationRunner._build_pathway() to make the UQ
                ensemble see the same composed extension pathways (GLP1,
                renal, cortisol, hepatic, ...) the point estimate does; see
                this module's PathwayFactory docstring for why this
                matters and why (params, meals) is the right signature.
                Must be deterministic given (params, meals), or
                random_seed reproducibility is silently voided.
                Mutually exclusive with extension_pathway_factory (below)
                -- pass at most one of the two.
            extension_distributions: Optional dict mapping DOTTED extension-
                parameter keys (e.g. "renal.beta") to distributions, sampled
                alongside `distributions` and passed to
                extension_pathway_factory as its third argument. Requires
                extension_pathway_factory to also be supplied (raises
                ValueError otherwise -- a distribution with nowhere to go
                would be sampled, recorded, and silently never reach any
                pathway, the exact "computed but discarded" bug class this
                project has repeatedly found and fixed elsewhere). Keys
                must not collide with `distributions`' own field names
                (raises ValueError) -- dotted keys make this structurally
                unlikely in practice, but both dicts are merged into
                UQResult.sampled_parameter_values/parameter_names, so a
                collision would silently overwrite one draw's record with
                the other's.
            extension_pathway_factory: Optional (params, meals,
                extension_overrides) -> BasePathway callable -- like
                pathway_factory, but ALSO receives this sample's own
                extension-parameter draws (a dict of dotted key -> sampled
                value, built from extension_distributions; empty dict if
                extension_distributions is None/empty). Mutually exclusive
                with pathway_factory -- pass at most one of the two. The
                factory is responsible for routing each dotted key to the
                right extension parameter object at CONSTRUCTION time (see
                SimulationRunner._build_pathway()'s extension_overrides
                parameter) -- never by mutating an already-constructed
                pathway's parameters in place, which would leave any
                y0/derived quantity computed from the OLD value stale (the
                same fixed-point-drift bug class as KNOWN BUGS FIXED #10-14
                in CLAUDE.md).

        Returns:
            UQResult with 95% confidence bands for all species
        """
        if pathway_factory is not None and extension_pathway_factory is not None:
            raise ValueError(
                "propagate() accepts at most one of pathway_factory/"
                "extension_pathway_factory, not both -- there is no defined "
                "precedence between them."
            )
        extension_distributions = extension_distributions or {}
        if extension_distributions and extension_pathway_factory is None:
            raise ValueError(
                "extension_distributions was supplied but extension_pathway_factory "
                "was not -- the sampled extension values would have nowhere to go "
                "and be silently discarded. Pass extension_pathway_factory too."
            )
        colliding_keys = set(distributions) & set(extension_distributions)
        if colliding_keys:
            raise ValueError(
                f"extension_distributions key(s) {sorted(colliding_keys)} collide "
                f"with a core distributions field name -- use a dotted extension "
                f"key (e.g. 'renal.beta') to avoid ambiguity in "
                f"UQResult.sampled_parameter_values/parameter_names."
            )
        if extension_pathway_factory is not None:
            # Fail once, at call entry, on a wrong-arity factory -- without
            # this, a caller-supplied factory that doesn't accept 3
            # arguments raises the same TypeError on every single sample
            # (n_samples times), burning the whole sample budget on a
            # single, individually-uninformative-looking failure message
            # per sample instead of one clear diagnostic up front. A dry
            # bind (not a real call) correctly handles *args/**kwargs/
            # defaults the same way an actual call would, without
            # hand-rolling arity-counting logic.
            try:
                inspect.signature(extension_pathway_factory).bind(None, None, {})
            except TypeError as e:
                raise ValueError(
                    f"extension_pathway_factory must accept 3 positional arguments "
                    f"(params, meals, extension_overrides) -- got a signature "
                    f"incompatible with that call: {e}"
                ) from e

        rng = np.random.default_rng(self.random_seed)
        meals = meals or []

        # Run N simulations
        all_trajectories: dict[str, list[np.ndarray]] = {}
        sampled_param_values: dict[str, list[float]] = {
            name: [] for name in list(distributions) + list(extension_distributions)
        }
        confidences: list[float] = []
        times_s: Optional[np.ndarray] = None
        n_failed = 0
        failure_reasons: list[str] = []
        species_pathway: dict[str, str] = {}

        for i in range(self.n_samples):
            if self.progress_callback:
                self.progress_callback(i + 1, self.n_samples)

            # Sample a parameter set
            sampled_params = self._sample_params(distributions, rng)
            # Deliberately NOT floored at positivity like the core sample
            # above -- an out-of-range extension draw should reach the
            # extension parameter object's own constructor validation
            # (e.g. RenalClearanceParameters rejects beta outside [0, 1))
            # and be counted as a real, reasoned failure in
            # failure_reasons, not silently mapped to "the drug does
            # nothing" the way the core floor does for a kinetic rate
            # constant.
            sampled_extension_values = self._sample_extension_values(extension_distributions, rng)

            # Build pathway with sampled parameters
            try:
                if extension_pathway_factory is not None:
                    pathway = extension_pathway_factory(
                        sampled_params, meals, sampled_extension_values
                    )
                elif pathway_factory is not None:
                    pathway = pathway_factory(sampled_params, meals)
                else:
                    pathway = GlucoseInsulinPathway(params=sampled_params, meals=meals)
                system = ODESystem(pathway=pathway)
                result = system.simulate(
                    duration_s=duration_s,
                    meals=meals,
                    n_output_points=n_output_points,
                    check_clamps=False,
                )
            except Exception as e:
                n_failed += 1
                failure_reasons.append(f"sample {i}: {type(e).__name__}: {e}")
                continue

            if not result.success:
                n_failed += 1
                failure_reasons.append(f"sample {i}: ODE solver failed: {result.message}")
                continue

            # Store time array (same for all — we use same t_eval)
            if times_s is None:
                times_s = result.times_s
                for name in result.species:
                    all_trajectories[name] = []
                species_pathway = _species_pathway_map(pathway)
            elif set(result.species) != set(all_trajectories):
                # A caller-supplied pathway_factory must be deterministic
                # given (params, meals) -- if it isn't (e.g. it varies
                # which extension pathways compose based on something
                # other than params/meals), later samples silently raise
                # a bare KeyError below with no indication why. Fail with
                # a clear, named cause instead.
                raise ValueError(
                    f"sample {i}'s pathway produced a different species set "
                    f"than the first successful sample: "
                    f"{sorted(set(result.species) ^ set(all_trajectories))} differ. "
                    f"pathway_factory must be deterministic given (params, meals)."
                )

            # Store species trajectories
            for name, values in result.species.items():
                all_trajectories[name].append(values)

            # Store the parameter draw that produced this trajectory --
            # BayesianCalibrator (calibration.py) needs both together to
            # resample parameter draws by posterior weight (Sampling-
            # Importance-Resampling). Only kept for successful samples,
            # same lifetime/ordering as all_trajectories above.
            for name in distributions:
                sampled_param_values[name].append(getattr(sampled_params, name))
            for name, value in sampled_extension_values.items():
                sampled_param_values[name].append(value)

            # Track watchdog confidence
            confidences.append(system.overall_confidence(result))

        if times_s is None or not all_trajectories:
            # All simulations failed
            return UQResult(
                times_s=np.array([]),
                bands={},
                n_samples=self.n_samples,
                n_failed=n_failed,
                overall_confidence=0.0,
                parameter_names=list(distributions.keys()) + list(extension_distributions.keys()),
                failure_reasons=failure_reasons,
            )

        # Compute confidence bands
        bands = {}
        raw_traj: dict[str, np.ndarray] = {}

        for name, traj_list in all_trajectories.items():
            traj = np.stack(traj_list, axis=0)  # shape (n_successful, n_times)
            raw_traj[name] = traj
            bands[name] = ConfidenceBand(
                lower=np.percentile(traj, 2.5, axis=0),
                median=np.percentile(traj, 50.0, axis=0),
                upper=np.percentile(traj, 97.5, axis=0),
            )

        overall_conf = float(np.mean(confidences)) if confidences else 0.0
        sampled_param_arrays = {
            name: np.array(values) for name, values in sampled_param_values.items()
        }

        return UQResult(
            times_s=times_s,
            bands=bands,
            n_samples=self.n_samples,
            n_failed=n_failed,
            overall_confidence=overall_conf,
            parameter_names=list(distributions.keys()) + list(extension_distributions.keys()),
            raw_trajectories=raw_traj if store_trajectories else None,
            sampled_parameter_values=sampled_param_arrays if store_trajectories else None,
            failure_reasons=failure_reasons,
            species_pathway=species_pathway,
        )

    def _sample_params(
        self,
        distributions: dict[str, ParameterDistribution],
        rng: np.random.Generator,
    ) -> PathwayParameters:
        """Sample one parameter set from the distributions."""
        overrides = {}
        for field_name, dist in distributions.items():
            sampled = dist.sample(rng)
            # Enforce positivity for kinetic parameters
            if sampled <= 0:
                sampled = dist.mean() * 1e-6  # floor at effectively zero
            overrides[field_name] = sampled
        return dataclasses.replace(self.base_params, **overrides)

    def _sample_extension_values(
        self,
        extension_distributions: dict[str, ParameterDistribution],
        rng: np.random.Generator,
    ) -> dict[str, float]:
        """
        Sample one extension-parameter draw per key, from the same rng
        stream as _sample_params() (called immediately before this, per
        sample, in propagate()'s loop) -- deliberately NOT applying
        _sample_params()'s own positivity floor. That floor exists to keep
        a bad kinetic-rate draw from crashing the ODE solver by silently
        mapping it to "effectively zero"; for an extension parameter, the
        analogous move would silently reinterpret an out-of-range draw as
        "the drug has no effect", which is not the same thing and would
        hide a genuinely miscalibrated distribution rather than surface
        it. An out-of-range draw is instead left to reach the extension
        parameter object's own constructor validation (e.g.
        RenalClearanceParameters rejects beta outside [0, 1)), which
        propagate()'s existing try/except turns into a real, reasoned
        entry in failure_reasons -- the honest outcome, whose fix is
        picking a properly bounded distribution, not silently truncating.
        """
        return {name: dist.sample(rng) for name, dist in extension_distributions.items()}
