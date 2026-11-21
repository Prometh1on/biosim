# SAVE TO: packages\data-seeds\src\biosim_seeds\session.py
"""
PlatformSession -- single developer entry point for the biosim platform.

Wires together: seeded IDRegistry + ParameterStore, literature parameter
loading, caching, simulation, reporting, SBML export, and validation.

Instead of constructing and connecting 8 packages manually, use:

    session = PlatformSession()
    result  = session.simulate("high_sugar_meal")
    report  = session.report(result)
    sbml    = session.export_sbml()

The session is stateful (holds the seeded stores and cache) but
each simulate() call is independent. Safe to reuse across many runs.
"""

from __future__ import annotations

from typing import Optional

from biosim_ode import GlucoseInsulinPathway, PathwayParameters
from biosim_infra import CachedRunner, PrivacyAuditor, create_bundle

from .loader import ParameterLoader
from .registry_seed import create_seeded_registry
from .parameter_seed import create_seeded_store


class PlatformSession:
    """
    Fully-wired biosim session.

    Args:
        use_literature_params: If True (default), loads ODE parameters from
                               the seeded ParameterStore. If False, uses
                               hardcoded PathwayParameters defaults.
        cache_size:            Number of simulation results to cache.
    """

    def __init__(
        self,
        use_literature_params: bool = True,
        cache_size: int = 500,
    ) -> None:
        from biosim_curation import CurationQueue

        self._registry        = create_seeded_registry()
        self._store           = create_seeded_store(self._registry)
        self._loader          = ParameterLoader(self._store, self._registry)
        self._curation_queue  = CurationQueue()
        self._base_params   = (
            self._loader.load_pathway_params()
            if use_literature_params
            else PathwayParameters()
        )
        # Regression: this used to bypass CachedRunner.__init__() entirely
        # via CachedRunner.__new__(CachedRunner) + a manual _cache
        # assignment -- unnecessary fragility, since CachedRunner.__init__
        # already accepts an optional `cache` param for exactly this.
        # base_params threaded through so simulate() (below) can route
        # through this cache without losing self._base_params -- same fix
        # as SimulationRunner.__init__, applied to CachedRunner too.
        from biosim_infra import SimulationCache
        self._cached_runner = CachedRunner(
            cache=SimulationCache(max_size=cache_size),
            base_params=self._base_params,
        )

    # ------------------------------------------------------------------
    # Core simulation
    # ------------------------------------------------------------------

    def simulate(
        self,
        scenario,
        profile=None,
        variants: Optional[list] = None,
        run_uncertainty: bool = False,
        n_mc_samples: int = 100,
    ):
        """
        Run a simulation from a scenario name or ScenarioDefinition.

        Args:
            scenario:        Scenario name (str) from the library, or a
                             ScenarioDefinition instance.
            profile:         Optional PhysiologicalProfile for personalisation.
            variants:        Optional list of Variant objects for genetic personalisation.
            run_uncertainty: Run Monte Carlo UQ alongside the main simulation.
            n_mc_samples:    Number of MC samples if run_uncertainty=True.

        Returns:
            SimulationOutput
        """
        import dataclasses
        from biosim_api import get_scenario, ScenarioDefinition, PersonalizationSpec

        if isinstance(scenario, str):
            scenario = get_scenario(scenario)

        # Inject personalisation
        spec = PersonalizationSpec(
            physiological_profile=profile,
            genetic_variants=variants or [],
        )
        scenario = dataclasses.replace(
            scenario,
            personalization=spec,
            run_uncertainty=run_uncertainty,
            n_mc_samples=n_mc_samples,
        )

        # Regression: this method used to branch on profile/variants -- the
        # no-profile branch built a GlucoseInsulinPathway/ODESystem from
        # self._base_params (the literature-loaded parameters), then
        # discarded both and ran the scenario through a plain
        # SimulationRunner() instead; the personalised branch called
        # SimulationRunner() the same bare way. SimulationRunner itself
        # hardcoded a fresh PathwayParameters() internally, so
        # self._base_params never actually reached the simulation either
        # way (only export_sbml() used it) -- use_literature_params=True
        # (the default) had zero effect on simulate()'s actual output.
        # Invisible only because today's seeded literature values happen to
        # equal the class defaults exactly. Fixed by threading
        # self._base_params through SimulationRunner's own base_params
        # parameter (see runner.py) -- both cases now go through the same
        # call, so the branch (and the pathway/system it built only to
        # discard) is gone entirely.
        #
        # Regression (separate from the above): this used to call a fresh
        # SimulationRunner() directly here, bypassing self._cached_runner
        # entirely -- PlatformSession's advertised caching (cache_size
        # constructor arg, cache_stats()) was real, correctly built, and
        # never actually used. Now routes through the cache instead.
        # Safe to do only as of the SimulationCache._key() collision fix
        # earlier this session -- before that fix, caching a
        # personalization-bearing scenario risked serving another
        # scenario's stale, wrong-personalization result; that's plausibly
        # why this was never wired in the first place.
        #
        # Sync (not just set-once at __init__) so a caller that replaces
        # self._base_params after construction still gets it honoured --
        # self._cached_runner otherwise holds its own independent copy
        # captured at construction time, which would silently drift out
        # of sync with self._base_params if the latter is ever reassigned.
        self._cached_runner._base_params = self._base_params
        return self._cached_runner.run(scenario)

    def compare(self, scenario_a, scenario_b, profile=None) -> object:
        """
        Run two scenarios and return a differential comparison.

        Args:
            scenario_a/b: Scenario names or ScenarioDefinition instances.
            profile:      Optional shared PhysiologicalProfile.

        Returns:
            ComparisonOutput
        """
        from biosim_api import ComparisonRunner, get_scenario, ScenarioDefinition, PersonalizationSpec
        import dataclasses

        def _resolve(s):
            if isinstance(s, str):
                return get_scenario(s)
            return s

        sa = _resolve(scenario_a)
        sb = _resolve(scenario_b)

        if profile is not None:
            spec = PersonalizationSpec(physiological_profile=profile)
            sa = dataclasses.replace(sa, personalization=spec)
            sb = dataclasses.replace(sb, personalization=spec)

        return ComparisonRunner().compare(sa, sb)

    # ------------------------------------------------------------------
    # Output layer
    # ------------------------------------------------------------------

    def report(self, output) -> object:
        """Generate a plain language report from a SimulationOutput."""
        from biosim_lay import generate_report
        return generate_report(output)

    def assess_risk(self, output) -> object:
        """Return a traffic-light risk profile for a SimulationOutput."""
        from biosim_lay import assess_risk
        return assess_risk(output)

    def validate(self) -> object:
        """
        Run the standard OGTT validation against clinical benchmarks.
        Returns a ValidationReport.
        """
        from biosim_validation import ModelValidator
        return ModelValidator().validate_ogtt()

    # ------------------------------------------------------------------
    # Export
    # ------------------------------------------------------------------

    def export_sbml(self, params: Optional[PathwayParameters] = None) -> str:
        """Export the ODE model as SBML XML."""
        from biosim_expert import SBMLExporter
        pathway = GlucoseInsulinPathway(params=params or self._base_params)
        return SBMLExporter().export(pathway)

    def export_csv(self, output) -> str:
        """Export simulation time series as CSV."""
        from biosim_expert import CSVExporter
        return CSVExporter().export_time_series(output)

    def export_metrics_csv(self, outputs: list) -> str:
        """Export a metrics comparison table as CSV."""
        from biosim_expert import CSVExporter
        return CSVExporter().export_metrics_table(outputs)

    def export_curation_queue_csv(self) -> str:
        """
        Export this session's pending curation items as CSV (MASTER_PLAN
        5.2d) -- the minimal consumer completing the producer -> queue ->
        consumer loop 5.2's own scoping identified as missing. Pending
        items already include everything 5.2a (ParameterStore conflicts)
        and 5.2b (sensitivity-driven prioritisation) enqueue.
        """
        from biosim_expert import CSVExporter
        return CSVExporter().export_curation_queue(self._curation_queue.pending())

    def bundle(self, output, random_seed: int = 42) -> object:
        """Create a reproducibility bundle for a SimulationOutput."""
        return create_bundle(output, random_seed=random_seed)

    def audit_privacy(self, output) -> object:
        """Run a privacy audit on a SimulationOutput."""
        return PrivacyAuditor().audit(output)

    # ------------------------------------------------------------------
    # Curation (MASTER_PLAN 5.2b)
    # ------------------------------------------------------------------

    def prioritize_curation_from_sensitivity(
        self,
        parameter_names: Optional[list] = None,
        meals: Optional[list] = None,
        duration_s: float = 14400.0,
        **prioritise_kwargs,
    ) -> object:
        """
        Run OAT sensitivity analysis against this session's own
        base_params and enqueue the highest-impact parameters into this
        session's curation queue -- closing the loop the biosim_uncertainty.
        curation_bridge.prioritise_from_sensitivity() bridge function
        already implemented but that nothing ever called from a real
        entry point (MASTER_PLAN 5.2's own scoping finding; this is the
        first real caller).

        This method OWNS the "run sensitivity analysis, then prioritise"
        pairing; it deliberately does NOT own sensitivity analysis
        itself (OATSensitivityAnalyzer stays a standalone, reusable
        class in biosim_uncertainty) -- callers who already have a
        SensitivityResult from elsewhere should call
        `prioritise_from_sensitivity()` directly instead of through here.

        Args:
            parameter_names: PathwayParameters fields to analyze.
                              Defaults to every parameter
                              curation_bridge.PARAM_BIOLOGY_MAP knows how
                              to map to a real biological entity (analyzing
                              a parameter with no known mapping would be
                              silently skipped by prioritise_from_sensitivity
                              anyway -- see its own `if bio is None: continue`).
            meals:           Meal inputs for the sensitivity simulations.
                              Defaults to OATSensitivityAnalyzer's own
                              default (no meal, fasted-state sensitivity)
                              -- this method wires the two steps together,
                              it doesn't impose a new opinion on which
                              scenario is "representative" for curation
                              purposes. CONCRETE CONSEQUENCE, confirmed by
                              direct execution: under this default, p1/p2/
                              p3/ke (insulin action/clearance -- the
                              parameters that govern the TRANSIENT response
                              to a meal) show EXACTLY zero sensitivity on
                              every standard metric, even at a threshold
                              of 0.01 -- with no meal perturbing the
                              system away from its fixed point, these
                              parameters have nothing to act on, so only
                              the basal/steady-state parameters (Gb, Ib,
                              EC50_mTOR) ever get prioritised. Pass a real
                              meal (e.g. via `biosim_api.get_scenario(...)
                              .meal_events`, converted to `MealInput`) if
                              transient-response parameters should also be
                              considered for curation.
            duration_s:      Simulation duration for each sensitivity run.
            **prioritise_kwargs: Forwarded to prioritise_from_sensitivity()
                              (high_sensitivity_threshold,
                              critical_sensitivity_threshold,
                              top_n_per_metric).

        Returns:
            biosim_uncertainty.curation_bridge.BridgeResult
        """
        from biosim_uncertainty import (
            OATSensitivityAnalyzer,
            PARAM_BIOLOGY_MAP,
            prioritise_from_sensitivity,
        )

        names = parameter_names or list(PARAM_BIOLOGY_MAP.keys())
        result = OATSensitivityAnalyzer(base_params=self._base_params).analyze(
            parameter_names=names, meals=meals, duration_s=duration_s,
        )
        return prioritise_from_sensitivity(
            result, self._curation_queue, self._store, self._registry,
            **prioritise_kwargs,
        )

    # ------------------------------------------------------------------
    # Introspection
    # ------------------------------------------------------------------

    def provenance(self) -> dict:
        """Return a report of which parameters came from literature."""
        return self._loader.provenance_report()

    def cache_stats(self) -> dict:
        """Return simulation cache hit/miss statistics."""
        return self._cached_runner.cache.stats()

    def store_stats(self) -> dict:
        """Return IDRegistry and ParameterStore statistics."""
        return self._loader.store_stats()

    @property
    def base_params(self) -> PathwayParameters:
        """The loaded base pathway parameters."""
        return self._base_params

    @property
    def registry(self):
        return self._registry

    @property
    def store(self):
        return self._store

    @property
    def curation_queue(self):
        return self._curation_queue
