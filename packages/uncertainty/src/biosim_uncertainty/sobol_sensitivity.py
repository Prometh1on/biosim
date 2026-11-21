# SAVE TO: packages\uncertainty\src\biosim_uncertainty\sobol_sensitivity.py
"""
Sobol global sensitivity analysis.

`sensitivity.py`'s own module docstring has, since this platform's
earliest sessions, named this as the explicit follow-on to one-at-a-time
(OAT) sensitivity: "Sobol global sensitivity (Sprint 6) addresses both
[OAT] limitations" -- local-only validity, and no captured parameter
interactions -- but no such implementation existed anywhere in the
codebase (confirmed via a repo-wide grep for "Sobol" before writing this
module: the only hit was that one docstring sentence).

Method: variance-based global sensitivity indices (Sobol 2001), computed
via Saltelli's (2002) sampling scheme and SALib's reference
implementation (`SALib.sample.sobol`/`SALib.analyze.sobol`) -- not a
hand-rolled reimplementation. This is a deliberate choice, following this
project's own "prefer a battle-tested library over hand-rolled code"
research-and-reuse rule: the pick-freeze estimator's bias-correction and
confidence-interval bootstrap are easy to get subtly wrong by hand, and
SALib is the standard, widely-cited reference implementation in the
global-sensitivity-analysis literature (Herman & Usher 2017, J. Open
Source Software, describes SALib itself).

Two indices per (parameter, metric) pair:
    S1 (first-order):  fraction of output variance explained by this
                        parameter ALONE (its main effect).
    ST (total-order):  fraction of output variance explained by this
                        parameter's main effect PLUS every interaction
                        it participates in with other parameters.
    ST > S1 for a parameter means it has real interaction effects OAT
    (which perturbs one parameter at a time, holding all others fixed)
    cannot detect at all.

Parameter distributions are NOT arbitrary bounds invented for this
analysis -- they reuse this platform's own already-established,
literature-anchored `LogNormalDistribution.from_confidence_tier()`
objects (the exact same distributions `MonteCarloEngine.propagate()`
already draws from for Monte Carlo UQ), via SALib's `dists=['lognorm']`
problem-definition option, which takes each parameter's (mu, sigma) in
the same log-space parameterisation `LogNormalDistribution` itself
already computes (`mu`, `sigma` public properties) -- not [min, max] bounds
that would need a separate, disconnected justification.

Order-sensitivity caveat, unlike `MonteCarloEngine.propagate()`: a
Sobol analysis is NOT simply a large batch of independent samples --
the pick-freeze estimator depends on the EXACT structure of the
Saltelli sample matrix (the A matrix, B matrix, and k A_i "one-column-
swapped" matrices). Dropping or replacing a single failed row, the way
`MonteCarloEngine.propagate()`/`OATSensitivityAnalyzer._run()` tolerate
individual sample failures, would silently corrupt every downstream
index, not just the one estimate that failed -- there is no valid way
to "skip a failure and continue" for this method. If any model
evaluation fails, `SobolSensitivityAnalyzer.analyze()` raises
`RuntimeError` listing every failure, rather than silently computing
indices over a matrix with a substituted or dropped value.
"""

from __future__ import annotations

import dataclasses
import inspect
from dataclasses import dataclass
from typing import Callable, Mapping, Optional
import numpy as np

from biosim_ode import ODESystem, GlucoseInsulinPathway, PathwayParameters, SimulationResult
from biosim_ode.meal import MealInput

from .distributions import LogNormalDistribution
from .propagation import ExtensionPathwayFactory, PathwayFactory
from .sensitivity import STANDARD_METRICS


@dataclass
class SobolResult:
    """
    Result of a Sobol global sensitivity analysis.

    first_order[metric][param]  = S1 (main-effect index)
    total_order[metric][param]  = ST (main effect + all interactions)
    *_conf holds each index's bootstrap 95% CI half-width (SALib's
    `_conf` convention) -- an index whose magnitude is smaller than its
    own confidence interval should be read as "not distinguishable from
    zero at this sample size," not as a confidently-small effect.
    """
    first_order:      dict[str, dict[str, float]]
    first_order_conf: dict[str, dict[str, float]]
    total_order:      dict[str, dict[str, float]]
    total_order_conf: dict[str, dict[str, float]]
    parameter_names:  list[str]
    metric_names:     list[str]
    n_base_samples:   int
    n_model_evaluations: int

    def top_parameters(self, metric: str, n: int = 5, by: str = "total_order") -> list[tuple[str, float]]:
        """
        Return the top n most influential parameters for a given metric,
        ranked by total-order index by default (captures interactions;
        first-order alone can understate a parameter's real influence).
        """
        table = self.total_order if by == "total_order" else self.first_order
        if metric not in table:
            raise KeyError(f"Metric '{metric}' not in result. Available: {self.metric_names}")
        ranked = sorted(table[metric].items(), key=lambda x: abs(x[1]), reverse=True)
        return ranked[:n]

    def interaction_strength(self, metric: str, param: str) -> float:
        """
        ST - S1 for one (metric, param) pair -- the portion of this
        parameter's total influence that comes from interactions with
        OTHER parameters, exactly the effect OAT sensitivity analysis
        cannot see (it perturbs one parameter at a time, holding every
        other parameter fixed at its base value).
        """
        return self.total_order[metric][param] - self.first_order[metric][param]

    def summary(self) -> dict:
        top = {}
        for metric in self.metric_names:
            top[metric] = {
                p: {"S1": round(self.first_order[metric][p], 4),
                    "ST": round(self.total_order[metric][p], 4)}
                for p, _ in self.top_parameters(metric, n=5)
            }
        return {
            "n_base_samples": self.n_base_samples,
            "n_model_evaluations": self.n_model_evaluations,
            "top_sensitivities": top,
        }


class SobolSensitivityAnalyzer:
    """
    Sobol global sensitivity analyzer.

    Total model evaluations = n_base_samples * (2 * k + 2), where k is
    the number of parameters being analyzed (Saltelli's extended
    sampling scheme, `calc_second_order=False` variant -- second-order
    interaction indices are not computed by this analyzer, since this
    platform's own ODE parameter set is small enough that the total-
    order index already captures "does this parameter interact with
    anything," which is the actionable question for the curation-
    prioritisation use case this feeds; a full second-order breakdown
    would double the sample cost for a question this project doesn't
    currently need answered).
    """

    def __init__(
        self,
        base_params: Optional[PathwayParameters] = None,
        pathway_factory: Optional[PathwayFactory] = None,
        extension_pathway_factory: Optional[ExtensionPathwayFactory] = None,
    ) -> None:
        """
        Args:
            base_params: base PathwayParameters; Saltelli-sampled draws
                override individual fields.
            pathway_factory: optional (params, meals) -> BasePathway
                callable, used in place of a bare GlucoseInsulinPathway
                for every sample -- same rationale and contract as
                OATSensitivityAnalyzer's own pathway_factory and
                propagation.PathwayFactory's docstring. Defaults to None
                (unchanged prior behaviour). Mutually exclusive with
                extension_pathway_factory.
            extension_pathway_factory (MASTER_PLAN 2.4): optional
                (params, meals, extension_overrides) -> BasePathway
                callable -- the identical additive treatment already
                given to MonteCarloEngine.propagate() and
                OATSensitivityAnalyzer (see
                propagation.ExtensionPathwayFactory's own docstring).
                Unlike OATSensitivityAnalyzer (which needs a separate
                extension_base_values to perturb FROM), Sobol needs no
                such thing -- analyze()'s own `distributions` dict
                already carries a distribution (mean+CV) for every
                sampled parameter, core or extension alike; a dotted
                name (e.g. "renal.beta") in that SAME dict is routed to
                extension_overrides at evaluation time, not through
                dataclasses.replace(self.base_params, ...) (which would
                raise -- PathwayParameters is a frozen dataclass with no
                such field). Mutually exclusive with pathway_factory.
        """
        if pathway_factory is not None and extension_pathway_factory is not None:
            raise ValueError(
                "SobolSensitivityAnalyzer accepts at most one of pathway_factory/"
                "extension_pathway_factory, not both -- there is no defined "
                "precedence between them."
            )
        if extension_pathway_factory is not None:
            # Same eager arity check as MonteCarloEngine.propagate()/
            # OATSensitivityAnalyzer -- fail once, at construction, rather
            # than identically on every one of n_base_samples*(2k+2) model
            # evaluations.
            try:
                inspect.signature(extension_pathway_factory).bind(None, None, {})
            except TypeError as e:
                raise ValueError(
                    f"extension_pathway_factory must accept 3 positional arguments "
                    f"(params, meals, extension_overrides) -- got a signature "
                    f"incompatible with that call: {e}"
                ) from e
        self.base_params = base_params or PathwayParameters()
        self.pathway_factory = pathway_factory
        self.extension_pathway_factory = extension_pathway_factory

    def analyze(
        self,
        distributions: dict[str, LogNormalDistribution],
        metrics: Optional[dict[str, Callable[[SimulationResult], float]]] = None,
        meals: Optional[list[MealInput]] = None,
        duration_s: float = 14400.0,
        n_base_samples: int = 32,
        seed: Optional[int] = None,
    ) -> SobolResult:
        """
        Run Sobol global sensitivity analysis.

        Args:
            distributions: parameter_name -> LogNormalDistribution,
                the same distributions MonteCarloEngine.propagate()
                already draws from (see this module's docstring for why
                these, not arbitrary [min, max] bounds).
            metrics: metric_name -> callable(SimulationResult) -> float.
                Defaults to all STANDARD_METRICS.
            n_base_samples: SALib's own convention (Saltelli 2002)
                recommends a power of 2 for the underlying Sobol
                sequence's low-discrepancy property to hold; not
                enforced here (SALib itself only warns), but respected
                by every caller in this codebase's own tests/usage.
            seed: for reproducibility. None means SALib's own default
                (numpy's global RNG state), matching this analyzer's own
                lack of any other reproducibility guarantee absent one.

        Returns:
            SobolResult with S1/ST indices and their confidence
            intervals per (metric, parameter) pair.

        Raises:
            RuntimeError: if any Saltelli sample fails to simulate --
                see this module's docstring for why a Sobol analysis
                cannot tolerate a dropped/substituted failure the way
                MonteCarloEngine.propagate() can.
        """
        from SALib.sample import sobol as sobol_sample
        from SALib.analyze import sobol as sobol_analyze

        metrics = metrics or STANDARD_METRICS
        meals = meals or []
        param_names = list(distributions.keys())
        k = len(param_names)
        if k == 0:
            raise ValueError("distributions must contain at least one parameter")
        # Code-review finding: a dotted extension name (e.g. "renal.beta")
        # in `distributions` with no extension_pathway_factory set used to
        # be silently dropped from the pathway entirely -- every sample's
        # zero-variance output on that column then crashed deep inside
        # SALib.analyze.sobol with an opaque TypeError, instead of a clear
        # diagnostic naming the actual mistake. Checked once here (not
        # per-row), mirroring propagate()'s own identical guard.
        extension_names = [name for name in param_names if "." in name]
        if extension_names and self.extension_pathway_factory is None:
            raise ValueError(
                f"distributions key(s) {sorted(extension_names)} look like extension-"
                f"pathway parameters (dotted names), but extension_pathway_factory was "
                f"not supplied -- they would have nowhere to go and be silently "
                f"discarded. Pass extension_pathway_factory too."
            )

        problem = {
            "num_vars": k,
            "names": param_names,
            "bounds": [
                [distributions[name].mu, distributions[name].sigma]
                for name in param_names
            ],
            "dists": ["lognorm"] * k,
        }

        param_matrix = sobol_sample.sample(
            problem, N=n_base_samples, calc_second_order=False, seed=seed,
        )
        n_evaluations = param_matrix.shape[0]

        outputs: dict[str, np.ndarray] = {
            name: np.full(n_evaluations, np.nan) for name in metrics
        }
        failure_reasons: list[str] = []

        for row_idx, row in enumerate(param_matrix):
            row_values = dict(zip(param_names, row))
            # Dotted names (e.g. "renal.beta") route to extension_overrides,
            # not dataclasses.replace() -- PathwayParameters is a frozen
            # dataclass with a fixed field set and has no such field; this
            # mirrors MonteCarloEngine.propagate()'s own dotted-key
            # convention (see ExtensionPathwayFactory's docstring).
            core_values = {k: v for k, v in row_values.items() if "." not in k}
            extension_values = {k: v for k, v in row_values.items() if "." in k}
            params = dataclasses.replace(self.base_params, **core_values)
            result, reason = self._run(params, meals, duration_s, extension_values)
            if result is None:
                failure_reasons.append(f"sample {row_idx} ({row_values}): {reason}")
                continue
            for metric_name, metric_fn in metrics.items():
                outputs[metric_name][row_idx] = metric_fn(result)

        if failure_reasons:
            raise RuntimeError(
                f"{len(failure_reasons)}/{n_evaluations} Sobol sample evaluations "
                "failed -- a Sobol analysis cannot tolerate a dropped or "
                "substituted sample without corrupting every downstream index "
                "(see sobol_sensitivity.py's module docstring). Failures:\n"
                + "\n".join(failure_reasons)
            )

        first_order: dict[str, dict[str, float]] = {}
        first_order_conf: dict[str, dict[str, float]] = {}
        total_order: dict[str, dict[str, float]] = {}
        total_order_conf: dict[str, dict[str, float]] = {}

        for metric_name in metrics:
            sobol_indices = sobol_analyze.analyze(
                problem, outputs[metric_name], calc_second_order=False, seed=seed,
            )
            first_order[metric_name] = dict(zip(param_names, sobol_indices["S1"]))
            first_order_conf[metric_name] = dict(zip(param_names, sobol_indices["S1_conf"]))
            total_order[metric_name] = dict(zip(param_names, sobol_indices["ST"]))
            total_order_conf[metric_name] = dict(zip(param_names, sobol_indices["ST_conf"]))

        return SobolResult(
            first_order=first_order,
            first_order_conf=first_order_conf,
            total_order=total_order,
            total_order_conf=total_order_conf,
            parameter_names=param_names,
            metric_names=list(metrics.keys()),
            n_base_samples=n_base_samples,
            n_model_evaluations=n_evaluations,
        )

    def _run(
        self,
        params: PathwayParameters,
        meals: list[MealInput],
        duration_s: float,
        extension_overrides: Optional[Mapping[str, float]] = None,
    ) -> tuple[Optional[SimulationResult], Optional[str]]:
        """Run one simulation. Returns (None, reason) if it fails -- same shape as OATSensitivityAnalyzer._run()."""
        try:
            if self.extension_pathway_factory is not None:
                pathway = self.extension_pathway_factory(params, meals, extension_overrides or {})
            elif self.pathway_factory is not None:
                pathway = self.pathway_factory(params, meals)
            else:
                pathway = GlucoseInsulinPathway(params=params, meals=meals)
            system = ODESystem(pathway=pathway)
            result = system.simulate(
                duration_s=duration_s, meals=meals, n_output_points=200,
                check_clamps=False,
            )
            if result.success:
                return result, None
            return None, f"ODE solver failed: {result.message}"
        except Exception as e:
            return None, f"{type(e).__name__}: {e}"
