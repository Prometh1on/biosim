# SAVE TO: packages\uncertainty\src\biosim_uncertainty\sensitivity.py
"""
One-at-a-time (OAT) sensitivity analysis.

Identifies which parameters most affect which simulation outputs.
This feeds back to the curation pipeline: parameters with high
sensitivity need tighter uncertainty (more curation effort first).

Method:
    For each parameter θ_i, compute the normalised sensitivity index:

        S_i = Δoutput/output ÷ Δθ_i/θ_i

    Using central differences:
        S_i ≈ (f(θ_i*(1+δ)) - f(θ_i*(1-δ))) / (2 * δ * f(θ_i))

    where δ = perturbation_fraction (default 0.10 = ±10%).

    S_i = 1.0 means a 10% change in θ_i causes a 10% change in output.
    S_i = 0.0 means the output is insensitive to θ_i.
    S_i < 0 means an inverse relationship (more θ_i → less output).

Limitations of OAT:
    - Local: valid only near the base parameter set
    - Does not capture parameter interactions
    - Sobol global sensitivity (Sprint 6) addresses both limitations

Common output metrics:
    - glucose_peak: maximum blood glucose after meal
    - insulin_peak: maximum plasma insulin after meal
    - glucose_auc: area under glucose curve (total glucose exposure)
    - time_to_return: time for glucose to return within 10% of baseline
"""

from __future__ import annotations

import dataclasses
import inspect
from dataclasses import dataclass, field
from typing import Callable, Mapping, Optional
import numpy as np

from biosim_ode import ODESystem, GlucoseInsulinPathway, PathwayParameters, SimulationResult
from biosim_ode.meal import MealInput

from .propagation import ExtensionPathwayFactory, PathwayFactory


# ---------------------------------------------------------------------------
# Standard output metrics — callable from SimulationResult → scalar
# ---------------------------------------------------------------------------

def glucose_peak_metric(result: SimulationResult) -> float:
    """Peak blood glucose [μM]."""
    return result.glucose_peak()[0]


def insulin_peak_metric(result: SimulationResult) -> float:
    """Peak plasma insulin [μM]."""
    return result.insulin_peak()[0]


def glucose_auc_metric(result: SimulationResult) -> float:
    """Area under the glucose curve [μM * s]."""
    return float(np.trapezoid(result.species["G"], result.times_s))


def mtor_peak_metric(result: SimulationResult) -> float:
    """Peak mTOR activity [dimensionless]."""
    return float(np.max(result.species["mTOR"]))


def ampk_nadir_metric(result: SimulationResult) -> float:
    """Minimum AMPK activity [dimensionless] — lower = more suppressed."""
    return float(np.min(result.species["AMPK"]))


STANDARD_METRICS = {
    "glucose_peak":  glucose_peak_metric,
    "insulin_peak":  insulin_peak_metric,
    "glucose_auc":   glucose_auc_metric,
    "mtor_peak":     mtor_peak_metric,
    "ampk_nadir":    ampk_nadir_metric,
}


@dataclass
class SensitivityResult:
    """
    Result of a one-at-a-time sensitivity analysis.

    sensitivity_indices[metric_name][param_name] = S_i
    where S_i is the normalised sensitivity index.
    """
    sensitivity_indices: dict[str, dict[str, float]]
    base_outputs: dict[str, float]
    parameter_names: list[str]
    metric_names: list[str]
    perturbation_fraction: float
    n_failed: int = 0
    failure_reasons: list[str] = field(default_factory=list)

    def top_parameters(self, metric: str, n: int = 5) -> list[tuple[str, float]]:
        """
        Return the top n most influential parameters for a given metric,
        sorted by absolute sensitivity (descending).
        """
        if metric not in self.sensitivity_indices:
            raise KeyError(f"Metric '{metric}' not in result. Available: {self.metric_names}")
        scores = self.sensitivity_indices[metric]
        ranked = sorted(scores.items(), key=lambda x: abs(x[1]), reverse=True)
        return ranked[:n]

    def summary(self) -> dict:
        top = {}
        for metric in self.metric_names:
            top[metric] = {
                p: round(s, 4) for p, s in self.top_parameters(metric, n=5)
            }
        return {
            "perturbation_fraction": self.perturbation_fraction,
            "n_failed": self.n_failed,
            "top_sensitivities": top,
        }


class OATSensitivityAnalyzer:
    """
    One-at-a-time sensitivity analyzer.

    For each parameter, runs two simulations (perturbed up and down)
    and computes the normalised sensitivity index for each metric.

    Total simulations = 2 * n_parameters + 1 (one baseline).
    """

    def __init__(
        self,
        perturbation_fraction: float = 0.10,
        base_params: Optional[PathwayParameters] = None,
        pathway_factory: Optional[PathwayFactory] = None,
        extension_pathway_factory: Optional[ExtensionPathwayFactory] = None,
        extension_base_values: Optional[dict[str, float]] = None,
    ) -> None:
        """
        Args:
            perturbation_fraction: fractional perturbation for each OAT step.
            base_params: base PathwayParameters; perturbations override
                individual fields.
            pathway_factory: optional (params, meals) -> BasePathway
                callable, used in place of a bare GlucoseInsulinPathway
                for every simulation. Defaults to None (unchanged prior
                behaviour) -- see propagation.PathwayFactory's docstring
                for why this matters (a composed scenario's sensitivity
                indices would otherwise be computed against a
                structurally different ODE system than the one the
                point estimate actually runs). Mutually exclusive with
                extension_pathway_factory.
            extension_pathway_factory (MASTER_PLAN 2.4): optional
                (params, meals, extension_overrides) -> BasePathway
                callable -- the identical additive treatment already
                given to MonteCarloEngine.propagate() (see
                propagation.ExtensionPathwayFactory's own docstring),
                letting analyze() also compute sensitivity indices for
                extension-pathway-specific parameters (e.g.
                RenalClearanceParameters.beta) that don't live on
                PathwayParameters at all. Pass a dotted name (e.g.
                "renal.beta") in analyze()'s parameter_names to analyze
                it -- routed via extension_base_values below, not
                base_params. Mutually exclusive with pathway_factory.
            extension_base_values: dotted extension-parameter name ->
                point-estimate value to perturb FROM (e.g.
                {"renal.beta": 0.40}) -- the OAT analogue of
                extension_distributions' own mean, since OAT perturbs a
                single center point rather than sampling a distribution.
                Held fixed (unperturbed) whenever a DIFFERENT parameter
                is being analyzed, exactly like every core PathwayParameters
                field not currently being perturbed stays at its own
                base_params value. Ignored if extension_pathway_factory
                is None.
        """
        if not (0.0 < perturbation_fraction < 1.0):
            raise ValueError(
                f"perturbation_fraction must be in (0, 1), got {perturbation_fraction}"
            )
        if pathway_factory is not None and extension_pathway_factory is not None:
            raise ValueError(
                "OATSensitivityAnalyzer accepts at most one of pathway_factory/"
                "extension_pathway_factory, not both -- there is no defined "
                "precedence between them."
            )
        if extension_pathway_factory is not None:
            # Same eager arity check as MonteCarloEngine.propagate() --
            # fail once, at construction, rather than identically on
            # every one of 2*n_parameters+1 simulations.
            try:
                inspect.signature(extension_pathway_factory).bind(None, None, {})
            except TypeError as e:
                raise ValueError(
                    f"extension_pathway_factory must accept 3 positional arguments "
                    f"(params, meals, extension_overrides) -- got a signature "
                    f"incompatible with that call: {e}"
                ) from e
        # Code-review finding: without this guard, extension_base_values
        # supplied with no extension_pathway_factory silently produced a
        # misleading "this parameter has zero influence" result -- the
        # perturbed value never actually reached the pathway (upper/lower
        # runs were identical), the exact "computed but discarded" bug
        # class this project has repeatedly found and fixed elsewhere.
        # Mirrors propagate()'s own identical guard.
        if extension_base_values and extension_pathway_factory is None:
            raise ValueError(
                "extension_base_values was supplied but extension_pathway_factory "
                "was not -- any perturbed extension value would have nowhere to go "
                "and be silently discarded. Pass extension_pathway_factory too."
            )
        self.perturbation_fraction = perturbation_fraction
        self.base_params = base_params or PathwayParameters()
        self.pathway_factory = pathway_factory
        self.extension_pathway_factory = extension_pathway_factory
        self.extension_base_values = extension_base_values or {}

    def analyze(
        self,
        parameter_names: list[str],
        metrics: Optional[dict[str, Callable[[SimulationResult], float]]] = None,
        meals: Optional[list[MealInput]] = None,
        duration_s: float = 14400.0,
    ) -> SensitivityResult:
        """
        Run OAT sensitivity analysis.

        Args:
            parameter_names: List of PathwayParameters field names to
                analyze, OR dotted extension-parameter names (e.g.
                "renal.beta") present in extension_base_values -- routed
                automatically by whether the name contains "." (mirrors
                MonteCarloEngine.propagate()'s own dotted-key convention).
            metrics: Dict of metric_name → callable(SimulationResult) → float.
                     Defaults to all STANDARD_METRICS.
            meals: Meal inputs. Defaults to no meal.
            duration_s: Simulation duration in seconds.

        Returns:
            SensitivityResult with normalised sensitivity indices
        """
        metrics = metrics or STANDARD_METRICS
        meals = meals or []
        n_failed = 0
        failure_reasons: list[str] = []

        # Run baseline -- held at extension_base_values for every extension
        # parameter (empty dict, a no-op, when extension_pathway_factory
        # is unused), the same "every OTHER parameter stays at its own
        # base value" convention core PathwayParameters fields already get
        # via self.base_params.
        base_result, base_reason = self._run(self.base_params, meals, duration_s, self.extension_base_values)
        if base_result is None:
            raise RuntimeError(
                f"Baseline simulation failed — cannot run sensitivity analysis: {base_reason}"
            )

        base_outputs = {
            name: metric(base_result)
            for name, metric in metrics.items()
        }

        # Run perturbed simulations
        sensitivity_indices: dict[str, dict[str, float]] = {m: {} for m in metrics}

        for param_name in parameter_names:
            is_extension = "." in param_name
            base_val = (
                self.extension_base_values.get(param_name) if is_extension
                else getattr(self.base_params, param_name, None)
            )
            if base_val is None:
                continue
            if not isinstance(base_val, (int, float)):
                continue
            if base_val == 0:
                continue  # Cannot compute relative sensitivity at zero

            delta = abs(base_val) * self.perturbation_fraction
            lower_val = max(base_val - delta, base_val * 1e-6)  # keep positive

            if is_extension:
                upper_result, upper_reason = self._run(
                    self.base_params, meals, duration_s,
                    {**self.extension_base_values, param_name: base_val + delta},
                )
                lower_result, lower_reason = self._run(
                    self.base_params, meals, duration_s,
                    {**self.extension_base_values, param_name: lower_val},
                )
            else:
                upper_result, upper_reason = self._run(
                    dataclasses.replace(self.base_params, **{param_name: base_val + delta}),
                    meals, duration_s, self.extension_base_values,
                )
                lower_result, lower_reason = self._run(
                    dataclasses.replace(self.base_params, **{param_name: lower_val}),
                    meals, duration_s, self.extension_base_values,
                )

            if upper_result is None or lower_result is None:
                n_failed += 1
                if upper_result is None:
                    failure_reasons.append(f"{param_name} (upper perturbation): {upper_reason}")
                if lower_result is None:
                    failure_reasons.append(f"{param_name} (lower perturbation): {lower_reason}")
                for metric_name in metrics:
                    sensitivity_indices[metric_name][param_name] = float('nan')
                continue

            for metric_name, metric_fn in metrics.items():
                f_upper = metric_fn(upper_result)
                f_lower = metric_fn(lower_result)
                f_base = base_outputs[metric_name]

                if f_base == 0 or delta == 0:
                    s = float('nan')
                else:
                    # Central difference normalised sensitivity
                    df = f_upper - f_lower
                    dparam = (base_val + delta) - lower_val
                    s = (df / f_base) / (dparam / base_val)

                sensitivity_indices[metric_name][param_name] = s

        return SensitivityResult(
            sensitivity_indices=sensitivity_indices,
            base_outputs=base_outputs,
            parameter_names=parameter_names,
            metric_names=list(metrics.keys()),
            perturbation_fraction=self.perturbation_fraction,
            n_failed=n_failed,
            failure_reasons=failure_reasons,
        )

    def _run(
        self,
        params: PathwayParameters,
        meals: list[MealInput],
        duration_s: float,
        extension_overrides: Optional[Mapping[str, float]] = None,
    ) -> tuple[Optional[SimulationResult], Optional[str]]:
        """
        Run one simulation. Returns (None, reason) if it fails.

        Regression: this used to return a bare None on failure -- a real
        solver bug and a benign edge-case parameter perturbation were
        indistinguishable in SensitivityResult.n_failed. Same bug class
        as CohortSummary.failure_reasons/UQResult.failure_reasons.
        """
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
