"""
ODE system — integrates a pathway model and returns a SimulationResult.

Uses scipy.integrate.solve_ivp with the Radau method (implicit Runge-Kutta,
A-stable, handles stiff systems). Biological ODE systems are almost always
stiff — using a non-stiff solver will give wrong or non-converging results.

SimulationResult stores:
    - Time array in SECONDS (platform canonical unit)
    - All species values in PLATFORM CANONICAL UNITS (μM)
    - Metadata: solver success, model version, simulation parameters

Watchdog:
    The Boolean watchdog runs in parallel with the ODE. When the Boolean
    model's discrete state disagrees with the ODE's continuous state,
    the confidence score decreases. This implements the architecture
    decision from Sprint 0: ODE primary, Boolean as confidence watchdog.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional
import math

import numpy as np
from scipy.integrate import solve_ivp

from .base_pathway import BasePathway
from .pathway import GlucoseInsulinPathway, SPECIES_NAMES, N_SPECIES
from .meal import MealInput


# ---------------------------------------------------------------------------
# Unit conversion: ODE internal → platform canonical
#
# Keyed by internal unit string (as declared by each pathway's
# species_units), not by species position -- this lets ODESystem integrate
# any BasePathway (bare GlucoseInsulinPathway or a multi-pathway
# CompositePathway including cortisol/hepatic/NEFA/GLP1 extensions) without
# knowing its species layout in advance.
# ---------------------------------------------------------------------------
UNIT_CONVERSIONS: dict[str, tuple[str, float]] = {
    "mM":            ("uM", 1000.0),      # concentration: mM → μM
    "pM":            ("uM", 1e-6),        # concentration: pM → μM
    "nM":            ("uM", 1e-3),        # concentration: nM → μM
    "uM":             ("uM", 1.0),         # already canonical concentration
    "min-1":         ("s-1", 1 / 60),     # rate: min⁻¹ → s⁻¹
    "mM/min":        ("uM/s", 1000.0 / 60.0),  # rate-concentration: mM/min → μM/s
    "dimensionless": ("dimensionless", 1.0),
}

# Retained for backwards compatibility with existing callers/tests that
# import the bare-GI-specific unit tables directly.
INTERNAL_TO_CANONICAL_CONCENTRATION = {u: f for u, (_, f) in UNIT_CONVERSIONS.items()}
SPECIES_INTERNAL_UNITS = ["mM", "pM", "min-1", "nM", "dimensionless", "dimensionless"]
SPECIES_CANONICAL_UNITS = ["uM", "uM", "s-1", "uM", "dimensionless", "dimensionless"]


def _convert_species_unit(internal_unit: str) -> tuple[str, float]:
    """Look up (canonical_unit, factor) for an internal unit, failing fast on unknowns."""
    if internal_unit not in UNIT_CONVERSIONS:
        raise ValueError(
            f"Unknown internal species unit '{internal_unit}' -- add it to "
            "UNIT_CONVERSIONS in biosim_ode.system before using it in a pathway."
        )
    return UNIT_CONVERSIONS[internal_unit]


def _find_species_index(species_names: list[str], base_name: str) -> Optional[int]:
    """
    Locate a species by name, tolerating CompositePathway's
    "<pathway_name>.<species>" prefixing (e.g. "glucose_insulin_v1.G").
    Returns None if not present (e.g. a composite without the GI pathway).
    """
    for i, name in enumerate(species_names):
        if name == base_name or name.endswith(f".{base_name}"):
            return i
    return None


@dataclass
class WatchdogResult:
    """
    Result of comparing ODE and Boolean model outputs.

    When all species agree → high confidence (score near 1.0).
    When species disagree → confidence drops (widening uncertainty).

    This is the architecture from Sprint 0:
    ODE is primary, Boolean is confidence watchdog.
    """
    confidence_score: float          # 0.0 to 1.0
    disagreements: list[str]         # Species that disagreed
    ode_states: dict[str, float]     # ODE values at this time point
    boolean_states: dict[str, bool]  # Boolean model states

    @property
    def high_confidence(self) -> bool:
        return self.confidence_score >= 0.8

    @property
    def uncertain(self) -> bool:
        return self.confidence_score < 0.5


@dataclass
class SimulationResult:
    """
    Full result of a simulation run.

    All values are in platform canonical units (μM, s) after conversion.
    Time is always in seconds.

    Species available:
        G     [μM] Blood glucose (baseline ~5000 μM = 5 mM)
        I     [μM] Plasma insulin (baseline ~66e-6 μM = 66 pM)
        X     [s⁻¹] Remote insulin
        IGF1  [μM] Plasma IGF-1 (baseline ~13e-3 μM = 13 nM)
        mTOR  [dimensionless 0-1] mTORC1 activity
        AMPK  [dimensionless 0-1] AMPK activity
    """
    times_s: np.ndarray                     # Time in seconds
    species: dict[str, np.ndarray]          # Species name → values (canonical units)
    species_units: dict[str, str]           # Species name → canonical unit
    success: bool
    message: str
    solver_nfev: int = 0                    # Number of ODE evaluations
    watchdog: Optional[list[WatchdogResult]] = None
    # Point-of-use clamps (pathway.EFFECTIVE_PARAM_BOUNDS) that engaged at
    # any output time point: {(field, "floor"|"ceiling"): most extreme raw
    # value}. Keyed by side too, so a field pushed to its floor early and
    # its ceiling later reports both. Sampled at the output points only, so
    # an excursion lasting less than one output interval can be missed.
    clamp_hits: dict[tuple[str, str], float] = field(default_factory=dict)

    @property
    def duration_s(self) -> float:
        return float(self.times_s[-1] - self.times_s[0]) if len(self.times_s) > 0 else 0.0

    @property
    def duration_min(self) -> float:
        return self.duration_s / 60.0

    def at_time(self, t_s: float) -> dict[str, float]:
        """Interpolate all species values at a specific time (seconds)."""
        result = {}
        for name, values in self.species.items():
            result[name] = float(np.interp(t_s, self.times_s, values))
        return result

    def glucose_peak(self) -> tuple[float, float]:
        """Return (peak_glucose_uM, time_of_peak_s)."""
        g = self.species["G"]
        idx = int(np.argmax(g))
        return float(g[idx]), float(self.times_s[idx])

    def insulin_peak(self) -> tuple[float, float]:
        """Return (peak_insulin_uM, time_of_peak_s)."""
        i = self.species["I"]
        idx = int(np.argmax(i))
        return float(i[idx]), float(self.times_s[idx])

    def time_to_return_to_baseline(
        self,
        species_name: str,
        baseline: float,
        threshold_fraction: float = 0.1,
    ) -> Optional[float]:
        """
        Find the time when a species returns within threshold_fraction of baseline.

        Returns time in seconds, or None if it never returns within threshold.
        """
        values = self.species[species_name]
        baseline_band = abs(baseline) * threshold_fraction
        peak_idx = int(np.argmax(np.abs(values - baseline)))
        for idx in range(peak_idx, len(values)):
            if abs(values[idx] - baseline) <= baseline_band:
                return float(self.times_s[idx])
        return None


class ODESystem:
    """
    Integrates a pathway (bare GlucoseInsulinPathway or a CompositePathway of
    GlucoseInsulinPathway + extensions) and returns a SimulationResult.

    Usage:
        system = ODESystem()

        # Simulate 4 hours after a high-sugar meal
        from biosim_ode.meal import high_sugar_meal
        result = system.simulate_meal(high_sugar_meal(), duration_s=14400)

        # Access blood glucose over time
        glucose = result.species["G"]   # in μM
        times = result.times_s          # in seconds

        # Also works with any BasePathway, e.g. a GI+GLP1 composite:
        composite = glp1_pathway.couple_with_gi(GlucoseInsulinPathway(), GLP1Pathway.semaglutide())
        result = ODESystem(pathway=composite).simulate(duration_s=14400, meals=[...])
    """

    MODEL_VERSION = "glucose_insulin_v1.0"
    DEFAULT_RTOL = 1e-6
    DEFAULT_ATOL = 1e-8

    def __init__(
        self,
        pathway: Optional[BasePathway] = None,
        rtol: float = DEFAULT_RTOL,
        atol: float = DEFAULT_ATOL,
    ) -> None:
        self.pathway = pathway or GlucoseInsulinPathway()
        self.rtol = rtol
        self.atol = atol

    def simulate(
        self,
        duration_s: float,
        meals: Optional[list[MealInput]] = None,
        n_output_points: int = 500,
        y0_override: Optional[list[float]] = None,
        check_clamps: bool = True,
    ) -> SimulationResult:
        """
        Run a simulation for duration_s seconds.

        Args:
            duration_s: Simulation duration in seconds
            meals: List of meal inputs (empty = fasted simulation)
            n_output_points: Number of time points in the output
            y0_override: Optional initial state vector (internal units) to use
                instead of the pathway's own analytic basal state -- e.g. an
                IV glucose bolus, which appears as an instantaneously elevated
                initial glucose concentration rather than a gut-absorption
                Ra(t) forcing term. Must be the same length as pathway.y0.
                Defaults to None (today's behaviour, unchanged).
            check_clamps: Sweep the output trajectory for point-of-use clamp
                hits (SimulationResult.clamp_hits). Adds roughly 20-25% to a
                composite run, so ensemble callers that discard the result
                (Monte Carlo, OAT and Sobol sensitivity) pass False.

        Returns:
            SimulationResult with all species in canonical units (μM, s)
        """
        duration_min = duration_s / 60.0

        # Build pathway with meal inputs (works for a bare GlucoseInsulinPathway
        # or any CompositePathway that contains one -- see with_meals()).
        if getattr(self.pathway, "accepts_meals", False):
            pathway = self.pathway.with_meals(meals or [])
        else:
            pathway = self.pathway

        y0 = list(y0_override) if y0_override is not None else pathway.y0
        t_eval = np.linspace(0.0, duration_min, n_output_points)

        # Solve ODE — Radau for stiff systems
        sol = solve_ivp(
            fun=pathway.dydt,
            t_span=(0.0, duration_min),
            y0=y0,
            method="Radau",
            t_eval=t_eval,
            rtol=self.rtol,
            atol=self.atol,
            dense_output=False,
        )

        if not sol.success:
            return SimulationResult(
                times_s=np.array([]),
                species={},
                species_units={},
                success=False,
                message=f"ODE solver failed: {sol.message}",
                solver_nfev=sol.nfev,
            )

        # Convert to canonical units and build result. Species names/units
        # come from the pathway itself (not a hardcoded GI-only table) so
        # this works identically for composite pathways.
        times_s = sol.t * 60.0  # minutes → seconds
        species = {}
        species_units = {}

        for idx, name in enumerate(pathway.species_names):
            # CompositePathway prefixes names as "<pathway_name>.<species>";
            # strip the prefix so sim.species["G"] etc. work identically
            # whether the pathway is bare or composed with extensions.
            # Species names never collide across a platform pathway set.
            bare_name = name.rsplit(".", 1)[-1]
            internal_unit = pathway.species_units[idx]
            canonical_unit, factor = _convert_species_unit(internal_unit)
            species[bare_name] = sol.y[idx] * factor
            species_units[bare_name] = canonical_unit

        # Run boolean watchdog (only meaningful when G/I/mTOR/AMPK are present)
        watchdog_results = self._run_watchdog(sol.t, sol.y, pathway)
        clamp_hits = _sweep_clamp_hits(pathway, sol.t, sol.y) if check_clamps else {}

        return SimulationResult(
            times_s=times_s,
            species=species,
            species_units=species_units,
            success=True,
            message=f"Solved with Radau ({sol.nfev} evaluations)",
            solver_nfev=sol.nfev,
            watchdog=watchdog_results,
            clamp_hits=clamp_hits,
        )

    def simulate_meal(
        self,
        meal: MealInput,
        duration_s: float = 14400,  # 4 hours
        n_output_points: int = 500,
    ) -> SimulationResult:
        """Convenience method: simulate a single meal."""
        return self.simulate(duration_s=duration_s, meals=[meal], n_output_points=n_output_points)

    def simulate_fasted(
        self,
        duration_s: float = 3600,
        n_output_points: int = 200,
    ) -> SimulationResult:
        """Simulate fasted state — should remain at steady state."""
        return self.simulate(duration_s=duration_s, meals=[], n_output_points=n_output_points)

    # ------------------------------------------------------------------
    # Boolean watchdog
    # ------------------------------------------------------------------

    def _run_watchdog(
        self,
        t_min: np.ndarray,
        y_internal: np.ndarray,
        pathway: BasePathway,
    ) -> list[WatchdogResult]:
        """
        Run the Boolean watchdog at each time point.

        The Boolean model uses physiological thresholds to assign each
        species a binary high/low state. When ODE and Boolean agree
        on direction (both say G is high after a meal), confidence is high.
        Disagreements widen the confidence interval.

        Thresholds (in internal units):
            G:    hyperglycaemic if > 8.0 mM (144 mg/dL)
            I:    elevated if > 150 pM
            mTOR: active if > 0.5
            AMPK: suppressed if < 0.4

        Species are located by name (tolerating CompositePathway's
        "<pathway_name>." prefix) rather than fixed index, so this works for
        both a bare GlucoseInsulinPathway and any composite that includes it.
        Returns an empty list if G/I/mTOR/AMPK are not all present (e.g. a
        composite built without the core GI pathway) -- ExpertMetrics then
        falls back to full confidence, matching today's no-watchdog behaviour.
        """
        names = pathway.species_names
        idx_G, idx_I = _find_species_index(names, "G"), _find_species_index(names, "I")
        idx_mTOR, idx_AMPK = _find_species_index(names, "mTOR"), _find_species_index(names, "AMPK")
        if None in (idx_G, idx_I, idx_mTOR, idx_AMPK):
            return []

        results = []

        for i in range(len(t_min)):
            G    = y_internal[idx_G, i]
            I    = y_internal[idx_I, i]
            mTOR = y_internal[idx_mTOR, i]
            AMPK = y_internal[idx_AMPK, i]

            # Boolean states
            bool_states = {
                "G_elevated":      G > 8.0,
                "I_elevated":      I > 150.0,
                "mTOR_active":     mTOR > 0.5,
                "AMPK_suppressed": AMPK < 0.4,
            }

            # ODE-derived binary states (using same thresholds)
            ode_states = {
                "G": G,
                "I": I,
                "mTOR": mTOR,
                "AMPK": AMPK,
            }

            # Check consistency: if I is elevated, mTOR should be active
            # If I is elevated, AMPK should be suppressed
            # If G is elevated, I should be elevated (broken beta-cell response otherwise)
            #
            # Regression: this used to omit the G_elevated -> I_elevated
            # check entirely (only the two I_elevated-driven rules
            # existed here), and weighted every violation with a flat
            # 0.3 penalty. boolean_pathway.py's WatchdogEngine/
            # GlucoseInsulinBooleanPathway is a complete, tested
            # implementation of all three rules with deliberately
            # differentiated per-rule severities (0.25/0.25/0.30) --
            # its own docstring says it "replaces the basic watchdog in
            # system.py" -- but was never actually wired in anywhere in
            # production (confirmed via grep: referenced only by its own
            # tests). The omitted rule is the highest-severity one and
            # the most directly relevant to this platform's own T2D/
            # impaired-secretion mission (elevated glucose that fails to
            # trigger insulin secretion), so every confidence score this
            # platform has ever reported was blind to exactly the
            # failure mode most worth catching. Species are still
            # located by name above (not WatchdogEngine's fixed
            # y_internal[:6] indexing) to keep this robust for any
            # CompositePathway, not just the bare GlucoseInsulinPathway.
            disagreements = []
            severity = 0.0
            if bool_states["I_elevated"] and not bool_states["mTOR_active"]:
                disagreements.append("mTOR below threshold despite elevated insulin")
                severity += 0.25
            if bool_states["I_elevated"] and not bool_states["AMPK_suppressed"]:
                disagreements.append("AMPK not suppressed despite elevated insulin")
                severity += 0.25
            if bool_states["G_elevated"] and not bool_states["I_elevated"]:
                disagreements.append("Insulin not elevated despite elevated glucose (possible broken beta-cell response)")
                severity += 0.30

            score = max(0.0, 1.0 - severity)

            results.append(WatchdogResult(
                confidence_score=score,
                disagreements=disagreements,
                ode_states=ode_states,
                boolean_states=bool_states,
            ))

        return results

    def overall_confidence(self, result: SimulationResult) -> float:
        """Mean watchdog confidence score across the simulation."""
        if not result.watchdog:
            return 1.0
        return float(np.mean([w.confidence_score for w in result.watchdog]))


def _sweep_clamp_hits(
    pathway: BasePathway, t_min: np.ndarray, y: np.ndarray,
) -> dict[tuple[str, str], float]:
    """
    Recompute, at every output time point, whether any sub-pathway's
    point-of-use clamps engaged, keeping the most extreme raw value per
    (field, side). Works for a bare pathway (signal-free, so time-invariant)
    or a CompositePathway (signals recomputed exactly as dydt() injects
    them). Nested composites are not supported by the platform (see
    SimulationRunner._build_pathway()), so only direct sub-pathways are
    checked.
    """
    from .coupling import CompositePathway

    if not isinstance(pathway, CompositePathway):
        check = getattr(pathway, "point_of_use_clamp_hits", None)
        if check is None:
            return {}
        return {(name, side): raw for name, (side, raw) in check(None).items()}

    checkable = [p for p in pathway.pathways if hasattr(p, "point_of_use_clamp_hits")]
    if not checkable:
        return {}
    worst: dict[tuple[str, str], float] = {}
    for i, t in enumerate(t_min):
        injected = pathway.injected_signals(float(t), list(y[:, i]))
        for sub in checkable:
            for name, (side, raw) in sub.point_of_use_clamp_hits(injected[sub.name]).items():
                key = (name, side)
                worst[key] = _more_extreme(side, worst.get(key), raw)
    return worst


def _more_extreme(side: str, current: Optional[float], raw: float) -> float:
    # NaN is the most extreme floor value: it carries no magnitude at all.
    if current is None or math.isnan(raw):
        return raw
    if math.isnan(current):
        return current
    return min(current, raw) if side == "floor" else max(current, raw)
