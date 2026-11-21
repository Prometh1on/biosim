# SAVE TO: packages\uncertainty\src\biosim_uncertainty\result.py
"""
Result types for uncertainty quantification.

A UQResult stores the full ensemble of simulation outputs and
derived confidence bands. Everything is in platform canonical
units (μM for concentrations, seconds for time).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional
import numpy as np


@dataclass
class ConfidenceBand:
    """
    95% confidence band for a single species across time.

    All arrays have the same length as times_s in the parent UQResult.
    Units match the species' canonical unit (μM, s⁻¹, or dimensionless).
    """
    lower:  np.ndarray   # 2.5th percentile
    median: np.ndarray   # 50th percentile
    upper:  np.ndarray   # 97.5th percentile

    @property
    def width(self) -> np.ndarray:
        """Width of the 95% CI at each time point."""
        return self.upper - self.lower

    @property
    def relative_width(self) -> np.ndarray:
        """Width as fraction of median (coefficient of variation over time)."""
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.where(self.median != 0, self.width / np.abs(self.median), 0.0)

    def at_time_index(self, idx: int) -> tuple[float, float, float]:
        """Return (lower, median, upper) at a specific time index."""
        return float(self.lower[idx]), float(self.median[idx]), float(self.upper[idx])

    def mean_relative_width(self, floor_fraction: float = 0.05) -> float:
        """
        Mean relative width across time points where the median isn't
        close to zero relative to this band's own peak magnitude.

        relative_width's median!=0 guard only special-cases EXACT zero --
        a median that merely passes NEAR zero (a decaying bookkeeping
        compartment, e.g. GlucoseInsulinPathway's X species, exactly 0 at
        t=0 and decaying back toward ~0) produces an enormous
        relative_width at those points from an ordinary-sized numerator
        over a tiny denominator, which then dominates the mean and made
        most_uncertain_species() answer "X" for ordinary scenarios -- a
        bookkeeping artifact, not a genuine uncertainty signal. Excludes
        any time point where |median| falls below floor_fraction times
        this band's own peak |median| before averaging the remainder.
        floor_fraction (default 0.05) is an engineering choice, same
        category as MIN_GB_THRESHOLD_MARGIN_MM -- not a literature-
        derived statistic.
        """
        rw = self.relative_width
        finite = np.isfinite(rw)
        peak = float(np.max(np.abs(self.median))) if len(self.median) else 0.0
        if peak > 0:
            finite = finite & (np.abs(self.median) >= floor_fraction * peak)
        values = rw[finite]
        return float(np.mean(values)) if len(values) > 0 else 0.0


@dataclass
class UQResult:
    """
    Full uncertainty quantification result for a simulation scenario.

    Contains:
        - times_s: time array in seconds
        - bands: 95% confidence bands for each species
        - raw_trajectories: all N simulation trajectories (optional)
        - n_samples: number of Monte Carlo samples attempted
        - n_failed: number that failed to solve
        - overall_confidence: mean watchdog confidence across the ensemble
    """
    times_s:           np.ndarray
    bands:             dict[str, ConfidenceBand]
    n_samples:         int
    n_failed:          int
    overall_confidence: float
    parameter_names:   list[str] = field(default_factory=list)
    raw_trajectories:  Optional[dict[str, np.ndarray]] = None   # shape (n_samples, n_times)
    # Each successful sample's own sampled PathwayParameters value, per
    # parameter name -- same opt-in lifetime as raw_trajectories (only
    # populated when propagate(store_trajectories=True)). Added for
    # BayesianCalibrator (calibration.py): Sampling-Importance-Resampling
    # needs both the trajectory AND the parameter draw that produced it
    # for every sample, to resample parameter draws by posterior weight.
    # Before this field existed, propagate() sampled each parameter set,
    # used it to run the ODE, then discarded it -- only the trajectory
    # survived.
    sampled_parameter_values: Optional[dict[str, np.ndarray]] = None   # shape (n_successful,) per param name
    # Regression: n_failed used to be the only signal a caller had --
    # a real solver bug failing most/all samples was indistinguishable
    # from a few benign edge-case parameter draws. Same bug class as
    # CohortSummary.failure_reasons (biosim_seeds.population) and
    # TimeSlot.failure_reason (biosim_seeds.circadian_optimizer).
    failure_reasons:   list[str] = field(default_factory=list)
    # Maps each band's bare species name (the same keys as `bands`) back
    # to the name of the pathway that produced it -- e.g. "GLP1" ->
    # "GLP1Pathway", "G" -> "glucose_insulin_v1". Purely presentational
    # (species names are already globally unique across this platform's
    # pathway set, so nothing depends on this for correctness); useful
    # for a UI grouping confidence bands by pathway rather than listing
    # every band flat. Populated by MonteCarloEngine.propagate() from
    # the composed pathway's own species_names (which CompositePathway
    # prefixes as "<pathway_name>.<species>", the same prefix
    # ODESystem.simulate() strips to build `bands`' own bare keys) --
    # empty for a UQResult built any other way (e.g. hand-constructed in
    # a test), same graceful-degradation convention as every other
    # optional field on this dataclass.
    species_pathway:   dict[str, str] = field(default_factory=dict)

    @property
    def n_successful(self) -> int:
        return self.n_samples - self.n_failed

    @property
    def success_rate(self) -> float:
        return self.n_successful / self.n_samples if self.n_samples > 0 else 0.0

    @property
    def is_reliable(self) -> bool:
        """True if ≥90% of samples succeeded and overall confidence is high."""
        return self.success_rate >= 0.90 and self.overall_confidence >= 0.7

    def pathway_for(self, species: str) -> Optional[str]:
        """
        Name of the pathway that produced a given species' band, if known.

        Returns None if `species_pathway` was never populated (e.g. a
        hand-constructed UQResult in a test) or doesn't include this
        species -- callers should treat that as "unknown", not an error.
        """
        return self.species_pathway.get(species)

    def band(self, species: str) -> ConfidenceBand:
        """Get the confidence band for a species. Raises KeyError if not found."""
        if species not in self.bands:
            raise KeyError(
                f"Species '{species}' not in UQ result. "
                f"Available: {list(self.bands.keys())}"
            )
        return self.bands[species]

    def at_time(self, t_s: float) -> dict[str, tuple[float, float, float]]:
        """
        Interpolate confidence bands at a specific time.
        Returns dict of species → (lower, median, upper).
        """
        result = {}
        for name, band in self.bands.items():
            lower  = float(np.interp(t_s, self.times_s, band.lower))
            median = float(np.interp(t_s, self.times_s, band.median))
            upper  = float(np.interp(t_s, self.times_s, band.upper))
            result[name] = (lower, median, upper)
        return result

    def most_uncertain_species(self, candidates: Optional[list[str]] = None) -> str:
        """
        Return the name of the species with the highest mean relative
        uncertainty, optionally restricted to a caller-supplied candidate
        list (e.g. the six core GlucoseInsulinPathway species, to exclude
        an extension pathway's own bookkeeping states -- meal{i}_remaining
        _mM decays monotonically toward 0 by construction, the same
        near-zero-median pattern mean_relative_width()'s floor already
        guards against for X, but a caller may still want to scope the
        comparison explicitly). Defaults to every band, as before.
        """
        names = candidates if candidates is not None else list(self.bands.keys())
        return max(names, key=lambda s: self.bands[s].mean_relative_width())

    def uncertainty_summary(self) -> dict:
        """Summary dict for reporting and monitoring."""
        return {
            "n_samples": self.n_samples,
            "n_successful": self.n_successful,
            "success_rate": round(self.success_rate, 3),
            "n_failed": self.n_failed,
            "overall_confidence": round(self.overall_confidence, 3),
            "is_reliable": self.is_reliable,
            "species_uncertainty": {
                name: round(band.mean_relative_width(), 3)
                for name, band in self.bands.items()
            },
            "species_pathway": dict(self.species_pathway),
        }
