# SAVE TO: packages\uncertainty\src\biosim_uncertainty\distributions.py
"""
Parameter distributions for uncertainty quantification.

Kinetic parameters are not known precisely — they are measured with
experimental error, vary across cell types, and differ between papers.
This module represents each parameter as a probability distribution
rather than a point estimate.

Why log-normal for kinetic parameters:
    - Kinetic parameters (Km, kcat, Vmax) are always positive
    - Their uncertainty is multiplicative, not additive
      ('the true value is within 2-fold' is more meaningful than '±N units')
    - Log-normal is the maximum entropy distribution for positive quantities
      with known geometric mean and geometric standard deviation
    - Empirically, kinetic parameters across species/conditions are
      log-normally distributed (Bar-Even et al. 2011, Nat Chem Biol)

Coefficient of variation (CV) defaults by confidence tier:
    GOLD:   CV = 0.10  (±10%)
    SILVER: CV = 0.30  (±30%)
    BRONZE: CV = 0.60  (±60%)

These three specific numbers are a documented EXPERT-JUDGMENT DEFAULT,
not a directly-cited empirical statistic -- unlike the log-normal
functional form above (which Bar-Even et al. 2011 does support), no
single-number citation exists for "typical measurement CV by tier" and
one was not fabricated to look like one. This isn't a gap specific to
this codebase: the broader kinetic-parameter literature itself
acknowledges the underlying quantity is poorly characterized (e.g.
Marucci L, et al. "Simulating Biological Systems using Bio-PEPA."
Systems Biology, 2011, and reviews on enzyme kinetics database curation
note that "the uncertainty of measured kinetic values is poorly
characterized" field-wide -- reproducibility studies exist for specific
narrow comparisons, e.g. cross-method CV as low as ~5% when reanalyzing
the SAME dataset by different curve-fitting methods, but that measures
something different from cross-STUDY variability, which is what these
tiers are meant to represent). Same category of honest, isolated
assumption as `HPA_FEEDBACK_GAIN` in `cortisol_pathway.py` -- replace if
a real distribution of same-parameter cross-study CVs is ever published
and shown broadly applicable across the parameter types this platform
uses (Km, kcat, EC50, rate constants), not before.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
import math
import numpy as np


class ParameterDistribution(ABC):
    """Abstract base class for parameter distributions."""

    @abstractmethod
    def sample(self, rng: np.random.Generator) -> float:
        """Draw a single sample from the distribution."""
        ...

    @abstractmethod
    def mean(self) -> float:
        """Expected value."""
        ...

    @abstractmethod
    def std(self) -> float:
        """Standard deviation."""
        ...

    @property
    def cv(self) -> float:
        """Coefficient of variation (std/mean). Dimensionless."""
        if self.mean() == 0:
            return float('inf')
        return self.std() / abs(self.mean())

    def sample_n(self, n: int, rng: np.random.Generator) -> np.ndarray:
        """Draw n samples."""
        return np.array([self.sample(rng) for _ in range(n)])

    def is_degenerate(self) -> bool:
        """True if this is a fixed value (zero variance)."""
        return self.std() == 0.0


@dataclass(frozen=True)
class LogNormalDistribution(ParameterDistribution):
    """
    Log-normal distribution parameterised by mean and coefficient of variation.

    If X ~ LogNormal(mu, sigma):
        E[X]   = exp(mu + sigma²/2) = mean
        Var[X] = (exp(sigma²) - 1) * exp(2*mu + sigma²)
        CV     = sqrt(exp(sigma²) - 1)

    Parameterised by (mean, cv) for interpretability:
        sigma = sqrt(log(1 + cv²))
        mu    = log(mean) - sigma²/2
    """
    _mean: float
    _cv: float

    def __post_init__(self) -> None:
        if self._mean <= 0:
            raise ValueError(f"LogNormal mean must be positive, got {self._mean}")
        if self._cv <= 0:
            raise ValueError(f"LogNormal CV must be positive, got {self._cv}")

    @classmethod
    def from_mean_cv(cls, mean: float, cv: float) -> "LogNormalDistribution":
        return cls(_mean=mean, _cv=cv)

    @classmethod
    def from_confidence_tier(
        cls, mean: float, tier: str, additional_cv: float = 0.0
    ) -> "LogNormalDistribution":
        """
        Construct from a parameter value and its confidence tier.

        Args:
            mean: The point estimate (e.g. from ParameterStore.best())
            tier: 'gold', 'silver', or 'bronze'
            additional_cv: Extra uncertainty to combine with the tier's base
                CV in quadrature (sqrt(tier_cv**2 + additional_cv**2)) --
                e.g. `ComposedEffect.genetic_cv` from `biosim_genotype`,
                per `composer.py`'s documented contract that "the genetic
                uncertainty CV is combined in quadrature with the base
                parameter uncertainty." Defaults to 0.0 (no change), which
                reduces to the prior tier-only behaviour.

        The CV values below are a documented expert-judgment default, not
        a directly-cited empirical statistic -- see this module's
        docstring for why no single-number citation exists for this
        specific quantity.
        """
        cv_map = {"gold": 0.10, "silver": 0.30, "bronze": 0.60}
        if tier not in cv_map:
            raise ValueError(f"Unknown tier '{tier}'. Must be one of {list(cv_map)}")
        combined_cv = math.sqrt(cv_map[tier] ** 2 + additional_cv ** 2)
        return cls(_mean=mean, _cv=combined_cv)

    @property
    def _sigma(self) -> float:
        return math.sqrt(math.log(1.0 + self._cv ** 2))

    @property
    def _mu(self) -> float:
        return math.log(self._mean) - 0.5 * self._sigma ** 2

    @property
    def mu(self) -> float:
        """
        Public accessor for the underlying log-space mu parameter --
        added when sobol_sensitivity.py became the first consumer of
        this value from outside this module (previously reached into
        the underscore-prefixed _mu directly, a real encapsulation gap
        for a second real caller). Same value as _mu; kept both rather
        than renaming _mu, since _mu/_sigma are still used internally
        by sample()/percentile() and renaming risks nothing but adds
        churn for no benefit.
        """
        return self._mu

    @property
    def sigma(self) -> float:
        """Public accessor for the underlying log-space sigma parameter -- see mu's docstring."""
        return self._sigma

    def sample(self, rng: np.random.Generator) -> float:
        return float(rng.lognormal(mean=self._mu, sigma=self._sigma))

    def mean(self) -> float:
        return self._mean

    def std(self) -> float:
        return self._mean * self._cv

    def percentile(self, p: float) -> float:
        """Return the p-th percentile (0-100)."""
        from scipy.stats import lognorm
        return float(lognorm.ppf(p / 100.0, s=self._sigma, scale=math.exp(self._mu)))

    def fold_range_95(self) -> float:
        """
        95% confidence fold range (p97.5 / p2.5).
        A value of 4.0 means the 95% CI spans a 4-fold range.
        """
        return self.percentile(97.5) / self.percentile(2.5)


@dataclass(frozen=True)
class NormalDistribution(ParameterDistribution):
    """
    Normal distribution. Use for parameters that can take negative values
    (e.g. fold changes, log-transformed quantities).

    For kinetic parameters (always positive), prefer LogNormalDistribution.
    """
    _mean: float
    _std: float

    def __post_init__(self) -> None:
        if self._std < 0:
            raise ValueError(f"Normal std must be non-negative, got {self._std}")

    @classmethod
    def from_mean_std(cls, mean: float, std: float) -> "NormalDistribution":
        return cls(_mean=mean, _std=std)

    def sample(self, rng: np.random.Generator) -> float:
        return float(rng.normal(loc=self._mean, scale=self._std))

    def mean(self) -> float:
        return self._mean

    def std(self) -> float:
        return self._std


@dataclass(frozen=True)
class UniformDistribution(ParameterDistribution):
    """
    Uniform distribution over [lower, upper].
    Use when only bounds are known (maximum entropy for bounded parameters).
    """
    lower: float
    upper: float

    def __post_init__(self) -> None:
        if self.lower >= self.upper:
            raise ValueError(f"Uniform lower ({self.lower}) must be < upper ({self.upper})")

    def sample(self, rng: np.random.Generator) -> float:
        return float(rng.uniform(self.lower, self.upper))

    def mean(self) -> float:
        return (self.lower + self.upper) / 2.0

    def std(self) -> float:
        return (self.upper - self.lower) / math.sqrt(12.0)


@dataclass(frozen=True)
class FixedValue(ParameterDistribution):
    """
    Degenerate distribution — always returns the same value.
    Use for parameters with zero uncertainty (e.g. stoichiometric ratios,
    or parameters whose uncertainty contribution is negligible).
    """
    value: float

    def sample(self, rng: np.random.Generator) -> float:
        return self.value

    def mean(self) -> float:
        return self.value

    def std(self) -> float:
        return 0.0

    def is_degenerate(self) -> bool:
        return True
