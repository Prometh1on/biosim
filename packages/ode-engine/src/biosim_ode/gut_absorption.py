# SAVE TO: packages\ode-engine\src\biosim_ode\gut_absorption.py
"""
Gut absorption model — bi-exponential glucose appearance rate Ra(t).

Two first-order compartments in series (stomach -> gut -> plasma), which
gives the standard two-compartment first-order absorption ("Bateman")
curve below. This is a deliberately simple shape: gastric emptying is
proportional to what is left in the stomach, with no maximum rate.

It is NOT the Lehmann-Deutsch model (Lehmann & Deutsch 1992, J Biomed Eng,
PMID:1588781), which earlier versions of this docstring credited. That
paper's abstract (fetched from PubMed 2026-11-12) describes gastric
emptying parameterised by a "maximal rate of gastric emptying", i.e. a
capacity-limited process, which a first-order curve cannot represent.
This is an ABSTRACT-ONLY check (re-checked 2026-11-14): the paper has no
PubMed Central full text, so its gastric-emptying equations were not read.
Dalla Man et al. 2007 (PMID:17926672) likewise uses a nonlinear,
multi-compartment gut model. Both are cited below only as background.

Bi-exponential Ra(t):
    Ra(t) = D * k_a * k_e / (k_a - k_e) * (exp(-k_e*t) - exp(-k_a*t))

where:
    D   = dose [mM·kg⁻¹]
    k_e = gastric emptying rate [min⁻¹]
    k_a = intestinal absorption rate [min⁻¹]
    t   = time since meal [min]

Properties:
    - Ra(0) = 0 (absorption starts at zero)
    - Ra peaks at t_peak = ln(k_a/k_e) / (k_a - k_e)
    - Total glucose dose = integral of Ra dt = D (mass conserved)

GI-adjusted emptying (UNVALIDATED PLACEHOLDER -- only GI=100 is calibrated):
    High-GI foods → faster gastric emptying → earlier, sharper peak
    Low-GI foods  → slower emptying → later, flatter peak
    The slope linking GI to emptying rate has no literature calibration;
    see GutAbsorptionModel.for_gi().

Background literature (models this code does NOT implement):
    Lehmann & Deutsch 1992, J Biomed Eng. PMID:1588781
    Dalla Man et al. 2007, IEEE Trans Biomed Eng. PMID:17926672
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional


# Shared wording for every user-facing output that relies on a GI-based
# absorption-timing adjustment (GutAbsorptionModel.for_gi(), and
# biosim_seeds.gi_meals.gi_adjusted_kabs()).
GI_SLOPE_PLACEHOLDER_NOTE = (
    "The glycaemic-index (GI) timing adjustment is an unvalidated "
    "placeholder: only GI=100 (pure glucose) is calibrated against clinical "
    "data, so meal timing for other GI values is illustrative only."
)


# Where the default meal's absorption rate constants come from. No paper
# was found that gives two-compartment (Bateman) rate constants for an
# ordinary mixed meal, so the default meal reuses the OGTT preset's rate
# constants (GutAbsorptionModel.ogtt_calibrated_glucose()). This is a
# modelling choice, not a measurement:
#   - those constants were fitted to a 75 g glucose DRINK, not a mixed meal;
#   - they were fitted to three of the OGTT glucose benchmark points, so any
#     benchmark that uses a default meal shares that partial circularity.
# verification_status uses the parameter database's VerificationStatus
# vocabulary ("UNVERIFIED" = no source cited for this value).
DEFAULT_MEAL_ABSORPTION_PROVENANCE: dict[str, str] = {
    "label": "MODEL ASSUMPTION (no mixed-meal source)",
    "verification_status": "UNVERIFIED",
    "source": "Reuses GutAbsorptionModel.ogtt_calibrated_glucose() "
              "(k_emptying=0.0249, k_absorption=0.0416 per min)",
    "reason": "No cited source gives two-compartment rate constants for a "
              "mixed meal; the OGTT preset was fitted to a 75 g glucose drink "
              "against this model's own OGTT glucose benchmark.",
}


@dataclass(frozen=True)
class GutAbsorptionModel:
    """
    Bi-exponential gut absorption model parameters.

    k_emptying:   Gastric emptying rate [min⁻¹] (controls when food leaves stomach)
    k_absorption: Intestinal absorption rate [min⁻¹] (controls uptake from gut)

    Constraint: k_absorption > k_emptying (absorption faster than emptying)
    """
    k_emptying:   float = 0.055    # min⁻¹  (~12.6 min half-life)
    k_absorption: float = 0.180    # min⁻¹  (~3.9 min half-life in gut lumen)

    def __post_init__(self) -> None:
        if self.k_emptying <= 0:
            raise ValueError(f"k_emptying must be positive, got {self.k_emptying}")
        if self.k_absorption <= self.k_emptying:
            raise ValueError(
                f"k_absorption ({self.k_absorption}) must be > k_emptying ({self.k_emptying})"
            )

    def Ra(self, t_min: float, dose_mM_per_kg: float) -> float:
        """
        Glucose appearance rate at time t_min after meal ingestion.

        Args:
            t_min:          Time since meal [minutes]
            dose_mM_per_kg: Glucose dose [mM per kg body weight]

        Returns:
            Ra(t) in mM/min
        """
        if t_min < 0 or dose_mM_per_kg <= 0:
            return 0.0
        ka, ke = self.k_absorption, self.k_emptying
        coeff  = dose_mM_per_kg * ka * ke / (ka - ke)
        return coeff * (math.exp(-ke * t_min) - math.exp(-ka * t_min))

    @property
    def peak_time_min(self) -> float:
        """Time of peak absorption rate."""
        ka, ke = self.k_absorption, self.k_emptying
        return math.log(ka / ke) / (ka - ke)

    @property
    def time_to_90pct_absorbed_min(self) -> float:
        """Approximate time when 90% of dose has been absorbed."""
        # Numerical estimate: 5 × (1/k_emptying)
        return 5.0 / self.k_emptying

    @classmethod
    def for_gi(cls, glycaemic_index: float) -> "GutAbsorptionModel":
        """
        Return a GutAbsorptionModel calibrated to a food's glycaemic index.

        Anchored at GI=100 to `ogtt_calibrated_glucose()`
        (k_emptying=0.0249, k_absorption=0.0416) -- the one point on this
        curve with real clinical validation. This fixes a confirmed,
        citation-independent internal-consistency bug: the old formula
        returned k_emptying=0.090/k_absorption=0.315 at GI=100, a
        3.6x/7.6x discontinuity against the immediately-adjacent validated
        preset that `MealEvent.to_meal_input()` routes to at GI>=99.5 (see
        `test_for_gi_100_matches_ogtt_calibrated_glucose_exactly` in
        `test_absorption_and_comparison.py`, formerly
        `test_ogtt_calibrated_glucose_much_slower_than_for_gi_100` before
        this fix). Fixing this anchor is independent of the open question
        below -- the discontinuity was wrong regardless of how GI relates
        to timing.

        The GI-to-k_emptying SLOPE (how much lower GI should slow
        absorption) is still an explicit, un-validated placeholder, NOT a
        literature-derived calibration. The proportional shape (k_emptying
        at GI=30 is 1/3 of its value at GI=100) is preserved unchanged from
        the pre-existing, never-cited formula, rather than re-derived from
        either study below -- extrapolating a 2-point, single-population
        dataset across the full GI=30-100 range would be curve-fitting
        dressed up as derivation, which this project explicitly avoids
        (see CLAUDE.md's repeated rejections of exactly that pattern, e.g.
        the KG_secretion≈200 and joint-fit p1≈0.008 leads). `SimulationRunner`
        emits a runtime warning whenever this path is used (see
        `runner.py`) precisely because this slope question remains open:
        Cocate et al. 2011 (PMID:21208446) found no GI-driven timing
        difference at all between high/low GI meals (only peak height
        differed), while Louie et al. 2012/2013 (PMID:23344248) found low
        GI peaked ~17 min later, in a different population. Two studies,
        two GI points each, opposite conclusions -- not enough to referee
        which is right, let alone recalibrate the slope against either
        alone. See CLAUDE.md VALIDATION STATUS for the full literature
        search.

        Args:
            glycaemic_index: Food GI value (0-100, glucose reference = 100)
        """
        gi = max(10.0, min(100.0, glycaemic_index))
        # Anchor: k_emptying/k_absorption at GI=100 match
        # ogtt_calibrated_glucose() exactly (validated). Shape: k_emptying
        # at GI=30 is 1/3 of the GI=100 value, preserving the original
        # formula's proportional range (0.030/0.090 = 1/3) -- see docstring
        # for why this ratio isn't re-derived from the conflicting studies.
        k_e_at_100 = 0.0249
        k_e_at_30 = k_e_at_100 / 3.0
        k_e = k_e_at_30 + (gi - 30.0) / 70.0 * (k_e_at_100 - k_e_at_30)
        k_a_to_k_e_ratio = 0.0416 / 0.0249  # from ogtt_calibrated_glucose()
        k_a = k_e * k_a_to_k_e_ratio
        return cls(k_emptying=k_e, k_absorption=k_a)

    @classmethod
    def liquid_glucose(cls) -> "GutAbsorptionModel":
        """Pure glucose drink — fastest absorption (reference OGTT)."""
        return cls(k_emptying=0.082, k_absorption=0.250)

    @classmethod
    def ogtt_calibrated_glucose(cls) -> "GutAbsorptionModel":
        """
        Pure glucose (75g OGTT), rate constants fit against the 3 glucose-
        timing benchmarks in `ogtt_75g_healthy_adults()`
        (`peak_glucose`, `time_to_peak_glucose_min`, `glucose_at_120min`) via
        least-squares (`scipy.optimize.minimize`, Nelder-Mead) over the full
        composed ODE system -- NOT hand-picked. Unlike `liquid_glucose()`
        above (peak absorption at ~6.6 min, chosen for correct *qualitative*
        ordering relative to `mixed_meal()`/`high_fat_meal()` but never
        checked against real OGTT timing), this preset reproduces
        `peak_glucose`/`time_to_peak_glucose_min` within ~4%/~2% of their
        clinical reference values (see CLAUDE.md VALIDATION STATUS for the
        fit methodology and full before/after numbers).

        Deliberately NOT fit against `peak_insulin`/`time_to_peak_insulin_min`
        -- those stayed just as far off target (~-40% to -45%) across every
        absorption-shape tried, including this one, which is strong evidence
        that insulin peak/timing is governed by the insulin secretion
        kinetics in `pathway.py` (F_secretion_basal, ke, KG_secretion), not
        by gut absorption shape -- a distinct, not-yet-investigated
        calibration question that this preset does not attempt to fix.
        """
        return cls(k_emptying=0.0249, k_absorption=0.0416)

    @classmethod
    def default_meal(cls) -> "GutAbsorptionModel":
        """
        Absorption curve for a meal with no glycaemic index and no kabs.

        MODEL ASSUMPTION: same rate constants as ogtt_calibrated_glucose()
        (absorption peaks ~31 min after eating). See
        DEFAULT_MEAL_ABSORPTION_PROVENANCE for why, and what that costs.
        """
        return cls.ogtt_calibrated_glucose()

    @classmethod
    def mixed_meal(cls) -> "GutAbsorptionModel":
        """Mixed solid meal — moderate absorption."""
        return cls(k_emptying=0.050, k_absorption=0.160)

    @classmethod
    def high_fat_meal(cls) -> "GutAbsorptionModel":
        """High fat/protein meal — slowest gastric emptying."""
        return cls(k_emptying=0.028, k_absorption=0.090)


def Ra_biexp(
    t_min: float,
    dose_mM_per_kg: float,
    model: Optional[GutAbsorptionModel] = None,
) -> float:
    """
    Convenience function: compute Ra(t) using bi-exponential model.

    Args:
        t_min:          Time since meal [minutes]
        dose_mM_per_kg: Glucose dose [mM/kg]
        model:          GutAbsorptionModel (default: mixed meal)

    Returns:
        Ra(t) [mM/min]
    """
    model = model or GutAbsorptionModel.mixed_meal()
    return model.Ra(t_min, dose_mM_per_kg)


def Ra_curve(
    duration_min: float,
    dose_mM_per_kg: float,
    model: Optional[GutAbsorptionModel] = None,
    n_points: int = 200,
) -> tuple[list[float], list[float]]:
    """
    Compute the full Ra(t) curve over a simulation window.

    Returns:
        (times_min, Ra_values) both as lists of length n_points
    """
    import numpy as np
    model = model or GutAbsorptionModel.mixed_meal()
    times = list(np.linspace(0, duration_min, n_points))
    values = [model.Ra(t, dose_mM_per_kg) for t in times]
    return times, values
