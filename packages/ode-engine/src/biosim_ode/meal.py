"""
Meal input — models glucose appearance in blood after eating.

This is the forcing function for the 'eating X' use case.
It converts a food item (specified in grams of carbohydrate) into a
time-varying glucose appearance rate Ra(t) that drives the ODE system.

Two absorption shapes are available:

1. Two-compartment ("Bateman") curve -- the DEFAULT for the meal presets
   below and for scenario meals (biosim_api MealEvent) with no
   glycaemic_index and no kabs (since 2026-11-15). Set via
   `absorption_model=GutAbsorptionModel...`; see gut_absorption.py.
   Ra is zero when you eat, peaks later (~31 min with the default rate
   constants) and integrates to D_eff.

2. Single exponential, used when `absorption_model` is None:
    Ra(t) = D_eff * kabs * exp(-kabs * (t - t_meal))  for t >= t_meal
    Ra(t) = 0                                          for t < t_meal
   This peaks the instant you eat, which is not physiological. It is kept
   because a caller can still ask for it by setting kabs explicitly (e.g.
   biosim_seeds gi_meal_scenario(), via gi_adjusted_kabs()), and because
   the acarbose/miglitol effect was calibrated on it (SimulationRunner
   applies that effect to a two-compartment default meal by slowing both
   of its rate constants, and warns that the timing does not match exactly).

Both shapes share:
    D_eff = effective glucose reaching systemic blood [mM]
          = CHO_mmol * (1 - hepatic_extraction) / Vg_L
    CHO_mmol = carbohydrate_g * 1000 / 180 g/mol
    Vg_L = body_weight_kg * Vg_L_per_kg

Time is in minutes (ODE engine internal units).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

from .gut_absorption import GutAbsorptionModel


# Default physiological parameters
DEFAULT_BODY_WEIGHT_KG = 70.0
# Effective glucose distribution volume. Was an uncited 0.286 L/kg; now the
# same literature-anchored value used for GlucoseInsulinPathway.bergman_ivgtt_y0()
# (see BERGMAN_IVGTT_VG_L_PER_KG in pathway.py) -- Mora-Gonzalez D, Moreno-Cabanas A,
# Alvarez-Jimenez L, Morales-Palomo F, Ortega JF, Mora-Rodriguez R. "Glucose
# volume of distribution affects insulin sensitivity measured by intravenous
# glucose tolerance test." Scand J Med Sci Sports. 2024. PMID:38389141
# (verified live via PubMedFetcher): healthy-subject VdGLUCOSE = 230 +/- 49 mL/kg.
# Verified before changing: this only rescales Ra(t)'s AMPLITUDE (total dose),
# not its timing shape (peak_glucose/glucose_at_120min stay within their
# NORMAL/25% tolerance bands; time_to_peak_glucose_min is unaffected -- see
# CLAUDE.md KNOWN BUGS FIXED #19 for the before/after numbers).
DEFAULT_VG_L_PER_KG = 0.230
DEFAULT_KABS_PER_MIN = 0.025        # Absorption rate constant (t1/2 ~28 min)

# Was an uncited 0.30 ("30% haircut on every meal's dose amplitude, no
# citation anywhere in the file" -- flagged but not fixed in an earlier
# pass this session, see CLAUDE.md's now-resolved "NEXT" entry). Fixed
# using the actual primary source, verified live via PubMedFetcher
# (fetch_by_pubmed_id('651653'), not a WebSearch summary -- an earlier
# WebSearch summary had misattributed this exact PMID/title to a
# different author group, "Ferrannini, Bjorkman, Wahren, DeFronzo,"
# caught only by insisting on the live fetch):
#   Radziuk J, McDonald TJ, Rubenstein D, Dupre J. "Initial splanchnic
#   extraction of ingested glucose in normal man." Metabolism.
#   1978;27:657-669. PMID:651653.
# The abstract directly states the fraction needed, using a dual-tracer
# method that distinguishes oral-origin glucose from endogenous
# production (so this is NOT confounded by suppressed hepatic output,
# the concern that blocked fixing this earlier): of 45+/-1g and 96+/-5g
# oral glucose loads, 40+/-3g and 86+/-4g respectively reached the
# systemic pool, "representing 90% +/-4% of the administered glucose."
# Both doses (which bracket this platform's typical 75g meal) agree to
# within 1 percentage point (88.9% vs 89.6%), supporting extrapolation
# across that range. Extraction = 1 - 0.90 = 0.10 (SEM range ~0.06-0.14).
# Checked against every registered clinical benchmark (ModelValidator.
# validate_all()) before accepting, not just one point, per this
# session's established practice for exactly this class of fix (see
# meal.py's Vg_L_per_kg fix history below) -- verdict stays PASS on all
# four datasets with this change; see CLAUDE.md KNOWN BUGS FIXED for the
# full before/after numbers, including the one already-known gap
# (semaglutide's WIDE-tier, informational-only overshoot) that widens
# further and is documented honestly rather than used to justify
# reverting a directly primary-source-verified value.
DEFAULT_HEPATIC_EXTRACTION = 0.10
GLUCOSE_MOLECULAR_WEIGHT = 180.16   # g/mol


@dataclass
class MealInput:
    """
    Specifies a single food intake event.

    Args:
        carbs_g: Grams of carbohydrate (determines glucose dose)
        time_min: Time of meal onset in minutes from simulation start
        body_weight_kg: User body weight for dose normalization
        kabs: Gut absorption rate constant [min⁻¹]
        hepatic_extraction: Fraction of absorbed glucose removed by liver
        Vg_L_per_kg: Glucose distribution volume [L/kg]
        label: Optional human-readable label (e.g. "Breakfast", "Caffeine + croissant")
        absorption_model: Optional GutAbsorptionModel (bi-exponential, Ra(0)=0,
            realistic delayed peak -- see gut_absorption.py). When set, this
            REPLACES the monoexponential kabs formula below for Ra(t); kabs is
            then unused. Defaults to None (today's monoexponential behaviour,
            unchanged at this low level). The presets below and scenario
            meals (MealEvent.to_meal_input) set it to
            GutAbsorptionModel.default_meal() unless a GI or kabs is given.
    """
    carbs_g: float
    time_min: float = 0.0
    body_weight_kg: float = DEFAULT_BODY_WEIGHT_KG
    kabs: float = DEFAULT_KABS_PER_MIN
    hepatic_extraction: float = DEFAULT_HEPATIC_EXTRACTION
    Vg_L_per_kg: float = DEFAULT_VG_L_PER_KG
    label: str = ""
    absorption_model: Optional[GutAbsorptionModel] = None
    # Grams of protein co-ingested with this meal -- default 0.0 (backward
    # compatible, no existing meal is affected). Drives
    # NutrientSecretionPathway's protein-triggered endogenous insulin
    # secretion (see nutrient_secretion_pathway.py) when > 0; carried here
    # rather than only on MealEvent so the ODE-level MealInput itself is
    # self-contained, matching carbs_g's own placement.
    protein_g: float = 0.0

    def __post_init__(self) -> None:
        # math.isnan() checks added alongside the ordinary range checks below
        # -- found during the meal-absorption mass-conservation fix's
        # security review: `carbs_g < 0` and `body_weight_kg <= 0` are both
        # silently False for NaN (every NaN comparison is False), so a NaN
        # input passed both checks and reached MealAbsorptionPathway's y0
        # unvalidated, then hit a separate clamp() defect (also fixed this
        # pass) that mapped NaN to +Infinity instead of a safe floor. kabs's
        # `not (0.0 < kabs < 1.0)` and hepatic_extraction's chained
        # comparison already reject NaN correctly (both sides of a chained
        # comparison are False for NaN, so `not(...)` is True) -- that was
        # incidental to the chained form, not a deliberate NaN guard, but
        # it means only carbs_g/body_weight_kg needed this explicit check.
        if math.isnan(self.carbs_g) or self.carbs_g < 0:
            raise ValueError(f"carbs_g cannot be negative or NaN, got {self.carbs_g}")
        if math.isnan(self.body_weight_kg) or self.body_weight_kg <= 0:
            raise ValueError(f"body_weight_kg must be positive and finite, got {self.body_weight_kg}")
        if not (0.0 < self.kabs < 1.0):
            raise ValueError(f"kabs must be between 0 and 1 min⁻¹, got {self.kabs}")
        if not (0.0 <= self.hepatic_extraction < 1.0):
            raise ValueError(f"hepatic_extraction must be in [0, 1), got {self.hepatic_extraction}")
        if math.isnan(self.protein_g) or self.protein_g < 0:
            raise ValueError(f"protein_g cannot be negative or NaN, got {self.protein_g}")

    @property
    def carbs_mmol(self) -> float:
        """Total carbohydrate in millimoles."""
        return self.carbs_g * 1000.0 / GLUCOSE_MOLECULAR_WEIGHT

    @property
    def Vg_L(self) -> float:
        """Glucose distribution volume in litres."""
        return self.body_weight_kg * self.Vg_L_per_kg

    @property
    def D_eff_mM(self) -> float:
        """
        Effective glucose concentration reaching systemic blood [mM].
        Accounts for hepatic first-pass extraction.
        """
        return self.carbs_mmol * (1.0 - self.hepatic_extraction) / self.Vg_L

    def Ra(self, t_min: float) -> float:
        """
        Glucose appearance rate at time t [mM/min].

        Args:
            t_min: Time in minutes from simulation start

        Returns:
            Rate of glucose appearance in blood [mM/min]
        """
        dt = t_min - self.time_min
        if dt < 0:
            return 0.0
        if self.absorption_model is not None:
            # Bi-exponential model is mass-conserving (integral over infinite
            # time == its dose argument, same invariant as D_eff_mM below),
            # so passing D_eff_mM through unchanged delivers the identical
            # total glucose dose as the monoexponential path -- only the
            # time-shape of Ra(t) changes (Ra(0)=0, delayed realistic peak,
            # vs. an instantaneous peak at t=0).
            return self.absorption_model.Ra(dt, self.D_eff_mM)
        return self.D_eff_mM * self.kabs * math.exp(-self.kabs * dt)

    def peak_Ra(self) -> float:
        """Peak of the SINGLE-EXPONENTIAL Ra [mM/min], at t = time_min.
        Ignores absorption_model."""
        return self.D_eff_mM * self.kabs

    def time_to_half_peak(self) -> float:
        """Single-exponential half-time [min]. Ignores absorption_model."""
        return math.log(2.0) / self.kabs

    def total_glucose_appearance(self) -> float:
        """
        Total glucose appearing in blood over infinite time [mM * min → mM].
        ∫₀^∞ Ra(t) dt = D_eff_mM
        """
        return self.D_eff_mM


def total_Ra(meals: list[MealInput], t_min: float) -> float:
    """
    Sum glucose appearance rate from multiple meals at time t.

    Args:
        meals: List of MealInput objects (can be empty)
        t_min: Current time in minutes

    Returns:
        Total glucose appearance rate [mM/min]
    """
    return sum(meal.Ra(t_min) for meal in meals)


# ---------------------------------------------------------------------------
# Common meal presets for quick simulation setup
# ---------------------------------------------------------------------------

def high_sugar_meal(
    time_min: float = 0.0,
    body_weight_kg: float = DEFAULT_BODY_WEIGHT_KG,
) -> MealInput:
    """High-sugar meal: ~75g carbohydrate (equivalent to standard OGTT)."""
    return MealInput(carbs_g=75.0, time_min=time_min,
                     body_weight_kg=body_weight_kg, label="High-sugar meal (75g CHO)",
                     absorption_model=GutAbsorptionModel.default_meal())


def standard_meal(
    time_min: float = 0.0,
    body_weight_kg: float = DEFAULT_BODY_WEIGHT_KG,
) -> MealInput:
    """Standard mixed meal: ~50g carbohydrate."""
    return MealInput(carbs_g=50.0, time_min=time_min,
                     body_weight_kg=body_weight_kg, label="Standard meal (50g CHO)",
                     absorption_model=GutAbsorptionModel.default_meal())


def light_snack(
    time_min: float = 0.0,
    body_weight_kg: float = DEFAULT_BODY_WEIGHT_KG,
) -> MealInput:
    """Light snack: ~20g carbohydrate."""
    return MealInput(carbs_g=20.0, time_min=time_min,
                     body_weight_kg=body_weight_kg, label="Light snack (20g CHO)",
                     absorption_model=GutAbsorptionModel.default_meal())
