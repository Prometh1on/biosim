# SAVE TO: packages\ode-engine\src\biosim_ode\nutrient_secretion_pathway.py
"""
Nutrient (dietary protein) secretion pathway -- protein co-ingested with a
meal directly stimulates endogenous insulin secretion, independent of the
glucose-driven secretion this platform already models. Found via a fresh
Track 1 scoping fan-out (MASTER_PLAN row 1.28) after every previously-named
Track 1 candidate was done or blocked -- this platform's food model was
carbohydrate-only (`MealEvent.carbs_g`), with no representation of protein
at all.

Architecturally the SAME shape as `exogenous_insulin_pathway.py`'s own
ExogenousInsulinPathway: STATELESS (zero species -- the whole secretion
curve has a closed-form analytic solution, no ODE integration state
needed) and TRANSIENT (exactly zero before/without protein, decaying back
to exactly zero afterward -- never perturbs the analytic basal steady
state, `GlucoseInsulinPathway().y0`, for any protein dose). Distinct
mechanism from exogenous_insulin_delta (this is ENDOGENOUS beta-cell
secretion triggered by amino acids, not an injected external source), so
it is kept as its own rule-11 signal key rather than merged with either
exogenous_insulin_delta or glucose_independent_secretion (rule 8's
constant floor for sulfonylureas).

Literature, both live-verified via PubMedFetcher before implementation:

    Nuttall FQ, Mooradian AD, Gannon MC, Billington C, Krezowski P. "Effect
    of protein ingestion on the glucose and insulin response to a
    standardized oral glucose load." Diabetes Care. 1984;7(5):465-70.
    PMID:6389060. Type II diabetic subjects, 50g protein / 50g glucose /
    both, as a single meal. Reports a real SYNERGISTIC effect in T2D
    (combined insulin AUC 247+/-33 microU.h/ml vs. an additive 97+83=180
    from the two arms alone) -- NOT the regime this implementation
    targets (see below).

    Krezowski PA, Nuttall FQ, Gannon MC, Bartosh NH. "The effect of
    protein ingestion on the metabolic response to oral glucose in normal
    individuals." Am J Clin Nutr. 1986;44(6):847-56. PMID:3538843. n=8
    NORMAL (non-diabetic) subjects, same 50g/50g design. Confirmed
    verbatim: "Protein stimulated only a modest insulin rise and the area
    above fasting baseline was only 28% of that after glucose... The sum
    of the serum insulin area following protein ingestion and that
    following glucose ingestion was 100.4% of the combination meal... the
    effect on insulin secretion is not synergistic." Also: "When various
    amounts of protein were given with 50g glucose, the insulin area
    response was essentially first order" (Nuttall 1984's own finding,
    corroborated) -- supports a linear-in-grams dose law.

    WHY THE NORMAL-SUBJECT (Krezowski), NOT THE T2D (Nuttall), RESULT
    GOVERNS THIS MECHANISM'S SHAPE: this platform's core ODE (Bergman
    1981) is calibrated on a normal/healthy cohort, matching the same
    "which regime does this function actually represent" resolution this
    project already used for metformin's own p3/p1 removal (see
    drug_interactions.py's metformin_modifiers() docstring). Krezowski's
    additive (not synergistic) result in normals is exactly the shape a
    SIMPLE ADDITIVE rule-11 signal can represent; Nuttall's T2D synergy
    would need a genuinely different (multiplicative, glucose-gated)
    mechanism this implementation does NOT attempt -- disclosed as an
    out-of-scope divergence for T2D scenarios, not silently ignored.

    NOT IMPLEMENTED, DISCLOSED AS A REAL LIMITATION: both papers report a
    protein-driven GLUCAGON rise (Krezowski: "Addition of glucose to a
    protein meal results in a delayed rise in AAN and glucagon"). Row 1.1
    (glucagon's own fasting/basal EGP effect) is permanently blocked --
    no citable dose-response gain exists at the small perturbations this
    platform would need. An insulin-only implementation will therefore
    OVER-PREDICT glucose lowering for protein-heavy or protein-alone
    meals, since the real net effect is insulin's glucose-lowering action
    partially offset by glucagon's glucose-raising one. This mechanism is
    scoped to protein CO-INGESTED WITH CARBOHYDRATE (the calibration
    target, and the clinically typical use case); a protein-only meal's
    predicted glucose trajectory should not be trusted as validated.

Calibration, disclosed as a genuine engineering back-solve (the 28%
insulin-AUC RATIO is literature-quoted; the specific pM/g scale needed to
reproduce that ratio through THIS platform's own composed ODE is not) --
see NUTRIENT_SECRETION_PM_PER_G's own comment below for the exact
two-step numeric procedure and the resulting value.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Optional

from .base_pathway import BasePathway, PathwayDescriptor, PathwayMaturity

if TYPE_CHECKING:
    from .coupling import CompositePathway
    from .pathway import GlucoseInsulinPathway


# Disclosed, UNANCHORED timing assumption -- neither Krezowski nor Nuttall
# reports a time-to-peak or any timing curve for protein-driven insulin
# secretion, only a total AUC over several hours. A single monoexponential
# decay for the secretion RATE (mirroring MealInput.Ra()'s own shape,
# D_eff*kabs*exp(-kabs*t)) has the property that its own total time-integral
# (S1_0 below) is EXACTLY independent of the rate constant chosen -- so this
# constant affects only the signal's peak/timing shape, never the
# literature-calibrated total magnitude. Reused directly from this
# platform's own default glucose absorption rate constant
# (MealInput.DEFAULT_KABS_PER_MIN) as the least-arbitrary available choice,
# not because protein absorption is claimed to match glucose's own rate.
NUTRIENT_SECRETION_KA_PER_MIN = 0.025

# Engineering back-solve (NOT a directly-cited number -- see this module's
# own docstring for the real, literature-quoted ratio this reproduces).
# Two-step numeric procedure, run directly against the real composed ODE:
#   1. Measured this platform's OWN 50g-carbohydrate-alone meal insulin
#      AUC-above-baseline over a 240-minute window (matching Krezowski's
#      own ~4h protocol): 28567.0 pM.min.
#   2. Target for a 50g-protein-alone meal (Krezowski's own 28% ratio):
#      0.28 * 28567.0 = 7998.8 pM.min.
#   3. Ran the ACTUAL composed ODE with a trial NUTRIENT_SECRETION_PM_PER_G,
#      measured the resulting insulin AUC-above-baseline, and rescaled
#      linearly (the perturbation is modest enough that the composed
#      system's response is very close to linear in this forcing term --
#      confirmed empirically: two trial runs at different scales rescaled
#      to the identical implied constant to 3 significant figures) to
#      reach the final value below, which reproduces the target to
#      within 0.02%.
#
# Additivity cross-check (Krezowski's OWN "not synergistic" finding,
# 100.4% additive): composing a 50g-carbs+50g-protein meal in this
# platform's own ODE gives combined_AUC/(glucose_AUC+protein_AUC) =
# 0.977, not exactly 1.0 -- a real, expected consequence of this being a
# genuinely coupled nonlinear system (protein's extra insulin blunts the
# glucose excursion, which in turn reduces glucose-driven secretion in
# the COMBINED run relative to the sum of the two "alone" runs), not a
# fabricated synergy or a calibration error. 0.977 is close to
# Krezowski's own 1.004 (both within ~3% of exact additivity), which is
# itself just a point estimate from n=8 subjects -- not forced closer,
# disclosed honestly instead.
NUTRIENT_SECRETION_PM_PER_G = 14.93

# Engineering plausibility ceiling (mandatory security-reviewer finding,
# HIGH -- fixed before finalizing this milestone), NOT a literature-
# derived clinical maximum -- same category and reasoning as
# ExogenousInsulinPathway's own MAX_PLAUSIBLE_SC_INSULIN_DOSE_UNITS.
# NUTRIENT_SECRETION_PM_PER_G has no saturating/Hill term (a single
# literature data point, Krezowski's own 50g anchor, cannot support
# fitting one) -- the dose-response is purely linear, so nothing in the
# model otherwise prevents an arbitrarily large, REST-legal protein_g
# from producing an implausible hyperinsulinemic spike. Verified
# directly against the real composed ODE before choosing this value:
# a 300g-protein-ALONE meal (the prior REST ceiling) pushed
# peak_insulin_pM to ~1136-1262 pM, roughly 3.2-3.6x this platform's
# own validated OGTT peak_insulin reference (350 pM, CLAUDE.md
# VALIDATION STATUS) -- headroom above any realistic single-meal
# protein intake (even a very large steak-and-shake meal rarely
# exceeds ~100-150g), not a hard physiological limit.
MAX_PLAUSIBLE_PROTEIN_G_PER_MEAL = 150.0


class NutrientSecretionParameters:
    """Protein co-ingestion parameters for a single meal event."""

    def __init__(
        self,
        protein_g: float,
        ka: float = NUTRIENT_SECRETION_KA_PER_MIN,
        pm_per_g: float = NUTRIENT_SECRETION_PM_PER_G,
        meal_time_min: float = 0.0,
    ) -> None:
        if not math.isfinite(protein_g):
            raise ValueError(f"protein_g must be finite, got {protein_g}")
        if protein_g < 0:
            raise ValueError(f"protein_g cannot be negative, got {protein_g}")
        if protein_g > MAX_PLAUSIBLE_PROTEIN_G_PER_MEAL:
            raise ValueError(
                f"protein_g={protein_g} exceeds the documented plausibility "
                f"ceiling of {MAX_PLAUSIBLE_PROTEIN_G_PER_MEAL}g for a single "
                "meal -- see MAX_PLAUSIBLE_PROTEIN_G_PER_MEAL's own comment "
                "for why (verified against the real composed ODE's insulin "
                "output, not an arbitrary round number)"
            )
        if not math.isfinite(ka) or ka <= 0:
            raise ValueError(f"ka must be finite and positive, got {ka}")
        if not math.isfinite(pm_per_g) or pm_per_g < 0:
            raise ValueError(f"pm_per_g must be finite and non-negative, got {pm_per_g}")
        if not math.isfinite(protein_g * pm_per_g):
            raise ValueError(
                f"protein_g * pm_per_g (S1_0) must be finite, got "
                f"protein_g={protein_g}, pm_per_g={pm_per_g}"
            )
        if not math.isfinite(meal_time_min):
            raise ValueError(f"meal_time_min must be finite, got {meal_time_min}")
        self.protein_g = protein_g
        self.ka = ka
        self.pm_per_g = pm_per_g
        self.meal_time_min = meal_time_min

    @property
    def S1_0(self) -> float:
        """Total AUC contribution (pM.min) of this meal's protein-driven secretion."""
        return self.protein_g * self.pm_per_g

    def delta_pM_per_min(self, t_minutes: float) -> float:
        """
        The nutrient_secretion_delta signal at time t_minutes, for THIS
        meal alone. Exactly 0.0 before meal_time_min and for protein_g=0.
        """
        t_since = t_minutes - self.meal_time_min
        if t_since <= 0.0 or self.S1_0 <= 0.0:
            return 0.0
        return self.ka * self.S1_0 * math.exp(-self.ka * t_since)


class NutrientSecretionPathway(BasePathway):
    """
    Stateless (zero-species) pathway injecting a time-varying
    nutrient_secretion_delta signal into GlucoseInsulinPathway, summing
    the contribution of every meal event that co-ingested protein.
    """

    NAME    = "nutrient_secretion_v1"
    VERSION = "1.0.0"

    def __init__(self, meals: Optional[list[NutrientSecretionParameters]] = None) -> None:
        self._meals = list(meals or [])

    @property
    def descriptor(self) -> PathwayDescriptor:
        return PathwayDescriptor(
            name=self.NAME,
            version=self.VERSION,
            maturity=PathwayMaturity.EXPERIMENTAL,
            description=(
                "Stateless dietary-protein-driven endogenous insulin "
                "secretion -- protein co-ingested with a meal stimulates "
                "beta-cell secretion independent of glucose, via a "
                "closed-form monoexponential-decay secretion-rate curve "
                "per protein-bearing meal."
            ),
            species_names=[],
            output_signals=["nutrient_secretion_delta"],
            input_signals=[],
            dependencies=[],
            tags=["insulin", "meal", "protein", "endogenous-secretion", "time-varying"],
            literature_refs=["PMID:6389060", "PMID:3538843"],
        )

    @property
    def species_names(self) -> list[str]:
        return []

    @property
    def species_units(self) -> list[str]:
        return []

    @property
    def y0(self) -> list[float]:
        return []

    def dydt(
        self,
        t: float,
        y: list[float],
        signals: Optional[dict] = None,
    ) -> list[float]:
        return []

    def output_signals(self, t: float, y: list[float]) -> dict[str, float]:
        total = sum(m.delta_pM_per_min(t) for m in self._meals)
        return {"nutrient_secretion_delta": total}


def nutrient_secretion_couplings(gi_name: str, nutrient_name: str) -> list:
    """
    Build the CouplingSpec connecting a NutrientSecretionPathway to
    GlucoseInsulinPathway. Single signal, maps straight through -- no
    transform needed, since nutrient_secretion_delta is already expressed
    in GI's own signal vocabulary (rule 11).
    """
    from .coupling import CouplingSpec

    return [
        CouplingSpec(
            source_pathway=nutrient_name, signal_name="nutrient_secretion_delta",
            target_pathway=gi_name, target_input="nutrient_secretion_delta",
        ),
    ]


def couple_with_gi(gi: "GlucoseInsulinPathway", nutrient: "NutrientSecretionPathway") -> "CompositePathway":
    """
    Build the CompositePathway connecting a NutrientSecretionPathway to
    GlucoseInsulinPathway (see nutrient_secretion_couplings() for the
    signal contract).
    """
    from .coupling import CompositePathway

    return CompositePathway([gi, nutrient], couplings=nutrient_secretion_couplings(gi.name, nutrient.name))
