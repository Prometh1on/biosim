# SAVE TO: packages\ode-engine\src\biosim_ode\meal_absorption_pathway.py
"""
Meal absorption pathway -- mass-conserving gastric-emptying coupling.

Fixes a real bug (see CLAUDE.md Track 1 "gastric_emptying_modifier
mass-conservation" entry, numerically confirmed before this pathway was
written): `GlucoseInsulinPathway.dydt()` used to compute
`Ra = total_Ra(self.meals, t) * gastric_emptying_modifier`, multiplying
the meal appearance-rate curve's OUTPUT by GLP1Pathway's modifier. The
analytic property `integral_0^inf(D*k*e^-kt dt) = D` holds for any FIXED
rate constant k, but scaling the curve's amplitude by any factor c != 1
(constant or time-varying) scales the total absorbed mass by c too --
`integral(Ra*c dt) = c*D`, not D. Confirmed numerically: a chronic
semaglutide scenario (1.0mg, 75g meal, 5h) lost ~42% of the meal's
carbohydrate mass; liraglutide ~34%; a DPP-4 inhibitor ~13% -- the
platform was silently deleting glucose, not delaying it, for every
GLP-1 class drug.

Mass-conserving redesign: instead of computing Ra(t) from a closed-form
monoexponential evaluated fresh at each t, track the REMAINING
(unabsorbed) glucose mass per meal as its own ODE state:

    remaining_i(t_meal) = D_eff_i                         (full dose, not yet absorbed)
    d(remaining_i)/dt   = -gate_i(t) * kabs_i * modifier(t) * remaining_i(t)
    Ra_i(t)             =  gate_i(t) * kabs_i * modifier(t) * remaining_i(t)

where gate_i(t) = 1 for t >= meal_i.time_min else 0 (absorption hasn't
started yet). This conserves mass for ANY modifier(t) trajectory,
constant or time-varying: integral_onset^inf Ra_i dt = remaining_i(onset)
- remaining_i(inf) = D_eff_i - 0 = D_eff_i, as long as the modifier never
holds absorption at exactly zero forever (GLP1Pathway's own modifier is
floored at 0.30, so this always holds for every currently-shipped drug
preset). When modifier(t) == 1 for all t, remaining_i(t) reduces to the
exact analytic solution D_eff_i * exp(-kabs_i * dt), and Ra_i(t) is
bit-for-bit (to solver tolerance) identical to today's
`MealInput.Ra()` -- this pathway is a pure generalization, not a
behavior change, when no gastric-emptying modifier is present.

Architecture, following the established extension-pathway pattern
(GLP1Pathway/RenalPathway) rather than adding new species to
GlucoseInsulinPathway itself (which would grow GI's species count and
leak meal-tracking internals into `SimulationResult.species` in a way
that -- unlike GLP1Pathway's own single "GLP1" key, already a precedented
extra key today -- would be one key PER MEAL, dynamically named and
sized per scenario):

    GLP1.gastric_emptying_modifier -> MealAbsorptionPathway (input)
    MealAbsorptionPathway.meal{i}_remaining_mM -> GI (input, per meal)

GI's own dydt() reads `meal{i}_remaining_mM` (this pathway's output) plus
its own already-existing `gastric_emptying_modifier` input (unchanged
coupling from GLP1) to compute Ra itself -- see pathway.py's dydt() Ra
block. When this pathway is NOT composed (no GLP-1 class drug active,
the overwhelming majority of scenarios), GI falls back to the original
closed-form `meal.Ra(t) * modifier` computation, so every scenario that
doesn't involve a gastric-emptying-modifying drug is completely
unaffected -- confined blast radius, not a platform-wide behavior change.

Scope: every meal is state-tracked. Two-compartment (`GutAbsorptionModel`)
meals -- since 2026-11-15 the default for scenario meals -- get a stomach
and a gut state each, with the modifier slowing stomach emptying only (see
the class docstring). This lifted the old guard that rejected a GLP-1-class
drug combined with a `glycaemic_index` meal: before, a DPP-4 inhibitor +
75 g GI=100 OGTT in T2D deleted ~38% of the meal's glucose, because the
closed-form two-compartment Ra(t) was multiplied by the modifier.
"""

from __future__ import annotations

from typing import Optional

from .base_pathway import BasePathway, PathwayDescriptor, PathwayMaturity
from .kinetics import clamp
from .meal import MealInput


def _is_two_compartment(meal: MealInput) -> bool:
    return meal.absorption_model is not None


def _species_for_meal(i: int, meal: MealInput) -> list[str]:
    if _is_two_compartment(meal):
        return [f"meal{i}_stomach_mM", f"meal{i}_gut_mM"]
    return [f"meal{i}_remaining_mM"]


class MealAbsorptionPathway(BasePathway):
    """
    Per-meal unabsorbed-glucose states, so a gastric-emptying modifier
    delays a meal's glucose instead of deleting it.

    Single-exponential meal (kabs): one state, "meal{i}_remaining_mM".
    Two-compartment meal (GutAbsorptionModel): two states,
    "meal{i}_stomach_mM" and "meal{i}_gut_mM":

        d(stomach)/dt = -gate * k_emptying * modifier * stomach
        d(gut)/dt     =  gate * k_emptying * modifier * stomach - k_absorption * gut
        Ra            =  k_absorption * gut

    The modifier slows only stomach emptying (what GLP-1 drugs do); every
    unit that leaves the stomach reaches the blood, so mass is conserved
    for any modifier(t). With modifier == 1 this is exactly the closed-form
    GutAbsorptionModel.Ra() curve.

    Species are labelled by the meal's ORIGINAL index in `meals`, so GI can
    look them up with the same index it iterates `self.meals` with.
    """

    NAME    = "meal_absorption_v1"
    VERSION = "1.1.0"

    def __init__(self, meals: Optional[list[MealInput]] = None) -> None:
        self._meals = meals or []
        self._species: list[str] = []
        self._y0: list[float] = []
        for i, meal in enumerate(self._meals):
            self._species.extend(_species_for_meal(i, meal))
            self._y0.append(meal.D_eff_mM)
            if _is_two_compartment(meal):
                self._y0.append(0.0)

    @property
    def descriptor(self) -> PathwayDescriptor:
        return PathwayDescriptor(
            name=self.NAME,
            version=self.VERSION,
            maturity=PathwayMaturity.BETA,
            description=(
                "Mass-conserving per-meal gastric-emptying coupling: tracks "
                "unabsorbed glucose per meal (one state for a single-exponential "
                "meal, stomach + gut for a two-compartment meal) so a "
                "time-varying gastric_emptying_modifier delays absorption "
                "instead of deleting carbohydrate mass."
            ),
            species_names=self.species_names,
            output_signals=list(self.species_names),
            input_signals=["gastric_emptying_modifier"],
            dependencies=[],
            tags=["meal", "absorption", "gastric-emptying", "GLP-1", "mass-conservation"],
            literature_refs=[],
        )

    @property
    def species_names(self) -> list[str]:
        return list(self._species)

    @property
    def species_units(self) -> list[str]:
        return ["mM"] * len(self._species)

    @property
    def y0(self) -> list[float]:
        return list(self._y0)

    def dydt(self, t: float, y: list[float], signals: Optional[dict] = None) -> list[float]:
        s = signals or {}
        modifier = s.get("gastric_emptying_modifier", 1.0)
        derivs: list[float] = []
        pos = 0
        for meal in self._meals:
            gate = 1.0 if t >= meal.time_min else 0.0
            first = clamp(y[pos], lower=0.0)
            if _is_two_compartment(meal):
                gut = clamp(y[pos + 1], lower=0.0)
                emptying = gate * meal.absorption_model.k_emptying * modifier * first
                derivs.append(-emptying)
                derivs.append(emptying - meal.absorption_model.k_absorption * gut)
                pos += 2
            else:
                derivs.append(-gate * meal.kabs * modifier * first)
                pos += 1
        return derivs

    def output_signals(self, t: float, y: list[float]) -> dict[str, float]:
        return {
            name: clamp(value, lower=0.0)
            for name, value in zip(self._species, y)
        }


def meal_absorption_couplings(
    glp1_name: str, absorption_name: str, gi_name: str, meals: list[MealInput],
) -> list:
    """
    Build the CouplingSpecs wiring a MealAbsorptionPathway between a
    GLP1Pathway (source of the time-varying gastric_emptying_modifier)
    and GlucoseInsulinPathway (consumer of each meal's state: remaining
    mass for a single-exponential meal, gut content for a two-compartment
    meal).

    Exposed standalone (matching glp1_couplings()/renal_couplings()'
    existing convention) so SimulationRunner can flat-compose this
    alongside GLP1 and GI in a single CompositePathway.
    """
    from .coupling import CouplingSpec

    couplings = [
        CouplingSpec(
            source_pathway=glp1_name, signal_name="gastric_emptying_modifier",
            target_pathway=absorption_name, target_input="gastric_emptying_modifier",
        ),
    ]
    for i, meal in enumerate(meals):
        signal_name = _species_for_meal(i, meal)[-1]
        couplings.append(CouplingSpec(
            source_pathway=absorption_name, signal_name=signal_name,
            target_pathway=gi_name, target_input=signal_name,
        ))
    return couplings
