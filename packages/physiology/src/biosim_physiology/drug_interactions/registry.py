# SAVE TO: packages\physiology\src\biosim_physiology\drug_interactions\registry.py
"""Drug dispatch table, drug-class sets and scenario combination checks."""

from __future__ import annotations

from ._core import Drug, DrugDose
from .alpha_glucosidase import _alpha_glucosidase_inhibitor_modifiers
from .classic import (
    alcohol_modifiers,
    caffeine_modifiers,
    glucocorticoid_modifiers,
    metformin_modifiers,
)
from .incretins import _glp1_agonist_modifiers
from .insulin_sensitivity import (
    atorvastatin_modifiers,
    growth_hormone_modifiers,
    hydrochlorothiazide_modifiers,
    hydroxychloroquine_modifiers,
    indinavir_modifiers,
    pioglitazone_modifiers,
)
from .ode_coupled import _exogenous_insulin_modifiers, _sglt2_inhibitor_modifiers
from .sulfonylureas import _sulfonylurea_modifiers

_DRUG_FUNCTIONS = {
    Drug.CAFFEINE:       caffeine_modifiers,
    Drug.METFORMIN:      metformin_modifiers,
    Drug.GLUCOCORTICOID: glucocorticoid_modifiers,
    Drug.ALCOHOL:        alcohol_modifiers,
    Drug.SEMAGLUTIDE:      _glp1_agonist_modifiers,
    Drug.LIRAGLUTIDE:      _glp1_agonist_modifiers,
    Drug.DPP4_INHIBITOR:   _glp1_agonist_modifiers,
    Drug.TIRZEPATIDE:      _glp1_agonist_modifiers,
    Drug.DAPAGLIFLOZIN:    _sglt2_inhibitor_modifiers,
    Drug.EMPAGLIFLOZIN:    _sglt2_inhibitor_modifiers,
    Drug.GLIBENCLAMIDE:    _sulfonylurea_modifiers,
    Drug.GLIPIZIDE:        _sulfonylurea_modifiers,
    Drug.GLIMEPIRIDE:      _sulfonylurea_modifiers,
    Drug.ACARBOSE:         _alpha_glucosidase_inhibitor_modifiers,
    Drug.MIGLITOL:         _alpha_glucosidase_inhibitor_modifiers,
    Drug.INSULIN_REGULAR:  _exogenous_insulin_modifiers,
    Drug.INSULIN_LISPRO:   _exogenous_insulin_modifiers,
    Drug.PIOGLITAZONE:     pioglitazone_modifiers,
    Drug.INDINAVIR:        indinavir_modifiers,
    Drug.HYDROCHLOROTHIAZIDE: hydrochlorothiazide_modifiers,
    Drug.GROWTH_HORMONE:   growth_hormone_modifiers,
    Drug.ATORVASTATIN:     atorvastatin_modifiers,
    Drug.HYDROXYCHLOROQUINE: hydroxychloroquine_modifiers,
}


# Drugs whose real effect is realised via an ODE-coupled composite pathway
# in SimulationRunner rather than (or in addition to) a PathwayParameters
# modifier. Dose is in mg, matching DrugEvent.dose_mg / DrugDose.dose_mg.
#
# Drug.TIRZEPATIDE is DELIBERATELY NOT a member of this set (see
# _glp1_agonist_modifiers()'s own docstring): its real, citation-anchored
# effect on this platform is a static Gb modifier, not an ODE-coupled
# GLP1Pathway -- composing GLP1Pathway for it would require an uncited
# incretin_effect magnitude this session could not find (both papers
# describing it remain paywalled).
GLP1_DRUGS = {Drug.SEMAGLUTIDE, Drug.LIRAGLUTIDE, Drug.DPP4_INHIBITOR}

# Same category as GLP1_DRUGS above (ODE-coupled, not a static modifier) --
# see _sglt2_inhibitor_modifiers()'s docstring for why this one is even
# stricter than GLP-1 about not also having a static component.
SGLT2_DRUGS = {Drug.DAPAGLIFLOZIN, Drug.EMPAGLIFLOZIN}

# Same category as GLP1_DRUGS/SGLT2_DRUGS above (ODE-coupled, not a
# static modifier) -- see _exogenous_insulin_modifiers()'s docstring.
# Dose (DrugEvent/DrugDose.dose_mg) is interpreted as INSULIN UNITS for
# these drugs specifically, not milligrams (see this module's own
# docstring for the disclosed reasoning). Both members calibrated from
# the SAME within-study comparison (Howey et al. 1994's own arms A vs.
# C, see exogenous_insulin_pathway.py's own module docstring), not two
# separately-sourced numbers.
EXOGENOUS_INSULIN_DRUGS = {Drug.INSULIN_REGULAR, Drug.INSULIN_LISPRO}

# Drugs whose real effect is a genuinely ADDITIVE PathwayParameters field
# (glucose_independent_secretion) rather than a multiplicative one -- see
# _sulfonylurea_modifiers()'s docstring for why this can't go through
# compose_all_modifiers()'s pipeline. Extended to glipizide (2026-08-12
# deep-search plan, item D1) via GLIPIZIDE_SECRETION_FRACTION_OF_KE_IB's
# own docstring citations, then to glimepiride via a fresh head-to-head-
# ratio search angle -- see GLIMEPIRIDE_SECRETION_FRACTION_OF_KE_IB's
# own docstring citations, resolving MASTER_PLAN.md row 1.3b's own
# remaining open item.
SULFONYLUREA_DRUGS = {Drug.GLIBENCLAMIDE, Drug.GLIPIZIDE, Drug.GLIMEPIRIDE}

# Drugs modelled as a FLAT, presence-gated factor (dose_mg is only a >0
# flag, not a scaling axis): each active event multiplies the FULL factor,
# so N duplicate events silently compound it (HCQ security-review HIGH: 20
# duplicate HCQ events reached the 1e-4 p3 ceiling). Exactly one active
# event per scenario is allowed. NOT applied to time-resolved, legitimately
# repeatable drugs (exogenous insulin, caffeine, alcohol, glucocorticoid).
SINGLE_INSTANCE_DRUGS = {Drug.GROWTH_HORMONE, Drug.HYDROXYCHLOROQUINE}

# Drugs whose modelled effect scales with dose but which cannot represent
# more than one event faithfully, so duplicates are rejected and one event
# must carry the total dose:
#   * static steady-state factors (metformin/pioglitazone/indinavir/HCTZ/
#     atorvastatin): each event applies the full factor for its OWN dose,
#     so events compound rather than add (2 x pioglitazone 15mg != 30mg;
#     2 x metformin 1000mg double-applies past the 2000mg ceiling);
#   * first-event-only dispatch (acarbose/miglitol kabs scaling in the
#     runner, semaglutide/liraglutide GLP1Pathway dose): later events are
#     silently dropped, so the result depends on event ORDER (security
#     review: acarbose 25+100mg peaked at 9.74 mM, 100+25mg at 7.67 mM).
# Not included: drugs whose effect is dose-independent here (sulfonylureas,
# DPP-4 inhibitor, tirzepatide, SGLT2 inhibitors -- duplicates are inert)
# and time-resolved, legitimately repeatable ones (insulin, caffeine,
# alcohol, glucocorticoid).
DOSE_SCALED_SINGLE_EVENT_DRUGS = {
    Drug.METFORMIN, Drug.PIOGLITAZONE, Drug.INDINAVIR,
    Drug.HYDROCHLOROTHIAZIDE, Drug.ATORVASTATIN,
    Drug.ACARBOSE, Drug.MIGLITOL, Drug.SEMAGLUTIDE, Drug.LIRAGLUTIDE,
}


def validate_drug_combination(doses) -> None:
    """
    Fail fast on drug-event combinations whose effects this model cannot
    represent faithfully. `doses` is any iterable of objects with `.drug`
    (a Drug or its string value) and `.dose_mg` -- so it serves both
    DrugDose (personalizer) and DrugEvent (ScenarioDefinition) without
    conversion. Zero-dose events apply no effect and are ignored.

    Raises ValueError if:
      * more than one DISTINCT sulfonylurea is active (the override applies
        only the first, so the result would be order-dependent -- and
        combining sulfonylureas is not a real regimen); or
      * a SINGLE_INSTANCE_DRUGS member appears in more than one active event; or
      * a DOSE_SCALED_SINGLE_EVENT_DRUGS member appears in more than one
        active event (use one event with the total dose instead).
    """
    active = [(Drug(d.drug), d.dose_mg) for d in doses if d.dose_mg > 0]
    sulfonylureas = {drug for drug, _ in active if drug in SULFONYLUREA_DRUGS}
    if len(sulfonylureas) > 1:
        names = ", ".join(sorted(d.value for d in sulfonylureas))
        raise ValueError(
            "at most one sulfonylurea may be active in a scenario "
            f"(got {names}) -- their effects are not additive in this model "
            "and the result would be order-dependent"
        )
    for drug in SINGLE_INSTANCE_DRUGS:
        if sum(1 for d, _ in active if d == drug) > 1:
            raise ValueError(
                f"at most one {drug.value} event may be active in a scenario "
                "-- its effect is a flat factor, so duplicate events would "
                "silently compound it; use a single event"
            )
    for drug in DOSE_SCALED_SINGLE_EVENT_DRUGS:
        if sum(1 for d, _ in active if d == drug) > 1:
            raise ValueError(
                f"at most one {drug.value} event may be active in a scenario "
                "-- this model cannot represent several doses of it (they "
                "would compound or be silently dropped); use a single event "
                "with the total dose"
            )

# Drugs whose real effect is on MealInput.kabs (meal absorption rate
# constant), applied directly by SimulationRunner at meal-construction
# time -- not a PathwayParameters modifier at all, unlike every set
# above. Both members have their own, separately-verified same-meal-
# challenge quantitative anchor and their own kabs-scale function
# (acarbose_kabs_scale()/miglitol_kabs_scale()) -- SimulationRunner
# dispatches to the correct one per drug identity, never reusing one
# drug's calibration for the other (see MASTER_PLAN.md row 1.4b(ii)).
ALPHA_GLUCOSIDASE_INHIBITOR_DRUGS = {Drug.ACARBOSE, Drug.MIGLITOL}


def modifiers_for_dose(dose: DrugDose) -> dict[str, float]:
    """A single dose's own modifier dict (not combined with any others)."""
    fn = _DRUG_FUNCTIONS.get(dose.drug)
    return fn(dose) if fn else {}


def drug_modifiers(drug_doses: list[DrugDose]) -> dict[str, float]:
    """
    Compute combined parameter modifiers from a list of drug doses.
    Multiple drugs compose multiplicatively per parameter.
    """
    from ..modifiers import compose_all_modifiers
    return compose_all_modifiers([modifiers_for_dose(dose) for dose in drug_doses])
