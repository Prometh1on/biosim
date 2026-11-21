# SAVE TO: packages\simulation-api\src\biosim_api\scenario_warnings.py
"""
Stateless disclosure warnings for SimulationRunner.run().

Moved out of runner.py (2026-10-20, MASTER_PLAN Section 8) when run() was
split into pipeline stages; each function only builds warning text, so
behaviour is unchanged.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Mapping, Optional

if TYPE_CHECKING:
    from .scenario import ScenarioDefinition


def clamp_hit_warnings(clamp_hits: Mapping[tuple[str, str], float]) -> list[str]:
    """
    Disclose dydt()'s point-of-use clamps (EFFECTIVE_PARAM_BOUNDS). They
    engage on ODE-coupled signals (cortisol, NEFA, exercise...) applied at
    every solver step, which the personalizer's static clamp_warnings
    never sees -- without this, a clamped run looked identical to a clean one.
    """
    from biosim_ode.pathway import EFFECTIVE_PARAM_BOUNDS

    warnings = []
    for (name, side), raw in sorted(clamp_hits.items()):
        lower, upper = EFFECTIVE_PARAM_BOUNDS[name]
        bound = lower if side == "floor" else upper
        warnings.append(
            f"{name} hit its point-of-use {side} ({bound:.3g}) during "
            f"integration: the coupled-pathway-modified value reached "
            f"{raw:.3g}, so the simulation used the clamped value -- the "
            "composed effect is stronger than the model can represent."
        )
    return warnings


def unverified_variant_warnings(variants) -> list[str]:
    """
    One note per input variant whose effect size has no verified source
    (UNVERIFIED_VARIANT_NOTES, e.g. the IGF1 promoter variant), only when the
    genotype carries at least one copy of the effect allele -- with zero
    copies the effect is never applied. The note never contains the rsID
    (rule 7).
    """
    from biosim_genotype.effect_db import UNVERIFIED_VARIANT_NOTES, VariantEffectDB

    db = VariantEffectDB.default()

    def carries_effect_allele(rsid: str) -> bool:
        return any(
            v.dosage_for_allele(record.risk_allele) > 0
            for v in variants if v.rsid == rsid
            for record in db.get(rsid)
        )

    return [
        note for rsid, note in UNVERIFIED_VARIANT_NOTES.items()
        if carries_effect_allele(rsid)
    ]


def experimental_drug_warnings(scenario: "ScenarioDefinition") -> list[str]:
    """
    One "EXPERIMENTAL drug module" note (or, for a validated drug, a
    "VALIDATED drug module" note with its caveat) per distinct drug the
    scenario applies
    (biosim_physiology.drug_interactions.validation_status), in the same style
    as unverified_variant_warnings(). A zero dose applies nothing. Exogenous
    insulin has its own event type (doses in Units). Disclosure only -- a drug
    run is never blocked.
    """
    from biosim_physiology.drug_interactions.validation_status import (
        drug_validation_notes,
    )
    applied = [e.drug for e in scenario.drug_events if e.dose_mg > 0]
    applied += [
        f"insulin_{e.insulin_type}" for e in scenario.insulin_dose_events
        if e.dose_units > 0
    ]
    return drug_validation_notes(applied)


def _diabetes_states(profile) -> list[str]:
    from biosim_physiology.profile import DiseaseState
    diabetic = {
        DiseaseState.TYPE_2_DIABETES, DiseaseState.TYPE_1_DIABETES,
        DiseaseState.GESTATIONAL_DIABETES, DiseaseState.GCK_MODY,
    }
    return sorted(d.value for d in profile.disease_states if d in diabetic)


def hcq_cohort_warning(profile, drug_events) -> Optional[str]:
    """
    Track 1 row 1.44: Wasko et al. 2015 (PMID:26197707) studied NON-diabetic
    overweight/obese adults (BMI >= 25). Disclose, don't block, any other
    population -- same convention as the cold-acclimation warning.
    """
    if not any(de.drug == "hydroxychloroquine" and de.dose_mg > 0 for de in drug_events):
        return None
    reasons = []
    diabetic = _diabetes_states(profile)
    if diabetic:
        reasons.append(f"the profile has {', '.join(diabetic)}")
    if profile.bmi < 25.0:
        reasons.append(f"its BMI is {profile.bmi:.1f} (< 25)")
    if not reasons:
        return None
    return (
        "hydroxychloroquine's p3 factor comes from non-diabetic overweight/"
        "obese adults (Wasko et al. 2015, PMID:26197707), but "
        + " and ".join(reasons)
        + " -- applying it here is an unverified extrapolation beyond that cohort."
    )


def nafld_cohort_warning(profile, model_nefa_dynamics: bool) -> Optional[str]:
    """
    Track 1 row 1.45: Bugianesi et al. 2005 (PMID:15747110) studied
    NON-diabetic NAFLD; and NAFLD's p3 reduction overlaps other insulin-
    resistance mechanisms, so stacking them probably double-counts.
    """
    from biosim_physiology.profile import DiseaseState
    if DiseaseState.NAFLD not in profile.disease_states:
        return None
    reasons = []
    diabetic = _diabetes_states(profile)
    if diabetic:
        reasons.append(
            f"its source cohort was non-diabetic, but the profile has {', '.join(diabetic)}"
        )
    overlapping = sorted(
        d.value for d in profile.disease_states
        if d in {DiseaseState.INSULIN_RESISTANT, DiseaseState.TYPE_2_DIABETES}
    )
    if model_nefa_dynamics:
        overlapping.append("NEFA dynamics (adipose insulin resistance)")
    if overlapping:
        reasons.append(
            "composing it with " + ", ".join(overlapping) + " probably "
            "double-counts part of the same insulin resistance"
        )
    if not reasons:
        return None
    return (
        "NAFLD's p3 factor (Bugianesi et al. 2005, PMID:15747110): "
        + "; ".join(reasons) + "."
    )


def meal_composition_warnings(meal_events) -> list[str]:
    warnings: list[str] = []
    # GutAbsorptionModel.for_gi() (GI set but < 99.5) is documented as
    # uncalibrated and, per CLAUDE.md VALIDATION STATUS, backed by two
    # real studies (Cocate 2011 PMID:21208446, Louie 2012/2013
    # PMID:23344248) that directly disagree on whether GI drives peak
    # timing at all. No production scenario used this path when that was
    # last confirmed, but if one ever does, that should be loud in the
    # output rather than silently trusted as calibrated.
    from biosim_ode.gut_absorption import GI_SLOPE_PLACEHOLDER_NOTE

    for event in meal_events:
        if event.glycaemic_index is not None and event.glycaemic_index < 99.5:
            warnings.append(
                f"Meal '{event.label or event.carbs_g}' uses GutAbsorptionModel.for_gi"
                f"({event.glycaemic_index}). {GI_SLOPE_PLACEHOLDER_NOTE} Clinical "
                "studies also conflict on whether GI drives peak timing at all (see "
                "CLAUDE.md VALIDATION STATUS)."
            )

    # NutrientSecretionPathway (protein co-ingestion, MASTER_PLAN row
    # 1.28) does not model glucagon's counter-regulatory rise -- both
    # of its own calibration papers (Krezowski 1986/Nuttall 1984)
    # report protein also raises glucagon, and row 1.1 (glucagon's
    # own fasting/basal EGP effect) is permanently blocked, so this
    # mechanism will OVER-PREDICT glucose lowering for a meal where
    # protein dominates over carbohydrate -- exactly the same
    # "loud, not silently trusted" convention as the for_gi() warning
    # above. Threshold (protein_g > carbs_g) is a disclosed
    # engineering heuristic for "protein-dominant," not a literature
    # cutoff -- a meal with more protein than carbohydrate is the
    # case this mechanism's own calibration (50g protein ALONE vs.
    # 50g glucose ALONE) least resembles a "co-ingested with a normal
    # carb-dominant meal" scenario.
    for event in meal_events:
        if event.protein_g > 0 and event.protein_g > event.carbs_g:
            warnings.append(
                f"Meal '{event.label or event.carbs_g}' has protein_g="
                f"{event.protein_g} exceeding carbs_g={event.carbs_g} -- "
                "NutrientSecretionPathway does not model glucagon's "
                "counter-regulatory rise and will OVER-PREDICT glucose "
                "lowering for protein-dominant meals (see "
                "nutrient_secretion_pathway.py's own docstring). Treat "
                "glucose results from this meal as illustrative, not "
                "clinically calibrated."
            )
    return warnings


def cushings_warning(scenario: "ScenarioDefinition") -> Optional[str]:
    # Security-reviewer finding on row 1.38 (Cushing's syndrome): the
    # module-level comment on for_cushings() already discloses that
    # Cushing's textbook-defining feature -- hypercortisolism-driven
    # fasting hyperglycemia -- is NOT represented (gb_modifier stays
    # neutral; only the smaller postprandial p3_modifier effect is
    # modeled), but a code comment is invisible to a REST consumer.
    # Every other comparable partial-implementation gap in
    # SimulationRunner.run() (the for_gi() slope, the protein-dominant-meal
    # glucagon gap, the GLP-1/disease-Gb interaction) is
    # surfaced via `warnings` -- this one now matches that convention.
    if not scenario.cushings_syndrome:
        return None
    return (
        "cushings_syndrome=True models Cushing's syndrome's postprandial "
        "insulin-action impairment (p3_modifier) only -- it does NOT model "
        "the disease's defining fasting hyperglycemia (Gb is left "
        "unchanged). Do not treat this scenario's fasting/baseline "
        "glucose as clinically representative of Cushing's syndrome."
    )


def cold_acclimation_warning(profile) -> Optional[str]:
    # Track 1 row 1.40 (cold acclimation): Hanssen 2015's own cohort
    # (PMID:26147760) was n=8 TYPE 2 DIABETES subjects with NO healthy
    # arm -- applying its p3 factor to a non-diabetic profile is an
    # unverified extrapolation, disclosed rather than silently
    # generalized, same convention as the Cushing's warning above and
    # SimulationRunner._glp1_basal_glucose_warning() (runner.py).
    from biosim_physiology.profile import DiseaseState
    if profile is None or not profile.cold_acclimatized \
            or DiseaseState.TYPE_2_DIABETES in profile.disease_states:
        return None
    return (
        "cold_acclimatized=True is applied, but its source citation "
        "(Hanssen et al. 2015, PMID:26147760) studied Type 2 Diabetes "
        "subjects only (no healthy-subject arm) -- applying this "
        "insulin-sensitising effect to a profile without "
        "DiseaseState.TYPE_2_DIABETES is an unverified extrapolation "
        "beyond this citation's own evidence base."
    )
