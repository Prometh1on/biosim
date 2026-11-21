# SAVE TO: packages\data-seeds\src\biosim_seeds\condition_scenarios.py
"""
Physiological-condition scenario builder.

Pre-built scenarios exercising the ODE-coupled extension pathways wired
into SimulationRunner: acute psychological stress (CortisolPathway) and
type 2 diabetes hepatic kinetics (HepGlucPathway), alongside the existing
static personalisation modifiers.

Scenarios included:
    stress_meal_scenario()             -- acute stress + meal (CortisolPathway)
    type2_diabetes_meal_scenario()     -- T2D + meal (HepGlucPathway + disease modifiers)

See biosim_ode.cortisol_pathway.cortisol_couplings() and
biosim_ode.hepatic.hepatic_couplings() for why these compose additively
with the existing physiology.modifiers system rather than duplicating it.
"""

from __future__ import annotations

from biosim_api import MealEvent, PersonalizationSpec, ScenarioDefinition
from biosim_physiology import PhysiologicalProfile, DiseaseState


def stress_meal_scenario(
    carbs_g: float = 75.0,
    stress_level: float = 0.7,
    profile=None,
) -> ScenarioDefinition:
    """
    Acute psychological/physical stress + meal.

    Effect: cortisol elevation from stress_level (0-1) raises fasting
    glucose and reduces insulin action via the real HPA-axis ODE
    (CortisolPathway), independent of and additive to the existing
    empirical time-of-day circadian modifier.

    Mechanism: CRH -> cortisol -> glucocorticoid receptor occupancy ->
    reduced insulin action (p3) and elevated glucose set point (Gb)
    Literature: Rizza et al. 1982, J Clin Endocrinol Metab. PMID:7033265;
                Dinneen et al. 1993, J Clin Invest. PMID:8227343
    """
    return ScenarioDefinition(
        name=f"stress_{stress_level:.1f}_with_meal",
        duration_s=14400,
        description=(
            f"Acute stress (level {stress_level:.1f}) + {carbs_g:.0f}g meal. "
            f"Expected: elevated glucose peak via HPA-axis cortisol response."
        ),
        meal_events=[MealEvent(time_s=0.0, carbs_g=carbs_g, label=f"{carbs_g:.0f}g meal")],
        stress_level=stress_level,
        personalization=PersonalizationSpec(physiological_profile=profile or PhysiologicalProfile()),
    )


def type2_diabetes_meal_scenario(
    carbs_g: float = 75.0,
    profile=None,
) -> ScenarioDefinition:
    """
    Type 2 diabetes + meal.

    Effect: substantially higher and more prolonged glucose excursion.
    Combines the existing static disease-state modifiers (elevated Gb,
    reduced p3/secretion_scale) with impaired hepatic EGP suppression
    kinetics during the postprandial period (HepGlucPathway), which the
    static modifiers alone cannot represent.

    Mechanism: beta-cell dysfunction + peripheral insulin resistance +
    impaired hepatic glucose output suppression
    Literature: DeFronzo 2004, Med Clin North Am. PMID:15308380
    """
    base_profile = profile or PhysiologicalProfile()
    existing = [s for s in base_profile.disease_states if s != DiseaseState.HEALTHY]
    new_states = existing if DiseaseState.TYPE_2_DIABETES in existing else existing + [DiseaseState.TYPE_2_DIABETES]
    import dataclasses
    t2d_profile = dataclasses.replace(base_profile, disease_states=new_states)
    return ScenarioDefinition(
        name="type2_diabetes_with_meal",
        duration_s=14400,
        description=(
            f"Type 2 diabetes + {carbs_g:.0f}g meal. "
            f"Expected: substantially elevated and prolonged glucose peak."
        ),
        meal_events=[MealEvent(time_s=0.0, carbs_g=carbs_g, label=f"{carbs_g:.0f}g meal")],
        personalization=PersonalizationSpec(physiological_profile=t2d_profile),
    )
