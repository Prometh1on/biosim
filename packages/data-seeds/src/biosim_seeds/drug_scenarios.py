# SAVE TO: packages\data-seeds\src\biosim_seeds\drug_scenarios.py
"""
Drug interaction scenario builder.

Pre-built scenarios combining meals with common drugs/substances.
Each scenario documents the expected effect and the biological mechanism.

Scenarios included:
    caffeine_before_meal()       -- caffeine taken before eating
    metformin_with_meal()        -- T2D medication, daily dosing
    steroid_meal()               -- glucocorticoid (e.g. prednisone) effect
    alcohol_before_meal()        -- alcohol 1h before eating
    morning_coffee_breakfast()   -- 200mg caffeine + medium breakfast
    semaglutide_scenario()       -- weekly GLP-1 agonist (Ozempic/Wegovy)
    liraglutide_scenario()       -- daily GLP-1 agonist (Victoza/Saxenda)
    dpp4_inhibitor_scenario()    -- DPP-4 inhibitor (sitagliptin/saxagliptin)

All scenarios are designed to be run through SimulationRunner directly
or through PlatformSession.simulate(scenario_definition).

The GLP-1 class scenarios are handled differently from the others:
SimulationRunner detects the semaglutide/liraglutide/dpp4_inhibitor drug
name and composes a real GLP1Pathway ODE with GlucoseInsulinPathway
(biosim_ode.glp1_pathway.couple_with_gi) rather than applying a static
parameter modifier -- see biosim_physiology.drug_interactions.GLP1_DRUGS.
"""

from __future__ import annotations
from typing import Optional

from biosim_api import DrugEvent, MealEvent, PersonalizationSpec, ScenarioDefinition
from biosim_physiology import PhysiologicalProfile


def caffeine_before_meal(
    caffeine_mg: float = 200.0,
    carbs_g: float = 60.0,
    profile=None,
) -> ScenarioDefinition:
    """
    Caffeine taken 45 minutes before a meal.

    Effect: reduces insulin sensitivity ~10-15%, leading to a slightly
    higher and more prolonged glucose peak after eating.

    Mechanism: adenosine receptor antagonism + catecholamine release
    Literature: Keijzers et al. 2002, Diabetes Care. PMID:11815511
    """
    return ScenarioDefinition(
        name=f"caffeine_{caffeine_mg:.0f}mg_before_meal",
        duration_s=14400,
        description=(
            f"{caffeine_mg:.0f}mg caffeine + {carbs_g:.0f}g carb meal. "
            f"Caffeine reduces insulin sensitivity ~10-15%."
        ),
        meal_events=[MealEvent(time_s=0.0, carbs_g=carbs_g, label=f"{carbs_g:.0f}g meal")],
        # Regression: hours_before_simulation was never set (defaulted to
        # 0.0), so DrugEvent.to_drug_dose()'s hours_since_ingestion was
        # always exactly 0 -- caffeine_modifiers()'s concentration_fraction
        # is `0.0 if t <= 0 else ...`, making this scenario's caffeine
        # effect exactly zero despite the docstring's explicit "45 minutes
        # before a meal" claim. Verified: this scenario's output used to be
        # bit-for-bit identical to the same meal with the drug event
        # removed entirely. 0.75h = 45 min, matching the docstring.
        drug_events=[DrugEvent(
            time_s=0.0, drug="caffeine", dose_mg=caffeine_mg, hours_before_simulation=0.75,
        )],
        personalization=PersonalizationSpec(physiological_profile=profile or PhysiologicalProfile()),
    )


def metformin_with_meal(
    carbs_g: float = 75.0,
    metformin_mg: float = 1000.0,
    profile=None,
) -> ScenarioDefinition:
    """
    Metformin (chronic daily dosing) + standard meal.

    Effect: reduces hepatic glucose output and improves insulin sensitivity.
    Should produce lower peak glucose vs the same meal without metformin.

    Mechanism: AMPK activation in liver, reduced gluconeogenesis
    Literature: Bailey & Turner 1996, NEJM. PMID:8569826
    """
    return ScenarioDefinition(
        name=f"metformin_{metformin_mg:.0f}mg_with_meal",
        duration_s=14400,
        description=(
            f"Metformin {metformin_mg:.0f}mg/day (steady state) + {carbs_g:.0f}g meal. "
            f"Expected: lower peak glucose, improved insulin sensitivity."
        ),
        meal_events=[MealEvent(time_s=0.0, carbs_g=carbs_g, label=f"{carbs_g:.0f}g meal")],
        drug_events=[DrugEvent(time_s=0.0, drug="metformin", dose_mg=metformin_mg)],
        personalization=PersonalizationSpec(physiological_profile=profile or PhysiologicalProfile()),
    )


def steroid_meal(
    carbs_g: float = 60.0,
    prednisone_mg: float = 20.0,
    profile=None,
    hours_before_meal: float = 6.0,
) -> ScenarioDefinition:
    """
    Glucocorticoid (prednisone) effect on postprandial glucose.

    Effect: substantially reduces insulin sensitivity and elevates
    fasting glucose. Simulates steroid-induced hyperglycaemia.

    Mechanism: reduced GLUT4 translocation, increased hepatic glucose output
    Literature: Rizza et al. 1982, J Clin Endocrinol Metab. PMID:7033265

    hours_before_meal: dose-to-meal interval, now genuinely load-bearing
    (see CLAUDE.md's glucocorticoid PK/PD milestone entry) -- the real
    glycaemic effect lags plasma PK by hours (a genomic mechanism), and
    onsets/peaks/resolves on a real timeline (glucocorticoid_pd_effect_
    fraction() in biosim_ode.cortisol_pathway). Default 6.0h places the
    meal at this curve's peak (matching this platform's own literature
    citations: Yuen et al. 2012's "midday" elevation for an ~08:00 dose,
    PMID:21973241), NOT the old always-simultaneous dose+meal -- a dose
    and meal at t=0 together would show almost no glycaemic effect yet
    (regression caught by code review: this used to make this scenario's
    own regression test pass by a physiologically meaningless margin,
    ~0.0007% of the peak value, since the meal's own glucose peak occurs
    at ~25 minutes, long before the effect has developed at all).
    """
    return ScenarioDefinition(
        name=f"steroid_{prednisone_mg:.0f}mg_with_meal",
        duration_s=14400,
        description=(
            f"Prednisone {prednisone_mg:.0f}mg ({hours_before_meal:.0f}h before) + "
            f"{carbs_g:.0f}g meal. Steroid-induced insulin resistance expected."
        ),
        meal_events=[MealEvent(time_s=0.0, carbs_g=carbs_g, label=f"{carbs_g:.0f}g meal")],
        drug_events=[DrugEvent(
            time_s=0.0, drug="glucocorticoid", dose_mg=prednisone_mg,
            hours_before_simulation=hours_before_meal,
        )],
        personalization=PersonalizationSpec(physiological_profile=profile or PhysiologicalProfile()),
    )


def alcohol_before_meal(
    carbs_g: float = 60.0,
    alcohol_g: float = 28.0,   # 2 standard drinks
    hours_before: float = 1.0,
    profile=None,
) -> ScenarioDefinition:
    """
    Alcohol consumed before eating.

    Effect: transient suppression of hepatic glucose OUTPUT (a flux),
    resolving as alcohol clears -- NOT a shift in fasting glucose
    CONCENTRATION (Gb), which the real cited evidence explicitly rules
    out. See biosim_ode.alcohol_hepatic_pathway's own module docstring
    for the full derivation.

    Note: alcohol_g maps to dose_mg as 1g alcohol = 1000mg in the drug model.

    Mechanism: alcohol transiently suppresses hepatic glucose production
    (modelled via AlcoholHepaticPathway's ODE-coupled hepatic_egp_delta
    signal, not a static PathwayParameters modifier).
    Literature: Siler et al. 1998 (PMID:9815011, 48g oral alcohol, -12%
                hepatic glucose output); Shelmet et al. 1988
                (PMID:3280601, IV ethanol, corroborating peripheral/
                hepatic flux effects) -- both live-verified. Resolves the
                prior "citation unverified" placeholder (see
                biosim_physiology.drug_interactions module docstring).
    """
    dose_mg = alcohol_g * 1000.0   # convert g to mg for DrugDose
    return ScenarioDefinition(
        name=f"alcohol_{alcohol_g:.0f}g_before_meal",
        duration_s=14400,
        description=(
            f"{alcohol_g:.0f}g alcohol ({alcohol_g/14:.1f} standard drinks), "
            f"{hours_before:.0f}h before {carbs_g:.0f}g meal."
        ),
        meal_events=[MealEvent(time_s=0.0, carbs_g=carbs_g, label=f"{carbs_g:.0f}g meal")],
        # Regression: hours_before was encoded into time_s
        # (max(0, hours_before*3600-3600)), but DrugEvent.to_drug_dose()
        # computes hours_since_ingestion from hours_before_simulation, NOT
        # from time_s -- for any simulation_start_s=0.0 call (every real
        # caller), time_s's contribution is always max(0, 0-time_s)/3600=0
        # regardless of its value. hours_before had zero effect: sweeping
        # it from 0.25 to 5.0 produced bit-for-bit identical glucose/insulin
        # output every time. Fixed by using hours_before_simulation
        # directly, the same pattern library.py's caffeine_and_meal_scenario()
        # already uses correctly.
        drug_events=[DrugEvent(
            time_s=0.0, drug="alcohol", dose_mg=dose_mg, hours_before_simulation=hours_before,
        )],
        personalization=PersonalizationSpec(physiological_profile=profile or PhysiologicalProfile()),
    )


def morning_coffee_breakfast(
    profile=None,
) -> ScenarioDefinition:
    """
    A realistic morning scenario: coffee + toast + OJ.
    200mg caffeine + 44g carbohydrate (2 slices white bread + OJ).

    Coffee is modelled as taken 30 minutes before eating (a typical
    "coffee while getting ready" gap) -- same hours_before_simulation
    convention as simulation_api.library.caffeine_and_meal_scenario()'s
    default. Regression: this used to leave hours_before_simulation
    unset (defaulting to 0.0), making caffeine_modifiers()'s effect
    exactly zero regardless of dose -- see caffeine_before_meal()'s
    fix note above for the full mechanism.
    """
    return ScenarioDefinition(
        name="morning_coffee_breakfast",
        duration_s=14400,
        description="200mg caffeine (1 strong coffee) + toast + OJ (44g carbs).",
        meal_events=[MealEvent(time_s=0.0, carbs_g=44.0, label="Toast + OJ")],
        drug_events=[DrugEvent(
            time_s=0.0, drug="caffeine", dose_mg=200.0, hours_before_simulation=0.5,
        )],
        personalization=PersonalizationSpec(physiological_profile=profile or PhysiologicalProfile()),
    )


def semaglutide_scenario(
    carbs_g: float = 75.0,
    dose_mg: float = 1.0,
    profile=None,
) -> ScenarioDefinition:
    """
    Semaglutide (Ozempic/Wegovy) -- weekly GLP-1 receptor agonist -- + meal.

    Effect: substantially blunts the postprandial glucose peak via three
    incretin mechanisms: amplified glucose-dependent insulin secretion,
    suppressed glucagon (implicit in the incretin effect term), and slowed
    gastric emptying. Modelled via the real ODE-coupled GLP1Pathway
    (biosim_ode.glp1_pathway.couple_with_gi), not a static modifier.

    dose_mg: steady-state weekly dose (0.25-2.0mg typical range).

    Mechanism: GLP-1 receptor agonism, sustained (~1 week half-life)
    Literature: Lau et al. 2015, J Med Chem. PMID:26308095
    """
    return ScenarioDefinition(
        name=f"semaglutide_{dose_mg:.2f}mg_with_meal",
        duration_s=14400,
        description=(
            f"Semaglutide {dose_mg:.2f}mg/week (steady state) + {carbs_g:.0f}g meal. "
            f"Expected: substantially lower glucose peak via incretin amplification "
            f"of insulin secretion and slowed gastric emptying."
        ),
        meal_events=[MealEvent(time_s=0.0, carbs_g=carbs_g, label=f"{carbs_g:.0f}g meal")],
        drug_events=[DrugEvent(time_s=0.0, drug="semaglutide", dose_mg=dose_mg)],
        personalization=PersonalizationSpec(physiological_profile=profile or PhysiologicalProfile()),
    )


def liraglutide_scenario(
    carbs_g: float = 75.0,
    dose_mg: float = 1.2,
    profile=None,
) -> ScenarioDefinition:
    """
    Liraglutide (Victoza/Saxenda) -- daily GLP-1 receptor agonist -- + meal.

    Effect: lowers glucose peak, less pronounced than semaglutide due to
    liraglutide's shorter half-life (~13h vs ~1 week).

    dose_mg: standard T2D dose 1.2-1.8 mg/day.

    Mechanism: GLP-1 receptor agonism, daily dosing
    Literature: Buse et al. 2009, Lancet (LEAD-6). PMID:19515413
    """
    return ScenarioDefinition(
        name=f"liraglutide_{dose_mg:.2f}mg_with_meal",
        duration_s=14400,
        description=(
            f"Liraglutide {dose_mg:.2f}mg/day (steady state) + {carbs_g:.0f}g meal. "
            f"Expected: lower glucose peak via incretin amplification of insulin "
            f"secretion and slowed gastric emptying."
        ),
        meal_events=[MealEvent(time_s=0.0, carbs_g=carbs_g, label=f"{carbs_g:.0f}g meal")],
        drug_events=[DrugEvent(time_s=0.0, drug="liraglutide", dose_mg=dose_mg)],
        personalization=PersonalizationSpec(physiological_profile=profile or PhysiologicalProfile()),
    )


def dpp4_inhibitor_scenario(
    carbs_g: float = 75.0,
    profile=None,
) -> ScenarioDefinition:
    """
    DPP-4 inhibitor (sitagliptin/saxagliptin) + meal.

    Effect: modest glucose peak reduction -- DPP-4 inhibitors prevent
    endogenous GLP-1 degradation (doubling its effective half-life) rather
    than pharmacologically over-driving the receptor like an injectable
    agonist, so the effect is real but smaller than semaglutide/liraglutide.

    Mechanism: DPP-4 enzyme inhibition -> prolonged endogenous GLP-1 activity
    Literature: Deacon 2011, Diabetes Obes Metab. PMID:21114598
    """
    return ScenarioDefinition(
        name="dpp4_inhibitor_with_meal",
        duration_s=14400,
        description=(
            f"DPP-4 inhibitor (steady state) + {carbs_g:.0f}g meal. "
            f"Expected: modest glucose peak reduction via prolonged endogenous "
            f"GLP-1 activity."
        ),
        meal_events=[MealEvent(time_s=0.0, carbs_g=carbs_g, label=f"{carbs_g:.0f}g meal")],
        drug_events=[DrugEvent(time_s=0.0, drug="dpp4_inhibitor", dose_mg=100.0)],
        personalization=PersonalizationSpec(physiological_profile=profile or PhysiologicalProfile()),
    )
