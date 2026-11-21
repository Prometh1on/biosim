# SAVE TO: packages\physiology\src\biosim_physiology\drug_interactions\__init__.py
"""
Drug interaction modifiers for the glucose/insulin pathway.

Supported drugs (v1):
    CAFFEINE        acute reduction in insulin sensitivity
    METFORMIN       chronic improvement in insulin sensitivity + hepatic glucose
    GLUCOCORTICOID  reduction in insulin sensitivity + elevated fasting glucose
    ALCOHOL         complex biphasic effects (simplified for v1)
    DAPAGLIFLOZIN,
    EMPAGLIFLOZIN   SGLT2 inhibitors -- insulin-independent renal glucose
                    clearance, modelled entirely via an ODE-coupled pathway
                    (biosim_ode.renal_pathway.RenalPathway), not a modifier
                    here -- see _sglt2_inhibitor_modifiers() below.
    GLIBENCLAMIDE,
    GLIPIZIDE,
    GLIMEPIRIDE     Sulfonylureas -- glucose-INDEPENDENT insulin secretagogues
                    (KATP-channel closure via SUR1 binding), modelled as a
                    constant addition to PathwayParameters
                    .glucose_independent_secretion, applied directly by
                    PhysiologicalPersonalizer (not via the multiplicative
                    net_factors pipeline below -- see
                    glucose_independent_secretion_pM_per_min()). Each drug
                    has its own fraction (glipizide's real effect is ~2x
                    glibenclamide's per a head-to-head study, capped by
                    this model's own safety ceiling -- see
                    GLIPIZIDE_SECRETION_FRACTION_OF_KE_IB's own docstring;
                    glimepiride's is derived from a head-to-head fasting-
                    insulin RATIO to glibenclamide -- see
                    GLIMEPIRIDE_SECRETION_FRACTION_OF_KE_IB's own
                    docstring).
    TIRZEPATIDE     Dual GIP/GLP-1 receptor agonist -- unlike the other
                    GLP-1-class drugs above, modelled as a static Gb
                    modifier (not an ODE-coupled GLP1Pathway), since its
                    real, citation-anchored effect on this platform is
                    the fasting-glucose-lowering mechanism specifically
                    (see tirzepatide_gb_factor()/TIRZEPATIDE_GB_DOSE_MG
                    below) -- its OTHER real mechanism (a stronger
                    incretin_effect) remains uncited and unimplemented,
                    a deliberate partial implementation, not an oversight.
    ACARBOSE        Alpha-glucosidase inhibitor -- slows intestinal
                    carbohydrate breakdown/absorption (and, secondarily,
                    gastric emptying). Modelled as a direct reduction of
                    MealInput.kabs (the absorption RATE CONSTANT, not a
                    scale on Ra(t)'s output -- reducing kabs preserves
                    total carbohydrate mass exactly, since integral(D*k*
                    exp(-k*t) dt) = D for any k>0), applied by
                    SimulationRunner at meal-construction time -- not a
                    PathwayParameters modifier at all, so
                    _alpha_glucosidase_inhibitor_modifiers() below returns
                    {} like the ODE-coupled drugs above (for a different
                    reason: this drug's effect isn't on PathwayParameters
                    or an ODE-coupled pathway, it's on the meal inputs
                    themselves). See acarbose_kabs_scale() below.
    INSULIN_REGULAR Exogenous (injected) regular human insulin -- raises
                    plasma insulin from an EXTERNAL source rather than
                    via endogenous secretion, modelled entirely via an
                    ODE-coupled pathway (biosim_ode.exogenous_insulin_
                    pathway.ExogenousInsulinPathway), not a modifier here
                    -- see _exogenous_insulin_modifiers() below. Dose is
                    carried via DrugEvent/DrugDose.dose_mg, but interpreted
                    as INSULIN UNITS for this one drug, not milligrams --
                    a disclosed unit-reuse hack (insulin dosing is
                    conventionally in Units, not mg, and this platform has
                    no dedicated insulin-dose event type yet -- see
                    MASTER_PLAN row 1.2b).
    INSULIN_LISPRO  Exogenous (injected) rapid-acting insulin lispro
                    ([Lys(B28),Pro(B29)]-human insulin, zinc-free, brand
                    name Humalog) -- same ODE-coupled representation and
                    same dose_mg-as-Units convention as INSULIN_REGULAR
                    above, using ExogenousInsulinParameters
                    .rapid_acting_insulin() instead of .regular_insulin().
                    Calibrated from the SAME within-study comparison as
                    INSULIN_REGULAR (Howey et al. 1994's own arms A vs.
                    C), not a separately-sourced number.
    PIOGLITAZONE    Thiazolidinedione (PPAR-gamma insulin sensitizer) --
                    mechanistically distinct from every drug above (not a
                    secretagogue, incretin, glucosuric agent, or
                    absorption blocker). Modelled as a static p3
                    (whole-body insulin sensitivity) modifier, piecewise-
                    linearly interpolated across 5 real dose points
                    including an empirically-confirmed NULL zone at
                    7.5/15mg -- see pioglitazone_p3_factor()/
                    PIOGLITAZONE_P3_DOSE_MG below. Deliberately does NOT
                    also model the same paper's hepatic-ISI or
                    insulinogenic-index findings (would repeat KNOWN BUG
                    #14's over-determination mistake).

Each drug produces a time-varying or steady-state modifier on
pathway parameters. DrugDose carries the dose and timing needed
to compute the effect at simulation time.

Literature:
    Caffeine/insulin:    Keijzers et al. 2002, Diabetes Care. PMID:11815511
    Metformin:           Bailey & Turner 1996, NEJM. PMID:8569826
    Glucocorticoids:     Rizza et al. 1982, J Clin Endocrinol Metab. PMID:7033265
    Alcohol/glucose:     RESOLVED (was: citation unverified -- the previously-cited
                         PMID:9630250 was checked against PubMed and found unrelated
                         (snoRNA processing); no matching "Emanuele 1998" alcohol/
                         glucose paper could be confirmed; see CLAUDE.md known-bugs).
                         A fresh Track 1 scoping fan-out found the OLD Gb-based
                         mechanism was not just uncited but actively WRONG -- see
                         alcohol_modifiers()'s own docstring below and
                         biosim_ode.alcohol_hepatic_pathway's module docstring for
                         the real citations (Siler 1998 PMID:9815011, Shelmet 1988
                         PMID:3280601) and the corrected ODE-coupled mechanism.
    Sulfonylurea:        Riefflin A, Ayyagari U, Manley SE, Holman RR, Levy JC.
                         "The effect of glibenclamide on insulin secretion at
                         normal glucose concentrations." Diabetologia.
                         2015;58(1):43-9. PMID:25297572
    Acarbose:            Chen Z, Fu X, Kuang J, Chen J, Chen H, Pei J, Yang H.
                         "Single-dose acarbose decreased glucose-dependent
                         insulinotropic peptide and glucagon levels in Chinese
                         patients with newly diagnosed type 2 diabetes
                         mellitus after a mixed meal." BMC Endocr Disord.
                         2016. PMID:27686734
    Exogenous insulin:   Howey DC, Bowsher RR, Brunelle RL, Woodworth JR.
                         "[Lys(B28), Pro(B29)]-human insulin. A rapidly
                         absorbed analogue of human insulin." Diabetes.
                         1994;43(3):396-402. PMID:8314011
    Glipizide:           Groop LC et al. Acta Diabetol. 1991;28(2):162-8.
                         PMID:1777653; Groop L et al. Diabetes.
                         1987;36(11):1320-8. PMID:3117609
    Tirzepatide:         Furihata K, Mimura H, Urva S, Oura T, Ohwaki K,
                         Imaoka T. "A phase 1 multiple-ascending dose
                         study of tirzepatide in Japanese participants
                         with type 2 diabetes." Diabetes Obes Metab.
                         2022;24(2):239-246. PMID:34647404
    Pioglitazone:        Miyazaki Y, Matsuda M, DeFronzo RA. "Dose-response
                         effect of pioglitazone on insulin sensitivity and
                         insulin secretion in type 2 diabetes." Diabetes
                         Care. 2002;25(3):517-23. PMID:11874940

Package layout (split from one 1659-line module, 2026-10-20; every name
is re-exported here, so `from biosim_physiology.drug_interactions import X`
keeps working):
    _core               Drug, DrugDose
    classic             caffeine, metformin, glucocorticoid, alcohol
    incretins           tirzepatide, GLP-1 class
    insulin_sensitivity pioglitazone, indinavir, HCTZ, growth hormone,
                        atorvastatin, hydroxychloroquine (static p3 factors)
    sulfonylureas       glucose-independent secretion
    alpha_glucosidase   acarbose/miglitol meal kabs scaling
    ode_coupled         SGLT2 inhibitors, exogenous insulin
    registry            dispatch table, drug-class sets, combination checks
"""

from __future__ import annotations

from ._core import (  # noqa: F401 -- re-exported, see module docstring
    Drug,
    DrugDose,
)
from .classic import (  # noqa: F401 -- re-exported, see module docstring
    caffeine_modifiers,
    metformin_modifiers,
    glucocorticoid_modifiers,
    alcohol_modifiers,
)
from .incretins import (  # noqa: F401 -- re-exported, see module docstring
    TIRZEPATIDE_GB_DOSE_MG,
    TIRZEPATIDE_GB_MULTIPLIER,
    tirzepatide_gb_factor,
    _glp1_agonist_modifiers,
)
from .insulin_sensitivity import (  # noqa: F401 -- re-exported, see module docstring
    PIOGLITAZONE_P3_DOSE_MG,
    PIOGLITAZONE_P3_MULTIPLIER,
    pioglitazone_p3_factor,
    pioglitazone_modifiers,
    INDINAVIR_TESTED_DOSE_MG,
    INDINAVIR_P3_FACTOR_AT_TESTED_DOSE,
    indinavir_p3_factor,
    indinavir_modifiers,
    THIAZIDE_TESTED_DOSE_MG,
    THIAZIDE_P3_FACTOR_AT_TESTED_DOSE,
    hydrochlorothiazide_p3_factor,
    hydrochlorothiazide_modifiers,
    GROWTH_HORMONE_P3_FACTOR,
    growth_hormone_modifiers,
    ATORVASTATIN_P3_DOSE_MG,
    ATORVASTATIN_P3_FACTOR_AT_DOSE,
    atorvastatin_p3_factor,
    atorvastatin_modifiers,
    HYDROXYCHLOROQUINE_P3_FACTOR,
    hydroxychloroquine_modifiers,
)
from .ode_coupled import (  # noqa: F401 -- re-exported, see module docstring
    _sglt2_inhibitor_modifiers,
    _exogenous_insulin_modifiers,
)
from .sulfonylureas import (  # noqa: F401 -- re-exported, see module docstring
    _sulfonylurea_modifiers,
    SULFONYLUREA_SECRETION_FRACTION_OF_KE_IB,
    GLIPIZIDE_SECRETION_FRACTION_OF_KE_IB,
    GLIMEPIRIDE_SECRETION_FRACTION_OF_KE_IB,
    _SULFONYLUREA_SECRETION_FRACTION,
    glucose_independent_secretion_pM_per_min,
)
from .alpha_glucosidase import (  # noqa: F401 -- re-exported, see module docstring
    _alpha_glucosidase_inhibitor_modifiers,
    ACARBOSE_100MG_KABS_SCALE,
    acarbose_kabs_scale,
    MIGLITOL_50MG_KABS_SCALE,
    miglitol_kabs_scale,
)
from .registry import (  # noqa: F401 -- re-exported, see module docstring
    _DRUG_FUNCTIONS,
    GLP1_DRUGS,
    SGLT2_DRUGS,
    EXOGENOUS_INSULIN_DRUGS,
    SULFONYLUREA_DRUGS,
    SINGLE_INSTANCE_DRUGS,
    DOSE_SCALED_SINGLE_EVENT_DRUGS,
    validate_drug_combination,
    ALPHA_GLUCOSIDASE_INHIBITOR_DRUGS,
    modifiers_for_dose,
    drug_modifiers,
)
