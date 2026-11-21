# SAVE TO: packages\data-seeds\src\biosim_seeds\context_seed.py
"""
Real non-molecular CONTEXT entity registrations for ContextRegistry.

Deliberately separate from registry_seed.py (which seeds
biosim_id_resolver.IDRegistry, molecular entities only) -- see
biosim_context_registry.types's own module docstring for the full
architectural reasoning (project decision, 2026-08-18).

Changelog:
    2026-08-18  Initial seed: chronic_kidney_disease (DISEASE_STATE),
                sleep_restriction (PHYSIOLOGICAL_AXIS) -- the subject
                entities for physiology/modifiers.py's own CKD (Track 1
                row 1.24) and sleep-restriction (row 1.26) constants,
                which had no molecular entity to attach a ParameterStore
                subject to (unlike every drug/hormone chunk in
                registry_seed.py/parameter_seed.py's own campaign).
    2026-08-18  Added gck_mody (DISEASE_STATE) -- MASTER_PLAN 5.1 chunk 13.
    2026-08-18  Added gestational_diabetes (DISEASE_STATE) -- MASTER_PLAN
                5.1 chunk 14.
    2026-08-18  Added type_1_diabetes (DISEASE_STATE) -- MASTER_PLAN 5.1
                chunk 15.
    2026-08-18  Added type_2_diabetes, insulin_resistant, hypothyroidism
                (DISEASE_STATE) -- MASTER_PLAN 5.1 chunk 16, all with a
                real, disclosed citation-magnitude gap (see each entity's
                own description below).
    2026-08-18  Added intermittent_hypoxia (PHYSIOLOGICAL_AXIS) --
                MASTER_PLAN 5.1 chunk 22, curating Track 1 row 1.32's
                own p1/p3 factors.
    2026-08-19  Added cold_acclimation (PHYSIOLOGICAL_AXIS) --
                MASTER_PLAN 5.1 chunk 29, curating Track 1 row 1.40's
                own p3 factor.
    2026-08-19  Added menstrual_cycle (PHYSIOLOGICAL_AXIS) --
                MASTER_PLAN 5.1 chunk 30, curating Track 1 row 1.36's
                own luteal-phase p3 factor.
    2026-08-20  Added pcos (DISEASE_STATE) -- MASTER_PLAN 5.1 chunk 31,
                curating Track 1 row 1.29's own p3 factor.
    2026-08-20  Added acromegaly (DISEASE_STATE) -- MASTER_PLAN 5.1
                chunk 32, curating Track 1 row 1.37's own Gb factor.
    2026-08-20  Added osa (DISEASE_STATE) -- MASTER_PLAN 5.1 chunk 33,
                curating Track 1 row 1.39's own p3 factor.
    2026-08-20  Added hyperthyroidism (DISEASE_STATE) -- MASTER_PLAN
                5.1 chunk 35, curating the EGPParameters.
                for_hyperthyroidism() hepatic egp_basal anchor.
    2026-10-19  Added nafld (DISEASE_STATE) -- MASTER_PLAN 5.1 chunk 38,
                curating Track 1 row 1.45's own p3 factor.
"""

from __future__ import annotations

from biosim_context_registry import ContextEntityType, ContextRegistry


def seed_context_registry(registry: ContextRegistry) -> ContextRegistry:
    """Populate a ContextRegistry with the platform's real context entities."""
    _seed_disease_states(registry)
    _seed_physiological_axes(registry)
    return registry


def create_seeded_context_registry() -> ContextRegistry:
    """Create and return a new ContextRegistry seeded with all v1 context entities."""
    return seed_context_registry(ContextRegistry())


def _seed_disease_states(registry: ContextRegistry) -> None:
    registry.register(
        "chronic_kidney_disease",
        ContextEntityType.DISEASE_STATE,
        description=(
            "Moderate-severe, non-diabetic chronic kidney disease "
            "(DiseaseState.CHRONIC_KIDNEY_DISEASE, MASTER_PLAN row 1.24) "
            "-- reduces insulin sensitivity (p3) and clearance (ke), "
            "confirmed no CKD-specific secretion deficit. See "
            "physiology/modifiers.py's own citation comment "
            "(de Boer et al. 2016, PMID:26823551)."
        ),
    )
    registry.register(
        "gck_mody",
        ContextEntityType.DISEASE_STATE,
        description=(
            "Glucokinase (GCK, MODY2) mutations (DiseaseState.GCK_MODY, "
            "MASTER_PLAN row 1.11) -- a pure beta-cell glucose-sensing "
            "threshold reset (elevated fasting glucose only), confirmed "
            "no effect on p1/p3 (peripheral insulin action) or "
            "secretion_scale (maximal secretory capacity). See "
            "physiology/modifiers.py's own citation comment "
            "(Guenat et al. 2000, PMID:11119017)."
        ),
    )
    registry.register(
        "gestational_diabetes",
        ContextEntityType.DISEASE_STATE,
        description=(
            "Gestational diabetes (DiseaseState.GESTATIONAL_DIABETES, "
            "MASTER_PLAN row 1.11) -- folds BOTH normal-pregnancy insulin "
            "resistance and GDM's own additional decompensation into one "
            "preset (this platform has no separate 'pregnant, not "
            "diabetic' baseline). Reduces p1/p3, raises secretion_scale "
            "(compensatory hypersecretion) and Gb. See physiology/"
            "modifiers.py's own citation comment (Kautzky-Willer et al. "
            "1997, PMID:9353615 for p1/p3/secretion_scale; IADPSG "
            "Consensus Panel 2010, PMID:20190296 for the Gb diagnostic "
            "threshold)."
        ),
    )
    registry.register(
        "type_1_diabetes",
        ContextEntityType.DISEASE_STATE,
        description=(
            "Autoimmune beta-cell destruction (DiseaseState.TYPE_1_DIABETES, "
            "Track 1 row 1.10) -- near-total, not reduced, endogenous "
            "insulin secretion loss. Ib/secretion_scale factors are "
            "disclosed ENGINEERING CHOICES (chosen to reach this "
            "platform's own numeric floors after clamping), motivated by "
            "but not directly equal to a real qualitative citation -- see "
            "physiology/modifiers.py's own citation comment (Simonson et "
            "al. 1987, PMID:3100365, undetectable C-peptide in confirmed "
            "IDDM patients). Gb is a THIRD, separate engineering choice "
            "with NO citation at all (no consensus fasting-glucose value "
            "found) -- registered anyway with pubmed_id='UNVERIFIED' "
            "rather than omitted, so this gap is visible in the "
            "provenance store itself."
        ),
    )
    registry.register(
        "type_2_diabetes",
        ContextEntityType.DISEASE_STATE,
        description=(
            "Established T2D (DiseaseState.TYPE_2_DIABETES) -- substantial "
            "beta-cell loss and insulin resistance. Stumvoll et al. 2005 "
            "(PMID:15823385, a Lancet narrative review) is a real, "
            "topically-relevant citation supporting the DIRECTION of this "
            "preset (progressive beta-cell decline, early insulin "
            "insensitivity) but reports NO specific numeric magnitude for "
            "secretion_scale/p3/Gb -- confirmed via a dedicated scoping "
            "pass (live-fetched abstract, 997 chars, zero quantitative "
            "figures). Registered BRONZE (a real, motivating citation "
            "exists, same category as T1D's Ib/secretion_scale factors) "
            "rather than UNVERIFIED (which is reserved for a subject with "
            "NO citation attempt at all, e.g. this platform's own "
            "insulin_resistant entity below). The p1 factor (added "
            "2026-10-23, MASTER_PLAN row 1.46 (a)) is separately GOLD: "
            "Basu et al. 2009 (PMID:19196896) reports SG 0.019 vs. 0.034 "
            "min-1 in T2D vs. matched controls."
        ),
    )
    registry.register(
        "insulin_resistant",
        ContextEntityType.DISEASE_STATE,
        description=(
            "Pre-diabetes/metabolic syndrome (DiseaseState.INSULIN_RESISTANT) "
            "-- reduced insulin sensitivity with compensatory "
            "hypersecretion. UNLIKE every other disease state in this "
            "campaign, this preset has NO citation anywhere in "
            "physiology/modifiers.py at all -- not even in the module's "
            "own top-of-file citation list, confirmed via a dedicated "
            "scoping pass. All three factors registered with "
            "pubmed_id='UNVERIFIED', matching this project's own "
            "established convention (see type_1_diabetes's Gb_factor) "
            "for a real, load-bearing constant with zero literature "
            "backing, so this gap is visible in the provenance store "
            "rather than only in the absence of a code comment."
        ),
    )
    registry.register(
        "hypothyroidism",
        ContextEntityType.DISEASE_STATE,
        description=(
            "Overt hypothyroidism (DiseaseState.HYPOTHYROIDISM) -- reduced "
            "peripheral glucose disposal/insulin action (p1, p3). A "
            "genuine, UNRESOLVED three-way citation conflict exists, "
            "disclosed here and in physiology/modifiers.py's own comment, "
            "not resolved by fiat: Rochon et al. 2003 (PMID:12519082) and "
            "Handisurya et al. 2008 (PMID:18429948) directionally support "
            "reduced insulin action, but Shah et al. 1975 (PMID:1175861) "
            "found NO insulin resistance in hypothyroidism at all -- a "
            "direct glucose-nadir measurement that, if anything, suggests "
            "INCREASED sensitivity. None of the three studies reports a "
            "specific magnitude for p1/p3 either way (confirmed via a "
            "dedicated scoping pass of all three abstracts) -- p1=0.85/"
            "p3=0.88 are registered BRONZE, citing the two supporting "
            "studies, with the Shah 1975 contradiction disclosed directly "
            "in each parameter's own notes."
        ),
    )
    registry.register(
        "hyperthyroidism",
        ContextEntityType.DISEASE_STATE,
        description=(
            "Hyperthyroidism (DiseaseState.HYPERTHYROIDISM) -- elevated "
            "basal hepatic glucose production (EGP), represented SOLELY "
            "via the ODE-coupled HepGlucPathway "
            "(EGPParameters.for_hyperthyroidism()), NOT a static Gb/p3 "
            "modifier -- physiology/modifiers.py's own disease_modifiers() "
            "deliberately returns {} for this state to avoid double-"
            "counting the identical mechanism (same discipline already "
            "applied to glucocorticoid/CortisolPathway). Peripheral "
            "insulin-stimulated glucose disposal is reported as normal-or-"
            "increased (Dimitriadis & Raptis 2001, PMID:11460573), so no "
            "p3 modifier is added either. See biosim_ode.hepatic's own "
            "EGPParameters.for_hyperthyroidism() citation comment "
            "(Wennlund et al. 1986, PMID:3510000)."
        ),
    )
    registry.register(
        "nafld",
        ContextEntityType.DISEASE_STATE,
        description=(
            "Non-alcoholic fatty liver disease, biopsy-proven, "
            "non-diabetic (DiseaseState.NAFLD, Track 1 row 1.45) -- "
            "reduced peripheral glucose disposal during clamp, modelled "
            "here as a p3 reduction (an inference, BRONZE); measured endogenous "
            "glucose production was NORMAL, so Gb/ke/secretion_scale and "
            "hepatic EGP are untouched (real hepatic insulin resistance "
            "reported in the full text is NOT modelled). The 0.70 anchor "
            "(the low-dose clamp step; the high-dose step gives 0.55) is a "
            "disclosed JUDGMENT CALL, not a conservative bound. Composed "
            "with INSULIN_RESISTANT/TYPE_2_DIABETES it is probably partial "
            "double-counting of one underlying metabolic-syndrome state. "
            "See physiology/modifiers.py's own citation comment "
            "(Bugianesi et al. 2005, PMID:15747110)."
        ),
    )
    registry.register(
        "pcos",
        ContextEntityType.DISEASE_STATE,
        description=(
            "Polycystic ovary syndrome (DiseaseState.PCOS, Track 1 row "
            "1.29) -- peripheral insulin resistance (p3), mechanistically "
            "distinct from beta-cell secretion loss or hepatic glucose "
            "overproduction (both explicitly confirmed NULL by the same "
            "sources). A GENUINE, UNRESOLVED literature conflict exists "
            "on whether the resistance is obesity-independent (Dunaif "
            "1989's own explicit conclusion) or obesity-gated (Ketel "
            "2011's own lean-vs-lean comparison finds no significant "
            "difference) -- disclosed here and in physiology/"
            "modifiers.py's own comment, not resolved by fiat. The "
            "registered factor is Ketel's own OBESE-stratum ratio, "
            "applied regardless of the composed profile's own BMI (a "
            "real, disclosed limitation -- this platform has no separate "
            "lean/obese PCOS split). ke/Gb/hepatic glucose production "
            "confirmed NULL. See physiology/modifiers.py's own citation "
            "comment (Ketel et al. 2011, PMID:21920943)."
        ),
    )
    registry.register(
        "acromegaly",
        ContextEntityType.DISEASE_STATE,
        description=(
            "Chronic growth-hormone/IGF-1 excess (DiseaseState.ACROMEGALY, "
            "Track 1 row 1.37) -- Gb-ONLY, an honest PARTIAL "
            "implementation, not an assumed-null claim: the source "
            "confirms qualitatively that active acromegaly causes "
            "substantial insulin resistance during a clamp, but reports "
            "no extractable glucose-infusion-rate/M-value number for "
            "p3/ke/secretion_scale, so those fields are deliberately "
            "left untouched rather than guessed at. Represents ACTIVE, "
            "untreated acromegaly only (the source's own pre-operative "
            "measurement) -- post-adenomectomy normalization is reported "
            "only qualitatively, not used. See physiology/modifiers.py's "
            "own citation comment (Moller et al. 1992, PMID:1569148)."
        ),
    )
    registry.register(
        "osa",
        ContextEntityType.DISEASE_STATE,
        description=(
            "Obstructive sleep apnea as a CHRONIC clinical diagnosis "
            "(DiseaseState.OSA, Track 1 row 1.39) -- reduces insulin "
            "sensitivity (p3). Deliberately distinct from "
            "intermittent_hypoxia (an ACUTE, wakefulness-only exposure "
            "at a specific measured desaturation rate) -- the two are "
            "mutually exclusive at the profile.py level since both lower "
            "p3 via the same real mechanism and would silently double-"
            "count if composed. A GENUINE, UNRESOLVED literature conflict "
            "exists, disclosed here and in physiology/modifiers.py's own "
            "comment, not resolved by fiat: Harsch et al. 2004 found CPAP "
            "treatment measurably improves insulin sensitivity (used "
            "here), but Smurra et al. 2001 found no change in "
            "glucose-insulin metabolism after CPAP in a smaller cohort -- "
            "both studies' own conclusions name OBESITY as a likely "
            "confound/effect-modifier, and this platform composes OSA's "
            "factor with bmi_modifiers() regardless (no obesity-gating "
            "mechanism exists for any disease state, same 'disclosed, not "
            "solved' precedent as PCOS's own obesity confound). ke/Gb/"
            "secretion_scale confirmed unreported by either source. See "
            "physiology/modifiers.py's own citation comment (Harsch et "
            "al. 2004, PMID:14512265)."
        ),
    )


def _seed_physiological_axes(registry: ContextRegistry) -> None:
    registry.register(
        "sleep_restriction",
        ContextEntityType.PHYSIOLOGICAL_AXIS,
        description=(
            "Acute sleep restriction (PhysiologicalProfile.sleep_hours, "
            "MASTER_PLAN row 1.26) -- reduces insulin sensitivity (p3), "
            "confirmed no effect on secretion. See physiology/"
            "modifiers.py's own citation comment (Buxton et al. 2010, "
            "PMID:20585000)."
        ),
    )
    registry.register(
        "intermittent_hypoxia",
        ContextEntityType.PHYSIOLOGICAL_AXIS,
        description=(
            "Acute intermittent hypoxia during wakefulness "
            "(PhysiologicalProfile.desaturation_events_per_hour, Track 1 "
            "row 1.32) -- reduces both insulin sensitivity (p3) and "
            "glucose effectiveness (p1), confirmed no effect on "
            "secretion. Genuinely orthogonal to sleep_restriction above "
            "(this source induces hypoxia DURING WAKEFULNESS, "
            "deliberately isolated from sleep loss/fragmentation) and to "
            "CortisolPathway/stress_level (directly measured cortisol "
            "null in the same study). See physiology/modifiers.py's own "
            "citation comment (Louis & Punjabi 2009, PMID:19265062)."
        ),
    )
    registry.register(
        "cold_acclimation",
        ContextEntityType.PHYSIOLOGICAL_AXIS,
        description=(
            "Completed chronic (~10-day) cold-acclimation protocol "
            "(PhysiologicalProfile.cold_acclimatized, Track 1 row 1.40) "
            "-- raises insulin sensitivity (p3). T2D-ONLY evidence base, "
            "explicitly disclosed, not silently generalized to healthy "
            "subjects (the source cohort had no healthy-subject arm; "
            "runner.py emits a SimulationOutput.warnings entry whenever "
            "cold_acclimatized is set without DiseaseState."
            "TYPE_2_DIABETES present). See physiology/modifiers.py's own "
            "citation comment (Hanssen et al. 2015, PMID:26147760)."
        ),
    )
    registry.register(
        "menstrual_cycle",
        ContextEntityType.PHYSIOLOGICAL_AXIS,
        description=(
            "Menstrual cycle phase (PhysiologicalProfile."
            "menstrual_cycle_phase, Track 1 row 1.36) -- reduces insulin "
            "sensitivity (p3) at LUTEAL only; FOLLICULAR and MIDCYCLE are "
            "both neutral (the source's own conclusion explicitly found "
            "no significant midcycle-vs-follicular difference, despite a "
            "numerically lower point estimate -- a categorical, not "
            "continuous, modifier by design). Confirmed no effect on "
            "glucose effectiveness (p1). Distinct from GESTATIONAL_"
            "DIABETES (a disease state, not normal cycle physiology) and "
            "from PCOS (a distinct pathological insulin-resistance "
            "mechanism, not cycle-phase-dependent). See physiology/"
            "modifiers.py's own citation comment (Valdes & Elkind-Hirsch "
            "1991, PMID:1997519)."
        ),
    )
