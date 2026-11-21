# SAVE TO: packages\data-seeds\src\biosim_seeds\parameter_seed\diseases.py
"""Non-diabetes disease states."""

from __future__ import annotations

from biosim_context_registry import ContextEntityType, ContextRegistry
from biosim_id_resolver import IDRegistry
from biosim_parameter_db import (
    ParameterKind,
    ParameterStore,
)

from ._helpers import SeedParameterBuilder, _resolve, _resolve_context

# ---------------------------------------------------------------------------
# Cushing's syndrome cortisol ratio (MASTER_PLAN 5.1, chunk 34; Track 1
# row 1.38) -- registers cortisol_pathway.py's own
# CUSHINGS_CORTISOL_RATIO against the already-registered "cortisol"
# subject entity (registry_seed.py, chunk 6). Same subject as
# _seed_cortisol_pk above, a third CUSTOM entry under it.
#
# This citation was previously stalled: its own abstract (confirmed
# live via PubMedFetcher) is entirely about mitochondrial/OXPHOS gene
# expression, never mentioning cortisol at all -- the claimed 822/527
# figures live in the paper's Table 1, only reachable via full text.
# Verified live via the newly-built biosim_literature.verify_citation()
# tool (MASTER_PLAN 10, item L1) specifically to unblock this chunk:
# `python -m biosim_literature.verification 31695455 822 527` confirms
# both values appear in the fetched full text at "Plasma cortisol at
# 08:00 hrs (nmol/L) 527+/-51 822+/-80* 286+/-71" (controls / active
# Cushing's / treated).
# ---------------------------------------------------------------------------

def _seed_cushings_cortisol_ratio(store: ParameterStore, registry: IDRegistry) -> None:
    """
    Registers the Cushing's syndrome cortisol ratio. GOLD -- a direct
    ratio computed from a single primary source's own two reported
    Table 1 morning-cortisol values (822/527), no unit conversion --
    cross-checked live against the real
    CortisolParameters.CUSHINGS_CORTISOL_RATIO constant, not a
    duplicated literal.
    """
    cortisol_id = _resolve(registry, "cortisol")

    from biosim_ode.cortisol_pathway import CortisolParameters

    # Jezkova J, Durovcova V, Wenchich L, et al. "The relationship of
    # mitochondrial dysfunction and the development of insulin
    # resistance in Cushing's syndrome." Diabetes Metab Syndr Obes.
    # 2019;12:1459-1471. PMID:31695455 -- re-verified live this chunk
    # via the new full-text-capable verify_citation() tool (this
    # paper's own abstract never mentions cortisol; the claim is only
    # confirmable via full text, which the tool now checks). Table 1,
    # morning (08:00) plasma cortisol, same RIA assay: controls
    # 527+/-51 nmol/L (n=24) vs. active Cushing's 822+/-80 nmol/L
    # (n=24, P<0.05).
    store.add(
        SeedParameterBuilder()
        .for_subject(cortisol_id)
        .kind(ParameterKind.CUSTOM, custom_name="cushings_cortisol_ratio")
        .value(
            CortisolParameters.CUSHINGS_CORTISOL_RATIO,
            unit="dimensionless (ratio of morning plasma cortisol, active Cushing's vs. controls)",
        )
        .from_pubmed("31695455", figure_or_table="Table 1")
        .in_human()
        .clinical()
        .gold()
        .notes("Direct ratio of the paper's own two reported Table 1 "
               "morning-cortisol values (822/527) -- computed live "
               "against the real CortisolParameters."
               "CUSHINGS_CORTISOL_RATIO constant, not a duplicated "
               "literal. Only reachable via the paper's full text, not "
               "its abstract (which discusses only mitochondrial/OXPHOS "
               "gene expression) -- verified via the new full-text-"
               "capable biosim_literature.verify_citation() tool, built "
               "specifically to unblock this chunk (MASTER_PLAN 10, "
               "item L1). This is the ACHIEVED/OBSERVED cortisol ratio "
               "in real patients (Jezkova's own directly-measured "
               "figure), used as the TARGET steady state this "
               "platform's coupled CRH/cortisol solver is inverse-"
               "solved to reach -- see cortisol_pathway.py's own "
               "for_cushings() classmethod docstring and "
               "solve_cortisol_basal_for_target_steady_state() for the "
               "full derivation. Independently confirmed by direct "
               "execution this chunk: for_cushings(stress_level=0.0)'s "
               "own ACHIEVED cortisol_ss/CORTISOL_BASAL_NM ratio is "
               "1.5598, matching this literature ratio (1.5598) "
               "essentially exactly, by construction of the inverse "
               "solve -- the ~1.287x figure elsewhere in "
               "cortisol_pathway.py's own comments is what the naive, "
               "explicitly REJECTED approach (setting cortisol_basal "
               "directly to CORTISOL_BASAL_NM*1.560, without the "
               "inverse solve) would achieve instead; it is NOT what "
               "the current for_cushings() implementation produces.")
        .build()
    )


# ---------------------------------------------------------------------------
# Hyperthyroidism hepatic EGP anchor (MASTER_PLAN 5.1, chunk 35) --
# registers biosim_ode.hepatic's own EGPParameters.for_hyperthyroidism()
# egp_basal anchor against a new "hyperthyroidism" CONTEXT entity
# (context_seed.py). Same category as chronic_kidney_disease (chunk
# 11), pcos (chunk 31) -- a disease state, not a molecule. Unlike
# every other Track 5 disease-state chunk, the underlying mechanism
# here is a dynamic ODE-coupled pathway (HepGlucPathway), not a
# static physiology/modifiers.py factor -- physiology/modifiers.py's
# own disease_modifiers() deliberately returns {} for HYPERTHYROIDISM
# specifically to avoid double-counting this same mechanism.
# ---------------------------------------------------------------------------

def _seed_hyperthyroidism_egp(store: ParameterStore, context_registry: ContextRegistry) -> None:
    """
    Registers the hyperthyroidism hepatic egp_basal anchor factor.
    GOLD -- a directly-reported within-study percentage elevation, no
    unit conversion -- cross-checked live against the real
    EGPParameters.for_hyperthyroidism() output relative to
    EGP_BASAL_MM_MIN, not a duplicated literal.
    """
    hyperthyroidism_id = _resolve_context(
        context_registry, "hyperthyroidism", ContextEntityType.DISEASE_STATE,
    )

    from biosim_ode.hepatic import EGP_BASAL_MM_MIN, EGPParameters

    # Wennlund A, Felig P, Hagenfeldt L, Wahren J. "Hepatic glucose
    # production and splanchnic glucose exchange in hyperthyroidism."
    # J Clin Endocrinol Metab. 1986. PMID:3510000 -- verified live via
    # NCBI efetch this chunk. n=8 hyperthyroid vs. 10 normal subjects,
    # hepatic venous catheterization + [3H]3-glucose isotope dilution.
    # Directly quoted: "Basal HGP... was increased by 20% in the
    # hyperthyroid patients [14.2 +/- 0.6 (SEM) mumol/kg.min] as
    # compared to normal subjects (11.9 +/- 0.6, P less than 0.02)."
    egp_basal_factor = EGPParameters.for_hyperthyroidism().egp_basal / EGP_BASAL_MM_MIN
    store.add(
        SeedParameterBuilder()
        .for_subject(hyperthyroidism_id)
        .kind(ParameterKind.CUSTOM, custom_name="egp_basal_factor")
        .value(egp_basal_factor, unit="dimensionless (multiplier on HepGlucPathway's EGP_BASAL_MM_MIN)")
        .from_pubmed("3510000", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .gold()
        .notes("Directly-reported +20% basal hepatic glucose production "
               "elevation -- computed live against the real "
               "EGPParameters.for_hyperthyroidism() output, not a "
               "duplicated literal. This is the SOLE representation of "
               "hyperthyroidism's glucose effect on this platform -- "
               "physiology/modifiers.py's own disease_modifiers() "
               "deliberately returns {} for DiseaseState.HYPERTHYROIDISM "
               "(no static Gb/p3 modifier), since composing both would "
               "double-count the identical mechanism (same discipline "
               "already applied to glucocorticoid/CortisolPathway). "
               "ic50_pM is re-derived from this anchor via "
               "_solve_egp_suppression_ic50(), not independently cited "
               "-- egp_max/t2d_egp_scale stay at healthy defaults, since "
               "Wennlund reports no maximum-suppression-capacity defect, "
               "only a dose-response shift.")
        .build()
    )


# ---------------------------------------------------------------------------
# NAFLD p3 constant (MASTER_PLAN 5.1, chunk 38; Track 1 row 1.45) --
# registers physiology/modifiers.py's own NAFLD disease-state p3 factor
# against a new "nafld" CONTEXT entity (context_seed.py). Same category as
# pcos (chunk 31) -- a disease state, not a molecule.
# ---------------------------------------------------------------------------

def _seed_nafld(store: ParameterStore, context_registry: ContextRegistry) -> None:
    """
    Registers the NAFLD p3 factor (0.70). BRONZE -- the source directly
    reports a peripheral glucose-disposal deficit at fixed clamp insulin
    doses, not a p3-equivalent insulin-sensitivity ratio, so mapping it to
    a p3 multiplier is an inference (and the choice of WHICH reported step
    to use is a further judgment call), both disclosed in the notes.
    Cross-checked live against the real disease_modifiers() function.
    """
    nafld_id = _resolve_context(
        context_registry, "nafld", ContextEntityType.DISEASE_STATE,
    )

    from biosim_physiology.modifiers import disease_modifiers
    from biosim_physiology.profile import DiseaseState

    # Bugianesi E, Gastaldelli A, Vanni E, Gambino R, Cassader M, Baldi S,
    # Ponti V, Pagano G, Ferrannini E, Rizzetto M. "Insulin resistance in
    # non-diabetic patients with non-alcoholic fatty liver disease: sites
    # and mechanisms." Diabetologia. 2005. PMID:15747110 -- verified live
    # via NCBI efetch. Two-step euglycaemic clamp, n=12 biopsy-proven NAFLD
    # vs. n=6 controls. Directly quoted: "endogenous glucose production
    # (basal and during the clamp) was normal; however, peripheral glucose
    # disposal was markedly decreased (by 30% and 45% at the low and high
    # insulin doses, respectively, p<0.0001)".
    reference_factors = disease_modifiers([DiseaseState.NAFLD])
    store.add(
        SeedParameterBuilder()
        .for_subject(nafld_id)
        .kind(ParameterKind.CUSTOM, custom_name="p3_factor")
        .value(reference_factors["p3"], unit="dimensionless (multiplier on PathwayParameters.p3)")
        .from_pubmed("15747110", figure_or_table="abstract")
        .in_human()
        .clinical()
        .bronze()
        .notes("TIER NOTE: BRONZE, not GOLD -- the paper directly reports "
               "a peripheral glucose-DISPOSAL deficit at a given clamp "
               "insulin dose, not a minimal-model insulin-sensitivity (p3) "
               "ratio; reading that deficit as a p3 multiplier is an "
               "inference (disposal at unequal insulin levels is not "
               "insulin action per unit insulin). "
               "Low-insulin-dose-step peripheral glucose-disposal deficit "
               "(30%) -- computed live against the real disease_modifiers() "
               "function, not a duplicated literal. The step choice is a "
               "JUDGMENT CALL, not a conservative bound: the high-dose "
               "step (45%, factor 0.55) is not used, and a reviewer's "
               "reading of the paper's tables (not re-verified here) "
               "suggests insulin-normalised sensitivity is unresolved "
               "between ~1.0 and ~0.37. EGP was measured NORMAL (hepatic "
               "EGP untouched), but real hepatic insulin resistance in the "
               "full text is NOT modelled. Probably partial double-counting "
               "with INSULIN_RESISTANT (uncited 0.70) / TYPE_2_DIABETES "
               "when composed.")
        .build()
    )


# ---------------------------------------------------------------------------
# PCOS p3 constant (MASTER_PLAN 5.1, chunk 31; Track 1 row 1.29) --
# registers physiology/modifiers.py's own PCOS disease-state p3 factor
# against a new "pcos" CONTEXT entity (context_seed.py). Same category
# as chronic_kidney_disease (chunk 11) -- a disease state, not a
# molecule.
# ---------------------------------------------------------------------------

def _seed_pcos(store: ParameterStore, context_registry: ContextRegistry) -> None:
    """
    Registers the PCOS p3 factor. GOLD -- a direct ratio computed from
    a single primary source's own two reported M/I clamp values
    (0.5/0.8), no unit conversion -- cross-checked live against the
    real disease_modifiers() function, not a duplicated literal.
    """
    pcos_id = _resolve_context(
        context_registry, "pcos", ContextEntityType.DISEASE_STATE,
    )

    from biosim_physiology.modifiers import disease_modifiers
    from biosim_physiology.profile import DiseaseState

    # Ketel IJ, Stehouwer CD, Serne EH, et al. "Insulin-induced
    # capillary recruitment is impaired in both lean and obese women
    # with PCOS." Hum Reprod. 2011;26(11):3130-7. PMID:21920943 --
    # verified live via NCBI efetch this chunk. n=36 PCOS (20 lean, 16
    # obese) vs. n=27 matched controls, isoglycemic-hyperinsulinemic
    # clamp. Directly quoted: obese PCOS vs. obese controls M/I "0.5
    # (0.2-1.1) versus 0.8 (0.3-1.4) ... P < 0.01" (ratio 0.625,
    # significant); lean PCOS vs. lean controls M/I "1.5 (0.5-2.6)
    # versus 1.7 (1.0-3.7) ... P = 0.17" (ratio 0.882, NOT significant).
    reference_factors = disease_modifiers([DiseaseState.PCOS])
    store.add(
        SeedParameterBuilder()
        .for_subject(pcos_id)
        .kind(ParameterKind.CUSTOM, custom_name="p3_factor")
        .value(reference_factors["p3"], unit="dimensionless (multiplier on PathwayParameters.p3)")
        .from_pubmed("21920943", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .gold()
        .notes("Obese-stratum M/I ratio (0.5/0.8) -- computed live "
               "against the real disease_modifiers() function, not a "
               "duplicated literal. A GENUINE, UNRESOLVED literature "
               "conflict exists on whether this resistance is obesity-"
               "independent (Dunaif 1989, PMID:2670645) or obesity-"
               "gated (this same source's own lean-vs-lean comparison "
               "finds no significant difference, P=0.17) -- disclosed, "
               "not resolved by fiat. Applied regardless of the "
               "composed profile's own BMI -- this platform has no "
               "separate lean/obese PCOS split, a real, disclosed "
               "limitation. ke/Gb/hepatic glucose production confirmed "
               "NULL by the same/adjacent sources, not registered.")
        .build()
    )


# ---------------------------------------------------------------------------
# Acromegaly Gb constant (MASTER_PLAN 5.1, chunk 32; Track 1 row 1.37)
# -- registers physiology/modifiers.py's own acromegaly disease-state
# Gb factor against a new "acromegaly" CONTEXT entity (context_seed.py).
# Same category as chronic_kidney_disease (chunk 11) and pcos (chunk
# 31) -- a disease state, not a molecule.
# ---------------------------------------------------------------------------

def _seed_acromegaly(store: ParameterStore, context_registry: ContextRegistry) -> None:
    """
    Registers the acromegaly Gb factor. GOLD -- a direct ratio
    computed from a single primary source's own two reported plasma
    glucose values (5.9/5.2), no unit conversion -- cross-checked live
    against the real disease_modifiers() function, not a duplicated
    literal. Gb-only, an honest partial implementation: p3/ke/
    secretion_scale have no extractable magnitude in the source and
    are deliberately left unregistered, not guessed at.
    """
    acromegaly_id = _resolve_context(
        context_registry, "acromegaly", ContextEntityType.DISEASE_STATE,
    )

    from biosim_physiology.modifiers import disease_modifiers
    from biosim_physiology.profile import DiseaseState

    # Moller N, Schmitz O, Jorgensen JO, Astrup J, Bak JF, Christensen
    # SE, Alberti KG, Weeke J. "Basal- and insulin-stimulated substrate
    # metabolism in patients with active acromegaly before and after
    # adenomectomy." J Clin Endocrinol Metab. 1992;74(5):1012-9.
    # PMID:1569148 -- verified live via NCBI efetch this chunk. n=6
    # active (pre-adenomectomy) acromegaly vs. n=7 controls. Directly
    # quoted: "increased plasma glucose (5.9 +/- 0.2 vs. 5.2 +/- 0.2
    # mmol/L)... (P<0.05)".
    reference_factors = disease_modifiers([DiseaseState.ACROMEGALY])
    store.add(
        SeedParameterBuilder()
        .for_subject(acromegaly_id)
        .kind(ParameterKind.CUSTOM, custom_name="Gb_factor")
        .value(reference_factors["Gb"], unit="dimensionless (multiplier on PathwayParameters.Gb)")
        .from_pubmed("1569148", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .gold()
        .notes("Within-study plasma glucose ratio (5.9/5.2) at active "
               "(pre-adenomectomy) acromegaly -- computed live against "
               "the real disease_modifiers() function, not a duplicated "
               "literal. Gb-ONLY -- the source confirms qualitatively "
               "that active acromegaly causes substantial insulin "
               "resistance during a clamp ('substantially resistant to "
               "the actions of insulin on both glucose and lipid "
               "metabolism'), but reports no extractable glucose-"
               "infusion-rate/M-value number, so p3/ke/secretion_scale "
               "are deliberately NOT registered here -- an honest "
               "partial implementation, not an assumed-null claim. "
               "Represents ACTIVE, untreated acromegaly only -- post-"
               "adenomectomy normalization is reported only "
               "qualitatively ('all of these metabolic abnormalities "
               "were abolished'), no extractable post-op number, not "
               "used.")
        .build()
    )


# ---------------------------------------------------------------------------
# OSA p3 constant (MASTER_PLAN 5.1, chunk 33; Track 1 row 1.39) --
# registers physiology/modifiers.py's own OSA disease-state p3 factor
# against a new "osa" CONTEXT entity (context_seed.py). Same category
# as chronic_kidney_disease (chunk 11), pcos (chunk 31), and
# acromegaly (chunk 32) -- a disease state, not a molecule.
# ---------------------------------------------------------------------------

def _seed_osa(store: ParameterStore, context_registry: ContextRegistry) -> None:
    """
    Registers the OSA p3 factor. GOLD -- a direct ratio computed from
    a single primary source's own two reported clamp values (5.75/6.79),
    no unit conversion -- cross-checked live against the real
    disease_modifiers() function, not a duplicated literal.
    """
    osa_id = _resolve_context(
        context_registry, "osa", ContextEntityType.DISEASE_STATE,
    )

    from biosim_physiology.modifiers import disease_modifiers
    from biosim_physiology.profile import DiseaseState

    # Harsch IA, Schahin SP, Radespiel-Troger M, et al. "Continuous
    # positive airway pressure treatment rapidly improves insulin
    # sensitivity in patients with obstructive sleep apnea syndrome."
    # Am J Respir Crit Care Med. 2004;169(2):156-62. PMID:14512265 --
    # verified live via NCBI efetch this chunk. n=40 (AHI>20),
    # hyperinsulinemic euglycemic clamp before and after CPAP. Directly
    # quoted: "Insulin sensitivity significantly increased after 2 days
    # (5.75 +/- 4.20 baseline versus 6.79 +/- 4.91 micromol/kg.min; p =
    # 0.003)."
    reference_factors = disease_modifiers([DiseaseState.OSA])
    store.add(
        SeedParameterBuilder()
        .for_subject(osa_id)
        .kind(ParameterKind.CUSTOM, custom_name="p3_factor")
        .value(reference_factors["p3"], unit="dimensionless (multiplier on PathwayParameters.p3)")
        .from_pubmed("14512265", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .gold()
        .notes("Untreated (pre-CPAP baseline) clamp-derived insulin-"
               "sensitivity ratio (5.75/6.79) -- computed live against "
               "the real disease_modifiers() function, not a duplicated "
               "literal. A GENUINE, UNRESOLVED literature conflict "
               "exists: Smurra et al. 2001 (PMID:11311683) found no "
               "change in glucose-insulin metabolism after CPAP in a "
               "smaller (n=6) cohort, directly contradicting Harsch --"
               "disclosed, not resolved by fiat, same discipline as "
               "PCOS's own unresolved conflict. Both studies name "
               "OBESITY as a likely confound; this platform composes "
               "OSA's factor with bmi_modifiers() regardless (no "
               "obesity-gating mechanism exists for any disease state, "
               "a real, disclosed limitation). Mutually exclusive with "
               "intermittent_hypoxia at the profile.py level -- both "
               "lower p3 via the same mechanism and would double-count "
               "if composed. ke/Gb/secretion_scale confirmed unreported "
               "by either source, not registered.")
        .build()
    )


# ---------------------------------------------------------------------------
# Chronic kidney disease p3/ke constants (MASTER_PLAN 5.1, chunk 11;
# Track 1 row 1.24) -- registers physiology/modifiers.py's own CKD
# disease-state factors against the newly-registered
# "chronic_kidney_disease" CONTEXT entity (context_seed.py). First
# chunk in this campaign to resolve against ContextRegistry, not
# IDRegistry -- CKD is a disease state, not a molecule.
# ---------------------------------------------------------------------------

def _seed_ckd_disease_state(store: ParameterStore, context_registry: ContextRegistry) -> None:
    """
    Registers the two CKD-specific factors (p3, ke) as a pair of GOLD
    entries under the SAME subject -- both are directly-reported
    covariate-adjusted values from one primary source, not a back-
    solved estimate, so GOLD (not BRONZE like the drug-dose-response
    constants above) is the correct tier here.
    """
    ckd_id = _resolve_context(
        context_registry, "chronic_kidney_disease", ContextEntityType.DISEASE_STATE,
    )

    # de Boer IH, Zelnick L, Afkarian M, Ayers E, Curtin L, Himmelfarb J,
    # Ikizler TA, Kahn SE, Kestenbaum B, Utzschneider K. "Impaired
    # Glucose and Insulin Homeostasis in Moderate-Severe CKD." J Am Soc
    # Nephrol. 2016;27(9):2861-71. PMID:26823551 -- verified live via
    # NCBI efetch this session (see physiology/modifiers.py's own
    # citation comment for the full derivation). n=59 nondiabetic CKD
    # vs. n=39 healthy controls: covariate-adjusted difference in
    # insulin sensitivity (-0.7 mg/min per uU/ml on a control mean of
    # 5.0 -> factor (5.0-0.7)/5.0=0.86) and insulin clearance (-85
    # ml/min on a control mean of 998 -> factor (998-85)/998=0.915).
    store.add(
        SeedParameterBuilder()
        .for_subject(ckd_id)
        .kind(ParameterKind.CUSTOM, custom_name="p3_factor")
        .value(0.86, unit="dimensionless (multiplier on PathwayParameters.p3)")
        .from_pubmed("26823551", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .gold()
        .notes("Covariate-adjusted insulin-sensitivity factor -- the "
               "paper's own adjustment for physical activity, diet, fat "
               "mass, and fat-free mass, not the crude (unadjusted) "
               "ratio, avoiding double-counting against this platform's "
               "own separate BMI/fitness modifiers. Disclosed "
               "uncertainty: the adjusted difference's own 95% CI upper "
               "bound touches exactly 0.0 (not fully statistically "
               "robust in this n=59 cohort) -- used anyway per this "
               "project's own established discipline, flagged honestly "
               "rather than presented as equally solid as a tighter CI.")
        .build()
    )
    store.add(
        SeedParameterBuilder()
        .for_subject(ckd_id)
        .kind(ParameterKind.CUSTOM, custom_name="ke_factor")
        .value(0.915, unit="dimensionless (multiplier on PathwayParameters.ke)")
        .from_pubmed("26823551", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .gold()
        .notes("Covariate-adjusted insulin-clearance factor, same "
               "source and adjustment methodology as p3_factor above. "
               "The first disease-state factor in this platform to "
               "touch ke (insulin clearance) rather than only "
               "p1/p3/secretion_scale/Gb.")
        .build()
    )


def _seed_hypothyroidism(store: ParameterStore, context_registry: ContextRegistry) -> None:
    """
    Registers hypothyroidism's p1/p3 factors -- BRONZE, citing the two
    directionally-supportive studies, with the genuine Shah 1975
    contradiction disclosed in each note rather than silently dropped.
    """
    hypo_id = _resolve_context(
        context_registry, "hypothyroidism", ContextEntityType.DISEASE_STATE,
    )

    from biosim_physiology.modifiers import disease_modifiers
    from biosim_physiology.profile import DiseaseState

    factors = disease_modifiers([DiseaseState.HYPOTHYROIDISM])

    conflict_note = (
        "Directionally supported by Rochon et al. 2003 (Clin Sci, "
        "PMID:12519082, n=6, basal glucose disposal increased after "
        "restoring euthyroidism) and Handisurya et al. 2008 (Clin "
        "Endocrinol, PMID:18429948, n=12, insulin sensitivity improved "
        "after T4 replacement) -- both live-verified this session, "
        "neither reports a specific magnitude. A GENUINE, UNRESOLVED "
        "conflict exists: Shah et al. 1975 (Diabetes, PMID:1175861, "
        "n=10, direct insulin-tolerance-test measurement) found NO "
        "insulin resistance in hypothyroidism at all -- glucose nadir "
        "36+/-4 mg/dl (hypothyroid) vs. 43+/-3 (normal), 'might in fact "
        "suggest increased sensitivity.' Not resolved here, per this "
        "project's standing rule against picking a side by fiat -- BRONZE "
        "reflects both the lack of a reported magnitude AND this "
        "unresolved conflict."
    )
    store.add(
        SeedParameterBuilder()
        .for_subject(hypo_id)
        .kind(ParameterKind.CUSTOM, custom_name="p1_factor")
        .value(factors["p1"], unit="dimensionless (multiplier on PathwayParameters.p1)")
        .from_pubmed("12519082", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .bronze()
        .notes(conflict_note)
        .build()
    )
    store.add(
        SeedParameterBuilder()
        .for_subject(hypo_id)
        .kind(ParameterKind.CUSTOM, custom_name="p3_factor")
        .value(factors["p3"], unit="dimensionless (multiplier on PathwayParameters.p3)")
        .from_pubmed("18429948", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .bronze()
        .notes(conflict_note)
        .build()
    )
