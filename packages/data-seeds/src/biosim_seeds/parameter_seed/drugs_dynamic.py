# SAVE TO: packages\data-seeds\src\biosim_seeds\parameter_seed\drugs_dynamic.py
"""Drugs with ODE-coupled / pharmacokinetic parameters or meal-absorption
effects: exogenous insulin, SGLT2 inhibitors, cortisol, GLP-1 class,
prednisone, alcohol, alpha-glucosidase inhibitors."""

from __future__ import annotations

from biosim_id_resolver import IDRegistry
from biosim_parameter_db import (
    ParameterKind,
    ParameterStore,
)

from ._helpers import SeedParameterBuilder, _resolve

# ---------------------------------------------------------------------------
# Exogenous insulin PK parameters (MASTER_PLAN 5.1 -- first entry in the
# systematic parameter-breadth curation campaign: this platform's own
# extension-pathway constants, e.g. exogenous_insulin_pathway.py's
# INSULIN_REGULAR_KA_PER_MIN, have real, live-verified literature citations
# in their own module docstrings but were never formally registered in
# ParameterStore -- every extension pathway (SGLT2i, GLP-1, cortisol,
# sulfonylurea, acarbose/miglitol, exercise, NEFA) has this same gap; this
# is the first of many small follow-on chunks, not the whole campaign)
# ---------------------------------------------------------------------------

def _seed_exogenous_insulin_pk(store: ParameterStore, registry: IDRegistry) -> None:
    """
    Registers `biosim_ode.exogenous_insulin_pathway`'s own 4 calibrated
    two-sequential-depot PK constants. Deliberately BRONZE, not GOLD --
    the underlying Howey et al. 1994 clinical measurement (Tmax/Cmax) IS
    a real, directly-reported primary-source number, but these specific
    constants are a disclosed engineering BACK-SOLVE against that
    measurement (see exogenous_insulin_pathway.py's own docstrings for
    the exact derivation), not themselves directly-reported values --
    matching this file's own corrected BRONZE definition ("an indirect
    estimate... or a disclosed engineering back-solve/calibration
    against a real measurement", see the module docstring).
    """
    insulin_id = _resolve(registry, "insulin")

    # Howey DC, Bowsher RR, Brunelle RL, Woodworth JR. "[Lys(B28),
    # Pro(B29)]-human insulin. A rapidly absorbed analogue of human
    # insulin." Diabetes. 1994;43(3):396-402. PMID:8314011 -- verified
    # live via NCBI efetch (see exogenous_insulin_pathway.py's own
    # module docstring). n=10 healthy men, 10-U SC dose, euglycaemic
    # clamp. Arm C = regular human insulin: Tmax=101min, Cmax=308pM
    # excess above baseline.
    store.add(
        SeedParameterBuilder()
        .for_subject(insulin_id)
        .kind(ParameterKind.CUSTOM, custom_name="exogenous_insulin_regular_ka_per_min")
        .value(0.0113, unit="min-1")
        .from_pubmed("8314011", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .bronze()
        .notes("Two-sequential-equal-rate-constant SC depot absorption rate "
               "for regular human insulin. Back-solved (not directly "
               "reported) so that, combined with this platform's own "
               "independently-anchored ke=0.093/min (Sherwin 1974), the "
               "resulting Tmax matches Howey 1994's observed 101min for a "
               "10U dose. See exogenous_insulin_pathway.py's own "
               "INSULIN_REGULAR_KA_PER_MIN docstring for the full "
               "derivation, including why Kobayashi et al. 1983's "
               "(PMID:6339306) own one-compartment ka/ke pair was NOT "
               "reused (a flip-flop-kinetics apparent rate, not this "
               "model's true clearance).")
        .build()
    )
    store.add(
        SeedParameterBuilder()
        .for_subject(insulin_id)
        .kind(ParameterKind.CUSTOM, custom_name="exogenous_insulin_regular_pm_per_unit")
        .value(695.1, unit="pM/Unit")
        .from_pubmed("8314011", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .bronze()
        .notes("Per-Unit plasma insulin scaling factor for regular human "
               "insulin. Back-solved (not derived from an independently-"
               "cited insulin distribution volume, which no source "
               "anchors) so that a 10-Unit dose reproduces Howey 1994's "
               "observed peak excess (308pM) exactly, at the ka above. "
               "Same 'isolated, disclosed engineering calibration' "
               "category as this platform's own HPA_FEEDBACK_GAIN and "
               "SGLT2i beta constants -- a real, if imprecise, scaling "
               "factor, not a physical volume of distribution.")
        .build()
    )
    store.add(
        SeedParameterBuilder()
        .for_subject(insulin_id)
        .kind(ParameterKind.CUSTOM, custom_name="exogenous_insulin_lispro_ka_per_min")
        .value(0.0353, unit="min-1")
        .from_pubmed("8314011", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .bronze()
        .notes("Two-sequential-equal-rate-constant SC depot absorption rate "
               "for rapid-acting insulin lispro ([Lys(B28),Pro(B29)]-human "
               "insulin, brand name Humalog). Back-solved so that a "
               "10-Unit dose reproduces Howey 1994's own rapid-acting Tmax "
               "(42min) through the full composed ODE. Same within-study, "
               "same-protocol, same-10U-dose comparison as the regular-"
               "insulin entries above (Howey 1994's own arms A vs. C) -- "
               "no cross-study confound.")
        .build()
    )
    store.add(
        SeedParameterBuilder()
        .for_subject(insulin_id)
        .kind(ParameterKind.CUSTOM, custom_name="exogenous_insulin_lispro_pm_per_unit")
        .value(543.1, unit="pM/Unit")
        .from_pubmed("8314011", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .bronze()
        .notes("Per-Unit plasma insulin scaling factor for insulin lispro. "
               "Back-solved so that a 10-Unit dose reproduces Howey 1994's "
               "own observed rapid-acting peak excess (698pM) exactly, at "
               "the ka above.")
        .build()
    )


# ---------------------------------------------------------------------------
# Renal SGLT2i parameters (MASTER_PLAN 5.1, chunk 2 -- unlike the exogenous-
# insulin entries above, these ARE directly-quoted primary clinical
# measurements, not back-solved calibrations; see each entry's own notes
# for the exact quoted text)
# ---------------------------------------------------------------------------

def _seed_renal_sglt2i(store: ParameterStore, registry: IDRegistry) -> None:
    """
    Registers `biosim_ode.renal_pathway`'s own 5 directly-cited GFR/TmG
    constants against the newly-registered "SGLT2" subject entity
    (registry_seed.py). GOLD -- real, precisely-quoted, DIRECTLY
    REPORTED clinical hyperglycaemic-clamp measurements (the paper's own
    headline finding for exactly this quantity), matching this file's
    own established GOLD standard (Bergman/Sherwin's
    foundational minimal-model kinetics) rather than a review/secondary
    source or an indirect estimate.

    CORRECTED (housekeeping-backlog fix -- see this file's own module
    docstring): originally tiered SILVER, reasoning from this file's
    then-current (and, it turned out, wrong) docstring bar of
    "replicated in >=10 independent studies" -- a bar even the existing
    GOLD entries never actually met. Re-tiered to GOLD after checking
    this file's REAL established practice (GOLD = a primary source's
    own directly-reported value, not a replication count) directly
    against every other tier assignment in this file.

    Confirmed this re-tier has zero simulation-behaviour impact --
    checked directly, not assumed. CORRECTION to an earlier draft of
    this note (found by mandatory code review): `ParameterStore
    ._index_key()` deliberately does NOT include `custom_kind_name`
    (see chunk 1's own conflict-detection fix above), so all 5 of
    these CUSTOM entries genuinely DO share one `query()`/`best()`
    index bucket under subject="SGLT2" -- they are NOT independent,
    non-competing keys the way a first read of this note claimed. The
    real reason this re-tier is safe: every entry sharing that bucket
    was re-tiered SILVER -> GOLD together, as one block, so their
    RELATIVE ranking (and therefore `best()`'s own tie-broken pick,
    verified directly to be unchanged before/after) is preserved
    exactly. Re-tiering only ONE member of this group in a future
    session WOULD be capable of changing `best()`'s selection -- do
    not assume these 5 are safely independent.
    """
    sglt2_id = _resolve(registry, "SGLT2")

    # DeFronzo RA, Hompesch M, Kasichayanula S, et al. "Characterization
    # of renal glucose reabsorption in response to dapagliflozin in
    # healthy subjects and subjects with type 2 diabetes." Diabetes
    # Care. 2013;36(10):3169-3176. PMID:23735727, PMCID:PMC3781504 --
    # verified live via full text this session (see renal_pathway.py's
    # own module docstring). Stepped hyperglycaemic clamp with an
    # explicit TmG/splay/threshold model.
    store.add(
        SeedParameterBuilder()
        .for_subject(sglt2_id)
        .kind(ParameterKind.CUSTOM, custom_name="gfr_reference_ml_per_min")
        .value(127.5, unit="mL/min")
        .from_pubmed("23735727", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .gold()
        .notes("Glomerular filtration rate, reference value. GFR measured "
               "\"similar\" between T2D (127+/-34) and healthy (128+/-36) "
               "subjects in this study -- a single shared reference value "
               "(the mean of the two group means, 127.5) is used rather "
               "than two near-identical presets, per renal_pathway.py's "
               "own GFR_REFERENCE_ML_PER_MIN docstring.")
        .build()
    )
    store.add(
        SeedParameterBuilder()
        .for_subject(sglt2_id)
        .kind(ParameterKind.CUSTOM, custom_name="dapagliflozin_tmg_healthy_mg_per_min")
        .value(317.0, unit="mg/min")
        .from_pubmed("23735727", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .gold()
        .notes("Baseline maximal renal glucose transport (TmG), healthy "
               "subjects. Directly quoted: \"420 vs. 317 mg/min, excluding "
               "outlier; P < 0.001\" (T2D vs. healthy).")
        .build()
    )
    store.add(
        SeedParameterBuilder()
        .for_subject(sglt2_id)
        .kind(ParameterKind.CUSTOM, custom_name="dapagliflozin_tmg_t2d_mg_per_min")
        .value(420.0, unit="mg/min")
        .from_pubmed("23735727", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .gold()
        .notes("Baseline maximal renal glucose transport (TmG), type-2-"
               "diabetic subjects -- 32% higher than healthy (see the "
               "healthy entry's own quoted text). Real, population-"
               "specific reabsorptive-capacity elevation in T2D, this "
               "paper's own words: \"elevated TmG, splay, and threshold "
               "compared with controls.\"")
        .build()
    )

    # Al-Jobori H, Daniele G, Adams J, et al. "Empagliflozin and kinetics
    # of renal glucose transport in healthy individuals and individuals
    # with type 2 diabetes." Diabetes. 2017;66(7):1999-2006.
    # PMID:28428225 -- verified live this session (abstract only; no
    # open-access full text located).
    store.add(
        SeedParameterBuilder()
        .for_subject(sglt2_id)
        .kind(ParameterKind.CUSTOM, custom_name="empagliflozin_tmg_healthy_mg_per_min")
        .value(337.0, unit="mg/min")
        .from_pubmed("28428225", figure_or_table="abstract")
        .in_human()
        .clinical()
        .gold()
        .notes("Baseline maximal renal glucose transport (TmG), healthy "
               "subjects. Directly quoted: \"459 +/- 53 vs. 337 +/- 25 "
               "mg/min; P < 0.05\" (T2D vs. healthy).")
        .build()
    )
    store.add(
        SeedParameterBuilder()
        .for_subject(sglt2_id)
        .kind(ParameterKind.CUSTOM, custom_name="empagliflozin_tmg_t2d_mg_per_min")
        .value(459.0, unit="mg/min")
        .from_pubmed("28428225", figure_or_table="abstract")
        .in_human()
        .clinical()
        .gold()
        .notes("Baseline maximal renal glucose transport (TmG), type-2-"
               "diabetic subjects -- this paper's own words: \"T2D "
               "subjects had significantly greater TmG at baseline.\"")
        .build()
    )


# ---------------------------------------------------------------------------
# Cortisol PK/potency constants (MASTER_PLAN 5.1, chunk 6) -- registers
# cortisol_pathway.py's own two exogenous-steroid-dosing constants
# (PREDNISONE_TO_CORTISOL_MG_RATIO, DAILY_CORTISOL_PRODUCTION_MG) against
# the newly-registered "cortisol" subject entity (registry_seed.py).
# Both PMIDs below were independently re-verified live BEFORE this chunk
# was written -- the two PMIDs originally sitting in cortisol_pathway.py's
# own comment (412539, 1955512) turned out to be mistranscribed (see
# CLAUDE.md "DONE — Fixed two mistranscribed citations in
# cortisol_pathway.py"); the values themselves were unaffected, only the
# citation numbers were wrong.
# ---------------------------------------------------------------------------

def _seed_cortisol_pk(store: ParameterStore, registry: IDRegistry) -> None:
    """
    Registers cortisol_pathway.py's PREDNISONE_TO_CORTISOL_MG_RATIO and
    DAILY_CORTISOL_PRODUCTION_MG.

    BRONZE for the potency ratio -- NOT because it's a single study
    (this file's corrected tier definition no longer uses replication
    count at all, see the module docstring), but because no abstract
    text was ever available to independently verify the exact ratio
    against the paper's own words -- cited on bibliographic (title/
    author/journal/year) match alone, a genuinely less-certain
    situation than a directly-quotable statistic, even though the
    underlying fact (a real, correctly-identified paper) isn't itself
    "inferred or estimated" in the strict sense.

    GOLD for the production rate (re-tiered from an earlier SILVER
    assignment, housekeeping-backlog fix): a directly-quoted, primary-
    source statistic (Esteban 1991's own measured FPR), matching this
    file's own corrected GOLD standard exactly -- same shape of
    evidence as the renal SGLT2i entries above, now also GOLD.

    Note on `query()`/`best()` sharing (correction, found by mandatory
    code review -- see the renal SGLT2i docstring above for the full
    explanation): these two constants DO share one CUSTOM-kind index
    bucket under subject="cortisol" (`ParameterStore._index_key()`
    excludes `custom_kind_name` by design), so they are NOT
    independent keys. `best()` picks `daily_cortisol_production_mg`
    (GOLD, verified directly) over `prednisone_to_cortisol_mg_ratio`
    (BRONZE) both before and after this re-tier -- only one of the two
    was ever re-tiered here (the other, BRONZE, was untouched), so
    their relative ranking couldn't have flipped either way.
    """
    cortisol_id = _resolve(registry, "cortisol")

    # Meikle AW, Tyler FH. "Potency and duration of action of
    # glucocorticoids. Effects of hydrocortisone, prednisone and
    # dexamethasone on human pituitary-adrenal function." Am J Med.
    # 1977;63(2):200-207. PMID:888843 -- verified live via NCBI efetch
    # this session (title/authors/journal/year confirmed; abstract
    # unavailable for this pre-1990s record, a common gap already seen
    # elsewhere in this codebase). Standard clinical glucocorticoid
    # dose-equivalency: 5mg prednisone ~ 20mg hydrocortisone (cortisol),
    # i.e. 1mg prednisone ~ 4mg cortisol-equivalent.
    store.add(
        SeedParameterBuilder()
        .for_subject(cortisol_id)
        .kind(ParameterKind.CUSTOM, custom_name="prednisone_to_cortisol_mg_ratio")
        .value(4.0, unit="mg cortisol-equivalent per mg prednisone")
        .from_pubmed("888843", figure_or_table="text")
        .in_human()
        .clinical()
        .bronze()
        .notes("Standard clinical glucocorticoid potency-equivalence "
               "ratio (5mg prednisone ~ 20mg hydrocortisone). Real paper, "
               "correctly identified by title/author/journal/year match "
               "-- no verbatim abstract quote available (efetch returns "
               "an empty abstract for this pre-1990s record), so cited "
               "on bibliographic match alone, same practice already used "
               "elsewhere in this codebase for other older papers with "
               "no extractable abstract text.")
        .build()
    )

    # Esteban NV, Loughlin T, Yergey AL, Zawadzki JK, Booth JD,
    # Winterer JC, Loriaux DL. "Daily cortisol production rate in man
    # determined by stable isotope dilution/mass spectrometry." J Clin
    # Endocrinol Metab. 1991;72(1):39-45. PMID:1986026 -- verified live
    # via NCBI efetch this session. n=12 normal volunteers, stable
    # isotope ([9,12,12-2H3]cortisol) tracer dilution.
    store.add(
        SeedParameterBuilder()
        .for_subject(cortisol_id)
        .kind(ParameterKind.CUSTOM, custom_name="daily_cortisol_production_mg")
        .value(10.0, unit="mg/day hydrocortisone-equivalent")
        .from_pubmed("1986026", figure_or_table="Results text")
        .in_human()
        .clinical()
        .gold()
        .notes("Physiological cortisol production rate, unstressed adult. "
               "Directly quoted: \"Results from 12 normal volunteers "
               "revealed a FPR of 27.3 +/- 7.5 mumol/day (9.9 +/- 2.7 "
               "mg/day)\" -- this platform's own DAILY_CORTISOL_"
               "PRODUCTION_MG=10.0 matches the reported 9.9 mg/day "
               "almost exactly.")
        .build()
    )


# ---------------------------------------------------------------------------
# Alcohol acute hepatic EGP suppression constant (MASTER_PLAN 5.1,
# chunk 36) -- registers alcohol_hepatic_pathway.py's own
# ALCOHOL_EGP_REDUCTION_FRACTION_AT_TESTED_DOSE against the newly-
# registered "alcohol" metabolite subject entity (registry_seed.py).
# ---------------------------------------------------------------------------

def _seed_alcohol_egp(store: ParameterStore, registry: IDRegistry) -> None:
    """
    Registers the alcohol acute hepatic-glucose-production-suppression
    fraction at its own 48g reference dose. GOLD -- a directly-reported
    within-study percentage reduction, no unit conversion -- cross-
    checked live against the real
    ALCOHOL_EGP_REDUCTION_FRACTION_AT_TESTED_DOSE constant, not a
    duplicated literal.
    """
    alcohol_id = _resolve(registry, "alcohol")

    from biosim_ode.alcohol_hepatic_pathway import (
        ALCOHOL_EGP_REDUCTION_FRACTION_AT_TESTED_DOSE,
        ALCOHOL_TESTED_DOSE_G,
    )

    # Siler SQ, Neese RA, Christiansen MP, Hellerstein MK. "The
    # inhibition of gluconeogenesis following alcohol in humans." Am J
    # Physiol. 1998;275(5):E897-907. PMID:9815011 -- verified live via
    # NCBI efetch this chunk. Normal men, overnight-fasted, 48g oral
    # alcohol vs. placebo, stable-isotope tracers. Directly quoted:
    # "Hepatic glucose output decreased 12% after alcohol consumption,
    # from 2.03 +/- 0.21 to 1.79 +/- 0.21 mg.kg-1.min-1 (P < 0.05 vs.
    # baseline), but did not change following the placebo."
    store.add(
        SeedParameterBuilder()
        .for_subject(alcohol_id)
        .kind(ParameterKind.CUSTOM, custom_name=f"egp_reduction_fraction_at_{ALCOHOL_TESTED_DOSE_G:.0f}g")
        .value(
            ALCOHOL_EGP_REDUCTION_FRACTION_AT_TESTED_DOSE,
            unit="dimensionless fraction (reduction in hepatic glucose output)",
        )
        .from_pubmed("9815011", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .gold()
        .notes("Directly-reported hepatic glucose output reduction at "
               "the study's own 48g oral alcohol dose -- computed live "
               "against the real "
               "ALCOHOL_EGP_REDUCTION_FRACTION_AT_TESTED_DOSE constant, "
               "not a duplicated literal. Route-matched anchor (oral, "
               "same as this platform's own DrugDose/ScenarioDefinition "
               "alcohol events) -- a second, IV-route source (Shelmet et "
               "al. 1988, PMID:3280601) corroborates a real hepatic/"
               "peripheral effect but was NOT used as the anchor here, "
               "since route (oral vs. IV) plausibly changes the "
               "magnitude via first-pass hepatic alcohol metabolism. "
               "Both sources explicitly confirm NO effect on plasma "
               "glucose CONCENTRATION itself (only on the underlying "
               "EGP flux) -- Gb is deliberately NOT touched by this "
               "mechanism, matching the pathway's own ODE-coupled "
               "hepatic_egp_delta design (not a static Gb modifier).")
        .build()
    )


# ---------------------------------------------------------------------------
# GLP-1 basal/meal-response PK constants (MASTER_PLAN 5.1, chunk 23) --
# registers glp1_pathway.py's own GLP1_BASAL_PM/GLP1_MEAL_RESPONSE
# against the newly-registered "GLP-1" subject entity (registry_seed.py).
# ---------------------------------------------------------------------------

def _seed_glp1_pk(store: ParameterStore, registry: IDRegistry) -> None:
    """
    Registers glp1_pathway.py's GLP1_BASAL_PM and GLP1_MEAL_RESPONSE.

    Only these two of the module's eight constants are registered here --
    the rest (GLP1_EC50_SECRETION, GLP1_EC50_INSULIN, GLP1_MAX_INCRETIN,
    GLP1_EMPTYING_IC50, GLP1_EMPTYING_MAX) are the module's own disclosed
    uncited engineering assumptions (see glp1_pathway.py's own inline
    comments), not literature-anchored values -- matching this campaign's
    established practice of curating only genuinely anchored constants,
    not every free parameter in a module.

    GOLD for GLP1_BASAL_PM -- a directly-quoted, primary-source fasting
    concentration (Orskov 1994's own measured value).

    BRONZE for GLP1_MEAL_RESPONSE -- same category as
    NUTRIENT_SECRETION_PM_PER_G/the acarbose-miglitol kabs scales: a
    real, disclosed engineering back-solve (numerically re-derived via
    bisection against the composed ODE to reproduce the SAME paper's own
    41 pM postprandial peak), not a directly-reported number itself.

    GLP1_HALF_LIFE_MIN (2.0 min) is deliberately NOT registered here --
    checked live via PubMedFetcher against all three of this module's
    other cited PMIDs (Holst 2007 PMID:17928588, Drucker & Nauck 2006
    PMID:17098089, Nauck et al. 1993 PMID:8405741) plus a targeted
    follow-up search (code-reviewer finding: the first pass only re-
    checked pre-existing citations, not a fresh search for a dedicated
    GLP-1 PK paper). **Vilsboll T, Agerso H, Krarup T, Holst JJ.
    "Similar elimination rates of glucagon-like peptide-1 in obese type
    2 diabetic patients and healthy subjects." J Clin Endocrinol Metab.
    2003;88(1):220-4. PMID:12519856** (live-verified) is a real, directly
    on-topic IV-bolus PK study -- but its own abstract reports clearance
    (4.0-9.0 L/min, dose-dependent, nonlinear/saturable) and volume of
    distribution (a single 9-26 L RANGE, not paired per-dose) separately,
    with no directly-stated half-life number and no clean way to compute
    one without picking an arbitrary CL/Vd pairing the abstract itself
    doesn't specify -- full text was not fetched to resolve this. Holst's
    own abstract states GLP-1 "is extremely rapidly metabolized" but,
    like the other three, gives no specific numeric half-life either --
    a genuine, now more thoroughly searched, gap, not force-cited to an
    approximately-matching source or a fabricated CL/Vd pairing.
    """
    glp1_id = _resolve(registry, "GLP-1")

    # Orskov C, Rabenhoj L, Wettergren A, Kofod H, Holst JJ. "Tissue and
    # plasma concentrations of amidated and glycine-extended glucagon-
    # like peptide I in humans." Diabetes. 1994;43(4):535-9. PMID:8138058
    # -- verified live via NCBI efetch this session (twice, matching this
    # module's own established re-verification discipline for load-
    # bearing citations). n=6. Directly quoted: "The mean concentrations
    # of amidated GLP-I ... in fasting plasma were 7 +/- 1 ... pM (n=6)."
    store.add(
        SeedParameterBuilder()
        .for_subject(glp1_id)
        .kind(ParameterKind.CUSTOM, custom_name="glp1_basal_pm")
        .value(7.0, unit="pM")
        .from_pubmed("8138058", figure_or_table="Results text")
        .in_human()
        .clinical()
        .gold()
        .notes("Fasting plasma GLP-1 concentration. Directly quoted: "
               "\"The mean concentrations of amidated GLP-I ... in "
               "fasting plasma were 7 +/- 1 ... pM (n=6).\" This "
               "platform's own GLP1_BASAL_PM=7.0 matches exactly.")
        .build()
    )

    # Same source's own postprandial peak, used as the target for a
    # numerical back-solve (not a directly-reported constant itself):
    # "In response to a breakfast meal, the concentration of amidated
    # GLP-I rose significantly amounting to 41 +/- 5 pM 90 min after the
    # meal ingestion." GLP1_MEAL_RESPONSE=4.0 (pM per mM glucose above
    # threshold) was numerically re-derived via bisection against the
    # real composed GlucoseInsulinPathway+GLP1Pathway ODE (standard 75g
    # reference meal) to reproduce this 41 pM peak MAGNITUDE -- NOT the
    # paper's own 90-minute timing, which this platform's own reference
    # meal (glucose-and-therefore-GLP1 peak ~t=25min) cannot jointly
    # match without a second, unavailable glucose-kinetics anchor from
    # the same source (same "magnitude-only, timing undetermined" gap
    # already disclosed in glp1_pathway.py's own module docstring).
    store.add(
        SeedParameterBuilder()
        .for_subject(glp1_id)
        .kind(ParameterKind.CUSTOM, custom_name="glp1_meal_response_pm_per_mm")
        .value(4.0, unit="pM per mM glucose above threshold")
        .from_pubmed("8138058", figure_or_table="Results text")
        .in_human()
        .clinical()
        .bronze()
        .notes("Disclosed engineering back-solve, not a directly-"
               "reported value: numerically re-derived (bisection "
               "against the real composed ODE) to reproduce the same "
               "source's own directly-quoted 41 pM postprandial peak "
               "magnitude, not its 90-minute timing.")
        .build()
    )


# ---------------------------------------------------------------------------
# Alpha-glucosidase inhibitor kabs-scale constants (MASTER_PLAN 5.1, chunk 8)
# -- registers drug_interactions.py's ACARBOSE_100MG_KABS_SCALE and
# MIGLITOL_50MG_KABS_SCALE against two newly-registered "acarbose"/
# "miglitol" metabolite subjects (registry_seed.py).
# ---------------------------------------------------------------------------

def _seed_alpha_glucosidase_inhibitor_kabs(store: ParameterStore, registry: IDRegistry) -> None:
    """
    Registers the two dose-referenced meal-absorption-rate (kabs)
    reduction factors behind this platform's acarbose/miglitol
    mechanism. Both BRONZE -- each anchored to a real, precisely-quoted
    clinical measurement, but the STORED value itself is a disclosed
    engineering back-solve against that measurement (a kabs multiplier
    that reproduces the paper's own reported peak-height/peak-time
    numbers through this platform's own ODE model), not a directly-
    reported number -- same category as the exogenous-insulin PK
    constants.
    """
    acarbose_id = _resolve(registry, "acarbose")
    miglitol_id = _resolve(registry, "miglitol")

    # Chen Z, Fu X, Kuang J, Chen J, Chen H, Pei J, Yang H. "Single-dose
    # acarbose decreased glucose-dependent insulinotropic peptide and
    # glucagon levels in Chinese patients with newly diagnosed type 2
    # diabetes mellitus after a mixed meal." BMC Endocr Disord. 2016.
    # PMID:27686734 -- verified live via NCBI efetch this session (see
    # drug_interactions.py's own module-level citation comment for the
    # full derivation). n=12 newly diagnosed T2D, single 100mg acarbose
    # dose with a mixed meal. Table 2: peak glucose 13.22 mmol/L at
    # 60min (no drug) -> 9.45 mmol/L at 90min (with drug) -- a dual
    # peak-height (-28.5%) AND peak-time (1.5x delay) signature.
    store.add(
        SeedParameterBuilder()
        .for_subject(acarbose_id)
        .kind(ParameterKind.CUSTOM, custom_name="kabs_scale_at_100mg_reference_dose")
        .value(0.35, unit="dimensionless (multiplier on MealInput.kabs)")
        .from_pubmed("27686734", figure_or_table="Table 2")
        .in_human()
        .clinical()
        .bronze()
        .notes("Meal-absorption-rate-constant (kabs) reduction at the "
               "100mg reference dose. Back-solved (not directly "
               "reported) so that, run through this platform's own ODE "
               "model for a 75g meal, it simultaneously reproduces both "
               "of Chen et al. 2016's own independent observations: "
               "-27.0% peak glucose height (vs. the real -28.5%) and a "
               "1.51x peak-time delay (vs. the real 1.5x) -- a genuine "
               "two-number consistency check, not a single-target curve "
               "fit. Chosen over scaling Ra(t)'s output directly because "
               "only a true kabs reduction reproduces BOTH the lower AND "
               "later peak; a pure amplitude scale would only lower the "
               "peak without moving it.")
        .build()
    )

    # Samad AH, Ty Willing TS, Alberti KG, Taylor R. "Effects of BAYm
    # 1099, new alpha-glucosidase inhibitor, on acute metabolic
    # responses and metabolic control in NIDDM over 1 mo." Diabetes
    # Care. 1988;11(2):100-9. PMID:3042310 -- verified live via NCBI
    # efetch this session (BAYm 1099 confirmed to be miglitol's own
    # original development code). n=12 diet-treated NIDDM (T2D)
    # patients, 4 weeks miglitol 50mg TID vs. 4 weeks placebo, 25g
    # sucrose test meal: peak blood glucose 11.3+/-1.0 mM (drug) vs.
    # 14.3+/-1.4 mM (placebo), P<0.001 -- a -21.0% peak-height
    # reduction (single-quantity anchor; no peak-time ratio reported,
    # unlike Chen et al. 2016's dual signature for acarbose).
    store.add(
        SeedParameterBuilder()
        .for_subject(miglitol_id)
        .kind(ParameterKind.CUSTOM, custom_name="kabs_scale_at_50mg_reference_dose")
        .value(0.1446, unit="dimensionless (multiplier on MealInput.kabs)")
        .from_pubmed("3042310", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .bronze()
        .notes("Meal-absorption-rate-constant (kabs) reduction at the "
               "50mg TID reference dose. Back-solved so that, run "
               "through this platform's own ODE model for a matching "
               "25g sucrose meal, it reproduces Samad et al. 1988's own "
               "reported peak-height reduction (-21.0%) exactly. Unlike "
               "the acarbose entry above, only ONE quantity (peak "
               "height) was reported for miglitol -- the resulting "
               "~2.05x peak-time delay this value produces is a "
               "mechanical consequence of the shared kabs-reduction "
               "mechanism, not independently verified against a real "
               "timing figure, disclosed as such rather than presented "
               "as a confirmed second data point.")
        .build()
    )


# ---------------------------------------------------------------------------
# Prednisone/glucocorticoid p3/Gb/secretion_scale dose-response
# constants (MASTER_PLAN 5.1, chunk 20) -- registers
# drug_interactions.py's own glucocorticoid_modifiers() constants
# against the newly-registered "prednisone" metabolite subject
# (registry_seed.py).
#
# Important disclosed caveat, carried over from glucocorticoid_
# modifiers()'s own docstring: SimulationRunner does NOT use these
# values for a real glucocorticoid drug event -- CortisolPathway (the
# genuine PK/PD-timing-aware representation) is the sole path used
# there. This function (and therefore this registration) is real,
# citation-anchored, and reachable via any OTHER caller composing
# DrugDose lists directly -- not a dead/unused code path -- but it is
# NOT this platform's primary glucocorticoid mechanism. Disclosed here
# so a future reader doesn't assume registering this constant means
# it drives the platform's main steroid-response behavior.
# ---------------------------------------------------------------------------

def _seed_prednisone_glucocorticoid(store: ParameterStore, registry: IDRegistry) -> None:
    """
    Registers all three glucocorticoid_modifiers() constants at the
    30mg-prednisone-equivalent dose cap (dose_factor==1.0 exactly) --
    BRONZE, since the exact 0.30/0.10/0.10 magnitudes are Rizza et al.
    1982's own dose-response figures applied via a linear dose-scaling
    assumption below the cap, not a single directly-reported ratio for
    every dose.
    """
    prednisone_id = _resolve(registry, "prednisone")

    from biosim_physiology.drug_interactions import (
        Drug,
        DrugDose,
        glucocorticoid_modifiers,
    )

    # Rizza RA, Mandarino LJ, Gerich JE. "Cortisol-induced insulin
    # resistance in man: impaired suppression of glucose production and
    # stimulation of glucose utilization due to a postreceptor detect of
    # insulin action." J Clin Endocrinol Metab. 1982;54(1):131-8.
    # PMID:7033265 -- verified live via NCBI efetch this session (see
    # drug_interactions.py's own glucocorticoid_modifiers() docstring,
    # which cites this same PMID for "reduces p3 ~30%, elevates Gb
    # ~10%"). Values cross-checked live at the 30mg-prednisone-
    # equivalent reference dose (dose.dose_mg=30.0), where
    # dose_factor==min(30/10,3.0)/3.0==1.0 exactly.
    at_max_dose = glucocorticoid_modifiers(DrugDose(drug=Drug.GLUCOCORTICOID, dose_mg=30.0))
    factors = {
        "p3_reduction_at_30mg_reference_dose": 1.0 - at_max_dose["p3"],
        "gb_increase_at_30mg_reference_dose": at_max_dose["Gb"] - 1.0,
        "secretion_scale_increase_at_30mg_reference_dose": at_max_dose["secretion_scale"] - 1.0,
    }
    for custom_name, value in factors.items():
        store.add(
            SeedParameterBuilder()
            .for_subject(prednisone_id)
            .kind(ParameterKind.CUSTOM, custom_name=custom_name)
            .value(value, unit="dimensionless (fractional change applied to the corresponding PathwayParameters field)")
            .from_pubmed("7033265", figure_or_table="text, Results")
            .in_human()
            .clinical()
            .bronze()
            .notes("Value at the 30mg-prednisone-equivalent dose cap "
                   "(dose_factor==1.0 exactly), cross-checked live "
                   "against the real glucocorticoid_modifiers() "
                   "function, not a duplicated literal. BRONZE, not "
                   "GOLD -- glucocorticoid_modifiers()'s own linear "
                   "dose-scaling below this cap is a disclosed "
                   "simplification, not itself derived from a reported "
                   "dose-response curve in the cited source. IMPORTANT: "
                   "SimulationRunner does NOT use this function's "
                   "values for a real glucocorticoid drug event -- "
                   "CortisolPathway is the sole PK/PD-timing-aware "
                   "representation there (see glucocorticoid_"
                   "modifiers()'s own docstring). This registration is "
                   "real and reachable via other DrugDose-list callers, "
                   "not this platform's primary steroid mechanism.")
            .build()
        )
