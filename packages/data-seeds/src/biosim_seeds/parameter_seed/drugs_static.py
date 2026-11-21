# SAVE TO: packages\data-seeds\src\biosim_seeds\parameter_seed\drugs_static.py
"""Drugs modelled as a static factor (p3, Gb or secretion)."""

from __future__ import annotations

from biosim_id_resolver import IDRegistry
from biosim_parameter_db import (
    ParameterKind,
    ParameterStore,
)

from ._helpers import SeedParameterBuilder, _resolve

# ---------------------------------------------------------------------------
# Hydroxychloroquine chronic p3 constant (MASTER_PLAN 5.1, chunk 37,
# following Track 1 row 1.44) -- registers drug_interactions.py's own
# HYDROXYCHLOROQUINE_P3_FACTOR against the newly-registered
# "hydroxychloroquine" metabolite subject entity (registry_seed.py).
# ---------------------------------------------------------------------------

def _seed_hydroxychloroquine_p3(store: ParameterStore, registry: IDRegistry) -> None:
    """
    Registers HYDROXYCHLOROQUINE_P3_FACTOR (1.20). GOLD -- a directly-
    reported within-arm percentage change, no unit conversion, cross-
    checked live against the real module constant, not a duplicated
    literal. The ANCHOR CHOICE (within-arm +20% rather than the placebo-
    referenced 38.3%) is disclosed in the notes, not hidden.
    """
    hcq_id = _resolve(registry, "hydroxychloroquine")

    from biosim_physiology.drug_interactions import HYDROXYCHLOROQUINE_P3_FACTOR

    # Wasko MC, McClure CK, Kelsey SF, Huber K, Orchard T, Toledo FG.
    # "Antidiabetogenic effects of hydroxychloroquine on insulin
    # sensitivity and beta cell function: a randomised trial."
    # Diabetologia. 2015;58(10):2336-43. PMID:26197707 -- verified live
    # (abstract via NCBI efetch; the authors' own Discussion sentence via
    # full text from NCBI PMC, using biosim_literature.verify_citation()).
    # Randomised double-blind parallel RCT, HCQ 400 mg/day vs placebo,
    # n=17/15, 13 weeks, IVGTT + minimal model. Directly quoted: "+20.0% +/-
    # 7.1% vs -18.4% +/- 7.9%, respectively; p < 0.01; difference: 38.3%".
    store.add(
        SeedParameterBuilder()
        .for_subject(hcq_id)
        .kind(ParameterKind.CUSTOM, custom_name="hydroxychloroquine_p3_factor")
        .value(HYDROXYCHLOROQUINE_P3_FACTOR, unit="unitless multiplier on p3")
        .from_pubmed("26197707", figure_or_table="Results text")
        .in_human()
        .clinical()
        .gold()
        .notes("Within-arm minimal-model insulin-sensitivity increase on "
               "HCQ 400 mg/day for 13 weeks (a mean of per-subject "
               "percentages, not a ratio of group means) -- computed live "
               "against the real HYDROXYCHLOROQUINE_P3_FACTOR, not a "
               "duplicated literal. ANCHOR CHOICE: the placebo-referenced "
               "difference (38.3%) is NOT used because the placebo arm's "
               "own -18.4% decline is unexplained -- the authors' Discussion "
               "says so and reads the result by \"the change within the HCQ "
               "group (20 +/- 7.1%)\". Primary analysis was ITT with "
               "baseline carried forward for dropouts (conservative). "
               "Cohort: non-diabetic, overweight/obese with >=1 insulin-"
               "resistance marker -- the authors state effects in T2D "
               "\"must be inferred\". Single tested dose; treated as the "
               "chronic regimen. Corroborated (clamp, +26% muscle IS, "
               "explicit nulls for hepatic IS and insulin clearance) by "
               "Toledo et al. 2021, PMID:33528893 -- shared author "
               "(Toledo: first author there, last author on Wasko 2015), "
               "so not independent; it corroborates the DIRECTION, not "
               "the +20% magnitude. Wasko's own Table 2 ratio of group "
               "means is closer to +15% (~0.08 on 0.52) than the "
               "per-subject-mean +20% used here -- same direction, "
               "smaller size. Not registered: the same dataset's +45.4% "
               "beta-cell-function finding (KNOWN BUG #14 over-"
               "determination).")
        .build()
    )


# ---------------------------------------------------------------------------
# Growth hormone acute p3 dose-response constant (MASTER_PLAN 5.1, chunk 25,
# following Track 1 row 1.35) -- registers drug_interactions.py's own
# GROWTH_HORMONE_P3_FACTOR against the ALREADY-EXISTING "growth hormone"
# subject entity (registry_seed.py's own original v1 seed, force_numeric_id
# =12 -- NOT a new registration; an initial draft of this chunk attempted
# to re-register it and was correctly rejected by IDRegistry's own
# cross-reference-collision check before anything was committed).
# ---------------------------------------------------------------------------

def _seed_growth_hormone_p3(store: ParameterStore, registry: IDRegistry) -> None:
    """
    Registers drug_interactions.py's GROWTH_HORMONE_P3_FACTOR against
    the pre-existing "growth hormone" IDRegistry entity (UniProt
    P01241, registered since the original v1 seed -- no new entity
    needed for this chunk, only the parameter itself).

    GOLD -- a directly-quoted, primary-source percentage, not a back-
    solve. Same category as atorvastatin/hydrochlorothiazide's own
    directly-reported dose-response percentages elsewhere in this file.
    """
    gh_id = _resolve(registry, "growth hormone")

    # Fowelin J, Attvall S, von Schenck H, Smith U, Lager I.
    # "Characterization of the insulin-antagonistic effect of growth
    # hormone in man." Diabetologia. 1991;34(7):500-6. PMID:1916055 --
    # verified live via NCBI efetch this session. n=11 healthy subjects,
    # GH infused at 3 rates (6/12/24 mU/kg/h) for 1h, euglycaemic clamp +
    # D-(3-3H)-glucose vs. saline control. Directly quoted: "The insulin-
    # antagonistic effect of growth hormone started after approximately
    # 2 h, was maximal after 4-5 h (approximately 39% inhibition of
    # glucose infusion rate between control and growth hormone
    # 24 mU.kg-1.h-1)..." GROWTH_HORMONE_P3_FACTOR = 1.0 - 0.39 = 0.61.
    store.add(
        SeedParameterBuilder()
        .for_subject(gh_id)
        .kind(ParameterKind.CUSTOM, custom_name="growth_hormone_p3_factor")
        .value(1.0 - 0.39, unit="unitless multiplier on p3")
        .from_pubmed("1916055", figure_or_table="Results text")
        .in_human()
        .clinical()
        .gold()
        .notes("Acute insulin-antagonistic effect at the highest tested "
               "infusion rate (24 mU/kg/h). Directly quoted: "
               "\"approximately 39% inhibition of glucose infusion rate "
               "between control and growth hormone 24 mU.kg-1.h-1.\" "
               "This platform's own GROWTH_HORMONE_P3_FACTOR=0.61 "
               "matches (1.0-0.39) exactly. Corroborated directionally "
               "(no extractable magnitude) by Krag et al. 2007, "
               "PMID:17132823 -- see drug_interactions.py's own inline "
               "comment for the full cross-check.")
        .build()
    )


# ---------------------------------------------------------------------------
# Hydrochlorothiazide chronic p3 dose-response constant (MASTER_PLAN 5.1,
# chunk 26, following Track 1 row 1.34) -- registers drug_interactions.py's
# own THIAZIDE_P3_FACTOR_AT_TESTED_DOSE against the newly-registered
# "hydrochlorothiazide" subject entity (registry_seed.py).
# ---------------------------------------------------------------------------

def _seed_hydrochlorothiazide_p3(store: ParameterStore, registry: IDRegistry) -> None:
    """
    Registers drug_interactions.py's THIAZIDE_P3_FACTOR_AT_TESTED_DOSE.

    GOLD -- a direct ratio of two primary-source-reported clamp values
    (5.7/6.4), not a back-solve. Same category as pioglitazone's own
    2.7/1.6 ratio (also GOLD, both whole-body clamp-derived numbers, no
    unit conversion).
    """
    hctz_id = _resolve(registry, "hydrochlorothiazide")

    # Pollare T, Lithell H, Berne C. "A comparison of the effects of
    # hydrochlorothiazide and captopril on glucose and lipid metabolism
    # in patients with hypertension." N Engl J Med. 1989;321(13):868-73.
    # PMID:2671740 -- verified live via NCBI efetch this session.
    # Randomized, double-blind, crossover, n=50 essential hypertension
    # patients, two 4-month treatment periods, euglycemic clamp,
    # hydrochlorothiazide 40 +/- 12 mg/day. Directly quoted:
    # "hydrochlorothiazide caused a decrease from 6.4 +/- 2.0 to
    # 5.7 +/- 1.9 (P less than 0.01)" mg/kg.min insulin-mediated glucose
    # disposal. THIAZIDE_P3_FACTOR_AT_TESTED_DOSE = 5.7/6.4 = 0.891.
    store.add(
        SeedParameterBuilder()
        .for_subject(hctz_id)
        .kind(ParameterKind.CUSTOM, custom_name="thiazide_p3_factor_at_tested_dose")
        .value(5.7 / 6.4, unit="unitless multiplier on p3 at 40mg/day")
        .from_pubmed("2671740", figure_or_table="Results text")
        .in_human()
        .clinical()
        .gold()
        .notes("Chronic (4-month) insulin-mediated glucose disposal "
               "reduction at the study's own mean dose (40mg/day). "
               "Directly quoted: \"hydrochlorothiazide caused a "
               "decrease from 6.4 +/- 2.0 to 5.7 +/- 1.9 (P less than "
               "0.01).\" This platform's own "
               "THIAZIDE_P3_FACTOR_AT_TESTED_DOSE=0.891 matches "
               "(5.7/6.4) exactly. Largest and best-powered study "
               "(n=50) found for this mechanism.")
        .build()
    )


# ---------------------------------------------------------------------------
# Indinavir acute p3 dose-response constant (MASTER_PLAN 5.1, chunk 27,
# following Track 1 row 1.33) -- registers drug_interactions.py's own
# INDINAVIR_P3_FACTOR_AT_TESTED_DOSE against the newly-registered
# "indinavir" subject entity (registry_seed.py).
# ---------------------------------------------------------------------------

def _seed_indinavir_p3(store: ParameterStore, registry: IDRegistry) -> None:
    """
    Registers drug_interactions.py's INDINAVIR_P3_FACTOR_AT_TESTED_DOSE.

    GOLD -- a directly-quoted, primary-source PAIRED "average decrease"
    statistic, not a back-solve or a raw ratio of group means. Same
    category as growth hormone/hydrochlorothiazide's own directly-
    reported values elsewhere in this file.
    """
    indinavir_id = _resolve(registry, "indinavir")

    # Noor MA, Seneviratne T, Aweeka FT, et al. "Indinavir acutely
    # inhibits insulin-stimulated glucose disposal in humans: a
    # randomized, placebo-controlled study." AIDS. 2002;16(5):F1-8.
    # PMID:11964551 -- verified live via NCBI efetch this session. n=6
    # healthy HIV-seronegative adult men, randomized double-blind
    # crossover, single 1200mg oral indinavir vs. placebo, 180-min
    # euglycemic hyperinsulinemic clamp. Directly quoted: "Insulin-
    # stimulated glucose disposal per unit of insulin (M/I) decreased in
    # all subjects from 14.1 +/- 1.2 to 9.2 +/- 0.8 mg/kg.min per
    # microUI/ml... on indinavir (average decrease, 34.1 +/- 9.2%)."
    # INDINAVIR_P3_FACTOR_AT_TESTED_DOSE = 1.0 - 0.341 = 0.659 -- the
    # paper's own quoted PAIRED "average decrease" statistic, not the
    # raw ratio of group means (9.2/14.1=0.652, which corroborates
    # closely but is not what this platform's own constant uses).
    store.add(
        SeedParameterBuilder()
        .for_subject(indinavir_id)
        .kind(ParameterKind.CUSTOM, custom_name="indinavir_p3_factor_at_tested_dose")
        .value(1.0 - 0.341, unit="unitless multiplier on p3 at 1200mg")
        .from_pubmed("11964551", figure_or_table="Results text")
        .in_human()
        .clinical()
        .gold()
        .notes("Acute insulin-stimulated glucose disposal reduction at "
               "the study's own single tested dose (1200mg). Directly "
               "quoted: \"...on indinavir (average decrease, "
               "34.1 +/- 9.2%).\" This platform's own "
               "INDINAVIR_P3_FACTOR_AT_TESTED_DOSE=0.659 matches "
               "(1.0-0.341) exactly -- the paper's own PAIRED statistic, "
               "not the raw group-mean ratio (9.2/14.1=0.652, which "
               "corroborates within ~1 percentage point but is a "
               "different, less appropriate summary for this crossover "
               "design).")
        .build()
    )


# ---------------------------------------------------------------------------
# Atorvastatin chronic p3 dose-response constants (MASTER_PLAN 5.1,
# chunk 28, following Track 1 row 1.41) -- registers drug_interactions.py's
# own ATORVASTATIN_P3_FACTOR_AT_DOSE against the newly-registered
# "atorvastatin" subject entity (registry_seed.py).
# ---------------------------------------------------------------------------

def _seed_atorvastatin_p3(store: ParameterStore, registry: IDRegistry) -> None:
    """
    Registers the 3 statistically-significant dose-response anchor
    points from ATORVASTATIN_P3_FACTOR_AT_DOSE (20/40/80mg) -- GOLD,
    each a directly-reported within-arm percentage change, cross-
    checked live against the real atorvastatin_p3_factor() function at
    each exact anchor dose, not duplicated literals. Same multi-point
    registration pattern as tirzepatide's own 3 Gb anchor points.

    The 10mg point (1% decrease, p=0.312 vs. baseline) is deliberately
    NOT registered -- it is EMPIRICALLY NULL, not a small real effect,
    matching this platform's own established convention for a non-
    significant low-dose point (same reasoning as pioglitazone's own
    7.5/15mg null zone).
    """
    atorvastatin_id = _resolve(registry, "atorvastatin")

    from biosim_physiology.drug_interactions import atorvastatin_p3_factor

    # Koh KK, Quon MJ, Han SH, Lee Y, Kim SJ, Shin EK. "Atorvastatin
    # causes insulin resistance and increases ambient glycemia in
    # hypercholesterolemic patients." J Am Coll Cardiol. 2010;55(12):
    # 1209-16. PMID:20298928 -- verified live via NCBI efetch this
    # session. Randomized, single-blind, placebo-controlled, PARALLEL-
    # group study, n=44 placebo / 42-44 per atorvastatin arm
    # (10/20/40/80mg), 2-month treatment. Directly quoted: "Atorvastatin
    # 10, 20, 40, and 80 mg decreased insulin sensitivity (1%, 3%, 3%,
    # and 4%, respectively) when compared with either baseline (p =
    # 0.312, p = 0.008, p < 0.001, and p = 0.008, respectively...)."
    # 20/40/80mg are all significant vs. baseline; their own directly-
    # reported percentages (3%/3%/4%) are used as-is.
    dose_to_custom_name = {
        20.0: "p3_factor_at_20mg",
        40.0: "p3_factor_at_40mg",
        80.0: "p3_factor_at_80mg",
    }
    for dose_mg, custom_name in dose_to_custom_name.items():
        store.add(
            SeedParameterBuilder()
            .for_subject(atorvastatin_id)
            .kind(ParameterKind.CUSTOM, custom_name=custom_name)
            .value(
                atorvastatin_p3_factor(dose_mg),
                unit="unitless multiplier on p3",
            )
            .from_pubmed("20298928", figure_or_table="Results text")
            .in_human()
            .clinical()
            .gold()
            .notes(f"Directly-reported within-arm insulin-sensitivity "
                   f"decrease at the {dose_mg:.0f}mg daily dose (2-month "
                   f"treatment), one of 3 real significant anchor points "
                   f"this platform's own piecewise-linear "
                   f"atorvastatin_p3_factor() interpolates between. Two "
                   f"disclosed caveats (not hidden): PARALLEL-group "
                   f"design (not crossover, more between-subject "
                   f"variance) and SINGLE-blind (not double-blind) "
                   f"outcome assessment. The empirically NULL 10mg point "
                   f"(p=0.312 vs. baseline) is deliberately not "
                   f"registered as a separate anchor.")
            .build()
        )


# ---------------------------------------------------------------------------
# Sulfonylurea glucose-independent secretion (MASTER_PLAN 5.1, chunk 7) --
# registers physiology/drug_interactions.py's own
# SULFONYLUREA_SECRETION_FRACTION_OF_KE_IB against the "insulin" subject
# (already registered, shared with the exogenous-insulin PK entries above).
# ---------------------------------------------------------------------------

def _seed_sulfonylurea_secretion(store: ParameterStore, registry: IDRegistry) -> None:
    """
    Registers SULFONYLUREA_SECRETION_FRACTION_OF_KE_IB=0.5. Deliberately
    BRONZE, same category as the exogenous-insulin PK constants -- a
    real, disclosed engineering back-solve (Riefflin et al. 2015's own
    absolute pmol/min C-peptide secretion-rate figures are a different
    physical quantity from this model's plasma insulin concentration
    rate constant, and no distribution-volume/clearance relationship
    exists to convert between them), not a directly-reported number.
    """
    insulin_id = _resolve(registry, "insulin")

    # Riefflin A, Ayyagari U, Manley SE, Holman RR, Levy JC. "The effect
    # of glibenclamide on insulin secretion at normal glucose
    # concentrations." Diabetologia. 2015;58(1):43-9. PMID:25297572 --
    # verified live via NCBI efetch (see drug_interactions.py's own
    # module-level citation comment for the full derivation). RCT,
    # double-blind, placebo-controlled crossover, n=19 T2D, glibenclamide
    # 2.5mg twice daily x14 days vs. placebo x14 days, 3-step
    # hyperinsulinaemic glucose clamp. C-peptide secretion rate rose by
    # an ABSOLUTE increment "independent of clamp glucose concentration
    # (p=0.54)" -- genuinely glucose-independent, the textbook
    # KATP-channel mechanism behind sulfonylureas' hypoglycaemia risk.
    store.add(
        SeedParameterBuilder()
        .for_subject(insulin_id)
        .kind(ParameterKind.CUSTOM, custom_name="sulfonylurea_secretion_fraction_of_ke_ib")
        .value(0.5, unit="dimensionless (fraction of ke*Ib)")
        .from_pubmed("25297572", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .bronze()
        .notes("Glucose-independent insulin-secretion increment for an "
               "active sulfonylurea dose (glibenclamide), expressed as a "
               "fraction of this model's own ke*Ib basal-secretion "
               "reference. Riefflin et al. 2015's own quoted finding is "
               "genuinely glucose-independent (2015;58(1):43-9, C-peptide "
               "secretion rose 1.7-2.8-fold at every tested glucose "
               "level, all p<0.001, independent of clamp glucose "
               "p=0.54) -- but its absolute pmol/min C-peptide "
               "secretion-rate units are a different physical quantity "
               "from this model's plasma-insulin-concentration rate "
               "constant, with no reported distribution-volume/"
               "clearance relationship to convert between them. 0.5 is "
               "an isolated, disclosed engineering calibration chosen so "
               "the resulting effect is clearly large enough to "
               "represent Riefflin's real, clinically significant "
               "effect while staying below this model's own hard safety "
               "bound (PathwayParameters.__post_init__ requires "
               "glucose_independent_secretion < ke*Ib) -- same category "
               "as HPA_FEEDBACK_GAIN and the exogenous-insulin PK "
               "constants elsewhere in this store, not a directly-"
               "reported number. Not dose-scaled: Riefflin tested a "
               "single clinical dose (2.5mg twice daily); no dose-"
               "response curve was found to scale against (a disclosed "
               "limitation, MASTER_PLAN.md row 1.3b).")
        .build()
    )


# ---------------------------------------------------------------------------
# Pioglitazone p3 dose-response constant (MASTER_PLAN 5.1, chunk 9;
# Track 1 row 1.25) -- registers drug_interactions.py's
# PIOGLITAZONE_P3_MULTIPLIER at its 45mg (maximum labeled) reference
# dose against the newly-registered "pioglitazone" metabolite subject
# (registry_seed.py).
# ---------------------------------------------------------------------------

def _seed_pioglitazone_p3(store: ParameterStore, registry: IDRegistry) -> None:
    """
    Registers the whole-body insulin sensitivity (p3) multiplier at
    pioglitazone's 45mg reference dose. GOLD -- re-tiered from an
    earlier BRONZE assignment (housekeeping backlog item, MASTER_PLAN
    Section 8): the STORED value is the DIRECT ratio of the paper's own
    two directly-reported whole-body ISI numbers (2.7/1.6), with no
    unit conversion and no back-solved engineering constant involved
    (unlike acarbose/miglitol's kabs scales, which back-solve a
    different-unit engineering constant FROM a clinical ratio) --
    structurally the same evidentiary category as CKD's own
    covariate-adjusted p3_factor/ke_factor (chunk 11) and tirzepatide's
    own gb_factor_at_*mg (chunk 19), both already GOLD. Chunk 19's own
    mandatory code-reviewer pass flagged this exact inconsistency and
    judged GOLD the better-reasoned tier for this category; this entry
    is the resolution.
    """
    pioglitazone_id = _resolve(registry, "pioglitazone")

    # Miyazaki Y, Matsuda M, DeFronzo RA. "Dose-response effect of
    # pioglitazone on insulin sensitivity and insulin secretion in type
    # 2 diabetes." Diabetes Care. 2002;25(3):517-23. PMID:11874940 --
    # verified live via NCBI efetch this session (see
    # drug_interactions.py's own module-level citation comment for the
    # full derivation, including the 7.5/15mg empirically-null anchor
    # points this single reference-dose entry doesn't capture). 5-arm,
    # 26-week, diet-treated T2D RCT (n=58): "Whole-body ISI... increased
    # significantly in patients treated with... 45 mg (1.6 +/- 0.2 to
    # 2.7 +/- 0.6, P<0.05) per day of pioglitazone."
    store.add(
        SeedParameterBuilder()
        .for_subject(pioglitazone_id)
        .kind(ParameterKind.CUSTOM, custom_name="p3_multiplier_at_45mg_reference_dose")
        .value(2.7 / 1.6, unit="dimensionless (multiplier on PathwayParameters.p3)")
        .from_pubmed("11874940", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .gold()
        .notes("Whole-body insulin sensitivity index (ISI) ratio at "
               "pioglitazone's 45mg/day reference dose -- the direct "
               "post/pre ratio of Miyazaki et al. 2002's own reported "
               "values (2.7/1.6), not a back-solved engineering "
               "constant. GOLD (re-tiered from an earlier BRONZE "
               "assignment, MASTER_PLAN Section 8 housekeeping item) -- "
               "structurally the same category as CKD's own "
               "covariate-adjusted p3_factor/ke_factor and "
               "tirzepatide's own gb_factor_at_*mg, both a ratio "
               "computed directly from a primary source's own two "
               "reported numbers with no unit conversion. Does NOT "
               "capture this same paper's own empirically-null "
               "7.5/15mg anchor points or its 30mg anchor (2.5/1.8) -- "
               "see drug_interactions.py's pioglitazone_p3_factor() for "
               "the full 5-point piecewise-linear dose-response this "
               "single reference-dose entry is only one point of. "
               "Deliberately does NOT also register the same paper's "
               "hepatic-ISI or insulinogenic-index findings (would "
               "repeat KNOWN BUG #14's over-determination mistake, per "
               "drug_interactions.py's own docstring).")
        .build()
    )


# ---------------------------------------------------------------------------
# Metformin Gb reduction constant (MASTER_PLAN 5.1, chunk 17) --
# registers drug_interactions.py's own metformin_modifiers() Gb
# (fasting-glucose) reduction fraction against the newly-registered
# "metformin" metabolite subject (registry_seed.py). First IDRegistry
# (molecular) chunk since chunk 10 (protein) -- returns to Track 5 row
# 5.4's own "still many more small chunks available" backlog after six
# consecutive ContextRegistry (non-molecular) chunks.
# ---------------------------------------------------------------------------

def _seed_metformin_gb(store: ParameterStore, registry: IDRegistry) -> None:
    """
    Registers metformin's Gb (fasting-glucose) reduction fraction --
    BRONZE, a disclosed flux-to-concentration approximation, not a
    directly-reported concentration value itself.
    """
    metformin_id = _resolve(registry, "metformin")

    # DeFronzo RA, Barzilai N, Simonson DC. "Mechanism of metformin
    # action in obese and lean noninsulin-dependent diabetic subjects."
    # J Clin Endocrinol Metab. 1991;73(6):1294-301. PMID:8923861 --
    # verified live via NCBI efetch this session (see
    # drug_interactions.py's own metformin_modifiers() docstring for
    # the full derivation and the mechanism-weighting investigation it
    # resolves). Placebo-controlled, chronic (15wk), T2D, tracer-
    # measured: basal hepatic glucose production fell 12.9 -> 11.0
    # umol/kg.min (-14.7%, P<0.01).
    store.add(
        SeedParameterBuilder()
        .for_subject(metformin_id)
        .kind(ParameterKind.CUSTOM, custom_name="gb_reduction_fraction_at_max_dose")
        .value(0.147, unit="dimensionless (fractional reduction applied to PathwayParameters.Gb, dose-scaled)")
        .from_pubmed("8923861", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .bronze()
        .notes("BRONZE, not GOLD -- this is a disclosed flux-to-"
               "concentration APPROXIMATION, not a directly-reported "
               "concentration value. The primary source measures hepatic "
               "glucose PRODUCTION (a flux, umol/kg.min, -14.7%), but "
               "this platform's Gb is a steady-state glucose "
               "CONCENTRATION (mM) -- absent a more principled "
               "translation, drug_interactions.py's own "
               "metformin_modifiers() assumes a given %-reduction in "
               "hepatic production produces approximately the same "
               "%-reduction in the resulting basal concentration "
               "(holding clearance kinetics fixed), flagged there as an "
               "approximation to revisit, not an exact derivation. This "
               "0.147 value is the maximum (dose_factor=1.0, i.e. at or "
               "above the 2000mg dose cap) fractional reduction -- the "
               "real function applies a log-linear dose scaling below "
               "that cap, not captured by this single reference-value "
               "entry (same convention as pioglitazone's own single "
               "reference-dose entry above).")
        .build()
    )


# ---------------------------------------------------------------------------
# Caffeine p3/Gb dose-response constants (MASTER_PLAN 5.1, chunk 18) --
# registers drug_interactions.py's own caffeine_modifiers() constants
# against the newly-registered "caffeine" metabolite subject
# (registry_seed.py).
#
# Unlike metformin/pioglitazone (whose own dose-response has a clean
# dose value at which the scaling factor equals exactly 1.0), caffeine's
# own `concentration_fraction` (a product of a rising absorption
# exponential and a decaying elimination exponential, two different
# rate constants) never reaches exactly 1.0 at any time point -- it
# peaks at ~0.6406 at t~2.20h. This does NOT mean no live cross-check
# is possible, though (a mandatory code-reviewer pass correctly caught
# an earlier draft's wrong claim that it wasn't): `concentration_fraction`
# is a fully known closed-form function of t (both half-lives are
# literals in caffeine_modifiers()'s own source), so the test below
# calls the REAL function at a fixed t and BACK-SOLVES the underlying
# 0.15/0.02 constants by dividing out the independently-computed
# concentration_fraction(t) -- a genuine drift-detection guard, not a
# duplicated-literal comparison, matching every other chunk's own
# anti-drift convention.
# ---------------------------------------------------------------------------

def _seed_caffeine_p3(store: ParameterStore, registry: IDRegistry) -> None:
    """
    Registers caffeine's p3 dose-response amplitude (BRONZE, cited). An
    uncited Gb increase amplitude used to be registered too; it was removed
    with the term itself on 2026-11-20.
    """
    caffeine_id = _resolve(registry, "caffeine")

    # Keijzers GB, De Galan BE, Tack CJ, Smits P. "Caffeine can decrease
    # insulin sensitivity in humans." Diabetes Care. 2002;25(2):364-9.
    # PMID:11815511 -- verified live via NCBI efetch this session (see
    # drug_interactions.py's own caffeine_modifiers() docstring, which
    # cites this same PMID for "~15% reduction at 200 mg"). The stored
    # value (0.15) is caffeine_modifiers()'s own inline literal
    # (`max_reduction = 0.15 * dose_factor`, with dose_factor==1.0
    # exactly at the 200mg reference dose used in the docstring's own
    # claim) -- cross-checked live below (test_seeds.py), not just
    # confirmed by direct code reading.
    store.add(
        SeedParameterBuilder()
        .for_subject(caffeine_id)
        .kind(ParameterKind.CUSTOM, custom_name="p3_reduction_at_200mg_reference_dose")
        .value(0.15, unit="dimensionless (fractional reduction applied to PathwayParameters.p3 at peak concentration)")
        .from_pubmed("11815511", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .bronze()
        .notes("caffeine_modifiers()'s own dose-response amplitude at "
               "its 200mg reference dose (where dose_factor==1.0 "
               "exactly) -- this is the maximum PEAK reduction; the "
               "real function additionally scales this by a "
               "pharmacokinetic concentration_fraction (0 to <1, "
               "absorption t1/2=0.75h, elimination t1/2=5h) that is "
               "NOT captured by this single reference-value entry, "
               "matching pioglitazone's/metformin's own single "
               "reference-dose convention.")
        .build()
    )
    # (The caffeine "gb_increase_amplitude" 0.02 entry was removed on
    # 2026-11-20 with the term itself; see caffeine_modifiers().)


# ---------------------------------------------------------------------------
# Tirzepatide Gb dose-response constants (MASTER_PLAN 5.1, chunk 19;
# Track 1 row 1.5) -- registers drug_interactions.py's own
# TIRZEPATIDE_GB_DOSE_MG/TIRZEPATIDE_GB_MULTIPLIER anchor points against
# the newly-registered "tirzepatide" metabolite subject
# (registry_seed.py). All three values are cross-checked live below via
# tirzepatide_gb_factor() called at each exact anchor dose (5/10/15mg),
# not duplicated literals -- if either the anchor values or the
# interpolation logic ever drifts, this breaks immediately.
# ---------------------------------------------------------------------------

def _seed_tirzepatide_gb(store: ParameterStore, registry: IDRegistry) -> None:
    """
    Registers all three Gb dose-response anchor points -- GOLD, each a
    directly-reported per-arm fractional FPG reduction from ONE primary
    source (same "computed ratio from the paper's own reported values"
    category as CKD's own p3/ke factors), not a back-solved estimate.
    """
    tirzepatide_id = _resolve(registry, "tirzepatide")

    from biosim_physiology.drug_interactions import tirzepatide_gb_factor

    # Furihata K, Mimura H, Urva S, Oura T, Ohwaki K, Imaoka T. "A phase
    # 1 multiple-ascending dose study of tirzepatide in Japanese
    # participants with type 2 diabetes." Diabetes Obes Metab.
    # 2022;24(2):239-46. PMID:34647404 -- verified live via NCBI efetch
    # this session (see drug_interactions.py's own citation comment for
    # the full derivation). Once-weekly subcutaneous tirzepatide, 3 dose
    # arms (5/10/15mg): each arm's own fractional fasting-plasma-glucose
    # reduction relative to its own baseline (0.680/0.619/0.596) --
    # DIRECTLY-REPORTED within-arm ratios, not a parametric dose-
    # response curve fit (the code's own docstring is explicit that a
    # piecewise-LINEAR interpolation across exactly these 3 points is
    # used specifically because no single one of them alone can
    # validate an assumed curve shape).
    dose_to_custom_name = {
        5.0:  "gb_factor_at_5mg",
        10.0: "gb_factor_at_10mg",
        15.0: "gb_factor_at_15mg",
    }
    for dose_mg, custom_name in dose_to_custom_name.items():
        store.add(
            SeedParameterBuilder()
            .for_subject(tirzepatide_id)
            .kind(ParameterKind.CUSTOM, custom_name=custom_name)
            .value(
                tirzepatide_gb_factor(dose_mg),
                unit="dimensionless (multiplier on PathwayParameters.Gb)",
            )
            .from_pubmed("34647404", figure_or_table="text, Results")
            .in_human()
            .clinical()
            .gold()
            .notes(f"Directly-reported per-arm fractional FPG reduction at "
                   f"the {dose_mg:.0f}mg once-weekly dose, one of 3 real "
                   f"anchor points this platform's own piecewise-linear "
                   f"tirzepatide_gb_factor() interpolates between -- this "
                   f"single entry does NOT capture the interpolated shape "
                   f"between anchors, the sub-5mg linear-to-neutral ramp, "
                   f"or the above-15mg clamp (see that function's own "
                   f"docstring for the full policy). Does NOT capture "
                   f"tirzepatide's OTHER real mechanism (a stronger "
                   f"incretin_effect via GLP1Pathway) -- disclosed there "
                   f"as genuinely uncited (Heise 2022/Mather 2024, both "
                   f"paywalled), a deliberate partial implementation.")
            .build()
        )


# ---------------------------------------------------------------------------
# Glipizide secretion-fraction constant (MASTER_PLAN 5.1, chunk 21;
# Track 1 row 1.3b) -- registers drug_interactions.py's own
# GLIPIZIDE_SECRETION_FRACTION_OF_KE_IB against the "insulin" subject.
#
# A genuine gap found while scoping chunk 21 (originally planned as a
# pure documentation-completeness chunk, registering glibenclamide/
# glipizide's own IDRegistry entities): glibenclamide's own secretion
# fraction (chunk 7) was already curated, but glipizide's own -- a
# real, cited constant, already actively used in production code since
# Track 1 row 1.3b -- had NEVER been registered in ParameterStore at
# all. Closed here, under the SAME "insulin" subject as glibenclamide's
# own entry (not under the new "glipizide" entity), for consistency --
# both fractions are expressed relative to insulin's own ke*Ib basal-
# secretion reference, not an inherent chemical property of either drug.
# ---------------------------------------------------------------------------

def _seed_glipizide_secretion(store: ParameterStore, registry: IDRegistry) -> None:
    """
    Registers glipizide's own glucose-independent secretion fraction --
    BRONZE, a disclosed engineering compromise (capped below this
    platform's own hard safety ceiling), motivated by but not equal to
    a real head-to-head relative-ratio finding, same category as
    glibenclamide's own chunk-7 entry.
    """
    insulin_id = _resolve(registry, "insulin")

    from biosim_physiology.drug_interactions import (
        GLIPIZIDE_SECRETION_FRACTION_OF_KE_IB,
    )

    # Groop L, Luzi L, Melander A, Groop PH, Ratheiser K, Simonson DC,
    # DeFronzo RA. "Different effects of glyburide and glipizide on
    # insulin secretion and hepatic glucose production in normal and
    # NIDDM subjects." Diabetes. 1987;36(11):1320-8. PMID:3117609 --
    # verified live via NCBI efetch this session (see
    # drug_interactions.py's own module-level citation comment for the
    # full derivation, including the corroborating Groop et al. 1991,
    # PMID:1777653). IV-infused glyburide (glibenclamide) vs. glipizide
    # at MATCHED plasma drug concentrations, n=10 healthy + 6 NIDDM:
    # "the acute phase (0-10 min) of plasma insulin and C-peptide
    # increase in both control and NIDDM subjects was twice as great
    # with GZ [glipizide] compared with GB [glyburide]" -- a real,
    # T2D-inclusive, head-to-head relative ratio (glipizide ~2x
    # glibenclamide), used to motivate (not exactly reproduce) the
    # stored value: doubling glibenclamide's own 0.5 fraction to 1.0
    # would collide with this platform's own hard safety bound
    # (glucose_independent_secretion must be STRICTLY < ke*Ib), so it
    # is capped at 0.9 -- as large as the safety ceiling allows while
    # still representing a meaningfully larger effect than
    # glibenclamide's own 0.5.
    store.add(
        SeedParameterBuilder()
        .for_subject(insulin_id)
        .kind(ParameterKind.CUSTOM, custom_name="glipizide_secretion_fraction_of_ke_ib")
        .value(GLIPIZIDE_SECRETION_FRACTION_OF_KE_IB, unit="dimensionless (fraction of ke*Ib)")
        .from_pubmed("3117609", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .bronze()
        .notes("Disclosed engineering compromise, NOT a directly-"
               "reported fraction -- motivated by Groop et al. 1987's "
               "own head-to-head '~2x glibenclamide' relative finding, "
               "but capped at 0.9 (not doubled to 1.0) to stay strictly "
               "below this platform's own hard safety bound on "
               "glucose_independent_secretion (must be < ke*Ib). Also "
               "corroborated by Groop et al. 1991 (PMID:1777653, a "
               "healthy-cohort dose-response study finding glipizide's "
               "own effect is glucose-independent, the same mechanism "
               "shape as glibenclamide's own chunk-7 entry). Cross-"
               "checked live against the real GLIPIZIDE_SECRETION_"
               "FRACTION_OF_KE_IB module constant, not a duplicated "
               "literal.")
        .build()
    )
