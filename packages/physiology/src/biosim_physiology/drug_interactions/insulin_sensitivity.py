# SAVE TO: packages\physiology\src\biosim_physiology\drug_interactions\insulin_sensitivity.py
"""Drugs modelled as a static insulin-sensitivity (p3) factor."""

from __future__ import annotations

import math

from ._core import DrugDose

# Thiazolidinedione (pioglitazone), a PPAR-gamma insulin sensitizer --
# mechanistically distinct from every other drug class this platform
# models (not a secretagogue, incretin, glucosuric agent, or absorption
# blocker). Found via a fresh Track 1 scoping fan-out (MASTER_PLAN row
# 1.25) after every previously-named Track 1 candidate was done or
# blocked; verified live via PubMedFetcher before use.
#
# Miyazaki Y, Matsuda M, DeFronzo RA. "Dose-response effect of
# pioglitazone on insulin sensitivity and insulin secretion in type 2
# diabetes." Diabetes Care. 2002;25(3):517-23. PMID:11874940. A 5-arm
# (placebo/7.5/15/30/45mg), 26-week, randomized, DIET-TREATED T2D cohort
# (n=58, no sulfonylurea background confound). Confirmed verbatim:
# "Patients treated with 7.5 or 15 mg/day of pioglitazone had NO CHANGE
# in fasting plasma glucose (FPG) and fasting plasma insulin (FPI)
# concentrations or in plasma glucose (PG) and insulin concentrations
# during the OGTT." And: "Whole-body ISI... increased significantly in
# patients treated with 30 mg (1.8 +/- 0.3 to 2.5 +/- 0.3, P<0.05) or
# 45 mg (1.6 +/- 0.2 to 2.7 +/- 0.6, P<0.05) per day of pioglitazone."
#
# Route chosen deliberately: whole-body ISI is a dimensionless ratio
# (needs no baseline conversion) and maps directly onto p3, the way
# every other insulin-sensitivity modifier in this file does. This is
# the ONLY mechanism implemented from this paper -- the same source
# also reports a hepatic-ISI improvement (0.13->0.21 at 30mg) and an
# insulinogenic-index rise (compensatory secretion), but stacking all
# three from one dataset would repeat KNOWN BUG #14's over-determination
# mistake (three independently-chosen "improvement" multipliers from a
# single study, none separately identifiable) -- deliberately deferred,
# not an oversight.
#
# A Gb-additive route was explicitly considered and REJECTED: applying
# the paper's own FPG deltas (-66/-97 mg/dl) additively to this
# platform's T2D preset (Gb=6.5mM, G_threshold=4.5mM) would have 54-68%
# of the calibrated effect swallowed by MIN_GB_THRESHOLD_MARGIN_MM --
# the same composition-boundary failure tirzepatide's row 1.5 hit before
# row 1.23's clamp_warnings fix. p3 is not Gb/G_threshold-adjacent, so
# no such clamp interaction applies to the route used here.
PIOGLITAZONE_P3_DOSE_MG        = (0.0, 7.5, 15.0, 30.0, 45.0)
PIOGLITAZONE_P3_MULTIPLIER     = (1.0, 1.0, 1.0, 2.5 / 1.8, 2.7 / 1.6)


def pioglitazone_p3_factor(dose_mg: float) -> float:
    """
    Multiplicative factor on PathwayParameters.p3 for a given daily
    pioglitazone dose (mg), piecewise-linearly interpolated across the
    five real Miyazaki et al. 2002 (PMID:11874940) anchor points (see
    this module's own citation comment directly above) -- including the
    two EMPIRICALLY NULL points (7.5mg, 15mg), not an assumed monotonic
    curve. Neutral (1.0) at dose_mg=0 through 15mg; ramps between 15mg
    and 45mg; clamped at the 45mg value above 45mg (the drug's own
    labeled maximum daily dose -- real clinical dosing does not exceed
    this, so this is a deliberate "don't extrapolate past tested data"
    cap, matching tirzepatide_gb_factor()'s own convention).
    """
    if not math.isfinite(dose_mg):
        raise ValueError(f"dose_mg must be finite, got {dose_mg}")
    if dose_mg < 0:
        raise ValueError(f"dose_mg cannot be negative, got {dose_mg}")
    doses = PIOGLITAZONE_P3_DOSE_MG
    mults = PIOGLITAZONE_P3_MULTIPLIER
    if dose_mg <= doses[0]:
        return mults[0]
    if dose_mg >= doses[-1]:
        return mults[-1]
    for i in range(len(doses) - 1):
        if doses[i] <= dose_mg <= doses[i + 1]:
            t = (dose_mg - doses[i]) / (doses[i + 1] - doses[i])
            return mults[i] + (mults[i + 1] - mults[i]) * t
    raise AssertionError("unreachable")  # pragma: no cover


def pioglitazone_modifiers(dose: DrugDose) -> dict[str, float]:
    """
    Pioglitazone (thiazolidinedione, PPAR-gamma insulin sensitizer) --
    see PIOGLITAZONE_P3_DOSE_MG's own docstring above for the full
    citation and derivation. Only p3 is modified (deliberately, not an
    oversight -- see that comment for why the hepatic/secretion pieces
    of the same paper are not also stacked here).
    """
    if dose.dose_mg <= 0:
        return {}
    return {"p3": pioglitazone_p3_factor(dose.dose_mg)}


# ---------------------------------------------------------------------------
# Indinavir (HIV protease inhibitor) -- Track 1 row 1.33. Scoped by a
# computational-biologist fan-out during a fresh Track 1 roadmap pass
# (MASTER_PLAN Section 2a, candidate A1), independently re-verified live
# via NCBI efetch before implementation (this session, not trusted from
# the scoping report alone).
#
# Noor MA, Seneviratne T, Aweeka FT, Lo JC, Schwarz JM, Mulligan K,
# Schambelan M, Grunfeld C. "Indinavir acutely inhibits insulin-stimulated
# glucose disposal in humans: a randomized, placebo-controlled study."
# AIDS. 2002;16(5):F1-8. PMID:11964551. Randomized, double-blind,
# crossover, n=6 healthy HIV-seronegative adult men, single 1200mg oral
# indinavir vs. placebo, 180-min euglycemic hyperinsulinemic clamp.
# Confirmed verbatim: "Insulin-stimulated glucose disposal per unit of
# insulin (M/I) decreased in all subjects from 14.1 +/- 1.2 to 9.2 +/- 0.8
# mg/kg.min per microUI/ml (95% confidence interval for change, 3.7-6.1;
# P < 0.001) on indinavir (average decrease, 34.1 +/- 9.2%)."
#
# INDINAVIR_P3_FACTOR anchored to the paper's own quoted PAIRED
# "average decrease" statistic (1 - 0.341 = 0.659), not the raw ratio of
# GROUP MEANS (9.2/14.1 = 0.652) -- the paired statistic is the more
# statistically appropriate summary of a within-subject crossover design,
# though the two independently corroborate each other closely (within
# ~1 percentage point).
#
# Best population match of any candidate found in the Section 2a sweep:
# healthy, non-obese adult men, exactly this platform's own calibration
# cohort (unlike CKD/GDM/T1D/T2D, all of which have a real population-
# mismatch caveat of their own). secretion_scale/NEFA both explicitly
# NOT touched -- the same source directly measured "Free fatty acid
# levels were not significantly different at baseline and were
# suppressed equally with insulin administration during both studies,"
# ruling out a NEFAPathway double-count by the study's own data (the
# same discipline row 1.32 used for its own cortisol null). Zero
# double-counting risk confirmed: no existing reference to indinavir/
# protease inhibitors anywhere in this codebase prior to this chunk
# (grepped `packages/` directly).
#
# Single ACUTE dose only -- no dose-response curve exists in this
# source (only one dose, 1200mg, was tested). Modeled as a piecewise-
# linear ramp from neutral (1.0) at dose_mg=0 to the measured factor at
# 1200mg, clamped (not extrapolated) above 1200mg -- matching this
# module's own established "don't extrapolate past tested data"
# convention (PIOGLITAZONE_P3_DOSE_MG above, sleep_modifiers() in
# modifiers.py). n=6 -- disclosed as a real evidence-strength caveat,
# not hidden.
INDINAVIR_TESTED_DOSE_MG = 1200.0
INDINAVIR_P3_FACTOR_AT_TESTED_DOSE = 1.0 - 0.341  # paired "average decrease" statistic


def indinavir_p3_factor(dose_mg: float) -> float:
    """
    Multiplicative factor on PathwayParameters.p3 for a single acute
    indinavir dose (mg), linearly interpolated between neutral (1.0) at
    dose_mg=0 and the measured factor (0.659) at the paper's own single
    tested dose (1200mg) -- see INDINAVIR_TESTED_DOSE_MG's own docstring
    above. Clamped (not extrapolated) above 1200mg.
    """
    if not math.isfinite(dose_mg):
        raise ValueError(f"dose_mg must be finite, got {dose_mg}")
    if dose_mg < 0:
        raise ValueError(f"dose_mg cannot be negative, got {dose_mg}")
    t = min(dose_mg, INDINAVIR_TESTED_DOSE_MG) / INDINAVIR_TESTED_DOSE_MG
    return 1.0 + (INDINAVIR_P3_FACTOR_AT_TESTED_DOSE - 1.0) * t


def indinavir_modifiers(dose: DrugDose) -> dict[str, float]:
    """
    Indinavir (HIV protease inhibitor) -- see INDINAVIR_TESTED_DOSE_MG's
    own docstring above for the full citation and derivation. Only p3 is
    modified; secretion_scale is deliberately NOT touched (the source's
    own directly-measured NEFA null rules out a NEFAPathway route, and
    the source reports no secretion-magnitude finding at all).
    """
    if dose.dose_mg <= 0:
        return {}
    return {"p3": indinavir_p3_factor(dose.dose_mg)}


# ---------------------------------------------------------------------------
# Hydrochlorothiazide (thiazide diuretic) -- Track 1 row 1.34. Scoped by the
# same computational-biologist fan-out as indinavir above (MASTER_PLAN
# Section 2a, candidate A3), independently re-verified live via NCBI
# efetch before implementation.
#
# Pollare T, Lithell H, Berne C. "A comparison of the effects of
# hydrochlorothiazide and captopril on glucose and lipid metabolism in
# patients with hypertension." N Engl J Med. 1989;321(13):868-73.
# PMID:2671740. Randomized, double-blind, CROSSOVER, n=50 patients with
# essential hypertension, two 4-month treatment periods, euglycemic
# clamp, hydrochlorothiazide dose 40 +/- 12 mg/day. Confirmed verbatim:
# "hydrochlorothiazide caused a decrease from 6.4 +/- 2.0 to 5.7 +/- 1.9
# (P less than 0.01)" mg/kg.min insulin-mediated glucose disposal --
# THIAZIDE_P3_FACTOR = 5.7/6.4 = 0.891. Largest and best-powered study
# (n=50) found across the entire Section 2a forward-scoping sweep.
#
# This platform currently has NO antihypertensive drug class at all --
# genuinely new territory, not an extension of an existing mechanism.
#
# secretion_scale explicitly NOT touched, DESPITE the same source
# reporting a real, positive finding there -- a deliberate scope
# decision, not an oversight: "Hydrochlorothiazide increased the basal
# insulin concentration and the late insulin response to glucose." This
# is almost certainly the COMPENSATORY hyperinsulinemic response to the
# measured insulin resistance itself, not an independent secretory
# mechanism -- modelling both the p3 reduction AND a secretion_scale
# increase from the same single finding would repeat KNOWN BUG #14's
# over-determination mistake (three independently-chosen multipliers
# from one dataset, none separately identifiable/causally distinct).
#
# Real caveats, disclosed rather than hidden:
#   - Population is essential HYPERTENSIVES, not healthy adults -- a
#     genuine but milder population mismatch than CKD's own (this
#     platform's core p1/p2/p3/ke calibration cohort is healthy/normal-
#     tolerant, per Bergman 1981/Pacini 1986). Corrected on review
#     (the abstract's own exact wording, verbatim): "Captopril increased
#     the insulin-mediated disposal of glucose, as compared with
#     placebo, from 5.7... to 6.3... (P less than 0.05), whereas
#     hydrochlorothiazide caused a decrease from 6.4... to 5.7... (P
#     less than 0.01)" -- "as compared with placebo" grammatically
#     attaches to the captopril clause; whether it also governs the
#     hydrochlorothiazide clause by parallel sentence structure is
#     genuinely ambiguous from the abstract alone, not a settled fact
#     -- do not overclaim a placebo-referenced design for HCTZ
#     specifically. The lipid-outcome sentence later in the same
#     abstract DOES explicitly compare hydrochlorothiazide "as compared
#     with placebo," which is suggestive but not dispositive for the
#     glucose-disposal figure this constant is actually derived from.
#   - The mechanism is very likely POTASSIUM-mediated, not a direct
#     thiazide effect on insulin action: Helderman JH, Elahi D, Andersen
#     DK, et al. "Prevention of the glucose intolerance of thiazide
#     diuretics by maintenance of body potassium." Diabetes.
#     1983;32(2):106-11. PMID:6337892 (title independently confirmed
#     live via NCBI esummary this session -- full abstract not fetched,
#     cited here only for its own well-established title/finding, not a
#     magnitude this platform borrows) directly implies the effect is
#     conditional on K+ depletion. This platform has no potassium state
#     to represent that conditionality, so this modifier represents the
#     UNSUPPLEMENTED clinical regime (the one Pollare et al. actually
#     studied) by default -- disclosed, not hidden.
#
# Zero double-counting risk confirmed: grepped `packages/` directly --
# no existing reference to hydrochlorothiazide/thiazide/diuretic
# anywhere in this codebase prior to this chunk.
#
# Single tested chronic maintenance dose only (40mg/day mean) -- no
# dose-response curve exists in this source. Piecewise-linear ramp from
# neutral (1.0) at dose_mg=0 to the measured factor at 40mg, clamped
# (not extrapolated) above 40mg -- matching this module's own
# established "don't extrapolate past tested data" convention
# (INDINAVIR_TESTED_DOSE_MG above, PIOGLITAZONE_P3_DOSE_MG).
THIAZIDE_TESTED_DOSE_MG = 40.0
THIAZIDE_P3_FACTOR_AT_TESTED_DOSE = 5.7 / 6.4


def hydrochlorothiazide_p3_factor(dose_mg: float) -> float:
    """
    Multiplicative factor on PathwayParameters.p3 for a given chronic
    hydrochlorothiazide daily dose (mg), linearly interpolated between
    neutral (1.0) at dose_mg=0 and the measured factor (0.891) at the
    paper's own single tested mean dose (40mg/day) -- see
    THIAZIDE_TESTED_DOSE_MG's own docstring above. Clamped (not
    extrapolated) above 40mg.
    """
    if not math.isfinite(dose_mg):
        raise ValueError(f"dose_mg must be finite, got {dose_mg}")
    if dose_mg < 0:
        raise ValueError(f"dose_mg cannot be negative, got {dose_mg}")
    t = min(dose_mg, THIAZIDE_TESTED_DOSE_MG) / THIAZIDE_TESTED_DOSE_MG
    return 1.0 + (THIAZIDE_P3_FACTOR_AT_TESTED_DOSE - 1.0) * t


def hydrochlorothiazide_modifiers(dose: DrugDose) -> dict[str, float]:
    """
    Hydrochlorothiazide (thiazide diuretic) -- see
    THIAZIDE_TESTED_DOSE_MG's own docstring above for the full citation
    and derivation. Only p3 is modified; secretion_scale is
    deliberately NOT touched despite a real positive finding in the
    same source -- see that comment for why (compensatory response,
    not an independent mechanism -- would repeat KNOWN BUG #14).
    """
    if dose.dose_mg <= 0:
        return {}
    return {"p3": hydrochlorothiazide_p3_factor(dose.dose_mg)}


# ---------------------------------------------------------------------------
# Growth hormone (somatropin) -- Track 1 row 1.35. Scoped by the same
# computational-biologist fan-out as indinavir/hydrochlorothiazide above
# (MASTER_PLAN Section 2a, candidate A4), independently re-verified live
# via NCBI efetch before implementation.
#
# Fowelin J, Attvall S, von Schenck H, Smith U, Lager I. "Characterization
# of the insulin-antagonistic effect of growth hormone in man."
# Diabetologia. 1991;34(7):500-6. PMID:1916055. n=11 healthy subjects, GH
# infused at 3 rates (6/12/24, see unit note below) for 1h, euglycemic
# clamp + D-(3-3H)-glucose, vs. a saline control. Confirmed verbatim:
# "The insulin-antagonistic effect of growth hormone started after
# approximately 2 h, was maximal after 4-5 h (approximately 39%
# inhibition of glucose infusion rate between control and growth hormone
# 24 mU.kg-1.h-1) and lasted for 6-7 h after peak levels... Growth
# hormone infusion of 12 mU.kg-1.h-1 induced a similar insulin-
# antagonistic effect as the higher infusion rate whereas 6 mU.kg-1.h-1
# induced a smaller response." GROWTH_HORMONE_P3_FACTOR = 1 - 0.39 = 0.61.
#
# Corroborated directionally (no extractable magnitude) by Krag MB,
# Jessen N, Lange KH, Vestergaard PF, Christiansen JS, Schmitz O,
# Jorgensen JO. "Growth hormone-induced insulin resistance is associated
# with increased intramyocellular triglyceride content but unaltered
# VLDL-triglyceride kinetics." Am J Physiol Endocrinol Metab.
# 2007;292(3):E920-7. PMID:17132823 (independently re-verified live,
# title/author/journal/year confirmed) -- n=9 healthy men, 8-day GH
# 2mg/day, randomized double-blind placebo-controlled crossover:
# "Insulin resistance was induced by GH irrespective of concomitant
# acipimox (P<0.001)." Used only as directional corroboration, not a
# second magnitude -- the abstract itself reports no extractable p3-
# comparable number.
#
# DELIBERATE DESIGN DECISION: NO dose-response scaling is attempted here
# at all -- unlike every other acute drug mechanism in this file
# (indinavir, hydrochlorothiazide), this is modelled as a single flat
# peak-effect factor with NO ramp function, applied whenever GH is
# present (dose_mg > 0) regardless of the value. Two independent reasons,
# either one alone sufficient:
#   1. The source's own abstract is internally UNIT-INCONSISTENT: the
#      Methods sentence states "three different infusion rates (6, 12 or
#      24 mU.kg-1.min-1)" while the Results sentence states "growth
#      hormone 6 mU.kg-1.h-1... 12 mU.kg-1.h-1... 24 mU.kg-1.h-1" for the
#      SAME three numbers -- a 60x unit discrepancy in the same abstract.
#      A dose-scaled implementation would have to silently pick one
#      reading, manufacturing false precision from an ambiguous source --
#      this project's own established discipline (see MASTER_PLAN row
#      1.30's identical "do not adopt either reading without the full
#      text" stance) rejects exactly this move.
#   2. Even setting the unit question aside, the reported shape is
#      SATURATING, not linear: "Growth hormone infusion of 12 mU.kg-1.h-1
#      induced a SIMILAR insulin-antagonistic effect as the higher
#      infusion rate" -- 12 and 24 produce the same effect, so a naive
#      linear ramp between 0 and a single "reference dose" (the
#      indinavir/hydrochlorothiazide convention) would misrepresent this
#      specific source's own reported shape regardless of units.
# This platform's own `DrugDose.dose_mg` field cannot honestly represent
# an ambiguous-unit infusion RATE as a dose in milligrams -- rather than
# force a mapping that would mislead a REST client sending a real mg
# value, `dose_mg` is used here purely as a PRESENCE flag (>0 = GH
# administered), not a scaling axis. Correction (mandatory code-reviewer
# finding): this applies the study's own MAXIMAL/24 mU.kg-1.h-1-derived
# effect (the only tier with an extractable magnitude) at ANY nonzero
# dose, including a trace amount -- NOT a conservative lower-bound
# effect. The 6 mU.kg-1.h-1 tier's own "smaller response" is qualitative
# only, no magnitude reported, so there is no real number available to
# use as a low-dose anchor instead -- inventing one would fabricate
# precision this source doesn't support (the same row-1.30 precedent
# cited above). If a genuinely dose-scaled implementation is ever
# wanted, it needs the full text of PMID:1916055 first (to resolve the
# unit question AND obtain the 6 mU.kg-1.h-1 tier's own magnitude), not
# a guess at either.
#
# secretion_scale/hepatic_egp_delta explicitly NOT touched: the source's
# own text attributes the resistance to "a less pronounced insulin
# effect both to inhibit glucose production and to stimulate glucose
# utilization" -- i.e. a real hepatic component exists too, but with NO
# separate extractable magnitude, so adding a second, uncited
# hepatic_egp_delta contribution from the same single finding would
# repeat KNOWN BUG #14's over-determination mistake. Model p3 only.
#
# Relationship to acromegaly (Track 1 row 1.37, chronic GH excess --
# structurally the SAME underlying mechanism at a different exposure
# duration): deliberately NON-OVERLAPPING by construction, not by
# coincidence -- this drug-dose mechanism touches ONLY p3 (no Gb
# magnitude was extractable from this acute-infusion source), while
# acromegaly's own citation (Moller et al. 1992, PMID:1569148) supplies
# ONLY a Gb magnitude (no p3 magnitude was extractable there either,
# only a qualitative "substantially resistant to insulin" finding this
# project declined to transfer to p3 without a real number -- see
# DiseaseState.ACROMEGALY's own reasoning in modifiers.py, now
# implemented). A scenario combining both composes p3 (from this drug
# event) and Gb (from the disease state) multiplicatively without
# double-counting either magnitude, since neither number is borrowed
# from the other study.
GROWTH_HORMONE_P3_FACTOR = 1.0 - 0.39


def growth_hormone_modifiers(dose: DrugDose) -> dict[str, float]:
    """
    Growth hormone (somatropin) -- see GROWTH_HORMONE_P3_FACTOR's own
    docstring above for the full citation, the unit-ambiguity finding,
    and why NO dose-response scaling is attempted (a flat peak-effect
    factor applied whenever dose_mg > 0, not a ramp). Only p3 is
    modified; secretion_scale/hepatic_egp_delta deliberately NOT touched
    (see that comment for why).
    """
    if dose.dose_mg <= 0:
        return {}
    return {"p3": GROWTH_HORMONE_P3_FACTOR}


# ---------------------------------------------------------------------------
# Atorvastatin (HMG-CoA reductase inhibitor) -- Track 1 row 1.41. Flagged
# by an earlier Track 1 scoping fan-out (the same one that surfaced
# alcohol/row 1.31) as priority candidate #3, but never given its own
# MASTER_PLAN row or independent re-verification until now -- the
# citation below was fetched live via NCBI efetch THIS session, not
# trusted from the earlier informal report.
#
# Koh KK, Quon MJ, Han SH, Lee Y, Kim SJ, Shin EK. "Atorvastatin causes
# insulin resistance and increases ambient glycemia in
# hypercholesterolemic patients." J Am Coll Cardiol. 2010;55(12):
# 1209-16. PMID:20298928. Randomized, single-blind, placebo-controlled,
# PARALLEL-group (not crossover, unlike most other citations in this
# file) study, n=44 placebo / 42/44/43/40 for atorvastatin 10/20/40/80mg
# respectively, 2-month treatment period. Confirmed verbatim:
# "Atorvastatin 10, 20, 40, and 80 mg decreased insulin sensitivity (1%,
# 3%, 3%, and 4%, respectively) when compared with either baseline (p =
# 0.312, p = 0.008, p < 0.001, and p = 0.008, respectively, by paired t
# test) or placebo (p = 0.033 by ANOVA)."
#
# ATORVASTATIN_P3_FACTOR_AT_DOSE: the 10mg point (1% decrease) is
# EMPIRICALLY NULL vs. baseline (p=0.312, not significant) -- treated as
# neutral (1.0), not a small nonzero effect, matching this file's own
# established convention for a genuinely non-significant low-dose point
# (PIOGLITAZONE_P3_MULTIPLIER's own 7.5/15mg null zone, from the same
# Miyazaki 2002 source's own non-significant low-dose arms). 20/40/80mg
# are all significant vs. baseline (p=0.008, p<0.001, p=0.008
# respectively) -- their own directly-reported percentages (3%, 3%, 4%)
# are used as-is.
#
# GOLD tier -- these are directly-reported percentage changes in an
# established insulin-sensitivity index, not a back-solved engineering
# compromise, matching this file's own recently-clarified GOLD
# definition (see parameter_seed.py's own module docstring, Section 8
# housekeeping). Two real caveats, disclosed not hidden (mandatory
# code-reviewer finding: the first was already noted, the second added
# on review):
#   - PARALLEL-group design (separate patients per dose arm), not a
#     within-subject crossover -- more between-subject variance than
#     most other citations in this file, though the quoted percentages
#     are each arm's own within-arm baseline-vs-2-month change, not a
#     between-arm-only comparison.
#   - SINGLE-blind, not double-blind -- unlike every other blinding-
#     labeled citation in this file (indinavir, growth hormone,
#     glibenclamide, miglitol), this study's own outcome assessment was
#     only single-blind, a real (if partial) risk of measurement bias
#     that double-blind designs specifically control for. This doesn't
#     change the GOLD tier assignment (which tracks "directly reported,"
#     per this file's own codified definition, not blinding rigor), but
#     is worth disclosing alongside the parallel-group caveat above
#     rather than presenting the evidence base as more complete than it
#     is.
#
# secretion_scale explicitly NOT touched, DESPITE the same source
# reporting a real, dose-dependent positive finding there ("significantly
# increased fasting plasma insulin... 25%, 42%, 31%, and 45%,
# respectively") -- the same "compensatory hyperinsulinemic response to
# the measured resistance, not an independent secretory mechanism"
# reasoning already applied to hydrochlorothiazide above (KNOWN BUG #14
# avoidance). The glycated-hemoglobin finding (2%, 5%, 5%, 5%) is a
# longer-term integrated-glycemia readout this platform's own acute/
# postprandial ODE state doesn't represent -- not modelled here either.
#
# Zero double-counting risk confirmed: grepped `packages/` directly --
# no existing reference to atorvastatin/statin anywhere in this codebase
# prior to this chunk.
ATORVASTATIN_P3_DOSE_MG        = (0.0, 10.0, 20.0, 40.0, 80.0)
ATORVASTATIN_P3_FACTOR_AT_DOSE = (1.0, 1.0, 0.97, 0.97, 0.96)


def atorvastatin_p3_factor(dose_mg: float) -> float:
    """
    Multiplicative factor on PathwayParameters.p3 for a given daily
    atorvastatin dose (mg), piecewise-linearly interpolated across the
    five real Koh et al. 2010 (PMID:20298928) anchor points -- including
    the EMPIRICALLY NULL 10mg point (not significant vs. baseline,
    p=0.312), not an assumed monotonic curve from 0. Neutral (1.0) at
    dose_mg=0 through 10mg; ramps between 10mg and 80mg; clamped at the
    80mg value above 80mg (the highest tested dose in this source --
    "don't extrapolate past tested data", matching
    pioglitazone_p3_factor()'s own convention).
    """
    if not math.isfinite(dose_mg):
        raise ValueError(f"dose_mg must be finite, got {dose_mg}")
    if dose_mg < 0:
        raise ValueError(f"dose_mg cannot be negative, got {dose_mg}")
    doses = ATORVASTATIN_P3_DOSE_MG
    factors = ATORVASTATIN_P3_FACTOR_AT_DOSE
    if dose_mg <= doses[0]:
        return factors[0]
    if dose_mg >= doses[-1]:
        return factors[-1]
    for i in range(len(doses) - 1):
        if doses[i] <= dose_mg <= doses[i + 1]:
            t = (dose_mg - doses[i]) / (doses[i + 1] - doses[i])
            return factors[i] + (factors[i + 1] - factors[i]) * t
    raise AssertionError("unreachable")  # pragma: no cover


def atorvastatin_modifiers(dose: DrugDose) -> dict[str, float]:
    """
    Atorvastatin (HMG-CoA reductase inhibitor) -- see
    ATORVASTATIN_P3_DOSE_MG's own docstring above for the full citation
    and derivation. Only p3 is modified; secretion_scale is deliberately
    NOT touched despite a real, dose-dependent positive finding in the
    same source -- see that comment for why (compensatory response, not
    an independent mechanism -- would repeat KNOWN BUG #14).
    """
    if dose.dose_mg <= 0:
        return {}
    return {"p3": atorvastatin_p3_factor(dose.dose_mg)}


# ---------------------------------------------------------------------------
# Hydroxychloroquine (antimalarial/DMARD) -- Track 1 row 1.44. Found via a
# fresh Track 1 scoping fan-out for a genuinely new drug class (a
# repo-wide grep confirmed zero prior reference to hydroxychloroquine/
# chloroquine/antimalarial/DMARD/adiponectin anywhere in this codebase).
#
# Wasko MC, McClure CK, Kelsey SF, Huber K, Orchard T, Toledo FG.
# "Antidiabetogenic effects of hydroxychloroquine on insulin sensitivity
# and beta cell function: a randomised trial." Diabetologia.
# 2015;58(10):2336-2343. PMID:26197707 -- verified live via NCBI efetch
# this session. Randomized, double-blind, PARALLEL-arm (HCQ 400mg/day
# vs. placebo), n=17/15 non-diabetic overweight/obese adults with >=1
# insulin-resistance marker, 13 weeks. Insulin sensitivity measured by
# IVGTT + Bergman MINIMAL MODEL -- the SAME estimator this platform's
# own p3 is directly anchored to (Bergman 1981), a genuine metric match,
# not just a topically-relevant clamp study. Directly quoted: "There was
# a positive change in insulin sensitivity with HCQ but not placebo
# (mean+/-SEM: +20.0%+/-7.1% vs -18.4%+/-7.9%, respectively; p<0.01;
# difference: 38.3%+/-10.6%; 95% CI: 17%, 60%)."
#
# ANCHOR CHOICE, disclosed not hidden: uses the WITHIN-ARM active-
# treatment change (+20.0%), not the larger placebo-referenced
# difference (38.3%). The placebo arm's own -18.4% decline over just 13
# weeks is unusually large for a non-diabetic, overweight/obese cohort
# with no intervention -- and the authors themselves say so. Their own
# Discussion (full text, PMC4575248 -- an NIH author manuscript, freely
# retrievable via NCBI's PMC efetch; verified this session via
# biosim_literature.verify_citation()): "The exact effect size on insulin
# sensitivity was difficult to determine because the placebo group
# experienced an unexplained decline. However, if one considers only the
# change within the HCQ group (20 +/- 7.1%), the magnitude is still
# clinically significant." Using +20% therefore follows the authors' own
# reading, and is the conservative choice: it does not attribute an
# unexplained placebo-arm swing to the drug.
#
# (CORRECTION recorded: an earlier draft of this comment claimed full
# text could not be obtained and that Europe PMC's inEPMC="Y" flag was
# "misleading". Both were wrong -- Europe PMC's REST fullTextXML endpoint
# does fail for this non-open-access-subset manuscript, but NCBI's PMC
# serves it in full. The fetcher only tried Europe PMC; it now falls back
# to NCBI PMC efetch. Found by the chunk's code-reviewer.)
#
# Further, disclosed: the primary analysis was intention-to-treat with the
# baseline value carried forward for participants without week-13 data
# (pulling the within-arm change toward zero -- also conservative), and
# the authors' own sensitivity analysis (imputing the placebo-arm mean
# change for HCQ participants lacking week-13 data) still reports a
# significant effect. The +20.0% is a mean of per-subject percentage
# changes, not a ratio of group means. APPLICABILITY: the cohort is
# non-diabetic, overweight/obese with >=1 insulin-resistance marker; the
# authors state "the effects in type 2 diabetes must be inferred". This
# factor is applied to any profile, and a single HCQ dose event is treated
# as the chronic, steady-state 13-week regimen the study measured
# (same single-tested-dose simplification as GROWTH_HORMONE_P3_FACTOR).
#
# Corroborated independently (same direction, different lab technique,
# genuinely different-in-part author group -- Toledo FG appears on
# both, so NOT a fully independent replication, disclosed not hidden):
# Toledo FGS, Miller RG, Helbling NL, Zhang Y, DeLany JP. "The effects
# of hydroxychloroquine on insulin sensitivity, insulin clearance and
# inflammation in insulin-resistant adults: A randomized trial."
# Diabetes Obes Metab. 2021;23(6):1252-1261. PMID:33528893 -- also verified live.
# Same drug/dose/duration, hyperinsulinaemic-euglycaemic clamp:
# "Compared with placebo, HCQ significantly improved skeletal muscle
# insulin sensitivity by 26% (p = .019)... By contrast, HCQ had no
# effect on hepatic insulin sensitivity. HCQ did not affect insulin
# clearance..." -- NOT used as the anchor (Wasko's own minimal-model SI
# is the metric-matched source for this platform's own p3), but its
# TWO explicit nulls constrain scope for free: hepatic_egp_delta and ke
# are deliberately NOT touched by this mechanism, matching those
# directly-measured nulls exactly (same discipline as indinavir's own
# NEFA null).
#
# secretion_scale deliberately NOT modeled from the same dataset's own
# "+45.4% beta cell function" finding -- would repeat KNOWN BUG #14's
# over-determination mistake (a shared disposition-index confound
# between SI and beta-cell-function estimates from the same minimal-
# model fit, not two independent measurements).
#
# Single tested dose (400mg/day) only -- no dose-response curve exists
# in either source. Flat peak-effect factor, presence-gated (dose_mg >
# 0), not dose-scaled -- same design as GROWTH_HORMONE_P3_FACTOR above,
# for the identical reason (a single clinical dose, no ramp to fit).
# Further disclosures (HCQ security review): (1) The flat factor ignores
# hours_since_ingestion/is_chronic -- the source effect is a 13-WEEK
# steady-state result (HCQ's half-life is weeks), so ANY positive-dose
# event, even a single "just taken" dose, receives the full chronic effect;
# a scenario should only include HCQ to represent an established regimen.
# (2) Stacking with pioglitazone (1.20 x 1.6875 = 2.025 on p3) is an
# uncited, unvalidated extrapolation -- different mechanisms, no
# combination data; bounded by the p3-ceiling clamp-and-warn, not
# evidence-supported. (3) The abstract also reports "modest treatment
# effects on fasting plasma glucose and HbA1c (p<0.05)" -- Gb is
# deliberately NOT modelled (no extractable magnitude, and a fasting-
# glucose change is plausibly the downstream consequence of the SI/beta-cell
# changes already represented). (4) Hypoglycaemia risk (an HCQ labelling
# warning, notably alongside insulin/sulfonylureas) is not modelled: HCQ
# here only raises insulin sensitivity. (5) Exactly ONE active HCQ event is
# allowed per scenario (SINGLE_INSTANCE_DRUGS / validate_drug_combination,
# enforced at ScenarioDefinition construction): duplicate flat-factor
# events each multiplied the full 1.20 (20 duplicates reached the 1e-4 p3
# ceiling). (6) No cohort-mismatch warning exists (cf. cold_acclimatized's
# T2D-only warning): the factor applies to any profile although the cohort
# was non-diabetic overweight/obese -- logged as a follow-up.
HYDROXYCHLOROQUINE_P3_FACTOR = 1.0 + 0.20


def hydroxychloroquine_modifiers(dose: DrugDose) -> dict[str, float]:
    """
    Hydroxychloroquine -- see HYDROXYCHLOROQUINE_P3_FACTOR's own
    docstring above for the full citation, the placebo-drift anchor-
    choice reasoning, and why NO dose-response scaling is attempted (a
    flat peak-effect factor applied whenever dose_mg > 0, not a ramp).
    Only p3 is modified; hepatic_egp_delta/ke/secretion_scale
    deliberately NOT touched (see that comment for why).
    """
    if dose.dose_mg <= 0:
        return {}
    return {"p3": HYDROXYCHLOROQUINE_P3_FACTOR}
