# SAVE TO: packages\physiology\src\biosim_physiology\drug_interactions\sulfonylureas.py
"""Sulfonylureas: additive, glucose-independent insulin secretion."""

from __future__ import annotations

from ._core import Drug, DrugDose


def _sulfonylurea_modifiers(dose: DrugDose) -> dict[str, float]:
    """
    Sulfonylureas (glibenclamide, glipizide) intentionally return no PathwayParameters
    modifier here -- unlike every other drug in this file, their effect is
    genuinely ADDITIVE (a constant, glucose-independent secretion floor,
    PathwayParameters.glucose_independent_secretion, base value 0.0), not
    multiplicative -- and compose_all_modifiers()'s pipeline (a factor
    multiplying a base value) can never move a 0.0 base value away from
    zero regardless of the factor. Applied directly by
    PhysiologicalPersonalizer.personalize() instead (bypassing the generic
    net_factors loop, the same way the Gb/G_threshold safety margin and
    Brunzell's first_phase_gain zeroing already do) -- see
    glucose_independent_secretion_pM_per_min() below for the magnitude and
    citation. That application is gated on dose_mg > 0 (a zero-dose event
    must not apply an effect) but deliberately NOT gated on
    hours_since_ingestion/hours_before_simulation -- Riefflin et al. 2015
    (PMID:25297572) measured a CHRONIC (14-day), steady-state effect, so
    (matching metformin_modifiers()'s own established "always chronic,
    timing irrelevant" precedent -- see that function's docstring) this
    is a deliberate design choice, tested explicitly in test_physiology.py,
    not an unexamined gap.
    """
    return {}


# Riefflin A, Ayyagari U, Manley SE, Holman RR, Levy JC. "The effect of
# glibenclamide on insulin secretion at normal glucose concentrations."
# Diabetologia. 2015;58(1):43-9. PMID:25297572 (verified live via NCBI
# efetch this session): RCT, double-blind, placebo-controlled crossover,
# n=19 T2D, glibenclamide 2.5mg twice daily x14 days vs. placebo x14 days,
# 3-step hyperinsulinaemic glucose clamp (4.0/8.0/12.0 mmol/L). C-peptide
# secretion rate rose by an ABSOLUTE increment "independent of clamp
# glucose concentration (p=0.54)" -- 140/126/158 pmol/min at 4.0/8.0/12.0
# mmol/L respectively (proportional increase 2.8-/1.8-/1.7-fold vs.
# placebo, all p<0.001) -- i.e. genuinely glucose-INDEPENDENT, unlike this
# model's existing delta_secretion (Michaelis-Menten glucose-gated) and
# first_phase_secretion (glucose-derivative-gated) terms. This is the
# textbook mechanism behind sulfonylureas' real hypoglycaemia risk
# (KATP-channel closure via SUR1 binding triggers secretion regardless of
# ambient glucose) -- MASTER_PLAN.md's own original row 1.3 assumed a
# secretion_scale-only implementation ("reuses the existing
# secretion_scale_modifier rule-11 slot"), which would have been
# mechanistically BACKWARDS: secretion_scale multiplies delta_secretion,
# which floors at exactly -ke*Ib below G_threshold and goes MORE negative
# below basal glucose -- amplifying it would make a sulfonylurea appear to
# SUPPRESS insulin harder during hypoglycaemia, the opposite of its
# defining clinical risk (verified numerically before rejecting that
# design; see CLAUDE.md's sulfonylurea milestone).
#
# Riefflin's absolute pmol/min figures are whole-body C-peptide secretion
# RATES (via two-compartment deconvolution) -- a different physical
# quantity from this model's plasma INSULIN CONCENTRATION rate constant
# (glucose_independent_secretion, pM/min), and converting between them
# needs an insulin distribution-volume/clearance relationship this paper
# doesn't report and this codebase doesn't otherwise anchor. Rather than
# invent an unverified conversion factor, this uses the paper's ROBUST,
# GLUCOSE-INDEPENDENT FINDING (not its absolute units) to justify a
# bounded, disclosed engineering calibration against this model's OWN
# internal basal-secretion reference (ke*Ib) -- same "isolated, honestly
# labelled assumption" category as HPA_FEEDBACK_GAIN and the
# confidence-tier CV defaults elsewhere in this codebase -- rather than
# presenting a directly-cited number that doesn't actually exist:
#   0.5  ->  the analytic insulin fixed point at zero glucose (see
#            PathwayParameters.basal_state()) becomes 1.5x this patient's
#            own basal Ib -- clearly large enough to represent Riefflin's
#            real, clinically significant, glucose-independent effect
#            (1.7-2.8-fold at every tested glucose level, all p<0.001)
#            while staying well below the hard safety bound enforced by
#            PathwayParameters.__post_init__ (must be < ke*Ib, i.e. I_ss
#            must stay below 2x basal Ib) and not grossly contradicting
#            Go et al. 2004 (PMID:15133754)'s real-world, ambulatory
#            (non-clamped) finding that fasting insulin was "not
#            significantly increased" by glipizide/glibenclamide in their
#            own T2D cohort -- a genuine, disclosed tension between the
#            two studies, not resolved by silently picking whichever
#            number is more convenient.
SULFONYLUREA_SECRETION_FRACTION_OF_KE_IB = 0.5

# GLIPIZIDE_SECRETION_FRACTION_OF_KE_IB (2026-08-12 deep-search plan,
# item D1 -- glimepiride/glipizide extension): two real, live-verified,
# genuinely NEW-angle citations found (6 prior queries this session had
# all been rejected -- see the row's own history) --
#   Groop LC, Ratheiser K, Luzi L, Melander A, Simonson DC, Petrides A,
#   Bonadonna RC, Widen E, DeFronzo RA. "Effect of sulphonylurea on
#   glucose-stimulated insulin secretion in healthy and non-insulin
#   dependent diabetic subjects: a dose-response study." Acta
#   Diabetologica. 1991;28(2):162-8. PMID:1777653. n=12 healthy + 6
#   NIDDM, stepped hyperglycaemic clamp (1/3/7/17 mmol/L above basal),
#   oral 5mg glipizide: "Glipizide alone increased the plasma insulin
#   response by 180 pmol/l. A similar increase... was observed at each
#   hyperglycaemic step, indicating that glipizide did not affect the
#   sensitivity of the beta-cell to glucose" -- the same "constant
#   increment across multiple glucose levels" evidentiary shape as
#   Riefflin 2015's own glibenclamide finding, methodologically clean,
#   but for a HEALTHY (not T2D) cohort specifically for this number.
#   Groop L, Luzi L, Melander A, Groop PH, Ratheiser K, Simonson DC,
#   DeFronzo RA. "Different effects of glyburide and glipizide on
#   insulin secretion and hepatic glucose production in normal and
#   NIDDM subjects." Diabetes. 1987;36(11):1320-8. PMID:3117609. n=10
#   healthy + 6 NIDDM, IV-infused glyburide (glibenclamide) vs.
#   glipizide at MATCHED plasma drug concentrations: "the acute phase
#   (0-10 min) of plasma insulin and C-peptide increase in both control
#   and NIDDM subjects was twice as great with GZ [glipizide] compared
#   with GB [glyburide]" -- a real, T2D-INCLUSIVE, head-to-head relative
#   ratio (glipizide ~2x glibenclamide), directly usable to scale
#   glibenclamide's own already-adopted fraction.
#
# Doubling glibenclamide's own 0.5 fraction to exactly 1.0 would collide
# with PathwayParameters.__post_init__'s pre-existing hard safety bound
# (glucose_independent_secretion must be STRICTLY < ke*Ib -- see that
# check's own docstring) -- a fraction of exactly 1.0 would be rejected
# outright. Capped at 0.9 instead: as large as this model's own,
# already-established engineering safety ceiling allows while still
# representing a MEANINGFULLY larger effect than glibenclamide's own
# 0.5 (motivated by, but not an exact reproduction of, the real 2x
# ratio finding) -- same "disclosed, bounded engineering compromise"
# category as glibenclamide's own 0.5 itself (which was ALSO chosen to
# stay safely under this same ceiling, not a precise literature value).
GLIPIZIDE_SECRETION_FRACTION_OF_KE_IB = 0.9

# GLIMEPIRIDE_SECRETION_FRACTION_OF_KE_IB -- unlocks row 1.3b's own
# remaining open item (no clean absolute-magnitude citation had been
# found for glimepiride after two prior search passes). Found via a
# genuinely different query angle: a head-to-head RATIO to
# glibenclamide's own already-adopted fraction, rather than a fresh
# absolute number.
#
# Draeger KE, Wernicke-Panten K, Lomp HJ, Schuler E, Rosskamp R.
# "Long-term treatment of type 2 diabetic patients with the new oral
# antidiabetic agent glimepiride (Amaryl): a double-blind comparison
# with glibenclamide." Horm Metab Res. 1996;28(9):419-25. PMID:8911976
# -- verified live via NCBI efetch this session. International,
# prospective, double-blind trial, n=524 glimepiride (1-8mg once
# daily) vs. n=520 glibenclamide (2.5-20mg), matched to equivalent
# metabolic control (HbA1c 8.4% vs. 8.3%, fasting glucose 174 vs. 168
# mg/dl -- "not considered clinically relevant"). At matched glycaemia,
# directly quoted: "Statistically significant lower fasting insulin
# and C-peptide values were observed in glimepiride patients compared
# with glibenclamide (differences: insulin, -0.92 microU/ml [p = 0.04];
# C-peptide, -0.14 ng/ml [p = 0.03])." Corroborated directionally (no
# extractable magnitude) by Massi-Benedetti M, Herz M, Pfeiffer C.
# Horm Metab Res. 1996;28(9):451-5. PMID:8911982: "a significant
# suppression of endogenous insulin secretion was observed for
# glimepiride only" (vs. glibenclamide) during exercise -- same
# direction, glimepiride's secretion drive is weaker/more suppressible
# than glibenclamide's.
#
# Derivation, using only this platform's own reference Ib=66.0 pM/
# ke=0.093 and the ALREADY-ESTABLISHED pM-per-uU/mL conversion factor
# (6.945, the same constant this codebase already uses for p3's own
# Bergman-1981 derivation -- see pathway.py's own p3 docstring; NOT a
# fresh 6.0 approximation):
#   glibenclamide's own fixed-point elevation = 0.5 * Ib = 33.0 pM
#                                              = 33.0/6.945 = 4.752 uU/mL
#   glimepiride's own elevation = 4.752 - 0.92 (Draeger's directly-
#                                 quoted difference) = 3.832 uU/mL
#                               = 3.832 * 6.945 = 26.61 pM
#   glimepiride fraction = 26.61 / 66.0 = 0.403 -> rounded to 0.40
#
# A disclosed engineering calibration (same category as glibenclamide's
# own 0.5 and glipizide's own 0.9) -- anchored to Draeger's real,
# directly-quoted ratio, but the conversion to a dimensionless fraction
# depends on this platform's own reference Ib, not a literature value
# in its own right. Well below the hard safety bound (must be < 1.0,
# i.e. glucose_independent_secretion < ke*Ib) -- same margin category
# as glibenclamide's 0.5, nowhere near glipizide's already-tight 0.9.
#
# Limitation, disclosed not hidden: Draeger reports CHRONIC fasting
# insulin under stable dosing, a weaker evidentiary shape than
# Riefflin's/Groop's own constant-across-glucose-steps clamp design --
# and the constant-fraction design means this ratio is anchored at
# REFERENCE Ib specifically; a T1D-preset's own much-lower Ib would
# scale the absolute pM difference down proportionally, not
# necessarily matching Draeger's own absolute microU/mL gap at that
# very different insulin level.
GLIMEPIRIDE_SECRETION_FRACTION_OF_KE_IB = 0.40

_SULFONYLUREA_SECRETION_FRACTION: dict = {
    Drug.GLIBENCLAMIDE: SULFONYLUREA_SECRETION_FRACTION_OF_KE_IB,
    Drug.GLIPIZIDE:     GLIPIZIDE_SECRETION_FRACTION_OF_KE_IB,
    Drug.GLIMEPIRIDE:   GLIMEPIRIDE_SECRETION_FRACTION_OF_KE_IB,
}


def glucose_independent_secretion_pM_per_min(ke: float, Ib: float, drug: "Drug") -> float:
    """
    The constant, glucose-independent addition to insulin secretion for an
    active sulfonylurea dose (see SULFONYLUREA_SECRETION_FRACTION_OF_KE_IB's,
    GLIPIZIDE_SECRETION_FRACTION_OF_KE_IB's, and GLIMEPIRIDE_SECRETION_
    FRACTION_OF_KE_IB's own docstrings for the citations and calibration
    reasoning, per drug). Not dose-scaled for any of the three drugs --
    Riefflin et al. 2015 tested a single clinical dose (2.5mg twice
    daily glibenclamide), Groop et al. 1991's glipizide number is
    likewise a single-dose (5mg) finding, and Draeger et al. 1996's
    glimepiride number is a chronic stable-dosing comparison (1mg vs.
    2.5mg glibenclamide) with no within-study dose-response curve
    either; no dose-response curve was found to scale against for any
    of the three (a real, disclosed limitation, not an oversight -- see
    MASTER_PLAN.md row 1.3b).
    """
    return _SULFONYLUREA_SECRETION_FRACTION[Drug(drug)] * ke * Ib
