# SAVE TO: packages\physiology\src\biosim_physiology\drug_interactions\alpha_glucosidase.py
"""Alpha-glucosidase inhibitors: meal absorption (kabs) scaling."""

from __future__ import annotations

import math

from ._core import DrugDose


def _alpha_glucosidase_inhibitor_modifiers(dose: DrugDose) -> dict[str, float]:
    """
    Alpha-glucosidase inhibitors (acarbose, miglitol) intentionally
    return no PathwayParameters modifier here -- their real effect is on
    MealInput.kabs (the meal absorption rate constant), applied directly
    by SimulationRunner at meal-construction time, not on any
    PathwayParameters field at all. See acarbose_kabs_scale()/
    miglitol_kabs_scale() below for the per-drug magnitude and citation.
    """
    return {}


# Chen Z, Fu X, Kuang J, Chen J, Chen H, Pei J, Yang H. "Single-dose
# acarbose decreased glucose-dependent insulinotropic peptide and
# glucagon levels in Chinese patients with newly diagnosed type 2
# diabetes mellitus after a mixed meal." BMC Endocr Disord. 2016.
# PMID:27686734 (verified live via NCBI efetch this session): n=12 newly
# diagnosed T2D, single 100mg acarbose dose with a mixed meal (370kcal,
# 74% carbohydrate), same meal given a week apart without and with the
# drug (a real, if unblinded/non-randomized-order, same-meal comparison
# -- disclosed limitation, not a randomized double-blind crossover like
# this project's other drug benchmarks). Table 2's own numbers: peak
# glucose 13.22 mmol/L at 60min (no drug) -> 9.45 mmol/L at 90min (with
# drug) -- a **-28.5% peak-height reduction combined with a 60->90min
# (1.5x) peak-TIME DELAY**, not just a lower peak. This dual signature
# (lower AND later peak) is exactly what a genuine absorption RATE-
# CONSTANT slowdown produces (as opposed to e.g. a pure amplitude
# scale on Ra(t), which would only lower the peak without moving it) --
# confirming a direct kabs reduction is the mechanistically correct
# design, not merely a convenient one (verified by reproducing both
# numbers simultaneously in this platform's own ODE model before
# picking a magnitude -- see ACARBOSE_100MG_KABS_SCALE below).
# Insulin moved the same direction as glucose (-40.6% peak, no time
# shift) -- consistent with insulin secretion simply tracking a smaller/
# later glucose stimulus, not a separate anti-secretory mechanism (also
# corroborated by Zhou et al. 2013, PMID:23631607, a head-to-head study
# where a true secretagogue's insulin ROSE while acarbose's fell) --
# supporting a pure-kabs-modifier design with no additional
# PathwayParameters effect.
#
# ACARBOSE_100MG_KABS_SCALE = 0.35: the kabs multiplier that, run
# through this platform's own ODE model for a 75g meal, simultaneously
# reproduces Chen et al.'s two independent observations closely
# (-27.0% peak height vs. the real -28.5%; 1.51x peak-time delay vs.
# the real 1.5x) -- a genuine two-number consistency check, not a
# single-target curve fit picked to match only one of them.
ACARBOSE_100MG_KABS_SCALE = 0.35

# No study was found giving a direct acarbose-dose-to-kabs-percentage
# relationship (a real, disclosed gap -- see MASTER_PLAN.md row 1.4b).
# Fischer et al. 1998 (PMID:9625287, a double-blind placebo-controlled
# 495-patient dose-ranging RCT) gives a real dose-response CURVE, but in
# chronic (24-week) postprandial-AUC-percentage terms, not an acute
# single-dose kabs magnitude -- a different metric this platform cannot
# directly translate into a kabs value. That paper's own reported SHAPE
# (a saturating/plateau relationship: 25mg -11.6%/-11.3% at 1h/2h,
# 50-100mg plateau ~-13%, 200mg -18.5%) motivates using a log-dose
# response capped at the one directly-measured reference dose (100mg,
# Chen et al.'s own dose) rather than extrapolating linearly to higher
# doses this platform has no acute-kabs data for -- same "log-dose,
# calibrated at one measured point, not independently validated per-
# dose" limitation already disclosed for caffeine_modifiers().
def acarbose_kabs_scale(dose_mg: float) -> float:
    """
    Multiplicative factor on MealInput.kabs for an active acarbose dose
    (see ACARBOSE_100MG_KABS_SCALE's docstring for the citation and
    calibration). Returns 1.0 (no effect) for dose_mg <= 0. Capped at the
    100mg reference dose -- see this function's own module-level
    docstring for why higher doses aren't extrapolated further.
    """
    # Found by the mandatory security-reviewer pass: this function takes
    # a raw float directly (unlike every other modifier function in this
    # file, which takes a DrugDose/DrugEvent already validated finite by
    # its own __post_init__ -- see DrugEvent's own comment about exactly
    # this NaN/inf class of bug), so a NaN here isn't caught upstream by
    # anything. Without this guard, `min(1.0, float('nan'))` silently
    # returns 1.0 (Python's min() takes the first argument on a NaN
    # comparison), which would make dose_factor=1.0 and silently treat an
    # invalid dose as the maximum calibrated one, rather than raising.
    # Not reachable via the real DrugEvent/SimulationRunner path today
    # (DrugEvent.__post_init__ already rejects non-finite dose_mg before
    # this function is ever called), but this function is public and
    # independently callable with a raw float, so it needs its own
    # guard rather than relying on that upstream validation alone.
    if not math.isfinite(dose_mg):
        raise ValueError(f"dose_mg must be finite, got {dose_mg}")
    if dose_mg <= 0:
        return 1.0
    max_reduction = 1.0 - ACARBOSE_100MG_KABS_SCALE
    dose_factor = min(1.0, math.log(1.0 + dose_mg / 100.0) / math.log(2.0))
    return 1.0 - max_reduction * dose_factor


# Samad AH, Ty Willing TS, Alberti KG, Taylor R. "Effects of BAYm 1099,
# new alpha-glucosidase inhibitor, on acute metabolic responses and
# metabolic control in NIDDM over 1 mo." Diabetes Care. 1988;11(2):
# 100-9. PMID:3042310 (verified live via NCBI efetch this session --
# BAYm 1099 confirmed to be miglitol's own original development code,
# cross-checked against PMID:9540024, "the alpha-glucosidase inhibitor
# miglitol (BAY m 1099)"). n=12 diet-treated NIDDM (T2D) patients
# completing a double-blind, placebo-controlled crossover: 4 weeks
# miglitol 50mg TID vs. 4 weeks placebo, with test meals given during
# the last week of each phase. Test meal 2 (TM2, 25g sucrose alone --
# the purest, most demanding substrate for an alpha-glucosidase
# inhibitor, since sucrose is a single, rapidly-cleaved disaccharide
# with no starch-digestion component at all): peak blood glucose
# **11.3 +/- 1.0 mM (drug) vs. 14.3 +/- 1.4 mM (placebo), P<0.001 -- a
# -21.0% peak-height reduction**. "Peak blood glucose levels were
# significantly delayed" is stated qualitatively but with no reported
# ratio -- unlike acarbose's own Chen et al. 2016 calibration (which
# had both a peak-height AND a peak-time number to cross-check), this
# is a SINGLE-quantity anchor. `MIGLITOL_50MG_KABS_SCALE` below is
# calibrated to reproduce ONLY the peak-height number exactly, through
# the same mass-conserving kabs-reduction mechanism already verified
# mechanistically correct for this drug class (acarbose_kabs_scale()'s
# own docstring); the resulting ~2.05x peak-time delay this produces is
# a mechanical CONSEQUENCE of that mechanism, not independently verified
# against a real timing figure the way acarbose's 1.51x-vs-1.5x match
# was -- disclosed explicitly, not presented as a confirmed second data
# point.
MIGLITOL_50MG_KABS_SCALE = 0.1446

# No miglitol dose-ranging study reporting an ACUTE, single-dose,
# peak-height-comparable magnitude was found this session (a real,
# disclosed gap, same category as acarbose's own missing dose-to-kabs
# curve) -- PMID:9540024 (Hispanic NIDDM, 50/100/150/200mg TID dose-
# ranging) reports only 1-year HbA1c/fasting-glucose endpoints, a
# chronic, different metric this platform cannot translate into an
# acute kabs magnitude either. Rather than extrapolate acarbose's own
# log-dose SHAPE onto a mechanistically-related but pharmacologically
# distinct drug (miglitol is systemically absorbed and has a different
# relative potency across glucosidase subtypes than acarbose, which is
# not absorbed at all) with no miglitol-specific curve to justify that
# choice, this function is capped at the single verified reference dose
# (50mg) with no extrapolation beyond it -- a narrower, more
# conservative scope than acarbose_kabs_scale()'s own log-dose curve.
def miglitol_kabs_scale(dose_mg: float) -> float:
    """
    Multiplicative factor on MealInput.kabs for an active miglitol dose
    (see MIGLITOL_50MG_KABS_SCALE's docstring for the citation and
    calibration). Returns 1.0 (no effect) for dose_mg <= 0. Any dose >= 50mg
    (the one verified reference dose) gets the full calibrated effect --
    no dose-response curve was found to extrapolate beyond or below it,
    so this is deliberately a step function, not acarbose_kabs_scale()'s
    continuous log-dose curve.
    """
    if not math.isfinite(dose_mg):
        raise ValueError(f"dose_mg must be finite, got {dose_mg}")
    if dose_mg <= 0:
        return 1.0
    if dose_mg < 50.0:
        # No sub-50mg data either -- linear interpolation from (0, no
        # effect) to (50mg, the verified reduction) is the least
        # speculative option available (matches no effect at zero dose
        # by construction), disclosed as an interpolation, not a
        # separately-measured point.
        return 1.0 - (1.0 - MIGLITOL_50MG_KABS_SCALE) * (dose_mg / 50.0)
    return MIGLITOL_50MG_KABS_SCALE
