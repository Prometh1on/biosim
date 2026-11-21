# SAVE TO: packages\physiology\src\biosim_physiology\drug_interactions\incretins.py
"""Incretin-class drugs: tirzepatide's static Gb factor and the GLP-1 class."""

from __future__ import annotations

import math

from ._core import Drug, DrugDose

# Furihata K, Mimura H, Urva S, Oura T, Ohwaki K, Imaoka T. "A phase 1
# multiple-ascending dose study of tirzepatide in Japanese participants
# with type 2 diabetes." Diabetes Obes Metab. 2022;24(2):239-246.
# PMID:34647404 (open access, PMC9299227; live-verified twice this
# session -- once for the week-8 vs.-placebo deltas, once specifically
# for the baseline values needed to convert them into a fractional
# reduction). n=48 T2D (9 placebo/11 5mg/12 10mg/16 15mg), tirzepatide
# once weekly for 8 weeks vs. placebo.
#
# Table 1 (baseline fasting serum glucose, mg/dL, mean): placebo 173.9,
# 5mg 164.5, 10mg 181.2, 15mg 170.7. Results 3.4 (LS-mean decrease vs.
# placebo at week 8, mg/dL): 52.7 (5mg), 69.1 (10mg), 68.9 (15mg), all
# P<.05 or P<.0001 vs. placebo. Converting each arm's own absolute
# reduction to a FRACTION of that arm's own baseline (a real, disclosed
# assumption -- fractional effects are more likely to transfer across a
# different reference population than an absolute mg/dL number would,
# and this platform's own default Gb is notably lower than this
# trial's own, poorly-controlled-T2D baseline):
#   5mg:  52.7/164.5 = 32.0% reduction (Gb multiplier 0.680)
#   10mg: 69.1/181.2 = 38.1% reduction (Gb multiplier 0.619)
#   15mg: 68.9/170.7 = 40.4% reduction (Gb multiplier 0.596)
# A real, genuinely saturating (not linear) dose-response -- the
# 10mg->15mg step gains far less than 5mg->10mg did. Represented via
# piecewise-linear interpolation across these THREE real anchor points
# (not an assumed parametric curve shape, which none of these 3 points
# alone can validate) -- see TIRZEPATIDE_GB_DOSE_MG/
# TIRZEPATIDE_GB_MULTIPLIER below and tirzepatide_gb_factor()'s own
# extrapolation policy (linear to neutral below 5mg, clamped at the
# 15mg value above 15mg -- real clinical dosing tops out at 15mg, so
# this is a documented "don't extrapolate past tested data" cap, not
# an oversight).
#
# Deliberately NOT composing an ODE-coupled GLP1Pathway for tirzepatide
# this pass (unlike semaglutide/liraglutide/dpp4_inhibitor) -- the
# original scoping (MASTER_PLAN row 1.5) decided tirzepatide needs BOTH
# a distinct fasting/basal Gb-lowering component (now real, via
# Furihata 2022 above) AND a reused, stronger GLP1Pathway incretin_
# effect for its separately-real ISR7.2/beta-cell-sensitivity
# enhancement -- but the SECOND piece's exact magnitude remains
# genuinely uncited (both papers describing it, Heise 2022/Mather 2024,
# are paywalled with no extractable numbers, confirmed earlier this
# session). Inventing an uncited "stronger than semaglutide" incretin
# magnitude would violate this project's own anti-fabrication
# discipline -- implementing only the piece with real data (the static
# Gb modifier) is a deliberate, disclosed partial implementation, not
# an oversight. `Drug.TIRZEPATIDE` is therefore NOT a member of
# GLP1_DRUGS (registry.py).
TIRZEPATIDE_GB_DOSE_MG        = (5.0, 10.0, 15.0)
TIRZEPATIDE_GB_MULTIPLIER     = (0.680, 0.619, 0.596)


def tirzepatide_gb_factor(dose_mg: float) -> float:
    """
    Multiplicative factor on PathwayParameters.Gb for a given weekly
    tirzepatide dose (mg), piecewise-linearly interpolated across the
    three real Furihata et al. 2022 (PMID:34647404) anchor points (see
    this module's own citation comment directly above). Neutral (1.0)
    at dose_mg=0; linearly ramps from neutral to the 5mg value between
    0 and 5mg (a documented, disclosed simplification -- no real data
    exists below 5mg); clamped at the 15mg value above 15mg (real
    clinical dosing does not exceed 15mg, so this is a deliberate
    "don't extrapolate past tested data" cap).
    """
    if not math.isfinite(dose_mg):
        raise ValueError(f"dose_mg must be finite, got {dose_mg}")
    if dose_mg < 0:
        raise ValueError(f"dose_mg cannot be negative, got {dose_mg}")
    if dose_mg <= 0.0:
        return 1.0
    doses = TIRZEPATIDE_GB_DOSE_MG
    mults = TIRZEPATIDE_GB_MULTIPLIER
    if dose_mg <= doses[0]:
        # Linear ramp from neutral (dose=0) to the first real anchor point.
        t = dose_mg / doses[0]
        return 1.0 + (mults[0] - 1.0) * t
    if dose_mg >= doses[-1]:
        return mults[-1]
    for i in range(len(doses) - 1):
        if doses[i] <= dose_mg <= doses[i + 1]:
            t = (dose_mg - doses[i]) / (doses[i + 1] - doses[i])
            return mults[i] + (mults[i + 1] - mults[i]) * t
    raise AssertionError("unreachable")  # pragma: no cover


def _glp1_agonist_modifiers(dose: DrugDose) -> dict[str, float]:
    """
    GLP-1 class drugs (semaglutide, liraglutide, DPP-4 inhibitors) intentionally
    return no PathwayParameters modifier here.

    Unlike caffeine/metformin/glucocorticoid/alcohol, their effect is not a
    static multiplicative tweak: it is modelled as a genuine ODE-coupled
    incretin pathway (biosim_ode.glp1_pathway.GLP1Pathway) composed with
    GlucoseInsulinPathway via glp1_pathway.couple_with_gi(). SimulationRunner
    detects these drugs and builds that composite instead of applying a
    modifier here -- doing both would double-model the same mechanism.

    Tirzepatide (dual GIP/GLP-1 agonist) is the one exception: its real,
    citation-anchored effect on this platform IS a static Gb modifier
    (see TIRZEPATIDE_GB_DOSE_MG's own docstring above for the full
    derivation and the disclosed reason its OTHER mechanism, a stronger
    incretin_effect, is not yet implemented).
    """
    if dose.drug == Drug.TIRZEPATIDE:
        if dose.dose_mg <= 0:
            return {}
        return {"Gb": tirzepatide_gb_factor(dose.dose_mg)}
    return {}
