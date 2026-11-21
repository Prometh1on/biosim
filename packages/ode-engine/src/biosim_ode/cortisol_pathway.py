# SAVE TO: packages\ode-engine\src\biosim_ode\cortisol_pathway.py
"""
Cortisol axis pathway extension.

Cortisol is the primary glucocorticoid counter-regulatory hormone.
It opposes insulin action and elevates blood glucose, making it
essential for modelling stress responses, Cushing's syndrome,
steroid medication effects, and circadian metabolic variation.

State vector (internal units):
    y[0]  Cortisol  [nM]   plasma cortisol concentration
    y[1]  CRH      [pM]   hypothalamic corticotropin-releasing hormone
    y[2]  GR_occ   [0-1]  glucocorticoid receptor fractional occupancy

Literature:
    Rizza et al. 1982, J Clin Endocrinol Metab. PMID:7033265
    Scheen/Van Cauter et al. 1996, Am J Physiol. PMID:8770019
    Dinneen et al. 1993, J Clin Invest. PMID:8227343
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Optional

from scipy.optimize import brentq

from .base_pathway import BasePathway, PathwayDescriptor, PathwayMaturity
from .kinetics import hill_activation, approach_to_target

if TYPE_CHECKING:
    from .coupling import CompositePathway
    from .pathway import GlucoseInsulinPathway

CORTISOL_BASAL_NM    = 0.28
CRH_BASAL_PM         = 5.0
CORTISOL_PEAK_HOUR   = 7.0
CORTISOL_AMPLITUDE   = 0.18
CORTISOL_HALF_LIFE   = 80.0
EC50_GR_NM           = 0.15
HILL_GR              = 1.5
GR_RELAXATION        = 0.05
CRH_HALF_LIFE        = 10.0
CORTISOL_GB_SENSITIVITY  = 0.08

# GR_OCC_BASAL must be the ANALYTIC fixed point of the GR relaxation ODE
# (dGR = GR_RELAXATION * (hill_activation(cortisol, EC50_GR_NM, HILL_GR) - GR_occ)),
# not an independently-chosen constant -- KNOWN BUG (see CLAUDE.md): this
# was previously hardcoded to 0.15 (coincidentally equal to EC50_GR_NM),
# which is NOT hill_activation(CORTISOL_BASAL_NM, EC50_GR_NM, HILL_GR)
# (~0.718). That mismatch meant y0 was not actually a fixed point: every
# simulation using this pathway, even at zero stress/steroid, silently
# drifted GR occupancy from 0.15 toward ~0.72 over about an hour, cutting
# p3 (insulin action) by ~14% with no stress or steroid present at all.
# Computing it here instead of hardcoding a duplicate literal guarantees it
# can never drift out of sync with EC50_GR_NM/HILL_GR/CORTISOL_BASAL_NM again.
GR_OCC_BASAL = hill_activation(CORTISOL_BASAL_NM, EC50_GR_NM, HILL_GR)

# Fixed population-average basal GR occupancy used as the p3_modifier
# reference point (same fixed-reference pattern as NEFA_IR_REFERENCE) --
# must equal GR_OCC_BASAL exactly so p3_modifier is neutral (1.0) at the
# true basal steady state, not at an arbitrary constant.
GR_REFERENCE = GR_OCC_BASAL

# CORTISOL_P3_SENSITIVITY -- derived, not hardcoded (Section 8 housekeeping
# fix, found during Track 1 row 1.38/Cushing's syndrome scoping). The
# previous value (0.25) was UNCITED, and had a structural consequence no
# one had checked for: because `excess_gr = gr_occ - GR_REFERENCE` is
# bounded above by `1.0 - GR_REFERENCE` (~0.2817, since gr_occ is a Hill-
# saturation output that asymptotically approaches but never reaches 1.0),
# `p3_modifier = 1.0 - CORTISOL_P3_SENSITIVITY * excess_gr` could NEVER
# fall below ~0.9296 (a ~7% max reduction) no matter how much cortisol was
# present -- independently confirmed via direct execution (a gr_occ sweep,
# and a real SimulationRunner glucocorticoid-dose sweep from 10mg to
# 500mg moving glucose peak only 10.5044->10.5058 mM). This directly
# contradicted this module's own sibling function's docstring claim
# (glucocorticoid_modifiers() in physiology/drug_interactions.py: "Acute
# effect: reduces p3 ~30%") and that function's own note that
# CortisolPathway is the SOLE path SimulationRunner uses for glucocorticoid
# p3 effects -- meaning the documented ~30% figure never actually reached
# simulated output via this pathway.
#
# Re-derived from a real primary source ALREADY in this module's own
# "Literature:" list above, not a new citation: **Rizza RA, Mandarino LJ,
# Gerich JE. "Cortisol-induced insulin resistance in man: impaired
# suppression of glucose production and stimulation of glucose utilization
# due to a postreceptor detect of insulin action." J Clin Endocrinol Metab.
# 1982;54(1):131-8. PMID:7033265** -- a 24h IV cortisol infusion (2
# ug/kg/min) in 6 healthy volunteers vs. saline, directly measuring
# before/after insulin dose-response curves -- the ideal design for this
# purpose (a known, controlled cortisol elevation, not an endogenous
# disease state needing a within-study ratio guess). Confirmed verbatim:
# "increased the plasma cortisol concentration approximately 4-fold (37+/-3
# vs. 14+/-1 ug/dL)" -- the abstract's OWN "~4-fold" phrasing does not
# match its own quoted numbers (37/14=2.643, not 4); the ACTUAL computed
# ratio is used below, not the abstract's rounded claim. The paper reports
# TWO separate ED50 (insulin sensitivity, not Vmax/responsiveness) shifts:
# "suppression of glucose production (half-maximal response at 81+/-19 vs.
# 31+/-5 uU/mL)" (HEPATIC) and "stimulation of glucose utilization (half-
# maximal response at 104+/-9 vs. 64+/-7 uU/mL)" (PERIPHERAL), "with
# preservation of normal maximal responses to insulin" (confirms this is
# genuinely a sensitivity/ED50 shift, not a Vmax/responsiveness one --
# the right SHAPE of quantity for Bergman p3, unlike Nosadini et al. 1983
# (JCEM, PMID:6348064, the original candidate anchor for row 1.38/
# Cushing's syndrome), whose own abstract explicitly attributes ITS 0.513
# ratio to a responsiveness/Vmax defect instead ("without any significant
# difference in the levels of insulin half-maximally effective") -- a
# genuine, unresolved literature conflict on mechanism SHAPE between the
# two papers, disclosed here rather than silently ignored; Rizza's number
# is used specifically because it IS a sensitivity-type measure matching
# what p3 represents, not because Nosadini's finding is wrong.
#
# Only the PERIPHERAL ED50 shift is used, not the hepatic one -- primarily
# a SHAPE-matching argument: Bergman p3 represents peripheral glucose
# disposal specifically, and the hepatic ED50 shift is simply the wrong
# physiological quantity for that slot, independent of anything else in
# this codebase. A secondary, more hypothetical reason reinforces the same
# choice: HepGlucPathway/EGPParameters (hepatic.py) independently models
# hepatic glucose production, and there is today no live coupling from
# CortisolPathway into it (confirmed by reading hepatic.py directly --
# HepGlucPathway.dydt() is driven only by signals["insulin_pM"], and
# cortisol_couplings() wires only gb_modifier/p3_modifier to GI, nothing
# to hepatic_egp_delta) -- so using the hepatic figure here would not be a
# literal code-level double-count today, only a redundant/inconsistent
# representation of the same real-world finding IF a future scenario ever
# composed both pathways for the same disease/stress state. Matches this
# project's established mechanism-matching discipline (e.g. NEFA's
# antilipolysis vs. hepatic glucose production kept separate for the same
# shape-matching reason).
#
# Solved algebraically, not hand-tuned: at Rizza's own measured cortisol
# ratio, the peripheral ED50 shift implies p3_modifier should equal
# 64/104=0.6154 -- CORTISOL_P3_SENSITIVITY is exactly the value that makes
# this formula produce that target at that ratio, computed here (not a
# duplicate hardcoded literal) so it can never drift out of sync with
# GR_REFERENCE/EC50_GR_NM/HILL_GR/CORTISOL_BASAL_NM, same discipline as
# GR_OCC_BASAL above. New theoretical floor as cortisol->infinity: ~0.453
# (a ~55% max reduction) -- still above the function's own hardcoded 0.3
# safety floor (unchanged, remains an unreached backstop, not this fix's
# concern).
#
# Known, disclosed limitation NOT addressed by this fix: the GR-occupancy
# Hill curve (EC50_GR_NM=0.15, far smaller than CORTISOL_BASAL_NM=0.28)
# saturates very quickly for any real exogenous glucocorticoid addition,
# so the corrected p3_modifier still varies only modestly across the
# clinical dose range (e.g. ~48% reduction at 10mg prednisone-equivalent
# vs. ~55% at 120mg) -- this fix corrects the MAGNITUDE (the unreachable-
# floor bug), not the fast SATURATION SPEED, which is a separate,
# pre-existing property of EC50_GR_NM's own scale, out of this fix's scope.
#
# Also disclosed, not reconciled here: glucocorticoid_modifiers()'s own
# static p3 value (physiology/drug_interactions.py, capped at a 30%
# reduction at its own 30mg reference dose) now differs in magnitude from
# what this corrected pathway delivers (~53% at the same dose) -- that
# function is explicitly a dead path for SimulationRunner's own real
# glucocorticoid handling (see its own docstring), so this is a disclosed
# inconsistency between a live and a dead code path, not a double-counting
# risk; worth reconciling in a future pass, not fixed here (out of this
# fix's own scope).
_RIZZA_1982_CORTISOL_RATIO = 37.0 / 14.0            # 2.643 -- the ACTUAL ratio, not the abstract's own "~4-fold"
_RIZZA_1982_TARGET_P3_MODIFIER = 64.0 / 104.0       # 0.6154 -- peripheral (glucose utilization) ED50 ratio
_rizza_1982_cortisol_nM = CORTISOL_BASAL_NM * _RIZZA_1982_CORTISOL_RATIO
_rizza_1982_excess_gr = hill_activation(_rizza_1982_cortisol_nM, EC50_GR_NM, HILL_GR) - GR_REFERENCE
CORTISOL_P3_SENSITIVITY = (1.0 - _RIZZA_1982_TARGET_P3_MODIFIER) / _rizza_1982_excess_gr

# Exogenous corticosteroid (prednisone etc.) -> cortisol-equivalent GR ligand.
# Standard clinical glucocorticoid dose equivalency: 5mg prednisone ~ 20mg
# hydrocortisone (cortisol), i.e. 1mg prednisone ~ 4mg cortisol-equivalent
# (Meikle AW, Tyler FH. "Potency and duration of action of glucocorticoids.
# Effects of hydrocortisone, prednisone and dexamethasone on human
# pituitary-adrenal function." Am J Med. 1977. PMID:888843 -- standard
# endocrinology potency-ratio reference). Physiological cortisol production
# is ~10mg/day hydrocortisone-equivalent in unstressed adults (Esteban NV,
# Loughlin T, Yergey AL, et al. "Daily cortisol production rate in man
# determined by stable isotope dilution/mass spectrometry." J Clin
# Endocrinol Metab. 1991;72(1):39-45. PMID:1986026 -- confirmed
# verbatim, live-fetched: "FPR of 27.3 +/- 7.5 mumol/day (9.9 +/- 2.7
# mg/day)" in 12 normal volunteers, matching this constant almost
# exactly). CORRECTION: this block previously cited PMID:412539 (an
# unrelated 1977 paper on astrocyte histoenzymology in experimental
# epilepsy) and PMID:1955512 (a real paper -- DeFronzo/Barzilai/Simonson
# 1991 JCEM on metformin mechanism, already correctly cited elsewhere in
# this codebase for a DIFFERENT fact -- but not this one) for these two
# facts. Found and fixed during a Track 5 curation-campaign pass that
# went to independently verify these citations before seeding them into
# ParameterStore; see CLAUDE.md's "DONE — Fixed two mistranscribed
# citations in cortisol_pathway.py" entry. Combining these gives a dose ->
# cortisol-equivalent-nM conversion without needing a distribution-volume
# assumption: a dose that matches one full day's endogenous production
# contributes one cortisol_basal-worth of GR ligand.
PREDNISONE_TO_CORTISOL_MG_RATIO = 4.0    # mg cortisol-equivalent per mg prednisone
DAILY_CORTISOL_PRODUCTION_MG    = 10.0   # mg/day hydrocortisone-equivalent, unstressed adult

# Same "30 mg prednisone equivalent" saturation reference glucocorticoid_
# modifiers() (physiology/drug_interactions.py) already uses to cap its own
# dose_factor -- reused here (via the same dose->nM conversion above) rather
# than an independent constant, so the two representations' saturation
# points stay consistent. See CORTISOL_GB_SATURATION_NM's use in
# output_signals() below.
_GB_SATURATION_REFERENCE_DOSE_MG = 30.0
CORTISOL_GB_SATURATION_NM = (
    (_GB_SATURATION_REFERENCE_DOSE_MG * PREDNISONE_TO_CORTISOL_MG_RATIO / DAILY_CORTISOL_PRODUCTION_MG)
    * CORTISOL_BASAL_NM
)

# ---------------------------------------------------------------------------
# Glycaemic-effect PD timing (onset/peak/resolution) -- distinct from plasma
# PK, and lagging it by hours.
#
# Real prednisone plasma PK (Cai et al. 2025, BMC Pharmacol Toxicol,
# PMID:39962604, n=48 healthy Chinese subjects, single 5mg dose, fasted --
# independently re-verified against this paper's own Table 3, full text,
# not just its abstract): the active prednisolone moiety's Tmax = 1.84h,
# t1/2 = 2.32h. But the GLYCAEMIC effect does not track plasma
# concentration -- glucocorticoids act genomically (GR binding -> gene
# transcription -> protein synthesis -> metabolic effect), a mechanism
# that lags plasma Tmax by hours, not minutes. Two independent human
# clinical studies confirm the real effect's time course (both re-verified
# against live full text this session):
#   Burt et al. 2011 (JCEM, PMID:21411550, n=60 COPD, continuous glucose
#   monitoring, prednisolone ~08:00 dose): glucose "not significantly
#   different" from untreated 0000-1200h, significantly elevated only
#   1200-2400h ("Prednisolone predominantly causes hyperglycemia in the
#   afternoon and evening").
#   Yuen et al. 2012 (Clin Endocrinol, PMID:21973241, n=7 T2DM + 8
#   at-risk + 5 NGT, prednisone 20mg each morning x3 days, 24h frequent
#   sampling): glucose elevation concentrated midday (P=0.001) to
#   midnight (P=0.02) in the T2DM group, "dissipates overnight"; fasting
#   glucose/insulin/C-peptide UNCHANGED in ALL groups.
# The curve below is a CALIBRATION against these two studies' observed
# onset/peak/resolution WINDOWS, not a measured PK or PD curve -- no
# quantitative human glycaemic-effect time-course dataset at finer
# resolution was found, and no rat-model corticosteroid PK/PD parameter
# set (the only quantitative ones located during scoping) is transferable
# to this human platform -- checked and rejected, see CLAUDE.md's
# glucocorticoid PK/PD milestone entry.
GC_PD_ONSET_T50_HOURS           = 4.0   # hours to half-maximal onset
GC_PD_ONSET_HILL_N              = 4.0   # onset steepness (dimensionless)
GC_PD_RESOLUTION_HALFLIFE_HOURS = 5.0   # post-peak decay half-life (hours)
# Peak of the raw onset*decay expression in glucocorticoid_pd_effect_
# fraction() at the constants above, found by numerical maximization
# (occurs at h ~= 5.68h) -- purely a normalization constant so the curve
# peaks at exactly 1.0, not itself a literature value.
GC_PD_PEAK_RAW = 0.36520


def glucocorticoid_pd_effect_fraction(hours_since_dose: float) -> float:
    """
    Fraction (0-1) of the full glycaemic effect active at a given time
    since a glucocorticoid dose -- see the module-level constants above
    for the full citation trail and calibration rationale.

    Shape: ~0 in the first ~2h, material by ~4h (end of Burt's flat
    0000-1200h window for an 08:00 dose), peaks ~5.7h (just under the
    clinically observed 6-10h peak window -- close enough to be
    considered within it given this is a calibration, not a precision
    fit), falls to <15% of peak by hour 20-24 ("dissipates overnight").
    Returns 0.0 for hours_since_dose <= 0 (not yet dosed).

    Chronic/repeated daily dosing: hours_since_dose is wrapped into a
    24h cycle (a no-op for any hours_since_dose < 24, i.e. every acute-
    dose scenario is completely unaffected) before evaluating the curve.
    Regression (found by code review): without this, a chronically-dosed
    scenario (DrugEvent.hours_before_simulation spanning multiple days --
    this platform's existing convention for "steady-state, ongoing
    therapy", e.g. metformin's `is_chronic`/large-hours-elapsed pattern)
    would evaluate this SINGLE-DOSE decay curve at its raw, cumulative
    elapsed time and silently fade to ~0 forever after the first day --
    contradicting Yuen et al. 2012 (PMID:21973241) itself, which dosed
    prednisone 20mg EACH morning for 3 days and observed the SAME
    onset/peak/"dissipates overnight" pattern recur daily, not a
    cumulative decay. Wrapping models the real repeating daily cycle
    instead.
    """
    if hours_since_dose <= 0.0:
        return 0.0
    hours_in_cycle = hours_since_dose % 24.0
    onset = hours_in_cycle ** GC_PD_ONSET_HILL_N / (
        hours_in_cycle ** GC_PD_ONSET_HILL_N + GC_PD_ONSET_T50_HOURS ** GC_PD_ONSET_HILL_N
    )
    decay = math.exp(-math.log(2) / GC_PD_RESOLUTION_HALFLIFE_HOURS * hours_in_cycle)
    return min(1.0, (onset * decay) / GC_PD_PEAK_RAW)


# Reference dose/cap for glucocorticoid_modifiers()'s (physiology/
# drug_interactions.py) own secretion_scale dose-response -- duplicated
# here (both are one-line formulas) rather than imported, to avoid adding
# a new ode-engine -> physiology package dependency; biosim_api's
# runner.py already imports both packages and is the correct place to
# compose them.
GC_SECRETION_SCALE_REFERENCE_DOSE_MG = 10.0
GC_SECRETION_SCALE_DOSE_CAP_MG       = 30.0
GC_SECRETION_SCALE_MAGNITUDE         = 0.10


def glucocorticoid_secretion_scale_factor(
    dose_mg: float,
    hours_since_ingestion: float,
    is_t2d: bool,
) -> float:
    """
    Multiplier on secretion_scale for a glucocorticoid dose, gated by the
    same PD-timing curve as the Gb/p3 effect (glucocorticoid_pd_effect_
    fraction() above) and sign-flipped for T2D.

    Direction: glucocorticoid_modifiers() always applies a +10%-at-
    reference-dose COMPENSATORY increase -- real physiology in a patient
    who can still mount a compensatory beta-cell response to
    steroid-induced insulin resistance. But Yuen et al. 2012
    (PMID:21973241) found the OPPOSITE direction specifically in T2DM:
    insulin (P=0.03) and C-peptide (P=0.04) DECREASED postbreakfast in
    the T2DM group (reduced beta-cell reserve cannot compensate), while
    the 'at risk'/pre-diabetic group showed a significant COMPENSATORY
    INCREASE in C-peptide from midday to midnight (P=0.03) -- the
    direction this platform's existing +10% already represents. Shah
    et al. 2020 (PLoS One, PMID:32369480, n=5 T2DM + 5 healthy, prednisone
    40mg/day x3d, re-verified against full text this session)
    independently confirms the same T2D-specific direction: insulinogenic
    index fell 69.2% in T2DM vs 32-34% in healthy (p=0.008), HOMA-%B fell
    50.4% in T2DM vs ~19% in healthy -- a second, mechanistically
    consistent finding, not a single-study artifact.

    Magnitude is deliberately NOT re-derived from IGI/HOMA-%B's own
    percentages -- those are different constructs than this platform's
    secretion_scale (an ODE mechanism scalar), and reusing their exact
    percentages would overclaim precision this parameter doesn't have.
    The T2D-specific fix keeps the SAME magnitude constant as the
    existing compensatory case, only flipping its sign: an honest,
    documented order-of-magnitude estimate, not a literature-measured
    percentage, consistent with this codebase's standing calibration
    policy (see e.g. drug_interactions.py's own Gb-percentage notes).

    Non-T2D behaviour is numerically UNCHANGED from glucocorticoid_
    modifiers()'s existing +10%-at-reference-dose, dose-capped at
    30mg-equivalent (same convention that function already uses) --
    only now also gated by the real PD-timing curve instead of applying
    unconditionally regardless of dose timing.
    """
    if dose_mg <= 0.0:
        return 1.0
    cap_ratio = GC_SECRETION_SCALE_DOSE_CAP_MG / GC_SECRETION_SCALE_REFERENCE_DOSE_MG
    relative_dose = dose_mg / GC_SECRETION_SCALE_REFERENCE_DOSE_MG
    dose_factor = min(relative_dose, cap_ratio) / cap_ratio
    effect_fraction = glucocorticoid_pd_effect_fraction(hours_since_ingestion)
    magnitude = GC_SECRETION_SCALE_MAGNITUDE * dose_factor * effect_fraction
    return 1.0 - magnitude if is_t2d else 1.0 + magnitude


# HPA-axis negative feedback: glucocorticoid receptor occupancy (whether
# from endogenous cortisol or exogenous steroid) suppresses CRH release --
# the well-documented mechanism behind chronic-steroid adrenal suppression.
# ASSUMPTION (isolated here, no single-paper rate constant found -- unlike
# the dose-equivalency constants above, this is a qualitative-behaviour
# calibration, not a literature-derived value): gain chosen so that
# high-dose exogenous steroid (which saturates GR occupancy near 1.0)
# produces ~70% CRH suppression, while stress-driven cortisol elevation
# alone (which raises GR occupancy far less) is only mildly damped -- real
# pulsatile endogenous HPA regulation resists self-suppression much more
# than continuous exogenous ligand does, a distinction this 3-state ODE
# cannot fully capture. Replace this constant if a literature value for
# glucocorticoid-negative-feedback gain on CRH becomes available.
HPA_FEEDBACK_GAIN = 2.5


def _crh_suppression(gr_occ: float) -> float:
    """Same formula as dydt()'s own inline computation -- factored out so
    the steady-state solver below and the ODE use an identical expression."""
    excess_gr_feedback = max(0.0, gr_occ - GR_REFERENCE)
    return max(0.0, 1.0 - HPA_FEEDBACK_GAIN * excess_gr_feedback)


def solve_cortisol_crh_steady_state(
    cortisol_basal: float, crh_basal: float, stress_level: float,
) -> tuple[float, float, float]:
    """
    Analytic simultaneous fixed point of the 3-state Cortisol/CRH/GR_occ
    system, given `cortisol_basal`/`crh_basal`/`stress_level` -- computed
    here, not assumed to equal the raw input parameters, for the same
    reason `GR_OCC_BASAL` is computed rather than hardcoded (KNOWN BUG
    #10).

    Structural bug this closes (Section 8 housekeeping, MASTER_PLAN,
    found during a second attempt at Track 1 row 1.38/Cushing's
    syndrome): the OLD `basal_state()` returned `[cortisol_basal,
    crh_basal, gr_occ_basal]` literally, which is a genuine fixed point
    of `dydt()` ONLY in the trivial case `stress_level=0` (where
    `gr_occ_basal` also happens to equal `GR_REFERENCE`) -- for ANY
    `stress_level > 0`, `dCRH`'s own formula (`crh_target = crh_basal *
    (1+stress_level*2) * suppression(gr_occ)`) pulls CRH away from
    `crh_basal` immediately, and `dCortisol`'s formula
    (`cortisol_production = ke*cortisol_basal*(crh/crh_basal)`) then
    drags Cortisol along with it -- confirmed via direct execution:
    `CortisolParameters(stress_level=1.0)` alone (no dose, no disease
    preset) failed `validate_y0()` (`CRH: d/dt=0.693147`), silently,
    since this pathway was first built (no existing test checked a
    long-trajectory CRH invariant, only the glucose peak shape).

    **A first fix attempt (normalizing `cortisol_production`'s own
    crh-scaling term against a fixed reference) was tried, verified
    WRONG via the full test suite before shipping, and reverted**: it
    made `cortisol_production` always reduce to exactly `cortisol_basal`
    at the fixed point regardless of how elevated CRH's own target was
    -- silently eliminating `stress_level`'s entire steady-state effect
    on cortisol, not just the unwanted transient drift (5 existing
    `stress_level` tests failed, e.g. `test_stress_raises_cortisol`,
    revealing the ONLY reason they ever passed was the original bug's
    own transient drift -- there was no other mechanism by which
    `stress_level` affected cortisol at all under that fix).

    **This fix instead solves the GENUINE coupled fixed point**,
    leaving `dydt()`'s own `cortisol_production` formula completely
    UNCHANGED (still `crh/crh_basal`, not a normalized ratio) -- so
    `stress_level`'s real elevating effect on cortisol is preserved,
    not cancelled. Derivation: at any true fixed point,
    `y0.crh/crh_basal` must equal BOTH (a) `crh_target(gr_occ_ss)/
    crh_basal = (1+stress_level*2)*suppression(gr_occ_ss)` (from
    `dCRH=0`) and (b) whatever ratio makes `dCortisol=0`, i.e.
    `cortisol_ss/cortisol_basal` (since `cortisol_production =
    ke*cortisol_basal*(crh_ss/crh_basal)` must equal `ke*cortisol_ss`
    at the fixed point) -- so `cortisol_ss = cortisol_basal *
    (1+stress_level*2)*suppression(gr_occ_ss)` and, substituting into
    `dGR=0`'s own `gr_occ_ss = hill_activation(cortisol_ss, EC50_GR_NM,
    HILL_GR)`, `crh_basal` cancels out ENTIRELY, leaving a single-
    variable equation in `gr_occ_ss` alone:

        gr_occ_ss = hill_activation(
            cortisol_basal * (1+stress_level*2) * suppression(gr_occ_ss),
            EC50_GR_NM, HILL_GR,
        )

    Solved via `brentq` (bracketed in [0, 1) -- `f(0) > 0` always since
    `hill_activation(...) >= 0`; `f(1-eps) < 0` for any physiologically
    realistic `cortisol_basal`/`stress_level`, confirmed by direct
    execution across the full `stress_level` range and a Cushing's-
    scale elevated `cortisol_basal`). Reduces to EXACTLY the pre-fix
    default behaviour for the trivial case (`stress_level=0`,
    `cortisol_basal=CORTISOL_BASAL_NM`): `gr_occ_ss==GR_REFERENCE`
    exactly, since `GR_REFERENCE` is itself defined as
    `hill_activation(CORTISOL_BASAL_NM, EC50_GR_NM, HILL_GR)` -- verified
    via direct execution (residual < 1e-9). Independently confirmed via
    a full round-trip check: plugging the solved `(cortisol_ss, crh_ss,
    gr_occ_ss)` into the UNMODIFIED `dydt()` gives derivatives at
    machine-epsilon (~1e-15 to 1e-17) for every `stress_level` in
    [0, 1], and `cortisol_ss` genuinely elevates with `stress_level`
    (e.g. 0.28 -> 0.484 at `stress_level=0.8`), confirming the real
    effect is preserved, not cancelled.

    Returns (cortisol_ss, crh_ss, gr_occ_ss).
    """
    if cortisol_basal <= 0.0:
        raise ValueError(f"cortisol_basal must be positive, got {cortisol_basal}")
    # Mandatory security-reviewer finding (MEDIUM): CortisolParameters
    # .__init__ clamps stress_level to [0,1] before it ever reaches here,
    # but this function is public and can be called directly (as this
    # module's own tests already do) -- an un-clamped negative value
    # previously crashed deep inside hill_activation()/_crh_suppression()
    # ("x must be non-negative"), a confusing error unrelated to the real
    # cause. Guarded explicitly here instead.
    if stress_level < 0.0:
        raise ValueError(f"stress_level must be non-negative, got {stress_level}")
    drive = cortisol_basal * (1.0 + stress_level * 2.0)

    def f(gr_occ: float) -> float:
        return hill_activation(drive * _crh_suppression(gr_occ), EC50_GR_NM, HILL_GR) - gr_occ

    # Mandatory code-reviewer finding (MEDIUM): brentq's own bracket
    # precondition (f(0)>0, f(1-eps)<0) can fail for a sufficiently large
    # cortisol_basal (confirmed via direct execution: cortisol_basal=1e6,
    # ~3.6 million times the default, breaks the bracket) -- brentq's own
    # raw "f(a) and f(b) must have different signs" is not a helpful error
    # for a caller who hit this via a unit-confusion bug (e.g. passing a
    # pM value where nM is expected) rather than a deliberately extreme
    # value. Re-raised here with a message naming the actual input that
    # caused it.
    try:
        gr_occ_ss = brentq(f, 0.0, 1.0 - 1e-9)
    except ValueError as exc:
        raise ValueError(
            f"solve_cortisol_crh_steady_state() could not find a fixed point for "
            f"cortisol_basal={cortisol_basal}, stress_level={stress_level} -- this "
            "usually means cortisol_basal is implausibly large (check for a unit-"
            "confusion bug, e.g. passing a pM/uM value where nM is expected)."
        ) from exc
    ratio = (1.0 + stress_level * 2.0) * _crh_suppression(gr_occ_ss)
    cortisol_ss = cortisol_basal * ratio
    crh_ss = crh_basal * ratio
    return cortisol_ss, crh_ss, gr_occ_ss


def solve_cortisol_basal_for_target_steady_state(
    target_cortisol_ss: float, crh_basal: float, stress_level: float,
) -> float:
    """
    Inverse of solve_cortisol_crh_steady_state(): given a DESIRED
    steady-state cortisol level (e.g. a real, measured clinical value
    like Cushing's syndrome's own literature ratio), finds the
    `cortisol_basal` INPUT parameter that achieves it once the coupled
    negative-feedback solve is applied.

    Why this is needed (found while implementing Track 1 row 1.38,
    Cushing's syndrome, on the fixed foundation): the coupled solve's
    own negative feedback means `cortisol_ss < cortisol_basal * (1+
    stress_level*2)` for any `gr_occ_ss > GR_REFERENCE` -- the achieved
    steady state is LESS than what a naive "just multiply cortisol_basal
    by the target ratio" approach would suggest. Confirmed by direct
    execution: passing `cortisol_basal = CORTISOL_BASAL_NM * 1.560`
    (Cushing's own literature ratio) directly into `CortisolParameters`
    achieves an ACTUAL steady-state ratio of only ~1.287, not 1.560 --
    the negative feedback damps roughly half the intended elevation.

    A clinical citation like Jezkova et al. 2019 (PMID:31695455, see
    `CortisolParameters.CUSHINGS_CORTISOL_RATIO`) reports the OBSERVED,
    ACHIEVED cortisol level in real patients -- which already reflects
    whatever feedback operates in real physiology -- so it is the
    ACHIEVED steady state (`cortisol_ss`) that should match the target
    ratio, not the raw `cortisol_basal` input parameter. This function
    solves for that input.

    Monotonicity (verified by direct execution across a wide sweep,
    input multiplier 1.0-3.0x default): `cortisol_ss` increases
    monotonically with `cortisol_basal`, so a single bracketed `brentq`
    solve is well-posed. Bracket: `[crh_basal-independent single default,
    a generous multiple of the target]` -- widened automatically in BOTH
    directions if the initial bracket doesn't contain a root, since the
    exact multiplier needed depends on `stress_level` (which shifts the
    feedback curve -- a large `stress_level` means even a small
    `cortisol_basal` input can already overshoot the target, so `lo`
    must be able to shrink, not just `hi` grow; found by code review
    while composing this with `stress_level` near its own upper bound,
    where the initial `lo = target*0.5` no longer undershoots).
    """
    if target_cortisol_ss <= 0.0:
        raise ValueError(f"target_cortisol_ss must be positive, got {target_cortisol_ss}")

    def f(cortisol_basal_input: float) -> float:
        cortisol_ss, _, _ = solve_cortisol_crh_steady_state(
            cortisol_basal_input, crh_basal, stress_level,
        )
        return cortisol_ss - target_cortisol_ss

    def _cannot_bracket() -> ValueError:
        return ValueError(
            f"solve_cortisol_basal_for_target_steady_state() could not bracket "
            f"a solution for target_cortisol_ss={target_cortisol_ss}, "
            f"stress_level={stress_level} -- check these inputs are physiologically "
            "plausible."
        )

    lo, hi = target_cortisol_ss * 0.5, target_cortisol_ss * 2.0
    while f(hi) < 0.0:
        hi *= 2.0
        if hi > target_cortisol_ss * 1000.0:
            raise _cannot_bracket()
    while f(lo) > 0.0:
        lo *= 0.5
        if lo < target_cortisol_ss * 1e-6:
            raise _cannot_bracket()
    return brentq(f, lo, hi)


class CortisolParameters:
    def __init__(
        self,
        cortisol_basal: float = CORTISOL_BASAL_NM,
        crh_basal:      float = CRH_BASAL_PM,
        gr_occ_basal:   float = GR_OCC_BASAL,
        circadian_hour: Optional[float] = None,
        stress_level:   float = 0.0,
        exogenous_gc_nM: float = 0.0,
        exogenous_gc_dose_time_min: Optional[float] = None,
    ) -> None:
        self.cortisol_basal  = cortisol_basal
        self.crh_basal       = crh_basal
        self.gr_occ_basal    = gr_occ_basal
        self.circadian_hour  = circadian_hour
        self.stress_level    = max(0.0, min(1.0, stress_level))
        # Peak (fully-onset) exogenous glucocorticoid contribution to GR
        # ligand, in cortisol-equivalent nM (see PREDNISONE_TO_CORTISOL_MG_
        # RATIO above) -- see exogenous_gc_nM_at() below for how this is
        # actually applied at a given ODE time.
        self.exogenous_gc_nM = max(0.0, exogenous_gc_nM)
        # ODE-integration time (minutes, same clock dydt()/output_signals()
        # receive as `t`) at which the exogenous dose was/will be taken --
        # negative for a dose already taken before the simulated window,
        # positive for one scheduled during it. None means "no real dosing
        # timing known" (a direct CortisolParameters(exogenous_gc_nM=...)
        # construction, or a zero dose): exogenous_gc_nM_at() then returns
        # exogenous_gc_nM UNCHANGED for every t, exactly this platform's
        # pre-existing behaviour before real PK/PD timing was added (see
        # CLAUDE.md) -- so any caller not using for_prednisone()'s new
        # dose_time_min parameter is completely unaffected by this feature.
        self.exogenous_gc_dose_time_min = exogenous_gc_dose_time_min

    @classmethod
    def for_prednisone(
        cls,
        dose_mg: float,
        dose_time_min: float = 0.0,
        **kwargs,
    ) -> "CortisolParameters":
        """
        Exogenous corticosteroid dosing (prednisone-equivalent).

        Converts a prednisone dose to its PEAK (fully-onset) contribution
        to the GR ligand pool, expressed in the same cortisol-equivalent
        nM scale as endogenous cortisol -- see the module-level conversion
        constants for the literature basis. Any other CortisolParameters
        kwarg (stress_level, circadian_hour, ...) can be composed with
        this via **kwargs, e.g. a stressed patient also on steroids.

        dose_time_min: ODE-integration time (minutes) at which the dose
            was/will be taken, relative to simulation start (t=0).
            Defaults to 0.0 (dose taken exactly at simulation start,
            matching DrugDose's own "0 = just taken" convention) -- the
            glycaemic EFFECT genuinely has not developed yet at that
            instant (see glucocorticoid_pd_effect_fraction() above),
            unlike this platform's other static drug modifiers.

        Known, disclosed limitation (code-reviewer finding, Section 8
        housekeeping's CRH/Cortisol coupled-fixed-point fix): `y0` is
        computed by `solve_cortisol_crh_steady_state()` from
        `cortisol_basal`/`stress_level` ONLY -- it does NOT fold in
        `exogenous_gc_nM` at t=0. For dose_time_min=0.0 (the default),
        this is exact (`exogenous_gc_nM_at(0)` is itself gated to 0.0 by
        `glucocorticoid_pd_effect_fraction()`, so nothing is missing).
        But composing a large negative `dose_time_min` (a dose already
        WELL into its onset window before simulation start) WITH
        `stress_level > 0` leaves a small residual GR_occ inconsistency
        (confirmed: `for_prednisone(30.0, dose_time_min=-1800.0,
        stress_level=0.5)` -> `validate_y0(tolerance=1e-8)` still flags
        `GR_occ: d/dt~=0.0085`) -- genuinely SMALLER than the pre-fix
        residual for the same scenario (which had BOTH a much larger
        CRH residual, ~0.347, and a larger GR_occ residual, ~0.0137),
        but not fully closed. Not fixed here -- extending the solver to
        also account for a NON-ZERO exogenous ligand at t=0 is a real,
        separate piece of work (the exogenous term is itself time-gated
        by a curve, not a simple additive constant), logged as a
        follow-on to this fix rather than attempted as a drive-by
        extension.
        """
        cortisol_equivalent_mg = dose_mg * PREDNISONE_TO_CORTISOL_MG_RATIO
        exogenous_gc_nM = (
            cortisol_equivalent_mg / DAILY_CORTISOL_PRODUCTION_MG
        ) * CORTISOL_BASAL_NM
        return cls(
            exogenous_gc_nM=exogenous_gc_nM,
            exogenous_gc_dose_time_min=dose_time_min,
            **kwargs,
        )

    # Track 1 row 1.38 (third implementation attempt, on the now-fixed
    # CRH/Cortisol coupled-fixed-point foundation, Section 8 housekeeping):
    # active (untreated) Cushing's syndrome -- chronic endogenous HPA-axis
    # hypercortisolism. Route (a) from this row's own scoping: reuse this
    # pathway's own GR_REFERENCE machinery rather than a static disease-
    # state p3 modifier (route (b), explicitly rejected as a double-
    # counting risk against this exact pathway).
    #
    # Cortisol ratio: **Jezkova J, Durovcova V, Wenchich L, et al. "The
    # relationship of mitochondrial dysfunction and the development of
    # insulin resistance in Cushing's syndrome." Diabetes Metab Syndr Obes.
    # 2019;12:1459-1471. PMID:31695455** (open access, fetched in full,
    # not just the abstract) -- Table 1, morning (08:00) plasma cortisol,
    # same RIA assay: controls 527+/-51 nmol/L (n=24) vs. active Cushing's
    # 822+/-80 nmol/L (n=24, P<0.05). Ratio 822/527=1.560, applied as a
    # multiplier on CORTISOL_BASAL_NM -- NOT a literal nM value
    # (CORTISOL_BASAL_NM is an internal relative scale, confirmed during
    # this row's own first scoping pass by reverse-engineering
    # for_prednisone()'s own conversion arithmetic; only within-study
    # ratios transfer safely onto it).
    #
    # Unlike the SECOND implementation attempt (reverted, see MASTER_PLAN
    # row 1.38's own history), this preset does NOT set gr_occ_basal
    # directly -- solve_cortisol_crh_steady_state() (the fixed foundation)
    # computes gr_occ_ss itself from cortisol_basal/stress_level alone,
    # so passing gr_occ_basal here would be silently ignored by
    # basal_state() (see that method's own comment).
    #
    # Naively setting cortisol_basal to CORTISOL_BASAL_NM*1.560 directly
    # is NOT correct, though -- confirmed by direct execution: the
    # coupled solve's own negative feedback means the ACHIEVED steady-
    # state cortisol_ss ends up at only ~1.287x default, not 1.560x
    # (the feedback damps roughly half the intended elevation). Jezkova
    # 2019's own 1.560 ratio is the OBSERVED/ACHIEVED cortisol level in
    # real patients (already reflecting whatever feedback operates in
    # real physiology), so it is `cortisol_ss` (the achieved value) that
    # must match this ratio, not the raw `cortisol_basal` input. Solved
    # via the new `solve_cortisol_basal_for_target_steady_state()`
    # inverse solver -- see that function's own docstring for the full
    # derivation.
    #
    # Composed magnitude, disclosed not just computed: at this 1.560
    # ACHIEVED ratio (stress_level=0, the default for a resting Cushing's
    # patient), p3_modifier computes to ~0.778 (a real ~22% reduction)
    # with the corrected, Rizza-1982-derived CORTISOL_P3_SENSITIVITY (see
    # that constant's own citation comment) -- milder than Nosadini et
    # al. 1983's own clamp-derived 0.513 ratio, but a genuine, no-longer-
    # trivial effect. The gap between the two numbers is not a bug to
    # chase: Nosadini's own abstract explicitly attributes its 0.513
    # ratio to a RESPONSIVENESS/Vmax defect, not the sensitivity/ED50
    # shape Bergman p3 represents (the SAME disclosed conflict already
    # documented on CORTISOL_P3_SENSITIVITY itself) -- this preset does
    # not attempt to force a match to Nosadini's own, differently-shaped
    # number.
    #
    # Gb (fasting glucose) is NOT touched by this preset -- verified by
    # direct execution, not assumed. gb_modifier's own formula computes
    # `excess_cortisol = max(0, cortisol - p.cortisol_basal) +
    # exogenous_gc_nM_at(t)`, relative to the INSTANCE's own
    # cortisol_basal (not a fixed population reference) -- so at this
    # preset's own basal steady state (cortisol == cortisol_basal
    # exactly, by construction of the coupled solve), excess_cortisol is
    # always 0 and gb_modifier stays neutral (1.0). This is a real,
    # disclosed gap, not an oversight or a false-completeness claim:
    # Cushing's well-documented fasting hyperglycemia is NOT represented
    # by this preset at all -- an honest partial implementation (insulin-
    # sensitivity/p3 only), matching the established precedent for other
    # Gb-or-p3-only partial presets in this codebase (GCK_MODY,
    # tirzepatide, acromegaly). A future chunk could add a real,
    # separately-anchored Gb magnitude -- not attempted here, out of this
    # row's own scope once the two foundational blockers were resolved.
    CUSHINGS_CORTISOL_RATIO = 822.0 / 527.0  # 1.560, Jezkova 2019 (PMID:31695455)

    @classmethod
    def for_cushings(cls, stress_level: float = 0.0, **kwargs) -> "CortisolParameters":
        """
        Active (untreated) Cushing's syndrome -- see this method's own
        preceding module-level comment (CUSHINGS_CORTISOL_RATIO) for the
        full citation, and solve_cortisol_basal_for_target_steady_state()
        for why an inverse solve (not a direct multiply) is required.
        Composable with stress_level/an exogenous glucocorticoid dose
        (a real patient can have Cushing's AND be acutely stressed AND
        on exogenous steroids; all three route through the SAME shared,
        saturating GR-occupancy pool by physiological design, not a
        double-count -- see cortisol_pathway.py's module-level
        discussion of the stress+steroid case, which composes the
        identical way). `stress_level` is accepted explicitly (not via
        **kwargs) because the inverse solve itself needs it -- the
        target ACHIEVED cortisol_ss must still be reached even when
        stress_level>0 shifts the feedback curve. Unlike the reverted
        second attempt, does NOT pass gr_occ_basal -- the fixed
        basal_state() computes the genuine coupled fixed point itself
        from cortisol_basal alone.
        """
        target_cortisol_ss = CORTISOL_BASAL_NM * cls.CUSHINGS_CORTISOL_RATIO
        cortisol_basal = solve_cortisol_basal_for_target_steady_state(
            target_cortisol_ss, kwargs.get("crh_basal", CRH_BASAL_PM), stress_level,
        )
        return cls(cortisol_basal=cortisol_basal, stress_level=stress_level, **kwargs)

    def exogenous_gc_nM_at(self, t_minutes: float) -> float:
        """
        Time-gated exogenous GC-receptor-ligand contribution at ODE
        integration time t_minutes.

        Returns exogenous_gc_nM UNCHANGED (constant) if
        exogenous_gc_dose_time_min is None or the dose is zero -- see
        that field's own docstring for why. Otherwise gates it by
        glucocorticoid_pd_effect_fraction()'s onset/peak/resolution
        curve.
        """
        if self.exogenous_gc_nM <= 0.0 or self.exogenous_gc_dose_time_min is None:
            return self.exogenous_gc_nM
        hours_since_dose = (t_minutes - self.exogenous_gc_dose_time_min) / 60.0
        return self.exogenous_gc_nM * glucocorticoid_pd_effect_fraction(hours_since_dose)

    def basal_state(self) -> list[float]:
        # Section 8 housekeeping (MASTER_PLAN): solve_cortisol_crh_steady_
        # state() replaces the old literal [cortisol_basal, crh_basal,
        # gr_occ_basal] -- that was a genuine fixed point of dydt() only
        # in the trivial default case (stress_level=0); ANY stress_level>0
        # (or a disease preset that shifts gr_occ_basal) drifted from t=0
        # onward, silently, since this pathway was first built. See that
        # function's own docstring for the full derivation and the failed
        # first fix attempt it corrects.
        #
        # This computation is INDEPENDENT of self.gr_occ_basal -- gr_occ_ss
        # is fully determined by cortisol_basal/stress_level alone (the
        # derivation shows crh_basal itself cancels out of that part of
        # the equation too). self.gr_occ_basal is therefore no longer
        # consulted here; kept as a stored attribute purely for
        # introspection/copy-construction (with_stress()/with_time_of_day()
        # below still pass it through) since no caller currently overrides
        # it independently -- if a future caller needs to, this is the
        # method to revisit.
        cortisol, crh, gr_occ = solve_cortisol_crh_steady_state(
            self.cortisol_basal, self.crh_basal, self.stress_level,
        )
        if self.circadian_hour is not None:
            # Circadian phase is a deliberate SNAPSHOT perturbation, not a
            # chronic steady state to solve for -- same "not itself a
            # fixed point, by design" category as for_metformin()'s own
            # intentional onset drift (KNOWN BUG #14's closing note) and
            # this pathway's own pre-existing, unchanged behaviour for
            # circadian_hour specifically; out of this fix's scope.
            phase = 2.0 * math.pi * (self.circadian_hour - CORTISOL_PEAK_HOUR) / 24.0
            cortisol = max(0.01, cortisol + CORTISOL_AMPLITUDE * math.cos(phase))
        return [cortisol, crh, gr_occ]


class CortisolPathway(BasePathway):
    """3-species cortisol axis ODE model."""

    NAME    = "cortisol_v1"
    VERSION = "1.0.0"

    def __init__(self, params: Optional[CortisolParameters] = None) -> None:
        self.params = params or CortisolParameters()
        self._y0 = self.params.basal_state()
        self._descriptor = PathwayDescriptor(
            name=self.NAME,
            version=self.VERSION,
            maturity=PathwayMaturity.BETA,
            description=(
                "3-species cortisol axis: plasma cortisol, CRH, GR occupancy. "
                "Circadian rhythm and stress response included."
            ),
            species_names=["Cortisol", "CRH", "GR_occ"],
            output_signals=["cortisol_nM", "gr_occupancy", "gb_modifier", "p3_modifier"],
            input_signals=[],
            dependencies=[],
            tags=["endocrine", "cortisol", "stress", "circadian", "HPA-axis"],
            literature_refs=["PMID:7033265", "PMID:8770019"],
        )

    @property
    def descriptor(self) -> PathwayDescriptor:
        return self._descriptor

    @property
    def species_names(self) -> list[str]:
        return list(self._descriptor.species_names)

    @property
    def species_units(self) -> list[str]:
        return ["nM", "pM", "dimensionless"]

    @property
    def y0(self) -> list[float]:
        return list(self._y0)

    def dydt(self, t: float, y: list[float], signals=None) -> list[float]:
        cortisol = max(0.0, y[0])
        crh      = max(0.0, y[1])
        gr_occ   = max(0.0, min(1.0, y[2]))
        p = self.params

        # HPA-axis negative feedback: GR occupancy above the fixed
        # population reference suppresses CRH release (real mechanism
        # behind chronic-steroid adrenal suppression -- see
        # HPA_FEEDBACK_GAIN above for the calibration caveat). Gated on
        # *excess* over GR_REFERENCE so it is exactly zero at the true
        # basal steady state (gr_occ == GR_REFERENCE there), preserving
        # the existing zero-stress/zero-steroid basal invariant.
        # Uses the SAME _crh_suppression() helper solve_cortisol_crh_
        # steady_state() solves against (Section 8 housekeeping) -- kept
        # as one shared function specifically so the ODE and its own
        # steady-state solver can never drift apart again.
        suppression = _crh_suppression(gr_occ)
        crh_target = p.crh_basal * (1.0 + p.stress_level * 2.0) * suppression
        dCRH = -(math.log(2) / CRH_HALF_LIFE) * (crh - crh_target)

        ke_cortisol = math.log(2) / CORTISOL_HALF_LIFE
        cortisol_production = ke_cortisol * p.cortisol_basal * (crh / p.crh_basal)
        dCortisol = cortisol_production - ke_cortisol * cortisol

        # GR occupancy is driven by total glucocorticoid receptor ligand --
        # endogenous cortisol plus any exogenous corticosteroid (prednisone
        # etc., see CortisolParameters.for_prednisone()), time-gated by
        # exogenous_gc_nM_at() (see its docstring -- constant/unchanged
        # behaviour unless real dose timing was supplied). Additive and
        # exactly zero contribution at the default exogenous_gc_nM=0.0, so
        # this cannot change behaviour for any existing (non-steroid)
        # scenario.
        total_gr_ligand = cortisol + p.exogenous_gc_nM_at(t)
        gr_target = hill_activation(total_gr_ligand, EC50_GR_NM, HILL_GR)
        dGR = approach_to_target(gr_occ, gr_target, GR_RELAXATION)

        return [dCortisol, dCRH, dGR]

    def output_signals(self, t: float, y: list[float]) -> dict[str, float]:
        cortisol = max(0.0, y[0])
        gr_occ   = max(0.0, min(1.0, y[2]))
        p = self.params

        # Exogenous glucocorticoid raises fasting glucose via the same
        # hepatic-gluconeogenesis mechanism as endogenous cortisol excess
        # (Rizza et al. 1982 demonstrated both the Gb and p3 effects) --
        # additive, zero at the default exogenous_gc_nM=0.0, and time-gated
        # by exogenous_gc_nM_at() (see its docstring). Soft-capped (rather
        # than plain-linear) via CORTISOL_GB_SATURATION_NM: KNOWN BUG (see
        # CLAUDE.md) -- the plain-linear version was unbounded in dose,
        # unlike p3_modifier below (which naturally saturates via GR
        # occupancy's own Hill kinetics). The cap engages at the same
        # ~30mg-prednisone-equivalent reference glucocorticoid_modifiers()
        # already treats as its own dose_factor ceiling, preserving the
        # ORIGINAL CORTISOL_GB_SENSITIVITY slope for small/typical excess
        # (denominator ~= 1 there) and only bending the curve at doses far
        # beyond clinical range.
        excess_cortisol = max(0.0, cortisol - p.cortisol_basal) + p.exogenous_gc_nM_at(t)
        gb_modifier = 1.0 + CORTISOL_GB_SENSITIVITY * excess_cortisol / (
            1.0 + excess_cortisol / CORTISOL_GB_SATURATION_NM
        )

        # Use fixed population reference for p3 — not instance-specific basal
        excess_gr = max(0.0, gr_occ - GR_REFERENCE)
        p3_modifier = max(0.3, 1.0 - CORTISOL_P3_SENSITIVITY * excess_gr)

        return {
            "cortisol_nM":  cortisol,
            "gr_occupancy": gr_occ,
            "gb_modifier":  gb_modifier,
            "p3_modifier":  p3_modifier,
        }

    def with_stress(self, stress_level: float) -> "CortisolPathway":
        new_params = CortisolParameters(
            cortisol_basal=self.params.cortisol_basal,
            crh_basal=self.params.crh_basal,
            gr_occ_basal=self.params.gr_occ_basal,
            circadian_hour=self.params.circadian_hour,
            stress_level=stress_level,
        )
        return CortisolPathway(params=new_params)

    def with_time_of_day(self, hour: float) -> "CortisolPathway":
        new_params = CortisolParameters(
            cortisol_basal=self.params.cortisol_basal,
            crh_basal=self.params.crh_basal,
            gr_occ_basal=self.params.gr_occ_basal,
            circadian_hour=hour,
            stress_level=self.params.stress_level,
        )
        return CortisolPathway(params=new_params)


def cortisol_couplings(gi_name: str, cortisol_name: str) -> list:
    """
    Build the CouplingSpecs connecting CortisolPathway to
    GlucoseInsulinPathway:

        cortisol.gb_modifier -> GI input "gb_modifier"  (elevates fasting glucose)
        cortisol.p3_modifier -> GI input "p3_modifier"   (reduces insulin action)

    No GI -> cortisol signal is needed: cortisol is self-driven by its own
    stress_level/circadian_hour parameters, not by GI state. Note this is
    deliberately NOT used for routine time-of-day meal scenarios -- those
    stay on the existing empirical circadian_modifiers() in
    biosim_physiology.chronobiology, which is already validated and in
    production. This coupling is for genuinely new capability: modelling
    acute/chronic psychological stress (CortisolParameters.stress_level)
    raising cortisol independent of time of day.
    """
    from .coupling import CouplingSpec

    return [
        CouplingSpec(
            source_pathway=cortisol_name, signal_name="gb_modifier",
            target_pathway=gi_name, target_input="gb_modifier",
        ),
        CouplingSpec(
            source_pathway=cortisol_name, signal_name="p3_modifier",
            target_pathway=gi_name, target_input="p3_modifier",
        ),
    ]


def couple_with_gi(gi: "GlucoseInsulinPathway", cortisol: "CortisolPathway") -> "CompositePathway":
    """
    Build the CompositePathway connecting a CortisolPathway to
    GlucoseInsulinPathway (see cortisol_couplings() for the signal contract).
    """
    from .coupling import CompositePathway

    return CompositePathway([gi, cortisol], couplings=cortisol_couplings(gi.name, cortisol.name))
