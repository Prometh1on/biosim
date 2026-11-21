# SAVE TO: packages\ode-engine\src\biosim_ode\renal_pathway.py
"""
Renal glucose clearance pathway -- SGLT2 inhibitors (dapagliflozin,
empagliflozin).

SGLT2 (sodium-glucose cotransporter 2) in the proximal renal tubule
reabsorbs ~80-90% of filtered glucose. SGLT2 inhibitors block this,
causing genuine glucosuria (urinary glucose loss) -- an insulin-
INDEPENDENT mechanism, unlike every other drug pathway this platform
models (metformin/hepatic, GLP-1, cortisol are all insulin- or
meal-timing-mediated).

Architecturally distinct from every existing extension pathway
(hepatic/GLP1/cortisol/NEFA): this is the first pathway whose exported
signal is a genuine function of GI's own instantaneous state (glucose),
not just of time or the pathway's own internal dynamics. CompositePathway
computes every pathway's output_signals(t, y_slice) from ONLY that
pathway's own state slice, BEFORE any coupling injection -- so a signal
that depends on glucose cannot be computed in a stateless pathway (no
GI signal is available in output_signals at all). This pathway resolves
that by carrying one real state variable (the renal clearance rate
itself) that relaxes toward a glucose-dependent target computed inside
dydt() (which DOES receive injected signals), exactly the same
"input signal -> internal target -> relaxing state -> output signal"
shape HepGlucPathway already uses for insulin -> EGP.

Mechanism and design (see CLAUDE.md Track 1 SGLT2i milestone for the
full literature-verification trail):

    Filtered load:  F(G) = GFR_L_per_min * G                    [mmol/min]
    Reabsorption:   R(F) = (1 - beta) * F * hill_inhibition(F, Tm, p) ** (1/p)
    Excretion:      E(F) = F - R(F)

    beta in [0, 1): concentration-INDEPENDENT fractional loss of
        reabsorptive capacity -- the parameter that represents the
        drug's dominant, verified effect: a collapse of the renal
        glucose threshold, active even at normal/low plasma glucose,
        NOT merely a rightward-shifted Tm-limited saturation curve (a
        hard Tm-limited model, beta=0, predicts ZERO glucosuria below
        Tm/GFR -- falsified by DeFronzo et al. 2013's own measured
        appearance threshold of 9.5-10.9 mM against a Tm/GFR-implied
        threshold of ~14.1 mM). beta=0 exactly reproduces the no-drug
        physiology (near-total reabsorption at ordinary glucose, Tm-
        limited saturation only at severe hyperglycaemia).
    Tm (mg/min), p (dimensionless splay index): the smooth-saturation
        shape of reabsorption as filtered load approaches capacity
        ("splay" in the renal physiology literature -- real nephrons
        transition gradually, not at a sharp cutoff). p is chosen (not
        measured) to simultaneously reproduce near-zero basal
        glucosuria and the measured 9.5 mM healthy appearance
        threshold; see SPLAY_INDEX_P below for the numeric trail.

The exported signal is *already* the delta (drug minus no-drug
excretion at the same instantaneous glucose, converted to mM/min) --
so composing this pathway with beta=0 (no active SGLT2i drug) yields
a signal of EXACTLY zero at every glucose level, not merely at basal --
a stronger neutrality guarantee than hepatic_egp_delta's (which is only
exactly zero at that pathway's own basal state). This is necessary, not
just tidy: unlike metformin's insulin-gated hepatic suppression (truly
inactive at basal insulin), SGLT2i's renal mechanism is continuously
active whenever plasma glucose exceeds the (drug-lowered) threshold --
and normal fasting glucose (~5 mM) is already above that lowered
threshold (<2 mM). A static Gb-lowering modifier alongside this pathway
(the two-part "static set-point shift + dynamic postprandial kinetics"
pattern used for metformin) would double-count exactly the way
CortisolPathway/glucocorticoid_modifiers() did (see CLAUDE.md KNOWN
BUGS): here the ODE pathway is the SOLE, complete representation, and
no static Gb modifier is added anywhere for SGLT2i drugs (see
biosim_physiology.drug_interactions._sglt2_inhibitor_modifiers(), which
deliberately returns {} for the same reason _glp1_agonist_modifiers()
does).

Baseline (non-SGLT2i) insulin-independent glucose handling is carried by
p1 alone: the old uncited `extra_g_clearance` "renal etc." term (0.018
min-1) double-counted p1 and is 0 since B7 (2026-11-16). The delta
formulation here does not touch either: at beta=0 this pathway's
contribution is exactly zero.

Literature (verified live against NCBI/Europe PMC full text this
session, not from memory):
    Al-Jobori H et al. 2017, Diabetes. PMID:28428225 -- empagliflozin,
        stepped hyperglycaemic clamp, TmG and glucosuria threshold.
    Polidori D et al. 2013, Diabetes Care. PMID:23412078 -- canagliflozin
        same-meal challenge (confounded by SGLT1 co-inhibition -- see
        biosim_validation's validate_sglt2i() docstring).
    DeFronzo RA et al. 2013, Diabetes Care 36(10):3169-76. PMID:23735727,
        PMCID:PMC3781504 -- dapagliflozin, pancreatic/stepped
        hyperglycaemic clamp with an explicit TmG/splay/threshold model;
        GFR, TmG, splay AUC, and appearance-threshold values below are
        read directly from this paper's full text (fetched and confirmed
        this session, not paraphrased from its abstract alone).
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Optional

from .base_pathway import BasePathway, PathwayDescriptor, PathwayMaturity
from .kinetics import hill_inhibition, approach_to_target
from .meal import GLUCOSE_MOLECULAR_WEIGHT, DEFAULT_VG_L_PER_KG, DEFAULT_BODY_WEIGHT_KG

if TYPE_CHECKING:
    from .coupling import CompositePathway
    from .pathway import GlucoseInsulinPathway


# Reference (population-average) fasting glucose used to compute this
# pathway's own y0 -- matches PathwayParameters.Gb's own default (5.0 mM),
# NOT an independent assumption.
GLUCOSE_REFERENCE_MM = 5.0

# GFR measured "similar" between T2D (127+/-34) and healthy (128+/-36)
# mL/min subjects in PMID:23735727 -- a single shared reference value is
# used rather than two near-identical presets.
GFR_REFERENCE_ML_PER_MIN = 127.5

# Baseline maximal renal glucose transport (TmG), healthy vs. type-2-
# diabetic subjects -- both re-fetched and re-verified live (full text
# for DeFronzo, abstract-only for Al-Jobori; no open-access full text
# exists for the latter) rather than trusted from an earlier session's
# note. DeFronzo et al. 2013 (dapagliflozin, PMID:23735727): "420 vs.
# 317 mg/min, excluding outlier; P < 0.001" (T2D 32% higher than
# healthy). Al-Jobori et al. 2017 (empagliflozin, PMID:28428225):
# "459 +/- 53 vs. 337 +/- 25 mg/min; P < 0.05" (T2D subjects had
# significantly greater TmG at baseline). Real, population-specific
# reabsorptive-capacity elevation in T2D -- both papers' own words,
# "elevated TmG, splay, and threshold compared with controls" --
# distinct from and NOT a substitute for splay (SPLAY_INDEX_P below,
# not split by population -- see that constant's own docstring for why).
DAPAGLIFLOZIN_TMG_HEALTHY_MG_PER_MIN = 317.0
DAPAGLIFLOZIN_TMG_T2D_MG_PER_MIN = 420.0
EMPAGLIFLOZIN_TMG_HEALTHY_MG_PER_MIN = 337.0
EMPAGLIFLOZIN_TMG_T2D_MG_PER_MIN = 459.0

# Splay index (dimensionless Hill-like shape exponent for the reabsorption
# saturation curve). NOT itself measured -- DeFronzo et al. 2013 reports
# splay as an AUC in mg^2/min^2 (13,598/27,328 mg^2/min^2 at baseline,
# healthy/T2D) whose defining integral is in that paper's Supplementary
# Data (not retrievable via Europe PMC's body XML this session). p=8 is
# chosen so that, at GFR_REFERENCE_ML_PER_MIN and a normal healthy Tm
# (317 mg/min), reabsorption(F) is BOTH (a) close to total at ordinary
# fasting glucose (5 mM) and (b) still shows genuine saturation onset
# near the measured healthy appearance threshold (9.5 mM, same paper) --
# p=1 (plain Michaelis-Menten) badly over-predicts basal glucosuria at
# beta=0 (>40 g/day at 5 mM, physiologically absurd for an undrugged
# kidney); p->infinity recovers a hard Tm cutoff, which DeFronzo's own
# measured threshold (9.5 mM, well below the Tm/GFR-implied 14.1 mM)
# rules out. This is an inferred shape parameter, not a citable number --
# treat accordingly.
SPLAY_INDEX_P = 8.0

# Fractional (concentration-independent) reabsorption-capacity loss under
# a therapeutic SGLT2i dose. Derived (not directly measured) from PMID:
# 23735727's own discussion-section statement: "the increase in 24-h
# urine glucose excretion (60-80 g/day) observed with all orally
# administered SGLT2 inhibitors currently in clinical trials represents
# only 40-50% of filtered glucose load" -- a class-wide claim, not
# specific to dapagliflozin's own 12+12 subjects. Approximating daily
# filtered load as GFR_REFERENCE_ML_PER_MIN held at a constant fasting
# glucose of 5 mM all day (the same simplification this paper itself
# uses for its own "~162 g/day" headline figure, confirmed to reproduce
# it exactly: 0.1275 L/min * 5 mmol/L * 180.16 g/mol * 1440 min/day =
# 165 g/day) and beta*F=E at low F/Tm (negligible saturation there),
# beta = UGE_g_per_day / filtered_load_g_per_day ~= 70/165 ~= 0.42. This
# is a documented approximation (a "what if glucose stayed at fasting
# all day" simplification, the same honesty standard as metformin's
# Gb-percentage flux-to-concentration approximation in drug_interactions.
# py) -- not read as more precise than it is, and NOT separately
# validated against a dedicated per-drug 24h UGE study.
SGLT2I_BETA_CLASS_DEFAULT = 0.40

# MASTER_PLAN 2.2: coefficient of variation for sampling beta in a Monte
# Carlo UQ ensemble (see SimulationRunner's UQ block) -- unlike most
# extension-pathway parameters in this codebase, PMID:23735727's own
# discussion-section statement above (the source SGLT2I_BETA_CLASS_DEFAULT
# is derived from) already contains a real, if imprecise, RANGE, not just
# a bare point estimate -- so this is a genuine derived uncertainty, not
# an invented one (contrast with GLP1Parameters.agonist_scale, which has
# NO reported dispersion anywhere in its own evidence base and is
# deliberately NOT assigned a CV for exactly that reason -- see this
# constant's own milestone entry in CLAUDE.md).
#
# The source states the SAME underlying fact two slightly different ways,
# and they don't agree exactly -- disclosed honestly rather than quietly
# picking whichever is more convenient:
#   (a) "60-80 g/day" UGE -> beta = UGE/165 (this constant's own approximation)
#       range [0.364, 0.485], width 0.121, SD (uniform approx) = width/sqrt(12)
#       = 0.0350 -> CV relative to the 0.40 anchor = 0.0350/0.40 ~= 0.088
#   (b) "40-50% of filtered load" -> beta range [0.40, 0.50] directly,
#       width 0.10, SD = 0.10/sqrt(12) = 0.0289 -> CV ~= 0.0289/0.40 ~= 0.072
# RENAL_BETA_CV = 0.08 is the average of the two (rounded), anchored to
# SGLT2I_BETA_CLASS_DEFAULT itself (0.40) rather than either phrasing's
# own slightly-higher range mean (0.4242 / 0.45) -- sampling around the
# ALREADY-SHIPPED point estimate keeps the UQ band's median matching the
# point-estimate simulation (the same reasoning pathway_factory's own
# fix exists to preserve), rather than silently re-centering the ensemble
# a few percent away from what the point estimate actually reports.
RENAL_BETA_CV = 0.08

# Numerical relaxation rate for the renal-clearance state -- fast
# relative to meal/insulin timescales (real tubular transport responds
# within minutes, far faster than the ~30-60 min postprandial glucose
# excursion this model resolves), chosen as a smoothing device so the
# state tracks its target without introducing solver stiffness, not
# itself a literature-measured transit-time constant.
RENAL_RELAXATION = 0.5   # min^-1 (~2 min time constant)

# Defensive upper bound on the glucose value this pathway's kinetics are
# evaluated at -- found necessary during self-verification (the mandatory
# review-agent pass hit a session-quota outage; this was checked directly
# instead, see CLAUDE.md): `hill_inhibition(F, tm, p) ** (1/p)` with
# Python floats raises an uncaught OverflowError for a sufficiently
# extreme filtered load (F**p overflowing float64), and glucose_mM=inf
# silently produces a NaN delta (inf * 0.0) with no error at all --
# either would propagate into the ODE solver from a single malformed or
# divergent upstream state, not just a deliberately adversarial input.
# 150 mM (~2700 mg/dL) is far beyond any physiologically real glucose
# value (the highest documented human case is ~147 mM) and well beyond
# where this model's excretion kinetics have any real resolution left
# (already deeply saturated long before this point) -- clamping here
# cannot affect any realistic simulated scenario.
GLUCOSE_UPPER_BOUND_MM = 150.0


def _filtered_load_mg_per_min(glucose_mM: float, gfr_ml_per_min: float) -> float:
    """Filtered glucose load, GFR * plasma glucose, in mg/min."""
    gfr_L_per_min = gfr_ml_per_min / 1000.0
    bounded_glucose = max(0.0, min(GLUCOSE_UPPER_BOUND_MM, glucose_mM))
    filtered_mmol_per_min = gfr_L_per_min * bounded_glucose
    return filtered_mmol_per_min * GLUCOSE_MOLECULAR_WEIGHT


def _excretion_mg_per_min(filtered_mg_per_min: float, tm_mg_per_min: float, p: float, beta: float) -> float:
    """
    Urinary glucose excretion rate (mg/min) given filtered load, using the
    smooth Tm-limited reabsorption model described in this module's
    docstring. beta=0 reproduces near-total reabsorption at ordinary
    glucose (the no-drug case); beta>0 adds a concentration-independent
    fractional loss representing SGLT2 blockade.
    """
    if filtered_mg_per_min <= 0.0:
        return 0.0
    reabsorbed = (1.0 - beta) * filtered_mg_per_min * hill_inhibition(filtered_mg_per_min, tm_mg_per_min, p) ** (1.0 / p)
    return filtered_mg_per_min - reabsorbed


class RenalClearanceParameters:
    """Parameters for the SGLT2-inhibitor renal glucose clearance model."""

    def __init__(
        self,
        gfr_ml_per_min: float = GFR_REFERENCE_ML_PER_MIN,
        tm_mg_per_min: float = 317.0,
        splay_p: float = SPLAY_INDEX_P,
        beta: float = 0.0,   # 0.0 = no drug (neutral)
        relaxation: float = RENAL_RELAXATION,
        body_weight_kg: float = DEFAULT_BODY_WEIGHT_KG,
        vg_L_per_kg: float = DEFAULT_VG_L_PER_KG,
    ) -> None:
        if gfr_ml_per_min <= 0:
            raise ValueError(f"gfr_ml_per_min must be positive, got {gfr_ml_per_min}")
        if tm_mg_per_min <= 0:
            raise ValueError(f"tm_mg_per_min must be positive, got {tm_mg_per_min}")
        if splay_p <= 0:
            raise ValueError(f"splay_p must be positive, got {splay_p}")
        if not 0.0 <= beta < 1.0:
            raise ValueError(f"beta must be in [0, 1), got {beta}")
        self.gfr_ml_per_min = gfr_ml_per_min
        self.tm_mg_per_min = tm_mg_per_min
        self.splay_p = splay_p
        self.beta = beta
        self.relaxation = relaxation
        self.body_weight_kg = body_weight_kg
        self.vg_L_per_kg = vg_L_per_kg

    @property
    def vg_L(self) -> float:
        return self.body_weight_kg * self.vg_L_per_kg

    def target_delta_mM_per_min(self, glucose_mM: float) -> float:
        """
        The glucose-clearance delta this pathway should relax toward at a
        given instantaneous glucose (mM): -(E_drug(G) - E_nodrug(G)) / Vg_L,
        converted from mg/min excretion to mM/min via GLUCOSE_MOLECULAR_
        WEIGHT and this instance's own glucose distribution volume. Exactly
        zero for every glucose value when beta=0.
        """
        filtered = _filtered_load_mg_per_min(glucose_mM, self.gfr_ml_per_min)
        e_drug = _excretion_mg_per_min(filtered, self.tm_mg_per_min, self.splay_p, self.beta)
        e_nodrug = _excretion_mg_per_min(filtered, self.tm_mg_per_min, self.splay_p, 0.0)
        extra_excretion_mg_per_min = e_drug - e_nodrug
        extra_excretion_mmol_per_min = extra_excretion_mg_per_min / GLUCOSE_MOLECULAR_WEIGHT
        return -extra_excretion_mmol_per_min / self.vg_L

    @classmethod
    def for_dapagliflozin(
        cls, t2d: bool = False, beta_class_default: float = SGLT2I_BETA_CLASS_DEFAULT,
    ) -> "RenalClearanceParameters":
        """
        Dapagliflozin 10 mg/day x 7 days -- DeFronzo et al. 2013,
        PMID:23735727. TmG and GFR are this study's own subject-group
        values (317 mg/min healthy, excluding the one reported outlier /
        420 mg/min type-2-diabetic; 127.5 mL/min GFR, shared across both
        groups since the paper reports GFR as "similar" between them --
        see GFR_REFERENCE_ML_PER_MIN docstring). beta uses the class-wide
        default (see SGLT2I_BETA_CLASS_DEFAULT) regardless of t2d -- this
        paper does not report a population-split 24h UGE figure, only the
        class-wide 60-80 g/day statement that default is derived from; a
        T2D patient's own elevated glucose (via the disease-state Gb
        modifier) is what raises filtered load and therefore excretion,
        not a different beta.

        t2d=True selects this paper's own measured, elevated type-2-
        diabetic reabsorptive capacity (see
        DAPAGLIFLOZIN_TMG_T2D_MG_PER_MIN's docstring) instead of always
        using the healthy value regardless of the patient's actual
        disease state (the prior behaviour of this classmethod). Splay
        (SPLAY_INDEX_P) is deliberately NOT adjusted for t2d even though
        the same paper also reports elevated splay in T2D -- a real,
        disclosed limitation, not an oversight: splay is reported only as
        an aggregate AUC (mg^2/min^2) whose defining integral is in this
        paper's Supplementary Data, not retrievable via Europe PMC's body
        XML, so there is no principled way to convert it into this
        model's dimensionless splay_p exponent (see SPLAY_INDEX_P's own
        docstring for the same reasoning applied to the healthy value).

        beta_class_default (MASTER_PLAN 2.2): overrides SGLT2I_BETA_
        CLASS_DEFAULT -- lets a Monte Carlo UQ ensemble sample the
        CLASS-WIDE default (see RENAL_BETA_CV) rather than this
        classmethod's own fixed point estimate, without re-deriving the
        per-drug construction logic. Defaults to the module constant, so
        omitting it reproduces this classmethod's exact prior behaviour.
        """
        return cls(
            gfr_ml_per_min=GFR_REFERENCE_ML_PER_MIN,
            tm_mg_per_min=(
                DAPAGLIFLOZIN_TMG_T2D_MG_PER_MIN if t2d else DAPAGLIFLOZIN_TMG_HEALTHY_MG_PER_MIN
            ),
            beta=beta_class_default,
        )

    @classmethod
    def for_empagliflozin(
        cls, days_of_treatment: float = 14.0, t2d: bool = False,
        beta_class_default: float = SGLT2I_BETA_CLASS_DEFAULT,
    ) -> "RenalClearanceParameters":
        """
        Empagliflozin -- Al-Jobori et al. 2017, PMID:28428225 (that
        paper's abstract does not itself state the empagliflozin dose or
        a GFR measurement, so GFR reuses GFR_REFERENCE_ML_PER_MIN as an
        assumption carried over from the dapagliflozin study -- both are
        SGLT2 inhibitors with a similar acute renal haemodynamic profile,
        but this specific number is NOT independently confirmed for
        empagliflozin; flagged here rather than silently reused).

        days_of_treatment scales beta by the ratio of THIS drug's own
        measured TmG-reduction fraction at the requested duration to its
        day-14 value -- a documented proxy: beta itself represents
        threshold collapse, not TmG reduction, but no duration-resolved
        UGE data exists for empagliflozin to scale beta against directly,
        and the two mechanisms are expected to mature together
        pharmacodynamically. Clamped to the measured 48h-to-day14 range;
        values outside it fall back to the nearest measured endpoint
        rather than extrapolating.

        t2d selects this paper's own type-2-diabetic cohort figures
        throughout, not just TmG: baseline TmG (459 vs. 337 mg/min,
        healthy) AND the day2/day14 TmG-reduction fractions the beta
        scaling above is derived from -- "reduced the TmG in both
        individuals with and without diabetes by 44 +/- 7 and 53 +/- 6%,
        respectively" (48h) and "further reduced... on day 14 (by
        65 +/- 5 and 75 +/- 3%, respectively)" -- both directly stated in
        the same abstract for both groups, re-verified live this session,
        not a healthy-cohort value reused for T2D by default (the prior
        behaviour of this classmethod).

        beta_class_default (MASTER_PLAN 2.2): overrides
        SGLT2I_BETA_CLASS_DEFAULT -- the value scaled by
        fraction_of_day14 below. Sampling THIS class-wide default (not
        the already-duration-scaled final beta) keeps the day2/day14
        onset shape intact under UQ; see for_dapagliflozin()'s own
        docstring for the same parameter.
        """
        # max()/min()-based clamping silently resolves NaN to whichever
        # bound it's compared against first (NaN comparisons are always
        # False) -- e.g. days_of_treatment=nan used to silently produce
        # the day-14 beta with no error, misrepresenting an invalid
        # duration as a confident, specific value. Reject non-finite
        # input explicitly instead (found during self-verification when
        # the mandatory review-agent pass hit a session-quota outage --
        # see CLAUDE.md).
        if not math.isfinite(days_of_treatment):
            raise ValueError(f"days_of_treatment must be finite, got {days_of_treatment}")
        if t2d:
            day2_reduction, day14_reduction = 0.44, 0.65   # T2D TmG reduction, PMID:28428225
            tm_mg_per_min = EMPAGLIFLOZIN_TMG_T2D_MG_PER_MIN
        else:
            day2_reduction, day14_reduction = 0.53, 0.75   # healthy TmG reduction, PMID:28428225
            tm_mg_per_min = EMPAGLIFLOZIN_TMG_HEALTHY_MG_PER_MIN
        days = max(2.0, min(14.0, days_of_treatment))
        fraction_of_day14 = (
            day2_reduction + (day14_reduction - day2_reduction) * (days - 2.0) / (14.0 - 2.0)
        ) / day14_reduction
        return cls(
            gfr_ml_per_min=GFR_REFERENCE_ML_PER_MIN,
            tm_mg_per_min=tm_mg_per_min,
            beta=beta_class_default * fraction_of_day14,
        )


class RenalPathway(BasePathway):
    """
    1-species SGLT2-inhibitor renal glucose clearance model.

    Designed to couple to GlucoseInsulinPathway via CompositePathway:
    reads glucose_mM as an input signal, exports a glucose-dependent
    clearance-delta signal (see this module's docstring for why a real
    state variable, not a stateless computation, is structurally
    required here).
    """

    NAME    = "renal_clearance_v1"
    VERSION = "1.0.0"

    def __init__(self, params: Optional[RenalClearanceParameters] = None) -> None:
        self.params = params or RenalClearanceParameters()
        self._y0 = [self.params.target_delta_mM_per_min(GLUCOSE_REFERENCE_MM)]
        self._descriptor = PathwayDescriptor(
            name=self.NAME,
            version=self.VERSION,
            maturity=PathwayMaturity.EXPERIMENTAL,
            description=(
                "1-species SGLT2-inhibitor renal glucose clearance model. "
                "Couples to GlucoseInsulinPathway as an additional, "
                "insulin-independent glucose clearance term."
            ),
            species_names=["RGC"],
            output_signals=["renal_clearance_delta_mM_per_min"],
            input_signals=["glucose_mM"],
            dependencies=[],
            tags=["renal", "sglt2", "glucosuria", "drug"],
            literature_refs=["PMID:28428225", "PMID:23412078", "PMID:23735727"],
        )

    @property
    def descriptor(self) -> PathwayDescriptor:
        return self._descriptor

    @property
    def species_names(self) -> list[str]:
        return ["RGC"]

    @property
    def species_units(self) -> list[str]:
        return ["mM/min"]

    @property
    def y0(self) -> list[float]:
        return list(self._y0)

    def dydt(
        self,
        t: float,
        y: list[float],
        signals: Optional[dict] = None,
    ) -> list[float]:
        rgc = y[0]
        glucose_mM = GLUCOSE_REFERENCE_MM
        if signals and "glucose_mM" in signals:
            glucose_mM = max(0.0, signals["glucose_mM"])

        target = self.params.target_delta_mM_per_min(glucose_mM)
        dRGC = approach_to_target(rgc, target, self.params.relaxation)
        return [dRGC]

    def output_signals(self, t: float, y: list[float]) -> dict[str, float]:
        return {"renal_clearance_delta_mM_per_min": y[0]}


def renal_couplings(gi_name: str, renal: "RenalPathway") -> list:
    """
    Build the CouplingSpecs connecting a RenalPathway to
    GlucoseInsulinPathway:

        GI.glucose_mM                          -> renal input "glucose_mM"
        renal.renal_clearance_delta_mM_per_min -> GI input "renal_glucose_clearance_delta"

    No transform/subtraction is needed on the renal->GI leg (unlike
    hepatic_couplings' egp-minus-basal transform): RenalPathway's own
    exported signal is already the delta (see
    RenalClearanceParameters.target_delta_mM_per_min's docstring), exactly
    zero at every glucose level when beta=0.
    """
    from .coupling import CouplingSpec

    return [
        CouplingSpec(
            source_pathway=gi_name, signal_name="glucose_mM",
            target_pathway=renal.name, target_input="glucose_mM",
        ),
        CouplingSpec(
            source_pathway=renal.name, signal_name="renal_clearance_delta_mM_per_min",
            target_pathway=gi_name, target_input="renal_glucose_clearance_delta",
        ),
    ]


def couple_with_gi(gi: "GlucoseInsulinPathway", renal: "RenalPathway") -> "CompositePathway":
    """
    Build the CompositePathway connecting a RenalPathway to
    GlucoseInsulinPathway (see renal_couplings() for the signal contract).
    """
    from .coupling import CompositePathway

    return CompositePathway([gi, renal], couplings=renal_couplings(gi.name, renal))
