"""
Glucose/Insulin/IGF-1/mTOR/AMPK pathway — v1 mechanistic model.

This is the core biological model: the Bergman minimal model (1981) --
three coupled ODEs for G, I and X -- plus three insulin-driven read-outs
(IGF-1, mTOR, AMPK). The read-outs depend on insulin but nothing in dG, dI
or dX depends on them, so they never change the glucose/insulin prediction.

State vector (internal units):
    y[0]  G     [mM]          Blood glucose
    y[1]  I     [pM]          Plasma insulin
    y[2]  X     [min⁻¹]      Remote insulin (interstitial compartment)
    y[3]  IGF1  [nM]          Plasma IGF-1
    y[4]  mTOR  [0–1]         Fractional mTORC1 activity
    y[5]  AMPK  [0–1]         Fractional AMPK activity

Time units: minutes (internal). Converted to seconds at the ODESystem boundary.

Literature sources for default parameters:
    - Bergman RN et al. J Clin Invest. 1981;68(6):1456-1467  (p1, p2, p3, ke, Gb, Ib)
    - Dalla Man C et al. IEEE Trans Biomed Eng. 2007 (Vg, kabs)
    - Sarbassov DD et al. Science. 2005 (mTOR/insulin coupling)
    - Hardie DG et al. Nat Rev Mol Cell Biol. 2012 (AMPK regulation)
    - Le Roith D. NEJM. 1997 (IGF-1 physiology)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional
import math

from .base_pathway import BasePathway, PathwayDescriptor, PathwayMaturity
from .kinetics import (
    hill_activation,
    hill_inhibition,
    michaelis_menten,
    approach_to_target,
    clamp,
)
from .meal import MealInput, GLUCOSE_MOLECULAR_WEIGHT, DEFAULT_VG_L_PER_KG


# bergman_ivgtt_y0() below used to default to its own hardcoded 0.286 L/kg,
# independently of (but coincidentally identical to) MealInput's old
# uncited default. Both have since been anchored to the same real
# citation -- Mora-Gonzalez D, Moreno-Cabanas A, Alvarez-Jimenez L,
# Morales-Palomo F, Ortega JF, Mora-Rodriguez R. "Glucose volume of
# distribution affects insulin sensitivity measured by intravenous
# glucose tolerance test." Scand J Med Sci Sports. 2024. PMID:38389141
# (verified live via PubMedFetcher): healthy-subject VdGLUCOSE = 230 +/-
# 49 mL/kg, within the paper's own cited "standard" range of
# 145-230 mL/kg -- and consolidated into ONE constant
# (`DEFAULT_VG_L_PER_KG`, defined in meal.py) rather than two separately-
# named ones that happen to share a value today but could silently drift
# apart later. IV and oral routes both dilute an appearing glucose dose
# into the same single glucose compartment (G) in this ODE, so there is
# no principled reason for them to use different distribution volumes.
BERGMAN_IVGTT_VG_L_PER_KG = DEFAULT_VG_L_PER_KG


# ---------------------------------------------------------------------------
# Species index constants — never use magic numbers in the ODE
# ---------------------------------------------------------------------------
G_IDX    = 0   # Blood glucose [mM]
I_IDX    = 1   # Plasma insulin [pM]
X_IDX    = 2   # Remote insulin effect [min⁻¹]
IGF1_IDX = 3   # Plasma IGF-1 [nM]
MTOR_IDX = 4   # mTORC1 activity [0-1]
AMPK_IDX = 5   # AMPK activity [0-1]

SPECIES_NAMES = ["G", "I", "X", "IGF1", "mTOR", "AMPK"]
SPECIES_UNITS = ["mM", "pM", "min-1", "nM", "dimensionless", "dimensionless"]
N_SPECIES = 6

# Mirrors biosim_physiology.personalizer._PARAM_BOUNDS["p3"] exactly (that
# module cannot be imported here -- biosim_ode is a lower-level package
# biosim_physiology depends on, not the reverse). Security-reviewer finding
# (Section 8 housekeeping, found while re-deriving CORTISOL_P3_SENSITIVITY):
# the STATIC p3 clamp in personalizer.py only ever bounds PathwayParameters.p3
# BEFORE integration starts -- it has no visibility into ODE-COUPLED
# multiplicative signals (p3_modifier from CortisolPathway,
# insulin_resistance_modifier from NEFAPathway), which are re-applied at
# EVERY solver step and can drive the real, integrated p3_eff below this
# floor with zero disclosure anywhere (confirmed by direct execution: an
# extreme realistic composition -- GDM+insulin_resistant+hypothyroidism+CKD+
# PCOS disease states, sleep/hypoxia/luteal-phase axes, indinavir/
# hydrochlorothiazide/growth_hormone/atorvastatin/glucocorticoid(120mg)
# drug events, stress_level=1.0, NEFA dynamics on -- drove the real,
# solver-integrated p3_eff down to ~4.55e-9, nearly triple the pre-fix
# breach of ~9.30e-9, since CORTISOL_P3_SENSITIVITY's own increase widens
# CortisolPathway's own p3_modifier range). Clamped here, at the actual
# point of use inside dydt(), so the ODE can never integrate with a p3
# value the platform's own design says is unphysiological, regardless of
# which modifier(s) composed it there.
P3_EFF_FLOOR = 1e-8
P3_EFF_CEILING = 1e-4

# Extends the same point-of-use clamp to the other three rule-11
# multiplicative signals with an established personalizer.py bound
# (Section 8 housekeeping follow-on to the p3_eff fix above -- the
# identical unclamped-at-point-of-use gap exists for every one of them,
# not just p3). Gb_eff/p1_eff/secretion_scale_eff mirror
# personalizer.py's own _PARAM_BOUNDS entries for "Gb"/"p1"/
# "secretion_scale" exactly. extra_g_clearance has NO personalizer.py
# entry; it is bounded to the same ceiling as p1_eff (both are first-order
# clearance-rate constants). Since B7 (2026-11-16) its default is 0 and its
# floor is 0, and exercise no longer scales it: exercise uptake is the
# separate additive exercise_glucose_uptake_delta (0.072 x intensity),
# which is not clamped here.
GB_EFF_FLOOR = 2.0
GB_EFF_CEILING = 15.0
P1_EFF_FLOOR = 0.001
P1_EFF_CEILING = 0.20
SECRETION_SCALE_EFF_FLOOR = 0.05
SECRETION_SCALE_EFF_CEILING = 3.0
EXTRA_G_CLEARANCE_EFF_FLOOR = 0.0   # default is 0 since B7 (2026-11-16)
EXTRA_G_CLEARANCE_EFF_CEILING = 0.20

# Single source of truth for the point-of-use clamps above, keyed by the
# PathwayParameters field each one bounds (also the key used by
# point_of_use_clamp_hits()).
EFFECTIVE_PARAM_BOUNDS: dict[str, tuple[float, float]] = {
    "p3": (P3_EFF_FLOOR, P3_EFF_CEILING),
    "Gb": (GB_EFF_FLOOR, GB_EFF_CEILING),
    "p1": (P1_EFF_FLOOR, P1_EFF_CEILING),
    "extra_g_clearance": (EXTRA_G_CLEARANCE_EFF_FLOOR, EXTRA_G_CLEARANCE_EFF_CEILING),
    "secretion_scale": (SECRETION_SCALE_EFF_FLOOR, SECRETION_SCALE_EFF_CEILING),
}


def raw_effective_params(p: "PathwayParameters", s: dict) -> dict[str, float]:
    """
    The coupled-signal-modified parameter values dydt() would use BEFORE
    its point-of-use clamps -- shared by dydt() and the clamp-disclosure
    check so the two can never compose a modifier differently.
    """
    return {
        "p3": p.p3 * s.get("p3_modifier", 1.0) * s.get("insulin_resistance_modifier", 1.0),
        "Gb": p.Gb * s.get("gb_modifier", 1.0),
        "p1": p.p1 * s.get("p1_modifier", 1.0),
        "extra_g_clearance": p.extra_g_clearance * s.get("extra_g_clearance_modifier", 1.0),
        "secretion_scale": p.secretion_scale * s.get("secretion_scale_modifier", 1.0),
    }


@dataclass
class PathwayParameters:
    """
    All kinetic parameters for the glucose/insulin pathway.

    Default values are from the literature (see module docstring).
    These can be overridden per-simulation to model genetic variants,
    disease states, or drug effects (Sprints 7 and 8).
    """
    # Glucose dynamics (Bergman 1981)
    Gb: float = 5.0          # mM, basal glucose (~90 mg/dL)
    p1: float = 0.0317       # min⁻¹, glucose effectiveness
    # min⁻¹, X compartment clearance. Credited to Pacini & Bergman 1986
    # (MINMOD, PMID:3640682), but VALUE UNVERIFIED (checked 2026-11-14): that
    # paper describes the fitting program, its abstract gives no p2 value,
    # and no full text was available to check. See docs/PROVENANCE_GAPS.md.
    p2: float = 0.0123
    # min⁻² per pM, insulin action on X. Bergman RN, Phillips LS, Cobelli C.
    # 1981, J Clin Invest, PMID:7033284, Fig 6 caption: SI = 5.1 +/- 0.7
    # x10⁻⁴ min⁻¹.uU⁻¹.mL for the combined normal-tolerant groups (n=11) --
    # NOT the parameter_seed.py-recorded 5.03e-6 (see CLAUDE.md KNOWN BUGS
    # FIXED #15: that seed entry's "from Table 1" citation could not be
    # reconciled with this directly-quoted prose value from the same paper
    # once Table 1 itself was actually read). p3 = (SI * p2) / (pM per
    # uU/mL, ~6.945) = (5.1e-4 * 0.0123) / 6.945.
    p3: float = 9.032e-7
    # min⁻¹, insulin-independent glucose clearance ON TOP OF p1. Zero since
    # 2026-11-16 (B7): it was an unsourced 0.018 ("renal etc."), but in the
    # minimal model p1 (S_G) already stands for ALL insulin-independent
    # disposal, so a resting extra term counted it twice (diagnosed in
    # docs/VALIDATION.md). Exercise uptake has its own additive signal,
    # exercise_glucose_uptake_delta. Kept as a field (default 0) so callers
    # and the extra_g_clearance_modifier signal keep working.
    extra_g_clearance: float = 0.0

    # Insulin dynamics
    Ib: float = 66.0              # pM, basal insulin (~11 μU/mL)
    # min⁻¹, insulin elimination (t½ ~7.5 min). Credited to Sherwin et al.
    # 1974 (PMID:4856884), but VALUE UNVERIFIED (checked 2026-11-14): the
    # abstract reports a three-compartment model and says plasma insulin
    # disappearance curves "are multiexponential", with no single rate
    # constant; the PMC scan has no text layer, so the tables were not read.
    # See docs/PROVENANCE_GAPS.md.
    ke: float = 0.093
    # MODEL ASSUMPTION (no source): mM, glucose level above which beta cells
    # secrete. No source range found; note Rebelato 2018 (PMID:30166558, rat
    # islets) saw secretion start at ">5.6 mM", which this value does not use.
    G_threshold: float = 4.5
    # mM, Michaelis constant for secretion (excess glucose above G_threshold).
    # Rebelato et al. 2018, Sci Rep, PMID:30166558: GSIS EC50 = 10.13 mM
    # absolute glucose in control (non-glucose-adapted) islets -> excess above
    # this model's G_threshold = 10.13 - 4.5 = 5.63 mM. Previously an uncited
    # 2.0 -- see CLAUDE.md VALIDATION STATUS for the calibration investigation
    # (peak_insulin deviation improved from -44.9% to -24.7% against the OGTT
    # benchmark with this literature-anchored value).
    KG_secretion: float = 5.63
    # pM per (mM/min), first-phase ("dynamic") insulin secretion gain.
    # Mari A, Schmitz O, Gastaldelli A, Oestergaard T, Nyholm B, Ferrannini E.
    # 2002, Am J Physiol Endocrinol Metab, PMID:12388151: insulin secretion
    # decomposes into a "static" component (function of glucose level --
    # delta_secretion below already is this) plus a "dynamic" component
    # proportional to the POSITIVE time-derivative of glucose (first-phase
    # release, fires only while glucose is rising). Needs no new ODE state
    # (dG is already computed earlier in the same dydt() call) and is
    # exactly zero at any basal steady state by construction (dG=0 there),
    # so it cannot perturb any existing steady-state invariant. The gain
    # itself (unlike the citation for the functional FORM) was fit against
    # the OGTT benchmarks, same methodology as GutAbsorptionModel
    # .ogtt_calibrated_glucose() -- confirmed safe for the Bergman IVGTT
    # benchmark too: glucose only ever falls after an IV bolus (no Ra
    # forcing term), so dG is never positive there and this term
    # contributes exactly zero throughout that entire simulation.
    # MODEL ASSUMPTION (magnitude): the functional form is from Mari 2002, but
    # 121.3 was fitted to the OGTT benchmark, not taken from a source.
    first_phase_gain: float = 121.3
    IGF1b: float = 13.0              # nM, basal IGF-1
    # min⁻¹, IGF-1 clearance = ln2 / IGF1_HALF_LIFE_MIN (810 min, 13.5 h).
    # The IGF1 state is TOTAL serum IGF-1 (IGF1b 13 nM ~ 100 ng/mL), most of
    # it in the IGFBP-3/ALS ternary (~150-200 kD) complex. Guler et al. 1989,
    # PMID:2558477: in the 50-kD and 200-kD complexes "apparent half-lives
    # were 20-30 min and 12-15 h", free IGF-I "10-12 min"; before infusion,
    # 92 of 152 ug/l was in the 200-kD complex. One compartment with the
    # ternary midpoint; the faster 50-kD and free pools are not modelled.
    # Was an uncited 0.007 (t½ ~99 min) until 2026-11-14.
    k_igf1_clear: float = math.log(2) / 810.0
    # MODEL ASSUMPTION (no source): fractional rise in IGF-1 production per
    # unit of relative insulin excess. No source range found.
    IGF1_ins_sensitivity: float = 0.1

    # mTOR dynamics (Sarbassov 2005)
    EC50_mTOR: float = 100.0  # pM, insulin EC50 for mTOR activation
    n_mTOR: float = 1.5       # Hill coefficient
    # MODEL ASSUMPTION (no source): min⁻¹, how fast mTOR activity follows
    # insulin. No source range found.
    k_mTOR: float = 0.05

    # AMPK dynamics (Hardie 2012)
    EC50_AMPK: float = 200.0  # pM, insulin IC50 for AMPK inhibition
    # MODEL ASSUMPTIONS (no source, no source range found): Hill steepness of
    # insulin's AMPK inhibition; min⁻¹ speed AMPK follows insulin; and the
    # 0-1 normalisation ceiling of AMPK activity.
    n_AMPK: float = 1.0
    k_AMPK: float = 0.02
    AMPK_max: float = 1.0

    # Modifier fields — genetic/physiological personalisation (neutral defaults)
    secretion_scale: float = 1.0   # multiplier on insulin secretion response (1.0 = unmodified)

    # pM/min, constant addition to insulin secretion, independent of ambient
    # glucose (neutral default 0.0). Represents sulfonylurea-class drugs
    # (glibenclamide/glimepiride/glipizide), which close the beta-cell KATP
    # channel via SUR1 binding and trigger secretion regardless of glucose
    # level -- the textbook mechanism behind their real hypoglycaemia risk.
    # Deliberately NOT a rule-11 CompositePathway signal or a
    # secretion_scale-style multiplicative modifier: secretion_scale
    # multiplies delta_secretion (pathway.py dydt(), below), which is
    # Michaelis-Menten glucose-gated and floors at exactly -ke*Ib below
    # G_threshold -- amplifying it would make a sulfonylurea appear to
    # SUPPRESS insulin harder during hypoglycaemia, the opposite of its
    # defining clinical risk (verified numerically before choosing this
    # design; see CLAUDE.md's sulfonylurea milestone). See
    # biosim_physiology.drug_interactions.glucose_independent_secretion_
    # pM_per_min() for the literature and calibration.
    glucose_independent_secretion: float = 0.0

    def __post_init__(self) -> None:
        errors = []
        if self.Gb <= 0: errors.append("Gb must be positive")
        if self.Ib <= 0: errors.append("Ib must be positive")
        if self.ke <= 0: errors.append("ke must be positive")
        if self.p1 <= 0: errors.append("p1 must be positive")
        if self.p2 <= 0: errors.append("p2 must be positive")
        if self.EC50_mTOR <= 0: errors.append("EC50_mTOR must be positive")
        if self.EC50_AMPK <= 0: errors.append("EC50_AMPK must be positive")
        if self.secretion_scale <= 0: errors.append("secretion_scale must be positive")
        if self.glucose_independent_secretion < 0:
            errors.append("glucose_independent_secretion must be >= 0")
        # Hard safety bound, not a literature-derived value: at or above
        # ke*Ib, the analytic insulin fixed point at zero glucose (see
        # basal_state() below: I_ss = Ib + glucose_independent_secretion/ke)
        # would reach or exceed 2x this instance's own basal Ib -- i.e. the
        # drug's glucose-independent floor alone, with no glucose-driven
        # secretion at all, would match or exceed a healthy person's entire
        # normal fasting insulin level. Kept strictly below that ceiling so
        # this remains a bounded, physiologically modest floor shift, not
        # an unbounded/implausible one.
        if self.glucose_independent_secretion >= self.ke * self.Ib:
            errors.append(
                "glucose_independent_secretion must be < ke * Ib "
                f"(got {self.glucose_independent_secretion}, ke*Ib={self.ke * self.Ib})"
            )
        if errors:
            raise ValueError("Invalid PathwayParameters:\n" + "\n".join(f"  - {e}" for e in errors))

    @property
    def F_secretion_basal(self) -> float:
        """
        Beta-cell secretion scale factor [pM/min].
        Set so that secretion exactly balances clearance at basal glucose.
        This guarantees the model is at steady state with no meal input.
        """
        dG_above = max(self.Gb - self.G_threshold, 0.0)
        if dG_above == 0:
            return self.ke * self.Ib  # fallback
        return self.ke * self.Ib * (self.KG_secretion + dG_above) / dG_above

    def basal_state(self) -> list[float]:
        """
        Compute the analytic steady-state values for all species.
        This is the initial condition for every simulation.
        """
        # Analytic insulin fixed point. At G=Gb (dG=0), delta_secretion and
        # first_phase_secretion are both exactly zero by F_secretion_basal's
        # own construction (dydt() below), so dI=0 requires
        # I_ss = Ib + glucose_independent_secretion/ke. Reduces to exactly
        # Ib when glucose_independent_secretion is at its neutral default
        # (0.0), preserving every existing scenario's basal state unchanged.
        I_ss = self.Ib + self.glucose_independent_secretion / self.ke

        # X (remote insulin) relaxes toward p3*max(I-Ib,0)/p2 (dX below) --
        # when I_ss != Ib, X's own analytic fixed point is correspondingly
        # nonzero too, not the previously-safe-to-hardcode 0.0. Getting this
        # wrong would repeat the exact "hardcoded basal value not derived
        # from its own relaxation formula" bug class already fixed several
        # times in this codebase (see CLAUDE.md KNOWN BUGS FIXED #10-14) --
        # X would silently drift away from 0 at t=0 with zero meal input.
        # Reduces to exactly 0.0 when I_ss == Ib.
        X_ss = self.p3 * max(I_ss - self.Ib, 0.0) / self.p2

        # IGF-1's own production term (dydt() below) is a function of
        # I/Ib, not a fixed constant -- its steady state must be re-derived
        # from I_ss for the same reason X's must. Reduces to exactly
        # IGF1b when I_ss == Ib (the sensitivity term becomes zero).
        IGF1_ss = max(0.0, self.IGF1b * (1.0 + self.IGF1_ins_sensitivity * (I_ss / self.Ib - 1.0)))

        # mTOR/AMPK steady states: relax toward a Hill function of insulin --
        # must use I_ss (the true fixed point), not the fixed Ib, for the
        # same reason. Reduces to exactly the prior formulas when I_ss == Ib.
        mTOR_ss = hill_activation(I_ss, self.EC50_mTOR, self.n_mTOR)
        AMPK_ss = self.AMPK_max * hill_inhibition(I_ss, self.EC50_AMPK, self.n_AMPK)

        return [
            self.Gb,      # G
            I_ss,         # I
            X_ss,         # X
            IGF1_ss,      # IGF1
            mTOR_ss,      # mTOR
            AMPK_ss,      # AMPK
        ]


class GlucoseInsulinPathway(BasePathway):
    """
    Bergman minimal model (coupled G, I, X) plus three insulin-driven
    read-outs (IGF-1, mTOR, AMPK) that do not feed back into G, I or X.

    This is the primary v1 pathway. It is designed to be called by ODESystem
    via its dydt() method, which scipy.integrate.solve_ivp will integrate.

    The model is parameterized from PathwayParameters, which in production
    are sourced from the ParameterStore (Sprint 1). Default parameters are
    literature-validated and produce physiologically realistic OGTT responses.
    """

    def __init__(
        self,
        params: Optional[PathwayParameters] = None,
        meals: Optional[list[MealInput]] = None,
    ) -> None:
        self.params = params or PathwayParameters()
        self.meals = meals or []
        self._y0 = self.params.basal_state()

    @property
    def y0(self) -> list[float]:
        """Initial conditions -- always start at analytic steady state."""
        return list(self._y0)

    @property
    def descriptor(self) -> PathwayDescriptor:
        return PathwayDescriptor(
            name="glucose_insulin_v1",
            version="1.0.0",
            maturity=PathwayMaturity.BETA,
            description=(
                "Bergman minimal model (1981): coupled glucose/insulin/remote "
                "insulin action, plus three insulin-driven read-outs "
                "(IGF-1, mTOR, AMPK) that do not feed back."
            ),
            species_names=list(SPECIES_NAMES),
            output_signals=["glucose_mM", "insulin_pM", "mTOR_activity", "AMPK_activity", "IGF1_nM"],
            # Regression: this used to list only 6 of the 9 signals
            # dydt() actually reads (see CLAUDE.md's "EXTENSION-PATHWAY
            # COUPLING CONTRACT" -- rule 11) -- p1_modifier,
            # extra_g_clearance_modifier, and secretion_scale_modifier
            # were missing, so this metadata undersold what the pathway
            # actually accepts. Found again during the exogenous-insulin
            # milestone's scoping pass: renal_glucose_clearance_delta (added
            # for the SGLT2i pathway) was ALSO missing here despite being
            # documented in dydt()'s own docstring below -- the exact same
            # "contract drifted from implementation" slip, a second time.
            # Fixed alongside adding exogenous_insulin_delta (the 10th key)
            # for this milestone. "meal*_remaining_mM" (added for
            # MASTER_PLAN 3.2, wiring CouplingSpec.validate() into
            # production) is a wildcard entry -- meal_absorption_pathway.py's
            # MealAbsorptionPathway exports a genuinely dynamic, one-per-
            # meal signal name (meal0_remaining_mM, meal1_remaining_mM,
            # ...) that no static list can enumerate; see
            # CouplingSpec.validate()'s own trailing-"*" prefix-match support.
            # "meal*_gut_mM" is the same idea for two-compartment meals.
            input_signals=[
                "p3_modifier", "gb_modifier", "insulin_resistance_modifier",
                "hepatic_egp_delta", "incretin_effect", "gastric_emptying_modifier",
                "p1_modifier", "extra_g_clearance_modifier", "secretion_scale_modifier",
                "renal_glucose_clearance_delta", "exogenous_insulin_delta",
                "nutrient_secretion_delta", "exercise_glucose_uptake_delta",
                "meal*_remaining_mM", "meal*_gut_mM",
            ],
            dependencies=[],
            tags=["metabolic", "glucose", "insulin", "mTOR", "AMPK"],
            literature_refs=["PMID:7033284", "PMID:17926672"],
        )

    def output_signals(self, t: float, y: list) -> dict:
        """Export named signals for inter-pathway coupling."""
        from .kinetics import clamp
        G    = clamp(y[0], lower=0.01)
        I    = clamp(y[1], lower=0.0)
        IGF1 = clamp(y[3], lower=0.0)
        mTOR = clamp(y[4], lower=0.0, upper=1.0)
        AMPK = clamp(y[5], lower=0.0, upper=1.0)
        return {
            "glucose_mM":    G,
            "insulin_pM":    I,
            "IGF1_nM":       IGF1,
            "mTOR_activity": mTOR,
            "AMPK_activity": AMPK,
        }

    @property
    def species_names(self) -> list[str]:
        return list(SPECIES_NAMES)

    @property
    def species_units(self) -> list[str]:
        return list(SPECIES_UNITS)

    def point_of_use_clamp_hits(
        self, signals: Optional[dict[str, float]],
    ) -> dict[str, tuple[str, float]]:
        """
        Which dydt() point-of-use clamps (EFFECTIVE_PARAM_BOUNDS) would
        engage for these coupled signals: {field: ("floor"|"ceiling", raw)}.
        Empty when every composed value is within bounds. A NaN raw value
        is reported as a floor hit, since clamp() maps NaN to the floor.
        """
        hits: dict[str, tuple[str, float]] = {}
        for name, raw in raw_effective_params(self.params, signals or {}).items():
            lower, upper = EFFECTIVE_PARAM_BOUNDS[name]
            if math.isnan(raw) or raw < lower:
                hits[name] = ("floor", raw)
            elif raw > upper:
                hits[name] = ("ceiling", raw)
        return hits

    def dydt(self, t: float, y: list[float], signals=None) -> list[float]:
        """
        Compute derivatives for all species at time t with state y.

        This is the function passed to scipy.integrate.solve_ivp.
        Time t is in minutes.

        Args:
            t: Current time [minutes]
            y: Current state vector [G, I, X, IGF1, mTOR, AMPK]
            signals: Optional modifiers injected by coupled extension pathways
                     (cortisol, hepatic, NEFA, GLP1) via CompositePathway.
                     All keys are optional; absent keys use the neutral
                     (no-effect) default so this pathway behaves identically
                     to today whether run standalone or composed:
                         p3_modifier                  mult. on p3   (default 1.0)
                         gb_modifier                   mult. on Gb   (default 1.0)
                         insulin_resistance_modifier   mult. on p3   (default 1.0)
                         hepatic_egp_delta             add. to dG    (default 0.0, mM/min)
                         incretin_effect               scales delta_secretion (default 0.0)
                         gastric_emptying_modifier      mult. on Ra   (default 1.0)
                         p1_modifier                    mult. on p1   (default 1.0)
                         extra_g_clearance_modifier     mult. on extra_g_clearance (default 1.0)
                         exercise_glucose_uptake_delta  add. to the clearance rate (default 0.0, min^-1)
                         secretion_scale_modifier       mult. on secretion_scale (default 1.0)
                         renal_glucose_clearance_delta   add. to dG    (default 0.0, mM/min --
                              SGLT2-inhibitor insulin-independent renal glucosuria,
                              see renal_pathway.py. Unlike hepatic_egp_delta, this is
                              exactly zero at EVERY glucose level (not just basal) when
                              no SGLT2i drug is active, since the mechanism itself is
                              continuously active whenever present, not gated by insulin.)
                         exogenous_insulin_delta        add. to dI    (default 0.0, pM/min --
                              injected insulin appearing in plasma from a subcutaneous
                              depot, see exogenous_insulin_pathway.py. Exactly zero
                              before/without an active dose and asymptotically zero
                              again afterward -- a genuinely TRANSIENT signal, unlike
                              glucose_independent_secretion [rule 8's PathwayParameters
                              field, sulfonylureas], which is a constant floor and does
                              NOT go through this signals dict at all.)
                         nutrient_secretion_delta        add. to dI    (default 0.0, pM/min --
                              protein-co-ingestion-driven ENDOGENOUS insulin secretion,
                              see nutrient_secretion_pathway.py. Same transient shape as
                              exogenous_insulin_delta [exactly zero before/without
                              protein, decaying back to exactly zero afterward] but a
                              distinct mechanism -- kept as its own key rather than
                              merged with either exogenous_insulin_delta or
                              glucose_independent_secretion.)

                     p1_modifier/extra_g_clearance_modifier exist for
                     genuinely time-varying couplings (e.g. an exercise
                     bout's during-vs-after transition, computed fresh at
                     the coupled pathway's own output_signals(t, y) on
                     every solver step -- CompositePathway.dydt() already
                     calls output_signals(t, ...) with the CURRENT
                     integration time, not just t=0, so this is safe to
                     vary continuously through a simulation, not only a
                     one-shot personalisation override).

        Returns:
            List of derivatives [dG/dt, dI/dt, dX/dt, dIGF1/dt, dmTOR/dt, dAMPK/dt]
        """
        p = self.params
        s = signals or {}

        # Unpack state (clamp to physical bounds)
        G    = clamp(y[G_IDX],    lower=0.01)   # mM
        I    = clamp(y[I_IDX],    lower=0.0)    # pM
        X    = clamp(y[X_IDX],    lower=0.0)    # min⁻¹
        IGF1 = clamp(y[IGF1_IDX], lower=0.0)   # nM
        mTOR = clamp(y[MTOR_IDX], lower=0.0, upper=1.0)
        AMPK = clamp(y[AMPK_IDX], lower=0.0, upper=1.0)

        # Effective parameters after coupled-pathway modifiers (neutral by default)
        raw = raw_effective_params(p, s)
        p3_eff = clamp(raw["p3"], lower=P3_EFF_FLOOR, upper=P3_EFF_CEILING)
        Gb_eff = clamp(raw["Gb"], lower=GB_EFF_FLOOR, upper=GB_EFF_CEILING)
        p1_eff = clamp(raw["p1"], lower=P1_EFF_FLOOR, upper=P1_EFF_CEILING)
        extra_g_clearance_eff = clamp(
            raw["extra_g_clearance"],
            lower=EXTRA_G_CLEARANCE_EFF_FLOOR, upper=EXTRA_G_CLEARANCE_EFF_CEILING,
        )
        secretion_scale_eff = clamp(
            raw["secretion_scale"],
            lower=SECRETION_SCALE_EFF_FLOOR, upper=SECRETION_SCALE_EFF_CEILING,
        )

        # -------------------------------------------------------------------
        # 1. Blood glucose: dG/dt
        # Ra(t): meal glucose appearance rate [mM/min], slowed by GLP-1's
        # gastric-emptying effect when present.
        # Clearance: insulin-independent (p1) + remote insulin (X) + extra,
        # relative to the (possibly cortisol-shifted) glucose set point.
        # hepatic_egp_delta carries the *excess* hepatic glucose output above
        # its own basal (delta formulation — see dI below) so an unmodified
        # hepatic pathway at its basal EGP contributes exactly zero here.
        #
        # Mass conservation: when a MealAbsorptionPathway is composed
        # (SimulationRunner does this whenever a gastric-emptying-modifying
        # drug, e.g. GLP-1, is active -- see meal_absorption_pathway.py),
        # its "meal{i}_remaining_mM" signal replaces the closed-form
        # meal.Ra(t) call for that meal, using this same modifier on the
        # RATE the remaining mass decays at rather than on Ra's output --
        # this conserves total carbohydrate mass for any modifier(t)
        # trajectory, whereas multiplying meal.Ra(t)'s output by a
        # modifier != 1 (constant or time-varying) permanently deletes
        # (1 - modifier) of the meal's mass instead of merely delaying it
        # (confirmed numerically: ~13-49% of meal CHO mass deleted across
        # the shipped GLP-1 drug presets before this fix). Falls back to
        # the closed form when no MealAbsorptionPathway is composed (the
        # overwhelming majority of scenarios). A two-compartment
        # (GutAbsorptionModel) meal is tracked as stomach + gut; its
        # "meal{i}_gut_mM" signal gives Ra = k_absorption * gut, with the
        # modifier already applied to stomach emptying inside that pathway.
        # -------------------------------------------------------------------
        modifier = s.get("gastric_emptying_modifier", 1.0)
        Ra = 0.0
        for i, meal in enumerate(self.meals):
            remaining_key = f"meal{i}_remaining_mM"
            gut_key = f"meal{i}_gut_mM"
            if meal.absorption_model is not None and gut_key in s:
                Ra += meal.absorption_model.k_absorption * s[gut_key]
            elif meal.absorption_model is None and remaining_key in s:
                gate = 1.0 if t >= meal.time_min else 0.0
                Ra += gate * meal.kabs * modifier * s[remaining_key]
            else:
                Ra += meal.Ra(t) * modifier
        total_clearance_rate = (
            p1_eff + X + extra_g_clearance_eff
            + s.get("exercise_glucose_uptake_delta", 0.0)
        )
        dG = (
            Ra - total_clearance_rate * (G - Gb_eff)
            + s.get("hepatic_egp_delta", 0.0)
            + s.get("renal_glucose_clearance_delta", 0.0)
        )

        # -------------------------------------------------------------------
        # 2. Plasma insulin: dI/dt
        # Delta-secretion formulation: at G=Gb the secretion term is exactly
        # zero, guaranteeing dI/dt=0 at the analytic basal steady state.
        #
        # delta_secretion = F_vmax * (f(G) - f(Gb))
        # where f(G) = max(G-thresh,0) / (KG + max(G-thresh,0))
        #
        # basal_frac is deliberately evaluated at the nominal (uncoupled) p.Gb,
        # not Gb_eff: it is the fixed calibration set point, not the current
        # physiological target. This preserves realistic beta-cell
        # compensation — e.g. under sustained cortisol-driven hyperglycaemia
        # (Gb_eff > p.Gb), G settles above the calibration point and
        # delta_secretion becomes persistently positive, compensating with
        # extra insulin, rather than silently re-zeroing at the new set point.
        # incretin_effect (GLP-1) amplifies this same delta term and is exactly
        # zero at basal glucose regardless of its magnitude, so it cannot
        # perturb the fasting steady state either.
        # -------------------------------------------------------------------
        Gb_above  = max(p.Gb - p.G_threshold, 0.0)
        dG_above  = max(G   - p.G_threshold, 0.0)
        basal_frac   = Gb_above / (p.KG_secretion + Gb_above) if Gb_above > 0 else 0.0
        current_frac = dG_above / (p.KG_secretion + dG_above) if dG_above > 0 else 0.0
        delta_secretion = (
            p.F_secretion_basal * (current_frac - basal_frac) * secretion_scale_eff
            * (1.0 + s.get("incretin_effect", 0.0))
        )
        # First-phase ("dynamic") secretion: fires only while glucose is
        # actively rising (dG > 0), reusing the dG already computed above --
        # zero at any basal steady state (dG=0 there) by construction, so it
        # cannot perturb the fasting steady state. See PathwayParameters
        # .first_phase_gain for the literature source (Mari et al. 2002).
        first_phase_secretion = p.first_phase_gain * max(dG, 0.0) * secretion_scale_eff
        # glucose_independent_secretion (sulfonylurea-class drugs) is added
        # unconditionally, unlike delta_secretion/first_phase_secretion --
        # it is deliberately NOT gated by glucose level (see its own
        # docstring on PathwayParameters). Neutral default 0.0 preserves
        # every existing scenario's dynamics unchanged.
        # exogenous_insulin_delta (injected insulin appearing in plasma from
        # a subcutaneous depot, see exogenous_insulin_pathway.py) is added
        # unconditionally too -- a genuinely TRANSIENT signal (exactly zero
        # before/without a dose and decaying back to exactly zero after),
        # unlike glucose_independent_secretion's constant floor above.
        # nutrient_secretion_delta (protein-co-ingestion-driven endogenous
        # insulin secretion, see nutrient_secretion_pathway.py) is the SAME
        # transient shape as exogenous_insulin_delta -- exactly zero before/
        # without protein and decaying back to exactly zero afterward --
        # but a DIFFERENT mechanism (endogenous beta-cell secretion
        # triggered by amino acids, not an injected external source), so it
        # is kept as its own signal key rather than folded into either
        # exogenous_insulin_delta or glucose_independent_secretion.
        dI = (
            -p.ke * (I - p.Ib) + delta_secretion + first_phase_secretion
            + p.glucose_independent_secretion + s.get("exogenous_insulin_delta", 0.0)
            + s.get("nutrient_secretion_delta", 0.0)
        )

        # -------------------------------------------------------------------
        # 3. Remote insulin (interstitial compartment): dX/dt
        # Bergman remote compartment — lags plasma insulin
        # -------------------------------------------------------------------
        dX = -p.p2 * X + p3_eff * max(I - p.Ib, 0.0)

        # -------------------------------------------------------------------
        # 4. Plasma IGF-1: dIGF1/dt
        # Production scales weakly with insulin (chronic insulin → IGF-1 axis)
        # Clearance: first-order
        # -------------------------------------------------------------------
        igf1_production = p.k_igf1_clear * p.IGF1b * (
            1.0 + p.IGF1_ins_sensitivity * (I / p.Ib - 1.0)
        )
        igf1_production = max(igf1_production, 0.0)
        dIGF1 = igf1_production - p.k_igf1_clear * IGF1

        # -------------------------------------------------------------------
        # 5. mTORC1 activity: dmTOR/dt
        # Activated by insulin via AKT (Hill kinetics)
        # Relaxes to insulin-dependent target
        # -------------------------------------------------------------------
        mTOR_target = hill_activation(I, p.EC50_mTOR, p.n_mTOR)
        dmTOR = approach_to_target(mTOR, mTOR_target, p.k_mTOR)

        # -------------------------------------------------------------------
        # 6. AMPK activity: dAMPK/dt
        # Inhibited by insulin/glucose signalling
        # Relaxes to insulin-dependent target
        # -------------------------------------------------------------------
        AMPK_target = p.AMPK_max * hill_inhibition(I, p.EC50_AMPK, p.n_AMPK)
        dAMPK = approach_to_target(AMPK, AMPK_target, p.k_AMPK)

        return [dG, dI, dX, dIGF1, dmTOR, dAMPK]

    def steady_state_derivatives(self) -> list[float]:
        """
        Compute derivatives at the analytic basal steady state.
        All derivatives should be ~0. Used for model validation.
        """
        return self.dydt(0.0, self.params.basal_state())

    def bergman_ivgtt_y0(
        self,
        dose_g: float = 20.0,
        body_weight_kg: float = 70.0,
        Vg_L_per_kg: float = BERGMAN_IVGTT_VG_L_PER_KG,
    ) -> list[float]:
        """
        Initial state for a Bergman-style IVGTT (intravenous glucose
        tolerance test): an instantaneous glucose bolus injected directly
        into plasma, bypassing gut absorption and hepatic first-pass
        extraction entirely (unlike an oral MealInput). Everything else
        (insulin, IGF1, mTOR, AMPK) starts at its normal analytic basal
        state -- only G is elevated at t=0.

        Default 20g dose in a 70kg subject matches Bergman RN et al. 1981,
        J Clin Invest, PMID:7033284, Fig 1 (the source for the
        `bergman_iv_glucose` validation benchmark). Pass this as
        `y0_override` to `ODESystem.simulate()` with `meals=[]` -- there is
        no Ra(t) forcing term for an IV bolus, so no MealInput is involved.

        `Vg_L_per_kg` defaults to `BERGMAN_IVGTT_VG_L_PER_KG` (module-level
        constant above, PMID:38389141) -- this used to default to the same
        uncited 0.286 L/kg as `MealInput.DEFAULT_VG_L_PER_KG`, which is
        higher than the literature-recommended range for a rapid-mixing
        IVGTT bolus (145-230 mL/kg) and over-dilutes the bolus, understating
        the true initial glucose spike. Correcting this moved
        `glucose_return_time_min` closer to the Bergman 1981 reference
        (46.0 -> 49.6 min vs. a 90 min reference) but does not close the
        gap entirely -- the residual deviation is consistent with the
        classic 2-compartment Bergman minimal model's well-documented
        structural simplification of long-tail glucose return-to-baseline
        dynamics beyond ~60-90 minutes, not a further-fixable calibration
        bug in this specific parameter. See CLAUDE.md VALIDATION STATUS for
        why `p1`/`p2` were deliberately NOT re-tuned to chase the rest of
        this gap (both are independently anchored to separate citations
        that also calibrate the already-passing OGTT benchmarks).

        Follow-up: `MealInput.DEFAULT_VG_L_PER_KG` (the oral-route
        constant) was initially left untouched here to avoid disturbing
        OGTT calibration without a dedicated re-validation pass -- that
        pass has since happened (see CLAUDE.md KNOWN BUGS FIXED #19) and
        both routes now share this same literature-anchored value.
        """
        carbs_mmol = dose_g * 1000.0 / GLUCOSE_MOLECULAR_WEIGHT
        Vg_L = body_weight_kg * Vg_L_per_kg
        dose_mM = carbs_mmol / Vg_L
        y0 = self.y0
        y0[G_IDX] = self.params.Gb + dose_mM
        return y0

    def with_meals(self, meals: list[MealInput]) -> "GlucoseInsulinPathway":
        """Return a new pathway instance with the given meal inputs."""
        return GlucoseInsulinPathway(params=self.params, meals=meals)

    @property
    def accepts_meals(self) -> bool:
        return True
