# SAVE TO: packages\ode-engine\src\biosim_ode\hepatic.py
"""
Hepatic glucose production (EGP) extension.

Endogenous glucose production (EGP) from the liver is the other half
of postprandial glucose regulation. During fasting, EGP maintains
glucose at ~5 mM. After a meal, rising insulin suppresses EGP.

The base GlucoseInsulinPathway uses a simplified glucose effectiveness
term (p1) that partially captures EGP suppression, but does not model
EGP as an independent compartment. This module extends the model with
an explicit EGP state variable for scenarios where hepatic dynamics matter:

    - Metformin simulation (primary action is EGP suppression)
    - Dawn phenomenon (early morning EGP surge)
    - Type 2 diabetes (impaired EGP suppression)
    - Glucagon response modelling (Sprint 16, future)

State vector extension:
    y[0]  EGP  [mM/min]  hepatic glucose production rate

The EGP state is designed to couple to GlucoseInsulinPathway via
CompositePathway, injecting an additional glucose production term.

Literature:
    Cherrington 1999, Diabetes. PMID:10331429
    Rizza 2010, Diabetes. PMID:20705776
    DeFronzo 2004, Med Clin North Am. PMID:15308380
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Optional

from .base_pathway import BasePathway, PathwayDescriptor, PathwayMaturity
from .kinetics import hill_inhibition, approach_to_target

if TYPE_CHECKING:
    from .coupling import CompositePathway
    from .pathway import GlucoseInsulinPathway

# Reference (population-average) basal insulin at which EGP_BASAL_MM_MIN is
# defined to hold -- matches PathwayParameters.Ib used throughout the model.
EGP_REFERENCE_INSULIN_PM = 66.0


def _solve_egp_suppression_ic50(
    egp_basal: float,
    egp_max: float,
    hill_n: float,
    insulin_basal_pM: float = EGP_REFERENCE_INSULIN_PM,
) -> float:
    """
    Solve for the insulin IC50 that makes egp_basal the exact analytic fixed
    point of `dEGP = relaxation * (egp_max * hill_inhibition(x, ic50, n) - EGP)`
    at insulin_basal_pM.

    hill_inhibition(x, ic50, n) = ic50**n / (ic50**n + x**n) = f
        => ic50 = x * (f / (1 - f)) ** (1/n),  f = egp_basal / egp_max

    egp_basal and egp_max are treated as the trustworthy literature anchors
    (measured fasting EGP rate vs. documented fasted/stressed maximum) --
    EGP_SUPPRESSION_IC50 (a dose-response shape parameter, not itself directly
    measured) is the derived quantity. Same anchor-vs-derived reasoning as
    CortisolPathway.GR_OCC_BASAL (see CLAUDE.md KNOWN BUGS FIXED #10).
    """
    f = egp_basal / egp_max
    if not 0.0 < f < 1.0:
        raise ValueError("egp_basal must be strictly between 0 and egp_max")
    return insulin_basal_pM * (f / (1.0 - f)) ** (1.0 / hill_n)


def _metformin_egp_reduction_fraction(dose_mg: float) -> float:
    """
    Dose-response fraction by which metformin reduces basal EGP.
    Shared by EGPParameters.for_metformin() and
    for_type2_diabetes_and_metformin() so the formula isn't duplicated.
    """
    reduction = 0.25 * math.log(1.0 + dose_mg / 500.0) / math.log(3.0)
    return min(reduction, 0.60)


# Default EGP parameters
EGP_BASAL_MM_MIN     = 0.032   # mM/min basal EGP (~1.9 mM/h for 70kg)
EGP_MAX_MM_MIN       = 0.080   # mM/min maximum EGP (fasted/stressed)
EGP_SUPPRESSION_N    = 1.2     # Hill coefficient
EGP_SUPPRESSION_IC50 = _solve_egp_suppression_ic50(
    EGP_BASAL_MM_MIN, EGP_MAX_MM_MIN, EGP_SUPPRESSION_N,
)   # pM insulin for 50% EGP suppression -- derived, not measured (see above)
EGP_RELAXATION       = 0.03    # min⁻¹ time constant for EGP changes


class EGPParameters:
    """Parameters for the hepatic glucose production model."""

    def __init__(
        self,
        egp_basal:         float = EGP_BASAL_MM_MIN,
        egp_max:           float = EGP_MAX_MM_MIN,
        ic50_pM:           float = EGP_SUPPRESSION_IC50,
        hill_n:            float = EGP_SUPPRESSION_N,
        relaxation:        float = EGP_RELAXATION,
        t2d_egp_scale:     float = 1.0,   # > 1.0 for impaired EGP suppression
        metformin_scale:   float = 1.0,   # < 1.0 for metformin effect
    ) -> None:
        self.egp_basal       = egp_basal
        self.egp_max         = egp_max
        self.ic50_pM         = ic50_pM
        self.hill_n          = hill_n
        self.relaxation      = relaxation
        self.t2d_egp_scale   = t2d_egp_scale
        self.metformin_scale = metformin_scale

    @classmethod
    def for_type2_diabetes(cls) -> "EGPParameters":
        """
        T2D hepatic glucose production, per Beck-Nielsen H, Hother-Nielsen O,
        Staehr P. "Is hepatic glucose production increased in type 2
        diabetes mellitus?" Curr Diab Rep. 2002. PMID:12643178: fasting HGP
        "may be increased by a mean of 12%" and "suppression of HGP by
        insulin shows a rightward shift of the dose response curve (reduced
        insulin sensitivity) but normal maximal suppression (no maximum
        velocity defect)."

        "No maximum velocity defect" means `egp_max` is unchanged from
        healthy, and there is no independent scale-up mechanism beyond the
        dose-response shift -- `t2d_egp_scale` stays at its default 1.0.
        The only two things that move are `egp_basal` (the anchored,
        measured +12% quantity) and `ic50_pM` (the "rightward shift",
        derived via the same analytic-fixed-point technique used for the
        default calibration -- see `_solve_egp_suppression_ic50()` -- so
        the elevated basal is *itself* the fixed point of the shifted
        dose-response curve, not an independent free parameter stacked on
        top of it as the previous *2.5 IC50 / *1.40 scale version was).
        """
        egp_basal_t2d = EGP_BASAL_MM_MIN * 1.12
        return cls(
            egp_basal=egp_basal_t2d,
            ic50_pM=_solve_egp_suppression_ic50(egp_basal_t2d, EGP_MAX_MM_MIN, EGP_SUPPRESSION_N),
        )

    @classmethod
    def for_hyperthyroidism(cls) -> "EGPParameters":
        """
        Hyperthyroidism hepatic glucose production, per Wennlund A, Felig P,
        Hagenfeldt L, Wahren J. "Hepatic glucose production and splanchnic
        glucose exchange in hyperthyroidism." J Clin Endocrinol Metab. 1986.
        PMID:3510000: n=8 hyperthyroid vs. 10 normal subjects, hepatic venous
        catheterization + [3H]3-glucose isotope dilution -- basal HGP
        "increased by 20%" (14.2 +/- 0.6 vs. 11.9 +/- 0.6 umol/kg.min,
        P<0.02), AND, separately, with a larger insulin rise (66% vs. 37%),
        HGP "was less markedly suppressed... suggesting the presence of
        hepatic resistance to insulin" (P<0.01).

        Exact same anchor-then-derive shape as for_type2_diabetes() above --
        `egp_basal` is the measured +20% anchor, `ic50_pM` is re-derived via
        `_solve_egp_suppression_ic50()` so the elevated basal is the genuine
        fixed point of the (rightward-shifted) suppression curve, not an
        independent parameter. This single derivation naturally reproduces
        BOTH of Wennlund's separately-measured findings (elevated basal AND
        reduced suppression by a given insulin rise) rather than needing two
        free parameters for what the source paper shows is one mechanism.
        `egp_max`/`t2d_egp_scale` are left at their healthy defaults --
        Wennlund reports no maximum-suppression-capacity defect, only a
        dose-response shift, same reasoning as for_type2_diabetes().

        This is the SOLE representation of hyperthyroidism's glucose effect
        (see DiseaseState.HYPERTHYROIDISM in physiology/modifiers.py, which
        deliberately returns no static Gb/p3 modifier) -- composing this
        preset via HepGlucPathway alongside a static modifier would
        double-count the identical mechanism.

        Investigated directly (not assumed) what this predicts for a
        postprandial meal, composed with GlucoseInsulinPathway: naive
        expectation was a meaningfully HIGHER glucose peak (impaired
        suppression -> more EGP during the meal response). Measured
        instead a "close to normal" postprandial trajectory (within ~1%
        of a default/neutral hepatic coupling's own AUC) -- the same
        "no maximum velocity defect, modest rightward IC50 shift" shape
        as for_type2_diabetes(), which produces a small effect by the
        same anchor-then-derive construction regardless of how large the
        basal anchor is (see test_hyperthyroidism_hepatic_coupling_
        postmeal_close_to_normal in test_hepatic.py for the numbers).
        This does not mean Wennlund's real, statistically significant
        finding is wrong -- it means this platform's specific pathway
        kinetics (relaxation rate, Hill shape) don't let it dominate over
        the meal-driven insulin excursion at this scale. No independent
        postprandial-glucose clinical benchmark for hyperthyroidism was
        found this session to check this prediction against directly.
        """
        egp_basal_hyper = EGP_BASAL_MM_MIN * 1.20
        return cls(
            egp_basal=egp_basal_hyper,
            ic50_pM=_solve_egp_suppression_ic50(egp_basal_hyper, EGP_MAX_MM_MIN, EGP_SUPPRESSION_N),
        )

    @classmethod
    def for_metformin(cls, dose_mg: float = 1000.0) -> "EGPParameters":
        """
        Metformin primarily reduces hepatic glucose output -- the fasting
        EGP rate itself, not merely the shape of the insulin-suppression
        curve. At 1000 mg/day: ~25% EGP reduction.

        Regression: this used to be represented as a `metformin_scale`
        multiplier applied only to the suppression TARGET
        (`egp_target = egp_max * suppression * metformin_scale`),
        leaving `egp_basal`/`ic50_pM` at their unmodified defaults.
        Unlike `for_type2_diabetes()` (which correctly re-derives
        `ic50_pM` from its own elevated `egp_basal` via
        `_solve_egp_suppression_ic50()`), this meant `egp_basal` was NOT
        the analytic fixed point of `dEGP` for a metformin-only preset --
        confirmed directly: `dEGP` at t=0/basal insulin was -0.00024
        mM/min, not 0. Composed with `GlucoseInsulinPathway` (whose own
        basal state already reflects the SEPARATE, documented static Gb
        modifier from `physiology/drug_interactions.py`'s
        `metformin_modifiers()` -- -14.7% as of the mechanism-weighting
        correction in CLAUDE.md, was -4% when this fix was originally
        written; the fixed-point property below doesn't depend on the
        modifier's exact magnitude), this meant even a
        fasted (zero-meal) metformin scenario drifted glucose downward
        over time purely as an artifact -- confirmed: -0.10 mM by 4h,
        double-counting on top of the intended static effect, not a
        genuine additional mechanism.

        Fixed by applying the exact same anchor-then-derive technique
        already used for T2D, just in the opposite direction: lower
        `egp_basal` itself by the dose-dependent reduction fraction, then
        re-derive `ic50_pM` so the new, lower basal EGP is the genuine
        fixed point. `metformin_scale` is retired (always 1.0 now) since
        the effect is represented directly via `egp_basal`, the same way
        T2D's effect is. The 0.60 reduction cap preserves the same
        practical ceiling as the previous `metformin_scale` floor of 0.40.
        """
        reduction = _metformin_egp_reduction_fraction(dose_mg)
        egp_basal_metformin = EGP_BASAL_MM_MIN * (1.0 - reduction)
        return cls(
            egp_basal=egp_basal_metformin,
            ic50_pM=_solve_egp_suppression_ic50(egp_basal_metformin, EGP_MAX_MM_MIN, EGP_SUPPRESSION_N),
        )

    @classmethod
    def for_type2_diabetes_and_metformin(cls, dose_mg: float = 1000.0) -> "EGPParameters":
        """
        A T2D patient on metformin: T2D raises basal EGP by 12%
        (`for_type2_diabetes()`); metformin lowers it by its own
        dose-dependent fraction (`for_metformin()`). Applies both
        multiplicative adjustments to the same EGP_BASAL_MM_MIN anchor
        and re-derives ic50_pM once for the combined basal, so the
        combined preset is a genuine joint fixed point -- exactly the
        same anchor-then-derive pattern each preset uses alone, rather
        than composing two independently-derived ic50_pM values (which
        would not be a fixed point for the combined basal).
        """
        reduction = _metformin_egp_reduction_fraction(dose_mg)
        combined_basal = EGP_BASAL_MM_MIN * 1.12 * (1.0 - reduction)
        return cls(
            egp_basal=combined_basal,
            ic50_pM=_solve_egp_suppression_ic50(combined_basal, EGP_MAX_MM_MIN, EGP_SUPPRESSION_N),
        )


class HepGlucPathway(BasePathway):
    """
    1-species hepatic glucose production ODE model.

    Designed to couple to GlucoseInsulinPathway via CompositePathway.
    EGP output signal is injected into the glucose ODE as an
    additional production term.

    Coupling spec:
        CouplingSpec(
            source_pathway="hep_gluc_v1",
            signal_name="egp_mM_per_min",
            target_pathway="glucose_insulin_v1",
            target_input="hepatic_egp_delta",
        )
    """

    NAME    = "hep_gluc_v1"
    VERSION = "1.0.0"

    def __init__(self, params: Optional[EGPParameters] = None) -> None:
        self.params = params or EGPParameters()
        self._y0 = [self.params.egp_basal]
        self._descriptor = PathwayDescriptor(
            name=self.NAME,
            version=self.VERSION,
            maturity=PathwayMaturity.BETA,
            description=(
                "1-species hepatic glucose production model. "
                "EGP suppressed by insulin via Hill kinetics. "
                "Couples to GlucoseInsulinPathway as additional glucose source."
            ),
            species_names=["EGP"],
            output_signals=["egp_mM_per_min", "egp_suppression_fraction"],
            input_signals=["insulin_pM"],
            dependencies=[],
            tags=["hepatic", "liver", "EGP", "gluconeogenesis", "glycogenolysis"],
            literature_refs=["PMID:10331429", "PMID:20705776"],
        )

    @property
    def descriptor(self) -> PathwayDescriptor:
        return self._descriptor

    @property
    def species_names(self) -> list[str]:
        return ["EGP"]

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
        egp = max(0.0, y[0])
        p   = self.params

        # Get insulin from coupled pathway (or use basal if standalone)
        insulin_pM = 66.0
        if signals and "insulin_pM" in signals:
            insulin_pM = max(0.0, signals["insulin_pM"])

        # EGP target: maximum * insulin suppression * disease/drug modifiers
        suppression = hill_inhibition(insulin_pM, p.ic50_pM, p.hill_n)
        egp_target  = p.egp_max * suppression * p.t2d_egp_scale * p.metformin_scale
        egp_target  = max(0.0, min(egp_target, p.egp_max))

        dEGP = approach_to_target(egp, egp_target, p.relaxation)
        return [dEGP]

    def output_signals(self, t: float, y: list[float]) -> dict[str, float]:
        egp = max(0.0, y[0])
        p   = self.params
        suppression_fraction = 1.0 - (egp / p.egp_max) if p.egp_max > 0 else 0.0
        return {
            "egp_mM_per_min":          egp,
            "egp_suppression_fraction": suppression_fraction,
        }

    def basal_egp(self) -> float:
        return self.params.egp_basal


def hepatic_couplings(gi_name: str, hep: "HepGlucPathway") -> list:
    """
    Build the CouplingSpecs connecting a HepGlucPathway to
    GlucoseInsulinPathway:

        GI.insulin_pM        -> hep input "insulin_pM"  (drives EGP suppression)
        hep.egp_mM_per_min    -> GI input "hepatic_egp_delta"

    The transform subtracts THIS INSTANCE's own egp_basal (not the module
    default) before injection, so the coupling contributes zero at t=0
    regardless of preset (default/T2D/metformin) -- the delta formulation
    from CLAUDE.md rule 2, generalised. This means an elevated T2D fasting
    glucose still needs its own static Gb modifier (physiology/modifiers.py
    disease_modifiers) for the set-point shift; this coupling adds the
    genuinely distinct, complementary effect of impaired/enhanced EGP
    *kinetics* during the postprandial insulin response, which the static
    modifier cannot represent and which does not double-count the set-point
    shift (verified: hepatic_egp_delta is exactly 0 at this pathway's own
    basal state for every existing preset).
    """
    from .coupling import CouplingSpec

    basal = hep.params.egp_basal

    return [
        CouplingSpec(
            source_pathway=gi_name, signal_name="insulin_pM",
            target_pathway=hep.name, target_input="insulin_pM",
        ),
        CouplingSpec(
            source_pathway=hep.name, signal_name="egp_mM_per_min",
            target_pathway=gi_name, target_input="hepatic_egp_delta",
            transform=lambda egp, _basal=basal: egp - _basal,
        ),
    ]


def couple_with_gi(gi: "GlucoseInsulinPathway", hep: "HepGlucPathway") -> "CompositePathway":
    """
    Build the CompositePathway connecting a HepGlucPathway to
    GlucoseInsulinPathway (see hepatic_couplings() for the signal contract).
    """
    from .coupling import CompositePathway

    return CompositePathway([gi, hep], couplings=hepatic_couplings(gi.name, hep))
