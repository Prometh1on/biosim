# SAVE TO: packages\ode-engine\src\biosim_ode\nefa_pathway.py
"""
Non-esterified fatty acid (NEFA/FFA) pathway.

Free fatty acids (NEFA) mediate the link between adipose tissue lipolysis
and insulin resistance. Elevated plasma NEFA:
    - Compete with glucose for oxidation (Randle cycle)
    - Activate PKC-theta, impairing insulin receptor signalling
    - Provide gluconeogenic substrate to the liver

This pathway is critical for modelling:
    - Obesity-associated insulin resistance (elevated basal NEFA -- but
      see NEFA_BASAL_UM/NEFA_OBESITY_UM's own comment below: two
      independent studies now contradict this preset's own implied
      magnitude, disclosed there, not resolved by a change here)
    - Fasting physiology (NEFA rise as glycogen depletes)
    - Exercise physiology (NEFA mobilised from adipose)

State vector:
    y[0]  NEFA  [uM]  plasma non-esterified fatty acid concentration

Output signals (for coupling to GlucoseInsulinPathway):
    "nefa_uM"           current NEFA
    "insulin_resistance_modifier"  p3 multiplier (< 1.0 at elevated NEFA)

Literature:
    Randle et al. 1963, Lancet. PMID:13990765  (original glucose-FFA cycle)
    Boden 1997, Diabetes. PMID:8971073          (NEFA and insulin resistance)
    Roden et al. 1996, J Clin Invest. PMID:8675698
    Heiling, Miles & Jensen 1991, Am J Physiol. PMID:1951681 (NEFA_HALF_LIFE_MIN)
    Hagenfeldt, Hagenfeldt & Wennmalm 1975, Horm Metab Res. PMID:1213651
        (corroborating NEFA_HALF_LIFE_MIN)
    Stumvoll et al. 2000, JCEM (NEFA_INSULIN_IC50_PM)
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Optional

from .base_pathway import BasePathway, PathwayDescriptor, PathwayMaturity
from .kinetics import hill_inhibition, hill_activation, approach_to_target

if TYPE_CHECKING:
    from .coupling import CompositePathway
    from .pathway import GlucoseInsulinPathway

# Reference (population-average) basal insulin at which each preset's
# nefa_basal is defined to hold -- matches PathwayParameters.Ib elsewhere.
NEFA_REFERENCE_INSULIN_PM = 66.0

# Default NEFA parameters
NEFA_BASAL_UM         = 350.0   # uM fasting (normal: 300-600 uM) -- still uncited
NEFA_OBESITY_UM       = 700.0   # uM in obese/insulin-resistant state -- still uncited
#
# NEFA_BASAL_UM/NEFA_OBESITY_UM: TWO independent search passes (2026-08-12
# deep-search plan item C1, then a fresh-angle retry 2026-08-18) have now
# each found real, live-verified evidence -- and BOTH point away from, not
# toward, this preset's own implied premise that obesity roughly DOUBLES
# fasting/basal plasma NEFA:
#
#   Curry SD, et al. 2011, PMID:21304064 (live-verified) -- lean/obese/
#   gastric-bypass adult cohort: "Fasting FFA were not significantly
#   different" between groups.
#
#   Kelsey MM, Forster JE, Van Pelt RE, Reusch JE, Nadeau KJ. "Adipose
#   tissue insulin resistance in adolescents with and without type 2
#   diabetes." Pediatr Obes. 2014;9(5):373-80. PMID:23861170 (live-
#   verified twice: once via PubMedFetcher for the abstract/bibliographic
#   details, once via a direct NCBI efetch of the PMC full text for the
#   actual table, since the abstract alone doesn't report the number).
#   n=47 adolescents (lean/obese/T2D), hyperinsulinaemic-euglycaemic
#   clamp. Table 1, "Baseline NEFA (umol/L)": All 558(254); Lean
#   521(317); Obese 583(155); T2D 561(273); p=0.86 -- i.e. NO significant
#   difference across groups, and the lean-vs-obese point estimate ratio
#   (583/521~1.12) is far below this preset's own implied ~2.0x.
#
# NOT changed to either study's own numbers: Curry reports no usable
# number at all; Kelsey's is a real, quotable number but from an
# ADOLESCENT (not adult) population, with small per-group n (10 lean, 12
# obese) and very wide SDs (317, 155) relative to the group means --
# substituting it as a precise adult default would trade one unverified
# number for a differently-flawed one (population mismatch, high
# variance, and a non-significant group difference presented as if it
# were a point estimate), the same false-precision failure mode this
# project has already rejected elsewhere (KG_secretion~200, joint-fit
# p1~0.008). Left unchanged, now formally PARKED per this project's own
# search-exhaustion rule (2+ independent passes with genuinely different
# angles, MASTER_PLAN row 1.19/Section 9) -- do not re-search without a
# specific new lead (a different primary source reporting both lean AND
# obese ADULT fasting NEFA/FFA in one low-variance cohort).
#
# NEFA_HALF_LIFE_MIN (2026-08-12 deep-search plan, item C2): re-derived
# from an uncited 90.0 to 3.7, anchored to Heiling VJ, Miles JM, Jensen
# MD. "How valid are isotopic measurements of fatty acid oxidation?"
# Am J Physiol. 1991;261(5 Pt 1):E572-7. PMID:1951681 (live-verified
# twice this session, not trusted from a single pass): "Studies in 11
# healthy volunteers... After stopping tracer infusions, the half-life
# of plasma [14C]FFA was 3.7 +/- 0.1 min." Independently corroborated
# by a SEPARATE research group using a different tracer/fatty-acid
# species: Hagenfeldt L, Hagenfeldt K, Wennmalm A. "Turnover of plasma
# free arachidonic and oleic acids in men and women." Horm Metab Res.
# 1975;7(6):467-71. PMID:1213651 -- fractional turnover rate 0.47-0.80
# min^-1 (men/women), converting via t1/2=ln(2)/k to an even FASTER
# ~0.87-1.5 min. Both sources agree on the same order of magnitude (a
# few minutes), consistent with the well-established physiological fact
# that the plasma NEFA pool is small and turns over extremely rapidly
# (unlike glucose/insulin) -- the prior 90.0 value was not merely
# uncited but a genuine, large (roughly 24-140x) error, not just an
# imprecise placeholder. Verified this pathway's own basal fixed-point
# invariant is mathematically independent of this constant (the target-
# based dNEFA formula's fixed point doesn't depend on the rate constant
# governing approach TO that fixed point) before changing it.
NEFA_HALF_LIFE_MIN    = 3.7     # min clearance half-life
#
# NEFA_INSULIN_IC50_PM (2026-08-12 deep-search plan, item C3):
# re-derived from an uncited 80.0 to 68.0, anchored to Stumvoll M,
# Jacob S, Wahl HG, Hauer B, Loblein K, Grauer P, Becker R, Nielsen M,
# Renn W, Haring H. "Suppression of systemic, intramuscular, and
# subcutaneous adipose tissue lipolysis by insulin in humans." J Clin
# Endocrinol Metab. 2000. PMID:11061533 (n=13 lean healthy volunteers,
# three-step hyperinsulinaemic-euglycaemic clamp, live-verified):
# "The insulin concentrations resulting in a half-maximal suppression
# (EC50) of systemic lipolysis, adipose tissue, and muscle lipolysis
# were 51, 68, and 44 pmol/L, respectively." Uses the ADIPOSE-TISSUE-
# specific value (68 pM) as the closest physiological match to "NEFA
# release from adipose tissue," this pathway's own stated mechanism --
# not the systemic (51 pM) or muscle-compartment (44 pM) values from
# the same paper.
NEFA_INSULIN_IC50_PM  = 68.0    # pM insulin for 50% NEFA suppression (antilipolysis)
NEFA_INSULIN_N        = 1.0
#
# NEFA_IR_EC50_UM / NEFA_IR_MAX -- RESOLVED (MASTER_PLAN row 1.22, a
# fresh-angle retry). The original pass (2026-08-12 deep-search plan,
# item C4, see git history) found real, on-topic dose-response
# citations -- Gormsen LC et al. 2007 (PMID:17341555), Boden G et al.
# 1994 (PMID:8200979) -- but rejected both: neither paper reports its
# own basal/reference FFA level alongside its dose-response curve, so
# no excess-above-basal EC50 could be computed without importing a
# mismatched external reference.
#
# Resolved by finding a paper that reports BOTH in the SAME cohort:
# Belfort R, Mandarino L, Kashyap S, Wirfel K, Pratipanawatr T, Berria
# R, Defronzo RA, Cusi K. "Dose-response effect of elevated plasma free
# fatty acid on insulin signaling." Diabetes. 2005;54(6):1640-1648.
# PMID:15919784 -- independently live-verified twice (once via the
# search pass, once directly, word-for-word identical). n=21 lean
# healthy normal-glucose-tolerant subjects, graded 4h Liposyn (lipid
# emulsion) infusion (30/60/90 ml/h) + 2h euglycaemic insulin clamp.
# Quoted directly: "Basal plasma FFA concentration (~440 micromol/l)
# was increased to 695, 1,251, and 1,688 micromol/l after 4 h of
# Liposyn infusion and resulted in a dose-dependent reduction in
# insulin-stimulated glucose disposal (Rd) by 22, 30, and 34%,
# respectively."
#
# This paper's OWN basal (440 uM) is a real physiological value (not an
# artificially near-zero suppressed level), so excess-above-basal is
# directly computable from its own three (excess, reduction) pairs --
# (255, 0.22), (811, 0.30), (1248, 0.34) -- without importing a second
# study's differently-defined reference. Fitting the exact Michaelis-
# Menten form this pathway's own formula uses
# (reduction = ir_max * excess/(ir_ec50 + excess)) via nonlinear
# least-squares (independently re-derived, not copied from the search
# agent's own report): ir_max ~ 0.385, ir_ec50 ~ 197 uM (residuals all
# <1 percentage point). With only 3 data points, pairwise-exact solves
# range ir_ec50 ~ 162-410 uM and ir_max ~ 0.360-0.452 -- a real,
# disclosed range, not a single precise answer.
#
# NEFA_IR_EC50_UM changed from the uncited 500.0 to 200.0 (a genuine,
# substantive correction -- 500 sits clearly OUTSIDE the fitted
# plausible range of 162-410). NEFA_IR_MAX left UNCHANGED at 0.40: the
# fit's own range (0.360-0.452) comfortably includes the prior uncited
# guess, so the new evidence corroborates rather than contradicts it --
# changing a well-corroborated number just to match a 3-point fit's own
# central estimate more precisely would overstate this fit's real
# precision (same "a coincidence isn't secretly evidence" honesty
# convention already used elsewhere in this codebase).
NEFA_IR_EC50_UM       = 200.0   # uM excess-above-basal NEFA for half-maximal IR contribution
NEFA_IR_MAX           = 0.40    # maximum p3 reduction at saturating NEFA (40%) -- corroborated, not changed
NEFA_IR_REFERENCE     = NEFA_BASAL_UM   # fixed population-average basal NEFA (350 uM)


class NEFAParameters:
    def __init__(
        self,
        nefa_basal:    float = NEFA_BASAL_UM,
        nefa_max:      float = NEFA_OBESITY_UM,
        insulin_ic50:  float = NEFA_INSULIN_IC50_PM,
        hill_n:        float = NEFA_INSULIN_N,
        half_life_min: float = NEFA_HALF_LIFE_MIN,
        ir_ec50:       float = NEFA_IR_EC50_UM,
        ir_max:        float = NEFA_IR_MAX,
    ) -> None:
        self.nefa_basal    = nefa_basal
        self.nefa_max      = nefa_max
        self.insulin_ic50  = insulin_ic50
        self.hill_n        = hill_n
        self.half_life_min = half_life_min
        self.ir_ec50       = ir_ec50
        self.ir_max        = ir_max

    @classmethod
    def for_obesity(cls) -> "NEFAParameters":
        """Elevated basal NEFA with impaired insulin suppression."""
        return cls(nefa_basal=650.0, nefa_max=950.0, insulin_ic50=150.0)

    @classmethod
    def for_exercise(cls) -> "NEFAParameters":
        """During exercise: NEFA rises dramatically (adipose lipolysis)."""
        return cls(nefa_basal=800.0, nefa_max=1800.0, insulin_ic50=50.0)


class NEFAPathway(BasePathway):
    """
    1-species plasma NEFA ODE model.

    NEFA is suppressed by insulin (antilipolysis in adipose tissue).
    At elevated NEFA, the output signal insulin_resistance_modifier
    reduces p3 (insulin action) in the coupled GlucoseInsulinPathway.
    """

    NAME    = "nefa_v1"
    VERSION = "1.0.0"

    def __init__(self, params: Optional[NEFAParameters] = None) -> None:
        self.params = params or NEFAParameters()
        self._y0 = [self.params.nefa_basal]
        self._descriptor = PathwayDescriptor(
            name=self.NAME,
            version=self.VERSION,
            maturity=PathwayMaturity.EXPERIMENTAL,
            description=(
                "1-species plasma NEFA model. "
                "NEFA suppressed by insulin (antilipolysis). "
                "Elevated NEFA reduces insulin sensitivity via Randle cycle."
            ),
            species_names=["NEFA"],
            output_signals=["nefa_uM", "insulin_resistance_modifier"],
            input_signals=["insulin_pM"],
            dependencies=[],
            tags=["lipid", "NEFA", "FFA", "Randle", "insulin-resistance"],
            literature_refs=["PMID:13990765", "PMID:8971073"],
        )

    @property
    def descriptor(self) -> PathwayDescriptor:
        return self._descriptor

    @property
    def species_names(self) -> list[str]:
        return ["NEFA"]

    @property
    def species_units(self) -> list[str]:
        return ["uM"]

    @property
    def y0(self) -> list[float]:
        return list(self._y0)

    def dydt(self, t: float, y: list[float], signals=None) -> list[float]:
        nefa = max(0.0, y[0])
        p    = self.params

        insulin_pM = 66.0
        if signals and "insulin_pM" in signals:
            insulin_pM = max(0.0, signals["insulin_pM"])

        # Insulin suppresses lipolysis -> lower NEFA target at high insulin.
        # antilipolysis is normalized against its own value at the reference
        # basal insulin (NEFA_REFERENCE_INSULIN_PM) so nefa_target equals
        # p.nefa_basal EXACTLY when insulin_pM == NEFA_REFERENCE_INSULIN_PM,
        # for any insulin_ic50/hill_n/nefa_max. Without this normalization,
        # nefa_target = nefa_max*(1-A) + nefa_basal*A is a convex combination
        # of nefa_basal and nefa_max (A = hill_activation, in [0,1)) that can
        # only equal nefa_basal if nefa_max == nefa_basal or in the
        # unreachable limit A -> 1 (insulin -> infinity) -- no finite
        # insulin_ic50 could make nefa_basal a true fixed point of the
        # un-normalized formula (unlike HepGlucPathway's egp_target, which
        # has no additive nefa_basal-like floor and so is solvable by
        # re-deriving its IC50 alone -- this pathway's convex-combination
        # structure needed the normalization itself, not just a re-solved
        # constant).
        antilipolysis        = hill_activation(insulin_pM, p.insulin_ic50, p.hill_n)
        antilipolysis_basal  = hill_activation(NEFA_REFERENCE_INSULIN_PM, p.insulin_ic50, p.hill_n)
        suppression_ratio    = antilipolysis / antilipolysis_basal
        nefa_target = max(0.0, p.nefa_max - (p.nefa_max - p.nefa_basal) * suppression_ratio)

        ke = math.log(2) / p.half_life_min
        dNEFA = ke * nefa_target - ke * nefa
        return [dNEFA]

    def output_signals(self, t: float, y: list[float]) -> dict[str, float]:
        nefa = max(0.0, y[0])
        p    = self.params

        # Randle cycle: NEFA impairs glucose oxidation → IR modifier on p3.
        # Use the fixed population-average NEFA reference (not this instance's
        # own basal) so the modifier is exactly neutral (1.0) at the NEFA level
        # the Bergman 1981 p3 calibration cohort was measured at. Using
        # instance-specific basal here would zero out real insulin resistance
        # for e.g. the obesity preset, whose elevated basal NEFA IS the
        # mechanism being modelled (same fixed-reference pattern as
        # cortisol_pathway.GR_REFERENCE).
        excess_nefa    = max(0.0, nefa - NEFA_IR_REFERENCE)
        nefa_fraction  = excess_nefa / (p.ir_ec50 + excess_nefa) if excess_nefa > 0 else 0.0
        ir_modifier    = max(0.3, 1.0 - p.ir_max * nefa_fraction)

        return {
            "nefa_uM":                      nefa,
            "insulin_resistance_modifier":   ir_modifier,
        }

    @classmethod
    def for_obesity(cls) -> "NEFAPathway":
        return cls(NEFAParameters.for_obesity())

    @classmethod
    def for_exercise(cls) -> "NEFAPathway":
        return cls(NEFAParameters.for_exercise())


def nefa_couplings(gi_name: str, nefa_name: str) -> list:
    """
    Build the CouplingSpecs connecting a NEFAPathway to
    GlucoseInsulinPathway:

        GI.insulin_pM                    -> nefa input "insulin_pM"
        nefa.insulin_resistance_modifier  -> GI input "insulin_resistance_modifier"

    Deliberately NOT wired to BMI/obesity in SimulationRunner: the fixed
    NEFA_IR_REFERENCE calibration point means NEFAPathway.for_obesity()'s
    modifier is NOT neutral at its own (elevated) basal -- composing it
    for every obese profile would double-count insulin resistance already
    modelled by biosim_physiology.modifiers.bmi_modifiers(). Use this
    coupling for scenarios where NEFA's own dynamics add real value that a
    static BMI modifier cannot -- primarily extended fasting and
    post-exercise NEFA mobilisation -- with the default (basal-neutral)
    NEFAParameters, not the obesity preset.
    """
    from .coupling import CouplingSpec

    return [
        CouplingSpec(
            source_pathway=gi_name, signal_name="insulin_pM",
            target_pathway=nefa_name, target_input="insulin_pM",
        ),
        CouplingSpec(
            source_pathway=nefa_name, signal_name="insulin_resistance_modifier",
            target_pathway=gi_name, target_input="insulin_resistance_modifier",
        ),
    ]


def couple_with_gi(gi: "GlucoseInsulinPathway", nefa: "NEFAPathway") -> "CompositePathway":
    """
    Build the CompositePathway connecting a NEFAPathway to
    GlucoseInsulinPathway (see nefa_couplings() for the signal contract).
    """
    from .coupling import CompositePathway

    return CompositePathway([gi, nefa], couplings=nefa_couplings(gi.name, nefa.name))
