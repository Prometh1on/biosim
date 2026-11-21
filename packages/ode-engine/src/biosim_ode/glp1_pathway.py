# SAVE TO: packages\ode-engine\src\biosim_ode\glp1_pathway.py
"""
GLP-1 (Glucagon-like peptide-1) incretin pathway.

GLP-1 is the gut hormone secreted by L-cells in response to nutrient
ingestion. It is the most clinically important incretin because:
    1. It amplifies glucose-dependent insulin secretion (the incretin effect)
    2. It suppresses glucagon secretion
    3. It slows gastric emptying (reduces Ra peak)
    4. It reduces appetite centrally

GLP-1 agonists (semaglutide/Ozempic, liraglutide/Victoza, tirzepatide)
are the dominant new drug class for T2D and obesity. This pathway allows
the platform to simulate their effects.

State vector:
    y[0]  GLP1  [pM]  plasma GLP-1 concentration

Physiology:
    Basal GLP-1: ~7 pM fasting (Orskov et al. 1994, PMID:8138058, n=6)
    Peak GLP-1:  ~41 pM postprandial (same source; peak reported at
                 90 min post-breakfast -- this platform's own reference
                 meal peaks the underlying glucose signal, and hence
                 GLP1's own fast-tracking response, closer to t=25min;
                 GLP1_MEAL_RESPONSE is calibrated to the reported PEAK
                 MAGNITUDE only, not this specific timing -- see the
                 constant's own inline comment above)
    Half-life:   ~2 min (rapidly degraded by DPP-4)

GLP-1 agonist simulation:
    Use GLP1Parameters.for_semaglutide() or for_liraglutide() to model
    pharmacological GLP-1 receptor agonism. These extend the effective
    half-life to hours/days and amplify the receptor effect.

Output signals (for coupling to GlucoseInsulinPathway):
    "glp1_pM"              plasma GLP-1
    "incretin_effect"      fractional amplification of insulin secretion (0-1)
    "gastric_emptying_modifier"  < 1 = slower stomach emptying (scales kabs
                                 for a single-exponential meal, k_emptying for
                                 a two-compartment meal)

Literature:
    Drucker & Nauck 2006, Lancet. PMID:17098089
    Nauck et al. 1993, Diabetologia. PMID:8405741
    Holst 2007, Physiol Rev. PMID:17928588
    Orskov et al. 1994, Diabetes. PMID:8138058 (endogenous meal-response
        basal/peak GLP-1 concentration anchor -- GLP1_BASAL_PM,
        GLP1_MEAL_RESPONSE)
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Optional

from .base_pathway import BasePathway, PathwayDescriptor, PathwayMaturity
from .kinetics import hill_activation, approach_to_target

if TYPE_CHECKING:
    from .coupling import CompositePathway
    from .pathway import GlucoseInsulinPathway


# Default GLP-1 parameters
#
# GLP1_BASAL_PM and GLP1_MEAL_RESPONSE (2026-08-12 deep-search plan, item
# A1): anchored to Orskov C, Rabenhoj L, Wettergren A, Kofod H, Holst JJ.
# "Tissue and plasma concentrations of amidated and glycine-extended
# glucagon-like peptide I in humans." Diabetes. 1994;43(4):535-9.
# PMID:8138058 (live-fetched and independently re-verified twice this
# session, not trusted from a single pass). Directly quoted: "The mean
# concentrations of amidated GLP-I ... in fasting plasma were 7 +/- 1 ...
# pM (n=6). In response to a breakfast meal, the concentration of
# amidated GLP-I rose significantly amounting to 41 +/- 5 pM 90 min
# after the meal ingestion." Basal set to the paper's own 7.0 pM.
# GLP1_MEAL_RESPONSE was previously an uncited guess (0.035) that this
# file's own earlier history had already flagged as producing a
# postprandial GLP1 peak only ~1.05x basal -- far below any real
# reported range; numerically re-derived (bisection against the real
# composed GlucoseInsulinPathway+GLP1Pathway ODE, standard 75g
# reference meal) to reproduce the paper's own 41 pM peak magnitude.
# Deliberately NOT claimed to match the paper's own 90-minute timing --
# this platform's standard 75g reference meal peaks glucose (and
# therefore GLP1, whose ~2min half-life makes it track glucose almost
# instantaneously) around t=25min, not 90min; the paper does not report
# its own breakfast meal's glucose kinetics, so a joint peak-height AND
# peak-time fit (the technique used elsewhere in this codebase, e.g.
# acarbose) isn't possible here -- only the magnitude is anchored.
# Also unresolved: the paper's own abstract (truncated at 250 words,
# pre-PMC era, no full text accessible) never explicitly states the
# n=6 subjects were healthy/normal-glucose-tolerant -- no disease-state
# MeSH heading is attached to the record, which is circumstantial
# support, not direct confirmation.
GLP1_BASAL_PM        = 7.0     # pM fasting GLP-1
GLP1_HALF_LIFE_MIN   = 2.0     # min (DPP-4 degradation)
GLP1_MEAL_RESPONSE   = 4.0     # pM per mM glucose above threshold
GLP1_EC50_SECRETION  = 6.0     # mM glucose for half-max GLP-1 secretion -- still an uncited threshold assumption (see deep-search plan A1 note in this module's docstring); Orskov 1994 anchors only the basal/peak GLP-1 concentrations, not a glucose threshold.
GLP1_EC50_INSULIN    = 15.0    # pM GLP-1 for half-max incretin effect
GLP1_MAX_INCRETIN    = 0.60    # maximum fractional amplification of insulin secretion
GLP1_EMPTYING_IC50   = 20.0    # pM GLP-1 for 50% gastric emptying reduction
GLP1_EMPTYING_MAX    = 0.50    # maximum fractional reduction in gastric emptying rate


class GLP1Parameters:
    """Parameters for the GLP-1 incretin pathway."""

    def __init__(
        self,
        glp1_basal:        float = GLP1_BASAL_PM,
        half_life_min:     float = GLP1_HALF_LIFE_MIN,
        meal_response:     float = GLP1_MEAL_RESPONSE,
        ec50_secretion:    float = GLP1_EC50_SECRETION,
        ec50_insulin:      float = GLP1_EC50_INSULIN,
        max_incretin:      float = GLP1_MAX_INCRETIN,
        ec50_emptying:     float = GLP1_EMPTYING_IC50,
        max_emptying_reduction: float = GLP1_EMPTYING_MAX,
        agonist_scale:     float = 1.0,   # > 1.0 for pharmacological agonists
    ) -> None:
        self.glp1_basal    = glp1_basal
        self.half_life_min = half_life_min
        self.meal_response = meal_response
        self.ec50_secretion = ec50_secretion
        self.ec50_insulin  = ec50_insulin
        self.max_incretin  = max_incretin
        self.ec50_emptying = ec50_emptying
        self.max_emptying_reduction = max_emptying_reduction
        self.agonist_scale = agonist_scale

    @classmethod
    def for_semaglutide(cls, dose_mg: float = 1.0) -> "GLP1Parameters":
        """
        Semaglutide (Ozempic/Wegovy) — weekly GLP-1 receptor agonist.

        Semaglutide has >94% albumin binding and a half-life of ~1 week.
        Modelled as sustained receptor occupancy proportional to dose.
        At 1mg/week steady state: ~80% receptor occupancy.

        Reference: Lau et al. 2015, J Med Chem. PMID:26308095
        """
        # Dose-response: 0.5mg=60%, 1mg=80%, 2mg=90% receptor occupancy
        occupancy = 1.0 - math.exp(-0.8 * dose_mg)
        agonist_scale = 1.0 + occupancy * 8.0   # up to 9x basal GLP-1 effect
        return cls(
            glp1_basal=GLP1_BASAL_PM * (1.0 + occupancy * 5.0),
            half_life_min=10080.0,   # 1 week in minutes
            max_incretin=min(0.95, GLP1_MAX_INCRETIN + occupancy * 0.35),
            max_emptying_reduction=min(0.70, GLP1_EMPTYING_MAX + occupancy * 0.20),
            agonist_scale=agonist_scale,
        )

    @classmethod
    def for_liraglutide(cls, dose_mg: float = 1.2) -> "GLP1Parameters":
        """
        Liraglutide (Victoza/Saxenda) — daily GLP-1 receptor agonist.

        Half-life ~13 hours. Standard T2D dose 1.2-1.8 mg/day.
        Reference: Buse et al. 2009, Lancet (LEAD-6). PMID:19515413
        """
        occupancy = 1.0 - math.exp(-0.6 * dose_mg)
        return cls(
            glp1_basal=GLP1_BASAL_PM * (1.0 + occupancy * 3.0),
            half_life_min=780.0,   # 13 hours
            max_incretin=min(0.85, GLP1_MAX_INCRETIN + occupancy * 0.25),
            max_emptying_reduction=min(0.60, GLP1_EMPTYING_MAX + occupancy * 0.15),
            agonist_scale=1.0 + occupancy * 4.0,
        )

    @classmethod
    def dpp4_inhibited(cls) -> "GLP1Parameters":
        """
        DPP-4 inhibitor effect (sitagliptin/saxagliptin).
        Prevents GLP-1 degradation → doubles effective GLP-1 half-life.
        Reference: Deacon 2011, Diabetes Obes Metab. PMID:21114598

        The half-life doubling is the only mechanism modelled here (this
        used to also set a `glp1_peak` field intended to represent direct
        peak elevation, but dydt()/output_signals() never read it -- it
        was a pure no-op, removed rather than wired in since no primary
        source anchors a specific peak-elevation magnitude).
        """
        return cls(half_life_min=GLP1_HALF_LIFE_MIN * 2.0)


class GLP1Pathway(BasePathway):
    """
    1-species GLP-1 incretin ODE model.

    Couples to GlucoseInsulinPathway via:
        - incretin_effect: amplifies secretion_scale during meal
        - gastric_emptying_modifier: slows stomach emptying (flattens glucose peak)
    """

    NAME    = "glp1_v1"
    VERSION = "1.0.0"

    def __init__(
        self,
        params: Optional[GLP1Parameters] = None,
        basal_glucose_mM: float = 5.0,
    ) -> None:
        self.params = params or GLP1Parameters()
        self._basal_glucose_mM = basal_glucose_mM
        # y0 must be the analytic fixed point of dydt() at the ACTUAL
        # basal glucose this pathway will be composed with -- not just
        # glp1_basal * agonist_scale, which is only correct when
        # glucose_stim=0 at that basal glucose (i.e. basal_glucose_mM <=
        # ec50_secretion). dydt() scales the ENTIRE secretion term (basal
        # replacement + meal response) by agonist_scale, so the true
        # steady state at basal glucose G0 is
        #   agonist_scale * (glp1_basal + meal_response/ke * max(0, G0 - ec50_secretion))
        # `basal_glucose_mM` defaults to 5.0 (this pathway's own
        # historical default and dydt()'s own internal fallback when no
        # `glucose_mM` signal is present) -- for every basal_glucose_mM
        # <= ec50_secretion (GLP1_EC50_SECRETION=6.0 by default,
        # including every healthy-Gb scenario), this reduces to exactly
        # the prior formula, glp1_basal * agonist_scale.
        #
        # 2026-08-12 deep-search plan follow-up (mandatory
        # security-reviewer finding on the Orskov 1994 meal-response
        # recalibration, CRITICAL): before this parameter existed, this
        # pathway's y0 ALWAYS assumed glucose_stim=0 regardless of what
        # GI it was actually composed with -- for a disease-elevated Gb
        # preset (T2D Gb~6.5mM, GCK-MODY Gb~7.17mM, both already shipped
        # elsewhere in this codebase), glucose_stim_basal is genuinely
        # nonzero, so the pathway's own "basal steady state" was not
        # actually a fixed point of its own ODE -- a real, reproducible
        # drift (confirmed directly: T2D+semaglutide 1mg, standalone
        # fasted, drifted 142.0 -> 3352.0 pM over a standard 5h window
        # before this fix). This was a dormant, minor issue at the
        # prior, much smaller GLP1_MEAL_RESPONSE=0.035 (a ~17% 5h drift)
        # -- the Orskov-anchored recalibration's ~114x larger
        # meal_response amplified it into a severe one (~2260%),
        # exactly the KNOWN BUGS FIXED #10/#11/#12/#14 bug class (a
        # preset's own y0 not actually being its ODE's fixed point), now
        # fixed the same way: make y0 a genuine function of the real
        # basal state it will run at, not a hardcoded assumption.
        ke = math.log(2) / max(self.params.half_life_min, 0.1)
        glucose_stim_basal = max(0.0, basal_glucose_mM - self.params.ec50_secretion)
        self._y0 = [
            self.params.agonist_scale
            * (self.params.glp1_basal + self.params.meal_response * glucose_stim_basal / ke)
        ]
        self._descriptor = PathwayDescriptor(
            name=self.NAME,
            version=self.VERSION,
            maturity=PathwayMaturity.BETA,
            description=(
                "1-species GLP-1 incretin model. "
                "Nutrient-stimulated GLP-1 amplifies insulin secretion "
                "and slows gastric emptying. "
                "Supports simulation of GLP-1 agonists (semaglutide, liraglutide) "
                "and DPP-4 inhibitors."
            ),
            species_names=["GLP1"],
            output_signals=["glp1_pM", "incretin_effect", "gastric_emptying_modifier"],
            input_signals=["glucose_mM"],
            dependencies=[],
            tags=["incretin", "GLP-1", "gut", "semaglutide", "liraglutide", "DPP4"],
            literature_refs=["PMID:17098089", "PMID:8405741", "PMID:17928588", "PMID:8138058"],
        )

    @property
    def descriptor(self) -> PathwayDescriptor:
        return self._descriptor

    @property
    def species_names(self) -> list[str]:
        return ["GLP1"]

    @property
    def species_units(self) -> list[str]:
        return ["pM"]

    @property
    def y0(self) -> list[float]:
        return list(self._y0)

    def dydt(self, t: float, y: list[float], signals=None) -> list[float]:
        glp1 = max(0.0, y[0])
        p    = self.params
        ke   = math.log(2) / max(p.half_life_min, 0.1)

        # Glucose-stimulated GLP-1 secretion from L-cells. Falls back to
        # this instance's own configured basal_glucose_mM (not a bare
        # 5.0) when no glucose_mM signal is present -- keeps
        # validate_y0()'s signals=None check (BasePathway) consistent
        # with whatever basal glucose y0 was actually derived for.
        glucose_mM = self._basal_glucose_mM
        if signals and "glucose_mM" in signals:
            glucose_mM = max(0.0, signals["glucose_mM"])

        # GLP-1 secretion rises with glucose above threshold
        glucose_stim = max(0.0, glucose_mM - p.ec50_secretion)
        secretion    = ke * p.glp1_basal + p.meal_response * glucose_stim

        # Scale by agonist (pharmacological GLP-1 receptor agonist effect)
        secretion *= p.agonist_scale

        dGLP1 = secretion - ke * glp1
        return [dGLP1]

    def output_signals(self, t: float, y: list[float]) -> dict[str, float]:
        glp1 = max(0.0, y[0])
        p    = self.params

        # Incretin effect: fractional amplification of insulin secretion
        incretin = p.max_incretin * hill_activation(glp1, p.ec50_insulin, 1.5)

        # Gastric emptying modifier: slows absorption at high GLP-1
        emptying_reduction = p.max_emptying_reduction * hill_activation(
            glp1, p.ec50_emptying, 1.2
        )
        gastric_modifier = max(0.30, 1.0 - emptying_reduction)

        return {
            "glp1_pM":                 glp1,
            "incretin_effect":         incretin,
            "gastric_emptying_modifier": gastric_modifier,
        }

    @classmethod
    def semaglutide(cls, dose_mg: float = 1.0, basal_glucose_mM: float = 5.0) -> "GLP1Pathway":
        return cls(GLP1Parameters.for_semaglutide(dose_mg), basal_glucose_mM=basal_glucose_mM)

    @classmethod
    def liraglutide(cls, dose_mg: float = 1.2, basal_glucose_mM: float = 5.0) -> "GLP1Pathway":
        return cls(GLP1Parameters.for_liraglutide(dose_mg), basal_glucose_mM=basal_glucose_mM)

    @classmethod
    def dpp4_inhibited(cls, basal_glucose_mM: float = 5.0) -> "GLP1Pathway":
        return cls(GLP1Parameters.dpp4_inhibited(), basal_glucose_mM=basal_glucose_mM)


def glp1_couplings(gi_name: str, glp1_name: str) -> list:
    """
    Build the CouplingSpecs connecting a GLP1Pathway to GlucoseInsulinPathway:

        GI.glucose_mM               -> glp1 input "glucose_mM"
                                        (glucose-stimulated GLP-1 secretion)
        glp1.incretin_effect        -> GI input "incretin_effect"
                                        (amplifies insulin secretion)
        glp1.gastric_emptying_modifier -> GI input "gastric_emptying_modifier"
                                        (slows meal glucose appearance)

    Exposed standalone (not just inside couple_with_gi) so SimulationRunner
    can flat-compose GI with multiple simultaneously-active extensions
    (e.g. GLP-1 + cortisol) in a single CompositePathway -- CouplingSpec
    injection only resolves direct siblings within one CompositePathway,
    so nesting composites would silently break cross-extension signals.
    """
    from .coupling import CouplingSpec

    return [
        CouplingSpec(
            source_pathway=gi_name, signal_name="glucose_mM",
            target_pathway=glp1_name, target_input="glucose_mM",
        ),
        CouplingSpec(
            source_pathway=glp1_name, signal_name="incretin_effect",
            target_pathway=gi_name, target_input="incretin_effect",
        ),
        CouplingSpec(
            source_pathway=glp1_name, signal_name="gastric_emptying_modifier",
            target_pathway=gi_name, target_input="gastric_emptying_modifier",
        ),
    ]


def couple_with_gi(gi: "GlucoseInsulinPathway", glp1: "GLP1Pathway") -> "CompositePathway":
    """
    Build the CompositePathway that actually connects a GLP1Pathway to
    GlucoseInsulinPathway (see glp1_couplings() for the signal contract).

    Without this, GLP1Pathway.semaglutide()/liraglutide()/dpp4_inhibited()
    integrate on their own but never influence blood glucose or insulin --
    the ODE-coupled effect requires running the returned composite through
    ODESystem instead of a bare GlucoseInsulinPathway.

    IMPORTANT: `glp1` should already have been constructed with
    `basal_glucose_mM=gi.params.Gb` if `gi`'s own basal glucose differs
    from GLP1Pathway's default (5.0mM, e.g. a T2D/GCK-MODY disease
    preset) -- this function does not re-derive or check that for you,
    since GLP1Pathway's own y0 (its ODE's fixed point) is fixed at
    construction time, before this function ever sees it. Getting this
    wrong reproduces the fixed-point drift this same module's own
    GLP1Pathway.__init__ docstring documents. SimulationRunner._build_
    pathway()/._detect_glp1() already do this correctly for every
    production scenario.
    """
    from .coupling import CompositePathway

    return CompositePathway([gi, glp1], couplings=glp1_couplings(gi.name, glp1.name))
