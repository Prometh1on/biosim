# SAVE TO: packages\ode-engine\src\biosim_ode\alcohol_hepatic_pathway.py
"""
Alcohol hepatic glucose production (EGP) suppression pathway.

Found and built via a fresh Track 1 scoping fan-out (2026-08-18) that
re-checked `drug_interactions.py`'s own explicitly-named uncited
placeholder (`alcohol_modifiers()`, "effect magnitudes are a placeholder
estimate pending a confirmed source"). The scoping pass found the
placeholder wasn't just uncited -- it was actively WRONG, contradicted
by both real primary sources located: alcohol's genuine, measured effect
is a TRANSIENT REDUCTION IN HEPATIC GLUCOSE PRODUCTION (a flux), not a
shift in the fasting glucose CONCENTRATION set-point (Gb), which both
sources explicitly rule out.

Literature, both live-verified via NCBI efetch before implementation:

    Siler SQ, Neese RA, Christiansen MP, Hellerstein MK. "The inhibition
    of gluconeogenesis following alcohol in humans." Am J Physiol.
    1998;275(5):E897-907. PMID:9815011. Normal men, overnight-fasted,
    48g oral alcohol vs. placebo, stable-isotope tracers. Confirmed
    verbatim: "Hepatic glucose output decreased 12% after alcohol
    consumption, from 2.03 +/- 0.21 to 1.79 +/- 0.21 mg.kg-1.min-1
    (P<0.05 vs. baseline), but did not change following the placebo,"
    and, critically for the fix below: "...a smaller effect on glucose
    production and NO EFFECT ON PLASMA GLUCOSE CONCENTRATIONS."

    Shelmet JJ, Reichard GA, Skutches CL, Hoeldtke RD, Owen OE, Boden G.
    "Ethanol causes acute inhibition of carbohydrate, fat, and protein
    oxidation and insulin resistance." J Clin Invest. 1988;81(4):1137-45.
    PMID:3280601. IV ethanol (0.75 g/kg), six normal males. Corroborates
    a real hepatic/peripheral effect (basal glucose appearance -30%,
    clamp glucose disposal -36%) but ALSO explicitly states "had no
    effect on glucose tolerance" -- independent confirmation that the
    fasting/ambient glucose CONCENTRATION itself is not what moves.

    ROUTE-MATCHED ANCHOR CHOICE: Siler (oral, 48g) is used for the
    calibration below, not Shelmet (IV, ~52g-equivalent), because this
    platform's `DrugDose`/`ScenarioDefinition` alcohol events are
    oral doses (`drug_scenarios.alcohol_before_meal()` converts grams
    directly to `dose_mg`) -- the same "which regime does this function
    actually represent" resolution already used for `metformin_
    modifiers()`'s own p3/p1 removal and for `NutrientSecretionPathway`'s
    Krezowski-vs-Nuttall choice. This does not mean Shelmet's IV numbers
    are wrong -- route (oral vs. IV) plausibly explains part of the
    magnitude difference (first-pass hepatic alcohol metabolism differs
    by route), so Siler is simply the better-matched anchor, not the
    only real one.

Architecturally the SAME shape as `nutrient_secretion_pathway.py`'s own
NutrientSecretionPathway/ExogenousInsulinPathway: STATELESS (zero
species -- a closed-form analytic curve, no ODE integration state
needed) and TRANSIENT (exactly zero before/without alcohol, decaying
back to exactly zero afterward). Injects into the ALREADY-EXISTING
`hepatic_egp_delta` rule-11 signal (an additive delta -- see
coupling.py's ADDITIVE_SIGNAL_KEYS) rather than inventing a new signal
key, since this is exactly the same physical quantity HepGlucPathway
already produces (excess hepatic glucose output above basal) --
`hepatic_egp_delta` correctly SUMS contributions from multiple sources
(confirmed via coupling.py's `_combine_signal_values()`), so this
pathway composes safely alongside T2D/metformin/hyperthyroidism's own
HepGlucPathway-based coupling if both are active in the same scenario.

Deliberately NOT modeled as an `EGPParameters`/`HepGlucPathway` preset
(unlike T2D/hyperthyroidism/metformin): those all represent CHRONIC,
steady-state hepatic states (a disease, or a daily medication) where
the elevated/suppressed EGP basal rate is itself the analytic fixed
point of the whole simulation. Alcohol is fundamentally an ACUTE,
self-resolving perturbation -- exactly the same category distinction
this platform already draws between `ExogenousInsulinPathway`
(acute, time-varying subcutaneous depot) and the chronic modifier-based
drugs. Forcing alcohol into `EGPParameters.for_metformin()`'s own
chronic-preset pattern would make the EGP suppression CONSTANT for the
whole simulation, which is not what either primary source describes
(alcohol clears; the effect resolves).

Timing shape: UNANCHORED, disclosed assumption (same category as
NutrientSecretionPathway's own disclosed timing assumption) -- neither
Siler nor Shelmet reports a time-resolved suppression curve, only a
single value over a multi-hour window (Siler: "in the 5 h after alcohol
ingestion"). A symmetric peak-then-decay curve (peak ~2h post-
ingestion, resolving by ~6h) is assumed, carried over from the prior
(uncited) alcohol_modifiers() implementation's own timing shape, which
is broadly consistent with typical alcohol absorption/clearance
timescales and Siler's own 5h measurement window -- but this specific
curve SHAPE is not itself derived from either citation, only the peak
AMPLITUDE is (Siler's -12% HGP figure). Flagged here, not hidden.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Optional

from .base_pathway import BasePathway, PathwayDescriptor, PathwayMaturity
from .hepatic import EGP_BASAL_MM_MIN

if TYPE_CHECKING:
    from .coupling import CompositePathway
    from .pathway import GlucoseInsulinPathway


# Siler et al. 1998 (PMID:9815011): "Hepatic glucose output decreased
# 12%... after alcohol consumption" at a 48g oral dose -- a directly-
# reported fractional reduction, applied to this platform's own
# EGP_BASAL_MM_MIN (the same "fractional change applied to this
# platform's own basal rate" convention already used for metformin's
# Gb factor and T2D/hyperthyroidism's egp_basal anchors).
ALCOHOL_EGP_REDUCTION_FRACTION_AT_TESTED_DOSE = 0.12
ALCOHOL_TESTED_DOSE_G = 48.0   # Siler's own tested dose -- do not extrapolate past this

# UNANCHORED timing shape (see module docstring) -- carried over from
# the prior (uncited) alcohol_modifiers() implementation's own assumed
# curve, not derived from either citation.
ALCOHOL_PEAK_TIME_MIN = 120.0     # 2h
ALCOHOL_DECAY_TAU_MIN = 90.0      # 1.5h
ALCOHOL_RESOLVES_BY_MIN = 360.0   # 6h -- hard cutoff, matching the prior implementation


class AlcoholHepaticParameters:
    """Alcohol hepatic-EGP-suppression parameters for a single dose event."""

    def __init__(
        self,
        dose_g: float,
        dose_time_min: float = 0.0,
        egp_basal_mM_min: float = EGP_BASAL_MM_MIN,
    ) -> None:
        if not math.isfinite(dose_g):
            raise ValueError(f"dose_g must be finite, got {dose_g}")
        if dose_g < 0:
            raise ValueError(f"dose_g cannot be negative, got {dose_g}")
        if not math.isfinite(dose_time_min):
            raise ValueError(f"dose_time_min must be finite, got {dose_time_min}")
        if not math.isfinite(egp_basal_mM_min) or egp_basal_mM_min < 0:
            raise ValueError(f"egp_basal_mM_min must be finite and non-negative, got {egp_basal_mM_min}")
        self.dose_g = dose_g
        self.dose_time_min = dose_time_min
        self.egp_basal_mM_min = egp_basal_mM_min

    @property
    def peak_delta_mM_per_min(self) -> float:
        """
        The maximum (t=peak, dose at or above the tested cap) EGP
        suppression this dose would produce -- a NEGATIVE value (excess
        below basal), dose-scaled linearly up to Siler's own tested dose
        and capped there (no data exists above 48g to extrapolate from).
        """
        dose_factor = min(self.dose_g / ALCOHOL_TESTED_DOSE_G, 1.0)
        return -ALCOHOL_EGP_REDUCTION_FRACTION_AT_TESTED_DOSE * self.egp_basal_mM_min * dose_factor

    def delta_mM_per_min(self, t_minutes: float) -> float:
        """
        The hepatic_egp_delta contribution from THIS dose alone, at
        t_minutes. Exactly 0.0 before the dose and beyond the resolution
        window; peaks (most negative) at ALCOHOL_PEAK_TIME_MIN after
        dosing.
        """
        t_since = t_minutes - self.dose_time_min
        if t_since <= 0.0 or t_since > ALCOHOL_RESOLVES_BY_MIN or self.dose_g <= 0.0:
            return 0.0
        decay = math.exp(-abs(t_since - ALCOHOL_PEAK_TIME_MIN) / ALCOHOL_DECAY_TAU_MIN)
        return self.peak_delta_mM_per_min * decay


class AlcoholHepaticPathway(BasePathway):
    """
    Stateless (zero-species) pathway injecting a time-varying, transient
    hepatic_egp_delta signal into GlucoseInsulinPathway (directly, or
    additively alongside HepGlucPathway's own contribution -- see module
    docstring), summing the contribution of every alcohol dose event.
    """

    NAME    = "alcohol_hepatic_v1"
    VERSION = "1.0.0"

    def __init__(self, doses: Optional[list[AlcoholHepaticParameters]] = None) -> None:
        self._doses = list(doses or [])

    @property
    def descriptor(self) -> PathwayDescriptor:
        return PathwayDescriptor(
            name=self.NAME,
            version=self.VERSION,
            maturity=PathwayMaturity.EXPERIMENTAL,
            description=(
                "Stateless, transient alcohol-driven hepatic glucose "
                "production suppression -- a closed-form peak-then-decay "
                "curve per dose event, additive to hepatic_egp_delta."
            ),
            species_names=[],
            output_signals=["hepatic_egp_delta"],
            input_signals=[],
            dependencies=[],
            tags=["hepatic", "alcohol", "EGP", "time-varying"],
            literature_refs=["PMID:9815011", "PMID:3280601"],
        )

    @property
    def species_names(self) -> list[str]:
        return []

    @property
    def species_units(self) -> list[str]:
        return []

    @property
    def y0(self) -> list[float]:
        return []

    def dydt(
        self,
        t: float,
        y: list[float],
        signals: Optional[dict] = None,
    ) -> list[float]:
        return []

    def output_signals(self, t: float, y: list[float]) -> dict[str, float]:
        total = sum(d.delta_mM_per_min(t) for d in self._doses)
        # Aggregate plausibility floor (mandatory security-reviewer
        # finding): each individual dose's own peak contribution is
        # already capped via ALCOHOL_TESTED_DOSE_G's dose_factor, but
        # summing MANY doses has no aggregate cap -- a REST-legal
        # scenario can submit up to 50 drug_events (schemas.py's own
        # ScenarioIn.drug_events max_length=50). Verified directly
        # against the real composed ODE: 50 near-simultaneous
        # alcohol doses (each at/above the per-dose saturation point)
        # drove simulated glucose to ~1.65 mM (severe-hypoglycemia
        # range) despite every individual dose being correctly bounded
        # -- the raw additive hepatic_egp_delta signal has no floor
        # tying the AGGREGATE to the physical fact that hepatic glucose
        # output cannot fall below zero (unlike HepGlucPathway's own
        # Hill-inhibition EGP, which is inherently bounded in
        # [0, egp_max]). Floored here at -egp_basal_mM_min (the largest
        # basal reference among the composed doses, defaulting to
        # EGP_BASAL_MM_MIN) -- the aggregate can never imply more than
        # 100% suppression of the reference basal EGP at any instant.
        if self._doses:
            max_egp_basal = max(d.egp_basal_mM_min for d in self._doses)
            total = max(total, -max_egp_basal)
        return {"hepatic_egp_delta": total}


def alcohol_hepatic_couplings(gi_name: str, alcohol_name: str) -> list:
    """
    Build the CouplingSpec connecting an AlcoholHepaticPathway to
    GlucoseInsulinPathway. Single signal, maps straight through -- no
    transform needed, since hepatic_egp_delta is already expressed in
    GI's own signal vocabulary (rule 11) and this pathway's own output
    IS the delta directly (no basal to subtract, unlike HepGlucPathway).
    """
    from .coupling import CouplingSpec

    return [
        CouplingSpec(
            source_pathway=alcohol_name, signal_name="hepatic_egp_delta",
            target_pathway=gi_name, target_input="hepatic_egp_delta",
        ),
    ]


def couple_with_gi(gi: "GlucoseInsulinPathway", alcohol: "AlcoholHepaticPathway") -> "CompositePathway":
    """
    Build the CompositePathway connecting an AlcoholHepaticPathway to
    GlucoseInsulinPathway (see alcohol_hepatic_couplings() for the
    signal contract).
    """
    from .coupling import CompositePathway

    return CompositePathway(
        [gi, alcohol], couplings=alcohol_hepatic_couplings(gi.name, alcohol.name),
    )
