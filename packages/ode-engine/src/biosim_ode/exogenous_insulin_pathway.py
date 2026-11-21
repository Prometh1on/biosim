# SAVE TO: packages\ode-engine\src\biosim_ode\exogenous_insulin_pathway.py
"""
Exogenous insulin pathway -- injected (subcutaneous) regular human insulin
appearing in plasma. The first mechanism this platform models that raises
plasma insulin from an EXTERNAL source, rather than via endogenous
secretion (delta_secretion/first_phase_secretion) or a secretagogue drug
(glucose_independent_secretion, sulfonylureas) -- the prerequisite for a
future Type 1 Diabetes disease state (MASTER_PLAN 1.10), which needs
exogenous insulin to be simulatable at all (T1D patients have little/no
endogenous secretion and rely entirely on injected insulin).

Architecturally distinct from every existing extension pathway in one
respect and similar in another:
    - Like ExercisePathway, this is STATELESS (zero species) -- the whole
      absorption/appearance curve has a closed-form analytic solution (see
      below), so no new ODE integration state is needed at all.
    - Like sulfonylurea's glucose_independent_secretion, its output is
      additive to dI -- but UNLIKE that field (a constant floor that
      shifts the analytic insulin fixed point), this signal is genuinely
      TRANSIENT: exactly zero before/without a dose, rising to a peak,
      then decaying back to exactly zero. It never changes the analytic
      basal steady state (GlucoseInsulinPathway().y0) at all, for any
      dose -- a strictly easier fixed-point situation than sulfonylurea's
      own (which needed a real re-derivation of basal_state()).

Mechanism and design (verified numerically before picking a structure --
see CLAUDE.md's exogenous-insulin milestone entry for the full trail):

    A SINGLE subcutaneous depot compartment feeding directly into this
    model's own plasma insulin state (which already has its own real
    elimination rate, ke=0.093 min-1, Sherwin 1974) CANNOT reproduce a
    realistic time-to-peak-plasma-insulin (Tmax) at all: solving for the
    single-depot absorption rate constant that would give Howey et al.
    1994's observed Tmax=101min (10U regular insulin, healthy subjects)
    requires a depot half-life of ~20 HOURS, which in turn makes the
    total duration of action unbounded (>900 min, never returning to
    baseline within any clinically real time frame) -- physiologically
    absurd for regular insulin's real ~5-8h duration of action. This is
    the classic "flip-flop kinetics" problem in the real PK literature
    (confirmed independently by Potocka et al. 2011, PMID:20940337: "A
    2-compartment pharmacokinetic model with... 2 sequential first-order
    absorption processes... was developed" for exactly this reason).

    TWO SEQUENTIAL depot compartments (S1 -> S2 -> plasma) resolve this.
    Using the SAME rate constant ka for both transitions (S1->S2 and
    S2->plasma) is a special, exactly-solvable case (a repeated
    eigenvalue) with a closed-form analytic solution -- no ODE state
    needed:
        S1(t) = S1_0 * exp(-ka*t)
        S2(t) = S1_0 * ka * t * exp(-ka*t)      (verified by direct
                                                   differentiation: this
                                                   satisfies dS2/dt =
                                                   ka*(S1-S2) exactly)
        exogenous_insulin_delta(t) = ka * S2(t)
                                    = ka^2 * S1_0 * t * exp(-ka*t)
    This is a Bateman-function shape (the same family used for two-
    compartment drug PK generally) -- confirmed numerically to match a
    genuine two-state ODE integration to ~1e-8 relative error before
    being adopted as the sole implementation (no ODE state = no solver
    stiffness risk from this pathway, and no risk of the closed-form
    diverging from "the real dynamics" since there IS no separate real
    dynamics here -- this pathway's whole model IS the closed form).

Calibration, disclosed as a genuine engineering back-solve, not a
directly-cited pair of numbers (see INSULIN_REGULAR_KA_PER_MIN and
INSULIN_REGULAR_PM_PER_UNIT's own docstrings immediately below for the
exact reasoning and citation).

Literature:
    Howey DC, Bowsher RR, Brunelle RL, Woodworth JR. "[Lys(B28),
    Pro(B29)]-human insulin. A rapidly absorbed analogue of human
    insulin." Diabetes. 1994;43(3):396-402. PMID:8314011 (verified live
    via NCBI efetch this session) -- n=10 healthy men, 10-U SC dose
    during a euglycaemic clamp, comparing four treatments in the SAME
    study/protocol: "Serum insulin concentrations peaked more than two
    times higher (Cmax, 698 vs. 308 pM, A vs. C) and in less than half
    the time (Tmax, 42 vs. 101 min, A vs. C)" -- arm A = zinc-free
    LYSPRO (rapid-acting insulin lispro, see INSULIN_LISPRO_KA_PER_MIN/
    INSULIN_LISPRO_PM_PER_UNIT below), arm C = regular human insulin
    (see INSULIN_REGULAR_KA_PER_MIN/INSULIN_REGULAR_PM_PER_UNIT above).
    Both analogs' calibrations come from this SAME within-study
    comparison, not two separately-sourced numbers.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Optional

from .base_pathway import BasePathway, PathwayDescriptor, PathwayMaturity

if TYPE_CHECKING:
    from .coupling import CompositePathway
    from .pathway import GlucoseInsulinPathway


# Two-sequential-equal-rate-constant depot absorption rate, calibrated
# (NOT itself a directly-cited literature number) so that, combined with
# this model's own existing plasma insulin elimination rate (ke=0.093
# min-1, Sherwin 1974 -- see PathwayParameters.ke), the resulting Tmax
# matches Howey et al. 1994's observed 101 min for a 10U regular-insulin
# dose. Kobayashi et al. 1983 (PMID:6339306) independently fits a
# ONE-compartment SC-insulin model in diabetic (IDDM/NIDDM) patients
# (ka=0.026/min, ke=0.013/min) -- that ke is NOT reused here: it is the
# APPARENT elimination rate of that paper's own one-compartment fit
# (dominated by the slow absorption phase itself, "flip-flop kinetics"),
# not this model's independently-anchored TRUE plasma clearance rate
# (0.093/min) -- reusing Kobayashi's ka directly against THIS model's
# much faster ke would give Tmax~19min (verified directly, far too
# early), which is exactly why a fresh two-compartment calibration
# against THIS model's own ke was needed rather than reusing either
# paper's fitted rate constants wholesale.
INSULIN_REGULAR_KA_PER_MIN = 0.0113

# Back-solved (not derived from an independently-cited insulin
# distribution volume, which no source found this session anchors --
# see CLAUDE.md's exogenous-insulin milestone entry) so that, at
# INSULIN_REGULAR_KA_PER_MIN above, a 10-Unit dose reproduces Howey et
# al. 1994's observed peak EXCESS above baseline (308 pM) exactly. Same
# "isolated, disclosed engineering calibration" category as
# HPA_FEEDBACK_GAIN and the SGLT2i beta constants elsewhere in this
# codebase -- treated as a real, if imprecise, per-Unit scaling factor,
# not a physical volume of distribution.
INSULIN_REGULAR_PM_PER_UNIT = 695.1

# Rapid-acting insulin analog (insulin lispro, [Lys(B28),Pro(B29)]-human
# insulin, brand name Humalog -- the zinc-free "LYSPRO" arm specifically,
# Howey et al. 1994's own "rapidly absorbed" arm) -- back-solved the same
# way as INSULIN_REGULAR_KA_PER_MIN/INSULIN_REGULAR_PM_PER_UNIT above, so
# that a 10-Unit dose reproduces Howey et al. 1994's own observed
# rapid-acting Tmax (42 min) and peak EXCESS above baseline (698 pM)
# exactly through the FULL composed ODE. This calibration is CLEANER than
# regular insulin's own (no cross-study confound to reason through, like
# regular insulin's need to explicitly NOT reuse Kobayashi et al. 1983's
# ke): both the regular-insulin and rapid-acting-insulin anchor points
# come from the SAME within-study, same-protocol, same-10U-dose
# comparison (Howey et al. 1994's own arms A vs. C -- "Serum insulin
# concentrations peaked more than two times higher (Cmax, 698 vs. 308 pM,
# A vs. C) and in less than half the time (Tmax, 42 vs. 101 min, A vs.
# C)"), verified live via NCBI efetch this session.
INSULIN_LISPRO_KA_PER_MIN = 0.0353
INSULIN_LISPRO_PM_PER_UNIT = 543.1

# Engineering plausibility ceiling, not a literature-derived clinical
# maximum -- added specifically because dose_units silently reinterprets
# DrugEvent.dose_mg as INSULIN UNITS for these drugs (see Drug
# .INSULIN_REGULAR's/Drug.INSULIN_LISPRO's own docstrings): a scenario
# author or REST caller copy-pasting a dose_mg value from a milligram-
# dosed drug (metformin ~500-2000, glucocorticoid ~20-60) into an insulin
# event would otherwise silently produce a many-fold-too-large insulin
# spike with no warning -- the same severity class as the previously-
# fixed Gb/G_threshold margin collapse (see CLAUDE.md KNOWN BUGS). Real
# single subcutaneous insulin boluses (regular OR rapid-acting) are
# clinically a few Units up to roughly 20-50 U even for insulin-resistant
# patients; this pathway models a single SC bolus only, not continuous IV
# infusion (a different route this platform doesn't represent). Chosen
# with headroom above that range, not tuned to any specific benchmark.
# Named generically (not "REGULAR"-specific) since it's enforced once, in
# the shared __init__ below, for whichever insulin analog is constructed.
MAX_PLAUSIBLE_SC_INSULIN_DOSE_UNITS = 100.0


class ExogenousInsulinParameters:
    """Parameters for a single subcutaneous insulin dose (regular or a rapid-acting analog)."""

    def __init__(
        self,
        dose_units: float,
        ka: float = INSULIN_REGULAR_KA_PER_MIN,
        pm_per_unit: float = INSULIN_REGULAR_PM_PER_UNIT,
        dose_time_min: float = 0.0,
    ) -> None:
        if not math.isfinite(dose_units):
            raise ValueError(f"dose_units must be finite, got {dose_units}")
        if dose_units < 0:
            raise ValueError(f"dose_units cannot be negative, got {dose_units}")
        if dose_units > MAX_PLAUSIBLE_SC_INSULIN_DOSE_UNITS:
            raise ValueError(
                f"dose_units={dose_units} exceeds the documented plausibility "
                f"ceiling of {MAX_PLAUSIBLE_SC_INSULIN_DOSE_UNITS} U for a "
                "single subcutaneous insulin bolus -- dose_units is "
                "measured in Units, not milligrams; check this value wasn't "
                "copy-pasted from a milligram-dosed drug"
            )
        if not math.isfinite(ka) or ka <= 0:
            raise ValueError(f"ka must be finite and positive, got {ka}")
        if not math.isfinite(pm_per_unit) or pm_per_unit < 0:
            raise ValueError(f"pm_per_unit must be finite and non-negative, got {pm_per_unit}")
        if not math.isfinite(dose_units * pm_per_unit):
            raise ValueError(
                f"dose_units * pm_per_unit (S1_0) must be finite, got "
                f"dose_units={dose_units}, pm_per_unit={pm_per_unit}"
            )
        if not math.isfinite(dose_time_min):
            raise ValueError(f"dose_time_min must be finite, got {dose_time_min}")
        self.dose_units = dose_units
        self.ka = ka
        self.pm_per_unit = pm_per_unit
        self.dose_time_min = dose_time_min

    @property
    def S1_0(self) -> float:
        """Initial depot content, in the same units delta_pM_per_min() uses internally."""
        return self.dose_units * self.pm_per_unit

    def delta_pM_per_min(self, t_minutes: float) -> float:
        """
        The exogenous_insulin_delta signal at time t_minutes (see this
        module's own docstring for the closed-form derivation). Exactly
        0.0 before dose_time_min (the dose hasn't been given yet) and for
        dose_units=0 (no dose at all, e.g. a scenario with no insulin
        event -- this pathway is simply never composed in that case, but
        this also holds if constructed directly with a zero dose).
        """
        t_since = t_minutes - self.dose_time_min
        if t_since <= 0.0 or self.S1_0 <= 0.0:
            return 0.0
        return (self.ka ** 2) * self.S1_0 * t_since * math.exp(-self.ka * t_since)

    @classmethod
    def regular_insulin(cls, dose_units: float, dose_time_min: float = 0.0) -> "ExogenousInsulinParameters":
        """
        A single subcutaneous dose of regular (short-acting) human
        insulin -- see INSULIN_REGULAR_KA_PER_MIN/INSULIN_REGULAR_PM_
        PER_UNIT's own docstrings for the citation and calibration.
        """
        return cls(dose_units=dose_units, dose_time_min=dose_time_min)

    @classmethod
    def rapid_acting_insulin(cls, dose_units: float, dose_time_min: float = 0.0) -> "ExogenousInsulinParameters":
        """
        A single subcutaneous dose of rapid-acting insulin lispro
        ([Lys(B28),Pro(B29)]-human insulin, zinc-free) -- see
        INSULIN_LISPRO_KA_PER_MIN/INSULIN_LISPRO_PM_PER_UNIT's own
        docstrings for the citation and calibration.
        """
        return cls(
            dose_units=dose_units, ka=INSULIN_LISPRO_KA_PER_MIN,
            pm_per_unit=INSULIN_LISPRO_PM_PER_UNIT, dose_time_min=dose_time_min,
        )


class ExogenousInsulinPathway(BasePathway):
    """
    Stateless (zero-species) pathway injecting a time-varying
    exogenous_insulin_delta signal into GlucoseInsulinPathway, following a
    closed-form two-sequential-compartment subcutaneous absorption curve
    (see this module's own docstring for the derivation).
    """

    NAME    = "exogenous_insulin_v1"
    VERSION = "1.0.0"

    def __init__(self, params: ExogenousInsulinParameters) -> None:
        self.params = params

    @property
    def descriptor(self) -> PathwayDescriptor:
        return PathwayDescriptor(
            name=self.NAME,
            version=self.VERSION,
            maturity=PathwayMaturity.EXPERIMENTAL,
            description=(
                "Stateless subcutaneous insulin depot -- injected regular "
                "or rapid-acting (lispro) insulin appearing in plasma via "
                "a closed-form two-sequential-compartment (equal rate "
                "constant) absorption model."
            ),
            species_names=[],
            output_signals=["exogenous_insulin_delta"],
            input_signals=[],
            dependencies=[],
            tags=["insulin", "drug", "exogenous", "time-varying"],
            literature_refs=["PMID:8314011"],
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
        return {"exogenous_insulin_delta": self.params.delta_pM_per_min(t)}


def exogenous_insulin_couplings(gi_name: str, insulin_name: str) -> list:
    """
    Build the CouplingSpec connecting an ExogenousInsulinPathway to
    GlucoseInsulinPathway. Single signal, maps straight through -- no
    transform needed, since exogenous_insulin_delta is already expressed
    in GI's own signal vocabulary (rule 11).
    """
    from .coupling import CouplingSpec

    return [
        CouplingSpec(
            source_pathway=insulin_name, signal_name="exogenous_insulin_delta",
            target_pathway=gi_name, target_input="exogenous_insulin_delta",
        ),
    ]


def couple_with_gi(gi: "GlucoseInsulinPathway", insulin: "ExogenousInsulinPathway") -> "CompositePathway":
    """
    Build the CompositePathway connecting an ExogenousInsulinPathway to
    GlucoseInsulinPathway (see exogenous_insulin_couplings() for the
    signal contract).
    """
    from .coupling import CompositePathway

    return CompositePathway([gi, insulin], couplings=exogenous_insulin_couplings(gi.name, insulin.name))
