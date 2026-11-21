# SAVE TO: packages\simulation-api\src\biosim_api\scenario.py
"""
Scenario definition language.

A ScenarioDefinition is the complete, self-contained description of
a simulation run. It captures everything needed to reproduce the
result exactly: interventions, personalization, and solver settings.

Design principles:
    - Immutable after construction (frozen dataclasses where possible)
    - All times in seconds (platform canonical unit)
    - Serialisable to/from dict for API transport and result archiving
    - Every scenario has a name for logging and comparison
"""

from __future__ import annotations

import dataclasses
import hashlib
import hmac
import math
import secrets
from dataclasses import dataclass, field
from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:
    pass


@dataclass(frozen=True)
class MealEvent:
    """
    A single food intake event within a scenario.

    time_s:           seconds from simulation start when eating begins
    carbs_g:          grams of carbohydrate
    label:            human-readable label e.g. 'Breakfast', 'High-sugar snack'
    glycaemic_index:  optional food glycaemic index (0-100, glucose = 100).
                       When set, the two-compartment GutAbsorptionModel's rate
                       constants follow the GI. When None (and kabs is None),
                       the meal uses GutAbsorptionModel.default_meal() -- also
                       two-compartment (Ra(0)=0, delayed peak), with rate
                       constants that are a MODEL ASSUMPTION (see
                       DEFAULT_MEAL_ABSORPTION_PROVENANCE in gut_absorption.py).
                       Only GI>=99.5 (pure glucose) is calibrated; any lower
                       GI uses an UNVALIDATED PLACEHOLDER slope and the run
                       carries a warning saying so.
    kabs:             optional override for the monoexponential absorption
                       rate constant [min-1] (see MealInput.kabs). Mutually
                       exclusive with glycaemic_index (they're two different
                       absorption-shape mechanisms -- pick one). Exists so
                       callers that already compute their own kabs (e.g.
                       biosim_seeds.gi_meals.gi_adjusted_kabs()) can apply it.
                       Setting kabs selects the single-exponential curve
                       (peaks the instant you eat); None = two-compartment.
    complex_carb_fraction: fraction (0.0-1.0) of this meal's carbohydrate
                       that requires enzymatic digestion (starches,
                       disaccharides like sucrose/maltose) before it can be
                       absorbed, as opposed to free monosaccharides
                       (glucose, fructose) that absorb directly with no
                       digestion step. Alpha-glucosidase inhibitors
                       (acarbose -- see biosim_physiology.drug_interactions
                       .acarbose_kabs_scale()) act specifically on the
                       enzymatic-digestion step, so a meal that's mostly
                       free sugar (e.g. a fruit juice or soda) should see a
                       correspondingly SMALLER acarbose effect than one
                       that's mostly starch (e.g. bread, rice), not the
                       drug's full modelled effect regardless of
                       composition. Defaults to 1.0 (fully
                       enzymatic-digestion-dependent -- MASTER_PLAN 1.4's
                       original acarbose milestone's own behaviour,
                       unchanged) so every existing scenario is
                       unaffected; a caller who knows a meal is
                       substantially free sugar can set this lower.
    protein_g:         grams of protein co-ingested with this meal.
                       Drives NutrientSecretionPathway's protein-triggered
                       endogenous insulin secretion (see
                       biosim_ode.nutrient_secretion_pathway, MASTER_PLAN
                       row 1.28), independent of this meal's own
                       carbs_g-driven glucose appearance. Defaults to 0.0
                       (no protein, no effect -- every existing scenario
                       unaffected).
    """
    time_s:          float
    carbs_g:         float
    label:           str = ""
    glycaemic_index: Optional[float] = None
    kabs:            Optional[float] = None
    complex_carb_fraction: float = 1.0
    protein_g:       float = 0.0

    def __post_init__(self) -> None:
        if self.time_s < 0:
            raise ValueError(f"MealEvent time_s cannot be negative, got {self.time_s}")
        if self.carbs_g < 0:
            raise ValueError(f"MealEvent carbs_g cannot be negative, got {self.carbs_g}")
        if self.glycaemic_index is not None and self.kabs is not None:
            raise ValueError(
                "MealEvent cannot set both glycaemic_index and kabs -- they are "
                "two different absorption-shape mechanisms, pick one"
            )
        if not (0.0 <= self.complex_carb_fraction <= 1.0):
            raise ValueError(
                f"MealEvent complex_carb_fraction must be in [0.0, 1.0], "
                f"got {self.complex_carb_fraction}"
            )
        # Mandatory security-reviewer finding (MEDIUM, fixed before
        # finalizing this milestone): `< 0` alone doesn't reject NaN
        # (every NaN comparison is False) or +inf, matching the same gap
        # DrugEvent.__post_init__/InsulinDoseEvent.__post_init__ already
        # guard against explicitly with math.isfinite() -- protein_g had
        # inherited MealEvent's own weaker carbs_g/time_s convention
        # instead of that stronger one.
        # Mandatory security-reviewer finding (MEDIUM, fixed before
        # finalizing this milestone): `< 0` alone doesn't reject NaN
        # (every NaN comparison is False) or +inf, matching the same gap
        # DrugEvent.__post_init__/InsulinDoseEvent.__post_init__ already
        # guard against explicitly with math.isfinite() -- protein_g had
        # inherited MealEvent's own weaker carbs_g/time_s convention
        # instead of that stronger one.
        if not math.isfinite(self.protein_g) or self.protein_g < 0:
            raise ValueError(f"MealEvent protein_g must be finite and non-negative, got {self.protein_g}")

    def to_meal_input(self, body_weight_kg: float = 70.0):
        from biosim_ode.meal import MealInput
        from biosim_ode.gut_absorption import GutAbsorptionModel

        absorption_model = None
        if self.glycaemic_index is not None:
            # GI=100 (pure glucose, e.g. an OGTT load) gets the preset
            # actually fit against real OGTT timing benchmarks -- see
            # GutAbsorptionModel.ogtt_calibrated_glucose(). Every other GI
            # value still falls back to the generic for_gi() interpolation,
            # which is NOT yet validated against clinical timing data (see
            # CLAUDE.md VALIDATION STATUS) -- don't assume it's calibrated.
            if self.glycaemic_index >= 99.5:
                absorption_model = GutAbsorptionModel.ogtt_calibrated_glucose()
            else:
                absorption_model = GutAbsorptionModel.for_gi(self.glycaemic_index)
        elif self.kabs is None:
            # No GI and no kabs: two-compartment default (MODEL ASSUMPTION
            # rate constants, see GutAbsorptionModel.default_meal()).
            absorption_model = GutAbsorptionModel.default_meal()
        kwargs = {"kabs": self.kabs} if self.kabs is not None else {}
        return MealInput(
            carbs_g=self.carbs_g,
            time_min=self.time_s / 60.0,
            body_weight_kg=body_weight_kg,
            label=self.label or f"{self.carbs_g:.0f}g carbs at t={self.time_s:.0f}s",
            absorption_model=absorption_model,
            protein_g=self.protein_g,
            **kwargs,
        )


@dataclass(frozen=True)
class DrugEvent:
    """
    A drug/substance intake event within a scenario.

    time_s:   seconds from simulation start when drug is taken. Must be
              >= 0 -- a drug event always sits within the simulated
              window (unlike ExerciseBout, there is no in-progress-during-
              simulation drug pathway; personalization is a static
              t=0 snapshot). Use hours_before_simulation for a dose
              already taken before t=0.
    drug:     Drug enum value
    dose_mg:  dose in milligrams for every drug EXCEPT Drug.INSULIN_REGULAR
              and Drug.INSULIN_LISPRO, for which this field is silently
              reinterpreted as INSULIN UNITS (see ExogenousInsulinPathway's
              own module docstring and Drug.INSULIN_REGULAR's/
              Drug.INSULIN_LISPRO's docstrings in drug_interactions.py)
              -- a 10.0 dose_mg for either insulin drug means 10 Units,
              not 10mg. Do not copy a dose_mg value from a milligram-dosed
              drug (metformin ~500-2000, glucocorticoid ~20-60) into an
              insulin event.
    hours_before_simulation: For a dose already taken before the
              simulated window starts (e.g. caffeine 30 min before a
              meal) -- time_s is still 0 in that case (it can't be
              negative), so this field carries the already-elapsed time
              directly, rather than losing it. Regression: this used to
              have no effect at all -- to_drug_dose()'s old formula,
              (simulation_start_s - time_s) with time_s>=0 enforced and
              simulation_start_s always called as 0.0, could only ever
              produce <= 0, clamped to exactly 0 -- every dose in the
              platform was silently modeled as just-ingested regardless
              of intended timing. See KNOWN BUGS FIXED in CLAUDE.md.
    """
    time_s:  float
    drug:    str   # Drug.value string — avoids circular import at definition time
    dose_mg: float
    hours_before_simulation: float = 0.0

    def __post_init__(self) -> None:
        # Regression (found by security review, glucocorticoid PK/PD
        # milestone -- see CLAUDE.md): "< 0" does not reject NaN (NaN
        # comparisons are always False) or +inf, so a NaN/inf dose_mg/
        # time_s/hours_before_simulation used to survive construction and
        # only crash much later, deep inside scipy's ODE solver
        # (`asarray_chkfinite`/`lu_factor` raising an opaque "array must
        # not contain infs or NaNs") -- confirmed reachable via
        # glucocorticoid_secretion_scale_factor() (biosim_ode.cortisol_
        # pathway), whose own `dose_mg <= 0.0` guard has the identical gap.
        # Fixed at this single construction-time choke point instead of
        # patching every downstream consumer individually, consistent with
        # "validate at system boundaries" -- protects every current and
        # future caller of DrugDose/DrugEvent, not just this one path.
        if not math.isfinite(self.time_s):
            raise ValueError(f"DrugEvent time_s must be finite, got {self.time_s}")
        if not math.isfinite(self.dose_mg):
            raise ValueError(f"DrugEvent dose_mg must be finite, got {self.dose_mg}")
        if not math.isfinite(self.hours_before_simulation):
            raise ValueError(
                f"DrugEvent hours_before_simulation must be finite, got {self.hours_before_simulation}"
            )
        if self.time_s < 0:
            raise ValueError(f"DrugEvent time_s cannot be negative")
        if self.dose_mg < 0:
            raise ValueError(f"DrugEvent dose_mg cannot be negative")
        if self.hours_before_simulation < 0:
            raise ValueError(f"DrugEvent hours_before_simulation cannot be negative")
        # Regression: an unrecognized drug string used to survive
        # construction and only fail later, deep inside to_drug_dose()
        # (an uncaught Drug(self.drug) ValueError) -- a confusing crash
        # far from the actual mistake, and only reachable via the
        # personalization path; runner.py's _detect_glp1()/
        # _detect_hepatic() separately caught the same ValueError and
        # silently skipped composing the GLP-1/hepatic pathway with no
        # warning at all. Fail fast here instead, at construction time,
        # with a clear message -- consistent with this method's existing
        # validation of every other field.
        from biosim_physiology.drug_interactions import Drug
        try:
            Drug(self.drug)
        except ValueError:
            valid = ", ".join(sorted(d.value for d in Drug))
            raise ValueError(
                f"DrugEvent drug={self.drug!r} is not a recognized drug. "
                f"Valid values: {valid}"
            ) from None

    def to_drug_dose(self, simulation_start_s: float = 0.0):
        from biosim_physiology.drug_interactions import Drug, DrugDose
        hours_since = self.hours_before_simulation + max(0.0, simulation_start_s - self.time_s) / 3600.0
        return DrugDose(
            drug=Drug(self.drug),
            dose_mg=self.dose_mg,
            hours_since_ingestion=hours_since,
        )


@dataclass(frozen=True)
class InsulinDoseEvent:
    """
    A dedicated exogenous-insulin dose event, with genuine Unit semantics.

    MASTER_PLAN row 1.2b (event-type): DrugEvent's own dose_mg field is
    silently reinterpreted as insulin UNITS for Drug.INSULIN_REGULAR/
    Drug.INSULIN_LISPRO (see DrugEvent.dose_mg's own docstring) -- a
    disclosed but honestly confusing unit-reuse hack, kept only because it
    was the fastest way to ship 1.2a/1.2b without a new event type. This
    is the honestly-named alternative: an insulin dose is expressed in its
    own real unit from the start, never overloading a milligram field.
    Confirmed via the rapid-acting-analog work (1.2b analog) that the old
    DrugEvent-based path has NOT caused a real production problem -- this
    type is purely additive, not a replacement; SimulationRunner still
    supports both (see _detect_exogenous_insulin()'s own docstring), but
    a scenario may not set both for the same effective dose (see
    ScenarioDefinition.__post_init__'s mutual-exclusivity check).

    time_s:    seconds from simulation start when the dose is given. Must
               be >= 0 -- same "no in-progress drug pathway" reasoning as
               DrugEvent.time_s; use hours_before_simulation for a dose
               already given before t=0.
    insulin_type: "regular" or "lispro" (rapid-acting analog, insulin
               lispro) -- see ExogenousInsulinParameters.regular_insulin()/
               .rapid_acting_insulin() in biosim_ode.exogenous_insulin_pathway.
    dose_units: dose in genuine insulin Units (NOT milligrams). Bounded by
               ExogenousInsulinParameters.MAX_PLAUSIBLE_SC_INSULIN_DOSE_UNITS
               at pathway-construction time (single source of truth --
               not re-validated here) -- see that class's own docstring
               for why (a realistic mg-value copy-paste footgun).
    hours_before_simulation: same meaning as DrugEvent's own field -- time
               already elapsed before t=0 for a dose given before the
               simulated window starts.
    """
    time_s:  float
    insulin_type: str
    dose_units: float
    hours_before_simulation: float = 0.0

    def __post_init__(self) -> None:
        # Same finite/non-negative discipline as DrugEvent.__post_init__
        # (see that method's own comment for why NaN/+inf must be
        # rejected at construction, not left to crash deep inside scipy).
        if not math.isfinite(self.time_s):
            raise ValueError(f"InsulinDoseEvent time_s must be finite, got {self.time_s}")
        if not math.isfinite(self.dose_units):
            raise ValueError(f"InsulinDoseEvent dose_units must be finite, got {self.dose_units}")
        if not math.isfinite(self.hours_before_simulation):
            raise ValueError(
                f"InsulinDoseEvent hours_before_simulation must be finite, got {self.hours_before_simulation}"
            )
        if self.time_s < 0:
            raise ValueError("InsulinDoseEvent time_s cannot be negative")
        if self.dose_units < 0:
            raise ValueError("InsulinDoseEvent dose_units cannot be negative")
        if self.hours_before_simulation < 0:
            raise ValueError("InsulinDoseEvent hours_before_simulation cannot be negative")
        if self.insulin_type not in ("regular", "lispro"):
            raise ValueError(
                f"InsulinDoseEvent insulin_type must be 'regular' or 'lispro', "
                f"got {self.insulin_type!r}"
            )


@dataclass
class PersonalizationSpec:
    """
    Specification of all personalisation to apply to a scenario.
    Combines physiological profile, genetic variants, and a prior
    exercise bout.
    """
    physiological_profile: Optional[object] = None   # PhysiologicalProfile
    genetic_variants:      list             = field(default_factory=list)  # list[Variant]
    exercise_bout:         Optional[object] = None   # ExerciseBout -- see runner.py

    @property
    def has_any(self) -> bool:
        return (
            self.physiological_profile is not None
            or len(self.genetic_variants) > 0
            or self.exercise_bout is not None
        )


def _exercise_bout_to_dict(bout) -> Optional[dict]:
    """dataclasses.asdict() is safe here -- ExerciseBout carries no PII."""
    return dataclasses.asdict(bout) if bout is not None else None


# Process-local, in-memory only -- generated fresh at import time, never
# read from an env var or config file, and never itself serialized.
# See _variant_fingerprint()'s docstring for why a *keyed* scheme is
# required here, not a bare hash.
_FINGERPRINT_KEY = secrets.token_bytes(32)


def _variant_fingerprint(variant) -> str:
    """
    A keyed fingerprint for one genetic Variant, for cache-key purposes only.

    Never embed a raw rsid in anything that ends up in to_dict() -- rule 7
    (CLAUDE.md): "PersonalizationResult must never contain raw rsIDs or
    allele data", and biosim_genotype.variant's own module docstring says
    the same.

    A bare (unkeyed) hash of `rsid:genotype` is NOT sufficient here and was
    rejected after security review: `rsid` is drawn from a small, publicly
    cataloged space (this platform's own VariantEffectDB has ~14 known
    SNPs x 3 Genotype values = <=42 real combinations), so an unkeyed
    digest is an exhaustive lookup table, not a one-way function --
    truncation length is irrelevant to that, a full untruncated SHA-256
    would be equally reversible. Keying with a secret (HMAC) closes this:
    without `_FINGERPRINT_KEY`, brute-forcing the same small space no
    longer recovers anything, because the attacker doesn't know the key.

    Trade-off, by design: `_FINGERPRINT_KEY` is regenerated every process
    start and lives only in memory, so this fingerprint is NOT a stable
    cross-process/cross-restart identifier -- two processes hashing the
    same variant get different fingerprints. That's fine for
    SimulationCache (in-memory, per-process already). It means this
    fingerprint must NOT be relied on for cross-process equality (e.g. if
    ReproducibilityBundle is ever wired to actually persist
    `ScenarioDefinition.to_dict()`, per its own docstring's stated intent
    -- see CLAUDE.md -- reproducing a genetically-personalized scenario
    from an archived bundle needs a different, deliberately-designed
    scheme, not this one, and re-deriving it from the *original* Variant
    input, not from this fingerprint).

    Note: `biosim_infra.PrivacyAuditor` audits `SimulationOutput`/
    `PersonalizationResult` (rule 7's stated scope), not
    `ScenarioDefinition.to_dict()` -- it has no visibility into this
    function at all, and a fingerprint wouldn't pattern-match its
    rsid-shaped regex check either way. This function is the only privacy
    control for this particular surface; don't assume the auditor also
    covers it.

    Regression (caught by code review before merge): the first version of
    this function hashed only `rsid:genotype.value`, missing
    `observed_alleles`/`ref_allele`/`alt_allele`. `EffectComposer.compose()`
    calls `variant.dosage_for_allele(risk_allele)`, which reads
    observed_alleles/ref_allele/alt_allele FIRST and only falls back to
    genotype alone when none of those are set (see that method's own
    docstring -- this is the same real-world "23andMe homozygous call"
    bug class already fixed elsewhere in this codebase, see KNOWN BUGS
    FIXED in CLAUDE.md). Two homozygous calls for the same SNP with
    opposite observed letters (e.g. protective "CC" vs. risk "TT", both
    `Genotype.HOMOZYGOUS_ALT`) have the same rsid+genotype but opposite
    `dosage_for_allele()` results (0 vs. 2) -- and produced the same
    fingerprint under the old version, silently colliding in the cache
    for a genuinely different simulated output. `chromosome`/`position`
    are NOT included below -- confirmed `EffectComposer.compose()` only
    ever reads `variant.rsid` and `variant.dosage_for_allele()`, so those
    two fields are pure provenance metadata that can't affect the
    simulated output.
    """
    raw = (
        f"{variant.rsid}:{variant.genotype.value}:"
        f"{variant.observed_alleles}:{variant.ref_allele}:{variant.alt_allele}"
    )
    return hmac.new(_FINGERPRINT_KEY, raw.encode(), hashlib.sha256).hexdigest()[:16]


def _personalization_to_dict(spec: "PersonalizationSpec") -> dict:
    return {
        "physiological_profile": (
            dataclasses.asdict(spec.physiological_profile)
            if spec.physiological_profile is not None else None
        ),
        "genetic_variant_fingerprints": sorted(
            _variant_fingerprint(v) for v in spec.genetic_variants
        ),
        "exercise_bout": _exercise_bout_to_dict(spec.exercise_bout),
    }


@dataclass
class ScenarioDefinition:
    """
    Complete, self-contained description of a simulation scenario.

    This is the primary input to SimulationRunner. Everything needed
    to reproduce the simulation exactly is captured here.

    Fields:
        name:              Human-readable scenario identifier
        duration_s:        How long to simulate (seconds)
        meal_events:       Food intake events (time-stamped)
        drug_events:       Drug/substance intake events
        personalization:   Physiological + genetic personalisation
        run_uncertainty:   Whether to run Monte Carlo UQ (slower)
        n_mc_samples:      Number of Monte Carlo samples if UQ enabled
        description:       Free-text description for documentation
    """
    name:              str
    duration_s:        float = 14400.0   # 4 hours
    meal_events:       list[MealEvent]  = field(default_factory=list)
    drug_events:       list[DrugEvent]  = field(default_factory=list)
    # Dedicated exogenous-insulin dose events with genuine Unit semantics
    # (MASTER_PLAN 1.2b event-type) -- additive alongside drug_events, not
    # a replacement. See InsulinDoseEvent's own docstring and
    # __post_init__'s mutual-exclusivity check below (a scenario may not
    # mix this with a DrugEvent for insulin_regular/insulin_lispro).
    insulin_dose_events: list[InsulinDoseEvent] = field(default_factory=list)
    personalization:   PersonalizationSpec = field(default_factory=PersonalizationSpec)
    run_uncertainty:   bool  = False
    n_mc_samples:      int   = 100
    description:       str   = ""
    stress_level:      float = 0.0   # 0-1, drives CortisolPathway composition (see runner.py)
    cushings_syndrome: bool = False  # Track 1 row 1.38 -- chronic HPA-axis
        # hypercortisolism, a genuinely different mechanism/timescale from
        # stress_level (acute/chronic psychological stress, 0-1 continuous)
        # despite composing through the SAME CortisolPathway instance
        # (_detect_cortisol() in runner.py) -- both are HPA-axis states, not
        # a static PhysiologicalProfile.disease_states entry, since the
        # entire point of route (a) is reusing CortisolPathway's own
        # GR_REFERENCE machinery rather than risking a double-count via a
        # separate static disease-state p3 modifier (see MASTER_PLAN row
        # 1.38's own scoping for why route (b), a static modifier, was
        # rejected). Composable with stress_level/a glucocorticoid DrugEvent
        # -- a real patient can have Cushing's AND be acutely stressed AND
        # on exogenous steroids; all three route through one shared,
        # saturating GR-occupancy pool by physiological design (confirmed
        # non-additive/non-double-counting for stress+steroid already, same
        # reasoning applies here). Implemented on the fixed CortisolPathway
        # foundation (Section 8 housekeeping's coupled-fixed-point fix) --
        # a second attempt at this field was reverted when that foundation
        # was found broken; see HISTORY.md for the full account.
    model_nefa_dynamics: bool = False  # opt-in: compose NEFAPathway (EXPERIMENTAL maturity)
    exercise_bout_in_progress: Optional[object] = None   # ExerciseBout -- see runner.py
    # ^ Distinct from personalization.exercise_bout: THIS bout starts at
    # t=0 and runs for its own duration WITHIN the simulated window (drives
    # ExercisePathway composition, a time-gated coupling), vs.
    # personalization.exercise_bout which already ended before t=0 (a
    # static, one-shot PhysiologicalPersonalizer modifier). Don't set both
    # for the same bout -- they model different things.

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("ScenarioDefinition name cannot be empty")
        if self.duration_s <= 0:
            raise ValueError(f"duration_s must be positive, got {self.duration_s}")
        # Found by the mandatory security-reviewer pass for the
        # n_mc_samples fix below -- CRITICAL, addressed in the same
        # pass rather than deferred. An unbounded duration_s forces the
        # Radau ODE solver to integrate across an arbitrarily large
        # time span, an unbounded time-cost lever independent of
        # n_mc_samples/run_uncertainty entirely. 7 days (604800s) is
        # far beyond this platform's own longest documented standalone
        # integration checks (~600-1000 minutes, i.e. under 1 day) while
        # remaining generous for a legitimate multi-day scenario.
        if self.duration_s > 604800.0:
            raise ValueError(f"duration_s must be <= 604800 (7 days), got {self.duration_s}")
        # Same pass: unbounded event-list length here, combined with
        # MealAbsorptionPathway allocating one ODE state variable PER
        # meal event (whenever a GLP-1-class drug is active) and
        # ODESystem.simulate()'s Radau solver estimating a DENSE
        # finite-difference Jacobian (no jac/jac_sparsity supplied),
        # made a very large meal_events list (needing only one trivial
        # drug event to trigger MealAbsorptionPathway) an O(N^2)-memory,
        # O(N^3)-per-step cost lever -- reachable with no
        # run_uncertainty involved at all. 50 is far beyond any
        # realistic scenario (this platform's own library scenarios use
        # 1-3 meals) while remaining generous for a legitimate
        # multi-day, multi-meal simulation.
        for field_name, events in (
            ("meal_events", self.meal_events),
            ("drug_events", self.drug_events),
            ("insulin_dose_events", self.insulin_dose_events),
        ):
            if len(events) > 50:
                raise ValueError(
                    f"{field_name} must have <= 50 entries, got {len(events)}"
                )
        if self.n_mc_samples < 10:
            raise ValueError(f"n_mc_samples must be >= 10, got {self.n_mc_samples}")
        # Housekeeping-backlog fix: n_mc_samples previously had no upper
        # bound at all -- an unauthenticated REST client requesting
        # run_uncertainty=True with an arbitrarily large n_mc_samples
        # was an unbounded, uninterruptible compute job (measured
        # directly: ~18ms/sample through the real production path, so
        # even a "modest"-looking n_mc_samples=1_000_000 would run for
        # ~5 hours on a single request). MASTER_PLAN 3.11's own mandatory
        # code-reviewer pass already fixed the 20x-worse version of this
        # gap for /population and /batch (by rejecting run_uncertainty
        # outright on those two routes); this is the 1x version, on
        # plain /simulate, closing the same gap at its actual source.
        # 1000 is a real ceiling, not an arbitrary one: 10x this
        # platform's own default (100), and ~18s worst-case at the
        # measured per-sample cost -- comfortably inside a normal HTTP
        # request/response cycle, matching the same "completes within a
        # plain request" reasoning /population's n_individuals=20 and
        # /batch's 20-scenario cap already use.
        if self.n_mc_samples > 1000:
            raise ValueError(f"n_mc_samples must be <= 1000, got {self.n_mc_samples}")
        if not (0.0 <= self.stress_level <= 1.0):
            raise ValueError(f"stress_level must be 0-1, got {self.stress_level}")
        # Regression: "Don't set both for the same bout" (see
        # exercise_bout_in_progress's own comment above) was documentation
        # only, with no runtime guard -- the same unguarded-docstring-claim
        # shape that let metformin's and glucocorticoid's static-modifier-
        # vs-ODE-pathway double-counting ship unnoticed (see CLAUDE.md).
        # Setting both fields is never a valid combination (they model two
        # mutually exclusive things: a bout already over before t=0 vs. one
        # in progress during the simulated window) and would double-count
        # the same bout's insulin-sensitising effect through both the
        # static PhysiologicalPersonalizer modifier and ExercisePathway's
        # ODE coupling -- confirmed via direct investigation to produce a
        # real, measurable extra glucose drop when both are set. Fail fast
        # rather than leave this exploitable only by accident.
        if self.personalization.exercise_bout is not None and self.exercise_bout_in_progress is not None:
            raise ValueError(
                "Cannot set both personalization.exercise_bout and "
                "exercise_bout_in_progress -- they model mutually exclusive "
                "scenarios (a bout already ended before t=0 vs. one happening "
                "during the simulated window) and setting both double-counts "
                "the same bout's effect. Use exactly one."
            )
        # Same double-representation hazard as the exercise_bout check
        # above, for exogenous insulin: insulin_dose_events and a
        # DrugEvent for Drug.INSULIN_REGULAR/Drug.INSULIN_LISPRO are two
        # different ways to express the SAME mechanism (see
        # InsulinDoseEvent's own docstring) -- setting both would make
        # _detect_exogenous_insulin()'s "which one wins" ambiguous rather
        # than double-count (it prefers insulin_dose_events unconditionally),
        # but silently picking a winner the caller didn't realize they were
        # choosing between is exactly the kind of accidental-exploit shape
        # the exercise_bout check above already exists to prevent. Fail
        # fast instead.
        if self.insulin_dose_events:
            from biosim_physiology.drug_interactions import Drug
            legacy_insulin_drugs = {Drug.INSULIN_REGULAR.value, Drug.INSULIN_LISPRO.value}
            if any(d.drug in legacy_insulin_drugs for d in self.drug_events):
                raise ValueError(
                    "Cannot set both insulin_dose_events and a DrugEvent for "
                    "insulin_regular/insulin_lispro -- they are two different "
                    "ways to represent the same exogenous-insulin dose. Use "
                    "exactly one (insulin_dose_events is the preferred, "
                    "honestly-Unit-named form)."
                )

        # Drug combinations the model can't represent faithfully (two
        # distinct sulfonylureas; duplicate flat-factor events such as
        # HCQ; duplicate dose-scaled events such as metformin/acarbose).
        # Validated HERE, at construction, so the ValueError reaches
        # REST's existing 400 handler -- raising it only inside
        # PhysiologicalPersonalizer was swallowed by the runner's catch-all
        # whenever a physiological_profile was present, silently degrading
        # to a healthy-baseline simulation.
        if self.drug_events:
            from biosim_physiology.drug_interactions import validate_drug_combination
            validate_drug_combination(self.drug_events)

            # (A GLP-1-class drug combined with a glycaemic_index meal used to
            # be rejected here because it deleted meal glucose. Since
            # 2026-11-15 MealAbsorptionPathway state-tracks two-compartment
            # meals too, so the combination conserves mass and is allowed.)

    @property
    def has_meals(self) -> bool:
        return len(self.meal_events) > 0

    @property
    def has_drugs(self) -> bool:
        return len(self.drug_events) > 0

    @property
    def is_personalised(self) -> bool:
        return self.personalization.has_any

    def to_dict(self) -> dict:
        """
        Serialise to a JSON-compatible dict for archiving, API transport,
        and (via SimulationCache._key()) cache-key hashing.

        Regression: this used to omit MealEvent.glycaemic_index/kabs,
        DrugEvent.hours_before_simulation, exercise_bout_in_progress, and
        `personalization` entirely -- SimulationCache._key() hashes this
        dict verbatim, so two scenarios differing ONLY in one of those
        fields (e.g. a healthy vs. a Type-2-diabetes-personalized run of
        the same meal) produced byte-identical cache keys and would
        silently serve each other's cached result. Every field that can
        change the simulated output must appear here -- including
        insulin_dose_events (MASTER_PLAN 1.2b event-type), added at the
        same time as the field itself specifically to avoid repeating
        that exact regression for a brand-new field.
        """
        return {
            "name":          self.name,
            "duration_s":    self.duration_s,
            "description":   self.description,
            "run_uncertainty": self.run_uncertainty,
            "n_mc_samples":  self.n_mc_samples,
            "stress_level":  self.stress_level,
            "cushings_syndrome": self.cushings_syndrome,
            "model_nefa_dynamics": self.model_nefa_dynamics,
            "meal_events": [
                {
                    "time_s": m.time_s, "carbs_g": m.carbs_g, "label": m.label,
                    "glycaemic_index": m.glycaemic_index, "kabs": m.kabs,
                    "complex_carb_fraction": m.complex_carb_fraction,
                    "protein_g": m.protein_g,
                }
                for m in self.meal_events
            ],
            "drug_events": [
                {
                    "time_s": d.time_s, "drug": d.drug, "dose_mg": d.dose_mg,
                    "hours_before_simulation": d.hours_before_simulation,
                }
                for d in self.drug_events
            ],
            "insulin_dose_events": [
                {
                    "time_s": i.time_s, "insulin_type": i.insulin_type,
                    "dose_units": i.dose_units,
                    "hours_before_simulation": i.hours_before_simulation,
                }
                for i in self.insulin_dose_events
            ],
            "exercise_bout_in_progress": _exercise_bout_to_dict(self.exercise_bout_in_progress),
            "personalization": _personalization_to_dict(self.personalization),
        }
