# SAVE TO: packages\rest-api\src\biosim_rest\schemas.py
"""
Request/response schemas for the Biosim REST API.

Kept deliberately separate from biosim_api's internal dataclasses
(ScenarioDefinition, SimulationOutput, ...) per the wire-format boundary:
these are what clients see, not what the simulation engine uses internally.
Mirrors ScenarioDefinition.to_dict()'s existing fields exactly -- that
dict is already the accepted wire format for a scenario, so the REST
input contract doesn't invent a second one.
"""

from __future__ import annotations

from typing import Optional

from pydantic import BaseModel, Field


class MealEventIn(BaseModel):
    time_s: float = Field(ge=0)
    carbs_g: float = Field(ge=0)
    label: str = ""
    # Regression: MealEvent.glycaemic_index/kabs (the GI-adjusted
    # absorption-shape mechanism -- see KNOWN BUGS FIXED #15,
    # gi_meal_scenario()'s computed-but-discarded kabs) had no wire-schema
    # equivalent at all -- a REST client could never submit a GI-adjusted
    # meal, silently always getting the default monoexponential
    # absorption regardless of the food's real glycaemic index. Same root
    # cause/shape as the DrugEventIn.hours_before_simulation fix above,
    # just for meals instead of drugs.
    glycaemic_index: Optional[float] = Field(default=None, ge=0, le=100)
    kabs: Optional[float] = Field(default=None, gt=0)
    # Fraction of this meal's carbohydrate requiring enzymatic digestion
    # (starch/disaccharides) vs. free monosaccharides that absorb
    # directly -- lets an alpha-glucosidase inhibitor (acarbose)'s
    # effect scale down for mostly-free-sugar meals instead of applying
    # its full modelled effect regardless of composition (MASTER_PLAN
    # row 1.4b(iii)). Defaults to 1.0, matching MealEvent's own default
    # and preserving every existing scenario's behaviour unchanged.
    complex_carb_fraction: float = Field(default=1.0, ge=0.0, le=1.0)
    # Grams of protein co-ingested with this meal -- drives
    # NutrientSecretionPathway's protein-triggered endogenous insulin
    # secretion (MASTER_PLAN row 1.28). Defaults to 0.0 (no protein, no
    # effect), matching MealEvent's own default so every existing
    # scenario is unaffected. Upper bound matches
    # MAX_PLAUSIBLE_PROTEIN_G_PER_MEAL in nutrient_secretion_pathway.py
    # (mandatory security-reviewer finding, HIGH -- fixed before
    # finalizing this milestone: the un-ceilinged dose-response is
    # purely linear with no saturating term, and a 300g-protein-alone
    # meal was verified to push peak_insulin_pM to ~3.2-3.6x this
    # platform's own validated OGTT reference). An engineering
    # plausibility ceiling, not a literature maximum -- same category as
    # this schema's own body_weight_kg/height_cm bounds.
    protein_g: float = Field(default=0.0, ge=0.0, le=150.0)


class DrugEventIn(BaseModel):
    time_s: float = Field(ge=0)
    drug: str
    # dose_mg is in milligrams for every drug EXCEPT "insulin_regular"
    # and "insulin_lispro", for which it is silently reinterpreted as
    # INSULIN UNITS (see ExogenousInsulinPathway's module docstring and
    # DrugEvent.dose_mg's own docstring in biosim_api.scenario) -- do not
    # send a dose_mg value copy-pasted from a milligram-dosed drug
    # (metformin ~500-2000, glucocorticoid ~20-60) for either insulin
    # event.
    dose_mg: float = Field(ge=0)
    # Regression: this field was missing from the wire schema entirely --
    # DrugEvent.to_drug_dose()'s hours_since_ingestion is computed purely
    # from hours_before_simulation (time_s never contributes for the
    # standard simulation_start_s=0.0 call every real caller uses), so
    # every time-dependent drug (caffeine, alcohol) submitted through the
    # REST API was silently modeled at hours_since_ingestion=0 regardless
    # of the client's intended timing -- for caffeine this made the drug
    # effect exactly zero (concentration_fraction's t<=0 floor), with no
    # way to express "N hours before the meal" via this API at all. Same
    # root cause as the data-seeds scenario-builder bugs already fixed
    # this session, just at the REST wire-format layer instead.
    hours_before_simulation: float = Field(default=0.0, ge=0)


class InsulinDoseEventIn(BaseModel):
    """
    A dedicated exogenous-insulin dose event with genuine Unit semantics
    (MASTER_PLAN 1.2b event-type) -- the preferred alternative to
    DrugEventIn's dose_mg-as-Units hack for "insulin_regular"/
    "insulin_lispro" (see that field's own docstring). Mirrors
    biosim_api.scenario.InsulinDoseEvent field-for-field; may not be
    combined with a DrugEventIn for either insulin drug in the same
    scenario (ScenarioDefinition.__post_init__ rejects that combination).
    """
    time_s: float = Field(ge=0)
    insulin_type: str = Field(pattern="^(regular|lispro)$")
    dose_units: float = Field(ge=0)
    hours_before_simulation: float = Field(default=0.0, ge=0)


class PhysiologicalProfileIn(BaseModel):
    """
    Mirrors biosim_physiology.profile.PhysiologicalProfile field-for-field
    (MASTER_PLAN row 3.14 / process-improvement item C2). Confirmed via
    direct grep this was previously entirely absent from ScenarioIn -- a
    REST client could submit drugs and meals but had no way to say WHO
    the scenario represents (age/sex/BMI/disease state/fitness/time of
    day). `biological_sex`/`fitness_level`/`disease_states` are left as
    plain strings, not Pydantic enums -- matching DrugEventIn.drug's own
    established convention of deferring enum validation to the internal
    dataclass's own __post_init__ (BiologicalSex(...)/FitnessLevel(...)/
    DiseaseState(...) all raise ValueError for an invalid value, already
    caught by /simulate's existing `except ValueError -> 400` handler,
    so an invalid value here needs no separate validation path).
    Deliberately excludes genetic_variants (rule 7's raw-rsID privacy
    discipline needs its own dedicated design pass, not a drive-by
    addition here) and exercise_bout (a separate, not-yet-scoped
    follow-on -- see MASTER_PLAN 3.14's own note).
    """
    age_years: float = Field(default=30.0, ge=0, le=120)
    biological_sex: str = "not_specified"
    # Upper bounds added (security-reviewer finding, C2's mandatory pass):
    # unlike age_years/body_fat_percent, these had only a lower bound --
    # traced downstream to MealInput.Vg_L = body_weight_kg * Vg_L_per_kg;
    # an extreme NEAR-ZERO (still gt=0) value can inflate Ra(t)/D_eff_mM
    # toward overflow inside the ODE. Not a crash risk today (solve_ivp
    # already degrades gracefully to success=False, SimulationRunner to
    # a zeroed ExpertMetrics), but bounded here for consistency and
    # defense-in-depth rather than relying on that fallback alone. 500kg/
    # 300cm are both far beyond any real human (tallest recorded ~272cm).
    body_weight_kg: float = Field(default=70.0, gt=0, le=500.0)
    height_cm: float = Field(default=170.0, gt=0, le=300.0)
    body_fat_percent: Optional[float] = Field(default=None, gt=0, lt=100)
    fitness_level: str = "moderately_active"
    # Bounded (Q2 checklist, extreme-value check): 13 real DiseaseState
    # values exist (as of NAFLD) -- matching this project's established
    # defense-in-depth convention for list fields (meal_events/
    # drug_events/insulin_dose_events at max_length=50), rather than
    # leaving this the one unbounded list on this endpoint. NOTE this cap
    # bounds request size, NOT the composed p3: a 10-state stack can still
    # reach the p3 floor, which the personalizer's clamp_warnings handles.
    disease_states: list[str] = Field(default_factory=lambda: ["healthy"], max_length=10)
    time_of_day_hours: Optional[float] = Field(default=None, ge=0, lt=24)
    sleep_hours: Optional[float] = Field(default=None, ge=0, le=24)
    # Bounded (Q2 checklist, extreme-value check) at 200 events/h -- far
    # beyond any clinically documented AHI/desaturation rate (severe OSA
    # is typically <100/h) and, more importantly, well beyond
    # intermittent_hypoxia_modifiers()'s own internal clamp at 24.3
    # events/h (Louis & Punjabi 2009's own single tested rate), so this
    # bound is pure defense-in-depth against a nonsensical input, not a
    # physiologically meaningful ceiling.
    desaturation_events_per_hour: Optional[float] = Field(default=None, ge=0, le=200)
    # Left as a plain string, not a Pydantic enum -- matching
    # biological_sex/fitness_level/disease_states's own established
    # convention of deferring enum validation to the internal dataclass's
    # own __post_init__ (MenstrualCyclePhase(...) raises ValueError for an
    # invalid value, already caught by /simulate's existing
    # `except ValueError -> 400` handler).
    menstrual_cycle_phase: Optional[str] = Field(default=None)
    cold_acclimatized: bool = Field(
        default=False,
        description=(
            "Completed a ~10-day chronic cold-acclimation protocol (Hanssen "
            "et al. 2015, PMID:26147760) -- increases insulin sensitivity. "
            "The source study evaluated Type 2 Diabetes subjects only (no "
            "healthy-subject arm); setting this without "
            "disease_states=['type_2_diabetes'] triggers a disclosure "
            "warning in the response's `warnings` list."
        ),
    )


class ScenarioIn(BaseModel):
    """Request body for POST /simulate -- mirrors ScenarioDefinition.to_dict()."""
    name: str = Field(min_length=1)
    # Upper bound added (found by the mandatory security-reviewer pass
    # for the n_mc_samples fix immediately below -- CRITICAL, addressed
    # in the same pass rather than deferred, per this project's own
    # security-response protocol). duration_s previously had no upper
    # bound -- an arbitrarily large value forces the Radau ODE solver
    # to integrate across an arbitrarily large time span, an unbounded
    # time-cost lever independent of n_mc_samples. 7 days (604800s) is
    # far beyond any real scenario this platform has ever run (the
    # longest documented standalone integration checks in this
    # codebase's own history are ~600-1000 minutes, i.e. well under 1
    # day) while still being generous enough for any legitimate
    # multi-day fasting/chronic-dosing scenario.
    duration_s: float = Field(default=14400.0, gt=0, le=604800.0)
    description: str = ""
    # Upper bounds added (same security-review pass): unbounded list
    # length here, combined with MealAbsorptionPathway allocating one
    # ODE state variable PER meal event (whenever a GLP-1-class drug is
    # active) and ODESystem.simulate()'s Radau solver estimating a
    # DENSE finite-difference Jacobian (no jac/jac_sparsity supplied),
    # made a single unauthenticated request with a very large
    # meal_events list (needing only one trivial drug_events entry to
    # trigger MealAbsorptionPathway) an O(N^2)-memory, O(N^3)-per-step
    # DoS vector -- reachable with no run_uncertainty, no /population,
    # no /batch involved at all. 50 is far beyond any realistic
    # scenario (this platform's own library scenarios use 1-3 meals)
    # while remaining generous for a legitimate multi-day, multi-meal
    # simulation.
    meal_events: list[MealEventIn] = Field(default_factory=list, max_length=50)
    drug_events: list[DrugEventIn] = Field(default_factory=list, max_length=50)
    insulin_dose_events: list[InsulinDoseEventIn] = Field(default_factory=list, max_length=50)
    run_uncertainty: bool = False
    # Upper bound added (housekeeping-backlog fix, MASTER_PLAN 3.13-
    # adjacent): previously unbounded, an unauthenticated-client-
    # reachable, uninterruptible compute-cost lever on this single-
    # simulation endpoint (~18ms/sample measured -- see
    # ScenarioDefinition.__post_init__'s own comment for the full
    # reasoning and the 1000 ceiling's derivation). This mirrors, at
    # 1x cost, the same gap MASTER_PLAN 3.11 already closed for
    # /population and /batch at their own 20x-worse multiplier.
    n_mc_samples: int = Field(default=100, ge=10, le=1000)
    stress_level: float = Field(default=0.0, ge=0.0, le=1.0)
    cushings_syndrome: bool = Field(
        default=False,
        description=(
            "Models Cushing's syndrome's postprandial insulin-action impairment "
            "only (via CortisolPathway's p3_modifier) -- does NOT model the "
            "disease's defining fasting hyperglycemia (basal glucose Gb is left "
            "unchanged). A matching disclosure warning is included in the "
            "response's `warnings` list whenever this is set."
        ),
    )
    model_nefa_dynamics: bool = False
    physiological_profile: Optional[PhysiologicalProfileIn] = None


class ExpertMetricsOut(BaseModel):
    glucose_peak_mM: float
    glucose_peak_time_s: float
    glucose_auc_uM_s: float
    glucose_baseline_uM: float
    insulin_peak_pM: float
    insulin_peak_time_s: float
    mtor_peak: float
    ampk_nadir: float
    time_to_glucose_return_s: Optional[float]
    overall_watchdog_confidence: float
    # Regression (MASTER_PLAN 3.11b): ExpertMetrics.glucose_nadir_mM was
    # added to feed the hypoglycaemia-risk lay-output indicator (see
    # biosim_lay.report/risk_score), but this wire schema was never
    # updated to carry it -- a REST client had no way to see the lowest
    # glucose value reached anywhere in the trajectory, only the peak.
    # Same "computed a value, missing wire exposure" bug class already
    # fixed for DrugEventIn.hours_before_simulation/
    # MealEventIn.glycaemic_index.
    glucose_nadir_mM: float


class LaySummaryOut(BaseModel):
    times_min: list[float]
    energy_level: list[float]
    hunger_level: list[float]
    glucose_narrative: str
    energy_narrative: str
    confidence_note: str
    is_high_confidence: bool


class RiskIndicatorOut(BaseModel):
    """Mirrors biosim_lay.report.RiskIndicator -- traffic-light assessment
    of one aspect of the simulation."""
    label: str
    value: float
    unit: str
    traffic_light: str
    description: str
    caveat: str


class RiskProfileOut(BaseModel):
    """Mirrors biosim_lay.report.RiskProfile (assess_risk())."""
    indicators: list[RiskIndicatorOut]
    overall: str
    summary: str


class MetabolicRiskIndicatorOut(BaseModel):
    """Mirrors biosim_lay.risk_score.RiskIndicatorScore."""
    name: str
    value: float
    unit: str
    band: str
    interpretation: str
    caveat: str


class MetabolicRiskProfileOut(BaseModel):
    """Mirrors biosim_lay.risk_score.MetabolicRiskProfile (compute_risk_profile())."""
    indicators: list[MetabolicRiskIndicatorOut]
    overall_band: str
    headline: str
    disclaimer: str
    n_elevated_or_high: int


class SimulateResponse(BaseModel):
    scenario_name: str
    scenario_description: str
    success: bool
    expert_metrics: ExpertMetricsOut
    lay_summary: LaySummaryOut
    # Regression (MASTER_PLAN 3.11b): neither assess_risk() nor
    # compute_risk_profile()'s output was exposed via /simulate at all --
    # a REST client could get raw numbers (expert_metrics) or a narrative
    # (lay_summary) but never the traffic-light/population-percentile
    # risk assessment this platform's own biosim_lay package computes.
    # Both are computed from the same SimulationOutput /simulate already
    # builds, unconditionally (both functions degrade gracefully, not
    # crash, on a failed simulation -- see _compute_expert_metrics()'s
    # own failure branch, which returns a valid zeroed ExpertMetrics).
    risk_profile: RiskProfileOut
    metabolic_risk_profile: MetabolicRiskProfileOut
    times_s: list[float]
    glucose_uM: list[float]
    insulin_uM: list[float]
    warnings: list[str]


class CompareRequest(BaseModel):
    """Request body for POST /compare -- names of two library scenarios."""
    scenario_a: str
    scenario_b: str


class CompareResponse(BaseModel):
    scenario_a: SimulateResponse
    scenario_b: SimulateResponse
    peak_glucose_diff_mM: float
    peak_insulin_diff_pM: float
    narrative: str


class ScenarioListItem(BaseModel):
    name: str
    description: str
    duration_s: float


class HealthResponse(BaseModel):
    status: str
    version: str


class BenchmarkResultOut(BaseModel):
    name: str
    tolerance: str
    is_scored: bool
    reference: float
    unit: str
    model_value: float
    deviation_percent: float
    passed: bool


class ValidationDatasetOut(BaseModel):
    dataset_name: str
    verdict: str
    n_pass: int
    n_fail: int
    results: list[BenchmarkResultOut]
    notes: list[str]
    # "" for non-drug benchmarks; for a drug benchmark, says the drug module
    # is experimental. Never changes the verdict.
    module_label: str = ""


class ValidationSummaryResponse(BaseModel):
    """Response body for GET /validation -- current status of every
    pre-registered clinical benchmark (see ModelValidator.validate_all())."""
    datasets: list[ValidationDatasetOut]


class ErrorResponse(BaseModel):
    detail: str


class PopulationRequest(BaseModel):
    """
    Request body for POST /population -- runs `scenario` through a
    synthetic demographic cohort (MASTER_PLAN 3.11).

    Bounded and always sequential by design: `PopulationRunner` itself
    supports an opt-in `use_multiprocessing` flag backed by its own
    `ProcessPoolExecutor`, but exposing that toggle over a synchronous
    HTTP route would let concurrent requests each spawn their own OS
    process pool with no shared limit -- exactly the nested-parallelism
    resource-exhaustion risk this row's own scoping pass flagged and
    left unfixed. This endpoint never exposes that flag (always runs
    the `PopulationRunner` default, sequential) and caps `n_individuals`
    to a size a single HTTP request/response cycle can carry.
    """
    scenario: ScenarioIn
    n_individuals: int = Field(default=20, ge=1, le=20)
    random_seed: int = 42


class CohortSummaryOut(BaseModel):
    """Mirrors biosim_seeds.population.CohortSummary (minus the
    individual_profiles/individual_peaks_mM fields, which carry
    PhysiologicalProfile objects with no wire representation yet)."""
    n_total: int
    n_succeeded: int
    n_failed: int
    glucose_peak_mean_mM: float
    glucose_peak_sd_mM: float
    glucose_peak_p5_mM: float
    glucose_peak_p95_mM: float
    insulin_peak_mean_pM: float
    mtor_peak_mean: float
    ampk_nadir_mean: float
    pct_above_7p8_mM: float
    pct_above_11p1_mM: float
    failure_reasons: list[str]
    # e.g. the "EXPERIMENTAL drug module" note when the scenario applies a drug.
    warnings: list[str] = []


class BatchRequest(BaseModel):
    """
    Request body for POST /batch -- runs each of `scenarios`
    independently and reports per-scenario success/failure
    (MASTER_PLAN 3.11).

    Bounded (max 20 scenarios) and always single-threaded, same
    reasoning as `PopulationRequest` above: `BatchRunner` defaults to a
    real `ProcessPoolExecutor`, which this endpoint deliberately never
    uses -- see the /batch route in app.py.
    """
    scenarios: list[ScenarioIn] = Field(min_length=1, max_length=20)


class BatchResponseItem(BaseModel):
    name: str
    succeeded: bool
    result: Optional[SimulateResponse] = None
    failure_reason: Optional[str] = None


class BatchResponse(BaseModel):
    n_succeeded: int
    n_failed: int
    results: list[BatchResponseItem]
