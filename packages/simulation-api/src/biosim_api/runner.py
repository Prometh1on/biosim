# SAVE TO: packages\simulation-api\src\biosim_api\runner.py
"""
SimulationRunner -- central orchestrator for the full simulation pipeline.

Pipeline (in order):
    1. Build base PathwayParameters (population average)
    2. Apply physiological personalisation (age, sex, BMI, fitness, disease, drugs, circadian)
    3. Apply genetic personalisation (variant effects on parameters)
    4. Assemble MealInputs from MealEvents (body-weight adjusted)
    5. Run ODE simulation
    6. Optionally run Monte Carlo uncertainty quantification
    7. Compute ExpertMetrics from simulation result
    8. Compute LaySummary (phenotype translation)
    9. Return SimulationOutput

Sprint 9 note -- microservice extraction:
    This runner is designed to be extracted to a microservice at Sprint 9.
    All inputs are JSON-serialisable (ScenarioDefinition.to_dict()) and
    the output is a plain Python dataclass. The runner holds no state
    between calls and is safe to parallelise.
"""

from __future__ import annotations

import dataclasses
from typing import Optional

from biosim_ode import (
    GlucoseInsulinPathway,
    ODESystem,
    PathwayParameters,
)

from .output import SimulationOutput
from .pathway_assembly import PathwayAssemblyMixin
from .reporting import (  # noqa: F401 -- constant re-exported for existing importers
    UQ_WIDE_GLUCOSE_RELATIVE_WIDTH,
    ReportingMixin,
)
from .scenario import ScenarioDefinition
from .scenario_warnings import (
    clamp_hit_warnings,
    cold_acclimation_warning,
    cushings_warning,
    hcq_cohort_warning,
    meal_composition_warnings,
    nafld_cohort_warning,
    unverified_variant_warnings,
    experimental_drug_warnings,
)


class SimulationRunner(PathwayAssemblyMixin, ReportingMixin):
    """
    Runs a ScenarioDefinition through the full biosim pipeline.

    Usage:
        runner = SimulationRunner()
        output = runner.run(scenario)
        print(output.expert_metrics.glucose_peak_mM)
        print(output.lay_summary.glucose_narrative)
    """

    def __init__(
        self,
        n_output_points: int = 400,
        base_params: Optional[PathwayParameters] = None,
    ) -> None:
        self.n_output_points = n_output_points
        # Regression: this used to be hardcoded to a fresh PathwayParameters()
        # inside run() (below), so a caller-supplied base (e.g.
        # PlatformSession's literature-loaded self._base_params) was never
        # actually used for simulation -- only for export_sbml(). Currently
        # invisible because the seeded literature values happen to equal
        # these class defaults exactly, but that's a coincidence of today's
        # data, not a guarantee (the parameter DB is append-only/evolving).
        # See CLAUDE.md.
        self._base_params = base_params if base_params is not None else PathwayParameters()

    def run(self, scenario: ScenarioDefinition) -> SimulationOutput:
        """
        Execute a full simulation from a scenario definition.

        Args:
            scenario: Complete scenario specification

        Returns:
            SimulationOutput with all results and summaries
        """
        # 1 + 2. Physiological personalisation
        params, physio_result = self._personalize_physiology(scenario, self._base_params)
        # Surface any safety-clamp-engagement warning from personalisation
        # (e.g. MIN_GB_THRESHOLD_MARGIN_MM overriding a drug/disease
        # modifier's own composed Gb value) -- see PhysiologicalPersonalizer
        # .clamp_warnings' own docstring (MASTER_PLAN.md 1.23). Read after
        # _personalize_physiology()'s glucocorticoid step so this covers
        # `physio_result` in its final, possibly-replaced form, not a stale
        # pre-replace copy.
        clamp_warnings = list(physio_result.clamp_warnings) if physio_result is not None else []

        # 3. Genetic personalisation
        params, genetic_result, genetic_warnings = self._personalize_genetics(scenario, params)

        # 4. Assemble meal inputs
        meals, meal_warnings = _assemble_meals(scenario)

        # 5. Run ODE simulation
        gi_pathway = GlucoseInsulinPathway(params=params, meals=meals)
        scenario_warnings = self._scenario_warnings(scenario, gi_pathway)
        sim_result, simulation_warnings = self._simulate(scenario, gi_pathway, meals)

        # 6. Optional Monte Carlo UQ
        uq_result, uq_warnings = self._run_uncertainty(
            scenario, params, meals, genetic_result, sim_result,
        )

        # 7. Compute expert metrics
        expert_metrics = self._compute_expert_metrics(sim_result, params)

        # 8. Compute lay summary
        lay_summary = self._compute_lay_summary(sim_result, expert_metrics, uq_result)

        return SimulationOutput(
            scenario_name=scenario.name,
            scenario_description=scenario.description,
            scenario=scenario,
            simulation=sim_result,
            expert_metrics=expert_metrics,
            lay_summary=lay_summary,
            uq=uq_result,
            physiological_personalization=physio_result,
            genetic_personalization=genetic_result,
            warnings=[
                *clamp_warnings, *genetic_warnings, *meal_warnings,
                *scenario_warnings, *simulation_warnings, *uq_warnings,
            ],
        )

    def _personalize_physiology(
        self, scenario: ScenarioDefinition, params: PathwayParameters,
    ) -> tuple[PathwayParameters, Optional[object]]:
        """
        Stages 1 + 2: physiological personalisation, then the glucocorticoid
        secretion_scale adjustment. Returns (params, physio_result).
        """
        # ------------------------------------------------------------------
        # Glucocorticoid drug events are excluded from the drug_doses passed
        # to PhysiologicalPersonalizer here -- _detect_cortisol() below
        # ALWAYS composes a real ODE-coupled CortisolPathway whenever a
        # glucocorticoid event exists, and CortisolPathway's gb_modifier/
        # p3_modifier are a full, standalone, self-sufficient representation
        # of the SAME steroid-induced Gb/p3 effect glucocorticoid_modifiers()
        # (physiology/drug_interactions.py) also computes -- confirmed via a
        # fasted, zero-meal SimulationRunner run: prednisone 20mg composed
        # BOTH ways drifted fasting glucose to an extra +18% (and near-doubled
        # steady-state insulin) beyond either representation alone, same
        # double-counting bug class as metformin's hepatic coupling (see
        # hepatic.py's for_metformin() fix) -- just resolved the other
        # direction here, since (unlike hepatic's genuinely complementary
        # postprandial EGP kinetics) CortisolPathway has no meal/insulin
        # input at all and is meant to be a complete alternative, not an
        # addition. glucocorticoid_modifiers()'s secretion_scale compensatory
        # bump (+10%, or a T2D-specific decrease -- see
        # glucocorticoid_secretion_scale_factor()) is NOT represented by
        # CortisolPathway at all, so it is applied separately below, after
        # personalisation.
        #
        # Regression, fixed together with the timing work below (see
        # CLAUDE.md): this used to read only the event's dose_mg
        # (_glucocorticoid_dose_mg()) and build a FRESH DrugDose with no
        # hours_since_ingestion (defaulting to 0.0) purely to extract that
        # secretion_scale bump -- silently discarding the scenario's real
        # DrugEvent timing. Harmless while glucocorticoid_modifiers() itself
        # ignored timing, but became a live bug the moment secretion_scale
        # became timing-gated (glucocorticoid_secretion_scale_factor() below
        # returns a no-op 1.0 at hours_since_ingestion=0, i.e. "just dosed,
        # effect not yet onset") -- confirmed via
        # test_glucocorticoid_secretion_scale_uses_real_dose_timing (RED
        # against the pre-fix code, GREEN after). _glucocorticoid_event()
        # now carries the event through so both this block and
        # _detect_cortisol() below share one real source of dosing timing.
        # ------------------------------------------------------------------
        exercise_bout = scenario.personalization.exercise_bout
        glucocorticoid_event = self._glucocorticoid_event(scenario)
        personalization_drug_events = [
            de for de in scenario.drug_events if de.drug != "glucocorticoid"
        ]
        physio_result = None
        if (personalization_drug_events or exercise_bout is not None) \
                and scenario.personalization.physiological_profile is None:
            from biosim_physiology import PhysiologicalProfile
            drug_doses=[de.to_drug_dose(0.0) for de in personalization_drug_events]
            from biosim_physiology import PhysiologicalPersonalizer
            r2=PhysiologicalPersonalizer().personalize(
                PhysiologicalProfile(), params, drug_doses, exercise_bout=exercise_bout,
            )
            params=r2.personalized_params
            # Regression: this branch used to update `params` (so the
            # simulation math was always correct) but never assigned
            # `physio_result`, unlike the sibling branch below -- so
            # SimulationOutput.physiological_personalization stayed None
            # even though real personalization was applied. Every REST
            # /simulate request with drug_events (service.py never sets
            # `personalization`, so it always takes this branch) silently
            # lost that provenance -- confirmed downstream in infra.py's
            # create_bundle(), which falls back to unmodified default
            # PathwayParameters() whenever physiological_personalization
            # is None, recording the wrong parameters in the
            # ReproducibilityBundle.
            physio_result=r2
        if scenario.personalization.physiological_profile is not None:
            # Deliberately NOT wrapped in a catch-all: one used to turn any
            # failure here into a warning and simulate with base parameters,
            # so e.g. a T2D profile silently ran as a healthy subject (REST
            # returned 200). A ValueError now reaches REST's 400 handler and
            # anything else surfaces as the bug it is -- same as the
            # drug-only branch above, which never had a catch-all.
            from biosim_physiology import PhysiologicalPersonalizer
            profile = scenario.personalization.physiological_profile

            # Build drug doses from drug events (at t=0 of simulation)
            drug_doses = [
                de.to_drug_dose(simulation_start_s=0.0)
                for de in personalization_drug_events
            ]

            physio_result = PhysiologicalPersonalizer().personalize(
                profile=profile,
                base_params=params,
                drug_doses=drug_doses,
                exercise_bout=exercise_bout,
            )
            params = physio_result.personalized_params

        if glucocorticoid_event is not None:
            params, physio_result = _apply_glucocorticoid_secretion_scale(
                scenario, glucocorticoid_event, params, physio_result,
            )
        return params, physio_result

    @staticmethod
    def _personalize_genetics(
        scenario: ScenarioDefinition, params: PathwayParameters,
    ) -> tuple[PathwayParameters, Optional[object], list[str]]:
        """Stage 3. Returns (params, genetic_result, warnings)."""
        warnings: list[str] = []
        genetic_result = None
        if scenario.personalization.genetic_variants:
            try:
                from biosim_genotype import GenotypePersonalizer
                genetic_result = GenotypePersonalizer().personalize(
                    variants=scenario.personalization.genetic_variants,
                    base_params=params,
                )
                params = genetic_result.personalized_params
                if genetic_result.any_warnings:
                    warnings.extend(genetic_result.population_warnings)
            except Exception as e:
                # Deliberately do NOT embed str(e) here --
                # this branch is reachable with real Variant/rsid data in
                # scope (see biosim_seeds.population's genetic-variant
                # sampling, the first caller to exercise this at cohort
                # scale), and an exception message could in principle
                # echo a raw rsid/genotype value. Only the exception
                # type is surfaced; the full message is not exported
                # anywhere from this warning, per rule 7 (no raw rsIDs
                # in anything that could be exported, e.g.
                # CohortSummary.failure_reasons/population_export.py).
                warnings.append(
                    f"Genetic personalisation failed ({type(e).__name__}). "
                    "Using base parameters."
                )
        if genetic_result is not None:
            warnings.extend(unverified_variant_warnings(
                scenario.personalization.genetic_variants
            ))
        return params, genetic_result, warnings

    def _scenario_warnings(self, scenario: ScenarioDefinition, gi_pathway) -> list[str]:
        """
        Stage 5 disclosures that depend on the personalised pathway or the
        profile, not on the solver result (GLP-1/Gb, Cushing's, cohorts).
        """
        warnings = [
            w for w in (
                self._glp1_basal_glucose_warning(scenario, gi_pathway),
                cushings_warning(scenario),
                cold_acclimation_warning(scenario.personalization.physiological_profile),
            )
            if w is not None
        ]

        # Rows 1.44/1.45: same disclosure convention for HCQ's and NAFLD's
        # source cohorts. A drug-only request is judged against the default
        # profile, which is what the personalizer applied to it.
        from biosim_physiology import PhysiologicalProfile as _Profile
        cohort_profile = scenario.personalization.physiological_profile or _Profile()
        for cohort_warning in (
            hcq_cohort_warning(cohort_profile, scenario.drug_events),
            nafld_cohort_warning(cohort_profile, scenario.model_nefa_dynamics),
        ):
            if cohort_warning is not None:
                warnings.append(cohort_warning)

        # Drug validation label, once per drug: EXPERIMENTAL for all but
        # semaglutide, which is VALIDATED with a caveat (2026-11-21).
        warnings.extend(experimental_drug_warnings(scenario))
        return warnings

    def _glp1_basal_glucose_warning(
        self, scenario: ScenarioDefinition, gi_pathway,
    ) -> Optional[str]:
        # Known, disclosed limitation (2026-08-12 deep-search plan
        # follow-up, mandatory security-reviewer finding on the Orskov
        # 1994 GLP-1 meal-response recalibration): GLP1Pathway.__init__
        # now derives a genuine, stable analytic fixed point for ANY
        # basal_glucose_mM (fixing the prior severe, unbounded drift for
        # a disease-elevated Gb composed with a GLP-1 agonist) -- but
        # that fixed point's own MAGNITUDE can be physiologically
        # implausible once basal glucose exceeds GLP1_EC50_SECRETION,
        # because a long agonist half-life (e.g. semaglutide's ~1 week)
        # divides a genuinely nonzero glucose_stim_basal by a very small
        # clearance rate, producing an enormous equilibrium concentration
        # (e.g. T2D Gb=6.5mM + semaglutide 1mg: ~157,000 pM, vs. a
        # healthy-Gb value of ~142 pM). This does not perturb the
        # SIMULATED glucose/insulin trajectory in practice (incretin_
        # effect/gastric_emptying_modifier are already Hill-saturated
        # well below either magnitude), but the exported `glp1_pM`
        # output signal itself is not a believable plasma concentration
        # in this combination -- flagged here rather than silently
        # redesigning the agonist_scale/glucose-response coupling this
        # session (a larger architectural question, matching this
        # project's own established practice of a scoped, disclosed
        # follow-up rather than an unplanned mid-fix redesign).
        # Mandatory security-reviewer finding on THIS fix itself (LOW,
        # both addressed): checking membership in GLP1_DRUGS plus the
        # hardcoded module-level GLP1_EC50_SECRETION constant is
        # imprecise on two counts -- (1) a DPP-4 inhibitor has
        # agonist_scale=1.0 (no long-half-life amplification), so its
        # own fixed point stays physiologically ordinary even at a
        # disease-elevated Gb (~18.5 pM for T2D, well within this
        # module's own documented "~7-15 pM basal" range) and should
        # NOT trigger this warning; (2) a future preset could customize
        # ec50_secretion, silently desyncing the warning's threshold
        # from the actual instance's own behaviour. Fixed by detecting
        # the real GLP1Pathway instance (the same one _build_pathway()
        # will independently reconstruct -- a second, cheap, pure
        # construction, not a behavioural risk) and reading its own
        # params directly, rather than a drug-identity/module-constant
        # proxy for them.
        glp1_for_warning = self._detect_glp1(
            scenario.drug_events, basal_glucose_mM=gi_pathway.params.Gb,
        )
        if (
            glp1_for_warning is not None
            and glp1_for_warning.params.agonist_scale > 1.0
            and gi_pathway.params.Gb > glp1_for_warning.params.ec50_secretion
        ):
            return (
                "A GLP-1 class drug is active with a disease-elevated basal "
                f"glucose (Gb={gi_pathway.params.Gb:.2f} mM > "
                f"ec50_secretion={glp1_for_warning.params.ec50_secretion:.1f} mM) -- the "
                "GLP1Pathway's own basal steady-state concentration is a "
                "stable fixed point (no longer drifts), but its magnitude "
                "is not a physiologically believable plasma GLP-1 "
                "concentration in this combination. Simulated glucose/"
                "insulin dynamics are not meaningfully affected (the "
                "downstream incretin/gastric-emptying signals are already "
                "saturated), but do not treat the reported glp1_pM value "
                "itself as a real plasma concentration for this scenario."
            )
        return None

    def _simulate(
        self, scenario: ScenarioDefinition, gi_pathway, meals: list,
    ) -> tuple[object, list[str]]:
        """
        Stage 5: integrate the composed pathway. Returns (sim_result, warnings).

        Some interventions are not modelled as PathwayParameters modifiers
        -- they compose a real coupled ODE pathway with GlucoseInsulinPathway
        (GLP-1 drugs, stress-driven cortisol, T2D/metformin hepatic EGP
        kinetics, opt-in NEFA dynamics). _build_pathway() detects all of
        them and flat-composes whichever are active into one CompositePathway
        alongside the bare pathway, since CouplingSpec injection only
        resolves direct siblings within a single CompositePathway (nesting
        composites would silently break cross-extension signal passing).
        """
        pathway = self._build_pathway(scenario, gi_pathway)
        system = ODESystem(pathway=pathway)
        sim_result = system.simulate(
            duration_s=scenario.duration_s,
            meals=meals,
            n_output_points=self.n_output_points,
        )

        warnings: list[str] = []
        if not sim_result.success:
            warnings.append(f"ODE solver failed: {sim_result.message}")
        warnings.extend(clamp_hit_warnings(sim_result.clamp_hits))
        return sim_result, warnings

    def _run_uncertainty(
        self,
        scenario: ScenarioDefinition,
        params: PathwayParameters,
        meals: list,
        genetic_result,
        sim_result,
    ) -> tuple[Optional[object], list[str]]:
        """Stage 6: optional Monte Carlo UQ. Returns (uq_result, warnings)."""
        uq_result = None
        warnings: list[str] = []
        if scenario.run_uncertainty and sim_result.success:
            try:
                from biosim_uncertainty import LogNormalDistribution, MonteCarloEngine

                # Genetic variants add measurement/effect-size uncertainty on
                # top of base parameter uncertainty -- composer.py's own
                # docstring documents this as "combined in quadrature with
                # the base parameter uncertainty", but until now nothing
                # ever read ComposedEffect.genetic_cv back out: a
                # genotype-personalized run and a population-average run
                # got identically-wide confidence bands even though the
                # point estimate had been shifted by the variant. Per-
                # parameter (not the aggregate genetic_uncertainty_cv,
                # which is itself already a quadrature sum across all
                # modified parameters and would double-count here).
                genetic_cv_by_param: dict[str, float] = {}
                if genetic_result is not None:
                    genetic_cv_by_param = {
                        e.parameter_field: e.genetic_cv
                        for e in genetic_result.applied_effects
                    }

                engine = MonteCarloEngine(
                    n_samples=scenario.n_mc_samples,
                    random_seed=42,
                    base_params=params,
                )
                distributions = {
                    name: LogNormalDistribution.from_confidence_tier(
                        value, "silver", additional_cv=genetic_cv_by_param.get(name, 0.0)
                    )
                    for name, value in (
                        ("p1", params.p1),
                        ("p2", params.p2),
                        ("p3", params.p3),
                        ("ke", params.ke),
                        ("EC50_mTOR", params.EC50_mTOR),
                    )
                }
                # MASTER_PLAN 2.2: extension-pathway-specific parameters
                # (e.g. RenalClearanceParameters.beta) were still FIXED at
                # their point-estimate value inside every UQ sample, even
                # once the pathway STRUCTURE fix below made the ensemble
                # see the right composed pathway -- a drug scenario's
                # confidence band correctly reflected uncertainty in
                # p1/p2/p3/ke/EC50_mTOR but silently treated the drug's
                # own calibrated effect size as if it were exactly known.
                # Only renal.beta is sampled today: it's the one extension
                # parameter with a real, literature-derived RANGE behind
                # its point estimate (see RENAL_BETA_CV's own derivation
                # comment) -- GLP1Parameters.agonist_scale and every other
                # extension constant checked this session have no reported
                # dispersion at all, and assigning one would fabricate
                # false precision (see CLAUDE.md's negative-result entry
                # for agonist_scale specifically). Empty when no SGLT2i
                # drug is active in this scenario, in which case
                # extension_pathway_factory behaves identically to the
                # pathway_factory it replaces.
                extension_distributions: dict = {}
                if self._detect_renal(scenario) is not None:
                    from biosim_ode.renal_pathway import (
                        RENAL_BETA_CV,
                        SGLT2I_BETA_CLASS_DEFAULT,
                    )
                    extension_distributions["renal.beta"] = LogNormalDistribution.from_mean_cv(
                        SGLT2I_BETA_CLASS_DEFAULT, cv=RENAL_BETA_CV,
                    )
                # Regression: propagate() always built a bare
                # GlucoseInsulinPathway, never the same CompositePathway
                # (GLP1, renal/SGLT2i, cortisol, hepatic, exercise, meal-
                # absorption) _build_pathway() composes for the point
                # estimate (_simulate()) -- so for a drug/disease
                # scenario the reported "confidence band" could describe
                # a structurally different ODE system than sim_result.
                # Provably so for GLP-1/SGLT2i specifically: both have
                # ZERO static PathwayParameters effect
                # (_glp1_agonist_modifiers()/_sglt2_inhibitor_modifiers()
                # both return {}), so the pre-fix band was bit-for-bit
                # identical to a drug-free ensemble's. Reusing
                # _build_pathway() here (same as _simulate()) makes the UQ
                # ensemble see exactly what the point estimate does, for
                # every scenario, at zero behaviour change for any
                # scenario with no active extension (_build_pathway()
                # returns the bare pathway unchanged in that case).
                # extension_pathway_factory (not pathway_factory) is used
                # unconditionally -- with extension_distributions empty,
                # its 3rd argument is always {}, which _build_pathway()'s
                # own extension_overrides default already treats as a
                # no-op, so this is behaviourally identical to the old
                # 2-arg pathway_factory for every scenario with no active
                # SGLT2i drug.
                uq_result = engine.propagate(
                    distributions=distributions,
                    meals=meals,
                    duration_s=scenario.duration_s,
                    n_output_points=self.n_output_points,
                    extension_distributions=extension_distributions,
                    extension_pathway_factory=lambda p, m, ov: self._build_pathway(
                        scenario, GlucoseInsulinPathway(params=p, meals=m), ov,
                    ),
                )
            except Exception as e:
                warnings.append(f"Uncertainty quantification failed: {e}")
        return uq_result, warnings


def _apply_glucocorticoid_secretion_scale(
    scenario: ScenarioDefinition, glucocorticoid_event, params: PathwayParameters, physio_result,
) -> tuple[PathwayParameters, Optional[object]]:
    """
    glucocorticoid_modifiers()'s secretion_scale bump, which CortisolPathway
    does not represent (see _personalize_physiology()). Returns
    (params, physio_result).
    """
    from biosim_ode.cortisol_pathway import (
        glucocorticoid_secretion_scale_factor,
    )
    from biosim_physiology.profile import DiseaseState

    profile = scenario.personalization.physiological_profile
    is_t2d = (
        profile is not None
        and DiseaseState.TYPE_2_DIABETES in profile.disease_states
    )
    # to_drug_dose(0.0) is the same "elapsed time as of simulation
    # start" snapshot personalization_drug_events already use in
    # _personalize_physiology() -- correct for a dose already taken before t=0
    # (hours_before_simulation) or taken exactly at t=0, same
    # platform-wide static-personalization limitation as every
    # other drug (see DrugEvent's own docstring: "personalization
    # is a static t=0 snapshot").
    hours_since_ingestion = glucocorticoid_event.to_drug_dose(
        simulation_start_s=0.0
    ).hours_since_ingestion
    secretion_scale_factor = glucocorticoid_secretion_scale_factor(
        glucocorticoid_event.dose_mg, hours_since_ingestion, is_t2d,
    )
    params = dataclasses.replace(
        params, secretion_scale=params.secretion_scale * secretion_scale_factor
    )
    # Keep the reported provenance (SimulationOutput.physiological_
    # personalization) consistent with what was actually simulated --
    # otherwise a caller inspecting personalized_params.secretion_scale
    # would see a stale, pre-adjustment value.
    if physio_result is not None:
        physio_result = dataclasses.replace(physio_result, personalized_params=params)
    return params, physio_result


def _assemble_meals(scenario: ScenarioDefinition) -> tuple[list, list[str]]:
    """Stage 4: MealInputs (body-weight adjusted) plus meal warnings."""
    body_weight = 70.0
    if scenario.personalization.physiological_profile is not None:
        body_weight = scenario.personalization.physiological_profile.body_weight_kg

    meals = [
        event.to_meal_input(body_weight_kg=body_weight)
        for event in scenario.meal_events
    ]
    meals, agi_warnings = _apply_alpha_glucosidase_inhibitor(scenario, meals)
    return meals, meal_composition_warnings(scenario.meal_events) + agi_warnings


def _apply_alpha_glucosidase_inhibitor(
    scenario: ScenarioDefinition, meals: list,
) -> tuple[list, list[str]]:
    # Alpha-glucosidase inhibitors (acarbose) slow intestinal
    # carbohydrate absorption -- modelled as a reduction of the meal's
    # absorption rate constant(s): kabs for a single-exponential meal, both
    # k_emptying and k_absorption for a two-compartment meal (see below).
    # Scaling a rate constant preserves total carbohydrate mass exactly
    # (integral(D*k*exp(-k*t) dt) = D for any k>0), unlike scaling
    # Ra(t)'s output (the GLP-1 mass-conservation bug fixed earlier
    # in this file's own history). See biosim_physiology
    # .drug_interactions.acarbose_kabs_scale()'s own docstring for
    # the citation and calibration.
    #
    # Two-compartment meals -- the default meal (since 2026-11-15) and
    # glycaemic_index meals (since 2026-11-17; before that they were left
    # unaffected with a warning) -- are slowed by time-stretching their
    # curve: BOTH rate constants scaled by the same factor, with a warning
    # that the factor was calibrated on the single-exponential curve.
    # (MASTER_PLAN row 1.4b(i) had rejected scaling k_absorption ALONE,
    # which cannot reach the measured effect; scaling both avoids that
    # k_absorption > k_emptying limit.) Free glucose needs no digestion,
    # so set complex_carb_fraction=0 for a pure-glucose drink.
    #
    # Per-meal scaled by MealEvent.complex_carb_fraction (default
    # 1.0, preserving the original milestone's unconditional-full-
    # effect behaviour exactly) -- acarbose inhibits the enzymatic
    # digestion step (alpha-glucosidase), which free monosaccharides
    # (glucose, fructose) don't need at all, so a meal that's mostly
    # free sugar should see a correspondingly smaller effect than one
    # that's mostly starch/disaccharide, not the drug's full
    # modelled effect regardless of composition (MASTER_PLAN row
    # 1.4b(iii)).
    from biosim_physiology.drug_interactions import (
        ALPHA_GLUCOSIDASE_INHIBITOR_DRUGS,
        Drug,
        acarbose_kabs_scale,
        miglitol_kabs_scale,
    )
    # Per-drug dispatch, not a single shared calibration -- acarbose
    # and miglitol each have their own, separately-verified kabs-scale
    # function (different reference studies, different magnitudes).
    # Reusing one drug's calibration for the other here would have
    # been a real, silent mismatch (found while adding miglitol,
    # MASTER_PLAN 1.4b(ii) -- this dispatch previously had only one
    # member so the bug was latent, never actually reachable).
    _AGI_KABS_SCALE_FUNCTIONS = {
        Drug.ACARBOSE: acarbose_kabs_scale,
        Drug.MIGLITOL: miglitol_kabs_scale,
    }
    # If more than one alpha-glucosidase-inhibitor event is present
    # (acarbose + miglitol co-prescribed -- unusual and clinically
    # implausible, these are not typically combined), the first one
    # listed wins -- same "first wins" convention already documented
    # for _detect_glp1()/_detect_renal()/_detect_exogenous_insulin()
    # in pathway_assembly.py. Verified via
    # test_acarbose_and_miglitol_together_only_the_first_listed_is_applied
    # (found missing by the mandatory security-reviewer pass on the
    # miglitol milestone, MASTER_PLAN 1.4b(ii)).
    acarbose_event = None
    agi_drug = None
    for event in scenario.drug_events:
        try:
            drug = Drug(event.drug)
        except ValueError:
            continue
        if drug in ALPHA_GLUCOSIDASE_INHIBITOR_DRUGS:
            acarbose_event = event
            agi_drug = drug
            break
    if acarbose_event is None:
        return meals, []

    warnings: list[str] = []
    dose_kabs_scale = _AGI_KABS_SCALE_FUNCTIONS[agi_drug](acarbose_event.dose_mg)
    new_meals = []
    for meal, event in zip(meals, scenario.meal_events):
        effective_scale = 1.0 - (1.0 - dose_kabs_scale) * event.complex_carb_fraction
        if meal.absorption_model is not None:
            # Two-compartment meal: the default meal (since 2026-11-15) or a
            # glycaemic_index meal (since 2026-11-17, previously left
            # unaffected with a warning). Slowing BOTH
            # rate constants by the drug's factor stretches the curve in time
            # (mass conserved, k_absorption > k_emptying kept), the same thing
            # scaling kabs does to the single-exponential curve. The factor was
            # calibrated on that single-exponential curve: checked 2026-11-15,
            # at 100 mg acarbose the glucose-rise reduction is identical
            # (ratio 0.730 on both curves) but the peak delay is larger
            # (x2.08 here vs x1.51, against ~x1.5 measured). Said in the warning.
            if effective_scale < 1.0:
                model = meal.absorption_model
                new_meals.append(dataclasses.replace(meal, absorption_model=type(model)(
                    k_emptying=model.k_emptying * effective_scale,
                    k_absorption=model.k_absorption * effective_scale,
                )))
                warnings.append(
                    f"Meal '{event.label or event.carbs_g}': {agi_drug.value}'s slowing "
                    "of absorption was calibrated on the single-exponential curve. On "
                    "this two-compartment meal it lowers the glucose rise by the same "
                    "amount but delays the peak more than measured (about 2x instead "
                    "of 1.5x at 100 mg acarbose)."
                )
            else:
                new_meals.append(meal)
        else:
            new_meals.append(dataclasses.replace(meal, kabs=meal.kabs * effective_scale))
    return new_meals, warnings


