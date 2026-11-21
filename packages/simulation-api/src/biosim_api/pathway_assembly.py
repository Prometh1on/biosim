# SAVE TO: packages\simulation-api\src\biosim_api\pathway_assembly.py
"""
Extension-pathway assembly and detection for SimulationRunner.

Moved verbatim out of runner.py (2026-10-20) to keep that module readable;
SimulationRunner inherits these methods, so behaviour and call sites
(self._detect_glp1(...), SimulationRunner._build_pathway, ...) are unchanged.
"""

from __future__ import annotations

from typing import Mapping, Optional

from biosim_ode import (
    GlucoseInsulinPathway,
)
from biosim_ode.alcohol_hepatic_pathway import (
    AlcoholHepaticParameters,
    AlcoholHepaticPathway,
    alcohol_hepatic_couplings,
)
from biosim_ode.cortisol_pathway import (
    CortisolParameters,
    CortisolPathway,
    cortisol_couplings,
)
from biosim_ode.coupling import CompositePathway
from biosim_ode.exercise_pathway import ExercisePathway, exercise_couplings
from biosim_ode.exogenous_insulin_pathway import (
    ExogenousInsulinParameters,
    ExogenousInsulinPathway,
    exogenous_insulin_couplings,
)
from biosim_ode.glp1_pathway import GLP1Pathway, glp1_couplings
from biosim_ode.hepatic import EGPParameters, HepGlucPathway, hepatic_couplings
from biosim_ode.meal_absorption_pathway import (
    MealAbsorptionPathway,
    meal_absorption_couplings,
)
from biosim_ode.nefa_pathway import NEFAPathway, nefa_couplings
from biosim_ode.nutrient_secretion_pathway import (
    NutrientSecretionParameters,
    NutrientSecretionPathway,
    nutrient_secretion_couplings,
)
from biosim_ode.renal_pathway import (
    RenalClearanceParameters,
    RenalPathway,
    renal_couplings,
)

from .scenario import ScenarioDefinition


class PathwayAssemblyMixin:
    """
    Builds the composite ODE pathway for a scenario and detects which
    extension pathways (GLP-1, exogenous insulin, alcohol, renal SGLT2i,
    cortisol, hepatic EGP) its drugs and profile require.
    """

    # ------------------------------------------------------------------
    # Extension pathway composition
    # ------------------------------------------------------------------

    def _build_pathway(
        self,
        scenario: ScenarioDefinition,
        gi_pathway: GlucoseInsulinPathway,
        extension_overrides: Optional[Mapping[str, float]] = None,
    ):
        """
        Detect every extension pathway this scenario activates and
        flat-compose them with gi_pathway into one CompositePathway.
        Returns the bare gi_pathway unchanged if none are active.

        extension_overrides (MASTER_PLAN 2.2): optional dotted-key ->
        sampled-value overrides (e.g. {"renal.beta": 0.42}), forwarded to
        whichever _detect_*() helper owns that key. Defaults to {} (every
        extension uses its own point-estimate default), so every existing
        caller of this method -- the point-estimate build two calls
        below, and MonteCarloEngine.propagate()'s pathway_factory path
        before extension_pathway_factory existed -- is unaffected.
        """
        extension_overrides = extension_overrides or {}
        extensions = []   # list of (pathway_instance, couplings)

        glp1 = self._detect_glp1(scenario.drug_events, basal_glucose_mM=gi_pathway.params.Gb)
        if glp1 is not None:
            extensions.append((glp1, glp1_couplings(gi_pathway.name, glp1.name)))
            # Mass-conserving gastric-emptying coupling (see
            # meal_absorption_pathway.py) -- composed only when a
            # gastric-emptying-modifying drug is actually active, so every
            # other scenario keeps using GI's original closed-form Ra(t)
            # unchanged. glp1's own gastric_emptying_modifier output feeds
            # both GI (unchanged coupling above) and this pathway.
            absorption = MealAbsorptionPathway(gi_pathway.meals)
            extensions.append((
                absorption,
                meal_absorption_couplings(
                    glp1.name, absorption.name, gi_pathway.name, gi_pathway.meals,
                ),
            ))

        cortisol = self._detect_cortisol(scenario)
        if cortisol is not None:
            extensions.append((cortisol, cortisol_couplings(gi_pathway.name, cortisol.name)))

        hepatic = self._detect_hepatic(scenario)
        if hepatic is not None:
            extensions.append((hepatic, hepatic_couplings(gi_pathway.name, hepatic)))

        if scenario.model_nefa_dynamics:
            nefa = NEFAPathway()
            extensions.append((nefa, nefa_couplings(gi_pathway.name, nefa.name)))

        if scenario.exercise_bout_in_progress is not None:
            exercise = ExercisePathway(scenario.exercise_bout_in_progress)
            extensions.append((exercise, exercise_couplings(gi_pathway.name, exercise.name)))

        renal = self._detect_renal(scenario, extension_overrides)
        if renal is not None:
            extensions.append((renal, renal_couplings(gi_pathway.name, renal)))

        exogenous_insulin = self._detect_exogenous_insulin(
            scenario.drug_events, scenario.insulin_dose_events,
        )
        if exogenous_insulin is not None:
            extensions.append((
                exogenous_insulin,
                exogenous_insulin_couplings(gi_pathway.name, exogenous_insulin.name),
            ))

        # Dietary protein co-ingestion (MASTER_PLAN row 1.28) -- composed
        # only when at least one meal actually has protein_g > 0, so every
        # existing scenario (no meal sets this field) is completely
        # unaffected, same "confined blast radius" pattern as every other
        # opt-in extension pathway above.
        protein_params = [
            NutrientSecretionParameters(protein_g=m.protein_g, meal_time_min=m.time_min)
            for m in gi_pathway.meals if m.protein_g > 0
        ]
        if protein_params:
            nutrient = NutrientSecretionPathway(protein_params)
            extensions.append((nutrient, nutrient_secretion_couplings(gi_pathway.name, nutrient.name)))

        alcohol_doses = self._detect_alcohol(scenario.drug_events)
        if alcohol_doses:
            alcohol = AlcoholHepaticPathway(alcohol_doses)
            extensions.append((alcohol, alcohol_hepatic_couplings(gi_pathway.name, alcohol.name)))

        if not extensions:
            return gi_pathway

        all_pathways = [gi_pathway] + [ext for ext, _ in extensions]
        all_couplings = [c for _, couplings in extensions for c in couplings]
        return CompositePathway(all_pathways, couplings=all_couplings)

    def _detect_glp1(
        self, drug_events: list, basal_glucose_mM: float = 5.0,
    ) -> Optional[GLP1Pathway]:
        """
        Detect a GLP-1 class drug event (semaglutide, liraglutide, a DPP-4
        inhibitor) and build the matching GLP1Pathway instance for ODE
        coupling. Returns None if no such drug is present in the scenario.

        If more than one GLP-1 drug event is present (an unusual and
        clinically implausible combination), the first one wins -- this
        platform does not yet model interactions between multiple GLP-1
        class drugs in the same scenario.

        basal_glucose_mM: the ACTUAL composed GI pathway's own basal
        glucose (gi_pathway.params.Gb) -- forwarded so GLP1Pathway's own
        y0 is the genuine fixed point of its ODE at the real basal state
        it will be composed with, not always assuming a healthy 5.0mM
        basal. Defaults to 5.0 (this pathway's own historical default)
        for any caller that doesn't have a specific GI to hand.
        Mandatory security-reviewer finding on the Orskov-1994
        meal-response recalibration: without this, a disease-elevated-Gb
        preset (T2D ~6.5mM, GCK-MODY ~7.17mM) composed with a GLP-1
        agonist drifted severely away from its own "basal" state (e.g.
        T2D+semaglutide 1mg: 142.0 -> 3352.0 pM over a standard 5h
        fasted window) since glucose_stim_basal is genuinely nonzero at
        those Gb values but GLP1Pathway's own y0 never accounted for it.
        """
        from biosim_physiology.drug_interactions import GLP1_DRUGS, Drug

        for event in drug_events:
            try:
                drug = Drug(event.drug)
            except ValueError:
                continue
            if drug not in GLP1_DRUGS:
                continue
            if drug is Drug.SEMAGLUTIDE:
                return GLP1Pathway.semaglutide(
                    dose_mg=event.dose_mg, basal_glucose_mM=basal_glucose_mM,
                )
            if drug is Drug.LIRAGLUTIDE:
                return GLP1Pathway.liraglutide(
                    dose_mg=event.dose_mg, basal_glucose_mM=basal_glucose_mM,
                )
            if drug is Drug.DPP4_INHIBITOR:
                return GLP1Pathway.dpp4_inhibited(basal_glucose_mM=basal_glucose_mM)
        return None

    def _detect_exogenous_insulin(
        self, drug_events: list, insulin_dose_events: Optional[list] = None,
    ) -> Optional[ExogenousInsulinPathway]:
        """
        Detect an exogenous-insulin dose (regular human insulin or the
        rapid-acting analog insulin lispro) and build the matching
        ExogenousInsulinPathway instance for ODE coupling. Returns None
        if none is present. No static PathwayParameters modifier is
        applied for either insulin type at all (see
        _exogenous_insulin_modifiers()'s docstring in
        biosim_physiology.drug_interactions) -- this is their sole
        representation.

        Two ways to express a dose, checked in this order (MASTER_PLAN
        1.2b): `insulin_dose_events` (InsulinDoseEvent, genuine Unit
        semantics, the preferred form) takes priority if present;
        otherwise falls back to the legacy `drug_events`-based path
        (a DrugEvent for Drug.INSULIN_REGULAR/Drug.INSULIN_LISPRO, whose
        dose_mg is silently reinterpreted as Units -- a disclosed
        unit-reuse hack kept only for backward compatibility, see
        DrugEvent.dose_mg's own docstring). `ScenarioDefinition
        .__post_init__` already rejects setting both for the same
        scenario, so this priority order only ever matters for which
        (empty) path is skipped, never for resolving an actual conflict.

        Unlike every other _detect_*() helper in this file, this pathway
        is genuinely TIME-VARYING within the simulated window (see
        ExogenousInsulinPathway's own docstring) -- dose_time_min combines
        BOTH event.time_s (a dose given partway through the simulated
        window) and event.hours_before_simulation (a dose already given
        before t=0) into one net time offset relative to t=0, rather than
        using DrugEvent.to_drug_dose()'s hours_since_ingestion (which
        discards time_s entirely for the standard simulation_start_s=0.0
        call every other, STATIC drug-modifier path in this file uses --
        appropriate for a one-shot personalisation snapshot, not for a
        signal that must know exactly when within the window the dose
        was actually given).

        If more than one exogenous-insulin dose is present (multiple
        injections within one scenario -- clinically real, e.g. a
        correction bolus), the first one wins for now -- this platform
        does not yet model multiple overlapping insulin doses in the same
        scenario (a real, disclosed limitation, not silently wrong: their
        combined effect would need to be superposed, not just the first
        one kept). This applies within each path independently (the
        first InsulinDoseEvent, or the first matching DrugEvent).
        """
        if insulin_dose_events:
            event = insulin_dose_events[0]
            dose_time_min = (event.time_s / 60.0) - (event.hours_before_simulation * 60.0)
            if event.insulin_type == "regular":
                return ExogenousInsulinPathway(
                    ExogenousInsulinParameters.regular_insulin(
                        dose_units=event.dose_units, dose_time_min=dose_time_min,
                    )
                )
            if event.insulin_type == "lispro":
                return ExogenousInsulinPathway(
                    ExogenousInsulinParameters.rapid_acting_insulin(
                        dose_units=event.dose_units, dose_time_min=dose_time_min,
                    )
                )
            return None  # unreachable -- InsulinDoseEvent's own __post_init__ validates insulin_type

        from biosim_physiology.drug_interactions import EXOGENOUS_INSULIN_DRUGS, Drug

        for event in drug_events:
            try:
                drug = Drug(event.drug)
            except ValueError:
                continue
            if drug not in EXOGENOUS_INSULIN_DRUGS:
                continue
            dose_time_min = (event.time_s / 60.0) - (event.hours_before_simulation * 60.0)
            if drug is Drug.INSULIN_REGULAR:
                return ExogenousInsulinPathway(
                    ExogenousInsulinParameters.regular_insulin(
                        dose_units=event.dose_mg, dose_time_min=dose_time_min,
                    )
                )
            if drug is Drug.INSULIN_LISPRO:
                return ExogenousInsulinPathway(
                    ExogenousInsulinParameters.rapid_acting_insulin(
                        dose_units=event.dose_mg, dose_time_min=dose_time_min,
                    )
                )
        return None

    def _detect_alcohol(self, drug_events: list) -> list[AlcoholHepaticParameters]:
        """
        Detect every Drug.ALCOHOL event and build the matching
        AlcoholHepaticParameters list for ODE coupling (see
        AlcoholHepaticPathway's own module docstring for why this
        replaced the previous static Gb modifier). Unlike GLP-1/renal/
        exogenous-insulin's "first event wins" convention, ALL alcohol
        events are collected -- AlcoholHepaticPathway sums each dose's
        own transient contribution (a real, plausible scenario: more
        than one drink at different times), matching
        NutrientSecretionPathway's own multi-event convention.

        dose_time_min combines BOTH event.time_s and
        event.hours_before_simulation into one net time offset relative
        to t=0, same technique as _detect_exogenous_insulin() -- alcohol
        is genuinely time-varying within the simulated window, so the
        static-modifier-path's hours_since_ingestion (which discards
        time_s) is not appropriate here.
        """
        from biosim_physiology.drug_interactions import Drug

        doses = []
        for event in drug_events:
            try:
                drug = Drug(event.drug)
            except ValueError:
                continue
            if drug is not Drug.ALCOHOL:
                continue
            dose_time_min = (event.time_s / 60.0) - (event.hours_before_simulation * 60.0)
            doses.append(
                AlcoholHepaticParameters(dose_g=event.dose_mg / 1000.0, dose_time_min=dose_time_min)
            )
        return doses

    def _detect_renal(
        self, scenario: ScenarioDefinition, extension_overrides: Optional[Mapping[str, float]] = None,
    ) -> Optional[RenalPathway]:
        """
        Detect an SGLT2-inhibitor drug event (dapagliflozin, empagliflozin)
        and build the matching RenalPathway instance for ODE coupling.
        Returns None if no such drug is present in the scenario. No static
        PathwayParameters modifier is applied for these drugs at all (see
        _sglt2_inhibitor_modifiers()'s docstring in
        biosim_physiology.drug_interactions) -- this is their sole
        representation.

        Also detects a type-2-diabetes disease state, the same way
        _detect_hepatic() does, and passes it through to each drug's
        for_dapagliflozin()/for_empagliflozin() classmethod so the
        pathway uses that drug's own cited T2D-cohort renal-transport
        values (elevated baseline TmG) rather than always defaulting to
        the healthy-cohort value regardless of the patient's actual
        disease state -- see both classmethods' own docstrings for the
        exact figures and citations (Track 1 MASTER_PLAN 1.8).

        extension_overrides (MASTER_PLAN 2.2): if the dotted key
        "renal.beta" is present, forwarded as both classmethods'
        beta_class_default -- lets a Monte Carlo UQ sample override the
        class-wide beta point estimate (see RENAL_BETA_CV) rather than
        always using SGLT2I_BETA_CLASS_DEFAULT regardless of context.
        Defaults to {} (no override, exact prior behaviour).

        If more than one SGLT2i drug event is present (clinically
        implausible -- these are not co-prescribed), the first one wins,
        same convention as _detect_glp1().
        """
        from biosim_physiology.drug_interactions import SGLT2_DRUGS, Drug
        from biosim_physiology.profile import DiseaseState

        extension_overrides = extension_overrides or {}
        beta_kwargs = {}
        if "renal.beta" in extension_overrides:
            beta_kwargs["beta_class_default"] = extension_overrides["renal.beta"]

        profile = scenario.personalization.physiological_profile
        disease_states = profile.disease_states if profile is not None else []
        is_t2d = DiseaseState.TYPE_2_DIABETES in disease_states

        for event in scenario.drug_events:
            try:
                drug = Drug(event.drug)
            except ValueError:
                continue
            if drug not in SGLT2_DRUGS:
                continue
            if drug is Drug.DAPAGLIFLOZIN:
                return RenalPathway(
                    RenalClearanceParameters.for_dapagliflozin(t2d=is_t2d, **beta_kwargs)
                )
            if drug is Drug.EMPAGLIFLOZIN:
                days_of_treatment = event.hours_before_simulation / 24.0
                return RenalPathway(
                    RenalClearanceParameters.for_empagliflozin(
                        days_of_treatment, t2d=is_t2d, **beta_kwargs
                    )
                )
        return None

    @staticmethod
    def _glucocorticoid_event(scenario: ScenarioDefinition):
        """The first Drug.GLUCOCORTICOID DrugEvent in the scenario, or None."""
        for event in scenario.drug_events:
            if event.drug == "glucocorticoid":
                return event
        return None

    def _detect_cortisol(self, scenario: ScenarioDefinition) -> Optional[CortisolPathway]:
        """
        Detect scenario.stress_level > 0 and/or a glucocorticoid (steroid)
        drug event, and build one combined CortisolPathway instance for
        both -- a patient can be stressed AND on prednisone simultaneously,
        and the HPA axis is a single physiological system, not two
        independent instances. Returns None if neither is present.

        Regression: this docstring used to claim the exogenous-steroid
        contribution was "additive to the existing static
        glucocorticoid_modifiers()... not a replacement for it." That was
        wrong and unverified -- confirmed via a fasted, zero-meal
        SimulationRunner run that composing both produced an extra +18%
        fasting-glucose rise (and near-doubled steady-state insulin) beyond
        either alone. CortisolPathway's gb_modifier/p3_modifier have no
        meal/insulin input at all (see cortisol_couplings()'s docstring:
        "cortisol is self-driven... not by GI state") and are a complete,
        standalone representation of the SAME Gb/p3 effect
        glucocorticoid_modifiers() computes -- genuinely a REPLACEMENT for
        it, not a complementary addition (unlike hepatic's postprandial EGP
        kinetics, which really is complementary to metformin's static
        effect, since it responds to meal-driven insulin). `run()` excludes
        glucocorticoid drug events from PhysiologicalPersonalizer's static
        drug_doses for exactly this reason, applying only
        glucocorticoid_secretion_scale_factor()'s compensatory/decompensatory
        secretion_scale bump separately (the one piece CortisolPathway
        doesn't represent at all).

        Real dosing timing: `event.time_s`/`event.hours_before_simulation`
        are converted to the ODE-integration-time coordinate (minutes)
        CortisolParameters.for_prednisone()'s dose_time_min expects, the
        same conversion DrugEvent.to_drug_dose() uses internally -- unlike
        that method (which collapses to a single "hours elapsed as of
        t=0" snapshot), this preserves the sign for a dose scheduled
        DURING the simulation (time_s > 0), since CortisolPathway's
        exogenous_gc_nM_at() is evaluated at the ODE solver's real,
        continuously-varying t, not just at t=0.

        Track 1 row 1.38: scenario.cushings_syndrome composes into this
        SAME combined instance too, for the same reason stress+steroid
        already do (one HPA axis, not independent instances) -- Cushing's
        raises the TRUE basal cortisol/GR-occupancy (via
        solve_cortisol_basal_for_target_steady_state(), the same inverse
        solve CortisolParameters.for_cushings() uses on its own), on top
        of which acute stress and/or an exogenous steroid dose can still
        stack, all sharing the one saturating GR-ligand pool.
        """
        steroid_event = self._glucocorticoid_event(scenario)

        if scenario.stress_level <= 0 and steroid_event is None \
                and not scenario.cushings_syndrome:
            return None

        cushings_kwargs = {}
        if scenario.cushings_syndrome:
            from biosim_ode.cortisol_pathway import (
                CORTISOL_BASAL_NM,
                CRH_BASAL_PM,
                solve_cortisol_basal_for_target_steady_state,
            )
            target_cortisol_ss = CORTISOL_BASAL_NM * CortisolParameters.CUSHINGS_CORTISOL_RATIO
            cushings_kwargs = {
                "cortisol_basal": solve_cortisol_basal_for_target_steady_state(
                    target_cortisol_ss, CRH_BASAL_PM, scenario.stress_level,
                ),
            }

        if steroid_event is not None:
            dose_time_min = (
                steroid_event.time_s / 60.0 - steroid_event.hours_before_simulation * 60.0
            )
            params = CortisolParameters.for_prednisone(
                steroid_event.dose_mg,
                dose_time_min=dose_time_min,
                stress_level=scenario.stress_level,
                **cushings_kwargs,
            )
        elif scenario.cushings_syndrome:
            # Reuses the cortisol_basal already solved above rather than
            # calling for_cushings() again, which would re-derive the same
            # target_cortisol_ss/brentq solve a second time (code-reviewer
            # LOW finding) -- for_cushings() itself is exactly this
            # construction (see its own docstring: "does NOT pass
            # gr_occ_basal -- basal_state() computes the genuine coupled
            # fixed point itself from cortisol_basal alone").
            params = CortisolParameters(
                cortisol_basal=cushings_kwargs["cortisol_basal"], stress_level=scenario.stress_level,
            )
        else:
            params = CortisolParameters(stress_level=scenario.stress_level)
        return CortisolPathway(params)

    def _detect_hepatic(self, scenario: ScenarioDefinition) -> Optional[HepGlucPathway]:
        """
        Detect type-2-diabetes disease state, hyperthyroidism disease state,
        and/or a metformin drug event, and build the matching HepGlucPathway
        for ODE coupling. Combines T2D + metformin if both present (a T2D
        patient on metformin) into one EGPParameters instance rather than
        picking only one classmethod's effect.

        T2D/metformin's hepatic effect is additive to (not a replacement
        for) the existing static Gb/p3/secretion_scale modifiers in
        physiology.modifiers -- see hepatic_couplings() docstring for why
        the two do not double-count. Hyperthyroidism is the OPPOSITE case
        -- it is the SOLE representation (see DiseaseState.HYPERTHYROIDISM
        in physiology/modifiers.py, which deliberately has no static
        component), so it is handled independently below, not folded into
        the T2D/metformin combination logic.

        Scope limitation, documented rather than silently dropped: a
        combined hyperthyroidism + T2D/metformin preset has not been built
        (no combined-condition literature was scoped for it this session,
        unlike T2D+metformin's for_type2_diabetes_and_metformin()) --
        hyperthyroidism is only applied when T2D and metformin are BOTH
        absent. A scenario with hyperthyroidism alongside T2D/metformin
        still runs successfully, but only reflects the T2D/metformin
        hepatic effect, not hyperthyroidism's. This is a real coverage
        gap to close in a future milestone if that combination is needed,
        not a bug being papered over.
        """
        from biosim_physiology.drug_interactions import Drug
        from biosim_physiology.profile import DiseaseState

        profile = scenario.personalization.physiological_profile
        disease_states = profile.disease_states if profile is not None else []
        is_t2d = DiseaseState.TYPE_2_DIABETES in disease_states
        is_hyperthyroid = DiseaseState.HYPERTHYROIDISM in disease_states

        metformin_dose = None
        for event in scenario.drug_events:
            try:
                drug = Drug(event.drug)
            except ValueError:
                continue
            if drug is Drug.METFORMIN:
                metformin_dose = event.dose_mg
                break

        if is_t2d and metformin_dose is not None:
            return HepGlucPathway(EGPParameters.for_type2_diabetes_and_metformin(metformin_dose))
        if is_t2d:
            return HepGlucPathway(EGPParameters.for_type2_diabetes())
        if metformin_dose is not None:
            return HepGlucPathway(EGPParameters.for_metformin(metformin_dose))
        if is_hyperthyroid:
            return HepGlucPathway(EGPParameters.for_hyperthyroidism())
        return None
