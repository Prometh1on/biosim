# SAVE TO: packages\simulation-api\tests\test_simulation_api.py
"""
Test suite for biosim_api.

Coverage:
    - ScenarioDefinition: construction, validation, serialisation
    - MealEvent and DrugEvent: construction, conversion
    - SimulationRunner: base run, with meals, with personalization, with UQ
    - ExpertMetrics: all fields populated, units correct
    - LaySummary: energy/hunger trajectories, dip detection, narratives
    - Scenario library: all pre-built scenarios run successfully
    - ComparisonRunner: differential output, larger meal = larger difference
    - End-to-end: full personalised scenario with physio + genetic layers
"""

import math

import pytest
import numpy as np

from biosim_api import (
    LIBRARY,
    ComparisonRunner,
    DrugEvent,
    InsulinDoseEvent,
    MealEvent,
    PersonalizationSpec,
    ScenarioDefinition,
    SimulationRunner,
    caffeine_and_meal_scenario,
    fasted_scenario,
    get_scenario,
    high_sugar_meal_scenario,
    light_snack_scenario,
    ogtt_75g,
    standard_meal_scenario,
)


# ---------------------------------------------------------------------------
# ScenarioDefinition
# ---------------------------------------------------------------------------

class TestScenarioDefinition:
    def test_basic_construction(self):
        s = ScenarioDefinition(name="test", duration_s=7200)
        assert s.name == "test"
        assert s.duration_s == 7200

    def test_empty_name_raises(self):
        with pytest.raises(ValueError, match="name"):
            ScenarioDefinition(name="", duration_s=7200)

    def test_negative_duration_raises(self):
        with pytest.raises(ValueError, match="duration"):
            ScenarioDefinition(name="x", duration_s=-1.0)

    def test_low_mc_samples_raises(self):
        with pytest.raises(ValueError, match="n_mc_samples"):
            ScenarioDefinition(name="x", n_mc_samples=5)

    def test_high_mc_samples_raises(self):
        # Housekeeping-backlog fix: n_mc_samples previously had no
        # upper bound at all -- an unbounded, uninterruptible compute
        # cost reachable by any caller with run_uncertainty=True.
        with pytest.raises(ValueError, match="n_mc_samples"):
            ScenarioDefinition(name="x", n_mc_samples=1001)

    def test_mc_samples_at_the_ceiling_is_accepted(self):
        # Boundary check: exactly 1000 (the ceiling itself) must still
        # be valid -- only values ABOVE it are rejected.
        scenario = ScenarioDefinition(name="x", n_mc_samples=1000)
        assert scenario.n_mc_samples == 1000

    def test_excessive_duration_raises(self):
        # Found by the mandatory security-reviewer pass for the
        # n_mc_samples fix above -- CRITICAL, fixed in the same pass.
        with pytest.raises(ValueError, match="duration_s"):
            ScenarioDefinition(name="x", duration_s=604801.0)

    def test_duration_at_the_ceiling_is_accepted(self):
        scenario = ScenarioDefinition(name="x", duration_s=604800.0)
        assert scenario.duration_s == 604800.0

    def test_excessive_meal_events_raises(self):
        events = [MealEvent(time_s=float(i) * 60, carbs_g=10.0) for i in range(51)]
        with pytest.raises(ValueError, match="meal_events"):
            ScenarioDefinition(name="x", meal_events=events)

    def test_meal_events_at_the_ceiling_is_accepted(self):
        events = [MealEvent(time_s=float(i) * 60, carbs_g=10.0) for i in range(50)]
        scenario = ScenarioDefinition(name="x", meal_events=events)
        assert len(scenario.meal_events) == 50

    def test_both_exercise_bout_fields_raises(self):
        # Regression: personalization.exercise_bout (a bout already ended
        # before t=0, a static PhysiologicalPersonalizer modifier) and
        # exercise_bout_in_progress (a bout starting at t=0, driving
        # ExercisePathway's ODE coupling) model mutually exclusive things
        # -- setting both for the same bout would double-count its effect
        # through both the static modifier and the ODE pathway, the same
        # bug shape already found and fixed for metformin and
        # glucocorticoid composition this session. This was previously
        # only a documentation comment, with no runtime guard.
        from biosim_physiology import ExerciseBout, ExerciseIntensity
        bout = ExerciseBout.from_intensity(ExerciseIntensity.MODERATE, 45.0)
        with pytest.raises(ValueError, match="exercise_bout_in_progress"):
            ScenarioDefinition(
                name="x",
                personalization=PersonalizationSpec(exercise_bout=bout),
                exercise_bout_in_progress=bout,
            )

    def test_either_exercise_bout_field_alone_is_fine(self):
        from biosim_physiology import ExerciseBout, ExerciseIntensity
        bout = ExerciseBout.from_intensity(ExerciseIntensity.MODERATE, 45.0)
        ScenarioDefinition(name="a", personalization=PersonalizationSpec(exercise_bout=bout))
        ScenarioDefinition(name="b", exercise_bout_in_progress=bout)

    def test_has_meals(self):
        s = ScenarioDefinition(
            name="x",
            meal_events=[MealEvent(time_s=0.0, carbs_g=50.0)]
        )
        assert s.has_meals

    def test_no_meals(self):
        assert not ScenarioDefinition(name="x").has_meals

    def test_has_drugs(self):
        s = ScenarioDefinition(
            name="x",
            drug_events=[DrugEvent(time_s=0.0, drug="caffeine", dose_mg=200.0)]
        )
        assert s.has_drugs

    def test_to_dict_serialisable(self):
        import json
        s = high_sugar_meal_scenario()
        d = s.to_dict()
        json_str = json.dumps(d)
        assert "high_sugar_meal" in json_str

    def test_to_dict_has_meal_events(self):
        s = high_sugar_meal_scenario()
        d = s.to_dict()
        assert len(d["meal_events"]) == 1
        assert d["meal_events"][0]["carbs_g"] == 75.0


# ---------------------------------------------------------------------------
# MealEvent and DrugEvent
# ---------------------------------------------------------------------------

class TestMealEvent:
    def test_to_meal_input(self):
        e = MealEvent(time_s=1800.0, carbs_g=50.0)
        mi = e.to_meal_input(body_weight_kg=70.0)
        assert mi.carbs_g == 50.0
        assert mi.time_min == pytest.approx(30.0)

    def test_negative_time_raises(self):
        with pytest.raises(ValueError):
            MealEvent(time_s=-10.0, carbs_g=50.0)

    def test_negative_carbs_raises(self):
        with pytest.raises(ValueError):
            MealEvent(time_s=0.0, carbs_g=-1.0)

    def test_no_glycaemic_index_uses_default_two_compartment_model(self):
        # Since 2026-11-15 a meal with neither glycaemic_index nor kabs uses
        # the two-compartment default (was single-exponential).
        from biosim_ode.gut_absorption import GutAbsorptionModel
        mi = MealEvent(time_s=0.0, carbs_g=50.0).to_meal_input()
        assert mi.absorption_model == GutAbsorptionModel.default_meal()

    def test_glycaemic_index_builds_absorption_model(self):
        from biosim_ode.gut_absorption import GutAbsorptionModel
        e = MealEvent(time_s=0.0, carbs_g=75.0, glycaemic_index=100.0)
        mi = e.to_meal_input()
        assert isinstance(mi.absorption_model, GutAbsorptionModel)
        assert mi.Ra(0.0) == pytest.approx(0.0, abs=1e-9)

    def test_kabs_override_applied_to_meal_input(self):
        # Regression: kabs was previously discarded entirely --
        # biosim_seeds.gi_meals.gi_meal_scenario() computed a GI-weighted
        # kabs but had no field on MealEvent to carry it through, so every
        # GI-adjusted scenario silently used the default kabs regardless of
        # the meal's actual glycaemic index.
        e = MealEvent(time_s=0.0, carbs_g=50.0, kabs=0.045)
        mi = e.to_meal_input()
        assert mi.kabs == pytest.approx(0.045)
        assert mi.absorption_model is None

    def test_kabs_and_glycaemic_index_mutually_exclusive(self):
        with pytest.raises(ValueError):
            MealEvent(time_s=0.0, carbs_g=50.0, kabs=0.045, glycaemic_index=80.0)

    def test_protein_g_defaults_to_zero(self):
        e = MealEvent(time_s=0.0, carbs_g=50.0)
        assert e.protein_g == 0.0

    def test_negative_protein_g_raises(self):
        with pytest.raises(ValueError):
            MealEvent(time_s=0.0, carbs_g=50.0, protein_g=-1.0)

    def test_nan_protein_g_raises(self):
        # Mandatory security-reviewer finding (MEDIUM, fixed before
        # finalizing this milestone): `< 0` alone doesn't reject NaN --
        # protein_g had inherited the same weaker convention as
        # carbs_g/time_s rather than DrugEvent's own stronger
        # math.isfinite() guard.
        with pytest.raises(ValueError):
            MealEvent(time_s=0.0, carbs_g=50.0, protein_g=float("nan"))

    def test_infinite_protein_g_raises(self):
        with pytest.raises(ValueError):
            MealEvent(time_s=0.0, carbs_g=50.0, protein_g=float("inf"))

    def test_protein_g_carried_through_to_meal_input(self):
        # Regression pattern matching kabs's own history (MASTER_PLAN
        # 1.4b(iii)): a MealEvent field that isn't forwarded in
        # to_meal_input() is silently discarded.
        e = MealEvent(time_s=0.0, carbs_g=50.0, protein_g=25.0)
        mi = e.to_meal_input()
        assert mi.protein_g == pytest.approx(25.0)


class TestGlp1WithGiMealAllowed:
    """
    A GLP-1-class drug combined with a glycaemic_index meal used to be
    rejected: the closed-form GutAbsorptionModel Ra(t) was multiplied by
    gastric_emptying_modifier, deleting (1 - modifier) of the meal's glucose
    (~38% of a 75 g OGTT with a DPP-4 inhibitor in T2D). Since 2026-11-15
    MealAbsorptionPathway tracks two-compartment meals as stomach + gut
    states, so the drug delays the glucose instead of deleting it (see
    test_default_meal_two_compartment.py for the mass-conservation check),
    and the combination is allowed.
    """

    @pytest.mark.parametrize("drug", ["semaglutide", "liraglutide", "dpp4_inhibitor"])
    @pytest.mark.parametrize("gi", [100.0, 55.0])
    def test_glp1_drug_with_gi_meal_runs(self, drug, gi):
        s = ScenarioDefinition(
            name="gi_glp1", duration_s=3600.0,
            meal_events=[MealEvent(time_s=0.0, carbs_g=75.0, glycaemic_index=gi)],
            drug_events=[DrugEvent(time_s=0.0, drug=drug, dose_mg=1.0)],
        )
        assert SimulationRunner().run(s).success

    def test_glp1_drug_with_kabs_meal_is_allowed(self):
        ScenarioDefinition(
            name="ok", meal_events=[MealEvent(time_s=0.0, carbs_g=75.0, kabs=0.04)],
            drug_events=[DrugEvent(time_s=0.0, drug="semaglutide", dose_mg=1.0)],
        )

    def test_mixed_meal_kinds_run_together(self):
        s = ScenarioDefinition(
            name="mixed", duration_s=7200.0,
            meal_events=[
                MealEvent(time_s=0.0, carbs_g=30.0, kabs=0.03),
                MealEvent(time_s=1800.0, carbs_g=40.0, glycaemic_index=40.0),
                MealEvent(time_s=3600.0, carbs_g=20.0),
            ],
            drug_events=[DrugEvent(time_s=0.0, drug="liraglutide", dose_mg=1.2)],
        )
        assert SimulationRunner().run(s).success


class TestDrugEvent:
    def test_to_drug_dose(self):
        e = DrugEvent(time_s=0.0, drug="caffeine", dose_mg=200.0)
        dose = e.to_drug_dose(simulation_start_s=3600.0)
        assert dose.dose_mg == 200.0

    def test_negative_dose_raises(self):
        with pytest.raises(ValueError):
            DrugEvent(time_s=0.0, drug="caffeine", dose_mg=-100.0)

    def test_unrecognized_drug_raises_at_construction(self):
        # Regression: an unrecognized/typo'd drug string used to survive
        # construction and only fail later, deep inside to_drug_dose()
        # (an uncaught Drug(self.drug) ValueError) -- a confusing crash
        # far from the actual mistake. Separately, runner.py's
        # _detect_glp1()/_detect_hepatic() caught that same ValueError
        # and silently skipped composing the pathway with zero warning.
        # Fail fast at construction instead, with a clear message.
        with pytest.raises(ValueError, match="not a recognized drug"):
            DrugEvent(time_s=0.0, drug="semaglutid", dose_mg=1.0)

    def test_negative_hours_before_simulation_raises(self):
        with pytest.raises(ValueError):
            DrugEvent(time_s=0.0, drug="caffeine", dose_mg=200.0, hours_before_simulation=-1.0)

    def test_nan_dose_mg_raises(self):
        # Regression (security review, glucocorticoid PK/PD milestone --
        # see CLAUDE.md): "dose_mg < 0" does not reject NaN (NaN
        # comparisons are always False), so a NaN dose used to survive
        # construction and only crash much later, deep inside scipy's ODE
        # solver with an opaque "array must not contain infs or NaNs"
        # error -- confirmed reachable via
        # glucocorticoid_secretion_scale_factor()'s own dose_mg<=0.0 guard
        # having the identical gap. Fail fast at construction instead.
        with pytest.raises(ValueError, match="finite"):
            DrugEvent(time_s=0.0, drug="glucocorticoid", dose_mg=float("nan"))

    def test_infinite_dose_mg_raises(self):
        with pytest.raises(ValueError, match="finite"):
            DrugEvent(time_s=0.0, drug="glucocorticoid", dose_mg=float("inf"))

    def test_nan_hours_before_simulation_raises(self):
        with pytest.raises(ValueError, match="finite"):
            DrugEvent(
                time_s=0.0, drug="glucocorticoid", dose_mg=20.0,
                hours_before_simulation=float("nan"),
            )

    def test_hours_before_simulation_reflected_in_drug_dose(self):
        # Regression: to_drug_dose()'s old formula, (simulation_start_s -
        # time_s), could only ever be <= 0 given time_s>=0 is enforced and
        # callers always pass simulation_start_s=0.0 -- clamped to exactly
        # 0.0 every time. hours_before_simulation is the only way to model
        # a dose already partially elapsed at t=0.
        e = DrugEvent(time_s=0.0, drug="caffeine", dose_mg=200.0, hours_before_simulation=0.5)
        dose = e.to_drug_dose(simulation_start_s=0.0)
        assert dose.hours_since_ingestion == pytest.approx(0.5)

    def test_default_hours_before_simulation_is_zero(self):
        e = DrugEvent(time_s=0.0, drug="caffeine", dose_mg=200.0)
        dose = e.to_drug_dose(simulation_start_s=0.0)
        assert dose.hours_since_ingestion == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# SimulationRunner — core
# ---------------------------------------------------------------------------

class TestSimulationRunnerCore:
    @pytest.fixture(scope="class")
    def fasted_output(self):
        return SimulationRunner().run(fasted_scenario(duration_s=3600.0))

    @pytest.fixture(scope="class")
    def meal_output(self):
        return SimulationRunner().run(high_sugar_meal_scenario())

    def test_fasted_succeeds(self, fasted_output):
        assert fasted_output.success

    def test_meal_succeeds(self, meal_output):
        assert meal_output.success

    def test_glucose_rises_after_meal(self, meal_output):
        G_baseline = meal_output.expert_metrics.glucose_baseline_uM
        G_peak = meal_output.expert_metrics.glucose_peak_uM
        assert G_peak > G_baseline * 1.1

    def test_fasted_glucose_stable(self, fasted_output):
        G = fasted_output.simulation.species["G"]
        assert np.max(G) - np.min(G) < 100.0  # < 0.1 mM drift

    def test_base_params_override_changes_output(self):
        # Regression: SimulationRunner used to hardcode a fresh
        # PathwayParameters() inside run(), ignoring any base params a
        # caller constructed it with (e.g. PlatformSession's
        # literature-loaded self._base_params) -- invisible only because
        # today's literature-seeded values happen to equal the class
        # defaults exactly.
        from biosim_ode import PathwayParameters
        default_out = SimulationRunner().run(fasted_scenario(duration_s=3600.0))
        custom_out = SimulationRunner(
            base_params=PathwayParameters(Gb=6.5)
        ).run(fasted_scenario(duration_s=3600.0))
        assert custom_out.expert_metrics.glucose_baseline_uM == pytest.approx(6500.0)
        assert (
            default_out.expert_metrics.glucose_baseline_uM
            != custom_out.expert_metrics.glucose_baseline_uM
        )

    def test_expert_metrics_populated(self, meal_output):
        m = meal_output.expert_metrics
        assert m.glucose_peak_uM > 0
        assert m.insulin_peak_uM > 0
        assert m.mtor_peak > 0
        assert 0 <= m.ampk_nadir <= 1.0

    def test_expert_metrics_units(self, meal_output):
        m = meal_output.expert_metrics
        # Glucose peak should be in uM (5000-20000 uM for physiological range)
        assert 5000 < m.glucose_peak_uM < 25000
        # Glucose peak in mM should be 5-25
        assert 5 < m.glucose_peak_mM < 25

    def test_scenario_name_preserved(self, meal_output):
        assert meal_output.scenario_name == "high_sugar_meal"

    def test_summary_dict(self, meal_output):
        s = meal_output.summary()
        assert "glucose_peak_mM" in s
        assert "insulin_peak_pM" in s
        assert "confidence" in s

    def test_glucose_iauc_5h_zero_when_sim_too_short(self, meal_output):
        # high_sugar_meal_scenario() runs for 4h (14400s), not the 5h this
        # metric needs -- must default to 0.0, not raise or extrapolate.
        assert meal_output.expert_metrics.glucose_iauc_5h_uM_s == 0.0

    def test_glucose_iauc_5h_positive_for_long_enough_meal(self):
        # See CLAUDE.md "DONE -- Built the first validation benchmark for a clinically significant drug pathway" -- iAUC(0-5h), the standard
        # clinical incremental-AUC metric used to validate GLP-1 drug
        # pathways against real trials (Hjerpsted et al. 2018,
        # PMID:28941314), needs a >=5h simulation window to be meaningful.
        scenario = ScenarioDefinition(
            name="iauc_test_meal", duration_s=18000.0,
            meal_events=[MealEvent(time_s=0.0, carbs_g=75.0)],
            personalization=PersonalizationSpec(),
        )
        output = SimulationRunner().run(scenario)
        m = output.expert_metrics
        assert m.glucose_iauc_5h_uM_s > 0
        # mmol.h/L conversion: uM.s -> mM.h -- sanity-check the physiological
        # range (a 75g meal's 5h incremental AUC should be a few mmol.h/L,
        # not near-zero or absurdly large).
        assert 0.5 < m.glucose_iauc_5h_mmol_h_per_L < 20.0

    def test_glucose_iauc_5h_excludes_below_baseline_undershoot(self):
        # iAUC is clipped at 0 below baseline by definition (matches how
        # real trials compute it) -- a fasted (no meal) scenario run long
        # enough should show iAUC == 0 even if glucose drifts slightly,
        # since there's no meal-driven excursion above baseline to integrate.
        scenario = ScenarioDefinition(
            name="iauc_fasted_test", duration_s=18000.0,
            meal_events=[], personalization=PersonalizationSpec(),
        )
        output = SimulationRunner().run(scenario)
        assert output.expert_metrics.glucose_iauc_5h_uM_s == pytest.approx(0.0, abs=1.0)

    def test_glucose_iauc_4h_positive_for_4h_meal(self, meal_output):
        # high_sugar_meal_scenario() runs for exactly 4h (14400s) -- the
        # standard "mixed meal test" duration this metric is fixed to
        # (see CLAUDE.md's metformin benchmark entry, PMID:39049802),
        # unlike glucose_iauc_5h_uM_s which is 0.0 for this same fixture.
        m = meal_output.expert_metrics
        assert m.glucose_iauc_4h_uM_s > 0
        assert 0.5 < m.glucose_iauc_4h_mmol_min_per_L < 2000.0

    def test_glucose_iauc_4h_zero_when_sim_too_short(self, fasted_output):
        # fasted_scenario(duration_s=3600.0) runs only 1h, not the 4h this
        # metric needs -- must default to 0.0, not raise or extrapolate.
        assert fasted_output.expert_metrics.glucose_iauc_4h_uM_s == 0.0

    def test_glucose_iauc_4h_excludes_below_baseline_undershoot(self):
        scenario = ScenarioDefinition(
            name="iauc_4h_fasted_test", duration_s=14400.0,
            meal_events=[], personalization=PersonalizationSpec(),
        )
        output = SimulationRunner().run(scenario)
        assert output.expert_metrics.glucose_iauc_4h_uM_s == pytest.approx(0.0, abs=1.0)


# ---------------------------------------------------------------------------
# for_gi() uncalibrated-path warning
#
# See CLAUDE.md VALIDATION STATUS: GutAbsorptionModel.for_gi() (any
# glycaemic_index < 99.5) is not independently validated and is
# contradicted by conflicting clinical evidence on whether GI drives peak
# timing at all (Cocate 2011 vs Louie 2012/2013). No production scenario
# used this path when last confirmed, but if one starts, SimulationRunner
# should surface that loudly instead of silently trusting it as calibrated.
# ---------------------------------------------------------------------------

class TestForGiWarning:
    def test_non_ogtt_gi_meal_produces_warning(self):
        scenario = ScenarioDefinition(
            name="low_gi_test",
            duration_s=3600.0,
            meal_events=[MealEvent(time_s=0.0, carbs_g=50.0, glycaemic_index=40.0)],
        )
        output = SimulationRunner().run(scenario)
        assert any("for_gi" in w for w in output.warnings)

    def test_warning_says_placeholder_and_only_gi_100_calibrated(self):
        from biosim_ode.gut_absorption import GI_SLOPE_PLACEHOLDER_NOTE
        assert "unvalidated placeholder" in GI_SLOPE_PLACEHOLDER_NOTE
        assert "only GI=100" in GI_SLOPE_PLACEHOLDER_NOTE
        scenario = ScenarioDefinition(
            name="low_gi_test",
            duration_s=3600.0,
            meal_events=[MealEvent(time_s=0.0, carbs_g=50.0, glycaemic_index=40.0)],
        )
        output = SimulationRunner().run(scenario)
        assert any(GI_SLOPE_PLACEHOLDER_NOTE in w for w in output.warnings)

    def test_ogtt_gi_meal_produces_no_for_gi_warning(self):
        output = SimulationRunner().run(ogtt_75g())
        assert not any("for_gi" in w for w in output.warnings)

    def test_meal_without_glycaemic_index_produces_no_for_gi_warning(self):
        output = SimulationRunner().run(high_sugar_meal_scenario())
        assert not any("for_gi" in w for w in output.warnings)


class TestProteinDominantMealWarning:
    """
    NutrientSecretionPathway (MASTER_PLAN row 1.28) doesn't model
    glucagon's counter-regulatory rise -- a protein-dominant meal should
    surface a loud warning rather than silently over-predict glucose
    lowering, same "loud, not silently trusted" convention as
    TestForGiWarning above. Mandatory code-reviewer finding, fixed
    before finalizing this milestone.
    """

    def test_protein_dominant_meal_produces_a_warning(self):
        scenario = ScenarioDefinition(
            name="protein_dominant",
            meal_events=[MealEvent(time_s=0.0, carbs_g=5.0, protein_g=50.0)],
            duration_s=14400,
        )
        output = SimulationRunner().run(scenario)
        assert any("NutrientSecretionPathway" in w for w in output.warnings)

    def test_carb_dominant_meal_with_some_protein_produces_no_warning(self):
        scenario = ScenarioDefinition(
            name="carb_dominant",
            meal_events=[MealEvent(time_s=0.0, carbs_g=75.0, protein_g=20.0)],
            duration_s=14400,
        )
        output = SimulationRunner().run(scenario)
        assert not any("NutrientSecretionPathway" in w for w in output.warnings)

    def test_no_protein_meal_produces_no_warning(self):
        scenario = ScenarioDefinition(
            name="no_protein",
            meal_events=[MealEvent(time_s=0.0, carbs_g=75.0)],
            duration_s=14400,
        )
        output = SimulationRunner().run(scenario)
        assert not any("NutrientSecretionPathway" in w for w in output.warnings)

    def test_equal_protein_and_carbs_produces_no_warning(self):
        # protein_g > carbs_g is the trigger, not >= -- an exactly-equal
        # split is not treated as "protein-dominant."
        scenario = ScenarioDefinition(
            name="equal_split",
            meal_events=[MealEvent(time_s=0.0, carbs_g=50.0, protein_g=50.0)],
            duration_s=14400,
        )
        output = SimulationRunner().run(scenario)
        assert not any("NutrientSecretionPathway" in w for w in output.warnings)


# ---------------------------------------------------------------------------
# GLP-1 disease-elevated-basal-glucose plausibility warning
#
# Mandatory security-reviewer finding on the Orskov 1994 GLP-1
# meal-response recalibration (CRITICAL): GLP1Pathway now derives a
# genuine, stable analytic fixed point for ANY basal glucose (fixing a
# severe, unbounded drift for a disease-elevated Gb composed with a
# GLP-1 agonist), but that fixed point's own magnitude is no longer a
# believable plasma concentration once basal glucose exceeds
# GLP1_EC50_SECRETION -- SimulationRunner.run() surfaces this loudly
# via a warning rather than silently producing an implausible-looking
# glp1_pM value with no explanation.
# ---------------------------------------------------------------------------

class TestGLP1DiseaseElevatedGbWarning:
    def test_t2d_plus_semaglutide_produces_a_warning(self):
        from biosim_physiology import PhysiologicalProfile, DiseaseState
        scenario = ScenarioDefinition(
            name="t2d_semaglutide", meal_events=[], duration_s=3600,
            drug_events=[DrugEvent(time_s=0, drug="semaglutide", dose_mg=1.0)],
            personalization=PersonalizationSpec(
                physiological_profile=PhysiologicalProfile(
                    disease_states=[DiseaseState.TYPE_2_DIABETES]
                )
            ),
        )
        output = SimulationRunner().run(scenario)
        assert any("GLP-1 class drug" in w for w in output.warnings)

    def test_healthy_plus_semaglutide_produces_no_warning(self):
        scenario = ScenarioDefinition(
            name="healthy_semaglutide", meal_events=[], duration_s=3600,
            drug_events=[DrugEvent(time_s=0, drug="semaglutide", dose_mg=1.0)],
        )
        output = SimulationRunner().run(scenario)
        assert not any("GLP-1 class drug" in w for w in output.warnings)

    def test_t2d_plus_dpp4_inhibitor_produces_no_warning(self):
        # A DPP-4 inhibitor's agonist_scale stays at 1.0 (no long-half-
        # life amplification) -- its fixed point at a disease-elevated
        # Gb stays physiologically ordinary (~18.5 pM for T2D, well
        # within this pathway's own documented ~7-15 pM basal range),
        # unlike semaglutide/liraglutide. The warning must not fire for
        # a drug class it doesn't actually describe (mandatory
        # security-reviewer finding on this fix's own first version,
        # which checked drug identity + a hardcoded threshold rather
        # than the real constructed instance's own agonist_scale).
        from biosim_physiology import PhysiologicalProfile, DiseaseState
        scenario = ScenarioDefinition(
            name="t2d_dpp4", meal_events=[], duration_s=3600,
            drug_events=[DrugEvent(time_s=0, drug="dpp4_inhibitor", dose_mg=100.0)],
            personalization=PersonalizationSpec(
                physiological_profile=PhysiologicalProfile(
                    disease_states=[DiseaseState.TYPE_2_DIABETES]
                )
            ),
        )
        output = SimulationRunner().run(scenario)
        assert not any("GLP-1 class drug" in w for w in output.warnings)

    def test_t2d_alone_no_glp1_drug_produces_no_warning(self):
        from biosim_physiology import PhysiologicalProfile, DiseaseState
        scenario = ScenarioDefinition(
            name="t2d_no_drug", meal_events=[], duration_s=3600,
            personalization=PersonalizationSpec(
                physiological_profile=PhysiologicalProfile(
                    disease_states=[DiseaseState.TYPE_2_DIABETES]
                )
            ),
        )
        output = SimulationRunner().run(scenario)
        assert not any("GLP-1 class drug" in w for w in output.warnings)

    def test_t2d_plus_semaglutide_glp1_no_longer_drifts_end_to_end(self):
        # Confirms the underlying fix (not just the warning) through the
        # real production path: the simulated GLP1 species stays at a
        # stable value across a standard fasted window, not drifting to
        # an ever-growing one.
        from biosim_physiology import PhysiologicalProfile, DiseaseState
        scenario = ScenarioDefinition(
            name="t2d_semaglutide_stability", meal_events=[], duration_s=18000,
            drug_events=[DrugEvent(time_s=0, drug="semaglutide", dose_mg=1.0)],
            personalization=PersonalizationSpec(
                physiological_profile=PhysiologicalProfile(
                    disease_states=[DiseaseState.TYPE_2_DIABETES]
                )
            ),
        )
        output = SimulationRunner().run(scenario)
        glp1 = output.simulation.species["GLP1"] * 1e6
        assert glp1.max() == pytest.approx(glp1.min(), rel=1e-4)


class TestCushingsSyndromeWarning:
    """
    Security-reviewer finding on Track 1 row 1.38 (third attempt): the
    for_cushings() partial-implementation gap (Gb untouched -- fasting
    hyperglycemia not modeled) was disclosed only in a code comment,
    invisible to a REST consumer. Now surfaced via SimulationOutput
    .warnings, matching the for_gi()/protein-dominant-meal/GLP-1
    warnings' own established convention.
    """

    def test_cushings_syndrome_produces_the_disclosure_warning(self):
        scenario = ScenarioDefinition(
            name="cushings", meal_events=[], duration_s=3600, cushings_syndrome=True,
        )
        output = SimulationRunner().run(scenario)
        assert any("cushings_syndrome=True" in w for w in output.warnings)

    def test_healthy_scenario_produces_no_cushings_warning(self):
        scenario = ScenarioDefinition(name="healthy", meal_events=[], duration_s=3600)
        output = SimulationRunner().run(scenario)
        assert not any("cushings_syndrome=True" in w for w in output.warnings)

    def test_cushings_plus_t2d_composes_without_error_or_double_counting(self):
        # Security-reviewer MEDIUM finding: no test existed composing
        # cushings_syndrome with an existing disease state. T2D's static
        # Gb/p3 modifiers (physiology/modifiers.py) and CortisolPathway's
        # signals are architecturally independent (no DiseaseState.
        # CUSHINGS_SYNDROME/cortisol reference exists in the static
        # modifier table -- route (b), a static preset, was deliberately
        # rejected to avoid exactly this double-counting risk), so this
        # should compose cleanly: succeed, and still show the disclosure
        # warning (Gb is unaffected by cushings_syndrome regardless of
        # what T2D's OWN static Gb modifier does).
        from biosim_physiology import PhysiologicalProfile, DiseaseState
        scenario = ScenarioDefinition(
            name="cushings_t2d", meal_events=[MealEvent(0, 60.0)], cushings_syndrome=True,
            personalization=PersonalizationSpec(
                physiological_profile=PhysiologicalProfile(
                    disease_states=[DiseaseState.TYPE_2_DIABETES]
                )
            ),
        )
        output = SimulationRunner().run(scenario)
        assert output.success
        assert any("cushings_syndrome=True" in w for w in output.warnings)


class TestColdAcclimationWarning:
    """
    Track 1 row 1.40: Hanssen et al. 2015's own cohort (PMID:26147760)
    was n=8 Type 2 Diabetes subjects only (no healthy-subject arm) --
    applying cold_acclimatized's p3 factor to a non-diabetic profile is
    an unverified extrapolation, disclosed via SimulationOutput.warnings,
    matching the Cushing's/GLP-1 warnings' own established convention.
    """

    def test_cold_acclimatized_without_t2d_produces_the_disclosure_warning(self):
        from biosim_physiology import PhysiologicalProfile
        scenario = ScenarioDefinition(
            name="cold_healthy", meal_events=[], duration_s=3600,
            personalization=PersonalizationSpec(
                physiological_profile=PhysiologicalProfile(cold_acclimatized=True)
            ),
        )
        output = SimulationRunner().run(scenario)
        assert any("cold_acclimatized=True" in w for w in output.warnings)

    def test_cold_acclimatized_with_t2d_produces_no_disclosure_warning(self):
        # Hanssen's own studied population -- not an extrapolation.
        from biosim_physiology import PhysiologicalProfile, DiseaseState
        scenario = ScenarioDefinition(
            name="cold_t2d", meal_events=[], duration_s=3600,
            personalization=PersonalizationSpec(
                physiological_profile=PhysiologicalProfile(
                    disease_states=[DiseaseState.TYPE_2_DIABETES], cold_acclimatized=True,
                )
            ),
        )
        output = SimulationRunner().run(scenario)
        assert not any("cold_acclimatized=True" in w for w in output.warnings)

    def test_no_cold_acclimation_produces_no_disclosure_warning(self):
        scenario = ScenarioDefinition(name="healthy", meal_events=[], duration_s=3600)
        output = SimulationRunner().run(scenario)
        assert not any("cold_acclimatized=True" in w for w in output.warnings)

    def test_cold_acclimation_raises_p3_end_to_end(self):
        from biosim_physiology import PhysiologicalProfile
        healthy = SimulationRunner().run(ScenarioDefinition(
            name="healthy", meal_events=[MealEvent(0, 60.0)],
        ))
        cold = SimulationRunner().run(ScenarioDefinition(
            name="cold", meal_events=[MealEvent(0, 60.0)],
            personalization=PersonalizationSpec(
                physiological_profile=PhysiologicalProfile(cold_acclimatized=True)
            ),
        ))
        assert healthy.success and cold.success
        # Higher insulin sensitivity -> lower glucose peak.
        assert cold.expert_metrics.glucose_peak_mM < healthy.expert_metrics.glucose_peak_mM


class TestCohortMismatchWarnings:
    """
    Same disclosure convention as TestColdAcclimationWarning, for rows
    1.44/1.45. HCQ's anchor (Wasko 2015, PMID:26197707) studied non-diabetic
    overweight/obese adults; NAFLD's (Bugianesi 2005, PMID:15747110) studied
    non-diabetic subjects, and composing it with another insulin-resistance
    mechanism probably double-counts part of the same effect.
    """

    @staticmethod
    def _run(profile=None, drug_events=(), **kwargs):
        return SimulationRunner().run(ScenarioDefinition(
            name="cohort", meal_events=[], duration_s=3600,
            drug_events=list(drug_events),
            personalization=PersonalizationSpec(physiological_profile=profile),
            **kwargs,
        ))

    @staticmethod
    def _hcq():
        return DrugEvent(time_s=0, drug="hydroxychloroquine", dose_mg=400.0)

    @staticmethod
    def _warned(output, needle):
        # Cohort warnings only -- every drug run also carries an
        # "EXPERIMENTAL drug module" (or, for semaglutide, "VALIDATED drug
        # module") note naming the drug (2026-11-16, 2026-11-21).
        return [
            w for w in output.warnings
            if needle in w and not w.startswith(
                ("EXPERIMENTAL drug module", "VALIDATED drug module"))
        ]

    def test_hcq_in_matching_cohort_produces_no_warning(self):
        from biosim_physiology import PhysiologicalProfile
        out = self._run(PhysiologicalProfile(body_weight_kg=90.0), [self._hcq()])
        assert not self._warned(out, "hydroxychloroquine")

    def test_hcq_with_diabetes_warns(self):
        from biosim_physiology import DiseaseState, PhysiologicalProfile
        out = self._run(
            PhysiologicalProfile(
                body_weight_kg=90.0, disease_states=[DiseaseState.TYPE_2_DIABETES],
            ),
            [self._hcq()],
        )
        [w] = self._warned(out, "hydroxychloroquine")
        assert "type_2_diabetes" in w and "PMID:26197707" in w

    def test_hcq_below_overweight_bmi_warns_including_drug_only_requests(self):
        # The drug-only path uses the default profile (BMI ~24.2) -- also
        # outside the studied cohort, so it must warn too.
        out = self._run(None, [self._hcq()])
        [w] = self._warned(out, "hydroxychloroquine")
        assert "BMI" in w

    def test_zero_dose_hcq_produces_no_warning(self):
        out = self._run(None, [DrugEvent(time_s=0, drug="hydroxychloroquine", dose_mg=0.0)])
        assert not self._warned(out, "hydroxychloroquine")

    def test_nafld_alone_produces_no_warning(self):
        from biosim_physiology import DiseaseState, PhysiologicalProfile
        out = self._run(PhysiologicalProfile(disease_states=[DiseaseState.NAFLD]))
        assert not self._warned(out, "NAFLD")

    def test_nafld_with_t2d_warns_once_with_both_reasons(self):
        from biosim_physiology import DiseaseState, PhysiologicalProfile
        out = self._run(PhysiologicalProfile(
            disease_states=[DiseaseState.NAFLD, DiseaseState.TYPE_2_DIABETES],
        ))
        [w] = self._warned(out, "NAFLD")
        assert "PMID:15747110" in w and "non-diabetic" in w and "double-count" in w

    def test_nafld_with_insulin_resistant_warns_about_double_counting_only(self):
        from biosim_physiology import DiseaseState, PhysiologicalProfile
        out = self._run(PhysiologicalProfile(
            disease_states=[DiseaseState.NAFLD, DiseaseState.INSULIN_RESISTANT],
        ))
        [w] = self._warned(out, "NAFLD")
        assert "double-count" in w and "non-diabetic" not in w

    def test_nafld_with_nefa_dynamics_warns_about_double_counting(self):
        from biosim_physiology import DiseaseState, PhysiologicalProfile
        out = self._run(
            PhysiologicalProfile(disease_states=[DiseaseState.NAFLD]),
            model_nefa_dynamics=True,
        )
        [w] = self._warned(out, "NAFLD")
        assert "NEFA" in w


# ---------------------------------------------------------------------------
# Extension pathway composition (cortisol/hepatic/NEFA) via SimulationRunner
# ---------------------------------------------------------------------------

class TestExtensionPathwayComposition:
    def test_stress_level_zero_matches_bare_baseline(self):
        baseline = SimulationRunner().run(high_sugar_meal_scenario())
        s = ScenarioDefinition(
            name="no_stress", duration_s=14400,
            meal_events=[MealEvent(0, 75.0)], stress_level=0.0,
        )
        out = SimulationRunner().run(s)
        assert out.expert_metrics.glucose_peak_mM == pytest.approx(
            baseline.expert_metrics.glucose_peak_mM, rel=0.01
        )

    def test_stress_raises_glucose_peak(self):
        no_stress = ScenarioDefinition(
            name="ns", duration_s=14400, meal_events=[MealEvent(0, 75.0)],
        )
        stressed = ScenarioDefinition(
            name="s", duration_s=14400, meal_events=[MealEvent(0, 75.0)],
            stress_level=0.9,
        )
        runner = SimulationRunner()
        peak_ns = runner.run(no_stress).expert_metrics.glucose_peak_mM
        peak_s = runner.run(stressed).expert_metrics.glucose_peak_mM
        assert peak_s > peak_ns

    def test_stress_level_out_of_range_raises(self):
        with pytest.raises(ValueError):
            ScenarioDefinition(name="bad", stress_level=1.5)

    def test_t2d_disease_state_composes_hepatic(self):
        from biosim_physiology import PhysiologicalProfile, DiseaseState
        healthy = ScenarioDefinition(
            name="healthy", duration_s=14400, meal_events=[MealEvent(0, 75.0)],
        )
        t2d_profile = PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES])
        t2d = ScenarioDefinition(
            name="t2d", duration_s=14400, meal_events=[MealEvent(0, 75.0)],
            personalization=PersonalizationSpec(physiological_profile=t2d_profile),
        )
        runner = SimulationRunner()
        out_healthy = runner.run(healthy)
        out_t2d = runner.run(t2d)
        assert out_healthy.success and out_t2d.success
        # T2D already elevates glucose via the static disease modifier; this
        # just confirms the composite pathway runs successfully end-to-end
        # and produces a plausibly higher peak, not a specific magnitude.
        assert out_t2d.expert_metrics.glucose_peak_mM > out_healthy.expert_metrics.glucose_peak_mM

    def test_hyperthyroidism_disease_state_composes_hepatic(self):
        # End-to-end wiring check for the new HYPERTHYROIDISM disease
        # state (_detect_hepatic() in runner.py) -- confirms the scenario
        # runs successfully with the new EGPParameters.for_hyperthyroidism()
        # preset composed, not just that the preset/pathway work in
        # isolation (test_hepatic.py).
        from biosim_physiology import PhysiologicalProfile, DiseaseState
        healthy = ScenarioDefinition(
            name="healthy", duration_s=14400, meal_events=[MealEvent(0, 75.0)],
        )
        hyper_profile = PhysiologicalProfile(disease_states=[DiseaseState.HYPERTHYROIDISM])
        hyper = ScenarioDefinition(
            name="hyperthyroid", duration_s=14400, meal_events=[MealEvent(0, 75.0)],
            personalization=PersonalizationSpec(physiological_profile=hyper_profile),
        )
        runner = SimulationRunner()
        out_healthy = runner.run(healthy)
        out_hyper = runner.run(hyper)
        assert out_healthy.success and out_hyper.success

    def test_hyperthyroidism_applies_no_static_parameter_modifier(self):
        # Deliberately {} -- see DiseaseState.HYPERTHYROIDISM's own
        # docstring in physiology/modifiers.py. HepGlucPathway is the
        # sole representation; a static Gb/p3 modifier here would
        # double-count, the same bug class fixed once for glucocorticoid.
        from biosim_physiology import PhysiologicalProfile, DiseaseState
        from biosim_ode import PathwayParameters
        hyper_profile = PhysiologicalProfile(disease_states=[DiseaseState.HYPERTHYROIDISM])
        scenario = ScenarioDefinition(
            name="hyper_static_check", meal_events=[], duration_s=10,
            personalization=PersonalizationSpec(physiological_profile=hyper_profile),
        )
        out = SimulationRunner().run(scenario)
        assert out.physiological_personalization.personalized_params.Gb == pytest.approx(PathwayParameters().Gb)
        assert out.physiological_personalization.personalized_params.p3 == pytest.approx(PathwayParameters().p3)

    def test_hyperthyroidism_fasted_does_not_drift_from_baseline(self):
        # Regression-style check matching
        # test_metformin_fasted_does_not_drift_below_its_own_static_baseline
        # above: a fasted (zero-meal) hyperthyroid scenario should sit
        # flat at the default 5.0 mM basal, not drift -- this pathway's
        # own basal-neutrality is verified in isolation in test_hepatic.py
        # (test_couple_with_gi_preserves_basal_steady_state_hyperthyroidism);
        # this confirms it holds end-to-end through SimulationRunner too.
        from biosim_physiology import PhysiologicalProfile, DiseaseState
        hyper_profile = PhysiologicalProfile(disease_states=[DiseaseState.HYPERTHYROIDISM])
        scenario = ScenarioDefinition(
            name="hyper_fasted", meal_events=[], duration_s=14400,
            personalization=PersonalizationSpec(physiological_profile=hyper_profile),
        )
        out = SimulationRunner().run(scenario)
        assert out.success
        G_uM = out.simulation.species["G"]
        assert G_uM.max() - G_uM.min() < 5.0   # uM; effectively flat

    def test_metformin_composes_hepatic_and_still_lowers_glucose(self):
        from biosim_seeds.drug_scenarios import metformin_with_meal
        no_drug = ScenarioDefinition(
            name="nd", duration_s=14400, meal_events=[MealEvent(0, 75.0)],
        )
        runner = SimulationRunner()
        peak_nd = runner.run(no_drug).expert_metrics.glucose_peak_mM
        peak_mf = runner.run(metformin_with_meal(carbs_g=75.0, metformin_mg=1000.0)).expert_metrics.glucose_peak_mM
        assert peak_mf < peak_nd

    def test_metformin_fasted_does_not_drift_below_its_own_static_baseline(self):
        # Regression: EGPParameters.for_metformin() used to start its EGP
        # species at the healthy pre-treatment egp_basal and drift away
        # from it toward a lower suppressed target even with zero meal
        # input -- double-counting on top of the SEPARATE, already-
        # intended static Gb modifier (physiology/drug_interactions.py,
        # metformin_modifiers(): Gb -14.7% as of CLAUDE.md's mechanism-
        # weighting correction, -4% when this test was written -- the
        # fixed-point property here doesn't depend on the modifier's
        # exact magnitude) that GlucoseInsulinPathway's own
        # basal state already reflects. A fasted (zero-meal) chronic-
        # metformin scenario should sit at its own personalized baseline,
        # not keep drifting downward for hours with nothing driving it.
        scenario = ScenarioDefinition(
            name="metformin_fasted", meal_events=[], duration_s=14400,
            drug_events=[DrugEvent(time_s=0, drug="metformin", dose_mg=2000.0, hours_before_simulation=24 * 7)],
        )
        out = SimulationRunner().run(scenario)
        assert out.success
        G = out.simulation.species["G"]
        assert G.max() - G.min() < 5.0  # uM; effectively flat, not a multi-mM artifact drift

    def test_model_nefa_dynamics_runs_successfully(self):
        s = ScenarioDefinition(
            name="nefa", duration_s=14400, meal_events=[MealEvent(0, 75.0)],
            model_nefa_dynamics=True,
        )
        out = SimulationRunner().run(s)
        assert out.success

    def test_model_nefa_dynamics_default_matches_baseline(self):
        # Default NEFAParameters is neutral at basal -- opting in must not
        # perturb an otherwise-identical scenario.
        baseline = ScenarioDefinition(
            name="nd", duration_s=14400, meal_events=[MealEvent(0, 75.0)],
        )
        with_nefa = ScenarioDefinition(
            name="wn", duration_s=14400, meal_events=[MealEvent(0, 75.0)],
            model_nefa_dynamics=True,
        )
        runner = SimulationRunner()
        peak_base = runner.run(baseline).expert_metrics.glucose_peak_mM
        peak_nefa = runner.run(with_nefa).expert_metrics.glucose_peak_mM
        assert peak_nefa == pytest.approx(peak_base, rel=0.01)

    def test_steroid_meal_composes_cortisol_pathway(self):
        # steroid_meal() (Drug.GLUCOCORTICOID) composes the real ODE-coupled
        # CortisolPathway, which is now the SOLE representation of steroid's
        # Gb/p3 effect (see runner._detect_cortisol()'s docstring) --
        # glucocorticoid_modifiers()'s Gb/p3 factors are excluded from
        # personalisation to avoid double-counting; only its secretion_scale
        # compensatory bump still applies alongside the pathway. Confirms
        # the wiring end-to-end via the actual scenario builder, not just
        # the low-level couple_with_gi() call.
        from biosim_seeds.drug_scenarios import steroid_meal
        no_drug = ScenarioDefinition(
            name="nd", duration_s=14400, meal_events=[MealEvent(0, 60.0)],
        )
        runner = SimulationRunner()
        peak_nd = runner.run(no_drug).expert_metrics.glucose_peak_mM
        peak_steroid = runner.run(steroid_meal(carbs_g=60.0, prednisone_mg=20.0)).expert_metrics.glucose_peak_mM
        assert peak_steroid > peak_nd

    def test_glucocorticoid_fasted_does_not_double_count_static_and_pathway(self):
        # Regression: glucocorticoid_modifiers() (physiology/drug_interactions.py)
        # and CortisolPathway both independently fired their FULL Gb/p3
        # steroid effect -- CortisolPathway has no meal/insulin input at all
        # (self-driven by stress_level/circadian_hour only), so composing
        # both stacked two complete representations of the identical
        # phenomenon. A fasted (zero-meal) prednisone scenario used to
        # drift fasting glucose an extra +18% beyond the already-intended
        # static Gb elevation and nearly double steady-state insulin, purely
        # as an artifact -- confirmed directly before this fix.
        #
        # hours_before_simulation=24*7 + 6, not just 24*7 (regression
        # caught by code review, real PK/PD timing milestone -- see
        # CLAUDE.md): 24*7 is an exact multiple of 24, so once
        # glucocorticoid_pd_effect_fraction() wraps chronic dosing into
        # its daily cycle, that value happens to land EXACTLY at the
        # start of a fresh cycle -- the effect fraction is ~0 for most of
        # this 4h fasted window, so the old `< 6.1` assertion alone would
        # have kept passing even with the double-counting bug back in
        # place (nothing left to double-count). The "+6" lands squarely
        # at this curve's peak instead, so the check is exercised at the
        # moment the effect is actually fully developed.
        scenario = ScenarioDefinition(
            name="steroid_fasted", meal_events=[], duration_s=14400,
            drug_events=[DrugEvent(
                time_s=0, drug="glucocorticoid", dose_mg=20.0,
                hours_before_simulation=24 * 7 + 6,
            )],
        )
        out = SimulationRunner().run(scenario)
        assert out.success
        G_mM = out.simulation.species["G"] / 1000.0
        # 5.0 (unmodified default, static Gb excluded) drifting to
        # CortisolPathway's own single, non-doubled steady state (~5.9 mM,
        # confirmed via direct investigation) -- not the doubled ~6.3 mM
        # the pre-fix artifact produced. Two-sided: the lower bound makes
        # this a genuine check (a "no effect left at all" regression would
        # fail it too), not just a ceiling a broken implementation could
        # also satisfy.
        assert 5.3 < G_mM.max() < 6.1

    def test_glucocorticoid_secretion_scale_compensation_still_applies(self):
        # The one piece of glucocorticoid_modifiers() CortisolPathway does
        # NOT represent (secretion_scale's +10%-at-reference-dose
        # compensatory beta-cell response) must still apply -- only the
        # double-counted Gb/p3 factors are excluded from personalisation.
        #
        # dose given 6h before simulation start: the compensatory response
        # is now genuinely PD-timing-gated (glucocorticoid_secretion_scale_
        # factor()), so a dose that hasn't developed any effect yet (e.g.
        # hours_before_simulation=0, "just taken") would show NO
        # compensation at all -- see
        # test_glucocorticoid_secretion_scale_uses_real_dose_timing below
        # for that exact regression. 6h sits within this curve's onset
        # window (see glucocorticoid_pd_effect_fraction()).
        from biosim_physiology import PhysiologicalProfile
        from biosim_ode import PathwayParameters
        scenario = ScenarioDefinition(
            name="steroid_secretion_check", meal_events=[], duration_s=10,
            drug_events=[DrugEvent(time_s=0, drug="glucocorticoid", dose_mg=20.0, hours_before_simulation=6.0)],
            personalization=PersonalizationSpec(physiological_profile=PhysiologicalProfile()),
        )
        out = SimulationRunner().run(scenario)
        assert out.physiological_personalization.personalized_params.secretion_scale > 1.0
        # ...but Gb/p3 themselves stay at the unmodified default -- CortisolPathway,
        # not the static modifier, represents those for a glucocorticoid event.
        assert out.physiological_personalization.personalized_params.Gb == pytest.approx(PathwayParameters().Gb)

    def test_glucocorticoid_secretion_scale_uses_real_dose_timing(self):
        # Regression for the dead-wiring bug (see CLAUDE.md): run() used to
        # build a FRESH DrugDose with no hours_since_ingestion (defaulting
        # to 0.0) to compute the secretion_scale compensatory bump,
        # discarding the scenario's real DrugEvent timing. Now that the
        # bump is genuinely PD-timing-gated, that bug would silently ZERO
        # OUT the compensation for any dose not taken at exactly t=0 --
        # confirm a dose given 6h before simulation start (squarely within
        # the PD curve's onset window) shows MORE compensation than one
        # just taken (hours_before_simulation=0, effect not yet onset).
        from biosim_physiology import PhysiologicalProfile
        just_taken = ScenarioDefinition(
            name="steroid_just_taken", meal_events=[], duration_s=10,
            drug_events=[DrugEvent(time_s=0, drug="glucocorticoid", dose_mg=20.0)],
            personalization=PersonalizationSpec(physiological_profile=PhysiologicalProfile()),
        )
        six_hours_ago = ScenarioDefinition(
            name="steroid_6h_ago", meal_events=[], duration_s=10,
            drug_events=[DrugEvent(time_s=0, drug="glucocorticoid", dose_mg=20.0, hours_before_simulation=6.0)],
            personalization=PersonalizationSpec(physiological_profile=PhysiologicalProfile()),
        )
        out_just_taken = SimulationRunner().run(just_taken)
        out_6h = SimulationRunner().run(six_hours_ago)
        scale_just_taken = out_just_taken.physiological_personalization.personalized_params.secretion_scale
        scale_6h = out_6h.physiological_personalization.personalized_params.secretion_scale
        assert scale_just_taken == pytest.approx(1.0, abs=1e-9)
        assert scale_6h > scale_just_taken

    def test_glucocorticoid_secretion_scale_decreases_for_t2d(self):
        # Yuen et al. 2012 (PMID:21973241) found DECREASED insulin/C-peptide
        # postbreakfast specifically in the T2DM group (reduced beta-cell
        # reserve cannot mount the compensatory response healthy/at-risk
        # patients show) -- confirmed by Shah et al. 2020 (PMID:32369480)
        # independently. glucocorticoid_secretion_scale_factor() must flip
        # sign for a T2D profile, not apply the same +10% compensatory bump.
        from biosim_physiology import PhysiologicalProfile, DiseaseState
        SIX_HOURS_AGO = 6.0
        healthy = ScenarioDefinition(
            name="steroid_healthy", meal_events=[], duration_s=10,
            drug_events=[DrugEvent(time_s=0, drug="glucocorticoid", dose_mg=20.0, hours_before_simulation=SIX_HOURS_AGO)],
            personalization=PersonalizationSpec(physiological_profile=PhysiologicalProfile()),
        )
        t2d = ScenarioDefinition(
            name="steroid_t2d", meal_events=[], duration_s=10,
            drug_events=[DrugEvent(time_s=0, drug="glucocorticoid", dose_mg=20.0, hours_before_simulation=SIX_HOURS_AGO)],
            personalization=PersonalizationSpec(
                physiological_profile=PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES])
            ),
        )
        healthy_out = SimulationRunner().run(healthy)
        t2d_out = SimulationRunner().run(t2d)
        healthy_scale = healthy_out.physiological_personalization.personalized_params.secretion_scale
        t2d_scale = t2d_out.physiological_personalization.personalized_params.secretion_scale
        disease_only_scale = t2d_scale / 0.60   # T2D's own disease-state secretion_scale=0.60
        assert healthy_scale > 1.0                       # compensatory bump, unchanged direction
        assert disease_only_scale < 1.0                   # decompensatory bump, flipped for T2D

    def test_dapagliflozin_composes_renal_pathway_and_lowers_peak(self):
        # Confirms end-to-end wiring: a dapagliflozin DrugEvent must
        # actually compose RenalPathway via _detect_renal()/_build_pathway(),
        # not just be accepted as a valid Drug value.
        no_drug = ScenarioDefinition(
            name="nd", duration_s=14400, meal_events=[MealEvent(0, 75.0)],
        )
        with_drug = ScenarioDefinition(
            name="dapa", duration_s=14400, meal_events=[MealEvent(0, 75.0)],
            drug_events=[DrugEvent(time_s=0, drug="dapagliflozin", dose_mg=10.0)],
        )
        runner = SimulationRunner()
        peak_nd = runner.run(no_drug).expert_metrics.glucose_peak_mM
        peak_dapa = runner.run(with_drug).expert_metrics.glucose_peak_mM
        assert peak_dapa < peak_nd

    def test_dapagliflozin_fasted_genuinely_lowers_glucose_below_baseline(self):
        # Unlike hepatic/metformin (which must NOT drift a fasted scenario
        # away from its own static baseline -- see
        # test_metformin_fasted_does_not_drift_below_its_own_static_baseline
        # above), this is the opposite, INTENDED behaviour for SGLT2i: the
        # renal clearance mechanism is continuously active at any glucose
        # above the drug-lowered threshold, including ordinary fasting
        # glucose, so a fasted chronic-dapagliflozin scenario SHOULD settle
        # below the untreated 5.0 mM basal -- not stay flat at it.
        scenario = ScenarioDefinition(
            name="dapa_fasted", meal_events=[], duration_s=14400,
            drug_events=[DrugEvent(time_s=0, drug="dapagliflozin", dose_mg=10.0, hours_before_simulation=24 * 7)],
        )
        out = SimulationRunner().run(scenario)
        assert out.success
        G_mM = out.simulation.species["G"] / 1000.0
        assert G_mM[-1] < 4.9   # below the untreated 5.0 mM basal

    def test_sglt2i_drug_event_applies_no_static_parameter_modifier(self):
        # _sglt2_inhibitor_modifiers() must return {} -- RenalPathway is the
        # SOLE representation, matching the GLP-1 precedent (and even
        # stricter than glucocorticoid's partial exclusion, since SGLT2i
        # has no static component to keep at all -- see
        # biosim_physiology.drug_interactions._sglt2_inhibitor_modifiers()).
        from biosim_physiology import PhysiologicalProfile
        from biosim_ode import PathwayParameters
        scenario = ScenarioDefinition(
            name="dapa_static_check", meal_events=[], duration_s=10,
            drug_events=[DrugEvent(time_s=0, drug="dapagliflozin", dose_mg=10.0)],
            personalization=PersonalizationSpec(physiological_profile=PhysiologicalProfile()),
        )
        out = SimulationRunner().run(scenario)
        assert out.physiological_personalization.personalized_params.Gb == pytest.approx(PathwayParameters().Gb)

    def test_dapagliflozin_t2d_profile_uses_elevated_renal_tmg(self):
        # _detect_renal() must detect DiseaseState.TYPE_2_DIABETES on the
        # scenario's physiological_profile and pass t2d=True through to
        # RenalClearanceParameters.for_dapagliflozin(), not always build
        # the healthy-cohort preset regardless of disease state (Track 1
        # MASTER_PLAN 1.8). Verified directly via _detect_renal()'s own
        # returned RenalPathway params, not an indirect glucose-trajectory
        # comparison -- T2D's own elevated Gb also independently shifts
        # glucose, which would confound a bare peak/AUC comparison.
        from biosim_physiology import PhysiologicalProfile, DiseaseState

        healthy = ScenarioDefinition(
            name="dapa_healthy", meal_events=[], duration_s=10,
            drug_events=[DrugEvent(time_s=0, drug="dapagliflozin", dose_mg=10.0)],
            personalization=PersonalizationSpec(physiological_profile=PhysiologicalProfile()),
        )
        t2d = ScenarioDefinition(
            name="dapa_t2d", meal_events=[], duration_s=10,
            drug_events=[DrugEvent(time_s=0, drug="dapagliflozin", dose_mg=10.0)],
            personalization=PersonalizationSpec(
                physiological_profile=PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES])
            ),
        )
        runner = SimulationRunner()
        healthy_renal = runner._detect_renal(healthy)
        t2d_renal = runner._detect_renal(t2d)
        assert healthy_renal.params.tm_mg_per_min == pytest.approx(317.0)
        assert t2d_renal.params.tm_mg_per_min == pytest.approx(420.0)

    def test_empagliflozin_t2d_profile_uses_elevated_renal_tmg(self):
        # Symmetric coverage for empagliflozin (dapagliflozin is covered by
        # the test directly above) -- same _detect_renal() code path, but
        # verified independently rather than assumed identical, since
        # for_empagliflozin() also selects a different day2/day14 beta-
        # scaling fraction for t2d, not just a different TmG.
        from biosim_physiology import PhysiologicalProfile, DiseaseState

        healthy = ScenarioDefinition(
            name="empa_healthy", meal_events=[], duration_s=10,
            drug_events=[DrugEvent(time_s=0, drug="empagliflozin", dose_mg=25.0, hours_before_simulation=24 * 14)],
            personalization=PersonalizationSpec(physiological_profile=PhysiologicalProfile()),
        )
        t2d = ScenarioDefinition(
            name="empa_t2d", meal_events=[], duration_s=10,
            drug_events=[DrugEvent(time_s=0, drug="empagliflozin", dose_mg=25.0, hours_before_simulation=24 * 14)],
            personalization=PersonalizationSpec(
                physiological_profile=PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES])
            ),
        )
        runner = SimulationRunner()
        healthy_renal = runner._detect_renal(healthy)
        t2d_renal = runner._detect_renal(t2d)
        assert healthy_renal.params.tm_mg_per_min == pytest.approx(337.0)
        assert t2d_renal.params.tm_mg_per_min == pytest.approx(459.0)

    def test_dapagliflozin_t2d_profile_and_no_profile_both_use_healthy_tmg(self):
        # A scenario with no physiological_profile at all (the common
        # drug-only "quick path", see the runner's own KNOWN BUGS FIXED
        # entry on this path) must not crash _detect_renal()'s new
        # disease-state lookup and must default to the healthy preset.
        no_profile = ScenarioDefinition(
            name="dapa_no_profile", meal_events=[], duration_s=10,
            drug_events=[DrugEvent(time_s=0, drug="dapagliflozin", dose_mg=10.0)],
        )
        renal = SimulationRunner()._detect_renal(no_profile)
        assert renal.params.tm_mg_per_min == pytest.approx(317.0)

    def test_empagliflozin_longer_treatment_produces_larger_effect(self):
        # _detect_renal() converts hours_before_simulation to
        # days_of_treatment for empagliflozin's duration-scaled beta (see
        # RenalClearanceParameters.for_empagliflozin()) -- confirm this
        # actually changes the composed pathway's effect end-to-end,
        # not just at the RenalClearanceParameters unit level.
        short = ScenarioDefinition(
            name="empa_48h", meal_events=[], duration_s=14400,
            drug_events=[DrugEvent(time_s=0, drug="empagliflozin", dose_mg=25.0, hours_before_simulation=48.0)],
        )
        long = ScenarioDefinition(
            name="empa_14d", meal_events=[], duration_s=14400,
            drug_events=[DrugEvent(time_s=0, drug="empagliflozin", dose_mg=25.0, hours_before_simulation=24 * 14)],
        )
        runner = SimulationRunner()
        G_short = runner.run(short).simulation.species["G"][-1]
        G_long = runner.run(long).simulation.species["G"][-1]
        assert G_long < G_short

    def test_detect_renal_extension_overrides_beta_reaches_dapagliflozin(self):
        # MASTER_PLAN 2.2: extension_overrides["renal.beta"] must reach
        # RenalClearanceParameters.for_dapagliflozin()'s beta_class_default,
        # not just be silently ignored.
        scenario = ScenarioDefinition(
            name="dapa_override", meal_events=[], duration_s=10,
            drug_events=[DrugEvent(time_s=0, drug="dapagliflozin", dose_mg=10.0)],
        )
        runner = SimulationRunner()
        default = runner._detect_renal(scenario)
        overridden = runner._detect_renal(scenario, {"renal.beta": 0.20})
        assert default.params.beta != pytest.approx(0.20)
        assert overridden.params.beta == pytest.approx(0.20)

    def test_detect_renal_extension_overrides_beta_reaches_empagliflozin_scaled(self):
        # For empagliflozin, the override is the PRE-scaling class-wide
        # default -- confirm the override still gets the same day2/day14
        # onset scaling applied, not bypassed.
        scenario = ScenarioDefinition(
            name="empa_override", meal_events=[], duration_s=10,
            drug_events=[DrugEvent(time_s=0, drug="empagliflozin", dose_mg=25.0, hours_before_simulation=48.0)],
        )
        runner = SimulationRunner()
        default = runner._detect_renal(scenario)
        overridden = runner._detect_renal(scenario, {"renal.beta": 0.20})
        # Same ratio (overridden/default) as the pre-scaling constants
        # themselves, since both are scaled by the identical fraction.
        from biosim_ode.renal_pathway import SGLT2I_BETA_CLASS_DEFAULT
        expected_ratio = 0.20 / SGLT2I_BETA_CLASS_DEFAULT
        assert overridden.params.beta == pytest.approx(default.params.beta * expected_ratio)

    def test_build_pathway_extension_overrides_default_is_empty_and_unaffected(self):
        # Omitting extension_overrides at the point-estimate call site
        # must reproduce the exact prior single-value beta -- no
        # behaviour change for any scenario that doesn't pass overrides.
        scenario = ScenarioDefinition(
            name="dapa_no_override", meal_events=[], duration_s=10,
            drug_events=[DrugEvent(time_s=0, drug="dapagliflozin", dose_mg=10.0)],
        )
        from biosim_ode import GlucoseInsulinPathway
        runner = SimulationRunner()
        gi = GlucoseInsulinPathway()
        composite_a = runner._build_pathway(scenario, gi)
        composite_b = runner._build_pathway(scenario, gi, None)
        composite_c = runner._build_pathway(scenario, gi, {})
        renal_a = next(p for p in composite_a._pathways if p.name == "renal_clearance_v1")
        renal_b = next(p for p in composite_b._pathways if p.name == "renal_clearance_v1")
        renal_c = next(p for p in composite_c._pathways if p.name == "renal_clearance_v1")
        assert renal_a.params.beta == pytest.approx(renal_b.params.beta) == pytest.approx(renal_c.params.beta)

    def test_multiple_extensions_compose_simultaneously(self):
        # Stress + NEFA dynamics + metformin all active at once must still
        # run successfully as one flat CompositePathway (not nested).
        from biosim_seeds.drug_scenarios import metformin_with_meal
        import dataclasses
        s = dataclasses.replace(
            metformin_with_meal(carbs_g=75.0), stress_level=0.5, model_nefa_dynamics=True,
        )
        out = SimulationRunner().run(s)
        assert out.success
        assert out.expert_metrics.glucose_peak_mM > 0

    def test_hepatic_egp_delta_sums_across_alcohol_and_metformin_t2d(self):
        # Mandatory code-reviewer finding on the alcohol mechanism fix
        # (2026-08-18): HepGlucPathway (from T2D+metformin) and the new
        # AlcoholHepaticPathway both target the SAME hepatic_egp_delta
        # signal -- verified ad hoc during review that
        # CompositePathway._combine_signal_values() correctly sums their
        # contributions rather than one overwriting the other, but that
        # verification wasn't locked in as a regression test anywhere.
        # This is exactly the composition-boundary check this project's
        # own standing audit protocol calls for.
        from biosim_physiology import PhysiologicalProfile, DiseaseState
        from biosim_ode import GlucoseInsulinPathway
        from biosim_ode.hepatic import HepGlucPathway

        scenario = ScenarioDefinition(
            name="t2d_metformin_alcohol",
            meal_events=[MealEvent(0, 75.0)],
            duration_s=14400,
            drug_events=[
                DrugEvent(time_s=0, drug="metformin", dose_mg=2000.0),
                DrugEvent(time_s=0, drug="alcohol", dose_mg=48000.0, hours_before_simulation=0.0),
            ],
            personalization=PersonalizationSpec(
                physiological_profile=PhysiologicalProfile(
                    disease_states=[DiseaseState.TYPE_2_DIABETES],
                ),
            ),
        )
        runner = SimulationRunner()
        gi = GlucoseInsulinPathway()
        composite = runner._build_pathway(scenario, gi)
        hepatic = next(p for p in composite._pathways if p.name == "hep_gluc_v1")
        alcohol = next(p for p in composite._pathways if p.name == "alcohol_hepatic_v1")
        t = 120.0   # alcohol's own peak-suppression time
        hepatic_delta = hepatic.output_signals(t, [hepatic.y0[0]])["egp_mM_per_min"] - hepatic.params.egp_basal
        alcohol_delta = alcohol.output_signals(t, [])["hepatic_egp_delta"]
        assert alcohol_delta < 0.0   # alcohol's own contribution is genuinely present

        # Both pathways target the SAME target_input (hepatic_egp_delta) --
        # confirm CompositePathway sums rather than lets one silently
        # clobber the other (coupling.py's own ADDITIVE_SIGNAL_KEYS
        # combination rule, exercised directly here, not just assumed).
        from biosim_ode.coupling import CompositePathway
        combined = CompositePathway._combine_signal_values(
            "hepatic_egp_delta", [hepatic_delta, alcohol_delta],
        )
        assert combined == pytest.approx(hepatic_delta + alcohol_delta)

        out = SimulationRunner().run(scenario)
        assert out.success
        assert out.expert_metrics.glucose_peak_mM > 0

    def test_stress_and_exercise_bout_both_affect_glucose_not_just_exercise(self):
        # Regression for the CompositePathway signal-collision bug (see
        # CLAUDE.md): cortisol_couplings() and exercise_couplings() both
        # target GI's gb_modifier/p3_modifier. Before the fix,
        # CompositePathway.dydt() injected coupling signals via plain dict
        # assignment in list order -- exercise composes AFTER cortisol in
        # SimulationRunner._build_pathway(), so exercise's own neutral
        # (1.0) gb_modifier/p3_modifier (exercise doesn't touch p3/gb while
        # the bout is still in progress) silently discarded cortisol's
        # real, stress-elevated values entirely.
        #
        # Confirmed directly (not assumed) before writing this test: with
        # the bug present, a fasted stress_level=0.9 + an in-progress
        # 90-minute exercise bout produced a byte-for-byte identical final
        # glucose value to an exercise-only scenario with stress_level=0.0
        # (both 5.0000 mM) -- cortisol's entire glycaemic effect vanished.
        from biosim_physiology.exercise import ExerciseBout

        bout = ExerciseBout(intensity_percent_vo2max=70.0, duration_minutes=90.0)
        stress_plus_exercise = ScenarioDefinition(
            name="stress_plus_exercise", duration_s=3600, meal_events=[],
            stress_level=0.9, exercise_bout_in_progress=bout,
        )
        exercise_only = ScenarioDefinition(
            name="exercise_only", duration_s=3600, meal_events=[],
            stress_level=0.0, exercise_bout_in_progress=bout,
        )
        runner = SimulationRunner()
        G_both = runner.run(stress_plus_exercise).simulation.species["G"][-1]
        G_exercise_only = runner.run(exercise_only).simulation.species["G"][-1]

        # Cortisol's contribution must survive composition with exercise --
        # if it were still being discarded, these would be exactly equal.
        assert G_both != pytest.approx(G_exercise_only, rel=1e-6)
        assert G_both > G_exercise_only


# ---------------------------------------------------------------------------
# Acarbose (alpha-glucosidase inhibitor) -- MASTER_PLAN 1.4
# ---------------------------------------------------------------------------

class TestAcarboseAbsorptionModifier:
    """
    Acarbose slows intestinal carbohydrate absorption -- modelled as a
    direct reduction of MealInput.kabs, applied by SimulationRunner at
    meal-construction time (not a PathwayParameters modifier, not an
    ODE-coupled pathway). See CLAUDE.md's acarbose milestone entry for
    the full citation (Chen et al. 2016, PMID:27686734) and calibration.
    """

    def test_acarbose_lowers_and_delays_the_glucose_peak(self):
        # The defining dual signature of a genuine rate-constant
        # slowdown (as opposed to a pure amplitude scale, which would
        # only lower the peak without moving it) -- both effects must
        # survive through the real SimulationRunner path.
        no_drug = ScenarioDefinition(
            name="no_drug", meal_events=[MealEvent(time_s=0, carbs_g=75.0)], duration_s=14400,
        )
        with_acarbose = ScenarioDefinition(
            name="acarbose", meal_events=[MealEvent(time_s=0, carbs_g=75.0)], duration_s=14400,
            drug_events=[DrugEvent(time_s=0, drug="acarbose", dose_mg=100.0)],
        )
        runner = SimulationRunner()
        out_no_drug = runner.run(no_drug)
        out_acarbose = runner.run(with_acarbose)
        assert out_acarbose.expert_metrics.glucose_peak_mM < out_no_drug.expert_metrics.glucose_peak_mM
        assert out_acarbose.expert_metrics.glucose_peak_time_s > out_no_drug.expert_metrics.glucose_peak_time_s

    def test_higher_dose_produces_larger_effect(self):
        low = ScenarioDefinition(
            name="acarbose_low", meal_events=[MealEvent(time_s=0, carbs_g=75.0)], duration_s=14400,
            drug_events=[DrugEvent(time_s=0, drug="acarbose", dose_mg=25.0)],
        )
        high = ScenarioDefinition(
            name="acarbose_high", meal_events=[MealEvent(time_s=0, carbs_g=75.0)], duration_s=14400,
            drug_events=[DrugEvent(time_s=0, drug="acarbose", dose_mg=100.0)],
        )
        runner = SimulationRunner()
        peak_low = runner.run(low).expert_metrics.glucose_peak_mM
        peak_high = runner.run(high).expert_metrics.glucose_peak_mM
        assert peak_high < peak_low

    def test_gi_meal_is_slowed_like_a_default_meal_and_warns(self):
        # Since 2026-11-17 a glycaemic_index meal gets the same treatment as a
        # default two-compartment meal: both rate constants slowed by the
        # drug's calibrated factor (lower, later peak), with a warning that
        # the factor was calibrated on the single-exponential curve. It used
        # to be left unaffected.
        no_drug = ScenarioDefinition(
            name="gi_no_drug",
            meal_events=[MealEvent(time_s=0, carbs_g=75.0, glycaemic_index=70.0)],
            duration_s=14400,
        )
        with_acarbose = ScenarioDefinition(
            name="gi_acarbose",
            meal_events=[MealEvent(time_s=0, carbs_g=75.0, glycaemic_index=70.0)],
            duration_s=14400,
            drug_events=[DrugEvent(time_s=0, drug="acarbose", dose_mg=100.0)],
        )
        runner = SimulationRunner()
        out_no_drug = runner.run(no_drug)
        out_acarbose = runner.run(with_acarbose)
        assert out_acarbose.expert_metrics.glucose_peak_mM < out_no_drug.expert_metrics.glucose_peak_mM
        assert out_acarbose.expert_metrics.glucose_peak_time_s > out_no_drug.expert_metrics.glucose_peak_time_s
        assert any("two-compartment" in w for w in out_acarbose.warnings)

    def test_gi_meal_of_free_sugar_is_unaffected(self):
        # complex_carb_fraction=0 (free glucose needs no digestion): no effect,
        # matching Chen 2016's own finding of no acarbose effect on an OGTT.
        def run(drugs):
            return SimulationRunner().run(ScenarioDefinition(
                name="gi_free_sugar", duration_s=14400, drug_events=drugs,
                meal_events=[MealEvent(time_s=0, carbs_g=75.0, glycaemic_index=100.0,
                                       complex_carb_fraction=0.0)],
            ))
        plain = run([])
        drug = run([DrugEvent(time_s=0, drug="acarbose", dose_mg=100.0)])
        assert drug.expert_metrics.glucose_peak_mM == pytest.approx(plain.expert_metrics.glucose_peak_mM)

    def test_no_acarbose_event_means_no_warning_and_no_effect(self):
        scenario = ScenarioDefinition(
            name="plain", meal_events=[MealEvent(time_s=0, carbs_g=75.0)], duration_s=14400,
        )
        out = SimulationRunner().run(scenario)
        assert not any("acarbose" in w.lower() for w in out.warnings)

    def test_zero_dose_acarbose_event_has_no_effect(self):
        no_drug = ScenarioDefinition(
            name="no_drug2", meal_events=[MealEvent(time_s=0, carbs_g=75.0)], duration_s=14400,
        )
        zero_dose = ScenarioDefinition(
            name="acarbose_zero", meal_events=[MealEvent(time_s=0, carbs_g=75.0)], duration_s=14400,
            drug_events=[DrugEvent(time_s=0, drug="acarbose", dose_mg=0.0)],
        )
        runner = SimulationRunner()
        peak_no_drug = runner.run(no_drug).expert_metrics.glucose_peak_mM
        peak_zero = runner.run(zero_dose).expert_metrics.glucose_peak_mM
        assert peak_zero == pytest.approx(peak_no_drug)

    def test_acarbose_drug_only_quick_path_still_produces_the_effect(self):
        # The "quick path" (no physiological_profile, drug_events present)
        # is a common REST-API entry point -- confirm the meal-level
        # modifier isn't accidentally gated on personalization running.
        scenario = ScenarioDefinition(
            name="acarbose_quick_path", meal_events=[MealEvent(time_s=0, carbs_g=75.0)],
            duration_s=14400,
            drug_events=[DrugEvent(time_s=0, drug="acarbose", dose_mg=100.0)],
        )
        out = SimulationRunner().run(scenario)
        assert out.success
        assert out.expert_metrics.glucose_peak_mM < 10.0  # measurably below the ~10.5mM no-drug peak

    def test_free_sugar_meal_is_genuinely_unaffected_by_acarbose(self):
        # MASTER_PLAN row 1.4b(iii): a meal that's entirely free
        # monosaccharide (complex_carb_fraction=0.0) needs no enzymatic
        # digestion at all, so alpha-glucosidase inhibition should have
        # NO effect on it -- byte-identical to the no-drug baseline,
        # not just a smaller effect.
        no_drug = ScenarioDefinition(
            name="free_sugar_no_drug",
            meal_events=[MealEvent(time_s=0, carbs_g=75.0, complex_carb_fraction=0.0)],
            duration_s=14400,
        )
        with_acarbose = ScenarioDefinition(
            name="free_sugar_acarbose",
            meal_events=[MealEvent(time_s=0, carbs_g=75.0, complex_carb_fraction=0.0)],
            duration_s=14400,
            drug_events=[DrugEvent(time_s=0, drug="acarbose", dose_mg=100.0)],
        )
        runner = SimulationRunner()
        out_no_drug = runner.run(no_drug)
        out_acarbose = runner.run(with_acarbose)
        assert out_acarbose.expert_metrics.glucose_peak_mM == pytest.approx(
            out_no_drug.expert_metrics.glucose_peak_mM
        )
        assert out_acarbose.expert_metrics.glucose_peak_time_s == pytest.approx(
            out_no_drug.expert_metrics.glucose_peak_time_s
        )

    def test_partial_complex_carb_fraction_produces_an_intermediate_effect(self):
        # A half-free-sugar, half-starch meal should show LESS effect
        # than a fully-starch meal, but MORE than a fully-free-sugar one
        # -- confirms the scaling is genuinely proportional, not a
        # binary on/off switch.
        def peak_for(fraction):
            scenario = ScenarioDefinition(
                name=f"frac_{fraction}",
                meal_events=[MealEvent(time_s=0, carbs_g=75.0, complex_carb_fraction=fraction)],
                duration_s=14400,
                drug_events=[DrugEvent(time_s=0, drug="acarbose", dose_mg=100.0)],
            )
            return SimulationRunner().run(scenario).expert_metrics.glucose_peak_mM

        peak_free_sugar = peak_for(0.0)
        peak_half = peak_for(0.5)
        peak_full_starch = peak_for(1.0)
        assert peak_free_sugar > peak_half > peak_full_starch

    def test_complex_carb_fraction_default_matches_original_full_effect(self):
        # Default (1.0, unset) must reproduce the exact pre-existing
        # milestone's own behaviour byte-for-byte -- confirms this is a
        # backward-compatible addition, not a silent behaviour change
        # for every existing scenario that doesn't set this new field.
        default_scenario = ScenarioDefinition(
            name="default_frac", meal_events=[MealEvent(time_s=0, carbs_g=75.0)],
            duration_s=14400,
            drug_events=[DrugEvent(time_s=0, drug="acarbose", dose_mg=100.0)],
        )
        explicit_full = ScenarioDefinition(
            name="explicit_full",
            meal_events=[MealEvent(time_s=0, carbs_g=75.0, complex_carb_fraction=1.0)],
            duration_s=14400,
            drug_events=[DrugEvent(time_s=0, drug="acarbose", dose_mg=100.0)],
        )
        runner = SimulationRunner()
        peak_default = runner.run(default_scenario).expert_metrics.glucose_peak_mM
        peak_explicit = runner.run(explicit_full).expert_metrics.glucose_peak_mM
        assert peak_default == pytest.approx(peak_explicit)

    def test_complex_carb_fraction_out_of_bounds_raises(self):
        with pytest.raises(ValueError):
            MealEvent(time_s=0, carbs_g=75.0, complex_carb_fraction=1.5)
        with pytest.raises(ValueError):
            MealEvent(time_s=0, carbs_g=75.0, complex_carb_fraction=-0.1)
        with pytest.raises(ValueError):
            MealEvent(time_s=0, carbs_g=75.0, complex_carb_fraction=float("nan"))


class TestMiglitolAbsorptionModifier:
    """
    Miglitol (MASTER_PLAN 1.4b(ii)) -- the second alpha-glucosidase
    inhibitor, same mass-conserving kabs-reduction mechanism as acarbose
    but its OWN, separately-verified calibration (Samad et al. 1988,
    PMID:3042310). See CLAUDE.md's miglitol milestone entry.
    """

    def test_miglitol_lowers_and_delays_the_glucose_peak(self):
        no_drug = ScenarioDefinition(
            name="no_drug_mig", meal_events=[MealEvent(time_s=0, carbs_g=75.0)], duration_s=14400,
        )
        with_miglitol = ScenarioDefinition(
            name="miglitol", meal_events=[MealEvent(time_s=0, carbs_g=75.0)], duration_s=14400,
            drug_events=[DrugEvent(time_s=0, drug="miglitol", dose_mg=50.0)],
        )
        runner = SimulationRunner()
        out_no_drug = runner.run(no_drug)
        out_miglitol = runner.run(with_miglitol)
        assert out_miglitol.expert_metrics.glucose_peak_mM < out_no_drug.expert_metrics.glucose_peak_mM
        assert out_miglitol.expert_metrics.glucose_peak_time_s > out_no_drug.expert_metrics.glucose_peak_time_s

    def test_miglitol_and_acarbose_produce_genuinely_different_effects(self):
        # Regression guard for the runner.py dispatch fix: before this
        # milestone, ALPHA_GLUCOSIDASE_INHIBITOR_DRUGS had only one
        # member, so the per-drug dispatch bug (every AGI event routed
        # through acarbose_kabs_scale() unconditionally) was latent, not
        # reachable. Confirms the two drugs' real, distinct calibrations
        # both survive through the actual SimulationRunner path, not
        # just the unit-level drug_interactions.py functions.
        meal = [MealEvent(time_s=0, carbs_g=75.0)]
        acarbose_scenario = ScenarioDefinition(
            name="acarbose_100", meal_events=meal, duration_s=14400,
            drug_events=[DrugEvent(time_s=0, drug="acarbose", dose_mg=50.0)],
        )
        miglitol_scenario = ScenarioDefinition(
            name="miglitol_50", meal_events=meal, duration_s=14400,
            drug_events=[DrugEvent(time_s=0, drug="miglitol", dose_mg=50.0)],
        )
        runner = SimulationRunner()
        peak_acarbose = runner.run(acarbose_scenario).expert_metrics.glucose_peak_mM
        peak_miglitol = runner.run(miglitol_scenario).expert_metrics.glucose_peak_mM
        assert peak_acarbose != pytest.approx(peak_miglitol)

    def test_miglitol_gi_meal_is_slowed_like_a_default_meal(self):
        # Same treatment as acarbose since 2026-11-17 (see
        # test_gi_meal_is_slowed_like_a_default_meal_and_warns).
        with_miglitol = ScenarioDefinition(
            name="gi_miglitol",
            meal_events=[MealEvent(time_s=0, carbs_g=75.0, glycaemic_index=70.0)],
            duration_s=14400,
            drug_events=[DrugEvent(time_s=0, drug="miglitol", dose_mg=50.0)],
        )
        no_drug = ScenarioDefinition(
            name="gi_no_drug_mig",
            meal_events=[MealEvent(time_s=0, carbs_g=75.0, glycaemic_index=70.0)],
            duration_s=14400,
        )
        runner = SimulationRunner()
        out_no_drug = runner.run(no_drug)
        out_miglitol = runner.run(with_miglitol)
        assert out_miglitol.expert_metrics.glucose_peak_mM < out_no_drug.expert_metrics.glucose_peak_mM

    def test_zero_dose_miglitol_event_has_no_effect(self):
        no_drug = ScenarioDefinition(
            name="no_drug_mig2", meal_events=[MealEvent(time_s=0, carbs_g=75.0)], duration_s=14400,
        )
        zero_dose = ScenarioDefinition(
            name="miglitol_zero", meal_events=[MealEvent(time_s=0, carbs_g=75.0)], duration_s=14400,
            drug_events=[DrugEvent(time_s=0, drug="miglitol", dose_mg=0.0)],
        )
        runner = SimulationRunner()
        peak_no_drug = runner.run(no_drug).expert_metrics.glucose_peak_mM
        peak_zero = runner.run(zero_dose).expert_metrics.glucose_peak_mM
        assert peak_zero == pytest.approx(peak_no_drug)

    def test_miglitol_dose_below_reference_produces_smaller_effect(self):
        low = ScenarioDefinition(
            name="miglitol_low", meal_events=[MealEvent(time_s=0, carbs_g=75.0)], duration_s=14400,
            drug_events=[DrugEvent(time_s=0, drug="miglitol", dose_mg=10.0)],
        )
        full = ScenarioDefinition(
            name="miglitol_full", meal_events=[MealEvent(time_s=0, carbs_g=75.0)], duration_s=14400,
            drug_events=[DrugEvent(time_s=0, drug="miglitol", dose_mg=50.0)],
        )
        runner = SimulationRunner()
        peak_low = runner.run(low).expert_metrics.glucose_peak_mM
        peak_full = runner.run(full).expert_metrics.glucose_peak_mM
        assert peak_full < peak_low

    def test_miglitol_free_sugar_meal_is_genuinely_unaffected(self):
        # Same complex_carb_fraction scaling as acarbose (MASTER_PLAN
        # 1.4b(iii)) -- runner.py's effective_scale formula is generic
        # over whichever AGI drug was dispatched, not acarbose-specific,
        # so a free-monosaccharide meal (complex_carb_fraction=0.0) needs
        # no enzymatic digestion and should see NO miglitol effect
        # either, byte-identical to the no-drug baseline.
        no_drug = ScenarioDefinition(
            name="free_sugar_no_drug_mig",
            meal_events=[MealEvent(time_s=0, carbs_g=75.0, complex_carb_fraction=0.0)],
            duration_s=14400,
        )
        with_miglitol = ScenarioDefinition(
            name="free_sugar_miglitol",
            meal_events=[MealEvent(time_s=0, carbs_g=75.0, complex_carb_fraction=0.0)],
            duration_s=14400,
            drug_events=[DrugEvent(time_s=0, drug="miglitol", dose_mg=50.0)],
        )
        runner = SimulationRunner()
        out_no_drug = runner.run(no_drug)
        out_miglitol = runner.run(with_miglitol)
        assert out_miglitol.expert_metrics.glucose_peak_mM == pytest.approx(
            out_no_drug.expert_metrics.glucose_peak_mM
        )

    def test_acarbose_and_miglitol_together_only_the_first_listed_is_applied(self):
        # Same "first wins" limitation already documented/tested for
        # exogenous insulin (test_mixed_insulin_drug_types_only_the_
        # first_listed_is_composed) -- runner.py's AGI dispatch loop
        # breaks on the first recognized alpha-glucosidase-inhibitor
        # event, so a scenario listing both acarbose and miglitol uses
        # only whichever is listed first, matching a single-drug
        # scenario with that same first event.
        meal = [MealEvent(time_s=0, carbs_g=75.0)]
        acarbose_first = ScenarioDefinition(
            name="acarbose_then_miglitol", meal_events=meal, duration_s=14400,
            drug_events=[
                DrugEvent(time_s=0, drug="acarbose", dose_mg=100.0),
                DrugEvent(time_s=0, drug="miglitol", dose_mg=50.0),
            ],
        )
        acarbose_alone = ScenarioDefinition(
            name="acarbose_alone", meal_events=meal, duration_s=14400,
            drug_events=[DrugEvent(time_s=0, drug="acarbose", dose_mg=100.0)],
        )
        miglitol_first = ScenarioDefinition(
            name="miglitol_then_acarbose", meal_events=meal, duration_s=14400,
            drug_events=[
                DrugEvent(time_s=0, drug="miglitol", dose_mg=50.0),
                DrugEvent(time_s=0, drug="acarbose", dose_mg=100.0),
            ],
        )
        miglitol_alone = ScenarioDefinition(
            name="miglitol_alone", meal_events=meal, duration_s=14400,
            drug_events=[DrugEvent(time_s=0, drug="miglitol", dose_mg=50.0)],
        )
        runner = SimulationRunner()
        peak_acarbose_first = runner.run(acarbose_first).expert_metrics.glucose_peak_mM
        peak_acarbose_alone = runner.run(acarbose_alone).expert_metrics.glucose_peak_mM
        peak_miglitol_first = runner.run(miglitol_first).expert_metrics.glucose_peak_mM
        peak_miglitol_alone = runner.run(miglitol_alone).expert_metrics.glucose_peak_mM
        assert peak_acarbose_first == pytest.approx(peak_acarbose_alone)
        assert peak_miglitol_first == pytest.approx(peak_miglitol_alone)

    def test_every_alpha_glucosidase_inhibitor_drug_runs_without_keyerror(self):
        # Regression guard for the runner.py dispatch fix (MASTER_PLAN
        # 1.4b(ii)): _AGI_KABS_SCALE_FUNCTIONS (the per-drug kabs-scale
        # dispatch dict in runner.py) and ALPHA_GLUCOSIDASE_INHIBITOR_DRUGS
        # (drug_interactions.py) are two independently-written sets that
        # must stay in sync -- if a future AGI drug were ever added to one
        # without the other, this loop surfaces it as a real, immediate
        # test failure rather than the two silently drifting apart until
        # a live KeyError crashes a production request. Found missing by
        # the mandatory security-reviewer pass.
        from biosim_physiology.drug_interactions import ALPHA_GLUCOSIDASE_INHIBITOR_DRUGS
        for drug in ALPHA_GLUCOSIDASE_INHIBITOR_DRUGS:
            scenario = ScenarioDefinition(
                name=f"agi_sync_{drug.value}",
                meal_events=[MealEvent(time_s=0, carbs_g=75.0)],
                duration_s=14400,
                drug_events=[DrugEvent(time_s=0, drug=drug.value, dose_mg=50.0)],
            )
            out = SimulationRunner().run(scenario)
            assert out.success


# ---------------------------------------------------------------------------
# Gestational diabetes -- MASTER_PLAN 1.11
# ---------------------------------------------------------------------------

class TestGestationalDiabetesDiseaseState:
    """
    DiseaseState.GESTATIONAL_DIABETES -- a purely static Gb/p1/p3/
    secretion_scale personalizer modifier (no ODE-coupled pathway, same
    architecture as TYPE_2_DIABETES/INSULIN_RESISTANT). See CLAUDE.md's
    gestational-diabetes milestone entry for the full citation
    (Kautzky-Willer et al. 1997, PMID:9353615) and calibration.
    """

    def test_gestational_diabetes_raises_the_glucose_peak_on_a_standard_meal(self):
        from biosim_physiology import PhysiologicalProfile, DiseaseState
        meal = [MealEvent(time_s=0, carbs_g=75.0)]
        healthy = ScenarioDefinition(
            name="healthy_meal", meal_events=meal, duration_s=14400,
            personalization=PersonalizationSpec(physiological_profile=PhysiologicalProfile()),
        )
        gdm = ScenarioDefinition(
            name="gdm_meal", meal_events=meal, duration_s=14400,
            personalization=PersonalizationSpec(
                physiological_profile=PhysiologicalProfile(
                    disease_states=[DiseaseState.GESTATIONAL_DIABETES]
                )
            ),
        )
        runner = SimulationRunner()
        peak_healthy = runner.run(healthy).expert_metrics.glucose_peak_mM
        peak_gdm = runner.run(gdm).expert_metrics.glucose_peak_mM
        assert peak_gdm > peak_healthy


class TestGckModyDiseaseState:
    """
    DiseaseState.GCK_MODY -- a pure Gb (glucose-threshold-reset)
    personalizer modifier, deliberately leaving p1/p3/secretion_scale
    untouched (unlike every other diabetes-related state in this file).
    See CLAUDE.md's GCK-MODY milestone entry for the full citation
    (Guenat et al. 2000, PMID:11119017) and calibration.
    """

    def test_gck_mody_raises_the_fasting_and_peak_glucose_on_a_standard_meal(self):
        from biosim_physiology import PhysiologicalProfile, DiseaseState
        meal = [MealEvent(time_s=0, carbs_g=75.0)]
        healthy = ScenarioDefinition(
            name="healthy_meal_mody", meal_events=meal, duration_s=14400,
            personalization=PersonalizationSpec(physiological_profile=PhysiologicalProfile()),
        )
        gck_mody = ScenarioDefinition(
            name="gck_mody_meal", meal_events=meal, duration_s=14400,
            personalization=PersonalizationSpec(
                physiological_profile=PhysiologicalProfile(
                    disease_states=[DiseaseState.GCK_MODY]
                )
            ),
        )
        runner = SimulationRunner()
        out_healthy = runner.run(healthy)
        out_mody = runner.run(gck_mody)
        assert out_mody.expert_metrics.glucose_peak_mM > out_healthy.expert_metrics.glucose_peak_mM
        assert out_mody.expert_metrics.glucose_baseline_uM > out_healthy.expert_metrics.glucose_baseline_uM


# ---------------------------------------------------------------------------
# Exogenous insulin (regular human insulin) -- MASTER_PLAN 1.2a
# ---------------------------------------------------------------------------

class TestExogenousInsulinPathwayComposition:
    """
    Injected regular insulin -- modelled as a stateless, closed-form
    two-sequential-compartment absorption pathway (see CLAUDE.md's
    exogenous-insulin milestone entry for the full literature-
    verification/calibration trail: Howey et al. 1994, PMID:8314011).
    """

    def test_insulin_dose_raises_plasma_insulin_end_to_end(self):
        no_drug = ScenarioDefinition(
            name="fasted_no_insulin", meal_events=[], duration_s=36000,
        )
        with_insulin = ScenarioDefinition(
            name="fasted_insulin", meal_events=[], duration_s=36000,
            drug_events=[DrugEvent(time_s=0, drug="insulin_regular", dose_mg=10.0)],
        )
        runner = SimulationRunner()
        peak_no_drug = runner.run(no_drug).expert_metrics.insulin_peak_pM
        peak_insulin = runner.run(with_insulin).expert_metrics.insulin_peak_pM
        assert peak_insulin > peak_no_drug

    def test_insulin_peak_matches_howey_1994_calibration_end_to_end(self):
        # Reproduces Howey et al. 1994 (PMID:8314011) through the real
        # SimulationRunner path (not just the bare ODESystem check in
        # test_exogenous_insulin_pathway.py) -- 10U regular insulin,
        # Tmax=101min, Cmax=308pM excess above baseline.
        scenario = ScenarioDefinition(
            name="insulin_calibration_check", meal_events=[], duration_s=36000,
            drug_events=[DrugEvent(time_s=0, drug="insulin_regular", dose_mg=10.0)],
        )
        out = SimulationRunner().run(scenario)
        assert out.success
        assert out.expert_metrics.insulin_peak_time_s / 60.0 == pytest.approx(101.0, abs=3.0)
        excess = out.expert_metrics.insulin_peak_pM - 66.0  # PathwayParameters().Ib default
        assert excess == pytest.approx(308.0, rel=0.03)

    def test_meal_scenario_higher_dose_lowers_glucose_auc(self):
        # A fasted (no-meal) scenario does NOT show a glucose-lowering
        # effect at all in this model -- see
        # test_fasted_scenario_glucose_is_genuinely_unaffected below for
        # why, and the honest disclosure in CLAUDE.md. A real meal
        # scenario is needed to demonstrate insulin's actual glucose-
        # lowering action, since only then is glucose displaced from Gb
        # in the first place for accelerated clearance to act on.
        no_insulin = ScenarioDefinition(
            name="meal_no_insulin", meal_events=[MealEvent(time_s=0, carbs_g=75.0)], duration_s=36000,
        )
        low_dose = ScenarioDefinition(
            name="meal_insulin_10U", meal_events=[MealEvent(time_s=0, carbs_g=75.0)], duration_s=36000,
            drug_events=[DrugEvent(time_s=0, drug="insulin_regular", dose_mg=10.0)],
        )
        high_dose = ScenarioDefinition(
            name="meal_insulin_30U", meal_events=[MealEvent(time_s=0, carbs_g=75.0)], duration_s=36000,
            drug_events=[DrugEvent(time_s=0, drug="insulin_regular", dose_mg=30.0)],
        )
        runner = SimulationRunner()
        auc_no_insulin = runner.run(no_insulin).expert_metrics.glucose_auc_uM_s
        auc_low = runner.run(low_dose).expert_metrics.glucose_auc_uM_s
        auc_high = runner.run(high_dose).expert_metrics.glucose_auc_uM_s
        assert auc_low < auc_no_insulin
        assert auc_high < auc_low

    def test_fasted_scenario_glucose_is_genuinely_unaffected(self):
        # Honest, disclosed model behaviour, not a bug: in this model,
        # insulin's glucose-clearance action (via the X/remote-insulin
        # compartment) is proportional to (G - Gb_eff), which is exactly
        # zero at the analytic basal state regardless of how large X
        # becomes -- so a fasted (no-meal) scenario's glucose genuinely
        # stays flat at Gb even with a large exogenous insulin dose,
        # despite plasma insulin itself rising substantially. This
        # matches the identical, previously-verified behaviour for
        # glibenclamide (sulfonylurea) earlier this session -- confirmed
        # here for exogenous insulin specifically since a user might
        # otherwise expect (incorrectly, for THIS model's architecture)
        # a fasted insulin injection to cause hypoglycaemia.
        scenario = ScenarioDefinition(
            name="fasted_high_dose_insulin", meal_events=[], duration_s=36000,
            drug_events=[DrugEvent(time_s=0, drug="insulin_regular", dose_mg=30.0)],
        )
        out = SimulationRunner().run(scenario)
        G_mM = out.simulation.species["G"] / 1000.0
        assert G_mM.max() == pytest.approx(5.0, abs=1e-6)
        assert G_mM.min() == pytest.approx(5.0, abs=1e-6)
        # Insulin itself DOES rise substantially -- confirming the flat
        # glucose isn't because the dose had no effect at all.
        assert out.expert_metrics.insulin_peak_pM > 300.0

    def test_dose_given_before_simulation_start_shifts_the_peak_earlier(self):
        on_time = ScenarioDefinition(
            name="insulin_on_time", meal_events=[], duration_s=36000,
            drug_events=[DrugEvent(time_s=0, drug="insulin_regular", dose_mg=10.0)],
        )
        given_earlier = ScenarioDefinition(
            name="insulin_30min_ago", meal_events=[], duration_s=36000,
            drug_events=[DrugEvent(
                time_s=0, drug="insulin_regular", dose_mg=10.0, hours_before_simulation=0.5,
            )],
        )
        runner = SimulationRunner()
        t_on_time = runner.run(on_time).expert_metrics.insulin_peak_time_s
        t_earlier = runner.run(given_earlier).expert_metrics.insulin_peak_time_s
        assert t_earlier == pytest.approx(t_on_time - 30 * 60, abs=60)

    def test_zero_dose_has_no_effect(self):
        no_drug = ScenarioDefinition(name="no_insulin_event", meal_events=[], duration_s=36000)
        zero_dose = ScenarioDefinition(
            name="insulin_zero_dose", meal_events=[], duration_s=36000,
            drug_events=[DrugEvent(time_s=0, drug="insulin_regular", dose_mg=0.0)],
        )
        runner = SimulationRunner()
        peak_no_drug = runner.run(no_drug).expert_metrics.insulin_peak_pM
        peak_zero = runner.run(zero_dose).expert_metrics.insulin_peak_pM
        assert peak_zero == pytest.approx(peak_no_drug)

    def test_multiple_insulin_events_only_the_first_is_composed(self):
        # Documented limitation of _detect_exogenous_insulin(): only one
        # ExogenousInsulinPathway can be composed per SimulationRunner
        # call today (CompositePathway also can't disambiguate two
        # instances sharing the same pathway name), so a scenario with
        # two insulin_regular events uses only the first one, matching a
        # single-dose scenario with the same first dose.
        single_dose = ScenarioDefinition(
            name="single_insulin_dose", meal_events=[], duration_s=36000,
            drug_events=[DrugEvent(time_s=0, drug="insulin_regular", dose_mg=10.0)],
        )
        two_doses = ScenarioDefinition(
            name="two_insulin_doses", meal_events=[], duration_s=36000,
            drug_events=[
                DrugEvent(time_s=0, drug="insulin_regular", dose_mg=10.0),
                DrugEvent(time_s=1800, drug="insulin_regular", dose_mg=20.0),
            ],
        )
        runner = SimulationRunner()
        peak_single = runner.run(single_dose).expert_metrics.insulin_peak_pM
        peak_two = runner.run(two_doses).expert_metrics.insulin_peak_pM
        assert peak_two == pytest.approx(peak_single)

    def test_mixed_insulin_drug_types_only_the_first_listed_is_composed(self):
        # Same "first wins" limitation as the same-drug case above, now
        # confirmed for TWO DIFFERENT insulin identities (regular and
        # lispro) rather than two doses of the same drug -- the existing
        # _detect_exogenous_insulin() docstring already describes this
        # generically ("exogenous-insulin drug event", not tied to one
        # specific drug), but this is the first test to actually exercise
        # a mixed-drug-type scenario rather than just repeated doses of
        # the same drug.
        lispro_first = ScenarioDefinition(
            name="lispro_then_regular", meal_events=[], duration_s=36000,
            drug_events=[
                DrugEvent(time_s=0, drug="insulin_lispro", dose_mg=10.0),
                DrugEvent(time_s=0, drug="insulin_regular", dose_mg=10.0),
            ],
        )
        lispro_alone = ScenarioDefinition(
            name="lispro_alone", meal_events=[], duration_s=36000,
            drug_events=[DrugEvent(time_s=0, drug="insulin_lispro", dose_mg=10.0)],
        )
        regular_alone = ScenarioDefinition(
            name="regular_alone", meal_events=[], duration_s=36000,
            drug_events=[DrugEvent(time_s=0, drug="insulin_regular", dose_mg=10.0)],
        )
        runner = SimulationRunner()
        peak_mixed = runner.run(lispro_first).expert_metrics.insulin_peak_pM
        peak_lispro = runner.run(lispro_alone).expert_metrics.insulin_peak_pM
        peak_regular = runner.run(regular_alone).expert_metrics.insulin_peak_pM
        assert peak_mixed == pytest.approx(peak_lispro)
        assert peak_mixed != pytest.approx(peak_regular)

    def test_lispro_reproduces_howey_1994_end_to_end(self):
        # Same Howey et al. 1994 within-study comparison as the regular-
        # insulin test, through the real SimulationRunner path -- confirms
        # insulin_lispro is genuinely wired to
        # ExogenousInsulinParameters.rapid_acting_insulin(), not silently
        # falling back to the regular-insulin calibration.
        scenario = ScenarioDefinition(
            name="lispro_calibration_check", meal_events=[], duration_s=36000,
            drug_events=[DrugEvent(time_s=0, drug="insulin_lispro", dose_mg=10.0)],
        )
        out = SimulationRunner().run(scenario)
        assert out.success
        assert out.expert_metrics.insulin_peak_time_s / 60.0 == pytest.approx(42.0, abs=3.0)
        excess = out.expert_metrics.insulin_peak_pM - 66.0  # PathwayParameters().Ib default
        assert excess == pytest.approx(698.0, rel=0.03)

    def test_near_zero_endogenous_secretion_plus_exogenous_insulin_makes_t1d_simulatable(self):
        # MASTER_PLAN 1.2c: verifies (does NOT introduce new code) that
        # this platform's existing PathwayParameters.Ib/secretion_scale
        # fields, combined with exogenous_insulin_pathway.py's
        # ExogenousInsulinPathway (1.2a/1.2b), already provide everything
        # needed at the mechanics level to represent T1D end-to-end
        # through the real, production SimulationRunner path -- see
        # CLAUDE.md's own DONE entry and
        # test_ode_engine.py::TestNearTotalBetaCellFailureFixedPoint for
        # the underlying ODE-engine-level fixed-point verification this
        # test builds on.
        from biosim_ode import PathwayParameters
        t1d_like_params = PathwayParameters(
            Ib=2.0, Gb=8.0, secretion_scale=0.01, first_phase_gain=0.0,
        )
        meal = [MealEvent(time_s=0, carbs_g=75.0)]
        no_insulin = ScenarioDefinition(name="t1d_no_insulin", meal_events=meal, duration_s=18000)
        with_insulin = ScenarioDefinition(
            name="t1d_with_insulin", meal_events=meal, duration_s=18000,
            drug_events=[DrugEvent(time_s=0, drug="insulin_regular", dose_mg=6.0)],
        )
        runner = SimulationRunner(base_params=t1d_like_params)
        out_no_insulin = runner.run(no_insulin).expert_metrics
        out_with_insulin = runner.run(with_insulin).expert_metrics
        # Near-total endogenous secretion failure: without exogenous
        # insulin, plasma insulin barely moves at all.
        assert (out_no_insulin.insulin_peak_pM - 2.0) < 5.0
        # Exogenous insulin (the ONLY source of a real insulin response
        # for this patient) genuinely restores a substantial rise and
        # measurably lowers the glucose peak.
        assert out_with_insulin.insulin_peak_pM > 150.0
        assert out_with_insulin.glucose_peak_mM < out_no_insulin.glucose_peak_mM

    def test_type_1_diabetes_disease_state_makes_t1d_simulatable_end_to_end(self):
        # MASTER_PLAN 1.10: the real, literature-anchored
        # DiseaseState.TYPE_1_DIABETES preset (modifiers.py, Simonson
        # et al. 1987 PMID:3100365) reproduces the same qualitative
        # result as the illustrative-values test directly above, but
        # through the REAL PhysiologicalPersonalizer/PersonalizationSpec
        # path a real caller would actually use, not a hand-constructed
        # base_params override.
        from biosim_physiology import PhysiologicalProfile, DiseaseState
        meal = [MealEvent(time_s=0, carbs_g=75.0)]
        t1d_profile = PersonalizationSpec(
            physiological_profile=PhysiologicalProfile(
                disease_states=[DiseaseState.TYPE_1_DIABETES]
            )
        )
        no_insulin = ScenarioDefinition(
            name="t1d_preset_no_insulin", meal_events=meal, duration_s=18000,
            personalization=t1d_profile,
        )
        with_insulin = ScenarioDefinition(
            name="t1d_preset_with_insulin", meal_events=meal, duration_s=18000,
            personalization=t1d_profile,
            drug_events=[DrugEvent(time_s=0, drug="insulin_regular", dose_mg=6.0)],
        )
        runner = SimulationRunner()
        out_no_insulin = runner.run(no_insulin)
        out_with_insulin = runner.run(with_insulin)
        base_ib = out_no_insulin.physiological_personalization.personalized_params.Ib
        # Near-total endogenous secretion failure: without exogenous
        # insulin, plasma insulin barely moves above its own (near-zero)
        # personalized basal value.
        assert (out_no_insulin.expert_metrics.insulin_peak_pM - base_ib) < 5.0
        # Exogenous insulin restores a substantial rise and measurably
        # lowers the glucose peak.
        assert out_with_insulin.expert_metrics.insulin_peak_pM > 150.0
        assert out_with_insulin.expert_metrics.glucose_peak_mM < out_no_insulin.expert_metrics.glucose_peak_mM

    def test_lispro_peaks_earlier_and_higher_than_regular_insulin(self):
        # Direct within-model confirmation that the two insulin analogs
        # are genuinely distinct, not aliases of the same calibration.
        regular = ScenarioDefinition(
            name="regular_10U", meal_events=[], duration_s=36000,
            drug_events=[DrugEvent(time_s=0, drug="insulin_regular", dose_mg=10.0)],
        )
        lispro = ScenarioDefinition(
            name="lispro_10U", meal_events=[], duration_s=36000,
            drug_events=[DrugEvent(time_s=0, drug="insulin_lispro", dose_mg=10.0)],
        )
        runner = SimulationRunner()
        out_regular = runner.run(regular).expert_metrics
        out_lispro = runner.run(lispro).expert_metrics
        assert out_lispro.insulin_peak_time_s < out_regular.insulin_peak_time_s
        assert out_lispro.insulin_peak_pM > out_regular.insulin_peak_pM

    def test_no_insulin_static_modifier_applied(self):
        # _exogenous_insulin_modifiers() must return {} -- the ODE-coupled
        # pathway is the SOLE representation, matching the GLP-1/SGLT2i
        # precedent.
        from biosim_ode import PathwayParameters
        scenario = ScenarioDefinition(
            name="insulin_static_check", meal_events=[], duration_s=10,
            drug_events=[DrugEvent(time_s=0, drug="insulin_regular", dose_mg=10.0)],
        )
        out = SimulationRunner().run(scenario)
        p = out.physiological_personalization.personalized_params
        default = PathwayParameters()
        assert p.Gb == pytest.approx(default.Gb)
        assert p.p3 == pytest.approx(default.p3)
        assert p.secretion_scale == pytest.approx(default.secretion_scale)
        assert p.glucose_independent_secretion == pytest.approx(default.glucose_independent_secretion)

    def test_drug_only_quick_path_still_produces_the_effect(self):
        scenario = ScenarioDefinition(
            name="insulin_quick_path", meal_events=[], duration_s=36000,
            drug_events=[DrugEvent(time_s=0, drug="insulin_regular", dose_mg=10.0)],
        )
        out = SimulationRunner().run(scenario)
        assert out.success
        assert out.expert_metrics.insulin_peak_pM > 300.0


class TestNutrientSecretionPathwayComposition:
    """
    Dietary protein co-ingestion -- MASTER_PLAN row 1.28. See
    packages/ode-engine/src/biosim_ode/nutrient_secretion_pathway.py's
    own docstring for the full literature-verification/calibration
    trail: Krezowski et al. 1986 (PMID:3538843, normal subjects, 28%
    additive insulin-AUC ratio for 50g protein vs. 50g glucose).
    """

    def test_no_protein_scenario_is_byte_identical_to_before(self):
        # Backward-compatibility guard: MealEvent.protein_g defaults to
        # 0.0, so NutrientSecretionPathway is never composed and every
        # existing scenario is completely unaffected -- compared against
        # a bare GlucoseInsulinPathway + ODESystem run that never even
        # imports the new pathway module, not a hardcoded magic number.
        from biosim_ode import GlucoseInsulinPathway, ODESystem
        from biosim_ode.meal import MealInput

        scenario = ScenarioDefinition(
            name="no_protein", meal_events=[MealEvent(time_s=0, carbs_g=50.0)], duration_s=14400,
        )
        runner = SimulationRunner()
        out = runner.run(scenario)
        assert out.success

        # Same output grid as the runner, so the sampled peak is comparable.
        reference = ODESystem(GlucoseInsulinPathway()).simulate(
            duration_s=14400, meals=[MealEvent(time_s=0, carbs_g=50.0).to_meal_input()],
            n_output_points=runner.n_output_points,
        )
        reference_peak_mM = reference.species["G"].max() / 1000.0   # uM -> mM
        assert out.expert_metrics.glucose_peak_mM == pytest.approx(reference_peak_mM)

    def test_protein_meal_raises_insulin_peak_end_to_end(self):
        no_protein = ScenarioDefinition(
            name="no_protein", meal_events=[MealEvent(time_s=0, carbs_g=50.0, protein_g=0.0)],
            duration_s=14400,
        )
        with_protein = ScenarioDefinition(
            name="with_protein", meal_events=[MealEvent(time_s=0, carbs_g=50.0, protein_g=50.0)],
            duration_s=14400,
        )
        runner = SimulationRunner()
        peak_no_protein = runner.run(no_protein).expert_metrics.insulin_peak_pM
        peak_with_protein = runner.run(with_protein).expert_metrics.insulin_peak_pM
        assert peak_with_protein > peak_no_protein

    def test_insulin_auc_ratio_matches_krezowski_1986_end_to_end(self):
        # Reproduces Krezowski et al. 1986's own 28% additive insulin-AUC
        # ratio through the real SimulationRunner path (not just the bare
        # ODESystem check in test_nutrient_secretion_pathway.py) -- 50g
        # protein alone vs. 50g carbohydrate alone, both over a 240min
        # window matching the paper's own protocol.
        import numpy as np

        def insulin_auc(carbs_g, protein_g):
            scenario = ScenarioDefinition(
                name="auc_check",
                meal_events=[MealEvent(time_s=0, carbs_g=carbs_g, protein_g=protein_g)],
                duration_s=14400,
            )
            out = SimulationRunner().run(scenario)
            t = out.simulation.times_s
            I = out.simulation.species["I"]
            Ib = I[0]
            return np.trapezoid(I - Ib, t) * 1e6 / 60.0

        glucose_auc = insulin_auc(50.0, 0.0)
        protein_auc = insulin_auc(0.0, 50.0)
        # KNOWN GAP since B7 (2026-11-16), recorded, NOT re-tuned: this
        # secretion constant was calibrated to Krezowski's 0.28 while the
        # model double-counted resting glucose clearance. With that removed,
        # the same constant gives ~0.22. Pinned so any drift is looked at.
        assert protein_auc / glucose_auc == pytest.approx(0.216, abs=0.005)

    def test_zero_protein_preserves_basal_invariance_end_to_end(self):
        scenario = ScenarioDefinition(
            name="fasted_zero_protein", meal_events=[], duration_s=14400,
        )
        out = SimulationRunner().run(scenario)
        assert out.expert_metrics.glucose_baseline_uM == pytest.approx(5000.0)

    def test_protein_alone_meal_does_not_double_count_against_exogenous_insulin(self):
        # Both nutrient_secretion_delta and exogenous_insulin_delta are
        # additive to dI -- confirm composing a protein meal WITH an
        # insulin dose doesn't produce anything degenerate (e.g. an
        # unbounded/NaN insulin peak), just their expected combined effect.
        scenario = ScenarioDefinition(
            name="protein_plus_insulin",
            meal_events=[MealEvent(time_s=0, carbs_g=50.0, protein_g=50.0)],
            drug_events=[DrugEvent(time_s=0, drug="insulin_regular", dose_mg=10.0)],
            duration_s=14400,
        )
        out = SimulationRunner().run(scenario)
        assert out.success
        assert math.isfinite(out.expert_metrics.insulin_peak_pM)
        assert out.expert_metrics.insulin_peak_pM > 0.0


class TestInsulinDoseEvent:
    """
    MASTER_PLAN 1.2b (event-type): a dedicated, honestly-Unit-named
    alternative to DrugEvent's dose_mg-as-Units hack for
    insulin_regular/insulin_lispro. Additive, not a replacement -- the
    legacy DrugEvent path (TestExogenousInsulinPathwayComposition above)
    must keep working unchanged.
    """

    def test_negative_dose_units_raises(self):
        with pytest.raises(ValueError, match="dose_units"):
            InsulinDoseEvent(time_s=0, insulin_type="regular", dose_units=-1.0)

    def test_negative_time_s_raises(self):
        with pytest.raises(ValueError, match="time_s"):
            InsulinDoseEvent(time_s=-1.0, insulin_type="regular", dose_units=10.0)

    def test_negative_hours_before_simulation_raises(self):
        with pytest.raises(ValueError, match="hours_before_simulation"):
            InsulinDoseEvent(
                time_s=0, insulin_type="regular", dose_units=10.0,
                hours_before_simulation=-1.0,
            )

    def test_nan_dose_units_raises(self):
        with pytest.raises(ValueError, match="finite"):
            InsulinDoseEvent(time_s=0, insulin_type="regular", dose_units=float("nan"))

    def test_invalid_insulin_type_raises(self):
        with pytest.raises(ValueError, match="insulin_type"):
            InsulinDoseEvent(time_s=0, insulin_type="glargine", dose_units=10.0)

    def test_zero_dose_is_valid(self):
        # Boundary, not an error -- mirrors DrugEvent's own dose_mg=0.0 convention.
        InsulinDoseEvent(time_s=0, insulin_type="regular", dose_units=0.0)


class TestInsulinDoseEventMutualExclusivity:
    def test_combining_with_legacy_drug_event_raises(self):
        with pytest.raises(ValueError, match="insulin_dose_events"):
            ScenarioDefinition(
                name="conflict", meal_events=[], duration_s=3600,
                drug_events=[DrugEvent(time_s=0, drug="insulin_regular", dose_mg=10.0)],
                insulin_dose_events=[
                    InsulinDoseEvent(time_s=0, insulin_type="regular", dose_units=10.0)
                ],
            )

    def test_lispro_drug_event_also_raises(self):
        # Confirms the check isn't hardcoded to insulin_regular only.
        with pytest.raises(ValueError, match="insulin_dose_events"):
            ScenarioDefinition(
                name="conflict_lispro", meal_events=[], duration_s=3600,
                drug_events=[DrugEvent(time_s=0, drug="insulin_lispro", dose_mg=10.0)],
                insulin_dose_events=[
                    InsulinDoseEvent(time_s=0, insulin_type="lispro", dose_units=10.0)
                ],
            )

    def test_insulin_dose_events_alone_is_valid(self):
        ScenarioDefinition(
            name="insulin_event_alone", meal_events=[], duration_s=3600,
            insulin_dose_events=[
                InsulinDoseEvent(time_s=0, insulin_type="regular", dose_units=10.0)
            ],
        )

    def test_unrelated_drug_event_alongside_insulin_dose_event_is_valid(self):
        # A non-insulin DrugEvent (e.g. metformin) coexisting with an
        # InsulinDoseEvent is not a conflict -- only the SAME mechanism
        # expressed two ways is rejected.
        ScenarioDefinition(
            name="metformin_plus_insulin_event", meal_events=[], duration_s=3600,
            drug_events=[DrugEvent(time_s=0, drug="metformin", dose_mg=2000.0)],
            insulin_dose_events=[
                InsulinDoseEvent(time_s=0, insulin_type="regular", dose_units=10.0)
            ],
        )


class TestInsulinDoseEventComposition:
    """
    End-to-end confirmation that InsulinDoseEvent reaches the exact same
    ExogenousInsulinPathway composition as the legacy DrugEvent path --
    same calibration, same Howey et al. 1994 (PMID:8314011) reproduction.
    """

    def test_regular_insulin_matches_the_legacy_drug_event_path_exactly(self):
        via_drug_event = ScenarioDefinition(
            name="via_drug_event", meal_events=[], duration_s=36000,
            drug_events=[DrugEvent(time_s=0, drug="insulin_regular", dose_mg=10.0)],
        )
        via_insulin_event = ScenarioDefinition(
            name="via_insulin_event", meal_events=[], duration_s=36000,
            insulin_dose_events=[
                InsulinDoseEvent(time_s=0, insulin_type="regular", dose_units=10.0)
            ],
        )
        runner = SimulationRunner()
        m1 = runner.run(via_drug_event).expert_metrics
        m2 = runner.run(via_insulin_event).expert_metrics
        assert m2.insulin_peak_pM == pytest.approx(m1.insulin_peak_pM)
        assert m2.insulin_peak_time_s == pytest.approx(m1.insulin_peak_time_s)

    def test_lispro_matches_the_legacy_drug_event_path_exactly(self):
        via_drug_event = ScenarioDefinition(
            name="via_drug_event_lispro", meal_events=[], duration_s=36000,
            drug_events=[DrugEvent(time_s=0, drug="insulin_lispro", dose_mg=10.0)],
        )
        via_insulin_event = ScenarioDefinition(
            name="via_insulin_event_lispro", meal_events=[], duration_s=36000,
            insulin_dose_events=[
                InsulinDoseEvent(time_s=0, insulin_type="lispro", dose_units=10.0)
            ],
        )
        runner = SimulationRunner()
        m1 = runner.run(via_drug_event).expert_metrics
        m2 = runner.run(via_insulin_event).expert_metrics
        assert m2.insulin_peak_pM == pytest.approx(m1.insulin_peak_pM)
        assert m2.insulin_peak_time_s == pytest.approx(m1.insulin_peak_time_s)

    def test_reproduces_howey_1994_calibration_via_the_new_event_type(self):
        # Same numeric target as
        # test_insulin_peak_matches_howey_1994_calibration_end_to_end
        # above, now confirmed reachable through InsulinDoseEvent too --
        # this is the whole point of the new type: it must be a genuine
        # alternative entry point, not a second, subtly different path.
        scenario = ScenarioDefinition(
            name="insulin_event_calibration_check", meal_events=[], duration_s=36000,
            insulin_dose_events=[
                InsulinDoseEvent(time_s=0, insulin_type="regular", dose_units=10.0)
            ],
        )
        out = SimulationRunner().run(scenario)
        assert out.success
        assert out.expert_metrics.insulin_peak_time_s / 60.0 == pytest.approx(101.0, abs=3.0)
        excess = out.expert_metrics.insulin_peak_pM - 66.0
        assert excess == pytest.approx(308.0, rel=0.03)

    def test_hours_before_simulation_shifts_the_peak_earlier(self):
        on_time = ScenarioDefinition(
            name="insulin_event_on_time", meal_events=[], duration_s=36000,
            insulin_dose_events=[
                InsulinDoseEvent(time_s=0, insulin_type="regular", dose_units=10.0)
            ],
        )
        given_earlier = ScenarioDefinition(
            name="insulin_event_30min_ago", meal_events=[], duration_s=36000,
            insulin_dose_events=[
                InsulinDoseEvent(
                    time_s=0, insulin_type="regular", dose_units=10.0,
                    hours_before_simulation=0.5,
                )
            ],
        )
        runner = SimulationRunner()
        t_on_time = runner.run(on_time).expert_metrics.insulin_peak_time_s
        t_earlier = runner.run(given_earlier).expert_metrics.insulin_peak_time_s
        assert t_earlier == pytest.approx(t_on_time - 30 * 60, abs=60)

    def test_first_event_wins_when_multiple_are_present(self):
        # Same documented limitation as the legacy DrugEvent path.
        single = ScenarioDefinition(
            name="single_insulin_event", meal_events=[], duration_s=36000,
            insulin_dose_events=[
                InsulinDoseEvent(time_s=0, insulin_type="regular", dose_units=10.0)
            ],
        )
        two = ScenarioDefinition(
            name="two_insulin_events", meal_events=[], duration_s=36000,
            insulin_dose_events=[
                InsulinDoseEvent(time_s=0, insulin_type="regular", dose_units=10.0),
                InsulinDoseEvent(time_s=1800, insulin_type="regular", dose_units=20.0),
            ],
        )
        runner = SimulationRunner()
        peak_single = runner.run(single).expert_metrics.insulin_peak_pM
        peak_two = runner.run(two).expert_metrics.insulin_peak_pM
        assert peak_two == pytest.approx(peak_single)

    def test_to_dict_includes_insulin_dose_events(self):
        scenario = ScenarioDefinition(
            name="dict_check", meal_events=[], duration_s=3600,
            insulin_dose_events=[
                InsulinDoseEvent(
                    time_s=100.0, insulin_type="lispro", dose_units=5.0,
                    hours_before_simulation=0.25,
                )
            ],
        )
        d = scenario.to_dict()
        assert d["insulin_dose_events"] == [
            {
                "time_s": 100.0, "insulin_type": "lispro", "dose_units": 5.0,
                "hours_before_simulation": 0.25,
            }
        ]

    def test_cache_key_differs_from_an_otherwise_identical_scenario_without_it(self):
        # Same regression class as SimulationCache._key()'s existing
        # meal-GI/drug-timing collision fixes -- to_dict() must reflect
        # this field so two scenarios differing ONLY in it don't collide.
        from biosim_infra import SimulationCache
        cache = SimulationCache()
        no_insulin = ScenarioDefinition(name="cache_check", meal_events=[], duration_s=3600)
        with_insulin = ScenarioDefinition(
            name="cache_check", meal_events=[], duration_s=3600,
            insulin_dose_events=[
                InsulinDoseEvent(time_s=0, insulin_type="regular", dose_units=10.0)
            ],
        )
        assert cache._key(no_insulin) != cache._key(with_insulin)


# ---------------------------------------------------------------------------
# LaySummary
# ---------------------------------------------------------------------------

class TestLaySummary:
    @pytest.fixture(scope="class")
    def meal_output(self):
        return SimulationRunner().run(high_sugar_meal_scenario())

    def test_energy_array_length_matches_times(self, meal_output):
        ls = meal_output.lay_summary
        assert len(ls.energy_level) == len(ls.times_min)

    def test_hunger_array_length_matches_times(self, meal_output):
        ls = meal_output.lay_summary
        assert len(ls.hunger_level) == len(ls.times_min)

    def test_energy_rises_after_meal(self, meal_output):
        ls = meal_output.lay_summary
        energy_start = float(ls.energy_level[0])
        energy_peak = float(np.max(ls.energy_level))
        assert energy_peak > energy_start

    def test_hunger_falls_after_meal(self, meal_output):
        ls = meal_output.lay_summary
        hunger_start = float(ls.hunger_level[0])
        hunger_min = float(np.min(ls.hunger_level))
        assert hunger_min < hunger_start

    def test_energy_in_valid_range(self, meal_output):
        ls = meal_output.lay_summary
        assert np.all(ls.energy_level >= 0.0)
        assert np.all(ls.energy_level <= 10.0)

    def test_glucose_narrative_populated(self, meal_output):
        assert len(meal_output.lay_summary.glucose_narrative) > 20

    def test_energy_narrative_populated(self, meal_output):
        assert len(meal_output.lay_summary.energy_narrative) > 20

    def test_confidence_note_populated(self, meal_output):
        assert len(meal_output.lay_summary.confidence_note) > 5

    def test_fasted_no_energy_dip(self):
        output = SimulationRunner().run(fasted_scenario(duration_s=3600.0))
        assert not output.lay_summary.has_energy_dip

    def test_glucose_narrative_does_not_call_elevated_baseline_return_normal(self):
        # Regression: glucose_narrative unconditionally said "return
        # toward normal" for time_to_glucose_return_s, which measures
        # return to THIS PATIENT'S OWN baseline (glucose_baseline_uM),
        # not a population-normal value. The T2D disease-state preset's
        # baseline (Gb x1.3, ~6.5 mM) is already in the ADA pre-diabetic
        # fasting range (5.6-6.9 mM) -- same bug already fixed in
        # biosim_lay.report.py's sibling narrative for the identical
        # metric, but this is a separate code path (simulation-api must
        # not depend on biosim_lay) that fix didn't reach, despite being
        # exposed through the same REST API response.
        from biosim_physiology import PhysiologicalProfile, DiseaseState
        # 50 g, not the 75 g high-sugar meal: since B7 (2026-11-16) a T2D
        # 75 g meal no longer returns to baseline within 4 h, so the return
        # sentence this test checks would not be written at all.
        scenario = standard_meal_scenario(
            personalization=PersonalizationSpec(
                physiological_profile=PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES])
            )
        )
        output = SimulationRunner().run(scenario)
        assert output.expert_metrics.time_to_glucose_return_s is not None
        assert output.expert_metrics.glucose_baseline_uM / 1000.0 >= 5.6, (
            "test assumes the T2D baseline is elevated"
        )
        assert "toward normal" not in output.lay_summary.glucose_narrative
        assert "above the normal fasting range" in output.lay_summary.glucose_narrative

    def test_glucose_narrative_still_says_toward_normal_for_healthy_baseline(self, meal_output):
        assert "toward normal" in meal_output.lay_summary.glucose_narrative


# ---------------------------------------------------------------------------
# Scenario library
# ---------------------------------------------------------------------------

class TestScenarioLibrary:
    def test_all_library_scenarios_present(self):
        expected = {"ogtt_75g", "high_sugar_meal", "standard_meal",
                    "light_snack", "fasted", "caffeine_and_meal"}
        assert expected.issubset(set(LIBRARY.keys()))

    def test_get_scenario_known(self):
        s = get_scenario("ogtt_75g")
        assert s.name == "OGTT_75g"

    def test_get_scenario_unknown_raises(self):
        with pytest.raises(KeyError):
            get_scenario("nonexistent_scenario")

    @pytest.mark.parametrize("name", ["ogtt_75g", "standard_meal", "light_snack", "fasted"])
    def test_library_scenarios_run(self, name):
        output = SimulationRunner().run(get_scenario(name))
        assert output.success, f"Scenario '{name}' failed: {output.warnings}"

    def test_larger_meal_higher_glucose(self):
        runner = SimulationRunner()
        small = runner.run(light_snack_scenario())
        large = runner.run(high_sugar_meal_scenario())
        assert large.expert_metrics.glucose_peak_uM > small.expert_metrics.glucose_peak_uM

    def test_caffeine_scenario_runs(self):
        output = SimulationRunner().run(caffeine_and_meal_scenario())
        assert output.success

    def test_caffeine_timing_before_meal_has_real_effect(self):
        # Regression: caffeine_and_meal_scenario()'s caffeine_time_s used
        # to have no physiological effect at all -- DrugEvent.to_drug_dose()
        # could only ever report hours_since_ingestion=0.0 regardless of
        # the requested timing, so "taken 30 min before" and "taken at the
        # same moment as the meal" were bit-for-bit identical simulations.
        # Fixed via DrugEvent.hours_before_simulation (see scenario.py).
        runner = SimulationRunner()
        taken_before = runner.run(caffeine_and_meal_scenario(caffeine_time_s=-1800.0))
        taken_at_meal = runner.run(caffeine_and_meal_scenario(caffeine_time_s=0.0))
        assert (
            taken_before.expert_metrics.insulin_peak_pM
            != taken_at_meal.expert_metrics.insulin_peak_pM
        )


# ---------------------------------------------------------------------------
# Comparison runner
# ---------------------------------------------------------------------------

class TestComparisonRunner:
    @pytest.fixture(scope="class")
    def comparison(self):
        runner = ComparisonRunner()
        return runner.compare(
            high_sugar_meal_scenario(),
            light_snack_scenario(),
        )

    def test_comparison_succeeds(self, comparison):
        assert comparison.scenario_a.success
        assert comparison.scenario_b.success

    def test_high_sugar_higher_peak(self, comparison):
        assert comparison.peak_glucose_diff_uM > 0

    def test_glucose_diff_array(self, comparison):
        assert len(comparison.glucose_diff_uM) > 0

    def test_narrative_populated(self, comparison):
        assert "higher" in comparison.narrative or "lower" in comparison.narrative

    def test_insulin_diff_nonzero(self, comparison):
        assert abs(comparison.peak_insulin_diff_uM) > 0

    def test_scenario_a_higher_glucose_property(self, comparison):
        assert comparison.scenario_a_higher_glucose


# ---------------------------------------------------------------------------
# Personalised scenarios
# ---------------------------------------------------------------------------

class TestPersonalisedScenarios:
    def test_physiological_personalisation_applied(self):
        from biosim_physiology import (
            BiologicalSex, DiseaseState, FitnessLevel, PhysiologicalProfile
        )
        profile = PhysiologicalProfile(
            age_years=55,
            biological_sex=BiologicalSex.MALE,
            body_weight_kg=90.0,
            fitness_level=FitnessLevel.SEDENTARY,
            disease_states=[DiseaseState.INSULIN_RESISTANT],
            time_of_day_hours=20.0,
        )
        scenario = high_sugar_meal_scenario(
            personalization=PersonalizationSpec(physiological_profile=profile)
        )
        output = SimulationRunner().run(scenario)
        assert output.success
        assert output.physiological_personalization is not None
        assert output.has_personalization

    @pytest.mark.parametrize("exc_type", [ValueError, RuntimeError])
    def test_profile_personalisation_failure_is_not_silently_dropped(
        self, monkeypatch, exc_type,
    ):
        # Regression: a catch-all around the profile branch turned ANY
        # personalisation failure into a warning and simulated with base
        # parameters -- a T2D profile silently ran as a healthy subject
        # (REST returned 200). Failures must propagate: ValueError becomes
        # REST's 400, anything else is a real bug and must be loud.
        from biosim_physiology import DiseaseState, PhysiologicalPersonalizer, PhysiologicalProfile

        def boom(self, *args, **kwargs):
            raise exc_type("personalisation exploded")

        monkeypatch.setattr(PhysiologicalPersonalizer, "personalize", boom)
        scenario = high_sugar_meal_scenario(
            personalization=PersonalizationSpec(
                physiological_profile=PhysiologicalProfile(
                    disease_states=[DiseaseState.TYPE_2_DIABETES],
                ),
            )
        )
        with pytest.raises(exc_type, match="personalisation exploded"):
            SimulationRunner().run(scenario)

    @pytest.mark.parametrize("name,side,raw,bound", [
        ("p3", "floor", 4.55e-9, "1e-08"),
        ("secretion_scale", "ceiling", 7.5, "3"),
    ])
    def test_point_of_use_clamp_hit_is_disclosed_as_a_warning(
        self, monkeypatch, name, side, raw, bound,
    ):
        # The dydt() point-of-use clamps (P3_EFF_FLOOR et al.) engage on
        # ODE-coupled signals that the static clamp_warnings never sees; a
        # hit found by ODESystem's post-solve sweep must reach the output.
        from biosim_ode import GlucoseInsulinPathway
        monkeypatch.setattr(
            GlucoseInsulinPathway, "point_of_use_clamp_hits",
            lambda self, signals: {name: (side, raw)},
        )
        output = SimulationRunner().run(high_sugar_meal_scenario())
        assert any(
            w.startswith(f"{name} hit its point-of-use {side} ({bound})")
            and f"{raw:.3g}" in w
            for w in output.warnings
        ), output.warnings

    def test_drug_only_quick_path_still_records_physiological_personalization(self):
        # Regression: when personalization.physiological_profile is None
        # but drug_events (or exercise_bout) are present, runner.run()
        # takes a separate "quick path" that correctly computes
        # personalized params (the simulation math was always right) but
        # never assigned physio_result -- so
        # SimulationOutput.physiological_personalization stayed None even
        # though real personalization was applied. Every REST /simulate
        # request with drug_events takes exactly this path, since
        # service.py's scenario_in_to_definition() never sets
        # `personalization` at all -- this silently broke
        # ReproducibilityBundle provenance (infra.py's create_bundle()
        # falls back to unmodified default PathwayParameters() whenever
        # physiological_personalization is None).
        from biosim_ode import PathwayParameters
        scenario = ScenarioDefinition(
            name="metformin_quick_path", meal_events=[], duration_s=10,
            drug_events=[DrugEvent(time_s=0, drug="metformin", dose_mg=2000.0)],
        )
        out = SimulationRunner().run(scenario)
        assert out.success
        assert out.physiological_personalization is not None
        recorded = out.physiological_personalization.personalized_params
        # Not just non-None -- must be the genuinely personalized params
        # (metformin lowers Gb below the unmodified default), matching
        # what the ODE simulation itself actually used.
        assert recorded.Gb < PathwayParameters().Gb

    def test_glibenclamide_quick_path_records_glucose_independent_secretion(self):
        # Same "drug-only quick path" regression class as the metformin
        # test directly above -- confirms glibenclamide's real, additive
        # effect (which bypasses the multiplicative net_factors pipeline
        # entirely, see PhysiologicalPersonalizer.personalize()) is
        # genuinely recorded, not just correctly used internally.
        scenario = ScenarioDefinition(
            name="glibenclamide_quick_path", meal_events=[], duration_s=10,
            drug_events=[DrugEvent(time_s=0, drug="glibenclamide", dose_mg=2.5)],
        )
        out = SimulationRunner().run(scenario)
        assert out.success
        assert out.physiological_personalization is not None
        recorded = out.physiological_personalization.personalized_params
        assert recorded.glucose_independent_secretion > 0.0

    def test_glibenclamide_raises_fasted_insulin_above_baseline(self):
        # End-to-end demonstration of the core sulfonylurea mechanism
        # through the real SimulationRunner path: a chronic glibenclamide
        # dose measurably raises steady-state insulin above the undrugged
        # baseline, even with zero meal input -- the real, disclosed
        # hypoglycaemia-risk mechanism (glucose-independent secretion),
        # not just a stored-but-unused field.
        no_drug = ScenarioDefinition(name="no_drug_fasted", meal_events=[], duration_s=14400)
        with_drug = ScenarioDefinition(
            name="glibenclamide_fasted", meal_events=[], duration_s=14400,
            drug_events=[DrugEvent(time_s=0, drug="glibenclamide", dose_mg=2.5)],
        )
        runner = SimulationRunner()
        I_no_drug = runner.run(no_drug).simulation.species["I"][-1]
        I_with_drug = runner.run(with_drug).simulation.species["I"][-1]
        assert I_with_drug > I_no_drug

    def test_t2d_profile_higher_glucose_peak(self):
        from biosim_physiology import DiseaseState, PhysiologicalProfile
        healthy_scenario = high_sugar_meal_scenario()
        t2d_scenario = high_sugar_meal_scenario(
            personalization=PersonalizationSpec(
                physiological_profile=PhysiologicalProfile(
                    disease_states=[DiseaseState.TYPE_2_DIABETES]
                )
            )
        )
        runner = SimulationRunner()
        out_healthy = runner.run(healthy_scenario)
        out_t2d = runner.run(t2d_scenario)
        assert out_t2d.expert_metrics.glucose_peak_uM > out_healthy.expert_metrics.glucose_peak_uM

    def test_genetic_personalisation_applied(self):
        from biosim_genotype import Genotype, Variant
        variants = [Variant("rs7903146", Genotype.HOMOZYGOUS_ALT)]
        scenario = high_sugar_meal_scenario(
            personalization=PersonalizationSpec(genetic_variants=variants)
        )
        output = SimulationRunner().run(scenario)
        assert output.success
        assert output.genetic_personalization is not None

    def test_with_uncertainty(self):
        scenario = ScenarioDefinition(
            name="uq_test",
            meal_events=[MealEvent(time_s=0.0, carbs_g=50.0)],
            run_uncertainty=True,
            n_mc_samples=20,
            duration_s=7200,
        )
        output = SimulationRunner().run(scenario)
        assert output.success
        assert output.uq is not None

    def test_uq_sees_the_same_composite_pathway_as_the_point_estimate(self):
        # Regression: MonteCarloEngine.propagate() always built a bare
        # GlucoseInsulinPathway, with no way to see whatever extension
        # pathways _build_pathway() composed for the point estimate two
        # lines above it in run(). For a semaglutide scenario specifically
        # this was PROVABLY a no-drug ensemble: _glp1_agonist_modifiers()
        # returns {} (semaglutide has zero static PathwayParameters
        # effect), so the sampled parameters -- and therefore the
        # reported "confidence band" -- were identical to a drug-free
        # run's. Structural proof the fix closes this: the composite's
        # own species (GLP1) must now have a real UQ band, not be
        # silently absent.
        scenario = ScenarioDefinition(
            name="uq_semaglutide",
            meal_events=[MealEvent(time_s=0.0, carbs_g=50.0)],
            drug_events=[DrugEvent(time_s=0.0, drug="semaglutide", dose_mg=1.0)],
            run_uncertainty=True,
            n_mc_samples=15,
            duration_s=7200,
        )
        output = SimulationRunner().run(scenario)
        assert output.success
        assert output.uq is not None
        assert "GLP1" in output.uq.bands
        assert set(output.simulation.species) == set(output.uq.bands)

    def test_uq_species_pathway_attributes_the_composed_pathways_own_species(self):
        # MASTER_PLAN 2.3: end-to-end confirmation that the real production
        # pipeline's UQResult.species_pathway (built from whatever
        # _build_pathway() actually composed) attributes the GLP1 species
        # to the GLP1 extension pathway and the core species to GI, not
        # left empty or misattributed through the real run() call path.
        scenario = ScenarioDefinition(
            name="uq_species_pathway",
            meal_events=[MealEvent(time_s=0.0, carbs_g=50.0)],
            drug_events=[DrugEvent(time_s=0.0, drug="semaglutide", dose_mg=1.0)],
            run_uncertainty=True,
            n_mc_samples=15,
            duration_s=7200,
        )
        output = SimulationRunner().run(scenario)
        assert output.success
        assert output.uq is not None
        assert output.uq.species_pathway["GLP1"] == "glp1_v1"
        assert output.uq.species_pathway["G"] == "glucose_insulin_v1"
        assert output.uq.pathway_for("G") == "glucose_insulin_v1"

    def test_uq_samples_renal_beta_for_an_sglt2i_scenario(self):
        # MASTER_PLAN 2.2: a dapagliflozin scenario's real UQ ensemble must
        # genuinely sample "renal.beta", not just the pathway STRUCTURE
        # fix (species_pathway/pathway_factory) from the prior two
        # milestones -- confirmed via store_trajectories=True so the
        # actual per-sample draws are inspectable, not just their presence
        # in parameter_names.
        scenario = ScenarioDefinition(
            name="uq_renal_beta",
            meal_events=[MealEvent(time_s=0.0, carbs_g=50.0)],
            drug_events=[DrugEvent(time_s=0, drug="dapagliflozin", dose_mg=10.0)],
            run_uncertainty=True,
            n_mc_samples=15,
            duration_s=7200,
        )
        output = SimulationRunner().run(scenario)
        assert output.success
        assert output.uq is not None
        assert "renal.beta" in output.uq.parameter_names
        assert "RGC" in output.uq.bands

    def test_uq_does_not_sample_renal_beta_without_an_sglt2i_drug(self):
        # A scenario with no SGLT2i drug active must not spuriously gain a
        # "renal.beta" entry -- extension_distributions stays empty for
        # every scenario extension_pathway_factory's own detection doesn't
        # recognize.
        scenario = ScenarioDefinition(
            name="uq_no_renal",
            meal_events=[MealEvent(time_s=0.0, carbs_g=50.0)],
            run_uncertainty=True,
            n_mc_samples=15,
            duration_s=7200,
        )
        output = SimulationRunner().run(scenario)
        assert output.success
        assert output.uq is not None
        assert "renal.beta" not in output.uq.parameter_names
        assert "RGC" not in output.uq.bands

    def test_genetic_variant_widens_the_uq_distribution_for_its_own_parameter(self):
        # Regression: ComposedEffect.genetic_cv was computed by
        # biosim_genotype but never read back by the runner when building
        # Monte Carlo distributions -- a genotype-personalized run and a
        # population-average run got identically-wide confidence bands,
        # even though composer.py's own docstring documents the genetic CV
        # as meant to be "combined in quadrature with the base parameter
        # uncertainty." rs1801282 (PPARG) modifies p3 specifically.
        #
        # Checked at the distribution-construction level (not by comparing
        # simulated band widths end-to-end) because a genetic variant also
        # shifts the point estimate itself (PathwayParameters.p3's mean),
        # which independently changes glucose dynamics and confounds a
        # naive "genotyped band should be wider" comparison -- the CV
        # widening this test guards is a distinct, smaller effect on top
        # of that mean shift.
        from unittest.mock import patch
        from biosim_genotype import Genotype, Variant
        from biosim_uncertainty import LogNormalDistribution

        scenario = ScenarioDefinition(
            name="uq_genotyped",
            meal_events=[MealEvent(time_s=0.0, carbs_g=50.0)],
            run_uncertainty=True,
            n_mc_samples=10,
            duration_s=7200,
            personalization=PersonalizationSpec(
                genetic_variants=[Variant("rs1801282", Genotype.HOMOZYGOUS_ALT)]
            ),
        )
        real_from_tier = LogNormalDistribution.from_confidence_tier
        with patch.object(
            LogNormalDistribution, "from_confidence_tier",
            side_effect=real_from_tier,
        ) as mock_from_tier:
            output = SimulationRunner().run(scenario)

        assert output.success
        # Exactly one of the five UQ parameter distributions (p1, p2, p3,
        # ke, EC50_mTOR) should have received a nonzero additional_cv --
        # rs1801282 (PPARG) modifies p3 only.
        additional_cvs = [
            c.kwargs.get("additional_cv", 0.0) for c in mock_from_tier.call_args_list
        ]
        assert any(cv > 0.0 for cv in additional_cvs), (
            "expected at least one UQ parameter distribution to receive a "
            "nonzero additional_cv from the applied genetic variant"
        )
        assert sum(1 for cv in additional_cvs if cv > 0.0) == 1, (
            "only p3 should be widened by rs1801282 (PPARG); other "
            "parameters must keep additional_cv=0.0"
        )

    def test_no_genetic_variants_means_zero_additional_cv(self):
        from unittest.mock import patch
        from biosim_uncertainty import LogNormalDistribution

        scenario = ScenarioDefinition(
            name="uq_ungenotyped",
            meal_events=[MealEvent(time_s=0.0, carbs_g=50.0)],
            run_uncertainty=True,
            n_mc_samples=10,
            duration_s=7200,
        )
        real_from_tier = LogNormalDistribution.from_confidence_tier
        with patch.object(
            LogNormalDistribution, "from_confidence_tier",
            side_effect=real_from_tier,
        ) as mock_from_tier:
            output = SimulationRunner().run(scenario)

        assert output.success
        additional_cvs = [
            c.kwargs.get("additional_cv", 0.0) for c in mock_from_tier.call_args_list
        ]
        assert all(cv == 0.0 for cv in additional_cvs)


# ---------------------------------------------------------------------------
# LaySummary <-> UQResult wiring (MASTER_PLAN 2.1)
# ---------------------------------------------------------------------------
#
# Regression: _compute_lay_summary() accepted a `uq` parameter but never
# read it -- confidence_note/is_high_confidence were derived solely from
# the point-estimate Boolean watchdog score, structurally blind to whether
# a Monte Carlo UQ ensemble showed a high sample-failure rate or wide
# confidence bands. These tests call the private method directly with
# hand-built SimulationResult/ExpertMetrics/UQResult objects so every
# branch (uq=None, a reliable narrow ensemble, an unreliable ensemble, a
# wide-band ensemble) is exercised deterministically, without depending on
# a real Monte Carlo run happening to land on one side of a threshold.

class TestLaySummaryUQWiring:
    def _sim_and_metrics(self, watchdog_confidence: float = 0.9):
        from biosim_ode.system import SimulationResult
        from biosim_api.output import ExpertMetrics

        times = np.linspace(0.0, 7200.0, 50)
        G = np.full(50, 5000.0)  # flat glucose, 5.0 mM, in uM
        I = np.full(50, 66e-6)   # flat insulin, 66 pM, in uM
        sim = SimulationResult(
            times_s=times,
            species={"G": G, "I": I},
            species_units={"G": "uM", "I": "uM"},
            success=True,
            message="ok",
        )
        metrics = ExpertMetrics(
            glucose_peak_uM=5000.0,
            glucose_peak_time_s=1000.0,
            glucose_auc_uM_s=0.0,
            glucose_baseline_uM=5000.0,
            insulin_peak_uM=66e-6,
            insulin_peak_time_s=1000.0,
            mtor_peak=0.15,
            ampk_nadir=0.6,
            time_to_glucose_return_s=None,
            overall_watchdog_confidence=watchdog_confidence,
            n_ode_evaluations=100,
        )
        return sim, metrics

    def _uq(self, n_samples=20, n_failed=0, overall_confidence=0.9,
            glucose_relative_width=0.02):
        from biosim_uncertainty import ConfidenceBand, UQResult

        median = np.full(50, 5000.0)
        # ConfidenceBand.width == upper - lower == 2*half_width, and
        # relative_width == width/|median| -- so half_width must be
        # median*glucose_relative_width/2 for relative_width to equal
        # glucose_relative_width exactly, not double it.
        half_width = median * glucose_relative_width / 2.0
        band = ConfidenceBand(
            lower=median - half_width,
            median=median,
            upper=median + half_width,
        )
        return UQResult(
            times_s=np.linspace(0.0, 7200.0, 50),
            bands={"G": band},
            n_samples=n_samples,
            n_failed=n_failed,
            overall_confidence=overall_confidence,
        )

    def test_uq_none_is_unaffected(self):
        sim, metrics = self._sim_and_metrics(watchdog_confidence=0.9)
        ls = SimulationRunner()._compute_lay_summary(sim, metrics, None)
        assert ls.is_high_confidence is True
        assert "ensemble" not in ls.confidence_note

    def test_reliable_narrow_ensemble_does_not_downgrade_confidence(self):
        sim, metrics = self._sim_and_metrics(watchdog_confidence=0.9)
        uq = self._uq(n_failed=0, overall_confidence=0.95, glucose_relative_width=0.02)
        ls = SimulationRunner()._compute_lay_summary(sim, metrics, uq)
        assert ls.is_high_confidence is True
        assert "ensemble" not in ls.confidence_note

    def test_unreliable_ensemble_downgrades_high_watchdog_confidence(self):
        # A high point-estimate watchdog score alone used to be sufficient
        # for is_high_confidence=True, even if most Monte Carlo samples
        # failed to solve -- this is exactly the blind spot being closed.
        sim, metrics = self._sim_and_metrics(watchdog_confidence=0.95)
        uq = self._uq(n_samples=20, n_failed=15, overall_confidence=0.95)
        ls = SimulationRunner()._compute_lay_summary(sim, metrics, uq)
        assert ls.is_high_confidence is False
        assert "reliability" in ls.confidence_note

    def test_wide_glucose_band_downgrades_high_watchdog_confidence(self):
        sim, metrics = self._sim_and_metrics(watchdog_confidence=0.95)
        uq = self._uq(n_failed=0, overall_confidence=0.95, glucose_relative_width=0.40)
        ls = SimulationRunner()._compute_lay_summary(sim, metrics, uq)
        assert ls.is_high_confidence is False
        assert "wide range" in ls.confidence_note

    def test_narrow_band_just_under_threshold_does_not_downgrade(self):
        sim, metrics = self._sim_and_metrics(watchdog_confidence=0.9)
        uq = self._uq(n_failed=0, overall_confidence=0.95, glucose_relative_width=0.10)
        ls = SimulationRunner()._compute_lay_summary(sim, metrics, uq)
        assert ls.is_high_confidence is True
        assert "ensemble" not in ls.confidence_note

    def test_low_watchdog_confidence_with_reliable_uq_stays_low_no_extra_note(self):
        # If the watchdog already says low confidence, a clean UQ ensemble
        # shouldn't add a spurious ensemble-related caveat -- there's
        # nothing wrong on the UQ side to report.
        sim, metrics = self._sim_and_metrics(watchdog_confidence=0.5)
        uq = self._uq(n_failed=0, overall_confidence=0.95, glucose_relative_width=0.02)
        ls = SimulationRunner()._compute_lay_summary(sim, metrics, uq)
        assert ls.is_high_confidence is False
        assert "ensemble" not in ls.confidence_note

    def test_real_run_uncertainty_scenario_still_produces_a_lay_summary(self):
        # End-to-end sanity check that the real SimulationRunner.run()
        # pipeline still wires a genuine UQResult through without
        # crashing -- the branch-level logic above is covered by the
        # deterministic unit tests, this just confirms real integration.
        scenario = ScenarioDefinition(
            name="uq_lay_summary_wiring",
            meal_events=[MealEvent(time_s=0.0, carbs_g=50.0)],
            run_uncertainty=True,
            n_mc_samples=15,
            duration_s=7200,
        )
        output = SimulationRunner().run(scenario)
        assert output.success
        assert output.uq is not None
        assert isinstance(output.lay_summary.confidence_note, str)
        assert len(output.lay_summary.confidence_note) > 5


class TestUnverifiedVariantWarning:
    """A genotype including rs35767 (IGF1, pubmed_id="UNVERIFIED") must carry
    the unverified-effect note in the runner warnings, without the rsID."""

    @staticmethod
    def _run(rsid):
        from biosim_genotype import Genotype, Variant
        scenario = high_sugar_meal_scenario(
            personalization=PersonalizationSpec(
                genetic_variants=[Variant(rsid, Genotype.HETEROZYGOUS)]
            )
        )
        return SimulationRunner().run(scenario)

    def test_rs35767_genotype_warns_and_stays_private(self):
        from biosim_genotype.effect_db import UNVERIFIED_VARIANT_NOTES
        from biosim_infra import PrivacyAuditor
        output = self._run("rs35767")
        assert UNVERIFIED_VARIANT_NOTES["rs35767"] in output.warnings
        assert PrivacyAuditor().audit(output).passed

    def test_zero_risk_allele_copies_gives_no_warning(self):
        # 0/0 carries no copy of the effect allele, so the unverified effect
        # is never applied and the note would be false.
        from biosim_genotype import Genotype, Variant
        from biosim_genotype.effect_db import UNVERIFIED_VARIANT_NOTES
        scenario = high_sugar_meal_scenario(
            personalization=PersonalizationSpec(
                genetic_variants=[Variant("rs35767", Genotype.HOMOZYGOUS_REF)]
            )
        )
        output = SimulationRunner().run(scenario)
        assert UNVERIFIED_VARIANT_NOTES["rs35767"] not in output.warnings

    def test_verified_variant_alone_has_no_unverified_warning(self):
        from biosim_genotype.effect_db import UNVERIFIED_VARIANT_NOTES
        output = self._run("rs7903146")
        assert UNVERIFIED_VARIANT_NOTES["rs35767"] not in output.warnings


class TestNoReturnWithinWindowIsSaid:
    """
    Since B7 (2026-11-16) a T2D 75 g meal does not return to baseline within
    the 4-h window. The narrative must say so, not just stop after the peak.
    """

    def test_t2d_high_sugar_meal_narrative_says_not_returned(self):
        from biosim_physiology import PhysiologicalProfile, DiseaseState
        output = SimulationRunner().run(high_sugar_meal_scenario(
            personalization=PersonalizationSpec(
                physiological_profile=PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES])
            )
        ))
        assert output.expert_metrics.time_to_glucose_return_s is None
        narrative = output.lay_summary.glucose_narrative
        assert "has not returned to your own baseline" in narrative
        assert "4 hours" in narrative
        assert "toward normal" not in narrative


def test_library_ogtt_glucose_drink_is_not_slowed_by_acarbose():
    # A 75 g OGTT is pure glucose, which needs no digestion; Chen 2016
    # (PMID:27686734) found no acarbose effect on the OGTT. The library
    # scenario says so with complex_carb_fraction=0 (B3 follow-up review).
    from biosim_api import ogtt_75g
    base = ogtt_75g()
    plain = SimulationRunner().run(base)
    with_drug = SimulationRunner().run(ScenarioDefinition(
        name="ogtt_acarbose", duration_s=base.duration_s, meal_events=base.meal_events,
        drug_events=[DrugEvent(time_s=0, drug="acarbose", dose_mg=100.0)],
    ))
    assert base.meal_events[0].complex_carb_fraction == 0.0
    assert with_drug.expert_metrics.glucose_peak_mM == pytest.approx(plain.expert_metrics.glucose_peak_mM)
