# SAVE TO: packages\data-seeds\tests\test_condition_scenarios.py
"""
Tests for condition-based scenarios (stress, type 2 diabetes) that exercise
the cortisol/hepatic ODE-coupled extension pathways via SimulationRunner.
"""
import pytest
from biosim_seeds.condition_scenarios import stress_meal_scenario, type2_diabetes_meal_scenario
from biosim_api import SimulationRunner, MealEvent, ScenarioDefinition


class TestStressMealScenario:
    def test_runs_successfully(self):
        out = SimulationRunner().run(stress_meal_scenario())
        assert out.success

    def test_has_stress_level_set(self):
        s = stress_meal_scenario(stress_level=0.7)
        assert s.stress_level == pytest.approx(0.7)

    def test_stress_raises_glucose_vs_no_stress(self):
        no_stress = ScenarioDefinition(
            name="nd", duration_s=14400, meal_events=[MealEvent(0, 75.0)],
        )
        runner = SimulationRunner()
        peak_ns = runner.run(no_stress).expert_metrics.glucose_peak_mM
        peak_s = runner.run(stress_meal_scenario(carbs_g=75.0, stress_level=0.7)).expert_metrics.glucose_peak_mM
        assert peak_s > peak_ns

    def test_higher_stress_larger_effect(self):
        runner = SimulationRunner()
        peak_low = runner.run(stress_meal_scenario(stress_level=0.2)).expert_metrics.glucose_peak_mM
        peak_high = runner.run(stress_meal_scenario(stress_level=0.9)).expert_metrics.glucose_peak_mM
        assert peak_high > peak_low


class TestType2DiabetesMealScenario:
    def test_runs_successfully(self):
        out = SimulationRunner().run(type2_diabetes_meal_scenario())
        assert out.success

    def test_has_t2d_disease_state(self):
        from biosim_physiology import DiseaseState
        s = type2_diabetes_meal_scenario()
        assert DiseaseState.TYPE_2_DIABETES in s.personalization.physiological_profile.disease_states

    def test_t2d_raises_glucose_vs_healthy(self):
        no_drug = ScenarioDefinition(
            name="nd", duration_s=14400, meal_events=[MealEvent(0, 75.0)],
        )
        runner = SimulationRunner()
        peak_healthy = runner.run(no_drug).expert_metrics.glucose_peak_mM
        peak_t2d = runner.run(type2_diabetes_meal_scenario(carbs_g=75.0)).expert_metrics.glucose_peak_mM
        assert peak_t2d > peak_healthy

    def test_does_not_duplicate_healthy_state(self):
        from biosim_physiology import DiseaseState
        s = type2_diabetes_meal_scenario()
        assert DiseaseState.HEALTHY not in s.personalization.physiological_profile.disease_states
