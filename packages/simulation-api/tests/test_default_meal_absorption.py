# SAVE TO: packages\simulation-api\tests\test_default_meal_absorption.py
"""
Scenario meals with no glycaemic_index and no kabs use the two-compartment
("Bateman") absorption curve: zero at the moment of eating, a later peak,
and the full meal's glucose delivered in total.
"""
import numpy as np
import pytest
from scipy.integrate import quad

from biosim_api import (
    DrugEvent,
    MealEvent,
    ScenarioDefinition,
    SimulationRunner,
    high_sugar_meal_scenario,
)
from biosim_ode.gut_absorption import GutAbsorptionModel


class TestDefaultMealEventAbsorption:
    def test_default_meal_uses_two_compartment_model(self):
        meal = MealEvent(time_s=0.0, carbs_g=75.0).to_meal_input()
        assert meal.absorption_model == GutAbsorptionModel.default_meal()

    def test_default_meal_rate_is_zero_when_eaten(self):
        meal = MealEvent(time_s=600.0, carbs_g=75.0).to_meal_input()
        assert meal.Ra(meal.time_min) == 0.0

    def test_default_meal_rate_peaks_later(self):
        meal = MealEvent(time_s=0.0, carbs_g=75.0).to_meal_input()
        times = np.linspace(0.0, 300.0, 3001)
        rates = [meal.Ra(t) for t in times]
        assert 10.0 < times[int(np.argmax(rates))] < 120.0

    def test_default_meal_delivers_its_glucose(self):
        meal = MealEvent(time_s=0.0, carbs_g=75.0).to_meal_input()
        total, _ = quad(meal.Ra, 0.0, 2000.0, limit=200)
        assert total == pytest.approx(meal.D_eff_mM, rel=0.01)

    def test_library_high_sugar_meal_is_two_compartment(self):
        scenario = high_sugar_meal_scenario()
        meal = scenario.meal_events[0].to_meal_input()
        assert meal.Ra(0.0) == 0.0

    def test_explicit_kabs_keeps_single_exponential(self):
        # A caller who sets kabs asked for the single-exponential curve.
        meal = MealEvent(time_s=0.0, carbs_g=75.0, kabs=0.03).to_meal_input()
        assert meal.absorption_model is None
        assert meal.Ra(0.0) > 0.0


class TestDefaultMealWithDrugs:
    def test_glp1_drug_with_default_meal_runs(self):
        scenario = ScenarioDefinition(
            name="sema_default_meal",
            meal_events=[MealEvent(time_s=0, carbs_g=75.0)],
            drug_events=[DrugEvent(time_s=0, drug="semaglutide", dose_mg=1.0)],
            duration_s=14400,
        )
        out = SimulationRunner().run(scenario)
        assert out.expert_metrics.glucose_peak_mM > 5.0

    def test_glp1_drug_with_gi_meal_is_now_allowed(self):
        # Two-compartment meals are now state-tracked (mass-conserving), so
        # the old "glycaemic_index + GLP-1 drug" rejection is no longer needed.
        scenario = ScenarioDefinition(
            name="gi_dpp4",
            meal_events=[MealEvent(time_s=0, carbs_g=75.0, glycaemic_index=100.0)],
            drug_events=[DrugEvent(time_s=0, drug="dpp4_inhibitor", dose_mg=100.0)],
            duration_s=14400,
        )
        out = SimulationRunner().run(scenario)
        assert out.expert_metrics.glucose_peak_mM > 5.0

    def test_acarbose_still_slows_a_default_meal_and_says_how(self):
        plain = ScenarioDefinition(
            name="plain", meal_events=[MealEvent(time_s=0, carbs_g=75.0)], duration_s=14400,
        )
        with_acarbose = ScenarioDefinition(
            name="acarbose", meal_events=[MealEvent(time_s=0, carbs_g=75.0)], duration_s=14400,
            drug_events=[DrugEvent(time_s=0, drug="acarbose", dose_mg=100.0)],
        )
        runner = SimulationRunner()
        out_plain = runner.run(plain)
        out_acarbose = runner.run(with_acarbose)
        assert out_acarbose.expert_metrics.glucose_peak_mM < out_plain.expert_metrics.glucose_peak_mM
        assert out_acarbose.expert_metrics.glucose_peak_time_s > out_plain.expert_metrics.glucose_peak_time_s
        assert any("two-compartment" in w for w in out_acarbose.warnings)
