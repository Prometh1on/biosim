# SAVE TO: packages\data-seeds\tests\test_drug_scenarios.py
"""
Tests for drug interaction scenarios.
Each test verifies the expected direction of effect from the mechanism.
"""
import pytest
from biosim_seeds.drug_scenarios import (
    alcohol_before_meal, caffeine_before_meal, dpp4_inhibitor_scenario,
    liraglutide_scenario, metformin_with_meal, semaglutide_scenario,
    steroid_meal, morning_coffee_breakfast,
)
from biosim_api import SimulationRunner, high_sugar_meal_scenario, MealEvent, ScenarioDefinition


@pytest.fixture(scope="module")
def baseline():
    return SimulationRunner().run(high_sugar_meal_scenario())


class TestCaffeineScenario:
    def test_runs_successfully(self):
        out = SimulationRunner().run(caffeine_before_meal())
        assert out.success

    def test_has_drug_event(self):
        s = caffeine_before_meal()
        assert len(s.drug_events) == 1
        assert s.drug_events[0].drug == "caffeine"

    def test_caffeine_raises_glucose_vs_no_caffeine(self, baseline):
        """Caffeine reduces insulin sensitivity → higher glucose peak."""
        out_caff = SimulationRunner().run(
            caffeine_before_meal(caffeine_mg=300.0, carbs_g=75.0)
        )
        # With caffeine: glucose peak should be >= baseline (same carbs, less insulin sensitivity)
        # The effect is modest so we just check it runs and produces a plausible result
        assert out_caff.expert_metrics.glucose_peak_mM > 5.0

    def test_caffeine_actually_differs_from_true_no_drug_baseline(self):
        # Regression: hours_before_simulation was never set on this
        # scenario's DrugEvent (defaulted to 0.0), so
        # DrugEvent.to_drug_dose()'s hours_since_ingestion was always
        # exactly 0 -- caffeine_modifiers()'s concentration_fraction is
        # `0.0 if t <= 0 else ...`, making the drug's effect exactly zero
        # regardless of dose. This scenario's output was bit-for-bit
        # identical to the same meal with the drug event stripped
        # entirely, despite the docstring's explicit "45 minutes before a
        # meal" claim. test_caffeine_raises_glucose_vs_no_caffeine above
        # didn't catch this -- it only asserted `> 5.0`, true even at
        # zero effect.
        no_drug = ScenarioDefinition(
            name="nd", duration_s=14400, meal_events=[MealEvent(0, 60.0)],
        )
        with_caffeine = caffeine_before_meal(caffeine_mg=200.0, carbs_g=60.0)
        peak_nd = SimulationRunner().run(no_drug).expert_metrics.glucose_peak_mM
        peak_caf = SimulationRunner().run(with_caffeine).expert_metrics.glucose_peak_mM
        assert peak_caf > peak_nd

    def test_higher_caffeine_dose_reduces_p3(self):
        from biosim_physiology.drug_interactions import Drug, DrugDose, drug_modifiers
        low  = drug_modifiers([DrugDose(Drug.CAFFEINE, 100.0, hours_since_ingestion=1.0)])
        high = drug_modifiers([DrugDose(Drug.CAFFEINE, 400.0, hours_since_ingestion=1.0)])
        assert high.get("p3", 1.0) <= low.get("p3", 1.0)


class TestAlcoholBeforeMealScenario:
    """
    Regression: hours_before was encoded into DrugEvent.time_s
    (max(0, hours_before*3600-3600)), but DrugEvent.to_drug_dose()
    computes hours_since_ingestion from hours_before_simulation, not
    time_s -- for any simulation_start_s=0.0 call (every real caller),
    time_s's contribution to hours_since_ingestion is always exactly 0
    regardless of its value. hours_before had no effect whatsoever:
    sweeping it produced bit-for-bit identical output every time. Zero
    test coverage existed for this function before this fix.
    """
    def test_runs_successfully(self):
        out = SimulationRunner().run(alcohol_before_meal())
        assert out.success

    def test_has_drug_event(self):
        s = alcohol_before_meal()
        assert len(s.drug_events) == 1
        assert s.drug_events[0].drug == "alcohol"

    def test_hours_before_actually_changes_the_effect(self):
        peak_1h = SimulationRunner().run(
            alcohol_before_meal(hours_before=1.0)
        ).expert_metrics.glucose_peak_mM
        peak_5h = SimulationRunner().run(
            alcohol_before_meal(hours_before=5.0)
        ).expert_metrics.glucose_peak_mM
        assert peak_1h != pytest.approx(peak_5h)

    def test_alcohol_effect_resolves_as_time_passes(self):
        # AlcoholHepaticPathway resolves by 6h post-ingestion (corrected
        # 2026-08-18 -- previously attributed to the now-dead
        # alcohol_modifiers(), which returns {} unconditionally; the real
        # mechanism is the ODE-coupled AlcoholHepaticPathway, sharing the
        # same 360-minute resolution window) -- longer hours_before
        # should mean less residual effect at meal time (t=0), so glucose
        # peak should move back toward the true no-drug baseline.
        no_drug = ScenarioDefinition(
            name="nd", duration_s=14400, meal_events=[MealEvent(0, 60.0)],
        )
        peak_nd = SimulationRunner().run(no_drug).expert_metrics.glucose_peak_mM
        peak_1h = SimulationRunner().run(
            alcohol_before_meal(hours_before=1.0)
        ).expert_metrics.glucose_peak_mM
        peak_5h = SimulationRunner().run(
            alcohol_before_meal(hours_before=5.0)
        ).expert_metrics.glucose_peak_mM
        assert abs(peak_5h - peak_nd) < abs(peak_1h - peak_nd)


class TestMetforminScenario:
    def test_runs_successfully(self):
        out = SimulationRunner().run(metformin_with_meal())
        assert out.success

    def test_has_drug_event(self):
        s = metformin_with_meal()
        assert s.drug_events[0].drug == "metformin"

    def test_metformin_lowers_glucose_vs_no_drug(self):
        """Metformin should reduce peak glucose."""
        from biosim_api import MealEvent, ScenarioDefinition
        no_drug = ScenarioDefinition(
            name="no_drug", duration_s=14400,
            meal_events=[MealEvent(0, 75.0)],
        )
        with_metformin = metformin_with_meal(carbs_g=75.0, metformin_mg=1000.0)
        runner = SimulationRunner()
        out_nd = runner.run(no_drug)
        out_mf = runner.run(with_metformin)
        assert out_mf.expert_metrics.glucose_peak_mM < out_nd.expert_metrics.glucose_peak_mM


class TestSteroidScenario:
    def test_runs_successfully(self):
        out = SimulationRunner().run(steroid_meal())
        assert out.success

    def test_has_drug_event(self):
        s = steroid_meal()
        assert s.drug_events[0].drug == "glucocorticoid"

    def test_steroid_raises_glucose_vs_no_drug(self):
        """Steroid reduces insulin sensitivity → higher glucose peak."""
        from biosim_api import MealEvent, ScenarioDefinition
        no_drug = ScenarioDefinition(
            name="nd", duration_s=14400,
            meal_events=[MealEvent(0, 60.0)],
        )
        with_steroid = steroid_meal(carbs_g=60.0, prednisone_mg=20.0)
        runner = SimulationRunner()
        peak_nd = runner.run(no_drug).expert_metrics.glucose_peak_mM
        peak_st = runner.run(with_steroid).expert_metrics.glucose_peak_mM
        assert peak_st > peak_nd


class TestSemaglutideScenario:
    def test_runs_successfully(self):
        out = SimulationRunner().run(semaglutide_scenario())
        assert out.success

    def test_has_drug_event(self):
        s = semaglutide_scenario()
        assert s.drug_events[0].drug == "semaglutide"

    def test_semaglutide_reduces_glucose_peak_by_more_than_1mM(self):
        """CLAUDE.md Priority 1 acceptance criterion, exercised at the
        scenario/SimulationRunner level (not just the low-level ODE coupling)."""
        no_drug = ScenarioDefinition(
            name="nd", duration_s=14400, meal_events=[MealEvent(0, 75.0)],
        )
        runner = SimulationRunner()
        peak_nd = runner.run(no_drug).expert_metrics.glucose_peak_mM
        peak_sema = runner.run(semaglutide_scenario(carbs_g=75.0, dose_mg=1.0)).expert_metrics.glucose_peak_mM
        assert peak_nd - peak_sema > 1.0

    def test_higher_dose_larger_effect(self):
        runner = SimulationRunner()
        peak_low  = runner.run(semaglutide_scenario(dose_mg=0.5)).expert_metrics.glucose_peak_mM
        peak_high = runner.run(semaglutide_scenario(dose_mg=2.0)).expert_metrics.glucose_peak_mM
        assert peak_high < peak_low


class TestLiraglutideScenario:
    def test_runs_successfully(self):
        out = SimulationRunner().run(liraglutide_scenario())
        assert out.success

    def test_has_drug_event(self):
        s = liraglutide_scenario()
        assert s.drug_events[0].drug == "liraglutide"

    def test_liraglutide_lowers_glucose_vs_no_drug(self):
        no_drug = ScenarioDefinition(
            name="nd", duration_s=14400, meal_events=[MealEvent(0, 75.0)],
        )
        runner = SimulationRunner()
        peak_nd = runner.run(no_drug).expert_metrics.glucose_peak_mM
        peak_lira = runner.run(liraglutide_scenario()).expert_metrics.glucose_peak_mM
        assert peak_lira < peak_nd


class TestDpp4InhibitorScenario:
    def test_runs_successfully(self):
        out = SimulationRunner().run(dpp4_inhibitor_scenario())
        assert out.success

    def test_has_drug_event(self):
        s = dpp4_inhibitor_scenario()
        assert s.drug_events[0].drug == "dpp4_inhibitor"

    def test_dpp4_inhibitor_lowers_glucose_vs_no_drug(self):
        no_drug = ScenarioDefinition(
            name="nd", duration_s=14400, meal_events=[MealEvent(0, 75.0)],
        )
        runner = SimulationRunner()
        peak_nd = runner.run(no_drug).expert_metrics.glucose_peak_mM
        peak_dpp4 = runner.run(dpp4_inhibitor_scenario()).expert_metrics.glucose_peak_mM
        assert peak_dpp4 < peak_nd

    def test_smaller_effect_than_semaglutide(self):
        """DPP-4 inhibition of endogenous GLP-1 should be a milder effect
        than pharmacological receptor agonism."""
        runner = SimulationRunner()
        peak_dpp4 = runner.run(dpp4_inhibitor_scenario()).expert_metrics.glucose_peak_mM
        peak_sema = runner.run(semaglutide_scenario()).expert_metrics.glucose_peak_mM
        assert peak_dpp4 > peak_sema


class TestMorningCoffeeBreakfast:
    def test_runs_successfully(self):
        out = SimulationRunner().run(morning_coffee_breakfast())
        assert out.success

    def test_has_both_meal_and_drug(self):
        s = morning_coffee_breakfast()
        assert len(s.meal_events) == 1
        assert len(s.drug_events) == 1
        assert s.drug_events[0].dose_mg == 200.0

    def test_glucose_rises_after_meal(self):
        out = SimulationRunner().run(morning_coffee_breakfast())
        assert out.expert_metrics.glucose_peak_mM > 5.5

    def test_differs_from_true_no_drug_baseline(self):
        # Regression: this scenario's DrugEvent never set
        # hours_before_simulation (defaulted to 0.0), so
        # DrugEvent.to_drug_dose()'s hours_since_ingestion was always
        # exactly 0 -- caffeine_modifiers()'s concentration_fraction is
        # `0.0 if t <= 0 else ...`, making caffeine's effect exactly zero
        # regardless of dose. test_glucose_rises_after_meal above didn't
        # catch this (true for the meal alone, with or without caffeine).
        no_drug = ScenarioDefinition(
            name="nd", duration_s=14400, meal_events=[MealEvent(0, 44.0, label="Toast + OJ")],
        )
        peak_nd = SimulationRunner().run(no_drug).expert_metrics.glucose_peak_mM
        peak_mc = SimulationRunner().run(morning_coffee_breakfast()).expert_metrics.glucose_peak_mM
        assert peak_mc > peak_nd
