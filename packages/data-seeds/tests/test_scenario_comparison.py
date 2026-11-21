# SAVE TO: packages\data-seeds\tests\test_scenario_comparison.py
"""Tests for the lay scenario comparison report."""
import pytest
from biosim_seeds.scenario_comparison_report import build_lay_comparison, LayComparisonReport
from biosim_api import SimulationRunner, high_sugar_meal_scenario, light_snack_scenario


@pytest.fixture(scope="module")
def outputs():
    runner = SimulationRunner()
    return runner.run(high_sugar_meal_scenario()), runner.run(light_snack_scenario())


class TestBuildLayComparison:
    def test_returns_report(self, outputs):
        a, b = outputs
        report = build_lay_comparison(a, b)
        assert isinstance(report, LayComparisonReport)

    def test_scenario_names_preserved(self, outputs):
        a, b = outputs
        report = build_lay_comparison(a, b)
        assert report.scenario_a_name == a.scenario_name
        assert report.scenario_b_name == b.scenario_name

    def test_better_scenario_is_lower_glucose(self, outputs):
        a, b = outputs
        report = build_lay_comparison(a, b)
        # high_sugar_meal has higher glucose, so light_snack should be "better"
        assert "light_snack" in report.better_scenario

    def test_headline_populated(self, outputs):
        a, b = outputs
        assert len(build_lay_comparison(a, b).headline) > 10

    def test_glucose_comparison_populated(self, outputs):
        a, b = outputs
        report = build_lay_comparison(a, b)
        assert len(report.glucose_comparison) > 10
        assert "mmol" in report.glucose_comparison

    def test_energy_comparison_populated(self, outputs):
        a, b = outputs
        assert len(build_lay_comparison(a, b).energy_comparison) > 10

    def test_recommendation_populated(self, outputs):
        a, b = outputs
        assert len(build_lay_comparison(a, b).recommendation) > 10

    def test_as_text_returns_string(self, outputs):
        a, b = outputs
        text = build_lay_comparison(a, b).as_text()
        assert isinstance(text, str)
        assert "COMPARISON" in text
        assert "RECOMMENDATION" in text

    def test_glucose_diff_sign(self, outputs):
        a, b = outputs
        report = build_lay_comparison(a, b)
        # high_sugar has higher peak than light_snack → diff > 0
        assert report.glucose_diff_mM > 0

    def test_similar_meals_neutral_headline(self):
        runner = SimulationRunner()
        # Use same scenario twice → should note similarity
        out_a = runner.run(high_sugar_meal_scenario())
        out_b = runner.run(high_sugar_meal_scenario())
        report = build_lay_comparison(out_a, out_b)
        assert "similar" in report.headline.lower() or abs(report.glucose_diff_mM) < 0.5
