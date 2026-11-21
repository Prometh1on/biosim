# SAVE TO: packages\data-seeds\tests\test_gi_meals.py
"""
Tests for GI-adjusted meal scenarios.
"""
import pytest
from biosim_seeds.gi_meals import (
    GIFoodItem, GI_VALUES, gi_adjusted_kabs,
    gi_meal_scenario, high_gi_breakfast, low_gi_breakfast,
    compare_gi_scenarios,
)
from biosim_api import SimulationRunner


class TestGIFoodItem:
    def test_construction(self):
        item = GIFoodItem("white_rice_200g", 56.0, 73.0, "White rice")
        assert item.carbs_g == 56.0
        assert item.gi == 73.0

    def test_gi_load(self):
        item = GIFoodItem("test", 50.0, 70.0)
        assert item.gi_load == pytest.approx(35.0)

    def test_gi_category_high(self):
        assert GIFoodItem("test", 50, 75).gi_category == "high"

    def test_gi_category_medium(self):
        assert GIFoodItem("test", 50, 60).gi_category == "medium"

    def test_gi_category_low(self):
        assert GIFoodItem("test", 50, 40).gi_category == "low"

    def test_gi_category_negligible(self):
        assert GIFoodItem("test", 0, 0).gi_category == "negligible"


class TestGIAdjustedKabs:
    def test_reference_gi_gives_base_kabs(self):
        assert gi_adjusted_kabs(55.0) == pytest.approx(0.060)

    def test_high_gi_faster(self):
        assert gi_adjusted_kabs(100.0) > gi_adjusted_kabs(55.0)

    def test_low_gi_slower(self):
        assert gi_adjusted_kabs(30.0) < gi_adjusted_kabs(55.0)

    def test_gi_100_scales_correctly(self):
        assert gi_adjusted_kabs(100.0) == pytest.approx(0.060 * 100.0/55.0, rel=0.01)


class TestGIMealScenario:
    def test_high_gi_breakfast_items(self):
        items = high_gi_breakfast()
        assert len(items) >= 1
        assert any(i.gi >= 70 for i in items)

    def test_low_gi_breakfast_items(self):
        items = low_gi_breakfast()
        assert len(items) >= 1
        assert all(i.gi < 70 for i in items)

    def test_scenario_construction(self):
        items = high_gi_breakfast()
        s = gi_meal_scenario(items, name="test_gi")
        assert s.name == "test_gi"
        assert len(s.meal_events) == 1

    def test_empty_meal_raises(self):
        with pytest.raises(ValueError):
            gi_meal_scenario([], name="empty")

    def test_compare_gi_scenarios(self):
        high, low = compare_gi_scenarios()
        assert high.name == "high_gi_breakfast"
        assert low.name == "low_gi_breakfast"

    def test_high_gi_has_more_carbs_or_similar(self):
        high_items = high_gi_breakfast()
        low_items  = low_gi_breakfast()
        high_carbs = sum(i.carbs_g for i in high_items)
        low_carbs  = sum(i.carbs_g for i in low_items)
        # Both should have meaningful carb content
        assert high_carbs > 10
        assert low_carbs  > 10

    def test_gi_values_database_not_empty(self):
        assert len(GI_VALUES) > 10
        assert "white_bread_slice" in GI_VALUES
        assert "apple_medium" in GI_VALUES

    def test_gi_scenario_runs(self):
        items = high_gi_breakfast()
        scenario = gi_meal_scenario(items, name="gi_test")
        output = SimulationRunner().run(scenario)
        assert output.success
        assert output.expert_metrics.glucose_peak_mM > 5.0

    def test_gi_meal_scenario_carries_computed_kabs(self):
        # Regression: gi_meal_scenario() computed a GI-weighted kabs via
        # gi_adjusted_kabs() but the MealEvent it built had no field to
        # carry it, so the value was silently discarded and every
        # GI-adjusted scenario used the default kabs regardless of the
        # meal's actual glycaemic index.
        items = high_gi_breakfast()
        scenario = gi_meal_scenario(items, name="gi_test")
        weighted_gi = sum(i.gi * i.carbs_g for i in items) / sum(i.carbs_g for i in items)
        assert scenario.meal_events[0].kabs == pytest.approx(gi_adjusted_kabs(weighted_gi))

    def test_high_gi_peaks_earlier_than_low_gi_same_carb_load(self):
        # The actual physiological claim gi_meal_scenario() exists to
        # demonstrate: a higher-GI meal should produce an earlier glucose
        # peak than a lower-GI meal, for comparable carb content -- this
        # was silently false before the kabs-passthrough fix above (both
        # scenarios peaked at the identical simulation timestep).
        high_scenario = gi_meal_scenario(
            [GIFoodItem("high", carbs_g=50.0, gi=90.0)], name="high"
        )
        low_scenario = gi_meal_scenario(
            [GIFoodItem("low", carbs_g=50.0, gi=30.0)], name="low"
        )
        out_high = SimulationRunner().run(high_scenario)
        out_low = SimulationRunner().run(low_scenario)
        assert out_high.expert_metrics.glucose_peak_time_s < out_low.expert_metrics.glucose_peak_time_s


class TestGICategoryConsistencyWithLayOutput:
    """
    biosim_lay.FOOD_DB (the lay/dietary-guidance layer) stores a
    hardcoded qualitative gi_category per food, completely independent
    of this module's GI_VALUES (the numeric table actually driving
    gi_adjusted_kabs()'s absorption-rate math) -- the same food's GI
    fact duplicated across two packages with no cross-check, the same
    "can silently drift apart" shape already caught this session for
    Vg_L_per_kg and MODEL_VERSION. It had: cola_330ml (GI=63, a real
    "medium") mislabeled "high"; orange_juice_200ml (GI=50, "low")
    mislabeled "high"; banana_medium (GI=51), wholegrain_bread (GI=53),
    and grapes_100g (GI=46, all real "low") mislabeled "medium". A user
    asking biosim_lay for low-GI food swaps would never see foods its
    own numeric data (used identically in the ODE simulation layer)
    says qualify -- the opposite of the guidance the platform's own
    math implies. This package already depends on biosim_lay (see
    gi_meals.py's `from biosim_lay.food_items import FOOD_DB`), so this
    lives here rather than in biosim_lay to avoid inverting that
    dependency direction just to share one table.
    """
    def test_food_db_categories_match_gi_values_numeric_ground_truth(self):
        from biosim_lay.food_items import FOOD_DB, gi_category_for_value

        mismatches = []
        for food_id, gi in GI_VALUES.items():
            if food_id not in FOOD_DB or FOOD_DB[food_id].gi_category == "negligible":
                continue   # negligible-carb foods have no meaningful GI category
            expected = gi_category_for_value(gi)
            actual = FOOD_DB[food_id].gi_category
            if expected != actual:
                mismatches.append((food_id, gi, expected, actual))
        assert mismatches == [], (
            f"biosim_lay.FOOD_DB gi_category out of sync with GI_VALUES: {mismatches}"
        )


class TestGIPlaceholderLabel:
    """
    gi_meal_scenario() turns GI into an absorption rate with its own
    uncalibrated linear formula (gi_adjusted_kabs), so the runner's
    for_gi() warning never fires for it. The placeholder status must still
    reach every user-facing output: the scenario description (expert/REST)
    and the plain-language lay report.
    """

    def test_gi_meal_scenario_description_says_placeholder(self):
        from biosim_ode.gut_absorption import GI_SLOPE_PLACEHOLDER_NOTE
        scenario = gi_meal_scenario(low_gi_breakfast(), name="low_gi")
        assert GI_SLOPE_PLACEHOLDER_NOTE in scenario.description

    def test_lay_report_for_gi_meal_says_placeholder(self):
        from biosim_lay import generate_report
        from biosim_ode.gut_absorption import GI_SLOPE_PLACEHOLDER_NOTE
        output = SimulationRunner().run(gi_meal_scenario(low_gi_breakfast(), name="low_gi"))
        report = generate_report(output)
        assert GI_SLOPE_PLACEHOLDER_NOTE in report.glucose_section
