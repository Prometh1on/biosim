# SAVE TO: packages\data-seeds\src\biosim_seeds\gi_meals.py
"""
Glycaemic index (GI) adjusted meal scenarios.

GI measures how quickly a food raises blood glucose relative to pure glucose.
High-GI foods (white bread, sugary drinks) cause rapid absorption and sharp peaks.
Low-GI foods (lentils, oats, pasta) cause slower absorption and flatter curves.

This module provides:
    - GI data for common foods (from Atkinson et al. 2008 GI database)
    - gi_adjusted_scenario(): builds a scenario with GI-adjusted absorption
    - GIMealComposer: composes meals from food items with GI-weighted kabs

The GI values are literature values; the GI-to-absorption-rate slope is an
UNVALIDATED PLACEHOLDER (only GI=100 is calibrated anywhere in the platform).

Literature:
    Atkinson et al. 2008, Diabetes Care. PMID:18835944 (International GI tables)
    Jenkins et al. 1981, Am J Clin Nutr. PMID:6259925 (original GI concept)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from biosim_api import MealEvent, PersonalizationSpec, ScenarioDefinition
from biosim_ode.gut_absorption import GI_SLOPE_PLACEHOLDER_NOTE


# ---------------------------------------------------------------------------
# GI database (glucose = 100, white bread = 75 reference scale)
# ---------------------------------------------------------------------------

GI_VALUES: dict[str, float] = {
    # High GI (> 70)
    "glucose_75g_ogtt":   100.0,
    "white_bread_slice":   75.0,
    "white_rice_200g":     73.0,
    "cola_330ml":          63.0,  # borderline high
    "potato_medium":       85.0,
    "sports_drink_500ml":  78.0,

    # Medium GI (55-70)
    "pasta_200g":          55.0,
    "brown_rice_200g":     68.0,
    "banana_medium":       51.0,  # borderline medium
    "mango_half":          60.0,
    "oats_40g":            55.0,
    "orange_juice_200ml":  50.0,  # borderline medium

    # Low GI (< 55)
    "apple_medium":        36.0,
    "orange_medium":       43.0,
    "milk_200ml":          37.0,
    "wholegrain_bread":    53.0,
    "grapes_100g":         46.0,

    # Negligible GI
    "egg_boiled":           0.0,
    "chicken_breast":       0.0,
    "salad_plain":          5.0,
    "water_any":            0.0,
}


@dataclass(frozen=True)
class GIFoodItem:
    """A food item with GI information."""
    food_id:    str
    carbs_g:    float
    gi:         float
    label:      str = ""

    @property
    def gi_load(self) -> float:
        """Glycaemic load = (GI × carbs_g) / 100."""
        return (self.gi * self.carbs_g) / 100.0

    @property
    def gi_category(self) -> str:
        if self.gi >= 70: return "high"
        if self.gi >= 55: return "medium"
        if self.gi > 0:   return "low"
        return "negligible"


def gi_adjusted_kabs(gi: float, base_kabs: float = 0.060) -> float:
    """
    Compute absorption rate from glycaemic index. GI 55 = base_kabs.

    UNVALIDATED PLACEHOLDER: the linear GI-to-kabs slope has no literature
    calibration (see GI_SLOPE_PLACEHOLDER_NOTE in biosim_ode.gut_absorption).
    """
    gi = max(5.0, min(100.0, gi))
    return base_kabs * (gi / 55.0)


def gi_meal_scenario(
    food_items: list[GIFoodItem],
    name: str = "gi_meal",
    duration_s: float = 14400.0,
    profile=None,
) -> ScenarioDefinition:
    """
    Build a scenario where absorption rate is weighted by the meal's GI.

    The meal's overall kabs is the carbohydrate-weighted mean GI of all foods,
    then converted to a kabs via gi_adjusted_kabs().

    Args:
        food_items: List of GIFoodItem composing the meal
        name:       Scenario name
        duration_s: Simulation duration in seconds
        profile:    Optional PhysiologicalProfile

    Returns:
        ScenarioDefinition with GI-adjusted MealInput
    """
    total_carbs = sum(f.carbs_g for f in food_items)
    if total_carbs <= 0:
        raise ValueError("Meal has no carbohydrates")

    # Weighted mean GI
    weighted_gi = sum(f.gi * f.carbs_g for f in food_items) / total_carbs
    kabs = gi_adjusted_kabs(weighted_gi)

    label = " + ".join(f.label or f.food_id for f in food_items)
    gi_label = f"GI={weighted_gi:.0f} GL={sum(f.gi_load for f in food_items):.0f}"

    # Create a MealEvent that carries GI-adjusted kabs
    meal_event = MealEvent(
        time_s=0.0,
        carbs_g=total_carbs,
        label=f"{label} ({gi_label})",
        kabs=kabs,
    )

    return ScenarioDefinition(
        name=name,
        duration_s=duration_s,
        description=(
            f"GI-adjusted meal: {label}. Weighted GI={weighted_gi:.0f}. "
            f"{GI_SLOPE_PLACEHOLDER_NOTE}"
        ),
        meal_events=[meal_event],
        personalization=PersonalizationSpec(
            physiological_profile=profile
        ) if profile else PersonalizationSpec(),
    )


# ---------------------------------------------------------------------------
# Pre-built GI meal comparisons
# ---------------------------------------------------------------------------

def high_gi_breakfast() -> list[GIFoodItem]:
    """White toast + OJ + sugary cereal — high GI breakfast."""
    from biosim_lay.food_items import FOOD_DB
    return [
        GIFoodItem("white_bread_slice", 14.0 * 2, GI_VALUES["white_bread_slice"], "2 slices white toast"),
        GIFoodItem("orange_juice_200ml", 20.0, GI_VALUES["orange_juice_200ml"], "OJ 200ml"),
    ]


def low_gi_breakfast() -> list[GIFoodItem]:
    """Oat porridge + whole milk + apple — low GI breakfast."""
    return [
        GIFoodItem("oats_40g", 27.0, GI_VALUES["oats_40g"], "Porridge oats"),
        GIFoodItem("milk_200ml", 10.0, GI_VALUES["milk_200ml"], "Whole milk"),
        GIFoodItem("apple_medium", 20.0, GI_VALUES["apple_medium"], "Apple"),
    ]


def compare_gi_scenarios() -> tuple[ScenarioDefinition, ScenarioDefinition]:
    """
    Returns (high_gi_scenario, low_gi_scenario) for the same total carb load.
    Useful for demonstrating the GI effect.
    """
    high = gi_meal_scenario(high_gi_breakfast(), name="high_gi_breakfast")
    low  = gi_meal_scenario(low_gi_breakfast(),  name="low_gi_breakfast")
    return high, low
