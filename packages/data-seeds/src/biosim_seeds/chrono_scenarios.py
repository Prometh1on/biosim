# SAVE TO: packages\data-seeds\src\biosim_seeds\chrono_scenarios.py
"""
Chronobiological meal scenarios.

The same meal at 08:00 produces a lower glucose peak than at 20:00
due to circadian variation in insulin sensitivity (~30% difference).
"""

from __future__ import annotations
from biosim_api import MealEvent, PersonalizationSpec, ScenarioDefinition
from biosim_physiology import PhysiologicalProfile


def morning_meal(carbs_g: float = 75.0, hour: float = 8.0) -> ScenarioDefinition:
    return ScenarioDefinition(
        name=f"morning_meal_{carbs_g:.0f}g",
        duration_s=14400,
        description=f"{carbs_g:.0f}g meal at {hour:.0f}:00 (high insulin sensitivity)",
        meal_events=[MealEvent(time_s=0.0, carbs_g=carbs_g, label=f"Morning meal {hour:.0f}:00")],
        personalization=PersonalizationSpec(
            physiological_profile=PhysiologicalProfile(time_of_day_hours=hour)
        ),
    )


def evening_meal(carbs_g: float = 75.0, hour: float = 20.0) -> ScenarioDefinition:
    return ScenarioDefinition(
        name=f"evening_meal_{carbs_g:.0f}g",
        duration_s=14400,
        description=f"{carbs_g:.0f}g meal at {hour:.0f}:00 (lower insulin sensitivity)",
        meal_events=[MealEvent(time_s=0.0, carbs_g=carbs_g, label=f"Evening meal {hour:.0f}:00")],
        personalization=PersonalizationSpec(
            physiological_profile=PhysiologicalProfile(time_of_day_hours=hour)
        ),
    )


def meal_at_hour(carbs_g: float, hour: float) -> ScenarioDefinition:
    return ScenarioDefinition(
        name=f"meal_{carbs_g:.0f}g_at_{hour:.0f}h",
        duration_s=14400,
        description=f"{carbs_g:.0f}g meal at {hour:.0f}:00",
        meal_events=[MealEvent(time_s=0.0, carbs_g=carbs_g, label=f"Meal at {hour:.0f}:00")],
        personalization=PersonalizationSpec(
            physiological_profile=PhysiologicalProfile(time_of_day_hours=hour)
        ),
    )


def full_day_scenarios(carbs_per_meal: float = 60.0) -> list[ScenarioDefinition]:
    return [
        meal_at_hour(carbs_per_meal, 7.5),
        meal_at_hour(carbs_per_meal, 12.5),
        meal_at_hour(carbs_per_meal, 19.0),
    ]


def shift_worker_late_meal(carbs_g: float = 60.0) -> ScenarioDefinition:
    return ScenarioDefinition(
        name=f"shift_worker_meal_{carbs_g:.0f}g",
        duration_s=14400,
        description=f"{carbs_g:.0f}g meal at 03:00 (minimum insulin sensitivity)",
        meal_events=[MealEvent(time_s=0.0, carbs_g=carbs_g, label="Night-shift meal")],
        personalization=PersonalizationSpec(
            physiological_profile=PhysiologicalProfile(time_of_day_hours=3.0)
        ),
    )
