# SAVE TO: packages\data-seeds\src\biosim_seeds\exercise_scenarios.py
"""
Exercise scenario builder.

Pre-built scenarios combining exercise bouts with meal simulations.
Shows how the same meal produces different glucose dynamics depending
on whether exercise preceded it.

Key insight:
    A 45-minute moderate run 1 hour before eating measurably reduces the
    postprandial glucose peak compared to a sedentary control, via
    enhanced GLUT4 expression and residual AMPK activation
    (ExerciseBout.post_exercise_modifiers()). NOTE: this scenario builder
    did not actually apply that modifier at all until it was wired up (see
    CLAUDE.md KNOWN BUGS FIXED) -- exercise and sedentary scenarios
    produced bit-for-bit identical output. Now fixed and verified to
    produce a real, directional effect, but the actual magnitude with
    default parameters (peak -1.52%, AUC -1.89% for this example --
    re-measured after the 2026-08-12 p3-coefficient recalibration below,
    essentially unchanged from before it) is far smaller than an
    earlier, never-validated "~15-25%" estimate in this docstring --
    corrected to avoid re-asserting an unverified number.
    post_exercise_modifiers()'s p3 coefficient was recalibrated
    (2026-08-12 deep-search plan, item B1) from an uncited 0.50 to 0.31,
    anchored to two real clamp studies (Mikines et al. 1988,
    PMID:3126668; Wang et al. 2013, PMID:23493357) -- see that method's
    own docstring for the full derivation and disclosed intensity-
    mapping assumptions. The p1/Gb coefficients (0.10/0.03) remain
    uncited. There is still no real clinical post-exercise OGTT
    benchmark in this codebase to fit the composed effect against
    (unlike the OGTT/IVGTT benchmarks in ModelValidator).

Scenarios:
    post_exercise_meal()     -- meal after a run (any intensity)
    sedentary_meal()         -- same meal, no exercise
    fasted_run()             -- exercise in fasted state
    compare_exercise_effect() -- returns (exercise, no_exercise) pair
"""

from __future__ import annotations

from typing import Optional

from biosim_api import MealEvent, PersonalizationSpec, ScenarioDefinition
from biosim_physiology import ExerciseBout, ExerciseIntensity, PhysiologicalProfile


def post_exercise_meal(
    carbs_g: float = 75.0,
    intensity: ExerciseIntensity = ExerciseIntensity.MODERATE,
    duration_minutes: float = 45.0,
    hours_before_meal: float = 1.0,
    profile: Optional[PhysiologicalProfile] = None,
) -> ScenarioDefinition:
    """
    Meal eaten after an exercise bout.

    The exercise bout occurred `hours_before_meal` hours before the simulation.
    Post-exercise insulin sensitisation is applied to the parameters.

    Args:
        carbs_g:           Carbohydrate content of the meal
        intensity:         Exercise intensity
        duration_minutes:  Exercise duration
        hours_before_meal: Hours between end of exercise and start of meal
        profile:           Optional physiological profile
    """
    bout = ExerciseBout.from_intensity(intensity, duration_minutes, hours_before_meal)
    label = f"{intensity.value} exercise {duration_minutes:.0f}min, {hours_before_meal:.0f}h before meal"

    return ScenarioDefinition(
        name=f"post_exercise_meal_{intensity.value}",
        duration_s=14400,
        description=f"{carbs_g:.0f}g meal after {label}. Enhanced insulin sensitivity expected.",
        meal_events=[MealEvent(time_s=0.0, carbs_g=carbs_g, label=f"{carbs_g:.0f}g post-exercise meal")],
        personalization=PersonalizationSpec(physiological_profile=profile, exercise_bout=bout),
    )


def sedentary_meal(
    carbs_g: float = 75.0,
    profile: Optional[PhysiologicalProfile] = None,
) -> ScenarioDefinition:
    """Control meal with no prior exercise."""
    return ScenarioDefinition(
        name="sedentary_meal",
        duration_s=14400,
        description=f"{carbs_g:.0f}g meal, no prior exercise (sedentary control).",
        meal_events=[MealEvent(time_s=0.0, carbs_g=carbs_g, label=f"{carbs_g:.0f}g sedentary meal")],
        personalization=PersonalizationSpec(physiological_profile=profile) if profile
        else PersonalizationSpec(),
    )


def fasted_run(
    intensity: ExerciseIntensity = ExerciseIntensity.MODERATE,
    duration_minutes: float = 45.0,
) -> ScenarioDefinition:
    """
    Fasted exercise simulation (no meal, just the exercise state).

    Uses `ScenarioDefinition.exercise_bout_in_progress` (NOT
    `PersonalizationSpec.exercise_bout`, which models a bout that already
    ended before t=0) -- this bout starts at t=0 and runs for
    `duration_minutes` within the simulated window, composing a real
    time-gated `ExercisePathway`: `during_exercise_modifiers()` (elevated
    clearance) while the bout is in progress, transitioning to a decaying
    `post_exercise_modifiers()`-equivalent signal afterward. See
    `biosim_ode.exercise_pathway` and CLAUDE.md rule 11.
    """
    bout = ExerciseBout.from_intensity(intensity, duration_minutes)
    return ScenarioDefinition(
        name=f"fasted_run_{intensity.value}",
        duration_s=7200,
        description=f"{duration_minutes:.0f}min {intensity.value} run, fasted state.",
        meal_events=[],
        personalization=PersonalizationSpec(),
        exercise_bout_in_progress=bout,
    )


def compare_exercise_effect(
    carbs_g: float = 75.0,
    intensity: ExerciseIntensity = ExerciseIntensity.MODERATE,
    duration_minutes: float = 45.0,
    hours_before_meal: float = 1.0,
) -> tuple[ScenarioDefinition, ScenarioDefinition]:
    """
    Returns (post_exercise_scenario, sedentary_scenario) for direct comparison.

    Both scenarios have the same meal and duration. The exercise scenario
    has enhanced insulin sensitivity from the prior bout.
    """
    return (
        post_exercise_meal(carbs_g, intensity, duration_minutes, hours_before_meal),
        sedentary_meal(carbs_g),
    )
