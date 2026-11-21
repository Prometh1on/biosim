# SAVE TO: packages\simulation-api\src\biosim_api\library.py
"""
Pre-built scenario library and comparison runner.

The scenario library provides validated, named scenarios that serve
as both convenient starting points and regression benchmarks.
Every scenario in the library has been verified to produce
physiologically plausible simulation outputs.

ComparisonRunner runs two scenarios and returns the differential,
which is the basis for the 'what if I had eaten X instead of Y?'
feature in the lay interface.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional
import numpy as np

from .output import SimulationOutput
from .runner import SimulationRunner
from .scenario import DrugEvent, MealEvent, PersonalizationSpec, ScenarioDefinition


# ---------------------------------------------------------------------------
# Pre-built scenarios
# ---------------------------------------------------------------------------

def ogtt_75g(
    personalization: Optional[PersonalizationSpec] = None,
) -> ScenarioDefinition:
    """Standard 75g oral glucose tolerance test."""
    return ScenarioDefinition(
        name="OGTT_75g",
        duration_s=14400,
        description="Standard 75g oral glucose tolerance test — clinical benchmark",
        # glycaemic_index=100 routes to GutAbsorptionModel.ogtt_calibrated_glucose()
        # (see MealEvent.to_meal_input()) -- fit against real OGTT timing
        # benchmarks, not the generic for_gi() interpolation. An earlier
        # attempt used the generic for_gi(100) preset and made validation
        # WORSE (see CLAUDE.md VALIDATION STATUS); this calibrated preset
        # instead improves peak_glucose/time_to_peak_glucose_min accuracy.
        # complex_carb_fraction=0: free glucose needs no digestion, so an
        # alpha-glucosidase inhibitor has nothing to slow (Chen 2016,
        # PMID:27686734, saw no acarbose effect on the OGTT). Changes nothing
        # unless such a drug is added.
        meal_events=[MealEvent(time_s=0.0, carbs_g=75.0, label="75g glucose load",
                               glycaemic_index=100.0, complex_carb_fraction=0.0)],
        personalization=personalization or PersonalizationSpec(),
    )


def high_sugar_meal_scenario(
    meal_time_s: float = 0.0,
    personalization: Optional[PersonalizationSpec] = None,
) -> ScenarioDefinition:
    """75g carbohydrate high-sugar meal."""
    return ScenarioDefinition(
        name="high_sugar_meal",
        duration_s=14400,
        description="High-sugar meal: 75g carbohydrate",
        meal_events=[MealEvent(time_s=meal_time_s, carbs_g=75.0, label="High-sugar meal")],
        personalization=personalization or PersonalizationSpec(),
    )


def standard_meal_scenario(
    personalization: Optional[PersonalizationSpec] = None,
) -> ScenarioDefinition:
    """50g carbohydrate standard meal."""
    return ScenarioDefinition(
        name="standard_meal",
        duration_s=14400,
        description="Standard mixed meal: 50g carbohydrate",
        meal_events=[MealEvent(time_s=0.0, carbs_g=50.0, label="Standard meal")],
        personalization=personalization or PersonalizationSpec(),
    )


def light_snack_scenario(
    personalization: Optional[PersonalizationSpec] = None,
) -> ScenarioDefinition:
    """20g carbohydrate snack."""
    return ScenarioDefinition(
        name="light_snack",
        duration_s=7200,
        description="Light snack: 20g carbohydrate",
        meal_events=[MealEvent(time_s=0.0, carbs_g=20.0, label="Light snack")],
        personalization=personalization or PersonalizationSpec(),
    )


def fasted_scenario(duration_s: float = 3600.0) -> ScenarioDefinition:
    """Fasted state — no meal. Used to verify steady state."""
    return ScenarioDefinition(
        name="fasted",
        duration_s=duration_s,
        description="Fasted state — no meal input",
        meal_events=[],
        personalization=PersonalizationSpec(),
    )


def caffeine_and_meal_scenario(
    caffeine_dose_mg: float = 200.0,
    caffeine_time_s: float = -1800.0,   # 30 min before meal (negative = before t=0)
    carbs_g: float = 50.0,
    personalization: Optional[PersonalizationSpec] = None,
) -> ScenarioDefinition:
    """
    Caffeine taken before a meal.
    Caffeine is administered before simulation start; its effect at t=0
    is computed from the hours_since_ingestion.

    Note: drug_time_s = 0 means taken at simulation start.
    A negative value (e.g. -1800) means 30 min before simulation start.
    The DrugEvent handles the timing correctly.
    """
    caffeine_sim_time = max(0.0, caffeine_time_s)  # drug events must be at t>=0
    hours_before_meal = abs(min(0.0, caffeine_time_s)) / 3600.0
    return ScenarioDefinition(
        name="caffeine_and_meal",
        duration_s=14400,
        description=(
            f"{caffeine_dose_mg:.0f}mg caffeine "
            f"{'%.0f' % (abs(caffeine_time_s)/60)} min before a {carbs_g:.0f}g meal"
        ),
        meal_events=[MealEvent(time_s=0.0, carbs_g=carbs_g, label=f"{carbs_g:.0f}g meal")],
        drug_events=[DrugEvent(
            time_s=caffeine_sim_time,
            drug="caffeine",
            dose_mg=caffeine_dose_mg,
            hours_before_simulation=hours_before_meal,
        )],
        personalization=personalization or PersonalizationSpec(),

    )


LIBRARY: dict[str, ScenarioDefinition] = {
    "ogtt_75g":           ogtt_75g(),
    "high_sugar_meal":    high_sugar_meal_scenario(),
    "standard_meal":      standard_meal_scenario(),
    "light_snack":        light_snack_scenario(),
    "fasted":             fasted_scenario(),
    "caffeine_and_meal":  caffeine_and_meal_scenario(),
}


def get_scenario(name: str) -> ScenarioDefinition:
    """Retrieve a named scenario from the library."""
    if name not in LIBRARY:
        raise KeyError(
            f"Scenario '{name}' not in library. "
            f"Available: {sorted(LIBRARY.keys())}"
        )
    return LIBRARY[name]


# ---------------------------------------------------------------------------
# Comparison runner
# ---------------------------------------------------------------------------

@dataclass
class ComparisonOutput:
    """
    Differential output from comparing two scenarios.

    glucose_diff:  scenario_a glucose minus scenario_b glucose (uM)
    Positive = scenario_a has higher glucose at that time.
    """
    scenario_a:             SimulationOutput
    scenario_b:             SimulationOutput
    times_s:                np.ndarray
    glucose_diff_uM:        np.ndarray
    insulin_diff_uM:        np.ndarray
    mtor_diff:              np.ndarray
    peak_glucose_diff_uM:   float
    peak_insulin_diff_uM:   float
    narrative:              str

    @property
    def scenario_a_higher_glucose(self) -> bool:
        return self.peak_glucose_diff_uM > 0


class ComparisonRunner:
    """
    Runs two scenarios and produces a differential output.

    The comparison uses the same time grid for both simulations
    so that element-wise subtraction is valid.
    """

    def __init__(self, n_output_points: int = 400) -> None:
        self._runner = SimulationRunner(n_output_points=n_output_points)

    def compare(
        self,
        scenario_a: ScenarioDefinition,
        scenario_b: ScenarioDefinition,
    ) -> ComparisonOutput:
        """
        Run both scenarios and compute their differential.

        Both scenarios are run to the same duration (max of their durations).
        """
        # Normalise duration to the larger of the two
        max_duration = max(scenario_a.duration_s, scenario_b.duration_s)
        if scenario_a.duration_s != max_duration:
            import dataclasses
            scenario_a = dataclasses.replace(scenario_a, duration_s=max_duration)
        if scenario_b.duration_s != max_duration:
            import dataclasses
            scenario_b = dataclasses.replace(scenario_b, duration_s=max_duration)

        out_a = self._runner.run(scenario_a)
        out_b = self._runner.run(scenario_b)

        # Use scenario_a's time grid as reference
        times_s = out_a.simulation.times_s

        # Interpolate scenario_b onto scenario_a's grid
        G_a = out_a.simulation.species["G"]
        G_b = np.interp(times_s, out_b.simulation.times_s, out_b.simulation.species["G"])
        I_a = out_a.simulation.species["I"]
        I_b = np.interp(times_s, out_b.simulation.times_s, out_b.simulation.species["I"])
        mTOR_a = out_a.simulation.species["mTOR"]
        mTOR_b = np.interp(times_s, out_b.simulation.times_s, out_b.simulation.species["mTOR"])

        g_diff = G_a - G_b
        i_diff = I_a - I_b
        m_diff = mTOR_a - mTOR_b

        peak_g_diff = float(out_a.expert_metrics.glucose_peak_uM - out_b.expert_metrics.glucose_peak_uM)
        peak_i_diff = float(out_a.expert_metrics.insulin_peak_uM - out_b.expert_metrics.insulin_peak_uM)

        direction = "higher" if peak_g_diff > 0 else "lower"
        narrative = (
            f"'{scenario_a.name}' produces {abs(peak_g_diff/1000):.1f} mM "
            f"{direction} peak blood glucose than '{scenario_b.name}'."
        )

        return ComparisonOutput(
            scenario_a=out_a,
            scenario_b=out_b,
            times_s=times_s,
            glucose_diff_uM=g_diff,
            insulin_diff_uM=i_diff,
            mtor_diff=m_diff,
            peak_glucose_diff_uM=peak_g_diff,
            peak_insulin_diff_uM=peak_i_diff,
            narrative=narrative,
        )
