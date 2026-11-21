# SAVE TO: packages\data-seeds\src\biosim_seeds\circadian_optimizer.py
"""
Circadian meal optimizer.

Uses the platform's circadian model to find the time of day that
minimises postprandial glucose excursion for a given meal and
physiological profile.

The optimizer runs the simulation at multiple times of day and
returns a ranked schedule. This is the "what time should I eat?"
feature.

Key finding from chronobiology literature:
    The same meal produces the lowest glucose peak at 10:00-12:00
    (peak insulin sensitivity) and the highest peak at 21:00-23:00
    (nadir of insulin sensitivity). The difference is ~15-20% in
    peak glucose for the same carbohydrate load.

Usage:
    optimizer = CircadianOptimizer()
    schedule = optimizer.optimize(carbs_g=60.0, n_times=8)
    print(schedule.best_time_description)
    print(schedule.ranked_times)

Reference:
    Morris et al. 2015, PNAS. PMID:25870289
    Saad et al. 2012, Diabetes. PMID:22751690
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional
import numpy as np


@dataclass
class TimeSlot:
    """One time-of-day evaluation in the optimizer."""
    hour:              float
    glucose_peak_mM:   float
    insulin_peak_pM:   float
    energy_score:      float   # from lay summary energy_level (0-10)
    has_energy_dip:    bool
    success:           bool
    # Regression: a failed slot's exception/failure reason used to be
    # silently discarded (bare `except Exception: ...` with no field to
    # even carry a reason) -- same bug class as CohortSummary
    # .failure_reasons' bandit B110 fix (population.py). None for
    # successful slots.
    failure_reason:    Optional[str] = None


@dataclass
class CircadianSchedule:
    """
    Result of circadian optimization — ranked eating times.
    """
    carbs_g:     float
    time_slots:  list[TimeSlot]
    best_hour:   float
    worst_hour:  float

    @property
    def ranked_slots(self) -> list[TimeSlot]:
        """Slots sorted by glucose peak ascending (best first)."""
        return sorted(
            [s for s in self.time_slots if s.success],
            key=lambda s: s.glucose_peak_mM,
        )

    @property
    def best_time_description(self) -> str:
        h = int(self.best_hour)
        m = int((self.best_hour - h) * 60)
        period = "AM" if h < 12 else "PM"
        h12 = h if h <= 12 else h - 12
        h12 = h12 or 12
        return f"{h12}:{m:02d} {period}"

    @property
    def glucose_range_mM(self) -> float:
        """Difference between best and worst glucose peaks."""
        good = [s.glucose_peak_mM for s in self.time_slots if s.success]
        if not good:
            return 0.0
        return max(good) - min(good)

    @property
    def failure_reasons(self) -> list[str]:
        """Reasons for any time slot that failed to simulate, if any."""
        return [
            f"{s.hour:.1f}h: {s.failure_reason}"
            for s in self.time_slots
            if not s.success and s.failure_reason
        ]

    def to_dict(self) -> dict:
        return {
            "carbs_g":         self.carbs_g,
            "best_hour":       self.best_hour,
            "worst_hour":      self.worst_hour,
            "glucose_range_mM": round(self.glucose_range_mM, 2),
            "best_time":       self.best_time_description,
            "failure_reasons": self.failure_reasons,
            "ranked_slots": [
                {
                    "hour":            s.hour,
                    "glucose_peak_mM": round(s.glucose_peak_mM, 2),
                    "energy_score":    round(s.energy_score, 2),
                    "has_energy_dip":  s.has_energy_dip,
                }
                for s in self.ranked_slots
            ],
        }

    def narrative(self) -> str:
        if not self.ranked_slots:
            return "Optimization failed — no successful simulations."
        best  = self.ranked_slots[0]
        worst = self.ranked_slots[-1]
        diff  = worst.glucose_peak_mM - best.glucose_peak_mM
        return (
            f"For a {self.carbs_g:.0f}g carbohydrate meal, the best eating time "
            f"is around {self.best_time_description} "
            f"(peak blood sugar: {best.glucose_peak_mM:.1f} mmol/L). "
            f"Eating at {int(self.worst_hour):02d}:00 produces a "
            f"{diff:.1f} mmol/L higher peak "
            f"({worst.glucose_peak_mM:.1f} mmol/L). "
            f"This reflects the natural circadian rhythm of insulin sensitivity."
        )


class CircadianOptimizer:
    """
    Finds the optimal eating time for a given meal and physiological profile.

    Sweeps through times of day, runs a simulation at each, and ranks
    by postprandial glucose peak.

    Args:
        n_times:  Number of time points to evaluate (default: 8 = every 3 hours)
        profile:  Optional PhysiologicalProfile (defaults to healthy 30-year-old)
    """

    def __init__(
        self,
        n_times: int = 8,
        profile=None,
    ) -> None:
        self.n_times = n_times
        self.profile = profile

    def optimize(
        self,
        carbs_g: float = 60.0,
        duration_s: float = 14400.0,
    ) -> CircadianSchedule:
        """
        Run the optimizer for a given carbohydrate load.

        Args:
            carbs_g:    Grams of carbohydrate in the meal
            duration_s: Simulation duration in seconds

        Returns:
            CircadianSchedule with ranked times and narrative
        """
        from biosim_api import MealEvent, PersonalizationSpec, ScenarioDefinition, SimulationRunner
        from biosim_physiology import PhysiologicalProfile

        hours = np.linspace(0, 24, self.n_times, endpoint=False)
        runner = SimulationRunner(n_output_points=200)
        slots: list[TimeSlot] = []

        for hour in hours:
            profile = self.profile or PhysiologicalProfile(time_of_day_hours=float(hour))
            # If a profile was provided, clone it with the new time
            if self.profile is not None:
                import dataclasses
                profile = dataclasses.replace(self.profile, time_of_day_hours=float(hour))

            scenario = ScenarioDefinition(
                name=f"optim_{hour:.1f}h",
                duration_s=duration_s,
                meal_events=[MealEvent(time_s=0.0, carbs_g=carbs_g)],
                personalization=PersonalizationSpec(physiological_profile=profile),
            )
            try:
                output = runner.run(scenario)
                if output.success:
                    m  = output.expert_metrics
                    ls = output.lay_summary
                    energy_peak = float(np.max(ls.energy_level))
                    slots.append(TimeSlot(
                        hour=float(hour),
                        glucose_peak_mM=m.glucose_peak_mM,
                        insulin_peak_pM=m.insulin_peak_pM,
                        energy_score=energy_peak,
                        has_energy_dip=ls.has_energy_dip,
                        success=True,
                    ))
                else:
                    slots.append(TimeSlot(
                        float(hour), 0, 0, 5, False, False,
                        failure_reason="; ".join(output.warnings or ["simulation reported failure, no warning detail"]),
                    ))
            except Exception as e:
                slots.append(TimeSlot(
                    float(hour), 0, 0, 5, False, False,
                    failure_reason=f"{type(e).__name__}: {e}",
                ))

        succeeded = [s for s in slots if s.success]
        if not succeeded:
            return CircadianSchedule(carbs_g=carbs_g, time_slots=slots,
                                     best_hour=10.0, worst_hour=22.0)

        best  = min(succeeded, key=lambda s: s.glucose_peak_mM)
        worst = max(succeeded, key=lambda s: s.glucose_peak_mM)

        return CircadianSchedule(
            carbs_g=carbs_g,
            time_slots=slots,
            best_hour=best.hour,
            worst_hour=worst.hour,
        )
