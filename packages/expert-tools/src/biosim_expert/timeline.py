# SAVE TO: packages\expert-tools\src\biosim_expert\timeline.py
"""
Molecular timeline analysis.

Extracts the mechanistic narrative from a simulation trajectory as a
sequence of annotated events. Instead of raw time-series, this gives
researchers a structured account of what happens and when:

    t=0:   Meal ingested — glucose absorption begins
    t=18:  Glucose crosses 7 mM — beta cells at peak secretion
    t=31:  Insulin peaks at 312 pM — mTOR fully activated (0.84)
    t=45:  mTOR reaches maximum — AMPK at nadir (0.18)
    t=89:  Glucose returns to within 10% of baseline
    t=180: All species within 5% of baseline — simulation settled

This is designed for the expert researcher interface and for automated
report generation in publications.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional
import numpy as np


@dataclass
class TimelineEvent:
    """A single annotated event in the simulation timeline."""
    time_min:    float
    event_type:  str    # "meal_start", "glucose_peak", "insulin_peak", etc.
    description: str
    values:      dict[str, float] = field(default_factory=dict)
    significance: str = "normal"   # "normal", "notable", "warning"


@dataclass
class MolecularTimeline:
    """
    Full annotated timeline of a simulation.

    events are sorted by time_min ascending.
    """
    events:        list[TimelineEvent]
    duration_min:  float
    scenario_name: str

    def events_in_window(self, start_min: float, end_min: float) -> list[TimelineEvent]:
        return [e for e in self.events if start_min <= e.time_min <= end_min]

    def events_of_type(self, event_type: str) -> list[TimelineEvent]:
        return [e for e in self.events if e.event_type == event_type]

    def to_text(self) -> str:
        lines = [
            f"MOLECULAR TIMELINE: {self.scenario_name}",
            f"Duration: {self.duration_min:.0f} min",
            "",
        ]
        for e in self.events:
            flag = "" if e.significance == "normal" else f"[{e.significance.upper()}] "
            lines.append(f"  t={e.time_min:6.1f} min  {flag}{e.description}")
            for k, v in e.values.items():
                lines.append(f"              {k}: {v:.3g}")
        return "\n".join(lines)

    def to_dict(self) -> dict:
        return {
            "scenario":    self.scenario_name,
            "duration_min": self.duration_min,
            "events": [
                {
                    "time_min":    e.time_min,
                    "event_type":  e.event_type,
                    "description": e.description,
                    "values":      e.values,
                    "significance": e.significance,
                }
                for e in self.events
            ],
        }


class TimelineExtractor:
    """
    Extracts annotated events from a SimulationOutput.
    """

    def extract(self, output) -> MolecularTimeline:
        """
        Build a MolecularTimeline from a SimulationOutput.

        Args:
            output: SimulationOutput from SimulationRunner

        Returns:
            MolecularTimeline with all detected events
        """
        sim     = output.simulation
        m       = output.expert_metrics
        events: list[TimelineEvent] = []

        if not sim.success or not sim.species:
            return MolecularTimeline(events=[], duration_min=0,
                                     scenario_name=output.scenario_name)

        times_min = sim.times_s / 60.0
        G    = sim.species["G"] / 1000.0   # uM → mM
        I_pM = sim.species["I"]
        mTOR = sim.species["mTOR"]
        AMPK = sim.species["AMPK"]

        baseline_G = float(G[0])
        baseline_I = float(I_pM[0])

        # Event 1: simulation start / meal
        #
        # Regression: this used to check hasattr(output.simulation,
        # "meals") -- SimulationResult (biosim_ode.system) has no
        # `meals` field at all (never did), so this branch was always
        # unreachable and "meal_start" could never actually fire, even
        # for a real meal scenario -- exactly the event this module's own
        # docstring uses as its headline example. Fixed using
        # output.scenario (added this session specifically so consumers
        # like this one can read the real ScenarioDefinition), which
        # carries real MealEvent.time_s/carbs_g/label.
        scenario = getattr(output, "scenario", None)
        meal_events = scenario.meal_events if scenario is not None else []
        if meal_events:
            for meal in meal_events:
                events.append(TimelineEvent(
                    time_min=meal.time_s / 60.0,
                    event_type="meal_start",
                    description=f"Meal ingested — glucose absorption begins",
                    values={"carbs_g": meal.carbs_g},
                ))
        else:
            events.append(TimelineEvent(
                time_min=0.0, event_type="simulation_start",
                description="Simulation start — fasted state",
                values={"G_mM": round(baseline_G, 2), "I_pM": round(baseline_I, 1)},
            ))

        # Event 2: glucose crosses 7 mM threshold
        crosses = np.where((G[:-1] < 7.0) & (G[1:] >= 7.0))[0]
        if len(crosses) > 0:
            t_cross = float(times_min[crosses[0]])
            events.append(TimelineEvent(
                time_min=t_cross, event_type="glucose_threshold_7mM",
                description="Glucose crosses 7 mmol/L — active insulin secretion phase",
                values={"G_mM": round(float(G[crosses[0]]), 2)},
                significance="notable" if float(G[crosses[0]]) > 7.8 else "normal",
            ))

        # Event 3: glucose peak
        peak_G_idx = int(np.argmax(G))
        events.append(TimelineEvent(
            time_min=float(times_min[peak_G_idx]),
            event_type="glucose_peak",
            description=f"Peak blood glucose: {float(G[peak_G_idx]):.1f} mmol/L",
            values={
                "G_mM":  round(float(G[peak_G_idx]), 2),
                "I_pM":  round(float(I_pM[peak_G_idx]), 1),
                "mTOR":  round(float(mTOR[peak_G_idx]), 3),
            },
            significance="warning" if float(G[peak_G_idx]) > 11.1 else
                         "notable" if float(G[peak_G_idx]) > 7.8 else "normal",
        ))

        # Event 4: insulin peak
        peak_I_idx = int(np.argmax(I_pM))
        events.append(TimelineEvent(
            time_min=float(times_min[peak_I_idx]),
            event_type="insulin_peak",
            description=f"Peak insulin: {float(I_pM[peak_I_idx]):.0f} pM",
            values={
                "I_pM":  round(float(I_pM[peak_I_idx]), 1),
                "mTOR":  round(float(mTOR[peak_I_idx]), 3),
                "AMPK":  round(float(AMPK[peak_I_idx]), 3),
            },
        ))

        # Event 5: mTOR peak
        peak_mTOR_idx = int(np.argmax(mTOR))
        events.append(TimelineEvent(
            time_min=float(times_min[peak_mTOR_idx]),
            event_type="mtor_peak",
            description=f"mTOR fully activated: {float(mTOR[peak_mTOR_idx]):.2f} (fractional)",
            values={"mTOR": round(float(mTOR[peak_mTOR_idx]), 3)},
            significance="notable" if float(mTOR[peak_mTOR_idx]) > 0.7 else "normal",
        ))

        # Event 6: AMPK nadir
        nadir_AMPK_idx = int(np.argmin(AMPK))
        events.append(TimelineEvent(
            time_min=float(times_min[nadir_AMPK_idx]),
            event_type="ampk_nadir",
            description=f"AMPK minimally suppressed: {float(AMPK[nadir_AMPK_idx]):.2f}",
            values={"AMPK": round(float(AMPK[nadir_AMPK_idx]), 3)},
            significance="notable" if float(AMPK[nadir_AMPK_idx]) < 0.25 else "normal",
        ))

        # Event 7: glucose returns to near-baseline
        #
        # Regression: this used to hand-roll its own return-to-baseline
        # scan (searching forward from event 3's glucose-peak index) --
        # a second, independent implementation of what
        # ExpertMetrics.time_to_glucose_return_s already computes via
        # SimulationResult.time_to_return_to_baseline() (which anchors on
        # the point of MAXIMUM deviation from baseline, not specifically
        # the glucose peak -- the two coincide for a typical postprandial
        # rise-then-fall, but aren't guaranteed to). `m` (assigned above)
        # was never actually read anywhere in this method until now.
        if m.time_to_glucose_return_s is not None:
            t_return_min = m.time_to_glucose_return_s / 60.0
            idx_return = int(np.argmin(np.abs(times_min - t_return_min)))
            events.append(TimelineEvent(
                time_min=t_return_min,
                event_type="glucose_return_baseline",
                description=f"Blood glucose returns to within 10% of baseline",
                values={"G_mM": round(float(G[idx_return]), 2)},
            ))

        # Event 8: settled state
        events.append(TimelineEvent(
            time_min=float(times_min[-1]),
            event_type="simulation_end",
            description="Simulation complete",
            values={
                "G_mM":  round(float(G[-1]), 2),
                "I_pM":  round(float(I_pM[-1]), 1),
                "mTOR":  round(float(mTOR[-1]), 3),
                "AMPK":  round(float(AMPK[-1]), 3),
            },
        ))

        events.sort(key=lambda e: e.time_min)

        return MolecularTimeline(
            events=events,
            duration_min=float(times_min[-1]),
            scenario_name=output.scenario_name,
        )
