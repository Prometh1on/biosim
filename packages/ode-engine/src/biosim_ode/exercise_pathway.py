# SAVE TO: packages\ode-engine\src\biosim_ode\exercise_pathway.py
"""
Exercise pathway -- time-gated coupling for an exercise bout occurring
DURING the simulated window (as opposed to before it).

biosim_physiology.exercise.ExerciseBout already computes the right
modifiers (during_exercise_modifiers(), post_exercise_modifiers()); what
was missing was somewhere to plug a TIME-VARYING version of them into
GlucoseInsulinPathway. This pathway is stateless (zero species) and exists
purely to compute output_signals(t, y) as a function of the current
simulation time relative to the bout's start/end -- CompositePathway.dydt()
already calls every coupled pathway's output_signals(t, y) with the
CURRENT integration time on every solver step (see coupling.py), so no new
ODE-solver machinery is needed, only the coupling itself.

Two distinct exercise semantics exist in this codebase, deliberately kept
separate rather than overloaded onto one field:
    - PersonalizationSpec.exercise_bout: a bout that already ENDED before
      t=0 (simulation start = some time after the bout) -- a static,
      one-shot personalisation applied via PhysiologicalPersonalizer,
      using post_exercise_modifiers() evaluated once.
    - ScenarioDefinition.exercise_bout_in_progress (this module): a bout
      that STARTS at t=0 and runs for its own duration_minutes within the
      simulated window -- a genuinely time-varying coupling, via this
      ExercisePathway.

Literature:
    Richter & Hargreaves 2013, Physiol Rev. PMID:23899560
    Dela & Kjaer 2006, Essays Biochem. PMID:17144881
"""

from __future__ import annotations

import dataclasses
from typing import TYPE_CHECKING, Optional

from .base_pathway import BasePathway, PathwayDescriptor, PathwayMaturity

if TYPE_CHECKING:
    from biosim_physiology.exercise import ExerciseBout


class ExercisePathway(BasePathway):
    """
    Stateless (zero-species) pathway injecting time-gated exercise
    signals into GlucoseInsulinPathway.

    While t <= bout.duration_minutes (the bout is in progress):
        exercise_glucose_uptake_delta <- during_exercise_modifiers()["glucose_uptake_delta"]
        secretion_scale_modifier   <- during_exercise_modifiers()["secretion_scale"]
    After the bout ends (t > bout.duration_minutes):
        p1_modifier, p3_modifier, gb_modifier <- post_exercise_modifiers(),
        with the bout's own hours_before_simulation REPLACED by the actual
        elapsed time since the bout ended (t - duration_minutes, in hours)
        -- a genuinely decaying signal, not a one-shot snapshot.
    """

    NAME    = "exercise_v1"
    VERSION = "1.0.0"

    def __init__(self, bout: "ExerciseBout") -> None:
        self.bout = bout

    @property
    def descriptor(self) -> PathwayDescriptor:
        return PathwayDescriptor(
            name=self.NAME,
            version=self.VERSION,
            maturity=PathwayMaturity.EXPERIMENTAL,
            description=(
                "Stateless time-gated exercise coupling: during-bout "
                "clearance boost, transitioning to a decaying post-bout "
                "insulin-sensitisation signal."
            ),
            species_names=[],
            output_signals=[
                "p1_modifier", "p3_modifier", "gb_modifier",
                "exercise_glucose_uptake_delta", "secretion_scale_modifier",
            ],
            input_signals=[],
            dependencies=[],
            tags=["exercise", "physiology", "time-varying"],
            literature_refs=["PMID:23899560", "PMID:17144881"],
        )

    @property
    def species_names(self) -> list[str]:
        return []

    @property
    def species_units(self) -> list[str]:
        return []

    @property
    def y0(self) -> list[float]:
        return []

    def dydt(self, t: float, y: list[float], signals: Optional[dict] = None) -> list[float]:
        return []

    def output_signals(self, t: float, y: list[float]) -> dict[str, float]:
        bout = self.bout
        if t <= bout.duration_minutes:
            during = bout.during_exercise_modifiers()
            return {
                "exercise_glucose_uptake_delta": during["glucose_uptake_delta"],
                "secretion_scale_modifier": during["secretion_scale"],
                "p1_modifier": 1.0,
                "p3_modifier": 1.0,
                "gb_modifier": 1.0,
            }
        elapsed_hours = (t - bout.duration_minutes) / 60.0
        decaying_bout = dataclasses.replace(bout, hours_before_simulation=elapsed_hours)
        post = decaying_bout.post_exercise_modifiers()
        return {
            "exercise_glucose_uptake_delta": 0.0,
            "secretion_scale_modifier": 1.0,
            "p1_modifier": post["p1"],
            "p3_modifier": post["p3"],
            "gb_modifier": post["Gb"],
        }


def exercise_couplings(gi_name: str, exercise_name: str) -> list:
    """
    Build the CouplingSpecs connecting an ExercisePathway to
    GlucoseInsulinPathway. All four output signals map straight through
    to the matching GI input of the same name -- no transform needed,
    since ExercisePathway already computes them in GI's own signal
    vocabulary (see PathwayParameters coupling contract, CLAUDE.md rule 11).
    """
    from .coupling import CouplingSpec

    return [
        CouplingSpec(
            source_pathway=exercise_name, signal_name=signal,
            target_pathway=gi_name, target_input=signal,
        )
        for signal in (
            "p1_modifier", "p3_modifier", "gb_modifier",
            "exercise_glucose_uptake_delta", "secretion_scale_modifier",
        )
    ]


def couple_with_gi(gi, exercise: "ExercisePathway"):
    """
    Build the CompositePathway connecting an ExercisePathway to
    GlucoseInsulinPathway (see exercise_couplings() for the signal contract).
    """
    from .coupling import CompositePathway

    return CompositePathway([gi, exercise], couplings=exercise_couplings(gi.name, exercise.name))
