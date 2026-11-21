# SAVE TO: packages\ode-engine\tests\test_no_double_counted_clearance.py
"""
B7 (2026-11-16): no resting glucose clearance on top of p1.

In the Bergman minimal model, p1 (glucose effectiveness, S_G) already stands
for ALL insulin-independent glucose disposal. The model used to add an
unsourced `extra_g_clearance` = 0.018/min on top of it at rest, counting the
same disposal twice (diagnosed in docs/VALIDATION.md, section 4). At rest that
term is now zero. Exercise keeps its own extra uptake, as a separate additive
signal of exactly the size it had before.
"""
import pytest

from biosim_ode.exercise_pathway import ExercisePathway, exercise_couplings
from biosim_ode.pathway import GlucoseInsulinPathway, PathwayParameters
from biosim_physiology.exercise import (
    EXERCISE_GLUCOSE_UPTAKE_PER_RELATIVE_INTENSITY,
    ExerciseBout,
    ExerciseIntensity,
)


def _dG_one_mM_above_basal(signals=None) -> float:
    gi = GlucoseInsulinPathway()
    y = list(gi.y0)
    y[0] = gi.params.Gb + 1.0   # G - Gb = 1 mM, X = 0, no meal
    return gi.dydt(0.0, y, signals=signals or {})[0]


class TestRestingClearanceIsP1Only:
    def test_default_extra_clearance_is_zero(self):
        assert PathwayParameters().extra_g_clearance == 0.0

    def test_resting_glucose_disposal_is_exactly_p1(self):
        assert _dG_one_mM_above_basal() == pytest.approx(-PathwayParameters().p1)

    def test_basal_state_is_still_a_fixed_point(self):
        assert GlucoseInsulinPathway().validate_y0(tolerance=1e-8) == []


class TestExerciseUptakeKeptAsItsOwnSignal:
    def test_uptake_per_intensity_equals_the_old_increment(self):
        # Old: 0.018 * (1 + 4 * intensity) during a bout, minus 0.018 at rest.
        assert EXERCISE_GLUCOSE_UPTAKE_PER_RELATIVE_INTENSITY == pytest.approx(0.018 * 4.0)

    def test_gi_adds_the_exercise_uptake_signal(self):
        delta = 0.05
        dG = _dG_one_mM_above_basal({"exercise_glucose_uptake_delta": delta})
        assert dG == pytest.approx(-(PathwayParameters().p1 + delta))

    def test_exercise_pathway_emits_uptake_during_bout_only(self):
        bout = ExerciseBout.from_intensity(ExerciseIntensity.MODERATE, 30.0)
        ex = ExercisePathway(bout)
        during = ex.output_signals(10.0, [])["exercise_glucose_uptake_delta"]
        after = ex.output_signals(60.0, [])["exercise_glucose_uptake_delta"]
        assert during == pytest.approx(
            EXERCISE_GLUCOSE_UPTAKE_PER_RELATIVE_INTENSITY * bout.relative_intensity
        )
        assert after == 0.0

    def test_couplings_carry_the_uptake_signal(self):
        names = {c.signal_name for c in exercise_couplings("gi", "ex")}
        assert "exercise_glucose_uptake_delta" in names
        assert "extra_g_clearance_modifier" not in names
