# SAVE TO: packages\ode-engine\tests\test_exercise_pathway.py
"""
Test suite for the time-gated exercise pathway extension.

See CLAUDE.md KNOWN BUGS FIXED #16: the exercise personalisation layer
was entirely disconnected for a whole session before being wired up;
this pathway specifically closes the remaining "exercise happening
DURING the simulated window" gap (fasted_run()), distinct from
PersonalizationSpec.exercise_bout's static pre-simulation treatment.
"""
import pytest
from biosim_ode import GlucoseInsulinPathway, ODESystem, ExercisePathway
from biosim_ode.exercise_pathway import couple_with_gi, exercise_couplings
from biosim_physiology import ExerciseBout, ExerciseIntensity


class TestExercisePathway:
    def test_zero_species(self):
        bout = ExerciseBout.from_intensity(ExerciseIntensity.MODERATE, 45.0)
        p = ExercisePathway(bout)
        assert p.species_names == []
        assert p.y0 == []

    def test_dydt_returns_empty(self):
        bout = ExerciseBout.from_intensity(ExerciseIntensity.MODERATE, 45.0)
        p = ExercisePathway(bout)
        assert p.dydt(0.0, []) == []

    def test_during_bout_signals(self):
        bout = ExerciseBout.from_intensity(ExerciseIntensity.MODERATE, 45.0)
        p = ExercisePathway(bout)
        signals = p.output_signals(20.0, [])   # t=20 < duration_minutes=45
        assert signals["exercise_glucose_uptake_delta"] > 0.0
        assert signals["secretion_scale_modifier"] < 1.0
        assert signals["p1_modifier"] == pytest.approx(1.0)
        assert signals["p3_modifier"] == pytest.approx(1.0)
        assert signals["gb_modifier"] == pytest.approx(1.0)

    def test_after_bout_signals(self):
        bout = ExerciseBout.from_intensity(ExerciseIntensity.MODERATE, 45.0)
        p = ExercisePathway(bout)
        signals = p.output_signals(100.0, [])   # t=100 > duration_minutes=45
        assert signals["exercise_glucose_uptake_delta"] == pytest.approx(0.0)
        assert signals["secretion_scale_modifier"] == pytest.approx(1.0)
        assert signals["p1_modifier"] > 1.0
        assert signals["p3_modifier"] > 1.0
        assert signals["gb_modifier"] < 1.0

    def test_post_bout_effect_decays_with_elapsed_time(self):
        # The whole point of making this time-gated rather than a static
        # snapshot: the post-bout effect should weaken as t advances
        # further past the bout's end, not stay constant.
        bout = ExerciseBout.from_intensity(ExerciseIntensity.VIGOROUS, 45.0)
        p = ExercisePathway(bout)
        just_after = p.output_signals(50.0, [])["p3_modifier"]
        much_later = p.output_signals(600.0, [])["p3_modifier"]
        assert just_after > much_later > 1.0

    def test_transition_is_at_duration_minutes_boundary(self):
        bout = ExerciseBout.from_intensity(ExerciseIntensity.MODERATE, 45.0)
        p = ExercisePathway(bout)
        at_boundary = p.output_signals(45.0, [])
        just_after = p.output_signals(45.001, [])
        assert at_boundary["exercise_glucose_uptake_delta"] > 0.0   # still "during" at t==duration
        assert just_after["exercise_glucose_uptake_delta"] == pytest.approx(0.0)

    def test_couplings_forward_all_five_signals(self):
        # Locks the exact set so a future signal addition here can't
        # silently forget to also add it to exercise_couplings() -- the
        # same drift this fix itself closed for secretion_scale_modifier.
        specs = exercise_couplings("gi", "exercise")
        forwarded = {spec.signal_name for spec in specs}
        assert forwarded == {
            "p1_modifier", "p3_modifier", "gb_modifier",
            "exercise_glucose_uptake_delta", "secretion_scale_modifier",
        }

    def test_implements_base_pathway(self):
        from biosim_ode.base_pathway import BasePathway
        bout = ExerciseBout.from_intensity(ExerciseIntensity.LIGHT, 20.0)
        assert isinstance(ExercisePathway(bout), BasePathway)


class TestExerciseCoupledWithGI:
    def test_composite_y0_length_unchanged(self):
        # Zero-species pathway must contribute nothing to the composite
        # state vector.
        bout = ExerciseBout.from_intensity(ExerciseIntensity.MODERATE, 45.0)
        composite = couple_with_gi(GlucoseInsulinPathway(), ExercisePathway(bout))
        assert len(composite.y0) == 6

    def test_couple_with_gi_preserves_basal_steady_state(self):
        bout = ExerciseBout.from_intensity(ExerciseIntensity.VIGOROUS, 45.0)
        composite = couple_with_gi(GlucoseInsulinPathway(), ExercisePathway(bout))
        derivs = composite.dydt(0.0, composite.y0)
        assert derivs == pytest.approx([0.0] * 6, abs=1e-9)

    def test_secretion_scale_modifier_wired_through_to_gi(self):
        # Regression: during_exercise_modifiers()["secretion_scale"] was
        # computed but ExercisePathway.output_signals() dropped it before
        # this fix -- only extra_g_clearance_modifier was forwarded.
        # Confirms the coupling actually lowers insulin secretion during a
        # meal eaten mid-bout, not just that the raw signal value looks
        # right in isolation (see test_pathway_expansion.py's dydt-level
        # tests for that). Deliberately reconstructs the pre-fix behaviour
        # by patching output_signals() to strip the key, rather than
        # asserting against a hardcoded expected trajectory.
        from unittest.mock import patch
        from biosim_ode.meal import MealInput

        bout = ExerciseBout.from_intensity(ExerciseIntensity.VIGOROUS, 45.0)
        meals = [MealInput(carbs_g=50.0, time_min=0.0)]

        def run(strip_secretion_signal: bool) -> float:
            exercise = ExercisePathway(bout)
            if strip_secretion_signal:
                original = ExercisePathway.output_signals
                def patched(self, t, y):
                    signals = original(self, t, y)
                    signals.pop("secretion_scale_modifier", None)
                    return signals
                with patch.object(ExercisePathway, "output_signals", patched):
                    composite = couple_with_gi(GlucoseInsulinPathway(meals=meals), exercise)
                    result = ODESystem(pathway=composite).simulate(duration_s=3600, meals=meals)
            else:
                composite = couple_with_gi(GlucoseInsulinPathway(meals=meals), exercise)
                result = ODESystem(pathway=composite).simulate(duration_s=3600, meals=meals)
            assert result.success
            return result.species["I"].max()

        peak_insulin_without_wiring = run(strip_secretion_signal=True)
        peak_insulin_with_wiring = run(strip_secretion_signal=False)
        assert peak_insulin_with_wiring < peak_insulin_without_wiring

    def test_fasted_exercise_lowers_glucose_below_basal(self):
        # The actual physiological claim this pathway exists to
        # demonstrate: enhanced clearance during exercise, then residual
        # insulin sensitisation afterward, should pull glucose below its
        # basal value even with no meal at all.
        bout = ExerciseBout.from_intensity(ExerciseIntensity.VIGOROUS, 45.0)
        composite = couple_with_gi(GlucoseInsulinPathway(), ExercisePathway(bout))
        result = ODESystem(pathway=composite).simulate(duration_s=7200, meals=[])
        assert result.success
        G = result.species["G"]
        assert G.min() < G[0]
