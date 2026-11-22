# SAVE TO: packages\data-seeds\tests\test_pathway_extensions.py
"""Tests for NEFA pathway, molecular timeline, and advanced extensions."""
import pytest
from biosim_ode.nefa_pathway import NEFAPathway, NEFAParameters
from biosim_ode.base_pathway import BasePathway, PathwayMaturity
from biosim_expert.timeline import TimelineExtractor, MolecularTimeline, TimelineEvent
from biosim_api import SimulationRunner, high_sugar_meal_scenario, fasted_scenario


# ---------------------------------------------------------------------------
# NEFA pathway
# ---------------------------------------------------------------------------

class TestNEFAPathway:
    def test_instantiation(self):
        assert NEFAPathway().name == "nefa_v1"

    def test_species_names(self):
        assert "NEFA" in NEFAPathway().species_names

    def test_y0_at_basal(self):
        p = NEFAPathway()
        assert p.y0[0] == pytest.approx(350.0)

    def test_dydt_one_value(self):
        p = NEFAPathway()
        assert len(p.dydt(0.0, p.y0)) == 1

    def test_high_insulin_suppresses_nefa(self):
        p = NEFAPathway()
        # Start at max, with high insulin → NEFA should fall
        d_high_i = p.dydt(0.0, [p.params.nefa_max], signals={"insulin_pM": 400.0})[0]
        assert d_high_i < 0

    def test_no_insulin_elevates_nefa(self):
        p = NEFAPathway()
        # Start at basal, with no insulin → NEFA should rise
        d_no_i = p.dydt(0.0, [p.params.nefa_basal], signals={"insulin_pM": 0.0})[0]
        assert d_no_i > 0

    def test_ir_modifier_below_one_at_elevated_nefa(self):
        p = NEFAPathway()
        signals = p.output_signals(0.0, [700.0])   # elevated NEFA
        assert signals["insulin_resistance_modifier"] < 1.0

    def test_ir_modifier_near_one_at_basal_nefa(self):
        p = NEFAPathway()
        signals = p.output_signals(0.0, [350.0])   # basal NEFA
        # At basal NEFA, effect should be small
        assert signals["insulin_resistance_modifier"] > 0.7

    def test_ir_modifier_exactly_neutral_at_reference(self):
        # Regression: modifier must be exactly 1.0 at the fixed population
        # reference (not just "near one"), otherwise composing NEFA with
        # GlucoseInsulinPathway silently double-counts insulin resistance
        # already baked into the Bergman 1981 p3 calibration.
        from biosim_ode.nefa_pathway import NEFA_IR_REFERENCE
        p = NEFAPathway()
        signals = p.output_signals(0.0, [NEFA_IR_REFERENCE])
        assert signals["insulin_resistance_modifier"] == pytest.approx(1.0)

    def test_ir_modifier_uses_fixed_reference_not_instance_basal(self):
        # Obesity preset has an elevated *instance* basal NEFA -- the modifier
        # must still reflect reduced insulin sensitivity relative to the fixed
        # population reference, not be renormalized to the obese basal.
        obese = NEFAPathway.for_obesity()
        signals = obese.output_signals(0.0, [obese.params.nefa_basal])
        assert signals["insulin_resistance_modifier"] < 1.0

    def test_y0_is_genuine_fixed_point(self):
        # Regression for the nefa_target calibration bug: y0 must be an
        # actual fixed point of dydt at basal insulin, not just a
        # plausible-looking literal. This was previously false for every
        # preset -- e.g. the default preset's NEFA silently drifted from
        # 350 uM toward ~540 uM over ~600 minutes even with no disease/
        # exercise modifiers applied. Structural root cause: the
        # un-normalized nefa_target was a convex combination of nefa_basal
        # and nefa_max that could only equal nefa_basal in the unreachable
        # limit of infinite insulin -- fixed by normalizing antilipolysis
        # against its own value at the reference basal insulin.
        assert NEFAPathway().validate_y0(tolerance=1e-8) == []
        assert NEFAPathway.for_obesity().validate_y0(tolerance=1e-8) == []
        assert NEFAPathway.for_exercise().validate_y0(tolerance=1e-8) == []

    def test_nefa_stable_over_long_trajectory_at_basal_insulin(self):
        # The practical consequence of the bug above: standalone integration
        # at basal insulin (66 pM) must not let NEFA drift over time, not
        # just at t=0 -- for every preset, since each defines its own basal.
        from scipy.integrate import solve_ivp

        for preset in (NEFAPathway(), NEFAPathway.for_obesity(), NEFAPathway.for_exercise()):
            def dydt_at_basal_insulin(t, y, preset=preset):
                return preset.dydt(t, y, signals={"insulin_pM": 66.0})

            sol = solve_ivp(dydt_at_basal_insulin, [0, 600], preset.y0, method="Radau", rtol=1e-8, atol=1e-10)
            assert sol.success
            assert sol.y[0, -1] == pytest.approx(preset.params.nefa_basal, abs=1e-6)

    def test_obesity_params_elevated_basal(self):
        normal = NEFAPathway()
        obese  = NEFAPathway.for_obesity()
        assert obese.params.nefa_basal > normal.params.nefa_basal

    def test_implements_base_pathway(self):
        assert isinstance(NEFAPathway(), BasePathway)

    def test_maturity_experimental(self):
        assert NEFAPathway().maturity == PathwayMaturity.EXPERIMENTAL

    def test_output_signals_keys(self):
        p = NEFAPathway()
        signals = p.output_signals(0.0, p.y0)
        assert "nefa_uM" in signals
        assert "insulin_resistance_modifier" in signals

    def test_composite_with_glucose_insulin(self):
        from biosim_ode import GlucoseInsulinPathway
        from biosim_ode.coupling import CompositePathway
        composite = CompositePathway([GlucoseInsulinPathway(), NEFAPathway()])
        assert len(composite.y0) == 7
        derivs = composite.dydt(0.0, composite.y0)
        assert len(derivs) == 7

    def test_couple_with_gi_default_preserves_basal_steady_state(self):
        # Default (non-obesity) NEFAParameters is neutral at its own basal
        # (matches NEFA_IR_REFERENCE), so composing it must not perturb G/I.
        from biosim_ode import GlucoseInsulinPathway
        from biosim_ode.nefa_pathway import couple_with_gi
        composite = couple_with_gi(GlucoseInsulinPathway(), NEFAPathway())
        derivs = composite.dydt(0.0, composite.y0)
        assert derivs[0] == pytest.approx(0.0, abs=1e-6)
        assert derivs[1] == pytest.approx(0.0, abs=1e-6)

    def test_couple_with_gi_exercise_preset_reduces_insulin_action(self):
        # for_exercise()'s basal (800 uM) is itself well above the fixed
        # NEFA_IR_REFERENCE (350 uM), so composing it correctly reduces
        # effective p3 from t=0 -- the intended dynamic value-add for
        # exercise scenarios (elevated NEFA competing with glucose uptake).
        from biosim_ode import GlucoseInsulinPathway
        from biosim_ode.nefa_pathway import couple_with_gi
        gi = GlucoseInsulinPathway()
        composite = couple_with_gi(gi, NEFAPathway.for_exercise())
        y = composite.y0
        y[1] = gi.params.Ib + 50.0   # elevated insulin -> nonzero dX response

        dX_composed = composite.dydt(0.0, y)[2]
        dX_unmodified = gi.dydt(0.0, y[:6])[2]
        assert dX_composed < dX_unmodified


# ---------------------------------------------------------------------------
# Molecular timeline
# ---------------------------------------------------------------------------

class TestMolecularTimeline:
    @pytest.fixture(scope="class")
    @classmethod
    def meal_output(cls):
        return SimulationRunner().run(high_sugar_meal_scenario())

    @pytest.fixture(scope="class")
    @classmethod
    def timeline(cls, meal_output):
        return TimelineExtractor().extract(meal_output)

    def test_returns_timeline(self, meal_output):
        tl = TimelineExtractor().extract(meal_output)
        assert isinstance(tl, MolecularTimeline)

    def test_has_events(self, timeline):
        assert len(timeline.events) > 0

    def test_events_sorted_by_time(self, timeline):
        times = [e.time_min for e in timeline.events]
        assert times == sorted(times)

    def test_glucose_peak_event_present(self, timeline):
        types = [e.event_type for e in timeline.events]
        assert "glucose_peak" in types

    def test_insulin_peak_event_present(self, timeline):
        types = [e.event_type for e in timeline.events]
        assert "insulin_peak" in types

    def test_mtor_peak_event_present(self, timeline):
        types = [e.event_type for e in timeline.events]
        assert "mtor_peak" in types

    def test_ampk_nadir_event_present(self, timeline):
        types = [e.event_type for e in timeline.events]
        assert "ampk_nadir" in types

    def test_to_text_returns_string(self, timeline):
        text = timeline.to_text()
        assert isinstance(text, str)
        assert len(text) > 50
        assert "glucose_peak" in text or "mmol" in text

    def test_to_dict_structure(self, timeline):
        d = timeline.to_dict()
        assert "scenario" in d
        assert "events" in d
        assert isinstance(d["events"], list)

    def test_events_in_window(self, timeline):
        early = timeline.events_in_window(0, 60)
        assert all(0 <= e.time_min <= 60 for e in early)

    def test_events_of_type(self, timeline):
        peaks = timeline.events_of_type("glucose_peak")
        assert len(peaks) == 1

    def test_scenario_name_preserved(self, timeline, meal_output):
        assert timeline.scenario_name == meal_output.scenario_name

    def test_fasted_timeline_no_glucose_peak_above_baseline(self):
        output = SimulationRunner().run(fasted_scenario(3600))
        tl = TimelineExtractor().extract(output)
        peak_events = tl.events_of_type("glucose_peak")
        if peak_events:
            # Fasted: glucose should stay near baseline
            assert peak_events[0].values.get("G_mM", 0) < 7.0

    def test_significance_notable_for_high_glucose(self, timeline):
        peak = timeline.events_of_type("glucose_peak")[0]
        # High sugar meal should produce a notable peak
        g_peak = peak.values.get("G_mM", 0)
        if g_peak > 7.8:
            assert peak.significance in ("notable", "warning")
