# SAVE TO: packages\expert-tools\tests\test_timeline.py
"""
Test suite for biosim_expert.timeline.

TimelineExtractor had zero test coverage before this file -- found while
auditing the module for the "computed but never wired" bug class already
fixed repeatedly elsewhere in this project this session: `extract()`
assigned `m = output.expert_metrics` and never read it, while
hand-rolling a second, independent return-to-baseline computation in
event 7 instead of reusing `m.time_to_glucose_return_s` (already
computed via the more general SimulationResult.time_to_return_to_baseline()).
"""

import pytest

from biosim_api import SimulationRunner, high_sugar_meal_scenario, fasted_scenario
from biosim_expert import MolecularTimeline, TimelineEvent, TimelineExtractor


@pytest.fixture(scope="module")
def meal_output():
    return SimulationRunner().run(high_sugar_meal_scenario())


@pytest.fixture(scope="module")
def timeline(meal_output):
    return TimelineExtractor().extract(meal_output)


class TestTimelineExtractor:
    def test_returns_molecular_timeline(self, timeline):
        assert isinstance(timeline, MolecularTimeline)

    def test_events_present_for_a_meal_scenario(self, timeline):
        types = {e.event_type for e in timeline.events}
        assert "meal_start" in types
        assert "glucose_peak" in types
        assert "insulin_peak" in types
        assert "mtor_peak" in types
        assert "ampk_nadir" in types
        assert "simulation_end" in types

    def test_fasted_scenario_emits_simulation_start_not_meal_start(self):
        output = SimulationRunner().run(fasted_scenario(duration_s=3600.0))
        timeline = TimelineExtractor().extract(output)
        types = {e.event_type for e in timeline.events}
        assert "simulation_start" in types
        assert "meal_start" not in types

    def test_failed_simulation_returns_empty_timeline(self):
        class _FakeSim:
            success = False
            species = {}

        class _FakeOutput:
            simulation = _FakeSim()
            scenario_name = "broken"
            expert_metrics = None

        timeline = TimelineExtractor().extract(_FakeOutput())
        assert timeline.events == []
        assert timeline.duration_min == 0

    def test_glucose_return_event_matches_expert_metrics(self, meal_output, timeline):
        # Regression: event 7 used to hand-roll its own return-to-baseline
        # scan instead of reading the already-computed
        # ExpertMetrics.time_to_glucose_return_s -- confirm it now matches
        # exactly (same authoritative computation, not a second one that
        # could silently drift out of sync).
        return_events = timeline.events_of_type("glucose_return_baseline")
        expected_s = meal_output.expert_metrics.time_to_glucose_return_s
        if expected_s is None:
            assert return_events == []
        else:
            assert len(return_events) == 1
            assert return_events[0].time_min == pytest.approx(expected_s / 60.0)

    def test_no_glucose_return_event_when_never_returns(self):
        # A meal scenario short enough that glucose is still rising at
        # the end of the window (never returns to baseline within it)
        # should not emit event 7 at all.
        from biosim_api import ScenarioDefinition, MealEvent
        scenario = ScenarioDefinition(
            name="short", duration_s=600.0,
            meal_events=[MealEvent(time_s=0.0, carbs_g=75.0)],
        )
        output = SimulationRunner().run(scenario)
        assert output.expert_metrics.time_to_glucose_return_s is None
        timeline = TimelineExtractor().extract(output)
        assert timeline.events_of_type("glucose_return_baseline") == []


class TestMolecularTimeline:
    def test_events_in_window_filters_by_time(self):
        events = [
            TimelineEvent(time_min=0.0, event_type="a", description="a"),
            TimelineEvent(time_min=50.0, event_type="b", description="b"),
            TimelineEvent(time_min=200.0, event_type="c", description="c"),
        ]
        tl = MolecularTimeline(events=events, duration_min=240.0, scenario_name="x")
        windowed = tl.events_in_window(10.0, 100.0)
        assert [e.event_type for e in windowed] == ["b"]

    def test_events_of_type_filters_correctly(self):
        events = [
            TimelineEvent(time_min=0.0, event_type="a", description="a"),
            TimelineEvent(time_min=1.0, event_type="b", description="b"),
            TimelineEvent(time_min=2.0, event_type="a", description="a again"),
        ]
        tl = MolecularTimeline(events=events, duration_min=10.0, scenario_name="x")
        assert len(tl.events_of_type("a")) == 2
        assert len(tl.events_of_type("b")) == 1

    def test_to_text_contains_scenario_name_and_events(self, timeline):
        text = timeline.to_text()
        assert timeline.scenario_name in text
        for e in timeline.events:
            assert e.description in text

    def test_to_dict_round_trips_event_fields(self, timeline):
        d = timeline.to_dict()
        assert d["scenario"] == timeline.scenario_name
        assert len(d["events"]) == len(timeline.events)
        assert d["events"][0]["event_type"] == timeline.events[0].event_type
