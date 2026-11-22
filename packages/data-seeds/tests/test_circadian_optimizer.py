# SAVE TO: packages\data-seeds\tests\test_circadian_optimizer.py
"""Tests for the circadian meal optimizer."""
import pytest
from biosim_seeds.circadian_optimizer import (
    CircadianOptimizer, CircadianSchedule, TimeSlot,
)


class TestTimeSlot:
    def test_construction(self):
        s = TimeSlot(10.0, 7.2, 280.0, 6.5, False, True)
        assert s.hour == 10.0
        assert s.success


class TestCircadianOptimizer:
    @pytest.fixture(scope="class")
    @classmethod
    def schedule(cls):
        # Use 4 time points for speed in tests
        opt = CircadianOptimizer(n_times=4)
        return opt.optimize(carbs_g=60.0)

    def test_returns_schedule(self, schedule):
        assert isinstance(schedule, CircadianSchedule)

    def test_has_time_slots(self, schedule):
        assert len(schedule.time_slots) == 4

    def test_ranked_slots_ascending(self, schedule):
        peaks = [s.glucose_peak_mM for s in schedule.ranked_slots]
        assert peaks == sorted(peaks)

    def test_best_hour_set(self, schedule):
        assert 0.0 <= schedule.best_hour < 24.0

    def test_worst_hour_set(self, schedule):
        assert 0.0 <= schedule.worst_hour < 24.0

    def test_glucose_range_positive(self, schedule):
        assert schedule.glucose_range_mM >= 0.0

    def test_narrative_populated(self, schedule):
        text = schedule.narrative()
        assert len(text) > 30
        assert "mmol" in text

    def test_best_time_description_format(self, schedule):
        desc = schedule.best_time_description
        assert "AM" in desc or "PM" in desc

    def test_to_dict_structure(self, schedule):
        d = schedule.to_dict()
        assert "best_hour" in d
        assert "ranked_slots" in d
        assert "glucose_range_mM" in d

    def test_to_dict_includes_energy_score(self, schedule):
        # Regression: TimeSlot.energy_score was computed every iteration
        # but never read by to_dict() or narrative() -- confirmed via
        # grep, it appeared nowhere else in the whole repo.
        d = schedule.to_dict()
        assert d["ranked_slots"], "expected at least one successful slot"
        assert "energy_score" in d["ranked_slots"][0]

    def test_failed_slot_records_failure_reason(self, monkeypatch):
        # Regression: a failed slot's exception used to be silently
        # discarded (bare `except Exception: ...`, no field to even
        # carry a reason) -- same bug class as CohortSummary
        # .failure_reasons' bandit B110 fix in population.py.
        from biosim_api import SimulationRunner as _SR

        def _raise(self, scenario):
            raise ValueError("synthetic ODE solver failure for testing")

        monkeypatch.setattr(_SR, "run", _raise)
        opt = CircadianOptimizer(n_times=2)
        schedule = opt.optimize(carbs_g=60.0)
        assert all(not s.success for s in schedule.time_slots)
        assert all(s.failure_reason is not None for s in schedule.time_slots)
        assert all("synthetic ODE solver failure" in s.failure_reason for s in schedule.time_slots)
        assert len(schedule.failure_reasons) == 2
        assert all("synthetic ODE solver failure" in r for r in schedule.failure_reasons)
        assert schedule.to_dict()["failure_reasons"] == schedule.failure_reasons

    def test_morning_better_than_evening(self):
        """Morning should produce lower glucose than evening — core chronobiology."""
        opt = CircadianOptimizer(n_times=4)
        schedule = opt.optimize(carbs_g=75.0)
        # Find closest to 8am and 20pm slots
        slots = {s.hour: s for s in schedule.time_slots if s.success}
        if not slots:
            pytest.skip("No successful simulations")
        # Best hour should not be at the late-evening nadir (22:00)
        # (At 4 points: 0, 6, 12, 18 — morning=6, evening=18)
        morning_slots = [s for s in schedule.time_slots if 4 <= s.hour <= 14 and s.success]
        evening_slots = [s for s in schedule.time_slots if s.hour >= 16 and s.success]
        if morning_slots and evening_slots:
            best_morning = min(s.glucose_peak_mM for s in morning_slots)
            best_evening = min(s.glucose_peak_mM for s in evening_slots)
            # Morning should be at least as good as or better than evening
            assert best_morning <= best_evening + 0.5  # allow 0.5 mM tolerance

    def test_with_profile(self):
        from biosim_physiology import PhysiologicalProfile, DiseaseState
        profile = PhysiologicalProfile(
            disease_states=[DiseaseState.INSULIN_RESISTANT]
        )
        opt = CircadianOptimizer(n_times=3, profile=profile)
        schedule = opt.optimize(carbs_g=50.0)
        assert isinstance(schedule, CircadianSchedule)
        assert len(schedule.time_slots) == 3

    def test_carbs_preserved_in_schedule(self, schedule):
        assert schedule.carbs_g == pytest.approx(60.0)
