# SAVE TO: packages\data-seeds\tests\test_export_and_exercise.py
"""Tests for population export and exercise scenario builder."""
import json
import pytest
from biosim_seeds.population import PopulationSpec, PopulationRunner
from biosim_seeds.population_export import (
    summary_csv, individuals_csv, histogram_csv, to_json, compare_cohorts_csv,
)
from biosim_seeds.exercise_scenarios import (
    post_exercise_meal, sedentary_meal, fasted_run, compare_exercise_effect,
)
from biosim_physiology import ExerciseIntensity
from biosim_api import SimulationRunner, get_scenario


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def small_summary():
    spec = PopulationSpec(n_individuals=6, random_seed=42,
                          t2d_fraction=0.0, insulin_resistant_fraction=0.0)
    return PopulationRunner(spec).run(get_scenario("standard_meal"))


# ---------------------------------------------------------------------------
# Population export
# ---------------------------------------------------------------------------

class TestSummaryCSV:
    def test_returns_string(self, small_summary):
        assert isinstance(summary_csv(small_summary), str)

    def test_has_header(self, small_summary):
        lines = summary_csv(small_summary).strip().split("\n")
        assert "n_total" in lines[0]

    def test_has_data_row(self, small_summary):
        lines = summary_csv(small_summary).strip().split("\n")
        assert len(lines) == 2   # header + 1 data row

    def test_n_total_in_csv(self, small_summary):
        csv_str = summary_csv(small_summary)
        assert str(small_summary.n_total) in csv_str

    def test_failure_reasons_actually_appear_in_csv(self, monkeypatch):
        # Regression: this export layer wrote n_failed (a count) but
        # never the actual failure_reasons text, added specifically so a
        # real bug failing part of a cohort would be diagnosable (bandit
        # B110 fix). A user exporting to CSV saw THAT individuals failed
        # but never WHY.
        def _raise(self, scenario):
            raise ValueError("synthetic ODE solver failure for testing")

        monkeypatch.setattr(SimulationRunner, "run", _raise)
        spec = PopulationSpec(n_individuals=2, random_seed=42)
        summary = PopulationRunner(spec).run(get_scenario("standard_meal"))
        csv_str = summary_csv(summary)
        assert "synthetic ODE solver failure" in csv_str


class TestIndividualsCSV:
    def test_returns_string(self, small_summary):
        assert isinstance(individuals_csv(small_summary), str)

    def test_row_count_matches_succeeded(self, small_summary):
        lines = [l for l in individuals_csv(small_summary).strip().split("\n") if l]
        assert len(lines) == small_summary.n_succeeded + 1  # +1 for header

    def test_has_age_column(self, small_summary):
        assert "age_years" in individuals_csv(small_summary)

    def test_has_glucose_peak_column(self, small_summary):
        assert "glucose_peak_mM" in individuals_csv(small_summary)


class TestHistogramCSV:
    def test_returns_string(self, small_summary):
        assert isinstance(histogram_csv(small_summary), str)

    def test_has_bins(self, small_summary):
        lines = [l for l in histogram_csv(small_summary).strip().split("\n") if l]
        assert len(lines) > 2   # header + at least some bins

    def test_pct_column_present(self, small_summary):
        assert "pct" in histogram_csv(small_summary)


class TestToJSON:
    def test_returns_string(self, small_summary):
        assert isinstance(to_json(small_summary), str)

    def test_valid_json(self, small_summary):
        data = json.loads(to_json(small_summary))
        assert "summary" in data
        assert "individuals" in data

    def test_individual_count_matches(self, small_summary):
        data = json.loads(to_json(small_summary))
        assert len(data["individuals"]) == small_summary.n_succeeded

    def test_summary_has_mean_peak(self, small_summary):
        data = json.loads(to_json(small_summary))
        assert "glucose_peak_mean_mM" in data["summary"]

    def test_individual_has_sex(self, small_summary):
        data = json.loads(to_json(small_summary))
        if data["individuals"]:
            assert "sex" in data["individuals"][0]

    def test_failure_reasons_actually_appear_in_json(self, monkeypatch):
        # Same regression as TestSummaryCSV.test_failure_reasons_actually_appear_in_csv,
        # for the JSON export path.
        def _raise(self, scenario):
            raise ValueError("synthetic ODE solver failure for testing")

        monkeypatch.setattr(SimulationRunner, "run", _raise)
        spec = PopulationSpec(n_individuals=2, random_seed=42)
        summary = PopulationRunner(spec).run(get_scenario("standard_meal"))
        data = json.loads(to_json(summary))
        assert data["summary"]["n_failed"] == 2
        assert len(data["summary"]["failure_reasons"]) == 2
        assert all("synthetic ODE solver failure" in r for r in data["summary"]["failure_reasons"])


class TestCompareCohortsCSV:
    def test_two_cohorts(self, small_summary):
        csv_str = compare_cohorts_csv([small_summary, small_summary], ["a", "b"])
        lines = [l for l in csv_str.strip().split("\n") if l]
        assert len(lines) == 3   # header + 2 cohorts

    def test_labels_in_csv(self, small_summary):
        csv_str = compare_cohorts_csv([small_summary], ["my_cohort"])
        assert "my_cohort" in csv_str

    def test_sanitizes_formula_injection_in_label(self, small_summary):
        # Found by mandatory security review of a sibling exporter
        # (biosim_expert.batch.CSVExporter's own CSV-formula-injection
        # fix, OWASP CWE-1236) -- compare_cohorts_csv()'s own caller-
        # supplied `labels` argument is bare, unvalidated free text.
        csv_str = compare_cohorts_csv([small_summary], ["=cmd|'/c calc'!A1"])
        import csv, io
        rows = list(csv.DictReader(io.StringIO(csv_str)))
        assert rows[0]["cohort"] == "'=cmd|'/c calc'!A1"


class TestSummaryCSVFormulaInjection:
    # Same fix applied to summary_csv()'s own failure_reasons column.
    # The REAL PopulationRunner failure-capture path always prefixes
    # "age=X sex=Y: " onto every failure_reasons entry (population.py),
    # so a genuinely attacker-controlled leading character can't reach
    # this column through that path today -- this test constructs a
    # CohortSummary directly (not through a real cohort run) to verify
    # the defensive sanitization itself, independent of whether
    # today's one real producer happens to already be safe.
    def test_sanitizes_formula_injection_in_failure_reasons(self):
        from biosim_seeds.population import CohortSummary

        summary = CohortSummary(
            n_total=1, n_succeeded=0,
            glucose_peak_mean_mM=0.0, glucose_peak_sd_mM=0.0,
            glucose_peak_p5_mM=0.0, glucose_peak_p95_mM=0.0,
            insulin_peak_mean_pM=0.0, mtor_peak_mean=0.0, ampk_nadir_mean=0.0,
            pct_above_7p8_mM=0.0, pct_above_11p1_mM=0.0,
            individual_profiles=[], individual_peaks_mM=[],
            failure_reasons=["=cmd|'/c calc'!A1"],
        )
        csv_str = summary_csv(summary)
        import csv, io
        rows = list(csv.DictReader(io.StringIO(csv_str)))
        assert rows[0]["failure_reasons"] == "'=cmd|'/c calc'!A1"


# ---------------------------------------------------------------------------
# Exercise scenarios
# ---------------------------------------------------------------------------

class TestExerciseScenarios:
    def test_post_exercise_meal_construction(self):
        s = post_exercise_meal(75.0, ExerciseIntensity.MODERATE, 45.0, 1.0)
        assert s.name.startswith("post_exercise")
        assert len(s.meal_events) == 1
        assert s.meal_events[0].carbs_g == 75.0

    def test_sedentary_meal_no_drugs(self):
        s = sedentary_meal(60.0)
        assert s.name == "sedentary_meal"
        assert len(s.drug_events) == 0

    def test_fasted_run_no_meals(self):
        s = fasted_run()
        assert len(s.meal_events) == 0

    def test_fasted_run_carries_exercise_bout_in_progress(self):
        # Regression: fasted_run() previously built no exercise-related
        # object at all -- see CLAUDE.md KNOWN BUGS FIXED #16's follow-up.
        s = fasted_run()
        assert s.exercise_bout_in_progress is not None

    def test_fasted_run_lowers_glucose_below_basal(self):
        # fasted_run()'s whole point: exercise (during-bout clearance,
        # then post-bout sensitisation) should measurably lower glucose
        # even with zero meal input -- previously silently absent.
        output = SimulationRunner().run(fasted_run(ExerciseIntensity.VIGOROUS, 45.0))
        G = output.simulation.species["G"]
        assert G.min() < G[0]

    def test_compare_returns_two_scenarios(self):
        ex, sed = compare_exercise_effect()
        assert ex.name != sed.name
        assert ex.meal_events[0].carbs_g == sed.meal_events[0].carbs_g

    def test_sedentary_meal_runs(self):
        output = SimulationRunner().run(sedentary_meal(60.0))
        assert output.success

    def test_fasted_run_runs(self):
        output = SimulationRunner().run(fasted_run())
        assert output.success

    def test_post_exercise_runs(self):
        output = SimulationRunner().run(
            post_exercise_meal(60.0, ExerciseIntensity.LIGHT, 30.0, 2.0)
        )
        assert output.success

    def test_vigorous_intensity_valid(self):
        s = post_exercise_meal(75.0, ExerciseIntensity.VIGOROUS, 30.0, 0.5)
        assert "vigorous" in s.name

    def test_post_exercise_scenario_carries_exercise_bout(self):
        # Regression: post_exercise_meal() used to build an ExerciseBout
        # and reference it only in the free-text description -- never
        # attached to the ScenarioDefinition/PersonalizationSpec, so it
        # had zero effect on the actual simulation.
        s = post_exercise_meal(75.0, ExerciseIntensity.MODERATE, 45.0, 1.0)
        assert s.personalization.exercise_bout is not None

    def test_post_exercise_meal_differs_from_sedentary(self):
        # The actual physiological claim compare_exercise_effect() exists
        # to demonstrate: prior exercise measurably changes the glucose
        # response to the same meal. This was silently false before the
        # exercise_bout-wiring fix (identical output, bit-for-bit).
        exercise, sedentary = compare_exercise_effect()
        out_ex = SimulationRunner().run(exercise)
        out_sed = SimulationRunner().run(sedentary)
        assert out_ex.expert_metrics.glucose_peak_mM != out_sed.expert_metrics.glucose_peak_mM

    def test_post_exercise_meal_lowers_glucose_peak(self):
        # Post-exercise insulin sensitisation should REDUCE the peak, not
        # just change it in some arbitrary direction.
        exercise, sedentary = compare_exercise_effect()
        out_ex = SimulationRunner().run(exercise)
        out_sed = SimulationRunner().run(sedentary)
        assert out_ex.expert_metrics.glucose_peak_mM < out_sed.expert_metrics.glucose_peak_mM
