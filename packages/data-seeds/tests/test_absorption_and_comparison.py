# SAVE TO: packages\data-seeds\tests\test_absorption_and_comparison.py
"""Tests for bi-exponential gut absorption model and scenario comparison report."""
import math
import pytest
from biosim_ode.gut_absorption import (
    GutAbsorptionModel, Ra_biexp, Ra_curve,
)
from biosim_api.comparison_report import ScenarioComparison, ComparisonReport


# ---------------------------------------------------------------------------
# Gut absorption model
# ---------------------------------------------------------------------------

class TestGutAbsorptionModel:
    def test_default_construction(self):
        m = GutAbsorptionModel()
        assert m.k_emptying < m.k_absorption

    def test_ra_zero_at_t0(self):
        m = GutAbsorptionModel()
        assert m.Ra(0.0, 1.0) == pytest.approx(0.0)

    def test_ra_positive_after_t0(self):
        m = GutAbsorptionModel()
        assert m.Ra(5.0, 1.0) > 0

    def test_ra_peaks_then_falls(self):
        m = GutAbsorptionModel()
        dose = 1.0
        t_peak = m.peak_time_min
        ra_before = m.Ra(t_peak - 2, dose)
        ra_at     = m.Ra(t_peak,     dose)
        ra_after  = m.Ra(t_peak + 20, dose)
        assert ra_at >= ra_before
        assert ra_at >= ra_after

    def test_ra_zero_negative_time(self):
        assert GutAbsorptionModel().Ra(-1.0, 1.0) == 0.0

    def test_ra_zero_negative_dose(self):
        assert GutAbsorptionModel().Ra(5.0, 0.0) == 0.0

    def test_liquid_glucose_faster_than_mixed(self):
        liquid = GutAbsorptionModel.liquid_glucose()
        mixed  = GutAbsorptionModel.mixed_meal()
        # Liquid glucose peaks earlier
        assert liquid.peak_time_min < mixed.peak_time_min

    def test_mixed_faster_than_high_fat(self):
        mixed    = GutAbsorptionModel.mixed_meal()
        high_fat = GutAbsorptionModel.high_fat_meal()
        assert mixed.peak_time_min < high_fat.peak_time_min

    def test_for_gi_high_faster_than_low(self):
        high = GutAbsorptionModel.for_gi(90)
        low  = GutAbsorptionModel.for_gi(30)
        assert high.peak_time_min < low.peak_time_min
        assert high.k_emptying > low.k_emptying

    def test_for_gi_reference_value(self):
        # Range updated alongside the GI=100 anchor fix -- see
        # test_for_gi_100_matches_ogtt_calibrated_glucose_exactly below.
        # k_emptying now spans ~0.0083 (GI=30) to 0.0249 (GI=100); GI=55
        # sits roughly a third of the way up that range.
        m = GutAbsorptionModel.for_gi(55)
        assert 0.005 < m.k_emptying < 0.030

    def test_ogtt_calibrated_glucose_valid_construction(self):
        m = GutAbsorptionModel.ogtt_calibrated_glucose()
        assert m.k_emptying < m.k_absorption
        assert m.Ra(0.0, 1.0) == pytest.approx(0.0)

    def test_for_gi_100_matches_ogtt_calibrated_glucose_exactly(self):
        # Regression: for_gi(100)'s rate constants used to be a factor of
        # 3.6x/7.6x faster than ogtt_calibrated_glucose() -- a
        # citation-independent internal-consistency bug (the immediately
        # adjacent, real-data-validated preset that MealEvent routes to at
        # GI>=99.5), fixed by anchoring for_gi()'s GI=100 endpoint to it
        # directly. This replaces the old
        # test_ogtt_calibrated_glucose_much_slower_than_for_gi_100, whose
        # entire premise (for_gi(100) being wrong) is what this fix
        # resolves -- see for_gi()'s docstring for what's fixed (the
        # anchor) vs what remains an open, un-validated question (the
        # GI-to-timing slope shape).
        calibrated = GutAbsorptionModel.ogtt_calibrated_glucose()
        generic_gi_100 = GutAbsorptionModel.for_gi(100)
        assert generic_gi_100.k_emptying == pytest.approx(calibrated.k_emptying)
        assert generic_gi_100.k_absorption == pytest.approx(calibrated.k_absorption)

    def test_for_gi_low_end_now_slower_than_the_validated_anchor(self):
        # This assertion previously ran the OPPOSITE direction (see git
        # history / CLAUDE.md): for_gi()'s entire range was confirmed
        # faster than the one validated point, including its own
        # "slowest" GI=30 endpoint. The anchor fix above (for_gi(100) now
        # equals the validated ogtt_calibrated_glucose() exactly) flips
        # this by construction -- since for_gi(30)'s k_emptying is defined
        # as 1/3 of for_gi(100)'s, it is now necessarily slower than the
        # validated point, not faster. This is the flip the old test's own
        # comment anticipated and asked to be updated, not deleted, when
        # it happened.
        calibrated = GutAbsorptionModel.ogtt_calibrated_glucose()
        slowest_generic = GutAbsorptionModel.for_gi(30)
        assert slowest_generic.peak_time_min > calibrated.peak_time_min

    def test_invalid_ke_raises(self):
        with pytest.raises(ValueError):
            GutAbsorptionModel(k_emptying=0.0)

    def test_ka_le_ke_raises(self):
        with pytest.raises(ValueError):
            GutAbsorptionModel(k_emptying=0.1, k_absorption=0.05)

    def test_peak_time_formula(self):
        m = GutAbsorptionModel(k_emptying=0.05, k_absorption=0.15)
        expected = math.log(0.15 / 0.05) / (0.15 - 0.05)
        assert m.peak_time_min == pytest.approx(expected)


class TestRaBiexp:
    def test_convenience_function(self):
        ra = Ra_biexp(10.0, 1.0)
        assert ra > 0

    def test_with_custom_model(self):
        model = GutAbsorptionModel.liquid_glucose()
        ra = Ra_biexp(5.0, 1.0, model=model)
        assert ra > 0

    def test_ra_curve_returns_two_lists(self):
        times, values = Ra_curve(240.0, 1.0, n_points=50)
        assert len(times) == 50
        assert len(values) == 50

    def test_ra_curve_starts_near_zero(self):
        times, values = Ra_curve(240.0, 1.0)
        assert values[0] == pytest.approx(0.0, abs=0.001)


# ---------------------------------------------------------------------------
# Scenario comparison report
# ---------------------------------------------------------------------------

class TestScenarioComparison:
    @pytest.fixture(scope="class")
    @classmethod
    def report(cls):
        comp = ScenarioComparison(["high_sugar_meal", "standard_meal", "light_snack"])
        return comp.run()

    def test_returns_report(self, report):
        assert isinstance(report, ComparisonReport)

    def test_three_results(self, report):
        assert len(report.scenario_results) == 3

    def test_all_succeeded(self, report):
        assert all(r.success for r in report.scenario_results)

    def test_ranking_lowest_to_highest(self, report):
        peaks = [r.glucose_peak_mM for r in
                 sorted(report.scenario_results, key=lambda r: r.glucose_peak_mM)]
        assert peaks == sorted(peaks)

    def test_best_is_light_snack(self, report):
        # Light snack should have lowest glucose peak
        assert "light_snack" in report.best_glucose_scenario

    def test_worst_is_high_sugar(self, report):
        assert "high_sugar" in report.worst_glucose_scenario

    def test_narrative_populated(self, report):
        assert len(report.narrative) > 20

    def test_as_text_contains_scenario_names(self, report):
        text = report.as_text()
        assert "high_sugar_meal" in text
        assert "standard_meal" in text

    def test_as_csv_parseable(self, report):
        import csv, io
        csv_str = report.as_csv()
        reader = list(csv.DictReader(io.StringIO(csv_str)))
        assert len(reader) == 3
        assert "glucose_peak_mM" in reader[0]

    def test_delta_vs_reference(self, report):
        deltas = report.delta_vs_reference()
        assert deltas is not None
        assert len(deltas) == 2   # reference excluded

    def test_delta_directions(self, report):
        deltas = report.delta_vs_reference()
        # light_snack should be lower than high_sugar_meal (reference)
        if "light_snack" in deltas:
            assert deltas["light_snack"]["direction"] == "lower"

    def test_single_scenario(self):
        comp = ScenarioComparison(["fasted"])
        report = comp.run()
        assert len(report.scenario_results) == 1

    def test_two_scenario_comparison(self):
        from biosim_seeds import morning_meal, evening_meal
        comp = ScenarioComparison([morning_meal(), evening_meal()])
        report = comp.run()
        assert len(report.scenario_results) == 2
        # Evening should have higher glucose
        morning_r = next(r for r in report.scenario_results if "morning" in r.name)
        evening_r = next(r for r in report.scenario_results if "evening" in r.name)
        assert evening_r.glucose_peak_mM > morning_r.glucose_peak_mM


class TestComparisonReportCSVInjection:
    # Found by mandatory security review of a sibling exporter
    # (biosim_expert.batch.CSVExporter's own CSV-formula-injection fix)
    # -- ComparisonReport.as_csv() wrote r.name as a bare, unsanitized
    # CSV cell. A real scenario name can't itself start with a
    # formula-trigger character (validated by biosim_api.scenario's own
    # rules), but ScenarioResult.name is a plain str field with no
    # constructor-level validation, so this is constructed directly
    # rather than through a real ScenarioComparison run.
    @pytest.fixture(scope="class")
    @classmethod
    def report(cls):
        comp = ScenarioComparison(["high_sugar_meal", "standard_meal", "light_snack"])
        return comp.run()

    def test_as_csv_sanitizes_malicious_scenario_name(self):
        from biosim_api.comparison_report import ComparisonReport, ScenarioResult

        report = ComparisonReport(
            scenario_results=[
                ScenarioResult(
                    name="=cmd|'/c calc'!A1",
                    glucose_peak_mM=8.0, glucose_peak_min=45.0,
                    glucose_auc=1000.0, insulin_peak_pM=300.0,
                    mtor_peak=0.5, ampk_nadir=0.5,
                    return_time_min=60.0, confidence=0.9, success=True,
                ),
            ],
            reference_name=None,
            narrative="",
            ranking_by_glucose=["=cmd|'/c calc'!A1"],
        )
        import csv, io
        rows = list(csv.DictReader(io.StringIO(report.as_csv())))
        assert rows[0]["scenario"] == "'=cmd|'/c calc'!A1"

    def test_as_csv_does_not_sanitize_ordinary_scenario_names(self, report):
        # No-false-positives guard, against a real report.
        import csv, io
        rows = list(csv.DictReader(io.StringIO(report.as_csv())))
        names = {r["scenario"] for r in rows}
        assert names == {"high_sugar_meal", "standard_meal", "light_snack"}
