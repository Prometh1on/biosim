# SAVE TO: packages\validation\tests\test_heldout.py
"""
Held-out validation (finish-plan item B2): the model, unchanged, against
published glucose/insulin data it was never fitted to.

Headline band: STRICT, +/-20% of the published mean. Verdicts use the
existing PASS/MISS/NO DATA logic. The pinned deviations below are the
honest out-of-sample result on 2026-11-16 (commit after cc56fd2/ee1fa77);
if one moves, look at why -- never loosen it to make it pass.
"""
import pytest

from biosim_validation import ModelValidator, Tolerance, ValidationVerdict
from biosim_validation.heldout import (
    CALIBRATION_SOURCE_PMIDS,
    HELDOUT_STRICT_FRACTION,
    heldout_datasets,
)

HELDOUT_NAMES = {
    "HeldOut_Navarro2025_OGTT_75g_healthy",
    "HeldOut_Navarro2025_liquid_meal_healthy",
    "HeldOut_Eshghi2017_OGTT_4h_healthy_vs_T2D",
}


@pytest.fixture(scope="module")
def reports():
    return ModelValidator().validate_heldout()


class TestHeldOutDefinitions:
    def test_three_named_datasets(self):
        assert {d.name for d in heldout_datasets()} == HELDOUT_NAMES

    def test_strict_band_is_20_percent(self):
        assert HELDOUT_STRICT_FRACTION == 0.20
        for ds in heldout_datasets():
            for p in ds.points:
                assert p.tolerance is Tolerance.STRICT
                assert p.tolerance_fraction == pytest.approx(0.20)

    def test_every_point_cites_fetched_text(self):
        for ds in heldout_datasets():
            for p in ds.points:
                assert p.pubmed_id in {"40564902", "28744255"}
                assert "Table" in p.notes and '"' in p.notes, p.name

    def test_held_out_sources_were_never_used_for_calibration(self):
        held_out = {p.pubmed_id for ds in heldout_datasets() for p in ds.points}
        assert held_out.isdisjoint(CALIBRATION_SOURCE_PMIDS)

    def test_fasting_values_are_not_scored(self):
        # Fasting glucose/insulin is the model's starting state, not a prediction.
        for ds in heldout_datasets():
            assert not any(p.name.startswith("fasting") for p in ds.points)


class TestHeldOutResults:
    def test_returns_one_report_per_dataset(self, reports):
        assert set(reports) == HELDOUT_NAMES

    def test_held_out_group_is_separate_from_validate_all(self):
        assert HELDOUT_NAMES.isdisjoint(ModelValidator().validate_all())

    # Current values, after B7 (2026-11-16). These datasets had been looked
    # at (B2) before B7, so these are "seen" results, not fresh held-out ones.
    # The true out-of-sample baseline (before B7: 7 of 12 within +/-20%) is
    # recorded in docs/VALIDATION.md.
    @pytest.mark.parametrize("name, point, deviation", [
        ("HeldOut_Navarro2025_OGTT_75g_healthy", "glucose_30min", -3.4),
        ("HeldOut_Navarro2025_OGTT_75g_healthy", "glucose_60min", -1.5),
        ("HeldOut_Navarro2025_OGTT_75g_healthy", "glucose_90min", -10.4),
        ("HeldOut_Navarro2025_OGTT_75g_healthy", "glucose_120min", -14.1),
        ("HeldOut_Navarro2025_liquid_meal_healthy", "glucose_60min", -10.4),
        ("HeldOut_Navarro2025_liquid_meal_healthy", "glucose_90min", -3.9),
        ("HeldOut_Navarro2025_liquid_meal_healthy", "glucose_120min", -1.1),
        ("HeldOut_Navarro2025_liquid_meal_healthy", "glucose_180min", 4.7),
        ("HeldOut_Eshghi2017_OGTT_4h_healthy_vs_T2D", "glucose_iauc_healthy", 95.4),
        ("HeldOut_Eshghi2017_OGTT_4h_healthy_vs_T2D", "glucose_iauc_t2d", 19.0),
        ("HeldOut_Eshghi2017_OGTT_4h_healthy_vs_T2D", "insulin_iauc_healthy", 1.7),
        ("HeldOut_Eshghi2017_OGTT_4h_healthy_vs_T2D", "insulin_iauc_t2d", -56.1),
    ])
    def test_pinned_out_of_sample_deviation(self, reports, name, point, deviation):
        result = next(r for r in reports[name].results if r.name == point)
        assert result.deviation_percent == pytest.approx(deviation, abs=1.0)

    def test_every_scored_point_is_pinned_and_10_of_12_pass(self, reports):
        results = [r for rep in reports.values() for r in rep.results]
        assert len(results) == 12
        assert sum(r.passed for r in results) == 10

    def test_verdicts(self, reports):
        verdicts = {n: r.verdict for n, r in reports.items()}
        assert verdicts == {
            "HeldOut_Navarro2025_OGTT_75g_healthy": ValidationVerdict.PASS,
            "HeldOut_Navarro2025_liquid_meal_healthy": ValidationVerdict.PASS,
            "HeldOut_Eshghi2017_OGTT_4h_healthy_vs_T2D": ValidationVerdict.MISS,
        }

    def test_reports_say_held_out(self, reports):
        for r in reports.values():
            assert "held-out" in r.module_label.lower()


def test_held_out_sources_are_not_any_seeded_parameter_or_benchmark_source():
    from biosim_seeds.parameter_seed.verification import VERIFICATION_RECORDS
    held_out = {p.pubmed_id for ds in heldout_datasets() for p in ds.points}
    seeded = {r["pmid"] for r in VERIFICATION_RECORDS}
    import biosim_validation.benchmark as bm
    existing = {
        p.pubmed_id
        for name in dir(bm)
        if callable(getattr(bm, name)) and getattr(getattr(bm, name), "__module__", "") == bm.__name__
        and name[0].islower() and not name.startswith("_")
        for ds in [getattr(bm, name)()] if isinstance(ds, bm.BenchmarkDataset)
        for p in ds.points
    }
    assert existing, "found no existing benchmark datasets to compare against"
    assert held_out.isdisjoint(seeded | existing)
