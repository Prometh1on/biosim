# SAVE TO: packages\validation\tests\test_first_contact.py
"""
H2 first-contact held-out test (pre-registered in docs/VALIDATION.md section 8,
commit d804d64; model frozen at 85bb29e; run once on 2026-11-21).

The deviations pinned below ARE the first-contact result: 1 of 11 points within
+/-20%. If one moves, find out why and report it as a post-hoc change -- never
loosen a pin, and never tune the model against these data.
"""
import pytest

from biosim_validation import ModelValidator, Tolerance, ValidationVerdict
from biosim_validation.first_contact import (
    FIRST_CONTACT_PMIDS,
    first_contact_datasets,
)
from biosim_validation.grading import Grade, grade_report, locked_datasets
from biosim_validation.heldout import CALIBRATION_SOURCE_PMIDS, heldout_datasets

# (dataset, point) -> (deviation % at the single pre-registered run, grade)
FIRST_CONTACT_RESULT = {
    ("FirstContact_Wickenberg2010_OGTT_75g_healthy", "glucose_iauc_0_120"): (67.7, Grade.MISS),
    ("FirstContact_Wickenberg2010_OGTT_75g_healthy", "insulin_iauc_0_30"): (62.4, Grade.MISS),
    ("FirstContact_Wickenberg2010_OGTT_75g_healthy", "insulin_iauc_0_120"): (115.7, Grade.MISS),
    ("FirstContact_Josic2010_bread_meal_healthy", "glucose_iauc_0_30"): (-14.8, Grade.GOLD),
    ("FirstContact_Josic2010_bread_meal_healthy", "glucose_iauc_0_120"): (65.4, Grade.MISS),
    ("FirstContact_Josic2010_bread_meal_healthy", "insulin_iauc_0_30"): (43.5, Grade.MISS),
    ("FirstContact_Josic2010_bread_meal_healthy", "insulin_iauc_0_120"): (66.8, Grade.MISS),
    ("FirstContact_Yoshida2022_OGTT_75g_young_adults", "glucose_auc_young_male"): (34.3, Grade.MISS),
    ("FirstContact_Yoshida2022_OGTT_75g_young_adults", "insulin_auc_young_male"): (21.4, Grade.GOLD),
    ("FirstContact_Yoshida2022_OGTT_75g_young_adults", "glucose_auc_young_female"): (63.5, Grade.MISS),
    ("FirstContact_Yoshida2022_OGTT_75g_young_adults", "insulin_auc_young_female"): (22.4, Grade.GOLD),
}


@pytest.fixture(scope="module")
def reports():
    return ModelValidator().validate_first_contact()


class TestPreRegisteredDefinitions:
    def test_eleven_points_in_three_datasets(self):
        keys = {(d.name, p.name) for d in first_contact_datasets() for p in d.points}
        assert keys == set(FIRST_CONTACT_RESULT)

    def test_every_point_is_strict_cited_and_tabulated(self):
        for ds in first_contact_datasets():
            for p in ds.points:
                assert p.tolerance is Tolerance.STRICT
                assert p.pubmed_id in FIRST_CONTACT_PMIDS
                assert p.is_scored
                assert p.spread_kind in ("sd", "sem")
                assert "Table" in p.source_location
                assert p.quote

    def test_sources_were_never_used_before(self):
        from biosim_validation.grading import development_datasets
        used = set(CALIBRATION_SOURCE_PMIDS)
        for ds in development_datasets() + heldout_datasets():
            used |= {p.pubmed_id for p in ds.points}
        assert not (FIRST_CONTACT_PMIDS & used)

    def test_not_added_to_the_frozen_locked_manifest(self):
        frozen = {ds.name for ds in locked_datasets()}
        assert not frozen & {ds.name for ds in first_contact_datasets()}


class TestFirstContactResult:
    def test_every_point_matches_the_single_run(self, reports):
        for (dataset, point), (deviation, _) in FIRST_CONTACT_RESULT.items():
            result = next(r for r in reports[dataset].results if r.name == point)
            assert result.deviation_percent == pytest.approx(deviation, abs=0.5), (dataset, point)

    def test_one_of_eleven_within_20_percent(self, reports):
        passed = [r.name for rep in reports.values() for r in rep.results if r.passed]
        assert passed == ["glucose_iauc_0_30"]

    def test_every_dataset_verdict_is_miss(self, reports):
        assert {rep.verdict for rep in reports.values()} == {ValidationVerdict.MISS}

    def test_frozen_grades(self, reports):
        for (dataset, point), (_, grade) in FIRST_CONTACT_RESULT.items():
            assert grade_report(reports[dataset])[point] is grade, (dataset, point)
