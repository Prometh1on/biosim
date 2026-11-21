# SAVE TO: packages\validation\tests\test_grading.py
"""
Gold-standard grading (S0, 2026-11-20; rules in docs/GOLD_STANDARD.md).

A grade is computed alongside the existing PASS/MISS verdict and never
changes it. The rules and every benchmark value are frozen by a hash, so a
change to either needs a deliberate, recorded "benchmark change".
"""
import math
from pathlib import Path

import pytest

from biosim_validation import BenchmarkPoint, ModelValidator, Tolerance
from biosim_validation.grading import (
    CIRCULAR_REGISTER,
    Grade,
    benchmark_manifest_sha256,
    development_datasets,
    grade_point,
    locked_datasets,
)
from biosim_validation.frozen_manifest import FROZEN_SHA256

ROOT = Path(__file__).resolve().parents[3]


def _point(**kw):
    base = dict(name="p", description="d", reference=10.0, unit="u",
                tolerance=Tolerance.WIDE, pubmed_id="123")
    base.update(kw)
    return BenchmarkPoint(**base)


class TestGradeRules:
    def test_ci_difference_gold_inside_ci(self):
        p = _point(spread_kind="ci_diff", ci_low=5.0, ci_high=20.0)
        assert grade_point("D", p.check(6.0)) is Grade.GOLD

    def test_ci_difference_miss_outside_ci(self):
        p = _point(spread_kind="ci_diff", ci_low=8.0, ci_high=12.0)
        assert grade_point("D", p.check(13.0)) is Grade.MISS

    def test_sd_gold_within_half_sd_silver_within_one(self):
        p = _point(spread_kind="sd", sd=2.0, n=20)
        assert grade_point("D", p.check(10.9)) is Grade.GOLD
        assert grade_point("D", p.check(11.9)) is Grade.SILVER
        assert grade_point("D", p.check(12.5)) is Grade.MISS

    def test_sem_is_converted_to_sd(self):
        p = _point(spread_kind="sem", sem=1.0, n=4)          # SD = 2.0
        assert grade_point("D", p.check(11.0)) is Grade.GOLD
        assert grade_point("D", p.check(11.5)) is Grade.SILVER

    def test_no_spread_is_capped_at_silver(self):
        p = _point(spread_kind="none")
        assert grade_point("D", p.check(10.0)) is Grade.SILVER   # exact match, still SILVER
        assert grade_point("D", p.check(20.0)) is Grade.MISS

    def test_figure_readout_and_unverified_are_not_assessable(self):
        assert grade_point("D", _point(spread_kind="figure").check(10.0)) is Grade.NOT_ASSESSABLE
        assert grade_point("D", _point(pubmed_id="UNVERIFIED").check(10.0)) is Grade.NOT_ASSESSABLE

    def test_circular_dataset_is_never_gold(self):
        dataset = next(iter(CIRCULAR_REGISTER))
        p = _point(spread_kind="ci_diff", ci_low=5.0, ci_high=20.0)
        assert grade_point(dataset, p.check(10.0)) is Grade.CIRCULAR

    def test_fasting_glucose_is_set_not_scored(self):
        p = _point(name="fasting_glucose", spread_kind="sd", sd=1.0)
        assert grade_point("OGTT_75g_healthy", p.check(10.0)) is Grade.SET

    def test_grade_never_changes_the_verdict(self):
        before = {n: r.verdict for n, r in ModelValidator().validate_all().items()}
        reports = ModelValidator().validate_all()
        for name, rep in reports.items():
            for r in rep.results:
                grade_point(name, r)
        assert {n: r.verdict for n, r in reports.items()} == before


class TestSpreadSources:
    def test_every_spread_value_has_a_quote_and_location(self):
        for ds in development_datasets() + locked_datasets():
            for p in ds.points:
                if p.spread_kind in ("sd", "sem", "ci_diff"):
                    assert p.quote.strip() and p.source_location.strip(), (ds.name, p.name)
                    assert p.n, (ds.name, p.name)

    def test_drug_ci_values_match_the_abstracts(self):
        pts = {p.name: p for ds in development_datasets() for p in ds.points}
        sema = pts["glucose_iauc_5h_reduction_mmol_h_per_L"]
        assert (sema.ci_low, sema.ci_high, sema.n) == (0.27, 2.42, 30)
        assert "-1.34 mmol h/L [-2.42, -0.27]" in sema.quote
        met = pts["glucose_iauc_4h_reduction_mmol_min_per_L_t2d"]
        assert (met.ci_low, met.ci_high, met.n) == (64.0, 307.0, 15)
        assert "MD 186 mmol/L × min [95% CI 64-307]" in met.quote


class TestDevelopmentAndLockedSets:
    def test_sets_are_disjoint_and_labelled(self):
        dev = {ds.name for ds in development_datasets()}
        locked = {ds.name for ds in locked_datasets()}
        assert dev.isdisjoint(locked)
        for ds in development_datasets():
            assert all(p.benchmark_set == "dev" for p in ds.points), ds.name
        for ds in locked_datasets():
            assert all(p.benchmark_set == "locked" for p in ds.points), ds.name

    def test_b2_held_out_sets_are_locked(self):
        locked = {ds.name for ds in locked_datasets()}
        assert {
            "HeldOut_Navarro2025_OGTT_75g_healthy",
            "HeldOut_Navarro2025_liquid_meal_healthy",
            "HeldOut_Eshghi2017_OGTT_4h_healthy_vs_T2D",
        } <= locked

    def test_sets_match_what_the_validator_runs(self):
        # Review (S0): a dataset added to validate_all() must also be graded
        # and frozen, and the locked set must be exactly validate_heldout().
        assert {ds.name for ds in development_datasets()} == set(ModelValidator().validate_all())
        assert {ds.name for ds in locked_datasets()} == set(ModelValidator().validate_heldout())

    def test_every_spread_has_the_fields_its_rule_needs(self):
        for ds in development_datasets() + locked_datasets():
            for p in ds.points:
                if p.spread_kind == "ci_diff":
                    assert p.ci_low is not None and p.ci_high is not None and p.ci_low < p.ci_high
                elif p.spread_kind == "sd":
                    assert p.sd and p.sd > 0
                elif p.spread_kind == "sem":
                    assert p.sem and p.sem > 0 and p.n and p.n > 1
                else:
                    assert p.spread_kind in ("none", "figure"), (ds.name, p.name)

    def test_development_validation_never_runs_locked_sets(self):
        locked = {ds.name for ds in locked_datasets()}
        assert locked.isdisjoint(ModelValidator().validate_all())


class TestFrozen:
    def test_benchmarks_and_grading_code_match_the_frozen_hash(self):
        # Changing a reference, a spread, a set label or the grading code
        # changes this hash. Update FROZEN_SHA256 only together with a
        # HISTORY.md entry "Benchmark change: <new hash>" giving the reason.
        assert benchmark_manifest_sha256() == FROZEN_SHA256

    def test_frozen_hash_is_recorded_in_history(self):
        history = (ROOT / "HISTORY.md").read_text(encoding="utf-8")
        assert f"Benchmark freeze: {FROZEN_SHA256}" in history or \
               f"Benchmark change: {FROZEN_SHA256}" in history
