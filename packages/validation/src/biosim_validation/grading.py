# SAVE TO: packages\validation\src\biosim_validation\grading.py
"""
Gold-standard grading (S0, frozen 2026-11-20). Plain-English rules:
docs/GOLD_STANDARD.md.

A grade sits ALONGSIDE the existing PASS/MISS verdict and never changes it.
It asks a stricter question: is the model inside the uncertainty the source
itself reports?

    SET             an initial condition (fasting glucose): not a prediction.
    NOT_ASSESSABLE  no cited number to compare with (UNVERIFIED reference, or
                    a value read off a figure).
    CIRCULAR        the dataset was used to fit or choose a model parameter
                    (CIRCULAR_REGISTER), so agreement proves little.
    GOLD / SILVER / MISS, by the spread the source reports:
      ci_diff  treatment difference with a 95% CI: GOLD inside the CI, else MISS.
      sd / sem single-arm mean (SEM -> SD = SEM * sqrt(n)): GOLD within 0.5 SD,
               SILVER within 1 SD, else MISS.
      none     no spread reported: SILVER if inside the old tolerance band,
               else MISS (never GOLD).

The rules, every benchmark value and the development/locked split are frozen
by benchmark_manifest_sha256() == frozen_manifest.FROZEN_SHA256 (a test).
"""

from __future__ import annotations

import dataclasses
import hashlib
import inspect
import json
import math
from enum import Enum

from .benchmark import BenchmarkCheckResult, BenchmarkDataset, BenchmarkPoint


class Grade(str, Enum):
    GOLD = "GOLD"
    SILVER = "SILVER"
    MISS = "MISS"
    SET = "SET"
    CIRCULAR = "CIRCULAR"
    NOT_ASSESSABLE = "NOT-ASSESSABLE"


# Initial conditions: the model starts there, so matching proves nothing.
SET_POINTS: frozenset[tuple[str, str]] = frozenset({
    ("OGTT_75g_healthy", "fasting_glucose"),
})

# Independence register: datasets used to fit or choose a model parameter.
CIRCULAR_REGISTER: dict[str, str] = {
    "OGTT_75g_healthy": (
        "GutAbsorptionModel.ogtt_calibrated_glucose() (also the default meal) "
        "and first_phase_gain were fitted to this benchmark's points."
    ),
    "T2D_vs_nondiabetic_meal_insulin": (
        "The T2D glucose-effectiveness factor (x0.56) was taken from the same "
        "Basu 2009 cohort (PMID:19196896)."
    ),
}


def grade_point(dataset_name: str, result: BenchmarkCheckResult) -> Grade:
    p: BenchmarkPoint = result.benchmark
    if (dataset_name, p.name) in SET_POINTS:
        return Grade.SET
    if not p.is_scored or p.spread_kind == "figure":
        return Grade.NOT_ASSESSABLE
    if dataset_name in CIRCULAR_REGISTER:
        return Grade.CIRCULAR
    model = result.model_value
    if p.spread_kind == "ci_diff":
        return Grade.GOLD if p.ci_low <= model <= p.ci_high else Grade.MISS
    if p.spread_kind in ("sd", "sem"):
        sd = p.sd if p.spread_kind == "sd" else p.sem * math.sqrt(p.n)
        gap = abs(model - p.reference)
        if gap <= 0.5 * sd:
            return Grade.GOLD
        return Grade.SILVER if gap <= sd else Grade.MISS
    return Grade.SILVER if result.passed else Grade.MISS


def grade_report(report) -> dict[str, Grade]:
    """Point name -> grade for one ValidationReport."""
    return {r.name: grade_point(report.dataset_name, r) for r in report.results}


def development_datasets() -> list[BenchmarkDataset]:
    """Benchmarks that may guide model changes (ModelValidator.validate_all())."""
    from . import benchmark as bm
    return [
        bm.ogtt_75g_healthy_adults(), bm.bergman_iv_glucose(),
        bm.semaglutide_postprandial_glucose(), bm.steroid_t2d_postprandial_glucose(),
        bm.metformin_postprandial_glucose_t2d(),
        bm.canagliflozin_postprandial_glucose_and_insulin(),
        bm.pea_protein_coingestion_50g_glucose(), bm.t2d_vs_nondiabetic_meal_insulin(),
        bm.caffeine_ogtt_insulin(),
    ]


def locked_datasets() -> list[BenchmarkDataset]:
    """The locked test set: never used to choose or adjust anything; reported
    once at the end of each stream (ModelValidator.validate_heldout())."""
    from .heldout import heldout_datasets
    return heldout_datasets()


def _point_record(p: BenchmarkPoint) -> dict:
    d = dataclasses.asdict(p)
    d["tolerance"] = p.tolerance.value
    return d


def benchmark_manifest_sha256() -> str:
    """Hash of every benchmark point (both sets) plus the grading and band code."""
    from . import benchmark as bm
    payload = {
        "dev": {ds.name: [_point_record(p) for p in ds.points] for ds in development_datasets()},
        "locked": {ds.name: [_point_record(p) for p in ds.points] for ds in locked_datasets()},
        "grading_source": inspect.getsource(grade_point),
        "set_points": sorted(SET_POINTS),
        "circular": sorted(CIRCULAR_REGISTER),
        "band_source": inspect.getsource(bm.BenchmarkPoint.tolerance_fraction.fget),
    }
    blob = json.dumps(payload, sort_keys=True, ensure_ascii=True, default=str)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()
