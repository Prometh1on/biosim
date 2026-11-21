# SAVE TO: packages\validation\src\biosim_validation\frozen_manifest.py
"""
Frozen hash of every benchmark value, the development/locked split and the
grading code (grading.benchmark_manifest_sha256()), set in S0 on 2026-11-20
BEFORE any model change. test_grading.py fails if they change.

To change a benchmark or a grading rule on purpose: update this hash AND add
a HISTORY.md line "Benchmark change: <new hash>" with the reason. Never do
it to make a result pass.
"""

FROZEN_SHA256 = "098e8cbbdfd7bec2d73d499f433702d9b246fb0ca2f3e6fd20fb50da5e1e3d8f"
