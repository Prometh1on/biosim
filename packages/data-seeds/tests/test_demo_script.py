# SAVE TO: packages\data-seeds\tests\test_demo_script.py
"""The README's one-command demo (scripts/demo.py) runs and writes its figure."""
import runpy
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]


def test_demo_writes_figure_and_variant_raises_glucose(tmp_path):
    demo = runpy.run_path(str(ROOT / "scripts" / "demo.py"))
    out = tmp_path / "demo.png"
    result = demo["main"](out)
    assert out.exists() and out.stat().st_size > 10_000
    # TCF7L2 rs7903146 TT: less insulin secretion, so a higher glucose peak.
    assert result["variant_peak_mM"] > result["baseline_peak_mM"]
    assert result["variant_warnings"] == []   # no rsID or raw allele in warnings
