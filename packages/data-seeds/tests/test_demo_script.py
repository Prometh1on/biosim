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


def test_demo_default_output_is_untracked():
    # Running the README demo must not modify a tracked file: it writes to
    # outputs/ (gitignored), not over the committed README figure.
    demo = runpy.run_path(str(ROOT / "scripts" / "demo.py"))
    default = Path(demo["DEFAULT_OUT"])
    assert default.parent == ROOT / "outputs"
    gitignore = (ROOT / ".gitignore").read_text(encoding="utf-8").splitlines()
    assert "outputs/" in gitignore, "outputs/ must be gitignored"
    assert default != ROOT / "docs" / "figures" / "demo_ogtt_variant.png"
