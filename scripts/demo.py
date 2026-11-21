r"""
One-command demo: a 75 g oral glucose tolerance test (OGTT), with and without
one genetic variant.

    .venv\Scripts\python.exe scripts\demo.py          (Windows)
    .venv/bin/python scripts/demo.py                  (macOS/Linux)

It simulates the OGTT for a default healthy adult, then again for someone
carrying two copies of the TCF7L2 rs7903146 T allele (the model lowers insulin
secretion by 12% per copy -- a modelling choice, see effect_db.py), and saves
docs/figures/demo_ogtt_variant.png. All data are simulated; no personal data.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import conftest  # noqa: E402,F401  (adds packages/*/src to sys.path)

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from biosim_api import PersonalizationSpec, SimulationRunner, ogtt_75g  # noqa: E402
from biosim_genotype import Genotype, Variant  # noqa: E402
from biosim_genotype.effect_db import VariantEffectDB  # noqa: E402

DEFAULT_OUT = ROOT / "docs" / "figures" / "demo_ogtt_variant.png"


def main(out: Path = DEFAULT_OUT) -> dict:
    runner = SimulationRunner()
    baseline = runner.run(ogtt_75g())
    variant = runner.run(ogtt_75g(personalization=PersonalizationSpec(
        genetic_variants=[Variant("rs7903146", Genotype.HOMOZYGOUS_ALT)],
    )))

    # Two copies x the model's per-allele effect (read, not hardcoded).
    per_allele = VariantEffectDB.default().get("rs7903146")[0].effect_per_allele
    variant_label = f"TCF7L2 rs7903146 TT (secretion {2 * per_allele:+.0%})"

    fig, (ax_g, ax_i) = plt.subplots(1, 2, figsize=(11, 4.2))
    for output, label, color in (
        (baseline, "No variant", "#1f5f99"),
        (variant, variant_label, "#b2432f"),
    ):
        sim = output.simulation
        minutes = sim.times_s / 60.0
        ax_g.plot(minutes, sim.species["G"] / 1000.0, color=color, label=label)
        ax_i.plot(minutes, sim.species["I"] * 1e6, color=color, label=label)
    ax_g.set(title="Blood glucose", xlabel="minutes after drinking 75 g glucose",
             ylabel="mmol/L", xlim=(0, 240))
    ax_i.set(title="Plasma insulin", xlabel="minutes after drinking 75 g glucose",
             ylabel="pmol/L", xlim=(0, 240))
    for ax in (ax_g, ax_i):
        ax.legend(fontsize=8)
    fig.suptitle("Biosim demo: simulated 75 g OGTT, with and without one genetic variant")
    fig.tight_layout()
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=130)
    plt.close(fig)

    result = {
        "baseline_peak_mM": baseline.expert_metrics.glucose_peak_mM,
        "variant_peak_mM": variant.expert_metrics.glucose_peak_mM,
        # Privacy rule 7: no rsID may appear in any output text.
        "variant_warnings": [w for w in variant.warnings if re.search(r"rs\d+", w)],
    }
    print(f"Peak glucose: {result['baseline_peak_mM']:.1f} mmol/L without the variant, "
          f"{result['variant_peak_mM']:.1f} mmol/L with it. Figure: {out}")
    return result


if __name__ == "__main__":
    main()
