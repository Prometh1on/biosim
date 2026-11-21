r"""
Regenerate docs/figures/heldout_validation.png: the unchanged model against the
three held-out datasets (packages/validation/src/biosim_validation/heldout.py).

Run from the repo root:  .venv\Scripts\python.exe scripts\plot_heldout.py
Error bars are the published spread: SD for Navarro-Masip 2025, SEM for
Eshghi 2017 (as each paper's table footnote states).
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import conftest  # noqa: E402,F401  (adds packages/*/src to sys.path)

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from biosim_api import MealEvent, ScenarioDefinition, SimulationRunner, ogtt_75g  # noqa: E402
from biosim_validation import ModelValidator  # noqa: E402
from biosim_validation.heldout import MG_DL_PER_MMOL_L as MG  # noqa: E402


def curve(output):
    s = output.simulation
    return s.times_s / 60.0, s.species["G"] / 1000.0


def main() -> None:
    runner = SimulationRunner()
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

    # A: OGTT, Navarro-Masip 2025 Table 2 (mean ± SD, mg/dL -> mM)
    t, g = curve(runner.run(ogtt_75g()))
    ax = axes[0]
    ax.plot(t, g, color="#1f5f99", label="Model (unchanged)")
    ax.errorbar([0, 30, 60, 90, 120],
                np.array([94.1, 167.5, 187.5, 176.0, 150.4]) / MG,
                yerr=np.array([7.6, 31.0, 44.9, 53.2, 46.0]) / MG,
                fmt="o", color="#b2432f", capsize=3, label="Navarro-Masip 2025, n=24 (±SD)")
    ax.set(title="A. 75 g OGTT, healthy", xlabel="minutes", ylabel="glucose (mmol/L)", xlim=(0, 180))

    # B: liquid meal, same table
    t, g = curve(runner.run(ScenarioDefinition(
        name="fortisip", duration_s=14400,
        meal_events=[MealEvent(time_s=0.0, carbs_g=18.4, protein_g=6.0)])))
    ax = axes[1]
    ax.plot(t, g, color="#1f5f99", label="Model (unchanged)")
    ax.errorbar([0, 60, 90, 120, 180],
                np.array([102.4, 129.4, 115.2, 104.3, 90.0]) / MG,
                yerr=np.array([10.2, 27.4, 28.5, 23.8, 10.8]) / MG,
                fmt="o", color="#b2432f", capsize=3, label="Navarro-Masip 2025, n=24 (±SD)")
    ax.set(title="B. Liquid meal (18.4 g CHO, 6 g protein)", xlabel="minutes",
           ylabel="glucose (mmol/L)", xlim=(0, 240))

    # C: Eshghi 2017 Table 5 iAUCs, model / published ratio
    report = ModelValidator().validate_heldout()["HeldOut_Eshghi2017_OGTT_4h_healthy_vs_T2D"]
    sem = {"glucose_iauc_healthy": 138, "glucose_iauc_t2d": 193,
           "insulin_iauc_healthy": 14486, "insulin_iauc_t2d": 5078}
    ax = axes[2]
    names = [r.name.replace("_iauc_", "\n") for r in report.results]
    ratios = [r.model_value / r.benchmark.reference for r in report.results]
    errs = [sem[r.name] / r.benchmark.reference for r in report.results]
    x = np.arange(len(ratios))
    ax.bar(x, ratios, color="#1f5f99", label="Model / published")
    ax.errorbar(x, np.ones_like(x, float), yerr=errs, fmt="_", color="#b2432f", capsize=4,
                label="Published (=1, ±SEM)")
    ax.axhspan(0.8, 1.2, color="#999999", alpha=0.18, label="±20% band")
    ax.set_xticks(x, names)
    ax.set(title="C. 4-h OGTT iAUC, Eshghi 2017 (n=6+6)", ylabel="model / published")

    for ax in axes:
        ax.legend(fontsize=8)
    fig.suptitle("Held-out validation: model unchanged vs data it was never fitted to")
    fig.tight_layout()
    out = ROOT / "docs" / "figures" / "heldout_validation.png"
    fig.savefig(out, dpi=130)
    print(out)


if __name__ == "__main__":
    main()
