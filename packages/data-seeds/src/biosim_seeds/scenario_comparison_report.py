# SAVE TO: packages\data-seeds\src\biosim_seeds\scenario_comparison_report.py
"""
Scenario comparison report builder for the lay interface.

Generates a structured plain-English comparison between two simulation
outputs, designed for the "what if I had eaten X instead of Y?" feature.

Different from the expert ComparisonReport in biosim_api — this is
written in non-technical language for lay users.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass
class LayComparisonReport:
    """Plain-English comparison between two simulations."""
    scenario_a_name:      str
    scenario_b_name:      str
    glucose_diff_mM:      float        # positive = A is higher
    insulin_diff_pM:      float
    better_scenario:      str          # name of scenario with lower glucose
    headline:             str
    glucose_comparison:   str
    energy_comparison:    str
    recommendation:       str
    disclaimer:           str = (
        "This comparison is based on a mathematical model. "
        "Individual results vary. For personal health guidance, "
        "consult a healthcare professional."
    )

    @property
    def difference_pct(self) -> float:
        """Percentage difference in glucose peak."""
        return abs(self.glucose_diff_mM) / max(
            self._baseline_peak(), 0.01
        ) * 100.0

    def _baseline_peak(self) -> float:
        # Approximate — can't access original output here
        return 8.0

    def as_text(self) -> str:
        lines = [
            f"COMPARISON: {self.scenario_a_name} vs {self.scenario_b_name}",
            "=" * 60,
            "",
            self.headline,
            "",
            "BLOOD SUGAR",
            self.glucose_comparison,
            "",
            "ENERGY",
            self.energy_comparison,
            "",
            "RECOMMENDATION",
            self.recommendation,
            "",
            "-" * 60,
            self.disclaimer,
        ]
        return "\n".join(lines)


def build_lay_comparison(output_a, output_b) -> LayComparisonReport:
    """
    Build a plain-English comparison from two SimulationOutputs.

    Args:
        output_a: First SimulationOutput (the "what I ate" scenario)
        output_b: Second SimulationOutput (the "alternative" scenario)

    Returns:
        LayComparisonReport with all sections populated
    """
    ma = output_a.expert_metrics
    mb = output_b.expert_metrics
    la = output_a.lay_summary
    lb = output_b.lay_summary

    peak_a = ma.glucose_peak_mM
    peak_b = mb.glucose_peak_mM
    diff   = peak_a - peak_b   # positive = A higher

    better = output_b.scenario_name if diff > 0 else output_a.scenario_name
    direction = "lower" if diff > 0 else "higher"
    abs_diff  = abs(diff)

    # Headline
    if abs_diff < 0.3:
        headline = (
            f"'{output_a.scenario_name}' and '{output_b.scenario_name}' "
            f"produce very similar blood sugar responses."
        )
    else:
        headline = (
            f"'{better}' produces a {abs_diff:.1f} mmol/L {direction} "
            f"blood sugar peak than '{output_a.scenario_name if better == output_b.scenario_name else output_b.scenario_name}'."
        )

    # Glucose comparison
    glucose_comparison = (
        f"{output_a.scenario_name}: peaks at {peak_a:.1f} mmol/L "
        f"({la.glucose_narrative.split('.')[0]}).\n"
        f"{output_b.scenario_name}: peaks at {peak_b:.1f} mmol/L "
        f"({lb.glucose_narrative.split('.')[0]})."
    )

    # Energy comparison
    dip_a = "an energy dip" if la.has_energy_dip else "a smooth return to normal energy"
    dip_b = "an energy dip" if lb.has_energy_dip else "a smooth return to normal energy"
    energy_comparison = (
        f"{output_a.scenario_name} is expected to produce {dip_a}.\n"
        f"{output_b.scenario_name} is expected to produce {dip_b}."
    )

    # Recommendation
    if abs_diff < 0.3:
        recommendation = (
            "Both options are similar in their metabolic impact. "
            "Other factors (nutrition, fibre, protein) may be more relevant."
        )
    elif peak_b < 7.8 and peak_a >= 7.8:
        recommendation = (
            f"If you want to keep blood sugar in the normal range (below 7.8 mmol/L), "
            f"'{better}' is the better choice."
        )
    else:
        recommendation = (
            f"'{better}' produces a lower blood sugar response. "
            f"For better metabolic control, it is generally preferable to choose "
            f"lower-carbohydrate or lower-GI options."
        )

    return LayComparisonReport(
        scenario_a_name=output_a.scenario_name,
        scenario_b_name=output_b.scenario_name,
        glucose_diff_mM=diff,
        insulin_diff_pM=ma.insulin_peak_pM - mb.insulin_peak_pM,
        better_scenario=better,
        headline=headline,
        glucose_comparison=glucose_comparison,
        energy_comparison=energy_comparison,
        recommendation=recommendation,
    )
