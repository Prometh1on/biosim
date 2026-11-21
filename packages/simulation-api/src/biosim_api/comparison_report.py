# SAVE TO: packages\simulation-api\src\biosim_api\comparison_report.py
"""
Multi-scenario comparison report.

Runs any number of scenarios and produces a structured comparison
of their metabolic outcomes. Designed for:
    - Comparing different meals (what's the better breakfast?)
    - Comparing the same meal at different times of day
    - Comparing personalised vs base-model predictions
    - Drug intervention comparisons

Usage:
    from biosim_api.comparison_report import ScenarioComparison

    comparison = ScenarioComparison([
        "high_sugar_meal",
        "standard_meal",
        "light_snack",
    ])
    report = comparison.run()
    print(report.as_text())
    csv = report.as_csv()
"""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

# Leading characters a spreadsheet application (Excel/LibreOffice/Google
# Sheets) interprets as the start of a formula when a CSV cell is opened
# -- the classic "CSV/formula injection" vector (OWASP CWE-1236). Found
# during a sibling-exporter security review (biosim_expert.batch
# .CSVExporter's own equivalent fix) -- duplicated here as a small,
# local helper rather than adding a new biosim_api -> biosim_expert
# dependency for a 3-line function (biosim_expert sits above
# biosim_api in this project's own package layering).
_CSV_FORMULA_TRIGGER_CHARS = ("=", "+", "-", "@", "\t", "\r")


def _csv_safe_cell(value: str) -> str:
    """Prefix a single quote onto a string cell that would otherwise be
    interpreted as a formula by a spreadsheet application -- see
    biosim_expert.batch._csv_safe_cell for the full rationale."""
    if isinstance(value, str) and value.startswith(_CSV_FORMULA_TRIGGER_CHARS):
        return "'" + value
    return value


@dataclass
class ScenarioResult:
    """Metrics from one scenario in a comparison."""
    name:               str
    glucose_peak_mM:    float
    glucose_peak_min:   float
    glucose_auc:        float
    insulin_peak_pM:    float
    mtor_peak:          float
    ampk_nadir:         float
    return_time_min:    Optional[float]
    confidence:         float
    success:            bool


@dataclass
class ComparisonReport:
    """Structured comparison of multiple scenario outcomes."""
    scenario_results:   list[ScenarioResult]
    reference_name:     Optional[str]          # first scenario is the reference
    narrative:          str
    ranking_by_glucose: list[str]              # scenario names, low→high peak glucose

    @property
    def best_glucose_scenario(self) -> str:
        """Scenario with the lowest peak glucose."""
        return self.ranking_by_glucose[0]

    @property
    def worst_glucose_scenario(self) -> str:
        """Scenario with the highest peak glucose."""
        return self.ranking_by_glucose[-1]

    def as_text(self) -> str:
        lines = [
            "SCENARIO COMPARISON REPORT",
            "=" * 50,
            "",
            self.narrative,
            "",
            f"{'Scenario':<30} {'Peak G (mM)':>12} {'Time (min)':>11} {'Insulin (pM)':>12} {'Confidence':>11}",
            "-" * 78,
        ]
        for r in sorted(self.scenario_results, key=lambda x: x.glucose_peak_mM):
            flag = " *" if r.name == self.best_glucose_scenario else ""
            lines.append(
                f"{r.name:<30} {r.glucose_peak_mM:>11.2f} {r.glucose_peak_min:>11.0f} "
                f"{r.insulin_peak_pM:>12.0f} {r.confidence:>10.2f}{flag}"
            )
        lines += ["", "* = lowest peak glucose"]
        return "\n".join(lines)

    def as_csv(self) -> str:
        buf = io.StringIO()
        writer = csv.writer(buf)
        writer.writerow([
            "scenario", "glucose_peak_mM", "glucose_peak_min",
            "glucose_auc", "insulin_peak_pM", "mtor_peak",
            "ampk_nadir", "return_time_min", "confidence",
        ])
        for r in self.scenario_results:
            writer.writerow([
                _csv_safe_cell(r.name),
                f"{r.glucose_peak_mM:.3f}",
                f"{r.glucose_peak_min:.1f}",
                f"{r.glucose_auc:.1f}",
                f"{r.insulin_peak_pM:.1f}",
                f"{r.mtor_peak:.4f}",
                f"{r.ampk_nadir:.4f}",
                f"{r.return_time_min:.1f}" if r.return_time_min else "N/A",
                f"{r.confidence:.3f}",
            ])
        return buf.getvalue()

    def delta_vs_reference(self) -> Optional[dict[str, dict]]:
        """
        For each scenario, compute the absolute and % difference in glucose peak
        vs the reference scenario. Returns None if no reference is set.
        """
        if not self.reference_name:
            return None
        ref = next((r for r in self.scenario_results if r.name == self.reference_name), None)
        if ref is None:
            return None
        deltas = {}
        for r in self.scenario_results:
            if r.name == self.reference_name:
                continue
            abs_diff = r.glucose_peak_mM - ref.glucose_peak_mM
            pct_diff = 100.0 * abs_diff / max(ref.glucose_peak_mM, 0.01)
            deltas[r.name] = {
                "abs_diff_mM":  round(abs_diff, 3),
                "pct_diff":     round(pct_diff, 1),
                "direction":    "higher" if abs_diff > 0 else "lower",
            }
        return deltas


class ScenarioComparison:
    """
    Runs multiple scenarios and produces a ComparisonReport.

    Args:
        scenarios: List of scenario names (from library) or ScenarioDefinition objects
        n_output_points: Time resolution for each simulation
    """

    def __init__(
        self,
        scenarios: list,
        n_output_points: int = 300,
    ) -> None:
        self.scenarios       = scenarios
        self.n_output_points = n_output_points

    def run(self, reference_index: int = 0) -> ComparisonReport:
        """
        Run all scenarios and compile the comparison report.

        Args:
            reference_index: Index in scenarios list to use as the reference

        Returns:
            ComparisonReport with all metrics and narrative
        """
        from biosim_api import SimulationRunner, get_scenario, ScenarioDefinition
        runner = SimulationRunner(n_output_points=self.n_output_points)
        results = []

        for s in self.scenarios:
            if isinstance(s, str):
                scenario = get_scenario(s)
            else:
                scenario = s
            try:
                output = runner.run(scenario)
                m = output.expert_metrics
                results.append(ScenarioResult(
                    name=scenario.name,
                    glucose_peak_mM=m.glucose_peak_mM,
                    glucose_peak_min=m.glucose_peak_time_s / 60.0,
                    glucose_auc=m.glucose_auc_uM_s,
                    insulin_peak_pM=m.insulin_peak_pM,
                    mtor_peak=m.mtor_peak,
                    ampk_nadir=m.ampk_nadir,
                    return_time_min=(m.time_to_glucose_return_s / 60.0
                                     if m.time_to_glucose_return_s else None),
                    confidence=m.overall_watchdog_confidence,
                    success=output.success,
                ))
            except Exception as e:
                results.append(ScenarioResult(
                    name=getattr(scenario, "name", str(s)),
                    glucose_peak_mM=0, glucose_peak_min=0, glucose_auc=0,
                    insulin_peak_pM=0, mtor_peak=0, ampk_nadir=1,
                    return_time_min=None, confidence=0, success=False,
                ))

        ranked = sorted(results, key=lambda r: r.glucose_peak_mM)
        ranking = [r.name for r in ranked]

        reference_name = results[reference_index].name if results else None
        narrative = self._build_narrative(results, ranking)

        return ComparisonReport(
            scenario_results=results,
            reference_name=reference_name,
            narrative=narrative,
            ranking_by_glucose=ranking,
        )

    def _build_narrative(
        self,
        results: list[ScenarioResult],
        ranking: list[str],
    ) -> str:
        if not results:
            return "No scenarios completed successfully."
        succeeded = [r for r in results if r.success]
        if not succeeded:
            return "All scenarios failed."

        best  = ranking[0]
        worst = ranking[-1]
        best_r  = next(r for r in succeeded if r.name == best)
        worst_r = next((r for r in succeeded if r.name == worst), best_r)

        n = len(succeeded)
        diff = worst_r.glucose_peak_mM - best_r.glucose_peak_mM

        narrative = (
            f"Compared {n} scenario{'s' if n > 1 else ''}. "
            f"'{best}' produces the lowest peak blood glucose ({best_r.glucose_peak_mM:.1f} mmol/L)."
        )
        if n > 1 and diff > 0.1:
            narrative += (
                f" '{worst}' produces the highest peak ({worst_r.glucose_peak_mM:.1f} mmol/L), "
                f"a difference of {diff:.1f} mmol/L."
            )
        return narrative
