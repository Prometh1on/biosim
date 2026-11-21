# SAVE TO: packages\expert-tools\src\biosim_expert\batch.py
"""
Batch simulation runner and CSV/tabular result exporter.

BatchRunner runs multiple scenarios in parallel and collects results.
CSVExporter converts SimulationOutput and comparison results into
tabular formats for downstream analysis in R or Python.

Parallelism uses ProcessPoolExecutor by default, NOT ThreadPoolExecutor
(the prior implementation, and this docstring, claimed "Parallelism
uses ThreadPoolExecutor — safe because each simulation is independent"
-- true about independence, but the parallelism claim itself was never
verified and turned out to be false). Each simulation is CPU-bound
scipy/numpy ODE solving in pure Python control flow, which holds the
GIL for the overwhelming majority of its wall-clock time -- confirmed
empirically before changing anything (Track 3 roadmap milestone,
"local multiprocessing... start with local multi-core"): benchmarking
16 scenarios via the OLD ThreadPoolExecutor implementation measured
0.66x the sequential runtime (i.e. genuinely SLOWER than not
parallelizing at all -- thread-scheduling/GIL-contention overhead with
zero real concurrency benefit), while ProcessPoolExecutor measured a
genuine 1.55x speedup at a more realistic batch size (n=200, 8
workers, 12 cores available). The crossover is real and scale-
dependent: at very small n (a handful of scenarios), process-spawn
overhead (Windows uses "spawn," not "fork" -- every worker re-imports
this package's full dependency tree) can make ProcessPoolExecutor
slower than sequential too, same as ThreadPoolExecutor was at n=16.
This trade-off is inherent to the workload, not a defect in the
implementation, and matches this project's own "batch"/"population
cohort" use cases (tens to hundreds of scenarios), not tiny 2-3-
scenario comparisons.
"""

from __future__ import annotations

import csv
import io
from concurrent.futures import Executor, ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from typing import Callable, Optional
import numpy as np


def _run_scenario(scenario):
    """
    Module-level, picklable worker entry point for BatchRunner's pool.

    Must be a plain module-level function, not a bound method or
    closure -- ProcessPoolExecutor pickles the submitted callable (and
    its arguments) to send to a worker process, and bound methods of
    arbitrary instances are not reliably picklable across process
    boundaries. Constructs its own SimulationRunner() inside the worker
    rather than pickling a shared instance from the parent process --
    SimulationRunner() construction is cheap (verified: no expensive
    setup in __init__), so this has no meaningful cost, and avoids any
    question of whether SimulationRunner itself is safely picklable/
    shareable across processes.
    """
    from biosim_api import SimulationRunner
    return SimulationRunner().run(scenario)


@dataclass
class BatchResult:
    """Results of a batch simulation run."""
    scenarios:    list[str]
    outputs:      list[object]   # list[SimulationOutput]
    n_succeeded:  int
    n_failed:     int
    failed_names: list[str]
    failure_reasons: list[str] = field(default_factory=list)

    @property
    def success_rate(self) -> float:
        total = self.n_succeeded + self.n_failed
        return self.n_succeeded / total if total > 0 else 0.0

    def get(self, scenario_name: str) -> Optional[object]:
        for o in self.outputs:
            if o.scenario_name == scenario_name:
                return o
        return None


class BatchRunner:
    """
    Runs multiple scenarios in parallel and aggregates results.

    Usage:
        runner = BatchRunner(max_workers=4)
        scenarios = [
            ScenarioDefinition(name=f"meal_{carbs}g",
                               meal_events=[MealEvent(0, carbs)])
            for carbs in [20, 40, 60, 80]
        ]
        results = runner.run_all(scenarios)

    `executor_cls` defaults to ProcessPoolExecutor (see module docstring
    for why) -- exposed as a constructor parameter rather than hardcoded
    so tests that need same-process state visibility (e.g.
    `monkeypatch.setattr(SimulationRunner, "run", ...)` to simulate a
    failure) can pass `ThreadPoolExecutor` explicitly: a genuine,
    documented trade-off of process-based parallelism is that a
    monkeypatch in the parent process is invisible to a spawned child
    process, which re-imports this package fresh rather than sharing
    the parent's patched in-memory state.
    """

    def __init__(
        self,
        max_workers: int = 4,
        executor_cls: type[Executor] = ProcessPoolExecutor,
    ) -> None:
        self.max_workers = max_workers
        self.executor_cls = executor_cls

    def run_all(self, scenarios: list) -> BatchResult:
        """
        Run all scenarios in parallel, collecting results and errors.

        Args:
            scenarios: List of ScenarioDefinition objects

        Returns:
            BatchResult with all outputs, including failures
        """
        outputs = []
        failed = []
        failure_reasons = []

        with self.executor_cls(max_workers=self.max_workers) as executor:
            future_to_name = {
                executor.submit(_run_scenario, scenario): scenario.name
                for scenario in scenarios
            }
            for future in as_completed(future_to_name):
                name = future_to_name[future]
                try:
                    result = future.result()
                    outputs.append(result)
                except Exception as e:
                    failed.append(name)
                    failure_reasons.append(f"{name}: {type(e).__name__}: {e}")

        outputs.sort(key=lambda o: o.scenario_name)

        return BatchResult(
            scenarios=[s.name for s in scenarios],
            outputs=outputs,
            n_succeeded=len(outputs),
            n_failed=len(failed),
            failed_names=failed,
            failure_reasons=failure_reasons,
        )

    def parameter_sweep(
        self,
        base_scenario,
        parameter_name: str,
        values: list[float],
        modifier_fn: Callable,
    ) -> BatchResult:
        """
        Run a parameter sweep — vary one parameter across a range of values.

        Args:
            base_scenario:  Base ScenarioDefinition to modify
            parameter_name: Name of the parameter being swept (for labelling)
            values:         List of parameter values to test
            modifier_fn:    Function(scenario, value) -> ScenarioDefinition

        Returns:
            BatchResult with one output per value
        """
        import dataclasses
        scenarios = []
        for v in values:
            s = modifier_fn(base_scenario, v)
            # Regression: this used to guard the rename with
            # `if not hasattr(s, '_swept')` -- ScenarioDefinition never
            # declares or sets a `_swept` attribute anywhere in this
            # codebase (confirmed via grep), so hasattr() was always
            # False and the guard always took the same branch. Removed
            # the dead condition; the rename it gated always ran anyway.
            s = dataclasses.replace(
                s, name=f"{s.name}_{parameter_name}_{v:.3g}"
            )
            scenarios.append(s)
        return self.run_all(scenarios)


# Leading characters Excel/LibreOffice/Google Sheets interpret as the
# start of a formula when a CSV cell is opened in a spreadsheet
# application (the classic "CSV/formula injection" vector, OWASP
# CWE-1236). Housekeeping-backlog fix: found during MASTER_PLAN 5.2d's
# mandatory code-reviewer pass (rated LOW at the time -- no free-text
# field in this exporter was reachable from a less-trusted input path
# yet), now closed after 5.2c wired IngestionPipeline's own
# extraction-derived text (subject names, failure messages -- genuinely
# less trusted than this platform's own generated strings, since it
# can originate from an LLM extractor or a paper's own free-text
# fields) into CurationQueue, which this exporter's
# export_curation_queue() reads from directly.
_CSV_FORMULA_TRIGGER_CHARS = ("=", "+", "-", "@", "\t", "\r")


def _csv_safe_cell(value: str) -> str:
    """
    Prefix a single quote onto any string cell that would otherwise be
    interpreted as a formula by a spreadsheet application opening this
    CSV -- the standard OWASP-recommended mitigation. A leading quote
    is visible but harmless in every spreadsheet application and in
    any programmatic CSV reader (csv.reader/pandas.read_csv both return
    the quote as a literal character, same as any other cell content).
    Non-string values (numbers, already-formatted strings) are returned
    unchanged -- this only guards free text that could carry a
    leading trigger character from a source this platform doesn't
    fully control (an LLM extractor's output, a paper's own title/
    subject-name text).
    """
    if isinstance(value, str) and value.startswith(_CSV_FORMULA_TRIGGER_CHARS):
        return "'" + value
    return value


class CSVExporter:
    """
    Exports simulation results to CSV format for downstream analysis.
    """

    def export_time_series(
        self,
        output,                   # SimulationOutput
        species: Optional[list[str]] = None,
        include_uncertainty: bool = False,
    ) -> str:
        """
        Export time-series data as CSV.

        Columns: time_s, time_min, G_uM, I_uM, X_s-1, IGF1_uM, mTOR, AMPK
        If uncertainty is available and requested, adds _lower/_upper columns.

        Args:
            output:              SimulationOutput from SimulationRunner
            species:             Subset of species to export (None = all)
            include_uncertainty: Include UQ confidence bands if available

        Returns:
            CSV string
        """
        sim = output.simulation
        if not sim.success:
            return "# Simulation failed — no data to export\n"

        all_species = species or list(sim.species.keys())
        buf = io.StringIO()
        writer = csv.writer(buf)

        # Header
        header = ["time_s", "time_min"]
        for s in all_species:
            unit = sim.species_units.get(s, "unknown")
            header.append(f"{s}_{unit}")

        if include_uncertainty and output.uq is not None:
            for s in all_species:
                if s in output.uq.bands:
                    header.append(f"{s}_lower")
                    header.append(f"{s}_upper")

        writer.writerow(header)

        # Data rows
        times_s = sim.times_s
        times_min = times_s / 60.0
        for i in range(len(times_s)):
            row = [f"{times_s[i]:.2f}", f"{times_min[i]:.4f}"]
            for s in all_species:
                row.append(f"{sim.species[s][i]:.6g}")
            if include_uncertainty and output.uq is not None:
                for s in all_species:
                    if s in output.uq.bands:
                        band = output.uq.bands[s]
                        row.append(f"{band.lower[i]:.6g}")
                        row.append(f"{band.upper[i]:.6g}")
            writer.writerow(row)

        return buf.getvalue()

    def export_metrics_table(self, outputs: list) -> str:
        """
        Export a summary metrics table — one row per scenario.

        Useful for comparing results across a parameter sweep or
        across different personalisation profiles.

        Returns:
            CSV string with columns: scenario, glucose_peak_mM, insulin_peak_pM,
            mtor_peak, ampk_nadir, auc_uM_s, confidence, n_warnings
        """
        buf = io.StringIO()
        writer = csv.writer(buf)
        writer.writerow([
            "scenario",
            "glucose_peak_mM",
            "glucose_peak_time_min",
            "glucose_excursion_mM",
            "insulin_peak_pM",
            "mtor_peak",
            "ampk_nadir",
            "glucose_auc_uM_s",
            "time_to_return_min",
            "watchdog_confidence",
            "n_warnings",
        ])
        for out in outputs:
            m = out.expert_metrics
            return_min = (
                round(m.time_to_glucose_return_s / 60.0, 1)
                if m.time_to_glucose_return_s else "N/A"
            )
            writer.writerow([
                _csv_safe_cell(out.scenario_name),
                f"{m.glucose_peak_mM:.3f}",
                f"{m.glucose_peak_time_s / 60.0:.1f}",
                f"{m.glucose_excursion_uM / 1000.0:.3f}",
                f"{m.insulin_peak_pM:.1f}",
                f"{m.mtor_peak:.4f}",
                f"{m.ampk_nadir:.4f}",
                f"{m.glucose_auc_uM_s:.1f}",
                return_min,
                f"{m.overall_watchdog_confidence:.3f}",
                len(out.warnings),
            ])
        return buf.getvalue()

    def export_curation_queue(self, items: list) -> str:
        """
        Export curation items (e.g. `CurationQueue.pending()`) as CSV --
        MASTER_PLAN 5.2d, the minimal consumer 5.2's own scoping found
        was the last missing piece: real producers now exist (5.2a's
        ParameterStore-conflict bridge, 5.2b's sensitivity-analysis
        bridge) and the queue's own data structure was always real, but
        nothing anywhere presented it anywhere a human curator would
        actually look. A curator can open this CSV in a spreadsheet to
        triage -- no new UI/dashboard needed, reusing this project's
        existing CSV-export precedent (`export_metrics_table`) rather
        than inventing a new presentation mechanism.

        Args:
            items: list of CurationItem (e.g. from
                   `CurationQueue.pending()` -- any subset/ordering the
                   caller wants is fine, this method doesn't query the
                   queue itself)

        Returns:
            CSV string, one row per item -- flagged_at oldest first is
            the caller's own choice (pending() already sorts this way
            by default), not re-sorted here.
        """
        buf = io.StringIO()
        writer = csv.writer(buf)
        writer.writerow([
            "item_id", "item_type", "priority", "subject_name",
            "parameter_kind", "value", "unit", "reason",
            "source_description", "flagged_at", "assigned_to",
        ])
        for item in items:
            writer.writerow([
                _csv_safe_cell(item.item_id),
                item.item_type.value,
                item.priority.name,
                _csv_safe_cell(item.subject_name),
                _csv_safe_cell(item.parameter_kind),
                item.value,
                _csv_safe_cell(item.unit),
                _csv_safe_cell(item.reason),
                _csv_safe_cell(item.source_description),
                item.flagged_at.isoformat(),
                _csv_safe_cell(item.assigned_to or ""),
            ])
        return buf.getvalue()

    def export_comparison(self, comparison) -> str:
        """
        Export a comparison result — differential time series.

        Columns: time_s, G_a, G_b, G_diff, I_a, I_b, I_diff, mTOR_a, mTOR_b, mTOR_diff
        """
        buf = io.StringIO()
        writer = csv.writer(buf)
        writer.writerow([
            "time_s", "time_min",
            f"G_{comparison.scenario_a.scenario_name}_uM",
            f"G_{comparison.scenario_b.scenario_name}_uM",
            "G_diff_uM",
            f"I_{comparison.scenario_a.scenario_name}_uM",
            f"I_{comparison.scenario_b.scenario_name}_uM",
            "I_diff_uM",
            "mTOR_diff",
        ])
        times_s = comparison.times_s
        G_a = comparison.scenario_a.simulation.species["G"]
        G_b_interp = np.interp(
            times_s,
            comparison.scenario_b.simulation.times_s,
            comparison.scenario_b.simulation.species["G"],
        )
        I_a = comparison.scenario_a.simulation.species["I"]
        I_b_interp = np.interp(
            times_s,
            comparison.scenario_b.simulation.times_s,
            comparison.scenario_b.simulation.species["I"],
        )
        for i in range(len(times_s)):
            writer.writerow([
                f"{times_s[i]:.2f}",
                f"{times_s[i]/60.0:.4f}",
                f"{G_a[i]:.4f}",
                f"{G_b_interp[i]:.4f}",
                f"{comparison.glucose_diff_uM[i]:.4f}",
                f"{I_a[i]:.6g}",
                f"{I_b_interp[i]:.6g}",
                f"{comparison.insulin_diff_uM[i]:.6g}",
                f"{comparison.mtor_diff[i]:.6g}",
            ])
        return buf.getvalue()
