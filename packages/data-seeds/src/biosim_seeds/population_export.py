# SAVE TO: packages\data-seeds\src\biosim_seeds\population_export.py
"""
Population simulation result export.

Converts CohortSummary results into CSV and JSON formats for
downstream analysis in R, Python, or Excel.

Outputs:
    summary_csv()       -- one-row aggregate statistics table
    individuals_csv()   -- one row per individual with profile + peak
    histogram_csv()     -- glucose peak distribution (for plotting)
    to_json()           -- full structured JSON for API responses
"""

from __future__ import annotations

import csv
import io
import json
from typing import Optional

# Leading characters a spreadsheet application (Excel/LibreOffice/Google
# Sheets) interprets as the start of a formula when a CSV cell is opened
# -- the classic "CSV/formula injection" vector (OWASP CWE-1236). Found
# during a sibling-exporter security review (biosim_expert.batch
# .CSVExporter's own equivalent fix) -- duplicated here as a small,
# local helper rather than adding a new biosim_seeds -> biosim_expert
# dependency for a 3-line function. summary_csv()'s failure_reasons
# text is built from real exception messages (f"{name}: {type(e)
# .__name__}: {e}") -- a lower-probability but real vector if an
# exception's own message ever echoes attacker-influenced content.
_CSV_FORMULA_TRIGGER_CHARS = ("=", "+", "-", "@", "\t", "\r")


def _csv_safe_cell(value: str) -> str:
    """Prefix a single quote onto a string cell that would otherwise be
    interpreted as a formula by a spreadsheet application -- see
    biosim_expert.batch._csv_safe_cell for the full rationale."""
    if isinstance(value, str) and value.startswith(_CSV_FORMULA_TRIGGER_CHARS):
        return "'" + value
    return value


def summary_csv(summary) -> str:
    """
    Export aggregate cohort statistics as a single-row CSV.

    Columns: n_total, n_succeeded, mean_peak_mM, sd_peak_mM,
             p5_peak_mM, p95_peak_mM, pct_above_7p8, pct_above_11p1,
             mean_insulin_pM, mean_mtor_peak, mean_ampk_nadir
    """
    # Regression: this used to write n_failed (a count) but never the
    # actual failure_reasons text -- CohortSummary.failure_reasons was
    # added specifically so a real bug failing part of a cohort run
    # would be diagnosable (see the bandit B110 fix in CLAUDE.md), but
    # this export layer, this package's own stated reporting layer for
    # CohortSummary, silently dropped it. A user exporting to CSV saw
    # THAT individuals failed but never WHY.
    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow([
        "n_total", "n_succeeded", "n_failed",
        "mean_peak_glucose_mM", "sd_peak_glucose_mM",
        "p5_peak_glucose_mM", "p95_peak_glucose_mM",
        "pct_above_7p8_mM", "pct_above_11p1_mM",
        "mean_insulin_peak_pM", "mean_mtor_peak", "mean_ampk_nadir",
        "failure_reasons",
    ])
    writer.writerow([
        summary.n_total, summary.n_succeeded, summary.n_failed,
        f"{summary.glucose_peak_mean_mM:.3f}",
        f"{summary.glucose_peak_sd_mM:.3f}",
        f"{summary.glucose_peak_p5_mM:.3f}",
        f"{summary.glucose_peak_p95_mM:.3f}",
        f"{summary.pct_above_7p8_mM:.1f}",
        f"{summary.pct_above_11p1_mM:.1f}",
        f"{summary.insulin_peak_mean_pM:.1f}",
        f"{summary.mtor_peak_mean:.4f}",
        f"{summary.ampk_nadir_mean:.4f}",
        _csv_safe_cell(" | ".join(summary.failure_reasons)),
    ])
    return buf.getvalue()


def individuals_csv(summary) -> str:
    """
    Export per-individual results as CSV.

    Columns: individual_id, age_years, sex, bmi, fitness,
             disease_states, glucose_peak_mM
    """
    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow([
        "individual_id", "age_years", "sex", "bmi", "fitness",
        "disease_states", "glucose_peak_mM",
    ])
    for i, (profile, peak) in enumerate(
        zip(summary.individual_profiles, summary.individual_peaks_mM)
    ):
        writer.writerow([
            i + 1,
            f"{profile.age_years:.1f}",
            profile.biological_sex.value,
            f"{profile.bmi:.1f}",
            profile.fitness_level.value,
            "+".join(d.value for d in profile.disease_states),
            f"{peak:.3f}",
        ])
    return buf.getvalue()


def histogram_csv(summary, n_bins: int = 20) -> str:
    """
    Export glucose peak distribution as histogram data (bin_min, bin_max, count).
    Suitable for bar charts in Excel or R.
    """
    import numpy as np
    if not summary.individual_peaks_mM:
        return "bin_min,bin_max,count\n"

    peaks  = summary.individual_peaks_mM
    counts, edges = np.histogram(peaks, bins=n_bins)
    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow(["bin_min_mM", "bin_max_mM", "count", "pct"])
    total = sum(counts)
    for i, count in enumerate(counts):
        writer.writerow([
            f"{edges[i]:.3f}",
            f"{edges[i+1]:.3f}",
            int(count),
            f"{100.0 * count / max(total, 1):.1f}",
        ])
    return buf.getvalue()


def to_json(summary, indent: int = 2) -> str:
    """
    Export full CohortSummary as JSON.
    Includes summary statistics and per-individual data.
    """
    individuals = []
    for profile, peak in zip(summary.individual_profiles, summary.individual_peaks_mM):
        individuals.append({
            "age_years":     round(profile.age_years, 1),
            "sex":           profile.biological_sex.value,
            "bmi":           round(profile.bmi, 1),
            "fitness":       profile.fitness_level.value,
            "disease_states": [d.value for d in profile.disease_states],
            "glucose_peak_mM": round(peak, 3),
        })

    data = {
        "summary": {
            "n_total":              summary.n_total,
            "n_succeeded":          summary.n_succeeded,
            "n_failed":             summary.n_failed,
            "failure_reasons":      list(summary.failure_reasons),
            "glucose_peak_mean_mM": round(summary.glucose_peak_mean_mM, 3),
            "glucose_peak_sd_mM":   round(summary.glucose_peak_sd_mM, 3),
            "glucose_peak_p5_mM":   round(summary.glucose_peak_p5_mM, 3),
            "glucose_peak_p95_mM":  round(summary.glucose_peak_p95_mM, 3),
            "pct_above_7p8_mM":     round(summary.pct_above_7p8_mM, 1),
            "pct_above_11p1_mM":    round(summary.pct_above_11p1_mM, 1),
            "mean_insulin_peak_pM": round(summary.insulin_peak_mean_pM, 1),
            "mtor_peak_mean":       round(summary.mtor_peak_mean, 4),
            "ampk_nadir_mean":      round(summary.ampk_nadir_mean, 4),
        },
        "individuals": individuals,
    }
    return json.dumps(data, indent=indent)


def compare_cohorts_csv(summaries: list, labels: Optional[list[str]] = None) -> str:
    """
    Compare multiple CohortSummary results side by side.
    One row per cohort, columns are summary statistics.
    """
    labels = labels or [f"cohort_{i}" for i in range(len(summaries))]
    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow([
        "cohort", "n_total", "n_succeeded",
        "mean_peak_mM", "sd_peak_mM", "p5_mM", "p95_mM",
        "pct_above_7p8", "pct_above_11p1",
    ])
    for label, s in zip(labels, summaries):
        writer.writerow([
            _csv_safe_cell(label), s.n_total, s.n_succeeded,
            f"{s.glucose_peak_mean_mM:.3f}", f"{s.glucose_peak_sd_mM:.3f}",
            f"{s.glucose_peak_p5_mM:.3f}",  f"{s.glucose_peak_p95_mM:.3f}",
            f"{s.pct_above_7p8_mM:.1f}",    f"{s.pct_above_11p1_mM:.1f}",
        ])
    return buf.getvalue()
