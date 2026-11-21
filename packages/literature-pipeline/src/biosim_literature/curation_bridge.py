# SAVE TO: packages\literature-pipeline\src\biosim_literature\curation_bridge.py
"""
Extraction-to-curation bridge (MASTER_PLAN 5.2c).

Closes the third and last of the real gaps found while scoping
`CurationQueue`'s own workflow (MASTER_PLAN 5.2): `IngestionPipeline`
already correctly flags every problematic `ExtractionOutcome` with
`needs_curation=True` (bad unit, bad kind, an implausible value, an
unresolved entity name, or an unexpected build error) -- but nothing
anywhere ever surfaced those flagged outcomes to the curation queue for
a human to actually review. This module is the missing wiring,
mirroring `biosim_parameter_db.curation_bridge`'s own shape (a
standalone bridge function callers invoke explicitly, not automatic/
hooked into `IngestionPipeline.ingest_paper()` -- extraction and
curation-queue surfacing are different concerns, and forcing every
ingestion call to also touch a `CurationQueue` would couple two
packages that don't otherwise depend on each other).

Flow:
    1. IngestionPipeline runs and produces an IngestionReport, each
       ExtractionOutcome already flagged needs_curation as appropriate
       (already existing, unaffected by this module).
    2. This bridge reads report.outcomes, filters to the flagged ones,
       and enqueues each as a CurationQueue extraction item.
    3. Curators review and resolve items from the queue (existing
       CurationQueue.record_decision()).

Per MASTER_PLAN.md row 5.2c: technically unblocked today (this bridge
works with any Extractor, including the free MockExtractor -- it does
not require a real ANTHROPIC_API_KEY), but of limited real value until
genuine LLM-based extraction is actually running in production, since
the project's "keep dev free" decision (5.3) means the pipeline itself
sees little real production use today. Built anyway for completeness,
matching this campaign's own "small, well-scoped chunk" sizing for
this row.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from biosim_curation import CurationItemType, CurationPriority, CurationQueue

from .pipeline import IngestionOutcome, IngestionReport

# Priority per outcome type -- UNRESOLVED/IMPLAUSIBLE are genuine
# scientific-judgment calls a curator needs to make (a real entity or
# value that just doesn't fit this platform's current registry/
# plausibility bounds); BAD_UNIT/BAD_KIND/ERROR are more often a
# pipeline/extractor formatting problem than a scientific one, still
# worth a human look but rated one notch lower by default.
#
# Found by mandatory code review, disclosed rather than re-engineered:
# IMPLAUSIBLE actually covers TWO distinct plausibility.py sub-cases --
# a PHYSICAL_LIMITS FAIL (e.g. a negative Km; that check's own message
# text says "almost certainly a data entry error", the same character
# as the MEDIUM bucket) and a PLAUSIBLE_RANGES WARNING (a value outside
# the expected range but not physically impossible, genuinely a
# judgment call). Both map to this single IngestionOutcome.IMPLAUSIBLE
# with no sub-status carried through, so both get HIGH priority here --
# a defensible, conservative choice (the curator still sees the
# distinguishing text in the enqueued item's own `reason`), just not a
# perfect match to this comment's own "scientific judgment call"
# framing for the FAIL sub-case specifically.
_OUTCOME_PRIORITY: dict[IngestionOutcome, CurationPriority] = {
    IngestionOutcome.UNRESOLVED: CurationPriority.HIGH,
    IngestionOutcome.IMPLAUSIBLE: CurationPriority.HIGH,
    IngestionOutcome.BAD_UNIT: CurationPriority.MEDIUM,
    IngestionOutcome.BAD_KIND: CurationPriority.MEDIUM,
    IngestionOutcome.ERROR: CurationPriority.MEDIUM,
}


@dataclass
class BridgeResult:
    item_ids_enqueued: list[str] = field(default_factory=list)

    @property
    def n_enqueued(self) -> int:
        return len(self.item_ids_enqueued)

    def summary(self) -> dict:
        return {"n_enqueued": self.n_enqueued}


def _paper_key(report: IngestionReport) -> str:
    """Stable identifier for the paper a report came from, preferring
    pubmed_id (the most common, most stable case) over doi over the raw
    title -- mirrors Provenance's own documented fallback order used
    elsewhere in this codebase (biosim_parameter_db.curation_bridge's
    own _source() helper)."""
    return report.pubmed_id or report.doi or report.paper_title or "unknown paper"


def enqueue_needs_curation_outcomes(report: IngestionReport, queue: CurationQueue) -> BridgeResult:
    """
    Enqueue every ExtractionOutcome in `report` flagged needs_curation
    into `queue`, one CurationQueue item per flagged extraction.

    Idempotent: each item_id is built from the paper's own stable
    identifier plus the extraction's own (subject_name, parameter_kind)
    pair -- not a list index, which could shift across re-runs if the
    extractor returns extractions in a different order. Re-running the
    SAME paper through the SAME extractor enqueues the same item_ids,
    which CurationQueue.enqueue()'s own already-in-queue dedup handles
    (bumping priority if the new priority is higher, never duplicating).
    Two DIFFERENT flagged extractions for the same (subject, kind) pair
    from the SAME paper (a real but rare case -- e.g. two different
    values reported for the same quantity in one paper) will collide on
    item_id and only the first is retained; this is a disclosed
    limitation, not silently incorrect, since CurationQueue.enqueue()
    itself has no mechanism to distinguish "the same item reported
    again" from "a genuinely different item that happens to share a
    key" -- the same ambiguity every other bridge in this codebase
    accepts (e.g. enqueue_conflict()'s own canonicalized-pair item_id).

    Args:
        report: IngestionReport to read outcomes from
        queue:  CurationQueue to enqueue into

    Returns:
        BridgeResult listing every item_id enqueued this call.
    """
    result = BridgeResult()
    paper_key = _paper_key(report)

    for outcome in report.outcomes:
        if not outcome.needs_curation:
            continue

        raw = outcome.raw
        item_id = f"extraction:{paper_key}:{raw.subject_name}:{raw.parameter_kind}"
        priority = _OUTCOME_PRIORITY.get(outcome.outcome, CurationPriority.MEDIUM)
        value = outcome.normalized_value if outcome.normalized_value is not None else raw.value
        unit = outcome.normalized_unit or raw.unit

        item_id = queue.enqueue(
            item_id=item_id,
            item_type=CurationItemType.EXTRACTION,
            reason=f"{outcome.outcome.value}: {outcome.message}",
            subject_name=raw.subject_name,
            parameter_kind=raw.parameter_kind,
            value=value,
            unit=unit,
            source_description=report.paper_title or paper_key,
            priority=priority,
        )
        result.item_ids_enqueued.append(item_id)

    return result
