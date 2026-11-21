# SAVE TO: packages\parameter-db\src\biosim_parameter_db\curation_bridge.py
"""
Conflict-to-curation bridge (MASTER_PLAN 5.2a).

Closes a real gap found while scoping `CurationQueue`'s own workflow
(MASTER_PLAN 5.2): `ParameterStore` already actively detects and
records genuine scientific conflicts between competing literature
values (`ConflictRecord`, via `unresolved_conflicts()`) completely
independently of any LLM/API key -- but nothing anywhere surfaced that
real data to the curation queue for a human to actually review. This
module is the missing wiring, mirroring
`biosim_uncertainty.curation_bridge`'s own shape (a standalone bridge
function callers invoke explicitly, not automatic/hooked into `add()`
-- conflict detection and curation-queue surfacing are different
concerns, and forcing every `add()` call to also touch a `CurationQueue`
would couple two packages that don't otherwise depend on each other).

Flow:
    1. ParameterStore.add() detects and records conflicts as it goes
       (already existing, unaffected by this module)
    2. This bridge reads store.unresolved_conflicts(), looks up both
       parameters' units/provenance and the subject's human-readable
       name, and enqueues each as a CurationQueue conflict item
    3. Curators review and resolve items from the queue (existing
       CurationQueue.record_decision(); still separate from
       ParameterStore.resolve_conflict() -- resolving the CURATION item
       does not itself resolve the STORE's own conflict record, and
       vice versa; a curator acting on either should update both,
       matching this project's "no half-finished workflow" preference
       over silently auto-syncing two independently-owned records)
"""

from __future__ import annotations

from dataclasses import dataclass, field


def _source(param) -> str:
    """Provenance display string, in Provenance's own documented
    preference order (pubmed_id preferred, doi/arxiv_id fallback,
    source_description last)."""
    p = param.provenance
    return p.pubmed_id or p.doi or p.arxiv_id or p.source_description or "unknown source"


@dataclass
class BridgeResult:
    conflict_ids_enqueued: list[str] = field(default_factory=list)
    # Mandatory code-reviewer finding (MEDIUM): ParameterStore's own
    # conflict detection (_detect_conflict()) never checks that the two
    # compared parameters share a unit -- reachable via direct
    # ParameterBuilder use (which performs no cross-parameter unit
    # validation), even though it isn't hit by any real production data
    # path today. Rather than silently enqueue a fold_difference/unit
    # pairing that could be scientifically meaningless, conflicts with
    # mismatched units are skipped here and counted separately -- a
    # curator should never see a same-unit-looking comparison that
    # isn't actually one.
    conflict_ids_skipped_unit_mismatch: list[str] = field(default_factory=list)

    @property
    def n_conflicts_enqueued(self) -> int:
        return len(self.conflict_ids_enqueued)

    def summary(self) -> dict:
        return {
            "n_conflicts_enqueued": self.n_conflicts_enqueued,
            "n_skipped_unit_mismatch": len(self.conflict_ids_skipped_unit_mismatch),
        }


def enqueue_unresolved_conflicts(store, queue, registry) -> BridgeResult:
    """
    Enqueue every currently-unresolved ParameterStore conflict into the
    given CurationQueue.

    Idempotent: CurationQueue.enqueue_conflict() canonicalizes each
    conflict's item_id from its sorted parameter-ID pair, so calling
    this repeatedly (e.g. once per session) never creates a duplicate
    entry for the same conflict. Note this only means "safe to re-run,"
    not "re-running can change anything" -- ConflictRecord.fold_difference
    is set once at construction and never mutated, so enqueue_conflict()'s
    computed priority is identical on every call for a still-unresolved
    conflict; CurationQueue.enqueue()'s own general priority-upgrade path
    exists for other callers, not something this specific bridge can
    trigger.

    Args:
        store:    ParameterStore whose unresolved_conflicts() to read
        queue:    CurationQueue to enqueue into
        registry: IDRegistry, used to resolve each conflict's
                  subject_canonical_id into a human-readable name for
                  the curator (falls back to the raw canonical ID
                  string if the entity isn't registered or the ID
                  string is malformed -- a conflict is never silently
                  dropped just because this one cosmetic lookup fails)

    Returns:
        BridgeResult listing every conflict_id enqueued this call, and
        separately every conflict_id skipped for a unit mismatch
        between its two parameters (see BridgeResult's own docstring).
    """
    from biosim_id_resolver import CanonicalID, ResolutionError

    result = BridgeResult()

    for conflict in store.unresolved_conflicts():
        existing_param = store.get(conflict.existing_parameter_id)
        incoming_param = store.get(conflict.incoming_parameter_id)

        if existing_param.value.unit != incoming_param.value.unit:
            result.conflict_ids_skipped_unit_mismatch.append(conflict.conflict_id)
            continue

        try:
            canonical = CanonicalID.from_string(conflict.subject_canonical_id)
            subject_name = registry.get_record(canonical).name
        except (ResolutionError, ValueError):
            subject_name = conflict.subject_canonical_id

        queue.enqueue_conflict(
            parameter_id_a=conflict.existing_parameter_id,
            parameter_id_b=conflict.incoming_parameter_id,
            subject_name=subject_name,
            parameter_kind=conflict.kind.value,
            value_a=conflict.existing_value,
            value_b=conflict.incoming_value,
            unit=existing_param.value.unit,
            fold_difference=conflict.fold_difference,
            source_a=_source(existing_param),
            source_b=_source(incoming_param),
        )
        result.conflict_ids_enqueued.append(conflict.conflict_id)

    return result
