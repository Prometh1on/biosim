"""
Curation queue — manages the backlog of parameters and extractions
awaiting expert review.

Design principles:
    - Items are never deleted — only resolved or deferred
    - Priority ordering: CRITICAL first, then by flagged_at (oldest first)
    - Assignments prevent duplicate review work
    - Full audit trail: every state change is logged
    - Stats feed the curation dashboard (Sprint 10)

Thread safety: all public methods acquire the internal lock.
"""

from __future__ import annotations

import threading
import uuid
from datetime import datetime, timezone
from typing import Optional

from .record import (
    CurationDecision,
    CurationItem,
    CurationItemType,
    CurationPriority,
    CurationRecord,
)


class CurationQueueError(Exception):
    pass


class CurationQueue:
    """
    Priority-ordered queue of items awaiting expert curation.

    Items enter the queue from two sources:
        1. Literature pipeline — flagged extractions (bad unit, implausible value,
           unresolved entity, low extractor confidence)
        2. Parameter store — conflicts between parameters from different papers

    Items leave the queue when a curator records an APPROVED, REJECTED,
    or REVISED decision. DEFERRED items stay in the queue.
    """

    def __init__(self) -> None:
        self._items: dict[str, CurationItem] = {}
        self._lock = threading.RLock()

    # ------------------------------------------------------------------
    # Enqueue
    # ------------------------------------------------------------------

    def enqueue(
        self,
        item_id: str,
        item_type: CurationItemType,
        reason: str,
        subject_name: str,
        parameter_kind: str,
        value: float,
        unit: str,
        source_description: str,
        priority: CurationPriority = CurationPriority.MEDIUM,
    ) -> str:
        """
        Add an item to the curation queue.

        Returns the item_id (same as input, for convenience).
        If the item is already in the queue, updates the priority if the
        new priority is higher (lower numeric value).
        """
        with self._lock:
            if item_id in self._items:
                existing = self._items[item_id]
                if priority.value < existing.priority.value:
                    existing.priority = priority
                    existing.reason = f"{existing.reason} | {reason}"
                return item_id

            item = CurationItem(
                item_id=item_id,
                item_type=item_type,
                priority=priority,
                reason=reason,
                subject_name=subject_name,
                parameter_kind=parameter_kind,
                value=value,
                unit=unit,
                source_description=source_description,
            )
            self._items[item_id] = item
            return item_id

    def enqueue_conflict(
        self,
        parameter_id_a: str,
        parameter_id_b: str,
        subject_name: str,
        parameter_kind: str,
        value_a: float,
        value_b: float,
        unit: str,
        fold_difference: float,
        source_a: str,
        source_b: str,
    ) -> str:
        """Convenience method for enqueueing a parameter conflict."""
        # Sort so the same conflict reported with its two parameter IDs
        # swapped (e.g. by a re-run that iterates pairs in a different
        # order) canonicalizes to the same item_id -- enqueue()'s
        # already-in-queue dedup keys on item_id, so an order-dependent
        # id here would let the identical conflict be enqueued twice
        # instead of being recognized as the same item.
        id_a, id_b = sorted([parameter_id_a, parameter_id_b])
        conflict_id = f"conflict:{id_a}:{id_b}"
        priority = (
            CurationPriority.CRITICAL if fold_difference >= 10.0
            else CurationPriority.HIGH if fold_difference >= 5.0
            else CurationPriority.MEDIUM
        )
        return self.enqueue(
            item_id=conflict_id,
            item_type=CurationItemType.PARAMETER,
            reason=(
                f"Conflict: {fold_difference:.1f}x difference between "
                f"value {value_a} ({source_a}) and {value_b} ({source_b})"
            ),
            subject_name=subject_name,
            parameter_kind=parameter_kind,
            value=value_a,
            unit=unit,
            source_description=f"{source_a} vs {source_b}",
            priority=priority,
        )

    # ------------------------------------------------------------------
    # Assignment
    # ------------------------------------------------------------------

    def assign(self, item_id: str, curator_id: str) -> None:
        """Assign an item to a specific curator."""
        with self._lock:
            item = self._get_or_raise(item_id)
            if item.is_resolved:
                raise CurationQueueError(
                    f"Cannot assign {item_id} — it is already resolved."
                )
            item.assigned_to = curator_id

    def unassign(self, item_id: str) -> None:
        """Remove the assignment from an item, returning it to the pool."""
        with self._lock:
            item = self._get_or_raise(item_id)
            item.assigned_to = None

    # ------------------------------------------------------------------
    # Decision recording
    # ------------------------------------------------------------------

    def record_decision(self, item_id: str, record: CurationRecord) -> None:
        """
        Record a curation decision for an item.

        Args:
            item_id: The item being decided on
            record: The CurationRecord capturing the decision

        Raises:
            CurationQueueError: If item not found or already resolved
        """
        with self._lock:
            item = self._get_or_raise(item_id)
            if item.is_resolved:
                raise CurationQueueError(
                    f"Item {item_id} is already resolved with decision "
                    f"'{item.latest_decision}'. "
                    f"To change a decision, contact the platform administrator."
                )
            item.decision_records.append(record)

    # ------------------------------------------------------------------
    # Querying
    # ------------------------------------------------------------------

    def get(self, item_id: str) -> CurationItem:
        """Get a specific item by ID."""
        with self._lock:
            return self._get_or_raise(item_id)

    def pending(
        self,
        curator_id: Optional[str] = None,
        item_type: Optional[CurationItemType] = None,
        min_priority: CurationPriority = CurationPriority.LOW,
        limit: Optional[int] = None,
    ) -> list[CurationItem]:
        """
        Return pending (unresolved) items sorted by priority then age.

        Args:
            curator_id: Filter to items assigned to this curator (None = all)
            item_type: Filter by item type (None = all)
            min_priority: Minimum priority level (CRITICAL=1 is highest)
            limit: Maximum number of items to return
        """
        with self._lock:
            items = [
                item for item in self._items.values()
                if not item.is_resolved
                and item.priority.value <= min_priority.value
                and (curator_id is None or item.assigned_to == curator_id)
                and (item_type is None or item.item_type == item_type)
            ]
            items.sort(key=lambda i: (i.priority.value, i.flagged_at))
            if limit is not None:
                items = items[:limit]
            return items

    def resolved(self) -> list[CurationItem]:
        """Return all resolved items."""
        with self._lock:
            return [i for i in self._items.values() if i.is_resolved]

    def deferred(self) -> list[CurationItem]:
        """Return items with a DEFERRED decision but not yet finally resolved."""
        with self._lock:
            return [
                i for i in self._items.values()
                if not i.is_resolved and i.latest_decision == CurationDecision.DEFERRED
            ]

    def stats(self) -> dict:
        """Return queue statistics for the curation dashboard."""
        with self._lock:
            all_items = list(self._items.values())
            pending = [i for i in all_items if not i.is_resolved]
            by_priority: dict[str, int] = {}
            for item in pending:
                k = item.priority.name
                by_priority[k] = by_priority.get(k, 0) + 1
            by_type: dict[str, int] = {}
            for item in pending:
                k = item.item_type.value
                by_type[k] = by_type.get(k, 0) + 1
            return {
                "total_items": len(all_items),
                "pending": len(pending),
                "resolved": len([i for i in all_items if i.is_resolved]),
                "deferred": len(self.deferred()),
                "assigned": len([i for i in pending if i.assigned_to]),
                "unassigned": len([i for i in pending if not i.assigned_to]),
                "by_priority": by_priority,
                "by_type": by_type,
            }

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    def _get_or_raise(self, item_id: str) -> CurationItem:
        if item_id not in self._items:
            raise CurationQueueError(f"Item '{item_id}' not found in curation queue.")
        return self._items[item_id]
