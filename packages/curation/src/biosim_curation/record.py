"""
Curation record types.

A CurationRecord is the permanent log of an expert's decision on a
parameter or extraction. Records are immutable once created.

Decisions:
    APPROVED  — parameter is valid; confidence upgraded to GOLD
    REJECTED  — parameter is invalid; marked inactive in the store
    REVISED   — curator corrected the value; a new parameter is created
    DEFERRED  — needs more information; stays in queue with a note

CurationItems are the items in the queue awaiting review.
They carry enough context for a curator to make a decision without
needing to look anything up.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
from typing import Optional
import uuid


class CurationDecision(str, Enum):
    APPROVED = "approved"
    REJECTED = "rejected"
    REVISED  = "revised"
    DEFERRED = "deferred"


class CurationPriority(int, Enum):
    CRITICAL = 1   # On the critical path of v1 model — block simulation
    HIGH     = 2   # Frequently used parameter — high impact on outputs
    MEDIUM   = 3   # Used in model but not on critical path
    LOW      = 4   # Edge case or rarely-used pathway


class CurationItemType(str, Enum):
    PARAMETER  = "parameter"   # A parameter in the ParameterStore
    EXTRACTION = "extraction"  # A RawExtraction from the literature pipeline


@dataclass(frozen=True)
class CurationRecord:
    """
    Permanent record of a curator's decision.
    Immutable — decisions cannot be changed, only superseded by new decisions.
    """
    record_id: str
    item_id: str
    item_type: CurationItemType
    curator_id: str
    decision: CurationDecision
    notes: str
    created_at: datetime

    revised_value: Optional[float] = None
    revised_unit: Optional[str] = None
    confidence_override: Optional[str] = None  # ConfidenceTier value if upgrading

    def __post_init__(self) -> None:
        if self.decision == CurationDecision.REVISED:
            if self.revised_value is None:
                raise ValueError(
                    "CurationRecord with decision=REVISED must provide revised_value. "
                    "A revision without a corrected value is meaningless."
                )
        if not self.notes.strip():
            raise ValueError(
                "CurationRecord requires non-empty notes. "
                "Every decision must be documented with a rationale."
            )
        if not self.curator_id.strip():
            raise ValueError("curator_id cannot be empty.")

    @classmethod
    def approve(
        cls,
        item_id: str,
        item_type: CurationItemType,
        curator_id: str,
        notes: str,
    ) -> "CurationRecord":
        return cls(
            record_id=str(uuid.uuid4()),
            item_id=item_id,
            item_type=item_type,
            curator_id=curator_id,
            decision=CurationDecision.APPROVED,
            notes=notes,
            created_at=datetime.now(timezone.utc),
            confidence_override="gold",
        )

    @classmethod
    def reject(
        cls,
        item_id: str,
        item_type: CurationItemType,
        curator_id: str,
        notes: str,
    ) -> "CurationRecord":
        return cls(
            record_id=str(uuid.uuid4()),
            item_id=item_id,
            item_type=item_type,
            curator_id=curator_id,
            decision=CurationDecision.REJECTED,
            notes=notes,
            created_at=datetime.now(timezone.utc),
        )

    @classmethod
    def revise(
        cls,
        item_id: str,
        item_type: CurationItemType,
        curator_id: str,
        notes: str,
        revised_value: float,
        revised_unit: str,
    ) -> "CurationRecord":
        return cls(
            record_id=str(uuid.uuid4()),
            item_id=item_id,
            item_type=item_type,
            curator_id=curator_id,
            decision=CurationDecision.REVISED,
            notes=notes,
            revised_value=revised_value,
            revised_unit=revised_unit,
            created_at=datetime.now(timezone.utc),
        )

    @classmethod
    def defer(
        cls,
        item_id: str,
        item_type: CurationItemType,
        curator_id: str,
        notes: str,
    ) -> "CurationRecord":
        return cls(
            record_id=str(uuid.uuid4()),
            item_id=item_id,
            item_type=item_type,
            curator_id=curator_id,
            decision=CurationDecision.DEFERRED,
            notes=notes,
            created_at=datetime.now(timezone.utc),
        )


@dataclass
class CurationItem:
    """
    An item in the curation queue.
    Carries enough context for a curator to make a decision.
    """
    item_id: str
    item_type: CurationItemType
    priority: CurationPriority
    reason: str                    # Why this was flagged for curation
    subject_name: str              # Human-readable entity name
    parameter_kind: str            # e.g. "Km", "kcat"
    value: float
    unit: str
    source_description: str        # Paper/source this came from
    flagged_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
    assigned_to: Optional[str] = None
    decision_records: list[CurationRecord] = field(default_factory=list)

    @property
    def is_resolved(self) -> bool:
        return any(
            r.decision in (CurationDecision.APPROVED, CurationDecision.REJECTED, CurationDecision.REVISED)
            for r in self.decision_records
        )

    @property
    def latest_decision(self) -> Optional[CurationDecision]:
        if not self.decision_records:
            return None
        return sorted(self.decision_records, key=lambda r: r.created_at)[-1].decision
