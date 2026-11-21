"""Tests for biosim_curation."""
import pytest
from biosim_curation import (
    CurationDecision, CurationItem, CurationItemType,
    CurationPriority, CurationQueue, CurationQueueError, CurationRecord,
)

PARAM_ID = "param-001"
CURATOR = "dr.jones"

def make_item_kwargs(**kw):
    defaults = dict(
        item_id=PARAM_ID, item_type=CurationItemType.PARAMETER,
        reason="Implausible value", subject_name="insulin receptor",
        parameter_kind="Km", value=0.45, unit="uM",
        source_description="PMID:12345678",
        priority=CurationPriority.HIGH,
    )
    defaults.update(kw)
    return defaults


class TestCurationRecord:
    def test_approve_factory(self):
        r = CurationRecord.approve(PARAM_ID, CurationItemType.PARAMETER, CURATOR, "Verified against original paper")
        assert r.decision == CurationDecision.APPROVED
        assert r.confidence_override == "gold"

    def test_reject_factory(self):
        r = CurationRecord.reject(PARAM_ID, CurationItemType.PARAMETER, CURATOR, "Wrong cell type")
        assert r.decision == CurationDecision.REJECTED

    def test_revise_factory(self):
        r = CurationRecord.revise(PARAM_ID, CurationItemType.PARAMETER, CURATOR, "Corrected unit", 0.045, "uM")
        assert r.decision == CurationDecision.REVISED
        assert r.revised_value == 0.045

    def test_defer_factory(self):
        r = CurationRecord.defer(PARAM_ID, CurationItemType.PARAMETER, CURATOR, "Need original paper")
        assert r.decision == CurationDecision.DEFERRED

    def test_revised_without_value_raises(self):
        with pytest.raises(ValueError, match="revised_value"):
            CurationRecord(
                record_id="r1", item_id=PARAM_ID, item_type=CurationItemType.PARAMETER,
                curator_id=CURATOR, decision=CurationDecision.REVISED,
                notes="Some note", created_at=__import__('datetime').datetime.now(),
            )

    def test_empty_notes_raises(self):
        with pytest.raises(ValueError, match="notes"):
            CurationRecord(
                record_id="r1", item_id=PARAM_ID, item_type=CurationItemType.PARAMETER,
                curator_id=CURATOR, decision=CurationDecision.APPROVED,
                notes="   ", created_at=__import__('datetime').datetime.now(),
            )

    def test_empty_curator_raises(self):
        with pytest.raises(ValueError, match="curator_id"):
            CurationRecord(
                record_id="r1", item_id=PARAM_ID, item_type=CurationItemType.PARAMETER,
                curator_id="  ", decision=CurationDecision.APPROVED,
                notes="Valid note", created_at=__import__('datetime').datetime.now(),
            )

    def test_record_is_immutable(self):
        r = CurationRecord.approve(PARAM_ID, CurationItemType.PARAMETER, CURATOR, "Looks good")
        with pytest.raises(Exception):
            r.decision = CurationDecision.REJECTED


class TestCurationQueue:
    def test_enqueue_and_get(self):
        q = CurationQueue()
        q.enqueue(**make_item_kwargs())
        item = q.get(PARAM_ID)
        assert item.item_id == PARAM_ID
        assert item.priority == CurationPriority.HIGH

    def test_get_unknown_raises(self):
        q = CurationQueue()
        with pytest.raises(CurationQueueError):
            q.get("nonexistent")

    def test_pending_returns_unresolved(self):
        q = CurationQueue()
        q.enqueue(**make_item_kwargs())
        assert len(q.pending()) == 1

    def test_pending_empty_after_approval(self):
        q = CurationQueue()
        q.enqueue(**make_item_kwargs())
        r = CurationRecord.approve(PARAM_ID, CurationItemType.PARAMETER, CURATOR, "Verified")
        q.record_decision(PARAM_ID, r)
        assert len(q.pending()) == 0

    def test_resolved_after_approval(self):
        q = CurationQueue()
        q.enqueue(**make_item_kwargs())
        r = CurationRecord.approve(PARAM_ID, CurationItemType.PARAMETER, CURATOR, "Verified")
        q.record_decision(PARAM_ID, r)
        assert len(q.resolved()) == 1

    def test_priority_ordering(self):
        q = CurationQueue()
        q.enqueue(**make_item_kwargs(item_id="low", priority=CurationPriority.LOW))
        q.enqueue(**make_item_kwargs(item_id="critical", priority=CurationPriority.CRITICAL))
        q.enqueue(**make_item_kwargs(item_id="medium", priority=CurationPriority.MEDIUM))
        items = q.pending()
        assert items[0].item_id == "critical"
        assert items[-1].item_id == "low"

    def test_duplicate_enqueue_upgrades_priority(self):
        q = CurationQueue()
        q.enqueue(**make_item_kwargs(priority=CurationPriority.LOW))
        q.enqueue(**make_item_kwargs(priority=CurationPriority.CRITICAL))
        assert q.get(PARAM_ID).priority == CurationPriority.CRITICAL

    def test_duplicate_enqueue_does_not_downgrade_priority(self):
        q = CurationQueue()
        q.enqueue(**make_item_kwargs(priority=CurationPriority.CRITICAL))
        q.enqueue(**make_item_kwargs(priority=CurationPriority.LOW))
        assert q.get(PARAM_ID).priority == CurationPriority.CRITICAL

    def test_double_resolve_raises(self):
        q = CurationQueue()
        q.enqueue(**make_item_kwargs())
        r = CurationRecord.approve(PARAM_ID, CurationItemType.PARAMETER, CURATOR, "OK")
        q.record_decision(PARAM_ID, r)
        with pytest.raises(CurationQueueError, match="already resolved"):
            q.record_decision(PARAM_ID, r)

    def test_assign_and_filter(self):
        q = CurationQueue()
        q.enqueue(**make_item_kwargs(item_id="p1"))
        q.enqueue(**make_item_kwargs(item_id="p2"))
        q.assign("p1", "dr.jones")
        results = q.pending(curator_id="dr.jones")
        assert len(results) == 1
        assert results[0].item_id == "p1"

    def test_assign_resolved_raises(self):
        q = CurationQueue()
        q.enqueue(**make_item_kwargs())
        r = CurationRecord.approve(PARAM_ID, CurationItemType.PARAMETER, CURATOR, "OK")
        q.record_decision(PARAM_ID, r)
        with pytest.raises(CurationQueueError):
            q.assign(PARAM_ID, "dr.smith")

    def test_unassign(self):
        q = CurationQueue()
        q.enqueue(**make_item_kwargs())
        q.assign(PARAM_ID, "dr.jones")
        q.unassign(PARAM_ID)
        assert q.get(PARAM_ID).assigned_to is None

    def test_deferred_stays_in_pending(self):
        q = CurationQueue()
        q.enqueue(**make_item_kwargs())
        r = CurationRecord.defer(PARAM_ID, CurationItemType.PARAMETER, CURATOR, "Need more info")
        q.record_decision(PARAM_ID, r)
        assert len(q.pending()) == 1
        assert len(q.deferred()) == 1

    def test_enqueue_conflict(self):
        q = CurationQueue()
        q.enqueue_conflict(
            "p1", "p2", "Insulin receptor", "Km",
            10.0, 100.0, "uM", 10.0, "PMID:111", "PMID:222"
        )
        items = q.pending()
        assert len(items) == 1
        assert items[0].priority == CurationPriority.CRITICAL

    def test_enqueue_conflict_same_pair_reversed_order_is_the_same_item(self):
        # Regression: conflict_id used to be built as
        # f"conflict:{parameter_id_a}:{parameter_id_b}" with no
        # canonicalization, so the identical conflict reported with its
        # two parameter IDs swapped (e.g. a re-run iterating pairs in a
        # different order) produced a DIFFERENT item_id and bypassed
        # enqueue()'s already-in-queue dedup entirely, creating a
        # duplicate curation item for the same underlying conflict.
        q = CurationQueue()
        first_id = q.enqueue_conflict(
            "p1", "p2", "Insulin receptor", "Km",
            10.0, 100.0, "uM", 10.0, "PMID:111", "PMID:222"
        )
        second_id = q.enqueue_conflict(
            "p2", "p1", "Insulin receptor", "Km",
            100.0, 10.0, "uM", 10.0, "PMID:222", "PMID:111"
        )
        assert first_id == second_id
        assert len(q.pending()) == 1

    def test_filter_by_type(self):
        q = CurationQueue()
        q.enqueue(**make_item_kwargs(item_id="p1", item_type=CurationItemType.PARAMETER))
        q.enqueue(**make_item_kwargs(item_id="e1", item_type=CurationItemType.EXTRACTION))
        params = q.pending(item_type=CurationItemType.PARAMETER)
        assert len(params) == 1

    def test_stats(self):
        q = CurationQueue()
        q.enqueue(**make_item_kwargs(item_id="p1", priority=CurationPriority.CRITICAL))
        q.enqueue(**make_item_kwargs(item_id="p2", priority=CurationPriority.HIGH))
        r = CurationRecord.approve("p1", CurationItemType.PARAMETER, CURATOR, "OK")
        q.record_decision("p1", r)
        stats = q.stats()
        assert stats["total_items"] == 2
        assert stats["pending"] == 1
        assert stats["resolved"] == 1
        assert stats["by_priority"]["HIGH"] == 1

    def test_limit_pending(self):
        q = CurationQueue()
        for i in range(10):
            q.enqueue(**make_item_kwargs(item_id=f"p{i}"))
        assert len(q.pending(limit=3)) == 3

    def test_is_resolved_false_initially(self):
        item = CurationItem(
            item_id="x", item_type=CurationItemType.PARAMETER,
            priority=CurationPriority.MEDIUM, reason="test",
            subject_name="test", parameter_kind="Km",
            value=1.0, unit="uM", source_description="test",
        )
        assert not item.is_resolved

    def test_latest_decision_none_when_empty(self):
        item = CurationItem(
            item_id="x", item_type=CurationItemType.PARAMETER,
            priority=CurationPriority.MEDIUM, reason="test",
            subject_name="test", parameter_kind="Km",
            value=1.0, unit="uM", source_description="test",
        )
        assert item.latest_decision is None
