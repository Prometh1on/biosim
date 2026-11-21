# SAVE TO: packages\parameter-db\tests\test_snapshot_validator.py
"""
Test suite for biosim_parameter_db.snapshot_validator (MASTER_PLAN 3.6).

Coverage:
    - A real, valid snapshot (built via ParameterStore.save(), not
      hand-crafted) reports zero errors and correct counts
    - Missing file, malformed JSON both report an error, never raise
    - Duplicate parameter_id, missing parameter_id/seq
    - Out-of-order seq, duplicate seq
    - Dangling supersedes_parameter_id / conflict references
    - Copy-then-read never mutates the original file
    - snapshot_sha256 matches the real file content
"""

import hashlib
import json

import pytest

from biosim_parameter_db import (
    ConfidenceTier,
    Organism,
    ParameterBuilder,
    ParameterKind,
    ParameterStore,
    validate_snapshot,
)

INSULIN_RECEPTOR_ID = "BIOSIM:PROTEIN:000001"


def make_km(value: float = 45.0, tier: ConfidenceTier = ConfidenceTier.SILVER):
    return (
        ParameterBuilder()
        .for_subject(INSULIN_RECEPTOR_ID)
        .kind(ParameterKind.KM)
        .value(value, unit="uM")
        .from_pubmed("12345678")
        .organism(Organism.HUMAN)
        .confidence(tier)
        .build()
    )


class TestValidSnapshot:
    def test_real_snapshot_has_no_errors(self, tmp_path):
        store = ParameterStore()
        store.add(make_km(value=45.0, tier=ConfidenceTier.GOLD))
        store.add(make_km(value=10.0))  # same subject/kind -> a real conflict
        old_id = store.add(make_km(value=1.0))
        store.supersede(old_id, make_km(value=2.0))
        path = tmp_path / "snapshot.json"
        store.save(path)

        report = validate_snapshot(path)
        assert report.is_valid, report.errors
        assert report.errors == []
        assert report.n_parameters == 4
        # Several of these values differ by more than the 2x fold
        # threshold, so genuine ConflictRecords are expected -- the exact
        # count isn't the point of this test, only that a real,
        # store-produced snapshot validates clean.
        assert report.n_conflicts >= 1

    def test_snapshot_sha256_matches_real_file_content(self, tmp_path):
        store = ParameterStore()
        store.add(make_km())
        path = tmp_path / "snapshot.json"
        store.save(path)

        report = validate_snapshot(path)
        expected = hashlib.sha256(path.read_bytes()).hexdigest()
        assert report.snapshot_sha256 == expected

    def test_copy_then_read_does_not_mutate_the_original_file(self, tmp_path):
        store = ParameterStore()
        store.add(make_km())
        path = tmp_path / "snapshot.json"
        store.save(path)
        before = path.read_bytes()

        validate_snapshot(path)

        after = path.read_bytes()
        assert before == after


class TestMissingOrMalformedFile:
    def test_missing_file_reports_error_not_raise(self, tmp_path):
        report = validate_snapshot(tmp_path / "does_not_exist.json")
        assert not report.is_valid
        assert any("does not exist" in e for e in report.errors)

    def test_malformed_json_reports_error_not_raise(self, tmp_path):
        path = tmp_path / "broken.json"
        path.write_text("{not valid json,,,")
        report = validate_snapshot(path)
        assert not report.is_valid
        assert any("not valid JSON" in e for e in report.errors)

    def test_non_array_parameters_field_reports_error(self, tmp_path):
        path = tmp_path / "snapshot.json"
        path.write_text(json.dumps({"parameters": "oops, a string", "conflicts": []}))
        report = validate_snapshot(path)
        assert not report.is_valid
        assert any("not a JSON array" in e for e in report.errors)

    def test_non_array_conflicts_field_reports_error(self, tmp_path):
        # Code-review gap: the symmetric "parameters not an array" case
        # above had a test; this one didn't.
        path = tmp_path / "snapshot.json"
        path.write_text(json.dumps({"parameters": [], "conflicts": "oops, a string"}))
        report = validate_snapshot(path)
        assert not report.is_valid
        assert any("not a JSON array" in e for e in report.errors)


class TestCorruptedParameterFields:
    def _write(self, tmp_path, parameters, conflicts=None):
        path = tmp_path / "snapshot.json"
        path.write_text(json.dumps({"parameters": parameters, "conflicts": conflicts or []}))
        return path

    def test_duplicate_parameter_id_detected(self, tmp_path):
        path = self._write(tmp_path, [
            {"parameter_id": "p1", "seq": 0},
            {"parameter_id": "p1", "seq": 1},
        ])
        report = validate_snapshot(path)
        assert not report.is_valid
        assert any("duplicate parameter_id" in e for e in report.errors)

    def test_missing_parameter_id_detected(self, tmp_path):
        path = self._write(tmp_path, [{"seq": 0}])
        report = validate_snapshot(path)
        assert not report.is_valid
        assert any("missing 'parameter_id'" in e for e in report.errors)

    def test_missing_seq_detected(self, tmp_path):
        path = self._write(tmp_path, [{"parameter_id": "p1"}])
        report = validate_snapshot(path)
        assert not report.is_valid
        assert any("missing 'seq'" in e for e in report.errors)

    def test_out_of_order_seq_detected(self, tmp_path):
        # Array order (p1 then p2) contradicts seq order (p2's seq=0 is
        # lower than p1's seq=5) -- exactly the "reordered snapshot"
        # scenario a single-threaded backfill over this array would get
        # wrong.
        path = self._write(tmp_path, [
            {"parameter_id": "p1", "seq": 5},
            {"parameter_id": "p2", "seq": 0},
        ])
        report = validate_snapshot(path)
        assert not report.is_valid
        assert any("array order does not match seq order" in e for e in report.errors)

    def test_duplicate_seq_detected(self, tmp_path):
        path = self._write(tmp_path, [
            {"parameter_id": "p1", "seq": 0},
            {"parameter_id": "p2", "seq": 0},
        ])
        report = validate_snapshot(path)
        assert not report.is_valid
        assert any("duplicate seq=0" in e for e in report.errors)

    def test_ascending_seq_with_gaps_is_not_flagged(self, tmp_path):
        # Gaps are fine -- only ORDER matters, not contiguity (a real
        # store's seq counter never has gaps today, but nothing in the
        # snapshot format guarantees that forever).
        path = self._write(tmp_path, [
            {"parameter_id": "p1", "seq": 0},
            {"parameter_id": "p2", "seq": 7},
            {"parameter_id": "p3", "seq": 100},
        ])
        report = validate_snapshot(path)
        assert report.is_valid, report.errors

    def test_dangling_supersedes_parameter_id_detected(self, tmp_path):
        path = self._write(tmp_path, [
            {"parameter_id": "p1", "seq": 0, "supersedes_parameter_id": "ghost"},
        ])
        report = validate_snapshot(path)
        assert not report.is_valid
        assert any("supersedes_parameter_id" in e and "ghost" in e for e in report.errors)

    def test_valid_supersedes_parameter_id_not_flagged(self, tmp_path):
        path = self._write(tmp_path, [
            {"parameter_id": "p1", "seq": 0},
            {"parameter_id": "p2", "seq": 1, "supersedes_parameter_id": "p1"},
        ])
        report = validate_snapshot(path)
        assert report.is_valid, report.errors

    def test_dangling_conflict_reference_detected(self, tmp_path):
        path = self._write(
            tmp_path,
            [{"parameter_id": "p1", "seq": 0}],
            conflicts=[{
                "conflict_id": "c1",
                "existing_parameter_id": "p1",
                "incoming_parameter_id": "ghost",
            }],
        )
        report = validate_snapshot(path)
        assert not report.is_valid
        assert any(
            "incoming_parameter_id='ghost'" in e and "does not resolve" in e
            for e in report.errors
        )

    def test_valid_conflict_reference_not_flagged(self, tmp_path):
        path = self._write(
            tmp_path,
            [
                {"parameter_id": "p1", "seq": 0},
                {"parameter_id": "p2", "seq": 1},
            ],
            conflicts=[{
                "conflict_id": "c1",
                "existing_parameter_id": "p1",
                "incoming_parameter_id": "p2",
            }],
        )
        report = validate_snapshot(path)
        assert report.is_valid, report.errors

    def test_conflict_missing_field_detected(self, tmp_path):
        # Code-review gap: the "missing field" branch (as opposed to
        # "field present but dangling") had no direct test.
        path = self._write(
            tmp_path,
            [{"parameter_id": "p1", "seq": 0}],
            conflicts=[{"conflict_id": "c1", "existing_parameter_id": "p1"}],
        )
        report = validate_snapshot(path)
        assert not report.is_valid
        assert any("missing 'incoming_parameter_id'" in e for e in report.errors)


class TestMalformedElementsDoNotCrash:
    # Code-review HIGH finding: the validator's own contract is "never
    # raises" -- these confirm that promise actually holds for realistic
    # malformed ELEMENTS within an otherwise-valid array (a non-dict
    # entry, or a non-hashable id/seq value), not just a malformed file
    # as a whole. Before the fix, each of these crashed with an unhandled
    # AttributeError or TypeError instead of returning a report.
    def _write(self, tmp_path, parameters, conflicts=None):
        path = tmp_path / "snapshot.json"
        path.write_text(json.dumps({"parameters": parameters, "conflicts": conflicts or []}))
        return path

    def test_non_dict_parameter_element_does_not_crash(self, tmp_path):
        path = self._write(tmp_path, ["not an object", 123, None])
        report = validate_snapshot(path)
        assert not report.is_valid
        assert any("not a JSON object" in e for e in report.errors)

    def test_non_scalar_parameter_id_does_not_crash(self, tmp_path):
        # A list value for parameter_id is unhashable -- `pid in
        # known_parameter_ids` used to raise TypeError here.
        path = self._write(tmp_path, [{"parameter_id": ["a", "b"], "seq": 0}])
        report = validate_snapshot(path)
        assert not report.is_valid
        assert any("'parameter_id' must be a plain scalar" in e for e in report.errors)

    def test_non_scalar_seq_does_not_crash(self, tmp_path):
        path = self._write(tmp_path, [{"parameter_id": "p1", "seq": {"nested": "dict"}}])
        report = validate_snapshot(path)
        assert not report.is_valid
        assert any("'seq' must be a plain scalar" in e for e in report.errors)

    def test_non_scalar_supersedes_parameter_id_does_not_crash(self, tmp_path):
        path = self._write(tmp_path, [
            {"parameter_id": "p1", "seq": 0, "supersedes_parameter_id": ["not", "hashable"]},
        ])
        report = validate_snapshot(path)
        assert not report.is_valid
        assert any("'supersedes_parameter_id' must be a plain scalar" in e for e in report.errors)

    def test_non_dict_conflict_element_does_not_crash(self, tmp_path):
        path = self._write(tmp_path, [{"parameter_id": "p1", "seq": 0}], conflicts=["not an object"])
        report = validate_snapshot(path)
        assert not report.is_valid
        assert any("not a JSON object" in e for e in report.errors)

    def test_non_scalar_conflict_reference_does_not_crash(self, tmp_path):
        path = self._write(
            tmp_path,
            [{"parameter_id": "p1", "seq": 0}],
            conflicts=[{
                "conflict_id": "c1",
                "existing_parameter_id": "p1",
                "incoming_parameter_id": {"nested": "dict"},
            }],
        )
        report = validate_snapshot(path)
        assert not report.is_valid
        assert any(
            "incoming_parameter_id' must be a plain scalar" in e for e in report.errors
        )

    def test_non_dict_top_level_does_not_crash(self, tmp_path):
        path = tmp_path / "snapshot.json"
        path.write_text(json.dumps(["just", "a", "list"]))
        report = validate_snapshot(path)
        assert not report.is_valid
        assert any("not a JSON object" in e for e in report.errors)

    def test_copy_failure_between_existence_check_and_read_does_not_crash(self, tmp_path, monkeypatch):
        # Code-review MEDIUM finding: path.is_file() and the actual copy
        # are not atomic -- a concurrent delete/permission change in that
        # window used to propagate an uncaught OSError.
        import biosim_parameter_db.snapshot_validator as module

        path = tmp_path / "snapshot.json"
        path.write_text(json.dumps({"parameters": [], "conflicts": []}))

        def _raise(*args, **kwargs):
            raise OSError("simulated race: file vanished mid-copy")

        monkeypatch.setattr(module.shutil, "copyfile", _raise)
        report = module.validate_snapshot(path)
        assert not report.is_valid
        assert any("could not read snapshot file" in e for e in report.errors)
