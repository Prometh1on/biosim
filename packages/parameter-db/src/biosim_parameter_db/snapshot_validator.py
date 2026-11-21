# SAVE TO: packages\parameter-db\src\biosim_parameter_db\snapshot_validator.py
"""
Pre-flight validation for a ParameterStore JSON snapshot, before it is
trusted as input to the Track 3 database-migration backfill
(infra/MIGRATION_PLAN.md Step 2, MASTER_PLAN row 3.6).

Deliberately narrow scope (MASTER_PLAN 3.6): this is the pure-Python
snapshot reader/validator only -- no SQL/Cypher emission, no live
database connection, testable now against real snapshots without
Docker. The actual INSERT/MERGE statement generation is a separate,
later chunk (row 3.7), informed by whatever this validator and the
eventual disposable Docker validation boot (row 3.5) surface.

Operational contract (infra/MIGRATION_PLAN.md's own resolved "who runs
the backfill, and when" open question):
    - NEVER opens a live ParameterStore object, in-process or otherwise
      -- takes a file PATH to a JSON snapshot as its only input, matching
      migration_audit.source_snapshot's own documented contract
      (infra/postgres/001_schema.sql: "path/hash of the ParameterStore
      JSON snapshot the backfill read").
    - Copy-then-read, not read-in-place. Disclosed honestly (code
      review): this is NOT a true atomicity guarantee -- shutil.copyfile()
      is an ordinary chunked read/write with no locking, exactly as
      exposed to a concurrent writer as a direct read would be, and
      ParameterStore.save() itself writes via a plain path.write_text()
      (truncate-then-write), not a write-temp-then-atomic-rename
      pattern, so there is no atomicity on the writer's side for this
      step to protect against in the first place. In PRACTICE a torn
      write from write_text() almost always yields empty or
      syntactically-broken JSON, which degrades to the (already
      correctly handled) "not valid JSON" error path rather than a
      silently-wrong validation -- but that safety comes from
      write_text()'s own truncate-then-write shape, not from this
      module's copy step. A real fix, if this ever matters in practice,
      is ParameterStore.save() writing to a temp file and os.replace()-ing
      it into place; not done here since it's a change to a different,
      widely-used module for a narrow, currently-hypothetical race.
    - Never writes back to the file-based store -- this module only reads.

Reads the RAW JSON directly, NOT via ParameterStore.load() -- that
loader applies graceful fallbacks for missing/malformed fields (e.g. a
missing "seq" silently falls back to the next available counter value,
a missing "detected_at" is silently omitted) specifically so a
long-lived in-memory store keeps working across incremental schema
additions. Those same fallbacks would mask exactly the kind of snapshot
corruption this validator exists to catch before a backfill trusts it.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


def _is_scalar_id(value: object) -> bool:
    """
    True if `value` is a plain JSON scalar (str/int/float/bool) -- safe
    to use as a dict/set key. A list or dict value (a real, foreseeable
    form of snapshot corruption) is NOT hashable and would crash a bare
    `x in some_set` check with TypeError -- found by code review.
    Callers check for None (a genuinely missing field) separately, with
    their own "missing field" message, before calling this.
    """
    return isinstance(value, (str, int, float, bool))


@dataclass
class SnapshotValidationReport:
    """
    Result of validating one ParameterStore JSON snapshot file.

    is_valid is True iff errors is empty. snapshot_sha256 is the SHA-256
    of the snapshot's exact byte content at the moment it was copied for
    reading -- suitable for migration_audit.source_snapshot's own
    documented "path/hash" format once the actual backfill (row 3.7)
    writes that record.
    """
    snapshot_path: str
    snapshot_sha256: str
    n_parameters: int
    n_conflicts: int
    errors: list[str] = field(default_factory=list)

    @property
    def is_valid(self) -> bool:
        return not self.errors


def validate_snapshot(path: Path) -> SnapshotValidationReport:
    """
    Validate a ParameterStore JSON snapshot before it is used as backfill
    input. Checks:

        1. The file exists and contains valid JSON with the expected
           top-level "parameters"/"conflicts" array shape.
        2. Every parameter has a "parameter_id" and a "seq", and no two
           parameters share either.
        3. The "parameters" array is already in strictly ascending "seq"
           order -- infra/MIGRATION_PLAN.md's own requirement that the
           backfill insert in EXACT snapshot array order, single-threaded,
           so a re-run after a partial failure (or a hand-edited/merged
           snapshot) can't silently diverge Postgres's `seq` assignment
           from the array order a naive re-read would otherwise trust.
        4. Every parameter's "supersedes_parameter_id" (when not null)
           resolves to a real parameter_id within this same snapshot.
        5. Every conflict's "existing_parameter_id"/"incoming_parameter_id"
           resolves to a real parameter_id within this same snapshot.

    Never raises -- on a missing file, an unreadable file, malformed
    JSON, or a malformed ELEMENT within an otherwise-valid array (a
    non-object entry, or an id/seq field that isn't a plain scalar),
    returns a report with errors populated instead of propagating an
    exception. This is a deliberate API guarantee (code review found and
    closed two real gaps in it: non-dict array elements and non-hashable
    id fields, both directly foreseeable forms of hand-edited/corrupted
    snapshot data, not just theoretical) -- callers can always inspect
    report.is_valid/report.errors uniformly, with no try/except needed
    around this call.
    """
    path = Path(path)
    if not path.is_file():
        return SnapshotValidationReport(
            snapshot_path=str(path),
            snapshot_sha256="",
            n_parameters=0,
            n_conflicts=0,
            errors=[f"snapshot file does not exist: {path}"],
        )

    try:
        with tempfile.TemporaryDirectory(prefix="biosim_snapshot_validate_") as tmp_dir:
            tmp_path = Path(tmp_dir) / path.name
            shutil.copyfile(path, tmp_path)
            raw_bytes = tmp_path.read_bytes()
    except OSError as e:
        # TOCTOU gap: path.is_file() above and the copy/read here are not
        # atomic -- the file could be deleted, renamed, or become
        # unreadable in between (code review finding). Converted to a
        # report entry rather than left to propagate, same as every
        # other failure mode this function documents as non-raising.
        return SnapshotValidationReport(
            snapshot_path=str(path),
            snapshot_sha256="",
            n_parameters=0,
            n_conflicts=0,
            errors=[f"could not read snapshot file: {type(e).__name__}: {e}"],
        )

    snapshot_sha256 = hashlib.sha256(raw_bytes).hexdigest()

    try:
        data = json.loads(raw_bytes)
    except json.JSONDecodeError as e:
        return SnapshotValidationReport(
            snapshot_path=str(path),
            snapshot_sha256=snapshot_sha256,
            n_parameters=0,
            n_conflicts=0,
            errors=[f"snapshot is not valid JSON: {e}"],
        )

    errors: list[str] = []

    parameters = data.get("parameters", []) if isinstance(data, dict) else []
    if not isinstance(data, dict):
        errors.append("top-level snapshot is not a JSON object")
    if not isinstance(parameters, list):
        errors.append("'parameters' is present but is not a JSON array")
        parameters = []

    conflicts = data.get("conflicts", []) if isinstance(data, dict) else []
    if not isinstance(conflicts, list):
        errors.append("'conflicts' is present but is not a JSON array")
        conflicts = []

    known_parameter_ids: set = set()
    seen_seqs: set = set()
    previous_seq: Optional[float] = None

    for idx, p in enumerate(parameters):
        if not isinstance(p, dict):
            errors.append(f"parameters[{idx}] is not a JSON object (got {type(p).__name__})")
            continue

        pid = p.get("parameter_id")
        if pid is None:
            errors.append(f"parameters[{idx}] is missing 'parameter_id'")
            continue
        if not _is_scalar_id(pid):
            errors.append(
                f"parameters[{idx}]: 'parameter_id' must be a plain scalar, "
                f"got {type(pid).__name__}"
            )
            continue
        if pid in known_parameter_ids:
            errors.append(f"parameters[{idx}] ('{pid}'): duplicate parameter_id")
        known_parameter_ids.add(pid)

        seq = p.get("seq")
        if seq is None:
            errors.append(f"parameters[{idx}] ('{pid}'): missing 'seq'")
            continue
        if not _is_scalar_id(seq):
            errors.append(
                f"parameters[{idx}] ('{pid}'): 'seq' must be a plain scalar, "
                f"got {type(seq).__name__}"
            )
            continue
        if seq in seen_seqs:
            errors.append(f"parameters[{idx}] ('{pid}'): duplicate seq={seq}")
        seen_seqs.add(seq)

        # Strict < (not <=) so an equal-seq pair is reported once, by the
        # duplicate-seq check above, not twice. Only compares when both
        # sides are genuinely orderable (seq is validated as a scalar
        # above, but a str seq -- while hashable and thus safe from the
        # crash this guards against -- can't be ordered against a prior
        # numeric seq either; skip the ordering check rather than risk a
        # TypeError from comparing incompatible scalar types).
        if previous_seq is not None and isinstance(seq, (int, float)) and not isinstance(seq, bool) \
                and isinstance(previous_seq, (int, float)) and not isinstance(previous_seq, bool) \
                and seq < previous_seq:
            errors.append(
                f"parameters[{idx}] ('{pid}'): seq={seq} is less than the "
                f"previous entry's seq={previous_seq} -- array order does not "
                f"match seq order; a single-threaded backfill over this array "
                f"would silently diverge from the original insertion/tie order"
            )
        previous_seq = seq

    for idx, p in enumerate(parameters):
        if not isinstance(p, dict):
            continue  # already reported above
        supersedes = p.get("supersedes_parameter_id")
        if supersedes is None:
            continue
        if not _is_scalar_id(supersedes):
            errors.append(
                f"parameters[{idx}] ('{p.get('parameter_id')}'): "
                f"'supersedes_parameter_id' must be a plain scalar, got "
                f"{type(supersedes).__name__}"
            )
            continue
        if supersedes not in known_parameter_ids:
            errors.append(
                f"parameters[{idx}] ('{p.get('parameter_id')}'): "
                f"supersedes_parameter_id '{supersedes}' does not resolve to "
                f"any parameter in this snapshot"
            )

    for idx, c in enumerate(conflicts):
        if not isinstance(c, dict):
            errors.append(f"conflicts[{idx}] is not a JSON object (got {type(c).__name__})")
            continue
        conflict_id = c.get("conflict_id")
        for field_name in ("existing_parameter_id", "incoming_parameter_id"):
            ref = c.get(field_name)
            if ref is None:
                errors.append(f"conflicts[{idx}]: missing '{field_name}'")
            elif not _is_scalar_id(ref):
                errors.append(
                    f"conflicts[{idx}] ('{conflict_id}'): '{field_name}' must be "
                    f"a plain scalar, got {type(ref).__name__}"
                )
            elif ref not in known_parameter_ids:
                errors.append(
                    f"conflicts[{idx}] ('{conflict_id}'): "
                    f"{field_name}='{ref}' does not resolve to any parameter "
                    f"in this snapshot"
                )

    return SnapshotValidationReport(
        snapshot_path=str(path),
        snapshot_sha256=snapshot_sha256,
        n_parameters=len(parameters),
        n_conflicts=len(conflicts),
        errors=errors,
    )
