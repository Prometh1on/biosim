"""
Parameter store — the central database for all kinetic parameters.

Architecture:
    - Parameters are stored in-memory with JSON persistence
    - Every write is append-only — no parameter is ever deleted
    - Conflicts (two parameters for the same entity/kind with different values)
      are stored and flagged, never silently resolved
    - Query results always return the highest-confidence, non-superseded
      parameter unless history is explicitly requested

Thread safety:
    All public methods acquire a read or write lock as appropriate.
    Read operations are concurrent-safe. Writes are serialized.
"""

from __future__ import annotations

import json
import threading
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from .types import (
    ConfidenceTier,
    ExperimentalContext,
    Organism,
    ExperimentalSystem,
    Parameter,
    ParameterKind,
    ParameterValue,
    Provenance,
)


CONFIDENCE_RANK = {
    ConfidenceTier.GOLD: 3,
    ConfidenceTier.SILVER: 2,
    ConfidenceTier.BRONZE: 1,
}

ORGANISM_RANK = {
    Organism.HUMAN: 5,
    Organism.RAT: 3,
    Organism.MOUSE: 3,
    Organism.YEAST: 1,
    Organism.ECOLI: 1,
    Organism.OTHER: 0,
    Organism.UNKNOWN: 0,
}

SYSTEM_RANK = {
    ExperimentalSystem.CLINICAL: 5,
    ExperimentalSystem.IN_VIVO: 4,
    ExperimentalSystem.EX_VIVO: 3,
    ExperimentalSystem.IN_VITRO: 2,
    ExperimentalSystem.IN_SILICO: 1,
    ExperimentalSystem.UNKNOWN: 0,
}


class ConflictRecord:
    """
    Records a conflict between two parameters for the same entity/kind.

    Conflicts are not errors — they represent genuine scientific disagreement.
    They must be resolved by expert curation, not by the software.
    The store continues to function with conflicts present, using the
    highest-priority parameter as the working value until resolved.
    """
    def __init__(
        self,
        existing_id: str,
        incoming_id: str,
        subject_canonical_id: str,
        kind: ParameterKind,
        existing_value: float,
        incoming_value: float,
        fold_difference: float,
    ) -> None:
        self.conflict_id = str(uuid.uuid4())
        self.existing_parameter_id = existing_id
        self.incoming_parameter_id = incoming_id
        self.subject_canonical_id = subject_canonical_id
        self.kind = kind
        self.existing_value = existing_value
        self.incoming_value = incoming_value
        self.fold_difference = fold_difference
        self.detected_at = datetime.now(timezone.utc)
        self.resolved = False
        self.resolution_notes = ""

    def to_dict(self) -> dict:
        return {
            "conflict_id": self.conflict_id,
            "existing_parameter_id": self.existing_parameter_id,
            "incoming_parameter_id": self.incoming_parameter_id,
            "subject_canonical_id": self.subject_canonical_id,
            "kind": self.kind.value,
            "existing_value": self.existing_value,
            "incoming_value": self.incoming_value,
            "fold_difference": self.fold_difference,
            "detected_at": self.detected_at.isoformat(),
            "resolved": self.resolved,
            "resolution_notes": self.resolution_notes,
        }


class ParameterStore:
    """
    Central store for all kinetic parameters in the biosim platform.

    Key guarantees:
        1. Every parameter has full provenance — no orphan values
        2. Conflicts are detected and flagged, never silently resolved
        3. Full version history is preserved — superseded parameters
           are retained for reproducibility
        4. Query priority: GOLD > SILVER > BRONZE, HUMAN > other organisms,
           CLINICAL/IN_VIVO > IN_VITRO > IN_SILICO
    """

    # Conflict threshold — flag as conflict if values differ by more than this factor
    CONFLICT_FOLD_THRESHOLD = 2.0

    def __init__(self) -> None:
        self._parameters: dict[str, Parameter] = {}
        self._conflicts: dict[str, ConflictRecord] = {}
        # Explicit monotonic insertion counter, deliberately mirroring
        # Postgres's `parameter.seq bigint GENERATED ALWAYS AS IDENTITY`
        # (infra/postgres/001_schema.sql) so the Python and SQL tiebreakers
        # are the same concept, not two things that happen to agree.
        # Assigned in add() and load(); used by get_all_for_subject(),
        # query(), and best() to break ties the same way
        # parameter_all_for_subject()'s `ORDER BY ... seq ASC` does.
        #
        # Regression (MASTER_PLAN 3.1): this comment used to claim
        # query()/best() didn't need this tiebreak at all, on the theory
        # that "within one index key the id list is already a subsequence
        # of insertion order." That's true for the in-memory add()-only
        # path, but false the moment a store is round-tripped through
        # save()/load() with a snapshot whose "parameters" array isn't in
        # true seq order (a hand-edited snapshot, a merge of two partial
        # exports, a future non-Python writer -- load()'s own docstring
        # already anticipates exactly this) -- load() appends to
        # _subject_kind_index in ARRAY order, not seq order, so a
        # reordered-but-seq-preserved snapshot silently desyncs the two.
        # Reproduced directly before fixing: saved two tied-priority
        # parameters (added A then B), reversed the JSON array while
        # leaving both parameters' own persisted "seq" untouched, reloaded
        # -- query()/best() returned B (wrong, followed array order) while
        # the already-fixed get_all_for_subject() still correctly returned
        # A (seq order). best() is "the primary method used by the ODE
        # engine to look up parameters" per its own docstring -- a
        # higher-stakes inconsistency than the one already fixed for
        # get_all_for_subject() alone.
        self._insertion_seq: dict[str, int] = {}
        self._next_seq = 0
        # Index: (subject_canonical_id, kind.value, object_canonical_id,
        # reaction_canonical_id) -> list of parameter_ids. Includes object/
        # reaction (both None for the common case of a subject-only
        # parameter like Km or basal abundance) so two parameters for the
        # SAME subject/kind but DIFFERENT objects -- e.g. the insulin
        # receptor's Kd for insulin vs. its Kd for IGF-1, real, clinically
        # relevant cross-reactivity this platform already tracks both
        # ligands for -- are indexed separately rather than treated as
        # rival measurements of the same quantity. Regression-guarded by
        # test_object_scoped_parameters_do_not_conflict_across_different_objects.
        self._subject_kind_index: dict[tuple[str, str, Optional[str], Optional[str]], list[str]] = {}
        self._lock = threading.RLock()

    # ------------------------------------------------------------------
    # Write operations
    # ------------------------------------------------------------------

    def add(self, parameter: Parameter) -> str:
        """
        Add a parameter to the store.

        If a parameter already exists for the same subject/kind combination,
        conflict detection runs automatically. Both parameters are stored;
        a ConflictRecord is created if values diverge beyond the threshold.

        Returns:
            The parameter_id of the stored parameter.
        """
        with self._lock:
            # Ensure parameter_id is set
            if not parameter.parameter_id:
                object.__setattr__(parameter, "parameter_id", str(uuid.uuid4()))

            self._parameters[parameter.parameter_id] = parameter
            self._insertion_seq[parameter.parameter_id] = self._next_seq
            self._next_seq += 1

            # Update subject/kind/object/reaction index
            key = self._index_key(parameter)
            self._subject_kind_index.setdefault(key, []).append(parameter.parameter_id)

            # Conflict detection
            #
            # Regression found via MASTER_PLAN 5.1's own curation work:
            # _index_key() (and _subject_kind_index generally) deliberately
            # does NOT include custom_kind_name -- kept that way so query()/
            # best()/get_all_for_subject() stay fully backward compatible
            # (ParameterLoader.load_pathway_params()'s own p3 mapping calls
            # store.best(insr_id, ParameterKind.CUSTOM) with no
            # custom_kind_name filter at all; making that an exact-match
            # parameter defaulting to None would have silently broken that
            # real production lookup, since p3's own custom_kind_name is
            # "insulin_action_p3", never None). But that same coarseness
            # meant CONFLICT DETECTION treated every CUSTOM-kind parameter
            # sharing a bare subject/object/reaction as a rival measurement
            # of the same quantity -- confirmed directly: registering 4
            # exogenous-insulin PK constants (all subject="insulin",
            # kind=CUSTOM, no object/reaction, but 4 DIFFERENT
            # custom_kind_names for 4 semantically unrelated physical
            # quantities) produced 5 nonsensical conflicts, fold differences
            # up to ~61,500x. Exactly the "two distinct physical quantities
            # treated as rival measurements of the same one" failure mode
            # query()'s own docstring already warns object_canonical_id/
            # reaction_canonical_id exist to prevent -- extended here to
            # the one scope dimension it hadn't been applied to yet, but
            # ONLY for conflict detection (query()'s own lookup semantics
            # are untouched).
            existing_ids = [
                pid for pid in self._subject_kind_index[key]
                if pid != parameter.parameter_id
                and not self._parameters[pid].is_superseded
                and (
                    parameter.kind != ParameterKind.CUSTOM
                    or self._parameters[pid].custom_kind_name == parameter.custom_kind_name
                )
            ]

            for existing_id in existing_ids:
                existing = self._parameters[existing_id]
                self._detect_conflict(existing, parameter)

            return parameter.parameter_id

    def supersede(
        self,
        old_parameter_id: str,
        new_parameter: Parameter,
    ) -> str:
        """
        Mark an existing parameter as superseded by a new one.

        Use this when a better measurement is published. The old parameter
        is retained in history but excluded from active queries.
        """
        with self._lock:
            if old_parameter_id not in self._parameters:
                raise KeyError(f"Parameter {old_parameter_id} not found.")

            old = self._parameters[old_parameter_id]
            old.supersede(new_parameter.parameter_id)

            # Link the new parameter back to the old
            object.__setattr__(
                new_parameter,
                "supersedes_parameter_id",
                old_parameter_id
            )

            return self.add(new_parameter)

    def resolve_conflict(self, conflict_id: str, notes: str) -> None:
        """Mark a conflict as resolved with curator notes."""
        with self._lock:
            if conflict_id not in self._conflicts:
                raise KeyError(f"Conflict {conflict_id} not found.")
            self._conflicts[conflict_id].resolved = True
            self._conflicts[conflict_id].resolution_notes = notes

    # ------------------------------------------------------------------
    # Query operations
    # ------------------------------------------------------------------

    def get(self, parameter_id: str) -> Parameter:
        """Retrieve a specific parameter by ID."""
        with self._lock:
            if parameter_id not in self._parameters:
                raise KeyError(f"Parameter {parameter_id} not found.")
            return self._parameters[parameter_id]

    def query(
        self,
        subject_canonical_id: str,
        kind: ParameterKind,
        min_confidence: ConfidenceTier = ConfidenceTier.BRONZE,
        organism: Optional[Organism] = None,
        include_superseded: bool = False,
        object_canonical_id: Optional[str] = None,
        reaction_canonical_id: Optional[str] = None,
    ) -> list[Parameter]:
        """
        Query parameters for a subject entity and kind.

        Results are sorted by priority:
            1. Confidence tier (GOLD > SILVER > BRONZE)
            2. Organism relevance (HUMAN > others)
            3. Experimental system (CLINICAL > IN_VIVO > IN_VITRO > IN_SILICO)

        Args:
            subject_canonical_id: The canonical ID string of the subject entity
            kind: The type of parameter to retrieve
            min_confidence: Minimum confidence tier to include
            organism: Filter to a specific organism (None = all)
            include_superseded: Whether to include superseded parameters
            object_canonical_id: Match only parameters scoped to this exact
                object entity (e.g. a binding partner for a KD). Defaults
                to None, matching only subject-only parameters (those with
                no object) -- NOT a wildcard over all objects, so a Kd
                stored `.for_object(insulin_id)` requires passing
                `object_canonical_id=insulin_id` explicitly to retrieve it.
                This is deliberate: silently mixing a subject's parameters
                across different objects (e.g. its Kd for insulin vs. its
                Kd for IGF-1) would treat two distinct physical quantities
                as rival measurements of the same one.
            reaction_canonical_id: Same matching behaviour as
                object_canonical_id, for reaction-scoped parameters.

        Returns:
            List of matching parameters, sorted by priority (best first)
        """
        with self._lock:
            key = (subject_canonical_id, kind.value, object_canonical_id, reaction_canonical_id)
            ids = self._subject_kind_index.get(key, [])

            results = []
            for pid in ids:
                p = self._parameters[pid]
                if not include_superseded and p.is_superseded:
                    continue
                if CONFIDENCE_RANK[p.confidence_tier] < CONFIDENCE_RANK[min_confidence]:
                    continue
                if organism and p.experimental_context.organism != organism:
                    continue
                results.append(p)

            # Explicit self._insertion_seq tiebreak, not just
            # `key=self._priority_score, reverse=True` -- see
            # self._insertion_seq's own docstring comment in __init__ for
            # the reproduced bug this fixes (a reordered-but-seq-
            # preserved snapshot silently desyncing query()/best() from
            # get_all_for_subject()). Same negate-and-sort-ascending
            # pattern get_all_for_subject() already uses, for the same
            # reason: it doesn't rely on any sort's tie-stability
            # semantics or on `ids`' own list order, only on the
            # explicit, always-correct seq value.
            return sorted(
                results,
                key=lambda p: (
                    -CONFIDENCE_RANK[p.confidence_tier],
                    -ORGANISM_RANK[p.experimental_context.organism],
                    -SYSTEM_RANK[p.experimental_context.system],
                    self._insertion_seq[p.parameter_id],
                ),
            )

    def best(
        self,
        subject_canonical_id: str,
        kind: ParameterKind,
        min_confidence: ConfidenceTier = ConfidenceTier.BRONZE,
        organism: Optional[Organism] = None,
        object_canonical_id: Optional[str] = None,
        reaction_canonical_id: Optional[str] = None,
    ) -> Optional[Parameter]:
        """
        Return the single best parameter for a subject/kind combination.

        Returns None if no parameter exists above the minimum confidence threshold.
        This is the primary method used by the ODE engine to look up parameters.
        See query()'s docstring for object_canonical_id/reaction_canonical_id's
        exact-match (not wildcard) semantics.
        """
        results = self.query(
            subject_canonical_id=subject_canonical_id,
            kind=kind,
            min_confidence=min_confidence,
            organism=organism,
            object_canonical_id=object_canonical_id,
            reaction_canonical_id=reaction_canonical_id,
        )
        return results[0] if results else None

    def get_all_for_subject(
        self,
        subject_canonical_id: str,
        include_superseded: bool = False,
    ) -> dict[str, list[Parameter]]:
        """
        Return all parameters for a subject entity, grouped by kind.
        Useful for building a complete kinetic profile for an entity.

        A subject+kind spanning multiple object/reaction scopes (e.g. the
        insulin receptor's Kd for insulin AND its Kd for IGF-1) returns
        ALL of them under one kind key -- this method is a profile/
        inspection view, not a value-selection path, so it deliberately
        does not disambiguate by scope the way query()/best() do. Callers
        needing to distinguish scopes read object_canonical_id/
        reaction_canonical_id off the returned Parameters directly.

        Within a kind, results are ordered by priority (confidence tier,
        then organism, then experimental system), with ties broken by
        global insertion order -- matching
        infra/postgres/002_query_api.sql's parameter_all_for_subject()
        `ORDER BY v.kind, ct.rank DESC, o.rank DESC, es.rank DESC,
        v.seq ASC` exactly, so the two are equivalent modulo shape (SQL
        returns a flat ordered SETOF; Python groups by kind while
        preserving that same per-kind order).
        """
        with self._lock:
            result: dict[str, list[Parameter]] = {}
            for (subj, kind_str, _obj, _reaction), ids in self._subject_kind_index.items():
                if subj != subject_canonical_id:
                    continue
                for pid in ids:
                    p = self._parameters[pid]
                    if not include_superseded and p.is_superseded:
                        continue
                    result.setdefault(kind_str, []).append(p)
            for kind_str, params in result.items():
                # NOT `key=self._priority_score, reverse=True` -- that alone
                # would not be wrong on tie-stability grounds (Python's sort
                # is stable regardless of reverse=True: ties keep their
                # relative INPUT order, verified directly). The real reason
                # an explicit seq tiebreak is required is upstream of the
                # sort: `result` above is built by iterating
                # self._subject_kind_index.items(), grouped by when each
                # distinct scope key was first created, with each key's own
                # id-list in its own append order -- not a flat merge in
                # true global add() order. Interleaved insertion into
                # different scopes (add scope A, then scope B, then scope A
                # again) would put both A entries ahead of B's, even though
                # B was added in between. Negating the rank components and
                # sorting ascending on the explicit self._insertion_seq
                # (not reliant on iteration order or sort-stability
                # semantics) restores true global insertion order regardless.
                result[kind_str] = sorted(
                    params,
                    key=lambda p: (
                        -CONFIDENCE_RANK[p.confidence_tier],
                        -ORGANISM_RANK[p.experimental_context.organism],
                        -SYSTEM_RANK[p.experimental_context.system],
                        self._insertion_seq[p.parameter_id],
                    ),
                )
            return result

    def unresolved_conflicts(self) -> list[ConflictRecord]:
        """Return all unresolved conflicts, sorted by fold difference (largest first)."""
        with self._lock:
            return sorted(
                [c for c in self._conflicts.values() if not c.resolved],
                key=lambda c: c.fold_difference,
                reverse=True,
            )

    def stats(self) -> dict:
        """Return store statistics for monitoring."""
        with self._lock:
            total = len(self._parameters)
            superseded = sum(1 for p in self._parameters.values() if p.is_superseded)
            by_tier = {}
            by_kind = {}
            for p in self._parameters.values():
                if p.is_superseded:
                    continue
                tier = p.confidence_tier.value
                kind = p.kind.value
                by_tier[tier] = by_tier.get(tier, 0) + 1
                by_kind[kind] = by_kind.get(kind, 0) + 1
            return {
                "total_parameters": total,
                "active_parameters": total - superseded,
                "superseded_parameters": superseded,
                "unresolved_conflicts": len([c for c in self._conflicts.values() if not c.resolved]),
                "total_conflicts": len(self._conflicts),
                "by_confidence_tier": by_tier,
                "by_parameter_kind": by_kind,
            }

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def save(self, path: Path) -> None:
        """Persist the full store to a JSON file."""
        with self._lock:
            data = {
                "parameters": [
                    self._serialize_parameter(p, self._insertion_seq[p.parameter_id])
                    for p in self._parameters.values()
                ],
                "conflicts": [c.to_dict() for c in self._conflicts.values()],
            }
            path.write_text(json.dumps(data, indent=2, default=str))

    @classmethod
    def load(cls, path: Path) -> "ParameterStore":
        """Load a store from a JSON file."""
        store = cls()
        data = json.loads(path.read_text())
        for p_dict in data.get("parameters", []):
            p = cls._deserialize_parameter(p_dict)
            store._parameters[p.parameter_id] = p
            # Prefer the persisted seq (robust to the parameters array
            # ever not being in true original-insertion order -- a
            # hand-edited snapshot, a merge of two partial exports, a
            # future non-Python writer); fall back to array order only
            # for pre-fix snapshots that never wrote "seq" at all.
            seq = p_dict.get("seq")
            if seq is None:
                seq = store._next_seq
            store._insertion_seq[p.parameter_id] = seq
            store._next_seq = max(store._next_seq, seq + 1)
            key = cls._index_key(p)
            store._subject_kind_index.setdefault(key, []).append(p.parameter_id)
        for c_dict in data.get("conflicts", []):
            c = ConflictRecord(
                existing_id=c_dict["existing_parameter_id"],
                incoming_id=c_dict["incoming_parameter_id"],
                subject_canonical_id=c_dict["subject_canonical_id"],
                kind=ParameterKind(c_dict["kind"]),
                existing_value=c_dict["existing_value"],
                incoming_value=c_dict["incoming_value"],
                fold_difference=c_dict["fold_difference"],
            )
            c.conflict_id = c_dict["conflict_id"]
            c.resolved = c_dict["resolved"]
            c.resolution_notes = c_dict["resolution_notes"]
            # Round-trip detected_at if present -- same
            # default_factory-regenerates-on-load bug as
            # Parameter.created_at/Provenance.extracted_at, just for
            # ConflictRecord's own timestamp; graceful fallback for
            # pre-fix snapshots that never wrote it.
            detected_at_str = c_dict.get("detected_at")
            if detected_at_str is not None:
                c.detected_at = datetime.fromisoformat(detected_at_str)
            store._conflicts[c.conflict_id] = c
        return store

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _index_key(p: Parameter) -> tuple[str, str, Optional[str], Optional[str]]:
        return (p.subject_canonical_id, p.kind.value, p.object_canonical_id, p.reaction_canonical_id)

    def _detect_conflict(self, existing: Parameter, incoming: Parameter) -> None:
        """Detect and register a conflict between two parameter values."""
        if existing.value.value == 0 or incoming.value.value == 0:
            return
        fold = max(
            existing.value.value / incoming.value.value,
            incoming.value.value / existing.value.value,
        )
        if fold >= self.CONFLICT_FOLD_THRESHOLD:
            conflict = ConflictRecord(
                existing_id=existing.parameter_id,
                incoming_id=incoming.parameter_id,
                subject_canonical_id=existing.subject_canonical_id,
                kind=existing.kind,
                existing_value=existing.value.value,
                incoming_value=incoming.value.value,
                fold_difference=fold,
            )
            self._conflicts[conflict.conflict_id] = conflict

    @staticmethod
    def _serialize_parameter(p: Parameter, seq: int) -> dict:
        return {
            "parameter_id": p.parameter_id,
            "seq": seq,
            "kind": p.kind.value,
            "custom_kind_name": p.custom_kind_name,
            "value": {
                "value": p.value.value,
                "unit": p.value.unit,
                "uncertainty_lower": p.value.uncertainty_lower,
                "uncertainty_upper": p.value.uncertainty_upper,
                "uncertainty_type": p.value.uncertainty_type,
            },
            "subject_canonical_id": p.subject_canonical_id,
            "object_canonical_id": p.object_canonical_id,
            "reaction_canonical_id": p.reaction_canonical_id,
            "confidence_tier": p.confidence_tier.value,
            "provenance": {
                "pubmed_id": p.provenance.pubmed_id,
                "doi": p.provenance.doi,
                "arxiv_id": p.provenance.arxiv_id,
                "source_description": p.provenance.source_description,
                "extracted_by": p.provenance.extracted_by,
                "extracted_at": p.provenance.extracted_at.isoformat(),
                "figure_or_table": p.provenance.figure_or_table,
                "is_primary_source": p.provenance.is_primary_source,
                "verification_status": p.provenance.verification_status,
                "verification_quote": p.provenance.verification_quote,
            },
            "experimental_context": {
                "organism": p.experimental_context.organism.value,
                "system": p.experimental_context.system.value,
                "cell_type": p.experimental_context.cell_type,
                "tissue": p.experimental_context.tissue,
                "temperature_celsius": p.experimental_context.temperature_celsius,
                "ph": p.experimental_context.ph,
                "notes": p.experimental_context.notes,
            },
            "supersedes_parameter_id": p.supersedes_parameter_id,
            "curator_notes": p.curator_notes,
            "created_at": p.created_at.isoformat(),
            "is_superseded": p.is_superseded,
        }

    @staticmethod
    def _deserialize_parameter(d: dict) -> Parameter:
        value = ParameterValue(
            value=d["value"]["value"],
            unit=d["value"]["unit"],
            uncertainty_lower=d["value"]["uncertainty_lower"],
            uncertainty_upper=d["value"]["uncertainty_upper"],
            uncertainty_type=d["value"]["uncertainty_type"],
        )
        provenance_kwargs = dict(
            pubmed_id=d["provenance"]["pubmed_id"],
            doi=d["provenance"]["doi"],
            arxiv_id=d["provenance"]["arxiv_id"],
            source_description=d["provenance"]["source_description"],
            extracted_by=d["provenance"]["extracted_by"],
            figure_or_table=d["provenance"]["figure_or_table"],
            verification_status=d["provenance"].get("verification_status", "UNCHECKED"),
            verification_quote=d["provenance"].get("verification_quote", ""),
            is_primary_source=d["provenance"]["is_primary_source"],
        )
        # Round-trip extracted_at if present (was silently regenerated via
        # default_factory=datetime.now(timezone.utc) before this fix --
        # graceful .get() fallback for pre-fix snapshots that never wrote it).
        extracted_at_str = d["provenance"].get("extracted_at")
        if extracted_at_str is not None:
            provenance_kwargs["extracted_at"] = datetime.fromisoformat(extracted_at_str)
        provenance = Provenance(**provenance_kwargs)
        context = ExperimentalContext(
            organism=Organism(d["experimental_context"]["organism"]),
            system=ExperimentalSystem(d["experimental_context"]["system"]),
            cell_type=d["experimental_context"]["cell_type"],
            tissue=d["experimental_context"]["tissue"],
            temperature_celsius=d["experimental_context"]["temperature_celsius"],
            ph=d["experimental_context"]["ph"],
            notes=d["experimental_context"]["notes"],
        )
        parameter_kwargs = dict(
            parameter_id=d["parameter_id"],
            kind=ParameterKind(d["kind"]),
            value=value,
            subject_canonical_id=d["subject_canonical_id"],
            provenance=provenance,
            experimental_context=context,
            confidence_tier=ConfidenceTier(d["confidence_tier"]),
            object_canonical_id=d.get("object_canonical_id"),
            reaction_canonical_id=d.get("reaction_canonical_id"),
            custom_kind_name=d.get("custom_kind_name"),
            supersedes_parameter_id=d.get("supersedes_parameter_id"),
            curator_notes=d.get("curator_notes", ""),
            is_superseded=d.get("is_superseded", False),
        )
        # Round-trip created_at if present -- same rationale as
        # provenance.extracted_at above (was silently regenerated via
        # default_factory=datetime.now(timezone.utc) before this fix).
        created_at_str = d.get("created_at")
        if created_at_str is not None:
            parameter_kwargs["created_at"] = datetime.fromisoformat(created_at_str)
        return Parameter(**parameter_kwargs)
