"""
Canonical ID resolver registry.

This is the single entry point for all ID resolution in the platform.
No other module is permitted to perform ontology lookups directly.

Resolution flow:
    ExternalID → lookup in cross-reference index → CanonicalID
    CanonicalID → lookup in registry → ResolutionRecord

Thread safety:
    The registry is read-only after initialization. Writes occur only
    during data ingestion pipelines, which acquire a write lock.
    All read operations are safe for concurrent use.

Usage:
    registry = IDRegistry.from_mapping_files(Path("data/mappings"))
    canonical = registry.resolve(ExternalID(OntologySource.UNIPROT, "P06213"))
    record = registry.get_record(canonical)
"""

from __future__ import annotations

import json
import threading
from pathlib import Path
from typing import Optional

from .types import (
    CanonicalID,
    EntityType,
    ExternalID,
    OntologySource,
    RESOLUTION_PRIORITY,
    ResolutionRecord,
)


class ResolutionError(Exception):
    """Raised when an external ID cannot be resolved to a canonical ID."""
    pass


class AmbiguousResolutionError(ResolutionError):
    """
    Raised when an external ID maps to multiple canonical IDs.

    This indicates a data integrity problem in the cross-reference tables
    and must be resolved by curation, not by the resolver silently picking one.
    """
    def __init__(self, external_id: ExternalID, candidates: list[CanonicalID]) -> None:
        self.external_id = external_id
        self.candidates = candidates
        super().__init__(
            f"Ambiguous resolution: {external_id} maps to multiple canonical IDs: "
            f"{[str(c) for c in candidates]}. "
            f"This requires manual curation to resolve."
        )


class IDRegistry:
    """
    Thread-safe registry mapping external ontology IDs to canonical platform IDs.

    Internal indexes:
        _canonical_index: CanonicalID → ResolutionRecord
        _xref_index: (OntologySource, accession) → list[CanonicalID]
        _deprecated_index: old CanonicalID string → current CanonicalID
        _name_index: lowercase name → list[CanonicalID] (fuzzy search support)
    """

    def __init__(self) -> None:
        self._canonical_index: dict[str, ResolutionRecord] = {}
        self._xref_index: dict[tuple[str, str], list[str]] = {}
        self._deprecated_index: dict[str, str] = {}
        self._name_index: dict[str, list[str]] = {}
        self._lock = threading.RLock()
        self._next_ids: dict[EntityType, int] = {t: 1 for t in EntityType}

    # ------------------------------------------------------------------
    # Registration (used by data ingestion pipelines only)
    # ------------------------------------------------------------------

    def register(
        self,
        name: str,
        entity_type: EntityType,
        description: str = "",
        cross_references: Optional[dict[OntologySource, list[str]]] = None,
        deprecated_canonical_ids: Optional[list[CanonicalID]] = None,
        force_numeric_id: Optional[int] = None,
    ) -> CanonicalID:
        """
        Register a new biological entity and return its canonical ID.

        Args:
            name: Human-readable name (e.g. "Insulin receptor", "D-Glucose")
            entity_type: The type of biological entity
            description: Optional longer description
            cross_references: Dict mapping ontology sources to lists of accessions
            deprecated_canonical_ids: Old canonical IDs that now map to this entity
            force_numeric_id: For deterministic seeding only — do not use in production

        Returns:
            The newly assigned CanonicalID

        Raises:
            ValueError: If a cross-reference accession is already registered
                        to a different canonical ID (data integrity violation)
        """
        with self._lock:
            xrefs = {k: list(v) for k, v in (cross_references or {}).items()}
            deprecated = deprecated_canonical_ids or []

            # Regression: a merge used to only re-point the cross-references
            # EXPLICITLY RESTATED in this call's `cross_references` argument
            # -- any accession the deprecated entity owned under an ontology
            # source the caller didn't happen to restate was silently left
            # pointing at the now-stale canonical ID forever. resolve() for
            # that accession would then return a deprecated ID with no
            # warning at all (resolve() doesn't consult _deprecated_index,
            # only get_record() does) -- confirmed directly: merge an entity
            # that had both HMDB and KEGG cross-refs into a new one that
            # only restates HMDB, and resolve()-ing the KEGG accession still
            # returns the deprecated canonical ID. Fixed by auto-inheriting
            # every cross-reference the deprecated entity(ies) actually had,
            # merged in before the collision check below -- callers no
            # longer need to manually re-list everything the old entity
            # owned just to keep it resolvable.
            for old_id in deprecated:
                old_record = self._canonical_index.get(str(old_id))
                if old_record is None:
                    continue
                for source, accessions in old_record.cross_references.items():
                    existing = xrefs.setdefault(source, [])
                    for acc in accessions:
                        if acc not in existing:
                            existing.append(acc)

            # Regression: the collision check below used to fire
            # unconditionally, with no exemption for accessions belonging
            # to an entity being deprecated BY THIS SAME CALL -- but
            # deprecated_canonical_ids is the documented mechanism for
            # exactly that workflow (merge two entities that legitimately
            # share a real-world identifier into one consolidated record).
            # Since register() is the only way to add data (no separate
            # update()/merge() method exists), this made the merge
            # workflow entirely non-functional: any attempt to register a
            # merged entity carrying a cross-reference already owned by
            # the entity it's deprecating raised ValueError unconditionally.
            # Confirmed live and reproduced with a direct register() call.
            # Every existing test deprecated a canonical ID that was never
            # actually registered, which is why this was never caught.
            deprecated_strs = {str(old_id) for old_id in deprecated}

            # Check for cross-reference collisions before registering
            # anything -- exempt accessions owned ENTIRELY by entities
            # being deprecated by this call; a genuine third-party owner
            # (anything not in deprecated_strs) still collides as before.
            for source, accessions in xrefs.items():
                for acc in accessions:
                    key = (source.value, acc.upper())
                    non_deprecated_owners = [
                        c for c in self._xref_index.get(key, []) if c not in deprecated_strs
                    ]
                    if non_deprecated_owners:
                        raise ValueError(
                            f"Cross-reference collision: {source.value}:{acc} is already "
                            f"registered to {non_deprecated_owners}. "
                            f"Resolve this in the curation layer before registering."
                        )

            # Assign canonical ID
            if force_numeric_id is not None:
                numeric_id = force_numeric_id
                if force_numeric_id >= self._next_ids[entity_type]:
                    self._next_ids[entity_type] = force_numeric_id + 1
            else:
                numeric_id = self._next_ids[entity_type]
                self._next_ids[entity_type] += 1

            canonical = CanonicalID(entity_type=entity_type, numeric_id=numeric_id)
            canonical_str = str(canonical)

            # Regression, code-reviewer finding: force_numeric_id was
            # never validated against ids already present, even though
            # both real deserialization entry points -- from_dict() and
            # from_mapping_files(), the class's own documented production
            # usage -- call register(force_numeric_id=...) unconditionally
            # for every loaded record, not just in tests. A corrupted or
            # hand-edited snapshot (a "never trust external data" boundary
            # per this project's own security rules) could produce a real
            # self-loop (deprecated_canonical_ids=[X], force_numeric_id=
            # X's own numeric_id -- confirmed reachable through the public
            # API) that _resolve_deprecation_chain()'s cycle guard exists
            # to catch, or -- more subtly -- reconstruct a canonical_str
            # that ALSO exists as a stale _deprecated_index key, which
            # get_record() checks before _canonical_index and would then
            # incorrectly report a live, just-loaded entity as deprecated.
            # Fail fast at registration instead of leaving either failure
            # mode to surface later, confusingly, at lookup time.
            if canonical_str in self._canonical_index:
                raise ValueError(
                    f"force_numeric_id collision: {canonical_str} is already a "
                    f"live, registered canonical ID. Cannot register it again."
                )
            if canonical_str in self._deprecated_index:
                raise ValueError(
                    f"force_numeric_id collision: {canonical_str} is already a "
                    f"deprecated canonical ID (superseded by "
                    f"{self._deprecated_index[canonical_str]}). Assigning it as a "
                    f"live id again would make it ambiguous which one get_record() "
                    f"should report -- and if this new registration also deprecates "
                    f"the SAME id via deprecated_canonical_ids, that is a self-loop."
                )

            # Build and store the resolution record
            record = ResolutionRecord(
                canonical_id=canonical,
                name=name,
                description=description,
                cross_references={k: list(v) for k, v in xrefs.items()},
                deprecated_canonical_ids=list(deprecated),
            )
            self._canonical_index[canonical_str] = record

            # Build cross-reference index -- re-point any accession
            # previously owned by a now-deprecated entity to THIS canonical
            # ID instead of appending a second, stale owner alongside it.
            # Without this, even after exempting the collision check above,
            # resolve() would return the deprecated entity as a second
            # candidate for the shared accession and start raising
            # AmbiguousResolutionError for what should now be an
            # unambiguous, correctly-merged lookup.
            for source, accessions in xrefs.items():
                for acc in accessions:
                    key = (source.value, acc.upper())
                    remaining = [
                        c for c in self._xref_index.get(key, []) if c not in deprecated_strs
                    ]
                    remaining.append(canonical_str)
                    self._xref_index[key] = remaining

            # Build deprecated ID index, and remove the deprecated entity's
            # own stale entry from _canonical_index -- get_record() already
            # checks _deprecated_index first and never falls through to a
            # deprecated entry here, so this was pure dead weight, but
            # stats() iterates _canonical_index directly and was silently
            # double-counting every merged-away entity as still separately
            # live (confirmed directly: merging one entity into another
            # left total_entities at 2, not 1).
            for old_id in deprecated:
                old_str = str(old_id)
                self._deprecated_index[old_str] = canonical_str
                self._canonical_index.pop(old_str, None)

                # Regression: _name_index was never touched on merge --
                # search_by_name() still returned the now-deprecated
                # old_id, and get_record() on that result then raised
                # ResolutionError, even though the entity itself is still
                # perfectly findable by name, just under a new canonical
                # id. Confirmed directly: registering "Old Name" then
                # merging it into a new entity left
                # search_by_name("old name") pointing at the dead id
                # forever. Repoint (not remove) every name-index entry
                # that referenced old_str to the new canonical_str, so a
                # search by the ENTITY'S OLD NAME still resolves --
                # deduping in case the new entity's own name registration
                # below would otherwise create a second identical entry.
                for ids in self._name_index.values():
                    if old_str in ids:
                        ids[:] = list(dict.fromkeys(
                            canonical_str if cid == old_str else cid for cid in ids
                        ))

            # Build name index (lowercase, for fuzzy search). Guarded
            # against a duplicate, not a bare append: if the merging
            # entity's own name happens to equal a deprecated old_id's
            # name, the repoint step above already added canonical_str
            # under this exact name_key.
            name_key = name.lower().strip()
            ids = self._name_index.setdefault(name_key, [])
            if canonical_str not in ids:
                ids.append(canonical_str)

            return canonical

    # ------------------------------------------------------------------
    # Resolution (the public API — all other modules use only these)
    # ------------------------------------------------------------------

    def resolve(self, external_id: ExternalID) -> CanonicalID:
        """
        Resolve an external ontology ID to a canonical platform ID.

        Args:
            external_id: An ExternalID from any supported ontology source

        Returns:
            The canonical ID for this entity

        Raises:
            ResolutionError: If the external ID is not found in the registry
            AmbiguousResolutionError: If the external ID maps to multiple
                                      canonical IDs (data integrity violation)
        """
        key = (external_id.source.value, external_id.accession.upper())
        candidates = self._xref_index.get(key, [])

        if not candidates:
            raise ResolutionError(
                f"Cannot resolve {external_id}. "
                f"This ID is not registered in the platform. "
                f"If this is a valid biological entity, add it via the curation pipeline."
            )

        if len(candidates) > 1:
            raise AmbiguousResolutionError(
                external_id,
                [CanonicalID.from_string(c) for c in candidates]
            )

        return CanonicalID.from_string(candidates[0])

    def resolve_any(self, accession: str) -> list[CanonicalID]:
        """
        Attempt to resolve an accession string without knowing its source.
        Searches all ontology sources in priority order.

        Returns all matching canonical IDs (may be from different sources).
        Empty list if not found anywhere.

        Use this only for search/discovery — prefer resolve() when the
        source is known.
        """
        results = []
        seen = set()
        for source in RESOLUTION_PRIORITY:
            key = (source.value, accession.upper())
            for canonical_str in self._xref_index.get(key, []):
                if canonical_str not in seen:
                    seen.add(canonical_str)
                    results.append(CanonicalID.from_string(canonical_str))
        return results

    def _resolve_deprecation_chain(self, canonical_str: str) -> str:
        """
        Follow _deprecated_index to the final, currently-live canonical ID.

        Regression: _deprecated_index[old] is set ONCE, at the moment
        `old` is deprecated, and is never updated if the entity it now
        points at is LATER deprecated too (a transitive merge chain,
        A->B->C). get_record(A) used to report "current canonical ID is
        B" as if that were final -- but B could itself be deprecated,
        landing a caller who trusted the hint on another dead id.
        Confirmed directly: registering B with deprecated=[A], then
        registering C with deprecated=[B], left _deprecated_index[A] ==
        "B" forever, never updated to "C". This walks the chain fully so
        the reported "current" id is always genuinely live.
        """
        seen = set()
        current = canonical_str
        while current in self._deprecated_index:
            if current in seen:
                raise ResolutionError(
                    f"Circular deprecation chain detected at {current} -- "
                    f"this indicates registry data corruption, not normal usage."
                )
            seen.add(current)
            current = self._deprecated_index[current]
        return current

    def get_record(self, canonical_id: CanonicalID) -> ResolutionRecord:
        """
        Retrieve the full resolution record for a canonical ID.

        Raises:
            ResolutionError: If the canonical ID is not found.
                             Also checks deprecated IDs and raises with
                             a migration hint if found there instead.
        """
        canonical_str = str(canonical_id)

        # Check deprecated index first and give a helpful error
        if canonical_str in self._deprecated_index:
            current = self._resolve_deprecation_chain(canonical_str)
            raise ResolutionError(
                f"{canonical_str} is a deprecated canonical ID. "
                f"The current canonical ID is {current}. "
                f"Update your references — deprecated IDs will be removed in a future release."
            )

        record = self._canonical_index.get(canonical_str)
        if record is None:
            raise ResolutionError(f"No record found for canonical ID {canonical_str}.")
        return record

    def cross_reference(
        self,
        canonical_id: CanonicalID,
        target_source: OntologySource,
    ) -> Optional[str]:
        """
        Get the accession for a canonical ID in a specific external ontology.

        Returns None if no cross-reference exists for this source.
        Use this to translate platform IDs back to researcher-facing ontology IDs.
        """
        record = self.get_record(canonical_id)
        return record.get_external_id(target_source)

    def search_by_name(self, name: str) -> list[CanonicalID]:
        """
        Search for entities by human-readable name (case-insensitive, exact match).
        For fuzzy/partial search, use the search service layer — not this method.
        """
        name_key = name.lower().strip()
        return [
            CanonicalID.from_string(c)
            for c in self._name_index.get(name_key, [])
        ]

    def is_registered(self, external_id: ExternalID) -> bool:
        """Return True if this external ID is in the registry."""
        key = (external_id.source.value, external_id.accession.upper())
        return key in self._xref_index

    # ------------------------------------------------------------------
    # Serialization (for persistence and migration)
    # ------------------------------------------------------------------

    def to_dict(self) -> dict:
        """Serialize the full registry to a JSON-compatible dict."""
        records = []
        for canonical_str, record in self._canonical_index.items():
            records.append({
                "canonical_id": canonical_str,
                "name": record.name,
                "description": record.description,
                "cross_references": {
                    source.value: accessions
                    for source, accessions in record.cross_references.items()
                },
                "deprecated_canonical_ids": [
                    str(d) for d in record.deprecated_canonical_ids
                ],
            })
        return {
            "records": records,
            "next_ids": {k.value: v for k, v in self._next_ids.items()},
            # Regression: a transitive merge chain (A deprecated into B,
            # then B itself later deprecated into C) used to lose the
            # A->B hop entirely on round-trip. Each record's own
            # deprecated_canonical_ids only captures the deprecations
            # applied AT THAT REGISTRATION -- C's record lists only
            # [B], never [A] -- and B's own record (which held [A]) is
            # removed from _canonical_index once merged into C, so
            # `records` above never re-derives it either. Serializing
            # the raw, flat _deprecated_index directly preserves every
            # hop regardless of chain depth, since it's appended to
            # (never overwritten) as merges happen -- see
            # _resolve_deprecation_chain(), which is what actually
            # walks it back to the live id.
            "deprecated_index": dict(self._deprecated_index),
        }

    @classmethod
    def from_dict(cls, data: dict) -> "IDRegistry":
        """Reconstruct a registry from a serialized dict."""
        registry = cls()
        for r in data.get("records", []):
            canonical = CanonicalID.from_string(r["canonical_id"])
            xrefs = {
                OntologySource(source): accessions
                for source, accessions in r.get("cross_references", {}).items()
            }
            deprecated = [
                CanonicalID.from_string(d)
                for d in r.get("deprecated_canonical_ids", [])
            ]
            registry.register(
                name=r["name"],
                entity_type=canonical.entity_type,
                description=r.get("description", ""),
                cross_references=xrefs,
                deprecated_canonical_ids=deprecated,
                force_numeric_id=canonical.numeric_id,
            )
        for entity_type_str, next_id in data.get("next_ids", {}).items():
            registry._next_ids[EntityType(entity_type_str)] = next_id
        # Restore the full, possibly multi-hop deprecation chain directly
        # -- register() above only reconstructs the DIRECT (one-hop)
        # deprecations each surviving record's own deprecated_canonical_ids
        # still lists, which silently drops earlier hops in a transitive
        # merge chain (see to_dict()'s comment). `.update()` rather than
        # assignment so this is additive on top of whatever register()
        # already populated, and `data.get(..., {})` defaults to `{}` for
        # data serialized before this fix, so old snapshots still load
        # (just without the deeper chain restored -- a graceful
        # degradation, not a crash).
        registry._deprecated_index.update(data.get("deprecated_index", {}))
        return registry

    @classmethod
    def from_mapping_files(cls, mapping_dir: Path) -> "IDRegistry":
        """
        Load a registry from JSON mapping files in a directory.

        Expects files named by entity type:
            proteins.json, metabolites.json, pathways.json, genes.json

        Each file is a list of registration records in the format
        produced by to_dict()["records"] -- and, per the same convention
        as to_dict()'s own top-level output, MAY optionally also carry a
        "deprecated_index" key for that file's entities (see the
        MASTER_PLAN 3.3 regression note below); omitted entirely for
        files that predate this fix or were produced by something
        writing only the bare records list, which still load correctly
        (just without any deeper merge-chain history restored, the same
        graceful degradation from_dict() already established).
        """
        registry = cls()
        for filename in sorted(mapping_dir.glob("*.json")):
            with open(filename) as f:
                data = json.load(f)
            for r in data.get("records", []):
                canonical = CanonicalID.from_string(r["canonical_id"])
                xrefs = {
                    OntologySource(source): accessions
                    for source, accessions in r.get("cross_references", {}).items()
                }
                deprecated = [
                    CanonicalID.from_string(d)
                    for d in r.get("deprecated_canonical_ids", [])
                ]
                registry.register(
                    name=r["name"],
                    entity_type=canonical.entity_type,
                    description=r.get("description", ""),
                    cross_references=xrefs,
                    deprecated_canonical_ids=deprecated,
                    force_numeric_id=canonical.numeric_id,
                )
            # Regression (MASTER_PLAN 3.3): this method's own docstring
            # already names it, alongside from_dict(), as the class's
            # documented PRODUCTION deserialization entry point -- but
            # unlike from_dict() (which restores the full, possibly
            # multi-hop _deprecated_index directly, see that method's
            # own comment), this one only ever reconstructed the ONE-HOP
            # deprecations each surviving record's own
            # deprecated_canonical_ids still lists, via register() above.
            # Same transitive-merge-chain loss to_dict()'s own comment
            # describes: in a chain A->B->C, C's record lists only [B],
            # never [A], so a merge history spanning more than one hop
            # was silently dropped every time a real caller used this
            # entry point instead of from_dict() -- found while auditing
            # this class's other real deserialization path, not because
            # any current production caller had actually hit it yet (no
            # real mapping-file producer emits a multi-hop chain today).
            # `.update()`, not assignment, so multiple files' own
            # deprecated_index entries accumulate additively across the
            # whole directory rather than the last file processed
            # silently discarding the others'.
            registry._deprecated_index.update(data.get("deprecated_index", {}))
        return registry

    # ------------------------------------------------------------------
    # Diagnostics
    # ------------------------------------------------------------------

    def stats(self) -> dict:
        """Return registry statistics for monitoring and health checks."""
        by_type: dict[str, int] = {}
        by_source: dict[str, int] = {}
        for record in self._canonical_index.values():
            t = record.canonical_id.entity_type.value
            by_type[t] = by_type.get(t, 0) + 1
            for source in record.cross_references:
                by_source[source.value] = by_source.get(source.value, 0) + 1
        return {
            "total_entities": len(self._canonical_index),
            "total_xrefs": len(self._xref_index),
            "deprecated_mappings": len(self._deprecated_index),
            "by_entity_type": by_type,
            "by_ontology_source": by_source,
        }
