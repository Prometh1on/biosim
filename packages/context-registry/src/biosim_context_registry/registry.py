# SAVE TO: packages\context-registry\src\biosim_context_registry\registry.py
"""
Registry for non-molecular biological CONTEXT entities (disease states,
physiological/behavioral axes) -- see types.py's own module docstring
for why this is deliberately separate from, and much simpler than,
biosim_id_resolver.IDRegistry.

Usage:
    registry = ContextRegistry()
    cid = registry.register("chronic_kidney_disease", ContextEntityType.DISEASE_STATE,
                             description="...")
    record = registry.get_record(cid)
    same_cid = registry.resolve_by_name("chronic_kidney_disease", ContextEntityType.DISEASE_STATE)
"""

from __future__ import annotations

from typing import Optional

from .types import ContextEntityType, ContextID, ContextRecord, normalize_name


class ContextResolutionError(Exception):
    """Raised when a context entity name cannot be resolved."""
    pass


class ContextRegistry:
    """
    A registry mapping (entity_type, normalized name) to a ContextID and
    its ContextRecord. No auto-incrementing IDs, no external cross-
    references, no deprecation/merge machinery -- see types.py's own
    docstring for why none of that is needed here (yet).
    """

    def __init__(self) -> None:
        self._records: dict[str, ContextRecord] = {}   # str(ContextID) -> record

    def register(
        self,
        name: str,
        entity_type: ContextEntityType,
        description: str = "",
    ) -> ContextID:
        """
        Register a new context entity and return its ContextID.

        Raises:
            ValueError: if an entity of this type and (normalized) name
                        is already registered -- this registry has no
                        merge/deprecation/rename machinery (unlike
                        IDRegistry, see types.py's own docstring), so a
                        duplicate registration fails loudly rather than
                        silently overwriting the existing record.
        """
        normalized = normalize_name(name)
        cid = ContextID(entity_type=entity_type, name=normalized)
        key = str(cid)
        if key in self._records:
            raise ValueError(
                f"Context entity already registered: {key}. "
                f"This registry has no rename/merge workflow -- if this is "
                f"a genuine duplicate, fix the caller; if it's an "
                f"intentional rename, construct a new ContextRegistry (or "
                f"edit context_seed.py directly) rather than mutating this "
                f"one in place."
            )
        self._records[key] = ContextRecord(
            context_id=cid, display_name=name, description=description,
        )
        return cid

    def resolve_by_name(
        self, name: str, entity_type: ContextEntityType,
    ) -> ContextID:
        """
        Resolve a context entity by its (case/whitespace-insensitive)
        name and entity type.

        Raises:
            ContextResolutionError: if no matching entity is registered.
        """
        cid = ContextID(entity_type=entity_type, name=normalize_name(name))
        if str(cid) not in self._records:
            raise ContextResolutionError(
                f"Cannot resolve context entity '{name}' ({entity_type.value}). "
                f"This entity is not registered -- add it via context_seed.py "
                f"if it's a real, needed subject."
            )
        return cid

    def get_record(self, context_id: ContextID) -> ContextRecord:
        """
        Retrieve the full record for a ContextID.

        Raises:
            ContextResolutionError: if the ContextID is not registered.
        """
        record = self._records.get(str(context_id))
        if record is None:
            raise ContextResolutionError(f"No record found for context ID {context_id}.")
        return record

    def is_registered(self, name: str, entity_type: ContextEntityType) -> bool:
        return str(ContextID(entity_type=entity_type, name=normalize_name(name))) in self._records

    def stats(self) -> dict:
        """Return registry statistics for monitoring and health checks."""
        by_type: dict[str, int] = {}
        for record in self._records.values():
            t = record.context_id.entity_type.value
            by_type[t] = by_type.get(t, 0) + 1
        return {
            "total_entities": len(self._records),
            "by_entity_type": by_type,
        }
