# SAVE TO: packages\context-registry\src\biosim_context_registry\types.py
"""
Canonical identifier types for non-molecular biological CONTEXT entities.

Deliberately separate from `biosim_id_resolver.IDRegistry` (project decision, 2026-
08-18): IDRegistry is specifically an identifier registry for MOLECULAR
entities (proteins/metabolites/pathways/genes/reactions/compartments),
each cross-referenced against a real external ontology (UniProt, ChEBI,
HMDB, Reactome, KEGG) that this platform independently verifies live
before registering anything. Disease states and physiological axes are
semantically different -- there is no equivalent external molecular
ontology this platform integrates with for them -- and folding them into
IDRegistry would muddle a registry whose whole design (cross-reference
collision detection, deprecation/merge chains) exists specifically to
manage competing claims about a single MOLECULE's identity across
databases. A context entity's identity is just its own name; there is
nothing external to cross-reference or reconcile.

Deliberately MUCH SIMPLER than IDRegistry, not a mirror of it -- no
auto-incrementing numeric IDs (a context entity's name IS its identity,
there's no risk of two different diseases coincidentally sharing a
name the way two different molecules might share an external accession),
no cross-reference index, no deprecation/merge machinery. Build that
complexity later if a real need for it ever materializes (e.g. this
platform starts cross-referencing ICD-10/SNOMED-CT/HPO codes) -- not
speculatively now, matching this project's own YAGNI discipline.

Canonical ID format:
    CONTEXT:<entity_type>:<normalized_name>
    e.g. CONTEXT:DISEASE_STATE:chronic_kidney_disease
         CONTEXT:PHYSIOLOGICAL_AXIS:sleep_restriction
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from enum import Enum

# Cheap defense-in-depth cap, not a load-bearing limit today: every current
# caller (context_seed.py) passes short hardcoded literals. Worth having
# before register() is ever reachable from external/high-volume input.
MAX_CONTEXT_NAME_LENGTH = 200

_VALID_NORMALIZED_NAME = re.compile(r"^[a-z0-9_]+$")


class ContextEntityType(str, Enum):
    """
    Only the two kinds this platform currently has a real need for.
    NOT pre-building INTERVENTION/PHENOTYPE speculatively -- add a new
    member here when a real chunk actually needs one, following this
    project's own established "don't build for hypothetical future
    requirements" discipline.
    """
    DISEASE_STATE = "DISEASE_STATE"
    PHYSIOLOGICAL_AXIS = "PHYSIOLOGICAL_AXIS"


def normalize_name(name: str) -> str:
    """
    Lowercase, strip, and collapse whitespace runs to a single underscore --
    the canonical name form.

    Raises:
        ValueError: if the result is empty, exceeds MAX_CONTEXT_NAME_LENGTH,
            or contains any character outside [a-z0-9_]. In particular ':'
            is rejected: ContextID's own string format
            (CONTEXT:<entity_type>:<name>) uses ':' as its delimiter, so a
            name containing one would produce an ID that
            ContextID.from_string() could never parse back.
    """
    normalized = re.sub(r"\s+", "_", name.strip().lower())
    if not normalized:
        raise ValueError("Context entity name cannot be empty")
    if len(normalized) > MAX_CONTEXT_NAME_LENGTH:
        raise ValueError(
            f"Context entity name exceeds {MAX_CONTEXT_NAME_LENGTH} characters: '{name}'"
        )
    if not _VALID_NORMALIZED_NAME.match(normalized):
        raise ValueError(
            f"Context entity name '{name}' normalizes to '{normalized}', which "
            f"contains characters outside [a-z0-9_] -- only letters, digits, "
            f"and underscores are allowed (':' in particular would break "
            f"ContextID.from_string()'s round-trip)."
        )
    return normalized


@dataclass(frozen=True)
class ContextID:
    """
    The platform's internal canonical identifier for a non-molecular
    context entity. Unlike IDRegistry's CanonicalID (an opaque
    auto-incrementing number), a ContextID is directly derived from the
    entity's own normalized name -- deterministic, human-readable, and
    requiring no registry state to construct or predict.
    """
    entity_type: ContextEntityType
    name: str   # already normalized (see normalize_name)

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("ContextID name cannot be empty")
        normalized = normalize_name(self.name)
        if self.name != normalized:
            raise ValueError(
                f"ContextID name must already be normalized -- got "
                f"'{self.name}', expected '{normalized}'. Use "
                f"ContextRegistry.register()/resolve_by_name(), which "
                f"normalize automatically, rather than constructing a "
                f"ContextID directly with a raw name."
            )

    def __str__(self) -> str:
        return f"CONTEXT:{self.entity_type.value}:{self.name}"

    @classmethod
    def from_string(cls, s: str) -> "ContextID":
        parts = s.strip().split(":")
        if len(parts) != 3 or parts[0] != "CONTEXT":
            raise ValueError(f"Invalid context ID format: '{s}'. Expected 'CONTEXT:<TYPE>:<name>'")
        try:
            entity_type = ContextEntityType(parts[1])
        except ValueError:
            raise ValueError(f"Unknown context entity type '{parts[1]}' in context ID '{s}'")
        return cls(entity_type=entity_type, name=parts[2])


@dataclass(frozen=True)
class ContextRecord:
    """The registered record for a context entity: its ID, display name, and description."""
    context_id: ContextID
    display_name: str
    description: str = ""
