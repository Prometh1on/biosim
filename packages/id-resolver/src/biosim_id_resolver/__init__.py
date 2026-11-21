"""
biosim_id_resolver — canonical biological entity ID resolution.

Public API (the only surface other modules may import from):

    from biosim_id_resolver import (
        CanonicalID,
        EntityType,
        ExternalID,
        IDRegistry,
        OntologySource,
        ResolutionError,
        AmbiguousResolutionError,
    )

Do not import from submodules directly. The internal module structure
is subject to change; only this public surface is stable.
"""

from .registry import AmbiguousResolutionError, IDRegistry, ResolutionError
from .types import (
    CanonicalID,
    EntityType,
    ExternalID,
    OntologySource,
    ResolutionRecord,
)

__all__ = [
    "CanonicalID",
    "EntityType",
    "ExternalID",
    "IDRegistry",
    "OntologySource",
    "ResolutionError",
    "AmbiguousResolutionError",
    "ResolutionRecord",
]

__version__ = "0.1.0"
