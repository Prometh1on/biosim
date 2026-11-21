# SAVE TO: packages\context-registry\src\biosim_context_registry\__init__.py
"""
biosim_context_registry -- registry for non-molecular biological CONTEXT
entities (disease states, physiological/behavioral axes), deliberately
separate from biosim_id_resolver.IDRegistry (see types.py's own module
docstring for the full architectural reasoning, project decision
2026-08-18).

Public API:
    from biosim_context_registry import (
        ContextRegistry, ContextID, ContextRecord, ContextEntityType,
        ContextResolutionError,
    )
"""
from .registry import ContextRegistry, ContextResolutionError
from .types import ContextEntityType, ContextID, ContextRecord, normalize_name

__all__ = [
    "ContextEntityType",
    "ContextID",
    "ContextRecord",
    "ContextRegistry",
    "ContextResolutionError",
    "normalize_name",
]
__version__ = "0.1.0"
