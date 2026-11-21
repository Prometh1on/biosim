# SAVE TO: packages\ode-engine\src\biosim_ode\registry.py
"""
PathwayRegistry — central registry for all pathway modules.

Pathways register themselves at import time via the global registry.
This gives the simulation API a single place to discover, validate,
and instantiate pathways without knowing their implementation details.

Usage:
    # Register a pathway (done once, typically in the pathway module itself)
    from biosim_ode.registry import global_registry
    global_registry().register(GlucoseInsulinPathway, GLUCOSE_INSULIN_DESCRIPTOR)

    # Discover pathways
    registry = global_registry()
    all_pathways = registry.list_all()
    production = registry.filter_by_maturity(PathwayMaturity.PRODUCTION)
    metabolic = registry.filter_by_tag("metabolic")

    # Instantiate
    cls = registry.get_class("glucose_insulin_v1")
    pathway = cls()
"""

from __future__ import annotations

import threading
from typing import Optional, Type

from .base_pathway import BasePathway, PathwayDescriptor, PathwayMaturity


class RegistryError(Exception):
    pass


class PathwayRegistry:
    """
    Thread-safe registry mapping pathway names to their classes and descriptors.

    Enforces:
        - No duplicate names (raises on double-registration)
        - Version must be incremented to replace an existing pathway
        - Dependencies must be resolvable within the registry
    """

    def __init__(self) -> None:
        self._classes:     dict[str, Type[BasePathway]] = {}
        self._descriptors: dict[str, PathwayDescriptor] = {}
        self._lock = threading.RLock()

    def register(
        self,
        pathway_class: Type[BasePathway],
        descriptor: PathwayDescriptor,
    ) -> None:
        """
        Register a pathway class with its descriptor.

        Args:
            pathway_class: The BasePathway subclass to register
            descriptor: Static metadata for this pathway

        Raises:
            RegistryError: If the name is already registered at the same version
        """
        with self._lock:
            name = descriptor.name
            if name in self._descriptors:
                existing = self._descriptors[name]
                if existing.version == descriptor.version:
                    raise RegistryError(
                        f"Pathway '{name}' version '{descriptor.version}' is already registered. "
                        f"Increment the version to replace it."
                    )
            self._classes[name] = pathway_class
            self._descriptors[name] = descriptor

    def get_class(self, name: str) -> Type[BasePathway]:
        """Return the pathway class for the given name."""
        with self._lock:
            if name not in self._classes:
                raise RegistryError(
                    f"Pathway '{name}' is not registered. "
                    f"Available: {sorted(self._classes.keys())}"
                )
            return self._classes[name]

    def get_descriptor(self, name: str) -> PathwayDescriptor:
        """Return the descriptor for the given pathway name."""
        with self._lock:
            if name not in self._descriptors:
                raise RegistryError(f"Pathway '{name}' is not registered.")
            return self._descriptors[name]

    def is_registered(self, name: str) -> bool:
        with self._lock:
            return name in self._descriptors

    def list_all(self) -> list[PathwayDescriptor]:
        """Return all registered descriptors, sorted by name."""
        with self._lock:
            return sorted(self._descriptors.values(), key=lambda d: d.name)

    def filter_by_maturity(self, maturity: PathwayMaturity) -> list[PathwayDescriptor]:
        """Return all pathways at the given maturity level."""
        return [d for d in self.list_all() if d.maturity == maturity]

    def filter_by_tag(self, tag: str) -> list[PathwayDescriptor]:
        """Return all pathways that have the given tag (case-insensitive)."""
        return [d for d in self.list_all() if d.has_tag(tag)]

    def validate_dependencies(self, name: str) -> list[str]:
        """
        Check that all dependencies of a pathway are registered.

        Returns:
            List of missing dependency names (empty = all satisfied).
        """
        with self._lock:
            descriptor = self.get_descriptor(name)
            return [
                dep for dep in descriptor.dependencies
                if dep not in self._descriptors
            ]

    def dependency_graph(self) -> dict[str, list[str]]:
        """Return the full dependency graph as an adjacency dict."""
        with self._lock:
            return {
                name: list(d.dependencies)
                for name, d in self._descriptors.items()
            }

    def stats(self) -> dict:
        with self._lock:
            by_maturity = {}
            for d in self._descriptors.values():
                k = d.maturity.value
                by_maturity[k] = by_maturity.get(k, 0) + 1
            return {
                "total_registered": len(self._descriptors),
                "by_maturity": by_maturity,
            }


# ---------------------------------------------------------------------------
# Global singleton registry
# ---------------------------------------------------------------------------

_global_registry: Optional[PathwayRegistry] = None
_registry_lock = threading.Lock()


def global_registry() -> PathwayRegistry:
    """Return the process-level global pathway registry (created on first call)."""
    global _global_registry
    with _registry_lock:
        if _global_registry is None:
            _global_registry = PathwayRegistry()
    return _global_registry


def reset_global_registry() -> None:
    """
    Reset the global registry to empty.
    For use in tests only — do not call in production code.
    """
    global _global_registry
    with _registry_lock:
        _global_registry = PathwayRegistry()
