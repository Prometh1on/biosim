# SAVE TO: packages\ode-engine\src\biosim_ode\base_pathway.py
"""
BasePathway abstract base class and pathway metadata types.

Every pathway in the platform must implement BasePathway.
This contract is what allows the CompositePathway to combine
multiple pathways without knowing their internals, and what
allows the PathwayRegistry to manage them uniformly.

Pathway maturity tiers:
    PRODUCTION   Validated against clinical data (Sprint 12).
                 Suitable for research publication. Default CI is narrow.
    BETA         Literature-parameterized, not yet clinically validated.
                 Suitable for exploration. CI is wider.
    EXPERIMENTAL Sparse data, many parameters estimated.
                 For development only. Wide CI, explicit uncertainty warnings.

Output signals:
    Pathways export named scalar values at each time step via output_signals().
    These are the values other pathways can read via CouplingSpec.
    Example: GlucoseInsulinPathway exports "insulin_pM", "glucose_mM", "mTOR_activity".
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import Enum
from typing import Optional


class PathwayMaturity(str, Enum):
    PRODUCTION   = "production"    # Validated against clinical data
    BETA         = "beta"          # Parameterized, not fully validated
    EXPERIMENTAL = "experimental"  # Sparse data, development only


@dataclass
class PathwayDescriptor:
    """
    Static metadata about a pathway module.
    Stored in the PathwayRegistry independently of pathway instances.
    """
    name: str                          # Unique identifier, e.g. "glucose_insulin_v1"
    version: str                       # Semantic version, e.g. "1.0.0"
    maturity: PathwayMaturity
    description: str
    species_names: list[str]           # All state variable names
    output_signals: list[str]          # Named values this pathway exports
    input_signals: list[str]           # Named values this pathway needs from others
    dependencies: list[str]            # Names of other pathways this depends on
    tags: list[str]                    # e.g. ["metabolic", "glucose", "insulin"]
    author: str = "biosim-core"
    literature_refs: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("PathwayDescriptor name cannot be empty")
        if not self.version:
            raise ValueError("PathwayDescriptor version cannot be empty")

    @property
    def n_species(self) -> int:
        return len(self.species_names)

    def has_tag(self, tag: str) -> bool:
        return tag.lower() in [t.lower() for t in self.tags]

    def has_dependency(self, name: str) -> bool:
        return name in self.dependencies


class BasePathway(ABC):
    """
    Abstract base class for all pathway modules.

    Every pathway must implement:
        descriptor   — static metadata (name, version, maturity, signals)
        species_names — names of all state variables
        species_units — units of all state variables (internal ODE units)
        y0           — initial conditions at basal steady state
        dydt         — ODE right-hand side
        output_signals — values exported to other pathways at each time step

    The dydt signature accepts an optional signals dict, which carries
    values injected from other pathways via CouplingSpec. At the single-
    pathway level, signals is None and the pathway uses its own parameters.
    """

    @property
    @abstractmethod
    def descriptor(self) -> PathwayDescriptor:
        """Static pathway metadata."""
        ...

    @property
    def name(self) -> str:
        return self.descriptor.name

    @property
    def maturity(self) -> PathwayMaturity:
        return self.descriptor.maturity

    @property
    def is_production(self) -> bool:
        return self.maturity == PathwayMaturity.PRODUCTION

    @property
    @abstractmethod
    def species_names(self) -> list[str]:
        """Names of all state variables in the same order as y0 and dydt."""
        ...

    @property
    @abstractmethod
    def species_units(self) -> list[str]:
        """Internal ODE units for each species (e.g. 'mM', 'pM', 'min-1')."""
        ...

    @property
    @abstractmethod
    def y0(self) -> list[float]:
        """
        Initial conditions at analytic basal steady state.
        All derivatives must be ~0 at this state with no external input.
        """
        ...

    @abstractmethod
    def dydt(
        self,
        t: float,
        y: list[float],
        signals: Optional[dict[str, float]] = None,
    ) -> list[float]:
        """
        ODE right-hand side.

        Args:
            t: Current time (internal units — minutes for this engine)
            y: Current state vector
            signals: Optional dict of values injected from other pathways.
                     None when running in standalone mode.

        Returns:
            List of derivatives in the same order as y and species_names.
        """
        ...

    @abstractmethod
    def output_signals(self, t: float, y: list[float]) -> dict[str, float]:
        """
        Compute the named scalar values this pathway exports at time t.
        These values are available to other pathways via CouplingSpec.

        Returns:
            Dict mapping signal name to current value.
        """
        ...

    def validate_y0(self, tolerance: float = 1e-4) -> list[str]:
        """
        Check that y0 is a valid steady state (all derivatives near zero).
        Returns a list of species names whose derivatives exceed the tolerance.
        """
        y0 = self.y0
        derivs = self.dydt(0.0, y0, signals=None)
        violations = []
        for i, (name, d) in enumerate(zip(self.species_names, derivs)):
            if abs(d) > tolerance:
                violations.append(f"{name}: d/dt = {d:.6g} (tolerance {tolerance})")
        return violations
