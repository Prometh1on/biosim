# SAVE TO: packages\ode-engine\src\biosim_ode\coupling.py
"""
Inter-pathway coupling — the architecture for multi-pathway models.

CouplingSpec declares a directed dependency between two pathways:
    source_pathway exports a named signal
    target_pathway reads that signal as an input modifier

CompositePathway combines multiple pathways into one ODE system by:
    1. Concatenating their state vectors
    2. At each time step, computing output signals from each pathway
    3. Passing those signals as inputs to dependent pathways
    4. Calling each pathway's dydt with its slice of the state vector

This is the extensibility mechanism that allows new pathways
(cortisol axis, lipid metabolism, gut microbiome metabolites) to
plug into the model without modifying existing pathway code.

For v1, CompositePathway supports up to 4 coupled pathways.
The architecture handles N pathways; the limit is practical, not structural.

Example (future Sprint 8):
    coupling = CouplingSpec(
        source_pathway="glucose_insulin_v1",
        signal_name="insulin_pM",
        target_pathway="cortisol_axis_v1",
        target_input="insulin_suppression_pM",
        transform=lambda x: x,  # identity — direct pass-through
    )
"""

from __future__ import annotations

import fnmatch
from dataclasses import dataclass, field
from typing import Callable, Optional
import numpy as np

from .base_pathway import BasePathway, PathwayDescriptor, PathwayMaturity

# EXTENSION-PATHWAY COUPLING CONTRACT (CLAUDE.md rule 11): every
# multiplicative modifier signal GlucoseInsulinPathway.dydt() reads has a
# neutral default of 1.0; every additive delta/effect signal has a neutral
# default of 0.0. When more than one CouplingSpec targets the SAME
# (target_pathway, target_input), CompositePathway.dydt() must combine their
# contributions accordingly (product for modifiers, sum for deltas) rather
# than letting the later coupling in the list silently overwrite the
# earlier one's value -- see KNOWN BUGS in CLAUDE.md (cortisol_couplings()
# and exercise_couplings() both target gb_modifier/p3_modifier on GI).
MULTIPLICATIVE_SIGNAL_KEYS = frozenset({
    "p3_modifier",
    "gb_modifier",
    "insulin_resistance_modifier",
    "gastric_emptying_modifier",
    "p1_modifier",
    "extra_g_clearance_modifier",
    "secretion_scale_modifier",
})
ADDITIVE_SIGNAL_KEYS = frozenset({
    "hepatic_egp_delta",
    "exercise_glucose_uptake_delta",
    "incretin_effect",
    "renal_glucose_clearance_delta",
    "exogenous_insulin_delta",
})


@dataclass
class CouplingSpec:
    """
    Declares a directed signal dependency between two pathways.

    source_pathway  → exports signal_name
    target_pathway  ← receives it as target_input
    transform       — optional function applied to the signal before injection
                      (default: identity). Use for unit conversion or scaling.
    """
    source_pathway: str
    signal_name: str
    target_pathway: str
    target_input: str
    transform: Callable[[float], float] = field(default=lambda x: x)

    def validate(
        self,
        available_pathways: list[str],
        source_signals: list[str],
        target_inputs: list[str],
    ) -> list[str]:
        """
        Validate this coupling against available pathway metadata.

        Returns:
            List of validation error strings (empty = valid)
        """
        errors = []
        if self.source_pathway not in available_pathways:
            errors.append(
                f"Source pathway '{self.source_pathway}' is not registered"
            )
        if self.target_pathway not in available_pathways:
            errors.append(
                f"Target pathway '{self.target_pathway}' is not registered"
            )
        if self.source_pathway == self.target_pathway:
            errors.append(
                f"A pathway cannot couple to itself: '{self.source_pathway}'"
            )
        if self.signal_name not in source_signals:
            errors.append(
                f"Signal '{self.signal_name}' is not exported by '{self.source_pathway}'"
            )
        # Regression: this method accepted target_inputs as a parameter
        # but never actually checked self.target_input against it -- a
        # typo'd target_input in any pathway-extension wiring would
        # never be caught here.
        #
        # Wildcard support (MASTER_PLAN 3.2, needed to wire this into
        # production): a target_inputs entry containing "*" is matched
        # via fnmatch (shell-glob semantics: "*" matches any substring,
        # anywhere in the pattern, not just a trailing prefix) -- needed
        # for GlucoseInsulinPathway's "meal*_remaining_mM" entry, which
        # represents the meal-absorption pathway's genuinely dynamic,
        # one-per-meal signal names (meal0_remaining_mM, meal1_remaining_mM,
        # ...). A pathway's static PathwayDescriptor can't enumerate an
        # unbounded, scenario-dependent key set any other way -- exact
        # membership would reject this real, already-shipped coupling.
        # A plain equality check first (not just relying on fnmatch on
        # every entry) keeps ordinary, non-wildcard names an exact match,
        # not accidentally glob-interpreted if one ever contains a
        # literal "*"/"?"/"[" character.
        #
        # fnmatchcase(), not fnmatch() -- found by the mandatory code-
        # reviewer pass: plain fnmatch.fnmatch() normalizes both operands
        # via os.path.normcase() first, which lowercases on Windows but
        # is a no-op on POSIX -- so a case-mismatched typo (e.g.
        # "Meal0_remaining_mM") would silently PASS this validation on a
        # Windows dev machine while correctly failing on Linux CI/
        # production, exactly the kind of typo this milestone exists to
        # catch. fnmatchcase() is always case-sensitive regardless of OS.
        if not any(
            self.target_input == t or ("*" in t and fnmatch.fnmatchcase(self.target_input, t))
            for t in target_inputs
        ):
            errors.append(
                f"Target input '{self.target_input}' is not accepted by '{self.target_pathway}'"
            )
        return errors


class CompositePathway(BasePathway):
    """
    Combines multiple pathways into a single ODE system.

    The composite state vector is the concatenation of each
    individual pathway's state vector:
        y = [y_pathway_0..., y_pathway_1..., y_pathway_2...]

    Species names are prefixed with the pathway name to avoid collisions:
        "glucose_insulin_v1.G", "glucose_insulin_v1.I", ...

    Signal passing:
        At each time step, output_signals() is called for each pathway.
        These are collected and injected into dependent pathways via
        the CouplingSpec list before calling their dydt().
    """

    def __init__(
        self,
        pathways: list[BasePathway],
        couplings: Optional[list[CouplingSpec]] = None,
    ) -> None:
        if not pathways:
            raise ValueError("CompositePathway requires at least one pathway")
        self._pathways = list(pathways)
        self._couplings = couplings or []
        # Wired in for MASTER_PLAN 3.2 -- CouplingSpec.validate() has
        # existed since its own regression fix (checking target_input
        # against target_inputs) but was never actually CALLED anywhere
        # in production, so a typo'd signal/target_input name in any
        # pathway-extension wiring would silently compose into an inert,
        # zero-effect coupling with no error at all -- exactly the bug
        # class in CLAUDE.md's KNOWN BUGS #7 (CompositePathway coupling
        # was architecturally inert), which cost a full session to find
        # the first time. Fail-fast here at construction, matching this
        # project's established convention for catching a programmer
        # error at the earliest possible point (DrugEvent.__post_init__,
        # PhysiologicalProfile's disease-state exclusivity checks, etc.)
        # rather than letting a silently-inert coupling ship. Verified
        # safe against every coupling function in this package before
        # wiring this in (systematic sweep, not just the one instance
        # this pass happened to add) -- the one real gap found
        # (MealAbsorptionPathway's dynamic meal{i}_remaining_mM signal
        # names) is handled by CouplingSpec.validate()'s new trailing-"*"
        # wildcard support and GlucoseInsulinPathway's own
        # "meal*_remaining_mM" input_signals entry, not by weakening this
        # check.
        self._validate_couplings()
        self._slice_indices = self._compute_slices()
        self._descriptor = self._build_descriptor()

    def _validate_couplings(self) -> None:
        pathways_by_name = {p.name: p for p in self._pathways}
        available_pathways = list(pathways_by_name.keys())
        errors: list[str] = []
        for coupling in self._couplings:
            source = pathways_by_name.get(coupling.source_pathway)
            target = pathways_by_name.get(coupling.target_pathway)
            source_signals = source.descriptor.output_signals if source is not None else []
            target_inputs = target.descriptor.input_signals if target is not None else []
            errors.extend(coupling.validate(available_pathways, source_signals, target_inputs))
        if errors:
            raise ValueError(
                "CompositePathway coupling validation failed:\n"
                + "\n".join(f"  - {e}" for e in errors)
            )

    # ------------------------------------------------------------------
    # BasePathway interface
    # ------------------------------------------------------------------

    @property
    def descriptor(self) -> PathwayDescriptor:
        return self._descriptor

    @property
    def species_names(self) -> list[str]:
        names = []
        for p in self._pathways:
            names.extend([f"{p.name}.{s}" for s in p.species_names])
        return names

    @property
    def species_units(self) -> list[str]:
        units = []
        for p in self._pathways:
            units.extend(p.species_units)
        return units

    @property
    def y0(self) -> list[float]:
        y0 = []
        for p in self._pathways:
            y0.extend(p.y0)
        return y0

    def dydt(
        self,
        t: float,
        y: list[float],
        signals: Optional[dict[str, float]] = None,
    ) -> list[float]:
        """
        Compute derivatives for all pathways.

        At each call:
            1. Split composite y into per-pathway slices
            2. Compute each pathway's output signals
            3. Build signal dicts for each pathway from couplings
            4. Call each pathway's dydt with its signals
            5. Concatenate and return all derivatives
        """
        y_slices = self._split_y(y)
        injected = self._injected_signals(t, y_slices)

        # Call each pathway's dydt with its signals
        all_derivs = []
        for i, p in enumerate(self._pathways):
            pathway_signals = injected[p.name] if injected[p.name] else None
            derivs = p.dydt(t, y_slices[i], signals=pathway_signals)
            all_derivs.extend(derivs)

        return all_derivs

    def injected_signals(self, t: float, y: list[float]) -> dict[str, dict[str, float]]:
        """
        The coupled signals each sub-pathway's dydt() receives at (t, y),
        keyed by pathway name -- exactly what dydt() injects, for post-solve
        inspection (e.g. ODESystem's point-of-use clamp disclosure).
        """
        return self._injected_signals(t, self._split_y(y))

    def _injected_signals(
        self, t: float, y_slices: list[list[float]],
    ) -> dict[str, dict[str, float]]:
        # Compute output signals for all pathways
        all_outputs: dict[str, dict[str, float]] = {}
        for i, p in enumerate(self._pathways):
            all_outputs[p.name] = p.output_signals(t, y_slices[i])

        # Build signal injection dicts from coupling specs. Contributions
        # are collected per (target_pathway, target_input) rather than
        # assigned directly, so a target_input written by more than one
        # coupling accumulates instead of the later coupling in the list
        # silently discarding the earlier one's value -- see
        # MULTIPLICATIVE_SIGNAL_KEYS/ADDITIVE_SIGNAL_KEYS above.
        contributions: dict[str, dict[str, list[float]]] = {p.name: {} for p in self._pathways}
        for coupling in self._couplings:
            source_outputs = all_outputs.get(coupling.source_pathway, {})
            if coupling.signal_name in source_outputs:
                raw_value = source_outputs[coupling.signal_name]
                transformed = coupling.transform(raw_value)
                contributions[coupling.target_pathway].setdefault(
                    coupling.target_input, []
                ).append(transformed)

        injected: dict[str, dict[str, float]] = {p.name: {} for p in self._pathways}
        for target_name, target_inputs in contributions.items():
            for target_input, values in target_inputs.items():
                injected[target_name][target_input] = self._combine_signal_values(
                    target_input, values
                )
        return injected

    @staticmethod
    def _combine_signal_values(target_input: str, values: list[float]) -> float:
        """
        Combine multiple couplings' contributions to the same target_input.

        A single contribution passes through unchanged (identical to the
        pre-fix behaviour for every currently-single-writer target_input,
        e.g. hepatic's "insulin_pM" or renal/GLP1's "glucose_mM"). Two or
        more contributions to a recognized rule-11 signal combine per its
        documented semantics: multiplicative modifiers multiply, additive
        deltas/effects sum. A collision on a target_input that is neither
        is genuinely ambiguous -- raise rather than silently guess.
        """
        if len(values) == 1:
            return values[0]
        if target_input in MULTIPLICATIVE_SIGNAL_KEYS:
            result = 1.0
            for v in values:
                result *= v
            return result
        if target_input in ADDITIVE_SIGNAL_KEYS:
            return sum(values)
        raise ValueError(
            f"Multiple couplings write to target_input '{target_input}', but "
            "it is not a recognized multiplicative/additive rule-11 signal "
            "(see MULTIPLICATIVE_SIGNAL_KEYS/ADDITIVE_SIGNAL_KEYS in "
            "coupling.py) -- cannot determine how to combine its values."
        )

    def output_signals(self, t: float, y: list[float]) -> dict[str, float]:
        """Return merged output signals from all pathways, namespaced by pathway."""
        y_slices = self._split_y(y)
        merged = {}
        for i, p in enumerate(self._pathways):
            for k, v in p.output_signals(t, y_slices[i]).items():
                merged[f"{p.name}.{k}"] = v
        return merged

    # ------------------------------------------------------------------
    # Composite-specific methods
    # ------------------------------------------------------------------

    @property
    def pathway_names(self) -> list[str]:
        return [p.name for p in self._pathways]

    @property
    def couplings(self) -> list[CouplingSpec]:
        return list(self._couplings)

    @property
    def pathways(self) -> tuple[BasePathway, ...]:
        return tuple(self._pathways)

    @property
    def accepts_meals(self) -> bool:
        """True if at least one sub-pathway can be driven by meal inputs."""
        return any(hasattr(p, "with_meals") for p in self._pathways)

    def with_meals(self, meals: list) -> "CompositePathway":
        """
        Return a new CompositePathway with meal inputs applied to whichever
        sub-pathway accepts them (currently only GlucoseInsulinPathway).

        This lets ODESystem drive a composite exactly like a bare
        GlucoseInsulinPathway: it can call with_meals(meals) without knowing
        whether the pathway is standalone or composed with extensions.
        """
        meal_capable = [p for p in self._pathways if hasattr(p, "with_meals")]
        if not meal_capable:
            raise ValueError(
                "with_meals() called on a CompositePathway containing no "
                "meal-capable pathway (none of "
                f"{[p.name for p in self._pathways]} implement with_meals())"
            )
        if len(meal_capable) > 1:
            raise ValueError(
                "Ambiguous with_meals(): more than one sub-pathway accepts "
                f"meals ({[p.name for p in meal_capable]})"
            )
        target_name = meal_capable[0].name
        new_pathways = [
            p.with_meals(meals) if p.name == target_name else p
            for p in self._pathways
        ]
        return CompositePathway(new_pathways, couplings=self._couplings)

    def get_pathway_slice(self, pathway_name: str) -> tuple[int, int]:
        """Return (start, end) indices for a pathway's species in the composite y."""
        for i, p in enumerate(self._pathways):
            if p.name == pathway_name:
                return self._slice_indices[i]
        raise KeyError(f"Pathway '{pathway_name}' not in composite")

    def extract_pathway_state(
        self,
        pathway_name: str,
        y: list[float],
    ) -> list[float]:
        """Extract the state slice for a specific pathway from the composite y."""
        start, end = self.get_pathway_slice(pathway_name)
        return y[start:end]

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _compute_slices(self) -> list[tuple[int, int]]:
        slices = []
        offset = 0
        for p in self._pathways:
            n = len(p.species_names)
            slices.append((offset, offset + n))
            offset += n
        return slices

    def _split_y(self, y: list[float]) -> list[list[float]]:
        return [list(y[s:e]) for s, e in self._slice_indices]

    def _build_descriptor(self) -> PathwayDescriptor:
        all_tags = list({tag for p in self._pathways for tag in p.descriptor.tags})
        all_deps = [p.name for p in self._pathways]
        worst_maturity = min(
            self._pathways,
            key=lambda p: ["production", "beta", "experimental"].index(p.maturity.value),
        ).maturity

        return PathwayDescriptor(
            name=f"composite({'_'.join(p.name for p in self._pathways)})",
            version="1.0.0",
            maturity=worst_maturity,
            description=f"Composite of: {', '.join(p.name for p in self._pathways)}",
            species_names=self.species_names,
            output_signals=[
                f"{p.name}.{s}"
                for p in self._pathways
                for s in p.descriptor.output_signals
            ],
            input_signals=[],
            dependencies=all_deps,
            tags=all_tags,
        )
