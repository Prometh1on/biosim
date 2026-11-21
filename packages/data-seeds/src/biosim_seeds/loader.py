# SAVE TO: packages\data-seeds\src\biosim_seeds\loader.py
"""
Seed loader — creates pre-seeded IDRegistry and ParameterStore instances
and connects them to the ODE pathway model.

This is the entry point for the real-data phase. Instead of using
hardcoded parameter defaults, the simulation can now pull values from
the literature-sourced ParameterStore.

ParameterLoader.load_pathway_params() returns a PathwayParameters instance
populated with the best available literature values. It falls back
gracefully to hardcoded defaults for any parameter not yet in the store.
"""

from __future__ import annotations

import dataclasses
from typing import Optional

from biosim_id_resolver import IDRegistry
from biosim_parameter_db import ConfidenceTier, Organism, ParameterKind, ParameterStore
from biosim_ode import PathwayParameters

from .parameter_seed import create_seeded_store
from .registry_seed import create_seeded_registry


# Mapping from PathwayParameters field name → (entity_name, ParameterKind,
# object_entity_name). Used by ParameterLoader to look up each ODE parameter
# in the ParameterStore. entity_name/object_entity_name must match names
# registered in the seeded IDRegistry. object_entity_name is None for
# subject-only parameters (Km, abundance, etc.); it must be set for any
# parameter seeded with ParameterBuilder.for_object() (e.g. EC50_mTOR is
# stored as mTOR's EC50 *for insulin receptor activation specifically* --
# see ParameterStore.query()'s docstring for why object_canonical_id is an
# exact-match filter, not a wildcard: omitting it here would silently miss
# the seeded value and fall back to the hardcoded default instead).
PARAMETER_MAPPINGS: list[tuple[str, str, ParameterKind, Optional[str]]] = [
    # (pathway_field,       entity_name,         parameter_kind,          object_entity_name)
    ("Gb",                 "D-Glucose",          ParameterKind.ABUNDANCE, None),
    ("Ib",                 "insulin",            ParameterKind.ABUNDANCE, None),
    ("ke",                 "insulin",            ParameterKind.RATE,      None),
    ("p2",                 "insulin receptor",   ParameterKind.RATE,      None),
    ("p3",                 "insulin receptor",   ParameterKind.CUSTOM,    None),
    ("EC50_mTOR",          "mTOR",               ParameterKind.EC50,      "insulin receptor"),
    ("IGF1b",              "IGF-1",              ParameterKind.ABUNDANCE, None),
    # MASTER_PLAN 5.1 chunk 3: Rebelato et al. 2018, PMID:30166558.
    # Confirmed empirically to make zero difference to any simulation
    # output -- PathwayParameters()'s own hardcoded default (5.63) already
    # equals this literature-derived value (applied directly in a prior
    # session); this mapping only makes it also resolvable via
    # ParameterStore/ParameterLoader, matching this campaign's own goal.
    # Mandatory code-reviewer finding (LOW): this entry is registered
    # under Organism.OTHER (rodent islet data -- see its own seed
    # comment), NOT Organism.HUMAN, so load_pathway_params()'s own
    # organism=Organism.HUMAN default query below does NOT find it --
    # it's retrieved only via this method's own "retry without organism
    # filter" fallback a few lines down. Confirmed this fallback path is
    # what actually fires here (not a coincidental match) before relying
    # on it; disclosed explicitly so a future reader doesn't assume a
    # rodent-tagged entry would surface through the human-first lookup
    # alone.
    ("KG_secretion",       "insulin",            ParameterKind.EC50,      None),
    # MASTER_PLAN 5.1 chunk 4: p1 (glucose effectiveness) was already
    # seeded (Bergman 1981 GOLD under D-Glucose/RATE; a wrongly-cited
    # DeFronzo 1979 replicate was removed in B4, 2026-11-16) but never
    # mapped. store.best() returns Bergman's 0.0317, identical to
    # PathwayParameters()'s default; the entry is .in_human(), matching
    # load_pathway_params()'s default query organism.
    ("p1",                 "D-Glucose",          ParameterKind.RATE,      None),
]

# Parameters not yet mapped because they require composite computation
# or have no direct literature analogue (k_mTOR, k_AMPK etc.)
# These will be added in future sprints as literature coverage expands.
UNMAPPED_PARAMETERS = [
    "k_mTOR",       # mTOR relaxation rate — no direct measurement available
    "k_AMPK",       # AMPK relaxation rate — no direct measurement available
    "EC50_AMPK",    # mapped as IC50 on AMPK
    "G_threshold",  # beta-cell glucose threshold — not yet mapped
    # Regression: "p3" and "KG_secretion" were both listed here despite
    # already being (or, for KG_secretion, now being) real entries in
    # PARAMETER_MAPPINGS above -- p3's own stale entry predates this
    # chunk (found in passing while adding KG_secretion's real mapping
    # to the same list); removed both rather than leaving a second,
    # contradictory "not yet mapped" claim sitting next to the mapping
    # that already handles it.
]


class ParameterLoader:
    """
    Loads ODE pathway parameters from a seeded ParameterStore.

    Falls back to PathwayParameters defaults for any parameter not found
    in the store. This ensures the model always has valid parameters
    even when the literature database is partially populated.

    Usage:
        loader = ParameterLoader()
        params = loader.load_pathway_params()
        pathway = GlucoseInsulinPathway(params=params)
    """

    def __init__(
        self,
        store:    Optional[ParameterStore] = None,
        registry: Optional[IDRegistry]     = None,
    ) -> None:
        self._registry = registry or create_seeded_registry()
        self._store    = store    or create_seeded_store(self._registry)

    def load_pathway_params(
        self,
        organism:       Organism       = Organism.HUMAN,
        min_confidence: ConfidenceTier = ConfidenceTier.BRONZE,
    ) -> PathwayParameters:
        """
        Build a PathwayParameters instance from literature values.

        Args:
            organism:       Prefer parameters measured in this organism
            min_confidence: Minimum confidence tier to accept

        Returns:
            PathwayParameters with literature values where available,
            falling back to hardcoded defaults elsewhere.
        """
        defaults   = PathwayParameters()
        overrides: dict[str, float] = {}
        provenance: dict[str, str]  = {}

        for field_name, entity_name, kind, object_entity_name in PARAMETER_MAPPINGS:
            entity_ids = self._registry.search_by_name(entity_name)
            if not entity_ids:
                continue
            entity_id = str(entity_ids[0])

            object_id: Optional[str] = None
            if object_entity_name is not None:
                object_ids = self._registry.search_by_name(object_entity_name)
                if not object_ids:
                    continue
                object_id = str(object_ids[0])

            param = self._store.best(
                subject_canonical_id=entity_id,
                kind=kind,
                min_confidence=min_confidence,
                organism=organism,
                object_canonical_id=object_id,
            )
            if param is None:
                # Try without organism filter
                param = self._store.best(
                    subject_canonical_id=entity_id,
                    kind=kind,
                    min_confidence=min_confidence,
                    object_canonical_id=object_id,
                )

            if param is not None:
                base_val = getattr(defaults, field_name, None)
                if base_val is not None and isinstance(base_val, (int, float)):
                    overrides[field_name]  = param.value.value
                    provenance[field_name] = param.provenance.pubmed_id or "unknown"

        loaded = dataclasses.replace(defaults, **overrides)
        return loaded

    def provenance_report(self) -> dict[str, str]:
        """
        Return a report of which parameters came from the store vs defaults.
        Maps parameter field name → PubMed ID (or 'default' if not loaded).
        """
        defaults   = PathwayParameters()
        loaded     = self.load_pathway_params()
        report: dict[str, str] = {}

        for f in dataclasses.fields(defaults):
            default_val = getattr(defaults, f.name)
            loaded_val  = getattr(loaded,   f.name)
            if isinstance(default_val, float):
                if abs(loaded_val - default_val) > 1e-12:
                    report[f.name] = "from_literature"
                else:
                    report[f.name] = "default"
        return report

    @property
    def store(self) -> ParameterStore:
        return self._store

    @property
    def registry(self) -> IDRegistry:
        return self._registry

    def store_stats(self) -> dict:
        return {
            "registry": self._registry.stats(),
            "store":    self._store.stats(),
        }
