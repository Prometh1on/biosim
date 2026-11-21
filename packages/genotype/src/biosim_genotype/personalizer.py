# SAVE TO: packages\genotype\src\biosim_genotype\personalizer.py
"""GenotypePersonalizer -- applies genetic variant effects to PathwayParameters."""
from __future__ import annotations
import dataclasses
import math
from dataclasses import dataclass, field
from typing import Optional
from biosim_ode import PathwayParameters
from .composer import ComposedEffect, EffectComposer
from .effect_db import Population, VariantEffectDB
from .variant import Variant

@dataclass
class PersonalizationResult:
    """Result of genetic personalisation. Contains no raw genetic data."""
    original_params:      PathwayParameters
    personalized_params:  PathwayParameters
    n_variants_input:     int
    n_variants_matched:   int
    n_variants_applied:   int
    n_parameters_modified: int
    genetic_uncertainty_cv: float
    population_warnings:  list[str]
    applied_effects:      list[ComposedEffect]

    @property
    def any_warnings(self) -> bool:
        return len(self.population_warnings) > 0

    @property
    def is_meaningfully_personalised(self) -> bool:
        return any(abs(e.percent_change) > 1.0 for e in self.applied_effects)

    @property
    def largest_effect(self):
        if not self.applied_effects:
            return None
        return max(self.applied_effects, key=lambda e: abs(e.fold_change - 1.0))

    def summary(self) -> dict:
        return {
            "n_variants_input":     self.n_variants_input,
            "n_variants_matched":   self.n_variants_matched,
            "n_variants_applied":   self.n_variants_applied,
            "n_parameters_modified": self.n_parameters_modified,
            "genetic_uncertainty_cv": round(self.genetic_uncertainty_cv, 4),
            "is_meaningfully_personalised": self.is_meaningfully_personalised,
            "any_population_warnings": self.any_warnings,
            "modified_parameters": {
                e.parameter_field: round(e.percent_change, 2)
                for e in self.applied_effects if e.is_modified
            },
        }

class GenotypePersonalizer:
    """Applies genetic variant effects to a PathwayParameters instance."""
    def __init__(self, effect_db=None, composer=None):
        self._effect_db = effect_db or VariantEffectDB.default()
        self._composer  = composer or EffectComposer()

    def personalize(self, variants, base_params=None, user_population=Population.ALL):
        base_params = base_params or PathwayParameters()
        n_input = len(variants)
        known_rsids = set(self._effect_db.all_rsids())
        relevant = [v for v in variants if v.rsid in known_rsids]
        n_matched = len(relevant)
        if not relevant:
            return PersonalizationResult(
                original_params=base_params, personalized_params=base_params,
                n_variants_input=n_input, n_variants_matched=0,
                n_variants_applied=0, n_parameters_modified=0,
                genetic_uncertainty_cv=0.0, population_warnings=[],
                applied_effects=[],
            )
        base_values = dataclasses.asdict(base_params)
        composed_effects = self._composer.compose(relevant, base_values, self._effect_db, user_population)
        # Counts distinct fingerprints, not raw rsids -- contributing_variants
        # holds privacy-safe fingerprints (see composer.py's ComposedEffect).
        applied_variant_fingerprints = {
            fp for e in composed_effects for fp, _ in e.contributing_variants
        }
        overrides = {e.parameter_field: e.final_value for e in composed_effects if e.is_modified}
        personalized = dataclasses.replace(base_params, **overrides)
        all_warnings = [w for e in composed_effects for w in e.population_warnings]
        combined_cv = math.sqrt(sum(e.genetic_cv ** 2 for e in composed_effects))
        return PersonalizationResult(
            original_params=base_params, personalized_params=personalized,
            n_variants_input=n_input, n_variants_matched=n_matched,
            n_variants_applied=len(applied_variant_fingerprints),
            n_parameters_modified=len([e for e in composed_effects if e.is_modified]),
            genetic_uncertainty_cv=combined_cv, population_warnings=all_warnings,
            applied_effects=composed_effects,
        )
