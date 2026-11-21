# SAVE TO: packages\genotype\src\biosim_genotype\__init__.py
"""
biosim_genotype -- genetic variant personalisation.

PRIVACY: All variant processing is in-memory. No raw genetic data
is written to disk or transmitted. PersonalizationResult contains
only derived parameter modifiers.

EPISTEMIC WARNING: Genetic effects explain a small fraction of
inter-individual variance. Results are for research exploration,
not clinical prediction.

Public API:
    from biosim_genotype import (
        Variant, Genotype,
        VariantEffectDB, EffectRecord, EffectType, Population,
        VariantParser, ParseResult,
        EffectComposer, ComposedEffect,
        GenotypePersonalizer, PersonalizationResult,
    )
"""
from .composer import ComposedEffect, EffectComposer
from .effect_db import EffectRecord, EffectType, Population, VariantEffectDB
from .parser import ParseResult, VariantParser
from .personalizer import GenotypePersonalizer, PersonalizationResult
from .variant import Genotype, Variant

__all__ = [
    "ComposedEffect",
    "EffectComposer",
    "EffectRecord",
    "EffectType",
    "GenotypePersonalizer",
    "Genotype",
    "ParseResult",
    "PersonalizationResult",
    "Population",
    "Variant",
    "VariantEffectDB",
    "VariantParser",
]
__version__ = "0.1.0"
