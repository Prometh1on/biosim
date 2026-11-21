# SAVE TO: packages\data-seeds\src\biosim_seeds\parameter_seed\_helpers.py
"""Name-to-ID resolution helpers and the seed ParameterBuilder."""

from __future__ import annotations

import dataclasses

from biosim_context_registry import ContextEntityType, ContextRegistry
from biosim_id_resolver import IDRegistry
from biosim_parameter_db import Parameter, ParameterBuilder

from .verification import lookup


class SeedParameterBuilder(ParameterBuilder):
    """ParameterBuilder that attaches the audited source-verification status
    (verification.VERIFICATION_RECORDS) to each seeded parameter. A value
    with no matching record stays UNCHECKED, which the seed tests reject."""

    def build(self) -> Parameter:
        parameter = super().build()
        label = parameter.custom_kind_name or parameter.kind.value
        record = lookup(parameter.provenance.pubmed_id, label, parameter.value.value)
        if record is None:
            return parameter
        provenance = dataclasses.replace(
            parameter.provenance,
            verification_status=record["status"],
            verification_quote=record["quote"],
        )
        return dataclasses.replace(parameter, provenance=provenance)


def _resolve(registry: IDRegistry, name: str) -> str:
    """Resolve a MOLECULAR entity name to its canonical ID string."""
    results = registry.search_by_name(name)
    if not results:
        raise ValueError(f"Entity '{name}' not found in registry. Run seed_registry() first.")
    return str(results[0])


def _resolve_context(
    context_registry: ContextRegistry, name: str, entity_type: ContextEntityType,
) -> str:
    """
    Resolve a non-molecular CONTEXT entity name to its canonical ID
    string (see biosim_context_registry -- deliberately separate from
    IDRegistry, so a separate resolve helper, not an overload of
    _resolve() above).
    """
    return str(context_registry.resolve_by_name(name, entity_type))
