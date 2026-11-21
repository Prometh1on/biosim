# SAVE TO: packages\genotype\src\biosim_genotype\composer.py
"""
Effect composer — combines multiple variants into pathway parameter modifiers.

Composition model (v1):
    Parameters with MULTIPLICATIVE effects from multiple variants
    are combined multiplicatively:
        modifier = product_i (1 + effect_i * dosage_i)

    Parameters with ADDITIVE_ABSOLUTE effects are summed:
        offset = sum_i (effect_i * dosage_i)

    Mixed types: multiplicative and additive effects on the same
    parameter are applied in sequence (multiplicative first, then additive).

Uncertainty:
    Each applied variant inflates the parameter's uncertainty.
    The genetic uncertainty CV is combined in quadrature with the
    base parameter uncertainty (provided by the UQ layer in Sprint 5).

Population warnings:
    If a user's population differs from the effect record's population,
    a warning is attached to the ComposedEffect but the effect is still
    applied. Population-specific effect sizes are a Sprint 8 concern.
"""

from __future__ import annotations

import hashlib
import hmac
import secrets
from dataclasses import dataclass, field
from typing import Optional
import math

from .variant import Variant
from .effect_db import EffectRecord, EffectType, Population, VariantEffectDB


# Process-local, in-memory only -- generated fresh at import time, never
# read from an env var or config file, and never itself serialized. Same
# rationale and fix shape as biosim_api.scenario's identical fingerprint
# (that one can't be imported here without an unwanted cross-package
# dependency -- biosim_genotype sits below biosim_api -- so this is an
# independent, equivalent implementation).
_RSID_FINGERPRINT_KEY = secrets.token_bytes(32)


def _rsid_fingerprint(rsid: str) -> str:
    """
    A keyed fingerprint for one rsID, for ComposedEffect.contributing_variants.

    Rule 7 (CLAUDE.md): "PersonalizationResult must never contain raw
    rsIDs or allele data" -- PersonalizationResult.applied_effects
    carries ComposedEffect objects, and this project's own known-variant
    catalog (VariantEffectDB) has only ~14 real rsIDs, so an unkeyed hash
    would be an exhaustive lookup table, not a one-way function (see the
    identical finding already fixed for biosim_api.scenario
    ._variant_fingerprint()). Keying with a process-random secret closes
    that: without `_RSID_FINGERPRINT_KEY`, the same small space no
    longer resolves to anything. Trade-off, by design: this makes the
    fingerprint process-local, not a stable cross-process identifier --
    fine for counting distinct contributing variants within one
    personalize() call (the only real use today), not for cross-run
    identity matching.
    """
    return hmac.new(_RSID_FINGERPRINT_KEY, rsid.encode(), hashlib.sha256).hexdigest()[:16]


@dataclass
class ComposedEffect:
    """
    The net effect of all relevant variants on one pathway parameter.

    parameter_field:      PathwayParameters field name
    base_value:           value before genetic modification
    final_value:           value after all variants applied
    contributing_variants: privacy-safe rsID fingerprints and their dosages
                          -- NOT raw rsIDs (rule 7, CLAUDE.md; see
                          _rsid_fingerprint() above for why a bare rsID
                          isn't safe even hashed)
    genetic_cv:           additional CV due to genetic uncertainty
    population_warnings:  flags when effect population differs from user
    """
    parameter_field:       str
    base_value:            float
    final_value:           float
    contributing_variants: list[tuple[str, int]]  # (rsid_fingerprint, dosage)
    genetic_cv:            float
    population_warnings:   list[str] = field(default_factory=list)

    @property
    def fold_change(self) -> float:
        """Ratio final / base. 1.0 = no change."""
        if self.base_value == 0:
            return 1.0
        return self.final_value / self.base_value

    @property
    def percent_change(self) -> float:
        """Percent change from base: (final - base) / base * 100."""
        if self.base_value == 0:
            return 0.0
        return (self.final_value - self.base_value) / self.base_value * 100.0

    @property
    def n_variants(self) -> int:
        return len(self.contributing_variants)

    @property
    def is_modified(self) -> bool:
        return abs(self.fold_change - 1.0) > 1e-9


class EffectComposer:
    """
    Combines genetic variant effects into composed parameter modifiers.

    The composer does NOT modify PathwayParameters directly.
    It returns ComposedEffect objects that the GenotypePersonalizer
    then applies to create a personalized PathwayParameters instance.
    """

    def __init__(self, max_fold_change: float = 5.0) -> None:
        """
        Args:
            max_fold_change: Maximum allowed fold change per parameter.
                             Effects that would exceed this are capped with a warning.
                             Prevents biologically implausible extreme modifications.
        """
        if max_fold_change <= 1.0:
            raise ValueError(f"max_fold_change must be > 1.0, got {max_fold_change}")
        self.max_fold_change = max_fold_change

    def compose(
        self,
        variants: list[Variant],
        base_param_values: dict[str, float],
        effect_db: VariantEffectDB,
        user_population: Population = Population.ALL,
    ) -> list[ComposedEffect]:
        """
        Compose all variant effects into per-parameter modifiers.

        Args:
            variants:          List of observed variants for this individual
            base_param_values: Dict of PathwayParameters field → current value
            effect_db:         The variant effect database to look up in
            user_population:   Individual's ancestry (for mismatch warnings)

        Returns:
            List of ComposedEffect — one per modified parameter.
            Parameters with no variants affecting them are not included.
        """
        # Group variants by parameter field
        by_param: dict[str, list[tuple[EffectRecord, int]]] = {}

        for variant in variants:
            records = effect_db.get(variant.rsid)
            for record in records:
                dosage = variant.dosage_for_allele(record.risk_allele)
                if dosage == 0:
                    continue  # no risk alleles — no effect
                by_param.setdefault(record.parameter_field, []).append(
                    (record, dosage)
                )

        composed = []
        for param_field, record_dosage_pairs in by_param.items():
            base = base_param_values.get(param_field)
            if base is None:
                continue  # parameter not in this model version

            effect = self._compose_parameter(
                param_field=param_field,
                base_value=base,
                record_dosage_pairs=record_dosage_pairs,
                user_population=user_population,
            )
            composed.append(effect)

        return composed

    def _compose_parameter(
        self,
        param_field: str,
        base_value: float,
        record_dosage_pairs: list[tuple[EffectRecord, int]],
        user_population: Population,
    ) -> ComposedEffect:
        """Compose all effects on a single parameter."""
        multiplicative_modifier = 1.0
        additive_offset = 0.0
        contributing: list[tuple[str, int]] = []
        population_warnings: list[str] = []
        genetic_cv_sq = 0.0

        for record, dosage in record_dosage_pairs:
            # Population mismatch warning
            #
            # Regression (privacy, rule 7): this used to embed the raw
            # rsid and gene name directly in the warning text --
            # PersonalizationResult.population_warnings is aggregated
            # straight into SimulationOutput.warnings (via runner.py),
            # which is fully exported/REST-serializable, so this was a
            # live raw-genetic-identifier leak sitting right alongside
            # the contributing_variants leak just fixed elsewhere in
            # this file -- same bug class, a different field, caught by
            # the mandatory security-reviewer pass on that fix, not by
            # this pass independently. Fixed by naming the affected
            # model parameter instead (not privacy-sensitive -- already
            # exposed via PersonalizationResult.summary()'s
            # modified_parameters), losing per-variant attribution when
            # multiple contributing records share a parameter, which is
            # the correct trade-off here.
            if (user_population != Population.ALL
                    and record.population != Population.ALL
                    and record.population != user_population):
                population_warnings.append(
                    f"{param_field}: effect size derived from "
                    f"{record.population.value} population; your ancestry is "
                    f"{user_population.value}. Effect may not generalise."
                )

            # Accumulate effect
            if record.effect_type == EffectType.MULTIPLICATIVE:
                multiplicative_modifier *= (1.0 + record.effect_per_allele * dosage)
            else:
                additive_offset += record.effect_per_allele * dosage

            contributing.append((_rsid_fingerprint(record.rsid), dosage))

            # Accumulate genetic uncertainty (in quadrature)
            uc = record.uncertainty_contribution * dosage
            genetic_cv_sq += uc ** 2

        # Apply combined effect
        final = base_value * multiplicative_modifier + additive_offset

        # Cap extreme fold changes -- symmetric (regression: this used to
        # only catch INCREASES. fold=abs(final/base_value) is always < 1.0
        # for a pure-decrease multiplicative modifier (a product of terms
        # each in (0,1)), so it could never exceed max_fold_change no
        # matter how many decreasing variants stacked. Confirmed with this
        # DB's own real GWAS data: stacking all 9 real secretion_scale-
        # lowering variants homozygous collapses secretion_scale to 0.211
        # (79% loss) with zero warning -- contradicts polygenic-inheritance
        # theory (common SNPs -> individually small effects) and this
        # module's own stated intent just below ("Extreme polygenic
        # combinations may not be physiologically valid"), which was never
        # actually enforced for the decrease direction).
        if base_value != 0:
            fold = final / base_value
            min_fold = 1.0 / self.max_fold_change
            if fold > self.max_fold_change:
                population_warnings.append(
                    f"{param_field}: combined effect would produce {fold:.2f}x change "
                    f"(capped at {self.max_fold_change:.1f}x). "
                    f"Extreme polygenic combinations may not be physiologically valid."
                )
                final = base_value * self.max_fold_change
            elif fold < min_fold:
                population_warnings.append(
                    f"{param_field}: combined effect would produce {fold:.2f}x change "
                    f"(capped at {min_fold:.2f}x). "
                    f"Extreme polygenic combinations may not be physiologically valid."
                )
                final = base_value * min_fold
            # Note: the fold-cap above now also guarantees final > 0 for
            # any base_value > 0 -- a separate "clamp non-positive final
            # to 1% of base" fallback used to live here, but fold < 0
            # (final <= 0) is always caught by the `fold < min_fold`
            # branch first (min_fold > 0 whenever max_fold_change > 1),
            # so it was dead code once the cap became symmetric. Removed
            # rather than left unreachable.

        return ComposedEffect(
            parameter_field=param_field,
            base_value=base_value,
            final_value=final,
            contributing_variants=contributing,
            genetic_cv=math.sqrt(genetic_cv_sq),
            population_warnings=population_warnings,
        )
