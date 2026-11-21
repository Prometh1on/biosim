"""
Unit normalizer for the literature ingestion pipeline.

All kinetic parameters entering the platform must be converted to canonical units:

    Concentrations       : micromolar (μM)
    Time                 : seconds (s)
    First-order rates    : per second (s⁻¹)
    Second-order rates   : per micromolar per second (μM⁻¹s⁻¹)
    Reaction rates       : micromolar per second (μM/s)
    Half-life            : seconds (s)
    Dimensionless        : no conversion (Hill coefficients, fold changes, ratios)

Unit conversion is explicit and exhaustive — unknown units raise NormalizationError
rather than passing through silently. A silent pass-through of an unknown unit
is a silent error in every downstream simulation.

Adding a new unit:
    Add an entry to the appropriate section of UNIT_CONVERSIONS below.
    The value is the multiplier to apply to convert TO the canonical unit.
    Document the addition in docs/unit-standards.md.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


class NormalizationError(Exception):
    """
    Raised when a unit cannot be normalized.

    Always raised rather than silently passing through — an unknown unit
    in a kinetic parameter is a data integrity problem, not an edge case.
    """
    pass


# ---------------------------------------------------------------------------
# Unit conversion tables
# Each entry: input_unit_string → (multiplier, canonical_unit_string)
# normalized_value = raw_value * multiplier
# ---------------------------------------------------------------------------

# Concentrations → μM
CONCENTRATION_UNITS: dict[str, tuple[float, str]] = {
    "M":   (1e6,  "uM"),
    "mM":  (1e3,  "uM"),
    "uM":  (1.0,  "uM"),
    "μM":  (1.0,  "uM"),
    "nM":  (1e-3, "uM"),
    "pM":  (1e-6, "uM"),
    "fM":  (1e-9, "uM"),
}

# Time → s
TIME_UNITS: dict[str, tuple[float, str]] = {
    "ms":  (1e-3,       "s"),
    "s":   (1.0,        "s"),
    "min": (60.0,       "s"),
    "h":   (3600.0,     "s"),
    "hr":  (3600.0,     "s"),
    "d":   (86400.0,    "s"),
    "day": (86400.0,    "s"),
}

# First-order rate constants → s⁻¹
FIRST_ORDER_RATE_UNITS: dict[str, tuple[float, str]] = {
    "s-1":    (1.0,        "s-1"),
    "s⁻¹":   (1.0,        "s-1"),
    "min-1":  (1/60,       "s-1"),
    "min⁻¹":  (1/60,       "s-1"),
    "h-1":    (1/3600,     "s-1"),
    "h⁻¹":   (1/3600,     "s-1"),
    "hr-1":   (1/3600,     "s-1"),
    "d-1":    (1/86400,    "s-1"),
}

# Second-order rate constants → μM⁻¹s⁻¹
SECOND_ORDER_RATE_UNITS: dict[str, tuple[float, str]] = {
    "M-1s-1":      (1e-6,       "uM-1s-1"),
    "M⁻¹s⁻¹":     (1e-6,       "uM-1s-1"),
    "mM-1s-1":     (1e-3,       "uM-1s-1"),
    "uM-1s-1":     (1.0,        "uM-1s-1"),
    "μM-1s-1":     (1.0,        "uM-1s-1"),
    "uM-1min-1":   (1/60,       "uM-1s-1"),
    "nM-1s-1":     (1e3,        "uM-1s-1"),
    "M-1min-1":    (1e-6/60,    "uM-1s-1"),
    "M-1h-1":      (1e-6/3600,  "uM-1s-1"),
}

# Reaction rates (Vmax, flux) → μM/s
REACTION_RATE_UNITS: dict[str, tuple[float, str]] = {
    "M/s":    (1e6,       "uM/s"),
    "mM/s":   (1e3,       "uM/s"),
    "uM/s":   (1.0,       "uM/s"),
    "μM/s":   (1.0,       "uM/s"),
    "nM/s":   (1e-3,      "uM/s"),
    "uM/min": (1/60,      "uM/s"),
    "mM/min": (1e3/60,    "uM/s"),
    "uM/h":   (1/3600,    "uM/s"),
    "mM/h":   (1e3/3600,  "uM/s"),
    "nM/min": (1e-3/60,   "uM/s"),
    "nM/h":   (1e-3/3600, "uM/s"),
}

# Dimensionless — no conversion, just validate and tag
DIMENSIONLESS_UNITS: set[str] = {
    "dimensionless",
    "fold_change",
    "fold-change",
    "ratio",
    "unitless",
    "",
}

# Half-life uses time units
HALF_LIFE_UNITS = TIME_UNITS

# Merged lookup for auto-detection
ALL_UNIT_TABLES = [
    CONCENTRATION_UNITS,
    TIME_UNITS,
    FIRST_ORDER_RATE_UNITS,
    SECOND_ORDER_RATE_UNITS,
    REACTION_RATE_UNITS,
]


@dataclass(frozen=True)
class NormalizedValue:
    """Result of a successful unit normalization."""
    value: float
    canonical_unit: str
    original_value: float
    original_unit: str
    conversion_factor: float


class UnitNormalizer:
    """
    Converts raw parameter values to canonical platform units.

    Usage:
        normalizer = UnitNormalizer()

        # Known unit type
        result = normalizer.normalize_concentration(45.0, "nM")
        # → NormalizedValue(value=0.045, canonical_unit="uM", ...)

        # Auto-detect unit type
        result = normalizer.normalize(45.0, "nM")
        # → NormalizedValue(value=0.045, canonical_unit="uM", ...)
    """

    def _convert(
        self,
        value: float,
        unit: str,
        table: dict[str, tuple[float, str]],
        unit_family: str,
    ) -> NormalizedValue:
        """Internal helper: look up unit in table and apply conversion."""
        unit_stripped = unit.strip()
        if unit_stripped not in table:
            raise NormalizationError(
                f"Unknown {unit_family} unit '{unit}'. "
                f"Known {unit_family} units: {sorted(table.keys())}. "
                f"To add it, update the appropriate table in normalizer.py "
                f"and document in docs/unit-standards.md."
            )
        factor, canonical = table[unit_stripped]
        return NormalizedValue(
            value=value * factor,
            canonical_unit=canonical,
            original_value=value,
            original_unit=unit_stripped,
            conversion_factor=factor,
        )

    def normalize_concentration(self, value: float, unit: str) -> NormalizedValue:
        return self._convert(value, unit, CONCENTRATION_UNITS, "concentration")

    def normalize_time(self, value: float, unit: str) -> NormalizedValue:
        return self._convert(value, unit, TIME_UNITS, "time")

    def normalize_first_order_rate(self, value: float, unit: str) -> NormalizedValue:
        return self._convert(value, unit, FIRST_ORDER_RATE_UNITS, "first-order rate")

    def normalize_second_order_rate(self, value: float, unit: str) -> NormalizedValue:
        return self._convert(value, unit, SECOND_ORDER_RATE_UNITS, "second-order rate")

    def normalize_reaction_rate(self, value: float, unit: str) -> NormalizedValue:
        return self._convert(value, unit, REACTION_RATE_UNITS, "reaction rate")

    def normalize_half_life(self, value: float, unit: str) -> NormalizedValue:
        return self._convert(value, unit, HALF_LIFE_UNITS, "half-life")

    def normalize_dimensionless(self, value: float, unit: str) -> NormalizedValue:
        if unit not in DIMENSIONLESS_UNITS:
            raise NormalizationError(
                f"Unit '{unit}' is not recognised as dimensionless. "
                f"Known dimensionless units: {sorted(DIMENSIONLESS_UNITS)}. "
                f"If this is correct, add it to DIMENSIONLESS_UNITS in normalizer.py "
                f"and document it in docs/unit-standards.md."
            )
        return NormalizedValue(
            value=value,
            canonical_unit="dimensionless",
            original_value=value,
            original_unit=unit,
            conversion_factor=1.0,
        )

    def normalize(self, value: float, unit: str) -> NormalizedValue:
        """
        Auto-detect the unit type and normalize.

        Searches all unit tables in order. Raises NormalizationError if the
        unit is not found in any table. Use the specific methods above when
        the unit type is known — auto-detection is for ingestion pipelines
        where the unit type may not be explicitly declared.
        """
        unit_stripped = unit.strip()

        if unit_stripped in DIMENSIONLESS_UNITS:
            return self.normalize_dimensionless(value, unit_stripped)

        for table in ALL_UNIT_TABLES:
            if unit_stripped in table:
                factor, canonical = table[unit_stripped]
                return NormalizedValue(
                    value=value * factor,
                    canonical_unit=canonical,
                    original_value=value,
                    original_unit=unit_stripped,
                    conversion_factor=factor,
                )

        raise NormalizationError(
            f"Unknown unit '{unit}'. "
            f"This unit is not in the platform's unit registry. "
            f"To add it: update the appropriate table in normalizer.py "
            f"and document the addition in docs/unit-standards.md. "
            f"Do NOT convert the value manually and use a known unit — "
            f"silent unit mismatches cause simulation errors."
        )

    def is_known_unit(self, unit: str) -> bool:
        """Return True if the unit is recognised by the normalizer."""
        unit_stripped = unit.strip()
        if unit_stripped in DIMENSIONLESS_UNITS:
            return True
        return any(unit_stripped in table for table in ALL_UNIT_TABLES)

    def canonical_unit_for(self, unit: str) -> Optional[str]:
        """Return the canonical unit for a given input unit, or None if unknown."""
        unit_stripped = unit.strip()
        if unit_stripped in DIMENSIONLESS_UNITS:
            return "dimensionless"
        for table in ALL_UNIT_TABLES:
            if unit_stripped in table:
                return table[unit_stripped][1]
        return None
