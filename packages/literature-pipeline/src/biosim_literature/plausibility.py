"""
Plausibility checker for extracted parameters.

Checks extracted parameter values against physiologically realistic ranges.
The goal is not to reject correct measurements — the ranges are intentionally
very wide. The goal is to catch unit conversion errors, transcription errors,
and order-of-magnitude mistakes before they corrupt the model.

A value that fails plausibility is not rejected outright. It is flagged with
a PlausibilityWarning and routes to the expert curation queue. Curators can
override the check with documented justification.

Adding new ranges:
    Add to PLAUSIBLE_RANGES_UM below.
    All values are in canonical units (post-normalization).
    Document the biological basis for the range in a comment.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Optional

from biosim_parameter_db import ParameterKind


class PlausibilityStatus(str, Enum):
    PASS    = "pass"
    WARNING = "warning"   # Outside expected range — flag for review
    FAIL    = "fail"      # Physically impossible (negative Km, Hill > 10, etc.)


@dataclass
class PlausibilityResult:
    status: PlausibilityStatus
    parameter_kind: ParameterKind
    value: float
    canonical_unit: str
    lower_bound: Optional[float]
    upper_bound: Optional[float]
    message: str

    @property
    def passed(self) -> bool:
        return self.status == PlausibilityStatus.PASS

    @property
    def needs_review(self) -> bool:
        return self.status in (PlausibilityStatus.WARNING, PlausibilityStatus.FAIL)


# ---------------------------------------------------------------------------
# Plausible ranges in canonical units (post-normalization)
# (lower_bound, upper_bound)
# These are VERY wide by design — we catch errors, not unusual biology.
# ---------------------------------------------------------------------------

PLAUSIBLE_RANGES: dict[ParameterKind, tuple[float, float]] = {
    # Km, Ki, Kd, EC50, IC50: 1 fM to 100 mM in μM
    # Outside this range almost certainly indicates a unit error
    ParameterKind.KM:          (1e-6, 1e8),
    ParameterKind.KI:          (1e-6, 1e8),
    ParameterKind.KD:          (1e-6, 1e8),
    ParameterKind.EC50:        (1e-6, 1e8),
    ParameterKind.IC50:        (1e-6, 1e8),

    # Vmax in μM/s: 1 fM/s to 1 mM/s
    ParameterKind.VMAX:        (1e-9, 1e3),

    # kcat in s⁻¹: 1e-6 to 1e6 (covers enzymes from very slow to diffusion-limited)
    ParameterKind.KCAT:        (1e-6, 1e6),

    # Hill coefficient: 0.1 to 10
    # Below 0.1 or above 10 is not impossible but is extremely unusual
    ParameterKind.HILL:        (0.1, 10.0),

    # Half-life in seconds: 1 ms to 100 years
    ParameterKind.HALF_LIFE:   (1e-3, 3.15e9),

    # Basal abundance in μM: 1 fM to 1 mM
    ParameterKind.ABUNDANCE:   (1e-9, 1e3),

    # Generic rate in s⁻¹
    ParameterKind.RATE:        (1e-9, 1e6),

    # Fold change: -1000 to +1000 (dimensionless, can be negative)
    ParameterKind.FOLD_CHANGE: (-1000.0, 1000.0),
}

# Hard physical limits — values outside these are impossible, not just unusual
# These trigger FAIL rather than WARNING
PHYSICAL_LIMITS: dict[ParameterKind, tuple[Optional[float], Optional[float]]] = {
    ParameterKind.KM:        (0.0, None),   # Must be positive
    ParameterKind.KI:        (0.0, None),
    ParameterKind.KD:        (0.0, None),
    ParameterKind.EC50:      (0.0, None),
    ParameterKind.IC50:      (0.0, None),
    ParameterKind.VMAX:      (0.0, None),
    ParameterKind.KCAT:      (0.0, None),
    ParameterKind.HILL:      (0.0, None),   # Hill coefficient must be positive
    ParameterKind.HALF_LIFE: (0.0, None),
    ParameterKind.ABUNDANCE: (0.0, None),
    ParameterKind.RATE:      (0.0, None),
}


class PlausibilityChecker:
    """
    Checks whether a normalized parameter value is physiologically plausible.

    Two-tier checking:
        1. Physical limits: impossible values → FAIL (e.g. negative Km)
        2. Plausible ranges: unusual but not impossible → WARNING

    FAIL and WARNING both route to curation queue.
    Only PASS values are automatically accepted.
    """

    def check(
        self,
        kind: ParameterKind,
        value: float,
        canonical_unit: str,
    ) -> PlausibilityResult:
        """
        Check a single parameter value for plausibility.

        Args:
            kind: The type of kinetic parameter
            value: The value in canonical units (post-normalization)
            canonical_unit: The canonical unit string (for reporting)

        Returns:
            PlausibilityResult with status, bounds, and human-readable message
        """
        # Check physical limits first
        if kind in PHYSICAL_LIMITS:
            lower_limit, upper_limit = PHYSICAL_LIMITS[kind]
            if lower_limit is not None and value <= lower_limit:
                return PlausibilityResult(
                    status=PlausibilityStatus.FAIL,
                    parameter_kind=kind,
                    value=value,
                    canonical_unit=canonical_unit,
                    lower_bound=lower_limit,
                    upper_bound=upper_limit,
                    message=(
                        f"{kind.value} value {value} {canonical_unit} is physically impossible "
                        f"(must be > {lower_limit}). This is almost certainly a data entry error."
                    ),
                )
            if upper_limit is not None and value >= upper_limit:
                return PlausibilityResult(
                    status=PlausibilityStatus.FAIL,
                    parameter_kind=kind,
                    value=value,
                    canonical_unit=canonical_unit,
                    lower_bound=lower_limit,
                    upper_bound=upper_limit,
                    message=(
                        f"{kind.value} value {value} {canonical_unit} exceeds physical limit "
                        f"of {upper_limit}. This is almost certainly a data entry error."
                    ),
                )

        # Check plausible ranges
        if kind in PLAUSIBLE_RANGES:
            lower, upper = PLAUSIBLE_RANGES[kind]
            if not (lower <= value <= upper):
                return PlausibilityResult(
                    status=PlausibilityStatus.WARNING,
                    parameter_kind=kind,
                    value=value,
                    canonical_unit=canonical_unit,
                    lower_bound=lower,
                    upper_bound=upper,
                    message=(
                        f"{kind.value} value {value} {canonical_unit} is outside the expected "
                        f"plausible range [{lower}, {upper}] {canonical_unit}. "
                        f"This may be correct for an unusual biological system, "
                        f"or may indicate a unit conversion error. "
                        f"Flagged for expert review."
                    ),
                )

        # For CUSTOM and unranged kinds, pass with a note
        return PlausibilityResult(
            status=PlausibilityStatus.PASS,
            parameter_kind=kind,
            value=value,
            canonical_unit=canonical_unit,
            lower_bound=PLAUSIBLE_RANGES.get(kind, (None, None))[0],
            upper_bound=PLAUSIBLE_RANGES.get(kind, (None, None))[1],
            message=f"{kind.value} value {value} {canonical_unit} passed plausibility check.",
        )

    def check_batch(
        self,
        items: list[tuple[ParameterKind, float, str]],
    ) -> list[PlausibilityResult]:
        """
        Check multiple parameters at once.

        Args:
            items: List of (kind, value, canonical_unit) tuples

        Returns:
            List of PlausibilityResult in the same order as input
        """
        return [self.check(kind, value, unit) for kind, value, unit in items]

    def summary(self, results: list[PlausibilityResult]) -> dict:
        """Return a summary dict for a batch of plausibility results."""
        return {
            "total": len(results),
            "passed": sum(1 for r in results if r.status == PlausibilityStatus.PASS),
            "warnings": sum(1 for r in results if r.status == PlausibilityStatus.WARNING),
            "failed": sum(1 for r in results if r.status == PlausibilityStatus.FAIL),
            "needs_review": [r for r in results if r.needs_review],
        }
