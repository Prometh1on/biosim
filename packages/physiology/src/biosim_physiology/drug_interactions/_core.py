# SAVE TO: packages\physiology\src\biosim_physiology\drug_interactions\_core.py
"""Drug enum and the DrugDose value object."""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum


class Drug(str, Enum):
    CAFFEINE        = "caffeine"
    METFORMIN       = "metformin"
    GLUCOCORTICOID  = "glucocorticoid"
    ALCOHOL         = "alcohol"
    SEMAGLUTIDE     = "semaglutide"
    LIRAGLUTIDE     = "liraglutide"
    DPP4_INHIBITOR  = "dpp4_inhibitor"
    DAPAGLIFLOZIN   = "dapagliflozin"
    EMPAGLIFLOZIN   = "empagliflozin"
    GLIBENCLAMIDE   = "glibenclamide"
    GLIPIZIDE       = "glipizide"
    GLIMEPIRIDE     = "glimepiride"
    TIRZEPATIDE     = "tirzepatide"
    ACARBOSE        = "acarbose"
    MIGLITOL        = "miglitol"
    INSULIN_REGULAR = "insulin_regular"
    INSULIN_LISPRO  = "insulin_lispro"
    PIOGLITAZONE    = "pioglitazone"
    INDINAVIR       = "indinavir"
    HYDROCHLOROTHIAZIDE = "hydrochlorothiazide"
    GROWTH_HORMONE  = "growth_hormone"
    ATORVASTATIN    = "atorvastatin"
    HYDROXYCHLOROQUINE = "hydroxychloroquine"

    @property
    def validation_status(self) -> "DrugValidationStatus":
        """EXPERIMENTAL or VALIDATED, from validation_status.DRUG_VALIDATION
        (the single source). KeyError if a new drug was not added there."""
        from .validation_status import DRUG_VALIDATION
        return DRUG_VALIDATION[self].status


@dataclass(frozen=True)
class DrugDose:
    """
    A single drug administration with dose and timing.

    drug:                   Which drug
    dose_mg:                Dose in milligrams
    hours_since_ingestion:  Time since ingestion at simulation start.
                            0 = just taken; for chronic drugs use a large number
                            to indicate steady-state dosing.
    is_chronic:             True for medications taken daily (metformin).
                            Steady-state effect applied, timing ignored.
    """
    drug:                  Drug
    dose_mg:               float
    hours_since_ingestion: float = 0.0
    is_chronic:            bool  = False

    def __post_init__(self) -> None:
        # Regression (found by security review, glucocorticoid PK/PD
        # milestone -- see CLAUDE.md): "< 0" does not reject NaN (NaN
        # comparisons are always False) or +inf -- fixed at the matching
        # choke point in DrugEvent.__post_init__ (biosim_api/scenario.py)
        # too, but DrugDose is independently constructable (not every
        # caller goes through DrugEvent.to_drug_dose()), so this needs its
        # own guard, not just the upstream one.
        if not math.isfinite(self.dose_mg):
            raise ValueError(f"dose_mg must be finite, got {self.dose_mg}")
        if self.dose_mg < 0:
            raise ValueError(f"dose_mg cannot be negative, got {self.dose_mg}")
        if not math.isfinite(self.hours_since_ingestion):
            raise ValueError(f"hours_since_ingestion must be finite, got {self.hours_since_ingestion}")
        if self.hours_since_ingestion < 0:
            raise ValueError(f"hours_since_ingestion cannot be negative")
