# SAVE TO: packages\physiology\src\biosim_physiology\drug_interactions\validation_status.py
"""
Validation status of every drug module (single source of truth).

Each Drug gets one record here; `Drug.validation_status` reads it, the
simulation runner turns it into a warning for any scenario that applies the
drug, and the validator labels drug benchmarks with it.

A drug module is VALIDATED only if it meets VALIDATED_CRITERION below AND the
result has been reviewed. Semaglutide is the one VALIDATED drug: since B7
(2026-11-16) its single point meets the numeric criterion (+0.7%), and on
2026-11-21 the user reviewed it and chose VALIDATED knowing the match is not
robust (docs/VALIDATION.md, "Semaglutide review"). That caveat is stored with
the record and shown wherever the label is. Every other drug fails the
criterion and stays EXPERIMENTAL. The % errors are copied from
ModelValidator.validate_all() after B7 (2026-11-16); a test in
packages/validation/tests compares them with the live validator, so they
cannot silently go stale.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Iterable

from ._core import Drug


class DrugValidationStatus(str, Enum):
    EXPERIMENTAL = "EXPERIMENTAL"
    VALIDATED = "VALIDATED"


VALIDATED_CRITERION = (
    "every point of a cited clinical benchmark for this drug is within "
    "±20% of the published value, that benchmark was not used to fit "
    "any part of the model, and a person has reviewed the result"
)


@dataclass(frozen=True)
class DrugValidation:
    """
    status:             EXPERIMENTAL or VALIDATED (see VALIDATED_CRITERION).
    benchmark:          ModelValidator dataset name, or "" if none exists.
    deviations_percent: % error of each benchmark point (model vs published).
    summary:            one plain-English line on the evidence.
    reviewed:           a person has reviewed the result (needed for VALIDATED).
    caveat:             limits of a VALIDATED result, shown with the label.
    """
    status: DrugValidationStatus
    benchmark: str = ""
    deviations_percent: tuple[float, ...] = ()
    summary: str = "no clinical benchmark is registered for it"
    reviewed: bool = False
    caveat: str = ""


_E = DrugValidationStatus.EXPERIMENTAL
_NO_BENCHMARK = DrugValidation(status=_E)

DRUG_VALIDATION: dict[Drug, DrugValidation] = {
    # User decision 2026-11-21 ("Validate anyway"), after the robustness
    # check in docs/VALIDATION.md. Label change only: no parameter was tuned.
    Drug.SEMAGLUTIDE: DrugValidation(
        DrugValidationStatus.VALIDATED, "Semaglutide_1mg_postprandial", (0.7,),
        "its one benchmark point (glucose iAUC reduction, Hjerpsted 2018, "
        "1 mg, n=30) is within 1% of the published value; reviewed 2026-11-21",
        reviewed=True,
        caveat=(
            "The match depends on population assumptions: it holds for the "
            "default healthy adult and default 75 g meal, but modelling the "
            "trial's obese population (Hjerpsted et al. 2018, PMID:28941314; "
            "BMI 30 / 35) gives -26% / -39%, and 50 g / 100 g meals give "
            "-34% / +34%. One point, on an OGTT-fitted meal curve."
        ),
    ),
    Drug.GLUCOCORTICOID: DrugValidation(
        _E, "Prednisone_40mg_T2D_postprandial", (-68.2,),
        "its benchmark (prednisone 40 mg in T2D) misses by 68%",
    ),
    Drug.METFORMIN: DrugValidation(
        _E, "Metformin_T2D_postprandial", (-30.6,),
        "its benchmark is 31% low (inside a wide ±50% band only)",
    ),
    Drug.DAPAGLIFLOZIN: DrugValidation(
        _E, "Canagliflozin_300mg_postprandial", (-84.9, -81.3),
        "the SGLT2-inhibitor benchmark (canagliflozin data, run with "
        "dapagliflozin) misses by 81-85%",
    ),
    Drug.EMPAGLIFLOZIN: DrugValidation(
        _E, summary="no benchmark of its own; the SGLT2-inhibitor class "
        "benchmark misses by 81-85%",
    ),
    Drug.CAFFEINE: DrugValidation(
        _E, "Caffeine_5mgkg_OGTT_insulin", (-97.8, -99.7),
        "its benchmark's insulin rise is almost entirely missing (about +1% "
        "modelled vs +25-60% measured)",
    ),
    Drug.ALCOHOL: _NO_BENCHMARK,
    Drug.LIRAGLUTIDE: _NO_BENCHMARK,
    Drug.DPP4_INHIBITOR: _NO_BENCHMARK,
    Drug.GLIBENCLAMIDE: _NO_BENCHMARK,
    Drug.GLIPIZIDE: _NO_BENCHMARK,
    Drug.GLIMEPIRIDE: _NO_BENCHMARK,
    Drug.TIRZEPATIDE: _NO_BENCHMARK,
    Drug.ACARBOSE: _NO_BENCHMARK,
    Drug.MIGLITOL: _NO_BENCHMARK,
    Drug.INSULIN_REGULAR: _NO_BENCHMARK,
    Drug.INSULIN_LISPRO: _NO_BENCHMARK,
    Drug.PIOGLITAZONE: _NO_BENCHMARK,
    Drug.INDINAVIR: _NO_BENCHMARK,
    Drug.HYDROCHLOROTHIAZIDE: _NO_BENCHMARK,
    Drug.GROWTH_HORMONE: _NO_BENCHMARK,
    Drug.ATORVASTATIN: _NO_BENCHMARK,
    Drug.HYDROXYCHLOROQUINE: _NO_BENCHMARK,
}


def _note(drug: Drug, record: DrugValidation) -> str:
    return (
        f"EXPERIMENTAL drug module '{drug.value}': not validated against "
        f"clinical data ({record.summary}). Treat its effect on glucose and "
        "insulin as illustrative, not a prediction."
    )


def _validated_note(drug: Drug, record: DrugValidation) -> str:
    return (
        f"VALIDATED drug module '{drug.value}' ({record.summary}). Caveat: "
        f"{record.caveat}"
    )


# User-facing disclosure per drug, in the same style as
# biosim_genotype.effect_db.UNVERIFIED_VARIANT_NOTES.
EXPERIMENTAL_DRUG_NOTES: dict[Drug, str] = {
    drug: _note(drug, record)
    for drug, record in DRUG_VALIDATION.items()
    if record.status is DrugValidationStatus.EXPERIMENTAL
}

# Validated drugs still carry a note, so the caveat travels with the label.
VALIDATED_DRUG_NOTES: dict[Drug, str] = {
    drug: _validated_note(drug, record)
    for drug, record in DRUG_VALIDATION.items()
    if record.status is DrugValidationStatus.VALIDATED
}


def _notes_for(drug_names: Iterable[str], table: dict[Drug, str]) -> list[str]:
    notes: list[str] = []
    for name in drug_names:
        note = table.get(Drug(name))
        if note is not None and note not in notes:
            notes.append(note)
    return notes


def experimental_drug_notes(drug_names: Iterable[str]) -> list[str]:
    """One note per distinct experimental drug, in first-seen order."""
    return _notes_for(drug_names, EXPERIMENTAL_DRUG_NOTES)


def drug_validation_notes(drug_names: Iterable[str]) -> list[str]:
    """One note per distinct drug, experimental or validated, first-seen order."""
    return _notes_for(drug_names, {**EXPERIMENTAL_DRUG_NOTES, **VALIDATED_DRUG_NOTES})
