# SAVE TO: packages\physiology\tests\test_drug_validation_status.py
"""
Every drug the model can apply carries a validation status, from one table.

Since 2026-11-16 every drug module was EXPERIMENTAL. On 2026-11-21 the user
reviewed semaglutide and promoted it to VALIDATED, with a recorded caveat:
its match depends on population assumptions (modelling the trial's obese
population gives -26% to -39%). Every other drug stays EXPERIMENTAL. A drug
added to the Drug enum without a status makes these tests fail.
"""
import pytest

from biosim_physiology.drug_interactions import Drug
from biosim_physiology.drug_interactions.validation_status import (
    DRUG_VALIDATION,
    EXPERIMENTAL_DRUG_NOTES,
    VALIDATED_CRITERION,
    VALIDATED_DRUG_NOTES,
    DrugValidationStatus,
    drug_validation_notes,
    experimental_drug_notes,
)


class TestEveryDrugHasAStatus:
    @pytest.mark.parametrize("drug", list(Drug))
    def test_drug_has_validation_status(self, drug):
        assert isinstance(drug.validation_status, DrugValidationStatus)

    def test_table_covers_exactly_the_drug_enum(self):
        assert set(DRUG_VALIDATION) == set(Drug)

    @pytest.mark.parametrize("drug", list(Drug))
    def test_status_comes_from_the_table(self, drug):
        assert drug.validation_status is DRUG_VALIDATION[drug].status


class TestValidatedCriterion:
    def test_criterion_is_stated_and_strict(self):
        assert "20%" in VALIDATED_CRITERION
        assert "not used" in VALIDATED_CRITERION

    @pytest.mark.parametrize("drug", list(Drug))
    def test_validated_only_if_every_point_is_within_20_percent(self, drug):
        record = DRUG_VALIDATION[drug]
        if record.status is DrugValidationStatus.VALIDATED:
            assert record.deviations_percent
            assert all(abs(d) <= 20.0 for d in record.deviations_percent)

    @pytest.mark.parametrize("drug", list(Drug))
    def test_validated_needs_a_review_and_a_caveat(self, drug):
        record = DRUG_VALIDATION[drug]
        if record.status is DrugValidationStatus.VALIDATED:
            assert record.reviewed
            assert record.caveat

    def test_only_semaglutide_is_validated(self):
        validated = {
            d for d, r in DRUG_VALIDATION.items()
            if r.status is DrugValidationStatus.VALIDATED
        }
        assert validated == {Drug.SEMAGLUTIDE}

    def test_semaglutide_caveat_is_recorded(self):
        # User decision 2026-11-21 ("Validate anyway") after seeing that
        # modelling Hjerpsted 2018's obese population gives -26% to -39%.
        record = DRUG_VALIDATION[Drug.SEMAGLUTIDE]
        assert record.status is DrugValidationStatus.VALIDATED
        assert record.deviations_percent == (0.7,)
        assert "28941314" in record.caveat
        assert "obese" in record.caveat
        assert "-26%" in record.caveat and "-39%" in record.caveat
        note = VALIDATED_DRUG_NOTES[Drug.SEMAGLUTIDE]
        assert "VALIDATED drug module 'semaglutide'" in note
        assert record.caveat in note


class TestExperimentalDrugNotes:
    @pytest.mark.parametrize("drug", list(Drug))
    def test_every_drug_has_exactly_one_note_naming_it(self, drug):
        if drug.validation_status is DrugValidationStatus.EXPERIMENTAL:
            note = EXPERIMENTAL_DRUG_NOTES[drug]
            assert "EXPERIMENTAL" in note
            assert drug not in VALIDATED_DRUG_NOTES
        else:
            note = VALIDATED_DRUG_NOTES[drug]
            assert drug not in EXPERIMENTAL_DRUG_NOTES
        assert drug.value in note

    def test_notes_for_listed_drugs_only_once_each(self):
        notes = experimental_drug_notes(["metformin", "metformin", "caffeine"])
        assert notes == [
            EXPERIMENTAL_DRUG_NOTES[Drug.METFORMIN],
            EXPERIMENTAL_DRUG_NOTES[Drug.CAFFEINE],
        ]

    def test_experimental_notes_skip_validated_drugs(self):
        assert experimental_drug_notes(["semaglutide"]) == []

    def test_validation_notes_cover_both_labels(self):
        notes = drug_validation_notes(["semaglutide", "metformin", "semaglutide"])
        assert notes == [
            VALIDATED_DRUG_NOTES[Drug.SEMAGLUTIDE],
            EXPERIMENTAL_DRUG_NOTES[Drug.METFORMIN],
        ]

    def test_no_drugs_no_notes(self):
        assert experimental_drug_notes([]) == []
        assert drug_validation_notes([]) == []


def test_caffeine_has_no_fasting_glucose_term():
    # Drug stream (2026-11-20): the uncited +2% Gb term was removed. In this
    # model a Gb change also shifts the beta-cell reference (a chronic set
    # point), so for an acute drug it LOWERED insulin -- a wrong-sign artefact.
    # Keijzers 2002 (PMID:11815511) reports only reduced insulin sensitivity.
    from biosim_physiology.drug_interactions import DrugDose, caffeine_modifiers
    mods = caffeine_modifiers(DrugDose(drug=Drug.CAFFEINE, dose_mg=350.0, hours_since_ingestion=1.0))
    assert set(mods) == {"p3"}
    assert mods["p3"] < 1.0
