# SAVE TO: packages\simulation-api\tests\test_experimental_drug_warnings.py
"""
A simulation that applies a drug says, in its warnings, that the drug module
is EXPERIMENTAL -- and still runs (the label never blocks a drug run). The one
VALIDATED drug (semaglutide, 2026-11-21) carries its validated note and caveat
instead.
"""
from biosim_api import DrugEvent, MealEvent, ScenarioDefinition, SimulationRunner
from biosim_physiology.drug_interactions import Drug
from biosim_physiology.drug_interactions.validation_status import (
    EXPERIMENTAL_DRUG_NOTES,
    VALIDATED_DRUG_NOTES,
)


def _run(drug_events):
    return SimulationRunner().run(ScenarioDefinition(
        name="drug_label", duration_s=7200,
        meal_events=[MealEvent(time_s=0, carbs_g=50.0)],
        drug_events=drug_events,
    ))


class TestExperimentalDrugWarnings:
    def test_drug_run_carries_its_experimental_note(self):
        out = _run([DrugEvent(time_s=0, drug="metformin", dose_mg=1000.0)])
        assert out.success
        assert EXPERIMENTAL_DRUG_NOTES[Drug.METFORMIN] in out.warnings

    def test_each_drug_gets_its_own_note_once(self):
        out = _run([
            DrugEvent(time_s=0, drug="caffeine", dose_mg=200.0),
            DrugEvent(time_s=3600, drug="caffeine", dose_mg=100.0),
            DrugEvent(time_s=0, drug="semaglutide", dose_mg=1.0),
        ])
        assert out.warnings.count(EXPERIMENTAL_DRUG_NOTES[Drug.CAFFEINE]) == 1
        assert VALIDATED_DRUG_NOTES[Drug.SEMAGLUTIDE] in out.warnings

    def test_validated_drug_carries_caveat_not_experimental_note(self):
        out = _run([DrugEvent(time_s=0, drug="semaglutide", dose_mg=1.0)])
        assert out.success
        assert VALIDATED_DRUG_NOTES[Drug.SEMAGLUTIDE] in out.warnings
        assert not any("EXPERIMENTAL drug module" in w for w in out.warnings)

    def test_zero_dose_applies_nothing_so_no_note(self):
        out = _run([DrugEvent(time_s=0, drug="metformin", dose_mg=0.0)])
        assert EXPERIMENTAL_DRUG_NOTES[Drug.METFORMIN] not in out.warnings

    def test_insulin_dose_event_carries_its_note(self):
        from biosim_api import InsulinDoseEvent
        out = SimulationRunner().run(ScenarioDefinition(
            name="lispro", duration_s=7200,
            meal_events=[MealEvent(time_s=0, carbs_g=50.0)],
            insulin_dose_events=[InsulinDoseEvent(time_s=0, insulin_type="lispro", dose_units=4.0)],
        ))
        assert EXPERIMENTAL_DRUG_NOTES[Drug.INSULIN_LISPRO] in out.warnings

    def test_no_drug_no_experimental_note(self):
        out = _run([])
        assert not any("EXPERIMENTAL drug module" in w for w in out.warnings)
        assert not any("VALIDATED drug module" in w for w in out.warnings)
