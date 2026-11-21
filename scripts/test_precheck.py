# SAVE TO: scripts\test_precheck.py
"""
Tests for scripts/precheck.py (MASTER_PLAN Section 10, item C4).

Uses hand-crafted temporary source snippets rather than the real
codebase for the unit-level checks (so a future real-code change can't
silently make these tests vacuous), plus one integration-level test
confirming the script runs clean against the actual repo today.
"""
import sys
import textwrap
from pathlib import Path

import pytest

SCRIPTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS_DIR))

import precheck  # noqa: E402


def _write(tmp_path: Path, name: str, source: str) -> Path:
    path = tmp_path / name
    path.write_text(textwrap.dedent(source), encoding="utf-8")
    return path


class TestDoseFunctionFinitenessGuard:
    def test_flags_a_raw_float_dose_param_with_no_guard(self, tmp_path):
        path = _write(tmp_path, "bad.py", """
            def some_new_drug_scale(dose_mg: float) -> float:
                if dose_mg <= 0:
                    return 1.0
                return 1.0 - 0.1 * dose_mg
        """)
        findings = precheck.check_dose_functions_have_finiteness_guards(path)
        assert len(findings) == 1
        assert "some_new_drug_scale" in findings[0]

    def test_does_not_flag_a_function_with_an_isfinite_guard(self, tmp_path):
        path = _write(tmp_path, "good.py", """
            import math

            def some_new_drug_scale(dose_mg: float) -> float:
                if not math.isfinite(dose_mg):
                    raise ValueError("bad dose")
                if dose_mg <= 0:
                    return 1.0
                return 1.0 - 0.1 * dose_mg
        """)
        assert precheck.check_dose_functions_have_finiteness_guards(path) == []

    def test_does_not_flag_an_already_validated_drugdose_parameter(self, tmp_path):
        # DrugDose/DrugEvent/InsulinDoseEvent already validate finiteness
        # in their own __post_init__ -- a function taking one of these,
        # not a raw float, needs no guard of its own.
        path = _write(tmp_path, "exempt.py", """
            def some_modifier(dose: DrugDose) -> dict:
                return {"Gb": 1.0}
        """)
        findings = precheck.check_dose_functions_have_finiteness_guards(path)
        # dose is not literally named dose_mg/dose_units, so this
        # particular example wouldn't trigger the check at all -- confirm
        # the exemption ALSO holds when the param is named correctly.
        assert findings == []

    def test_does_not_flag_a_dose_mg_param_typed_as_drugdose(self, tmp_path):
        path = _write(tmp_path, "exempt2.py", """
            def some_modifier(dose_mg: DrugDose) -> dict:
                return {"Gb": 1.0}
        """)
        assert precheck.check_dose_functions_have_finiteness_guards(path) == []

    def test_does_not_flag_an_isnan_style_guard(self, tmp_path):
        path = _write(tmp_path, "isnan.py", """
            import math

            def some_new_drug_scale(dose_mg: float) -> float:
                if math.isnan(dose_mg):
                    raise ValueError("bad dose")
                return 1.0
        """)
        assert precheck.check_dose_functions_have_finiteness_guards(path) == []


class TestDuplicateCustomParameterRegistrations:
    def test_flags_a_genuine_duplicate_add_registration(self, tmp_path):
        path = _write(tmp_path, "dupe.py", """
            store.add(
                ParameterBuilder()
                .for_subject(insr_id)
                .kind(ParameterKind.CUSTOM, custom_name="insulin_action_p3")
                .value(1.0)
                .build()
            )
            store.add(
                ParameterBuilder()
                .for_subject(insr_id)
                .kind(ParameterKind.CUSTOM, custom_name="insulin_action_p3")
                .value(2.0)
                .build()
            )
        """)
        findings = precheck.check_duplicate_custom_parameter_registrations(path)
        assert len(findings) == 1
        assert "insr_id" in findings[0]
        assert "insulin_action_p3" in findings[0]

    def test_does_not_flag_a_legitimate_supersede_pair(self, tmp_path):
        # Real, minimal reproduction of parameter_seed.py's own p3
        # supersede pattern -- this is the exact false positive found
        # and fixed while building this check against the real codebase.
        path = _write(tmp_path, "supersede.py", """
            old_id = store.add(
                ParameterBuilder()
                .for_subject(insr_id)
                .kind(ParameterKind.CUSTOM, custom_name="insulin_action_p3")
                .value(1.0)
                .build()
            )
            store.supersede(
                old_id,
                ParameterBuilder()
                .for_subject(insr_id)
                .kind(ParameterKind.CUSTOM, custom_name="insulin_action_p3")
                .value(2.0)
                .build()
            )
        """)
        assert precheck.check_duplicate_custom_parameter_registrations(path) == []

    def test_does_not_flag_distinct_custom_names_under_the_same_subject(self, tmp_path):
        path = _write(tmp_path, "distinct.py", """
            store.add(
                ParameterBuilder()
                .for_subject(insulin_id)
                .kind(ParameterKind.CUSTOM, custom_name="exogenous_insulin_regular_ka_per_min")
                .value(1.0)
                .build()
            )
            store.add(
                ParameterBuilder()
                .for_subject(insulin_id)
                .kind(ParameterKind.CUSTOM, custom_name="exogenous_insulin_regular_pm_per_unit")
                .value(2.0)
                .build()
            )
        """)
        assert precheck.check_duplicate_custom_parameter_registrations(path) == []

    def test_does_not_flag_the_same_custom_name_under_different_subjects(self, tmp_path):
        path = _write(tmp_path, "different_subjects.py", """
            store.add(
                ParameterBuilder()
                .for_subject(subject_a)
                .kind(ParameterKind.CUSTOM, custom_name="some_shared_name")
                .value(1.0)
                .build()
            )
            store.add(
                ParameterBuilder()
                .for_subject(subject_b)
                .kind(ParameterKind.CUSTOM, custom_name="some_shared_name")
                .value(2.0)
                .build()
            )
        """)
        assert precheck.check_duplicate_custom_parameter_registrations(path) == []

    def test_flags_a_duplicate_split_across_two_files(self, tmp_path):
        # parameter_seed is a package of several modules (2026-10-20 split),
        # so a duplicate can sit in two different files.
        block = """
            store.add(
                ParameterBuilder()
                .for_subject(insr_id)
                .kind(ParameterKind.CUSTOM, custom_name="insulin_action_p3")
                .value(1.0)
                .build()
            )
        """
        first = _write(tmp_path, "a.py", block)
        second = _write(tmp_path, "b.py", block)
        findings = precheck.check_duplicate_custom_parameter_registrations([first, second])
        assert len(findings) == 1
        assert "b.py" in findings[0] and "a.py" in findings[0]


def _fake_repo(tmp_path: Path, monkeypatch, seed_files: dict[str, str]) -> None:
    """A minimal repo layout with the physiology/data-seeds packages split into subpackages."""
    physiology = tmp_path / "packages" / "physiology" / "src" / "biosim_physiology"
    (physiology / "drug_interactions").mkdir(parents=True)
    (physiology / "modifiers").mkdir()
    _write(physiology / "drug_interactions", "new_drug.py", """
        def new_drug_factor(dose_mg: float) -> float:
            return 1.0
    """)
    seeds = tmp_path / "packages" / "data-seeds" / "src" / "biosim_seeds" / "parameter_seed"
    seeds.mkdir(parents=True)
    for name, source in seed_files.items():
        _write(seeds, name, source)
    (tmp_path / "MASTER_PLAN.md").write_text("", encoding="utf-8")
    monkeypatch.setattr(precheck, "REPO_ROOT", tmp_path)


class TestMainCoversSplitPackages:
    _DUPE = """
        store.add(
            ParameterBuilder()
            .for_subject(insr_id)
            .kind(ParameterKind.CUSTOM, custom_name="insulin_action_p3")
            .build()
        )
    """

    def test_checks_files_inside_subpackages(self, tmp_path, monkeypatch, capsys):
        _fake_repo(tmp_path, monkeypatch, {"a.py": self._DUPE, "b.py": self._DUPE})
        assert precheck.main() == 1
        out = capsys.readouterr().out
        assert "new_drug_factor" in out  # nested physiology file was scanned
        assert "duplicate CUSTOM registration" in out  # cross-file seed duplicate caught

    def test_fails_loudly_when_no_seed_files_are_found(self, tmp_path, monkeypatch, capsys):
        # A moved/renamed parameter_seed must not silently skip check 2.
        _fake_repo(tmp_path, monkeypatch, {})
        assert precheck.main() == 1
        assert "parameter_seed" in capsys.readouterr().out


class TestGbModifierCoverageMatrixReminder:
    def test_warns_on_a_gb_modifier_not_mentioned_in_master_plan(self, tmp_path, monkeypatch):
        source_path = _write(tmp_path, "new_drug.py", """
            def brand_new_drug_modifiers(dose):
                return {"Gb": 0.9}
        """)
        master_plan = tmp_path / "MASTER_PLAN.md"
        master_plan.write_text("nothing relevant here", encoding="utf-8")
        monkeypatch.setattr(precheck, "REPO_ROOT", tmp_path)
        warnings = precheck.check_gb_modifiers_are_in_the_coverage_matrix([source_path])
        assert len(warnings) == 1
        assert "brand_new_drug_modifiers" in warnings[0]

    def test_does_not_warn_when_the_function_name_is_already_mentioned(self, tmp_path, monkeypatch):
        source_path = _write(tmp_path, "known_drug.py", """
            def already_documented_modifiers(dose):
                return {"Gb": 0.9}
        """)
        master_plan = tmp_path / "MASTER_PLAN.md"
        master_plan.write_text("See already_documented_modifiers() in Section 11.", encoding="utf-8")
        monkeypatch.setattr(precheck, "REPO_ROOT", tmp_path)
        assert precheck.check_gb_modifiers_are_in_the_coverage_matrix([source_path]) == []

    def test_does_not_warn_on_a_function_with_no_gb_modifier(self, tmp_path, monkeypatch):
        source_path = _write(tmp_path, "unrelated.py", """
            def some_other_function(x):
                return {"p3": 0.9}
        """)
        master_plan = tmp_path / "MASTER_PLAN.md"
        master_plan.write_text("nothing relevant here", encoding="utf-8")
        monkeypatch.setattr(precheck, "REPO_ROOT", tmp_path)
        assert precheck.check_gb_modifiers_are_in_the_coverage_matrix([source_path]) == []


class TestRunsCleanAgainstTheRealRepo:
    def test_main_exits_zero_against_the_actual_codebase(self):
        # Integration check: the real repo has already been swept and
        # the one genuine finding (the p3 supersede false positive) was
        # fixed before this test was written -- confirms the script
        # itself stays runnable and clean against real, live source,
        # not just the synthetic fixtures above.
        assert precheck.main() == 0
