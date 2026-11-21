# SAVE TO: packages\expert-tools\tests\test_expert_tools.py
"""
Test suite for biosim_expert.

Coverage:
    - SBMLExporter: valid XML, all species present, all parameters present
    - BatchRunner: parallel execution, result collection, failure handling
    - CSVExporter: time series export, metrics table, comparison export
"""

import csv
import io
import xml.etree.ElementTree as ET
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import pytest

from biosim_api import (
    MealEvent, ScenarioDefinition, SimulationRunner,
    high_sugar_meal_scenario, light_snack_scenario, fasted_scenario,
    ComparisonRunner,
)
from biosim_expert import BatchResult, BatchRunner, CSVExporter, SBMLExporter
from biosim_ode import GlucoseInsulinPathway


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def meal_output():
    return SimulationRunner().run(high_sugar_meal_scenario())


@pytest.fixture(scope="module")
def fasted_output():
    return SimulationRunner().run(fasted_scenario(duration_s=3600))


@pytest.fixture(scope="module")
def comparison():
    return ComparisonRunner().compare(high_sugar_meal_scenario(), light_snack_scenario())


# ---------------------------------------------------------------------------
# SBMLExporter
# ---------------------------------------------------------------------------

class TestSBMLExporter:
    def test_export_returns_string(self):
        xml = SBMLExporter().export(GlucoseInsulinPathway())
        assert isinstance(xml, str)
        assert len(xml) > 100

    def test_valid_xml(self):
        xml = SBMLExporter().export(GlucoseInsulinPathway())
        # Should parse without errors
        root = ET.fromstring(xml)
        assert root is not None

    def test_sbml_namespace(self):
        xml = SBMLExporter().export(GlucoseInsulinPathway())
        assert "sbml.org/sbml" in xml

    def test_all_species_present(self):
        xml = SBMLExporter().export(GlucoseInsulinPathway())
        for species in ["G", "I", "X", "IGF1", "mTOR", "AMPK"]:
            assert f'id="{species}"' in xml

    def test_all_parameters_present(self):
        xml = SBMLExporter().export(GlucoseInsulinPathway())
        for param in ["p1", "p2", "p3", "ke", "Gb", "Ib", "EC50_mTOR"]:
            assert f'id="{param}"' in xml

    def test_rate_rules_present(self):
        xml = SBMLExporter().export(GlucoseInsulinPathway())
        assert "rateRule" in xml
        assert 'variable="G"' in xml
        assert 'variable="I"' in xml

    def test_initial_conditions_embedded(self):
        p = GlucoseInsulinPathway()
        xml = SBMLExporter().export(p)
        # Basal glucose is 5.0 mM
        assert "5" in xml

    def test_custom_model_id(self):
        xml = SBMLExporter().export(
            GlucoseInsulinPathway(), model_id="my_custom_model"
        )
        assert 'id="my_custom_model"' in xml

    def test_unit_definitions_present(self):
        xml = SBMLExporter().export(GlucoseInsulinPathway())
        assert "listOfUnitDefinitions" in xml
        assert "per_minute" in xml

    def test_personalised_params_exported(self):
        from biosim_ode import PathwayParameters
        import dataclasses
        custom_params = dataclasses.replace(PathwayParameters(), p1=0.05)
        pathway = GlucoseInsulinPathway(params=custom_params)
        xml = SBMLExporter().export(pathway)
        assert "0.05" in xml

    def test_two_compartment_meals_export(self):
        # Since 2026-11-17 the exporter encodes the two-compartment curve
        # (default meals and GI meals); it used to raise.
        from biosim_api import MealEvent
        from biosim_ode import high_sugar_meal
        gi_meal = MealEvent(time_s=0.0, carbs_g=50.0, glycaemic_index=80.0).to_meal_input()
        xml = SBMLExporter().export(GlucoseInsulinPathway(meals=[high_sugar_meal(), gi_meal]))
        assert 'id="meal0_k_emptying"' in xml and 'id="meal1_k_absorption"' in xml


class TestSBMLExportMatchesRealODE:
    """
    Regression for a previously-undocumented correctness gap: _rate_rules()
    encoded only -ke*(I-Ib) for dI/dt (no glucose-dependent secretion at
    all) and no meal Ra(t) forcing term for dG/dt whatsoever -- an
    external tool (COPASI, CellDesigner, etc.) loading this export and
    running any meal challenge would see glucose simply decay toward Gb
    with zero insulin response, a qualitatively different, non-functional
    model rather than a faithful export of this platform's real dynamics.
    Existing tests (above) only checked XML structure (species/parameter
    IDs present, valid XML) -- none ever simulated the exported model or
    numerically compared it against GlucoseInsulinPathway.dydt(), the
    exact blind spot this class closes.

    Verification approach: parse the exported MathML with a small,
    purpose-built evaluator (supporting exactly the operators _rate_rules()
    emits) and compare its output against dydt()'s real return value at
    several (t, y) points -- not a full SBML/libsbml simulation (not
    available in this environment), but a direct, rigorous check that the
    emitted math is actually equivalent to the real ODE, not just
    structurally present.
    """

    @staticmethod
    def _tag(el) -> str:
        return el.tag.rsplit("}", 1)[-1]

    @classmethod
    def _eval(cls, el, context: dict):
        import math as _math
        tag = cls._tag(el)
        if tag == "math":
            return cls._eval(list(el)[0], context)
        if tag == "cn":
            return float(el.text.strip())
        if tag == "ci":
            return context[el.text.strip()]
        if tag == "csymbol":
            return context["t"]
        if tag == "piecewise":
            otherwise = None
            for child in el:
                ctag = cls._tag(child)
                if ctag == "piece":
                    value_el, cond_el = list(child)
                    if cls._eval(cond_el, context):
                        return cls._eval(value_el, context)
                elif ctag == "otherwise":
                    otherwise = list(child)[0]
            if otherwise is not None:
                return cls._eval(otherwise, context)
            raise ValueError("piecewise matched no piece and had no otherwise")
        if tag == "apply":
            children = list(el)
            op = cls._tag(children[0])
            args = children[1:]
            vals = [cls._eval(a, context) for a in args]
            if op == "plus":
                return sum(vals)
            if op == "times":
                result = 1.0
                for v in vals:
                    result *= v
                return result
            if op == "minus":
                return -vals[0] if len(vals) == 1 else vals[0] - vals[1]
            if op == "divide":
                return vals[0] / vals[1]
            if op == "power":
                return vals[0] ** vals[1]
            if op == "exp":
                return _math.exp(vals[0])
            if op == "geq":
                return vals[0] >= vals[1]
            if op == "gt":
                return vals[0] > vals[1]
            raise ValueError(f"unsupported MathML operator: {op}")
        raise ValueError(f"unsupported MathML tag: {tag}")

    @classmethod
    def _rate_rule_values(cls, xml_str: str, context: dict) -> dict:
        root = ET.fromstring(xml_str)
        result = {}
        for rr in root.iter():
            if cls._tag(rr) == "rateRule":
                variable = rr.attrib["variable"]
                math_el = next(c for c in rr if cls._tag(c) == "math")
                result[variable] = cls._eval(math_el, context)
        return result

    @staticmethod
    def _param_context(xml_str: str) -> dict:
        root = ET.fromstring(xml_str)
        context = {}
        for el in root.iter():
            if TestSBMLExportMatchesRealODE._tag(el) == "parameter":
                context[el.attrib["id"]] = float(el.attrib["value"])
        return context

    def _assert_matches_dydt(self, pathway, t: float, y: list):
        # abs=1e-5, not 1e-9: exported parameters are formatted with 6
        # significant figures (e.g. F_secretion_basal, meal D_eff), so a
        # tiny rounding difference vs. the live Python float is expected
        # and real -- still tight enough that a missing/wrong term (which
        # this test class exists to catch) shows up as an error many
        # orders of magnitude larger than this, not a false negative.
        xml_str = SBMLExporter().export(pathway)
        context = self._param_context(xml_str)
        context.update({"G": y[0], "I": y[1], "X": y[2], "IGF1": y[3], "mTOR": y[4], "AMPK": y[5], "t": t})
        mathml_derivs = self._rate_rule_values(xml_str, context)
        real_derivs = pathway.dydt(t, y)
        species = ["G", "I", "X", "IGF1", "mTOR", "AMPK"]
        for i, name in enumerate(species):
            assert mathml_derivs[name] == pytest.approx(real_derivs[i], abs=1e-5), (
                f"{name}: MathML gave {mathml_derivs[name]}, dydt() gave {real_derivs[i]}"
            )

    def test_matches_at_basal_steady_state(self):
        pathway = GlucoseInsulinPathway()
        self._assert_matches_dydt(pathway, t=0.0, y=pathway.y0)

    def test_matches_with_a_meal_in_progress(self):
        from biosim_ode.meal import MealInput
        pathway = GlucoseInsulinPathway(meals=[MealInput(carbs_g=50.0, time_min=0.0)])
        # X high enough that glucose is falling (dG < 0): the exporter does
        # not encode the dG-driven first-phase secretion term (documented in
        # sbml_export.py), which is zero only when dG <= 0. Before B7
        # (2026-11-16) X=0.02 sufficed because the double-counted resting
        # clearance made dG negative here.
        y = [pathway.params.Gb + 3.0, pathway.params.Ib + 40.0, 0.05,
             pathway.params.IGF1b, 0.2, 0.5]
        assert pathway.dydt(30.0, y)[0] <= 0.0
        self._assert_matches_dydt(pathway, t=30.0, y=y)

    def test_matches_before_meal_onset(self):
        from biosim_ode.meal import MealInput
        pathway = GlucoseInsulinPathway(meals=[MealInput(carbs_g=50.0, time_min=60.0)])
        self._assert_matches_dydt(pathway, t=10.0, y=pathway.y0)

    def test_matches_with_insulin_below_basal(self):
        # Exercises dX's I<Ib clamp (the fix -- dX previously used an
        # unclamped I-Ib, wrong whenever insulin dips below basal, e.g.
        # fasting/exercise). Deliberately keeps G exactly at Gb so both
        # dG=0 (no meal, no clearance imbalance) and delta_secretion=0
        # (current_frac==basal_frac at G==Gb) -- isolates the dX clamp
        # from the documented first-phase-secretion gap below.
        pathway = GlucoseInsulinPathway()
        y = [pathway.params.Gb, pathway.params.Ib - 20.0, 0.05,
             pathway.params.IGF1b, 0.15, 0.6]
        self._assert_matches_dydt(pathway, t=5.0, y=y)

    def test_matches_below_threshold_except_for_the_documented_first_phase_gap(self):
        # G below G_threshold exercises dG_above's zero-branch clamp, but
        # in a no-meal scenario this also necessarily makes dG>0 (the
        # clearance term restores G toward Gb), which fires
        # first_phase_secretion -- deliberately not encoded in the
        # exported SBML (see this module's docstring: it depends on the
        # instantaneous derivative of G, not G's value, which plain SBML
        # rate rules can't reference without duplicating the entire dG
        # expression). Confirms every OTHER species matches exactly, and
        # that I's discrepancy is EXACTLY the expected first-phase
        # contribution -- not some other, unexplained bug hiding behind
        # the documented one.
        pathway = GlucoseInsulinPathway()
        y = [pathway.params.Gb - 1.0, pathway.params.Ib - 20.0, 0.0,
             pathway.params.IGF1b, 0.15, 0.6]
        t = 5.0
        real_derivs = pathway.dydt(t, y)
        xml_str = SBMLExporter().export(pathway)
        context = self._param_context(xml_str)
        context.update({"G": y[0], "I": y[1], "X": y[2], "IGF1": y[3], "mTOR": y[4], "AMPK": y[5], "t": t})
        mathml_derivs = self._rate_rule_values(xml_str, context)

        indices = {"G": 0, "X": 2, "IGF1": 3, "mTOR": 4, "AMPK": 5}
        for name, idx in indices.items():
            assert mathml_derivs[name] == pytest.approx(real_derivs[idx], abs=1e-5)

        dG_real = real_derivs[0]
        expected_first_phase = (
            pathway.params.first_phase_gain * max(dG_real, 0.0) * pathway.params.secretion_scale
        )
        assert dG_real > 0.0, "test assumes dG>0 here so first-phase actually fires"
        assert real_derivs[1] - mathml_derivs["I"] == pytest.approx(expected_first_phase, rel=1e-3)

    def test_matches_with_multiple_meals(self):
        # Verifies multi-meal Ra(t) summation specifically (dG only) --
        # with a second meal still actively absorbing at t=150, dG is
        # positive here, which fires the documented first-phase-secretion
        # gap for I (see test_matches_below_threshold_..._first_phase_gap
        # above); asserting the full derivative vector would conflate the
        # two concerns, so this checks G directly instead.
        from biosim_ode.meal import MealInput
        pathway = GlucoseInsulinPathway(meals=[
            MealInput(carbs_g=30.0, time_min=0.0),
            MealInput(carbs_g=40.0, time_min=120.0),
        ])
        y = [pathway.params.Gb + 1.5, pathway.params.Ib + 15.0, 0.01,
             pathway.params.IGF1b, 0.18, 0.55]
        t = 150.0
        real_derivs = pathway.dydt(t, y)
        xml_str = SBMLExporter().export(pathway)
        context = self._param_context(xml_str)
        context.update({"G": y[0], "I": y[1], "X": y[2], "IGF1": y[3], "mTOR": y[4], "AMPK": y[5], "t": t})
        mathml_derivs = self._rate_rule_values(xml_str, context)
        assert mathml_derivs["G"] == pytest.approx(real_derivs[0], abs=1e-5)

    def _assert_glucose_matches(self, pathway, t: float):
        y = [pathway.params.Gb + 1.5, pathway.params.Ib + 15.0, 0.01,
             pathway.params.IGF1b, 0.18, 0.55]
        real_derivs = pathway.dydt(t, y)
        xml_str = SBMLExporter().export(pathway)
        context = self._param_context(xml_str)
        context.update({"G": y[0], "I": y[1], "X": y[2], "IGF1": y[3], "mTOR": y[4], "AMPK": y[5], "t": t})
        mathml_derivs = self._rate_rule_values(xml_str, context)
        assert mathml_derivs["G"] == pytest.approx(real_derivs[0], abs=1e-5)

    @pytest.mark.parametrize("t", [0.0, 15.0, 45.0, 200.0])
    def test_default_two_compartment_meal_matches(self, t):
        from biosim_ode import high_sugar_meal
        self._assert_glucose_matches(GlucoseInsulinPathway(meals=[high_sugar_meal()]), t)

    def test_gi_meal_and_kabs_meal_together_match(self):
        from biosim_api import MealEvent
        from biosim_ode.meal import MealInput
        meals = [
            MealInput(carbs_g=30.0, time_min=0.0),
            MealEvent(time_s=3600.0, carbs_g=40.0, glycaemic_index=55.0).to_meal_input(),
        ]
        for t in (30.0, 59.0, 90.0):
            self._assert_glucose_matches(GlucoseInsulinPathway(meals=meals), t)


# ---------------------------------------------------------------------------
# BatchRunner
# ---------------------------------------------------------------------------

class TestBatchRunner:
    def test_runs_all_scenarios(self):
        scenarios = [high_sugar_meal_scenario(), light_snack_scenario(), fasted_scenario(3600)]
        result = BatchRunner(max_workers=2).run_all(scenarios)
        assert result.n_succeeded == 3
        assert result.n_failed == 0

    def test_result_count_matches(self):
        scenarios = [high_sugar_meal_scenario(), fasted_scenario(3600)]
        result = BatchRunner().run_all(scenarios)
        assert len(result.outputs) == 2

    def test_get_by_name(self):
        scenarios = [high_sugar_meal_scenario(), fasted_scenario(3600)]
        result = BatchRunner().run_all(scenarios)
        out = result.get("high_sugar_meal")
        assert out is not None
        assert out.scenario_name == "high_sugar_meal"

    def test_get_unknown_returns_none(self):
        result = BatchRunner().run_all([fasted_scenario(3600)])
        assert result.get("nonexistent") is None

    def test_success_rate_all_pass(self):
        result = BatchRunner().run_all([fasted_scenario(3600)])
        assert result.success_rate == pytest.approx(1.0)

    def test_batch_result_scenarios_list(self):
        scenarios = [high_sugar_meal_scenario(), light_snack_scenario()]
        result = BatchRunner().run_all(scenarios)
        assert "high_sugar_meal" in result.scenarios
        assert "light_snack" in result.scenarios

    def test_failure_reason_is_recorded_not_swallowed(self, monkeypatch):
        # Regression: run_all()'s `except Exception as e: failed.append(name)`
        # discarded str(e) entirely -- a scenario name showing up in
        # failed_names gave no way to tell a real bug (affecting every
        # scenario) from an unrelated, expected one-off failure. Mirrors
        # the same fix already applied to PopulationRunner
        # (test_population.py::TestPopulationRunnerFailureHandling).
        #
        # Uses ThreadPoolExecutor explicitly, not the (now-default)
        # ProcessPoolExecutor: monkeypatch.setattr mutates SimulationRunner
        # in THIS process's memory only -- a spawned worker process
        # re-imports biosim_api fresh and never sees the patch, so this
        # failure-injection technique only works with a same-process
        # (thread-based) executor. Real production failures (a genuine
        # exception raised inside a worker process) are still correctly
        # captured via future.result() regardless of executor type --
        # see test_real_failure_across_process_boundary_is_captured below,
        # which proves that with the actual default executor.
        def _raise(self, scenario):
            raise ValueError("synthetic ODE solver failure for testing")

        monkeypatch.setattr(SimulationRunner, "run", _raise)
        scenarios = [high_sugar_meal_scenario(), fasted_scenario(3600)]
        result = BatchRunner(executor_cls=ThreadPoolExecutor).run_all(scenarios)

        assert result.n_failed == 2
        assert len(result.failure_reasons) == 2
        assert all(
            "ValueError" in reason and "synthetic ODE solver failure" in reason
            for reason in result.failure_reasons
        )

    def test_real_failure_across_process_boundary_is_captured(self):
        # Confirms a genuine exception raised INSIDE a worker process
        # (not monkeypatched -- a real, missing-attribute failure)
        # propagates back through future.result() and is recorded
        # correctly, using the real default executor
        # (ProcessPoolExecutor), not the thread-based one the test
        # above uses for monkeypatch compatibility. types.SimpleNamespace
        # (not a full ScenarioDefinition) has the `.name` attribute
        # run_all() itself needs for bookkeeping, but none of the other
        # attributes SimulationRunner.run() needs -- so it constructs
        # fine (no __post_init__ validation to dodge) and fails
        # genuinely, deterministically, once actually simulated inside
        # the worker. Also must be picklable across the process
        # boundary -- SimpleNamespace is a plain stdlib type, confirmed
        # directly before writing this test.
        import types
        bad_scenario = types.SimpleNamespace(name="genuinely_invalid")
        result = BatchRunner().run_all([high_sugar_meal_scenario(), bad_scenario])
        assert result.n_succeeded == 1
        assert result.n_failed == 1
        assert len(result.failure_reasons) == 1
        assert "AttributeError" in result.failure_reasons[0]

    def test_default_executor_is_process_pool_not_thread_pool(self):
        # Regression guard for the fix itself: ThreadPoolExecutor was
        # measured providing ZERO real parallelism for this CPU-bound
        # workload (0.66x sequential -- see batch.py's module
        # docstring) -- confirms the default wasn't silently reverted.
        assert BatchRunner().executor_cls is ProcessPoolExecutor

    def test_process_pool_executor_actually_runs_scenarios_correctly(self):
        # End-to-end confirmation that real multiprocessing produces the
        # SAME results as before (not just "doesn't crash") -- pickling
        # ScenarioDefinition/SimulationOutput across the process
        # boundary must not silently corrupt or drop data.
        scenarios = [high_sugar_meal_scenario(), light_snack_scenario(), fasted_scenario(3600)]
        result = BatchRunner(max_workers=2).run_all(scenarios)
        assert result.n_succeeded == 3
        out = result.get("high_sugar_meal")
        assert out is not None
        assert out.expert_metrics.glucose_peak_mM > 5.0

    def test_parameter_sweep_renames_each_scenario_uniquely(self):
        # parameter_sweep() had zero test coverage before this. Also
        # regression: the rename used to be gated by an always-False
        # `hasattr(s, '_swept')` check (dead code, removed) -- confirm
        # the rename it gated still happens for every value.
        import dataclasses

        def _modifier(scenario, value):
            return dataclasses.replace(scenario, stress_level=value)

        result = BatchRunner().parameter_sweep(
            base_scenario=fasted_scenario(3600),
            parameter_name="stress_level",
            values=[0.0, 0.5, 1.0],
            modifier_fn=_modifier,
        )
        assert result.n_succeeded == 3
        assert len(set(result.scenarios)) == 3
        for value in [0.0, 0.5, 1.0]:
            assert any(f"stress_level_{value:.3g}" in name for name in result.scenarios)


# ---------------------------------------------------------------------------
# CSVExporter
# ---------------------------------------------------------------------------

class TestCSVExporter:
    def test_time_series_is_string(self, meal_output):
        csv_str = CSVExporter().export_time_series(meal_output)
        assert isinstance(csv_str, str)
        assert len(csv_str) > 100

    def test_time_series_has_header(self, meal_output):
        csv_str = CSVExporter().export_time_series(meal_output)
        lines = csv_str.strip().split("\n")
        header = lines[0]
        assert "time_s" in header
        assert "time_min" in header

    def test_time_series_has_glucose_column(self, meal_output):
        csv_str = CSVExporter().export_time_series(meal_output)
        assert "G_" in csv_str

    def test_time_series_row_count(self, meal_output):
        csv_str = CSVExporter().export_time_series(meal_output)
        lines = [l for l in csv_str.strip().split("\n") if l]
        # header + data rows
        assert len(lines) > 100

    def test_time_series_species_subset(self, meal_output):
        csv_str = CSVExporter().export_time_series(meal_output, species=["G", "I"])
        assert "G_" in csv_str
        assert "I_" in csv_str
        assert "mTOR" not in csv_str

    def test_metrics_table_is_string(self, meal_output, fasted_output):
        csv_str = CSVExporter().export_metrics_table([meal_output, fasted_output])
        assert isinstance(csv_str, str)

    def test_metrics_table_has_header(self, meal_output):
        csv_str = CSVExporter().export_metrics_table([meal_output])
        first_line = csv_str.strip().split("\n")[0]
        assert "glucose_peak_mM" in first_line
        assert "scenario" in first_line

    def test_metrics_table_row_per_output(self, meal_output, fasted_output):
        csv_str = CSVExporter().export_metrics_table([meal_output, fasted_output])
        lines = [l for l in csv_str.strip().split("\n") if l]
        assert len(lines) == 3  # header + 2 data rows

    def test_metrics_table_parseable(self, meal_output):
        csv_str = CSVExporter().export_metrics_table([meal_output])
        reader = csv.DictReader(io.StringIO(csv_str))
        rows = list(reader)
        assert len(rows) == 1
        assert float(rows[0]["glucose_peak_mM"]) > 5.0

    def test_comparison_export(self, comparison):
        csv_str = CSVExporter().export_comparison(comparison)
        assert "G_diff_uM" in csv_str
        assert "I_diff_uM" in csv_str

    def test_curation_queue_export_is_string(self):
        from biosim_curation import CurationItemType, CurationPriority, CurationQueue

        queue = CurationQueue()
        queue.enqueue(
            item_id="test1", item_type=CurationItemType.PARAMETER,
            reason="test reason", subject_name="insulin receptor",
            parameter_kind="Km", value=45.0, unit="uM",
            source_description="PubMed:12345678", priority=CurationPriority.HIGH,
        )
        csv_str = CSVExporter().export_curation_queue(queue.pending())
        assert isinstance(csv_str, str)

    def test_curation_queue_export_has_header(self):
        from biosim_curation import CurationQueue

        csv_str = CSVExporter().export_curation_queue(CurationQueue().pending())
        header = csv_str.strip().split("\n")[0]
        assert "item_id" in header
        assert "priority" in header
        assert "subject_name" in header

    def test_curation_queue_export_row_per_item(self):
        from biosim_curation import CurationItemType, CurationPriority, CurationQueue

        queue = CurationQueue()
        queue.enqueue(
            item_id="a", item_type=CurationItemType.PARAMETER, reason="r1",
            subject_name="s1", parameter_kind="Km", value=1.0, unit="uM",
            source_description="src1",
        )
        queue.enqueue(
            item_id="b", item_type=CurationItemType.PARAMETER, reason="r2",
            subject_name="s2", parameter_kind="Vmax", value=2.0, unit="uM/s",
            source_description="src2",
        )
        csv_str = CSVExporter().export_curation_queue(queue.pending())
        lines = [l for l in csv_str.strip().split("\n") if l]
        assert len(lines) == 3  # header + 2 items

    def test_curation_queue_export_parseable_and_values_correct(self):
        from biosim_curation import CurationItemType, CurationPriority, CurationQueue

        queue = CurationQueue()
        queue.enqueue(
            item_id="parseable_test", item_type=CurationItemType.PARAMETER,
            reason="unit conflict", subject_name="D-Glucose",
            parameter_kind="abundance", value=5.5, unit="mM",
            source_description="PubMed:99999999", priority=CurationPriority.CRITICAL,
        )
        csv_str = CSVExporter().export_curation_queue(queue.pending())
        reader = csv.DictReader(io.StringIO(csv_str))
        rows = list(reader)
        assert len(rows) == 1
        assert rows[0]["item_id"] == "parseable_test"
        assert rows[0]["priority"] == "CRITICAL"
        assert rows[0]["subject_name"] == "D-Glucose"
        assert float(rows[0]["value"]) == 5.5
        assert rows[0]["unit"] == "mM"
        assert rows[0]["assigned_to"] == ""

    def test_curation_queue_export_empty_queue(self):
        from biosim_curation import CurationQueue

        csv_str = CSVExporter().export_curation_queue(CurationQueue().pending())
        lines = [l for l in csv_str.strip().split("\n") if l]
        assert len(lines) == 1  # header only, no crash on empty input

    def test_curation_queue_export_handles_embedded_commas_and_quotes(self):
        # Mandatory code-reviewer finding (MEDIUM): real reason/
        # source_description strings this project actually generates
        # (see enqueue_conflict()'s own f-string, which embeds
        # parenthesized source values) commonly contain commas -- a
        # naive manual-string-join "simplification" of this method
        # would silently corrupt the CSV for exactly this real-world
        # shape. Round-trips a value with an embedded comma AND a
        # double-quote through csv.writer/csv.DictReader and confirms
        # it comes back byte-identical.
        from biosim_curation import CurationItemType, CurationPriority, CurationQueue

        tricky_reason = 'Conflict: 5.2x diff between value 1.0 (PubMed:123, Smith) and 2.0 (Jones "et al")'
        tricky_source = 'PubMed:123, Smith 2020 vs PubMed:456, "Jones" 2021'
        queue = CurationQueue()
        queue.enqueue(
            item_id="tricky_test", item_type=CurationItemType.PARAMETER,
            reason=tricky_reason, subject_name="insulin receptor",
            parameter_kind="Km", value=10.0, unit="uM",
            source_description=tricky_source, priority=CurationPriority.MEDIUM,
        )
        csv_str = CSVExporter().export_curation_queue(queue.pending())
        rows = list(csv.DictReader(io.StringIO(csv_str)))
        assert len(rows) == 1
        assert rows[0]["reason"] == tricky_reason
        assert rows[0]["source_description"] == tricky_source

    def test_curation_queue_export_sanitizes_formula_injection_in_free_text(self):
        # Housekeeping-backlog fix (was disclosed-not-fixed at MASTER_PLAN
        # 5.2d's own mandatory code-reviewer pass; closed now that 5.2c
        # wires genuinely-less-trusted extraction-derived text into this
        # same export path). A subject_name/reason/source_description
        # starting with =, +, -, @, tab, or CR is the classic Excel/
        # LibreOffice/Google Sheets "CSV formula injection" vector
        # (OWASP CWE-1236) -- confirms every free-text column gets a
        # defensive leading single-quote, and that the quote survives a
        # real csv round-trip rather than being stripped or corrupting
        # the row.
        from biosim_curation import CurationItemType, CurationPriority, CurationQueue

        queue = CurationQueue()
        queue.enqueue(
            item_id="=cmd|'/c calc'!A1",
            item_type=CurationItemType.PARAMETER,
            reason="=HYPERLINK(\"http://evil.example/\",\"click me\")",
            subject_name="+1+1",
            parameter_kind="-1-1",
            value=1.0,
            unit="@SUM(A1:A9)",
            source_description="@import(bad)",
            priority=CurationPriority.MEDIUM,
        )
        csv_str = CSVExporter().export_curation_queue(queue.pending())
        rows = list(csv.DictReader(io.StringIO(csv_str)))
        assert len(rows) == 1
        row = rows[0]
        for column in (
            "item_id", "subject_name", "parameter_kind",
            "unit", "reason", "source_description",
        ):
            assert row[column].startswith("'"), (
                f"{column}={row[column]!r} should be quote-prefixed"
            )
            # The original value (minus the defensive prefix) must
            # still be fully recoverable -- this is a display-safety
            # measure, not silent data loss.
            assert row[column][1:] in (
                "=cmd|'/c calc'!A1",
                "+1+1",
                "-1-1",
                "@SUM(A1:A9)",
                "=HYPERLINK(\"http://evil.example/\",\"click me\")",
                "@import(bad)",
            )

    def test_curation_queue_export_does_not_sanitize_ordinary_text(self):
        # A no-false-positives guard alongside the injection test above --
        # ordinary text (including a real hyphenated/at-sign-containing
        # value that does NOT start with a trigger character) must pass
        # through completely unchanged.
        from biosim_curation import CurationItemType, CurationPriority, CurationQueue

        queue = CurationQueue()
        queue.enqueue(
            item_id="ordinary-id",
            item_type=CurationItemType.PARAMETER,
            reason="insulin receptor Kd, 3-5x outlier",
            subject_name="insulin receptor",
            parameter_kind="Kd",
            value=1.0,
            unit="uM",
            source_description="curator@example.com noted this",
            priority=CurationPriority.MEDIUM,
        )
        csv_str = CSVExporter().export_curation_queue(queue.pending())
        rows = list(csv.DictReader(io.StringIO(csv_str)))
        assert rows[0]["item_id"] == "ordinary-id"
        assert rows[0]["reason"] == "insulin receptor Kd, 3-5x outlier"
        assert rows[0]["source_description"] == "curator@example.com noted this"

    def test_curation_queue_export_sanitizes_assigned_to(self):
        # Found by mandatory security review: the injection test above
        # covers 6 of the 7 sanitized fields but not `assigned_to`,
        # which is set via CurationQueue.assign(item_id, curator_id) --
        # a different call path than enqueue()'s own keyword arguments.
        from biosim_curation import CurationItemType, CurationPriority, CurationQueue

        queue = CurationQueue()
        queue.enqueue(
            item_id="assign_test", item_type=CurationItemType.PARAMETER,
            reason="r", subject_name="s", parameter_kind="Km",
            value=1.0, unit="uM", source_description="src",
            priority=CurationPriority.MEDIUM,
        )
        queue.assign("assign_test", curator_id="=cmd|'/c calc'!A1")
        csv_str = CSVExporter().export_curation_queue(queue.pending())
        rows = list(csv.DictReader(io.StringIO(csv_str)))
        assert rows[0]["assigned_to"] == "'=cmd|'/c calc'!A1"

    def test_metrics_table_sanitizes_formula_injection_in_scenario_name(self):
        # Same fix, the other real free-text field in this exporter --
        # SimulationOutput.scenario_name flows in from whatever created
        # the scenario (a REST client's own ScenarioIn.name, in
        # production), a less-trusted source than this file's own
        # internally-generated strings.
        from biosim_api.output import ExpertMetrics, SimulationOutput
        from biosim_ode import ODESystem

        import dataclasses

        output = SimulationRunner().run(fasted_scenario(1))
        malicious_output = dataclasses.replace(output, scenario_name="=cmd|'/c calc'!A1")
        csv_str = CSVExporter().export_metrics_table([malicious_output])
        reader = csv.DictReader(io.StringIO(csv_str))
        rows = list(reader)
        assert rows[0]["scenario"] == "'=cmd|'/c calc'!A1"

    def test_failed_simulation_export(self):
        from biosim_api.output import SimulationOutput, ExpertMetrics, LaySummary
        import numpy as np
        # Create a minimal failed output
        from biosim_ode import ODESystem
        # Just check it doesn't crash with a failed sim
        output = SimulationRunner().run(fasted_scenario(1))
        csv_str = CSVExporter().export_time_series(output)
        assert isinstance(csv_str, str)
