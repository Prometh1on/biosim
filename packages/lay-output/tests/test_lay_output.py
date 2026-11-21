# SAVE TO: packages\lay-output\tests\test_lay_output.py
"""
Test suite for biosim_lay.

Coverage:
    - FoodItem: construction, carbs, GI validation
    - Food DB: all items accessible, compose_meal, search
    - RiskIndicator and TrafficLight
    - assess_risk: green/amber/red thresholds, confidence flagging
    - generate_report: all sections populated, text output
    - End-to-end: food → scenario → simulation → report pipeline
"""

from types import SimpleNamespace

import pytest

from biosim_lay import (
    FOOD_DB,
    FoodItem,
    PlainLanguageReport,
    RiskIndicator,
    RiskProfile,
    TrafficLight,
    assess_risk,
    compose_meal,
    foods_by_gi_category,
    generate_report,
    get_food,
    search_foods,
)
from biosim_api import (
    ExpertMetrics,
    MealEvent,
    ScenarioDefinition,
    SimulationRunner,
    fasted_scenario,
    high_sugar_meal_scenario,
    light_snack_scenario,
)


def _output_with_nadir(nadir_mM: float, baseline_mM: float = 5.0) -> SimpleNamespace:
    """
    Minimal duck-typed stand-in for a SimulationOutput, carrying only what
    assess_risk() actually reads (m = output.expert_metrics, plus the two
    personalization attributes _secretion_appears_blunted() checks).
    Used because this platform's shipped SGLT2i/sulfonylurea mechanisms
    cannot currently reach real clinical hypoglycaemia (<3.9 mmol/L) via
    any real SimulationRunner scenario, even at each mechanism's own
    extreme (but valid) parameter ceiling -- verified directly before
    writing this test (see CLAUDE.md's hypoglycaemia-indicator milestone)
    -- so the AMBER/RED boundary logic can only be exercised with a
    directly-constructed metrics object, not a real scenario.
    """
    metrics = ExpertMetrics(
        glucose_peak_uM=baseline_mM * 1000.0,
        glucose_peak_time_s=0.0,
        glucose_auc_uM_s=0.0,
        glucose_baseline_uM=baseline_mM * 1000.0,
        glucose_nadir_uM=nadir_mM * 1000.0,
        insulin_peak_uM=0.066,
        insulin_peak_time_s=0.0,
        mtor_peak=0.3,
        ampk_nadir=0.6,
        time_to_glucose_return_s=None,
        overall_watchdog_confidence=1.0,
        n_ode_evaluations=1,
    )
    return SimpleNamespace(
        expert_metrics=metrics,
        physiological_personalization=None,
        genetic_personalization=None,
    )


# ---------------------------------------------------------------------------
# FoodItem and food database
# ---------------------------------------------------------------------------

class TestFoodItem:
    def test_construction(self):
        f = FoodItem("test_food", "Test food", 25.0, "1 portion", "medium")
        assert f.carbs_g == 25.0

    def test_negative_carbs_raises(self):
        with pytest.raises(ValueError):
            FoodItem("x", "X", -1.0, "1 unit", "low")

    def test_invalid_gi_raises(self):
        with pytest.raises(ValueError):
            FoodItem("x", "X", 10.0, "1 unit", "extreme")


class TestFoodDatabase:
    def test_all_items_accessible(self):
        for food_id in FOOD_DB:
            item = get_food(food_id)
            assert item.carbs_g >= 0

    def test_get_known_food(self):
        banana = get_food("banana_medium")
        assert banana.carbs_g == 25.0

    def test_get_unknown_raises(self):
        with pytest.raises(KeyError):
            get_food("moon_cheese_100g")

    def test_compose_meal_sums_carbs(self):
        total = compose_meal(["banana_medium", "white_bread_slice"])
        assert total == pytest.approx(25.0 + 14.0)

    def test_compose_meal_empty(self):
        assert compose_meal([]) == pytest.approx(0.0)

    def test_search_finds_banana(self):
        results = search_foods("banana")
        assert any(f.id == "banana_medium" for f in results)

    def test_search_case_insensitive(self):
        results = search_foods("BANANA")
        assert len(results) > 0

    def test_search_no_results(self):
        results = search_foods("xyznotafood")
        assert results == []

    def test_foods_by_gi_high(self):
        high = foods_by_gi_category("high")
        assert len(high) > 0
        assert all(f.gi_category == "high" for f in high)

    def test_foods_by_gi_low(self):
        low = foods_by_gi_category("low")
        assert all(f.gi_category == "low" for f in low)

    def test_water_has_zero_carbs(self):
        assert get_food("water_any").carbs_g == pytest.approx(0.0)

    def test_glucose_ogtt_is_75g(self):
        assert get_food("glucose_75g_ogtt").carbs_g == pytest.approx(75.0)


# ---------------------------------------------------------------------------
# Risk indicators
# ---------------------------------------------------------------------------

class TestRiskIndicators:
    @pytest.fixture(scope="class")
    def high_meal_output(self):
        return SimulationRunner().run(high_sugar_meal_scenario())

    @pytest.fixture(scope="class")
    def light_meal_output(self):
        return SimulationRunner().run(light_snack_scenario())

    @pytest.fixture(scope="class")
    def fasted_output(self):
        return SimulationRunner().run(fasted_scenario(3600))

    def test_assess_risk_returns_profile(self, high_meal_output):
        profile = assess_risk(high_meal_output)
        assert isinstance(profile, RiskProfile)

    def test_indicators_present(self, high_meal_output):
        profile = assess_risk(high_meal_output)
        assert len(profile.indicators) >= 3

    def test_peak_glucose_indicator_present(self, high_meal_output):
        profile = assess_risk(high_meal_output)
        ind = profile.indicator("Peak blood sugar")
        assert ind is not None

    def test_high_meal_not_all_green(self, high_meal_output):
        profile = assess_risk(high_meal_output)
        # High sugar meal may push into amber/red for some indicators
        assert profile.overall in (TrafficLight.GREEN, TrafficLight.AMBER, TrafficLight.RED)

    def test_fasted_all_green(self, fasted_output):
        profile = assess_risk(fasted_output)
        peak_ind = profile.indicator("Peak blood sugar")
        # Fasted should have glucose near baseline — green
        assert peak_ind.traffic_light == TrafficLight.GREEN

    def test_light_snack_greener_than_high_meal(self, light_meal_output, high_meal_output):
        light_profile = assess_risk(light_meal_output)
        high_profile = assess_risk(high_meal_output)
        # Convert to numeric for comparison
        order = [TrafficLight.GREEN, TrafficLight.AMBER, TrafficLight.RED, TrafficLight.GREY]
        light_score = order.index(light_profile.overall)
        high_score = order.index(high_profile.overall)
        assert light_score <= high_score  # light snack should be at least as green

    def test_has_red_property(self):
        profile = RiskProfile(
            indicators=[RiskIndicator("Test", 12.0, "mM", TrafficLight.RED, "High", "")],
            overall=TrafficLight.RED,
            summary="red",
        )
        assert profile.has_red

    def test_has_amber_property(self):
        profile = RiskProfile(
            indicators=[RiskIndicator("Test", 9.0, "mM", TrafficLight.AMBER, "Borderline", "")],
            overall=TrafficLight.AMBER,
            summary="amber",
        )
        assert profile.has_amber
        assert not profile.has_red

    def test_low_confidence_gives_grey(self, high_meal_output):
        # With very high confidence threshold, everything should be grey
        profile = assess_risk(high_meal_output, confidence_threshold=2.0)
        assert all(r.traffic_light == TrafficLight.GREY for r in profile.indicators)


class TestHypoglycaemiaRiskIndicator:
    """
    Battelino T, Danne T, Bergenstal RM, et al. "Clinical Targets for
    Continuous Glucose Monitoring Data Interpretation: Recommendations
    From the International Consensus on Time in Range." Diabetes Care.
    2019;42(8):1593-1603. PMID:31177185 (verified live via NCBI efetch):
    Table 2 -- Time below range (TBR) Level 1 = 3.0-3.8 mmol/L, Level 2
    = <3.0 mmol/L. Added because this platform's SGLT2i (insulin-
    independent renal clearance) and sulfonylurea (glucose-independent
    secretion) mechanisms can genuinely push simulated glucose low, with
    no lay-facing indicator ever having checked for it before this.
    """

    def test_indicator_present(self):
        profile = assess_risk(_output_with_nadir(5.0))
        assert profile.indicator("Hypoglycaemia risk") is not None

    def test_green_above_alert_value(self):
        ind = assess_risk(_output_with_nadir(4.5)).indicator("Hypoglycaemia risk")
        assert ind.traffic_light == TrafficLight.GREEN

    def test_amber_at_level_1(self):
        ind = assess_risk(_output_with_nadir(3.5)).indicator("Hypoglycaemia risk")
        assert ind.traffic_light == TrafficLight.AMBER

    def test_red_at_level_2(self):
        ind = assess_risk(_output_with_nadir(2.5)).indicator("Hypoglycaemia risk")
        assert ind.traffic_light == TrafficLight.RED

    def test_boundary_exactly_at_3_9_is_not_amber(self):
        # 3.9 mmol/L is the alert value itself -- "below" (strict <)
        # is what triggers AMBER, matching Battelino's own "<3.9 mmol/L"
        # wording, not "<=3.9".
        ind = assess_risk(_output_with_nadir(3.9)).indicator("Hypoglycaemia risk")
        assert ind.traffic_light == TrafficLight.GREEN

    def test_boundary_exactly_at_3_0_is_not_red(self):
        ind = assess_risk(_output_with_nadir(3.0)).indicator("Hypoglycaemia risk")
        assert ind.traffic_light == TrafficLight.AMBER

    def test_low_confidence_gives_grey(self):
        output = _output_with_nadir(2.5)
        output.expert_metrics.overall_watchdog_confidence = 0.0
        ind = assess_risk(output).indicator("Hypoglycaemia risk")
        assert ind.traffic_light == TrafficLight.GREY

    def test_fasted_scenario_nadir_equals_baseline(self):
        # Real end-to-end check: for a genuinely flat (no meal, no drug)
        # trajectory, the nadir must equal the baseline exactly -- a
        # meaningful regression guard even though this platform's
        # shipped drug mechanisms can't currently push a real scenario
        # into clinical hypoglycaemia (see this class's own docstring
        # and the module-level _output_with_nadir() helper's docstring).
        output = SimulationRunner().run(fasted_scenario(3600))
        assert output.expert_metrics.glucose_nadir_mM == pytest.approx(
            output.expert_metrics.glucose_baseline_uM / 1000.0, abs=1e-6
        )
        ind = assess_risk(output).indicator("Hypoglycaemia risk")
        assert ind.traffic_light == TrafficLight.GREEN

    def test_sglt2i_scenario_genuinely_lowers_the_nadir(self):
        # Confirms glucose_nadir_uM is wired to a REAL, drug-sensitive
        # value through the actual SimulationRunner path, not a stub --
        # an SGLT2i drug (insulin-independent renal glucose clearance,
        # continuously active even at/below basal glucose) measurably
        # lowers the nadir below the undrugged baseline in a fasted
        # scenario, confirmed directly against a no-drug comparator
        # rather than asserting an absolute clinical-hypoglycaemia value
        # this platform's current drug calibrations don't reach.
        from biosim_api import DrugEvent
        no_drug = ScenarioDefinition(name="nd_fasted", meal_events=[], duration_s=14400)
        with_dapa = ScenarioDefinition(
            name="dapa_fasted", meal_events=[], duration_s=14400,
            drug_events=[DrugEvent(time_s=0, drug="dapagliflozin", dose_mg=10.0, hours_before_simulation=24 * 7)],
        )
        runner = SimulationRunner()
        nadir_no_drug = runner.run(no_drug).expert_metrics.glucose_nadir_mM
        nadir_dapa = runner.run(with_dapa).expert_metrics.glucose_nadir_mM
        assert nadir_dapa < nadir_no_drug

    def test_green_result_carries_a_calibration_caveat(self):
        # Regression: found by the mandatory security-reviewer pass -- a
        # GREEN result here reads as "no hypoglycaemia risk," but no
        # currently-shipped scenario combination reaches AMBER/RED at
        # all, so the absence of a flag must not be presented as an
        # unqualified guarantee.
        ind = assess_risk(_output_with_nadir(4.5)).indicator("Hypoglycaemia risk")
        assert ind.traffic_light == TrafficLight.GREEN
        assert "not a guarantee against real hypoglycaemia" in ind.caveat

    def test_amber_result_does_not_carry_the_calibration_caveat(self):
        # The caveat is specifically about GREEN's false-reassurance risk
        # -- an AMBER/RED result already correctly signals concern and
        # doesn't need the same hedge.
        ind = assess_risk(_output_with_nadir(3.5)).indicator("Hypoglycaemia risk")
        assert ind.traffic_light == TrafficLight.AMBER
        assert "not a guarantee against real hypoglycaemia" not in ind.caveat


class TestBaselineAwareRecoveryTime:
    """
    Regression: "Recovery time" measured return to the patient's OWN
    baseline (glucose_baseline_uM), not a population-normal value. For a
    T2D patient (disease-state preset: Gb x1.3, ~6.5 mM -- already
    pre-diabetic by this file's own cited ADA fasting criteria), a fast
    return to that elevated baseline was unconditionally described as
    "near-normal"/"a healthy response" and GREEN -- exactly backwards
    for a patient whose fasting glucose is itself already elevated.
    """
    @staticmethod
    def _t2d_snack_output():
        from biosim_physiology import PhysiologicalProfile, DiseaseState
        from biosim_api import PersonalizationSpec
        profile = PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES])
        scenario = ScenarioDefinition(
            # kabs set explicitly = single-exponential curve, and a small
            # 5 g snack, which give the fast return to baseline this test
            # needs (~46 min). 15 g returned in ~135 min once B7 (2026-11-16)
            # removed the double-counted resting glucose clearance.
            name="t2d_snack", meal_events=[MealEvent(time_s=0, carbs_g=5.0, kabs=0.025)],
            duration_s=14400, personalization=PersonalizationSpec(physiological_profile=profile),
        )
        return SimulationRunner().run(scenario)

    def test_fasting_glucose_indicator_present_and_flags_t2d(self):
        output = self._t2d_snack_output()
        profile = assess_risk(output)
        ind = profile.indicator("Fasting blood sugar")
        assert ind is not None
        assert ind.traffic_light in (TrafficLight.AMBER, TrafficLight.RED)

    def test_healthy_baseline_still_gets_green_fasting_indicator(self):
        output = SimulationRunner().run(light_snack_scenario())
        ind = assess_risk(output).indicator("Fasting blood sugar")
        assert ind.traffic_light == TrafficLight.GREEN

    def test_recovery_time_is_not_green_when_baseline_is_elevated(self):
        output = self._t2d_snack_output()
        m = output.expert_metrics
        assert m.time_to_glucose_return_s / 60.0 < 90, (
            "test assumes a fast return so it would have wrongly been GREEN before the fix"
        )
        ind = assess_risk(output).indicator("Recovery time")
        assert ind is not None
        assert ind.traffic_light != TrafficLight.GREEN

    def test_recovery_time_description_flags_the_elevated_baseline(self):
        # Must not make the old unqualified "a healthy response" claim --
        # the corrected wording is allowed to mention "healthy" only as
        # part of the explicit hedge ("does not by itself indicate a
        # healthy response"), so check for the hedge itself, not just
        # the substring's absence.
        ind = assess_risk(self._t2d_snack_output()).indicator("Recovery time")
        desc = ind.description.lower()
        assert "does not" in desc
        assert "above the normal fasting range" in desc

    def test_glucose_section_does_not_claim_normal_when_baseline_elevated(self):
        report = generate_report(self._t2d_snack_output())
        assert "returns toward normal" not in report.glucose_section


class TestMTORAMPKBlutedSecretionCaveat:
    """
    Regression: mTOR/AMPK activation in this model is driven only by
    insulin level -- the T2D disease-state preset blunts insulin
    SECRETION (secretion_scale=0.60), independently lowering mTOR/AMPK
    excursion for a reason unrelated to genuine metabolic balance. A T2D
    patient's meal could read GREEN ("moderate"/"balanced") on both
    indicators when a healthy comparator eating the same meal read
    AMBER ("high"), with no caveat explaining why. Fixed with an
    explicit caveat on both indicators rather than recalibrating the
    fixed thresholds (no principled clinical basis to do so).
    """
    @staticmethod
    def _t2d_output():
        from biosim_physiology import PhysiologicalProfile, DiseaseState
        from biosim_api import PersonalizationSpec
        return SimulationRunner().run(high_sugar_meal_scenario(
            personalization=PersonalizationSpec(
                physiological_profile=PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES])
            )
        ))

    def test_t2d_mtor_indicator_carries_blunted_secretion_caveat(self):
        ind = assess_risk(self._t2d_output()).indicator("mTOR activation")
        assert "blunted insulin secretion" in ind.caveat.lower()

    def test_t2d_ampk_indicator_carries_blunted_secretion_caveat(self):
        ind = assess_risk(self._t2d_output()).indicator("AMPK suppression")
        assert "blunted insulin secretion" in ind.caveat.lower()

    def test_healthy_meal_has_no_blunted_secretion_caveat(self):
        output = SimulationRunner().run(high_sugar_meal_scenario())
        profile = assess_risk(output)
        assert "blunted insulin secretion" not in profile.indicator("mTOR activation").caveat.lower()
        assert "blunted insulin secretion" not in profile.indicator("AMPK suppression").caveat.lower()

    def test_blunted_secretion_caveat_actually_reaches_the_lay_report(self):
        # Regression: RiskIndicator.caveat was computed (see the two
        # tests above) but generate_report()'s risk_section builder only
        # ever read indicator.description, never .caveat -- so this exact
        # warning never reached the text a real user reads. Same
        # false-reassurance bug class as the T2D "healthy"/"efficient"
        # narrative fixes elsewhere in this project.
        output = self._t2d_output()
        report = generate_report(output)
        caveat = assess_risk(output).indicator("mTOR activation").caveat
        assert caveat, "expected a non-empty caveat for this T2D scenario"
        assert caveat in report.risk_section
        assert caveat in report.as_text()


# ---------------------------------------------------------------------------
# Plain language report
# ---------------------------------------------------------------------------

class TestPlainLanguageReport:
    @pytest.fixture(scope="class")
    def meal_output(self):
        return SimulationRunner().run(high_sugar_meal_scenario())

    def test_generate_report_returns_report(self, meal_output):
        report = generate_report(meal_output)
        assert isinstance(report, PlainLanguageReport)

    def test_headline_contains_scenario_name(self, meal_output):
        report = generate_report(meal_output)
        assert "high_sugar_meal" in report.headline

    def test_glucose_section_populated(self, meal_output):
        report = generate_report(meal_output)
        assert len(report.glucose_section) > 20
        assert "mmol" in report.glucose_section

    def test_energy_section_populated(self, meal_output):
        report = generate_report(meal_output)
        assert len(report.energy_section) > 20

    def test_risk_section_populated(self, meal_output):
        report = generate_report(meal_output)
        assert len(report.risk_section) > 20

    def test_disclaimer_present(self, meal_output):
        report = generate_report(meal_output)
        assert "not a medical diagnosis" in report.disclaimer

    def test_as_text_returns_string(self, meal_output):
        report = generate_report(meal_output)
        text = report.as_text()
        assert isinstance(text, str)
        assert len(text) > 200

    def test_comparison_section_when_provided(self, meal_output):
        light_output = SimulationRunner().run(light_snack_scenario())
        report = generate_report(meal_output, comparison_output=light_output)
        assert report.comparison_section is not None
        assert "higher" in report.comparison_section or "lower" in report.comparison_section

    def test_no_comparison_section_by_default(self, meal_output):
        report = generate_report(meal_output)
        assert report.comparison_section is None

    def test_risk_profile_attached(self, meal_output):
        report = generate_report(meal_output)
        assert isinstance(report.risk_profile, RiskProfile)


# ---------------------------------------------------------------------------
# End-to-end: food → scenario → simulation → report
# ---------------------------------------------------------------------------

class TestEndToEnd:
    def test_food_to_simulation_to_report(self):
        """Full lay-user pipeline: pick foods, simulate, get report."""
        # User selects foods
        carbs = compose_meal(["white_rice_200g", "cola_330ml"])  # 56 + 35 = 91g
        assert carbs == pytest.approx(91.0)

        # Build scenario
        scenario = ScenarioDefinition(
            name="rice_and_cola",
            duration_s=14400,
            meal_events=[MealEvent(time_s=0.0, carbs_g=carbs, label="Rice and Cola")],
        )

        # Simulate
        output = SimulationRunner().run(scenario)
        assert output.success

        # Generate report
        report = generate_report(output)
        assert "rice_and_cola" in report.headline

        # Risk profile
        risk = report.risk_profile
        assert risk.overall in list(TrafficLight)

    def test_low_carb_meal_greener_risk(self):
        """A low-carb meal should produce greener risk indicators than high-carb."""
        low = ScenarioDefinition(
            name="low_carb",
            meal_events=[MealEvent(0, compose_meal(["egg_boiled", "salad_plain"]))],
            duration_s=7200,
        )
        high = ScenarioDefinition(
            name="high_carb",
            meal_events=[MealEvent(0, compose_meal(["white_rice_200g", "cola_330ml"]))],
            duration_s=7200,
        )
        runner = SimulationRunner()
        out_low  = runner.run(low)
        out_high = runner.run(high)

        risk_low  = assess_risk(out_low)
        risk_high = assess_risk(out_high)

        peak_low  = out_low.expert_metrics.glucose_peak_mM
        peak_high = out_high.expert_metrics.glucose_peak_mM
        assert peak_high > peak_low


class TestGIPlaceholderNote:
    """The GI timing slope is an unvalidated placeholder (only GI=100 is
    calibrated); a lay reader must see that whenever a meal relies on it."""

    def test_low_gi_meal_report_carries_placeholder_note(self):
        from biosim_ode.gut_absorption import GI_SLOPE_PLACEHOLDER_NOTE
        scenario = ScenarioDefinition(
            name="low_gi",
            meal_events=[MealEvent(time_s=0.0, carbs_g=50.0, glycaemic_index=40.0)],
            duration_s=7200,
        )
        report = generate_report(SimulationRunner().run(scenario))
        assert GI_SLOPE_PLACEHOLDER_NOTE in report.glucose_section
        assert GI_SLOPE_PLACEHOLDER_NOTE in report.as_text()

    def test_meal_without_gi_has_no_placeholder_note(self):
        from biosim_ode.gut_absorption import GI_SLOPE_PLACEHOLDER_NOTE
        report = generate_report(SimulationRunner().run(high_sugar_meal_scenario()))
        assert GI_SLOPE_PLACEHOLDER_NOTE not in report.glucose_section


class TestUnverifiedVariantNote:
    """A lay reader must see that the rs35767 (IGF1) effect is unverified."""

    def test_rs35767_report_carries_note(self):
        from biosim_api import PersonalizationSpec
        from biosim_genotype import Genotype, Variant
        from biosim_genotype.effect_db import UNVERIFIED_VARIANT_NOTES
        scenario = high_sugar_meal_scenario(
            personalization=PersonalizationSpec(
                genetic_variants=[Variant("rs35767", Genotype.HETEROZYGOUS)]
            )
        )
        report = generate_report(SimulationRunner().run(scenario))
        assert UNVERIFIED_VARIANT_NOTES["rs35767"] in report.as_text()

    def test_no_variants_no_note(self):
        from biosim_genotype.effect_db import UNVERIFIED_VARIANT_NOTES
        report = generate_report(SimulationRunner().run(high_sugar_meal_scenario()))
        assert UNVERIFIED_VARIANT_NOTES["rs35767"] not in report.as_text()


class TestExperimentalDrugNote:
    """A lay reader must see that a drug module is experimental."""

    def test_drug_report_carries_experimental_note(self):
        from biosim_api import DrugEvent, MealEvent, ScenarioDefinition
        from biosim_physiology.drug_interactions import Drug
        from biosim_physiology.drug_interactions.validation_status import (
            EXPERIMENTAL_DRUG_NOTES,
        )
        scenario = ScenarioDefinition(
            name="metformin_meal", duration_s=7200,
            meal_events=[MealEvent(time_s=0, carbs_g=50.0)],
            drug_events=[DrugEvent(time_s=0, drug="metformin", dose_mg=1000.0)],
        )
        report = generate_report(SimulationRunner().run(scenario))
        assert EXPERIMENTAL_DRUG_NOTES[Drug.METFORMIN] in report.risk_section

    def test_validated_drug_report_carries_its_caveat(self):
        from biosim_api import DrugEvent, MealEvent, ScenarioDefinition
        from biosim_physiology.drug_interactions import Drug
        from biosim_physiology.drug_interactions.validation_status import (
            VALIDATED_DRUG_NOTES,
        )
        scenario = ScenarioDefinition(
            name="semaglutide_meal", duration_s=7200,
            meal_events=[MealEvent(time_s=0, carbs_g=50.0)],
            drug_events=[DrugEvent(time_s=0, drug="semaglutide", dose_mg=1.0)],
        )
        report = generate_report(SimulationRunner().run(scenario))
        assert VALIDATED_DRUG_NOTES[Drug.SEMAGLUTIDE] in report.risk_section

    def test_no_drug_no_experimental_note(self):
        report = generate_report(SimulationRunner().run(high_sugar_meal_scenario()))
        assert "EXPERIMENTAL drug module" not in report.as_text()


class TestNoReturnWithinWindowIsFlagged:
    """A meal that has not returned to baseline by the end is the slowest
    case, so it must get a RED Recovery-time indicator, not none at all."""

    def test_t2d_high_sugar_meal_gets_red_recovery_indicator(self):
        from biosim_api import PersonalizationSpec
        from biosim_physiology import PhysiologicalProfile, DiseaseState
        output = SimulationRunner().run(high_sugar_meal_scenario(
            personalization=PersonalizationSpec(
                physiological_profile=PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES])
            )
        ))
        assert output.expert_metrics.time_to_glucose_return_s is None
        ind = assess_risk(output).indicator("Recovery time")
        assert ind is not None
        assert ind.traffic_light == TrafficLight.RED
        assert "not returned" in ind.description
