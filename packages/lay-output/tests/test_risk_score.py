# SAVE TO: packages\lay-output\tests\test_risk_score.py
"""Tests for the metabolic risk score module."""
from types import SimpleNamespace

import pytest
from biosim_lay.risk_score import (
    RiskBand, RiskIndicatorScore, MetabolicRiskProfile, compute_risk_profile,
)
from biosim_api import (
    DrugEvent,
    ExpertMetrics,
    ScenarioDefinition,
    SimulationRunner,
    high_sugar_meal_scenario,
    fasted_scenario,
    light_snack_scenario,
)


def _output_with_nadir(nadir_mM: float, baseline_mM: float = 5.0) -> SimpleNamespace:
    """
    Minimal duck-typed stand-in for a SimulationOutput -- see the
    identical helper's docstring in test_lay_output.py for why a real
    SimulationRunner scenario can't reach the AMBER/RED boundary this
    tests (this platform's shipped SGLT2i/sulfonylurea mechanisms don't
    currently push simulated glucose below the clinical hypoglycaemia
    alert value, even at each mechanism's own extreme valid parameter).
    compute_risk_profile() gates on output.simulation.success now (see the
    "Simulation status"/RiskBand.UNKNOWN regression fix) -- success=True
    here so this double isn't itself treated as a failed simulation; the
    "time above 7.8" block is separately skipped because this SimpleNamespace
    has no `species` attribute at all (that block's own `hasattr` guard).
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
        simulation=SimpleNamespace(success=True),
        physiological_personalization=None,
        genetic_personalization=None,
    )


@pytest.fixture(scope="module")
def high_meal_output():
    return SimulationRunner().run(high_sugar_meal_scenario())


@pytest.fixture(scope="module")
def fasted_output():
    return SimulationRunner().run(fasted_scenario(duration_s=3600))


@pytest.fixture(scope="module")
def light_output():
    return SimulationRunner().run(light_snack_scenario())


class TestRiskBand:
    def test_enum_values(self):
        assert RiskBand.LOW.value == "low"
        assert RiskBand.HIGH.value == "high"


class TestRiskIndicatorScore:
    def test_construction(self):
        s = RiskIndicatorScore("Test", 7.5, "mmol/L", RiskBand.MODERATE, "Moderate")
        assert s.name == "Test"
        assert s.band == RiskBand.MODERATE


class TestComputeRiskProfile:
    def test_returns_profile(self, high_meal_output):
        profile = compute_risk_profile(high_meal_output)
        assert isinstance(profile, MetabolicRiskProfile)

    def test_has_indicators(self, high_meal_output):
        profile = compute_risk_profile(high_meal_output)
        assert len(profile.indicators) >= 3

    def test_peak_glucose_indicator_present(self, high_meal_output):
        profile = compute_risk_profile(high_meal_output)
        names = [i.name for i in profile.indicators]
        assert "Peak glucose" in names

    def test_peak_glucose_never_calls_a_who_borderline_peak_healthy(self):
        # Regression: LOW/MODERATE bands used to be gated purely on
        # excursion (peak - the patient's OWN baseline), never on the
        # absolute peak_mM -- for a T2D patient whose own baseline is
        # already elevated, a small further rise landed the absolute
        # peak in the WHO borderline range (7.8-11.1 mM) while still
        # being called "a healthy response" / "typical for this meal".
        from biosim_seeds.condition_scenarios import type2_diabetes_meal_scenario
        output = SimulationRunner().run(type2_diabetes_meal_scenario(carbs_g=25))
        assert output.expert_metrics.glucose_peak_mM >= 7.8, (
            "test setup assumption failed -- this scenario should reach "
            "the WHO borderline range for the regression to be meaningful"
        )
        profile = compute_risk_profile(output)
        peak_ind = next(i for i in profile.indicators if i.name == "Peak glucose")
        assert peak_ind.band != RiskBand.LOW
        assert "healthy response" not in peak_ind.interpretation
        assert "typical for this meal" not in peak_ind.interpretation

    def test_mtor_ampk_indicator_present(self, high_meal_output):
        profile = compute_risk_profile(high_meal_output)
        names = [i.name for i in profile.indicators]
        assert any("mTOR" in n for n in names)

    def test_fasted_all_low(self, fasted_output):
        profile = compute_risk_profile(fasted_output)
        peak_ind = next(i for i in profile.indicators if i.name == "Peak glucose")
        assert peak_ind.band in (RiskBand.LOW, RiskBand.MODERATE)

    def test_high_meal_worse_than_light(self, high_meal_output, light_output):
        band_order = [RiskBand.LOW, RiskBand.MODERATE, RiskBand.ELEVATED, RiskBand.HIGH]
        high_score  = band_order.index(compute_risk_profile(high_meal_output).overall_band)
        light_score = band_order.index(compute_risk_profile(light_output).overall_band)
        assert high_score >= light_score

    def test_headline_populated(self, high_meal_output):
        profile = compute_risk_profile(high_meal_output)
        assert len(profile.headline) > 10

    def test_to_dict_structure(self, high_meal_output):
        d = compute_risk_profile(high_meal_output).to_dict()
        assert "overall_band" in d
        assert "headline" in d
        assert "indicators" in d
        assert isinstance(d["indicators"], list)

    def test_to_dict_includes_disclaimer_and_indicator_caveats(self, high_meal_output):
        # Regression: to_dict() used to omit both `disclaimer` (the
        # class-level "this is a research tool, not a medical test" text)
        # and each indicator's `caveat` entirely -- real, meaningful
        # content silently dropped from every dict/JSON consumer of this
        # profile.
        profile = compute_risk_profile(high_meal_output)
        d = profile.to_dict()
        assert d["disclaimer"] == profile.disclaimer
        assert d["disclaimer"] != ""
        mtor_ind = next(i for i in d["indicators"] if "mTOR" in i["name"])
        assert mtor_ind["caveat"] != ""

    def test_n_elevated_count(self, high_meal_output):
        profile = compute_risk_profile(high_meal_output)
        count = sum(1 for i in profile.indicators
                    if i.band in (RiskBand.ELEVATED, RiskBand.HIGH))
        assert profile.n_elevated_or_high == count

    def test_disclaimer_present(self, high_meal_output):
        profile = compute_risk_profile(high_meal_output)
        assert len(profile.disclaimer) > 20

    def test_t2d_profile_elevated_risk(self):
        from biosim_physiology import PhysiologicalProfile, DiseaseState
        from biosim_api import PersonalizationSpec
        import dataclasses
        profile_t2d = PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES])
        scenario = high_sugar_meal_scenario(
            personalization=PersonalizationSpec(physiological_profile=profile_t2d)
        )
        output = SimulationRunner().run(scenario)
        risk = compute_risk_profile(output)
        band_order = [RiskBand.LOW, RiskBand.MODERATE, RiskBand.ELEVATED, RiskBand.HIGH]
        healthy_risk = compute_risk_profile(
            SimulationRunner().run(high_sugar_meal_scenario())
        )
        t2d_score     = band_order.index(risk.overall_band)
        healthy_score = band_order.index(healthy_risk.overall_band)
        assert t2d_score >= healthy_score


class TestComputeRiskProfileOnFailedSimulation:
    """
    Regression, found by the mandatory code-reviewer pass on the REST
    exposure of this function (MASTER_PLAN 3.11b): unlike report.py's
    assess_risk() (gated on overall_watchdog_confidence, turning GREY on
    a failed simulation), compute_risk_profile() never checked
    output.simulation.success at all. Reproduced directly against
    SimulationRunner._compute_expert_metrics()'s own documented failure-
    branch ExpertMetrics values: it reported overall_band=LOW and
    headline "Metabolic indicators are all in a healthy range for this
    meal" for a simulation that produced no real data -- a materially
    wrong, more-confident-sounding result than assess_risk() gives the
    identical failed SimulationOutput, now reachable over the public
    REST API for the first time via this same milestone's own change.
    """
    @staticmethod
    def _failed_output():
        # Exactly SimulationRunner._compute_expert_metrics()'s own
        # documented failure-branch values (runner.py), not an
        # arbitrary/synthetic set of zeros.
        metrics = ExpertMetrics(
            glucose_peak_uM=0, glucose_peak_time_s=0, glucose_auc_uM_s=0,
            glucose_baseline_uM=5000.0, glucose_nadir_uM=5000.0,
            insulin_peak_uM=0, insulin_peak_time_s=0,
            mtor_peak=0, ampk_nadir=1.0,
            time_to_glucose_return_s=None,
            overall_watchdog_confidence=0.0,
            n_ode_evaluations=0,
        )
        return SimpleNamespace(
            expert_metrics=metrics,
            simulation=SimpleNamespace(success=False),
            physiological_personalization=None,
            genetic_personalization=None,
        )

    def test_overall_band_is_unknown_not_low(self):
        profile = compute_risk_profile(self._failed_output())
        assert profile.overall_band == RiskBand.UNKNOWN

    def test_headline_does_not_claim_healthy(self):
        profile = compute_risk_profile(self._failed_output())
        assert "healthy" not in profile.headline.lower()
        assert "failed" in profile.headline.lower()

    def test_no_indicator_reports_a_low_confident_band(self):
        # Every indicator must be UNKNOWN, not a real LOW/MODERATE/
        # ELEVATED/HIGH band computed from meaningless failure-branch
        # zeros -- this is the exact bug: LOW ("healthy") was reachable
        # from data that was never really simulated.
        profile = compute_risk_profile(self._failed_output())
        assert all(i.band == RiskBand.UNKNOWN for i in profile.indicators)

    def test_successful_simulation_is_unaffected(self, high_meal_output):
        # The gate must not fire for a genuinely successful simulation.
        profile = compute_risk_profile(high_meal_output)
        assert profile.overall_band != RiskBand.UNKNOWN


class TestInsulinEfficiencyProxy:
    """
    Regression: this indicator divided glucose_auc_uM_s by 1e6 while
    claiming to convert uM.s -> mM.s (the correct divisor is 1e3) --
    produced a value 1000x too small, permanently pinned near the LOW
    band for every disease-state combination this platform can simulate,
    including T2D. Invisible to test_t2d_profile_elevated_risk above
    because overall_band is the WORST of all four indicators, and the
    other three already differentiated T2D from healthy on their own --
    this class checks the indicator itself, closing that blind spot.
    """
    @staticmethod
    def _indicator(output):
        profile = compute_risk_profile(output)
        return next(i for i in profile.indicators if i.name == "Insulin efficiency proxy")

    def test_unit_conversion_is_uM_s_to_mM_s(self, high_meal_output):
        # Directly pins the conversion factor so a future edit can't
        # silently reintroduce the /1e6 bug: value == (glucose_auc_uM_s
        # / 1e3) / insulin_peak_pM, not /1e6.
        m = high_meal_output.expert_metrics
        indicator = self._indicator(high_meal_output)
        expected = (m.glucose_auc_uM_s / 1e3) / m.insulin_peak_pM
        assert indicator.value == pytest.approx(round(expected, 2))

    def test_t2d_scores_worse_than_healthy_on_this_indicator_specifically(self):
        from biosim_physiology import PhysiologicalProfile, DiseaseState
        from biosim_api import PersonalizationSpec

        healthy = SimulationRunner().run(high_sugar_meal_scenario())
        t2d = SimulationRunner().run(high_sugar_meal_scenario(
            personalization=PersonalizationSpec(
                physiological_profile=PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES])
            )
        ))
        band_order = [RiskBand.LOW, RiskBand.MODERATE, RiskBand.ELEVATED, RiskBand.HIGH]
        healthy_band = self._indicator(healthy).band
        t2d_band = self._indicator(t2d).band
        assert band_order.index(t2d_band) > band_order.index(healthy_band), (
            "T2D must score strictly worse than healthy on the insulin "
            "efficiency proxy itself, not just on the overall (worst-of-4) band"
        )

    def test_t2d_is_not_low_band(self):
        from biosim_physiology import PhysiologicalProfile, DiseaseState
        from biosim_api import PersonalizationSpec

        t2d = SimulationRunner().run(high_sugar_meal_scenario(
            personalization=PersonalizationSpec(
                physiological_profile=PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES])
            )
        ))
        assert self._indicator(t2d).band != RiskBand.LOW

    def test_healthy_meals_stay_low(self, high_meal_output, light_output):
        assert self._indicator(high_meal_output).band == RiskBand.LOW
        assert self._indicator(light_output).band == RiskBand.LOW


class TestMTORAMPKBlutedSecretionCaveat:
    """
    Regression: mTOR/AMPK are driven only by insulin level in this model
    -- the T2D disease-state preset blunts insulin SECRETION
    (secretion_scale=0.60), which independently lowers mTOR/AMPK
    excursion for a reason unrelated to genuine metabolic balance. Before
    this fix, a T2D patient scored "Balanced anabolic/catabolic
    signalling" (LOW) on the mTOR/AMPK ratio at a meal that scored a
    healthy comparator MODERATE, with no caveat explaining why. Fixed
    with an explicit caveat rather than recalibrating the band
    thresholds (no principled clinical basis to do so).
    """
    @staticmethod
    def _ratio_indicator(output):
        profile = compute_risk_profile(output)
        return next(i for i in profile.indicators if i.name == "mTOR/AMPK ratio")

    def test_t2d_ratio_carries_blunted_secretion_caveat(self):
        from biosim_physiology import PhysiologicalProfile, DiseaseState
        from biosim_api import PersonalizationSpec
        t2d = SimulationRunner().run(high_sugar_meal_scenario(
            personalization=PersonalizationSpec(
                physiological_profile=PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES])
            )
        ))
        caveat = self._ratio_indicator(t2d).caveat.lower()
        assert "blunted insulin secretion" in caveat

    def test_healthy_ratio_has_no_blunted_secretion_caveat(self, high_meal_output):
        caveat = self._ratio_indicator(high_meal_output).caveat.lower()
        assert "blunted insulin secretion" not in caveat


class TestHypoglycaemiaRiskIndicator:
    """
    Battelino et al. 2019, PMID:31177185 (Diabetes Care) -- see report.py's
    _HYPO_LEVEL1_MM/_HYPO_LEVEL2_MM docstring for the full citation.
    Same milestone as report.py's own "Hypoglycaemia risk" indicator --
    this platform's SGLT2i/sulfonylurea mechanisms can genuinely push
    simulated glucose low, with nothing in this file's own output ever
    checking for it before.
    """

    @staticmethod
    def _indicator(output):
        profile = compute_risk_profile(output)
        return next(i for i in profile.indicators if i.name == "Hypoglycaemia risk")

    def test_indicator_present(self):
        assert self._indicator(_output_with_nadir(5.0)) is not None

    def test_low_above_alert_value(self):
        assert self._indicator(_output_with_nadir(4.5)).band == RiskBand.LOW

    def test_elevated_at_level_1(self):
        assert self._indicator(_output_with_nadir(3.5)).band == RiskBand.ELEVATED

    def test_high_at_level_2(self):
        assert self._indicator(_output_with_nadir(2.5)).band == RiskBand.HIGH

    def test_boundary_exactly_at_3_9_is_low(self):
        assert self._indicator(_output_with_nadir(3.9)).band == RiskBand.LOW

    def test_boundary_exactly_at_3_0_is_elevated(self):
        assert self._indicator(_output_with_nadir(3.0)).band == RiskBand.ELEVATED

    def test_fasted_scenario_nadir_equals_baseline(self, fasted_output):
        assert fasted_output.expert_metrics.glucose_nadir_mM == pytest.approx(
            fasted_output.expert_metrics.glucose_baseline_uM / 1000.0, abs=1e-6
        )
        assert self._indicator(fasted_output).band == RiskBand.LOW

    def test_sglt2i_scenario_genuinely_lowers_the_nadir(self):
        # Same real-wiring check as test_lay_output.py's identically-named
        # test -- confirms glucose_nadir_uM is a real, drug-sensitive
        # value through SimulationRunner, not asserting an absolute
        # clinical-hypoglycaemia value this platform's current drug
        # calibrations don't reach.
        no_drug = ScenarioDefinition(name="nd_fasted", meal_events=[], duration_s=14400)
        with_dapa = ScenarioDefinition(
            name="dapa_fasted", meal_events=[], duration_s=14400,
            drug_events=[DrugEvent(time_s=0, drug="dapagliflozin", dose_mg=10.0, hours_before_simulation=24 * 7)],
        )
        runner = SimulationRunner()
        nadir_no_drug = runner.run(no_drug).expert_metrics.glucose_nadir_mM
        nadir_dapa = runner.run(with_dapa).expert_metrics.glucose_nadir_mM
        assert nadir_dapa < nadir_no_drug

    def test_low_result_carries_a_calibration_caveat(self):
        # Same regression as test_lay_output.py's identically-named test.
        ind = self._indicator(_output_with_nadir(4.5))
        assert ind.band == RiskBand.LOW
        assert "not a guarantee against real hypoglycaemia" in ind.caveat

    def test_elevated_result_does_not_carry_the_calibration_caveat(self):
        ind = self._indicator(_output_with_nadir(3.5))
        assert ind.band == RiskBand.ELEVATED
        assert "not a guarantee against real hypoglycaemia" not in ind.caveat
