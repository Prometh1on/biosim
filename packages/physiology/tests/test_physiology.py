# SAVE TO: packages\physiology\tests\test_physiology.py
"""
Test suite for biosim_physiology.

Coverage:
    - PhysiologicalProfile: construction, BMI, validation
    - Age modifiers: sensitivity and secretion decline with age
    - Sex modifiers: sex-specific parameter differences
    - BMI modifiers: increasing adiposity reduces insulin sensitivity
    - Fitness modifiers: athletes more sensitive than sedentary
    - Disease state modifiers: T2D, insulin resistance, hypothyroidism
    - Chronobiology: morning vs evening sensitivity difference
    - Drug interactions: caffeine, metformin, glucocorticoid
    - PhysiologicalPersonalizer: full pipeline, ODE still runs, validate_y0
    - Composition: multiple modifiers combine correctly
"""

import math
import pytest

from biosim_physiology import (
    BiologicalSex,
    DiseaseState,
    Drug,
    DrugDose,
    FitnessLevel,
    PhysiologicalPersonalizationResult,
    PhysiologicalPersonalizer,
    PhysiologicalProfile,
    age_modifiers,
    bmi_modifiers,
    circadian_factor,
    circadian_modifiers,
    circadian_sensitivity_percent,
    cold_acclimation_modifiers,
    compose_all_modifiers,
    disease_modifiers,
    drug_modifiers,
    fitness_modifiers,
    intermittent_hypoxia_modifiers,
    menstrual_cycle_modifiers,
    MenstrualCyclePhase,
    sex_modifiers,
    sleep_modifiers,
)
from biosim_ode import GlucoseInsulinPathway, ODESystem, PathwayParameters, high_sugar_meal


# ---------------------------------------------------------------------------
# PhysiologicalProfile
# ---------------------------------------------------------------------------

class TestPhysiologicalProfile:
    def test_default_construction(self):
        p = PhysiologicalProfile()
        assert p.age_years == 30.0
        assert p.bmi == pytest.approx(70.0 / (1.70 ** 2), rel=0.001)

    def test_bmi_calculation(self):
        p = PhysiologicalProfile(body_weight_kg=90.0, height_cm=180.0)
        assert p.bmi == pytest.approx(27.78, rel=0.01)

    def test_bmi_category_normal(self):
        assert PhysiologicalProfile(body_weight_kg=70, height_cm=175).bmi_category == "normal"

    def test_bmi_category_overweight(self):
        assert PhysiologicalProfile(body_weight_kg=90, height_cm=175).bmi_category == "overweight"

    def test_bmi_category_obese(self):
        assert PhysiologicalProfile(body_weight_kg=98, height_cm=175).bmi_category == "obese"

    def test_lean_mass_when_body_fat_known(self):
        p = PhysiologicalProfile(body_weight_kg=80.0, body_fat_percent=20.0)
        assert p.lean_mass_kg == pytest.approx(64.0)

    def test_lean_mass_none_when_unknown(self):
        p = PhysiologicalProfile()
        assert p.lean_mass_kg is None

    def test_is_diabetic(self):
        p = PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES])
        assert p.is_diabetic

    def test_is_not_diabetic(self):
        p = PhysiologicalProfile()
        assert not p.is_diabetic

    def test_invalid_age_raises(self):
        with pytest.raises(ValueError):
            PhysiologicalProfile(age_years=150)

    def test_invalid_weight_raises(self):
        with pytest.raises(ValueError):
            PhysiologicalProfile(body_weight_kg=0)

    def test_invalid_time_of_day_raises(self):
        with pytest.raises(ValueError):
            PhysiologicalProfile(time_of_day_hours=25.0)

    def test_invalid_sleep_hours_raises(self):
        with pytest.raises(ValueError):
            PhysiologicalProfile(sleep_hours=25.0)

    def test_negative_sleep_hours_raises(self):
        with pytest.raises(ValueError):
            PhysiologicalProfile(sleep_hours=-1.0)

    def test_sleep_hours_defaults_to_none(self):
        assert PhysiologicalProfile().sleep_hours is None

    def test_negative_desaturation_events_per_hour_raises(self):
        with pytest.raises(ValueError):
            PhysiologicalProfile(desaturation_events_per_hour=-1.0)

    def test_nan_desaturation_events_per_hour_raises(self):
        # Code-reviewer LOW finding: a plain `< 0` check silently accepts
        # NaN (float('nan') < 0 is False in Python) -- unlike sleep_hours's
        # own chained-range check, which incidentally rejects NaN. Fixed
        # with an explicit math.isfinite() check for parity.
        with pytest.raises(ValueError):
            PhysiologicalProfile(desaturation_events_per_hour=float("nan"))

    def test_desaturation_events_per_hour_defaults_to_none(self):
        assert PhysiologicalProfile().desaturation_events_per_hour is None

    def test_menstrual_cycle_phase_defaults_to_none(self):
        assert PhysiologicalProfile().menstrual_cycle_phase is None

    def test_menstrual_cycle_phase_accepts_a_real_enum_member(self):
        p = PhysiologicalProfile(menstrual_cycle_phase=MenstrualCyclePhase.LUTEAL)
        assert p.menstrual_cycle_phase == MenstrualCyclePhase.LUTEAL

    def test_cold_acclimatized_defaults_to_false(self):
        assert PhysiologicalProfile().cold_acclimatized is False

    def test_with_time(self):
        p = PhysiologicalProfile()
        p2 = p.with_time(14.0)
        assert p2.time_of_day_hours == 14.0
        assert p.time_of_day_hours is None  # original unchanged

    def test_duplicate_disease_state_raises(self):
        # Found during a security-reviewer pass on the HYPERTHYROIDISM
        # addition: an unvalidated duplicate silently double-applies that
        # state's modifier (e.g. T2D listed twice squares Gb/p3/
        # secretion_scale -- a 2x-compounded, undocumented severity).
        with pytest.raises(ValueError):
            PhysiologicalProfile(disease_states=[
                DiseaseState.TYPE_2_DIABETES, DiseaseState.TYPE_2_DIABETES,
            ])

    def test_hypo_and_hyperthyroidism_together_raises(self):
        # A patient cannot genuinely be both hypo- and hyperthyroid.
        # Worse than a generic contradiction: HYPERTHYROIDISM's static
        # modifier is deliberately {} (see modifiers.py), so without this
        # check the combination would be silently indistinguishable from
        # HYPOTHYROIDISM alone at the static-modifier layer while STILL
        # composing for_hyperthyroidism()'s ODE pathway in
        # SimulationRunner -- a physiologically nonsensical mix of both
        # conditions' effects with no error at all.
        with pytest.raises(ValueError):
            PhysiologicalProfile(disease_states=[
                DiseaseState.HYPOTHYROIDISM, DiseaseState.HYPERTHYROIDISM,
            ])

    def test_t1d_and_t2d_together_raises(self):
        # T1D (near-total secretion loss) and T2D (resistance-driven,
        # partially-preserved secretion) are opposite beta-cell
        # mechanisms -- combining them would multiply TYPE_1_DIABETES's
        # near-zero Ib/secretion_scale factors by TYPE_2_DIABETES's own,
        # a physiologically nonsensical "both at once" state.
        with pytest.raises(ValueError):
            PhysiologicalProfile(disease_states=[
                DiseaseState.TYPE_1_DIABETES, DiseaseState.TYPE_2_DIABETES,
            ])

    def test_t1d_and_insulin_resistant_together_raises(self):
        with pytest.raises(ValueError):
            PhysiologicalProfile(disease_states=[
                DiseaseState.TYPE_1_DIABETES, DiseaseState.INSULIN_RESISTANT,
            ])

    def test_t1d_alone_is_valid(self):
        # Sanity check that the new mutual-exclusivity guard doesn't
        # false-positive on the legitimate, single-state case.
        profile = PhysiologicalProfile(disease_states=[DiseaseState.TYPE_1_DIABETES])
        assert profile.disease_states == [DiseaseState.TYPE_1_DIABETES]

    def test_gestational_diabetes_and_t1d_together_raises(self):
        # GDM is, by clinical definition, diabetes first recognized
        # during pregnancy in someone WITHOUT pre-existing diabetes --
        # pre-existing T1D/T2D in pregnancy is "pregestational diabetes",
        # a distinct entity this platform does not model.
        with pytest.raises(ValueError):
            PhysiologicalProfile(disease_states=[
                DiseaseState.GESTATIONAL_DIABETES, DiseaseState.TYPE_1_DIABETES,
            ])

    def test_gestational_diabetes_and_t2d_together_raises(self):
        with pytest.raises(ValueError):
            PhysiologicalProfile(disease_states=[
                DiseaseState.GESTATIONAL_DIABETES, DiseaseState.TYPE_2_DIABETES,
            ])

    def test_gestational_diabetes_alone_is_valid(self):
        profile = PhysiologicalProfile(disease_states=[DiseaseState.GESTATIONAL_DIABETES])
        assert profile.disease_states == [DiseaseState.GESTATIONAL_DIABETES]

    def test_gestational_diabetes_and_insulin_resistant_together_is_valid(self):
        # Deliberately NOT excluded (unlike T1D/T2D above) -- pre-
        # existing insulin resistance/pre-diabetes is a real, common GDM
        # risk factor, not a contradictory diagnosis.
        profile = PhysiologicalProfile(disease_states=[
            DiseaseState.GESTATIONAL_DIABETES, DiseaseState.INSULIN_RESISTANT,
        ])
        assert set(profile.disease_states) == {
            DiseaseState.GESTATIONAL_DIABETES, DiseaseState.INSULIN_RESISTANT,
        }

    def test_gck_mody_and_t1d_together_raises(self):
        with pytest.raises(ValueError):
            PhysiologicalProfile(disease_states=[
                DiseaseState.GCK_MODY, DiseaseState.TYPE_1_DIABETES,
            ])

    def test_gck_mody_and_t2d_together_raises(self):
        with pytest.raises(ValueError):
            PhysiologicalProfile(disease_states=[
                DiseaseState.GCK_MODY, DiseaseState.TYPE_2_DIABETES,
            ])

    def test_gck_mody_alone_is_valid(self):
        profile = PhysiologicalProfile(disease_states=[DiseaseState.GCK_MODY])
        assert profile.disease_states == [DiseaseState.GCK_MODY]

    def test_gck_mody_and_insulin_resistant_together_is_valid(self):
        # GCK-MODY's threshold-reset mechanism is independent of
        # insulin resistance -- not a contradictory combination.
        profile = PhysiologicalProfile(disease_states=[
            DiseaseState.GCK_MODY, DiseaseState.INSULIN_RESISTANT,
        ])
        assert set(profile.disease_states) == {
            DiseaseState.GCK_MODY, DiseaseState.INSULIN_RESISTANT,
        }

    def test_osa_and_desaturation_events_per_hour_together_raises(self):
        # Track 1 row 1.39: the first disease-state x profile-field
        # mutual exclusion -- OSA (chronic diagnosis) and
        # desaturation_events_per_hour (acute measured exposure rate)
        # represent the same real oxyhemoglobin-desaturation mechanism
        # via two independent p3-lowering routes.
        with pytest.raises(ValueError):
            PhysiologicalProfile(
                disease_states=[DiseaseState.OSA],
                desaturation_events_per_hour=24.3,
            )

    def test_osa_alone_is_valid(self):
        profile = PhysiologicalProfile(disease_states=[DiseaseState.OSA])
        assert profile.disease_states == [DiseaseState.OSA]

    def test_desaturation_events_per_hour_alone_is_valid(self):
        profile = PhysiologicalProfile(
            disease_states=[DiseaseState.HEALTHY], desaturation_events_per_hour=24.3,
        )
        assert profile.desaturation_events_per_hour == pytest.approx(24.3)

    def test_osa_and_type_2_diabetes_together_is_valid(self):
        # OSA-related insulin resistance is a real, common T2D
        # comorbidity -- deliberately not mutually exclusive with any
        # disease state, only with the profile-field it double-counts.
        profile = PhysiologicalProfile(disease_states=[
            DiseaseState.OSA, DiseaseState.TYPE_2_DIABETES,
        ])
        assert set(profile.disease_states) == {
            DiseaseState.OSA, DiseaseState.TYPE_2_DIABETES,
        }

    def test_osa_and_zero_desaturation_events_per_hour_together_is_valid(self):
        # Code-reviewer finding: intermittent_hypoxia_modifiers() itself
        # treats exactly 0 events/hour as a neutral no-op (returns {}),
        # so composing OSA with a MEASURED rate of 0.0 double-counts
        # nothing -- only a genuinely nonzero rate is excluded.
        profile = PhysiologicalProfile(
            disease_states=[DiseaseState.OSA], desaturation_events_per_hour=0.0,
        )
        assert profile.desaturation_events_per_hour == 0.0

    def test_osa_and_nonzero_desaturation_events_per_hour_together_raises(self):
        with pytest.raises(ValueError):
            PhysiologicalProfile(
                disease_states=[DiseaseState.OSA], desaturation_events_per_hour=0.1,
            )


# ---------------------------------------------------------------------------
# Age modifiers
# ---------------------------------------------------------------------------

class TestAgeModifiers:
    def test_reference_age_neutral(self):
        m = age_modifiers(30.0)
        assert m["p3"] == pytest.approx(1.0)
        assert m["secretion_scale"] == pytest.approx(1.0)

    def test_older_lower_p3(self):
        m30 = age_modifiers(30)
        m60 = age_modifiers(60)
        assert m60["p3"] < m30["p3"]

    def test_older_lower_secretion(self):
        m30 = age_modifiers(30)
        m70 = age_modifiers(70)
        assert m70["secretion_scale"] < m30["secretion_scale"]

    def test_floor_not_exceeded(self):
        m = age_modifiers(120)
        assert m["p3"] >= 0.30
        assert m["secretion_scale"] >= 0.50

    def test_young_no_penalty(self):
        m = age_modifiers(20)
        assert m["p3"] == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# Sex modifiers
# ---------------------------------------------------------------------------

class TestSexModifiers:
    def test_male_p1_higher(self):
        male = sex_modifiers(BiologicalSex.MALE)
        female = sex_modifiers(BiologicalSex.FEMALE)
        assert male["p1"] > female["p1"]

    def test_female_p3_higher(self):
        male = sex_modifiers(BiologicalSex.MALE)
        female = sex_modifiers(BiologicalSex.FEMALE)
        assert female["p3"] > male["p3"]

    def test_not_specified_no_modifier(self):
        m = sex_modifiers(BiologicalSex.NOT_SPECIFIED)
        assert m == {}


# ---------------------------------------------------------------------------
# BMI modifiers
# ---------------------------------------------------------------------------

class TestBMIModifiers:
    def test_normal_bmi_no_change(self):
        m = bmi_modifiers(22.0)
        assert m == {}

    def test_overweight_reduces_p3(self):
        m = bmi_modifiers(27.0)
        assert m["p3"] < 1.0

    def test_obese_more_reduction(self):
        m_ow = bmi_modifiers(27.0)
        m_ob = bmi_modifiers(32.0)
        assert m_ob["p3"] < m_ow["p3"]

    def test_severely_obese_largest_reduction(self):
        m = bmi_modifiers(40.0)
        assert m["p3"] == pytest.approx(0.50)


# ---------------------------------------------------------------------------
# Fitness modifiers
# ---------------------------------------------------------------------------

class TestFitnessModifiers:
    def test_moderate_reference(self):
        m = fitness_modifiers(FitnessLevel.MODERATELY_ACTIVE)
        assert m["p3"] == pytest.approx(1.0)

    def test_athlete_higher_than_sedentary(self):
        ath = fitness_modifiers(FitnessLevel.ATHLETE)
        sed = fitness_modifiers(FitnessLevel.SEDENTARY)
        assert ath["p3"] > sed["p3"]

    def test_sedentary_below_1(self):
        m = fitness_modifiers(FitnessLevel.SEDENTARY)
        assert m["p3"] < 1.0

    def test_athlete_above_1(self):
        m = fitness_modifiers(FitnessLevel.ATHLETE)
        assert m["p3"] > 1.0


# ---------------------------------------------------------------------------
# Sleep restriction modifiers
# ---------------------------------------------------------------------------

class TestSleepModifiers:
    """
    Fresh Track 1 scoping (2026-08-17/18), row 1.26: sleep-restriction
    p3 modifier, anchored to Buxton et al. 2010 (PMID:20585000,
    independently live-verified by the orchestrating session before
    implementation). See sleep_modifiers()'s own docstring for the full
    derivation.
    """

    def test_none_returns_no_adjustment(self):
        assert sleep_modifiers(None) == {}

    def test_ten_hours_is_neutral(self):
        assert sleep_modifiers(10.0) == {}

    def test_above_ten_hours_is_also_neutral_not_extrapolated(self):
        # No data above the paper's own sleep-replete condition (10h) --
        # not modelled as beneficial, just genuinely unmodified.
        assert sleep_modifiers(12.0) == {}

    def test_five_hours_matches_buxton_2010_clamp_derived_reduction(self):
        m = sleep_modifiers(5.0)
        assert m == {"p3": pytest.approx(0.89)}

    def test_below_five_hours_is_clamped_not_extrapolated(self):
        assert sleep_modifiers(2.0) == sleep_modifiers(5.0)
        assert sleep_modifiers(0.0) == sleep_modifiers(5.0)

    def test_intermediate_duration_interpolates_between_anchors(self):
        f6, f7, f8 = (sleep_modifiers(h)["p3"] for h in (6.0, 7.0, 8.0))
        assert f6 < f7 < f8 < 1.0

    def test_negative_sleep_hours_raises(self):
        with pytest.raises(ValueError):
            sleep_modifiers(-1.0)

    def test_nan_sleep_hours_raises(self):
        with pytest.raises(ValueError):
            sleep_modifiers(float("nan"))

    def test_does_not_touch_secretion(self):
        # Buxton 2010's own confirmed null: "without significant
        # alterations in the insulin secretory response."
        m = sleep_modifiers(5.0)
        assert "secretion_scale" not in m


# ---------------------------------------------------------------------------
# Intermittent hypoxia modifiers
# ---------------------------------------------------------------------------

class TestIntermittentHypoxiaModifiers:
    """
    Track 1 row 1.32: acute intermittent-hypoxia p3/p1 modifier, anchored
    to Louis & Punjabi 2009 (PMID:19265062, independently live-verified
    twice by two separate parties before implementation). See
    intermittent_hypoxia_modifiers()'s own docstring for the full
    derivation.
    """

    def test_none_returns_no_adjustment(self):
        assert intermittent_hypoxia_modifiers(None) == {}

    def test_zero_events_is_neutral(self):
        assert intermittent_hypoxia_modifiers(0.0) == {}

    def test_reference_rate_matches_louis_punjabi_2009_ratios(self):
        m = intermittent_hypoxia_modifiers(24.3)
        assert m["p3"] == pytest.approx(3.4 / 4.1)
        assert m["p1"] == pytest.approx(1.3 / 1.9)

    def test_above_reference_rate_is_clamped_not_extrapolated(self):
        # No dose-response data beyond the paper's own single tested rate.
        assert intermittent_hypoxia_modifiers(50.0) == intermittent_hypoxia_modifiers(24.3)
        assert intermittent_hypoxia_modifiers(200.0) == intermittent_hypoxia_modifiers(24.3)

    def test_intermediate_rate_interpolates_between_anchors(self):
        f6, f12, f18 = (
            intermittent_hypoxia_modifiers(r)["p3"] for r in (6.0, 12.0, 18.0)
        )
        assert 1.0 > f6 > f12 > f18 > intermittent_hypoxia_modifiers(24.3)["p3"]
        g6, g12, g18 = (
            intermittent_hypoxia_modifiers(r)["p1"] for r in (6.0, 12.0, 18.0)
        )
        assert 1.0 > g6 > g12 > g18 > intermittent_hypoxia_modifiers(24.3)["p1"]

    def test_negative_events_per_hour_raises(self):
        with pytest.raises(ValueError):
            intermittent_hypoxia_modifiers(-1.0)

    def test_nan_events_per_hour_raises(self):
        with pytest.raises(ValueError):
            intermittent_hypoxia_modifiers(float("nan"))

    def test_does_not_touch_secretion(self):
        # Louis & Punjabi 2009's own confirmed Results null: "pancreatic
        # insulin secretion was comparable between the two conditions."
        m = intermittent_hypoxia_modifiers(24.3)
        assert "secretion_scale" not in m

    def test_p1_effect_is_larger_than_p3_effect_at_reference_rate(self):
        # The paper's own S(G) (p1) effect is larger than its S(I) (p3)
        # effect at the tested rate -- 31.6% vs. 17.1% reduction -- a
        # genuine asymmetry worth pinning as a regression, not an
        # implementation artifact.
        m = intermittent_hypoxia_modifiers(24.3)
        assert (1.0 - m["p1"]) > (1.0 - m["p3"])


# ---------------------------------------------------------------------------
# Menstrual cycle phase modifiers
# ---------------------------------------------------------------------------

class TestMenstrualCycleModifiers:
    """
    Track 1 row 1.36: menstrual cycle phase p3 modifier, anchored to
    Valdes & Elkind-Hirsch 1991 (PMID:1997519, independently
    live-verified before implementation). See
    menstrual_cycle_modifiers()'s own docstring for the full derivation,
    including the correction that MIDCYCLE is neutral (not a smaller
    real effect) per the paper's own explicit non-significance claim.
    """

    def test_none_returns_no_adjustment(self):
        assert menstrual_cycle_modifiers(None) == {}

    def test_follicular_is_neutral(self):
        assert menstrual_cycle_modifiers(MenstrualCyclePhase.FOLLICULAR) == {}

    def test_midcycle_is_neutral_not_a_smaller_effect(self):
        # The paper's own conclusion: "no significant change at
        # midcycle" -- this is a correction to the original candidate's
        # own scoping, which had (incorrectly) proposed a 3-point ramp.
        assert menstrual_cycle_modifiers(MenstrualCyclePhase.MIDCYCLE) == {}

    def test_luteal_matches_valdes_1991(self):
        m = menstrual_cycle_modifiers(MenstrualCyclePhase.LUTEAL)
        assert m == {"p3": pytest.approx(3.20 / 6.20)}

    def test_luteal_reduces_p3_below_neutral(self):
        m = menstrual_cycle_modifiers(MenstrualCyclePhase.LUTEAL)
        assert m["p3"] < 1.0

    def test_does_not_touch_p1(self):
        # Valdes & Elkind-Hirsch's own confirmed null: "No change in
        # glucose effectiveness occurred as a function of the
        # menstrual cycle."
        m = menstrual_cycle_modifiers(MenstrualCyclePhase.LUTEAL)
        assert "p1" not in m


class TestColdAcclimationModifiers:
    """
    Track 1 row 1.40: this platform's first p3-INCREASING environmental
    axis, anchored to Hanssen et al. 2015 (PMID:26147760, independently
    live-verified before implementation). Deliberately scoped to ONLY
    the chronic ~10-day acclimation regime -- NOT acute cold/shivering
    exposure, which Sellers et al. 2021 (PMID:33090911) found moves in
    the OPPOSITE direction in a genuinely different regime.
    """

    def test_false_returns_no_adjustment(self):
        assert cold_acclimation_modifiers(False) == {}

    def test_true_matches_hanssen_2015(self):
        m = cold_acclimation_modifiers(True)
        assert m == {"p3": pytest.approx(1.43)}

    def test_true_increases_p3_above_neutral(self):
        m = cold_acclimation_modifiers(True)
        assert m["p3"] > 1.0


# ---------------------------------------------------------------------------
# Disease state modifiers
# ---------------------------------------------------------------------------

class TestDiseaseModifiers:
    def test_healthy_no_change(self):
        m = disease_modifiers([DiseaseState.HEALTHY])
        assert m == {}

    def test_t2d_reduces_secretion(self):
        m = disease_modifiers([DiseaseState.TYPE_2_DIABETES])
        assert m["secretion_scale"] < 1.0

    def test_t2d_reduces_p3(self):
        m = disease_modifiers([DiseaseState.TYPE_2_DIABETES])
        assert m["p3"] < 1.0

    def test_t2d_elevates_gb(self):
        m = disease_modifiers([DiseaseState.TYPE_2_DIABETES])
        assert m["Gb"] > 1.0

    def test_t2d_reduces_p1_by_basu_2009_glucose_effectiveness_ratio(self):
        # MASTER_PLAN row 1.46 (a): Basu 2009 (PMID:19196896), same
        # 14-T2D vs. 11-matched-control cohort as validate_t2d_meal_insulin():
        # oral-minimal-model SG 0.019 vs. 0.034 min-1.
        m = disease_modifiers([DiseaseState.TYPE_2_DIABETES])
        assert m["p1"] == pytest.approx(0.019 / 0.034, abs=0.005)
        assert m["p1"] == 0.56

    def test_t2d_p1_factor_reaches_personalized_parameters_unclamped(self):
        result = PhysiologicalPersonalizer().personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES])
        )
        assert result.personalized_params.p1 == pytest.approx(PathwayParameters().p1 * 0.56)

    def test_worst_reachable_p1_stack_keeps_headroom_above_the_floor(self):
        # Security-review sweep (2026-10-23): the lowest reachable p1 stack is
        # T2D + hypothyroidism + female + saturated hypoxia + 22:00 + age 85,
        # p1 ~0.0088 (~9x the 0.001 floor). Guards the headroom against a
        # future p1-lowering mechanism (GDM's 0.67 cannot combine with T2D).
        result = PhysiologicalPersonalizer().personalize(
            PhysiologicalProfile(
                age_years=85,
                biological_sex=BiologicalSex.FEMALE,
                disease_states=[DiseaseState.TYPE_2_DIABETES, DiseaseState.HYPOTHYROIDISM],
                desaturation_events_per_hour=200,
                time_of_day_hours=22.0,
            )
        )
        assert result.personalized_params.p1 >= 0.005
        assert not any("p1" in w for w in result.clamp_warnings)

    def test_insulin_resistant_reduces_p3(self):
        m = disease_modifiers([DiseaseState.INSULIN_RESISTANT])
        assert m["p3"] < 1.0

    def test_hypothyroid_reduces_p1(self):
        m = disease_modifiers([DiseaseState.HYPOTHYROIDISM])
        assert m["p1"] < 1.0

    def test_hypothyroid_no_longer_elevates_fasting_glucose(self):
        # Regression: Gb=1.04 was removed -- uncited, and actively
        # contradicted by its own cited source (Brenta 2011, PMID:
        # 21941681, verified against its own full text: "a poor
        # utilization of glucose in hypothyroidism may be offset by a
        # reduced release to circulation maintaining a balance") plus
        # two independent human clamp studies both finding fasting
        # glucose unchanged after correcting hypothyroidism (Handisurya
        # 2008, PMID:18429948; Rochon 2003, PMID:12519082). No source
        # supports an elevated fasting glucose in hypothyroidism.
        m = disease_modifiers([DiseaseState.HYPOTHYROIDISM])
        assert "Gb" not in m

    def test_hyperthyroidism_has_no_static_modifier(self):
        # Deliberately {} -- hyperthyroidism's real glucose effect is
        # hepatic and insulin-gated (Wennlund et al. 1986, PMID:3510000),
        # represented via EGPParameters.for_hyperthyroidism()/
        # HepGlucPathway instead. A static modifier here would
        # double-count, the same bug class already fixed for
        # glucocorticoid/CortisolPathway.
        m = disease_modifiers([DiseaseState.HYPERTHYROIDISM])
        assert m == {}

    def test_comorbidity_compounds(self):
        m_ir = disease_modifiers([DiseaseState.INSULIN_RESISTANT])
        m_both = disease_modifiers([DiseaseState.INSULIN_RESISTANT, DiseaseState.HYPOTHYROIDISM])
        assert m_both["p3"] < m_ir["p3"]

    def test_t1d_reduces_ib_and_secretion_scale(self):
        m = disease_modifiers([DiseaseState.TYPE_1_DIABETES])
        assert m["Ib"] < 1.0
        assert m["secretion_scale"] < 1.0

    def test_t1d_elevates_gb_more_than_t2d(self):
        # T1D represents near-total insulin deficiency, a more severe
        # mechanism than T2D's resistance-with-partial-secretion -- its
        # disclosed engineering-choice Gb factor should exceed T2D's own
        # literature-anchored 1.30.
        m_t1d = disease_modifiers([DiseaseState.TYPE_1_DIABETES])
        m_t2d = disease_modifiers([DiseaseState.TYPE_2_DIABETES])
        assert m_t1d["Gb"] > m_t2d["Gb"]

    def test_gestational_diabetes_reduces_p1_and_p3(self):
        m = disease_modifiers([DiseaseState.GESTATIONAL_DIABETES])
        assert m["p1"] < 1.0
        assert m["p3"] < 1.0

    def test_gestational_diabetes_compensatory_hypersecretion(self):
        # Real, measured (not assumed) compensatory increase -- Kautzky-
        # Willer et al. 1997 (PMID:9353615): "insulin secretion was 30%
        # higher... in subjects with GDM than in pregnant women with
        # NGT or in nonpregnant lean women."
        m = disease_modifiers([DiseaseState.GESTATIONAL_DIABETES])
        assert m["secretion_scale"] > 1.0

    def test_gestational_diabetes_gb_matches_iadpsg_threshold(self):
        # Metzger et al. 2010 (PMID:20190296, IADPSG consensus, verified
        # against the primary source's own full text): "fasting plasma
        # glucose >=5.1 mmol/l (92 mg/dl)... diagnose as GDM."
        m = disease_modifiers([DiseaseState.GESTATIONAL_DIABETES])
        assert m["Gb"] * 5.0 == pytest.approx(5.1, abs=0.01)

    def test_gestational_diabetes_p3_more_severely_blunted_than_t2d(self):
        # 84% blunting (GDM, vs. a true nonpregnant reference) vs. T2D's
        # own 55% reduction -- a real, quantified difference from an
        # independently-anchored source, not an assumed ordering.
        m_gdm = disease_modifiers([DiseaseState.GESTATIONAL_DIABETES])
        m_t2d = disease_modifiers([DiseaseState.TYPE_2_DIABETES])
        assert m_gdm["p3"] < m_t2d["p3"]

    def test_gck_mody_elevates_gb(self):
        m = disease_modifiers([DiseaseState.GCK_MODY])
        assert m["Gb"] > 1.0

    def test_gck_mody_gb_matches_guenat_2000_ratio(self):
        # Guenat et al. 2000 (PMID:11119017): MODY2 basal glucose
        # 7.6 mmol/l vs. healthy 5.3 mmol/l in the SAME study --
        # applied as a within-study ratio to this platform's own
        # default (5.0*1.434=7.17 mM).
        m = disease_modifiers([DiseaseState.GCK_MODY])
        assert m["Gb"] * 5.0 == pytest.approx(7.6 / 5.3 * 5.0, abs=0.01)

    def test_gck_mody_does_not_touch_insulin_sensitivity_or_secretion(self):
        # The defining, mechanistically distinguishing fact for this
        # preset (Guenat et al. 2000, same citation): "Insulin-mediated
        # glucose disposal... were unchanged in MODY2 patients" -- unlike
        # every other diabetes-related state in this file, p1/p3/
        # secretion_scale must be entirely ABSENT from the returned
        # dict, not merely set to a neutral 1.0 (which would still be a
        # claim about the parameter -- omission is the correct way to
        # express "no effect", per compose_all_modifiers()'s own
        # documented convention, same precedent as the earlier p1/p3
        # removals from metformin_modifiers()).
        m = disease_modifiers([DiseaseState.GCK_MODY])
        assert "p1" not in m
        assert "p3" not in m
        assert "secretion_scale" not in m

    def test_ckd_reduces_p3_and_ke_matching_de_boer_2016_adjusted_ratios(self):
        # de Boer et al. 2016 (PMID:26823551): adjusted difference of
        # -0.7 mg/min per uU/ml on a control mean of 5.0 -> p3 factor
        # (5.0-0.7)/5.0=0.86; adjusted difference of -85 ml/min on a
        # control mean of 998 -> ke factor (998-85)/998=0.915.
        m = disease_modifiers([DiseaseState.CHRONIC_KIDNEY_DISEASE])
        assert m["p3"] == pytest.approx(0.86)
        assert m["ke"] == pytest.approx(0.915)

    def test_ckd_does_not_touch_secretion_or_fasting_levels(self):
        # The same source explicitly confirms no CKD-specific secretion
        # deficit ("We did not observe a CKD-specific deficit in insulin
        # secretion"), and reports no absolute Gb/Ib anchor at all --
        # omission is the correct way to express "not modified here",
        # same convention as GCK_MODY's own p1/p3/secretion_scale omission
        # above.
        #
        # Asserts the FULL dict, not just the absent keys (a mandatory
        # code-reviewer pass correctly found the absent-keys-only version
        # of this test passed identically against a no-op {} modifier --
        # not RED-proof for CKD's own p3/ke values specifically). Exact
        # equality discriminates both "these keys are absent" and "the
        # modifier isn't accidentally empty" in one assertion.
        m = disease_modifiers([DiseaseState.CHRONIC_KIDNEY_DISEASE])
        assert m == {"p3": pytest.approx(0.86), "ke": pytest.approx(0.915)}

    def test_pcos_reduces_p3_matching_ketel_2011_obese_stratum(self):
        # Ketel et al. 2011 (PMID:21920943): obese PCOS vs. obese
        # controls M/I ratio 0.5/0.8=0.625, the significant (P<0.01)
        # comparison -- the anchor used, disclosed as obese-specific.
        m = disease_modifiers([DiseaseState.PCOS])
        assert m["p3"] == pytest.approx(0.625)

    def test_pcos_does_not_touch_secretion_ke_or_gb(self):
        # ke confirmed unaffected by Dunaif 1989 ("the metabolic
        # clearance rate of insulin did not differ"); Gb/hepatic glucose
        # production confirmed unaffected by Lewy 2001 ("Fasting glucose
        # and hepatic glucose production were comparable"); no source
        # reports an independent secretion-capacity magnitude.
        #
        # Asserts the FULL dict, not just the absent keys -- same
        # discriminating-test discipline already established for CKD
        # above (a mandatory code-reviewer finding on that chunk).
        m = disease_modifiers([DiseaseState.PCOS])
        assert m == {"p3": pytest.approx(0.625)}

    def test_acromegaly_raises_gb_matching_moller_1992(self):
        # Moller et al. 1992 (PMID:1569148): active acromegaly vs.
        # controls fasting glucose ratio 5.9/5.2=1.1346.
        m = disease_modifiers([DiseaseState.ACROMEGALY])
        assert m["Gb"] == pytest.approx(5.9 / 5.2)

    def test_acromegaly_does_not_touch_p3_ke_or_secretion(self):
        # The source reports insulin resistance only qualitatively (no
        # M-value/GIR number) -- an honest partial implementation, not a
        # guessed magnitude. Asserts the FULL dict, same discriminating-
        # test discipline as CKD/PCOS above.
        m = disease_modifiers([DiseaseState.ACROMEGALY])
        assert m == {"Gb": pytest.approx(5.9 / 5.2)}

    def test_osa_reduces_p3_matching_harsch_2004_baseline(self):
        # Harsch et al. 2004 (PMID:14512265): pre-CPAP (untreated OSA)
        # vs. 2-days-post-CPAP clamp M-value ratio 5.75/6.79=0.8469.
        m = disease_modifiers([DiseaseState.OSA])
        assert m["p3"] == pytest.approx(5.75 / 6.79)

    def test_osa_does_not_touch_secretion_ke_or_gb(self):
        # Neither Harsch nor Smurra reports a fasting-glucose or
        # secretion-capacity number -- an honest partial implementation.
        # Asserts the FULL dict, same discriminating-test discipline as
        # CKD/PCOS/acromegaly above.
        m = disease_modifiers([DiseaseState.OSA])
        assert m == {"p3": pytest.approx(5.75 / 6.79)}

    def test_nafld_reduces_p3_matching_bugianesi_2005_low_dose_step(self):
        # Bugianesi et al. 2005 (PMID:15747110): peripheral glucose
        # disposal decreased 30% at the LOW clamp insulin dose -> 0.70.
        # The 45% high-dose figure (0.55) is deliberately NOT the anchor.
        m = disease_modifiers([DiseaseState.NAFLD])
        assert m["p3"] == pytest.approx(0.70)
        assert m["p3"] != pytest.approx(0.55)

    def test_nafld_does_not_touch_gb_ke_secretion_or_hepatic_production(self):
        # The source directly measured endogenous glucose production as
        # NORMAL and reports no Gb/ke/secretion magnitude -- full-dict
        # assertion, same discriminating-test discipline as CKD/PCOS/OSA.
        m = disease_modifiers([DiseaseState.NAFLD])
        assert m == {"p3": pytest.approx(0.70)}


# ---------------------------------------------------------------------------
# Chronobiology
# ---------------------------------------------------------------------------

class TestChronobiology:
    def test_peak_at_10am(self):
        f = circadian_factor(10.0, amplitude=0.15)
        assert f == pytest.approx(1.15)

    def test_trough_at_22pm(self):
        f = circadian_factor(22.0, amplitude=0.15)
        assert f == pytest.approx(0.85)

    def test_morning_more_sensitive_than_evening(self):
        m_morning = circadian_modifiers(10.0)
        m_evening = circadian_modifiers(22.0)
        assert m_morning["p3"] > m_evening["p3"]

    def test_none_returns_empty(self):
        assert circadian_modifiers(None) == {}

    def test_sensitivity_percent_range(self):
        for h in range(24):
            pct = circadian_sensitivity_percent(float(h))
            assert 70 <= pct <= 130  # within ±30% of mean


# ---------------------------------------------------------------------------
# Drug interactions
# ---------------------------------------------------------------------------

class TestDrugInteractions:
    def test_caffeine_reduces_p3_at_peak(self):
        dose = DrugDose(Drug.CAFFEINE, dose_mg=200.0, hours_since_ingestion=0.75)
        m = drug_modifiers([dose])
        assert m["p3"] < 1.0

    def test_caffeine_no_effect_before_ingestion(self):
        dose = DrugDose(Drug.CAFFEINE, dose_mg=200.0, hours_since_ingestion=0.0)
        m = drug_modifiers([dose])
        # At t=0, caffeine hasn't been absorbed yet
        assert m.get("p3", 1.0) == pytest.approx(1.0, abs=0.01)

    def test_caffeine_higher_dose_larger_effect(self):
        low  = DrugDose(Drug.CAFFEINE, dose_mg=100.0,  hours_since_ingestion=1.0)
        high = DrugDose(Drug.CAFFEINE, dose_mg=400.0,  hours_since_ingestion=1.0)
        m_low  = drug_modifiers([low])
        m_high = drug_modifiers([high])
        assert m_high.get("p3", 1.0) < m_low.get("p3", 1.0)

    def test_metformin_has_no_p3_effect(self):
        # Regression: this used to include p3: 1.0 + 0.28*dose_factor,
        # left as "STILL DISPUTED" pending resolution of a genuine
        # conflict between a chronic-T2D tracer/clamp study (PMID:8923861,
        # no peripheral effect) and an acute-obese-IGT study (PMID:8582131,
        # real SI increase). Resolved by checking which regime this
        # function actually represents: metformin_modifiers() has no
        # is_chronic gating -- every dose is the chronic/T2D effect, the
        # same regime as PMID:8923861 (and this platform's own metformin
        # benchmark, Hansen et al. 2024) -- not PMID:8582131's acute/
        # pre-diabetic population. Removed entirely rather than guessed.
        dose = DrugDose(Drug.METFORMIN, dose_mg=1000.0, is_chronic=True)
        m = drug_modifiers([dose])
        assert m.get("p3", 1.0) == 1.0

    def test_metformin_reduces_gb(self):
        dose = DrugDose(Drug.METFORMIN, dose_mg=1000.0, is_chronic=True)
        m = drug_modifiers([dose])
        assert m["Gb"] < 1.0

    def test_metformin_has_no_p1_effect(self):
        # Regression: this used to include p1: 1.0 + 0.10*dose_factor,
        # cited to a narrative review (PMID:8569826) that doesn't
        # actually contain this number (confirmed via PubMedFetcher --
        # empty abstract, no indexed data). Two independent
        # placebo-controlled studies specifically measuring Bergman
        # minimal-model glucose effectiveness (S_G, this model's p1
        # analog) found metformin produces no significant change:
        # PMID:8582131 ("SG were similar in both tests") and PMID:16968813
        # ("In contrast to metformin, pioglitazone improved... glucose
        # effectiveness" -- i.e. metformin did not). Removed the p1 entry
        # entirely rather than leaving an uncited number in place.
        dose = DrugDose(Drug.METFORMIN, dose_mg=1000.0, is_chronic=True)
        m = drug_modifiers([dose])
        assert m.get("p1", 1.0) == 1.0

    def test_glucocorticoid_reduces_p3(self):
        dose = DrugDose(Drug.GLUCOCORTICOID, dose_mg=20.0)
        m = drug_modifiers([dose])
        assert m["p3"] < 1.0

    def test_glibenclamide_produces_no_multiplicative_modifier(self):
        # Regression/design guard: sulfonylurea's real effect is genuinely
        # additive (glucose_independent_secretion, base value 0.0) -- the
        # multiplicative net_factors pipeline can never move a 0.0 base
        # value, so _sulfonylurea_modifiers() must return {} here (matching
        # the GLP-1/SGLT2i precedent), with the real effect applied
        # separately by PhysiologicalPersonalizer (see TestPhysiologicalPersonalizer
        # below), not via drug_modifiers()/compose_all_modifiers() at all.
        dose = DrugDose(Drug.GLIBENCLAMIDE, dose_mg=2.5)
        m = drug_modifiers([dose])
        assert m == {}

    def test_acarbose_produces_no_pathwayparameters_modifier(self):
        # Design guard: acarbose's real effect is on MealInput.kabs
        # (applied directly by SimulationRunner at meal-construction
        # time), not on any PathwayParameters field at all -- matching
        # the GLP-1/SGLT2i/sulfonylurea precedent of returning {} here.
        dose = DrugDose(Drug.ACARBOSE, dose_mg=100.0)
        m = drug_modifiers([dose])
        assert m == {}

    def test_acarbose_kabs_scale_zero_dose_is_neutral(self):
        from biosim_physiology.drug_interactions import acarbose_kabs_scale
        assert acarbose_kabs_scale(0.0) == pytest.approx(1.0)
        assert acarbose_kabs_scale(-5.0) == pytest.approx(1.0)

    def test_acarbose_kabs_scale_rejects_non_finite_dose(self):
        # Regression: found by the mandatory security-reviewer pass --
        # this function takes a raw float (unlike every other modifier
        # function in this file, which takes an already-validated
        # DrugDose/DrugEvent), so a NaN dose used to silently return 1.0
        # via min(1.0, nan) (Python's min() takes the first argument on
        # a NaN comparison), misrepresenting an invalid dose as the
        # maximum calibrated one instead of raising.
        from biosim_physiology.drug_interactions import acarbose_kabs_scale
        with pytest.raises(ValueError):
            acarbose_kabs_scale(float("nan"))
        with pytest.raises(ValueError):
            acarbose_kabs_scale(float("inf"))

    def test_acarbose_kabs_scale_matches_100mg_calibration(self):
        # Chen et al. 2016 (PMID:27686734): a 100mg dose is this
        # function's own calibration reference (see
        # ACARBOSE_100MG_KABS_SCALE's docstring).
        from biosim_physiology.drug_interactions import (
            ACARBOSE_100MG_KABS_SCALE, acarbose_kabs_scale,
        )
        assert acarbose_kabs_scale(100.0) == pytest.approx(ACARBOSE_100MG_KABS_SCALE)

    def test_acarbose_kabs_scale_lower_dose_smaller_effect(self):
        from biosim_physiology.drug_interactions import acarbose_kabs_scale
        scale_50 = acarbose_kabs_scale(50.0)
        scale_100 = acarbose_kabs_scale(100.0)
        # A "smaller effect" means a scale CLOSER to 1.0 (less reduction).
        assert scale_50 > scale_100
        assert 0.0 < scale_50 < 1.0

    def test_acarbose_kabs_scale_capped_beyond_100mg(self):
        # No acute-dose kabs data exists beyond the 100mg reference dose
        # (see acarbose_kabs_scale()'s own module-level docstring) --
        # doses above 100mg must not extrapolate further.
        from biosim_physiology.drug_interactions import acarbose_kabs_scale
        assert acarbose_kabs_scale(200.0) == pytest.approx(acarbose_kabs_scale(100.0))
        assert acarbose_kabs_scale(1000.0) == pytest.approx(acarbose_kabs_scale(100.0))

    def test_acarbose_kabs_scale_never_reaches_zero_or_negative(self):
        # MealInput.kabs must stay strictly positive (0 < kabs < 1) --
        # confirm the scale factor can never produce a zero/negative
        # kabs regardless of dose.
        from biosim_physiology.drug_interactions import acarbose_kabs_scale
        for dose in [1.0, 50.0, 100.0, 500.0, 1e6]:
            assert acarbose_kabs_scale(dose) > 0.0

    def test_miglitol_produces_no_pathwayparameters_modifier(self):
        # Same design guard as acarbose -- miglitol's real effect is on
        # MealInput.kabs, not any PathwayParameters field.
        dose = DrugDose(Drug.MIGLITOL, dose_mg=50.0)
        m = drug_modifiers([dose])
        assert m == {}

    def test_miglitol_kabs_scale_zero_dose_is_neutral(self):
        from biosim_physiology.drug_interactions import miglitol_kabs_scale
        assert miglitol_kabs_scale(0.0) == pytest.approx(1.0)
        assert miglitol_kabs_scale(-5.0) == pytest.approx(1.0)

    def test_miglitol_kabs_scale_rejects_non_finite_dose(self):
        from biosim_physiology.drug_interactions import miglitol_kabs_scale
        with pytest.raises(ValueError):
            miglitol_kabs_scale(float("nan"))
        with pytest.raises(ValueError):
            miglitol_kabs_scale(float("inf"))

    def test_miglitol_kabs_scale_matches_50mg_calibration(self):
        # Samad et al. 1988 (PMID:3042310): 50mg TID is this function's
        # own calibration reference dose (see MIGLITOL_50MG_KABS_SCALE's
        # docstring).
        from biosim_physiology.drug_interactions import (
            MIGLITOL_50MG_KABS_SCALE, miglitol_kabs_scale,
        )
        assert miglitol_kabs_scale(50.0) == pytest.approx(MIGLITOL_50MG_KABS_SCALE)

    def test_miglitol_kabs_scale_capped_beyond_50mg(self):
        # No acute-dose data exists beyond the 50mg reference dose --
        # doses above it must not extrapolate further (a step function,
        # unlike acarbose_kabs_scale()'s continuous log-dose curve --
        # see miglitol_kabs_scale()'s own docstring for why).
        from biosim_physiology.drug_interactions import miglitol_kabs_scale
        assert miglitol_kabs_scale(100.0) == pytest.approx(miglitol_kabs_scale(50.0))
        assert miglitol_kabs_scale(1000.0) == pytest.approx(miglitol_kabs_scale(50.0))

    def test_miglitol_kabs_scale_lower_dose_smaller_effect(self):
        from biosim_physiology.drug_interactions import miglitol_kabs_scale
        scale_25 = miglitol_kabs_scale(25.0)
        scale_50 = miglitol_kabs_scale(50.0)
        assert scale_25 > scale_50
        assert 0.0 < scale_25 < 1.0

    def test_miglitol_kabs_scale_never_reaches_zero_or_negative(self):
        from biosim_physiology.drug_interactions import miglitol_kabs_scale
        for dose in [1.0, 25.0, 50.0, 500.0, 1e6]:
            assert miglitol_kabs_scale(dose) > 0.0

    def test_miglitol_and_acarbose_have_independent_calibrations(self):
        # Regression guard: runner.py previously dispatched every
        # alpha-glucosidase-inhibitor event through acarbose_kabs_scale()
        # unconditionally (latent until miglitol was added to
        # ALPHA_GLUCOSIDASE_INHIBITOR_DRUGS, since that set had only one
        # member before) -- confirms the two drugs' calibrations are
        # genuinely different, so a dispatch bug reusing one for the
        # other would be numerically detectable.
        from biosim_physiology.drug_interactions import (
            acarbose_kabs_scale, miglitol_kabs_scale,
        )
        assert acarbose_kabs_scale(50.0) != pytest.approx(miglitol_kabs_scale(50.0))

    def test_glucocorticoid_elevates_gb(self):
        dose = DrugDose(Drug.GLUCOCORTICOID, dose_mg=20.0)
        m = drug_modifiers([dose])
        assert m["Gb"] > 1.0

    def test_negative_dose_raises(self):
        with pytest.raises(ValueError):
            DrugDose(Drug.CAFFEINE, dose_mg=-100.0)

    def test_nan_dose_raises(self):
        # Regression (security review, glucocorticoid PK/PD milestone --
        # see CLAUDE.md): "dose_mg < 0" does not reject NaN (NaN
        # comparisons are always False) -- fixed at the matching choke
        # point in DrugEvent.__post_init__ (biosim_api/scenario.py) too,
        # but DrugDose is independently constructable, so it needs its
        # own guard.
        with pytest.raises(ValueError, match="finite"):
            DrugDose(Drug.GLUCOCORTICOID, dose_mg=float("nan"))

    def test_nan_hours_since_ingestion_raises(self):
        with pytest.raises(ValueError, match="finite"):
            DrugDose(Drug.GLUCOCORTICOID, dose_mg=20.0, hours_since_ingestion=float("nan"))

    def test_sglt2_inhibitors_have_no_static_modifier(self):
        # SGLT2i's effect is realised entirely via the ODE-coupled
        # RenalPathway (biosim_ode.renal_pathway), composed by
        # SimulationRunner -- a static modifier here would double-count
        # exactly the way glucocorticoid's did before that bug was fixed
        # (see CLAUDE.md KNOWN BUGS). Matches _glp1_agonist_modifiers()'s
        # empty-dict precedent.
        for drug in (Drug.DAPAGLIFLOZIN, Drug.EMPAGLIFLOZIN):
            dose = DrugDose(drug, dose_mg=10.0, is_chronic=True)
            m = drug_modifiers([dose])
            assert m == {}


class TestTirzepatideGbFactor:
    """
    2026-08-12 deep-search plan, item A6: tirzepatide's fasting-Gb
    modifier, anchored to Furihata et al. 2022 (PMID:34647404,
    independently live-verified twice this session -- once for the
    week-8 vs.-placebo deltas, once specifically for the baseline
    values needed to convert them into a fraction). See
    tirzepatide_gb_factor()'s own docstring for the full derivation.
    """

    def test_zero_dose_is_neutral(self):
        from biosim_physiology.drug_interactions import tirzepatide_gb_factor
        assert tirzepatide_gb_factor(0.0) == pytest.approx(1.0)

    def test_negative_dose_raises(self):
        from biosim_physiology.drug_interactions import tirzepatide_gb_factor
        with pytest.raises(ValueError):
            tirzepatide_gb_factor(-1.0)

    def test_nan_dose_raises(self):
        from biosim_physiology.drug_interactions import tirzepatide_gb_factor
        with pytest.raises(ValueError):
            tirzepatide_gb_factor(float("nan"))

    def test_5mg_matches_furihata_2022(self):
        # 52.7/164.5 = 0.32036... reduction -> 0.6796... multiplier
        from biosim_physiology.drug_interactions import tirzepatide_gb_factor
        assert tirzepatide_gb_factor(5.0) == pytest.approx(0.680, abs=1e-3)

    def test_10mg_matches_furihata_2022(self):
        # 69.1/181.2 = 0.38134... reduction -> 0.6187... multiplier
        from biosim_physiology.drug_interactions import tirzepatide_gb_factor
        assert tirzepatide_gb_factor(10.0) == pytest.approx(0.619, abs=1e-3)

    def test_15mg_matches_furihata_2022(self):
        # 68.9/170.7 = 0.40364... reduction -> 0.5964... multiplier
        from biosim_physiology.drug_interactions import tirzepatide_gb_factor
        assert tirzepatide_gb_factor(15.0) == pytest.approx(0.596, abs=1e-3)

    def test_dose_response_is_monotonically_decreasing_but_saturating(self):
        # Real, disclosed finding: 10mg->15mg gains far less than
        # 5mg->10mg did -- confirm the interpolation preserves this
        # genuine saturation, not a linear extrapolation.
        from biosim_physiology.drug_interactions import tirzepatide_gb_factor
        f5, f10, f15 = (tirzepatide_gb_factor(d) for d in (5.0, 10.0, 15.0))
        assert f5 > f10 > f15
        assert (f5 - f10) > (f10 - f15)

    def test_dose_above_15mg_is_clamped_not_extrapolated(self):
        from biosim_physiology.drug_interactions import tirzepatide_gb_factor
        assert tirzepatide_gb_factor(20.0) == tirzepatide_gb_factor(15.0)
        assert tirzepatide_gb_factor(100.0) == tirzepatide_gb_factor(15.0)

    def test_intermediate_dose_interpolates_between_anchors(self):
        from biosim_physiology.drug_interactions import tirzepatide_gb_factor
        f5, f7_5, f10 = (tirzepatide_gb_factor(d) for d in (5.0, 7.5, 10.0))
        assert f5 > f7_5 > f10

    def test_low_dose_ramps_linearly_from_neutral(self):
        from biosim_physiology.drug_interactions import tirzepatide_gb_factor
        f0, f2_5, f5 = (tirzepatide_gb_factor(d) for d in (0.0, 2.5, 5.0))
        assert f0 > f2_5 > f5
        assert f2_5 == pytest.approx((f0 + f5) / 2.0)

    def test_no_static_modifier_from_other_glp1_class_drugs(self):
        # Confirms tirzepatide's dispatch didn't accidentally give
        # semaglutide/liraglutide/dpp4_inhibitor a static effect too --
        # their real mechanism stays exclusively in GLP1Pathway.
        for drug in (Drug.SEMAGLUTIDE, Drug.LIRAGLUTIDE, Drug.DPP4_INHIBITOR):
            dose = DrugDose(drug, dose_mg=1.0)
            assert drug_modifiers([dose]) == {}

    def test_zero_dose_tirzepatide_event_produces_no_modifier(self):
        dose = DrugDose(Drug.TIRZEPATIDE, dose_mg=0.0)
        assert drug_modifiers([dose]) == {}

    def test_nonzero_dose_tirzepatide_lowers_gb(self):
        dose = DrugDose(Drug.TIRZEPATIDE, dose_mg=5.0)
        m = drug_modifiers([dose])
        assert m.get("Gb", 1.0) < 1.0

    def test_tirzepatide_not_in_glp1_drugs(self):
        # Deliberate scope decision -- see GLP1_DRUGS's own docstring.
        from biosim_physiology.drug_interactions import GLP1_DRUGS
        assert Drug.TIRZEPATIDE not in GLP1_DRUGS


class TestPioglitazoneP3Factor:
    """
    Fresh Track 1 scoping (2026-08-17/18), row 1.25: pioglitazone's
    insulin-sensitizing p3 modifier, anchored to Miyazaki et al. 2002
    (PMID:11874940, independently live-verified by the orchestrating
    session before implementation). See pioglitazone_p3_factor()'s own
    docstring for the full derivation.
    """

    def test_zero_dose_is_neutral(self):
        from biosim_physiology.drug_interactions import pioglitazone_p3_factor
        assert pioglitazone_p3_factor(0.0) == pytest.approx(1.0)

    def test_negative_dose_raises(self):
        from biosim_physiology.drug_interactions import pioglitazone_p3_factor
        with pytest.raises(ValueError):
            pioglitazone_p3_factor(-1.0)

    def test_nan_dose_raises(self):
        from biosim_physiology.drug_interactions import pioglitazone_p3_factor
        with pytest.raises(ValueError):
            pioglitazone_p3_factor(float("nan"))

    def test_7_5mg_and_15mg_are_empirically_neutral(self):
        # Miyazaki 2002's own confirmed null result -- these are NOT
        # small nonzero effects, they are exactly neutral, matching
        # "no change in FPG and FPI... or PG and insulin during OGTT."
        from biosim_physiology.drug_interactions import pioglitazone_p3_factor
        assert pioglitazone_p3_factor(7.5) == pytest.approx(1.0)
        assert pioglitazone_p3_factor(15.0) == pytest.approx(1.0)

    def test_30mg_matches_miyazaki_2002_whole_body_isi_ratio(self):
        # Whole-body ISI: 1.8 -> 2.5 at 30mg/day -> ratio 2.5/1.8.
        from biosim_physiology.drug_interactions import pioglitazone_p3_factor
        assert pioglitazone_p3_factor(30.0) == pytest.approx(2.5 / 1.8, abs=1e-3)

    def test_45mg_matches_miyazaki_2002_whole_body_isi_ratio(self):
        # Whole-body ISI: 1.6 -> 2.7 at 45mg/day -> ratio 2.7/1.6.
        from biosim_physiology.drug_interactions import pioglitazone_p3_factor
        assert pioglitazone_p3_factor(45.0) == pytest.approx(2.7 / 1.6, abs=1e-3)

    def test_dose_response_increases_only_above_the_null_zone(self):
        from biosim_physiology.drug_interactions import pioglitazone_p3_factor
        f0, f7_5, f15, f30, f45 = (
            pioglitazone_p3_factor(d) for d in (0.0, 7.5, 15.0, 30.0, 45.0)
        )
        assert f0 == f7_5 == f15
        assert f15 < f30 < f45

    def test_dose_above_45mg_is_clamped_not_extrapolated(self):
        from biosim_physiology.drug_interactions import pioglitazone_p3_factor
        assert pioglitazone_p3_factor(60.0) == pioglitazone_p3_factor(45.0)
        assert pioglitazone_p3_factor(200.0) == pioglitazone_p3_factor(45.0)

    def test_intermediate_dose_interpolates_between_anchors(self):
        from biosim_physiology.drug_interactions import pioglitazone_p3_factor
        f15, f22_5, f30 = (pioglitazone_p3_factor(d) for d in (15.0, 22.5, 30.0))
        assert f15 < f22_5 < f30

    def test_zero_dose_pioglitazone_event_produces_no_modifier(self):
        dose = DrugDose(Drug.PIOGLITAZONE, dose_mg=0.0)
        assert drug_modifiers([dose]) == {}

    def test_low_dose_pioglitazone_event_produces_a_neutral_modifier(self):
        # Genuinely null at 7.5/15mg, not just "small" -- unlike
        # dose_mg<=0 (which short-circuits to {} in pioglitazone_
        # modifiers()), a positive dose in the null zone still returns
        # an explicit {"p3": 1.0}, since pioglitazone_p3_factor() itself
        # is exactly 1.0 there. Asserts the exact dict rather than a
        # permissive .get() default, per a mandatory code-reviewer
        # finding that the prior version of this test would have passed
        # identically whether the key was present-and-1.0 or absent.
        dose = DrugDose(Drug.PIOGLITAZONE, dose_mg=7.5)
        m = drug_modifiers([dose])
        assert m == {"p3": pytest.approx(1.0)}

    def test_nonzero_effective_dose_increases_p3(self):
        dose = DrugDose(Drug.PIOGLITAZONE, dose_mg=30.0)
        m = drug_modifiers([dose])
        assert m["p3"] > 1.0

    def test_pioglitazone_only_modifies_p3(self):
        # Deliberate scope decision (see pioglitazone_p3_factor()'s own
        # docstring) -- the hepatic-ISI/insulinogenic-index pieces of
        # the same paper are NOT stacked here, to avoid KNOWN BUG #14's
        # over-determination mistake.
        dose = DrugDose(Drug.PIOGLITAZONE, dose_mg=45.0)
        m = drug_modifiers([dose])
        assert set(m.keys()) == {"p3"}


# ---------------------------------------------------------------------------
# Indinavir (HIV protease inhibitor) p3 factor -- Track 1 row 1.33
# ---------------------------------------------------------------------------

class TestIndinavirP3Factor:
    """
    Track 1 row 1.33: acute indinavir p3 reduction, anchored to Noor
    et al. 2002 (PMID:11964551, independently live-verified before
    implementation). See indinavir_p3_factor()'s own docstring for the
    full derivation.
    """

    def test_zero_dose_is_neutral(self):
        from biosim_physiology.drug_interactions import indinavir_p3_factor
        assert indinavir_p3_factor(0.0) == pytest.approx(1.0)

    def test_negative_dose_raises(self):
        from biosim_physiology.drug_interactions import indinavir_p3_factor
        with pytest.raises(ValueError):
            indinavir_p3_factor(-1.0)

    def test_nan_dose_raises(self):
        from biosim_physiology.drug_interactions import indinavir_p3_factor
        with pytest.raises(ValueError):
            indinavir_p3_factor(float("nan"))

    def test_tested_dose_matches_noor_2002_paired_statistic(self):
        from biosim_physiology.drug_interactions import indinavir_p3_factor
        assert indinavir_p3_factor(1200.0) == pytest.approx(1.0 - 0.341)

    def test_above_tested_dose_is_clamped_not_extrapolated(self):
        from biosim_physiology.drug_interactions import indinavir_p3_factor
        assert indinavir_p3_factor(2400.0) == indinavir_p3_factor(1200.0)
        assert indinavir_p3_factor(10000.0) == indinavir_p3_factor(1200.0)

    def test_intermediate_dose_interpolates_between_anchors(self):
        from biosim_physiology.drug_interactions import indinavir_p3_factor
        f300, f600, f900 = (indinavir_p3_factor(d) for d in (300.0, 600.0, 900.0))
        assert 1.0 > f300 > f600 > f900 > indinavir_p3_factor(1200.0)

    def test_zero_dose_event_produces_no_modifier(self):
        dose = DrugDose(Drug.INDINAVIR, dose_mg=0.0)
        assert drug_modifiers([dose]) == {}

    def test_tested_dose_reduces_p3(self):
        dose = DrugDose(Drug.INDINAVIR, dose_mg=1200.0)
        m = drug_modifiers([dose])
        assert m["p3"] < 1.0

    def test_only_modifies_p3(self):
        # secretion_scale deliberately NOT touched -- the source's own
        # directly-measured NEFA null rules out a NEFAPathway route, and
        # the source reports no secretion-magnitude finding at all.
        dose = DrugDose(Drug.INDINAVIR, dose_mg=1200.0)
        m = drug_modifiers([dose])
        assert set(m.keys()) == {"p3"}


# ---------------------------------------------------------------------------
# Hydrochlorothiazide (thiazide diuretic) p3 factor -- Track 1 row 1.34
# ---------------------------------------------------------------------------

class TestHydrochlorothiazideP3Factor:
    """
    Track 1 row 1.34: chronic hydrochlorothiazide p3 reduction, anchored
    to Pollare, Lithell & Berne 1989 (PMID:2671740, independently
    live-verified before implementation). See
    hydrochlorothiazide_p3_factor()'s own docstring for the full
    derivation.
    """

    def test_zero_dose_is_neutral(self):
        from biosim_physiology.drug_interactions import hydrochlorothiazide_p3_factor
        assert hydrochlorothiazide_p3_factor(0.0) == pytest.approx(1.0)

    def test_negative_dose_raises(self):
        from biosim_physiology.drug_interactions import hydrochlorothiazide_p3_factor
        with pytest.raises(ValueError):
            hydrochlorothiazide_p3_factor(-1.0)

    def test_nan_dose_raises(self):
        from biosim_physiology.drug_interactions import hydrochlorothiazide_p3_factor
        with pytest.raises(ValueError):
            hydrochlorothiazide_p3_factor(float("nan"))

    def test_tested_dose_matches_pollare_1989(self):
        from biosim_physiology.drug_interactions import hydrochlorothiazide_p3_factor
        assert hydrochlorothiazide_p3_factor(40.0) == pytest.approx(5.7 / 6.4)

    def test_above_tested_dose_is_clamped_not_extrapolated(self):
        from biosim_physiology.drug_interactions import hydrochlorothiazide_p3_factor
        assert hydrochlorothiazide_p3_factor(80.0) == hydrochlorothiazide_p3_factor(40.0)
        assert hydrochlorothiazide_p3_factor(1000.0) == hydrochlorothiazide_p3_factor(40.0)

    def test_intermediate_dose_interpolates_between_anchors(self):
        from biosim_physiology.drug_interactions import hydrochlorothiazide_p3_factor
        f10, f20, f30 = (hydrochlorothiazide_p3_factor(d) for d in (10.0, 20.0, 30.0))
        assert 1.0 > f10 > f20 > f30 > hydrochlorothiazide_p3_factor(40.0)

    def test_zero_dose_event_produces_no_modifier(self):
        dose = DrugDose(Drug.HYDROCHLOROTHIAZIDE, dose_mg=0.0)
        assert drug_modifiers([dose]) == {}

    def test_tested_dose_reduces_p3(self):
        dose = DrugDose(Drug.HYDROCHLOROTHIAZIDE, dose_mg=40.0)
        m = drug_modifiers([dose])
        assert m["p3"] < 1.0

    def test_only_modifies_p3(self):
        # secretion_scale deliberately NOT touched despite a real
        # positive finding in the same source ("increased the basal
        # insulin concentration and the late insulin response") --
        # treated as a compensatory response to the measured resistance,
        # not an independent mechanism (would repeat KNOWN BUG #14's
        # over-determination mistake).
        dose = DrugDose(Drug.HYDROCHLOROTHIAZIDE, dose_mg=40.0)
        m = drug_modifiers([dose])
        assert set(m.keys()) == {"p3"}


# ---------------------------------------------------------------------------
# Growth hormone (somatropin) p3 factor -- Track 1 row 1.35
# ---------------------------------------------------------------------------

class TestGrowthHormoneModifiers:
    """
    Track 1 row 1.35: acute growth hormone p3 reduction, anchored to
    Fowelin et al. 1991 (PMID:1916055, independently live-verified
    before implementation). Deliberately NOT dose-scaled -- see
    GROWTH_HORMONE_P3_FACTOR's own docstring for why (a unit-
    inconsistent source, and a saturating not linear reported shape).
    """

    def test_zero_dose_produces_no_modifier(self):
        dose = DrugDose(Drug.GROWTH_HORMONE, dose_mg=0.0)
        assert drug_modifiers([dose]) == {}

    def test_nonzero_dose_matches_fowelin_1991_peak_inhibition(self):
        dose = DrugDose(Drug.GROWTH_HORMONE, dose_mg=1.0)
        m = drug_modifiers([dose])
        assert m == {"p3": pytest.approx(1.0 - 0.39)}

    def test_different_nonzero_doses_produce_the_identical_factor(self):
        # Deliberately NOT a dose-response ramp -- a flat peak-effect
        # factor whenever GH is present at all, regardless of the
        # numeric dose_mg value (see the module docstring's own
        # unit-ambiguity/saturating-shape reasoning for why no ramp is
        # attempted). Confirms this is genuinely flat, not an
        # accidentally-still-linear function that happens to look flat
        # at small values.
        low = drug_modifiers([DrugDose(Drug.GROWTH_HORMONE, dose_mg=0.001)])
        high = drug_modifiers([DrugDose(Drug.GROWTH_HORMONE, dose_mg=1000.0)])
        assert low == high == {"p3": pytest.approx(1.0 - 0.39)}

    def test_only_modifies_p3(self):
        # secretion_scale/hepatic_egp_delta deliberately NOT touched --
        # the source's own hepatic-component mention has no separate
        # extractable magnitude (would repeat KNOWN BUG #14's
        # over-determination mistake).
        dose = DrugDose(Drug.GROWTH_HORMONE, dose_mg=1.0)
        m = drug_modifiers([dose])
        assert set(m.keys()) == {"p3"}

    def test_reduces_p3_below_neutral(self):
        dose = DrugDose(Drug.GROWTH_HORMONE, dose_mg=1.0)
        m = drug_modifiers([dose])
        assert m["p3"] < 1.0


# ---------------------------------------------------------------------------
# Body-size validation (found by the NAFLD security review, pre-existing)
# ---------------------------------------------------------------------------

class TestBodySizeValidation:
    """
    height_cm=1e-300 (REST schema only required gt=0) made profile.bmi
    divide by zero INSIDE personalize(); the runner's catch-all then dropped
    the WHOLE profile and returned a healthy-baseline 200 -- a requested
    disease state silently simulated as healthy. NaN/inf weight or height
    passed __post_init__ too (bmi_modifiers(NaN) fell through to the
    severely-obese branch). Validate at construction instead.
    """

    @pytest.mark.parametrize("kwargs", [
        {"height_cm": 1e-300},                       # bmi() underflow -> ZeroDivisionError
        {"height_cm": 1e-100},                       # finite but absurd BMI (~7e205)
        {"height_cm": 1e300},                        # (h/100)**2 overflow -> OverflowError, must still be ValueError
        {"height_cm": float("nan")},
        {"height_cm": float("inf")},
        {"body_weight_kg": float("nan")},
        {"body_weight_kg": float("inf")},
        {"body_weight_kg": 1e-9},                    # absurdly low BMI
        {"body_weight_kg": 500.0, "height_cm": 100.0},  # BMI 500
    ])
    def test_implausible_or_non_finite_body_size_is_rejected(self, kwargs):
        with pytest.raises(ValueError, match="finite and positive|implausible BMI"):
            PhysiologicalProfile(**kwargs)

    def test_bmi_bounds_are_inclusive_and_just_outside_is_rejected(self):
        lo, hi = PhysiologicalProfile.MIN_PLAUSIBLE_BMI, PhysiologicalProfile.MAX_PLAUSIBLE_BMI
        assert math.isfinite(PhysiologicalProfile(body_weight_kg=lo * 1.7 ** 2, height_cm=170.0).bmi)
        assert math.isfinite(PhysiologicalProfile(body_weight_kg=hi * 1.7 ** 2, height_cm=170.0).bmi)
        with pytest.raises(ValueError, match="implausible BMI"):
            PhysiologicalProfile(body_weight_kg=(lo - 0.01) * 1.7 ** 2, height_cm=170.0)
        with pytest.raises(ValueError, match="implausible BMI"):
            PhysiologicalProfile(body_weight_kg=(hi + 0.01) * 1.7 ** 2, height_cm=170.0)

    def test_bmi_bound_constants_are_not_dataclass_fields(self):
        import dataclasses
        names = {f.name for f in dataclasses.fields(PhysiologicalProfile)}
        assert "MIN_PLAUSIBLE_BMI" not in names and "MAX_PLAUSIBLE_BMI" not in names
        assert PhysiologicalProfile() == PhysiologicalProfile()

    @pytest.mark.parametrize("weight,height", [
        (70.0, 170.0), (45.0, 170.0), (120.0, 170.0),  # normal / underweight / severe obesity
        (200.0, 170.0),                                 # BMI ~69, extreme but real
        (40.0, 150.0),
        (14.5, 170.0), (720.0, 170.0),                  # near the BMI 5 / 250 edges
        (3.5, 50.0),                                    # newborn-sized, BMI 14
    ])
    def test_real_human_body_sizes_are_accepted(self, weight, height):
        p = PhysiologicalProfile(body_weight_kg=weight, height_cm=height)
        assert math.isfinite(p.bmi)


# ---------------------------------------------------------------------------
# Drug-combination validation (shared by ScenarioDefinition + personalizer)
# ---------------------------------------------------------------------------

class TestValidateDrugCombination:
    """
    Found by the HCQ security review (HIGH) and the NAFLD/guard code
    review (MEDIUM): duplicate events of a FLAT, presence-gated drug each
    multiply the full factor (20 HCQ events hit the 1e-4 p3 ceiling), and
    a guard living only in personalize() is swallowed by the runner's
    try/except when a profile is present. One shared validator, enforced
    at scenario construction.
    """

    def test_duplicate_hcq_events_are_rejected(self):
        from biosim_physiology.drug_interactions import validate_drug_combination
        with pytest.raises(ValueError, match="hydroxychloroquine"):
            validate_drug_combination([
                DrugDose(Drug.HYDROXYCHLOROQUINE, dose_mg=200.0),
                DrugDose(Drug.HYDROXYCHLOROQUINE, dose_mg=200.0),
            ])

    def test_duplicate_growth_hormone_events_are_rejected(self):
        from biosim_physiology.drug_interactions import validate_drug_combination
        with pytest.raises(ValueError, match="growth_hormone"):
            validate_drug_combination([
                DrugDose(Drug.GROWTH_HORMONE, dose_mg=1.0),
                DrugDose(Drug.GROWTH_HORMONE, dose_mg=1.0),
            ])

    def test_single_hcq_and_zero_dose_duplicate_are_accepted(self):
        from biosim_physiology.drug_interactions import validate_drug_combination
        validate_drug_combination([DrugDose(Drug.HYDROXYCHLOROQUINE, dose_mg=400.0)])
        validate_drug_combination([
            DrugDose(Drug.HYDROXYCHLOROQUINE, dose_mg=400.0),
            DrugDose(Drug.HYDROXYCHLOROQUINE, dose_mg=0.0),  # applies no effect
        ])

    def test_distinct_sulfonylureas_are_rejected_but_same_one_twice_is_not(self):
        from biosim_physiology.drug_interactions import validate_drug_combination
        with pytest.raises(ValueError, match="sulfonylurea"):
            validate_drug_combination([
                DrugDose(Drug.GLIBENCLAMIDE, dose_mg=2.5),
                DrugDose(Drug.GLIPIZIDE, dose_mg=5.0),
            ])
        validate_drug_combination([
            DrugDose(Drug.GLIBENCLAMIDE, dose_mg=2.5),
            DrugDose(Drug.GLIBENCLAMIDE, dose_mg=2.5),
        ])

    @pytest.mark.parametrize("drug,dose_mg", [
        (Drug.METFORMIN, 1000.0),
        (Drug.PIOGLITAZONE, 15.0),
        (Drug.INDINAVIR, 600.0),
        (Drug.HYDROCHLOROTHIAZIDE, 20.0),
        (Drug.ATORVASTATIN, 10.0),
        (Drug.ACARBOSE, 25.0),
        (Drug.MIGLITOL, 25.0),
        (Drug.SEMAGLUTIDE, 0.5),
        (Drug.LIRAGLUTIDE, 0.6),
    ])
    def test_duplicate_dose_scaled_events_are_rejected(self, drug, dose_mg):
        # Static factors compound per event (2 x pioglitazone 15mg != 30mg);
        # first-event-only drugs silently drop later events, making the
        # result order-dependent. One event with the total dose is the only
        # faithful representation.
        from biosim_physiology.drug_interactions import validate_drug_combination
        with pytest.raises(ValueError, match=f"{drug.value}.*total dose"):
            validate_drug_combination([
                DrugDose(drug, dose_mg=dose_mg),
                DrugDose(drug, dose_mg=dose_mg),
            ])

    def test_duplicate_dose_scaled_events_would_otherwise_compound(self):
        # Documents why the rejection exists: without it, two 15mg
        # pioglitazone events produce a larger p3 factor than one 30mg event.
        from biosim_physiology.drug_interactions import drug_modifiers
        twice = drug_modifiers([DrugDose(Drug.PIOGLITAZONE, 15.0)] * 2)["p3"]
        once = drug_modifiers([DrugDose(Drug.PIOGLITAZONE, 30.0)])["p3"]
        assert twice != pytest.approx(once)

    def test_single_dose_scaled_event_and_zero_dose_duplicate_are_accepted(self):
        from biosim_physiology.drug_interactions import validate_drug_combination
        validate_drug_combination([
            DrugDose(Drug.METFORMIN, dose_mg=2000.0),
            DrugDose(Drug.METFORMIN, dose_mg=0.0),  # applies no effect
            DrugDose(Drug.PIOGLITAZONE, dose_mg=30.0),
            DrugDose(Drug.ATORVASTATIN, dose_mg=40.0),
        ])

    def test_dose_independent_drugs_may_repeat(self):
        # Their modelled effect does not depend on dose, so a duplicate is
        # inert -- rejecting it would only break harmless requests.
        from biosim_physiology.drug_interactions import validate_drug_combination
        for drug in (Drug.DPP4_INHIBITOR, Drug.TIRZEPATIDE, Drug.DAPAGLIFLOZIN):
            validate_drug_combination([DrugDose(drug, 5.0), DrugDose(drug, 10.0)])

    def test_legitimately_repeatable_drugs_are_not_restricted(self):
        # Exogenous insulin / caffeine are time-resolved, repeatable
        # events -- the validator must not touch them.
        from biosim_physiology.drug_interactions import validate_drug_combination
        validate_drug_combination([
            DrugDose(Drug.INSULIN_REGULAR, dose_mg=5.0),
            DrugDose(Drug.INSULIN_REGULAR, dose_mg=5.0),
            DrugDose(Drug.CAFFEINE, dose_mg=100.0),
            DrugDose(Drug.CAFFEINE, dose_mg=100.0),
        ])

    def test_personalizer_still_enforces_it_directly(self):
        # Defence in depth for direct (non-scenario) callers.
        personalizer = PhysiologicalPersonalizer()
        with pytest.raises(ValueError, match="hydroxychloroquine"):
            personalizer.personalize(
                PhysiologicalProfile(),
                drug_doses=[
                    DrugDose(Drug.HYDROXYCHLOROQUINE, dose_mg=200.0),
                    DrugDose(Drug.HYDROXYCHLOROQUINE, dose_mg=200.0),
                ],
            )


# ---------------------------------------------------------------------------
# Hydroxychloroquine p3 factor -- Track 1 row 1.44
# ---------------------------------------------------------------------------

class TestHydroxychloroquineModifiers:
    """
    Track 1 row 1.44: a genuinely new drug class (antimalarial/DMARD),
    found via a fresh Track 1 scoping fan-out. Anchored to Wasko et al.
    2015 (PMID:26197707, independently live-verified before
    implementation) -- a within-arm minimal-model-SI increase (+20.0%),
    the conservative choice over the larger, placebo-referenced 38.3%
    difference (see HYDROXYCHLOROQUINE_P3_FACTOR's own docstring for
    the full anchor-choice reasoning). Deliberately NOT dose-scaled --
    a single tested dose (400mg/day), same design as growth hormone
    above.
    """

    def test_zero_dose_produces_no_modifier(self):
        dose = DrugDose(Drug.HYDROXYCHLOROQUINE, dose_mg=0.0)
        assert drug_modifiers([dose]) == {}

    def test_nonzero_dose_matches_wasko_2015_within_arm_increase(self):
        dose = DrugDose(Drug.HYDROXYCHLOROQUINE, dose_mg=400.0)
        m = drug_modifiers([dose])
        assert m == {"p3": pytest.approx(1.0 + 0.20)}

    def test_different_nonzero_doses_produce_the_identical_factor(self):
        # Deliberately NOT a dose-response ramp -- a flat peak-effect
        # factor whenever HCQ is present at all, regardless of the
        # numeric dose_mg value (a single tested dose, no dose-response
        # curve in either source).
        low = drug_modifiers([DrugDose(Drug.HYDROXYCHLOROQUINE, dose_mg=0.001)])
        high = drug_modifiers([DrugDose(Drug.HYDROXYCHLOROQUINE, dose_mg=1000.0)])
        assert low == high == {"p3": pytest.approx(1.0 + 0.20)}

    def test_only_modifies_p3(self):
        # hepatic_egp_delta/ke/secretion_scale deliberately NOT touched
        # -- Toledo et al. 2021 (PMID:33528893) directly measured both
        # hepatic insulin sensitivity and insulin clearance as NULL, and
        # secretion_scale is deliberately not modeled from the same
        # dataset's own beta-cell-function finding (would repeat KNOWN
        # BUG #14's over-determination mistake).
        dose = DrugDose(Drug.HYDROXYCHLOROQUINE, dose_mg=400.0)
        m = drug_modifiers([dose])
        assert set(m.keys()) == {"p3"}

    def test_raises_p3_above_neutral(self):
        # Unlike every other drug in this file (all p3-lowering), HCQ
        # is p3-RAISING (a genuine insulin-sensitizing effect) --
        # confirms the direction is correct, not accidentally inverted.
        dose = DrugDose(Drug.HYDROXYCHLOROQUINE, dose_mg=400.0)
        m = drug_modifiers([dose])
        assert m["p3"] > 1.0


# ---------------------------------------------------------------------------
# Atorvastatin p3 factor -- Track 1 row 1.41
# ---------------------------------------------------------------------------

class TestAtorvastatinP3Factor:
    """
    Track 1 row 1.41: atorvastatin p3 reduction, anchored to Koh et al.
    2010 (PMID:20298928, independently live-verified before
    implementation -- flagged informally in an earlier session, never
    re-checked or given its own row until now). See
    atorvastatin_p3_factor()'s own docstring for the full derivation.
    """

    def test_zero_dose_is_neutral(self):
        from biosim_physiology.drug_interactions import atorvastatin_p3_factor
        assert atorvastatin_p3_factor(0.0) == pytest.approx(1.0)

    def test_negative_dose_raises(self):
        from biosim_physiology.drug_interactions import atorvastatin_p3_factor
        with pytest.raises(ValueError):
            atorvastatin_p3_factor(-1.0)

    def test_nan_dose_raises(self):
        from biosim_physiology.drug_interactions import atorvastatin_p3_factor
        with pytest.raises(ValueError):
            atorvastatin_p3_factor(float("nan"))

    def test_ten_mg_is_empirically_null_not_a_small_effect(self):
        # Koh et al. 2010's own 10mg arm: 1% decrease vs. baseline, but
        # NOT statistically significant (p=0.312) -- genuinely null,
        # not just small, matching pioglitazone_p3_factor()'s own
        # established null-zone convention for a non-significant
        # low-dose point from the same style of source.
        from biosim_physiology.drug_interactions import atorvastatin_p3_factor
        assert atorvastatin_p3_factor(10.0) == pytest.approx(1.0)

    def test_twenty_mg_matches_koh_2010(self):
        from biosim_physiology.drug_interactions import atorvastatin_p3_factor
        assert atorvastatin_p3_factor(20.0) == pytest.approx(0.97)

    def test_eighty_mg_matches_koh_2010(self):
        from biosim_physiology.drug_interactions import atorvastatin_p3_factor
        assert atorvastatin_p3_factor(80.0) == pytest.approx(0.96)

    def test_above_eighty_mg_is_clamped_not_extrapolated(self):
        from biosim_physiology.drug_interactions import atorvastatin_p3_factor
        assert atorvastatin_p3_factor(160.0) == atorvastatin_p3_factor(80.0)
        assert atorvastatin_p3_factor(1000.0) == atorvastatin_p3_factor(80.0)

    def test_intermediate_dose_interpolates_between_anchors(self):
        from biosim_physiology.drug_interactions import atorvastatin_p3_factor
        f10, f15, f20 = (atorvastatin_p3_factor(d) for d in (10.0, 15.0, 20.0))
        assert f10 > f15 > f20

    def test_zero_dose_event_produces_no_modifier(self):
        dose = DrugDose(Drug.ATORVASTATIN, dose_mg=0.0)
        assert drug_modifiers([dose]) == {}

    def test_ten_mg_event_produces_an_explicit_neutral_modifier(self):
        # Genuinely null at 10mg, not just "small" -- a positive dose in
        # the null zone still returns an explicit {"p3": 1.0}, matching
        # the same code-reviewer-driven precedent already established
        # for pioglitazone's own low-dose null zone (asserting the
        # exact dict, not a permissive .get() default).
        dose = DrugDose(Drug.ATORVASTATIN, dose_mg=10.0)
        m = drug_modifiers([dose])
        assert m == {"p3": pytest.approx(1.0)}

    def test_effective_dose_reduces_p3(self):
        dose = DrugDose(Drug.ATORVASTATIN, dose_mg=40.0)
        m = drug_modifiers([dose])
        assert m["p3"] < 1.0

    def test_only_modifies_p3(self):
        # secretion_scale deliberately NOT touched despite a real,
        # dose-dependent positive finding in the same source
        # ("significantly increased fasting plasma insulin") --
        # compensatory response, not an independent mechanism (would
        # repeat KNOWN BUG #14's over-determination mistake).
        dose = DrugDose(Drug.ATORVASTATIN, dose_mg=40.0)
        m = drug_modifiers([dose])
        assert set(m.keys()) == {"p3"}


# ---------------------------------------------------------------------------
# PhysiologicalPersonalizer — full pipeline
# ---------------------------------------------------------------------------

class TestPhysiologicalPersonalizer:
    def test_default_profile_minimal_change(self):
        p = PhysiologicalPersonalizer()
        result = p.personalize(PhysiologicalProfile())
        # Default profile = reference: most modifiers should be 1.0
        assert result.n_parameters_modified <= 2

    def test_older_diabetic_reduces_sensitivity(self):
        personalizer = PhysiologicalPersonalizer()
        profile = PhysiologicalProfile(
            age_years=65,
            disease_states=[DiseaseState.TYPE_2_DIABETES],
        )
        result = personalizer.personalize(profile)
        base_p3 = PathwayParameters().p3
        assert result.personalized_params.p3 < base_p3

    def test_applied_modifiers_bmi_lists_every_parameter_it_actually_touches(self):
        # Regression: AppliedModifier used to hardcode parameter="p3" for
        # BMI regardless of which parameters bmi_modifiers() actually
        # returns -- for overweight/obese/severely-obese BMI, that's both
        # p3 AND Gb (see modifiers.py), so the old single-parameter shape
        # silently dropped Gb from the record entirely.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(PhysiologicalProfile(body_weight_kg=100.0, height_cm=170.0))
        bmi_mod = next(m for m in result.applied_modifiers if m.source.startswith("bmi_"))
        assert set(bmi_mod.parameters.keys()) == {"p3", "Gb"}

    def test_applied_modifiers_t2d_lists_real_parameters_not_the_word_multiple(self):
        # Regression: disease modifiers used to be recorded as
        # parameter="multiple", factor=1.0 -- a placeholder that named
        # nothing and always claimed a neutral factor even though T2D
        # genuinely modifies secretion_scale/p3/Gb.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES])
        )
        disease_mod = next(m for m in result.applied_modifiers if m.source.startswith("disease_"))
        assert set(disease_mod.parameters.keys()) == {"secretion_scale", "p3", "Gb", "p1"}
        assert disease_mod.parameters["Gb"] > 1.0  # T2D elevates fasting glucose
        assert disease_mod.parameters["p3"] < 1.0  # T2D reduces insulin sensitivity

    def test_applied_modifiers_metformin_lists_gb_only(self):
        # Regression: every drug's AppliedModifier used to report
        # factor=net_factors.get("p3", 1.0) regardless of which
        # parameters that specific drug modifies -- metformin's real Gb
        # effect (now -14.7%, see CLAUDE.md) was invisible here entirely.
        # Metformin's p3 effect was itself later removed (see
        # test_metformin_has_no_p3_effect) -- Gb is now the only
        # parameter metformin_modifiers() actually touches.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(),
            drug_doses=[DrugDose(Drug.METFORMIN, dose_mg=2000.0, is_chronic=True)],
        )
        drug_mod = next(m for m in result.applied_modifiers if m.source.startswith("drug_metformin"))
        assert set(drug_mod.parameters.keys()) == {"Gb"}
        assert drug_mod.parameters["Gb"] < 1.0

    def test_glibenclamide_sets_glucose_independent_secretion(self):
        # The core wiring this milestone adds: since drug_modifiers()/
        # compose_all_modifiers() structurally cannot express an additive
        # effect on a 0.0 base value, PhysiologicalPersonalizer must apply
        # it directly (bypassing net_factors), the same special-cased
        # pattern already used for the Gb/G_threshold margin clamp and
        # Brunzell's first_phase_gain zeroing.
        from biosim_physiology.drug_interactions import glucose_independent_secretion_pM_per_min
        personalizer = PhysiologicalPersonalizer()
        base = PathwayParameters()
        result = personalizer.personalize(
            PhysiologicalProfile(),
            drug_doses=[DrugDose(Drug.GLIBENCLAMIDE, dose_mg=2.5)],
        )
        expected = glucose_independent_secretion_pM_per_min(base.ke, base.Ib, Drug.GLIBENCLAMIDE)
        assert result.personalized_params.glucose_independent_secretion == pytest.approx(expected)
        assert expected > 0.0

    def test_no_sulfonylurea_leaves_glucose_independent_secretion_at_zero(self):
        # Regression safety: every other drug/profile combination must
        # leave the new field at its neutral default.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES]),
            drug_doses=[DrugDose(Drug.METFORMIN, dose_mg=2000.0, is_chronic=True)],
        )
        assert result.personalized_params.glucose_independent_secretion == pytest.approx(0.0)

    def test_glibenclamide_zero_dose_applies_no_effect(self):
        # Regression: found by the mandatory security-reviewer pass -- the
        # special-cased override block only checked drug identity, not
        # dose_mg, so a dose_mg=0.0 event previously still triggered the
        # full effect (every other dose-scaled modifier function in this
        # file already guards `dose.dose_mg <= 0: return {}`; this
        # bypasses those functions entirely so needed its own copy of the
        # same guard).
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(),
            drug_doses=[DrugDose(Drug.GLIBENCLAMIDE, dose_mg=0.0)],
        )
        assert result.personalized_params.glucose_independent_secretion == pytest.approx(0.0)
        assert not any(m.source.startswith("drug_glibenclamide") for m in result.applied_modifiers)

    def test_glibenclamide_applies_no_effect_for_type_1_diabetes(self):
        # Regression: found by the mandatory security-reviewer pass on
        # the T1D milestone -- glucose_independent_secretion_pM_per_min()
        # depends only on ke/Ib, not secretion_scale (deliberately, since
        # the KATP-channel mechanism is glucose-INDEPENDENT, not
        # beta-cell-mass-independent), so T1D's near-zero secretion_scale
        # (beta cells floored at the model's 5% engineering minimum) had
        # no bearing on whether this override applied -- producing a
        # real, nonzero, ~50%-of-floor insulin rise from a drug that
        # clinically requires substantial functional beta-cell mass to
        # work at all and is not prescribed in confirmed T1D. Fixed by
        # skipping the override (and its audit-trail entry) entirely
        # when TYPE_1_DIABETES is present.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.TYPE_1_DIABETES]),
            drug_doses=[DrugDose(Drug.GLIBENCLAMIDE, dose_mg=2.5)],
        )
        assert result.personalized_params.glucose_independent_secretion == pytest.approx(0.0)
        assert not any(m.source.startswith("drug_glibenclamide") for m in result.applied_modifiers)

    def test_glibenclamide_still_applies_normally_for_type_2_diabetes(self):
        # Guard against the T1D fix above accidentally becoming a
        # blanket "any diabetes disease state" exclusion -- T2D is the
        # actual population Riefflin et al. 2015 studied, and its
        # already-calibrated effect must be completely unaffected.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES]),
            drug_doses=[DrugDose(Drug.GLIBENCLAMIDE, dose_mg=2.5)],
        )
        assert result.personalized_params.glucose_independent_secretion > 0.0
        assert any(m.source.startswith("drug_glibenclamide") for m in result.applied_modifiers)

    def test_glipizide_sets_glucose_independent_secretion(self):
        # 2026-08-12 deep-search plan, item D1: glipizide extension,
        # anchored to Groop et al. 1991 (PMID:1777653) and Groop et al.
        # 1987 (PMID:3117609) -- see GLIPIZIDE_SECRETION_FRACTION_OF_KE_IB's
        # own docstring for the full derivation.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES]),
            drug_doses=[DrugDose(Drug.GLIPIZIDE, dose_mg=5.0)],
        )
        assert result.personalized_params.glucose_independent_secretion > 0.0
        assert any(m.source.startswith("drug_glipizide") for m in result.applied_modifiers)

    def test_glipizide_effect_is_meaningfully_larger_than_glibenclamide(self):
        # Regression guard: glipizide's own fraction (0.9) must stay
        # distinctly larger than glibenclamide's (0.5), motivated by the
        # real ~2x head-to-head ratio finding (Groop et al. 1987), not
        # silently collapsed to the same shared constant.
        personalizer = PhysiologicalPersonalizer()
        glib = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES]),
            drug_doses=[DrugDose(Drug.GLIBENCLAMIDE, dose_mg=2.5)],
        )
        glip = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES]),
            drug_doses=[DrugDose(Drug.GLIPIZIDE, dose_mg=5.0)],
        )
        assert glip.personalized_params.glucose_independent_secretion \
            > glib.personalized_params.glucose_independent_secretion

    def test_glipizide_zero_dose_applies_no_effect(self):
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(),
            drug_doses=[DrugDose(Drug.GLIPIZIDE, dose_mg=0.0)],
        )
        assert result.personalized_params.glucose_independent_secretion == pytest.approx(0.0)
        assert not any(m.source.startswith("drug_glipizide") for m in result.applied_modifiers)

    def test_glipizide_applies_no_effect_for_type_1_diabetes(self):
        # Same rationale as glibenclamide's own T1D exclusion (the KATP-
        # channel mechanism requires substantial functional beta-cell mass).
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.TYPE_1_DIABETES]),
            drug_doses=[DrugDose(Drug.GLIPIZIDE, dose_mg=5.0)],
        )
        assert result.personalized_params.glucose_independent_secretion == pytest.approx(0.0)

    def test_glipizide_stays_below_the_hard_safety_ceiling(self):
        # GLIPIZIDE_SECRETION_FRACTION_OF_KE_IB=0.9 was deliberately capped
        # below 1.0 (the exact 2x-glibenclamide value would have hit
        # PathwayParameters.__post_init__'s hard safety bound). Confirms
        # this directly: constructing personalized params must not raise.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES]),
            drug_doses=[DrugDose(Drug.GLIPIZIDE, dose_mg=5.0)],
        )
        p = result.personalized_params
        assert p.glucose_independent_secretion < p.ke * p.Ib

    def test_glimepiride_sets_glucose_independent_secretion(self):
        # Resolves MASTER_PLAN.md row 1.3b's own remaining open item --
        # anchored to Draeger et al. 1996 (PMID:8911976), a head-to-head
        # ratio to glibenclamide's own already-adopted fraction, see
        # GLIMEPIRIDE_SECRETION_FRACTION_OF_KE_IB's own docstring for the
        # full derivation.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES]),
            drug_doses=[DrugDose(Drug.GLIMEPIRIDE, dose_mg=1.0)],
        )
        assert result.personalized_params.glucose_independent_secretion > 0.0
        assert any(m.source.startswith("drug_glimepiride") for m in result.applied_modifiers)

    def test_glimepiride_effect_is_smaller_than_glibenclamide(self):
        # Regression guard: glimepiride's own fraction (0.40) must stay
        # distinctly SMALLER than glibenclamide's (0.5) -- motivated by
        # Draeger et al. 1996's own directly-quoted finding that
        # glimepiride patients had LOWER fasting insulin than matched
        # glibenclamide patients at equivalent glycaemic control, not
        # silently collapsed to the same shared constant or (incorrectly)
        # larger like glipizide's own real 2x-glibenclamide finding.
        personalizer = PhysiologicalPersonalizer()
        glib = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES]),
            drug_doses=[DrugDose(Drug.GLIBENCLAMIDE, dose_mg=2.5)],
        )
        glim = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES]),
            drug_doses=[DrugDose(Drug.GLIMEPIRIDE, dose_mg=1.0)],
        )
        assert glim.personalized_params.glucose_independent_secretion \
            < glib.personalized_params.glucose_independent_secretion

    def test_glimepiride_zero_dose_applies_no_effect(self):
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(),
            drug_doses=[DrugDose(Drug.GLIMEPIRIDE, dose_mg=0.0)],
        )
        assert result.personalized_params.glucose_independent_secretion == pytest.approx(0.0)
        assert not any(m.source.startswith("drug_glimepiride") for m in result.applied_modifiers)

    def test_glimepiride_applies_no_effect_for_type_1_diabetes(self):
        # Same rationale as glibenclamide/glipizide's own T1D exclusion
        # (the KATP-channel mechanism requires substantial functional
        # beta-cell mass) -- inherited automatically via SULFONYLUREA_DRUGS
        # membership, no separate gating logic needed.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.TYPE_1_DIABETES]),
            drug_doses=[DrugDose(Drug.GLIMEPIRIDE, dose_mg=1.0)],
        )
        assert result.personalized_params.glucose_independent_secretion == pytest.approx(0.0)

    def test_glimepiride_stays_below_the_hard_safety_ceiling(self):
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES]),
            drug_doses=[DrugDose(Drug.GLIMEPIRIDE, dose_mg=1.0)],
        )
        p = result.personalized_params
        assert p.glucose_independent_secretion < p.ke * p.Ib

    def test_multiple_distinct_sulfonylureas_are_rejected(self):
        # Security-review finding (glimepiride chunk): the override loop
        # applies only the FIRST sulfonylurea dose (break on first match)
        # while the audit trail listed every one, so combining two gave an
        # order-dependent result and a misleading audit trail. Clinically
        # nonsensical anyway -- fail fast, matching this project's
        # mutual-exclusion convention.
        personalizer = PhysiologicalPersonalizer()
        with pytest.raises(ValueError, match="sulfonylurea"):
            personalizer.personalize(
                PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES]),
                drug_doses=[
                    DrugDose(Drug.GLIBENCLAMIDE, dose_mg=2.5),
                    DrugDose(Drug.GLIMEPIRIDE, dose_mg=1.0),
                ],
            )

    def test_zero_dose_second_sulfonylurea_is_not_rejected(self):
        # A zero-dose event applies no effect (dose_mg > 0 gate), so it
        # is not a real combination.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES]),
            drug_doses=[
                DrugDose(Drug.GLIBENCLAMIDE, dose_mg=2.5),
                DrugDose(Drug.GLIMEPIRIDE, dose_mg=0.0),
            ],
        )
        assert result.personalized_params.glucose_independent_secretion > 0.0

    def test_glimepiride_result_is_a_genuine_fixed_point(self):
        # Q2 checklist item 2 (security-review gap): end-to-end basal
        # fixed-point check, same as glibenclamide's own test.
        personalizer = PhysiologicalPersonalizer()
        for states in ([], [DiseaseState.TYPE_2_DIABETES]):
            result = personalizer.personalize(
                PhysiologicalProfile(disease_states=states),
                drug_doses=[DrugDose(Drug.GLIMEPIRIDE, dose_mg=1.0)],
            )
            pathway = GlucoseInsulinPathway(params=result.personalized_params)
            assert pathway.validate_y0(tolerance=1e-8) == []

    def test_glimepiride_factor_matches_the_real_module_constant(self):
        from biosim_physiology.drug_interactions import GLIMEPIRIDE_SECRETION_FRACTION_OF_KE_IB
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES]),
            drug_doses=[DrugDose(Drug.GLIMEPIRIDE, dose_mg=1.0)],
        )
        p = result.personalized_params
        assert p.glucose_independent_secretion == pytest.approx(
            GLIMEPIRIDE_SECRETION_FRACTION_OF_KE_IB * p.ke * p.Ib
        )

    def test_glibenclamide_effect_is_deliberately_timing_independent(self):
        # Design decision, not an oversight (see _sulfonylurea_modifiers()'s
        # own docstring): Riefflin et al. 2015 (PMID:25297572) measured a
        # CHRONIC, steady-state effect, so this is modeled as immediately
        # and fully active once any nonzero-dose glibenclamide event is
        # present, regardless of hours_since_ingestion/hours_before_
        # simulation -- matching metformin_modifiers()'s own established
        # "always chronic, timing irrelevant" precedent. Confirmed
        # explicitly here (rather than left as an untested assumption) so
        # a future session doesn't mistake this for an unfixed bug.
        personalizer = PhysiologicalPersonalizer()
        just_taken = personalizer.personalize(
            PhysiologicalProfile(),
            drug_doses=[DrugDose(Drug.GLIBENCLAMIDE, dose_mg=2.5, hours_since_ingestion=0.0)],
        )
        long_ago = personalizer.personalize(
            PhysiologicalProfile(),
            drug_doses=[DrugDose(Drug.GLIBENCLAMIDE, dose_mg=2.5, hours_since_ingestion=500.0)],
        )
        assert (
            just_taken.personalized_params.glucose_independent_secretion
            == pytest.approx(long_ago.personalized_params.glucose_independent_secretion)
        )
        assert just_taken.personalized_params.glucose_independent_secretion > 0.0

    def test_applied_modifiers_glibenclamide_records_the_real_additive_effect(self):
        # Regression guard for the audit-trail gap this milestone would
        # otherwise reproduce: modifiers_for_dose(glibenclamide_dose)
        # returns {} (by design), so the generic per-drug loop in
        # _describe_modifiers() would otherwise silently skip recording
        # this drug's real effect entirely (the same "AppliedModifier
        # hardcoding" bug class already fixed for BMI/T2D/metformin above).
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(),
            drug_doses=[DrugDose(Drug.GLIBENCLAMIDE, dose_mg=2.5)],
        )
        drug_mod = next(m for m in result.applied_modifiers if m.source.startswith("drug_glibenclamide"))
        assert set(drug_mod.parameters.keys()) == {"glucose_independent_secretion"}
        assert drug_mod.parameters["glucose_independent_secretion"] > 0.0

    def test_glibenclamide_result_is_a_genuine_fixed_point(self):
        # End-to-end fixed-point safety check through the real
        # personalization pipeline (not just the bare PathwayParameters
        # unit tests in test_ode_engine.py) -- confirms overrides.get("ke"/
        # "Ib") correctly reflects any prior modifier's effect on those
        # fields (none currently touch ke/Ib, but the lookup must not
        # silently use stale base values if that ever changes).
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(),
            drug_doses=[DrugDose(Drug.GLIBENCLAMIDE, dose_mg=2.5)],
        )
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        assert pathway.validate_y0(tolerance=1e-8) == []

    def test_t1d_ib_and_secretion_scale_clamp_to_engineering_floor(self):
        # disease_modifiers()'s TYPE_1_DIABETES factors (0.05 for Ib,
        # 0.02 for secretion_scale) are deliberately chosen to clamp to
        # _PARAM_BOUNDS' own floors (10.0 pM / 0.05) after multiplying
        # the population-average defaults -- verifying the actual
        # clamped result, not just the raw factor.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.TYPE_1_DIABETES])
        )
        assert result.personalized_params.Ib == pytest.approx(10.0)
        assert result.personalized_params.secretion_scale == pytest.approx(0.05)

    def test_t1d_gb_crosses_brunzell_threshold_and_zeros_first_phase_gain(self):
        # T1D's disclosed Gb=11.0 mM exceeds BRUNZELL_AIR_ABSENT_THRESHOLD_MM
        # (~6.39 mM) -- confirms the EXISTING Brunzell gate (added for T2D)
        # correctly generalizes to this new disease state with no further
        # special-casing, exactly as its own docstring claims.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.TYPE_1_DIABETES])
        )
        assert result.personalized_params.Gb == pytest.approx(11.0)
        assert result.personalized_params.first_phase_gain == 0.0

    def test_t1d_result_is_a_genuine_fixed_point(self):
        # Confirms basal_state()'s already-verified generalization
        # (CLAUDE.md "Confirmed the ODE-engine already supports
        # near-total beta-cell failure") holds for the REAL, named T1D
        # preset, not just the illustrative values used in that
        # investigation.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.TYPE_1_DIABETES])
        )
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        assert pathway.validate_y0(tolerance=1e-8) == []

    def test_gestational_diabetes_result_is_a_genuine_fixed_point(self):
        # A purely multiplicative Gb/p1/p3/secretion_scale preset (like
        # T2D/INSULIN_RESISTANT), so this should hold by the same
        # analytic-fixed-point machinery those already rely on -- worth
        # confirming explicitly rather than assuming, per this project's
        # own established discipline for every new disease-state preset.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.GESTATIONAL_DIABETES])
        )
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        assert pathway.validate_y0(tolerance=1e-8) == []

    def test_gestational_diabetes_and_insulin_resistant_composed_is_a_genuine_fixed_point(self):
        # GESTATIONAL_DIABETES is deliberately allowed to combine with
        # INSULIN_RESISTANT (a real, common comorbidity, unlike the
        # TYPE_1_DIABETES exclusion) -- found missing by the mandatory
        # code-reviewer pass: only the "alone" case had a fixed-point
        # regression test. Confirms the compounded 0.16*0.70=0.112 p3
        # factor and 1.02*1.03=1.0506 Gb factor (both computed and
        # verified against _PARAM_BOUNDS during review -- neither
        # clamps) still produce a genuine analytic fixed point, not just
        # a mathematically-valid-but-unverified one.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[
                DiseaseState.GESTATIONAL_DIABETES, DiseaseState.INSULIN_RESISTANT,
            ])
        )
        assert result.personalized_params.Gb == pytest.approx(5.0 * 1.02 * 1.03)
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        assert pathway.validate_y0(tolerance=1e-8) == []

    def test_gestational_diabetes_gb_stays_below_brunzell_threshold(self):
        # Unlike T1D/T2D (whose elevated Gb crosses BRUNZELL_AIR_ABSENT_
        # THRESHOLD_MM and zeroes first_phase_gain), GDM's modest 5.1 mM
        # Gb should NOT cross it -- consistent with GDM being a milder,
        # largely-compensated state that retains a first-phase response.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.GESTATIONAL_DIABETES])
        )
        assert result.personalized_params.first_phase_gain > 0.0

    def test_gck_mody_result_is_a_genuine_fixed_point(self):
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.GCK_MODY])
        )
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        assert pathway.validate_y0(tolerance=1e-8) == []

    def test_gck_mody_leaves_insulin_sensitivity_at_the_population_default(self):
        # End-to-end confirmation (not just the raw modifier dict) that
        # p3/p1/secretion_scale genuinely stay at their unmodified
        # defaults through the real personalizer pipeline -- the defining
        # mechanistic fact this preset represents (Guenat et al. 2000,
        # PMID:11119017: peripheral insulin action "unchanged in MODY2").
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.GCK_MODY])
        )
        default = PathwayParameters()
        assert result.personalized_params.p1 == pytest.approx(default.p1)
        assert result.personalized_params.p3 == pytest.approx(default.p3)
        assert result.personalized_params.secretion_scale == pytest.approx(default.secretion_scale)

    def test_gck_mody_gb_crosses_brunzell_threshold(self):
        # Documents the known, disclosed uncertainty from modifiers.py's
        # own citation comment: GCK-MODY's 7.17 mM Gb exceeds
        # BRUNZELL_AIR_ABSENT_THRESHOLD_MM (~6.39 mM), so the existing,
        # generic T2D-era gate fires for this preset too, even though no
        # citation confirms this generalizes correctly to GCK-MODY's own
        # threshold-reset mechanism specifically. Locking in the CURRENT
        # behaviour (not asserting it's correct) so a future session
        # revisiting this sees a real, deliberate test rather than
        # silently-untested behaviour.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.GCK_MODY])
        )
        assert result.personalized_params.first_phase_gain == 0.0

    def test_ckd_result_is_a_genuine_fixed_point(self):
        # Confirms the fixed point holds specifically for CKD's OWN
        # p3/ke combination, not just that the composition machinery
        # preserves the invariant in general (already covered by other
        # disease-state tests) -- a mandatory code-reviewer pass
        # correctly found the original version of this test passed
        # identically under a no-op {} modifier, so it wasn't actually
        # discriminating CKD's own fix. The default-value assertions
        # below make it fail if CKD's modifier is ever accidentally
        # reverted to {}.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.CHRONIC_KIDNEY_DISEASE])
        )
        default = PathwayParameters()
        assert result.personalized_params.p3 != pytest.approx(default.p3)
        assert result.personalized_params.ke != pytest.approx(default.ke)
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        assert pathway.validate_y0(tolerance=1e-8) == []

    def test_ckd_reduces_p3_and_ke_end_to_end(self):
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.CHRONIC_KIDNEY_DISEASE])
        )
        default = PathwayParameters()
        assert result.personalized_params.p3 == pytest.approx(default.p3 * 0.86)
        assert result.personalized_params.ke == pytest.approx(default.ke * 0.915)
        assert result.personalized_params.secretion_scale == pytest.approx(default.secretion_scale)

    def test_ckd_composes_with_type_2_diabetes(self):
        # CKD is a common real comorbidity with T2D, deliberately not
        # mutually exclusive (unlike the T1D/T2D pairing) -- confirm the
        # two compose multiplicatively rather than raising or silently
        # dropping one.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[
                DiseaseState.TYPE_2_DIABETES, DiseaseState.CHRONIC_KIDNEY_DISEASE,
            ])
        )
        default = PathwayParameters()
        assert result.personalized_params.p3 == pytest.approx(default.p3 * 0.55 * 0.86)
        assert result.personalized_params.ke == pytest.approx(default.ke * 0.915)
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        assert pathway.validate_y0(tolerance=1e-8) == []

    def test_pcos_result_is_a_genuine_fixed_point(self):
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.PCOS])
        )
        default = PathwayParameters()
        assert result.personalized_params.p3 != pytest.approx(default.p3)
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        assert pathway.validate_y0(tolerance=1e-8) == []

    def test_pcos_reduces_p3_end_to_end(self):
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.PCOS])
        )
        default = PathwayParameters()
        assert result.personalized_params.p3 == pytest.approx(default.p3 * 0.625)
        assert result.personalized_params.ke == pytest.approx(default.ke)
        assert result.personalized_params.secretion_scale == pytest.approx(default.secretion_scale)

    def test_pcos_composes_with_type_2_diabetes(self):
        # PCOS is a common real comorbidity with T2D, deliberately not
        # mutually exclusive (same non-exclusive category as CKD above).
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[
                DiseaseState.TYPE_2_DIABETES, DiseaseState.PCOS,
            ])
        )
        default = PathwayParameters()
        assert result.personalized_params.p3 == pytest.approx(default.p3 * 0.55 * 0.625)
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        assert pathway.validate_y0(tolerance=1e-8) == []

    def test_acromegaly_result_is_a_genuine_fixed_point(self):
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.ACROMEGALY])
        )
        default = PathwayParameters()
        assert result.personalized_params.Gb != pytest.approx(default.Gb)
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        assert pathway.validate_y0(tolerance=1e-8) == []

    def test_acromegaly_raises_gb_end_to_end(self):
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.ACROMEGALY])
        )
        default = PathwayParameters()
        assert result.personalized_params.Gb == pytest.approx(default.Gb * (5.9 / 5.2))
        assert result.personalized_params.p3 == pytest.approx(default.p3)
        assert result.personalized_params.ke == pytest.approx(default.ke)
        assert result.personalized_params.secretion_scale == pytest.approx(default.secretion_scale)

    def test_acromegaly_composes_with_type_2_diabetes(self):
        # Acromegaly-induced diabetes is a real, recognized clinical
        # entity -- deliberately not mutually exclusive (same
        # non-exclusive category as CKD/PCOS above).
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[
                DiseaseState.TYPE_2_DIABETES, DiseaseState.ACROMEGALY,
            ])
        )
        default = PathwayParameters()
        assert result.personalized_params.Gb == pytest.approx(default.Gb * 1.30 * (5.9 / 5.2))
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        assert pathway.validate_y0(tolerance=1e-8) == []

    def test_acromegaly_composed_with_max_metformin_stays_near_but_above_gb_floor(self):
        # G1 "survives composition" check, computed directly rather than
        # assumed: acromegaly's own Gb-raising factor (1.1346) composed
        # with metformin's own strongest chronic Gb-lowering factor
        # (2000mg) lands the RAW (pre-floor-clamp) value at ~4.84 mM --
        # just above the MIN_GB_THRESHOLD_MARGIN_MM floor (G_threshold
        # 4.5 + 0.3 = 4.8 mM), not clamped, but close enough that this is
        # pinned as a regression rather than left to silently drift
        # across a floor if either factor changes later.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.ACROMEGALY]),
            drug_doses=[DrugDose(Drug.METFORMIN, dose_mg=2000.0, is_chronic=True)],
        )
        assert result.personalized_params.Gb > 4.8
        assert result.personalized_params.Gb == pytest.approx(4.839, abs=0.01)
        assert result.clamp_warnings == []

    def test_acromegaly_composed_with_max_tirzepatide_hits_gb_floor_and_discloses_it(self):
        # The mirror case: tirzepatide at its max dose (15mg) ALREADY
        # clamps to the Gb floor with zero acromegaly present (raw
        # ~2.98 mM). Adding acromegaly's own raising factor cannot make
        # this worse -- the composed raw value moves up (to ~3.38 mM)
        # but is still below the floor, so the pre-existing
        # clamp_warnings mechanism (not a new one) correctly fires.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.ACROMEGALY]),
            drug_doses=[DrugDose(Drug.TIRZEPATIDE, dose_mg=15.0)],
        )
        assert result.personalized_params.Gb == pytest.approx(4.8)
        assert any("Gb" in w or "fasting glucose" in w for w in result.clamp_warnings)

    def test_osa_result_is_a_genuine_fixed_point(self):
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.OSA])
        )
        default = PathwayParameters()
        assert result.personalized_params.p3 != pytest.approx(default.p3)
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        assert pathway.validate_y0(tolerance=1e-8) == []

    def test_osa_reduces_p3_end_to_end(self):
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.OSA])
        )
        default = PathwayParameters()
        assert result.personalized_params.p3 == pytest.approx(default.p3 * (5.75 / 6.79))
        assert result.personalized_params.ke == pytest.approx(default.ke)
        assert result.personalized_params.secretion_scale == pytest.approx(default.secretion_scale)

    def test_osa_composes_with_type_2_diabetes(self):
        # OSA is a common real comorbidity with T2D, deliberately not
        # mutually exclusive (same non-exclusive category as CKD/PCOS/
        # acromegaly above) -- only excluded against
        # desaturation_events_per_hour, its own double-counting risk.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[
                DiseaseState.TYPE_2_DIABETES, DiseaseState.OSA,
            ])
        )
        default = PathwayParameters()
        assert result.personalized_params.p3 == pytest.approx(default.p3 * 0.55 * (5.75 / 6.79))
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        assert pathway.validate_y0(tolerance=1e-8) == []

    def test_osa_composed_with_pcos_and_t2d_stays_well_above_the_p3_floor(self):
        # G1 "survives composition" check, computed directly rather than
        # assumed: OSA's own p3 factor (0.8469) composed with two other
        # p3-lowering states (T2D 0.55, PCOS 0.625) -- the most p3-
        # lowering realistic stack this platform can currently produce
        # from static disease modifiers alone -- to confirm the personalizer's
        # own p3 floor (_PARAM_BOUNDS, 1e-8) is nowhere close to binding.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[
                DiseaseState.TYPE_2_DIABETES, DiseaseState.PCOS, DiseaseState.OSA,
            ])
        )
        default = PathwayParameters()
        expected = default.p3 * 0.55 * 0.625 * (5.75 / 6.79)
        assert result.personalized_params.p3 == pytest.approx(expected)
        # code-reviewer finding: this stack lands at ~2.63e-7, ~26x (not
        # "two orders of magnitude") above the 1e-8 floor -- assert the
        # precise margin actually observed, not an approximate one.
        assert result.personalized_params.p3 > 1e-7  # one order above the 1e-8 floor
        assert result.clamp_warnings == []

    def test_nafld_result_is_a_genuine_fixed_point(self):
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.NAFLD])
        )
        default = PathwayParameters()
        assert result.personalized_params.p3 != pytest.approx(default.p3)
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        assert pathway.validate_y0(tolerance=1e-8) == []

    def test_nafld_reduces_only_p3_end_to_end(self):
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.NAFLD])
        )
        default = PathwayParameters()
        p = result.personalized_params
        assert p.p3 == pytest.approx(default.p3 * 0.70)
        assert p.Gb == pytest.approx(default.Gb)
        assert p.ke == pytest.approx(default.ke)
        assert p.secretion_scale == pytest.approx(default.secretion_scale)

    def test_nafld_composes_with_type_2_diabetes(self):
        # A real, common comorbidity -- deliberately not mutually
        # exclusive (same category as CKD/PCOS/OSA).
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[
                DiseaseState.TYPE_2_DIABETES, DiseaseState.NAFLD,
            ])
        )
        default = PathwayParameters()
        assert result.personalized_params.p3 == pytest.approx(default.p3 * 0.55 * 0.70)
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        assert pathway.validate_y0(tolerance=1e-8) == []

    def test_nafld_stacked_with_t2d_pcos_osa_stays_above_the_p3_floor(self):
        # G1 "survives composition": a T2D-branch stack (T2D 0.55 x PCOS
        # 0.625 x OSA 0.847 x NAFLD 0.70), computed directly: ~1.84e-7,
        # ~18x above the 1e-8 p3 floor, no clamp warning. NOT the legal
        # worst case (security review): GDM (p3 0.16) composes with these
        # and, with age/BMI/sleep/drug axes added, the legal minimum
        # reaches ~6e-10 -- covered by the next test.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[
                DiseaseState.TYPE_2_DIABETES, DiseaseState.PCOS,
                DiseaseState.OSA, DiseaseState.NAFLD,
            ])
        )
        default = PathwayParameters()
        expected = default.p3 * 0.55 * 0.625 * (5.75 / 6.79) * 0.70
        assert result.personalized_params.p3 == pytest.approx(expected)
        assert result.personalized_params.p3 > 1e-7
        assert result.clamp_warnings == []

    def test_gdm_branch_worst_case_stack_trips_the_p3_clamp_warning(self):
        # Security review (NAFLD): GDM-branch disease stacks go far lower
        # than the T2D branch. Confirm the personalizer's own generic
        # clamp-and-warn fires as designed rather than silently flooring.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(
                disease_states=[
                    DiseaseState.GESTATIONAL_DIABETES, DiseaseState.INSULIN_RESISTANT,
                    DiseaseState.PCOS, DiseaseState.NAFLD,
                    DiseaseState.CHRONIC_KIDNEY_DISEASE, DiseaseState.HYPOTHYROIDISM,
                ],
                age_years=100.0, sleep_hours=0.0,
                body_weight_kg=120.0, height_cm=170.0,  # BMI ~41
                fitness_level=FitnessLevel.SEDENTARY, time_of_day_hours=22.0,
            ),
        )
        # Raw composition is ~5.9e-9 (below the 1e-8 floor), verified
        # directly; the floor clamps it and the warning fires.
        assert result.personalized_params.p3 == pytest.approx(1e-8)  # clamped at the floor
        assert any("p3" in w for w in result.clamp_warnings)

    def test_cold_acclimation_result_is_a_genuine_fixed_point(self):
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(cold_acclimatized=True)
        )
        default = PathwayParameters()
        assert result.personalized_params.p3 != pytest.approx(default.p3)
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        assert pathway.validate_y0(tolerance=1e-8) == []

    def test_cold_acclimation_increases_p3_end_to_end(self):
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(cold_acclimatized=True)
        )
        default = PathwayParameters()
        assert result.personalized_params.p3 == pytest.approx(default.p3 * 1.43)
        assert result.personalized_params.ke == pytest.approx(default.ke)
        assert result.personalized_params.secretion_scale == pytest.approx(default.secretion_scale)

    def test_cold_acclimation_composes_with_type_2_diabetes(self):
        # Hanssen's own cohort IS type 2 diabetes -- this is the
        # citation's own studied population, not an extrapolation.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[
                DiseaseState.TYPE_2_DIABETES,
            ], cold_acclimatized=True)
        )
        default = PathwayParameters()
        assert result.personalized_params.p3 == pytest.approx(default.p3 * 0.55 * 1.43)
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        assert pathway.validate_y0(tolerance=1e-8) == []

    def test_cold_acclimation_composed_with_largest_p3_raising_factors_stays_well_below_ceiling(self):
        # D1 audit note: this is the platform's first p3-INCREASING
        # environmental axis, still touches the same _PARAM_BOUNDS
        # ceiling (1e-4) every other p3-touching mechanism does. This test
        # went through two rounds of correction, both from the same
        # mistake: claiming to have found "the largest"/"true worst case"
        # p3-raising stack without actually enumerating every p3-raising
        # mechanism in the codebase first. Round 1 (code-reviewer): missed
        # pioglitazone (45mg, 2.7/1.6=1.6875), which is larger than
        # fitness/sex/exercise individually. Round 2 (security-reviewer,
        # independently, same mistake caught a second time): missed
        # circadian time-of-day (up to 1.15x at time_of_day_hours=10.0)
        # and underweight BMI (1.05x). Composed with fitness (1.55), sex
        # (1.10), a vigorous just-completed post-exercise bout (1.62),
        # pioglitazone at 45mg (1.6875), circadian at 10:00 (verified
        # directly, not assumed), and underweight BMI (verified via a
        # 45kg/170cm profile, BMI~15.6) -- lands at ~7.27e-6, ~13.8x below
        # the ceiling. Deliberately NOT claimed as exhaustive this time --
        # the lesson from both rounds is that this codebase has enough
        # independent p3-raising mechanisms that a hand-picked "largest N"
        # list is easy to get wrong; this composition is offered as a
        # real, directly-verified data point with a real margin, not a
        # proof of the mathematical maximum.
        from biosim_physiology import ExerciseBout
        from biosim_physiology.drug_interactions import Drug, DrugDose
        personalizer = PhysiologicalPersonalizer()
        bout = ExerciseBout(
            intensity_percent_vo2max=100.0, duration_minutes=120, hours_before_simulation=0.0,
        )
        dose = DrugDose(Drug.PIOGLITAZONE, dose_mg=45.0, is_chronic=True)
        result = personalizer.personalize(
            PhysiologicalProfile(
                fitness_level=FitnessLevel.ATHLETE,
                biological_sex=BiologicalSex.FEMALE,
                cold_acclimatized=True,
                time_of_day_hours=10.0,
                body_weight_kg=45.0,
                height_cm=170.0,
            ),
            drug_doses=[dose],
            exercise_bout=bout,
        )
        assert result.personalized_params.p3 < 1e-5  # well below the 1e-4 ceiling
        # The pre-existing Gb clamp (unrelated to this test's own p3
        # concern -- exercise_bout also suppresses Gb, which this
        # specific composition happens to push below the safety floor)
        # is expected and disclosed here, not a p3-ceiling finding.
        assert any("Gb" in w or "fasting glucose" in w for w in result.clamp_warnings)

    def test_hydroxychloroquine_composed_with_largest_p3_raising_factors_stays_well_below_ceiling(self):
        # Same "not claimed exhaustive" discipline as the cold-
        # acclimation test above -- HCQ is a NEW p3-RAISING mechanism
        # (1.20x), so it must be checked against the same _PARAM_BOUNDS
        # ceiling (1e-4) directly, not assumed safe by extrapolation.
        # Composed here with the identical worst-case stack the cold-
        # acclimation test uses (fitness/sex/cold-acclimation/circadian/
        # underweight-BMI/pioglitazone), PLUS HCQ on top -- directly
        # computed (not hand-multiplied) at 8.72e-6, still ~11.5x below
        # the ceiling.
        from biosim_physiology import ExerciseBout
        from biosim_physiology.drug_interactions import Drug, DrugDose
        personalizer = PhysiologicalPersonalizer()
        bout = ExerciseBout(
            intensity_percent_vo2max=100.0, duration_minutes=120, hours_before_simulation=0.0,
        )
        pioglitazone = DrugDose(Drug.PIOGLITAZONE, dose_mg=45.0, is_chronic=True)
        hcq = DrugDose(Drug.HYDROXYCHLOROQUINE, dose_mg=400.0)
        result = personalizer.personalize(
            PhysiologicalProfile(
                fitness_level=FitnessLevel.ATHLETE,
                biological_sex=BiologicalSex.FEMALE,
                cold_acclimatized=True,
                time_of_day_hours=10.0,
                body_weight_kg=45.0,
                height_cm=170.0,
            ),
            drug_doses=[pioglitazone, hcq],
            exercise_bout=bout,
        )
        assert result.personalized_params.p3 < 1e-5  # well below the 1e-4 ceiling
        # Pin that HCQ actually participated (code-reviewer finding: the
        # ceiling assertion alone also passes at 7.27e-6 without HCQ).
        assert any(m.source.startswith("drug_hydroxychloroquine") for m in result.applied_modifiers)
        # Same pre-existing, unrelated Gb clamp as the cold-acclimation
        # test above (exercise_bout's own Gb suppression, not an HCQ
        # finding).
        assert any("Gb" in w or "fasting glucose" in w for w in result.clamp_warnings)

    def test_acromegaly_composes_with_growth_hormone_without_double_counting(self):
        # Mandatory code-reviewer finding: profile.py's and
        # drug_interactions.py's own non-overlap reasoning (ACROMEGALY
        # touches only Gb, Drug.GROWTH_HORMONE touches only p3 -- the
        # reason no mutual exclusion was added between the chronic
        # disease state and the acute drug event) was previously
        # unverified by any test. Confirms both move independently to
        # exactly their own isolated values with zero clamp interaction.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.ACROMEGALY]),
            drug_doses=[DrugDose(Drug.GROWTH_HORMONE, dose_mg=1.0)],
        )
        default = PathwayParameters()
        assert result.personalized_params.Gb == pytest.approx(default.Gb * (5.9 / 5.2))
        assert result.personalized_params.p3 == pytest.approx(default.p3 * (1.0 - 0.39))
        assert result.clamp_warnings == []
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        assert pathway.validate_y0(tolerance=1e-8) == []

    def test_sleep_restriction_reduces_p3_end_to_end(self):
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(sleep_hours=5.0)
        )
        default = PathwayParameters()
        assert result.personalized_params.p3 == pytest.approx(default.p3 * 0.89)
        assert result.personalized_params.secretion_scale == pytest.approx(default.secretion_scale)

    def test_full_sleep_produces_no_modifier_end_to_end(self):
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(sleep_hours=10.0)
        )
        default = PathwayParameters()
        assert result.personalized_params.p3 == pytest.approx(default.p3)

    def test_sleep_restriction_result_is_a_genuine_fixed_point(self):
        # RED-proofed per this session's own established convention
        # (see the CKD entry's own code-reviewer finding) -- confirms
        # the fixed point holds for sleep's OWN p3 value specifically.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(sleep_hours=5.0)
        )
        default = PathwayParameters()
        assert result.personalized_params.p3 != pytest.approx(default.p3)
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        assert pathway.validate_y0(tolerance=1e-8) == []

    def test_sleep_restriction_composes_with_type_2_diabetes(self):
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(
                disease_states=[DiseaseState.TYPE_2_DIABETES],
                sleep_hours=5.0,
            )
        )
        default = PathwayParameters()
        assert result.personalized_params.p3 == pytest.approx(default.p3 * 0.55 * 0.89)
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        assert pathway.validate_y0(tolerance=1e-8) == []

    def test_intermittent_hypoxia_reduces_p3_and_p1_end_to_end(self):
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(desaturation_events_per_hour=24.3)
        )
        default = PathwayParameters()
        assert result.personalized_params.p3 == pytest.approx(default.p3 * (3.4 / 4.1))
        assert result.personalized_params.p1 == pytest.approx(default.p1 * (1.3 / 1.9))
        assert result.personalized_params.secretion_scale == pytest.approx(
            default.secretion_scale
        )

    def test_no_hypoxia_exposure_produces_no_modifier_end_to_end(self):
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(desaturation_events_per_hour=0.0)
        )
        default = PathwayParameters()
        assert result.personalized_params.p3 == pytest.approx(default.p3)
        assert result.personalized_params.p1 == pytest.approx(default.p1)

    def test_intermittent_hypoxia_result_is_a_genuine_fixed_point(self):
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(desaturation_events_per_hour=24.3)
        )
        default = PathwayParameters()
        assert result.personalized_params.p3 != pytest.approx(default.p3)
        assert result.personalized_params.p1 != pytest.approx(default.p1)
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        assert pathway.validate_y0(tolerance=1e-8) == []

    def test_intermittent_hypoxia_composes_with_sleep_restriction(self):
        # Genuinely orthogonal axes (the source study induces hypoxia
        # DURING WAKEFULNESS, deliberately isolated from sleep loss) --
        # confirms both p3 factors multiply together rather than one
        # silently overwriting the other.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(sleep_hours=5.0, desaturation_events_per_hour=24.3)
        )
        default = PathwayParameters()
        assert result.personalized_params.p3 == pytest.approx(
            default.p3 * 0.89 * (3.4 / 4.1)
        )
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        assert pathway.validate_y0(tolerance=1e-8) == []

    def test_luteal_phase_reduces_p3_end_to_end(self):
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(menstrual_cycle_phase=MenstrualCyclePhase.LUTEAL)
        )
        default = PathwayParameters()
        assert result.personalized_params.p3 == pytest.approx(default.p3 * (3.20 / 6.20))
        assert result.personalized_params.p1 == pytest.approx(default.p1)

    def test_follicular_and_midcycle_produce_no_modifier_end_to_end(self):
        personalizer = PhysiologicalPersonalizer()
        default = PathwayParameters()
        for phase in (MenstrualCyclePhase.FOLLICULAR, MenstrualCyclePhase.MIDCYCLE):
            result = personalizer.personalize(
                PhysiologicalProfile(menstrual_cycle_phase=phase)
            )
            assert result.personalized_params.p3 == pytest.approx(default.p3)

    def test_luteal_phase_result_is_a_genuine_fixed_point(self):
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(menstrual_cycle_phase=MenstrualCyclePhase.LUTEAL)
        )
        default = PathwayParameters()
        assert result.personalized_params.p3 != pytest.approx(default.p3)
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        assert pathway.validate_y0(tolerance=1e-8) == []

    def test_luteal_phase_composes_with_sleep_restriction(self):
        # Genuinely orthogonal axes (distinct physiological mechanisms,
        # each with its own independent citation) -- confirms both p3
        # factors multiply together rather than one overwriting the
        # other.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(
                sleep_hours=5.0, menstrual_cycle_phase=MenstrualCyclePhase.LUTEAL,
            )
        )
        default = PathwayParameters()
        assert result.personalized_params.p3 == pytest.approx(
            default.p3 * 0.89 * (3.20 / 6.20)
        )
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        assert pathway.validate_y0(tolerance=1e-8) == []

    def test_athlete_increases_p3(self):
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(fitness_level=FitnessLevel.ATHLETE)
        )
        base_p3 = PathwayParameters().p3
        assert result.personalized_params.p3 > base_p3

    def test_morning_vs_evening_p3(self):
        personalizer = PhysiologicalPersonalizer()
        r_morning = personalizer.personalize(PhysiologicalProfile(time_of_day_hours=10.0))
        r_evening = personalizer.personalize(PhysiologicalProfile(time_of_day_hours=22.0))
        assert r_morning.personalized_params.p3 > r_evening.personalized_params.p3

    def test_no_exercise_bout_leaves_params_unchanged(self):
        # Backward compatibility: omitting exercise_bout (the default)
        # must behave exactly as before this parameter existed.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(PhysiologicalProfile())
        assert result.personalized_params.p3 == pytest.approx(PathwayParameters().p3)

    def test_exercise_bout_increases_p3(self):
        from biosim_physiology import ExerciseBout, ExerciseIntensity
        personalizer = PhysiologicalPersonalizer()
        bout = ExerciseBout.from_intensity(ExerciseIntensity.MODERATE, 45.0, 1.0)
        result = personalizer.personalize(PhysiologicalProfile(), exercise_bout=bout)
        assert result.personalized_params.p3 > PathwayParameters().p3

    def test_gb_threshold_margin_never_collapses_below_safety_bound(self):
        # Regression: PathwayParameters.F_secretion_basal is inversely
        # proportional to (Gb - G_threshold), by design, to hit an exact
        # analytic fixed point at basal glucose. Gb and G_threshold were
        # bounded independently with no joint constraint, so an ordinary,
        # realistic Gb-reducing modifier (metformin, at its own max
        # chronic dose) alone shrinks the margin to a few hundredths of a
        # mM, blowing F_secretion_basal (and the resulting insulin
        # response to a meal) up by an order of magnitude or more -- a
        # numerical artifact of the fixed-point formula, not a modeled
        # mechanism.
        #
        # Correction (fresh Track 1 scoping fan-out, 2026-08-18): this
        # test originally stacked metformin WITH an alcohol dose, on the
        # (then-correct) premise that alcohol also reduced Gb. Alcohol's
        # mechanism was found to be actively wrong (contradicted by both
        # available primary sources) and fixed to a transient, ODE-
        # coupled hepatic_egp_delta contribution instead -- it no longer
        # touches Gb or this static personalization pipeline AT ALL
        # (confirmed directly: adding an alcohol DrugDose here produces
        # a byte-identical Gb to metformin alone). The alcohol dose is
        # removed rather than left in as now-inert dead weight; metformin
        # alone already exercises this clamp fully.
        personalizer = PhysiologicalPersonalizer()
        doses = [
            DrugDose(Drug.METFORMIN, dose_mg=2000.0, is_chronic=True, hours_since_ingestion=24 * 7),
        ]
        result = personalizer.personalize(PhysiologicalProfile(), drug_doses=doses)
        p = result.personalized_params
        assert p.Gb - p.G_threshold >= 0.3 - 1e-9

    def test_gb_threshold_margin_clamp_bounds_the_insulin_spike(self):
        # See test_gb_threshold_margin_never_collapses_below_safety_bound's
        # own comment above -- the alcohol dose was removed (2026-08-18)
        # since alcohol no longer touches Gb at all; metformin alone
        # already engages this clamp.
        personalizer = PhysiologicalPersonalizer()
        doses = [
            DrugDose(Drug.METFORMIN, dose_mg=2000.0, is_chronic=True, hours_since_ingestion=24 * 7),
        ]
        result = personalizer.personalize(PhysiologicalProfile(), drug_doses=doses)
        meals = [high_sugar_meal()]
        sim = ODESystem(pathway=GlucoseInsulinPathway(params=result.personalized_params, meals=meals)) \
            .simulate(duration_s=14400, meals=meals)
        base_sim = ODESystem(pathway=GlucoseInsulinPathway(params=PathwayParameters(), meals=meals)) \
            .simulate(duration_s=14400, meals=meals)
        # Before the original fix this was ~4.7x; a real insulin
        # sensitiser (metformin) alone should not multiply peak insulin
        # by anywhere near that -- bound it generously at 2x.
        assert sim.species["I"].max() < base_sim.species["I"].max() * 2.0

    def test_clamp_warning_fires_when_the_safety_clamp_engages(self):
        # MASTER_PLAN.md 1.23 (2026-08-17 self-evaluation): the clamp
        # above silently overrode a modifier's own composed Gb value
        # with no signal anywhere -- a /verify pass found all three
        # tirzepatide doses (5/10/15mg) produce byte-identical output
        # for this platform's own default T2D preset because the clamp
        # floors every one of them at the same value. This test
        # confirms the clamp now records a warning naming both the
        # requested and clamped value when it actually engages.
        personalizer = PhysiologicalPersonalizer()
        doses = [DrugDose(Drug.TIRZEPATIDE, dose_mg=15.0)]
        result = personalizer.personalize(PhysiologicalProfile(), drug_doses=doses)
        assert len(result.clamp_warnings) == 1
        warning = result.clamp_warnings[0]
        assert "MIN_GB_THRESHOLD_MARGIN_MM" in warning
        assert "2.9" in warning  # requested Gb (5.0 * 0.596 = 2.98)
        assert "4.8" in warning  # clamped Gb (G_threshold 4.5 + 0.3 margin)
        # Section 8 housekeeping (part c): Gb also has its own point-of-use
        # clamp inside dydt() -- the joint clamp's own warning should say so.
        assert "point-of-use clamp inside the ODE solver" in warning

    def test_clamp_warning_absent_when_the_clamp_never_engages(self):
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(PhysiologicalProfile())
        assert result.clamp_warnings == []

    def test_clamp_warning_absent_for_a_dose_that_stays_clear_of_the_margin(self):
        # A mild, single Gb-lowering modifier that doesn't approach the
        # margin at all should not produce a warning -- confirms the
        # check is genuinely conditional, not firing unconditionally
        # whenever any Gb-lowering drug is present.
        #
        # Correction (mandatory security-reviewer finding on the alcohol
        # mechanism fix, 2026-08-18): this test originally used a mild
        # Drug.ALCOHOL dose. Alcohol's mechanism was found to be wrong
        # (see test_gb_threshold_margin_never_collapses_below_safety_
        # bound's own comment) and fixed to a transient, ODE-coupled
        # hepatic_egp_delta contribution that does NOT touch Gb at all --
        # a "mild" alcohol dose is now a ZERO-effect dose, not a small
        # one, so this test no longer exercised "a modifier that stays
        # clear of the margin" as its own comment claims (a Gb-lowering
        # modifier that happens to lower Gb by exactly 0.0 is a
        # degenerate case, not a genuine near-miss). Replaced with a low
        # metformin dose, which genuinely still reduces Gb (a real,
        # non-zero effect) but stays comfortably clear of the margin.
        personalizer = PhysiologicalPersonalizer()
        doses = [DrugDose(Drug.METFORMIN, dose_mg=100.0, is_chronic=True)]
        result = personalizer.personalize(PhysiologicalProfile(), drug_doses=doses)
        assert result.personalized_params.Gb < PathwayParameters().Gb   # a real, non-zero effect
        assert result.personalized_params.Gb - result.personalized_params.G_threshold > 0.3
        assert result.clamp_warnings == []

    def test_generic_bound_clamp_warning_fires_for_t1d_ib_and_secretion_scale(self):
        # D1 Composition-Boundary Behavioral Audit (2026-08-18): every
        # _PARAM_BOUNDS-clamped field except Gb/G_threshold silently
        # clamped with no warning at all. T1D's own Ib_factor=0.05/
        # secretion_scale_factor=0.02 were ALREADY, since row 1.10,
        # deliberately chosen to clamp to _PARAM_BOUNDS' own floors (see
        # test_t1d_ib_and_secretion_scale_clamp_to_engineering_floor
        # above) -- meaning this exact, already-shipped scenario has been
        # silently clamping with zero visibility since T1D was added.
        # This is the first regression test to actually check for a
        # warning here, not just the clamped numeric result.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.TYPE_1_DIABETES])
        )
        assert len(result.clamp_warnings) == 2
        joined = " ".join(result.clamp_warnings)
        assert "secretion_scale" in joined and "0.02" in joined and "0.05" in joined
        assert "Ib" in joined and "3.3" in joined and "10" in joined

    def test_clamp_warning_caveat_present_only_for_ode_coupled_params(self):
        # Section 8 housekeeping (part c, found during the
        # CORTISOL_P3_SENSITIVITY security review): a warning saying a
        # value "was clamped to X" is potentially misleading for the four
        # parameters (Gb/p1/p3/secretion_scale) that ALSO have their own
        # point-of-use clamp inside GlucoseInsulinPathway.dydt() -- an
        # ODE-coupled pathway can move the value further during
        # simulation, past this static clamp, with no other disclosure.
        # secretion_scale's warning (from the SAME T1D scenario as the
        # test above) should now name this explicitly; Ib's should NOT,
        # since Ib has no corresponding rule-11 signal/point-of-use clamp
        # at all -- the static clamp is fully accurate for it.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.TYPE_1_DIABETES])
        )
        secretion_scale_warning = next(
            w for w in result.clamp_warnings if "secretion_scale" in w
        )
        ib_warning = next(w for w in result.clamp_warnings if "Ib" in w and "secretion_scale" not in w)
        assert "point-of-use clamp inside the ODE solver" in secretion_scale_warning
        assert "point-of-use clamp inside the ODE solver" not in ib_warning

    def test_generic_bound_clamp_warning_absent_when_no_generic_bound_engages(self):
        # The Gb/G_threshold joint clamp (row 1.23) and the generic
        # per-parameter bound clamp (this fix) are deliberately
        # independent -- a tirzepatide dose that engages the FORMER
        # must not also spuriously produce a generic-bound warning for
        # any OTHER parameter it doesn't touch.
        personalizer = PhysiologicalPersonalizer()
        doses = [DrugDose(Drug.TIRZEPATIDE, dose_mg=15.0)]
        result = personalizer.personalize(PhysiologicalProfile(), drug_doses=doses)
        assert len(result.clamp_warnings) == 1  # only the Gb/G_threshold one
        assert "MIN_GB_THRESHOLD_MARGIN_MM" in result.clamp_warnings[0]

    def test_generic_bound_clamp_warning_excludes_gb_and_g_threshold_lower_bound_only(self):
        # Gb/G_threshold's LOWER bound is deliberately excluded from the
        # generic warning even if _PARAM_BOUNDS itself would have clamped
        # them -- the joint clamp below is more specific for that
        # direction and already produces its own warning; a second,
        # independent warning for the same field would be confusing, not
        # more informative. (Renamed from "..._excludes_gb_and_g_threshold"
        # -- a security-reviewer finding on row 1.37/ACROMEGALY found the
        # UPPER bound direction was ALSO silently excluded, with no other
        # disclosure mechanism at all; fixed below, see the two new tests
        # directly following this one.)
        personalizer = PhysiologicalPersonalizer()
        doses = [DrugDose(Drug.TIRZEPATIDE, dose_mg=15.0)]
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.TYPE_1_DIABETES]),
            drug_doses=doses,
        )
        for warning in result.clamp_warnings:
            if "Gb" in warning or "G_threshold" in warning:
                assert "MIN_GB_THRESHOLD_MARGIN_MM" in warning

    def test_generic_bound_clamp_warning_now_covers_gb_upper_bound(self):
        # Security-reviewer finding on row 1.37 (ACROMEGALY): the
        # previous blanket exclusion of Gb/G_threshold from the generic
        # warning silently covered BOTH directions, but the joint clamp
        # below only ever discloses the LOWER direction -- the UPPER
        # bound (_PARAM_BOUNDS['Gb'][1] == 15.0) had no disclosure
        # mechanism at all. No currently shippable, mutual-exclusion-
        # respecting combination reaches 15.0 today (confirmed: the
        # widest realistic stack found by the security review tops out
        # around 14.4 mM) -- this test uses a synthetic elevated
        # base_params purely to exercise the code path at the ceiling,
        # not as a claim that this base value is realistic.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.ACROMEGALY]),
            base_params=PathwayParameters(Gb=14.0),
        )
        assert result.personalized_params.Gb == pytest.approx(15.0)
        assert any(
            "Gb" in w and "15" in w for w in result.clamp_warnings
        )

    def test_generic_bound_clamp_warning_still_absent_for_gb_lower_bound(self):
        # The fix above must not resurrect a duplicate warning for the
        # LOWER direction -- that stays the joint clamp's own job.
        personalizer = PhysiologicalPersonalizer()
        doses = [DrugDose(Drug.TIRZEPATIDE, dose_mg=15.0)]
        result = personalizer.personalize(PhysiologicalProfile(), drug_doses=doses)
        assert len(result.clamp_warnings) == 1
        assert "MIN_GB_THRESHOLD_MARGIN_MM" in result.clamp_warnings[0]

    def test_joint_clamp_warning_reports_true_raw_value_not_intermediate_clamped_value(self):
        # Security-reviewer MEDIUM finding on row 1.37: the joint clamp's
        # warning message used to read the ALREADY-generic-bound-clamped
        # Gb value (e.g. floored at 2.0), not the true raw composed
        # value, understating how far off the calibrated result actually
        # was in an extreme-enough scenario. Synthetic low base_params
        # used purely to exercise the code path (no shipped mechanism
        # alone produces a starting Gb this low).
        personalizer = PhysiologicalPersonalizer()
        doses = [DrugDose(Drug.METFORMIN, dose_mg=2000.0, is_chronic=True)]
        result = personalizer.personalize(
            PhysiologicalProfile(),
            base_params=PathwayParameters(Gb=1.0),
            drug_doses=doses,
        )
        joined = " ".join(result.clamp_warnings)
        # True raw value: 1.0 * (5.7/6.4) = 0.853 mM -- NOT the
        # generic-bound-floored intermediate (2.0 mM).
        assert "0.853" in joined
        assert "2.000" not in joined

    def test_generic_bound_clamp_warning_reaches_simulation_via_indinavir_stacking(self):
        # Track 1 row 1.33 (indinavir) stacked with the same extreme
        # disease/physiological combination the D1 audit itself swept
        # (GDM+insulin_resistant+hypothyroidism+CKD+sleep_hours=5+
        # desaturation_events_per_hour=24.3) does NOT yet breach the p3
        # floor (confirmed directly: ratio ~0.041, still ~3.7x above
        # it) -- this pins that finding as a regression so a future
        # change to any of these factors is caught if it silently closes
        # (or further widens) the margin, without asserting a clamp that
        # doesn't currently occur.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(
                disease_states=[
                    DiseaseState.GESTATIONAL_DIABETES, DiseaseState.INSULIN_RESISTANT,
                    DiseaseState.HYPOTHYROIDISM, DiseaseState.CHRONIC_KIDNEY_DISEASE,
                ],
                sleep_hours=5.0, desaturation_events_per_hour=24.3,
            ),
            drug_doses=[DrugDose(Drug.INDINAVIR, dose_mg=1200.0)],
        )
        default_p3 = PathwayParameters().p3
        ratio = result.personalized_params.p3 / default_p3
        assert 0.01 < ratio < 0.10  # not yet clamped, but within an order of magnitude of it
        assert not any("p3" in w for w in result.clamp_warnings)

    def test_p3_floor_is_now_breached_by_a_realistic_extreme_stack_and_correctly_disclosed(self):
        # Mandatory security-reviewer finding on row 1.29 (PCOS): the
        # margin the test above pins (indinavir alone, ~3.7x above the
        # floor) has since collapsed further -- stacking PCOS on top of
        # the same extreme disease/physiological combination PLUS the
        # other new drugs added this session (hydrochlorothiazide,
        # growth hormone, atorvastatin) and menstrual_cycle_phase=LUTEAL
        # genuinely DOES breach the p3 floor now. This is NOT a silently
        # discovered bug -- the generic clamp_warnings mechanism (row
        # 1.33) was already built specifically for this scenario and is
        # confirmed here to fire correctly, not silently clamp.
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(
            PhysiologicalProfile(
                disease_states=[
                    DiseaseState.GESTATIONAL_DIABETES, DiseaseState.INSULIN_RESISTANT,
                    DiseaseState.HYPOTHYROIDISM, DiseaseState.CHRONIC_KIDNEY_DISEASE,
                    DiseaseState.PCOS,
                ],
                sleep_hours=5.0, desaturation_events_per_hour=24.3,
                menstrual_cycle_phase=MenstrualCyclePhase.LUTEAL,
            ),
            drug_doses=[
                DrugDose(Drug.INDINAVIR, dose_mg=1200.0),
                DrugDose(Drug.HYDROCHLOROTHIAZIDE, dose_mg=40.0, is_chronic=True),
                DrugDose(Drug.GROWTH_HORMONE, dose_mg=1.0),
                DrugDose(Drug.ATORVASTATIN, dose_mg=80.0, is_chronic=True),
            ],
        )
        assert result.personalized_params.p3 == pytest.approx(1e-8)  # clamped to the floor
        assert any("p3" in w and "1e-08" in w for w in result.clamp_warnings)

    def test_clamp_warning_reaches_simulation_output_through_the_real_runner(self):
        # End-to-end: confirms the warning survives the glucocorticoid
        # block's dataclasses.replace() and actually reaches
        # SimulationOutput.warnings, not just the personalizer's own
        # return value.
        from biosim_api.runner import SimulationRunner
        from biosim_api.scenario import DrugEvent, ScenarioDefinition

        scenario = ScenarioDefinition(
            name="tirzepatide_clamp_warning",
            duration_s=14400,
            meal_events=[],
            drug_events=[DrugEvent(time_s=0, drug="tirzepatide", dose_mg=15.0)],
        )
        output = SimulationRunner().run(scenario)
        assert any("MIN_GB_THRESHOLD_MARGIN_MM" in w for w in output.warnings)

    def test_clamp_warning_is_genuinely_mechanism_agnostic_not_tirzepatide_specific(self):
        # MASTER_PLAN.md Section 11 (D2 coverage matrix): computing the
        # real crossover thresholds found that metformin's OWN
        # calibrated effect (Hansen et al. 2024, -14.7% Gb reduction at
        # max 2000mg dose) also engages this same clamp for 3 of 5
        # shipped disease presets (insulin-resistant, gestational
        # diabetes, and the healthy default), not just tirzepatide --
        # and correctly stays silent for T2D, metformin's own
        # calibration population, which keeps it clear of the margin.
        # Locks this in so a future "simplify to tirzepatide only"
        # edit doesn't silently narrow the fix's real scope.
        personalizer = PhysiologicalPersonalizer()
        doses = [DrugDose(Drug.METFORMIN, dose_mg=2000.0, is_chronic=True, hours_since_ingestion=24 * 7)]

        for profile in (
            PhysiologicalProfile(),
            PhysiologicalProfile(disease_states=[DiseaseState.INSULIN_RESISTANT]),
            PhysiologicalProfile(disease_states=[DiseaseState.GESTATIONAL_DIABETES]),
        ):
            result = personalizer.personalize(profile, drug_doses=doses)
            assert len(result.clamp_warnings) == 1, profile.disease_states

        t2d_result = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES]), drug_doses=doses
        )
        assert t2d_result.clamp_warnings == []

    def test_first_phase_gain_abolished_above_brunzell_threshold(self):
        # Regression: first_phase_gain (the model's only representation
        # of acute/first-phase insulin release, Mari et al. 2002,
        # PMID:12388151) was a bare constant with no consumer that could
        # ever change it -- T2D's disease modifier scaled it identically
        # to second-phase secretion via the shared secretion_scale
        # factor, contradicting Brunzell et al. 1976 (JCEM, PMID:1262429,
        # n=66, IVGTT, verified live): "The acute insulin response was
        # present in subjects with fasting glucose levels below 115 mg/dl
        # but was absent above this level." 115 mg/dL = 115/18.0 ~= 6.39
        # mM. The T2D preset's own Gb (5.0 * 1.30 = 6.5 mM) sits above
        # that threshold, yet retained 55-60% of first-phase gain before
        # this fix. Gated on the FINAL effective Gb (post-composition,
        # post-safety-clamp), not on "is T2D present" -- this directly
        # implements Brunzell's own glucose-level-based finding rather
        # than a diagnosis-based one, so it generalises correctly to any
        # combination of modifiers that happens to cross the threshold.
        personalizer = PhysiologicalPersonalizer()
        profile = PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES])
        result = personalizer.personalize(profile)
        assert result.personalized_params.Gb >= 115.0 / 18.0
        assert result.personalized_params.first_phase_gain == pytest.approx(0.0)

    def test_first_phase_gain_preserved_below_brunzell_threshold(self):
        # INSULIN_RESISTANT's own Gb factor (1.03) keeps fasting glucose
        # well below the Brunzell threshold (5.0 * 1.03 = 5.15 mM) --
        # Brunzell's own finding is that the acute response is PRESENT
        # below 115 mg/dL, so it must not be touched here. This is the
        # gate correctly discriminating pre-diabetic from diabetic
        # fasting glucose, not a hardcoded per-disease-state special case.
        personalizer = PhysiologicalPersonalizer()
        profile = PhysiologicalProfile(disease_states=[DiseaseState.INSULIN_RESISTANT])
        result = personalizer.personalize(profile)
        assert result.personalized_params.Gb < 115.0 / 18.0
        assert result.personalized_params.first_phase_gain == pytest.approx(
            PathwayParameters().first_phase_gain
        )

    def test_first_phase_gain_preserved_for_healthy_profile(self):
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(PhysiologicalProfile())
        assert result.personalized_params.first_phase_gain == pytest.approx(
            PathwayParameters().first_phase_gain
        )

    def test_first_phase_gain_abolition_does_not_perturb_basal_steady_state(self):
        # first_phase_secretion = first_phase_gain * max(dG, 0) *
        # secretion_scale -- exactly zero at any basal steady state
        # (dG=0 there) for ANY gain value, so abolishing it cannot
        # perturb the analytic fixed point. Verify this directly rather
        # than trust the algebra alone, matching this project's own
        # established discipline for every fixed-point-adjacent change.
        personalizer = PhysiologicalPersonalizer()
        profile = PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES])
        result = personalizer.personalize(profile)
        assert result.personalized_params.first_phase_gain == pytest.approx(0.0)
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        violations = pathway.validate_y0(tolerance=1e-8)
        assert violations == [], violations

    def test_personalized_params_runs_ode(self):
        personalizer = PhysiologicalPersonalizer()
        profile = PhysiologicalProfile(
            age_years=50,
            biological_sex=BiologicalSex.FEMALE,
            body_weight_kg=80,
            fitness_level=FitnessLevel.LIGHTLY_ACTIVE,
            disease_states=[DiseaseState.INSULIN_RESISTANT],
            time_of_day_hours=20.0,
        )
        result = personalizer.personalize(profile)
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        system = ODESystem(pathway=pathway)
        sim = system.simulate_meal(high_sugar_meal(), duration_s=7200)
        assert sim.success

    def test_personalized_validate_y0_passes(self):
        personalizer = PhysiologicalPersonalizer()
        profile = PhysiologicalProfile(age_years=55, disease_states=[DiseaseState.TYPE_2_DIABETES])
        result = personalizer.personalize(profile)
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        violations = pathway.validate_y0(tolerance=1e-3)
        assert violations == [], f"y0 not at steady state: {violations}"

    def test_with_drugs(self):
        personalizer = PhysiologicalPersonalizer()
        profile = PhysiologicalProfile()
        drugs = [
            DrugDose(Drug.CAFFEINE, dose_mg=200.0, hours_since_ingestion=1.0),
        ]
        result = personalizer.personalize(profile, drug_doses=drugs)
        base_p3 = PathwayParameters().p3
        # Caffeine should reduce p3
        assert result.personalized_params.p3 < base_p3

    def test_summary_structure(self):
        personalizer = PhysiologicalPersonalizer()
        result = personalizer.personalize(PhysiologicalProfile(age_years=45))
        s = result.summary()
        assert "age" in s
        assert "n_modified" in s
        assert "net_factors" in s

    def test_t2d_higher_glucose_peak(self):
        personalizer = PhysiologicalPersonalizer()
        healthy = personalizer.personalize(PhysiologicalProfile())
        diabetic = personalizer.personalize(
            PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES])
        )
        sys_h = ODESystem(GlucoseInsulinPathway(params=healthy.personalized_params))
        sys_d = ODESystem(GlucoseInsulinPathway(params=diabetic.personalized_params))
        peak_h = sys_h.simulate_meal(high_sugar_meal(), duration_s=7200).glucose_peak()[0]
        peak_d = sys_d.simulate_meal(high_sugar_meal(), duration_s=7200).glucose_peak()[0]
        assert peak_d > peak_h


class TestPostExerciseP3CoefficientCalibration:
    """
    2026-08-12 deep-search plan, item B1: post_exercise_modifiers()'s
    max_p3_boost coefficient recalibrated from an uncited 0.50 to 0.31,
    anchored to Mikines et al. 1988 (PMID:3126668) and Wang et al. 2013
    (PMID:23493357), both live-verified this session. See
    post_exercise_modifiers()'s own docstring for the full derivation
    and disclosed intensity-mapping assumptions.
    """

    def test_max_p3_boost_matches_the_citation_derived_coefficient(self):
        from biosim_physiology import ExerciseBout

        # Immediately post-exercise (hours_before_simulation=0), decay=1.0,
        # so p3_multiplier == 1.0 + 0.31*volume exactly.
        bout = ExerciseBout(intensity_percent_vo2max=100.0, duration_minutes=60.0)
        volume = bout.relative_intensity * bout.duration_hours   # 1.0 * 1.0 = 1.0
        expected = 1.0 + 0.31 * volume
        assert bout.post_exercise_modifiers()["p3"] == pytest.approx(expected)

    def test_recalibration_reduced_the_boost_from_the_old_uncited_value(self):
        # Regression guard against silently reverting to the old,
        # uncited 0.50 coefficient.
        from biosim_physiology import ExerciseBout

        bout = ExerciseBout(intensity_percent_vo2max=100.0, duration_minutes=60.0)
        volume = bout.relative_intensity * bout.duration_hours
        old_boost = 1.0 + 0.50 * volume
        new_boost = bout.post_exercise_modifiers()["p3"]
        assert new_boost < old_boost


class TestPostExerciseP1AndGbCitationDisclosure:
    """
    Fresh Track 1 scoping (2026-08-18): p1_boost's 0.10 and
    gb_suppression's 0.03 were previously uncited. A literature search
    (live-verified) found genuinely DIVERGENT evidence for p1 (real
    citations exist, but two studies using this platform's own SG
    convention are null in exactly the matched populations) and no
    usable quantitative anchor at all for Gb. Deliberately NOT upgraded
    to derived point estimates -- see post_exercise_modifiers()'s own
    docstring for the full reasoning. These tests lock in that the
    disclosure is present and that the underlying VALUES are unchanged
    (a documentation-only pass, not a recalibration).
    """

    def test_p1_boost_value_is_unchanged_by_the_disclosure(self):
        from biosim_physiology import ExerciseBout

        bout = ExerciseBout(intensity_percent_vo2max=100.0, duration_minutes=60.0)
        volume = bout.relative_intensity * bout.duration_hours
        expected = 1.0 + 0.10 * volume
        assert bout.post_exercise_modifiers()["p1"] == pytest.approx(expected)

    def test_gb_suppression_value_is_unchanged_by_the_disclosure(self):
        from biosim_physiology import ExerciseBout

        bout = ExerciseBout(intensity_percent_vo2max=100.0, duration_minutes=60.0)
        volume = bout.relative_intensity * bout.duration_hours
        expected = 1.0 - 0.03 * volume
        assert bout.post_exercise_modifiers()["Gb"] == pytest.approx(expected)

    def test_docstring_discloses_the_p1_literature_search_and_its_divergence(self):
        from biosim_physiology import ExerciseBout

        doc = ExerciseBout.post_exercise_modifiers.__doc__
        assert "15840752" in doc  # Hayashi 2005
        assert "33905292" in doc  # Farkas 2022 (null in able-bodied controls)
        assert "18448376" in doc  # Bordenave 2008 (null in T2D)
        assert "not derivable" in doc
        assert "convergent evidence" in doc

    def test_docstring_discloses_gb_remains_an_open_gap(self):
        from biosim_physiology import ExerciseBout

        doc = ExerciseBout.post_exercise_modifiers.__doc__
        assert "29679309" in doc  # Nakayama 2019, blocked full text
        assert "open literature gap" in doc
