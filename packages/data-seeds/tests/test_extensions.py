# SAVE TO: packages\data-seeds\tests\test_extensions.py
"""Tests for exercise, cortisol pathway, and sensitivity-curation bridge."""
import math
import pytest
from biosim_physiology.exercise import (
    ExerciseBout, ExerciseIntensity, light_walk, moderate_run, hiit_session,
)
from biosim_ode.cortisol_pathway import CortisolPathway, CortisolParameters
from biosim_ode import ODESystem, GlucoseInsulinPathway, high_sugar_meal
from biosim_uncertainty.curation_bridge import (
    prioritise_from_sensitivity, PARAM_BIOLOGY_MAP, BridgeResult,
)


# ---------------------------------------------------------------------------
# Exercise
# ---------------------------------------------------------------------------

class TestExerciseBout:
    def test_construction(self):
        e = ExerciseBout(intensity_percent_vo2max=60.0, duration_minutes=45.0)
        assert e.intensity_percent_vo2max == 60.0

    def test_invalid_intensity_raises(self):
        with pytest.raises(ValueError):
            ExerciseBout(intensity_percent_vo2max=110.0, duration_minutes=30.0)

    def test_invalid_duration_raises(self):
        with pytest.raises(ValueError):
            ExerciseBout(intensity_percent_vo2max=60.0, duration_minutes=0.0)

    def test_from_intensity(self):
        e = ExerciseBout.from_intensity(ExerciseIntensity.MODERATE, 45.0)
        assert e.intensity_percent_vo2max == 50.0

    def test_light_walk(self):
        e = light_walk()
        assert e.intensity_percent_vo2max == ExerciseIntensity.LIGHT.percent_vo2max

    def test_post_exercise_increases_p3(self):
        e = ExerciseBout(70.0, 60.0, hours_before_simulation=0.0)
        assert e.post_exercise_modifiers()["p3"] > 1.0

    def test_effect_decays_with_time(self):
        e0  = ExerciseBout(70.0, 60.0, hours_before_simulation=0.0)
        e24 = ExerciseBout(70.0, 60.0, hours_before_simulation=24.0)
        assert e0.post_exercise_modifiers()["p3"] > e24.post_exercise_modifiers()["p3"]

    def test_higher_intensity_larger_effect(self):
        low  = ExerciseBout(30.0, 60.0)
        high = ExerciseBout(80.0, 60.0)
        assert high.post_exercise_modifiers()["p3"] > low.post_exercise_modifiers()["p3"]

    def test_during_exercise_increases_clearance(self):
        assert ExerciseBout(70.0, 45.0).during_exercise_modifiers()["glucose_uptake_delta"] > 0.0

    def test_during_exercise_suppresses_secretion(self):
        assert ExerciseBout(70.0, 45.0).during_exercise_modifiers()["secretion_scale"] < 1.0

    def test_mets(self):
        assert ExerciseBout(60.0, 30.0).mets == pytest.approx(6.0)

    def test_post_exercise_meal_lower_peak(self):
        import dataclasses
        from biosim_ode import PathwayParameters
        base_params = PathwayParameters()
        e = moderate_run(45.0, hours_before=1.0)
        mods = e.post_exercise_modifiers()
        ex_params = dataclasses.replace(base_params, **{
            k: getattr(base_params, k) * v
            for k, v in mods.items() if hasattr(base_params, k)
        })
        peak_base, _ = ODESystem(GlucoseInsulinPathway(params=base_params)).simulate_meal(high_sugar_meal(), 7200).glucose_peak()
        peak_ex,   _ = ODESystem(GlucoseInsulinPathway(params=ex_params)).simulate_meal(high_sugar_meal(), 7200).glucose_peak()
        assert peak_ex < peak_base


# ---------------------------------------------------------------------------
# Cortisol pathway
# ---------------------------------------------------------------------------

class TestCortisolPathway:
    def test_instantiation(self):
        assert CortisolPathway().name == "cortisol_v1"

    def test_species_names(self):
        p = CortisolPathway()
        assert "Cortisol" in p.species_names
        assert "GR_occ" in p.species_names

    def test_y0_length(self):
        assert len(CortisolPathway().y0) == 3

    def test_y0_positive(self):
        assert all(v >= 0 for v in CortisolPathway().y0)

    def test_dydt_returns_3_values(self):
        p = CortisolPathway()
        assert len(p.dydt(0.0, p.y0)) == 3

    def test_output_signals_keys(self):
        p = CortisolPathway()
        signals = p.output_signals(0.0, p.y0)
        assert "cortisol_nM" in signals
        assert "gr_occupancy" in signals
        assert "p3_modifier" in signals

    def test_basal_gb_modifier_near_one(self):
        signals = CortisolPathway().output_signals(0.0, CortisolPathway().y0)
        assert 0.95 <= signals["gb_modifier"] <= 1.10

    def test_stress_raises_cortisol(self):
        from scipy.integrate import solve_ivp
        def run_cortisol(pathway):
            sol = solve_ivp(pathway.dydt, [0, 600], pathway.y0, method="Radau")
            return sol.y[0, -1]
        assert run_cortisol(CortisolPathway().with_stress(0.8)) > run_cortisol(CortisolPathway())

    def test_morning_higher_than_midnight(self):
        p_m = CortisolPathway().with_time_of_day(7.0)
        p_n = CortisolPathway().with_time_of_day(2.0)
        assert p_m.y0[0] > p_n.y0[0]

    def test_high_cortisol_reduces_p3_signal(self):
        # Use GR occupancy well above the analytic reference point
        # (GR_REFERENCE = hill_activation(CORTISOL_BASAL_NM, EC50_GR_NM,
        # HILL_GR) ~= 0.718 -- see the KNOWN BUG note on GR_OCC_BASAL).
        from biosim_ode.cortisol_pathway import GR_REFERENCE
        p = CortisolPathway()
        high_y = [0.80, 5.0, 0.95]   # cortisol, CRH, GR_occ=0.95 >> GR_REFERENCE
        assert 0.95 > GR_REFERENCE   # sanity: the test value must exceed the reference
        signals = p.output_signals(0.0, high_y)
        assert signals["p3_modifier"] < 0.95

    def test_y0_is_genuine_fixed_point(self):
        # Regression for the GR_OCC_BASAL calibration bug: y0 must be an
        # actual fixed point of dydt (all derivatives ~0), not just a
        # plausible-looking initial guess. This was previously false --
        # GR_occ silently drifted from 0.15 toward ~0.72 over about an
        # hour even at zero stress/steroid.
        p = CortisolPathway()
        assert p.validate_y0(tolerance=1e-6) == []

    def test_p3_modifier_neutral_at_true_basal_over_time(self):
        # The practical consequence of the bug above: composing default
        # (zero stress/steroid) cortisol must not silently change p3 as
        # the simulation progresses. Integrate standalone for an hour and
        # confirm p3_modifier stays at 1.0 throughout, not just at t=0.
        from scipy.integrate import solve_ivp
        p = CortisolPathway()
        sol = solve_ivp(p.dydt, [0, 60], p.y0, method="Radau", rtol=1e-8, atol=1e-10)
        signals_end = p.output_signals(60.0, sol.y[:, -1])
        assert signals_end["p3_modifier"] == pytest.approx(1.0, abs=1e-4)

    def test_gb_and_p3_modifier_neutral_over_a_real_600_minute_fasted_run(self):
        # Item 6 of the glucocorticoid PK/PD milestone (see CLAUDE.md):
        # this project's own established technique for this exact bug
        # class (KNOWN BUGS FIXED #10-#14) is a full-trajectory
        # integration, not just a t=0 check -- a static y0 check alone
        # cannot detect drift that only manifests after the solver has
        # actually moved forward. Re-run at 600+ minutes (not just the 60
        # minutes above) specifically because this session changed BOTH
        # gb_modifier's formula (soft-capped via CORTISOL_GB_SATURATION_NM)
        # and exogenous_gc_nM's evaluation (now time-gated via
        # exogenous_gc_nM_at(t) instead of a bare attribute read) -- either
        # change could in principle have reintroduced basal drift that a
        # short 60-minute window wouldn't reveal.
        from scipy.integrate import solve_ivp
        p = CortisolPathway()
        sol = solve_ivp(p.dydt, [0, 600], p.y0, method="Radau", rtol=1e-8, atol=1e-10)
        signals_end = p.output_signals(600.0, sol.y[:, -1])
        assert signals_end["p3_modifier"] == pytest.approx(1.0, abs=1e-4)
        assert signals_end["gb_modifier"] == pytest.approx(1.0, abs=1e-4)

    def test_gb_modifier_time_gated_prednisone_neutral_before_dose(self):
        # Companion check for the SAME invariant with a real prednisone
        # dose scheduled to occur LATER in the simulated window (dose_time_
        # min > 0, e.g. an in-progress scenario where the dose hasn't been
        # taken yet at t=0) -- gb_modifier/p3_modifier must stay neutral
        # for the entire pre-dose interval, not just at t=0, since
        # exogenous_gc_nM_at() gates on (t - dose_time_min), which is
        # negative (not yet dosed) throughout this whole window.
        from scipy.integrate import solve_ivp
        DOSE_AT_MINUTE = 700.0   # dose scheduled well after this 600-min window
        p = CortisolPathway(CortisolParameters.for_prednisone(20.0, dose_time_min=DOSE_AT_MINUTE))
        sol = solve_ivp(p.dydt, [0, 600], p.y0, method="Radau", rtol=1e-8, atol=1e-10)
        signals_end = p.output_signals(600.0, sol.y[:, -1])
        assert signals_end["p3_modifier"] == pytest.approx(1.0, abs=1e-4)
        assert signals_end["gb_modifier"] == pytest.approx(1.0, abs=1e-4)

    def test_hpa_feedback_zero_at_basal_state(self):
        # excess_gr_feedback must be exactly 0 at y0 (gr_occ == GR_REFERENCE
        # there), so the feedback term cannot perturb the basal invariant
        # already covered by test_y0_is_genuine_fixed_point.
        p = CortisolPathway()
        dCRH_default = p.dydt(0.0, p.y0)[1]
        assert dCRH_default == pytest.approx(0.0, abs=1e-9)

    def test_chronic_steroid_suppresses_endogenous_cortisol(self):
        # The clinically well-known mechanism this feedback exists to
        # capture: sustained high-dose exogenous steroid suppresses the
        # HPA axis, driving *endogenous* cortisol production below its
        # own basal level (adrenal suppression). Integrate long enough
        # for the fast GR/CRH dynamics to settle.
        from scipy.integrate import solve_ivp
        p = CortisolPathway(CortisolParameters.for_prednisone(20.0))
        sol = solve_ivp(p.dydt, [0, 1000], p.y0, method="Radau", rtol=1e-8, atol=1e-10)
        endogenous_cortisol_final = sol.y[0, -1]
        assert endogenous_cortisol_final < p.params.cortisol_basal

    def test_stress_still_raises_cortisol_with_feedback_present(self):
        # Regression: adding HPA feedback must not eliminate the stress
        # response -- stress-driven GR occupancy elevation is far smaller
        # than steroid-saturated occupancy, so suppression stays mild.
        from scipy.integrate import solve_ivp
        default_final = solve_ivp(
            CortisolPathway().dydt, [0, 600], CortisolPathway().y0, method="Radau"
        ).y[0, -1]
        stressed = CortisolPathway(CortisolParameters(stress_level=0.8))
        stressed_final = solve_ivp(stressed.dydt, [0, 600], stressed.y0, method="Radau").y[0, -1]
        assert stressed_final > default_final

    def test_implements_base_pathway(self):
        from biosim_ode.base_pathway import BasePathway
        assert isinstance(CortisolPathway(), BasePathway)

    def test_composite_with_glucose_insulin(self):
        from biosim_ode.coupling import CompositePathway
        composite = CompositePathway([GlucoseInsulinPathway(), CortisolPathway()])
        assert len(composite.y0) == 9

    def test_composite_dydt_runs(self):
        from biosim_ode.coupling import CompositePathway
        composite = CompositePathway([GlucoseInsulinPathway(), CortisolPathway()])
        assert len(composite.dydt(0.0, composite.y0)) == 9

    def test_couple_with_gi_preserves_basal_steady_state(self):
        from biosim_ode.cortisol_pathway import couple_with_gi
        composite = couple_with_gi(GlucoseInsulinPathway(), CortisolPathway())
        derivs = composite.dydt(0.0, composite.y0)
        assert derivs[0] == pytest.approx(0.0, abs=1e-6)   # dG
        assert derivs[1] == pytest.approx(0.0, abs=1e-6)   # dI

    def test_couple_with_gi_prednisone_raises_glucose_peak(self):
        from biosim_ode.cortisol_pathway import couple_with_gi
        from biosim_ode import ODESystem, high_sugar_meal

        baseline = ODESystem(pathway=GlucoseInsulinPathway()) \
            .simulate_meal(high_sugar_meal(), duration_s=14400)

        steroid_composite = couple_with_gi(
            GlucoseInsulinPathway(), CortisolPathway(CortisolParameters.for_prednisone(20.0))
        )
        treated = ODESystem(pathway=steroid_composite).simulate_meal(high_sugar_meal(), duration_s=14400)

        assert baseline.success and treated.success
        peak_base, _ = baseline.glucose_peak()
        peak_treat, _ = treated.glucose_peak()
        assert peak_treat > peak_base

    def test_couple_with_gi_stress_raises_glucose_peak(self):
        from biosim_ode.cortisol_pathway import couple_with_gi
        from biosim_ode import ODESystem, high_sugar_meal

        baseline = ODESystem(pathway=GlucoseInsulinPathway()) \
            .simulate_meal(high_sugar_meal(), duration_s=14400)

        stressed_composite = couple_with_gi(
            GlucoseInsulinPathway(), CortisolPathway(CortisolParameters(stress_level=0.9))
        )
        stressed = ODESystem(pathway=stressed_composite).simulate_meal(high_sugar_meal(), duration_s=14400)

        assert baseline.success and stressed.success
        peak_base, _ = baseline.glucose_peak()
        peak_stress, _ = stressed.glucose_peak()
        assert peak_stress > peak_base

    def test_for_prednisone_default_matches_no_steroid(self):
        default = CortisolParameters()
        steroid = CortisolParameters.for_prednisone(0.0)
        assert steroid.exogenous_gc_nM == pytest.approx(default.exogenous_gc_nM)

    def test_for_prednisone_no_immediate_effect_at_t0(self):
        # Regression (real PK/PD timing, see CLAUDE.md): the glycaemic
        # effect genuinely has NOT developed yet at the instant of dosing
        # -- glucocorticoids act genomically (GR binding -> transcription
        # -> protein synthesis -> effect), a mechanism that lags plasma
        # concentration by hours (Burt 2011 PMID:21411550, Yuen 2012
        # PMID:21973241). for_prednisone()'s default dose_time_min=0.0
        # ("just taken", matching DrugDose's own convention) must therefore
        # show NO drive difference from an undosed patient at t=0 --
        # unlike the old always-on-immediately behaviour this replaces.
        y0 = CortisolPathway().y0
        dGR_default = CortisolPathway().dydt(0.0, y0)[2]
        dGR_steroid = CortisolPathway(CortisolParameters.for_prednisone(20.0)).dydt(0.0, y0)[2]
        assert dGR_steroid == pytest.approx(dGR_default, abs=1e-9)

    def test_for_prednisone_raises_gr_target_within_onset_window(self):
        # At the PD curve's ~5.7h peak (see glucocorticoid_pd_effect_
        # fraction()), the drive IS stronger with prednisone present --
        # dGR is a relaxation toward a target, not an instant jump, so
        # verify at a time when the effect has actually developed rather
        # than at t=0 (see test_for_prednisone_no_immediate_effect_at_t0
        # for why t=0 must NOT show this).
        PEAK_MIN = 342.0   # ~5.7h, this curve's peak (minutes)
        y0 = CortisolPathway().y0
        dGR_default = CortisolPathway().dydt(PEAK_MIN, y0)[2]
        dGR_steroid = CortisolPathway(CortisolParameters.for_prednisone(20.0)).dydt(PEAK_MIN, y0)[2]
        assert dGR_steroid > dGR_default

    def test_for_prednisone_gb_modifier_raised_within_onset_window(self):
        # Regression: gb_modifier is now ALSO time-gated (exogenous_gc_nM_
        # at()), not "responds immediately, unlike GR occupancy" as this
        # test used to assert -- see CLAUDE.md's glucocorticoid PK/PD
        # milestone entry. Evaluate at the curve's peak, matching the GR
        # occupancy test above.
        PEAK_MIN = 342.0
        p = CortisolPathway(CortisolParameters.for_prednisone(20.0))
        default = CortisolPathway()
        signals = p.output_signals(PEAK_MIN, p.y0)
        default_signals = default.output_signals(PEAK_MIN, default.y0)
        assert signals["gb_modifier"] > default_signals["gb_modifier"]

    def test_for_prednisone_gb_modifier_not_yet_elevated_at_t0(self):
        # Direct companion to test_for_prednisone_no_immediate_effect_at_t0
        # for the Gb signal specifically: fasting glucose must not jump the
        # instant the dose is taken -- Yuen et al. 2012 (PMID:21973241)
        # measured fasting glucose UNCHANGED in all groups despite 3 days
        # of morning prednisone dosing.
        p = CortisolPathway(CortisolParameters.for_prednisone(20.0))
        default = CortisolPathway()
        signals = p.output_signals(0.0, p.y0)
        default_signals = default.output_signals(0.0, default.y0)
        assert signals["gb_modifier"] == pytest.approx(default_signals["gb_modifier"], abs=1e-9)

    def test_for_prednisone_higher_dose_larger_effect(self):
        PEAK_MIN = 342.0
        y0 = CortisolPathway().y0
        dGR_low = CortisolPathway(CortisolParameters.for_prednisone(5.0)).dydt(PEAK_MIN, y0)[2]
        dGR_high = CortisolPathway(CortisolParameters.for_prednisone(40.0)).dydt(PEAK_MIN, y0)[2]
        assert dGR_high > dGR_low

    def test_gb_modifier_bounded_at_extreme_dose(self):
        # Regression (KNOWN BUG, see CLAUDE.md): gb_modifier used to be
        # plain-linear and UNBOUNDED in dose -- a 5000mg "dose" would
        # produce an absurd, unbounded fasting-glucose multiplier. Now
        # soft-capped via CORTISOL_GB_SATURATION_NM; confirm a huge dose
        # doesn't blow past a sane ceiling, while a moderate dose still
        # shows real, distinguishable elevation (not clamped to the same
        # value).
        PEAK_MIN = 342.0
        y0 = CortisolPathway().y0
        moderate = CortisolPathway(CortisolParameters.for_prednisone(20.0))
        extreme = CortisolPathway(CortisolParameters.for_prednisone(5000.0))
        gb_moderate = moderate.output_signals(PEAK_MIN, y0)["gb_modifier"]
        gb_extreme = extreme.output_signals(PEAK_MIN, y0)["gb_modifier"]
        assert gb_extreme > gb_moderate   # still monotonic
        assert gb_extreme < 2.0           # but genuinely bounded, not blown up

    def test_pd_effect_fraction_shape_matches_calibration_windows(self):
        # Direct check of the calibration windows glucocorticoid_pd_effect_
        # fraction() is documented to match: Burt et al. 2011
        # (PMID:21411550) found glucose unchanged 0000-1200h (i.e. flat
        # through ~4h post an 08:00 dose) and Yuen et al. 2012
        # (PMID:21973241) found the elevation concentrated midday-to-
        # midnight, "dissipating overnight".
        from biosim_ode.cortisol_pathway import glucocorticoid_pd_effect_fraction as f
        assert f(0.0) == 0.0
        assert f(-1.0) == 0.0
        assert f(2.0) < 0.2          # negligible early (Burt's flat window)
        assert f(6.0) > 0.9          # near-peak within the 6-10h window
        assert f(24.0) < 0.15        # dissipated by the next morning

    def test_pd_effect_fraction_repeats_daily_for_chronic_dosing(self):
        # Regression (found by code review, glucocorticoid PK/PD milestone
        # -- see CLAUDE.md): a genuinely chronically-dosed patient
        # (DrugEvent.hours_before_simulation spanning multiple days, this
        # platform's existing convention for ongoing therapy) must show
        # the SAME within-day onset/peak/resolution pattern recurring each
        # day, not a cumulative decay to permanent zero -- Yuen et al.
        # 2012 (PMID:21973241) dosed prednisone 20mg EACH morning for 3
        # days and observed the identical daily pattern recur, not a
        # fading effect. Any exact multiple of 24h must equal f(0)=0 (a
        # fresh cycle start); day 2/3/8 at the same within-day offset must
        # equal day 1 at that offset -- confirming the curve is genuinely
        # periodic, not a single decaying pulse evaluated at a huge
        # cumulative time (which is what it silently did before this fix).
        from biosim_ode.cortisol_pathway import glucocorticoid_pd_effect_fraction as f
        assert f(24.0) == pytest.approx(f(0.0), abs=1e-9)   # both 0.0
        assert f(6.0) == pytest.approx(f(24.0 + 6.0), abs=1e-9)
        assert f(6.0) == pytest.approx(f(24.0 * 7 + 6.0), abs=1e-9)
        assert f(6.0) == pytest.approx(f(24.0 * 2 + 6.0), abs=1e-9)
        # Acute (< 24h) dosing is never wrapped (23.9 % 24 == 23.9) --
        # only genuinely chronic (>= 24h) dosing engages the daily cycle.
        assert f(23.9) > 0.0

    def test_for_prednisone_composes_with_stress(self):
        # Verify **kwargs forwarding: a stressed patient also on steroids
        p = CortisolParameters.for_prednisone(20.0, stress_level=0.5)
        assert p.exogenous_gc_nM > 0
        assert p.stress_level == pytest.approx(0.5)

    def test_exogenous_gc_preserved_at_zero_matches_baseline(self):
        # Regression: additive exogenous_gc_nM must be exactly neutral at
        # its 0.0 default (proves the extension can't affect any existing,
        # non-steroid cortisol scenario).
        from biosim_ode.coupling import CompositePathway
        composite_default = CompositePathway([GlucoseInsulinPathway(), CortisolPathway()])
        composite_explicit_zero = CompositePathway(
            [GlucoseInsulinPathway(), CortisolPathway(CortisolParameters(exogenous_gc_nM=0.0))]
        )
        assert composite_default.dydt(0.0, composite_default.y0) == pytest.approx(
            composite_explicit_zero.dydt(0.0, composite_explicit_zero.y0)
        )

    def test_couple_with_gi_no_stress_matches_baseline_closely(self):
        # Default (unstressed) cortisol should barely perturb the meal response
        # -- it must not silently duplicate the existing circadian modifier's
        # effect just by being composed.
        from biosim_ode.cortisol_pathway import couple_with_gi
        from biosim_ode import ODESystem, high_sugar_meal

        baseline = ODESystem(pathway=GlucoseInsulinPathway()) \
            .simulate_meal(high_sugar_meal(), duration_s=14400)
        composite = couple_with_gi(GlucoseInsulinPathway(), CortisolPathway())
        composed = ODESystem(pathway=composite).simulate_meal(high_sugar_meal(), duration_s=14400)

        peak_base, _ = baseline.glucose_peak()
        peak_composed, _ = composed.glucose_peak()
        assert peak_composed == pytest.approx(peak_base, rel=0.01)


# ---------------------------------------------------------------------------
# Curation bridge
# ---------------------------------------------------------------------------

class TestCurationBridge:
    @pytest.fixture(scope="class")
    @classmethod
    def sensitivity_result(cls):
        from biosim_uncertainty import OATSensitivityAnalyzer, STANDARD_METRICS
        return OATSensitivityAnalyzer(0.10).analyze(
            parameter_names=["p1", "ke", "p3"],
            metrics={"glucose_peak": STANDARD_METRICS["glucose_peak"]},
            meals=[high_sugar_meal()], duration_s=7200,
        )

    @pytest.fixture(scope="class")
    @classmethod
    def seeded_deps(cls):
        from biosim_seeds import create_seeded_registry, create_seeded_store
        r = create_seeded_registry()
        return r, create_seeded_store(r)

    @pytest.fixture
    def queue(self):
        from biosim_curation import CurationQueue
        return CurationQueue()

    def test_bridge_returns_result(self, sensitivity_result, queue, seeded_deps):
        r, s = seeded_deps
        result = prioritise_from_sensitivity(sensitivity_result, queue, s, r)
        assert isinstance(result, BridgeResult)

    def test_bridge_result_summary(self, sensitivity_result, queue, seeded_deps):
        r, s = seeded_deps
        result = prioritise_from_sensitivity(sensitivity_result, queue, s, r)
        summary = result.summary()
        assert "n_actions" in summary
        assert "n_items_enqueued" in summary

    def test_param_biology_map_not_empty(self):
        assert len(PARAM_BIOLOGY_MAP) > 0
        assert "p1" in PARAM_BIOLOGY_MAP
        assert "ke" in PARAM_BIOLOGY_MAP

    def test_very_sensitive_param_gets_enqueued(self, seeded_deps):
        from biosim_uncertainty import OATSensitivityAnalyzer, STANDARD_METRICS
        from biosim_curation import CurationQueue, CurationPriority
        r, s = seeded_deps
        q = CurationQueue()
        sr = OATSensitivityAnalyzer(0.15).analyze(
            parameter_names=["p1"],
            metrics={"glucose_peak": STANDARD_METRICS["glucose_peak"]},
            meals=[high_sugar_meal()], duration_s=7200,
        )
        result = prioritise_from_sensitivity(sr, q, s, r,
            high_sensitivity_threshold=0.01,
            critical_sensitivity_threshold=0.5,
        )
        pending = q.pending()
        if pending:
            assert any(item.priority in (CurationPriority.CRITICAL, CurationPriority.HIGH)
                       for item in pending)
