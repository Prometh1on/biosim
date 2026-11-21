# SAVE TO: packages\ode-engine\tests\test_glp1.py
"""
Test suite for the GLP-1 incretin pathway extension and its coupling to
GlucoseInsulinPathway.

Until this file, GLP1Pathway had zero test coverage anywhere in the repo
and was never imported outside its own module -- semaglutide/liraglutide/
DPP-4 inhibitor simulation was effectively unreachable and unverified.
"""
import pytest
import numpy as np
from scipy.integrate import solve_ivp

from biosim_ode.glp1_pathway import GLP1Pathway, GLP1Parameters, couple_with_gi
from biosim_ode.base_pathway import BasePathway, PathwayMaturity
from biosim_ode import GlucoseInsulinPathway, ODESystem, high_sugar_meal, PathwayParameters
from biosim_ode.meal import MealInput


class TestGLP1Parameters:
    def test_default_construction(self):
        p = GLP1Parameters()
        assert p.agonist_scale == pytest.approx(1.0)

    def test_semaglutide_raises_agonist_scale(self):
        default = GLP1Parameters()
        sema = GLP1Parameters.for_semaglutide(1.0)
        assert sema.agonist_scale > default.agonist_scale

    def test_semaglutide_higher_dose_more_occupancy(self):
        low  = GLP1Parameters.for_semaglutide(0.5)
        high = GLP1Parameters.for_semaglutide(2.0)
        assert high.agonist_scale > low.agonist_scale
        assert high.max_incretin > low.max_incretin

    def test_semaglutide_extends_half_life(self):
        p = GLP1Parameters.for_semaglutide(1.0)
        assert p.half_life_min == pytest.approx(10080.0)   # 1 week

    def test_liraglutide_extends_half_life(self):
        p = GLP1Parameters.for_liraglutide(1.2)
        assert p.half_life_min == pytest.approx(780.0)   # 13 hours

    def test_liraglutide_smaller_effect_than_semaglutide(self):
        lira = GLP1Parameters.for_liraglutide(1.2)
        sema = GLP1Parameters.for_semaglutide(1.0)
        assert lira.agonist_scale < sema.agonist_scale

    def test_dpp4_inhibited_doubles_half_life(self):
        default = GLP1Parameters()
        dpp4 = GLP1Parameters.dpp4_inhibited()
        assert dpp4.half_life_min == pytest.approx(default.half_life_min * 2.0)

    def test_dpp4_inhibited_raises_simulated_glp1_peak(self):
        # Regression: dpp4_inhibited() used to also set a `glp1_peak`
        # field, asserted on directly here without ever running a
        # simulation -- but dydt()/output_signals() never read that
        # field at all, so it was a pure no-op (the exact "discarded
        # computation, only asserted on the bare params object" blind
        # spot this session's CLAUDE.md documents catching repeatedly
        # elsewhere: gi_meals kabs, exercise secretion_scale). The field
        # has been removed; this now checks the one mechanism that
        # actually is wired in -- half-life extension -- against the
        # real simulated GLP1 trajectory through the coupled pathway.
        composite = couple_with_gi(GlucoseInsulinPathway(), GLP1Pathway())
        baseline = ODESystem(pathway=composite).simulate_meal(high_sugar_meal(), duration_s=14400)

        dpp4_composite = couple_with_gi(GlucoseInsulinPathway(), GLP1Pathway.dpp4_inhibited())
        treated = ODESystem(pathway=dpp4_composite).simulate_meal(high_sugar_meal(), duration_s=14400)

        assert max(treated.species["GLP1"]) > max(baseline.species["GLP1"])


class TestGLP1Pathway:
    def test_instantiation(self):
        assert GLP1Pathway().name == "glp1_v1"

    def test_implements_base_pathway(self):
        assert isinstance(GLP1Pathway(), BasePathway)

    def test_maturity_beta(self):
        assert GLP1Pathway().maturity == PathwayMaturity.BETA

    def test_species_names(self):
        assert GLP1Pathway().species_names == ["GLP1"]

    def test_y0_at_basal(self):
        # 7.0 pM, Orskov et al. 1994 (PMID:8138058, live-verified) --
        # "mean concentrations of amidated GLP-I ... in fasting plasma
        # were 7 +/- 1 ... pM (n=6)".
        assert GLP1Pathway().y0[0] == pytest.approx(7.0)

    @pytest.mark.parametrize("params", [
        GLP1Parameters(),
        GLP1Parameters.for_semaglutide(1.0),
        GLP1Parameters.for_semaglutide(2.0),
        GLP1Parameters.for_liraglutide(1.2),
        GLP1Parameters.dpp4_inhibited(),
    ])
    def test_y0_is_genuine_fixed_point(self, params):
        # Same bug class as KNOWN BUGS FIXED #10/#11/#12/#14: a preset's
        # y0 must be the actual analytic fixed point of its own ODE, not
        # a separately-hardcoded value that happens to look plausible.
        # dydt() scales the ENTIRE secretion term (basal replacement +
        # meal response) by agonist_scale, so the true fixed point at
        # basal glucose (glucose_stim=0) is glp1_basal * agonist_scale,
        # not glp1_basal alone.
        p = GLP1Pathway(params)
        violations = p.validate_y0(tolerance=1e-8)
        assert violations == [], violations

    @pytest.mark.parametrize("params,label", [
        (GLP1Parameters.for_semaglutide(1.0), "semaglutide_1mg"),
        (GLP1Parameters.for_liraglutide(1.2), "liraglutide_1.2mg"),
    ])
    def test_glp1_stable_over_long_trajectory_at_basal_glucose(self, params, label):
        # Before the fix: semaglutide's y0 (30.03 pM) was far below its
        # own ODE's analytic fixed point (162.3 pM, since agonist_scale
        # scales the basal-replacement term too) -- a real, reachable
        # drift within a standard 5h meal simulation window (+9% at 300
        # min for semaglutide, +48% for liraglutide, confirmed by direct
        # standalone integration). Mirrors
        # test_egp_stable_over_long_trajectory_at_basal_insulin's
        # convention (KNOWN BUGS FIXED #11).
        p = GLP1Pathway(params)
        y = p.y0[0]
        dt = 0.1
        for _ in range(6000):   # 600 minutes
            y += p.dydt(0.0, [y], signals={"glucose_mM": 5.0})[0] * dt
        assert y == pytest.approx(p.y0[0], rel=1e-3), label

    def test_dydt_returns_one_value(self):
        p = GLP1Pathway()
        assert len(p.dydt(0.0, p.y0)) == 1

    def test_glucose_stimulates_secretion(self):
        p = GLP1Pathway()
        d_low  = p.dydt(0.0, p.y0, signals={"glucose_mM": 5.0})[0]
        d_high = p.dydt(0.0, p.y0, signals={"glucose_mM": 15.0})[0]
        assert d_high > d_low

    def test_output_signals_keys(self):
        p = GLP1Pathway()
        signals = p.output_signals(0.0, p.y0)
        assert "glp1_pM" in signals
        assert "incretin_effect" in signals
        assert "gastric_emptying_modifier" in signals

    def test_incretin_effect_rises_with_glp1(self):
        p = GLP1Pathway()
        low  = p.output_signals(0.0, [8.0])["incretin_effect"]
        high = p.output_signals(0.0, [60.0])["incretin_effect"]
        assert high > low

    def test_gastric_emptying_modifier_falls_with_glp1(self):
        p = GLP1Pathway()
        low_glp1_modifier  = p.output_signals(0.0, [8.0])["gastric_emptying_modifier"]
        high_glp1_modifier = p.output_signals(0.0, [60.0])["gastric_emptying_modifier"]
        assert high_glp1_modifier < low_glp1_modifier
        assert high_glp1_modifier >= 0.30   # floored per implementation

    def test_semaglutide_classmethod(self):
        p = GLP1Pathway.semaglutide(1.0)
        assert p.params.agonist_scale > 1.0

    def test_liraglutide_classmethod(self):
        p = GLP1Pathway.liraglutide(1.2)
        assert p.params.agonist_scale > 1.0

    def test_dpp4_inhibited_classmethod(self):
        default = GLP1Pathway()
        dpp4 = GLP1Pathway.dpp4_inhibited()
        assert dpp4.params.half_life_min > default.params.half_life_min

    def test_standalone_integration_stable(self):
        p = GLP1Pathway()
        sol = solve_ivp(
            lambda t, y: p.dydt(t, y, signals={"glucose_mM": 5.0}),
            [0, 60], p.y0, method="Radau",
        )
        assert sol.success
        assert sol.y[0, -1] > 0


class TestGLP1CoupledWithGI:
    def test_composite_shape(self):
        composite = couple_with_gi(GlucoseInsulinPathway(), GLP1Pathway())
        assert len(composite.y0) == 7   # 6 GI species + 1 GLP1 species

    def test_composite_dydt_runs(self):
        composite = couple_with_gi(GlucoseInsulinPathway(), GLP1Pathway())
        derivs = composite.dydt(0.0, composite.y0)
        assert len(derivs) == 7


class TestBasalGlucoseAwareFixedPoint:
    """
    Mandatory security-reviewer finding on the Orskov 1994 recalibration
    (CRITICAL): GLP1Pathway.y0 previously always assumed glucose_stim=0
    (healthy basal glucose <= GLP1_EC50_SECRETION), so composing it with
    a disease-elevated-Gb GlucoseInsulinPathway (T2D ~6.5mM, GCK-MODY
    ~7.17mM, both already shipped elsewhere in this codebase) was NOT
    actually a fixed point -- reproduced directly before this fix:
    T2D (Gb=6.5) + semaglutide 1mg drifted 142.0 -> 3352.0 pM over a
    standard 5h fasted window (the much smaller pre-recalibration
    GLP1_MEAL_RESPONSE=0.035 made this a minor, ~17% 5h drift; the
    Orskov-anchored ~114x larger meal_response amplified it into a
    severe one). Fixed by adding basal_glucose_mM to GLP1Pathway.__init__
    (and threading it through the semaglutide()/liraglutide()/
    dpp4_inhibited() classmethods and SimulationRunner._detect_glp1()),
    so y0 is the genuine analytic fixed point for whatever basal glucose
    the pathway is actually composed with, not always 5.0mM.
    """

    def test_default_basal_glucose_mM_preserves_healthy_behavior(self):
        # basal_glucose_mM defaults to 5.0 (this pathway's own historical
        # default, and dydt()'s own internal fallback) -- every existing
        # healthy-Gb caller is unaffected.
        assert GLP1Pathway().y0[0] == pytest.approx(7.0)
        assert GLP1Pathway.semaglutide(1.0).y0[0] == pytest.approx(
            GLP1Pathway(GLP1Parameters.for_semaglutide(1.0), basal_glucose_mM=5.0).y0[0]
        )

    @pytest.mark.parametrize("basal_glucose_mM,label", [
        (6.5, "t2d_like"),
        (7.17, "gck_mody_like"),
    ])
    def test_disease_elevated_basal_glucose_is_a_genuine_fixed_point(
        self, basal_glucose_mM, label,
    ):
        gi = GlucoseInsulinPathway(params=PathwayParameters(Gb=basal_glucose_mM))
        for glp1_params in [
            GLP1Parameters.for_semaglutide(1.0),
            GLP1Parameters.for_liraglutide(1.2),
            GLP1Parameters.dpp4_inhibited(),
        ]:
            glp1 = GLP1Pathway(glp1_params, basal_glucose_mM=basal_glucose_mM)
            composite = couple_with_gi(gi, glp1)
            derivs = composite.dydt(0.0, composite.y0)
            assert abs(derivs[-1]) < 1e-6, (
                f"{label}: GLP1 derivative at y0 is {derivs[-1]:.6g}, "
                f"not a genuine fixed point"
            )

    def test_disease_elevated_basal_glucose_no_longer_drifts_over_a_long_trajectory(self):
        # RED-proofs the exact reproduction from the security review:
        # T2D (Gb=6.5) + semaglutide 1mg, fasted, over a standard 5h
        # window -- before this fix, drifted 142.0 -> 3352.0 pM.
        gi = GlucoseInsulinPathway(params=PathwayParameters(Gb=6.5))
        glp1 = GLP1Pathway(GLP1Parameters.for_semaglutide(1.0), basal_glucose_mM=6.5)
        composite = couple_with_gi(gi, glp1)
        system = ODESystem(composite)
        result = system.simulate(duration_s=300 * 60, meals=[])
        glp1_pM = result.species["GLP1"] * 1e6
        assert glp1_pM.max() == pytest.approx(glp1_pM.min(), rel=1e-6), (
            f"GLP1 drifted from {glp1_pM[0]:.1f} to {glp1_pM[-1]:.1f} pM "
            f"over a standard 5h fasted window -- basal steady state is "
            f"not actually stable"
        )

    def test_healthy_basal_glucose_still_unaffected_by_the_new_parameter(self):
        # Confirms this fix's default doesn't perturb the already-
        # registered clinical benchmark's own healthy-Gb comparator.
        gi = GlucoseInsulinPathway()   # default healthy Gb=5.0
        glp1 = GLP1Pathway(GLP1Parameters.for_semaglutide(1.0))   # default basal_glucose_mM=5.0
        composite = couple_with_gi(gi, glp1)
        assert composite.validate_y0(tolerance=1e-8) == []


class TestOrskov1994MealResponseCalibration:
    """
    2026-08-12 deep-search plan, item A1: GLP1_MEAL_RESPONSE (and
    GLP1_BASAL_PM) are now anchored to Orskov et al. 1994 (PMID:8138058)
    -- "fasting plasma [GLP-1] were 7 +/- 1 pM ... rose ... amounting to
    41 +/- 5 pM 90 min after the meal ingestion" -- instead of the prior
    uncited GLP1_MEAL_RESPONSE=0.035, which this file's own history
    already flagged as producing a peak only ~1.05x basal (far below
    any real reported range). Peak MAGNITUDE is matched through the
    real composed ODE using this platform's own standard 75g reference
    meal; peak TIMING is deliberately not matched (the paper's own
    breakfast-meal glucose kinetics aren't reported, so there's nothing
    to jointly fit against, unlike acarbose's dual peak-height/peak-time
    calibration elsewhere in this codebase).
    """

    def test_basal_matches_orskov_1994(self):
        assert GLP1Pathway().y0[0] == pytest.approx(7.0)

    @staticmethod
    def _peak_glp1_pM(meal) -> float:
        gi = GlucoseInsulinPathway(meals=[meal])
        glp1 = GLP1Pathway(GLP1Parameters())
        composite = couple_with_gi(gi, glp1)
        result = ODESystem(composite).simulate(duration_s=300 * 60, meals=[meal])
        return float((result.species["GLP1"] * 1e6).max())

    def test_calibration_setup_now_overshoots_orskov_1994_peak(self):
        # GLP1_MEAL_RESPONSE was calibrated to Orskov's 41 +/- 5 pM with the
        # single-exponential 75 g meal AND the double-counted resting glucose
        # clearance. Since B7 (2026-11-16) removed that, the same setup gives
        # ~50 pM. KNOWN GAP, recorded, NOT re-tuned. Pinned to catch drift.
        meal = MealInput(carbs_g=75.0)
        assert self._peak_glp1_pM(meal) == pytest.approx(50.4, abs=1.0)

    def test_default_two_compartment_meal_vs_orskov_peak(self):
        # With the default two-compartment meal the peak was ~30 pM (B3,
        # 2026-11-15); since B7 (2026-11-16) it is ~39 pM, inside Orskov's
        # reported 41 +/- 5 pM. A side effect of B7, not a re-calibration.
        # Pinned so any drift is looked at.
        assert self._peak_glp1_pM(high_sugar_meal()) == pytest.approx(38.9, abs=1.0)

    def test_composed_meal_response_basal_matches_orskov_1994(self):
        gi = GlucoseInsulinPathway(meals=[high_sugar_meal()])
        glp1 = GLP1Pathway(GLP1Parameters())
        composite = couple_with_gi(gi, glp1)
        system = ODESystem(composite)
        result = system.simulate(duration_s=300 * 60, meals=[high_sugar_meal()])
        glp1_pM = result.species["GLP1"] * 1e6
        assert glp1_pM[0] == pytest.approx(7.0)

    def test_old_uncited_meal_response_would_have_produced_a_far_smaller_peak(self):
        # Regression guard against silently reverting to the old,
        # uncited GLP1_MEAL_RESPONSE=0.035 -- confirms this test suite
        # actually discriminates between the two calibrations, not just
        # a coincidence of the new default.
        gi = GlucoseInsulinPathway(meals=[high_sugar_meal()])
        old_params = GLP1Parameters(glp1_basal=7.0, meal_response=0.035)
        glp1 = GLP1Pathway(old_params)
        composite = couple_with_gi(gi, glp1)
        system = ODESystem(composite)
        result = system.simulate(duration_s=300 * 60, meals=[high_sugar_meal()])
        glp1_pM = result.species["GLP1"] * 1e6
        assert glp1_pM.max() < 10.0   # nowhere near the real 41 pM peak

    def test_composite_preserves_basal_steady_state(self):
        # A neutral (non-agonist) GLP1Pathway at its own basal must not
        # perturb the GI basal steady state (dG/dI must remain ~0).
        composite = couple_with_gi(GlucoseInsulinPathway(), GLP1Pathway())
        derivs = composite.dydt(0.0, composite.y0)
        assert derivs[0] == pytest.approx(0.0, abs=1e-6)   # dG
        assert derivs[1] == pytest.approx(0.0, abs=1e-6)   # dI

    def test_semaglutide_reduces_glucose_peak_vs_baseline(self):
        """
        CLAUDE.md Priority 1 acceptance criterion: semaglutide must reduce
        the postprandial glucose peak by > 1 mM vs an unmedicated baseline,
        via the real ODE-coupled incretin pathway (not a shortcut modifier).
        """
        baseline = ODESystem(pathway=GlucoseInsulinPathway()) \
            .simulate_meal(high_sugar_meal(), duration_s=14400)

        composite = couple_with_gi(GlucoseInsulinPathway(), GLP1Pathway.semaglutide(1.0))
        treated = ODESystem(pathway=composite).simulate_meal(high_sugar_meal(), duration_s=14400)

        assert baseline.success and treated.success
        peak_base_mM, _ = baseline.glucose_peak()
        peak_treat_mM, _ = treated.glucose_peak()
        peak_base_mM /= 1000.0    # uM -> mM
        peak_treat_mM /= 1000.0
        assert peak_base_mM - peak_treat_mM > 1.0

    def test_liraglutide_reduces_glucose_peak_vs_baseline(self):
        baseline = ODESystem(pathway=GlucoseInsulinPathway()) \
            .simulate_meal(high_sugar_meal(), duration_s=14400)
        composite = couple_with_gi(GlucoseInsulinPathway(), GLP1Pathway.liraglutide(1.2))
        treated = ODESystem(pathway=composite).simulate_meal(high_sugar_meal(), duration_s=14400)

        peak_base, _ = baseline.glucose_peak()
        peak_treat, _ = treated.glucose_peak()
        assert peak_treat < peak_base

    def test_dpp4_inhibitor_reduces_glucose_peak_vs_baseline(self):
        baseline = ODESystem(pathway=GlucoseInsulinPathway()) \
            .simulate_meal(high_sugar_meal(), duration_s=14400)
        composite = couple_with_gi(GlucoseInsulinPathway(), GLP1Pathway.dpp4_inhibited())
        treated = ODESystem(pathway=composite).simulate_meal(high_sugar_meal(), duration_s=14400)

        peak_base, _ = baseline.glucose_peak()
        peak_treat, _ = treated.glucose_peak()
        assert peak_treat < peak_base

    def test_semaglutide_stronger_effect_than_liraglutide(self):
        base_meal = high_sugar_meal()
        sema_composite = couple_with_gi(GlucoseInsulinPathway(), GLP1Pathway.semaglutide(1.0))
        lira_composite = couple_with_gi(GlucoseInsulinPathway(), GLP1Pathway.liraglutide(1.2))

        peak_sema, _ = ODESystem(pathway=sema_composite).simulate_meal(base_meal, 14400).glucose_peak()
        peak_lira, _ = ODESystem(pathway=lira_composite).simulate_meal(base_meal, 14400).glucose_peak()
        assert peak_sema < peak_lira

    def test_watchdog_still_runs_for_glp1_composite(self):
        composite = couple_with_gi(GlucoseInsulinPathway(), GLP1Pathway.semaglutide(1.0))
        result = ODESystem(pathway=composite).simulate_meal(high_sugar_meal(), duration_s=14400)
        assert result.watchdog and len(result.watchdog) > 0
