# SAVE TO: packages\ode-engine\tests\test_meal_absorption_pathway.py
"""
Test suite for MealAbsorptionPathway -- the mass-conserving fix for
gastric_emptying_modifier (see CLAUDE.md "DONE -- Fixed
gastric_emptying_modifier mass-non-conservation").

Numerically confirmed before this fix: a real semaglutide scenario
(1.0mg, 75g meal, 5h) deleted ~42% of the meal's carbohydrate mass
because GlucoseInsulinPathway.dydt() multiplied total_Ra(meals, t)'s
OUTPUT by gastric_emptying_modifier directly. The integral
integral_0^inf(D*k*e^-kt dt) = D holds for any FIXED k, but scaling
the curve's output by any c != 1 scales the total integral by c too
-- this holds even for a constant modifier, not just a time-varying
one. MealAbsorptionPathway fixes this by tracking remaining unabsorbed
glucose mass per meal as its own ODE state, so gastric_emptying_modifier
scales the RATE the remaining mass decays at, not Ra's output --
mass-conserving for any modifier(t) trajectory.
"""
import numpy as np
import pytest

from biosim_ode.meal import MealInput
from biosim_ode.meal_absorption_pathway import MealAbsorptionPathway, meal_absorption_couplings
from biosim_ode.glp1_pathway import GLP1Pathway, GLP1Parameters
from biosim_ode.pathway import GlucoseInsulinPathway
from biosim_ode.coupling import CompositePathway
from biosim_ode.base_pathway import BasePathway


def _build_composite(meals, glp1_params):
    """Mirror SimulationRunner._build_pathway()'s 3-way GLP1 composition."""
    from biosim_ode.glp1_pathway import glp1_couplings

    gi = GlucoseInsulinPathway(meals=meals)
    glp1 = GLP1Pathway(glp1_params)
    absorption = MealAbsorptionPathway(meals)
    all_couplings = (
        glp1_couplings(gi.name, glp1.name)
        + meal_absorption_couplings(glp1.name, absorption.name, gi.name, meals)
    )
    return CompositePathway([gi, glp1, absorption], couplings=all_couplings)


class TestMealAbsorptionPathwayUnit:
    def test_implements_base_pathway(self):
        assert isinstance(MealAbsorptionPathway([]), BasePathway)

    def test_zero_meals_is_stateless(self):
        p = MealAbsorptionPathway([])
        assert p.species_names == []
        assert p.y0 == []
        assert p.dydt(0.0, []) == []
        assert p.output_signals(0.0, []) == {}

    def test_one_monoexponential_meal_one_species(self):
        meal = MealInput(carbs_g=75.0)
        p = MealAbsorptionPathway([meal])
        assert p.species_names == ["meal0_remaining_mM"]
        assert p.y0 == pytest.approx([meal.D_eff_mM])

    def test_bi_exponential_meal_tracked_as_stomach_and_gut(self):
        from biosim_ode.gut_absorption import GutAbsorptionModel

        mono = MealInput(carbs_g=75.0, time_min=0.0)
        bi = MealInput(carbs_g=50.0, time_min=30.0, absorption_model=GutAbsorptionModel.ogtt_calibrated_glucose())
        p = MealAbsorptionPathway([mono, bi])
        # Since 2026-11-15 the bi-exponential meal is tracked too, as two states.
        assert p.species_names == ["meal0_remaining_mM", "meal1_stomach_mM", "meal1_gut_mM"]

    def test_two_monoexponential_meals_two_species_correctly_ordered(self):
        meal0 = MealInput(carbs_g=75.0, time_min=0.0)
        meal1 = MealInput(carbs_g=50.0, time_min=120.0)
        p = MealAbsorptionPathway([meal0, meal1])
        assert p.species_names == ["meal0_remaining_mM", "meal1_remaining_mM"]
        assert p.y0 == pytest.approx([meal0.D_eff_mM, meal1.D_eff_mM])

    def test_remaining_pinned_before_meal_onset(self):
        meal = MealInput(carbs_g=75.0, time_min=60.0)
        p = MealAbsorptionPathway([meal])
        derivs = p.dydt(30.0, [meal.D_eff_mM], signals={"gastric_emptying_modifier": 0.5})
        assert derivs == pytest.approx([0.0])

    def test_remaining_decays_after_meal_onset(self):
        meal = MealInput(carbs_g=75.0, time_min=0.0)
        p = MealAbsorptionPathway([meal])
        derivs = p.dydt(10.0, [meal.D_eff_mM], signals={"gastric_emptying_modifier": 1.0})
        expected_Ra = meal.kabs * meal.D_eff_mM
        assert derivs == pytest.approx([-expected_Ra])

    def test_modifier_slows_decay_rate(self):
        meal = MealInput(carbs_g=75.0, time_min=0.0)
        p = MealAbsorptionPathway([meal])
        deriv_full = p.dydt(10.0, [meal.D_eff_mM], signals={"gastric_emptying_modifier": 1.0})[0]
        deriv_half = p.dydt(10.0, [meal.D_eff_mM], signals={"gastric_emptying_modifier": 0.5})[0]
        assert abs(deriv_half) == pytest.approx(abs(deriv_full) * 0.5)

    def test_default_modifier_is_neutral(self):
        meal = MealInput(carbs_g=75.0, time_min=0.0)
        p = MealAbsorptionPathway([meal])
        derivs_no_signals = p.dydt(10.0, [meal.D_eff_mM])
        derivs_neutral = p.dydt(10.0, [meal.D_eff_mM], signals={"gastric_emptying_modifier": 1.0})
        assert derivs_no_signals == pytest.approx(derivs_neutral)

    def test_output_signals_are_clamped_nonnegative(self):
        meal = MealInput(carbs_g=75.0, time_min=0.0)
        p = MealAbsorptionPathway([meal])
        signals = p.output_signals(10.0, [-5.0])
        assert signals["meal0_remaining_mM"] == pytest.approx(0.0)

    def test_output_signals_key_matches_original_meal_index(self):
        from biosim_ode.gut_absorption import GutAbsorptionModel

        mono = MealInput(carbs_g=75.0, time_min=0.0)
        bi = MealInput(carbs_g=50.0, time_min=30.0, absorption_model=GutAbsorptionModel.ogtt_calibrated_glucose())
        mono2 = MealInput(carbs_g=20.0, time_min=90.0)
        p = MealAbsorptionPathway([mono, bi, mono2])
        signals = p.output_signals(0.0, p.y0)
        assert set(signals.keys()) == {
            "meal0_remaining_mM", "meal1_stomach_mM", "meal1_gut_mM", "meal2_remaining_mM",
        }


class TestMealAbsorptionCouplings:
    def test_bi_exponential_meal_sends_its_gut_state(self):
        from biosim_ode.gut_absorption import GutAbsorptionModel

        mono = MealInput(carbs_g=75.0, time_min=0.0)
        bi = MealInput(carbs_g=50.0, time_min=30.0, absorption_model=GutAbsorptionModel.ogtt_calibrated_glucose())
        couplings = meal_absorption_couplings("glp1", "absorption", "gi", [mono, bi])
        target_signals = {c.signal_name for c in couplings if c.source_pathway == "absorption"}
        assert target_signals == {"meal0_remaining_mM", "meal1_gut_mM"}

    def test_includes_modifier_coupling_to_absorption_pathway(self):
        meal = MealInput(carbs_g=75.0)
        couplings = meal_absorption_couplings("glp1", "absorption", "gi", [meal])
        modifier_couplings = [
            c for c in couplings
            if c.source_pathway == "glp1" and c.target_pathway == "absorption"
        ]
        assert len(modifier_couplings) == 1
        assert modifier_couplings[0].signal_name == "gastric_emptying_modifier"


class TestMassConservationEndToEnd:
    """
    The core regression: this class RED-proofs the original bug (verified
    via `git stash` -- these tests fail against the pre-fix pathway.py,
    which multiplied Ra's output by the modifier directly).
    """

    def _integrate_remaining(self, meals, glp1_params, duration_min=300.0, n=2000):
        composite = _build_composite(meals, glp1_params)
        t_eval = np.linspace(0.0, duration_min, n)
        from scipy.integrate import solve_ivp
        sol = solve_ivp(
            fun=composite.dydt, t_span=(0.0, duration_min), y0=composite.y0,
            method="Radau", t_eval=t_eval, rtol=1e-6, atol=1e-8,
        )
        assert sol.success
        # MealAbsorptionPathway's species is the last one in the composite.
        remaining = sol.y[-1]
        return t_eval, remaining

    def test_semaglutide_conserves_meal_mass(self):
        meal = MealInput(carbs_g=75.0, time_min=0.0)
        glp1_params = GLP1Parameters.for_semaglutide(dose_mg=1.0)
        t_eval, remaining = self._integrate_remaining([meal], glp1_params)

        total_absorbed = remaining[0] - remaining[-1]
        pct_absorbed = 100.0 * total_absorbed / meal.D_eff_mM
        # Must be close to fully absorbed over 5h (small residual is
        # correctly-delayed, still-absorbing tail -- not deleted mass).
        #
        # Threshold lowered from 95.0 to 90.0 (2026-08-12 deep-search
        # plan, item A1): the real, citation-anchored GLP1_MEAL_RESPONSE
        # calibration (Orskov et al. 1994, PMID:8138058) makes endogenous
        # GLP-1 secretion genuinely stronger under a meal than the old,
        # uncited 0.035 constant did -- semaglutide's agonist_scale
        # multiplies this larger meal-response term too, so plasma GLP-1
        # (and therefore gastric-emptying suppression) is now genuinely
        # higher throughout the postprandial window, extending the
        # absorption tail slightly past the 5h window (94.8% at 5h, down
        # from >95% before this fix). Confirmed via direct integration
        # to a longer horizon that this is delay, not deletion: 99.7% at
        # 10h, 99.999% at 20h, ~100% at 40h -- see
        # test_semaglutide_meal_mass_eventually_reaches_completion below.
        assert pct_absorbed > 90.0, (
            f"only {pct_absorbed:.1f}% of meal mass absorbed over 5h -- "
            f"mass is being deleted, not delayed"
        )

    def test_semaglutide_meal_mass_eventually_reaches_completion(self):
        # Companion to the relaxed 5h threshold above: confirms directly
        # that the residual mass at 5h is genuinely a delayed tail, not
        # deleted mass -- integrating further shows it asymptotes to
        # 100%, not plateauing short of it.
        meal = MealInput(carbs_g=75.0, time_min=0.0)
        glp1_params = GLP1Parameters.for_semaglutide(dose_mg=1.0)
        t_eval, remaining = self._integrate_remaining(
            [meal], glp1_params, duration_min=1200.0, n=2000,
        )
        total_absorbed = remaining[0] - remaining[-1]
        pct_absorbed = 100.0 * total_absorbed / meal.D_eff_mM
        assert pct_absorbed > 99.9

    def test_liraglutide_conserves_meal_mass(self):
        meal = MealInput(carbs_g=75.0, time_min=0.0)
        glp1_params = GLP1Parameters.for_liraglutide(dose_mg=1.2)
        t_eval, remaining = self._integrate_remaining([meal], glp1_params)
        total_absorbed = remaining[0] - remaining[-1]
        pct_absorbed = 100.0 * total_absorbed / meal.D_eff_mM
        assert pct_absorbed > 95.0

    def test_dpp4_inhibitor_conserves_meal_mass(self):
        meal = MealInput(carbs_g=75.0, time_min=0.0)
        glp1_params = GLP1Parameters.dpp4_inhibited()
        t_eval, remaining = self._integrate_remaining([meal], glp1_params)
        total_absorbed = remaining[0] - remaining[-1]
        pct_absorbed = 100.0 * total_absorbed / meal.D_eff_mM
        assert pct_absorbed > 95.0

    def test_neutral_glp1_matches_legacy_closed_form_exactly(self):
        """
        When gastric_emptying_modifier is held at exactly 1.0 throughout
        (a constant, not time-varying, injected signal -- NOT a default
        GLP1Pathway, whose basal 7pM GLP-1 already produces a small
        nonzero emptying effect even without any agonist; that's a real,
        separate physiological detail, not something this test needs to
        reproduce), the state-tracked Ra must reduce to EXACTLY
        (solver-tolerance) the same trajectory as the original closed-form
        MealInput.Ra() -- proving this is a pure generalization, not a
        behavior change, for the modifier's neutral value.
        """
        meal = MealInput(carbs_g=75.0, time_min=0.0)

        gi = GlucoseInsulinPathway(meals=[meal])
        absorption = MealAbsorptionPathway([meal])

        def _composite_dydt_forced_neutral(t, y):
            gi_y, abs_y = y[:6], y[6:]
            abs_signals = {"gastric_emptying_modifier": 1.0}
            abs_derivs = absorption.dydt(t, abs_y, signals=abs_signals)
            abs_outputs = absorption.output_signals(t, abs_y)
            gi_signals = {"gastric_emptying_modifier": 1.0, **abs_outputs}
            gi_derivs = gi.dydt(t, gi_y, signals=gi_signals)
            return list(gi_derivs) + list(abs_derivs)

        from scipy.integrate import solve_ivp
        t_eval = np.linspace(0.0, 120.0, 500)
        y0 = list(gi.y0) + list(absorption.y0)
        sol = solve_ivp(
            fun=_composite_dydt_forced_neutral, t_span=(0.0, 120.0), y0=y0,
            method="Radau", t_eval=t_eval, rtol=1e-9, atol=1e-11,
        )
        assert sol.success
        G_state_tracked = sol.y[0]

        bare_gi = GlucoseInsulinPathway(meals=[meal])
        sol_bare = solve_ivp(
            fun=bare_gi.dydt, t_span=(0.0, 120.0), y0=bare_gi.y0,
            method="Radau", t_eval=t_eval, rtol=1e-9, atol=1e-11,
        )
        assert sol_bare.success
        G_closed_form = sol_bare.y[0]

        np.testing.assert_allclose(G_state_tracked, G_closed_form, rtol=1e-4)


class TestGIUsesStateTrackedRaWhenComposed:
    """
    Surgical regression test on the exact line changed in
    GlucoseInsulinPathway.dydt(): the earlier TestMassConservationEndToEnd
    tests above (correctly) proved MealAbsorptionPathway's OWN state
    conserves mass, but that alone is NOT sufficient proof GI actually
    reads and uses it -- confirmed the hard way: those tests kept passing
    even against a `git stash`-reverted, pre-fix pathway.py (the state
    pathway's decay is entirely self-contained and doesn't depend on GI at
    all). This class targets GI's dG computation directly and DOES fail
    against the pre-fix code (verified via git stash before this file was
    finalized).
    """

    def test_dG_reflects_injected_remaining_state_not_closed_form(self):
        meal = MealInput(carbs_g=75.0, time_min=0.0)
        gi = GlucoseInsulinPathway(meals=[meal])
        y = gi.y0
        t = 30.0

        # An injected remaining value deliberately DIFFERENT from what
        # meal.Ra(t)'s closed form implies at this t (a real, mass-
        # conserving trajectory would have decayed further under a
        # sub-1.0 modifier than the unmodified closed form does).
        injected_remaining = meal.D_eff_mM * 0.5
        signals = {
            "gastric_emptying_modifier": 0.7,
            "meal0_remaining_mM": injected_remaining,
        }
        derivs = gi.dydt(t, y, signals=signals)

        expected_Ra_state_based = meal.kabs * 0.7 * injected_remaining
        # dG at y0 (G==Gb): clearance term is 0, so dG == Ra exactly.
        assert derivs[0] == pytest.approx(expected_Ra_state_based)

        # And this must NOT equal the old closed-form-times-modifier
        # value the pre-fix code would have used instead (proves the
        # test actually discriminates the two formulas, not a
        # coincidental match).
        old_buggy_Ra = meal.Ra(t) * 0.7
        assert derivs[0] != pytest.approx(old_buggy_Ra)

    def test_dG_ignores_injected_state_for_meal_with_absorption_model(self):
        """
        A bi-exponential (GutAbsorptionModel) meal is tracked by its
        "meal{i}_gut_mM" state, never "meal{i}_remaining_mM" -- an injected
        remaining-key must not be read for it. With no gut key present, GI
        falls back to the closed-form Ra().
        """
        from biosim_ode.gut_absorption import GutAbsorptionModel

        meal = MealInput(
            carbs_g=75.0, time_min=0.0,
            absorption_model=GutAbsorptionModel.ogtt_calibrated_glucose(),
        )
        gi = GlucoseInsulinPathway(meals=[meal])
        y = gi.y0
        t = 30.0
        signals = {
            "gastric_emptying_modifier": 0.7,
            "meal0_remaining_mM": 999.0,  # a deliberately absurd value
        }
        derivs = gi.dydt(t, y, signals=signals)
        expected = meal.Ra(t) * 0.7
        assert derivs[0] == pytest.approx(expected)


class TestBackwardCompatibilityWithoutMealAbsorptionPathway:
    def test_gi_dydt_falls_back_to_closed_form_when_not_composed(self):
        """
        GlucoseInsulinPathway.dydt() must produce IDENTICAL results to
        before this fix whenever no MealAbsorptionPathway is composed
        (the overwhelming majority of scenarios -- confined blast radius).
        """
        meal = MealInput(carbs_g=75.0, time_min=0.0)
        gi = GlucoseInsulinPathway(meals=[meal])
        y = gi.y0
        derivs_no_signals = gi.dydt(10.0, y)
        # An unrelated signal (no gastric_emptying_modifier, no
        # meal0_remaining_mM key) must not change the Ra computation at all.
        derivs_with_unrelated_signal = gi.dydt(10.0, y, signals={"p3_modifier": 1.2})
        assert derivs_no_signals[0] == pytest.approx(derivs_with_unrelated_signal[0])

        # And the dG contribution from Ra must match the original closed
        # form exactly: total_clearance*(G-Gb) is 0 at y0=basal, so dG at
        # y0 is Ra alone.
        assert derivs_no_signals[0] == pytest.approx(meal.Ra(10.0))
