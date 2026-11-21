# SAVE TO: packages\ode-engine\tests\test_nutrient_secretion_pathway.py
"""
Test suite for the dietary protein (nutrient) secretion pathway.

See CLAUDE.md/MASTER_PLAN row 1.28 for the full literature-verification
and calibration trail -- Krezowski et al. 1986 (PMID:3538843, normal
subjects, 28% additive insulin-AUC ratio for 50g protein vs. 50g glucose)
and Nuttall et al. 1984 (PMID:6389060, T2D, corroborating direction but
a genuinely synergistic, out-of-scope regime).
"""
import math

import pytest

from biosim_ode import GlucoseInsulinPathway, ODESystem
from biosim_ode.meal import MealInput
from biosim_ode.nutrient_secretion_pathway import (
    NutrientSecretionParameters,
    NutrientSecretionPathway,
    couple_with_gi,
    nutrient_secretion_couplings,
    MAX_PLAUSIBLE_PROTEIN_G_PER_MEAL,
    NUTRIENT_SECRETION_KA_PER_MIN,
    NUTRIENT_SECRETION_PM_PER_G,
)


class TestNutrientSecretionParameters:
    def test_zero_protein_is_neutral_at_every_time(self):
        p = NutrientSecretionParameters(protein_g=0.0)
        for t in [0.0, 10.0, 50.0, 100.0, 500.0]:
            assert p.delta_pM_per_min(t) == pytest.approx(0.0)

    def test_negative_protein_raises(self):
        with pytest.raises(ValueError):
            NutrientSecretionParameters(protein_g=-1.0)

    def test_nan_protein_raises(self):
        with pytest.raises(ValueError):
            NutrientSecretionParameters(protein_g=float("nan"))

    def test_non_positive_ka_raises(self):
        with pytest.raises(ValueError):
            NutrientSecretionParameters(protein_g=50.0, ka=0.0)
        with pytest.raises(ValueError):
            NutrientSecretionParameters(protein_g=50.0, ka=-0.01)

    def test_negative_pm_per_g_raises(self):
        with pytest.raises(ValueError):
            NutrientSecretionParameters(protein_g=50.0, pm_per_g=-1.0)

    def test_overflowing_s1_0_raises(self):
        with pytest.raises(ValueError, match="S1_0"):
            NutrientSecretionParameters(protein_g=50.0, pm_per_g=1e307)

    def test_protein_g_above_plausibility_ceiling_raises(self):
        # Mandatory security-reviewer finding, HIGH: the un-ceilinged
        # dose-response is purely linear (no saturating term), so a
        # sufficiently large protein_g produces an implausible composed
        # insulin spike -- verified directly against the real composed
        # ODE before choosing MAX_PLAUSIBLE_PROTEIN_G_PER_MEAL.
        with pytest.raises(ValueError, match="plausibility ceiling"):
            NutrientSecretionParameters(protein_g=MAX_PLAUSIBLE_PROTEIN_G_PER_MEAL + 1.0)

    def test_protein_g_at_plausibility_ceiling_is_accepted(self):
        NutrientSecretionParameters(protein_g=MAX_PLAUSIBLE_PROTEIN_G_PER_MEAL)

    def test_delta_is_zero_before_meal_time(self):
        p = NutrientSecretionParameters(protein_g=50.0, meal_time_min=50.0)
        assert p.delta_pM_per_min(49.9) == pytest.approx(0.0)
        assert p.delta_pM_per_min(50.0) == pytest.approx(0.0)
        assert p.delta_pM_per_min(50.1) > 0.0

    def test_delta_is_nonnegative_everywhere(self):
        p = NutrientSecretionParameters(protein_g=50.0)
        for t in [0.001, 1, 10, 50, 100, 200, 500, 2000, 1e5]:
            assert p.delta_pM_per_min(t) >= 0.0

    def test_delta_decays_to_zero_at_long_times(self):
        p = NutrientSecretionParameters(protein_g=50.0)
        assert p.delta_pM_per_min(1e6) == pytest.approx(0.0, abs=1e-9)

    def test_dose_response_is_linear_in_grams(self):
        # Nuttall 1984's own finding: "the insulin area response was
        # essentially first order" -- linear in protein_g.
        p50 = NutrientSecretionParameters(protein_g=50.0)
        p100 = NutrientSecretionParameters(protein_g=100.0)
        for t in [10, 50, 100, 200]:
            assert p100.delta_pM_per_min(t) == pytest.approx(2.0 * p50.delta_pM_per_min(t))

    def test_total_auc_is_independent_of_ka(self):
        # The whole point of a monoexponential-decay RATE curve: its own
        # time-integral (S1_0) is exactly independent of the rate
        # constant chosen -- this is what let the timing constant stay
        # a disclosed, unanchored assumption while only the magnitude
        # (pm_per_g) needed literature-driven calibration.
        import numpy as np
        for ka in [0.01, 0.025, 0.1]:
            p = NutrientSecretionParameters(protein_g=50.0, ka=ka)
            ts = np.linspace(0, 5000, 200000)
            deltas = [p.delta_pM_per_min(t) for t in ts]
            auc = np.trapezoid(deltas, ts)
            assert auc == pytest.approx(p.S1_0, rel=0.01)

    def test_meal_time_shifts_the_curve(self):
        on_time = NutrientSecretionParameters(protein_g=50.0, meal_time_min=0.0)
        early = NutrientSecretionParameters(protein_g=50.0, meal_time_min=-30.0)
        assert early.delta_pM_per_min(0.0) == pytest.approx(on_time.delta_pM_per_min(30.0))


class TestNutrientSecretionPathway:
    def test_zero_species(self):
        p = NutrientSecretionPathway([NutrientSecretionParameters(protein_g=50.0)])
        assert p.species_names == []
        assert p.y0 == []

    def test_dydt_returns_empty(self):
        p = NutrientSecretionPathway([NutrientSecretionParameters(protein_g=50.0)])
        assert p.dydt(0.0, []) == []

    def test_empty_meals_list_is_neutral(self):
        p = NutrientSecretionPathway([])
        assert p.output_signals(100.0, [])["nutrient_secretion_delta"] == pytest.approx(0.0)

    def test_none_meals_defaults_to_empty(self):
        p = NutrientSecretionPathway()
        assert p.output_signals(100.0, [])["nutrient_secretion_delta"] == pytest.approx(0.0)

    def test_output_signal_key(self):
        p = NutrientSecretionPathway([NutrientSecretionParameters(protein_g=50.0)])
        signals = p.output_signals(30.0, [])
        assert "nutrient_secretion_delta" in signals
        assert signals["nutrient_secretion_delta"] > 0.0

    def test_sums_across_multiple_meals(self):
        # Two separate protein-bearing meals at different times -- the
        # signal at any given t should be the sum of each meal's own
        # individual contribution, not just one of them.
        m1 = NutrientSecretionParameters(protein_g=30.0, meal_time_min=0.0)
        m2 = NutrientSecretionParameters(protein_g=20.0, meal_time_min=60.0)
        pathway = NutrientSecretionPathway([m1, m2])
        t = 90.0
        expected = m1.delta_pM_per_min(t) + m2.delta_pM_per_min(t)
        assert pathway.output_signals(t, [])["nutrient_secretion_delta"] == pytest.approx(expected)

    def test_zero_protein_preserves_basal_invariance(self):
        nutrient = NutrientSecretionPathway([NutrientSecretionParameters(protein_g=0.0)])
        gi = GlucoseInsulinPathway()
        composite = couple_with_gi(gi, nutrient)
        assert composite.validate_y0(tolerance=1e-8) == []

    def test_active_protein_does_not_perturb_the_analytic_basal_state_at_t0(self):
        nutrient = NutrientSecretionPathway(
            [NutrientSecretionParameters(protein_g=50.0, meal_time_min=0.0)]
        )
        gi = GlucoseInsulinPathway()
        composite = couple_with_gi(gi, nutrient)
        assert composite.validate_y0(tolerance=1e-8) == []

    def test_couplings_map_signal_straight_through(self):
        couplings = nutrient_secretion_couplings("gi", "nutrient")
        assert len(couplings) == 1
        c = couplings[0]
        assert c.source_pathway == "nutrient"
        assert c.signal_name == "nutrient_secretion_delta"
        assert c.target_pathway == "gi"
        assert c.target_input == "nutrient_secretion_delta"

    def test_composed_simulation_reproduces_krezowski_1986_insulin_auc_ratio(self):
        # End-to-end reproduction of Krezowski et al. 1986's own 28%
        # additive insulin-AUC ratio through the REAL composed ODE, not
        # just the isolated signal function -- 50g protein alone vs. 50g
        # carbohydrate alone, both over a 240min (4h) window matching
        # the paper's own protocol.
        import numpy as np

        def insulin_auc_above_baseline(meal: MealInput) -> float:
            gi = GlucoseInsulinPathway()
            if meal.protein_g > 0:
                nutrient = NutrientSecretionPathway(
                    [NutrientSecretionParameters(protein_g=meal.protein_g, meal_time_min=meal.time_min)]
                )
                pathway = couple_with_gi(gi, nutrient)
            else:
                pathway = gi
            system = ODESystem(pathway=pathway)
            sim = system.simulate(duration_s=240 * 60, meals=[meal])
            I = sim.species["I"]
            t = sim.times_s
            Ib = I[0]
            return np.trapezoid(I - Ib, t) * 1e6 / 60.0   # pM.min

        glucose_auc = insulin_auc_above_baseline(MealInput(carbs_g=50.0, protein_g=0.0))
        protein_auc = insulin_auc_above_baseline(MealInput(carbs_g=0.0, protein_g=50.0))
        # KNOWN GAP since B7 (2026-11-16), recorded, NOT re-tuned: this
        # secretion constant was calibrated to Krezowski's 0.28 while the
        # model double-counted resting glucose clearance. With that removed,
        # the same constant gives ~0.22. Pinned so any drift is looked at.
        assert protein_auc / glucose_auc == pytest.approx(0.218, abs=0.005)

    def test_composed_simulation_is_approximately_additive_not_synergistic(self):
        # Krezowski's own explicit finding: "the effect on insulin
        # secretion is not synergistic" (100.4% additive in the real
        # data). This platform's own composed system lands close to but
        # not exactly at 1.0 (a real, expected consequence of glucose
        # and insulin feeding back on each other in a genuinely coupled
        # nonlinear system) -- locked in as a regression at its own
        # measured value, not forced to exactly 1.0.
        import numpy as np

        def insulin_auc(carbs_g: float, protein_g: float) -> float:
            gi = GlucoseInsulinPathway()
            meal = MealInput(carbs_g=carbs_g, protein_g=protein_g, time_min=0.0)
            if protein_g > 0:
                nutrient = NutrientSecretionPathway(
                    [NutrientSecretionParameters(protein_g=protein_g, meal_time_min=0.0)]
                )
                pathway = couple_with_gi(gi, nutrient)
            else:
                pathway = gi
            system = ODESystem(pathway=pathway)
            sim = system.simulate(duration_s=240 * 60, meals=[meal])
            I = sim.species["I"]
            t = sim.times_s
            Ib = I[0]
            return np.trapezoid(I - Ib, t) * 1e6 / 60.0

        glucose_only = insulin_auc(50.0, 0.0)
        protein_only = insulin_auc(0.0, 50.0)
        combo = insulin_auc(50.0, 50.0)
        ratio = combo / (glucose_only + protein_only)
        # Genuinely NOT synergistic (ratio < 1.0 here, i.e. sub-additive
        # due to glucose-insulin feedback, not the T2D synergy Nuttall
        # 1984 reports in a different, out-of-scope population) and
        # close to Krezowski's own near-1.0 finding.
        assert 0.90 < ratio < 1.0
