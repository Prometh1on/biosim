"""
Test suite for biosim_ode.

Coverage:
    - Kinetic functions: Hill activation/inhibition, Michaelis-Menten, edge cases
    - MealInput: Ra(t) shape, total glucose, presets
    - PathwayParameters: validation, basal state computation
    - GlucoseInsulinPathway: steady state (zero derivatives), physiological response
    - ODESystem: fasted steady state, post-meal glucose/insulin spike,
                 mTOR increase, AMPK decrease, return to baseline
    - SimulationResult: unit conversion, interpolation, peak detection
    - WatchdogResult: confidence scoring
"""

import math
import pytest
import numpy as np

from biosim_ode import (
    GlucoseInsulinPathway,
    MealInput,
    ODESystem,
    PathwayParameters,
    SimulationResult,
    SPECIES_NAMES,
    high_sugar_meal,
    hill_activation,
    hill_inhibition,
    light_snack,
    michaelis_menten,
    standard_meal,
    total_Ra,
)
from biosim_ode.kinetics import (
    approach_to_target,
    clamp,
    first_order,
    michaelis_menten_inhibited,
)
from biosim_ode.meal import DEFAULT_BODY_WEIGHT_KG


# ---------------------------------------------------------------------------
# Kinetics
# ---------------------------------------------------------------------------

class TestHillActivation:
    def test_at_ec50_gives_half(self):
        assert hill_activation(100.0, 100.0, 1.0) == pytest.approx(0.5)

    def test_at_zero_gives_zero(self):
        assert hill_activation(0.0, 100.0, 1.0) == pytest.approx(0.0)

    def test_at_very_high_approaches_one(self):
        assert hill_activation(1e9, 100.0, 1.0) == pytest.approx(1.0, abs=1e-5)

    def test_hill_2_is_more_sigmoidal(self):
        h1 = hill_activation(50.0, 100.0, 1.0)
        h2 = hill_activation(50.0, 100.0, 2.0)
        assert h2 < h1  # n=2 is more sigmoidal — less activation at sub-EC50

    def test_negative_activator_raises(self):
        with pytest.raises(ValueError):
            hill_activation(-1.0, 100.0, 1.0)

    def test_zero_ec50_raises(self):
        with pytest.raises(ValueError):
            hill_activation(50.0, 0.0, 1.0)

    def test_negative_ec50_raises(self):
        with pytest.raises(ValueError):
            hill_activation(50.0, -10.0, 1.0)

    def test_zero_hill_raises(self):
        with pytest.raises(ValueError):
            hill_activation(50.0, 100.0, 0.0)


class TestHillInhibition:
    def test_at_ic50_gives_half(self):
        assert hill_inhibition(100.0, 100.0, 1.0) == pytest.approx(0.5)

    def test_at_zero_inhibitor_gives_one(self):
        assert hill_inhibition(0.0, 100.0, 1.0) == pytest.approx(1.0)

    def test_at_high_inhibitor_approaches_zero(self):
        assert hill_inhibition(1e9, 100.0, 1.0) == pytest.approx(0.0, abs=1e-5)

    def test_complementary_to_activation(self):
        for conc in [10.0, 50.0, 100.0, 200.0]:
            assert hill_activation(conc, 100.0, 1.5) + hill_inhibition(conc, 100.0, 1.5) == pytest.approx(1.0)


class TestMichaelisMenten:
    def test_at_km_gives_half_vmax(self):
        assert michaelis_menten(100.0, 10.0, 100.0) == pytest.approx(5.0)

    def test_zero_substrate_gives_zero(self):
        assert michaelis_menten(0.0, 10.0, 5.0) == pytest.approx(0.0)

    def test_high_substrate_approaches_vmax(self):
        assert michaelis_menten(1e9, 10.0, 5.0) == pytest.approx(10.0, rel=1e-5)

    def test_negative_substrate_raises(self):
        with pytest.raises(ValueError):
            michaelis_menten(-1.0, 10.0, 5.0)

    def test_zero_km_raises(self):
        with pytest.raises(ValueError):
            michaelis_menten(5.0, 10.0, 0.0)

    def test_inhibited_slower_than_uninhibited(self):
        uninhibited = michaelis_menten(100.0, 10.0, 50.0)
        inhibited = michaelis_menten_inhibited(100.0, 10.0, 50.0, inhibitor=20.0, ki=10.0)
        assert inhibited < uninhibited


class TestMiscKinetics:
    def test_approach_to_target_positive(self):
        d = approach_to_target(current=0.2, target=0.8, rate=0.1)
        assert d == pytest.approx(0.06)  # 0.1 * (0.8 - 0.2)

    def test_approach_to_target_negative(self):
        d = approach_to_target(current=0.9, target=0.2, rate=0.1)
        assert d == pytest.approx(-0.07)

    def test_clamp_lower(self):
        assert clamp(-5.0, lower=0.0) == pytest.approx(0.0)

    def test_clamp_upper(self):
        assert clamp(1.5, lower=0.0, upper=1.0) == pytest.approx(1.0)

    def test_clamp_passthrough(self):
        assert clamp(0.5, lower=0.0, upper=1.0) == pytest.approx(0.5)

    def test_clamp_nan_maps_to_lower_not_upper(self):
        # Found during the meal-absorption mass-conservation fix's
        # security review: max(lower, min(upper, nan)) is silently
        # `upper` (NaN comparisons are always False), the opposite of a
        # safe floor -- every pathway calling clamp(y[i], lower=0.0) to
        # guard a diverging/malformed state relied on this incorrectly.
        import math
        assert clamp(float("nan")) == pytest.approx(0.0)
        assert clamp(float("nan"), lower=0.0, upper=1.0) == pytest.approx(0.0)
        assert not math.isnan(clamp(float("nan")))


# ---------------------------------------------------------------------------
# MealInput
# ---------------------------------------------------------------------------

class TestMealInput:
    def test_ra_zero_before_meal(self):
        meal = MealInput(carbs_g=75.0, time_min=30.0)
        assert meal.Ra(29.9) == pytest.approx(0.0)

    def test_ra_positive_at_meal_time(self):
        meal = MealInput(carbs_g=75.0, time_min=0.0)
        assert meal.Ra(0.0) > 0.0

    def test_ra_decays_exponentially(self):
        meal = MealInput(carbs_g=75.0, time_min=0.0, kabs=0.025)
        r0 = meal.Ra(0.0)
        r1 = meal.Ra(1.0)
        expected_ratio = math.exp(-0.025)
        assert r1 / r0 == pytest.approx(expected_ratio, rel=1e-4)

    def test_total_appearance_equals_d_eff(self):
        meal = MealInput(carbs_g=75.0, time_min=0.0, kabs=0.025)
        assert meal.total_glucose_appearance() == pytest.approx(meal.D_eff_mM)

    def test_d_eff_scales_with_carbs(self):
        m1 = MealInput(carbs_g=50.0)
        m2 = MealInput(carbs_g=100.0)
        assert m2.D_eff_mM == pytest.approx(m1.D_eff_mM * 2.0, rel=1e-4)

    def test_peak_ra(self):
        meal = MealInput(carbs_g=75.0, time_min=0.0, kabs=0.025)
        assert meal.peak_Ra() == pytest.approx(meal.D_eff_mM * 0.025)

    def test_negative_carbs_raises(self):
        with pytest.raises(ValueError, match="carbs_g"):
            MealInput(carbs_g=-10.0)

    def test_nan_carbs_raises(self):
        # Found during the meal-absorption mass-conservation fix's security
        # review: `carbs_g < 0` is silently False for NaN, so a NaN input
        # passed validation and reached D_eff_mM/MealAbsorptionPathway's y0
        # unvalidated.
        with pytest.raises(ValueError, match="carbs_g"):
            MealInput(carbs_g=float("nan"))

    def test_nan_body_weight_raises(self):
        with pytest.raises(ValueError, match="body_weight_kg"):
            MealInput(carbs_g=75.0, body_weight_kg=float("nan"))

    def test_protein_g_defaults_to_zero(self):
        meal = MealInput(carbs_g=75.0)
        assert meal.protein_g == 0.0

    def test_negative_protein_g_raises(self):
        with pytest.raises(ValueError, match="protein_g"):
            MealInput(carbs_g=75.0, protein_g=-1.0)

    def test_nan_protein_g_raises(self):
        with pytest.raises(ValueError, match="protein_g"):
            MealInput(carbs_g=75.0, protein_g=float("nan"))

    def test_total_appearance_equals_d_eff_at_reduced_kabs(self):
        # Regression guard for the acarbose milestone: reducing kabs
        # directly (as opposed to scaling Ra(t)'s OUTPUT, the GLP-1
        # mass-conservation bug fixed earlier in this codebase) must
        # preserve total carbohydrate mass exactly, for ANY kabs value --
        # confirmed here at the actual acarbose-calibrated 100mg scale
        # factor (0.35x) and at a smaller reduction (0.6x), not just the
        # unmodified default.
        from biosim_physiology.drug_interactions import ACARBOSE_100MG_KABS_SCALE
        default_meal = MealInput(carbs_g=75.0, time_min=0.0)
        for scale in (ACARBOSE_100MG_KABS_SCALE, 0.6, 0.1):
            reduced = MealInput(carbs_g=75.0, time_min=0.0, kabs=default_meal.kabs * scale)
            assert reduced.total_glucose_appearance() == pytest.approx(reduced.D_eff_mM)
            assert reduced.D_eff_mM == pytest.approx(default_meal.D_eff_mM)

    def test_absorption_model_defaults_to_monoexponential(self):
        # Backward compatibility: absorption_model=None (the default) must
        # leave Ra() behaviour completely unchanged from before this field
        # existed.
        meal = MealInput(carbs_g=75.0, time_min=0.0, kabs=0.025)
        assert meal.absorption_model is None
        assert meal.Ra(10.0) == pytest.approx(meal.D_eff_mM * 0.025 * math.exp(-0.025 * 10.0))

    def test_absorption_model_overrides_monoexponential(self):
        from biosim_ode.gut_absorption import GutAbsorptionModel
        model = GutAbsorptionModel.mixed_meal()
        meal = MealInput(carbs_g=75.0, time_min=0.0, absorption_model=model)
        # Ra(0) == 0 is the defining bi-exponential property (Ra(0)=0),
        # unlike the monoexponential model's instantaneous peak at t=0.
        assert meal.Ra(0.0) == pytest.approx(0.0, abs=1e-9)
        assert meal.Ra(5.0) == pytest.approx(model.Ra(5.0, meal.D_eff_mM))

    def test_absorption_model_before_meal_still_zero(self):
        from biosim_ode.gut_absorption import GutAbsorptionModel
        meal = MealInput(carbs_g=75.0, time_min=30.0, absorption_model=GutAbsorptionModel.mixed_meal())
        assert meal.Ra(10.0) == 0.0

    def test_absorption_model_mass_conserving(self):
        # Switching absorption shape must not change the total glucose dose
        # delivered -- both models integrate to D_eff_mM over infinite time.
        from biosim_ode.gut_absorption import GutAbsorptionModel
        model = GutAbsorptionModel.mixed_meal()
        meal = MealInput(carbs_g=75.0, time_min=0.0, absorption_model=model)
        import numpy as np
        times = np.linspace(0, 2000, 200_000)
        total = np.trapezoid([meal.Ra(t) for t in times], times)
        assert total == pytest.approx(meal.D_eff_mM, rel=1e-3)

    def test_total_ra_sums_meals(self):
        m1 = MealInput(carbs_g=50.0, time_min=0.0)
        m2 = MealInput(carbs_g=50.0, time_min=0.0)
        assert total_Ra([m1, m2], 0.0) == pytest.approx(m1.Ra(0.0) + m2.Ra(0.0))

    def test_total_ra_empty_is_zero(self):
        assert total_Ra([], 100.0) == pytest.approx(0.0)

    def test_high_sugar_meal_preset(self):
        meal = high_sugar_meal()
        assert meal.carbs_g == 75.0
        assert meal.time_min == 0.0

    def test_standard_meal_preset(self):
        meal = standard_meal()
        assert meal.carbs_g == 50.0

    def test_light_snack_preset(self):
        meal = light_snack()
        assert meal.carbs_g == 20.0


# ---------------------------------------------------------------------------
# PathwayParameters
# ---------------------------------------------------------------------------

class TestPathwayParameters:
    def test_default_construction(self):
        p = PathwayParameters()
        assert p.Gb == pytest.approx(5.0)
        assert p.Ib == pytest.approx(66.0)

    def test_invalid_gb_raises(self):
        with pytest.raises(ValueError):
            PathwayParameters(Gb=0.0)

    def test_invalid_ib_raises(self):
        with pytest.raises(ValueError):
            PathwayParameters(Ib=-1.0)

    def test_basal_state_length(self):
        p = PathwayParameters()
        y0 = p.basal_state()
        assert len(y0) == 6

    def test_basal_state_glucose_equals_gb(self):
        p = PathwayParameters()
        y0 = p.basal_state()
        assert y0[0] == pytest.approx(p.Gb)

    def test_basal_state_insulin_equals_ib(self):
        p = PathwayParameters()
        y0 = p.basal_state()
        assert y0[1] == pytest.approx(p.Ib)

    def test_basal_state_x_is_zero(self):
        p = PathwayParameters()
        y0 = p.basal_state()
        assert y0[2] == pytest.approx(0.0)

    def test_basal_state_mtor_in_range(self):
        p = PathwayParameters()
        y0 = p.basal_state()
        assert 0.0 < y0[4] < 1.0

    def test_basal_state_ampk_in_range(self):
        p = PathwayParameters()
        y0 = p.basal_state()
        assert 0.0 < y0[5] < 1.0

    def test_f_secretion_basal_positive(self):
        p = PathwayParameters()
        assert p.F_secretion_basal > 0.0

    def test_glucose_independent_secretion_default_is_neutral(self):
        assert PathwayParameters().glucose_independent_secretion == pytest.approx(0.0)

    def test_glucose_independent_secretion_negative_raises(self):
        with pytest.raises(ValueError):
            PathwayParameters(glucose_independent_secretion=-1.0)

    def test_glucose_independent_secretion_at_or_above_ke_ib_raises(self):
        # Hard safety bound: at ke*Ib, the analytic insulin fixed point at
        # zero glucose would reach 2x this instance's own basal Ib.
        p = PathwayParameters()
        with pytest.raises(ValueError):
            PathwayParameters(glucose_independent_secretion=p.ke * p.Ib)
        with pytest.raises(ValueError):
            PathwayParameters(glucose_independent_secretion=p.ke * p.Ib * 1.5)

    def test_glucose_independent_secretion_just_below_bound_is_accepted(self):
        p = PathwayParameters()
        PathwayParameters(glucose_independent_secretion=p.ke * p.Ib * 0.999)  # must not raise

    def test_kg_secretion_matches_literature_gsis_ec50(self):
        # Rebelato et al. 2018, Sci Rep, PMID:30166558: GSIS EC50 = 10.13 mM
        # absolute glucose in control (non-glucose-adapted) islets. This
        # model's KG_secretion is a Michaelis constant on EXCESS glucose
        # above G_threshold, so the literature-anchored value is
        # EC50 - G_threshold, not EC50 itself.
        p = PathwayParameters()
        assert p.KG_secretion == pytest.approx(10.13 - 4.5, abs=1e-6)

    def test_p3_matches_bergman_1981_own_reported_si(self):
        # Bergman RN, Phillips LS, Cobelli C. 1981, J Clin Invest,
        # PMID:7033284, Fig 6 caption: SI = 5.1 +/- 0.7 x10-4 min-1.uU-1.mL
        # for the combined normal-tolerant groups (n=11) -- read directly
        # from the paper (pypdf text extraction), not a later replication
        # study. p3 = SI * p2 / (pM per uU/mL, 6.945).
        p = PathwayParameters()
        si_per_uUmL = (p.p3 / p.p2) * 6.945
        assert si_per_uUmL * 1e4 == pytest.approx(5.1, rel=0.01)

    def test_first_phase_gain_positive(self):
        assert PathwayParameters().first_phase_gain > 0.0


class TestGlucoseIndependentSecretionFixedPoint:
    """
    Regression coverage for the sulfonylurea milestone: adding
    glucose_independent_secretion changes the analytic basal fixed point
    (I, X, IGF1, mTOR, AMPK all shift), which basal_state() must re-derive
    correctly -- not just leave hardcoded at the old, no-longer-correct
    values (the exact bug class in CLAUDE.md KNOWN BUGS FIXED #10-14).
    """

    def test_default_field_preserves_exact_prior_basal_state(self):
        # Byte-for-byte regression guard: every existing scenario in the
        # codebase relies on glucose_independent_secretion=0.0 producing
        # the EXACT same basal_state() as before this field existed.
        p = PathwayParameters()
        y0 = p.basal_state()
        assert y0[1] == pytest.approx(p.Ib, abs=1e-12)      # I
        assert y0[2] == pytest.approx(0.0, abs=1e-12)        # X
        assert y0[3] == pytest.approx(p.IGF1b, abs=1e-12)    # IGF1

    def test_nonzero_value_shifts_insulin_fixed_point_by_exact_formula(self):
        p = PathwayParameters()
        gis = 0.5 * p.ke * p.Ib
        p_su = PathwayParameters(glucose_independent_secretion=gis)
        y0 = p_su.basal_state()
        assert y0[1] == pytest.approx(p.Ib + gis / p.ke)

    def test_y0_is_a_genuine_fixed_point_when_nonzero(self):
        # The core regression this milestone's math had to get right:
        # not just I, but X/IGF1/mTOR/AMPK's own basal values must all be
        # re-derived from the new I_ss, or the composite y0 would not be
        # a genuine joint fixed point -- validate_y0() (used throughout
        # this codebase's own fixed-point discipline) is the standard,
        # tight-tolerance check for exactly this.
        gis = 0.5 * PathwayParameters().ke * PathwayParameters().Ib
        pathway = GlucoseInsulinPathway(params=PathwayParameters(glucose_independent_secretion=gis))
        violations = pathway.validate_y0(tolerance=1e-8)
        assert violations == []

    def test_stable_over_long_trajectory_at_basal_glucose(self):
        # A static y0 check at t=0 cannot by itself rule out a subtler
        # fixed-point bug (KNOWN BUGS FIXED #10-14's own lesson) --
        # integrate forward with zero meal input and confirm the state
        # doesn't drift.
        gis = 0.5 * PathwayParameters().ke * PathwayParameters().Ib
        pathway = GlucoseInsulinPathway(params=PathwayParameters(glucose_independent_secretion=gis))
        y = pathway.y0
        dt = 1.0
        for _ in range(600):
            d = pathway.dydt(0.0, y)
            y = [yi + di * dt for yi, di in zip(y, d)]
        y0 = pathway.y0
        for i, (state, initial) in enumerate(zip(y, y0)):
            assert state == pytest.approx(initial, rel=1e-3, abs=1e-6), f"species index {i} drifted"

    def test_addition_is_glucose_independent_unlike_delta_secretion(self):
        # The whole point of this mechanism: unlike delta_secretion (which
        # floors at -ke*Ib below G_threshold), the added contribution to
        # dI must be EXACTLY the same constant at low glucose as at basal
        # glucose -- proving the drug's secretagogue effect isn't silently
        # re-gated by the model's existing glucose-sensing machinery.
        gis = 0.5 * PathwayParameters().ke * PathwayParameters().Ib
        p0 = PathwayParameters()
        p1 = PathwayParameters(glucose_independent_secretion=gis)
        pathway0 = GlucoseInsulinPathway(params=p0)
        pathway1 = GlucoseInsulinPathway(params=p1)
        y_low = [2.0, p0.Ib, 0.0, p0.IGF1b, 0.0, 0.0]  # severe hypoglycaemia
        d0 = pathway0.dydt(0.0, y_low)
        d1 = pathway1.dydt(0.0, y_low)
        assert d1[1] - d0[1] == pytest.approx(gis)


class TestNearTotalBetaCellFailureFixedPoint:
    """
    MASTER_PLAN 1.2c: T1D (row 1.10) fundamentally needs a way to
    represent near-total beta-cell failure (near-zero endogenous
    insulin, with exogenous insulin -- see exogenous_insulin_pathway.py,
    MASTER_PLAN 1.2a/1.2b -- as the ONLY source of plasma insulin). This
    class verifies a real, previously-untested finding: no new
    PathwayParameters field or ODE mechanism is needed at all --
    `Ib` and `secretion_scale` are ALREADY free, independently-
    personalisable fields, and `basal_state()`'s existing fixed-point
    derivation (built during the sulfonylurea milestone specifically to
    generalize beyond a hardcoded Ib -- see
    TestGlucoseIndependentSecretionFixedPoint above) ALREADY correctly
    re-derives X/IGF1/mTOR/AMPK from whatever I_ss results, for ANY Ib
    value, not just population-average ones.

    first_phase_gain is set to 0.0 explicitly in every test here --
    unlike Ib/secretion_scale, it is a fixed constant with NO Ib
    dependence (Mari et al. 2002's dynamic secretion term is
    Ib-independent by construction), so a near-zero-Ib patient would
    otherwise still show a full-strength first-phase response, which is
    physiologically wrong for near-total beta-cell destruction -- this
    is a caller-responsibility today (the existing Brunzell-threshold
    personalizer gate zeroes first_phase_gain based on GLUCOSE LEVEL,
    a distinct T2D-specific phenomenon, not on beta-cell-function
    identity, so it must not be relied on here).

    This is a verification-only finding, not a new mechanism -- no
    PathwayParameters/pathway.py code was changed to make this pass. A
    literature-anchored T1D disease-state preset (specific Ib/
    secretion_scale values, wired into physiology/modifiers.py's
    disease-state system) remains MASTER_PLAN 1.10's own, separate,
    not-yet-started scope -- this class only confirms the ODE-engine
    mechanics that preset will need are already correct.
    """

    def test_low_ib_is_a_genuine_fixed_point(self):
        p = PathwayParameters(Ib=2.0, secretion_scale=0.01, first_phase_gain=0.0)
        pathway = GlucoseInsulinPathway(params=p)
        assert pathway.validate_y0(tolerance=1e-8) == []

    def test_basal_state_derives_every_species_from_the_low_ib_fixed_point(self):
        # Not just I -- X/IGF1/mTOR/AMPK must all be re-derived from the
        # actual (low) I_ss, not left at population-average defaults.
        p_healthy = PathwayParameters()
        p_t1d = PathwayParameters(Ib=2.0, secretion_scale=0.01, first_phase_gain=0.0)
        y0_healthy = p_healthy.basal_state()
        y0_t1d = p_t1d.basal_state()
        assert y0_t1d[1] == pytest.approx(2.0)          # I == Ib exactly (no secretagogue)
        assert y0_t1d[2] == pytest.approx(0.0, abs=1e-9)  # X: exactly 0 (I_ss == this instance's own Ib)
        # IGF1_ss's formula is RELATIVE to this instance's own Ib
        # (self.IGF1b * (1 + sensitivity*(I_ss/self.Ib - 1))) -- since T1D
        # has no secretagogue, I_ss == Ib exactly for BOTH instances, so
        # IGF1_ss stays at each instance's own IGF1b regardless of what Ib
        # actually is. Not a bug: X/mTOR/AMPK use I_ss's ABSOLUTE value
        # (via p3/EC50_mTOR/EC50_AMPK, fixed pM thresholds), so THOSE
        # correctly differ; IGF1 alone is intentionally relative-to-own-
        # basal in this model.
        assert y0_t1d[3] == pytest.approx(y0_healthy[3])  # IGF1: same (both at own IGF1b)
        assert y0_t1d[4] < y0_healthy[4]                 # mTOR: lower (absolute Hill function of I)
        assert y0_t1d[5] > y0_healthy[5]                 # AMPK: higher (absolute Hill-inhibited by I)

    def test_stable_over_long_fasted_trajectory(self):
        # Same KNOWN-BUGS-#10-14 discipline as the sulfonylurea class
        # above -- a static y0 check alone cannot rule out post-t=0
        # drift.
        p = PathwayParameters(Ib=2.0, secretion_scale=0.01, first_phase_gain=0.0)
        pathway = GlucoseInsulinPathway(params=p)
        y = pathway.y0
        dt = 1.0
        for _ in range(600):
            d = pathway.dydt(0.0, y)
            y = [yi + di * dt for yi, di in zip(y, d)]
        for i, (state, initial) in enumerate(zip(y, pathway.y0)):
            assert state == pytest.approx(initial, rel=1e-3, abs=1e-6), f"species index {i} drifted"

    def test_postprandial_endogenous_response_is_near_totally_absent(self):
        # The other half of the claim: not just a stable basal state, but
        # a genuinely near-absent RESPONSE to a meal (F_secretion_basal
        # scales proportionally with Ib, and secretion_scale/
        # first_phase_gain independently gate the dynamic terms) --
        # confirms this isn't just "a low but still meal-responsive"
        # basal insulin level.
        p_t1d = PathwayParameters(Ib=2.0, Gb=8.0, secretion_scale=0.01, first_phase_gain=0.0)
        p_healthy = PathwayParameters(Gb=8.0)
        meal = MealInput(carbs_g=75.0, time_min=0.0)
        pathway_t1d = GlucoseInsulinPathway(params=p_t1d)
        pathway_healthy = GlucoseInsulinPathway(params=p_healthy)
        system_t1d = ODESystem(pathway=pathway_t1d)
        system_healthy = ODESystem(pathway=pathway_healthy)
        # ODESystem.simulate()'s own meals= parameter overrides the
        # pathway's meals (defaults to none) -- must be passed explicitly
        # here, not just set on the pathway's own constructor.
        sim_t1d = system_t1d.simulate(duration_s=300 * 60, meals=[meal])
        sim_healthy = system_healthy.simulate(duration_s=300 * 60, meals=[meal])
        I_t1d = sim_t1d.species["I"] * 1e6
        I_healthy = sim_healthy.species["I"] * 1e6
        t1d_excursion = I_t1d.max() - I_t1d[0]
        healthy_excursion = I_healthy.max() - I_healthy[0]
        assert t1d_excursion < 0.05 * healthy_excursion


class TestFirstPhaseSecretion:
    """
    Mari et al. 2002, PMID:12388151: insulin secretion has a "dynamic"
    component proportional to the positive derivative of glucose
    (first-phase release), on top of the existing "static" delta_secretion
    (a function of glucose level). See PathwayParameters.first_phase_gain.
    """

    def test_zero_at_basal_steady_state(self):
        # dG=0 at basal (no meal, at analytic steady state) -> first-phase
        # term must contribute exactly zero, preserving every existing
        # basal-invariant test in this file.
        p = GlucoseInsulinPathway()
        derivs = p.steady_state_derivatives()
        assert derivs[1] == pytest.approx(0.0, abs=1e-9)   # dI

    def test_fires_when_glucose_rising(self):
        p = GlucoseInsulinPathway(meals=[MealInput(carbs_g=75.0, time_min=0.0)])
        y = p.y0
        y[0] = p.params.Gb + 1.0   # elevated but still rising (Ra > clearance early on)
        dI_with_meal = p.dydt(5.0, y)[1]

        p_no_gain = GlucoseInsulinPathway(
            params=PathwayParameters(first_phase_gain=0.0),
            meals=[MealInput(carbs_g=75.0, time_min=0.0)],
        )
        dI_without_first_phase = p_no_gain.dydt(5.0, y)[1]
        assert dI_with_meal > dI_without_first_phase

    def test_zero_when_glucose_falling(self):
        # An IV-bolus-like state: glucose elevated, no meal (no Ra), so
        # dG < 0 (pure clearance) -- first-phase must contribute nothing,
        # confirming it's mathematically inert during the whole Bergman
        # IVGTT trajectory (no Ra forcing term there at all).
        p = GlucoseInsulinPathway()
        y = p.y0
        y[0] = p.params.Gb + 5.0
        dI_with_gain = p.dydt(0.0, y)[1]

        p_no_gain = GlucoseInsulinPathway(params=PathwayParameters(first_phase_gain=0.0))
        dI_without = p_no_gain.dydt(0.0, y)[1]
        assert dI_with_gain == pytest.approx(dI_without)


class TestBergmanIVGTT:
    def test_y0_length(self):
        p = GlucoseInsulinPathway()
        y0 = p.bergman_ivgtt_y0()
        assert len(y0) == 6

    def test_only_glucose_elevated(self):
        # Only G should differ from the normal analytic basal state --
        # insulin, IGF1, mTOR, AMPK all start at their usual basal values
        # (the bolus is instantaneous and hasn't triggered any response yet).
        p = GlucoseInsulinPathway()
        y0 = p.bergman_ivgtt_y0()
        basal = p.y0
        assert y0[0] > basal[0]
        assert y0[1:] == pytest.approx(basal[1:])

    def test_larger_dose_higher_initial_glucose(self):
        p = GlucoseInsulinPathway()
        y0_small = p.bergman_ivgtt_y0(dose_g=10.0)
        y0_large = p.bergman_ivgtt_y0(dose_g=40.0)
        assert y0_large[0] > y0_small[0]

    def test_no_hepatic_extraction_applied(self):
        # IV glucose bypasses hepatic first-pass extraction entirely, unlike
        # MealInput.D_eff_mM (oral route) -- the full carbs_mmol/Vg_L bump
        # should appear, not a (1 - hepatic_extraction)-scaled-down value.
        from biosim_ode.meal import GLUCOSE_MOLECULAR_WEIGHT
        from biosim_ode.pathway import BERGMAN_IVGTT_VG_L_PER_KG
        p = GlucoseInsulinPathway()
        y0 = p.bergman_ivgtt_y0(dose_g=20.0, body_weight_kg=70.0)
        expected_bump = (20.0 * 1000.0 / GLUCOSE_MOLECULAR_WEIGHT) / (70.0 * BERGMAN_IVGTT_VG_L_PER_KG)
        assert y0[0] == pytest.approx(p.params.Gb + expected_bump)

    def test_ivgtt_vg_uses_literature_anchored_constant(self):
        # Regression: bergman_ivgtt_y0() used to default to its own
        # hardcoded, uncited 0.286 L/kg. Now uses BERGMAN_IVGTT_VG_L_PER_KG
        # (Mora-Gonzalez et al. 2024, PMID:38389141, healthy-subject IVGTT
        # VdGLUCOSE mean). This was initially kept as a separately-named
        # constant from MealInput.DEFAULT_VG_L_PER_KG specifically to avoid
        # touching oral-meal calibration without a dedicated re-validation
        # pass -- that pass has since happened (see CLAUDE.md KNOWN BUGS
        # FIXED #19) and both are now deliberately consolidated to the
        # SAME value (IV and oral routes dilute into the same single
        # glucose compartment, so there's no principled reason for them to
        # differ) via `BERGMAN_IVGTT_VG_L_PER_KG = DEFAULT_VG_L_PER_KG` in
        # pathway.py. If this ever needs to diverge again, that's a
        # deliberate future decision, not an accidental one -- update this
        # test's assertion, don't just delete it.
        from biosim_ode.pathway import BERGMAN_IVGTT_VG_L_PER_KG
        from biosim_ode.meal import DEFAULT_VG_L_PER_KG
        assert BERGMAN_IVGTT_VG_L_PER_KG == pytest.approx(0.230)
        assert BERGMAN_IVGTT_VG_L_PER_KG == DEFAULT_VG_L_PER_KG

    def test_y0_override_used_instead_of_pathway_y0(self):
        p = GlucoseInsulinPathway()
        y0 = p.bergman_ivgtt_y0()
        sys = ODESystem(pathway=p)
        result = sys.simulate(duration_s=60, meals=[], y0_override=y0, n_output_points=5)
        assert result.success
        assert result.species["G"][0] == pytest.approx(y0[0] * 1000.0)  # mM -> uM canonical

    def test_y0_override_none_preserves_default_behaviour(self):
        # Backward compatibility: omitting y0_override must behave exactly
        # as before this parameter existed.
        p = GlucoseInsulinPathway()
        sys = ODESystem(pathway=p)
        result = sys.simulate(duration_s=60, meals=[], n_output_points=5)
        assert result.success
        assert result.species["G"][0] == pytest.approx(p.params.Gb * 1000.0)


# ---------------------------------------------------------------------------
# Pathway steady state
# ---------------------------------------------------------------------------

class TestPathwaySteadyState:
    """
    At the analytic basal state with no meal, all derivatives must be ~0.
    This is the fundamental validation that the model is self-consistent.
    """

    def test_dG_dt_near_zero_at_basal(self):
        pathway = GlucoseInsulinPathway()
        derivs = pathway.steady_state_derivatives()
        assert abs(derivs[0]) < 1e-6, f"dG/dt = {derivs[0]} — not at steady state"

    def test_dI_dt_near_zero_at_basal(self):
        pathway = GlucoseInsulinPathway()
        derivs = pathway.steady_state_derivatives()
        assert abs(derivs[1]) < 1e-4, f"dI/dt = {derivs[1]} — not at steady state"

    def test_dX_dt_near_zero_at_basal(self):
        pathway = GlucoseInsulinPathway()
        derivs = pathway.steady_state_derivatives()
        assert abs(derivs[2]) < 1e-8, f"dX/dt = {derivs[2]} — not at steady state"

    def test_dIGF1_dt_near_zero_at_basal(self):
        pathway = GlucoseInsulinPathway()
        derivs = pathway.steady_state_derivatives()
        assert abs(derivs[3]) < 1e-6, f"dIGF1/dt = {derivs[3]} — not at steady state"

    def test_dmTOR_dt_near_zero_at_basal(self):
        pathway = GlucoseInsulinPathway()
        derivs = pathway.steady_state_derivatives()
        assert abs(derivs[4]) < 1e-6, f"dmTOR/dt = {derivs[4]} — not at steady state"

    def test_dAMPK_dt_near_zero_at_basal(self):
        pathway = GlucoseInsulinPathway()
        derivs = pathway.steady_state_derivatives()
        assert abs(derivs[5]) < 1e-6, f"dAMPK/dt = {derivs[5]} — not at steady state"


# ---------------------------------------------------------------------------
# ODESystem — fasted simulation (must stay near steady state)
# ---------------------------------------------------------------------------

class TestODESystemFasted:
    def test_fasted_simulation_succeeds(self):
        system = ODESystem()
        result = system.simulate_fasted(duration_s=3600)
        assert result.success

    def test_fasted_glucose_stays_near_baseline(self):
        system = ODESystem()
        result = system.simulate_fasted(duration_s=3600)
        G = result.species["G"]
        G_baseline_uM = 5.0 * 1000.0  # 5 mM → 5000 μM
        assert np.max(np.abs(G - G_baseline_uM)) < 50.0, \
            "Fasted glucose drifted more than 50 μM from baseline"

    def test_fasted_insulin_stays_near_baseline(self):
        system = ODESystem()
        result = system.simulate_fasted(duration_s=3600)
        I = result.species["I"]
        I_baseline_uM = 66.0 * 1e-6  # 66 pM → μM
        assert np.max(np.abs(I - I_baseline_uM)) < I_baseline_uM * 0.05, \
            "Fasted insulin drifted more than 5% from baseline"

    def test_times_in_seconds(self):
        system = ODESystem()
        result = system.simulate_fasted(duration_s=3600)
        assert result.times_s[-1] == pytest.approx(3600.0, rel=0.01)

    def test_species_all_present(self):
        system = ODESystem()
        result = system.simulate_fasted(duration_s=3600)
        for name in SPECIES_NAMES:
            assert name in result.species


# ---------------------------------------------------------------------------
# ODESystem — postprandial simulation (physiological validation)
# ---------------------------------------------------------------------------

class TestODESystemMeal:
    """
    After a high-sugar meal, the model must show:
    1. Glucose rises above baseline
    2. Insulin rises above baseline
    3. mTOR activity increases
    4. AMPK activity decreases
    5. Glucose returns toward baseline within 4 hours

    These are the physiological benchmarks from Sprint 0.
    Exact magnitudes are validated in Sprint 12 (clinical validation).
    """

    @pytest.fixture(scope="class")
    def meal_result(self):
        system = ODESystem()
        return system.simulate_meal(high_sugar_meal(), duration_s=14400)

    def test_simulation_succeeds(self, meal_result):
        assert meal_result.success

    def test_glucose_rises_above_baseline(self, meal_result):
        G = meal_result.species["G"]
        G_baseline_uM = 5.0 * 1000.0
        assert np.max(G) > G_baseline_uM * 1.1, \
            "Glucose did not rise 10% above baseline after high-sugar meal"

    def test_glucose_peak_is_physiological(self, meal_result):
        peak_G, peak_t = meal_result.glucose_peak()
        # Physiological peak: 6.5-12 mM = 6500-12000 μM for 75g OGTT
        assert 5500 < peak_G < 20000, \
            f"Glucose peak {peak_G:.0f} μM outside physiological range"

    def test_glucose_peak_timing(self, meal_result):
        _, peak_t_s = meal_result.glucose_peak()
        # Peak should be within 20-90 minutes
        assert 20 * 60 < peak_t_s < 90 * 60, \
            f"Glucose peak at {peak_t_s/60:.1f} min — outside expected 20-90 min"

    def test_insulin_rises_above_baseline(self, meal_result):
        I = meal_result.species["I"]
        I_baseline_uM = 66.0 * 1e-6
        assert np.max(I) > I_baseline_uM * 1.5, \
            "Insulin did not rise 50% above baseline after high-sugar meal"

    def test_mtor_increases_post_meal(self, meal_result):
        mTOR = meal_result.species["mTOR"]
        mTOR_initial = mTOR[0]
        assert np.max(mTOR) > mTOR_initial * 1.1, \
            "mTOR did not increase after high-sugar meal"

    def test_ampk_decreases_post_meal(self, meal_result):
        AMPK = meal_result.species["AMPK"]
        AMPK_initial = AMPK[0]
        assert np.min(AMPK) < AMPK_initial * 0.9, \
            "AMPK did not decrease after high-sugar meal"

    def test_glucose_returns_toward_baseline(self, meal_result):
        G = meal_result.species["G"]
        G_baseline_uM = 5.0 * 1000.0
        G_final = G[-1]
        G_peak = np.max(G)
        excursion = G_peak - G_baseline_uM
        remaining_excursion = G_final - G_baseline_uM
        # Should return more than halfway to baseline in 4 hours
        assert remaining_excursion < excursion * 0.5, \
            "Glucose did not return 50% toward baseline within 4 hours"

    def test_larger_meal_higher_glucose_peak(self):
        system = ODESystem()
        r_small = system.simulate_meal(light_snack(), duration_s=14400)
        r_large = system.simulate_meal(high_sugar_meal(), duration_s=14400)
        peak_small, _ = r_small.glucose_peak()
        peak_large, _ = r_large.glucose_peak()
        assert peak_large > peak_small, \
            "Larger meal should produce higher glucose peak"

    def test_delayed_meal_shifts_peak_timing(self):
        system = ODESystem()
        r_immediate = system.simulate_meal(high_sugar_meal(time_min=0.0), duration_s=14400)
        r_delayed = system.simulate_meal(high_sugar_meal(time_min=30.0), duration_s=14400)
        _, t_immediate = r_immediate.glucose_peak()
        _, t_delayed = r_delayed.glucose_peak()
        assert t_delayed > t_immediate, \
            "Delayed meal should shift peak timing later"


# ---------------------------------------------------------------------------
# SimulationResult utility methods
# ---------------------------------------------------------------------------

class TestSimulationResult:
    @pytest.fixture
    def result(self):
        system = ODESystem()
        return system.simulate_meal(high_sugar_meal(), duration_s=7200)

    def test_duration_s(self, result):
        assert result.duration_s == pytest.approx(7200.0, rel=0.01)

    def test_duration_min(self, result):
        assert result.duration_min == pytest.approx(120.0, rel=0.01)

    def test_at_time_returns_dict(self, result):
        state = result.at_time(1800.0)
        assert "G" in state
        assert "I" in state
        assert "mTOR" in state

    def test_glucose_peak_returns_tuple(self, result):
        peak_G, peak_t = result.glucose_peak()
        assert peak_G > 0
        assert peak_t > 0

    def test_insulin_peak_returns_tuple(self, result):
        peak_I, peak_t = result.insulin_peak()
        assert peak_I > 0
        assert peak_t > 0


# ---------------------------------------------------------------------------
# Watchdog
# ---------------------------------------------------------------------------

class TestWatchdog:
    def test_watchdog_present_in_result(self):
        system = ODESystem()
        result = system.simulate_meal(high_sugar_meal(), duration_s=7200)
        assert result.watchdog is not None
        assert len(result.watchdog) > 0

    def test_fasted_confidence_is_high(self):
        system = ODESystem()
        result = system.simulate_fasted(duration_s=3600)
        confidence = system.overall_confidence(result)
        assert confidence >= 0.8, \
            f"Fasted simulation confidence {confidence:.2f} should be high"

    def test_watchdog_score_between_0_and_1(self):
        system = ODESystem()
        result = system.simulate_fasted(duration_s=3600)
        for w in result.watchdog:
            assert 0.0 <= w.confidence_score <= 1.0

    def test_elevated_glucose_without_insulin_response_flagged(self):
        # Regression: _run_watchdog() used to only check two rules driven
        # by elevated insulin (I_elevated -> mTOR_active,
        # I_elevated -> AMPK_suppressed) -- it never checked the reverse
        # direction, G_elevated -> I_elevated (a broken beta-cell
        # response: glucose is high but insulin secretion never rose to
        # meet it). boolean_pathway.py's WatchdogEngine implements this
        # third rule (severity 0.30, the highest of the three) but was
        # never wired into production. Synthetic state: G > 8.0 mM,
        # I < 150 pM, mTOR/AMPK at a consistent baseline (no violation
        # from the other two rules) -- isolates this one rule.
        system = ODESystem()
        pathway = GlucoseInsulinPathway()
        t_min = np.array([0.0])
        # Species order: G, I, X, IGF1, mTOR, AMPK
        y_internal = np.array([[9.0], [50.0], [0.0], [10.0], [0.1], [0.9]])
        results = system._run_watchdog(t_min, y_internal, pathway)
        assert len(results) == 1
        assert any(
            "beta-cell" in d.lower() or "insulin not elevated" in d.lower()
            for d in results[0].disagreements
        )
        assert results[0].confidence_score == pytest.approx(0.70)


class TestPublicAPICompleteness:
    """
    Regression guard for the housekeeping gap found during this
    project's coverage audit: __init__.py's __all__ used to omit
    CortisolPathway/HepGlucPathway/NEFAPathway/GLP1Pathway (and, once
    checked more closely, every extension pathway's own Parameters
    class and couple_with_gi/*_couplings convenience functions) despite
    all of them being imported at module level. Cosmetic (only affects
    `from biosim_ode import *`, which nothing in this codebase's own
    41 import sites actually uses), but real -- fixed by completing
    __all__ to match every name genuinely imported at module scope.
    """

    def test_every_name_in_all_actually_resolves(self):
        import biosim_ode
        missing = [name for name in biosim_ode.__all__ if not hasattr(biosim_ode, name)]
        assert missing == []

    def test_star_import_exposes_every_extension_pathway_class(self):
        ns: dict = {}
        exec("from biosim_ode import *", ns)
        for name in (
            "CortisolPathway", "CortisolParameters",
            "HepGlucPathway", "EGPParameters",
            "NEFAPathway", "NEFAParameters",
            "GLP1Pathway", "GLP1Parameters",
            "RenalPathway", "RenalClearanceParameters",
            "ExercisePathway", "ExogenousInsulinPathway", "ExogenousInsulinParameters",
        ):
            assert name in ns, f"{name} missing from `from biosim_ode import *`"

    def test_star_import_exposes_every_couple_with_gi_convenience_function(self):
        ns: dict = {}
        exec("from biosim_ode import *", ns)
        for name in (
            "couple_cortisol_with_gi", "cortisol_couplings",
            "couple_hepatic_with_gi", "hepatic_couplings",
            "couple_nefa_with_gi", "nefa_couplings",
            "couple_glp1_with_gi", "glp1_couplings",
            "couple_exercise_with_gi", "exercise_couplings",
            "couple_renal_with_gi", "renal_couplings",
            "couple_exogenous_insulin_with_gi", "exogenous_insulin_couplings",
        ):
            assert name in ns, f"{name} missing from `from biosim_ode import *`"


class TestReadOutsDoNotFeedBack:
    """IGF-1, mTOR and AMPK are insulin-driven read-outs: changing them must
    not change dG, dI or dX (README/CLAUDE.md, CV wording)."""

    def test_changing_readouts_leaves_glucose_insulin_x_unchanged(self):
        from biosim_ode import GlucoseInsulinPathway
        pathway = GlucoseInsulinPathway()
        y = list(pathway.y0)
        y[0], y[1], y[2] = 8.0, 200.0, 0.002   # a post-meal-like state
        perturbed = list(y)
        perturbed[3], perturbed[4], perturbed[5] = 40.0, 0.95, 0.05
        base = pathway.dydt(10.0, y)
        moved = pathway.dydt(10.0, perturbed)
        assert moved[:3] == base[:3]
        assert moved[3:] != base[3:]   # the read-outs themselves do respond
