# SAVE TO: packages\ode-engine\tests\test_renal_pathway.py
"""
Test suite for the SGLT2-inhibitor renal glucose clearance pathway.
"""
import pytest
from biosim_ode.renal_pathway import (
    RenalClearanceParameters,
    RenalPathway,
    GLUCOSE_REFERENCE_MM,
    SGLT2I_BETA_CLASS_DEFAULT,
)


class TestRenalClearanceParameters:
    def test_default_construction_is_neutral(self):
        p = RenalClearanceParameters()
        assert p.beta == pytest.approx(0.0)

    def test_beta_out_of_range_raises(self):
        with pytest.raises(ValueError):
            RenalClearanceParameters(beta=1.0)
        with pytest.raises(ValueError):
            RenalClearanceParameters(beta=-0.1)

    def test_non_positive_gfr_raises(self):
        with pytest.raises(ValueError):
            RenalClearanceParameters(gfr_ml_per_min=0.0)

    def test_non_positive_tm_raises(self):
        with pytest.raises(ValueError):
            RenalClearanceParameters(tm_mg_per_min=0.0)

    def test_non_positive_splay_raises(self):
        with pytest.raises(ValueError):
            RenalClearanceParameters(splay_p=0.0)

    def test_no_drug_delta_is_exactly_zero_across_glucose_range(self):
        # The core neutrality guarantee this pathway is designed around
        # (stronger than hepatic_egp_delta's basal-only zero): with beta=0
        # (no SGLT2i drug), the exported delta must be exactly zero at
        # EVERY glucose level, not just at the reference/basal value --
        # since the real mechanism is continuously active whenever
        # present, a basal-only guarantee would not be enough to prove no
        # double-counting against baseline glucose clearance (p1; the old
        # extra_g_clearance term is 0 since B7).
        p = RenalClearanceParameters(beta=0.0)
        for g in [0.5, 2.0, 5.0, 8.0, 10.0, 15.0, 20.0, 30.0]:
            assert p.target_delta_mM_per_min(g) == pytest.approx(0.0, abs=1e-12)

    def test_drug_delta_is_nonzero_and_negative_at_basal_glucose(self):
        # Unlike hepatic_egp_delta, the whole point of this mechanism is
        # that it is genuinely non-neutral even at ordinary/basal glucose
        # once a drug is active -- a real physiological fasting-glucose-
        # lowering effect, not merely postprandial kinetics.
        p = RenalClearanceParameters.for_dapagliflozin()
        delta = p.target_delta_mM_per_min(GLUCOSE_REFERENCE_MM)
        assert delta < 0.0

    def test_drug_delta_more_negative_at_higher_glucose(self):
        # More filtered glucose load -> more excess excretion under the
        # drug -> a larger (more negative) clearance delta.
        p = RenalClearanceParameters.for_dapagliflozin()
        delta_low  = p.target_delta_mM_per_min(5.0)
        delta_high = p.target_delta_mM_per_min(15.0)
        assert delta_high < delta_low

    def test_extreme_glucose_values_do_not_crash_or_produce_nan(self):
        # Regression: hill_inhibition(F, tm, p) ** (1/p) with plain Python
        # floats raised an uncaught OverflowError for a sufficiently
        # extreme filtered load, and glucose_mM=inf silently produced a
        # NaN delta (inf * 0.0) with no error at all -- found during
        # self-verification when the mandatory review-agent pass hit a
        # session-quota outage (see CLAUDE.md). Neither a divergent
        # upstream ODE state nor a malformed/adversarial input should be
        # able to crash this pathway or silently inject a NaN into the
        # solver.
        import math
        p = RenalClearanceParameters.for_dapagliflozin()
        for g in (1e5, 1e10, 1e100, float("inf"), float("nan")):
            delta = p.target_delta_mM_per_min(g)
            assert math.isfinite(delta), f"G={g!r} produced non-finite delta {delta!r}"

    def test_dapagliflozin_implied_daily_uge_matches_literature_range(self):
        # Sanity-check the whole unit-conversion chain (GFR, Tm, beta,
        # GLUCOSE_MOLECULAR_WEIGHT, Vg_L) against PMID:23735727's own
        # quoted class-wide range: "24-h urine glucose excretion (60-80
        # g/day)". Using the same "glucose held at fasting all day"
        # simplification the module's own beta derivation uses (see
        # SGLT2I_BETA_CLASS_DEFAULT docstring) -- this is the value beta
        # was picked to reproduce, so this test mostly guards against a
        # unit-conversion regression (e.g. a dropped /1000 or wrong MW),
        # not an independent literature check.
        from biosim_ode.meal import GLUCOSE_MOLECULAR_WEIGHT

        p = RenalClearanceParameters.for_dapagliflozin()
        delta = p.target_delta_mM_per_min(GLUCOSE_REFERENCE_MM)
        mg_per_min = -delta * p.vg_L * GLUCOSE_MOLECULAR_WEIGHT
        g_per_day = mg_per_min * 1440.0 / 1000.0
        assert 55.0 < g_per_day < 85.0

    def test_dapagliflozin_uses_beta_class_default(self):
        p = RenalClearanceParameters.for_dapagliflozin()
        assert p.beta == pytest.approx(SGLT2I_BETA_CLASS_DEFAULT)
        assert p.tm_mg_per_min == pytest.approx(317.0)

    def test_empagliflozin_defaults_to_day14_beta(self):
        p = RenalClearanceParameters.for_empagliflozin()
        assert p.beta == pytest.approx(SGLT2I_BETA_CLASS_DEFAULT)
        assert p.tm_mg_per_min == pytest.approx(337.0)

    def test_empagliflozin_beta_increases_with_treatment_duration(self):
        p_early = RenalClearanceParameters.for_empagliflozin(days_of_treatment=2.0)
        p_late  = RenalClearanceParameters.for_empagliflozin(days_of_treatment=14.0)
        assert p_early.beta < p_late.beta
        assert p_late.beta == pytest.approx(SGLT2I_BETA_CLASS_DEFAULT)

    def test_empagliflozin_duration_clamped_to_measured_range(self):
        # No extrapolation below the earliest measured point (48h) or
        # above the latest (day 14) -- both fall back to the nearest
        # measured endpoint's beta.
        below = RenalClearanceParameters.for_empagliflozin(days_of_treatment=0.5)
        at_2  = RenalClearanceParameters.for_empagliflozin(days_of_treatment=2.0)
        above = RenalClearanceParameters.for_empagliflozin(days_of_treatment=30.0)
        at_14 = RenalClearanceParameters.for_empagliflozin(days_of_treatment=14.0)
        assert below.beta == pytest.approx(at_2.beta)
        assert above.beta == pytest.approx(at_14.beta)

    def test_dapagliflozin_t2d_uses_elevated_tmg_not_the_healthy_value(self):
        # DeFronzo et al. 2013 (PMID:23735727) measured a genuinely
        # elevated baseline TmG in its type-2-diabetic cohort (420 vs.
        # 317 mg/min, healthy) -- t2d=True must select that value, not
        # silently default to the healthy figure regardless of disease
        # state (the prior behaviour of this classmethod, MASTER_PLAN 1.8).
        healthy = RenalClearanceParameters.for_dapagliflozin(t2d=False)
        t2d = RenalClearanceParameters.for_dapagliflozin(t2d=True)
        assert healthy.tm_mg_per_min == pytest.approx(317.0)
        assert t2d.tm_mg_per_min == pytest.approx(420.0)
        # beta is unchanged by t2d -- the paper gives no population-split
        # 24h UGE figure to derive a different value from.
        assert t2d.beta == pytest.approx(healthy.beta)

    def test_dapagliflozin_default_is_healthy_not_t2d(self):
        # t2d defaults to False -- byte-identical to the pre-existing,
        # already-tested no-argument call.
        default = RenalClearanceParameters.for_dapagliflozin()
        explicit_healthy = RenalClearanceParameters.for_dapagliflozin(t2d=False)
        assert default.tm_mg_per_min == pytest.approx(explicit_healthy.tm_mg_per_min)

    def test_dapagliflozin_t2d_elevated_tmg_yields_less_no_drug_excretion(self):
        # A higher baseline reabsorptive capacity means the undrugged
        # kidney reabsorbs more before saturating -- at a given moderate-
        # to-high glucose level, the T2D (no-drug) excretion rate should
        # be lower than the healthy (no-drug) excretion rate for the same
        # filtered load, confirming the elevated Tm actually changes the
        # pathway's behaviour, not just a stored-but-unused field.
        from biosim_ode.renal_pathway import _filtered_load_mg_per_min, _excretion_mg_per_min
        healthy = RenalClearanceParameters.for_dapagliflozin(t2d=False)
        t2d = RenalClearanceParameters.for_dapagliflozin(t2d=True)
        glucose_mM = 10.0
        filtered = _filtered_load_mg_per_min(glucose_mM, healthy.gfr_ml_per_min)
        e_healthy = _excretion_mg_per_min(filtered, healthy.tm_mg_per_min, healthy.splay_p, 0.0)
        e_t2d = _excretion_mg_per_min(filtered, t2d.tm_mg_per_min, t2d.splay_p, 0.0)
        assert e_t2d < e_healthy

    def test_empagliflozin_t2d_uses_elevated_tmg_and_its_own_reduction_fractions(self):
        # Al-Jobori et al. 2017 (PMID:28428225) reports BOTH a genuinely
        # elevated baseline TmG for T2D (459 vs. 337 mg/min) AND distinct
        # day2/day14 TmG-reduction fractions for the two groups (44%/65%
        # T2D vs. 53%/75% healthy) -- t2d=True must use the T2D cohort's
        # own figures throughout, not mix a T2D TmG with a healthy-cohort
        # beta-scaling fraction.
        healthy_day2 = RenalClearanceParameters.for_empagliflozin(days_of_treatment=2.0, t2d=False)
        t2d_day2 = RenalClearanceParameters.for_empagliflozin(days_of_treatment=2.0, t2d=True)
        assert healthy_day2.tm_mg_per_min == pytest.approx(337.0)
        assert t2d_day2.tm_mg_per_min == pytest.approx(459.0)
        # At day 2 the two groups' reduction fractions (44%/65% vs.
        # 53%/75%) genuinely differ -- 0.44/0.65 != 0.53/0.75 -- so the
        # resulting betas must differ too, not just the TmG.
        assert t2d_day2.beta != pytest.approx(healthy_day2.beta)

    def test_empagliflozin_t2d_default_duration_is_healthy_not_t2d(self):
        default = RenalClearanceParameters.for_empagliflozin()
        explicit_healthy = RenalClearanceParameters.for_empagliflozin(t2d=False)
        assert default.tm_mg_per_min == pytest.approx(explicit_healthy.tm_mg_per_min)
        assert default.beta == pytest.approx(explicit_healthy.beta)

    def test_empagliflozin_rejects_non_finite_duration(self):
        # Regression: max()/min()-based clamping silently resolved
        # days_of_treatment=nan to the day-14 beta with no error (NaN
        # comparisons are always False, so min(14.0, nan) silently
        # returns 14.0) -- found during self-verification when the
        # mandatory review-agent pass hit a session-quota outage. A
        # non-finite duration must raise, not silently misrepresent an
        # invalid input as a confident, specific value.
        with pytest.raises(ValueError):
            RenalClearanceParameters.for_empagliflozin(days_of_treatment=float("nan"))
        with pytest.raises(ValueError):
            RenalClearanceParameters.for_empagliflozin(days_of_treatment=float("inf"))

    def test_dapagliflozin_beta_class_default_override(self):
        # MASTER_PLAN 2.2: lets a Monte Carlo UQ ensemble sample beta
        # rather than always using the fixed SGLT2I_BETA_CLASS_DEFAULT.
        default = RenalClearanceParameters.for_dapagliflozin()
        overridden = RenalClearanceParameters.for_dapagliflozin(beta_class_default=0.30)
        assert default.beta == pytest.approx(SGLT2I_BETA_CLASS_DEFAULT)
        assert overridden.beta == pytest.approx(0.30)

    def test_dapagliflozin_beta_class_default_omitted_matches_prior_behaviour_exactly(self):
        omitted = RenalClearanceParameters.for_dapagliflozin(t2d=True)
        explicit = RenalClearanceParameters.for_dapagliflozin(
            t2d=True, beta_class_default=SGLT2I_BETA_CLASS_DEFAULT,
        )
        assert omitted.beta == pytest.approx(explicit.beta)
        assert omitted.tm_mg_per_min == pytest.approx(explicit.tm_mg_per_min)

    def test_empagliflozin_beta_class_default_override_scales_with_fraction_of_day14(self):
        # beta_class_default must be the pre-scaling class-wide default,
        # not the already-duration-scaled final beta -- overriding it
        # should still be scaled by the same day2/day14 onset fraction.
        default_day2 = RenalClearanceParameters.for_empagliflozin(days_of_treatment=2.0)
        overridden_day2 = RenalClearanceParameters.for_empagliflozin(
            days_of_treatment=2.0, beta_class_default=0.20,
        )
        ratio = default_day2.beta / SGLT2I_BETA_CLASS_DEFAULT
        assert overridden_day2.beta == pytest.approx(0.20 * ratio)

    def test_empagliflozin_beta_class_default_omitted_matches_prior_behaviour_exactly(self):
        omitted = RenalClearanceParameters.for_empagliflozin(days_of_treatment=8.0, t2d=True)
        explicit = RenalClearanceParameters.for_empagliflozin(
            days_of_treatment=8.0, t2d=True, beta_class_default=SGLT2I_BETA_CLASS_DEFAULT,
        )
        assert omitted.beta == pytest.approx(explicit.beta)


class TestRenalPathway:
    def test_instantiation(self):
        p = RenalPathway()
        assert p.name == "renal_clearance_v1"

    def test_species_names(self):
        assert RenalPathway().species_names == ["RGC"]

    def test_y0_length(self):
        assert len(RenalPathway().y0) == 1

    def test_dydt_returns_one_value(self):
        p = RenalPathway()
        derivs = p.dydt(0.0, p.y0)
        assert len(derivs) == 1

    def test_implements_base_pathway(self):
        from biosim_ode.base_pathway import BasePathway
        assert isinstance(RenalPathway(), BasePathway)

    def test_output_signals_key(self):
        p = RenalPathway()
        signals = p.output_signals(0.0, p.y0)
        assert "renal_clearance_delta_mM_per_min" in signals

    def test_no_drug_y0_is_genuine_fixed_point(self):
        p = RenalPathway(RenalClearanceParameters(beta=0.0))
        assert p.validate_y0(tolerance=1e-8) == []

    def test_drug_y0_is_genuine_fixed_point(self):
        # validate_y0() calls dydt(0, y0, signals=None) -- with no injected
        # glucose_mM signal, dydt must default to GLUCOSE_REFERENCE_MM,
        # exactly matching what y0 was computed against at construction.
        p = RenalPathway(RenalClearanceParameters.for_dapagliflozin())
        assert p.validate_y0(tolerance=1e-8) == []

    def test_no_drug_signal_stays_zero_standalone_across_glucose_values(self):
        # The pathway's own relaxing state, not just target_delta_mM_per_min
        # directly, must also report zero -- confirms the state genuinely
        # tracks its (always-zero) target rather than merely starting there.
        p = RenalPathway(RenalClearanceParameters(beta=0.0))
        for g in [1.0, 5.0, 12.0, 20.0]:
            derivs = p.dydt(0.0, p.y0, signals={"glucose_mM": g})
            assert derivs[0] == pytest.approx(0.0, abs=1e-12)

    def test_missing_glucose_signal_defaults_to_reference(self):
        p = RenalPathway(RenalClearanceParameters.for_dapagliflozin())
        with_signal    = p.dydt(0.0, p.y0, signals={"glucose_mM": GLUCOSE_REFERENCE_MM})
        without_signal = p.dydt(0.0, p.y0, signals=None)
        assert with_signal[0] == pytest.approx(without_signal[0])

    def test_higher_glucose_drives_state_toward_more_clearance(self):
        p = RenalPathway(RenalClearanceParameters.for_dapagliflozin())
        d_low  = p.dydt(0.0, p.y0, signals={"glucose_mM": 5.0})[0]
        d_high = p.dydt(0.0, p.y0, signals={"glucose_mM": 15.0})[0]
        # Starting from the fasting fixed point, a higher glucose signal
        # should push the clearance-delta state further negative (more
        # excretion), i.e. a more negative derivative.
        assert d_high < d_low


class TestRenalCoupledWithGI:
    def test_couple_with_gi_preserves_basal_steady_state_no_drug(self):
        from biosim_ode import GlucoseInsulinPathway
        from biosim_ode.renal_pathway import couple_with_gi
        composite = couple_with_gi(GlucoseInsulinPathway(), RenalPathway(RenalClearanceParameters(beta=0.0)))
        derivs = composite.dydt(0.0, composite.y0)
        assert derivs[0] == pytest.approx(0.0, abs=1e-6)   # dG
        assert derivs[1] == pytest.approx(0.0, abs=1e-6)   # dI

    def test_couple_with_gi_drug_active_perturbs_basal_glucose_downward(self):
        # Deliberately the OPPOSITE assertion from hepatic's basal-
        # preservation tests: this is the architectural point of the
        # whole design (see renal_pathway.py's module docstring) -- an
        # active SGLT2i drug genuinely lowers dG even at GI's own
        # untreated basal, because the mechanism is continuously active,
        # not insulin-gated. A test that asserted dG==0 here would be
        # asserting the drug has no effect at all.
        from biosim_ode import GlucoseInsulinPathway
        from biosim_ode.renal_pathway import couple_with_gi
        composite = couple_with_gi(GlucoseInsulinPathway(), RenalPathway(RenalClearanceParameters.for_dapagliflozin()))
        derivs = composite.dydt(0.0, composite.y0)
        assert derivs[0] < -1e-6   # dG, genuinely negative

    def test_fasted_chronic_dapagliflozin_settles_to_lower_glucose(self):
        # The emergent-fixed-point prediction: with no meals, glucose
        # should relax DOWN from its default basal to a new, lower
        # steady state under chronic dapagliflozin -- not drift
        # unboundedly, and not stay at the untreated basal.
        from scipy.integrate import solve_ivp
        from biosim_ode import GlucoseInsulinPathway
        from biosim_ode.renal_pathway import couple_with_gi

        composite = couple_with_gi(GlucoseInsulinPathway(), RenalPathway(RenalClearanceParameters.for_dapagliflozin()))

        def dydt(t, y):
            return composite.dydt(t, y)

        sol = solve_ivp(dydt, [0, 1440], composite.y0, method="Radau", rtol=1e-8, atol=1e-10)
        assert sol.success
        final_G = sol.y[0, -1]
        assert final_G < GlucoseInsulinPathway().y0[0]   # below untreated basal (5.0 mM)
        assert final_G > 0.0   # physically sane, not a runaway

    def test_dapagliflozin_coupling_lowers_glucose_peak_after_meal(self):
        from biosim_ode import GlucoseInsulinPathway, ODESystem, high_sugar_meal
        from biosim_ode.renal_pathway import couple_with_gi

        baseline = ODESystem(pathway=GlucoseInsulinPathway()) \
            .simulate_meal(high_sugar_meal(), duration_s=14400)
        composite = couple_with_gi(GlucoseInsulinPathway(), RenalPathway(RenalClearanceParameters.for_dapagliflozin()))
        treated = ODESystem(pathway=composite).simulate_meal(high_sugar_meal(), duration_s=14400)

        peak_base, _ = baseline.glucose_peak()
        peak_treat, _ = treated.glucose_peak()
        assert peak_treat < peak_base

    def test_no_drug_composition_does_not_perturb_meal_response_at_all(self):
        # Confirms zero double-counting against the existing (uncited)
        # baseline clearance (p1): composing this pathway with beta=0
        # must reproduce the same postprandial trajectory as not composing
        # it at all (up to adaptive-solver step-size noise from the extra
        # composed state dimension -- the signal itself is verified
        # analytically exact-zero in TestRenalClearanceParameters/
        # TestRenalPathway above; a real double-count would be orders of
        # magnitude larger than solver noise, not comparable to it).
        import numpy as np
        from biosim_ode import GlucoseInsulinPathway, ODESystem, high_sugar_meal
        from biosim_ode.renal_pathway import couple_with_gi

        baseline = ODESystem(pathway=GlucoseInsulinPathway()) \
            .simulate_meal(high_sugar_meal(), duration_s=14400)
        composite = couple_with_gi(GlucoseInsulinPathway(), RenalPathway(RenalClearanceParameters(beta=0.0)))
        no_drug = ODESystem(pathway=composite).simulate_meal(high_sugar_meal(), duration_s=14400)

        np.testing.assert_allclose(no_drug.species["G"], baseline.species["G"], rtol=1e-4)
