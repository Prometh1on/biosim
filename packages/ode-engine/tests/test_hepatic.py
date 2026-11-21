# SAVE TO: packages\ode-engine\tests\test_hepatic.py
"""
Test suite for the hepatic glucose production (EGP) pathway extension.
"""
import pytest
import numpy as np
from biosim_ode.hepatic import HepGlucPathway, EGPParameters, EGP_BASAL_MM_MIN


class TestEGPParameters:
    def test_default_construction(self):
        p = EGPParameters()
        assert p.egp_basal == pytest.approx(EGP_BASAL_MM_MIN)
        assert p.metformin_scale == pytest.approx(1.0)
        assert p.t2d_egp_scale == pytest.approx(1.0)

    def test_t2d_params_elevated_basal(self):
        p = EGPParameters.for_type2_diabetes()
        assert p.egp_basal > EGP_BASAL_MM_MIN

    def test_t2d_params_higher_ic50(self):
        default = EGPParameters()
        t2d = EGPParameters.for_type2_diabetes()
        assert t2d.ic50_pM > default.ic50_pM

    def test_t2d_no_maximum_velocity_defect(self):
        # Beck-Nielsen et al. 2002, PMID:12643178: T2D shows "normal maximal
        # suppression (no maximum velocity defect)" -- egp_max and
        # t2d_egp_scale must be unchanged from healthy default, since the
        # impairment is purely a dose-response (IC50) shift plus the
        # separately-measured +12% basal elevation, not an independent
        # scale-up mechanism.
        default = EGPParameters()
        t2d = EGPParameters.for_type2_diabetes()
        assert t2d.egp_max == pytest.approx(default.egp_max)
        assert t2d.t2d_egp_scale == pytest.approx(1.0)

    def test_t2d_basal_matches_literature_12_percent_elevation(self):
        # PMID:12643178: fasting HGP "may be increased by a mean of 12%"
        # -- not the previous, unanchored 30%.
        t2d = EGPParameters.for_type2_diabetes()
        assert t2d.egp_basal == pytest.approx(EGP_BASAL_MM_MIN * 1.12)

    def test_t2d_y0_is_genuine_fixed_point(self):
        # Regression: the T2D preset previously stacked THREE independent
        # "impairment" multipliers (egp_basal*1.30, ic50*2.5, scale*1.40)
        # that were jointly over-determined -- no combination could make
        # the elevated egp_basal an actual fixed point of its own ODE (see
        # CLAUDE.md KNOWN BUGS FIXED #14). Fixed by deriving ic50_pM from
        # egp_basal via the same analytic technique as the default preset.
        p = HepGlucPathway(EGPParameters.for_type2_diabetes())
        assert p.validate_y0(tolerance=1e-8) == []

    def test_metformin_reduces_basal_egp(self):
        # Regression: metformin used to be represented as a
        # `metformin_scale` multiplier on the suppression TARGET only,
        # leaving egp_basal at its unmodified default -- see
        # test_metformin_y0_is_genuine_fixed_point below for why that was
        # a real double-counting bug, not the "intentional onset dynamic"
        # an earlier version of this test claimed. Metformin's effect is
        # now represented directly via a lowered egp_basal, the same way
        # T2D's effect is represented via an elevated one.
        p = EGPParameters.for_metformin(1000.0)
        assert p.egp_basal < EGP_BASAL_MM_MIN
        assert p.egp_basal > 0.0
        assert p.metformin_scale == pytest.approx(1.0)   # retired; effect is in egp_basal now

    def test_higher_metformin_more_reduction(self):
        p_low  = EGPParameters.for_metformin(500.0)
        p_high = EGPParameters.for_metformin(2000.0)
        assert p_high.egp_basal < p_low.egp_basal

    def test_hyperthyroidism_params_elevated_basal(self):
        p = EGPParameters.for_hyperthyroidism()
        assert p.egp_basal > EGP_BASAL_MM_MIN

    def test_hyperthyroidism_basal_matches_literature_20_percent_elevation(self):
        # Wennlund et al. 1986, PMID:3510000: basal HGP "increased by 20%"
        # in hyperthyroid vs. normal subjects (P<0.02).
        p = EGPParameters.for_hyperthyroidism()
        assert p.egp_basal == pytest.approx(EGP_BASAL_MM_MIN * 1.20)

    def test_hyperthyroidism_higher_ic50(self):
        # The same rightward dose-response shift Wennlund measured
        # directly (less suppression despite a larger insulin rise).
        default = EGPParameters()
        hyper = EGPParameters.for_hyperthyroidism()
        assert hyper.ic50_pM > default.ic50_pM

    def test_hyperthyroidism_no_maximum_velocity_defect(self):
        # Wennlund reports a dose-response shift, not a suppression-
        # capacity defect -- same reasoning as T2D's analogous test.
        default = EGPParameters()
        hyper = EGPParameters.for_hyperthyroidism()
        assert hyper.egp_max == pytest.approx(default.egp_max)
        assert hyper.t2d_egp_scale == pytest.approx(1.0)

    def test_hyperthyroidism_y0_is_genuine_fixed_point(self):
        p = HepGlucPathway(EGPParameters.for_hyperthyroidism())
        assert p.validate_y0(tolerance=1e-8) == []

    def test_metformin_y0_is_genuine_fixed_point(self):
        # Regression: for_metformin()'s y0 used to NOT be a fixed point of
        # its own ODE -- a previous version of this test asserted that was
        # "intentional" (modeling "the onset of drug therapy"), but nothing
        # in the codebase ever actually distinguished onset-of-therapy from
        # chronic/established use (DrugDose.is_chronic is documented but
        # never read anywhere) -- the drift fired unconditionally, even for
        # a DrugDose explicitly marked is_chronic=True and dosed a week
        # ago. Composed with GlucoseInsulinPathway (whose own basal state
        # already reflects the SEPARATE, documented static Gb modifier
        # from physiology/drug_interactions.py's metformin_modifiers()
        # (-14.7% as of CLAUDE.md's mechanism-weighting correction, -4%
        # when this test was written -- the fixed-point property here
        # doesn't depend on the modifier's exact magnitude), this meant
        # even a fasted
        # (zero-meal) chronic-metformin scenario drifted glucose downward
        # over time purely as an artifact -- confirmed directly: -0.10 mM
        # by 4h, on top of the already-intended static effect. Fixed via
        # the same anchor-then-derive technique as the T2D preset
        # (test_t2d_y0_is_genuine_fixed_point above).
        p = HepGlucPathway(EGPParameters.for_metformin(1000.0))
        assert p.validate_y0(tolerance=1e-8) == []


class TestHepGlucPathway:
    def test_instantiation(self):
        p = HepGlucPathway()
        assert p.name == "hep_gluc_v1"

    def test_species_names(self):
        assert "EGP" in HepGlucPathway().species_names

    def test_y0_length(self):
        assert len(HepGlucPathway().y0) == 1

    def test_y0_positive(self):
        assert HepGlucPathway().y0[0] > 0

    def test_dydt_returns_one_value(self):
        p = HepGlucPathway()
        derivs = p.dydt(0.0, p.y0)
        assert len(derivs) == 1

    def test_high_insulin_suppresses_egp(self):
        p = HepGlucPathway()
        high_i = {"insulin_pM": 500.0}
        low_i  = {"insulin_pM": 10.0}
        dEGP_high = p.dydt(0.0, [p.params.egp_max], signals=high_i)[0]
        dEGP_low  = p.dydt(0.0, [p.params.egp_max], signals=low_i)[0]
        # High insulin should push EGP down faster
        assert dEGP_high < dEGP_low

    def test_high_insulin_negative_derivative(self):
        p = HepGlucPathway()
        derivs = p.dydt(0.0, [p.params.egp_max], signals={"insulin_pM": 300.0})
        assert derivs[0] < 0

    def test_output_signals_keys(self):
        p = HepGlucPathway()
        signals = p.output_signals(0.0, p.y0)
        assert "egp_mM_per_min" in signals
        assert "egp_suppression_fraction" in signals

    def test_implements_base_pathway(self):
        from biosim_ode.base_pathway import BasePathway
        assert isinstance(HepGlucPathway(), BasePathway)

    def test_basal_egp_accessor(self):
        p = HepGlucPathway()
        assert p.basal_egp() == pytest.approx(EGP_BASAL_MM_MIN)

    def test_composite_with_glucose_insulin(self):
        from biosim_ode import GlucoseInsulinPathway
        from biosim_ode.coupling import CompositePathway
        composite = CompositePathway([GlucoseInsulinPathway(), HepGlucPathway()])
        assert len(composite.y0) == 7

    def test_t2d_higher_egp_baseline(self):
        healthy = HepGlucPathway()
        t2d     = HepGlucPathway(EGPParameters.for_type2_diabetes())
        assert t2d.y0[0] > healthy.y0[0]

    def test_y0_is_genuine_fixed_point(self):
        # Regression for the EGP_SUPPRESSION_IC50 calibration bug: y0 must be
        # an actual fixed point of dydt at basal insulin, not just a
        # plausible-looking literal. This was previously false -- EGP
        # silently drifted from 0.032 toward ~0.0446 over ~600 minutes even
        # with the module's own default parameters and no disease/drug
        # modifiers applied.
        p = HepGlucPathway()
        assert p.validate_y0(tolerance=1e-8) == []

    def test_egp_stable_over_long_trajectory_at_basal_insulin(self):
        # The practical consequence of the bug above: standalone integration
        # at basal insulin (66 pM, matching PathwayParameters.Ib) must not
        # let EGP drift over time, not just at t=0.
        from scipy.integrate import solve_ivp
        p = HepGlucPathway()

        def dydt_at_basal_insulin(t, y):
            return p.dydt(t, y, signals={"insulin_pM": 66.0})

        sol = solve_ivp(dydt_at_basal_insulin, [0, 600], p.y0, method="Radau", rtol=1e-8, atol=1e-10)
        assert sol.success
        assert sol.y[0, -1] == pytest.approx(EGP_BASAL_MM_MIN, abs=1e-6)

    def test_simulate_egp_falls_with_insulin(self):
        from scipy.integrate import solve_ivp
        p = HepGlucPathway()
        # Start at max EGP, inject high insulin signal
        def dydt_with_insulin(t, y):
            return p.dydt(t, y, signals={"insulin_pM": 300.0})
        sol = solve_ivp(dydt_with_insulin, [0, 120], [p.params.egp_max], method="Radau")
        assert sol.success
        # EGP should fall from max toward suppressed level
        assert sol.y[0, -1] < sol.y[0, 0]


class TestHepaticCoupledWithGI:
    def test_couple_with_gi_preserves_basal_steady_state_default(self):
        from biosim_ode import GlucoseInsulinPathway
        from biosim_ode.hepatic import couple_with_gi
        composite = couple_with_gi(GlucoseInsulinPathway(), HepGlucPathway())
        derivs = composite.dydt(0.0, composite.y0)
        assert derivs[0] == pytest.approx(0.0, abs=1e-6)   # dG
        assert derivs[1] == pytest.approx(0.0, abs=1e-6)   # dI

    def test_couple_with_gi_preserves_basal_steady_state_t2d(self):
        # Regression: the coupling must be zero at THIS instance's own basal
        # even for a preset with an elevated egp_basal (T2D), not just the
        # module default -- otherwise composing T2D would spuriously perturb
        # the basal state on top of the separate disease-state Gb modifier.
        from biosim_ode import GlucoseInsulinPathway
        from biosim_ode.hepatic import couple_with_gi
        composite = couple_with_gi(GlucoseInsulinPathway(), HepGlucPathway(EGPParameters.for_type2_diabetes()))
        derivs = composite.dydt(0.0, composite.y0)
        assert derivs[0] == pytest.approx(0.0, abs=1e-6)
        assert derivs[1] == pytest.approx(0.0, abs=1e-6)

    def test_couple_with_gi_preserves_basal_steady_state_hyperthyroidism(self):
        # Same fixed-point property as T2D's analogous test above -- this
        # is the SOLE representation of hyperthyroidism's glucose effect
        # (no separate static Gb modifier exists, see
        # physiology/modifiers.py's DiseaseState.HYPERTHYROIDISM), so it's
        # even more important here that the coupling contributes nothing
        # at basal: a nonzero derivative here would mean hyperthyroidism
        # spuriously perturbs FASTING glucose, which the literature this
        # preset is anchored on (Wennlund et al. 1986, PMID:3510000)
        # never measured -- only postprandial/dynamic kinetics.
        from biosim_ode import GlucoseInsulinPathway
        from biosim_ode.hepatic import couple_with_gi
        composite = couple_with_gi(GlucoseInsulinPathway(), HepGlucPathway(EGPParameters.for_hyperthyroidism()))
        derivs = composite.dydt(0.0, composite.y0)
        assert derivs[0] == pytest.approx(0.0, abs=1e-6)
        assert derivs[1] == pytest.approx(0.0, abs=1e-6)

    def test_hyperthyroidism_hepatic_coupling_postmeal_close_to_normal(self):
        # Same comparison basis as
        # test_t2d_hepatic_coupling_postmeal_suppression_close_to_normal
        # below -- default-hepatic-COMPOSED vs. hyperthyroidism-hepatic-
        # COMPOSED (not bare GI with no hepatic pathway at all, which
        # measures "does having any hepatic coupling change dynamics",
        # a different and much noisier question). Investigated directly
        # before writing this assertion: naively expected the impaired
        # (rightward-shifted) suppression to RAISE the postprandial
        # peak, matching metformin's mirror-image lowering -- it does
        # not. Composed vs. composed, the difference is small (-0.4%
        # AUC, -0.2% peak, hyperthyroidism if anything very slightly
        # LOWER) -- the same "no maximum velocity defect, modest
        # rightward IC50 shift" shape as the T2D preset, which produces
        # a "close to normal" postprandial response by the same
        # anchor-then-derive construction, not a dramatic one. This is
        # NOT claiming hyperthyroidism doesn't affect real postprandial
        # glucose -- Wennlund 1986 measures a real, statistically
        # significant hepatic resistance to insulin -- it's that its
        # magnitude, run through this platform's specific pathway
        # kinetics (relaxation rate, Hill shape), doesn't dominate over
        # the meal-driven insulin excursion at this scale. Worth
        # revisiting if a real postprandial-glucose clinical benchmark
        # for hyperthyroidism is ever found (none was this session).
        from biosim_ode import GlucoseInsulinPathway, ODESystem, high_sugar_meal
        from biosim_ode.hepatic import couple_with_gi
        import numpy as np

        healthy = ODESystem(pathway=couple_with_gi(GlucoseInsulinPathway(), HepGlucPathway())) \
            .simulate_meal(high_sugar_meal(), duration_s=14400)
        hyper = ODESystem(pathway=couple_with_gi(GlucoseInsulinPathway(), HepGlucPathway(EGPParameters.for_hyperthyroidism()))) \
            .simulate_meal(high_sugar_meal(), duration_s=14400)

        auc_healthy = np.trapezoid(healthy.species["G"], healthy.times_s)
        auc_hyper = np.trapezoid(hyper.species["G"], hyper.times_s)
        assert auc_hyper == pytest.approx(auc_healthy, rel=0.01)

    def test_metformin_hepatic_coupling_lowers_glucose_peak(self):
        from biosim_ode import GlucoseInsulinPathway, ODESystem, high_sugar_meal
        from biosim_ode.hepatic import couple_with_gi

        baseline = ODESystem(pathway=GlucoseInsulinPathway()) \
            .simulate_meal(high_sugar_meal(), duration_s=14400)
        composite = couple_with_gi(GlucoseInsulinPathway(), HepGlucPathway(EGPParameters.for_metformin(1000.0)))
        treated = ODESystem(pathway=composite).simulate_meal(high_sugar_meal(), duration_s=14400)

        peak_base, _ = baseline.glucose_peak()
        peak_treat, _ = treated.glucose_peak()
        assert peak_treat < peak_base

    def test_t2d_hepatic_coupling_postmeal_suppression_close_to_normal(self):
        # Regression for KNOWN BUGS FIXED #14 (T2D EGP preset corrected to
        # match Beck-Nielsen et al. 2002, PMID:12643178): postprandial HGP
        # suppression in T2D is "close to normal" per that literature --
        # the T2D preset's impairment is a MODEST rightward IC50 shift plus
        # a +12% elevated basal (both captured in the delta-formulation's
        # own y0), NOT an independent suppression-magnitude defect. So the
        # coupling's own AUC contribution (isolated from the static Gb
        # modifier, both runs use the same Gb) should be close to healthy,
        # not meaningfully higher -- a big gap here would mean the T2D
        # preset drifted back toward the old, literature-contradicted
        # "no maximum velocity defect" violation (previously *2.5 IC50 and
        # *1.40 scale, both now known to overstate the real impairment).
        from biosim_ode import GlucoseInsulinPathway, ODESystem, high_sugar_meal
        from biosim_ode.hepatic import couple_with_gi
        import numpy as np

        healthy = ODESystem(pathway=couple_with_gi(GlucoseInsulinPathway(), HepGlucPathway())) \
            .simulate_meal(high_sugar_meal(), duration_s=14400)
        t2d = ODESystem(pathway=couple_with_gi(GlucoseInsulinPathway(), HepGlucPathway(EGPParameters.for_type2_diabetes()))) \
            .simulate_meal(high_sugar_meal(), duration_s=14400)

        auc_healthy = np.trapezoid(healthy.species["G"], healthy.times_s)
        auc_t2d = np.trapezoid(t2d.species["G"], t2d.times_s)
        assert auc_t2d == pytest.approx(auc_healthy, rel=0.01)
