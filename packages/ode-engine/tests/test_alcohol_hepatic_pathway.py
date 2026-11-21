# SAVE TO: packages\ode-engine\tests\test_alcohol_hepatic_pathway.py
"""
Test suite for the alcohol hepatic glucose production (EGP) suppression
pathway.

Fixes a real, confirmed bug found by a fresh Track 1 scoping fan-out:
the previous alcohol_modifiers() (biosim_physiology.drug_interactions)
reduced Gb (fasting glucose CONCENTRATION), which both available primary
sources -- Siler et al. 1998 (PMID:9815011) and Shelmet et al. 1988
(PMID:3280601) -- explicitly state does NOT happen. The real, cited
effect is a transient reduction in hepatic glucose PRODUCTION (a flux),
modelled here instead as an additive, time-varying hepatic_egp_delta
contribution -- see AlcoholHepaticPathway's own module docstring for
the full citation/derivation detail.
"""
import math

import pytest

from biosim_ode import GlucoseInsulinPathway, ODESystem
from biosim_ode.hepatic import EGP_BASAL_MM_MIN
from biosim_ode.alcohol_hepatic_pathway import (
    AlcoholHepaticParameters,
    AlcoholHepaticPathway,
    couple_with_gi,
    alcohol_hepatic_couplings,
    ALCOHOL_EGP_REDUCTION_FRACTION_AT_TESTED_DOSE,
    ALCOHOL_TESTED_DOSE_G,
    ALCOHOL_PEAK_TIME_MIN,
    ALCOHOL_RESOLVES_BY_MIN,
)


class TestAlcoholHepaticParameters:
    def test_zero_dose_is_neutral_at_every_time(self):
        p = AlcoholHepaticParameters(dose_g=0.0)
        for t in [0.0, 10.0, 120.0, 300.0, 500.0]:
            assert p.delta_mM_per_min(t) == pytest.approx(0.0)

    def test_negative_dose_raises(self):
        with pytest.raises(ValueError):
            AlcoholHepaticParameters(dose_g=-1.0)

    def test_nan_dose_raises(self):
        with pytest.raises(ValueError):
            AlcoholHepaticParameters(dose_g=float("nan"))

    def test_non_finite_dose_time_min_raises(self):
        with pytest.raises(ValueError):
            AlcoholHepaticParameters(dose_g=28.0, dose_time_min=float("inf"))

    def test_negative_egp_basal_raises(self):
        with pytest.raises(ValueError):
            AlcoholHepaticParameters(dose_g=28.0, egp_basal_mM_min=-0.01)

    def test_before_dose_time_is_exactly_zero(self):
        p = AlcoholHepaticParameters(dose_g=48.0, dose_time_min=100.0)
        assert p.delta_mM_per_min(100.0) == pytest.approx(0.0)
        assert p.delta_mM_per_min(50.0) == pytest.approx(0.0)

    def test_after_resolution_window_is_exactly_zero(self):
        p = AlcoholHepaticParameters(dose_g=48.0, dose_time_min=0.0)
        assert p.delta_mM_per_min(ALCOHOL_RESOLVES_BY_MIN + 0.01) == pytest.approx(0.0)

    def test_effect_is_negative_within_the_active_window(self):
        # A real hepatic glucose PRODUCTION suppression -- the delta must
        # be negative (excess below basal), never positive.
        p = AlcoholHepaticParameters(dose_g=48.0, dose_time_min=0.0)
        for t in [1.0, 60.0, ALCOHOL_PEAK_TIME_MIN, 200.0, 359.0]:
            assert p.delta_mM_per_min(t) < 0.0

    def test_peak_occurs_at_the_documented_peak_time(self):
        p = AlcoholHepaticParameters(dose_g=48.0, dose_time_min=0.0)
        peak_value = p.delta_mM_per_min(ALCOHOL_PEAK_TIME_MIN)
        # Sampling nearby times should all be less suppressive (closer to
        # zero) than the documented peak.
        for offset in (-30.0, -10.0, 10.0, 30.0):
            nearby = p.delta_mM_per_min(ALCOHOL_PEAK_TIME_MIN + offset)
            assert nearby > peak_value

    def test_dose_scales_linearly_below_the_tested_cap(self):
        p_half = AlcoholHepaticParameters(dose_g=ALCOHOL_TESTED_DOSE_G / 2.0, dose_time_min=0.0)
        p_full = AlcoholHepaticParameters(dose_g=ALCOHOL_TESTED_DOSE_G, dose_time_min=0.0)
        assert p_half.peak_delta_mM_per_min == pytest.approx(p_full.peak_delta_mM_per_min / 2.0)

    def test_dose_caps_at_the_tested_dose_no_extrapolation(self):
        p_at_cap = AlcoholHepaticParameters(dose_g=ALCOHOL_TESTED_DOSE_G, dose_time_min=0.0)
        p_above_cap = AlcoholHepaticParameters(dose_g=ALCOHOL_TESTED_DOSE_G * 3.0, dose_time_min=0.0)
        assert p_above_cap.peak_delta_mM_per_min == pytest.approx(p_at_cap.peak_delta_mM_per_min)

    def test_peak_delta_matches_the_cited_fraction_at_the_tested_dose(self):
        # Directly reproduces Siler et al. 1998's own -12% HGP figure,
        # applied to this platform's own EGP_BASAL_MM_MIN.
        p = AlcoholHepaticParameters(dose_g=ALCOHOL_TESTED_DOSE_G, dose_time_min=0.0)
        expected = -ALCOHOL_EGP_REDUCTION_FRACTION_AT_TESTED_DOSE * EGP_BASAL_MM_MIN
        assert p.peak_delta_mM_per_min == pytest.approx(expected)


class TestAlcoholHepaticPathway:
    def test_no_doses_is_neutral(self):
        pathway = AlcoholHepaticPathway([])
        assert pathway.output_signals(100.0, [])["hepatic_egp_delta"] == pytest.approx(0.0)

    def test_multiple_doses_sum(self):
        # Two separate drinks at different times -- their contributions
        # should sum, not overwrite.
        doses = [
            AlcoholHepaticParameters(dose_g=28.0, dose_time_min=0.0),
            AlcoholHepaticParameters(dose_g=28.0, dose_time_min=60.0),
        ]
        pathway = AlcoholHepaticPathway(doses)
        combined = pathway.output_signals(120.0, [])["hepatic_egp_delta"]
        individual_sum = sum(d.delta_mM_per_min(120.0) for d in doses)
        assert combined == pytest.approx(individual_sum)
        assert combined < doses[0].delta_mM_per_min(120.0)   # strictly more negative than either alone

    def test_many_overlapping_doses_are_floored_not_unbounded(self):
        # Mandatory security-reviewer finding: a REST-legal scenario can
        # submit up to 50 drug_events (schemas.py's own
        # ScenarioIn.drug_events max_length=50). Without an aggregate
        # floor, 50 near-simultaneous saturating-dose alcohol events
        # summed to an implausible ~6x-basal-EGP suppression, driving
        # simulated glucose to severe-hypoglycemia range. Confirms the
        # aggregate is floored at -egp_basal_mM_min regardless of how
        # many doses are stacked.
        many_doses = [
            AlcoholHepaticParameters(dose_g=ALCOHOL_TESTED_DOSE_G, dose_time_min=0.0)
            for _ in range(50)
        ]
        pathway = AlcoholHepaticPathway(many_doses)
        floored = pathway.output_signals(ALCOHOL_PEAK_TIME_MIN, [])["hepatic_egp_delta"]
        naive_sum = sum(d.delta_mM_per_min(ALCOHOL_PEAK_TIME_MIN) for d in many_doses)
        assert naive_sum < -EGP_BASAL_MM_MIN   # confirms the naive sum WOULD exceed the floor
        assert floored == pytest.approx(-EGP_BASAL_MM_MIN)

    def test_stateless_zero_species(self):
        pathway = AlcoholHepaticPathway()
        assert pathway.species_names == []
        assert pathway.species_units == []
        assert pathway.y0 == []
        assert pathway.dydt(0.0, []) == []

    def test_descriptor_reports_the_correct_signal(self):
        pathway = AlcoholHepaticPathway()
        assert pathway.descriptor.output_signals == ["hepatic_egp_delta"]
        assert pathway.descriptor.input_signals == []


class TestAlcoholHepaticCouplings:
    def test_couplings_map_signal_straight_through(self):
        couplings = alcohol_hepatic_couplings("gi", "alcohol")
        assert len(couplings) == 1
        c = couplings[0]
        assert c.source_pathway == "alcohol"
        assert c.signal_name == "hepatic_egp_delta"
        assert c.target_pathway == "gi"
        assert c.target_input == "hepatic_egp_delta"

    def test_zero_dose_preserves_basal_invariance(self):
        alcohol = AlcoholHepaticPathway([AlcoholHepaticParameters(dose_g=0.0)])
        gi = GlucoseInsulinPathway()
        composite = couple_with_gi(gi, alcohol)
        assert composite.validate_y0(tolerance=1e-8) == []

    def test_active_dose_does_not_perturb_the_analytic_basal_state_at_t0(self):
        # At t=0 with dose_time_min=0, delta_mM_per_min(0) is exactly 0
        # (t_since<=0 branch) -- the composite's own y0 must still be a
        # genuine fixed point at the instant of dosing.
        alcohol = AlcoholHepaticPathway([AlcoholHepaticParameters(dose_g=48.0, dose_time_min=0.0)])
        gi = GlucoseInsulinPathway()
        composite = couple_with_gi(gi, alcohol)
        assert composite.validate_y0(tolerance=1e-8) == []

    def test_fasted_trajectory_returns_close_to_basal_after_resolution(self):
        # G1 "survives composition" criterion -- a real integration run,
        # not just the t=0 fixed-point check. A fasted (zero-meal)
        # scenario dosed with alcohol should show a real transient dip
        # (the flux perturbation) but return CLOSE to the original Gb
        # by the time the effect has resolved -- confirming no
        # persistent set-point shift, matching Siler/Shelmet's own "no
        # effect on plasma glucose concentrations/tolerance" findings.
        alcohol = AlcoholHepaticPathway([AlcoholHepaticParameters(dose_g=48.0, dose_time_min=0.0)])
        gi = GlucoseInsulinPathway()
        composite = couple_with_gi(gi, alcohol)
        system = ODESystem(composite)
        result = system.simulate(duration_s=(ALCOHOL_RESOLVES_BY_MIN + 60.0) * 60.0)
        g = result.species["G"]
        gb_uM = gi.params.Gb * 1000.0
        assert min(g) < gb_uM * 0.999   # a real, measurable transient dip occurred
        assert g[-1] == pytest.approx(gb_uM, rel=0.01)   # returns close to basal, no persistent shift
