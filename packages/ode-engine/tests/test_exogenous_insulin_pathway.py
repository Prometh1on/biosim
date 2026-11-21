# SAVE TO: packages\ode-engine\tests\test_exogenous_insulin_pathway.py
"""
Test suite for the exogenous (injected) insulin pathway.

See CLAUDE.md's exogenous-insulin milestone entry (MASTER_PLAN 1.2a) for
the full literature-verification and calibration trail -- Howey et al.
1994, PMID:8314011 (10U regular insulin, healthy subjects: Tmax=101min,
Cmax=308pM excess above baseline).
"""
import math

import pytest

from biosim_ode import GlucoseInsulinPathway, ODESystem
from biosim_ode.exogenous_insulin_pathway import (
    ExogenousInsulinParameters,
    ExogenousInsulinPathway,
    couple_with_gi,
    exogenous_insulin_couplings,
    INSULIN_LISPRO_KA_PER_MIN,
    INSULIN_LISPRO_PM_PER_UNIT,
    INSULIN_REGULAR_KA_PER_MIN,
    INSULIN_REGULAR_PM_PER_UNIT,
)


class TestExogenousInsulinParameters:
    def test_zero_dose_is_neutral_at_every_time(self):
        p = ExogenousInsulinParameters.regular_insulin(dose_units=0.0)
        for t in [0.0, 10.0, 50.0, 100.0, 500.0]:
            assert p.delta_pM_per_min(t) == pytest.approx(0.0)

    def test_negative_dose_raises(self):
        with pytest.raises(ValueError):
            ExogenousInsulinParameters(dose_units=-1.0)

    def test_nan_dose_raises(self):
        with pytest.raises(ValueError):
            ExogenousInsulinParameters(dose_units=float("nan"))

    def test_non_positive_ka_raises(self):
        with pytest.raises(ValueError):
            ExogenousInsulinParameters(dose_units=10.0, ka=0.0)
        with pytest.raises(ValueError):
            ExogenousInsulinParameters(dose_units=10.0, ka=-0.01)

    def test_dose_above_plausibility_ceiling_raises(self):
        # Guards against the realistic footgun of copy-pasting a
        # dose_mg value from a milligram-dosed drug (metformin ~2000)
        # into an insulin_regular event, where dose_units means Units,
        # not mg -- see MAX_PLAUSIBLE_SC_INSULIN_DOSE_UNITS.
        with pytest.raises(ValueError, match="plausibility ceiling"):
            ExogenousInsulinParameters(dose_units=2000.0)

    def test_dose_at_plausibility_ceiling_is_accepted(self):
        ExogenousInsulinParameters(dose_units=100.0)

    def test_overflowing_s1_0_raises(self):
        # dose_units=20.0 is well within the plausibility ceiling on its
        # own -- this targets the S1_0 = dose_units * pm_per_unit overflow
        # guard specifically, via a caller-supplied pm_per_unit (not the
        # module default) large enough that the PRODUCT overflows float64
        # even though each factor individually is finite.
        with pytest.raises(ValueError, match="S1_0"):
            ExogenousInsulinParameters(dose_units=20.0, pm_per_unit=1e307)

    def test_delta_is_zero_before_dose_time(self):
        p = ExogenousInsulinParameters.regular_insulin(dose_units=10.0, dose_time_min=50.0)
        assert p.delta_pM_per_min(49.9) == pytest.approx(0.0)
        assert p.delta_pM_per_min(50.0) == pytest.approx(0.0)
        assert p.delta_pM_per_min(50.1) > 0.0

    def test_delta_is_nonnegative_everywhere(self):
        # The Bateman-type closed form (ka^2 * S1_0 * t * exp(-ka*t)) is
        # mathematically nonnegative for t>=0 by construction -- confirm
        # this holds numerically across a wide time range, not just near
        # the peak.
        p = ExogenousInsulinParameters.regular_insulin(dose_units=10.0)
        for t in [0.001, 1, 10, 50, 100, 101, 200, 500, 2000, 1e5]:
            assert p.delta_pM_per_min(t) >= 0.0

    def test_delta_decays_to_zero_at_long_times(self):
        p = ExogenousInsulinParameters.regular_insulin(dose_units=10.0)
        assert p.delta_pM_per_min(1e6) == pytest.approx(0.0, abs=1e-9)

    def test_dose_response_is_linear(self):
        # The whole system is linear in dose_units by construction (S1_0
        # scales with dose_units, and delta_pM_per_min is linear in S1_0).
        p10 = ExogenousInsulinParameters.regular_insulin(dose_units=10.0)
        p20 = ExogenousInsulinParameters.regular_insulin(dose_units=20.0)
        for t in [10, 50, 100, 200]:
            assert p20.delta_pM_per_min(t) == pytest.approx(2.0 * p10.delta_pM_per_min(t))

    def test_dose_time_shifts_the_curve(self):
        # A dose given 30 minutes before t=0 (dose_time_min=-30) at t=0
        # must match an on-time dose's curve at t=30.
        on_time = ExogenousInsulinParameters.regular_insulin(dose_units=10.0, dose_time_min=0.0)
        early = ExogenousInsulinParameters.regular_insulin(dose_units=10.0, dose_time_min=-30.0)
        assert early.delta_pM_per_min(0.0) == pytest.approx(on_time.delta_pM_per_min(30.0))

    def test_peak_matches_howey_1994_calibration(self):
        # Howey et al. 1994 (PMID:8314011): 10U regular insulin, Tmax=101
        # min, Cmax=308 pM excess above baseline -- this is the isolated-
        # signal-level check (ka^2*S1_0*t*exp(-ka*t)'s own peak); the
        # end-to-end check through the real composed ODE is in
        # test_simulation_api.py.
        p = ExogenousInsulinParameters.regular_insulin(dose_units=10.0)
        ts = [t / 10.0 for t in range(0, 3000)]
        deltas = [p.delta_pM_per_min(t) for t in ts]
        peak_idx = max(range(len(deltas)), key=lambda i: deltas[i])
        # The signal's own peak time is ka^-1 (standard Bateman-function
        # result for two equal-rate-constant compartments), not exactly
        # Howey's 101min Tmax for PLASMA INSULIN -- that's a downstream
        # consequence of this signal integrated against the model's own
        # ke, verified separately in test_simulation_api.py. Here we only
        # confirm the signal's own peak lands at the expected 1/ka.
        assert ts[peak_idx] == pytest.approx(1.0 / INSULIN_REGULAR_KA_PER_MIN, rel=0.02)


class TestExogenousInsulinPathway:
    def test_zero_species(self):
        p = ExogenousInsulinPathway(ExogenousInsulinParameters.regular_insulin(dose_units=10.0))
        assert p.species_names == []
        assert p.y0 == []

    def test_dydt_returns_empty(self):
        p = ExogenousInsulinPathway(ExogenousInsulinParameters.regular_insulin(dose_units=10.0))
        assert p.dydt(0.0, []) == []

    def test_output_signal_key(self):
        p = ExogenousInsulinPathway(ExogenousInsulinParameters.regular_insulin(dose_units=10.0))
        signals = p.output_signals(101.0, [])
        assert "exogenous_insulin_delta" in signals
        assert signals["exogenous_insulin_delta"] > 0.0

    def test_zero_dose_preserves_basal_invariance(self):
        # Rule 11's own requirement: a new/changed signal must keep dG=0
        # and dI=0 at the analytic basal state when the extension itself
        # is at its own basal (here: no dose, dose_units=0.0).
        insulin = ExogenousInsulinPathway(ExogenousInsulinParameters.regular_insulin(dose_units=0.0))
        gi = GlucoseInsulinPathway()
        composite = couple_with_gi(gi, insulin)
        assert composite.validate_y0(tolerance=1e-8) == []

    def test_active_dose_does_not_perturb_the_analytic_basal_state_at_t0(self):
        # Unlike glucose_independent_secretion (sulfonylurea), this signal
        # is exactly zero AT t=0 for a dose given at or after t=0 (it only
        # becomes nonzero strictly after the dose time) -- so the ANALYTIC
        # basal derivatives, evaluated at t=0 with y=y0, must still be
        # exactly zero even with an active dose scheduled at t=0.
        insulin = ExogenousInsulinPathway(ExogenousInsulinParameters.regular_insulin(dose_units=10.0, dose_time_min=0.0))
        gi = GlucoseInsulinPathway()
        composite = couple_with_gi(gi, insulin)
        assert composite.validate_y0(tolerance=1e-8) == []

    def test_couplings_map_signal_straight_through(self):
        couplings = exogenous_insulin_couplings("gi", "insulin")
        assert len(couplings) == 1
        c = couplings[0]
        assert c.source_pathway == "insulin"
        assert c.signal_name == "exogenous_insulin_delta"
        assert c.target_pathway == "gi"
        assert c.target_input == "exogenous_insulin_delta"

    def test_composed_simulation_reproduces_howey_1994_tmax_and_cmax(self):
        # End-to-end reproduction of Howey et al. 1994 (PMID:8314011)
        # through the REAL composed ODE (GlucoseInsulinPathway + this
        # pathway + ODESystem), not just the isolated signal function --
        # confirms the interaction with this model's own ke=0.093/min
        # (Sherwin 1974) produces the calibrated Tmax/Cmax, not only the
        # signal's own peak time (see
        # test_peak_matches_howey_1994_calibration above for that).
        params = ExogenousInsulinParameters.regular_insulin(dose_units=10.0)
        gi = GlucoseInsulinPathway()
        insulin = ExogenousInsulinPathway(params)
        composite = couple_with_gi(gi, insulin)
        system = ODESystem(pathway=composite)
        sim = system.simulate(duration_s=600 * 60)
        I_pM = sim.species["I"] * 1e6
        t_min = sim.times_s / 60.0
        baseline = I_pM[0]
        peak_idx = I_pM.argmax()
        assert t_min[peak_idx] == pytest.approx(101.0, abs=2.0)
        assert (I_pM[peak_idx] - baseline) == pytest.approx(308.0, rel=0.02)

    def test_composed_simulation_reproduces_howey_1994_lispro_tmax_and_cmax(self):
        # Same within-study comparison as the regular-insulin test above
        # (Howey et al. 1994, PMID:8314011) -- arm A (zinc-free LYSPRO,
        # rapid-acting) vs. arm C (regular insulin): "Serum insulin
        # concentrations peaked more than two times higher (Cmax, 698 vs.
        # 308 pM, A vs. C) and in less than half the time (Tmax, 42 vs.
        # 101 min, A vs. C)".
        params = ExogenousInsulinParameters.rapid_acting_insulin(dose_units=10.0)
        gi = GlucoseInsulinPathway()
        insulin = ExogenousInsulinPathway(params)
        composite = couple_with_gi(gi, insulin)
        system = ODESystem(pathway=composite)
        sim = system.simulate(duration_s=600 * 60)
        I_pM = sim.species["I"] * 1e6
        t_min = sim.times_s / 60.0
        baseline = I_pM[0]
        peak_idx = I_pM.argmax()
        assert t_min[peak_idx] == pytest.approx(42.0, abs=2.0)
        assert (I_pM[peak_idx] - baseline) == pytest.approx(698.0, rel=0.02)

    def test_rapid_acting_insulin_uses_distinct_ka_and_pm_per_unit(self):
        assert INSULIN_LISPRO_KA_PER_MIN != INSULIN_REGULAR_KA_PER_MIN
        assert INSULIN_LISPRO_PM_PER_UNIT != INSULIN_REGULAR_PM_PER_UNIT
        p = ExogenousInsulinParameters.rapid_acting_insulin(dose_units=10.0)
        assert p.ka == INSULIN_LISPRO_KA_PER_MIN
        assert p.pm_per_unit == INSULIN_LISPRO_PM_PER_UNIT

    def test_stable_long_trajectory_returns_to_baseline(self):
        # A static y0 check at t=0 alone cannot rule out a subtler drift
        # bug (CLAUDE.md KNOWN BUGS FIXED #10-14's own lesson) -- confirm
        # the composed system, integrated over a long enough window,
        # genuinely returns to the same basal state it started from once
        # the dose has fully decayed.
        params = ExogenousInsulinParameters.regular_insulin(dose_units=10.0)
        gi = GlucoseInsulinPathway()
        insulin = ExogenousInsulinPathway(params)
        composite = couple_with_gi(gi, insulin)
        system = ODESystem(pathway=composite)
        sim = system.simulate(duration_s=3000 * 60)
        I_pM = sim.species["I"] * 1e6
        assert I_pM[-1] == pytest.approx(I_pM[0], abs=1.0)
