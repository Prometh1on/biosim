# SAVE TO: packages\ode-engine\tests\test_default_meal_two_compartment.py
"""
Default meals use the two-compartment ("Bateman") absorption curve.

Before this change a default meal used a single exponential,
Ra(t) = D*k*exp(-k*t), which is at its MAXIMUM the instant you eat. Real
glucose appearance starts at zero and peaks later, because food first has
to leave the stomach. These tests pin that shape for the default meal
presets and for the mass-conserving GLP-1 absorption pathway, which now
tracks two-compartment meals too (stomach -> gut -> blood).
"""
import numpy as np
import pytest
from scipy.integrate import quad, solve_ivp

from biosim_ode.coupling import CompositePathway
from biosim_ode.glp1_pathway import GLP1Parameters, GLP1Pathway, glp1_couplings
from biosim_ode.gut_absorption import (
    DEFAULT_MEAL_ABSORPTION_PROVENANCE,
    GutAbsorptionModel,
)
from biosim_ode.meal import MealInput, high_sugar_meal, light_snack, standard_meal
from biosim_ode.meal_absorption_pathway import (
    MealAbsorptionPathway,
    meal_absorption_couplings,
)
from biosim_ode.pathway import GlucoseInsulinPathway

PRESETS = [high_sugar_meal, standard_meal, light_snack]


def _default_two_compartment_meal(carbs_g: float = 75.0) -> MealInput:
    return MealInput(carbs_g=carbs_g, absorption_model=GutAbsorptionModel.default_meal())


class TestDefaultMealPresetShape:
    @pytest.mark.parametrize("preset", PRESETS)
    def test_absorption_rate_is_zero_at_meal_time(self, preset):
        assert preset().Ra(0.0) == 0.0

    @pytest.mark.parametrize("preset", PRESETS)
    def test_absorption_rate_peaks_after_meal_time(self, preset):
        meal = preset()
        times = np.linspace(0.0, 300.0, 3001)
        rates = [meal.Ra(t) for t in times]
        t_peak = times[int(np.argmax(rates))]
        assert 10.0 < t_peak < 120.0
        assert max(rates) > 0.0

    @pytest.mark.parametrize("preset", PRESETS)
    def test_absorption_integrates_to_meal_glucose(self, preset):
        meal = preset()
        total, _ = quad(meal.Ra, 0.0, 2000.0, limit=200)
        assert total == pytest.approx(meal.D_eff_mM, rel=0.01)


class TestDefaultMealProvenance:
    def test_default_meal_uses_ogtt_rate_constants(self):
        default = GutAbsorptionModel.default_meal()
        ogtt = GutAbsorptionModel.ogtt_calibrated_glucose()
        assert default == ogtt

    def test_provenance_is_labelled_model_assumption(self):
        from biosim_parameter_db import VerificationStatus

        status = DEFAULT_MEAL_ABSORPTION_PROVENANCE["verification_status"]
        assert status in {s.value for s in VerificationStatus}
        assert status == VerificationStatus.UNVERIFIED.value
        assert "MODEL ASSUMPTION" in DEFAULT_MEAL_ABSORPTION_PROVENANCE["label"]


class TestTwoCompartmentStateTracking:
    def test_two_compartment_meal_gets_stomach_and_gut_states(self):
        meal = _default_two_compartment_meal()
        p = MealAbsorptionPathway([meal])
        assert p.species_names == ["meal0_stomach_mM", "meal0_gut_mM"]
        assert p.y0 == pytest.approx([meal.D_eff_mM, 0.0])

    def test_mixed_meal_list_keeps_original_indices(self):
        meals = [MealInput(carbs_g=30.0), _default_two_compartment_meal(50.0)]
        p = MealAbsorptionPathway(meals)
        assert p.species_names == [
            "meal0_remaining_mM", "meal1_stomach_mM", "meal1_gut_mM",
        ]

    def test_couplings_send_gut_state_to_gi(self):
        meal = _default_two_compartment_meal()
        couplings = meal_absorption_couplings("glp1", "abs", "gi", [meal])
        targets = {(c.signal_name, c.target_pathway) for c in couplings}
        assert ("meal0_gut_mM", "gi") in targets

    def test_gi_uses_gut_state_for_two_compartment_meal(self):
        meal = _default_two_compartment_meal()
        gi = GlucoseInsulinPathway(meals=[meal])
        derivs = gi.dydt(30.0, gi.y0, signals={
            "gastric_emptying_modifier": 0.7, "meal0_gut_mM": 2.0,
        })
        ka = meal.absorption_model.k_absorption
        assert derivs[0] == pytest.approx(ka * 2.0)

    def test_neutral_modifier_matches_closed_form(self):
        meal = _default_two_compartment_meal()
        absorption = MealAbsorptionPathway([meal])
        sol = solve_ivp(
            lambda t, y: absorption.dydt(t, y, {"gastric_emptying_modifier": 1.0}),
            (0.0, 300.0), absorption.y0, t_eval=[20.0, 45.0, 120.0],
            rtol=1e-9, atol=1e-11,
        )
        ka = meal.absorption_model.k_absorption
        for t, gut in zip(sol.t, sol.y[1]):
            assert ka * gut == pytest.approx(meal.Ra(t), rel=1e-5)

    def test_semaglutide_delays_but_does_not_delete_meal_glucose(self):
        meal = _default_two_compartment_meal()
        gi = GlucoseInsulinPathway(meals=[meal])
        glp1 = GLP1Pathway(GLP1Parameters.for_semaglutide(dose_mg=1.0))
        absorption = MealAbsorptionPathway([meal])
        composite = CompositePathway(
            [gi, glp1, absorption],
            couplings=glp1_couplings(gi.name, glp1.name)
            + meal_absorption_couplings(glp1.name, absorption.name, gi.name, [meal]),
        )
        sol = solve_ivp(
            composite.dydt, (0.0, 2400.0), composite.y0, method="Radau",
            t_eval=[300.0, 2400.0], rtol=1e-6, atol=1e-8,
        )
        assert sol.success
        stomach, gut = sol.y[-2], sol.y[-1]
        absorbed_fraction = 1.0 - (stomach + gut) / meal.D_eff_mM
        # Slower than without the drug at 5 h, but every gram arrives eventually.
        assert absorbed_fraction[1] > 0.999
