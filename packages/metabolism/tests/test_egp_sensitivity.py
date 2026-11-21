# SAVE TO: packages\metabolism\tests\test_egp_sensitivity.py
"""Tests for biosim_metabolism.egp_sensitivity (Track 4, MASTER_PLAN 4.3c)."""
import pytest

from biosim_metabolism.egp_sensitivity import run_egp_scenarios


class TestRunEgpScenarios:
    def test_returns_exactly_the_four_named_scenarios(self):
        results = run_egp_scenarios()
        names = {r.name for r in results}
        assert names == {
            "fully_unbounded",
            "glycogen_depleted__energy_unbounded",
            "glycogen_depleted__energy_wholebody",
            "glycogen_depleted__energy_hepaticspecific",
        }

    def test_fully_unbounded_scenario_hits_the_artificial_cap(self):
        # Reproduces the original, already-documented finding: with
        # everything permissive, the network's answer is vacuous (the
        # 1000.0 cap itself, not a real stoichiometric result).
        results = {r.name: r for r in run_egp_scenarios()}
        assert results["fully_unbounded"].max_g6pase_umol_kg_min == pytest.approx(1000.0)

    def test_glycogen_depleted_carbon_only_scenario_matches_hand_computed_sum(self):
        # 4.0/2 (lactate) + 2.5/2 (glycerol) + 2.5/2 (alanine) = 4.5 --
        # the exact carbon-precursor-bound sum, confirming the carbon
        # side alone (energy still fully permissive here) determines
        # this scenario's answer.
        results = {r.name: r for r in run_egp_scenarios()}
        assert results["glycogen_depleted__energy_unbounded"].max_g6pase_umol_kg_min == pytest.approx(4.5)

    def test_literature_anchored_energy_bounds_do_not_further_restrict_the_carbon_ceiling(self):
        # The single most important finding from this chunk: even the
        # TIGHTER of the two real, literature-anchored energy scenarios
        # (hepatic-specific O2) does not end up binding below what the
        # (separately, less confidently anchored) carbon-precursor
        # bounds already allow -- all three glycogen-depleted scenarios
        # agree exactly, regardless of how tightly energy is bounded.
        results = {r.name: r for r in run_egp_scenarios()}
        carbon_only = results["glycogen_depleted__energy_unbounded"].max_g6pase_umol_kg_min
        wholebody = results["glycogen_depleted__energy_wholebody"].max_g6pase_umol_kg_min
        hepatic = results["glycogen_depleted__energy_hepaticspecific"].max_g6pase_umol_kg_min
        assert wholebody == pytest.approx(carbon_only)
        assert hepatic == pytest.approx(carbon_only)

    def test_bounded_scenarios_convert_to_a_plausible_mM_per_min_range(self):
        results = {r.name: r for r in run_egp_scenarios()}
        # Hand-computed: 4.5 umol/kg/min * 70kg / 1000 / (0.230*70)
        # = 0.315 mmol/min / 16.1 L = 0.019565... mM/min
        converted = results["glycogen_depleted__energy_wholebody"].max_g6pase_mM_per_min()
        assert converted == pytest.approx(0.019565, abs=1e-5)

    def test_bounded_scenario_ceiling_sits_below_this_platforms_own_basal_and_max_egp(self):
        # A genuinely disclosed comparison, not a forced conclusion: in
        # this specific glycogen-depleted scenario, even the most
        # generous of the three bounded scenarios' ceiling sits BELOW
        # both EGP_BASAL_MM_MIN and EGP_MAX_MM_MIN -- consistent with
        # this scenario representing an extreme, late-fasting,
        # glycogen-exhausted state rather than the shorter-fast/basal
        # state the platform's own constants and the Adkins/Consoli
        # flux benchmarks (packages/metabolism/tests/test_flux_benchmarks.py)
        # actually represent -- these are not directly comparable
        # without accounting for that difference, which this test does
        # not attempt to resolve.
        from biosim_ode.hepatic import EGP_BASAL_MM_MIN, EGP_MAX_MM_MIN
        results = {r.name: r for r in run_egp_scenarios()}
        converted = results["glycogen_depleted__energy_wholebody"].max_g6pase_mM_per_min()
        assert converted < EGP_BASAL_MM_MIN
        assert converted < EGP_MAX_MM_MIN

    def test_each_scenario_has_a_non_empty_description(self):
        for r in run_egp_scenarios():
            assert len(r.description) > 20
