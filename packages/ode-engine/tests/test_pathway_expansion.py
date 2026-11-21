# SAVE TO: packages\ode-engine\tests\test_pathway_expansion.py
"""
Test suite for Sprint 6 -- pathway expansion framework.

Coverage:
    - PathwayDescriptor: construction, tag filtering, validation
    - PathwayRegistry: register, get, duplicate detection, filters
    - GlucoseInsulinPathway: implements BasePathway correctly
    - BooleanThresholds and BooleanState
    - GlucoseInsulinBooleanPathway: state computation, consistency rules
    - WatchdogEngine: fasted confidence high, violations detected correctly
    - CouplingSpec: validation logic
    - CompositePathway: single-pathway wrapping, species naming, slice extraction
    - Regression: existing ODE results unchanged after BasePathway refactor
"""

from typing import Optional

import pytest
import numpy as np

from biosim_ode import (
    GlucoseInsulinPathway,
    ODESystem,
    PathwayParameters,
    SPECIES_NAMES,
    high_sugar_meal,
)
from biosim_ode.base_pathway import BasePathway, PathwayDescriptor, PathwayMaturity
from biosim_ode.boolean_pathway import (
    BooleanState,
    BooleanThresholds,
    ConsistencyRule,
    GlucoseInsulinBooleanPathway,
    WatchdogEngine,
)
from biosim_ode.coupling import CompositePathway, CouplingSpec
from biosim_ode.registry import PathwayRegistry, RegistryError, reset_global_registry


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_descriptor(**kw) -> PathwayDescriptor:
    defaults = dict(
        name="test_pathway_v1",
        version="1.0.0",
        maturity=PathwayMaturity.BETA,
        description="Test pathway",
        species_names=["A", "B"],
        output_signals=["sig_a"],
        input_signals=[],
        dependencies=[],
        tags=["test", "metabolic"],
    )
    defaults.update(kw)
    return PathwayDescriptor(**defaults)


def basal_y() -> list[float]:
    return GlucoseInsulinPathway().params.basal_state()


def high_glucose_y() -> list[float]:
    y = basal_y()
    y[0] = 10.0   # G = 10 mM (high)
    y[1] = 300.0  # I = 300 pM (high)
    y[4] = 0.8    # mTOR = 0.8 (active)
    y[5] = 0.2    # AMPK = 0.2 (suppressed)
    return y


# ---------------------------------------------------------------------------
# PathwayDescriptor
# ---------------------------------------------------------------------------

class TestPathwayDescriptor:
    def test_construction(self):
        d = make_descriptor()
        assert d.name == "test_pathway_v1"
        assert d.n_species == 2

    def test_has_tag_case_insensitive(self):
        d = make_descriptor(tags=["Metabolic", "Glucose"])
        assert d.has_tag("metabolic")
        assert d.has_tag("GLUCOSE")

    def test_has_tag_false(self):
        d = make_descriptor(tags=["metabolic"])
        assert not d.has_tag("cortisol")

    def test_has_dependency(self):
        d = make_descriptor(dependencies=["glucose_insulin_v1"])
        assert d.has_dependency("glucose_insulin_v1")
        assert not d.has_dependency("cortisol_v1")

    def test_empty_name_raises(self):
        with pytest.raises(ValueError):
            make_descriptor(name="")

    def test_empty_version_raises(self):
        with pytest.raises(ValueError):
            make_descriptor(version="")


# ---------------------------------------------------------------------------
# PathwayRegistry
# ---------------------------------------------------------------------------

class TestPathwayRegistry:
    @pytest.fixture(autouse=True)
    def fresh_registry(self):
        reset_global_registry()
        yield
        reset_global_registry()

    def test_register_and_get(self):
        r = PathwayRegistry()
        r.register(GlucoseInsulinPathway, make_descriptor())
        cls = r.get_class("test_pathway_v1")
        assert cls is GlucoseInsulinPathway

    def test_duplicate_same_version_raises(self):
        r = PathwayRegistry()
        r.register(GlucoseInsulinPathway, make_descriptor())
        with pytest.raises(RegistryError, match="already registered"):
            r.register(GlucoseInsulinPathway, make_descriptor())

    def test_different_version_allowed(self):
        r = PathwayRegistry()
        r.register(GlucoseInsulinPathway, make_descriptor(version="1.0.0"))
        r.register(GlucoseInsulinPathway, make_descriptor(version="1.1.0"))
        assert r.get_descriptor("test_pathway_v1").version == "1.1.0"

    def test_get_unknown_raises(self):
        r = PathwayRegistry()
        with pytest.raises(RegistryError):
            r.get_class("nonexistent")

    def test_is_registered(self):
        r = PathwayRegistry()
        assert not r.is_registered("test_pathway_v1")
        r.register(GlucoseInsulinPathway, make_descriptor())
        assert r.is_registered("test_pathway_v1")

    def test_list_all_sorted(self):
        r = PathwayRegistry()
        r.register(GlucoseInsulinPathway, make_descriptor(name="z_pathway"))
        r.register(GlucoseInsulinPathway, make_descriptor(name="a_pathway", version="1.0.0"))
        names = [d.name for d in r.list_all()]
        assert names == sorted(names)

    def test_filter_by_maturity(self):
        r = PathwayRegistry()
        r.register(GlucoseInsulinPathway, make_descriptor(
            name="prod_p", maturity=PathwayMaturity.PRODUCTION))
        r.register(GlucoseInsulinPathway, make_descriptor(
            name="beta_p", version="1.0.0", maturity=PathwayMaturity.BETA))
        prod = r.filter_by_maturity(PathwayMaturity.PRODUCTION)
        assert len(prod) == 1
        assert prod[0].name == "prod_p"

    def test_filter_by_tag(self):
        r = PathwayRegistry()
        r.register(GlucoseInsulinPathway, make_descriptor(
            name="p1", tags=["metabolic", "glucose"]))
        r.register(GlucoseInsulinPathway, make_descriptor(
            name="p2", version="1.0.0", tags=["hormonal"]))
        metabolic = r.filter_by_tag("metabolic")
        assert len(metabolic) == 1
        assert metabolic[0].name == "p1"

    def test_validate_dependencies_missing(self):
        r = PathwayRegistry()
        r.register(GlucoseInsulinPathway, make_descriptor(
            dependencies=["missing_pathway"]))
        missing = r.validate_dependencies("test_pathway_v1")
        assert "missing_pathway" in missing

    def test_validate_dependencies_satisfied(self):
        r = PathwayRegistry()
        r.register(GlucoseInsulinPathway, make_descriptor(name="dep_p"))
        r.register(GlucoseInsulinPathway, make_descriptor(
            name="main_p", version="1.0.0", dependencies=["dep_p"]))
        missing = r.validate_dependencies("main_p")
        assert missing == []

    def test_stats(self):
        r = PathwayRegistry()
        r.register(GlucoseInsulinPathway, make_descriptor(
            maturity=PathwayMaturity.PRODUCTION))
        stats = r.stats()
        assert stats["total_registered"] == 1
        assert stats["by_maturity"]["production"] == 1


# ---------------------------------------------------------------------------
# GlucoseInsulinPathway implements BasePathway
# ---------------------------------------------------------------------------

class TestBasePathwayInterface:
    def test_isinstance_base_pathway(self):
        p = GlucoseInsulinPathway()
        assert isinstance(p, BasePathway)

    def test_descriptor_name(self):
        p = GlucoseInsulinPathway()
        assert p.descriptor.name == "glucose_insulin_v1"

    def test_descriptor_version(self):
        p = GlucoseInsulinPathway()
        assert p.descriptor.version == "1.0.0"

    def test_descriptor_maturity(self):
        p = GlucoseInsulinPathway()
        assert p.maturity == PathwayMaturity.BETA

    def test_descriptor_input_signals_matches_coupling_contract(self):
        # Regression: descriptor.input_signals used to list only 6 of the
        # 9 signals dydt() actually reads via signals.get(...) -- missing
        # p1_modifier, extra_g_clearance_modifier, secretion_scale_modifier
        # (see CLAUDE.md's "EXTENSION-PATHWAY COUPLING CONTRACT", rule 11).
        # Found AGAIN during the exogenous-insulin milestone's scoping pass:
        # renal_glucose_clearance_delta (added for the SGLT2i pathway) was
        # ALSO missing here despite being documented in dydt()'s own
        # docstring -- the same contract-drift slip a second time. Fixed
        # alongside adding exogenous_insulin_delta (the newest signal).
        # "meal*_remaining_mM" (MASTER_PLAN 3.2) added alongside wiring
        # CouplingSpec.validate() into CompositePathway.__init__ -- a
        # wildcard entry representing MealAbsorptionPathway's genuinely
        # dynamic, one-per-meal signal names, which no static list can
        # enumerate (see CouplingSpec.validate()'s own fnmatch support).
        p = GlucoseInsulinPathway()
        expected = {
            "p3_modifier", "gb_modifier", "insulin_resistance_modifier",
            "hepatic_egp_delta", "incretin_effect", "gastric_emptying_modifier",
            "p1_modifier", "extra_g_clearance_modifier", "secretion_scale_modifier",
            "renal_glucose_clearance_delta", "exogenous_insulin_delta",
            "nutrient_secretion_delta", "exercise_glucose_uptake_delta",
            "meal*_remaining_mM", "meal*_gut_mM",
        }
        assert set(p.descriptor.input_signals) == expected

    def test_is_not_production(self):
        p = GlucoseInsulinPathway()
        assert not p.is_production

    def test_output_signals_keys(self):
        p = GlucoseInsulinPathway()
        signals = p.output_signals(0.0, p.y0)
        assert "glucose_mM" in signals
        assert "insulin_pM" in signals
        assert "mTOR_activity" in signals
        assert "AMPK_activity" in signals

    def test_output_signals_at_basal(self):
        p = GlucoseInsulinPathway()
        signals = p.output_signals(0.0, p.y0)
        assert signals["glucose_mM"] == pytest.approx(5.0)
        assert signals["insulin_pM"] == pytest.approx(66.0)

    def test_validate_y0_passes(self):
        p = GlucoseInsulinPathway()
        violations = p.validate_y0(tolerance=1e-3)
        assert violations == [], f"y0 not at steady state: {violations}"

    def test_dydt_accepts_signals_kwarg(self):
        p = GlucoseInsulinPathway()
        derivs = p.dydt(0.0, p.y0, signals={"some_signal": 1.0})
        assert len(derivs) == 6

    def test_unknown_signal_keys_ignored(self):
        # Unrecognised signal keys must not raise -- coupling contracts grow
        # over time and pathways should tolerate signals meant for others.
        p = GlucoseInsulinPathway()
        derivs = p.dydt(0.0, p.y0, signals={"totally_unused_key": 42.0})
        assert derivs == pytest.approx(p.dydt(0.0, p.y0, signals=None))

    def test_neutral_signals_match_no_signals(self):
        p = GlucoseInsulinPathway()
        neutral = {
            "p3_modifier": 1.0, "gb_modifier": 1.0,
            "insulin_resistance_modifier": 1.0, "hepatic_egp_delta": 0.0,
            "incretin_effect": 0.0, "gastric_emptying_modifier": 1.0,
            "p1_modifier": 1.0, "extra_g_clearance_modifier": 1.0,
            "secretion_scale_modifier": 1.0,
        }
        assert p.dydt(0.0, p.y0, signals=neutral) == pytest.approx(
            p.dydt(0.0, p.y0, signals=None)
        )

    def test_p1_and_extra_g_clearance_modifiers_preserve_basal_steady_state(self):
        # Regression for KNOWN BUGS FIXED #16's follow-up: these two new
        # coupling signal keys (added for time-varying exercise couplings)
        # must not perturb the basal steady state even when non-neutral --
        # both are gated on (G - Gb_eff), which is zero at y0.
        p = GlucoseInsulinPathway()
        derivs = p.dydt(0.0, p.y0, signals={
            "p1_modifier": 1.5, "extra_g_clearance_modifier": 2.0,
        })
        assert derivs[0] == pytest.approx(0.0, abs=1e-9)   # dG

    def test_extra_g_clearance_modifier_increases_clearance_away_from_basal(self):
        # Default extra_g_clearance is 0 since B7; set one explicitly so the
        # modifier has something to scale.
        p = GlucoseInsulinPathway(params=PathwayParameters(extra_g_clearance=0.018))
        y = p.y0.copy()
        y[0] = p.params.Gb + 2.0   # elevated glucose -> nonzero clearance term
        dG_full = p.dydt(0.0, y, signals=None)[0]
        dG_boosted = p.dydt(0.0, y, signals={"extra_g_clearance_modifier": 2.0})[0]
        assert dG_boosted < dG_full   # faster clearance -> more negative dG

    def test_p1_modifier_increases_clearance_away_from_basal(self):
        p = GlucoseInsulinPathway()
        y = p.y0.copy()
        y[0] = p.params.Gb + 2.0
        dG_full = p.dydt(0.0, y, signals=None)[0]
        dG_boosted = p.dydt(0.0, y, signals={"p1_modifier": 2.0})[0]
        assert dG_boosted < dG_full

    def test_basal_steady_state_preserved_with_signals_present(self):
        # Even non-neutral incretin_effect / p3_modifier must not perturb the
        # analytic basal steady state -- both terms are gated on deviation
        # from Gb/G_threshold, which is zero at y0.
        p = GlucoseInsulinPathway()
        derivs = p.dydt(0.0, p.y0, signals={
            "p3_modifier": 0.5, "insulin_resistance_modifier": 0.7,
            "incretin_effect": 0.9,
        })
        assert derivs[0] == pytest.approx(0.0, abs=1e-9)   # dG
        assert derivs[1] == pytest.approx(0.0, abs=1e-9)   # dI

    def test_renal_glucose_clearance_delta_defaults_to_neutral(self):
        # renal_glucose_clearance_delta (rule 11, added for the SGLT2-
        # inhibitor pathway -- see renal_pathway.py) must default to 0.0
        # like every other signal key, so GlucoseInsulinPathway behaves
        # identically whether or not RenalPathway is composed.
        p = GlucoseInsulinPathway()
        dG_absent = p.dydt(0.0, p.y0, signals=None)[0]
        dG_present_but_neutral = p.dydt(0.0, p.y0, signals={"renal_glucose_clearance_delta": 0.0})[0]
        assert dG_absent == pytest.approx(dG_present_but_neutral)
        assert dG_absent == pytest.approx(0.0, abs=1e-9)

    def test_renal_glucose_clearance_delta_is_additive_to_dG(self):
        p = GlucoseInsulinPathway()
        dG_no_delta = p.dydt(0.0, p.y0, signals=None)[0]
        dG_with_delta = p.dydt(0.0, p.y0, signals={"renal_glucose_clearance_delta": -0.02})[0]
        assert dG_with_delta == pytest.approx(dG_no_delta - 0.02)

    def test_p3_modifier_reduces_X_response_to_insulin_excess(self):
        p = GlucoseInsulinPathway()
        y = p.y0.copy()
        y[1] = p.params.Ib + 50.0   # elevated insulin -> nonzero X response
        dX_full = p.dydt(0.0, y, signals=None)[2]
        dX_half = p.dydt(0.0, y, signals={"p3_modifier": 0.5})[2]
        assert dX_half == pytest.approx(dX_full * 0.5, rel=1e-9)

    def test_p3_eff_is_clamped_at_the_point_of_use_inside_dydt(self):
        # Security-reviewer finding (Section 8 housekeeping, found while
        # re-deriving CORTISOL_P3_SENSITIVITY): personalizer.py's static
        # _PARAM_BOUNDS/clamp_warnings mechanism only ever bounds
        # PathwayParameters.p3 BEFORE integration starts -- it has no
        # visibility into ODE-COUPLED multiplicative signals (p3_modifier,
        # insulin_resistance_modifier), which are re-applied at every
        # solver step and could silently drive the real, integrated p3_eff
        # below the platform's own documented floor (1e-8) with zero
        # disclosure. Confirmed by direct execution before this fix: an
        # extreme realistic composition drove real p3_eff down to ~4.55e-9.
        # Fixed by clamping p3_eff to (P3_EFF_FLOOR, P3_EFF_CEILING) at the
        # point of use inside dydt() -- this test uses an extreme p3_modifier
        # directly (0.005, well below anything any single shipped mechanism
        # produces alone) purely to exercise the clamp itself, matching this
        # project's own "synthetic value to exercise the code path" convention.
        p = GlucoseInsulinPathway()
        y = p.y0.copy()
        y[1] = p.params.Ib + 1000.0   # large elevated insulin excess
        raw_p3_eff = p.params.p3 * 0.005
        assert raw_p3_eff < 1e-8  # confirm this test scenario would breach the floor
        dX = p.dydt(0.0, y, signals={"p3_modifier": 0.005})[2]
        # dX = -p2*X + p3_eff*(I-Ib); X=0 at y0, so dX == p3_eff*(I-Ib) exactly.
        expected_dX_if_clamped = 1e-8 * 1000.0
        expected_dX_if_unclamped = raw_p3_eff * 1000.0
        assert dX == pytest.approx(expected_dX_if_clamped, rel=1e-6)
        assert dX != pytest.approx(expected_dX_if_unclamped, rel=1e-3)

    def test_gb_eff_is_clamped_at_the_point_of_use_inside_dydt(self):
        # Section 8 housekeeping follow-on to the p3_eff clamp above --
        # the identical unclamped-at-point-of-use gap exists for the other
        # three rule-11 multiplicative signals. dG is linear in Gb_eff
        # (coefficient -total_clearance_rate), so an extreme modifier's
        # dydt() output must equal the SAME call with the signal set
        # directly to the ceiling value if the clamp engages correctly.
        p = GlucoseInsulinPathway()
        y = p.y0.copy()
        raw_gb_eff = p.params.Gb * 100.0
        assert raw_gb_eff > 15.0  # confirm this scenario would breach GB_EFF_CEILING
        dG_extreme = p.dydt(0.0, y, signals={"gb_modifier": 100.0})[0]
        dG_at_ceiling = p.dydt(0.0, y, signals={"gb_modifier": 15.0 / p.params.Gb})[0]
        assert dG_extreme == pytest.approx(dG_at_ceiling, rel=1e-9)

    def test_p1_eff_is_clamped_at_the_point_of_use_inside_dydt(self):
        p = GlucoseInsulinPathway()
        y = p.y0.copy()
        y[0] = 10.0  # G elevated above Gb so p1_eff's term is non-zero
        raw_p1_eff = p.params.p1 * 1000.0
        assert raw_p1_eff > 0.20  # confirm this scenario would breach P1_EFF_CEILING
        # extra_g_clearance_modifier=0.0 here floors to EXTRA_G_CLEARANCE_EFF_FLOOR
        # (0.001), not true zero -- fine for this test, since the SAME floored
        # constant appears identically in both calls below and cancels out;
        # it isolates the comparison's SENSITIVITY to p1_modifier alone, not a
        # literal zeroing of the other term.
        signals_extreme = {"p1_modifier": 1000.0, "extra_g_clearance_modifier": 0.0}
        signals_at_ceiling = {"p1_modifier": 0.20 / p.params.p1, "extra_g_clearance_modifier": 0.0}
        dG_extreme = p.dydt(0.0, y, signals=signals_extreme)[0]
        dG_at_ceiling = p.dydt(0.0, y, signals=signals_at_ceiling)[0]
        assert dG_extreme == pytest.approx(dG_at_ceiling, rel=1e-9)

    def test_extra_g_clearance_eff_is_clamped_at_the_point_of_use_inside_dydt(self):
        p = GlucoseInsulinPathway(params=PathwayParameters(extra_g_clearance=0.018))
        y = p.y0.copy()
        y[0] = 10.0  # G elevated above Gb so extra_g_clearance_eff's term is non-zero
        raw_eff = p.params.extra_g_clearance * 100.0
        assert raw_eff > 0.20  # confirm this scenario would breach the ceiling
        # p1_modifier=0.0 here floors to P1_EFF_FLOOR (0.001), not true zero --
        # fine since the same floored constant appears identically in both
        # calls below and cancels out; isolates sensitivity to
        # extra_g_clearance_modifier alone, not a literal zeroing of p1_eff.
        signals_extreme = {"extra_g_clearance_modifier": 100.0, "p1_modifier": 0.0}
        signals_at_ceiling = {
            "extra_g_clearance_modifier": 0.20 / p.params.extra_g_clearance,
            "p1_modifier": 0.0,
        }
        dG_extreme = p.dydt(0.0, y, signals=signals_extreme)[0]
        dG_at_ceiling = p.dydt(0.0, y, signals=signals_at_ceiling)[0]
        assert dG_extreme == pytest.approx(dG_at_ceiling, rel=1e-9)

    def test_secretion_scale_eff_is_clamped_at_the_point_of_use_inside_dydt(self):
        p = GlucoseInsulinPathway()
        y = p.y0.copy()
        y[0] = 10.0  # G elevated above Gb so delta_secretion is non-zero
        raw_eff = p.params.secretion_scale * 100.0
        assert raw_eff > 3.0  # confirm this scenario would breach SECRETION_SCALE_EFF_CEILING
        dI_extreme = p.dydt(0.0, y, signals={"secretion_scale_modifier": 100.0})[1]
        dI_at_ceiling = p.dydt(0.0, y, signals={"secretion_scale_modifier": 3.0 / p.params.secretion_scale})[1]
        assert dI_extreme == pytest.approx(dI_at_ceiling, rel=1e-9)

    def test_secretion_scale_modifier_preserves_basal_steady_state(self):
        # Regression: during_exercise_modifiers()["secretion_scale"] was
        # computed by biosim_physiology.exercise but had no signal slot to
        # plug into -- ExercisePathway.output_signals() dropped it
        # entirely, so an in-progress exercise bout's catecholamine-driven
        # insulin secretion suppression never reached the ODE. Both terms
        # this modifier scales (delta_secretion, first_phase_secretion)
        # are already zero at the basal steady state by construction
        # (current_frac == basal_frac and dG <= 0 respectively), so a
        # non-neutral value must not perturb it.
        p = GlucoseInsulinPathway()
        derivs = p.dydt(0.0, p.y0, signals={"secretion_scale_modifier": 0.5})
        assert derivs[1] == pytest.approx(0.0, abs=1e-9)   # dI

    def test_secretion_scale_modifier_suppresses_secretion_away_from_basal(self):
        p = GlucoseInsulinPathway()
        y = p.y0.copy()
        y[0] = p.params.Gb + 2.0   # elevated glucose -> nonzero delta_secretion
        dI_full = p.dydt(0.0, y, signals=None)[1]
        dI_suppressed = p.dydt(0.0, y, signals={"secretion_scale_modifier": 0.5})[1]
        assert dI_suppressed < dI_full

    def test_gastric_emptying_modifier_slows_glucose_appearance(self):
        from biosim_ode.meal import high_sugar_meal
        p = GlucoseInsulinPathway(meals=[high_sugar_meal()])
        y = p.y0
        dG_full = p.dydt(5.0, y, signals=None)[0]
        dG_slowed = p.dydt(5.0, y, signals={"gastric_emptying_modifier": 0.5})[0]
        assert dG_slowed < dG_full

    def test_hepatic_egp_delta_adds_directly_to_dG(self):
        p = GlucoseInsulinPathway()
        dG_base = p.dydt(0.0, p.y0, signals=None)[0]
        dG_egp  = p.dydt(0.0, p.y0, signals={"hepatic_egp_delta": 0.02})[0]
        assert dG_egp == pytest.approx(dG_base + 0.02)

    def test_gb_modifier_shifts_glucose_clearance_target(self):
        p = GlucoseInsulinPathway()
        # At y0 (G == nominal Gb), raising the effective Gb should make dG
        # rise (clearance term now pulls toward a higher target).
        dG_base  = p.dydt(0.0, p.y0, signals=None)[0]
        dG_raised = p.dydt(0.0, p.y0, signals={"gb_modifier": 1.2})[0]
        assert dG_raised > dG_base

    def test_name_property(self):
        p = GlucoseInsulinPathway()
        assert p.name == p.descriptor.name


# ---------------------------------------------------------------------------
# Boolean pathway
# ---------------------------------------------------------------------------

class TestBooleanPathway:
    def test_basal_state_not_high(self):
        bp = GlucoseInsulinBooleanPathway()
        state = bp.compute_state(basal_y())
        assert not state.G_high
        assert not state.I_high

    def test_high_glucose_state(self):
        bp = GlucoseInsulinBooleanPathway()
        state = bp.compute_state(high_glucose_y())
        assert state.G_high
        assert state.I_high
        assert state.mTOR_active
        assert state.AMPK_suppressed

    def test_as_dict_keys(self):
        bp = GlucoseInsulinBooleanPathway()
        d = bp.compute_state(basal_y()).as_dict()
        assert "G_high" in d
        assert "I_high" in d
        assert "mTOR_active" in d
        assert "AMPK_suppressed" in d

    def test_no_violations_at_consistent_state(self):
        bp = GlucoseInsulinBooleanPathway()
        violations = bp.check_consistency(bp.compute_state(high_glucose_y()))
        assert violations == []

    def test_violation_detected_broken_mtor(self):
        bp = GlucoseInsulinBooleanPathway()
        y = high_glucose_y()
        y[4] = 0.1  # mTOR very low despite high insulin — violation
        state = bp.compute_state(y)
        violations = bp.check_consistency(state)
        descs = [v.description for v in violations]
        assert any("mTOR" in d for d in descs)

    def test_confidence_score_1_no_violations(self):
        bp = GlucoseInsulinBooleanPathway()
        assert bp.confidence_score([]) == pytest.approx(1.0)

    def test_confidence_score_decreases_with_violations(self):
        bp = GlucoseInsulinBooleanPathway()
        rule = ConsistencyRule("a", "b", "test", severity=0.25)
        assert bp.confidence_score([rule]) == pytest.approx(0.75)

    def test_confidence_score_floored_at_zero(self):
        bp = GlucoseInsulinBooleanPathway()
        rules = [ConsistencyRule("a", "b", f"rule{i}", severity=0.5) for i in range(5)]
        assert bp.confidence_score(rules) == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# WatchdogEngine
# ---------------------------------------------------------------------------

class TestWatchdogEngine:
    def test_fasted_confidence_high(self):
        engine = WatchdogEngine()
        result = engine.evaluate(basal_y(), t_min=0.0)
        assert result.confidence_score == pytest.approx(1.0)

    def test_consistent_high_state_confidence_high(self):
        engine = WatchdogEngine()
        result = engine.evaluate(high_glucose_y(), t_min=30.0)
        assert result.confidence_score == pytest.approx(1.0)

    def test_broken_mtor_lowers_confidence(self):
        engine = WatchdogEngine()
        y = high_glucose_y()
        y[4] = 0.1  # mTOR off despite high insulin
        result = engine.evaluate(y, t_min=30.0)
        assert result.confidence_score < 1.0

    def test_disagreements_list(self):
        engine = WatchdogEngine()
        y = high_glucose_y()
        y[4] = 0.1
        result = engine.evaluate(y, t_min=30.0)
        assert len(result.disagreements) > 0

    def test_evaluate_trajectory(self):
        engine = WatchdogEngine()
        p = GlucoseInsulinPathway()
        t = np.linspace(0, 60, 10)
        y_arr = np.tile(np.array(p.y0)[:, np.newaxis], (1, 10))
        results = engine.evaluate_trajectory(t, y_arr)
        assert len(results) == 10

    def test_overall_confidence_fasted(self):
        engine = WatchdogEngine()
        p = GlucoseInsulinPathway()
        t = np.linspace(0, 60, 20)
        y_arr = np.tile(np.array(p.y0)[:, np.newaxis], (1, 20))
        results = engine.evaluate_trajectory(t, y_arr)
        conf = engine.overall_confidence(results)
        assert conf == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# CouplingSpec
# ---------------------------------------------------------------------------

class TestCouplingSpec:
    def test_valid_coupling(self):
        c = CouplingSpec(
            source_pathway="glucose_insulin_v1",
            signal_name="insulin_pM",
            target_pathway="cortisol_v1",
            target_input="insulin_suppression",
        )
        errors = c.validate(
            available_pathways=["glucose_insulin_v1", "cortisol_v1"],
            source_signals=["insulin_pM", "glucose_mM"],
            target_inputs=["insulin_suppression"],
        )
        assert errors == []

    def test_self_coupling_invalid(self):
        c = CouplingSpec(
            source_pathway="p1",
            signal_name="sig",
            target_pathway="p1",
            target_input="inp",
        )
        errors = c.validate(["p1"], ["sig"], ["inp"])
        assert any("cannot couple to itself" in e for e in errors)

    def test_missing_source_pathway(self):
        c = CouplingSpec("missing", "sig", "p2", "inp")
        errors = c.validate(["p2"], ["sig"], ["inp"])
        assert any("missing" in e for e in errors)

    def test_missing_target_input_is_caught(self):
        # Regression: validate() used to accept target_inputs as a
        # parameter but never actually check self.target_input against
        # it -- a typo'd target_input in any pathway-extension wiring
        # would never be caught.
        c = CouplingSpec("p1", "sig", "p2", "typo_d_input")
        errors = c.validate(["p1", "p2"], ["sig"], ["correct_input"])
        assert any("typo_d_input" in e for e in errors)

    def test_transform_applied(self):
        c = CouplingSpec("p1", "sig", "p2", "inp", transform=lambda x: x * 2.0)
        assert c.transform(5.0) == pytest.approx(10.0)


# ---------------------------------------------------------------------------
# CompositePathway
# ---------------------------------------------------------------------------

class TestCompositePathway:
    def test_single_pathway_composite(self):
        p = GlucoseInsulinPathway()
        cp = CompositePathway([p])
        assert len(cp.y0) == len(p.y0)

    def test_species_names_prefixed(self):
        p = GlucoseInsulinPathway()
        cp = CompositePathway([p])
        for name in cp.species_names:
            assert name.startswith(p.name + ".")

    def test_single_pathway_dydt_equivalent(self):
        p = GlucoseInsulinPathway()
        cp = CompositePathway([p])
        y = p.y0
        direct = p.dydt(0.0, y)
        composite = cp.dydt(0.0, y)
        assert direct == pytest.approx(composite, rel=1e-9)

    def test_slice_extraction(self):
        p = GlucoseInsulinPathway()
        cp = CompositePathway([p])
        y = p.y0
        extracted = cp.extract_pathway_state(p.name, y)
        assert extracted == pytest.approx(y)

    def test_composite_pathway_names(self):
        p = GlucoseInsulinPathway()
        cp = CompositePathway([p])
        assert p.name in cp.pathway_names

    def test_empty_pathways_raises(self):
        with pytest.raises(ValueError):
            CompositePathway([])

    def test_output_signals_namespaced(self):
        p = GlucoseInsulinPathway()
        cp = CompositePathway([p])
        signals = cp.output_signals(0.0, p.y0)
        for key in signals:
            assert key.startswith(p.name + ".")

    def test_worst_maturity_propagates(self):
        p = GlucoseInsulinPathway()
        cp = CompositePathway([p])
        assert cp.maturity == PathwayMaturity.BETA

    def test_with_meals_applies_to_gi_subpathway(self):
        cp = CompositePathway([GlucoseInsulinPathway()])
        cp2 = cp.with_meals([high_sugar_meal()])
        # Default meals are two-compartment (Ra=0 at the moment of eating),
        # so look 30 min in.
        assert cp2.dydt(30.0, cp2.y0)[0] != pytest.approx(0.0)

    def test_with_meals_preserves_couplings(self):
        from biosim_ode.cortisol_pathway import CortisolPathway
        couplings = [
            CouplingSpec("cortisol_v1", "p3_modifier", "glucose_insulin_v1", "p3_modifier")
        ]
        cp = CompositePathway([GlucoseInsulinPathway(), CortisolPathway()], couplings=couplings)
        cp2 = cp.with_meals([high_sugar_meal()])
        assert cp2.couplings == couplings

    def test_with_meals_no_meal_capable_pathway_raises(self):
        from biosim_ode.cortisol_pathway import CortisolPathway
        cp = CompositePathway([CortisolPathway()])
        with pytest.raises(ValueError, match="no meal-capable pathway"):
            cp.with_meals([high_sugar_meal()])


class TestCortisolP3SensitivityFix:
    """
    Section 8 housekeeping fix (found during Track 1 row 1.38/Cushing's
    syndrome scoping): CORTISOL_P3_SENSITIVITY was uncited (0.25) and had
    an unreachable-floor consequence no one had checked for -- p3_modifier
    could never fall below ~0.9296 regardless of cortisol level, because
    excess_gr is mathematically bounded above by (1 - GR_REFERENCE) via
    the Hill-saturation shape of gr_occ. Re-derived from Rizza et al. 1982
    (PMID:7033265, already cited in this module) -- see cortisol_pathway.py's
    own citation comment for the full derivation.
    """

    def test_p3_modifier_floor_is_no_longer_capped_near_the_old_ceiling(self):
        from biosim_ode.cortisol_pathway import CortisolPathway, GR_REFERENCE
        pathway = CortisolPathway()
        # Manually elevated GR occupancy approaching its own asymptote (1.0)
        # -- the old, uncited CORTISOL_P3_SENSITIVITY=0.25 could never push
        # p3_modifier below ~0.9296 here; the fix should reach well below it.
        y = [1.0, 5.0, 0.99999]
        sig = pathway.output_signals(0.0, y)
        assert sig["p3_modifier"] < 0.5
        # Sanity: still respects the function's own hardcoded 0.3 safety
        # floor as an absolute backstop (unchanged by this fix).
        assert sig["p3_modifier"] >= 0.3

    def test_p3_modifier_matches_rizza_1982_anchor_point(self):
        # At Rizza et al. 1982's own measured cortisol ratio (37/14=2.643,
        # the ACTUAL computed ratio, not the abstract's own imprecise
        # "~4-fold" phrasing), composing through this pathway should
        # reproduce the paper's own peripheral ED50 ratio (64/104=0.6154)
        # -- this is the anchor CORTISOL_P3_SENSITIVITY was solved FOR, so
        # this test pins that derivation as a permanent regression.
        from biosim_ode.cortisol_pathway import CortisolPathway, CortisolParameters, CORTISOL_BASAL_NM
        rizza_ratio = 37.0 / 14.0
        params = CortisolParameters(cortisol_basal=CORTISOL_BASAL_NM)
        pathway = CortisolPathway(params=params)
        elevated_cortisol = CORTISOL_BASAL_NM * rizza_ratio
        # gr_occ at this cortisol level, computed the same way output_signals() does
        from biosim_ode.kinetics import hill_activation
        from biosim_ode.cortisol_pathway import EC50_GR_NM, HILL_GR
        gr_occ = hill_activation(elevated_cortisol, EC50_GR_NM, HILL_GR)
        sig = pathway.output_signals(0.0, [elevated_cortisol, 5.0, gr_occ])
        assert sig["p3_modifier"] == pytest.approx(64.0 / 104.0, abs=1e-9)

    def test_basal_steady_state_still_exactly_neutral(self):
        # The fix must not disturb the pathway's own fixed point -- excess_gr
        # is exactly 0.0 at basal regardless of CORTISOL_P3_SENSITIVITY's
        # value, so p3_modifier must still be exactly 1.0 at y0.
        from biosim_ode.cortisol_pathway import CortisolPathway
        pathway = CortisolPathway()
        sig = pathway.output_signals(0.0, pathway.y0)
        assert sig["p3_modifier"] == pytest.approx(1.0)
        assert pathway.validate_y0(tolerance=1e-8) == []

    def test_prednisone_dose_now_produces_a_real_p3_effect_not_a_trivial_one(self):
        # Before the fix, ANY prednisone dose produced at most a ~7%
        # p3 reduction (the old unreachable-floor bug) -- now a clinically
        # realistic dose produces a real, substantial reduction.
        from biosim_ode.cortisol_pathway import CortisolPathway, CortisolParameters
        from biosim_ode.kinetics import hill_activation
        from biosim_ode.cortisol_pathway import EC50_GR_NM, HILL_GR, GR_REFERENCE
        params = CortisolParameters.for_prednisone(dose_mg=30.0)
        cortisol_effective = params.cortisol_basal + params.exogenous_gc_nM
        gr_occ = hill_activation(cortisol_effective, EC50_GR_NM, HILL_GR)
        pathway = CortisolPathway(params=params)
        sig = pathway.output_signals(0.0, [cortisol_effective, 5.0, gr_occ])
        assert sig["p3_modifier"] < 0.6

    def test_fix_reaches_a_real_composed_ode_integration_not_just_the_bare_function(self):
        # Mandatory code-reviewer finding (MEDIUM): the other tests in this
        # class all call output_signals() directly on an isolated
        # CortisolPathway -- none exercised the fix through an actual ODE
        # solve. This one does: builds the SAME couple_with_gi() composition
        # SimulationRunner._detect_cortisol()/_build_pathway() uses in
        # production, runs a real ODESystem integration with a dose already
        # 6h into its PD-onset window (matching the fully-developed-effect
        # convention test_glucocorticoid_fasted_does_not_double_count_
        # static_and_pathway uses), and confirms the ACTUAL trajectory's own
        # Cortisol/GR_occ state -- not a hand-constructed one -- produces a
        # p3_modifier well below the old ~0.93 unreachable-floor ceiling.
        from biosim_ode.cortisol_pathway import CortisolPathway, CortisolParameters, couple_with_gi
        params = CortisolParameters.for_prednisone(dose_mg=30.0, dose_time_min=-360.0)
        composite = couple_with_gi(GlucoseInsulinPathway(), CortisolPathway(params=params))
        result = ODESystem(pathway=composite).simulate(duration_s=4 * 3600, meals=[])
        cortisol_traj = result.species["Cortisol"]
        gr_occ_traj = result.species["GR_occ"]
        final_signals = CortisolPathway(params=params).output_signals(
            4 * 3600 / 60.0, [cortisol_traj[-1], 5.0, gr_occ_traj[-1]],
        )
        assert final_signals["p3_modifier"] < 0.6

class TestCortisolCrhCoupledFixedPointFix:
    """
    Section 8 housekeeping (found during a second attempt at Track 1 row
    1.38/Cushing's syndrome): CortisolPathway's OLD basal_state() returned
    [cortisol_basal, crh_basal, gr_occ_basal] literally, a genuine fixed
    point of dydt() only in the trivial default case (stress_level=0).
    Any stress_level>0 drifted from t=0 onward, silently, since this
    pathway was first built. Fixed by solve_cortisol_crh_steady_state()
    -- see that function's own docstring for the full derivation and the
    first (wrong) fix attempt it corrects.
    """

    def test_stress_level_alone_is_now_a_genuine_fixed_point(self):
        from biosim_ode.cortisol_pathway import CortisolParameters, CortisolPathway
        for stress_level in [0.0, 0.3, 0.8, 1.0]:
            pathway = CortisolPathway(CortisolParameters(stress_level=stress_level))
            assert pathway.validate_y0(tolerance=1e-8) == [], (
                f"stress_level={stress_level} is not a genuine fixed point"
            )

    def test_default_case_is_byte_identical_to_before_the_fix(self):
        # Backward-compatibility check: the fix must not change ANY
        # existing (non-stress) scenario's y0 at all.
        from biosim_ode.cortisol_pathway import (
            CortisolParameters, CORTISOL_BASAL_NM, CRH_BASAL_PM, GR_REFERENCE,
        )
        params = CortisolParameters()
        y0 = params.basal_state()
        assert y0[0] == pytest.approx(CORTISOL_BASAL_NM, abs=1e-12)
        assert y0[1] == pytest.approx(CRH_BASAL_PM, abs=1e-12)
        assert y0[2] == pytest.approx(GR_REFERENCE, abs=1e-9)

    def test_stress_still_genuinely_elevates_cortisol_not_cancelled_out(self):
        # The exact regression the first (reverted) fix attempt
        # introduced: normalizing cortisol_production against crh's own
        # target made stress_level's steady-state effect on cortisol
        # vanish entirely. This confirms the CURRENT fix does not repeat
        # that mistake -- cortisol genuinely rises with stress_level.
        from biosim_ode.cortisol_pathway import CortisolParameters
        default_cortisol = CortisolParameters().basal_state()[0]
        stressed_cortisol = CortisolParameters(stress_level=0.8).basal_state()[0]
        assert stressed_cortisol > default_cortisol * 1.5  # a real, substantial elevation

    def test_stress_level_fixed_point_is_stable_over_a_long_trajectory(self):
        # Closes the exact blind spot that let the original bug go
        # undetected: every existing test only checked the glucose peak
        # shape, never a long-trajectory invariant for stress_level alone.
        from biosim_ode.cortisol_pathway import CortisolParameters, CortisolPathway, couple_with_gi
        composite = couple_with_gi(GlucoseInsulinPathway(), CortisolPathway(CortisolParameters(stress_level=0.8)))
        result = ODESystem(pathway=composite).simulate(duration_s=7 * 24 * 3600, meals=[])
        cortisol = result.species["Cortisol"]
        assert cortisol.max() == pytest.approx(cortisol.min(), rel=1e-6)
        assert cortisol[0] == pytest.approx(cortisol[-1], rel=1e-9)

    def test_prednisone_alone_still_unaffected_by_this_fix(self):
        from biosim_ode.cortisol_pathway import CortisolParameters, CortisolPathway
        pathway = CortisolPathway(CortisolParameters.for_prednisone(30.0))
        assert pathway.validate_y0(tolerance=1e-8) == []

    def test_nonpositive_cortisol_basal_raises_a_clear_error(self):
        # Security-reviewer finding: this guard existed but had zero
        # regression-test coverage anywhere in this class -- a future
        # refactor could silently remove it and the full suite would
        # still pass.
        from biosim_ode.cortisol_pathway import solve_cortisol_crh_steady_state
        with pytest.raises(ValueError, match="cortisol_basal"):
            solve_cortisol_crh_steady_state(0.0, 5.0, 0.0)
        with pytest.raises(ValueError, match="cortisol_basal"):
            solve_cortisol_crh_steady_state(-0.1, 5.0, 0.0)

    def test_negative_stress_level_raises_a_clear_error_not_a_confusing_one(self):
        # Security-reviewer finding (MEDIUM): CortisolParameters.__init__
        # clamps stress_level to [0,1] before it ever reaches this
        # function, but solve_cortisol_crh_steady_state() is public and
        # callable directly (as this test does) -- an un-clamped negative
        # value used to crash deep inside hill_activation() with an
        # unrelated-looking "x must be non-negative" error.
        from biosim_ode.cortisol_pathway import solve_cortisol_crh_steady_state
        with pytest.raises(ValueError, match="stress_level"):
            solve_cortisol_crh_steady_state(0.28, 5.0, -0.5)

    def test_extreme_cortisol_basal_raises_a_clear_error_not_a_raw_scipy_one(self):
        # Mandatory code-reviewer finding (MEDIUM): brentq's own bracket
        # precondition can fail for an implausibly large cortisol_basal
        # (e.g. a unit-confusion bug passing a pM/uM value where nM is
        # expected) -- confirmed this used to surface as scipy's own raw
        # "f(a) and f(b) must have different signs", not a helpful
        # domain error naming the actual input.
        from biosim_ode.cortisol_pathway import CortisolParameters
        with pytest.raises(ValueError, match="cortisol_basal"):
            CortisolParameters(cortisol_basal=1e6).basal_state()

    def test_stress_plus_chronic_dose_residual_is_smaller_than_pre_fix_and_disclosed(self):
        # Mandatory code-reviewer finding (MEDIUM): solve_cortisol_crh_
        # steady_state() doesn't fold in exogenous_gc_nM at t=0, so
        # composing a large-negative dose_time_min (a dose already well
        # into its onset window) WITH stress_level>0 still leaves a
        # small residual -- a genuine, disclosed, NOT-fully-closed
        # limitation (see for_prednisone()'s own docstring). This pins
        # the CURRENT (improved, not perfect) behavior as a regression:
        # the residual must stay small and must NOT regress back to the
        # much larger pre-fix magnitude (CRH~=0.347, GR_occ~=0.0137).
        from biosim_ode.cortisol_pathway import CortisolParameters, CortisolPathway
        params = CortisolParameters.for_prednisone(30.0, dose_time_min=-1800.0, stress_level=0.5)
        pathway = CortisolPathway(params)
        violations = pathway.validate_y0(tolerance=1e-8)
        assert len(violations) == 1
        assert "GR_occ" in violations[0]
        assert "CRH" not in " ".join(violations)  # the CRH residual IS now fully closed
        dydt = pathway.dydt(0.0, pathway.y0)
        assert abs(dydt[2]) < 0.01  # GR_occ residual stays small (was ~0.0137 pre-fix)


class TestCushingsSyndrome:
    """
    Track 1 row 1.38, third implementation attempt, on the now-fixed
    CRH/Cortisol coupled-fixed-point foundation. See
    CortisolParameters.for_cushings()/CUSHINGS_CORTISOL_RATIO and
    solve_cortisol_basal_for_target_steady_state() for the full
    citation and derivation.
    """

    def test_achieved_cortisol_ratio_matches_the_literature_target(self):
        # The whole point of the inverse solve: a naive direct multiply
        # (cortisol_basal = CORTISOL_BASAL_NM*1.560) only ACHIEVES a
        # ~1.287x steady state due to the coupled negative feedback --
        # this confirms the inverse solve genuinely closes that gap.
        # NOTE: cortisol_basal itself is the SOLVED INPUT, not the achieved
        # steady state -- check the achieved steady state via basal_state().
        from biosim_ode.cortisol_pathway import CortisolParameters, CORTISOL_BASAL_NM
        params = CortisolParameters.for_cushings()
        y0 = params.basal_state()
        achieved_cortisol_ss_ratio = y0[0] / CORTISOL_BASAL_NM
        assert achieved_cortisol_ss_ratio == pytest.approx(
            CortisolParameters.CUSHINGS_CORTISOL_RATIO, rel=1e-6,
        )

    def test_result_is_a_genuine_fixed_point(self):
        from biosim_ode.cortisol_pathway import CortisolParameters, CortisolPathway
        pathway = CortisolPathway(CortisolParameters.for_cushings())
        assert pathway.validate_y0(tolerance=1e-8) == []

    def test_stable_over_a_long_trajectory(self):
        # The exact check that caught the second (reverted) attempt's
        # failure -- confirms this one doesn't repeat it.
        from biosim_ode.cortisol_pathway import CortisolParameters, CortisolPathway, couple_with_gi
        composite = couple_with_gi(GlucoseInsulinPathway(), CortisolPathway(CortisolParameters.for_cushings()))
        result = ODESystem(pathway=composite).simulate(duration_s=7 * 24 * 3600, meals=[])
        cortisol = result.species["Cortisol"]
        assert cortisol.max() == pytest.approx(cortisol.min(), rel=1e-6)

    def test_p3_modifier_matches_the_derived_magnitude(self):
        from biosim_ode.cortisol_pathway import CortisolParameters, CortisolPathway
        pathway = CortisolPathway(CortisolParameters.for_cushings())
        sig = pathway.output_signals(0.0, pathway.y0)
        assert sig["p3_modifier"] == pytest.approx(0.7784, abs=1e-3)

    def test_gb_is_not_touched_an_honest_partial_implementation(self):
        from biosim_ode.cortisol_pathway import CortisolParameters, CortisolPathway
        pathway = CortisolPathway(CortisolParameters.for_cushings())
        sig = pathway.output_signals(0.0, pathway.y0)
        assert sig["gb_modifier"] == pytest.approx(1.0)

    def test_composes_with_stress_level_still_a_fixed_point(self):
        from biosim_ode.cortisol_pathway import CortisolParameters, CortisolPathway
        pathway = CortisolPathway(CortisolParameters.for_cushings(stress_level=0.5))
        assert pathway.validate_y0(tolerance=1e-8) == []

    def test_composes_with_high_stress_level_without_a_bracket_failure(self):
        # code-reviewer finding: the inverse solver's bracket only widened
        # `hi` upward, never `lo` downward -- at high stress_level even a
        # small cortisol_basal input already overshoots the Cushing's
        # target, so f(lo) > 0 too and brentq raised a raw, unwrapped
        # ValueError. Confirmed to crash at stress_level>=0.90 before the
        # fix (now widens lo downward symmetrically). stress_level is
        # schema-valid up to 1.0 (ScenarioDefinition/REST both enforce
        # [0.0, 1.0]), so the full range must work, not just low stress.
        from biosim_ode.cortisol_pathway import CortisolParameters, CortisolPathway, CORTISOL_BASAL_NM
        for stress_level in [0.85, 0.90, 0.95, 1.0]:
            params = CortisolParameters.for_cushings(stress_level=stress_level)
            pathway = CortisolPathway(params)
            assert pathway.validate_y0(tolerance=1e-8) == []
            achieved_ratio = params.basal_state()[0] / CORTISOL_BASAL_NM
            assert achieved_ratio == pytest.approx(
                CortisolParameters.CUSHINGS_CORTISOL_RATIO, rel=1e-6,
            )

    def test_solve_cortisol_basal_for_target_steady_state_is_monotonic_and_accurate(self):
        from biosim_ode.cortisol_pathway import (
            solve_cortisol_basal_for_target_steady_state,
            solve_cortisol_crh_steady_state,
            CORTISOL_BASAL_NM,
        )
        for target_ratio in [1.1, 1.56, 2.0]:
            target = CORTISOL_BASAL_NM * target_ratio
            cortisol_basal = solve_cortisol_basal_for_target_steady_state(target, 5.0, 0.0)
            achieved, _, _ = solve_cortisol_crh_steady_state(cortisol_basal, 5.0, 0.0)
            assert achieved == pytest.approx(target, rel=1e-6)

    def test_runner_composes_this_end_to_end_through_the_real_pipeline(self):
        from biosim_api import ScenarioDefinition, MealEvent, SimulationRunner
        runner = SimulationRunner()
        healthy = runner.run(ScenarioDefinition(name="healthy", meal_events=[MealEvent(0, 60.0)]))
        cushings = runner.run(ScenarioDefinition(
            name="cushings", meal_events=[MealEvent(0, 60.0)], cushings_syndrome=True,
        ))
        assert healthy.success and cushings.success
        assert cushings.expert_metrics.glucose_peak_mM > healthy.expert_metrics.glucose_peak_mM


# ---------------------------------------------------------------------------
# Signal-collision accumulation
#
# Regression: CompositePathway.dydt() used to inject each coupling's
# transformed value via plain dict assignment (`injected[target][input] =
# value`), iterating self._couplings in list order -- so if two different
# couplings targeted the SAME (target_pathway, target_input), the later one
# in the list silently discarded the earlier one's contribution entirely,
# with no warning. Real-world collision: cortisol_couplings() and
# exercise_couplings() both target GlucoseInsulinPathway's gb_modifier and
# p3_modifier (see CLAUDE.md). Fixed to accumulate: multiplicative rule-11
# signals (the "_modifier" keys) multiply together; additive rule-11
# signals (hepatic_egp_delta, incretin_effect, renal_glucose_clearance_delta)
# sum. A target_input written by exactly one coupling behaves identically
# to before (single-writer keys like hepatic's "insulin_pM" or renal/GLP1's
# "glucose_mM" are unaffected).
# ---------------------------------------------------------------------------

class _DummySourcePathway(BasePathway):
    """Test double: stateless, exports one fixed signal value."""

    def __init__(self, name: str, signal_name: str, value: float) -> None:
        self._name = name
        self._signal_name = signal_name
        self._value = value

    @property
    def descriptor(self) -> PathwayDescriptor:
        return PathwayDescriptor(
            name=self._name,
            version="1.0.0",
            maturity=PathwayMaturity.EXPERIMENTAL,
            description="test double source pathway",
            species_names=[],
            output_signals=[self._signal_name],
            input_signals=[],
            dependencies=[],
            tags=[],
        )

    @property
    def species_names(self) -> list[str]:
        return []

    @property
    def species_units(self) -> list[str]:
        return []

    @property
    def y0(self) -> list[float]:
        return []

    def dydt(self, t, y, signals=None) -> list[float]:
        return []

    def output_signals(self, t, y) -> dict[str, float]:
        return {self._signal_name: self._value}


class _RecordingTargetPathway(BasePathway):
    """Test double: records whatever signals dict it receives from dydt()."""

    def __init__(self, name: str) -> None:
        self._name = name
        self.received_signals: Optional[dict] = None

    @property
    def descriptor(self) -> PathwayDescriptor:
        return PathwayDescriptor(
            name=self._name,
            version="1.0.0",
            maturity=PathwayMaturity.EXPERIMENTAL,
            description="test double recording target pathway",
            species_names=["dummy"],
            output_signals=[],
            # Wildcard ("*", matched via CouplingSpec.validate()'s
            # fnmatch support) rather than an empty list -- this test
            # double deliberately accepts ANY signal (its whole purpose
            # is recording whatever it's given to test
            # CompositePathway's own accumulation mechanism generically,
            # including intentionally-unrecognized keys like
            # "not_a_rule_11_signal"), not real pathway-contract
            # validation, which MASTER_PLAN 3.2 wired into
            # CompositePathway.__init__ for actual production pathways.
            input_signals=["*"],
            dependencies=[],
            tags=[],
        )

    @property
    def species_names(self) -> list[str]:
        return ["dummy"]

    @property
    def species_units(self) -> list[str]:
        return ["dimensionless"]

    @property
    def y0(self) -> list[float]:
        return [0.0]

    def dydt(self, t, y, signals=None) -> list[float]:
        self.received_signals = signals
        return [0.0]

    def output_signals(self, t, y) -> dict[str, float]:
        return {}


class TestCompositePathwaySignalAccumulation:
    def test_two_sources_targeting_same_modifier_key_multiply(self):
        target = _RecordingTargetPathway("target")
        src_a = _DummySourcePathway("src_a", "sig_a", 2.0)
        src_b = _DummySourcePathway("src_b", "sig_b", 3.0)
        couplings = [
            CouplingSpec("src_a", "sig_a", "target", "p3_modifier"),
            CouplingSpec("src_b", "sig_b", "target", "p3_modifier"),
        ]
        cp = CompositePathway([target, src_a, src_b], couplings=couplings)
        cp.dydt(0.0, cp.y0)
        assert target.received_signals == pytest.approx({"p3_modifier": 6.0})

    def test_two_sources_targeting_same_delta_key_sum(self):
        target = _RecordingTargetPathway("target")
        src_a = _DummySourcePathway("src_a", "sig_a", 0.01)
        src_b = _DummySourcePathway("src_b", "sig_b", 0.02)
        couplings = [
            CouplingSpec("src_a", "sig_a", "target", "hepatic_egp_delta"),
            CouplingSpec("src_b", "sig_b", "target", "hepatic_egp_delta"),
        ]
        cp = CompositePathway([target, src_a, src_b], couplings=couplings)
        cp.dydt(0.0, cp.y0)
        assert target.received_signals == pytest.approx({"hepatic_egp_delta": 0.03})

    def test_single_writer_key_unaffected(self):
        # Regression safety net: hepatic's "insulin_pM" / renal's and GLP1's
        # "glucose_mM" style single-writer target_inputs must behave exactly
        # as before -- no accumulation logic should perturb a lone value.
        target = _RecordingTargetPathway("target")
        src = _DummySourcePathway("src", "insulin_pM", 66.0)
        couplings = [CouplingSpec("src", "insulin_pM", "target", "insulin_pM")]
        cp = CompositePathway([target, src], couplings=couplings)
        cp.dydt(0.0, cp.y0)
        assert target.received_signals == pytest.approx({"insulin_pM": 66.0})

    def test_unrecognized_key_with_multiple_writers_raises(self):
        # A collision on a target_input that isn't one of rule 11's known
        # multiplicative/additive signals is genuinely ambiguous -- fail
        # loudly rather than silently guessing a combination rule.
        target = _RecordingTargetPathway("target")
        src_a = _DummySourcePathway("src_a", "sig_a", 1.0)
        src_b = _DummySourcePathway("src_b", "sig_b", 2.0)
        couplings = [
            CouplingSpec("src_a", "sig_a", "target", "not_a_rule_11_signal"),
            CouplingSpec("src_b", "sig_b", "target", "not_a_rule_11_signal"),
        ]
        cp = CompositePathway([target, src_a, src_b], couplings=couplings)
        with pytest.raises(ValueError, match="not_a_rule_11_signal"):
            cp.dydt(0.0, cp.y0)

    def test_cortisol_and_exercise_both_targeting_gi_accumulate(self):
        """
        Real-world reproduction: a stressed/steroid patient with a
        concurrent in-progress exercise bout composes CortisolPathway and
        ExercisePathway together, both targeting GI's gb_modifier and
        p3_modifier. Before the fix, exercise's neutral (1.0) values
        (exercise doesn't touch p1/p3/gb while the bout is still in
        progress) silently discarded cortisol's real, non-neutral values
        entirely -- cortisol's entire glycaemic effect vanished. After the
        fix, both contributions are multiplied together.
        """
        from biosim_ode.cortisol_pathway import (
            CortisolPathway, CortisolParameters, cortisol_couplings,
        )
        from biosim_ode.exercise_pathway import ExercisePathway, exercise_couplings
        from biosim_physiology.exercise import ExerciseBout

        gi = GlucoseInsulinPathway()
        cortisol = CortisolPathway(CortisolParameters.for_prednisone(dose_mg=20.0, stress_level=0.5))
        bout = ExerciseBout(intensity_percent_vo2max=70.0, duration_minutes=30.0)
        exercise = ExercisePathway(bout)

        couplings = (
            cortisol_couplings(gi.name, cortisol.name)
            + exercise_couplings(gi.name, exercise.name)
        )
        composite = CompositePathway([gi, cortisol, exercise], couplings=couplings)

        # Manually drive cortisol's own state away from its basal fixed
        # point (elevated Cortisol/CRH/GR_occ) -- deterministic, no ODE
        # integration needed -- to a state where its gb_modifier/p3_modifier
        # are genuinely non-neutral, matching a real elevated-cortisol state
        # during a concurrent exercise bout.
        cortisol_y = [1.2, 8.0, 0.85]
        y = list(gi.y0) + cortisol_y + []
        t = 10.0  # within the 30-minute bout: exercise's own p3/gb are 1.0 here

        cortisol_signals = cortisol.output_signals(t, cortisol_y)
        exercise_signals = exercise.output_signals(t, [])
        # Sanity-check the premise: cortisol has real, non-neutral signals,
        # and exercise's overlapping keys are exactly neutral at this t.
        assert cortisol_signals["gb_modifier"] != pytest.approx(1.0)
        assert cortisol_signals["p3_modifier"] != pytest.approx(1.0)
        assert exercise_signals["gb_modifier"] == pytest.approx(1.0)
        assert exercise_signals["p3_modifier"] == pytest.approx(1.0)

        derivs = composite.dydt(t, y)

        expected_signals = {
            "gb_modifier": cortisol_signals["gb_modifier"] * exercise_signals["gb_modifier"],
            "p3_modifier": cortisol_signals["p3_modifier"] * exercise_signals["p3_modifier"],
            "p1_modifier": exercise_signals["p1_modifier"],
            "exercise_glucose_uptake_delta": exercise_signals["exercise_glucose_uptake_delta"],
            "secretion_scale_modifier": exercise_signals["secretion_scale_modifier"],
        }
        expected_gi_derivs = gi.dydt(t, gi.y0, signals=expected_signals)
        n_gi = len(gi.y0)
        assert derivs[:n_gi] == pytest.approx(expected_gi_derivs, rel=1e-9)


class TestCompositePathwayCouplingValidation:
    """
    MASTER_PLAN 3.2: CouplingSpec.validate() existed (with its own
    regression fix checking target_input against target_inputs) but was
    never actually called anywhere in production -- a typo'd signal_name/
    target_input in any pathway-extension wiring silently composed into
    an inert, zero-effect coupling with no error at all, exactly the bug
    class in CLAUDE.md's KNOWN BUGS #7 (an entire session lost to a
    silently-inert coupling). Now wired into CompositePathway.__init__,
    failing fast at construction time.
    """

    def test_typo_d_target_input_raises_at_construction(self):
        target = _RecordingTargetPathway("target")
        src = _DummySourcePathway("src", "sig_a", 1.0)
        # "target" only declares input_signals=["*"] (accepts anything,
        # per the test double's own design) -- use a real production
        # pathway instead, whose descriptor genuinely restricts what it
        # accepts, so this test discriminates a real typo, not a
        # permissive test double.
        gi = GlucoseInsulinPathway()
        couplings = [CouplingSpec("src", "sig_a", gi.name, "p3_modifer")]  # typo: modifer
        with pytest.raises(ValueError, match="p3_modifer"):
            CompositePathway([gi, src], couplings=couplings)

    def test_typo_d_signal_name_raises_at_construction(self):
        gi = GlucoseInsulinPathway()
        src = _DummySourcePathway("src", "sig_a", 1.0)
        couplings = [CouplingSpec("src", "sig_typo", gi.name, "p3_modifier")]
        with pytest.raises(ValueError, match="sig_typo"):
            CompositePathway([gi, src], couplings=couplings)

    def test_unregistered_pathway_name_raises_at_construction(self):
        gi = GlucoseInsulinPathway()
        src = _DummySourcePathway("src", "sig_a", 1.0)
        couplings = [CouplingSpec("src", "sig_a", "nonexistent_pathway", "p3_modifier")]
        with pytest.raises(ValueError, match="nonexistent_pathway"):
            CompositePathway([gi, src], couplings=couplings)

    def test_meal_absorption_dynamic_signal_wildcard_does_not_raise(self):
        # The one real, deliberate gap the static input_signals metadata
        # can't enumerate exactly -- confirms GlucoseInsulinPathway's
        # "meal*_remaining_mM" wildcard entry and CouplingSpec.validate()'s
        # fnmatch support correctly accept this genuinely dynamic,
        # one-per-meal signal name rather than rejecting a real,
        # already-shipped coupling.
        from biosim_ode.glp1_pathway import GLP1Pathway
        from biosim_ode.meal_absorption_pathway import (
            MealAbsorptionPathway, meal_absorption_couplings,
        )
        from biosim_ode.meal import MealInput

        gi = GlucoseInsulinPathway()
        glp1 = GLP1Pathway.semaglutide(dose_mg=1.0)
        meals = [MealInput(carbs_g=75.0)]
        absorption = MealAbsorptionPathway(meals)
        couplings = meal_absorption_couplings(glp1.name, absorption.name, gi.name, meals)
        # Must not raise.
        CompositePathway([gi, glp1, absorption], couplings=couplings)

    def test_wildcard_target_input_match_is_case_sensitive(self):
        # Regression: found by the mandatory code-reviewer pass -- plain
        # fnmatch.fnmatch() normalizes both operands via os.path
        # .normcase() first, which lowercases on Windows but is a no-op
        # on POSIX, so a case-mismatched typo against a wildcard entry
        # (e.g. "Meal0_remaining_mM" instead of "meal0_remaining_mM")
        # would silently PASS validation on Windows while correctly
        # failing on Linux CI/production -- exactly the class of typo
        # this milestone exists to catch, made platform-dependent
        # instead of deterministic. Fixed with fnmatch.fnmatchcase()
        # (always case-sensitive regardless of OS).
        gi = GlucoseInsulinPathway()
        src = _DummySourcePathway("src", "sig_a", 1.0)
        couplings = [CouplingSpec("src", "sig_a", gi.name, "Meal0_remaining_mM")]
        with pytest.raises(ValueError, match="Meal0_remaining_mM"):
            CompositePathway([gi, src], couplings=couplings)

    def test_every_production_extension_pathway_composes_without_raising(self):
        # Systematic confirmation that this fix doesn't break any of the
        # real, already-shipped extension pathways -- each composed with
        # GI via its own *_couplings() function, the exact composition
        # SimulationRunner._build_pathway() performs in production.
        from biosim_ode.cortisol_pathway import CortisolPathway, cortisol_couplings
        from biosim_ode.hepatic import HepGlucPathway, hepatic_couplings
        from biosim_ode.nefa_pathway import NEFAPathway, nefa_couplings
        from biosim_ode.glp1_pathway import GLP1Pathway, glp1_couplings
        from biosim_ode.exercise_pathway import ExercisePathway, exercise_couplings
        from biosim_ode.renal_pathway import RenalClearanceParameters, RenalPathway, renal_couplings
        from biosim_ode.exogenous_insulin_pathway import (
            ExogenousInsulinParameters, ExogenousInsulinPathway, exogenous_insulin_couplings,
        )
        from biosim_physiology.exercise import ExerciseBout

        gi = GlucoseInsulinPathway()
        cortisol = CortisolPathway()
        CompositePathway([gi, cortisol], couplings=cortisol_couplings(gi.name, cortisol.name))

        gi = GlucoseInsulinPathway()
        hepatic = HepGlucPathway()
        CompositePathway([gi, hepatic], couplings=hepatic_couplings(gi.name, hepatic))

        gi = GlucoseInsulinPathway()
        nefa = NEFAPathway()
        CompositePathway([gi, nefa], couplings=nefa_couplings(gi.name, nefa.name))

        gi = GlucoseInsulinPathway()
        glp1 = GLP1Pathway.semaglutide(dose_mg=1.0)
        CompositePathway([gi, glp1], couplings=glp1_couplings(gi.name, glp1.name))

        gi = GlucoseInsulinPathway()
        exercise = ExercisePathway(ExerciseBout(intensity_percent_vo2max=70.0, duration_minutes=30.0))
        CompositePathway([gi, exercise], couplings=exercise_couplings(gi.name, exercise.name))

        gi = GlucoseInsulinPathway()
        renal = RenalPathway(RenalClearanceParameters.for_dapagliflozin())
        CompositePathway([gi, renal], couplings=renal_couplings(gi.name, renal))

        gi = GlucoseInsulinPathway()
        insulin = ExogenousInsulinPathway(ExogenousInsulinParameters.regular_insulin(dose_units=10.0))
        CompositePathway(
            [gi, insulin], couplings=exogenous_insulin_couplings(gi.name, insulin.name)
        )


# ---------------------------------------------------------------------------
# NEFA pathway recalibration (2026-08-12 deep-search plan, items C2/C3)
# ---------------------------------------------------------------------------

class TestNEFACitationRecalibration:
    """
    NEFA_HALF_LIFE_MIN (90.0 -> 3.7, Heiling et al. 1991 PMID:1951681,
    corroborated by Hagenfeldt et al. 1975 PMID:1213651) and
    NEFA_INSULIN_IC50_PM (80.0 -> 68.0, Stumvoll et al. 2000
    PMID:11061533), both live-verified twice this session. See
    nefa_pathway.py's own module-level constant comments for the full
    derivation.
    """

    def test_half_life_matches_heiling_1991(self):
        from biosim_ode.nefa_pathway import NEFA_HALF_LIFE_MIN
        assert NEFA_HALF_LIFE_MIN == pytest.approx(3.7)

    def test_insulin_ic50_matches_stumvoll_2000_adipose_value(self):
        from biosim_ode.nefa_pathway import NEFA_INSULIN_IC50_PM
        assert NEFA_INSULIN_IC50_PM == pytest.approx(68.0)

    def test_basal_fixed_point_independent_of_half_life(self):
        # The target-based dNEFA formula's fixed point does not depend
        # on the rate constant governing approach to it -- confirmed for
        # both the new (fast) and a hypothetical slow half-life.
        from biosim_ode.nefa_pathway import NEFAPathway, NEFAParameters
        for half_life in [3.7, 90.0, 1.5]:
            p = NEFAPathway(NEFAParameters(half_life_min=half_life))
            assert p.validate_y0(tolerance=1e-8) == []

    def test_composed_basal_fixed_point_stable_over_long_fasted_trajectory(self):
        from biosim_ode.nefa_pathway import NEFAPathway, couple_with_gi
        composite = couple_with_gi(GlucoseInsulinPathway(), NEFAPathway())
        assert composite.validate_y0(tolerance=1e-8) == []
        result = ODESystem(pathway=composite).simulate(duration_s=600 * 60, meals=[])
        nefa = result.species["NEFA"]
        assert nefa.max() == pytest.approx(nefa.min(), abs=1e-9)
        assert nefa[0] == pytest.approx(350.0)

    def test_new_half_life_makes_nefa_respond_much_faster_to_a_meal(self):
        # Regression guard against silently reverting to the old,
        # uncited 90-minute half-life -- with the new, much faster
        # clearance, a real postprandial insulin rise should measurably
        # suppress NEFA within the first ~30 minutes (impossible with a
        # 90-minute half-life, whose own time constant alone exceeds
        # that window).
        from biosim_ode.nefa_pathway import NEFAPathway, couple_with_gi
        meal = high_sugar_meal()
        composite = couple_with_gi(
            GlucoseInsulinPathway(meals=[meal]), NEFAPathway(),
        )
        result = ODESystem(pathway=composite).simulate(duration_s=300 * 60, meals=[meal])
        t_min = result.times_s / 60.0
        nefa = result.species["NEFA"]
        idx_30min = (abs(t_min - 30.0)).argmin()
        assert nefa[idx_30min] < 0.5 * nefa[0]


class TestNEFAInsulinResistanceEC50Recalibration:
    """
    NEFA_IR_EC50_UM (500.0 -> 200.0, MASTER_PLAN row 1.22, a fresh-angle
    retry). Anchored to Belfort et al. 2005 (PMID:15919784) -- the first
    citation found reporting BOTH a graded FFA dose-response AND its own
    basal/reference FFA level in the same cohort, resolving the excess-
    above-basal ambiguity that blocked the original Gormsen/Boden search
    pass. NEFA_IR_MAX left unchanged (0.40) -- the new fit's own range
    (0.360-0.452) comfortably includes it. See nefa_pathway.py's own
    module-level constant comment for the full derivation.
    """

    def test_ec50_matches_belfort_2005_fit(self):
        from biosim_ode.nefa_pathway import NEFA_IR_EC50_UM
        assert NEFA_IR_EC50_UM == pytest.approx(200.0)

    def test_max_unchanged_and_within_the_new_fits_own_range(self):
        # The fit's own pairwise-exact-solve range (0.360-0.452) should
        # still bracket the unchanged 0.40 -- a regression guard against
        # silently drifting the two apart without updating this comment.
        from biosim_ode.nefa_pathway import NEFA_IR_MAX
        assert NEFA_IR_MAX == pytest.approx(0.40)
        assert 0.360 <= NEFA_IR_MAX <= 0.452

    def test_belfort_2005_own_three_points_are_correctly_fit(self):
        # Independently re-derive the fit directly from the paper's own
        # quoted numbers (basal ~440; doses 695/1251/1688; reductions
        # 22/30/34%) using this pathway's own exact Michaelis-Menten
        # formula, and confirm the adopted constants reproduce those
        # three real data points within a small, disclosed tolerance --
        # not just "some fit was done," but that THIS fit matches the
        # real, quoted numbers reasonably well.
        from biosim_ode.nefa_pathway import NEFA_IR_EC50_UM, NEFA_IR_MAX
        basal = 440.0
        doses = [695.0, 1251.0, 1688.0]
        reported_reductions = [0.22, 0.30, 0.34]
        for dose, reported in zip(doses, reported_reductions):
            excess = dose - basal
            predicted = NEFA_IR_MAX * excess / (NEFA_IR_EC50_UM + excess)
            assert predicted == pytest.approx(reported, abs=0.03)

    def test_basal_fixed_point_unaffected_by_the_ec50_change(self):
        # insulin_resistance_modifier is a pure output signal (rule 11)
        # -- confirms changing its dose-response shape doesn't perturb
        # NEFAPathway's own basal fixed point.
        from biosim_ode.nefa_pathway import NEFAPathway, couple_with_gi
        composite = couple_with_gi(GlucoseInsulinPathway(), NEFAPathway())
        assert composite.validate_y0(tolerance=1e-8) == []

    def test_ec50_change_genuinely_shifts_the_modifier_at_a_realistic_excess(self):
        # Regression guard against silently reverting to the old,
        # uncited 500.0 -- at a realistic elevated-NEFA excess (e.g. the
        # obesity preset), the new, much lower EC50 should produce a
        # measurably stronger (lower) insulin_resistance_modifier than
        # the old value would have.
        from biosim_ode.nefa_pathway import NEFAPathway, NEFAParameters
        pathway_new = NEFAPathway(NEFAParameters(ir_ec50=200.0, ir_max=0.40))
        pathway_old = NEFAPathway(NEFAParameters(ir_ec50=500.0, ir_max=0.40))
        # A realistic elevated NEFA level (e.g. obesity preset ~700 uM)
        signals_new = pathway_new.output_signals(0.0, [700.0])
        signals_old = pathway_old.output_signals(0.0, [700.0])
        assert signals_new["insulin_resistance_modifier"] < signals_old["insulin_resistance_modifier"]


class TestNEFABasalObesityDisclosure:
    """
    Fresh Track 1 re-scoping (2026-08-18), a second independent search
    pass for MASTER_PLAN row 1.19 (NEFA_BASAL_UM/NEFA_OBESITY_UM). Both
    this pass (Kelsey et al. 2014, PMID:23861170) and the original pass
    (Curry et al. 2011, PMID:21304064) found real evidence CONTRADICTING
    this preset's own implied ~2x obesity-elevation premise. Deliberately
    NOT recalibrated to either study's own numbers (population mismatch,
    high variance, non-significant group difference) -- documentation-
    only, this row is now formally parked. See nefa_pathway.py's own
    module-level comment for the full derivation.
    """

    def test_values_are_genuinely_unchanged_by_this_pass(self):
        from biosim_ode.nefa_pathway import NEFA_BASAL_UM, NEFA_OBESITY_UM
        assert NEFA_BASAL_UM == pytest.approx(350.0)
        assert NEFA_OBESITY_UM == pytest.approx(700.0)

    def test_source_discloses_both_counter_evidence_studies(self):
        import biosim_ode.nefa_pathway as nefa_module
        import inspect
        source = inspect.getsource(nefa_module)
        assert "21304064" in source   # Curry 2011
        assert "23861170" in source   # Kelsey 2014
        assert "formally PARKED" in source

    def test_source_discloses_the_real_kelsey_2014_table_numbers(self):
        import biosim_ode.nefa_pathway as nefa_module
        import inspect
        source = inspect.getsource(nefa_module)
        assert "521" in source
        assert "583" in source
        assert "0.86" in source


# ---------------------------------------------------------------------------
# ODESystem generalised to run any BasePathway (composite included)
# ---------------------------------------------------------------------------

class TestODESystemComposite:
    def test_composite_species_use_bare_names(self):
        from biosim_ode.cortisol_pathway import CortisolPathway
        composite = CompositePathway([GlucoseInsulinPathway(), CortisolPathway()])
        result = ODESystem(pathway=composite).simulate_meal(high_sugar_meal(), duration_s=3600)
        assert result.success
        assert "G" in result.species and "Cortisol" in result.species
        # No leftover "glucose_insulin_v1.G"-style prefixed keys
        assert not any("." in k for k in result.species)

    def test_composite_units_converted_correctly(self):
        from biosim_ode.cortisol_pathway import CortisolPathway
        composite = CompositePathway([GlucoseInsulinPathway(), CortisolPathway()])
        result = ODESystem(pathway=composite).simulate_fasted(duration_s=600)
        assert result.species_units["G"] == "uM"
        assert result.species_units["Cortisol"] == "uM"   # nM -> uM canonical
        assert result.species["Cortisol"][0] == pytest.approx(0.28 * 1e-3, rel=1e-6)

    def test_composite_watchdog_still_runs(self):
        from biosim_ode.cortisol_pathway import CortisolPathway
        composite = CompositePathway([GlucoseInsulinPathway(), CortisolPathway()])
        result = ODESystem(pathway=composite).simulate_meal(high_sugar_meal(), duration_s=3600)
        assert result.watchdog and len(result.watchdog) > 0

    def test_composite_without_gi_has_no_watchdog(self):
        from biosim_ode.cortisol_pathway import CortisolPathway
        from biosim_ode.nefa_pathway import NEFAPathway
        composite = CompositePathway([CortisolPathway(), NEFAPathway()])
        result = ODESystem(pathway=composite).simulate(duration_s=600, meals=[])
        assert result.success
        assert not result.watchdog

    def test_bare_gi_unaffected_by_generalisation(self):
        # Direct proof (not just "other tests still pass") that the
        # generalised species-name/unit lookup reproduces the exact same
        # values as before for the everyday bare-GI case.
        result = ODESystem().simulate_meal(high_sugar_meal(), duration_s=7200)
        assert result.species_units["G"] == "uM"
        assert result.species["G"][0] == pytest.approx(5000.0, rel=1e-6)


# ---------------------------------------------------------------------------
# Regression: existing ODE results unchanged after refactor
# ---------------------------------------------------------------------------

class TestRegressionAfterRefactor:
    """
    Confirm that wrapping GlucoseInsulinPathway as BasePathway
    does not change any simulation outputs.
    """

    def test_fasted_glucose_unchanged(self):
        system = ODESystem()
        result = system.simulate_fasted(duration_s=3600)
        G = result.species["G"]
        G_baseline = 5.0 * 1000.0
        assert np.max(np.abs(G - G_baseline)) < 50.0

    def test_meal_glucose_peak_unchanged(self):
        system = ODESystem()
        result = system.simulate_meal(high_sugar_meal(), duration_s=7200)
        peak_G, _ = result.glucose_peak()
        assert peak_G > 5500

    def test_validate_y0_still_passes(self):
        p = GlucoseInsulinPathway()
        assert p.validate_y0() == []


# ---------------------------------------------------------------------------
# Disclosure of the point-of-use clamps (P3_EFF_FLOOR et al.). personalizer.py's
# static clamp_warnings only sees PathwayParameters before integration; the
# dydt() clamps above engage on ODE-COUPLED signals with no disclosure at all.
# point_of_use_clamp_hits() recomputes the raw composed values, and
# ODESystem.simulate() sweeps it across the output trajectory.
# ---------------------------------------------------------------------------

class TestPointOfUseClampDisclosure:
    def test_no_hits_at_neutral_signals(self):
        assert GlucoseInsulinPathway().point_of_use_clamp_hits({}) == {}
        assert GlucoseInsulinPathway().point_of_use_clamp_hits(None) == {}

    def test_p3_floor_hit_reports_side_and_raw_value(self):
        p = GlucoseInsulinPathway()
        hits = p.point_of_use_clamp_hits({"p3_modifier": 0.005})
        assert set(hits) == {"p3"}
        side, raw = hits["p3"]
        assert side == "floor"
        assert raw == pytest.approx(p.params.p3 * 0.005)

    def test_insulin_resistance_modifier_composes_into_p3(self):
        # p3_eff = p3 * p3_modifier * insulin_resistance_modifier -- neither
        # alone breaches the floor here, their product does.
        p = GlucoseInsulinPathway()
        assert p.point_of_use_clamp_hits({"p3_modifier": 0.1}) == {}
        hits = p.point_of_use_clamp_hits(
            {"p3_modifier": 0.1, "insulin_resistance_modifier": 0.05}
        )
        assert hits["p3"][0] == "floor"

    def test_ceiling_hits_for_each_other_clamped_parameter(self):
        p = GlucoseInsulinPathway(params=PathwayParameters(extra_g_clearance=0.018))
        hits = p.point_of_use_clamp_hits({
            "gb_modifier": 100.0,
            "p1_modifier": 1000.0,
            "extra_g_clearance_modifier": 100.0,
            "secretion_scale_modifier": 100.0,
        })
        assert {k: v[0] for k, v in hits.items()} == {
            "Gb": "ceiling", "p1": "ceiling",
            "extra_g_clearance": "ceiling", "secretion_scale": "ceiling",
        }

    def test_nan_signal_is_reported_as_a_floor_hit(self):
        # clamp() maps NaN to the floor, so dydt() really does run at the
        # floor -- the disclosure must not skip it as "not below".
        hits = GlucoseInsulinPathway().point_of_use_clamp_hits({"p3_modifier": float("nan")})
        side, raw = hits["p3"]
        assert side == "floor" and np.isnan(raw)

    def test_sweep_keeps_the_most_extreme_raw_value(self):
        from biosim_ode.system import _more_extreme
        worst = None
        for raw in (5e-9, 2e-9, 7e-9):
            worst = _more_extreme("floor", worst, raw)
        assert worst == 2e-9
        worst = None
        for raw in (20.0, 40.0, 30.0):
            worst = _more_extreme("ceiling", worst, raw)
        assert worst == 40.0

    def test_floor_and_ceiling_hits_on_one_field_are_both_reported(self):
        # Code-review finding: a field pushed to its floor early and its
        # ceiling later (e.g. a time-gated coupling) must report both sides.
        from biosim_ode.system import _sweep_clamp_hits

        class _FlippingSource(_DummySourcePathway):
            def output_signals(self, t, y):
                return {self._signal_name: 0.005 if t < 30.0 else 500.0}

        gi = GlucoseInsulinPathway()
        cp = CompositePathway(
            [gi, _FlippingSource("src", "p3_mod", 0.0)],
            couplings=[CouplingSpec("src", "p3_mod", gi.name, "p3_modifier")],
        )
        t = np.array([0.0, 10.0, 40.0, 60.0])
        y = np.tile(np.array(cp.y0, dtype=float)[:, None], (1, len(t)))
        hits = _sweep_clamp_hits(cp, t, y)
        assert set(hits) == {("p3", "floor"), ("p3", "ceiling")}
        assert hits[("p3", "floor")] == pytest.approx(gi.params.p3 * 0.005)
        assert hits[("p3", "ceiling")] == pytest.approx(gi.params.p3 * 500.0)

    def test_simulation_reports_a_hit_driven_by_a_coupled_signal(self):
        gi = GlucoseInsulinPathway()
        src = _DummySourcePathway("src", "p3_mod", 0.005)
        cp = CompositePathway(
            [gi, src], couplings=[CouplingSpec("src", "p3_mod", gi.name, "p3_modifier")],
        )
        result = ODESystem(pathway=cp).simulate(duration_s=3600, n_output_points=20)
        assert result.success
        assert set(result.clamp_hits) == {("p3", "floor")}
        assert result.clamp_hits[("p3", "floor")] == pytest.approx(gi.params.p3 * 0.005)

    def test_check_clamps_false_skips_the_sweep(self):
        # Ensemble callers (UQ/OAT/Sobol) discard the result and opt out.
        gi = GlucoseInsulinPathway()
        cp = CompositePathway(
            [gi, _DummySourcePathway("src", "p3_mod", 0.005)],
            couplings=[CouplingSpec("src", "p3_mod", gi.name, "p3_modifier")],
        )
        result = ODESystem(pathway=cp).simulate(
            duration_s=3600, n_output_points=20, check_clamps=False,
        )
        assert result.success and result.clamp_hits == {}

    def test_simulation_has_no_hits_by_default(self):
        result = ODESystem(pathway=GlucoseInsulinPathway()).simulate_meal(
            high_sugar_meal(), duration_s=3600, n_output_points=20,
        )
        assert result.success
        assert result.clamp_hits == {}
