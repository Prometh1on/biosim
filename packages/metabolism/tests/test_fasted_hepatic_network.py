# SAVE TO: packages\metabolism\tests\test_fasted_hepatic_network.py
"""
Tests for biosim_metabolism.fasted_hepatic_network (Track 4, MASTER_PLAN
4.3b).

These tests are the "permanent structural checks" a hand-curated
stoichiometric network needs (per this chunk's own scoping pass):
elemental/moiety balance, closed-exchange loop-freedom, and canonical
precursor-to-glucose yield reproduction. Every numeric result below was
independently re-derived and verified by direct execution before being
written down here, not copied from the scoping pass's own prototype.
"""
import numpy as np
import pytest

from biosim_metabolism import (
    InfeasibleModelError,
    StoichiometricModel,
    flux_variability_analysis,
)
from biosim_metabolism.fasted_hepatic_network import (
    METABOLITES,
    _INTERNAL_REACTIONS,
    build_fasted_hepatic_network,
)

# Real conserved quantities tracked across this network (H and O are NOT
# tracked, per the module's own documented scope). "adenylate" is the
# adenine+ribose skeleton (ATP/ADP/AMP each = 1 unit); "nadpool"/
# "fadpool" the nicotinamide/flavin skeleton (NAD+/NADH, FAD/FADH2 each
# = 1 unit); "coapool" the CoA (pantetheine-ADP) moiety, conserved
# across free CoA AND acetyl-CoA (which still carries the full CoA
# backbone, plus a separately-tracked 2-carbon acetyl group).
_ELEMENTS = {
    "glycogen": {"C": 6, "N": 0, "P": 0},
    "g1p":      {"C": 6, "N": 0, "P": 1},
    "g6p":      {"C": 6, "N": 0, "P": 1},
    "f6p":      {"C": 6, "N": 0, "P": 1},
    "fdp":      {"C": 6, "N": 0, "P": 2},
    "t3p":      {"C": 3, "N": 0, "P": 1},
    "pep":      {"C": 3, "N": 0, "P": 1},
    "oaa":      {"C": 4, "N": 0, "P": 0},
    "pyr":      {"C": 3, "N": 0, "P": 0},
    "lac":      {"C": 3, "N": 0, "P": 0},
    "accoa":    {"C": 2, "N": 0, "P": 0, "coapool": 1},
    "akg":      {"C": 5, "N": 0, "P": 0},
    "glu":      {"C": 5, "N": 1, "P": 0},
    "ala":      {"C": 3, "N": 1, "P": 0},
    "nh4":      {"C": 0, "N": 1, "P": 0},
    "glyc":     {"C": 3, "N": 0, "P": 0},
    "glyc3p":   {"C": 3, "N": 0, "P": 1},
    "fa":       {"C": 16, "N": 0, "P": 0},
    "co2":      {"C": 1, "N": 0, "P": 0},
    "o2":       {"C": 0, "N": 0, "P": 0},
    "glc":      {"C": 6, "N": 0, "P": 0},
    "atp":      {"C": 0, "N": 0, "P": 3, "adenylate": 1},
    "adp":      {"C": 0, "N": 0, "P": 2, "adenylate": 1},
    "amp":      {"C": 0, "N": 0, "P": 1, "adenylate": 1},
    "pi":       {"C": 0, "N": 0, "P": 1},
    "nad":      {"nadpool": 1},
    "nadh":     {"nadpool": 1},
    "fad":      {"fadpool": 1},
    "fadh2":    {"fadpool": 1},
    "coa":      {"coapool": 1},
}

# The only two reactions permitted to be unbalanced -- deliberate,
# documented untracked-export "demand" reactions (see the module's own
# docstring).
_ALLOWED_UNBALANCED = {"UREA_SINK", "KETO_SINK"}


def _reaction_balance(stoich: dict) -> dict:
    balance = {}
    for met, coeff in stoich.items():
        for element, count in _ELEMENTS.get(met, {}).items():
            balance[element] = balance.get(element, 0) + coeff * count
    return {k: v for k, v in balance.items() if abs(v) > 1e-9}


class TestNetworkStructure:
    def test_network_has_the_expected_shape(self):
        model = build_fasted_hepatic_network()
        assert model.n_metabolites == 30
        assert model.n_reactions == 32

    def test_elements_table_covers_every_metabolite(self):
        # Regression guard, found by code review: _ELEMENTS's
        # get(met, {}) fallback silently treats an uncovered metabolite
        # as contributing zero to every element -- if a future
        # metabolite is added to the network without a matching
        # _ELEMENTS entry, the balance check below would never catch a
        # real imbalance involving it. This test makes that impossible.
        assert set(_ELEMENTS.keys()) == set(METABOLITES)

    def test_every_internal_reaction_is_elementally_balanced_except_the_two_flagged_sinks(self):
        for rxn, stoich in _INTERNAL_REACTIONS.items():
            imbalance = _reaction_balance(stoich)
            if rxn in _ALLOWED_UNBALANCED:
                assert imbalance, f"{rxn} was expected to be a deliberate untracked sink, but balances exactly"
            else:
                assert not imbalance, f"{rxn} is not elementally balanced: {imbalance}"

    def test_urea_sink_imbalance_is_exactly_one_urea_worth_of_carbon_and_nitrogen(self):
        imbalance = _reaction_balance(_INTERNAL_REACTIONS["UREA_SINK"])
        assert imbalance == {"C": -1, "N": -2}

    def test_keto_sink_imbalance_is_exactly_one_acetoacetate_worth_of_carbon(self):
        # 2 acetyl-CoA (2x2=4 tracked acetyl carbons) combine into
        # exactly 1 acetoacetate via the HMG-CoA pathway -- not "2x
        # acetoacetate" (a wording error found and fixed by code review;
        # the -4 C numeric assertion itself was always correct).
        imbalance = _reaction_balance(_INTERNAL_REACTIONS["KETO_SINK"])
        assert imbalance == {"C": -4}

    def test_closing_every_exchange_forces_every_internal_reaction_to_zero_flux(self):
        # No hidden, stoichiometrically-balanced internal loop that
        # could inflate an FVA range on its own -- the single most
        # important structural guarantee for a hand-curated network.
        # Regression, found by code review: EX_co2's bound used to be
        # hardcoded (not a real parameter), so this test could not
        # actually close it -- ex_co2_bound=0.0 now genuinely closes it
        # to exactly [0, 0], so EVERY exchange is closed here, not
        # all-but-one.
        model = build_fasted_hepatic_network(
            ex_glycogen_bound=0.0, ex_lac_bound=0.0, ex_glyc_bound=0.0,
            ex_ala_bound=0.0, ex_fa_bound=0.0, ex_o2_bound=0.0,
            ex_co2_bound=0.0, egp_upper_bound=0.0,
        )
        results = flux_variability_analysis(model)
        for r in results:
            assert r.min_flux == pytest.approx(0.0, abs=1e-9), r.reaction
            assert r.max_flux == pytest.approx(0.0, abs=1e-9), r.reaction

    def test_adk_is_structurally_required_for_fao_and_urea_sink_to_carry_any_flux(self):
        # FAO and UREA_SINK are the only two AMP-producing reactions;
        # ADK is the only AMP-consuming one. Removing ADK must force
        # both to exactly zero flux under steady state.
        #
        # Regression, found by code review: an earlier version of this
        # test manually rebuilt the model from _INTERNAL_REACTIONS,
        # duplicating (and getting wrong) the real builder's
        # reversibility bounds -- a real drift risk, since a future
        # change to build_fasted_hepatic_network()'s own bounds logic
        # wouldn't be reflected here. Fixed by building from the REAL
        # builder's own output and stripping just the ADK column,
        # so this test can never silently diverge from the real network.
        base_model = build_fasted_hepatic_network()
        adk_idx = base_model.reactions.index("ADK")
        reaction_names = [r for i, r in enumerate(base_model.reactions) if i != adk_idx]
        S = np.delete(np.asarray(base_model.stoichiometric_matrix), adk_idx, axis=1)
        lb = np.delete(np.asarray(base_model.lower_bounds), adk_idx)
        ub = np.delete(np.asarray(base_model.upper_bounds), adk_idx)
        model = StoichiometricModel(
            metabolites=list(base_model.metabolites), reactions=reaction_names,
            stoichiometric_matrix=S, lower_bounds=lb, upper_bounds=ub,
        )
        results = {r.reaction: r for r in flux_variability_analysis(model)}
        assert results["FAO"].max_flux == pytest.approx(0.0, abs=1e-9)
        assert results["UREA_SINK"].max_flux == pytest.approx(0.0, abs=1e-9)


class TestCanonicalYields:
    """
    Each real precursor, in isolation, with unlimited energy (fat/O2
    open), must reproduce its textbook glucose yield ratio exactly:
    1 glucose per 1 glycogen unit / 2 lactate / 2 glycerol / 2 alanine.
    """

    @staticmethod
    def _max_egp(**kwargs):
        model = build_fasted_hepatic_network(ex_fa_bound=1000.0, ex_o2_bound=1000.0, **kwargs)
        results = {r.reaction: r for r in flux_variability_analysis(model)}
        return results["G6PASE"].max_flux

    def test_glycogen_yield(self):
        assert self._max_egp(
            ex_glycogen_bound=1.0, ex_lac_bound=0, ex_glyc_bound=0, ex_ala_bound=0,
        ) == pytest.approx(1.0, abs=1e-6)

    def test_lactate_yield(self):
        assert self._max_egp(
            ex_glycogen_bound=0, ex_lac_bound=2.0, ex_glyc_bound=0, ex_ala_bound=0,
        ) == pytest.approx(1.0, abs=1e-6)

    def test_glycerol_yield(self):
        assert self._max_egp(
            ex_glycogen_bound=0, ex_lac_bound=0, ex_glyc_bound=2.0, ex_ala_bound=0,
        ) == pytest.approx(1.0, abs=1e-6)

    def test_alanine_yield(self):
        assert self._max_egp(
            ex_glycogen_bound=0, ex_lac_bound=0, ex_glyc_bound=0, ex_ala_bound=2.0,
        ) == pytest.approx(1.0, abs=1e-6)


class TestEnergyDependence:
    """
    Regression-locks the scoping pass's own key finding: converting
    lactate or alanine to glucose at the exact textbook 2:1 stoichiometric
    ratio (no surplus precursor left over to fund the reaction's own ATP
    cost via oxidation) is INFEASIBLE without fat availability -- glycerol
    is the one precursor that funds its own conversion (net NADH surplus),
    matching real biochemistry. This is the structural reason the whole
    network's FVA answer is currently energy-bound-driven, not purely
    stoichiometry-driven, once fat/O2/ATPM are meaningfully constrained.
    """

    @staticmethod
    def _fixed_ratio_feasible(precursor_ex: str, precursor_amount: float, glc_amount: float, fa_bound: float) -> bool:
        model = build_fasted_hepatic_network(
            ex_glycogen_bound=0, ex_lac_bound=0, ex_glyc_bound=0, ex_ala_bound=0,
            ex_fa_bound=fa_bound, ex_o2_bound=1000.0, egp_upper_bound=1000.0,
        )
        reactions = model.reactions
        lb = np.array(model.lower_bounds, dtype=float)
        ub = np.array(model.upper_bounds, dtype=float)
        i_p = reactions.index(precursor_ex)
        i_glc = reactions.index("EX_glc")
        lb[i_p] = ub[i_p] = precursor_amount
        lb[i_glc] = ub[i_glc] = glc_amount
        fixed_model = StoichiometricModel(
            metabolites=model.metabolites, reactions=reactions,
            stoichiometric_matrix=model.stoichiometric_matrix,
            lower_bounds=lb, upper_bounds=ub,
        )
        try:
            flux_variability_analysis(fixed_model)
            return True
        except InfeasibleModelError:
            return False

    def test_lactate_at_exact_ratio_is_infeasible_without_fat(self):
        assert not self._fixed_ratio_feasible("EX_lac", 2.0, 1.0, fa_bound=0.0)

    def test_lactate_at_exact_ratio_is_feasible_with_fat(self):
        assert self._fixed_ratio_feasible("EX_lac", 2.0, 1.0, fa_bound=1000.0)

    def test_alanine_at_exact_ratio_is_infeasible_without_fat(self):
        assert not self._fixed_ratio_feasible("EX_ala", 2.0, 1.0, fa_bound=0.0)

    def test_alanine_at_exact_ratio_is_feasible_with_fat(self):
        assert self._fixed_ratio_feasible("EX_ala", 2.0, 1.0, fa_bound=1000.0)

    def test_glycerol_at_exact_ratio_is_feasible_without_fat(self):
        # Glycerol's own dehydrogenation step (GPD1) produces a net NADH
        # surplus, funding its own conversion -- unlike lactate/alanine.
        assert self._fixed_ratio_feasible("EX_glyc", 2.0, 1.0, fa_bound=0.0)


class TestBuildFastedHepaticNetworkDefaults:
    def test_all_bounds_default_to_permissive_1000(self):
        model = build_fasted_hepatic_network()
        for name in ["EX_glycogen", "EX_lac", "EX_glyc", "EX_ala", "EX_fa", "EX_o2", "EX_glc"]:
            idx = model.reactions.index(name)
            assert model.upper_bounds[idx] == pytest.approx(1000.0)

    def test_explicit_bound_overrides_are_respected(self):
        model = build_fasted_hepatic_network(ex_lac_bound=3.5)
        idx = model.reactions.index("EX_lac")
        assert model.upper_bounds[idx] == pytest.approx(3.5)

    def test_ex_co2_is_reversible(self):
        model = build_fasted_hepatic_network()
        idx = model.reactions.index("EX_co2")
        assert model.lower_bounds[idx] < 0.0
