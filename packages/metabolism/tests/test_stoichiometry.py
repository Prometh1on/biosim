# SAVE TO: packages\metabolism\tests\test_stoichiometry.py
"""Tests for biosim_metabolism.stoichiometry (Track 4, MASTER_PLAN 4.3a)."""
import numpy as np
import pytest

from biosim_metabolism import (
    FVAResult,
    InfeasibleModelError,
    StoichiometricModel,
    flux_variability_analysis,
)


class TestStoichiometricModelValidation:
    def test_valid_construction(self):
        m = StoichiometricModel(
            metabolites=["M1"],
            reactions=["v1", "v2"],
            stoichiometric_matrix=np.array([[1.0, -1.0]]),
            lower_bounds=np.array([0.0, 0.0]),
            upper_bounds=np.array([10.0, 10.0]),
        )
        assert m.n_metabolites == 1
        assert m.n_reactions == 2

    def test_empty_metabolites_raises(self):
        with pytest.raises(ValueError, match="metabolite"):
            StoichiometricModel(
                metabolites=[], reactions=["v1"],
                stoichiometric_matrix=np.zeros((0, 1)),
                lower_bounds=np.array([0.0]), upper_bounds=np.array([1.0]),
            )

    def test_empty_reactions_raises(self):
        with pytest.raises(ValueError, match="reaction"):
            StoichiometricModel(
                metabolites=["M1"], reactions=[],
                stoichiometric_matrix=np.zeros((1, 0)),
                lower_bounds=np.array([]), upper_bounds=np.array([]),
            )

    def test_duplicate_metabolite_ids_raises(self):
        with pytest.raises(ValueError, match="unique"):
            StoichiometricModel(
                metabolites=["M1", "M1"], reactions=["v1"],
                stoichiometric_matrix=np.array([[1.0], [1.0]]),
                lower_bounds=np.array([0.0]), upper_bounds=np.array([1.0]),
            )

    def test_duplicate_reaction_ids_raises(self):
        with pytest.raises(ValueError, match="unique"):
            StoichiometricModel(
                metabolites=["M1"], reactions=["v1", "v1"],
                stoichiometric_matrix=np.array([[1.0, -1.0]]),
                lower_bounds=np.array([0.0, 0.0]), upper_bounds=np.array([1.0, 1.0]),
            )

    def test_wrong_matrix_shape_raises(self):
        with pytest.raises(ValueError, match="shape"):
            StoichiometricModel(
                metabolites=["M1", "M2"], reactions=["v1"],
                stoichiometric_matrix=np.array([[1.0]]),   # shape (1,1), expected (2,1)
                lower_bounds=np.array([0.0]), upper_bounds=np.array([1.0]),
            )

    def test_wrong_lower_bounds_shape_raises(self):
        with pytest.raises(ValueError, match="lower_bounds"):
            StoichiometricModel(
                metabolites=["M1"], reactions=["v1", "v2"],
                stoichiometric_matrix=np.array([[1.0, -1.0]]),
                lower_bounds=np.array([0.0]), upper_bounds=np.array([1.0, 1.0]),
            )

    def test_wrong_upper_bounds_shape_raises(self):
        with pytest.raises(ValueError, match="upper_bounds"):
            StoichiometricModel(
                metabolites=["M1"], reactions=["v1", "v2"],
                stoichiometric_matrix=np.array([[1.0, -1.0]]),
                lower_bounds=np.array([0.0, 0.0]), upper_bounds=np.array([1.0]),
            )

    def test_nan_in_matrix_raises(self):
        with pytest.raises(ValueError, match="NaN/Inf"):
            StoichiometricModel(
                metabolites=["M1"], reactions=["v1"],
                stoichiometric_matrix=np.array([[float("nan")]]),
                lower_bounds=np.array([0.0]), upper_bounds=np.array([1.0]),
            )

    def test_inf_in_lower_bounds_raises(self):
        with pytest.raises(ValueError, match="lower_bounds"):
            StoichiometricModel(
                metabolites=["M1"], reactions=["v1"],
                stoichiometric_matrix=np.array([[1.0]]),
                lower_bounds=np.array([float("-inf")]), upper_bounds=np.array([1.0]),
            )

    def test_nan_in_upper_bounds_raises(self):
        with pytest.raises(ValueError, match="upper_bounds"):
            StoichiometricModel(
                metabolites=["M1"], reactions=["v1"],
                stoichiometric_matrix=np.array([[1.0]]),
                lower_bounds=np.array([0.0]), upper_bounds=np.array([float("nan")]),
            )

    def test_lower_exceeds_upper_raises(self):
        with pytest.raises(ValueError, match="lower_bounds must not exceed"):
            StoichiometricModel(
                metabolites=["M1"], reactions=["v1"],
                stoichiometric_matrix=np.array([[1.0]]),
                lower_bounds=np.array([5.0]), upper_bounds=np.array([1.0]),
            )

    def test_non_numeric_matrix_raises_clean_value_error(self):
        # Same "same class of bug already found in renal_pathway.py/
        # meal.py" this milestone exists to guard against -- a non-numeric
        # input must fail clearly at construction, not silently coerce
        # into garbage or crash confusingly deep inside the LP solver.
        with pytest.raises(ValueError):
            StoichiometricModel(
                metabolites=["M1"], reactions=["v1"],
                stoichiometric_matrix=[["not_a_number"]],
                lower_bounds=np.array([0.0]), upper_bounds=np.array([1.0]),
            )

    def test_arrays_are_not_mutable_after_construction(self):
        # Regression, found by code review: frozen=True only prevents
        # rebinding the attributes, not mutating the numpy arrays they
        # point to. Construction must now store defensive, read-only
        # copies so an already-validated model can't be silently
        # invalidated afterward.
        m = StoichiometricModel(
            metabolites=["M1"], reactions=["v1", "v2"],
            stoichiometric_matrix=np.array([[1.0, -1.0]]),
            lower_bounds=np.array([0.0, 0.0]), upper_bounds=np.array([10.0, 10.0]),
        )
        with pytest.raises(ValueError, match="read-only"):
            m.upper_bounds[0] = float("nan")
        with pytest.raises(ValueError, match="read-only"):
            m.lower_bounds[0] = -1.0
        with pytest.raises(ValueError, match="read-only"):
            m.stoichiometric_matrix[0, 0] = 99.0

    def test_construction_copies_the_caller_s_own_arrays(self):
        # The stored arrays must be independent copies -- mutating the
        # ORIGINAL array the caller passed in (still writable, since it's
        # the caller's own object, not the dataclass's stored copy) must
        # not retroactively change an already-constructed, already-
        # validated model.
        original_ub = np.array([10.0, 10.0])
        m = StoichiometricModel(
            metabolites=["M1"], reactions=["v1", "v2"],
            stoichiometric_matrix=np.array([[1.0, -1.0]]),
            lower_bounds=np.array([0.0, 0.0]), upper_bounds=original_ub,
        )
        original_ub[0] = 999.0
        assert m.upper_bounds[0] == pytest.approx(10.0)


class TestFVAResult:
    def test_range_width(self):
        r = FVAResult(reaction="v1", min_flux=2.0, max_flux=5.0)
        assert r.range_width == pytest.approx(3.0)

    def test_is_fixed_true_when_min_equals_max(self):
        r = FVAResult(reaction="v1", min_flux=3.0, max_flux=3.0)
        assert r.is_fixed is True

    def test_is_fixed_false_when_min_differs_from_max(self):
        r = FVAResult(reaction="v1", min_flux=0.0, max_flux=5.0)
        assert r.is_fixed is False


class TestFluxVariabilityAnalysis:
    def test_linear_chain_forces_all_reactions_to_the_same_range(self):
        # A -> M1 -> M2 -> (out): v1, v2, v3 in strict series, no branching --
        # steady state forces v1 == v2 == v3 exactly, so FVA's answer for
        # every reaction must be the INTERSECTION of all three bounds.
        model = StoichiometricModel(
            metabolites=["M1", "M2"],
            reactions=["v1_in", "v2_mid", "v3_out"],
            stoichiometric_matrix=np.array([
                [1.0, -1.0, 0.0],
                [0.0, 1.0, -1.0],
            ]),
            lower_bounds=np.array([0.0, 2.0, 0.0]),
            upper_bounds=np.array([10.0, 5.0, 10.0]),
        )
        results = {r.reaction: r for r in flux_variability_analysis(model)}
        for name in ["v1_in", "v2_mid", "v3_out"]:
            assert results[name].min_flux == pytest.approx(2.0, abs=1e-6)
            assert results[name].max_flux == pytest.approx(5.0, abs=1e-6)

    def test_branch_point_narrows_range_below_the_reaction_s_own_bound(self):
        # v1_uptake -> M1, then M1 splits into v2_branchA and v3_branchB.
        # Steady state: v1 == v2 + v3. v1's own bound is [0,5] (the binding
        # constraint), v2/v3's own bounds are [0,10] each -- but stoichiometry
        # forces each branch's actual achievable max down to 5 (all of v1's
        # capacity, with the other branch at zero), well below its own
        # individually-declared upper bound of 10. This is the genuine,
        # non-trivial FVA result the row's own milestone exists to compute.
        model = StoichiometricModel(
            metabolites=["M1"],
            reactions=["v1_uptake", "v2_branchA", "v3_branchB"],
            stoichiometric_matrix=np.array([[1.0, -1.0, -1.0]]),
            lower_bounds=np.array([0.0, 0.0, 0.0]),
            upper_bounds=np.array([5.0, 10.0, 10.0]),
        )
        results = {r.reaction: r for r in flux_variability_analysis(model)}
        assert results["v1_uptake"].min_flux == pytest.approx(0.0, abs=1e-6)
        assert results["v1_uptake"].max_flux == pytest.approx(5.0, abs=1e-6)
        assert results["v2_branchA"].min_flux == pytest.approx(0.0, abs=1e-6)
        assert results["v2_branchA"].max_flux == pytest.approx(5.0, abs=1e-6)
        assert results["v3_branchB"].min_flux == pytest.approx(0.0, abs=1e-6)
        assert results["v3_branchB"].max_flux == pytest.approx(5.0, abs=1e-6)
        # The whole point: each branch's computed max (5) is strictly below
        # its own declared upper_bound (10) -- confirms FVA is genuinely
        # reading the stoichiometric coupling, not just echoing bounds.
        assert results["v2_branchA"].max_flux < 10.0
        assert results["v3_branchB"].max_flux < 10.0

    def test_single_unconstrained_reaction_matches_its_own_bounds_exactly(self):
        # No metabolites at all coupling reactions together -- one isolated
        # reaction should get exactly its own declared bounds back.
        model = StoichiometricModel(
            metabolites=["M1"],
            reactions=["v1"],
            stoichiometric_matrix=np.array([[0.0]]),  # v1 doesn't touch M1's balance
            lower_bounds=np.array([-3.0]), upper_bounds=np.array([7.0]),
        )
        results = flux_variability_analysis(model)
        assert len(results) == 1
        assert results[0].min_flux == pytest.approx(-3.0)
        assert results[0].max_flux == pytest.approx(7.0)

    def test_infeasible_model_raises_infeasible_model_error(self):
        # v1 -> M1 -> v2 (out): steady state forces v1 == v2, but their
        # bounds ([5,10] and [0,3]) don't overlap at all -- no feasible flux.
        model = StoichiometricModel(
            metabolites=["M1"],
            reactions=["v1_in", "v2_out"],
            stoichiometric_matrix=np.array([[1.0, -1.0]]),
            lower_bounds=np.array([5.0, 0.0]),
            upper_bounds=np.array([10.0, 3.0]),
        )
        with pytest.raises(InfeasibleModelError):
            flux_variability_analysis(model)

    def test_is_fixed_true_for_a_reaction_genuinely_pinned_by_a_real_solve(self):
        # Regression-guard for the is_fixed tolerance itself, per code
        # review: a real HiGHS solve (not a hand-constructed FVAResult)
        # pinning a reaction to a single value must still report
        # is_fixed=True despite ordinary solver floating-point noise.
        # A strict linear chain forces v1==v2==v3 exactly, so with
        # identical bounds on all three, every reaction is fully pinned.
        model = StoichiometricModel(
            metabolites=["M1", "M2"],
            reactions=["v1_in", "v2_mid", "v3_out"],
            stoichiometric_matrix=np.array([
                [1.0, -1.0, 0.0],
                [0.0, 1.0, -1.0],
            ]),
            lower_bounds=np.array([3.0, 3.0, 3.0]),
            upper_bounds=np.array([3.0, 3.0, 3.0]),
        )
        results = {r.reaction: r for r in flux_variability_analysis(model)}
        for name in ["v1_in", "v2_mid", "v3_out"]:
            assert results[name].is_fixed is True

    def test_result_count_matches_reaction_count(self):
        model = StoichiometricModel(
            metabolites=["M1", "M2"],
            reactions=["v1", "v2", "v3"],
            stoichiometric_matrix=np.array([
                [1.0, -1.0, 0.0],
                [0.0, 1.0, -1.0],
            ]),
            lower_bounds=np.array([0.0, 0.0, 0.0]),
            upper_bounds=np.array([10.0, 10.0, 10.0]),
        )
        results = flux_variability_analysis(model)
        assert len(results) == 3
        assert {r.reaction for r in results} == {"v1", "v2", "v3"}
