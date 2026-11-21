# SAVE TO: packages\metabolism\src\biosim_metabolism\stoichiometry.py
"""
Stoichiometric metabolic models and Flux Variability Analysis (FVA).

Track 4 -- substrate/tissue-resolved metabolism (MASTER_PLAN row 4.3a),
built on the reframing the project adopted 2026-08-10 (renamed from its
solver, "hybrid ODE/FBA", to its biology). This module is deliberately
objective-free: the Track 4 scoping pass (see CLAUDE.md's "First real
Track 4 scoping pass" entry) found that a biomass objective function --
the conventional FBA default -- is biologically inappropriate for this
platform's postprandial-hours timescale (biomass/growth objectives
target proliferating microbes/cancer lines; hepatocytes/myocytes in a
2-5h window are post-mitotic), and that Flux Variability Analysis
(feasible-range bounds, objective-free) is the architecturally and
scientifically correct first framing -- not a point-estimate objective
choice. `flux_variability_analysis()` therefore computes, for every
reaction independently, the full range of steady-state flux values
consistent with the model's stoichiometry and bounds alone.

Uses `scipy.optimize.linprog(method="highs")` -- HiGHS ships with the
already-installed scipy (confirmed: scipy 1.17.1), so this introduces
zero new dependencies, matching the platform's free/local-only
constraint (CLAUDE.md decision 2).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import linprog


@dataclass(frozen=True)
class StoichiometricModel:
    """
    A stoichiometric metabolic network for Flux Variability Analysis.

    metabolites:            metabolite IDs (rows of the stoichiometric matrix)
    reactions:               reaction IDs (columns of the stoichiometric matrix,
                             one flux variable per reaction)
    stoichiometric_matrix:   S, shape (len(metabolites), len(reactions)) --
                             S[i, j] is the stoichiometric coefficient of
                             metabolite i in reaction j (negative =
                             consumed, positive = produced, zero = not
                             involved).
    lower_bounds:            per-reaction flux lower bound, shape (len(reactions),)
    upper_bounds:            per-reaction flux upper bound, shape (len(reactions),)

    The steady-state assumption (S @ v == 0, i.e. every metabolite's net
    production rate is zero) is the standard FBA/FVA constraint -- this
    class does not model transient/dynamic metabolite accumulation (that
    is the ODE engine's job elsewhere in this platform).
    """

    metabolites: list
    reactions: list
    stoichiometric_matrix: np.ndarray
    lower_bounds: np.ndarray
    upper_bounds: np.ndarray

    def __post_init__(self) -> None:
        # Same NaN/shape/bounds validation discipline already established
        # for renal_pathway.py/meal.py in this codebase -- fail fast at
        # construction, not deep inside the LP solver.
        if len(self.metabolites) == 0:
            raise ValueError("StoichiometricModel requires at least one metabolite")
        if len(self.reactions) == 0:
            raise ValueError("StoichiometricModel requires at least one reaction")
        if len(set(self.metabolites)) != len(self.metabolites):
            raise ValueError("StoichiometricModel metabolite IDs must be unique")
        if len(set(self.reactions)) != len(self.reactions):
            raise ValueError("StoichiometricModel reaction IDs must be unique")

        n_met, n_rxn = len(self.metabolites), len(self.reactions)
        S = np.asarray(self.stoichiometric_matrix, dtype=float)
        lb = np.asarray(self.lower_bounds, dtype=float)
        ub = np.asarray(self.upper_bounds, dtype=float)

        if S.shape != (n_met, n_rxn):
            raise ValueError(
                f"stoichiometric_matrix shape {S.shape} does not match "
                f"(n_metabolites={n_met}, n_reactions={n_rxn})"
            )
        if lb.shape != (n_rxn,):
            raise ValueError(f"lower_bounds shape {lb.shape} must be ({n_rxn},)")
        if ub.shape != (n_rxn,):
            raise ValueError(f"upper_bounds shape {ub.shape} must be ({n_rxn},)")

        if not np.all(np.isfinite(S)):
            raise ValueError("stoichiometric_matrix must not contain NaN/Inf")
        if not np.all(np.isfinite(lb)):
            raise ValueError("lower_bounds must not contain NaN/Inf")
        if not np.all(np.isfinite(ub)):
            raise ValueError("upper_bounds must not contain NaN/Inf")
        if np.any(lb > ub):
            bad = self.reactions[int(np.argmax(lb > ub))]
            raise ValueError(
                f"lower_bounds must not exceed upper_bounds (violated at "
                f"reaction {bad!r})"
            )

        # Regression, found by code review: @dataclass(frozen=True) only
        # prevents REBINDING these three attributes -- it does not copy or
        # freeze the numpy arrays themselves, which stay fully mutable
        # (e.g. `model.upper_bounds[0] = float("nan")` succeeds silently).
        # That would invalidate every check just performed above without
        # flux_variability_analysis() ever re-validating, potentially
        # surfacing later as a genuinely misleading "this should not
        # happen" InfeasibleModelError for what is actually a caller-side
        # mutation, not a solver anomaly. Store defensive copies (not the
        # caller's own arrays) and mark them read-only so any further
        # mutation attempt raises immediately, matching this codebase's
        # own "ALWAYS create new objects, NEVER mutate existing ones" rule.
        for name, arr in (("stoichiometric_matrix", S), ("lower_bounds", lb), ("upper_bounds", ub)):
            arr = arr.copy()
            arr.flags.writeable = False
            object.__setattr__(self, name, arr)

    @property
    def n_metabolites(self) -> int:
        return len(self.metabolites)

    @property
    def n_reactions(self) -> int:
        return len(self.reactions)


@dataclass(frozen=True)
class FVAResult:
    """One reaction's feasible steady-state flux range."""

    reaction: str
    min_flux: float
    max_flux: float

    @property
    def range_width(self) -> float:
        return self.max_flux - self.min_flux

    @property
    def is_fixed(self) -> bool:
        """True if this reaction's flux is pinned to a single value (min == max)."""
        # abs_tol=1e-6 (not 1e-9): found by code review to be the correct
        # scale -- HiGHS's own default feasibility/optimality tolerances
        # are ~1e-7, and this module's own test suite already compares
        # solved fluxes with abs=1e-6 (see test_stoichiometry.py). A
        # tighter 1e-9 default would risk a false negative (a genuinely
        # pinned reaction reported as NOT fixed) purely from ordinary
        # solver noise, not a real difference between min_flux/max_flux.
        return math.isclose(self.min_flux, self.max_flux, abs_tol=1e-6)


class InfeasibleModelError(RuntimeError):
    """Raised when the model's steady-state constraints admit no feasible flux at all."""


def flux_variability_analysis(model: StoichiometricModel) -> list[FVAResult]:
    """
    Compute the feasible steady-state flux range for every reaction.

    For each reaction j, solves two LPs (minimize v_j, maximize v_j)
    subject to S @ v == 0 (steady state) and lb <= v <= ub -- the
    standard FVA definition, with no objective function at all (see
    this module's own docstring for why that's the deliberate choice
    here, not merely an omission).

    Raises InfeasibleModelError if the model's own steady-state
    constraints (S @ v == 0 within the given bounds) admit no feasible
    flux vector at all -- checked once, up front, rather than letting
    every one of the 2*n_reactions per-reaction solves fail individually
    with a less informative message.
    """
    n_rxn = model.n_reactions
    S = np.asarray(model.stoichiometric_matrix, dtype=float)
    lb = np.asarray(model.lower_bounds, dtype=float)
    ub = np.asarray(model.upper_bounds, dtype=float)
    bounds = list(zip(lb.tolist(), ub.tolist()))

    # Feasibility pre-check: an arbitrary feasible point (minimize 0).
    feasibility_probe = linprog(
        c=np.zeros(n_rxn), A_eq=S, b_eq=np.zeros(model.n_metabolites),
        bounds=bounds, method="highs",
    )
    if not feasibility_probe.success:
        raise InfeasibleModelError(
            f"StoichiometricModel has no feasible steady-state flux vector "
            f"within its own bounds (HiGHS status: {feasibility_probe.message})"
        )

    results = []
    for j, reaction in enumerate(model.reactions):
        c_min = np.zeros(n_rxn)
        c_min[j] = 1.0
        res_min = linprog(c=c_min, A_eq=S, b_eq=np.zeros(model.n_metabolites), bounds=bounds, method="highs")
        if not res_min.success:
            raise InfeasibleModelError(
                f"Minimizing flux through reaction {reaction!r} failed "
                f"(HiGHS status: {res_min.message}) despite the model "
                f"passing its own feasibility pre-check -- this should not "
                f"happen for a linear program with a bounded feasible region"
            )

        c_max = np.zeros(n_rxn)
        c_max[j] = -1.0   # linprog only minimizes -- maximize v_j by minimizing -v_j
        res_max = linprog(c=c_max, A_eq=S, b_eq=np.zeros(model.n_metabolites), bounds=bounds, method="highs")
        if not res_max.success:
            raise InfeasibleModelError(
                f"Maximizing flux through reaction {reaction!r} failed "
                f"(HiGHS status: {res_max.message}) despite the model "
                f"passing its own feasibility pre-check -- this should not "
                f"happen for a linear program with a bounded feasible region"
            )

        results.append(FVAResult(
            reaction=reaction,
            min_flux=float(res_min.fun),
            max_flux=float(-res_max.fun),
        ))
    return results
