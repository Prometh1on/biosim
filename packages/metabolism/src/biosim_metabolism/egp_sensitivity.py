# SAVE TO: packages\metabolism\src\biosim_metabolism\egp_sensitivity.py
"""
FVA sensitivity analysis for bounding EGP_MAX_MM_MIN (Track 4,
MASTER_PLAN 4.3c).

Deliberately reports a SENSITIVITY TABLE across named scenarios, never
a single headline number -- per `fasted_hepatic_network.py`'s own
scoping finding (a fully-unbounded network's answer is a trivial
restatement of the exchange bounds, not a stoichiometric result) and
per this project's own repeated discipline against reporting false
precision when the underlying literature genuinely doesn't support one.

**The honest, load-bearing finding this analysis produced** (see
CLAUDE.md's "DONE -- Bounded EGP_MAX_MM_MIN via a sensitivity table,
not a point estimate" entry for the full account): a dedicated
follow-up literature search DID find real, usable energy-side bounds
(fatty acid availability, O2 consumption -- see `energy_bounds.py`) that
were missing when `fasted_hepatic_network.py` was first built. But
running the network with those bounds applied revealed they are NOT
actually the binding constraint once glycogen availability is
meaningfully restricted (glycogenolysis costs zero ATP/redox, so
unlimited glycogen alone saturates the artificial 1000 cap regardless
of how tightly energy is bounded) -- the STILL-UNRESOLVED carbon-side
bound (specifically glycogen depletion RATE, not just its fractional
contribution to total glucose production, which several real studies
DO report) remains the actual bottleneck. Acquiring the energy bounds
was genuinely valuable (it definitively narrows down WHERE the
remaining uncertainty lives), even though it did not, by itself,
produce a trustworthy point comparison against `EGP_MAX_MM_MIN`.

Every scenario below deliberately EXCLUDES glycogen (`ex_glycogen_bound
=0.0`) rather than guessing at a bound for it -- representing a
glycogen-DEPLETED fasting state (matching, e.g., Landau et al. 1996's
own 42h-fast timepoint, where gluconeogenesis accounts for ~93-96% of
glucose production and glycogen stores are largely exhausted) rather
than the earlier, explicitly-rejected, partly-circular attempt to
derive a glycogen bound FROM the very EGP number this analysis is
trying to independently check.
"""

from __future__ import annotations

from dataclasses import dataclass

from .energy_bounds import (
    ATPM_HEPATIC_SPECIFIC_DERIVED_UMOL_KG_MIN,
    ATPM_WHOLE_BODY_DERIVED_UMOL_KG_MIN,
    EX_ALA_MEYER2002_UMOL_KG_MIN,
    EX_FA_WHOLE_BODY_UMOL_KG_MIN,
    EX_GLYC_MEYER2002_UMOL_KG_MIN,
    EX_LAC_MEYER2002_UMOL_KG_MIN,
    EX_O2_HEPATIC_SPECIFIC_UMOL_KG_MIN,
    EX_O2_WHOLE_BODY_UMOL_KG_MIN,
)
from .flux_benchmarks import umol_per_kg_min_to_mM_per_min
from .fasted_hepatic_network import build_fasted_hepatic_network
from .stoichiometry import flux_variability_analysis


@dataclass(frozen=True)
class EGPScenarioResult:
    """One named FVA scenario's feasible glucose-output (G6PASE) range."""

    name: str
    description: str
    min_g6pase_umol_kg_min: float
    max_g6pase_umol_kg_min: float

    def max_g6pase_mM_per_min(self) -> float:
        return umol_per_kg_min_to_mM_per_min(self.max_g6pase_umol_kg_min)

    def min_g6pase_mM_per_min(self) -> float:
        return umol_per_kg_min_to_mM_per_min(self.min_g6pase_umol_kg_min)


# Named scenarios, carbon-precursor bound crossed with energy-bound
# assumption. Glycogen is EXCLUDED (=0.0) in every scenario except the
# first (which exists only to reproduce the original, already-documented
# "everything permissive" finding for comparison) -- see this module's
# own docstring for why.
_SCENARIOS: dict = {
    "fully_unbounded": (
        "Every exchange fully permissive (the original finding this "
        "whole analysis started from) -- vacuous by construction, "
        "included only for comparison against the bounded scenarios below.",
        dict(),
    ),
    "glycogen_depleted__energy_unbounded": (
        "Glycogen depleted (fasted state, matching Landau et al. 1996's "
        "42h-fast timepoint); lactate/glycerol/alanine bounded per Meyer "
        "et al. 2002; energy side (fat/O2/ATPM) left fully permissive -- "
        "isolates what the carbon-precursor bounds alone imply.",
        dict(
            ex_glycogen_bound=0.0, ex_lac_bound=EX_LAC_MEYER2002_UMOL_KG_MIN,
            ex_glyc_bound=EX_GLYC_MEYER2002_UMOL_KG_MIN, ex_ala_bound=EX_ALA_MEYER2002_UMOL_KG_MIN,
        ),
    ),
    "glycogen_depleted__energy_wholebody": (
        "Same carbon bounds as above, PLUS whole-body literature-anchored "
        "energy bounds (Byrne et al. 2005 VO2, Bisschop et al. 2002 fatty "
        "acid Ra, a derived ATPM).",
        dict(
            ex_glycogen_bound=0.0, ex_lac_bound=EX_LAC_MEYER2002_UMOL_KG_MIN,
            ex_glyc_bound=EX_GLYC_MEYER2002_UMOL_KG_MIN, ex_ala_bound=EX_ALA_MEYER2002_UMOL_KG_MIN,
            ex_o2_bound=EX_O2_WHOLE_BODY_UMOL_KG_MIN, ex_fa_bound=EX_FA_WHOLE_BODY_UMOL_KG_MIN,
            atpm_lower_bound=0.0,
        ),
    ),
    "glycogen_depleted__energy_hepaticspecific": (
        "Same carbon bounds, but with the more anatomically appropriate "
        "hepatic/splanchnic-specific energy bounds (Brundin & Wahren "
        "1994's hepatic-vein-catheterization O2 measurement) instead of "
        "the whole-body ones -- the tighter of the two energy scenarios.",
        dict(
            ex_glycogen_bound=0.0, ex_lac_bound=EX_LAC_MEYER2002_UMOL_KG_MIN,
            ex_glyc_bound=EX_GLYC_MEYER2002_UMOL_KG_MIN, ex_ala_bound=EX_ALA_MEYER2002_UMOL_KG_MIN,
            ex_o2_bound=EX_O2_HEPATIC_SPECIFIC_UMOL_KG_MIN, ex_fa_bound=EX_FA_WHOLE_BODY_UMOL_KG_MIN,
            atpm_lower_bound=0.0,
        ),
    ),
}


def run_egp_scenarios() -> list[EGPScenarioResult]:
    """
    Run every named scenario in `_SCENARIOS` and return each one's
    feasible G6PASE (glucose-6-phosphatase, the network's EGP-flux
    reaction) range.
    """
    results = []
    for name, (description, kwargs) in _SCENARIOS.items():
        model = build_fasted_hepatic_network(**kwargs)
        fva = {r.reaction: r for r in flux_variability_analysis(model)}
        g6pase = fva["G6PASE"]
        results.append(EGPScenarioResult(
            name=name, description=description,
            min_g6pase_umol_kg_min=g6pase.min_flux,
            max_g6pase_umol_kg_min=g6pase.max_flux,
        ))
    return results
