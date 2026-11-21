# SAVE TO: packages\metabolism\src\biosim_metabolism\fasted_hepatic_network.py
"""
A small, hand-curated FASTED-STATE hepatic core reaction network (Track
4, MASTER_PLAN row 4.3b), for Flux Variability Analysis via
`stoichiometry.flux_variability_analysis()`.

Scope, stated up front: this is a fasted/postabsorptive-state model
ONLY. Glucokinase, PFK-1, pyruvate kinase, and glycogen synthase are
DELIBERATELY EXCLUDED -- each is physiologically suppressed in the
fasted state every flux benchmark in this package was measured under,
and each forms a classic ATP-dissipating futile cycle with a reaction
that IS included here (glucose-6-phosphatase, FBP1, PEPCK/PC, glycogen
phosphorylase respectively); including both halves of a futile cycle
would widen every FVA range without adding physiological meaning for a
fasted-state question. This choice is fully reversible (each excluded
enzyme is one more reaction row) but means this network CANNOT
represent the postprandial liver -- a separate variant would be needed
for that.

Reactions (32) and metabolites (30) below, verified two ways before
being written down: (1) every reaction's stoichiometry independently
re-derived from standard biochemistry (glycolysis/gluconeogenesis/TCA/
beta-oxidation/urea cycle), not copied from a single source; (2)
verified programmatically -- see test_fasted_hepatic_network.py's
elemental-balance check (every internal reaction conserves C/N/P/the
adenylate moiety/the NAD and FAD pools/CoA, except the two deliberate,
explicitly-flagged untracked-export "sink" reactions), closed-exchange
loop-freedom check (clamping every exchange to zero forces every
internal reaction's feasible flux to exactly zero -- no hidden
stoichiometrically-balanced cycle that could inflate an FVA range), and
four canonical-yield checks (1 glucose per 1 glycogen unit / 2 lactate /
2 glycerol / 2 alanine, each reproducing the textbook conversion ratio
exactly when run in isolation).

H and O (water, protons) are NOT tracked -- standard simplification for
a core model at this scale; every reaction below is balanced with
respect to every OTHER element/moiety tracked (C, N, P, adenylate,
NAD/NADH, FAD/FADH2, CoA).

Two reactions are DELIBERATELY, KNOWINGLY unbalanced -- both represent
an untracked export ("demand" reaction, standard FBA convention) rather
than a real missing metabolite:
    UREA_SINK  -- exports urea untracked (-1 C, -2 N per firing)
    KETO_SINK  -- exports 1 acetoacetate (from 2 acetyl-CoA, via the
                  HMG-CoA pathway) untracked (-4 C per firing; the
                  reaction's own coefficients, `accoa: -2, coa: +2`,
                  were correct all along -- this is a wording fix only,
                  found by code review: the acetyl-group carbon count
                  (2x2=4C) was previously mislabeled as "2x
                  acetoacetate" instead of the correct "1 acetoacetate")

CRITICAL structural finding from this chunk's own scoping pass (see
CLAUDE.md's "DONE -- Curated a fasted hepatic core reaction network"
entry for the full account): with only the carbon-precursor exchange
bounds set, FVA's answer for glucose-6-phosphatase flux is a TRIVIAL
ARITHMETIC SUM of the precursor exchange upper bounds -- the network's
own stoichiometry never binds, because unlimited fatty acid/O2
availability lets every ATP/redox demand be met for free. A
MEANINGFUL, stoichiometry-driven FVA answer requires the ENERGY side
(EX_fa, EX_o2, ATPM) to be genuinely bounded too -- and no literature
anchor for those three bounds was found this session. This module
therefore exposes every exchange/maintenance bound as an explicit
keyword argument with a permissive (1000.0) default, rather than baking
in any single confident-looking "the" bound -- MASTER_PLAN 4.3c (a
separate, later chunk) is responsible for deciding how to handle this
gap, not this one.

Every metabolite ID is a plain string (`"g6p"`, not a canonical
IDRegistry ID) -- `StoichiometricModel` itself has zero dependency on
`biosim_id_resolver`/`biosim_seeds`, and this module preserves that:
registering the ~22 metabolites here that AREN'T already in
`registry_seed.py` (8 of the 30 already are: glucose, G6P, pyruvate,
lactate, acetyl-CoA, ATP, AMP, glycogen) is a separate, independent
chunk, not a prerequisite for FVA to run.

Citations for the reaction-bound-relevant literature (not every
reaction needs one -- well-established glycolysis/TCA/beta-oxidation/
urea-cycle stoichiometry is standard biochemistry, verified here by
elemental balance and canonical-yield reproduction rather than a
citation per reaction):
    Consoli A, Nurjhan N, Capani F, Gerich J. "Predominant role of
    gluconeogenesis in increased hepatic glucose production in NIDDM."
    Diabetes. 1989;38(5):550-557. PMID:2653926 -- verified live via
    PubMedFetcher this session. Nondiabetic hepatic glucose output
    12.0 +/- 0.6 umol/kg/min; NIDDM 22.7 +/- 1.0. Independently
    corroborates ADKINS_2003_BASAL_EGP (flux_benchmarks.py) -- both
    land close to 1.7-1.8x this platform's own EGP_BASAL_MM_MIN.
    Separately: EGP_MAX_MM_MIN (0.080 mM/min) converts to 18.4
    umol/kg/min -- BELOW this paper's own measured NIDDM BASAL rate
    (22.7), a real discrepancy worth its own look, independent of
    anything in this module (see CLAUDE.md entry).
    Meyer C, Stumvoll M, Dostou J, Welle S, Haymond M, Gerich J.
    "Renal substrate exchange and gluconeogenesis in normal
    postabsorptive humans." Am J Physiol Endocrinol Metab.
    2002;282(2):E428-434. PMID:11788376 -- verified live via
    PubMedFetcher this session. Systemic lactate gluconeogenesis
    1.97 +/- 0.12 umol/kg/min, "about threefold greater than that from
    glycerol, glutamine, and alanine" (each individually). Units
    ambiguity flagged explicitly below (glucose-equivalent vs.
    substrate-mass units are not disambiguated in the abstract).
    Landau BR, Wahren J, Chandramouli V, Schumann WC, Ekberg K,
    Kalhan SC. "Contributions of gluconeogenesis to glucose production
    in the fasted state." J Clin Invest. 1996;98(2):378-385.
    PMID:8755648 -- verified live via PubMedFetcher this session.
    Gluconeogenesis contributes ~47% (14h fast), ~67% (22h), ~93% (42h)
    of total glucose production -- a genuinely DIFFERENT number from
    Rothman et al. 1991's 13C-MRS-based 64+/-5% (averaged over the same
    0-22h window), and from Meyer et al. 2002's 34%/43% (uncorrected/
    TCA-exchange-corrected) precursor-sum estimate. These measure
    different things (total tracer-detectable GNG-derived carbon vs.
    net uptake of circulating precursors specifically) and this genuine
    literature conflict is NOT resolved here -- do not force one number
    by picking whichever is convenient.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from .stoichiometry import StoichiometricModel

# All 30 metabolites, in a fixed order (used to build the stoichiometric
# matrix's row order consistently).
METABOLITES = [
    "glycogen", "g1p", "g6p", "f6p", "fdp", "t3p", "pep", "oaa", "pyr",
    "lac", "accoa", "akg", "glu", "ala", "nh4", "glyc", "glyc3p", "fa",
    "co2", "o2", "glc",
    "atp", "adp", "amp", "pi",
    "nad", "nadh", "fad", "fadh2", "coa",
]

# Every INTERNAL reaction's stoichiometry (metabolite -> coefficient;
# negative = consumed, positive = produced). Exchange reactions are
# built separately below (their bounds are the tunable parameters this
# module's own docstring explains are NOT confidently literature-
# anchored for the energy side).
_INTERNAL_REACTIONS: dict = {
    # --- A. Glycogen -> glucose ---
    "PYGL":     {"glycogen": -1, "pi": -1, "g1p": 1},
    "PGMT":     {"g1p": -1, "g6p": 1},
    "G6PASE":   {"g6p": -1, "glc": 1, "pi": 1},
    # --- B. Gluconeogenic backbone ---
    "PGI":      {"g6p": -1, "f6p": 1},
    "FBP1":     {"fdp": -1, "f6p": 1, "pi": 1},
    "FBA_TPI":  {"fdp": -1, "t3p": 2},
    "GAPD_SEG": {"t3p": -1, "nad": -1, "pi": -1, "adp": -1, "pep": 1, "nadh": 1, "atp": 1},
    # PEPCK's real cofactor is GTP, not ATP -- lumped as ATP here (a
    # standard, thermodynamically harmless simplification: GTP/ATP
    # interconvert 1:1 via nucleoside diphosphate kinase in vivo, and
    # this doesn't distort the network's own energy-cost conclusions,
    # e.g. lactate/alanine gluconeogenesis's total ATP-equivalent cost
    # is unaffected). Disclosed explicitly per code review -- every
    # other simplification in this module already was, this one wasn't.
    "PEPCK":    {"oaa": -1, "atp": -1, "pep": 1, "co2": 1, "adp": 1},
    "PC":       {"pyr": -1, "co2": -1, "atp": -1, "oaa": 1, "adp": 1, "pi": 1},
    # --- C. Precursor entry ---
    "LDH":      {"lac": -1, "nad": -1, "pyr": 1, "nadh": 1},
    "GLYK":     {"glyc": -1, "atp": -1, "glyc3p": 1, "adp": 1},
    "GPD1":     {"glyc3p": -1, "nad": -1, "t3p": 1, "nadh": 1},
    "ALT":      {"ala": -1, "akg": -1, "pyr": 1, "glu": 1},
    "GLUD1":    {"glu": -1, "nad": -1, "akg": 1, "nh4": 1, "nadh": 1},
    "UREA_SINK": {"nh4": -2, "co2": -1, "atp": -3, "adp": 2, "amp": 1, "pi": 4},
    # --- D. Acetyl-CoA / TCA / energy ---
    "PDH":      {"pyr": -1, "nad": -1, "coa": -1, "accoa": 1, "co2": 1, "nadh": 1},
    "CS_IDH":   {"accoa": -1, "oaa": -1, "nad": -1, "akg": 1, "co2": 1, "coa": 1, "nadh": 1},
    "AKG_OAA":  {"akg": -1, "nad": -2, "fad": -1, "adp": -1, "pi": -1,
                 "oaa": 1, "co2": 1, "nadh": 2, "fadh2": 1, "atp": 1},
    "FAO":      {"fa": -1, "coa": -8, "atp": -1, "fad": -7, "nad": -7,
                 "accoa": 8, "amp": 1, "pi": 2, "fadh2": 7, "nadh": 7},
    "KETO_SINK": {"accoa": -2, "coa": 2},
    "ETC_NADH": {"nadh": -1, "o2": -0.5, "adp": -2.5, "pi": -2.5, "nad": 1, "atp": 2.5},
    "ETC_FADH2": {"fadh2": -1, "o2": -0.5, "adp": -1.5, "pi": -1.5, "fad": 1, "atp": 1.5},
    "ATPM":     {"atp": -1, "adp": 1, "pi": 1},
    # Modeled irreversible ([0, 1000], not in _REVERSIBLE below) even
    # though adenylate kinase is a genuinely near-equilibrium,
    # physiologically reversible enzyme -- the only direction this
    # network actually needs is ATP+AMP -> 2 ADP (recycling AMP from
    # FAO/UREA_SINK), so this doesn't block anything here, but it's a
    # real, undisclosed-until-now simplification, per code review.
    "ADK":      {"atp": -1, "amp": -1, "adp": 2},
}

# Reversible internal reactions get bounds (-1000, 1000); irreversible
# ones get (0, 1000) -- both loose/permissive by default, matching this
# module's own stated policy of not baking a specific bound into the
# network builder for anything whose real physiological range wasn't
# independently verified this session.
_REVERSIBLE = {"PGMT", "PGI", "FBA_TPI", "GAPD_SEG", "LDH", "GPD1", "ALT", "GLUD1"}

_EXCHANGE_METABOLITES = {
    "EX_glycogen": ("glycogen", 1),   # source: produces glycogen
    "EX_lac":      ("lac", 1),        # uptake: produces lac
    "EX_glyc":     ("glyc", 1),
    "EX_ala":      ("ala", 1),
    "EX_fa":       ("fa", 1),
    "EX_o2":       ("o2", 1),
    "EX_co2":      ("co2", -1),       # efflux: consumes co2
    "EX_glc":      ("glc", -1),       # efflux: consumes glc (the EGP product)
}


def build_fasted_hepatic_network(
    ex_glycogen_bound: float = 1000.0,
    ex_lac_bound: float = 1000.0,
    ex_glyc_bound: float = 1000.0,
    ex_ala_bound: float = 1000.0,
    ex_fa_bound: float = 1000.0,
    ex_o2_bound: float = 1000.0,
    ex_co2_bound: float = 1000.0,
    atpm_lower_bound: float = 0.0,
    egp_upper_bound: float = 1000.0,
) -> StoichiometricModel:
    """
    Build the 32-reaction fasted hepatic core network as a
    StoichiometricModel, ready for flux_variability_analysis().

    Every exchange/maintenance bound is an explicit keyword argument
    with a permissive (1000.0) default -- see this module's own
    docstring for exactly why none of these are baked in as a single
    confident "the literature says X" value: the ATPM/EX_fa/EX_o2
    energy-side bounds have NO literature anchor found this session,
    and even the carbon-precursor bounds (EX_lac/EX_glyc/EX_ala) rest on
    "about threefold" approximations from a single paper (Meyer et al.
    2002, PMID:11788376) with an explicitly-flagged, unresolved
    units ambiguity. `egp_upper_bound` bounds EX_glc (the measured
    output, and by mass balance at steady state, also G6PASE's own
    flux) -- default loose since this reaction's own feasible range is
    exactly the quantity a later chunk (4.3c) wants to compute, not
    something to pre-constrain here.

    Units: umol per kg body weight per min (matching
    FluxBenchmark.value_umol_per_kg_min's own convention exactly, so a
    result from this network can be compared against
    flux_benchmarks.ADKINS_2003_BASAL_EGP/DROUIN_1998_... via the SAME
    FluxBenchmark.value_mM_per_min() conversion, with no separate unit
    translation needed).

    `ex_co2_bound` sets EX_co2's UPPER bound only (its lower bound is
    always fixed at -1000.0 -- EX_co2 is always reversible, see below);
    exposed as a real parameter, not a hardcoded constant, specifically
    so a test can close it along with every other exchange to verify
    genuine closed-system loop-freedom (regression found by code
    review: the prior version hardcoded EX_co2's bound, so the
    "clamp every exchange to zero" loop-freedom test could not actually
    close EX_co2, undermining its own claimed guarantee).
    """
    reactions = dict(_INTERNAL_REACTIONS)
    exchange_bounds = {
        "EX_glycogen": ex_glycogen_bound,
        "EX_lac": ex_lac_bound,
        "EX_glyc": ex_glyc_bound,
        "EX_ala": ex_ala_bound,
        "EX_fa": ex_fa_bound,
        "EX_o2": ex_o2_bound,
        "EX_co2": ex_co2_bound,
        "EX_glc": egp_upper_bound,
    }
    for ex_name, (metabolite, sign) in _EXCHANGE_METABOLITES.items():
        reactions[ex_name] = {metabolite: sign}

    reaction_names = list(reactions.keys())
    n_met, n_rxn = len(METABOLITES), len(reaction_names)
    met_index = {m: i for i, m in enumerate(METABOLITES)}

    S = np.zeros((n_met, n_rxn))
    for j, rxn in enumerate(reaction_names):
        for met, coeff in reactions[rxn].items():
            S[met_index[met], j] = coeff

    lower_bounds = np.zeros(n_rxn)
    upper_bounds = np.zeros(n_rxn)
    for j, rxn in enumerate(reaction_names):
        if rxn in exchange_bounds:
            lower_bounds[j] = 0.0
            upper_bounds[j] = exchange_bounds[rxn]
        elif rxn == "ATPM":
            lower_bounds[j] = atpm_lower_bound
            upper_bounds[j] = 1000.0
        elif rxn in _REVERSIBLE:
            lower_bounds[j] = -1000.0
            upper_bounds[j] = 1000.0
        else:
            lower_bounds[j] = 0.0
            upper_bounds[j] = 1000.0
    # EX_co2 is reversible (blood bicarbonate is not a limiting
    # substrate for PC) -- confirmed by the scoping pass this doesn't
    # change any scenario's result, kept as the defensively correct
    # choice. Symmetric with ex_co2_bound (not a hardcoded -1000) so
    # passing ex_co2_bound=0.0 genuinely closes this exchange to
    # exactly [0, 0] -- required for the loop-freedom test to actually
    # close every exchange, not all-but-one.
    lower_bounds[reaction_names.index("EX_co2")] = -ex_co2_bound

    return StoichiometricModel(
        metabolites=list(METABOLITES),
        reactions=reaction_names,
        stoichiometric_matrix=S,
        lower_bounds=lower_bounds,
        upper_bounds=upper_bounds,
    )
