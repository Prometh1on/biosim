# SAVE TO: packages\data-seeds\src\biosim_seeds\parameter_seed\core_dynamics.py
"""Core pathway kinetics: glucose, insulin, IGF-1, mTOR, AMPK, beta-cell
secretion Km, NEFA and protein-stimulated secretion."""

from __future__ import annotations

from biosim_id_resolver import IDRegistry
from biosim_parameter_db import (
    Organism,
    ParameterKind,
    ParameterStore,
)

from ._helpers import SeedParameterBuilder, _resolve

# ---------------------------------------------------------------------------
# Glucose dynamics parameters
# ---------------------------------------------------------------------------

def _seed_glucose_dynamics(store: ParameterStore, registry: IDRegistry) -> None:
    glucose_id = _resolve(registry, "D-Glucose")

    # p1 — glucose effectiveness (insulin-independent clearance)
    # Bergman et al. 1981, J Clin Invest 68:1456-1467
    store.add(
        SeedParameterBuilder()
        .for_subject(glucose_id)
        .kind(ParameterKind.RATE)
        .value(0.0317, unit="min-1")
        .with_uncertainty(0.023, 0.041, "sd")
        .from_pubmed("7033284", figure_or_table="Table 1")
        .in_human()
        .clinical()
        .gold()
        .notes("p1: glucose effectiveness. Mean of 14 healthy subjects. "
               "Adaptation from original μU/mL insulin units to pM.")
        .build()
    )

    # (A "DeFronzo 1979 replicate" for p1, 0.028, was removed in B4,
    # 2026-11-16: PMID:382871 is the glucose-clamp methods paper and does
    # not estimate glucose effectiveness. See docs/PROVENANCE_GAPS.md.)

    # Basal glucose (Gb)
    store.add(
        SeedParameterBuilder()
        .for_subject(glucose_id)
        .kind(ParameterKind.ABUNDANCE)
        .value(5.0, unit="mM")
        .with_uncertainty(4.5, 5.5, "sd")
        .from_pubmed("7033284")
        .in_human()
        .clinical()
        .gold()
        .notes("Gb: fasting plasma glucose. Population average healthy adults. "
               "WHO normal range 3.9-5.6 mM.")
        .build()
    )


# ---------------------------------------------------------------------------
# Insulin dynamics parameters
# ---------------------------------------------------------------------------

def _seed_insulin_dynamics(store: ParameterStore, registry: IDRegistry) -> None:
    insulin_id = _resolve(registry, "insulin")
    insr_id    = _resolve(registry, "insulin receptor")

    # ke — insulin elimination rate
    # Sherwin et al. 1974, J Clin Invest 53:1481-1492
    store.add(
        SeedParameterBuilder()
        .for_subject(insulin_id)
        .kind(ParameterKind.RATE)
        .value(0.093, unit="min-1")
        .with_uncertainty(0.065, 0.121, "sd")
        .from_pubmed("4856884", figure_or_table="Table 1")
        .in_human()
        .clinical()
        .gold()
        .notes("ke: insulin elimination rate constant. t½ ≈ 7.5 min. "
               "Measured by IV insulin infusion in healthy subjects. "
               "VALUE UNVERIFIED (2026-11-14): Sherwin 1974's abstract "
               "reports a three-compartment model (plasma 4.5% body wt) "
               "and says disappearance curves 'are multiexponential'; "
               "0.093/min was not found in fetched text (PMC scan has no "
               "text layer). See docs/PROVENANCE_GAPS.md.")
        .build()
    )

    # Basal insulin (Ib)
    store.add(
        SeedParameterBuilder()
        .for_subject(insulin_id)
        .kind(ParameterKind.ABUNDANCE)
        .value(66.0, unit="pM")
        .with_uncertainty(40.0, 90.0, "range")
        .from_pubmed("7033284")
        .in_human()
        .clinical()
        .gold()
        .notes("Ib: fasting plasma insulin. ~11 μU/mL × 6 pM per μU/mL. "
               "Wide inter-individual variation (40-100 pM is normal range).")
        .build()
    )

    # p2 — X compartment clearance
    # Bergman 1981 / Pacini & Bergman 1986
    store.add(
        SeedParameterBuilder()
        .for_subject(insr_id)
        .kind(ParameterKind.RATE)
        .value(0.0123, unit="min-1")
        .with_uncertainty(0.0083, 0.0163, "sd")
        .from_pubmed("3640682", figure_or_table="Table 2")
        .in_human()
        .clinical()
        .gold()
        .notes("p2: clearance of remote insulin compartment. "
               "Refined estimate from Pacini & Bergman 1986. "
               "VALUE UNVERIFIED (2026-11-14): PMID:3640682 is the MINMOD "
               "software paper; its abstract gives no p2 value and no full "
               "text was available. See docs/PROVENANCE_GAPS.md.")
        .build()
    )

    # p3 — insulin action parameter
    # Original seed entry claimed a "Table 1" transcription of Bergman
    # 1981 (5.03e-6) that could not be reconciled with the same paper's
    # own directly-quoted prose SI value once the actual PDF was read this
    # session (see CLAUDE.md KNOWN BUGS FIXED #15) -- superseded, not
    # deleted, per the append-only store's own history-preserving design.
    old_p3_id = store.add(
        SeedParameterBuilder()
        .for_subject(insr_id)
        .kind(ParameterKind.CUSTOM, custom_name="insulin_action_p3")
        .value(5.03e-6, unit="min-1")
        .with_uncertainty(3.0e-6, 7.5e-6, "sd")
        .from_pubmed("7033284", figure_or_table="Table 1")
        .in_human()
        .clinical()
        .gold()
        .notes("p3: insulin action. Units: min⁻² per pM. "
               "Converted from original μU/mL (1 μU/mL ≈ 6 pM). "
               "Represents rate of increase of X per unit insulin above basal. "
               "SUPERSEDED: this Table 1 transcription could not be reconciled "
               "with the same paper's own Fig 6 caption (SI = 5.1 +/- 0.7 x10-4 "
               "for the combined normal groups) -- see the superseding entry.")
        .build()
    )
    store.supersede(
        old_p3_id,
        SeedParameterBuilder()
        .for_subject(insr_id)
        .kind(ParameterKind.CUSTOM, custom_name="insulin_action_p3")
        .value(9.032e-7, unit="min-1")
        .with_uncertainty(7.3e-7, 1.09e-6, "sd")   # from SI = 5.1 +/- 0.7 x10-4
        .from_pubmed("7033284", figure_or_table="Fig. 6 caption")
        .in_human()
        .clinical()
        .gold()
        .notes("p3: insulin action. Units: min⁻² per pM. Derived from "
               "Bergman 1981's own directly-quoted SI = 5.1 +/- 0.7 x10⁻⁴ "
               "min⁻¹.uU⁻¹.mL (Fig 6 caption, groups I+II+III combined, "
               "n=11, normal-tolerant subjects): p3 = SI * p2 / (pM per "
               "uU/mL, 6.945). Read directly from the paper's own PDF text "
               "(pypdf extraction), not a later replication study or a "
               "garbled OCR table -- the most reliable value obtainable "
               "from this source without full manual table transcription.")
        .build()
    )

    # Insulin receptor KM for insulin
    # De Meyts & Whittaker 2002, Nat Rev Drug Discov
    store.add(
        SeedParameterBuilder()
        .for_subject(insr_id)
        .for_object(insulin_id)
        .kind(ParameterKind.KD)
        .value(0.1, unit="nM")  # 0.1 nM = 100 pM
        .with_uncertainty(0.05, 0.5, "range")
        .from_pubmed("12360255")
        .in_human()
        .in_vitro()
        .silver()
        .notes("Kd: insulin receptor affinity for insulin. ~0.1-0.5 nM. "
               "Varies with receptor state and experimental conditions.")
        .build()
    )


# ---------------------------------------------------------------------------
# IGF-1 dynamics parameters
# ---------------------------------------------------------------------------

def _seed_igf1_dynamics(store: ParameterStore, registry: IDRegistry) -> None:
    igf1_id = _resolve(registry, "IGF-1")

    # IGF-1 basal concentration
    # Le Roith 1997, NEJM 336:633-640 (review)
    store.add(
        SeedParameterBuilder()
        .for_subject(igf1_id)
        .kind(ParameterKind.ABUNDANCE)
        .value(13.0, unit="nM")
        .with_uncertainty(8.0, 22.0, "range")
        .from_pubmed("9032050")
        .in_human()
        .clinical()
        .silver()
        .notes("IGF1b: fasting plasma IGF-1 in healthy adults age 20-40. "
               "~100 ng/mL ÷ 7.65 kDa = 13 nM. Declines with age.")
        .build()
    )

    # IGF-1 half-life -- must equal ln2 / PathwayParameters.k_igf1_clear
    # (test_seeds.py TestIgf1HalfLifeConsistency). Was 600 min credited to
    # the Le Roith 1997 review with no quoted number, while the runtime used
    # an uncited ~99 min; both replaced 2026-11-14.
    store.add(
        SeedParameterBuilder()
        .for_subject(igf1_id)
        .kind(ParameterKind.HALF_LIFE)
        .value(810.0, unit="min")  # 13.5 h, midpoint of 12-15 h
        .with_uncertainty(720.0, 900.0, "range")
        .from_pubmed("2558477", figure_or_table="abstract")
        .in_human()
        .clinical()
        .silver()
        .notes("Total serum IGF-1, dominated by the IGFBP-3/ALS ternary "
               "(200-kD) complex. Guler et al. 1989: in the 50-kD and "
               "200-kD complexes 'apparent half-lives were 20-30 min and "
               "12-15 h'; free IGF-I '10-12 min'. 92 of 152 ug/l was in "
               "the 200-kD complex. SILVER: tracer study in two adults.")
        .build()
    )


# ---------------------------------------------------------------------------
# mTOR dynamics parameters
# ---------------------------------------------------------------------------

def _seed_mtor_dynamics(store: ParameterStore, registry: IDRegistry) -> None:
    mtor_id   = _resolve(registry, "mTOR")
    insr_id   = _resolve(registry, "insulin receptor")

    # EC50 for insulin activation of mTORC1 (via AKT)
    # Sarbassov et al. 2005, Science 307:1098-1101
    store.add(
        SeedParameterBuilder()
        .for_subject(mtor_id)
        .for_object(insr_id)
        .kind(ParameterKind.EC50)
        .value(100.0, unit="pM")
        .with_uncertainty(50.0, 200.0, "range")
        .from_pubmed("15718470")
        .in_human()
        .in_vitro()
        .cell_type("HEK293")
        .bronze()
        .notes("EC50_mTOR: insulin concentration for half-maximal mTORC1 activation. "
               "Estimated from pAKT/mTOR phosphorylation dose-response curves. "
               "Wide uncertainty — indirect measurement from Western blots.")
        .build()
    )

    # Hill coefficient for mTOR activation
    store.add(
        SeedParameterBuilder()
        .for_subject(mtor_id)
        .kind(ParameterKind.HILL)
        .value(1.5, unit="dimensionless")
        .with_uncertainty(1.0, 2.5, "range")
        .from_pubmed("15718470")
        .in_human()
        .in_vitro()
        .bronze()
        .notes("Hill coefficient for mTOR/insulin dose-response. "
               "Estimated from sigmoidal shape of activation curves.")
        .build()
    )


# ---------------------------------------------------------------------------
# AMPK dynamics parameters
# ---------------------------------------------------------------------------

def _seed_ampk_dynamics(store: ParameterStore, registry: IDRegistry) -> None:
    ampk_id = _resolve(registry, "AMPK")

    # EC50 for insulin inhibition of AMPK
    # Hardie et al. 2012, Nat Rev Mol Cell Biol 13:251-262 (review)
    store.add(
        SeedParameterBuilder()
        .for_subject(ampk_id)
        .kind(ParameterKind.IC50)
        .value(200.0, unit="pM")
        .with_uncertainty(100.0, 400.0, "range")
        .from_pubmed("22436748")
        .in_human()
        .in_vitro()
        .bronze()
        .notes("IC50_AMPK: insulin concentration for half-maximal AMPK inhibition. "
               "Indirect estimate from mTOR/AMPK cross-talk literature. "
               "Direct measurement not available — bronze confidence.")
        .build()
    )


# ---------------------------------------------------------------------------
# Beta-cell GSIS EC50 (MASTER_PLAN 5.1, chunk 3 -- unlike chunks 1/2, this IS
# a CORE PathwayParameters field, not an extension-pathway constant, and
# WAS already listed in loader.py's own UNMAPPED_PARAMETERS -- this chunk
# closes that specific gap by registering it here AND wiring it into
# PARAMETER_MAPPINGS, so it is genuinely consumed by the real simulation
# path, not just stored)
# ---------------------------------------------------------------------------

def _seed_beta_cell_secretion_km(store: ParameterStore, registry: IDRegistry) -> None:
    """
    Registers `PathwayParameters.KG_secretion` (5.63 mM) -- confirmed
    empirically to make ZERO difference to any simulation output:
    `PathwayParameters()`'s own hardcoded default ALREADY equals this
    literature-derived value (the citation was already applied directly
    to the class default in a prior session; this chunk only makes it
    ALSO resolvable via ParameterStore/ParameterLoader, matching this
    campaign's own goal, not changing what value is actually used).

    GOLD (re-tiered from an earlier SILVER assignment, housekeeping-
    backlog fix -- see this file's own module docstring): Rebelato
    2018's own quoted "(EC50 10.13 mM)" is a directly-reported primary-
    source statistic, not an estimate read off a graph (unlike, e.g.,
    Sarbassov 2005's BRONZE mTOR EC50 below, which genuinely IS "estimated
    from... dose-response curves") -- matching this file's own GOLD
    standard once corrected to what it actually tracks. The only new
    arithmetic is subtracting this model's own G_threshold (4.5, exact,
    not itself uncertain), not an independent estimation step. This is
    the ONLY EC50-kind entry under subject="insulin", so re-tiering it
    has zero effect on any `best()`/`query()` selection.
    """
    insulin_id = _resolve(registry, "insulin")

    # Rebelato E, Santos LR, Carpinelli AR, Rorsman P, Abdulkader F.
    # "Short-term high glucose culture potentiates pancreatic beta cell
    # function." Sci Rep. 2018;8(1):13492. PMID:30166558 -- verified
    # live via Europe PMC full text this session (PMCID:PMC6117280,
    # open access): "In control islets insulin secretion was initiated
    # at glucose concentrations >5.6 mM (EC50 10.13 mM)." Rodent
    # pancreatic islets, non-glucose-adapted (5.6 mM control culture).
    store.add(
        SeedParameterBuilder()
        .for_subject(insulin_id)
        .kind(ParameterKind.EC50)
        .value(5.63, unit="mM")
        .from_pubmed("30166558", figure_or_table="Results text")
        .organism(Organism.OTHER)  # "Rodent" per the paper's own abstract -- not further specified as rat vs. mouse, so not claimed as either
        .in_vitro()
        .gold()
        .notes("KG_secretion: Michaelis constant for glucose-stimulated "
               "insulin secretion, expressed as EXCESS glucose above "
               "this model's own G_threshold (4.5 mM) -- Rebelato 2018's "
               "own quoted absolute-glucose EC50 (10.13 mM) minus "
               "G_threshold: 10.13 - 4.5 = 5.63. RODENT islet data (this "
               "paper's own abstract says only \"rodent,\" not rat or "
               "mouse specifically), NOT human -- disclosed, not "
               "silently presented as human-derived; used here for lack "
               "of a human-specific GSIS EC50 anchor. Previously an "
               "uncited default of 2.0 -- this literature-anchored value "
               "improved peak_insulin deviation from -44.9% to -24.7% "
               "against the OGTT clinical benchmark (see CLAUDE.md "
               "VALIDATION STATUS).")
        .build()
    )


# ---------------------------------------------------------------------------
# NEFA basal-clearance/antilipolysis-potency PK constants (MASTER_PLAN 5.1,
# chunk 24) -- registers nefa_pathway.py's own NEFA_HALF_LIFE_MIN/
# NEFA_INSULIN_IC50_PM against the newly-registered "NEFA" subject entity
# (registry_seed.py).
# ---------------------------------------------------------------------------

def _seed_nefa_pk(store: ParameterStore, registry: IDRegistry) -> None:
    """
    Registers nefa_pathway.py's NEFA_HALF_LIFE_MIN and
    NEFA_INSULIN_IC50_PM.

    Both GOLD -- unlike most prior chunks in this campaign, neither is a
    back-solve: both are directly-quoted primary-source measurements.
    NEFA_BASAL_UM (the pathway's own third real constant) is deliberately
    NOT registered here -- the module's own inline comment marks it
    "still uncited" despite two independent search passes (2026-08-12
    deep-search plan), an honest, already-disclosed gap this chunk does
    not attempt to close.
    """
    nefa_id = _resolve(registry, "NEFA")

    # Heiling VJ, Miles JM, Jensen MD. "How valid are isotopic
    # measurements of fatty acid oxidation?" Am J Physiol.
    # 1991;261(5 Pt 1):E572-7. PMID:1951681 -- verified live via NCBI
    # efetch this session. n=11 healthy volunteers. Directly quoted:
    # "After stopping tracer infusions, the half-life of plasma [14C]FFA
    # was 3.7 +/- 0.1 min." Independently corroborated (order of
    # magnitude, not exact value) by Hagenfeldt, Hagenfeldt & Wennmalm
    # 1975 (PMID:1213651, a separate research group/tracer/fatty-acid
    # species -- see nefa_pathway.py's own inline comment for the full
    # cross-check), not separately registered here since it doesn't
    # report the same directly-comparable statistic.
    store.add(
        SeedParameterBuilder()
        .for_subject(nefa_id)
        .kind(ParameterKind.CUSTOM, custom_name="nefa_half_life_min")
        .value(3.7, unit="min")
        .from_pubmed("1951681", figure_or_table="Results text")
        .in_human()
        .clinical()
        .gold()
        .notes("Plasma NEFA/FFA clearance half-life. Directly quoted: "
               "\"After stopping tracer infusions, the half-life of "
               "plasma [14C]FFA was 3.7 +/- 0.1 min.\" This platform's "
               "own NEFA_HALF_LIFE_MIN=3.7 matches exactly.")
        .build()
    )

    # Stumvoll M, Jacob S, Wahl HG, et al. "Suppression of systemic,
    # intramuscular, and subcutaneous adipose tissue lipolysis by
    # insulin in humans." J Clin Endocrinol Metab. 2000;85(11):4257-64.
    # PMID:11061533 -- verified live via NCBI efetch this session. n=13
    # lean healthy volunteers, three-step hyperinsulinaemic-euglycaemic
    # clamp. Directly quoted: "The insulin concentrations resulting in a
    # half-maximal suppression (EC50) of systemic lipolysis, adipose
    # tissue, and muscle lipolysis were 51, 68, and 44 pmol/L,
    # respectively." Uses the ADIPOSE-TISSUE-specific value (68 pM),
    # the closest physiological match to "NEFA release from adipose
    # tissue" (this pathway's own stated mechanism) -- not the systemic
    # (51 pM) or muscle-compartment (44 pM) values from the same paper.
    store.add(
        SeedParameterBuilder()
        .for_subject(nefa_id)
        .kind(ParameterKind.CUSTOM, custom_name="nefa_insulin_ic50_pm")
        .value(68.0, unit="pM")
        .from_pubmed("11061533", figure_or_table="Results text")
        .in_human()
        .clinical()
        .gold()
        .notes("Insulin concentration for half-maximal suppression of "
               "adipose-tissue lipolysis (antilipolysis EC50). Directly "
               "quoted: \"...adipose tissue... lipolysis were 51, 68, "
               "and 44 pmol/L, respectively\" -- the adipose-tissue "
               "figure (68) used, matching this platform's own "
               "NEFA_INSULIN_IC50_PM=68.0 exactly. NOT the systemic "
               "(51) or muscle (44) values from the same paper.")
        .build()
    )


# ---------------------------------------------------------------------------
# Dietary protein co-ingestion insulin-secretion constant (MASTER_PLAN
# 5.1, chunk 10; Track 1 row 1.28) -- registers
# nutrient_secretion_pathway.py's own NUTRIENT_SECRETION_PM_PER_G
# against the newly-registered generic "protein" metabolite subject
# (registry_seed.py).
# ---------------------------------------------------------------------------

def _seed_protein_secretion(store: ParameterStore, registry: IDRegistry) -> None:
    """
    Registers the protein-co-ingestion-driven insulin-secretion scale
    constant. BRONZE -- same category as the acarbose/miglitol kabs-
    scale constants: anchored to a real, precisely-quoted clinical
    ratio, but the STORED value is a disclosed engineering back-solve
    (a pM-per-gram scale reproducing that ratio through this platform's
    own composed ODE), not a directly-reported number.
    """
    protein_id = _resolve(registry, "protein")

    # Krezowski PA, Nuttall FQ, Gannon MC, Bartosh NH. "The effect of
    # protein ingestion on the metabolic response to oral glucose in
    # normal individuals." Am J Clin Nutr. 1986;44(6):847-56.
    # PMID:3538843 -- verified live via NCBI efetch this session (see
    # nutrient_secretion_pathway.py's own module-level citation comment
    # for the full derivation). n=8 normal subjects: "the area [insulin
    # AUC above fasting] was only 28% of that after glucose" for 50g
    # protein vs. 50g glucose alone.
    store.add(
        SeedParameterBuilder()
        .for_subject(protein_id)
        .kind(ParameterKind.CUSTOM, custom_name="secretion_pm_per_gram")
        .value(14.93, unit="pM per gram protein (scale on nutrient_secretion_delta)")
        .from_pubmed("3538843", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .bronze()
        .notes("Engineering back-solve, not a directly-reported number "
               "-- the literature-quoted quantity is the 28% insulin-AUC "
               "RATIO (50g protein alone vs. 50g glucose alone), not a "
               "pM/gram scale. Derived by measuring this platform's own "
               "50g-carbohydrate-alone meal insulin AUC-above-baseline "
               "(28567.0 pM.min over 240min, matching Krezowski's own "
               "~4h protocol), targeting 0.28x that for a 50g-protein-"
               "alone meal (7998.8 pM.min), then numerically rescaling "
               "against the REAL composed ODE until the target was "
               "reproduced to within 0.02% -- see "
               "nutrient_secretion_pathway.py's own module-level comment "
               "for the full two-step procedure. Does NOT capture the "
               "companion Nuttall et al. 1984 T2D finding (PMID:6389060, "
               "a genuinely SYNERGISTIC, not additive, combined effect) "
               "-- deliberately out of scope, since this platform's core "
               "ODE is calibrated on the healthy/normal cohort Krezowski "
               "represents, not T2D.")
        .build()
    )
