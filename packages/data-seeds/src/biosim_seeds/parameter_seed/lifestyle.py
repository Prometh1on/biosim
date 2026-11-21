# SAVE TO: packages\data-seeds\src\biosim_seeds\parameter_seed\lifestyle.py
"""Lifestyle, environmental and cyclical axes."""

from __future__ import annotations

from biosim_context_registry import ContextEntityType, ContextRegistry
from biosim_parameter_db import (
    ParameterKind,
    ParameterStore,
)

from ._helpers import SeedParameterBuilder, _resolve_context

# ---------------------------------------------------------------------------
# Cold acclimation p3 constant (MASTER_PLAN 5.1, chunk 29; Track 1 row
# 1.40) -- registers physiology/modifiers.py's own chronic cold-
# acclimation p3 factor against a new "cold_acclimation" CONTEXT entity
# (context_seed.py). Same category as sleep_restriction (chunk 12) and
# intermittent_hypoxia (chunk 22) -- a physiological/behavioral axis,
# not a molecule.
# ---------------------------------------------------------------------------

def _seed_cold_acclimation(store: ParameterStore, context_registry: ContextRegistry) -> None:
    """
    Registers the chronic cold-acclimation p3 factor at its own
    completed-protocol reference point. GOLD -- a directly-reported
    within-group percentage change, cross-checked live against the
    real cold_acclimation_modifiers() function, not a duplicated
    literal.
    """
    cold_id = _resolve_context(
        context_registry, "cold_acclimation", ContextEntityType.PHYSIOLOGICAL_AXIS,
    )

    from biosim_physiology.modifiers import cold_acclimation_modifiers

    # Hanssen MJ, Hoeks J, Brans B, et al. "Short-term cold acclimation
    # improves insulin sensitivity in patients with type 2 diabetes
    # mellitus." Nat Med. 2015;21(8):863-865. PMID:26147760 -- verified
    # live via NCBI efetch this chunk. n=8 type 2 diabetes subjects, 10
    # days of cold acclimation (14-15 degC): directly quoted "increased
    # peripheral insulin sensitivity by ~43%". T2D-ONLY evidence base --
    # the cohort had no healthy-subject arm (see physiology/
    # modifiers.py's own citation comment for the full disclosure).
    reference_factor = cold_acclimation_modifiers(True)["p3"]
    store.add(
        SeedParameterBuilder()
        .for_subject(cold_id)
        .kind(ParameterKind.CUSTOM, custom_name="p3_factor_at_reference")
        .value(reference_factor, unit="dimensionless (multiplier on PathwayParameters.p3)")
        .from_pubmed("26147760", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .gold()
        .notes("Directly-reported ~43% peripheral-insulin-sensitivity "
               "increase (1.43x) after a completed 10-day chronic "
               "cold-acclimation protocol -- computed live against the "
               "real cold_acclimation_modifiers() function, not a "
               "duplicated literal. T2D-ONLY (n=8, no healthy-subject "
               "arm) -- this platform discloses that gap via a "
               "SimulationOutput.warnings entry (runner.py) whenever "
               "cold_acclimatized is set without DiseaseState."
               "TYPE_2_DIABETES present, not silently generalized.")
        .build()
    )


# ---------------------------------------------------------------------------
# Menstrual cycle luteal-phase p3 constant (MASTER_PLAN 5.1, chunk 30;
# Track 1 row 1.36) -- registers physiology/modifiers.py's own luteal-
# phase p3 factor against a new "menstrual_cycle" CONTEXT entity
# (context_seed.py). Same category as sleep_restriction (chunk 12),
# intermittent_hypoxia (chunk 22), and cold_acclimation (chunk 29) -- a
# physiological/behavioral axis, not a molecule.
# ---------------------------------------------------------------------------

def _seed_menstrual_cycle(store: ParameterStore, context_registry: ContextRegistry) -> None:
    """
    Registers the luteal-phase p3 factor. GOLD -- a direct ratio
    computed from a single primary source's own two reported minimal-
    model SI numbers (3.20/6.20), no unit conversion -- cross-checked
    live against the real menstrual_cycle_modifiers() function, not a
    duplicated literal. FOLLICULAR/MIDCYCLE are both neutral (1.0) and
    deliberately NOT registered as separate anchors -- the source's own
    conclusion explicitly found no significant midcycle-vs-follicular
    difference, matching this campaign's established convention for an
    empirically-null point (same reasoning as pioglitazone's/
    atorvastatin's own null low-dose zones).
    """
    from biosim_physiology.profile import MenstrualCyclePhase

    cycle_id = _resolve_context(
        context_registry, "menstrual_cycle", ContextEntityType.PHYSIOLOGICAL_AXIS,
    )

    from biosim_physiology.modifiers import menstrual_cycle_modifiers

    # Valdes CT, Elkind-Hirsch KE. "Intravenous glucose tolerance test-
    # derived insulin sensitivity changes during the menstrual cycle."
    # J Clin Endocrinol Metab. 1991;72(3):642-6. PMID:1997519 --
    # verified live via NCBI efetch this chunk. n=8 cycling women,
    # tolbutamide-modified IVGTT analyzed via the Bergman minimal model.
    # Directly quoted: "The mean SI (x10(-4)/min.microU/mL) decreased in
    # a stepwise fashion from the follicular level of 6.20 +/- 0.91 ...
    # and was lowest in the luteal phase (3.20 +/- 0.25; P less than
    # 0.007)."
    reference_factors = menstrual_cycle_modifiers(MenstrualCyclePhase.LUTEAL)
    store.add(
        SeedParameterBuilder()
        .for_subject(cycle_id)
        .kind(ParameterKind.CUSTOM, custom_name="p3_factor_at_luteal")
        .value(reference_factors["p3"], unit="dimensionless (multiplier on PathwayParameters.p3)")
        .from_pubmed("1997519", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .gold()
        .notes("Direct ratio of the paper's own two reported minimal-"
               "model SI values (3.20/6.20) at the luteal vs. follicular "
               "phase -- computed live against the real "
               "menstrual_cycle_modifiers() function, not a duplicated "
               "literal. FOLLICULAR and MIDCYCLE deliberately NOT "
               "registered as separate anchors -- the source's own "
               "conclusion explicitly states 'no significant change at "
               "midcycle', despite a numerically lower point estimate "
               "(4.95), matching this campaign's own null-low-point "
               "exclusion convention. Evidence-strength caveat, "
               "disclosed not hidden: single small study (n=8), no "
               "independent human minimal-model replication located.")
        .build()
    )


# ---------------------------------------------------------------------------
# Sleep restriction p3 constant (MASTER_PLAN 5.1, chunk 12; Track 1 row
# 1.26) -- registers physiology/modifiers.py's own sleep-restriction p3
# factor against the newly-registered "sleep_restriction" CONTEXT
# entity (context_seed.py). Second chunk resolving against
# ContextRegistry -- sleep restriction is a physiological/behavioral
# axis, not a molecule.
# ---------------------------------------------------------------------------

def _seed_sleep_restriction(store: ParameterStore, context_registry: ContextRegistry) -> None:
    """
    Registers the sleep-restriction p3 factor at its own 5h/night
    (fully restricted) reference point. GOLD -- a directly-reported
    clamp-derived value, not a back-solved estimate.
    """
    sleep_id = _resolve_context(
        context_registry, "sleep_restriction", ContextEntityType.PHYSIOLOGICAL_AXIS,
    )

    # Buxton OM, Pavlova M, Reid EW, Wang W, Simonson DC, Adler GK.
    # "Sleep restriction for 1 week reduces insulin sensitivity in
    # healthy men." Diabetes. 2010;59(9):2126-33. PMID:20585000 --
    # verified live via NCBI efetch this session (see physiology/
    # modifiers.py's own citation comment for the full derivation). n=20
    # healthy men, 10h/night (replete) vs. 5h/night x7 nights
    # (restricted): "insulin sensitivity assessed by clamp was reduced
    # by 11+/-5.5% (P<0.04) after sleep restriction."
    store.add(
        SeedParameterBuilder()
        .for_subject(sleep_id)
        .kind(ParameterKind.CUSTOM, custom_name="p3_factor_at_5h_reference")
        .value(0.89, unit="dimensionless (multiplier on PathwayParameters.p3)")
        .from_pubmed("20585000", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .gold()
        .notes("Clamp-derived insulin-sensitivity factor at 5h/night "
               "(the paper's own fully-restricted condition) -- chosen "
               "over the same paper's noisier IVGTT-derived figure "
               "(-20+/-24%, wider SD). Only one point of this "
               "platform's own piecewise-linear interpolation between "
               "5h and 10h (see physiology/modifiers.py's "
               "sleep_modifiers() for the full curve) -- this single "
               "reference-point entry does NOT capture the interpolated "
               "shape between the two anchors.")
        .build()
    )


# ---------------------------------------------------------------------------
# Intermittent hypoxia p1/p3 constants (MASTER_PLAN 5.1, chunk 22;
# Track 1 row 1.32) -- registers physiology/modifiers.py's own
# intermittent-hypoxia factors against a new "intermittent_hypoxia"
# CONTEXT entity (context_seed.py). Same category as sleep_restriction
# (chunk 12) -- a physiological/behavioral axis, not a molecule.
# ---------------------------------------------------------------------------

def _seed_intermittent_hypoxia(store: ParameterStore, context_registry: ContextRegistry) -> None:
    """
    Registers the two intermittent-hypoxia factors (p3, p1) at their own
    24.3 events/h (the paper's own single tested rate) reference point,
    as a pair of GOLD entries under the SAME subject -- both are ratios
    computed directly from a single primary source's own two reported
    numbers, no unit conversion, no back-solved estimate (the same
    evidentiary category as CKD's p3_factor/ke_factor and tirzepatide's
    gb_factor_at_*mg, per this file's own module-level GOLD definition).
    """
    hypoxia_id = _resolve_context(
        context_registry, "intermittent_hypoxia", ContextEntityType.PHYSIOLOGICAL_AXIS,
    )

    from biosim_physiology.modifiers import intermittent_hypoxia_modifiers

    # Louis M, Punjabi NM. "Effects of acute intermittent hypoxia on
    # glucose metabolism in awake healthy volunteers." J Appl Physiol
    # (1985). 2009;106(5):1538-44. PMID:19265062 -- verified live via
    # NCBI efetch twice this session, by two separate parties (see
    # physiology/modifiers.py's own citation comment for the full
    # derivation). n=13 healthy volunteers, 5h intermittent hypoxia or
    # normoxia during wakefulness: "intermittent hypoxia was associated
    # with a decrease in S(I) [4.1 vs. 3.4 (mU/l)-1.min-1; P=0.0179] and
    # S(G) (1.9 vs. 1.3 min-1x10-2, P=0.0065)."
    reference_factors = intermittent_hypoxia_modifiers(24.3)
    store.add(
        SeedParameterBuilder()
        .for_subject(hypoxia_id)
        .kind(ParameterKind.CUSTOM, custom_name="p3_factor_at_reference_rate")
        .value(reference_factors["p3"], unit="dimensionless (multiplier on PathwayParameters.p3)")
        .from_pubmed("19265062", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .gold()
        .notes("S(I) ratio (3.4/4.1) at the paper's own tested 24.3 "
               "events/h intermittent-hypoxia rate -- computed live "
               "against the real intermittent_hypoxia_modifiers() "
               "function, not a duplicated literal. Only one point of "
               "this platform's own piecewise-linear interpolation "
               "between 0 and 24.3 events/h (see physiology/"
               "modifiers.py's intermittent_hypoxia_modifiers() for the "
               "full curve) -- this single reference-point entry does "
               "NOT capture the interpolated shape.")
        .build()
    )
    store.add(
        SeedParameterBuilder()
        .for_subject(hypoxia_id)
        .kind(ParameterKind.CUSTOM, custom_name="p1_factor_at_reference_rate")
        .value(reference_factors["p1"], unit="dimensionless (multiplier on PathwayParameters.p1)")
        .from_pubmed("19265062", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .gold()
        .notes("S(G) ratio (1.3/1.9) at the same 24.3 events/h "
               "reference rate as p3_factor_at_reference_rate above -- "
               "the LARGER, more statistically robust of the paper's two "
               "own findings (-31.6%, P=0.0065 vs. p3's -17.1%, "
               "P=0.0179). Computed live against the same real function, "
               "not a duplicated literal.")
        .build()
    )
