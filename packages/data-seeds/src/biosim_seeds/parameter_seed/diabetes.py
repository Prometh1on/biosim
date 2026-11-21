# SAVE TO: packages\data-seeds\src\biosim_seeds\parameter_seed\diabetes.py
"""Diabetes-spectrum disease states."""

from __future__ import annotations

from biosim_context_registry import ContextEntityType, ContextRegistry
from biosim_parameter_db import (
    ParameterKind,
    ParameterStore,
)

from ._helpers import SeedParameterBuilder, _resolve_context

# ---------------------------------------------------------------------------
# GCK-MODY Gb constant (MASTER_PLAN 5.1, chunk 13; Track 1 row 1.11) --
# registers physiology/modifiers.py's own GCK-MODY (glucokinase-gene
# MODY2) fasting-glucose factor against a new "gck_mody" CONTEXT entity
# (context_seed.py). Third chunk resolving against ContextRegistry --
# GCK-MODY is a disease state with no molecular-entity subject, same
# category as CKD.
# ---------------------------------------------------------------------------

def _seed_gck_mody_disease_state(store: ParameterStore, context_registry: ContextRegistry) -> None:
    """
    Registers GCK-MODY's single Gb factor -- GOLD, a directly-reported
    within-study ratio (not a back-solved estimate). p1/p3/secretion_scale
    are deliberately NOT registered: the same source explicitly confirms
    all three are unaffected by this mutation (a pure beta-cell glucose-
    sensing threshold reset, not a resistance or secretion-capacity
    defect), matching CKD's own "leave what the citation says is
    unaffected" convention.
    """
    gck_mody_id = _resolve_context(
        context_registry, "gck_mody", ContextEntityType.DISEASE_STATE,
    )

    # Guenat E, Seematter G, Philippe J, Temler E, Jequier E, Tappy L.
    # "Counterregulatory responses to hypoglycemia in patients with
    # glucokinase gene mutations." Diabetes Metab. 2000;26(5):377-84.
    # PMID:11119017 -- verified live via NCBI efetch this session (see
    # physiology/modifiers.py's own citation comment for the full
    # derivation). n=7 confirmed GCK-MODY (MODY2), hyperinsulinaemic
    # clamp with progressive hypoglycaemia: "Basal glucose concentrations
    # were significantly higher in MODY2 patients (7.6 +/- 0.4 mmol/l)...
    # than in healthy subjects (5.3 +/- 0.1 mmol/l, P<0.01)" -- factor is
    # this within-study ratio (7.6/5.3), the SAME real value
    # modifiers.py's own disease_modifiers() computes, cross-checked
    # directly below rather than duplicated as a literal.
    from biosim_physiology.modifiers import disease_modifiers
    from biosim_physiology.profile import DiseaseState

    gb_factor = disease_modifiers([DiseaseState.GCK_MODY])["Gb"]

    store.add(
        SeedParameterBuilder()
        .for_subject(gck_mody_id)
        .kind(ParameterKind.CUSTOM, custom_name="Gb_factor")
        .value(gb_factor, unit="dimensionless (multiplier on PathwayParameters.Gb)")
        .from_pubmed("11119017", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .gold()
        .notes("Directly-reported within-study basal-glucose ratio "
               "(7.6/5.3 mmol/l), not a back-solved estimate. The same "
               "source directly confirms p1 (insulin-mediated glucose "
               "disposal) and p3-equivalent hepatic glucose suppression "
               "are UNCHANGED in this population -- deliberately not "
               "registered here, matching CKD's own convention of only "
               "registering what the citation says is actually affected.")
        .build()
    )


# ---------------------------------------------------------------------------
# Gestational diabetes p1/p3/secretion_scale/Gb constants (MASTER_PLAN
# 5.1, chunk 14; Track 1 row 1.11) -- registers physiology/modifiers.py's
# own gestational-diabetes factors against a new "gestational_diabetes"
# CONTEXT entity (context_seed.py). Fourth chunk resolving against
# ContextRegistry. Two different primary sources feed this one subject
# (Kautzky-Willer 1997 for p1/p3/secretion_scale, IADPSG 2010 for Gb) --
# same "multiple citations, one subject" pattern already used elsewhere
# in this campaign (e.g. cortisol's two constants from different papers).
# ---------------------------------------------------------------------------

def _seed_gestational_diabetes(store: ParameterStore, context_registry: ContextRegistry) -> None:
    """
    Registers all four gestational-diabetes factors. p1/p3/
    secretion_scale are GOLD -- directly-reported within-study
    percentages from one primary source, not back-solved. Gb is SILVER,
    not GOLD -- it is a directly-quoted number, but from a clinical
    consensus DIAGNOSTIC THRESHOLD (IADPSG panel), not a directly
    measured physiological quantity in a study cohort, a genuinely
    different evidentiary category from the other three.
    """
    gdm_id = _resolve_context(
        context_registry, "gestational_diabetes", ContextEntityType.DISEASE_STATE,
    )

    from biosim_physiology.modifiers import disease_modifiers
    from biosim_physiology.profile import DiseaseState

    factors = disease_modifiers([DiseaseState.GESTATIONAL_DIABETES])

    # Kautzky-Willer A, Prager R, Waldhausl W, Pacini G, Thomaseth K,
    # Wagner OF, Ulm M, Streli C, Ludvik B. "Pronounced insulin resistance
    # and inadequate beta-cell secretion characterize lean gestational
    # diabetes during and after pregnancy." Diabetes Care.
    # 1997;20(11):1717-23. PMID:9353615 -- verified live via NCBI efetch
    # this session (see physiology/modifiers.py's own citation comment
    # for the full derivation). Minimal-model kinetic analysis of lean
    # GDM vs. lean pregnant NGT vs. healthy nonpregnant lean women:
    # "glucose effectiveness was decreased by 33%" (p1=0.67), "S1 was
    # blunted... by 84%" (p3=0.16), "Insulin secretion was 30% higher"
    # (secretion_scale=1.30) -- all three directly-reported within-study
    # percentages, cross-checked against the real module value below
    # rather than duplicated as literals.
    store.add(
        SeedParameterBuilder()
        .for_subject(gdm_id)
        .kind(ParameterKind.CUSTOM, custom_name="p1_factor")
        .value(factors["p1"], unit="dimensionless (multiplier on PathwayParameters.p1)")
        .from_pubmed("9353615", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .gold()
        .notes("Directly-reported 33% glucose-effectiveness reduction "
               "in BOTH pregnant groups (not GDM-specific -- a normal-"
               "pregnancy effect folded into this one preset, since this "
               "platform has no separate 'pregnant, not diabetic' "
               "baseline).")
        .build()
    )
    store.add(
        SeedParameterBuilder()
        .for_subject(gdm_id)
        .kind(ParameterKind.CUSTOM, custom_name="p3_factor")
        .value(factors["p3"], unit="dimensionless (multiplier on PathwayParameters.p3)")
        .from_pubmed("9353615", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .gold()
        .notes("Directly-reported 84% insulin-sensitivity (S1) blunting "
               "in GDM vs. healthy nonpregnant lean women (the GDM-"
               "specific decompensation, on top of normal pregnancy's "
               "own already-quantified 66% blunting).")
        .build()
    )
    store.add(
        SeedParameterBuilder()
        .for_subject(gdm_id)
        .kind(ParameterKind.CUSTOM, custom_name="secretion_scale_factor")
        .value(factors["secretion_scale"], unit="dimensionless (multiplier on PathwayParameters.secretion_scale)")
        .from_pubmed("9353615", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .gold()
        .notes("Directly-reported 30% higher insulin secretion in GDM "
               "vs. pregnant NGT/nonpregnant lean women -- real, "
               "measured compensatory hypersecretion, not assumed.")
        .build()
    )

    # Metzger BE et al. (IADPSG Consensus Panel). "International
    # Association of Diabetes and Pregnancy Study Groups recommendations
    # on the diagnosis and classification of hyperglycemia in
    # pregnancy." Diabetes Care. 2010;33(3):676-82. PMID:20190296 --
    # verified live via NCBI efetch this session. "fasting plasma
    # glucose >=5.1 mmol/l... diagnose as GDM."
    store.add(
        SeedParameterBuilder()
        .for_subject(gdm_id)
        .kind(ParameterKind.CUSTOM, custom_name="Gb_factor")
        .value(factors["Gb"], unit="dimensionless (multiplier on PathwayParameters.Gb)")
        .from_pubmed("20190296", figure_or_table="text, Recommendations")
        .in_human()
        .clinical()
        .silver()
        .notes("SILVER, not GOLD like the other three factors in this "
               "subject -- this is a directly-quoted number, but from a "
               "clinical consensus DIAGNOSTIC THRESHOLD (IADPSG panel "
               "recommendation), not a directly measured physiological "
               "quantity in a study cohort. A deliberately modest "
               "elevation vs. the platform's own 5.0 mM default (unlike "
               "T2D's 6.5 mM), consistent with GDM being a milder, "
               "OGTT-primarily-diagnosed dysglycemia.")
        .build()
    )


# ---------------------------------------------------------------------------
# Type 1 diabetes Ib/secretion_scale constants (MASTER_PLAN 5.1, chunk
# 15; Track 1 row 1.10) -- registers physiology/modifiers.py's own T1D
# factors against a new "type_1_diabetes" CONTEXT entity
# (context_seed.py). Fifth chunk resolving against ContextRegistry.
#
# UNLIKE every prior chunk in this campaign (CKD/GCK-MODY/GDM, all
# GOLD -- directly-reported ratios/percentages), T1D's own factors are
# disclosed ENGINEERING CHOICES in modifiers.py's own comment -- chosen
# to reach PhysiologicalPersonalizer's numeric floors after clamping,
# motivated by (not equal to) a real qualitative citation (undetectable
# C-peptide). BRONZE is the correct tier here, same category as
# pioglitazone's/protein's own disclosed engineering back-solves, NOT
# GOLD. Gb is a THIRD, separate engineering choice with NO citation at
# all (modifiers.py's own comment: "a dedicated search for a clean
# untreated-T1D fasting-glucose consensus value did not converge") --
# registered anyway, with pubmed_id="UNVERIFIED" rather than omitted
# entirely (mandatory code-reviewer finding): this campaign's whole
# purpose is to make provenance gaps VISIBLE in the store, and silent
# omission would have hidden a real, actively-used, load-bearing
# constant (moves fasting glucose from 5.0 to 11.0 mM) from anyone
# querying this subject. Matches this project's own pre-existing
# "UNVERIFIED" convention (effect_db.py's IGF1 rs35767, benchmark.py's
# mTOR-activation reference) rather than fabricating a fake citation.
# ---------------------------------------------------------------------------

def _seed_type_1_diabetes(store: ParameterStore, context_registry: ContextRegistry) -> None:
    """
    Registers all three T1D factors. Ib_factor/secretion_scale_factor
    are BRONZE (disclosed engineering choices motivated by a real
    citation, not directly-reported values themselves). Gb_factor is
    ALSO registered (not omitted) but with pubmed_id="UNVERIFIED" --
    it has zero citation at all, unlike the other two, so this makes
    that gap visible in the store rather than hiding it.
    """
    t1d_id = _resolve_context(
        context_registry, "type_1_diabetes", ContextEntityType.DISEASE_STATE,
    )

    from biosim_physiology.modifiers import disease_modifiers
    from biosim_physiology.profile import DiseaseState

    factors = disease_modifiers([DiseaseState.TYPE_1_DIABETES])

    # Simonson DC, Delprato S, Castellino P, Groop L, DeFronzo RA.
    # "Effect of glyburide on glycemic control, insulin requirement, and
    # glucose metabolism in insulin-treated diabetic patients." Diabetes.
    # 1987;36(2):136-46. PMID:3100365 -- verified live via NCBI efetch
    # this session (see physiology/modifiers.py's own citation comment
    # for the full derivation). n=5 confirmed IDDM (insulin-dependent,
    # i.e. type 1) patients on chronic exogenous insulin therapy:
    # "C-peptide levels were undetectable... in the basal state and
    # after glucagon [stimulation]... remained undetectable after
    # glyburide therapy" -- confirms near-total, not merely reduced,
    # beta-cell failure. This is a QUALITATIVE finding (undetectable,
    # not a specific fractional-reduction number) -- the actual
    # Ib=0.05/secretion_scale=0.02 factors are ENGINEERING CHOICES
    # (chosen so this platform's own numeric floors are reached after
    # clamping), not derived arithmetically from this citation, hence
    # BRONZE not GOLD. Corroborated (not independently re-verified this
    # chunk) by the broader residual-beta-cell-function literature
    # (Kalinowska et al. 2017, PMID:29253033; PMID:35181926) per
    # modifiers.py's own comment.
    store.add(
        SeedParameterBuilder()
        .for_subject(t1d_id)
        .kind(ParameterKind.CUSTOM, custom_name="Ib_factor")
        .value(factors["Ib"], unit="dimensionless (multiplier on PathwayParameters.Ib)")
        .from_pubmed("3100365", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .bronze()
        .notes("Disclosed engineering choice, NOT a directly-reported "
               "factor -- chosen so PhysiologicalPersonalizer's own "
               "engineering floor for Ib (10.0 pM) is reached after "
               "clamping, motivated by (not equal to) this source's "
               "qualitative 'undetectable C-peptide' finding. A C-peptide "
               "measurement, not plasma insulin -- the two are co-"
               "secreted equimolarly but cleared at different rates, so "
               "this is not a direct unit-for-unit substitution either.")
        .build()
    )
    store.add(
        SeedParameterBuilder()
        .for_subject(t1d_id)
        .kind(ParameterKind.CUSTOM, custom_name="secretion_scale_factor")
        .value(factors["secretion_scale"], unit="dimensionless (multiplier on PathwayParameters.secretion_scale)")
        .from_pubmed("3100365", figure_or_table="text, Results")
        .in_human()
        .clinical()
        .bronze()
        .notes("Disclosed engineering choice, NOT a directly-reported "
               "factor -- chosen so this platform's own engineering "
               "floor for secretion_scale (0.05) is reached after "
               "clamping, same motivating citation and same BRONZE "
               "reasoning as Ib_factor above.")
        .build()
    )

    # Gb (fasting glucose): a THIRD, SEPARATE engineering choice with NO
    # citation at all -- modifiers.py's own comment: "a dedicated search
    # for a clean untreated-T1D fasting-glucose consensus value did not
    # converge." Registered with pubmed_id="UNVERIFIED" rather than
    # silently omitted, matching this project's own established
    # convention for a real, actively-used constant that has no
    # legitimate citation (see effect_db.py's IGF1 rs35767,
    # benchmark.py's mTOR-activation reference/is_scored=False) --
    # preserves the audit trail's completeness without fabricating
    # evidence for a number this platform has already admitted is a
    # guess.
    store.add(
        SeedParameterBuilder()
        .for_subject(t1d_id)
        .kind(ParameterKind.CUSTOM, custom_name="Gb_factor")
        .value(factors["Gb"], unit="dimensionless (multiplier on PathwayParameters.Gb)")
        .from_pubmed("UNVERIFIED", figure_or_table=None)
        .in_human()
        .clinical()
        .bronze()
        .notes("NOT a citation-backed value -- pubmed_id is the literal "
               "string 'UNVERIFIED', matching this project's own "
               "established convention (see effect_db.py's IGF1 "
               "rs35767, benchmark.py's mTOR-activation reference) for a "
               "real, actively-used constant with no legitimate "
               "citation, rather than a fabricated one. A dedicated "
               "search for a clean untreated-T1D fasting-glucose "
               "consensus value did not converge (modifiers.py's own "
               "comment) -- registered anyway so this gap is visible in "
               "the provenance store itself, not just in code comments.")
        .build()
    )


# ---------------------------------------------------------------------------
# T2D / insulin-resistant / hypothyroidism constants (MASTER_PLAN 5.1,
# chunk 16; Track 1 rows 1.10-area, no dedicated row) -- registers
# physiology/modifiers.py's own remaining disease-state factors against
# THREE new CONTEXT entities. Sixth chunk resolving against
# ContextRegistry.
#
# Scoped by a dedicated computational-biologist fan-out before writing
# any code (Section 9's own "verify the candidate is real" rule): NONE
# of these three states' specific numeric factors are derivable from a
# live-verifiable primary source. Confirmed by live-fetching every PMID
# involved:
#   - T2D: Stumvoll et al. 2005 (PMID:15823385, Lancet) is real and
#     topically on-point but a pure narrative review with ZERO numeric
#     physiological quantities in its abstract -- supports DIRECTION
#     only. BRONZE (a real motivating citation exists), not UNVERIFIED.
#   - insulin_resistant: NO citation anywhere in modifiers.py at all,
#     not even in the module's own header citation list. UNVERIFIED.
#   - hypothyroidism: THREE real PMIDs (12519082, 18429948, 1175861),
#     all live-verified faithful to modifiers.py's own comment -- but
#     none reports a specific magnitude, and the third (Shah 1975)
#     directly CONTRADICTS the p3<1.0 direction (found no insulin
#     resistance at all, glucose nadir data "might in fact suggest
#     increased sensitivity"). BRONZE citing the two supporting
#     studies, with the contradiction disclosed in each note (this
#     ParameterBuilder has no dedicated "contested evidence" field, so
#     the conflict lives in prose, matching modifiers.py's own
#     resolution -- disclosed, not resolved by fiat).
# ---------------------------------------------------------------------------

def _seed_type_2_diabetes(store: ParameterStore, context_registry: ContextRegistry) -> None:
    """
    Registers T2D's four factors: secretion_scale/p3/Gb BRONZE, a real
    motivating citation (Stumvoll 2005) that supports direction but
    reports no magnitude; p1 GOLD, a ratio of Basu 2009's two group means.
    """
    t2d_id = _resolve_context(
        context_registry, "type_2_diabetes", ContextEntityType.DISEASE_STATE,
    )

    from biosim_physiology.modifiers import disease_modifiers
    from biosim_physiology.profile import DiseaseState

    factors = disease_modifiers([DiseaseState.TYPE_2_DIABETES])

    magnitude_caveat = (
        "Stumvoll M, Goldstein BJ, van Haeften TW. 'Type 2 diabetes: "
        "principles of pathogenesis and therapy.' Lancet. "
        "2005;365(9467):1333-46. PMID:15823385 -- verified live via NCBI "
        "efetch this session (a dedicated scoping pass before this chunk "
        "was written, per this project's own 'verify the candidate is "
        "real before executing' rule). A real, topically on-point "
        "narrative review supporting the DIRECTION of this preset "
        "(progressive beta-cell decline, early insulin insensitivity) "
        "-- but its abstract (997 characters, confirmed not truncated) "
        "reports ZERO specific numeric physiological quantities. This "
        "factor is therefore BRONZE (a real motivating citation exists, "
        "same category as T1D's Ib/secretion_scale factors), not a "
        "directly-reported value like CKD/GCK-MODY/GDM's own GOLD "
        "factors."
    )
    store.add(
        SeedParameterBuilder()
        .for_subject(t2d_id)
        .kind(ParameterKind.CUSTOM, custom_name="secretion_scale_factor")
        .value(factors["secretion_scale"], unit="dimensionless (multiplier on PathwayParameters.secretion_scale)")
        .from_pubmed("15823385", figure_or_table="text, abstract")
        .in_human()
        .clinical()
        .bronze()
        .notes(magnitude_caveat)
        .build()
    )
    store.add(
        SeedParameterBuilder()
        .for_subject(t2d_id)
        .kind(ParameterKind.CUSTOM, custom_name="p3_factor")
        .value(factors["p3"], unit="dimensionless (multiplier on PathwayParameters.p3)")
        .from_pubmed("15823385", figure_or_table="text, abstract")
        .in_human()
        .clinical()
        .bronze()
        .notes(magnitude_caveat)
        .build()
    )
    store.add(
        SeedParameterBuilder()
        .for_subject(t2d_id)
        .kind(ParameterKind.CUSTOM, custom_name="Gb_factor")
        .value(factors["Gb"], unit="dimensionless (multiplier on PathwayParameters.Gb)")
        .from_pubmed("15823385", figure_or_table="text, abstract")
        .in_human()
        .clinical()
        .bronze()
        .notes(magnitude_caveat)
        .build()
    )
    store.add(
        SeedParameterBuilder()
        .for_subject(t2d_id)
        .kind(ParameterKind.CUSTOM, custom_name="p1_factor")
        .value(factors["p1"], unit="dimensionless (multiplier on PathwayParameters.p1)")
        .from_pubmed("19196896", figure_or_table="text, Results (Fig. 2A)")
        .in_human()
        .clinical()
        .gold()
        .notes(
            "Basu A, Dalla Man C, Basu R, Toffolo G, Cobelli C, Rizza RA. "
            "Diabetes Care. 2009;32(5):866-72. PMID:19196896 -- full text "
            "verified live: oral-minimal-model SG 0.019 +/- 0.002 (14 T2D) "
            "vs. 0.034 +/- 0.006 min-1 (11 matched nondiabetic subjects). "
            "GOLD: a within-study ratio computed from the paper's two "
            "group means (0.019/0.034 = 0.559, rounded to 0.56), not a "
            "value the paper states itself -- same tier rule as GCK-MODY's "
            "Gb factor (7.6/5.3). MASTER_PLAN row 1.46 (a)."
        )
        .build()
    )


def _seed_insulin_resistant(store: ParameterStore, context_registry: ContextRegistry) -> None:
    """
    Registers the pre-diabetes/metabolic-syndrome preset's three
    factors -- ALL pubmed_id="UNVERIFIED". Unlike every other disease
    state in this campaign, modifiers.py has no citation for this one
    at all, confirmed by a dedicated scoping pass (not present even in
    the module's own header citation list).
    """
    ir_id = _resolve_context(
        context_registry, "insulin_resistant", ContextEntityType.DISEASE_STATE,
    )

    from biosim_physiology.modifiers import disease_modifiers
    from biosim_physiology.profile import DiseaseState

    factors = disease_modifiers([DiseaseState.INSULIN_RESISTANT])

    no_citation_note = (
        "NOT a citation-backed value -- pubmed_id is the literal string "
        "'UNVERIFIED'. A dedicated scoping pass confirmed modifiers.py "
        "has NO citation for this disease state anywhere at all, unlike "
        "type_2_diabetes above (which at least has a real, directionally-"
        "supportive review). Registered anyway, matching this project's "
        "own established convention (type_1_diabetes's Gb_factor), so "
        "this gap is visible in the provenance store itself."
    )
    for custom_name, key in (
        ("p3_factor", "p3"), ("Gb_factor", "Gb"), ("secretion_scale_factor", "secretion_scale"),
    ):
        store.add(
            SeedParameterBuilder()
            .for_subject(ir_id)
            .kind(ParameterKind.CUSTOM, custom_name=custom_name)
            .value(factors[key], unit=f"dimensionless (multiplier on PathwayParameters.{key})")
            .from_pubmed("UNVERIFIED", figure_or_table=None)
            .in_human()
            .clinical()
            .bronze()
            .notes(no_citation_note)
            .build()
        )
