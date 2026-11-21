# SAVE TO: packages\physiology\src\biosim_physiology\modifiers\disease.py
"""Disease-state modifiers."""

from __future__ import annotations

import functools
from typing import TYPE_CHECKING, Callable

if TYPE_CHECKING:
    from ..profile import (
        DiseaseState,
    )


# ---------------------------------------------------------------------------
# Disease state modifiers
# ---------------------------------------------------------------------------

def disease_modifiers(disease_states: list["DiseaseState"]) -> dict[str, float]:
    """
    Combine all active disease state modifiers multiplicatively.
    Multiple comorbidities compound.
    """
    combined: dict[str, float] = {}

    for state in disease_states:
        m = _single_disease_modifier(state)
        for param, factor in m.items():
            combined[param] = combined.get(param, 1.0) * factor

    return combined


def _healthy_modifiers() -> dict[str, float]:
    return {}


def _insulin_resistant_modifiers() -> dict[str, float]:
    # Pre-diabetic: ~30% reduction in insulin sensitivity,
    # slightly elevated fasting glucose, mild secretion compensation
    return {
        "p3":              0.70,
        "Gb":              1.03,
        "secretion_scale": 1.05,  # compensatory hypersecretion
    }


def _type_2_diabetes_modifiers() -> dict[str, float]:
    # Established T2D: substantial beta-cell loss and insulin resistance
    # Secretion: ~60% of normal (progressive loss varies widely)
    # Insulin sensitivity: ~55% of normal
    # Fasting glucose: elevated ~1.5 mM above normal
    #
    # Glucose effectiveness (p1), MASTER_PLAN row 1.46 (a): Basu A, Dalla
    # Man C, Basu R, Toffolo G, Cobelli C, Rizza RA. Diabetes Care.
    # 2009;32(5):866-72. PMID:19196896 (full text verified live): 14 T2D
    # vs. 11 age/weight/BMI-matched nondiabetic subjects, mixed meal,
    # oral minimal model -- "S G (0.019 +/- 0.002 vs. 0.034 +/- 0.006
    # min-1) were lower (P < 0.03) in diabetic than in nondiabetic
    # subjects". 0.019/0.034 = 0.559 -> 0.56. The control SG (0.034)
    # matches this platform's own healthy p1 (0.0317), and Basu et al.
    # 1997 (PMID:9032098) found impaired basal glucose effectiveness in
    # an independent NIDDM cohort. Same cohort as validate_t2d_meal_insulin().
    # Fixed-point neutral: p1 multiplies (G - Gb), which is zero at basal.
    # Does NOT close row 1.46's insulin gap (structural, see MASTER_PLAN).
    return {
        "secretion_scale": 0.60,
        "p3":              0.55,
        "Gb":              1.30,  # ~6.5 mM fasting glucose
        "p1":              0.56,
    }


def _type_1_diabetes_modifiers() -> dict[str, float]:
    # Autoimmune beta-cell destruction -- near-total loss of
    # ENDOGENOUS insulin secretion, mechanistically opposite to
    # TYPE_2_DIABETES/INSULIN_RESISTANT above (which represent
    # resistance with partially-preserved, often compensatorily
    # elevated, secretion). This is the disease-state prerequisite
    # for T1D scenarios that also compose the already-built
    # ExogenousInsulinPathway (exogenous_insulin_pathway.py) --
    # confirmed via a dedicated fixed-point/trajectory investigation
    # (CLAUDE.md "Confirmed the ODE-engine already supports
    # near-total beta-cell failure") that basal_state() already
    # derives every species correctly from ANY (Ib, secretion_scale)
    # pair, including near-zero ones, with no further ODE-engine
    # change required -- this entry supplies the literature-anchored
    # numeric preset that investigation deliberately left open.
    #
    # Ib and secretion_scale factors, verified live via PubMedFetcher:
    # Simonson DC, Delprato S, Castellino P, Groop L, DeFronzo RA.
    # "Effect of glyburide on glycemic control, insulin requirement,
    # and glucose metabolism in insulin-treated diabetic patients."
    # Diabetes. 1987;36(2):136-46. PMID:3100365. n=5 confirmed IDDM
    # (insulin-dependent, i.e. type 1) patients on chronic exogenous
    # insulin therapy (54+/-8 U/day): "C-peptide levels were
    # undetectable (less than 0.01 pmol/ml) in the basal state and
    # after glucagon [stimulation]... remained undetectable after
    # glyburide therapy" -- i.e. genuinely unresponsive to a
    # secretagogue, confirming near-total, not merely reduced, beta-
    # cell failure. Corroborated by the broader residual-beta-cell-
    # function literature (Kalinowska et al. 2017, PMID:29253033;
    # PMID:35181926; both independently verified live this session):
    # most confirmed long-duration T1D patients have
    # C-peptide at or below assay detection limits, though a
    # minority retain measurable residual secretion -- this platform
    # represents the modal ("near-total, not absolute zero") case.
    #
    # Ib's <0.01 pmol/ml (<10 pM) C-peptide finding is a C-PEPTIDE
    # measurement, not plasma INSULIN -- the two are co-secreted
    # equimolarly but cleared at different rates, so this is not a
    # direct unit-for-unit substitution. Rather than invent an
    # unverified C-peptide-to-insulin clearance conversion, the
    # factor below is chosen so PhysiologicalPersonalizer's own
    # engineering floor for Ib (_PARAM_BOUNDS, 10.0 pM) is reached
    # after clamping -- treating this platform's own numerical floor
    # as the operational proxy for "at or below what this model can
    # represent as distinct from zero", consistent in spirit with
    # the primary source's own "undetectable" (i.e. at-or-below-
    # assay-floor) finding, not presented as a precise measured
    # plasma insulin concentration.
    #
    # Gb (fasting glucose) is NOT derived from either citation above
    # -- deliberately disclosed as an ENGINEERING CHOICE, not a
    # literature-anchored magnitude (a dedicated search for a clean
    # untreated-T1D fasting-glucose consensus value did not converge,
    # see MASTER_PLAN.md row 1.10's own search history). Left
    # untouched, this preset's near-zero Ib/secretion_scale would
    # still form a mathematically valid analytic fixed point (rule
    # 11 F_secretion_basal solves for ANY (Gb, Ib) pair) but at an
    # implausible, self-contradictory one: normal fasting glucose
    # maintained despite near-zero insulin, with nothing else in
    # this model regulating it. Set high enough to clearly exceed
    # Brunzell's own AIR-absent threshold (~6.39 mM,
    # personalizer.py's BRUNZELL_AIR_ABSENT_THRESHOLD_MM -- so
    # first_phase_gain is correctly zeroed by that EXISTING,
    # already-tested mechanism, with no new special-casing needed)
    # and well above TYPE_2_DIABETES's own 6.5 mM (reflecting more
    # severe, total insulin deficiency), while staying safely within
    # _PARAM_BOUNDS' 15.0 mM ceiling. Weak corroborating context,
    # not a precise anchor: Simonson et al.'s own T1D cohort ran a
    # mean 24h glucose of ~9.83 mM (177 mg/dl) even ON substantial
    # (54 U/day) exogenous insulin -- a genuinely untreated preset
    # baseline should plausibly exceed that treated value.
    return {
        "Ib":              0.05,   # -> 66.0*0.05=3.3, clamped to floor 10.0 pM
        "secretion_scale": 0.02,   # -> 1.0*0.02=0.02, clamped to floor 0.05
        "Gb":              2.2,    # -> 5.0*2.2=11.0 mM
    }


def _hypothyroidism_modifiers() -> dict[str, float]:
    # Reduced peripheral glucose disposal/insulin action (p1, p3) --
    # directionally supported by Rochon et al. 2003 (Clin Sci,
    # PMID:12519082, n=6, clamp+tracer: basal glucose disposal
    # increased after restoring euthyroidism) and Handisurya et al.
    # 2008 (Clin Endocrinol, PMID:18429948, n=12 overt hypothyroid,
    # euglycaemic-hyperinsulinaemic clamp: insulin sensitivity
    # improved after T4 replacement). A genuine, unresolved conflict
    # exists (Shah et al. 1975, Diabetes, PMID:1175861, n=10, direct
    # insulin-tolerance-test measurement, found NO insulin
    # resistance) -- not resolved here, flagged per this project's
    # standing rule against picking a side by fiat.
    #
    # Regression: Gb (fasting glucose) previously REMOVED here --
    # was an uncited 1.04, and is actively contradicted by evidence:
    # this function's own cited source, Brenta 2011 (J Thyroid Res,
    # PMID:21941681), states verbatim (verified against its own full
    # text) that "a poor utilization of glucose in hypothyroidism
    # may be offset by a reduced release to circulation maintaining
    # a balance of the glucose metabolism" -- i.e. reduced hepatic
    # glucose output (which pushes glucose down) offsets reduced
    # peripheral disposal (which pushes glucose up). Two independent
    # human clamp studies directly confirm the net: Handisurya 2008
    # ("Fasting and dynamic glycaemia as assessed by the AUC(Glucose)
    # remained unaltered following substitution therapy") and Rochon
    # 2003 ("basal plasma insulin and glucose levels... were similar
    # in each case"). No source was found supporting an elevated
    # fasting glucose in hypothyroidism -- removed rather than left
    # as a plausible-looking but unsupported, wrong-signed number.
    return {
        "p1":  0.85,
        "p3":  0.88,
    }


def _hyperthyroidism_modifiers() -> dict[str, float]:
    # Deliberately returns {} -- NOT a missing case. Hyperthyroidism's
    # real, literature-supported effect on fasting/postprandial
    # glucose is hepatic in origin and genuinely insulin-gated (see
    # Wennlund et al. 1986, JCEM, PMID:3510000, n=8 hyperthyroid vs.
    # 10 normal, hepatic venous catheterization + isotope dilution:
    # basal hepatic glucose production +20% AND, separately, LESS
    # suppressed by a larger insulin rise -- a rightward-shifted
    # insulin dose-response, not a static offset), which is exactly
    # what EGPParameters.for_hyperthyroidism()/HepGlucPathway (see
    # hepatic.py) represents. A static Gb modifier here alongside
    # that coupling would double-count the identical mechanism, the
    # same bug class already fixed once for glucocorticoid/
    # CortisolPathway (see CLAUDE.md KNOWN BUGS) -- the ODE pathway
    # is the sole representation, matching
    # _glp1_agonist_modifiers()'s/_sglt2_inhibitor_modifiers()'s
    # existing empty-dict precedent for the same reason. Peripheral
    # insulin-stimulated glucose disposal in hyperthyroidism is
    # reported as normal-or-increased, not reduced (Dimitriadis &
    # Raptis 2001, Exp Clin Endocrinol Diabetes, PMID:11460573:
    # "rates of insulin-stimulated glucose disposal in peripheral
    # tissues in hyperthyroidism have been found, in general, either
    # normal or increased"), so no p3 modifier is added either rather
    # than guessing a magnitude for an increase.
    return {}


def _gestational_diabetes_modifiers() -> dict[str, float]:
    # Pregnancy-specific insulin resistance + inadequate compensatory
    # secretion, diagnosed in the 2nd/3rd trimester. Unlike every
    # other disease state above, this platform has no "pregnant, not
    # diabetic" baseline to model AS A SEPARATE LAYER -- both the
    # normal-pregnancy insulin-resistance shift and GDM's additional
    # decompensation are folded directly into ONE preset here (the
    # same "one disease state, not a stacked baseline" design already
    # used by every other entry in this module), anchored to a
    # single study that measured GDM against a TRUE non-pregnant
    # reference, not just a matched-pregnant comparator.
    #
    # Verified live via PubMedFetcher: **Kautzky-Willer A, Prager R,
    # Waldhausl W, Pacini G, Thomaseth K, Wagner OF, Ulm M, Streli C,
    # Ludvik B. "Pronounced insulin resistance and inadequate
    # beta-cell secretion characterize lean gestational diabetes
    # during and after pregnancy." Diabetes Care. 1997;20(11):
    # 1717-23. PMID:9353615.** Mathematical modelling (minimal-model
    # kinetic analysis of OGTT/IVGTT glucose-insulin-C-peptide data --
    # the same modelling tradition Pacini's own co-authorship already
    # anchors this platform's p2/KG_secretion-adjacent parameters to)
    # comparing lean GDM, lean pregnant women with normal glucose
    # tolerance (NGT), and healthy nonpregnant lean women:
    # - "S1 [insulin sensitivity] was blunted in pregnant women with
    #   GDM by 84%... and in those with NGT by 66% compared with
    #   lean nonpregnant women" -- confirms normal (non-diabetic)
    #   pregnancy ALREADY causes major insulin resistance (66%) on
    #   its own, with GDM adding a real, further, quantified
    #   decompensation on top (84% vs. 66%).
    # - "glucose effectiveness was decreased by 33% in both pregnant
    #   groups... vs. healthy nonpregnant lean control subjects" --
    #   a pregnancy effect, not GDM-specific (applies equally to
    #   both pregnant groups), so folded into this preset too since
    #   there is no separate "pregnant, NGT" state to carry it.
    # - "Insulin secretion was 30% higher... in subjects with GDM
    #   than in pregnant women with NGT or in nonpregnant lean
    #   women" -- real, measured compensatory hypersecretion (same
    #   "insufficient but real" pattern as INSULIN_RESISTANT's own
    #   secretion_scale=1.05 above, just independently anchored here
    #   rather than assumed).
    #
    # Gb: NOT derived from this paper (it reports S1/glucose-
    # effectiveness/secretion percentages, not a fasting glucose
    # value) -- anchored instead to the internationally-adopted
    # clinical diagnostic threshold itself, verified live via full
    # text: **Metzger BE et al. (IADPSG Consensus Panel).
    # "International Association of Diabetes and Pregnancy Study
    # Groups recommendations on the diagnosis and classification of
    # hyperglycemia in pregnancy." Diabetes Care. 2010;33(3):
    # 676-82. PMID:20190296.** -- "fasting plasma glucose >=5.1
    # mmol/l (92 mg/dl)... diagnose as GDM." A deliberately modest
    # elevation vs. the default 5.0 mM (unlike T2D's 6.5 mM) --
    # consistent with GDM being a milder, OGTT-primarily-diagnosed
    # dysglycemia, not the same order of fasting hyperglycemia as
    # established T2D.
    #
    # Not touched: beta-cell "sensitivity to glucose" (decreased
    # 40-50% in GDM vs. the other two groups per the same paper) --
    # a genuinely different quantity from secretion_scale (a glucose-
    # SENSING parameter, closer in spirit to KG_secretion, which no
    # other disease state in this module touches either) --
    # deliberately not introduced as a new modifier axis for one
    # disease state alone; secretion_scale already carries this
    # platform's general "beta-cell function" representation for
    # every other entry above, and GDM's REAL measured impact on it
    # (a net increase, not a decrease, per the actual absolute
    # secretion figure) is what's represented here.
    return {
        "p1":              0.67,
        "p3":              0.16,
        "secretion_scale": 1.30,
        "Gb":              1.02,   # -> 5.0*1.02=5.1 mM (IADPSG threshold)
    }


def _gck_mody_modifiers() -> dict[str, float]:
    # Glucokinase (GCK, MODY2) mutations directly impair the beta
    # cell's own glucose SENSOR -- glucokinase is the rate-limiting
    # enzyme that lets the beta cell "measure" ambient glucose, so a
    # heterozygous inactivating mutation resets the glucose set-point
    # upward, WITHOUT reducing maximal secretory capacity or
    # peripheral insulin sensitivity. Mechanistically the OPPOSITE
    # shape from every other diabetes-related state above (T1D =
    # secretion loss, T2D/GDM/INSULIN_RESISTANT = resistance with
    # partial/compensatory secretion) -- this is a pure threshold
    # reset with everything else intact, matching this condition's
    # own well-established clinical reputation as mild, stable, and
    # rarely requiring treatment.
    #
    # Verified live via PubMedFetcher: **Guenat E, Seematter G,
    # Philippe J, Temler E, Jequier E, Tappy L. "Counterregulatory
    # responses to hypoglycemia in patients with glucokinase gene
    # mutations." Diabetes Metab. 2000;26(5):377-84. PMID:11119017.**
    # n=7 confirmed GCK-MODY (MODY2) patients, n=5 T2D, n=13 healthy
    # controls, hyperinsulinaemic clamp with progressive hypoglycaemia
    # (a small but real, controlled, directly-measured cohort --
    # GCK-MODY is rare, so this is the standard evidence base for
    # this population, same category as the small "lean GDM" cohort
    # already cited above). Confirmed verbatim: "Basal glucose
    # concentrations were significantly higher in MODY2 patients
    # (7.6 +/- 0.4 mmol/l) and type 2 diabetic patients (12.4 +/- 2.3
    # mmol/l) than in healthy subjects (5.3 +/- 0.1 mmol/l, P<0.01)"
    # -- Gb factor below is this within-study ratio (7.6/5.3=1.434),
    # applied to this platform's own default, same convention as
    # every other disease-state preset in this module.
    #
    # Critically, the SAME paper directly confirms p1/p3 should stay
    # UNCHANGED, not guessed at zero effect by omission: "Insulin-
    # mediated glucose disposal and suppression of endogenous glucose
    # production at euglycemia were unchanged in MODY2 patients, but
    # were blunted in type 2 diabetes" -- a real, measured null
    # result for peripheral insulin action, directly contrasted
    # against T2D in the SAME cohort. secretion_scale is also left
    # untouched -- no citation found (or expected, given the
    # mechanism) for a genuine maximal-secretory-capacity change;
    # PMID:36342518 (Henquin/dorzagliatin GCK-MODY study, verified
    # live) independently corroborates the "rightward-shifted
    # response curve, not a reduced one" shape ("upward and leftward
    # shift in ISR-glucose response" when a glucokinase activator
    # DOES improve it -- implying the untreated baseline is shifted
    # in the opposite direction, consistent with a threshold reset
    # rather than reduced capacity).
    #
    # Known, disclosed uncertainty (not resolved here): this Gb value
    # (7.17 mM) exceeds BRUNZELL_AIR_ABSENT_THRESHOLD_MM (~6.39 mM,
    # personalizer.py), so the existing, generic, glucose-level-based
    # first_phase_gain-zeroing gate (Brunzell et al. 1976, added for
    # T2D) fires automatically for this preset too. Whether that
    # generalizes correctly to GCK-MODY's specific threshold-reset
    # mechanism (as opposed to T1D's genuine secretion loss, the
    # population the gate is best-supported for) is NOT independently
    # confirmed by any citation found this pass -- the source above
    # addresses hypoglycaemic counter-regulation, not hyperglycaemic
    # first-phase secretion, in this population. Left as-is (letting
    # the existing, real, glucose-level-keyed rule apply uniformly)
    # rather than inventing an unsupported carve-out for one disease
    # state -- flagged honestly rather than silently accepted or
    # silently overridden either way.
    return {
        "Gb": 7.6 / 5.3,   # -> 5.0*1.4340=7.17 mM
    }


def _chronic_kidney_disease_modifiers() -> dict[str, float]:
    # Moderate-severe, NON-DIABETIC CKD reduces both insulin
    # SENSITIVITY (p3) and insulin CLEARANCE (ke) -- a genuinely
    # different shape from every disease state above, none of which
    # touch ke. Mechanistically distinct from T2D/INSULIN_RESISTANT
    # (resistance with partially-preserved secretion): here secretion
    # itself is explicitly confirmed UNCHANGED by the same source, so
    # secretion_scale is deliberately left untouched (return dict has
    # no secretion_scale key), matching GCK_MODY's own "leave what
    # the citation says is unaffected" convention above.
    #
    # Verified live via PubMedFetcher: **de Boer IH, Zelnick L,
    # Afkarian M, Ayers E, Curtin L, Himmelfarb J, Ikizler TA, Kahn
    # SE, Kestenbaum B, Utzschneider K. "Impaired Glucose and Insulin
    # Homeostasis in Moderate-Severe CKD." J Am Soc Nephrol.
    # 2016;27(9):2861-71. PMID:26823551.** n=59 nondiabetic CKD (mean
    # eGFR 37.6 ml/min/1.73m2) vs. n=39 healthy controls,
    # hyperinsulinaemic-euglycaemic clamp + IVGTT/OGTT. Confirmed
    # verbatim: "Participants with CKD had lower insulin sensitivity
    # than participants without CKD (mean[SD] 3.9[2.0] versus 5.0
    # [2.0] mg/min per uU/ml; P<0.01)... Insulin clearance... was
    # also lower in participants with CKD than controls (876 [226]
    # versus 998 [212] ml/min; P<0.01)."
    #
    # The crude ratios (3.9/5.0=0.78, 876/998=0.878) are NOT what's
    # used below -- this platform already applies bmi_modifiers()/
    # fitness modifiers separately, and the same paper reports the
    # CKD-specific residual after adjusting for exactly those
    # confounds, avoiding a double-counting guess: "Adjustment for
    # physical activity, diet, fat mass, and fat-free mass in
    # addition to demographics and smoking partially attenuated
    # associations of CKD with insulin sensitivity (adjusted
    # difference, -0.7... mg/min per uU/ml) and insulin clearance
    # (adjusted difference, -85... ml/min)." Applying these deltas to
    # the paper's own control-group means gives the CKD-specific
    # factors used here: p3 = (5.0-0.7)/5.0 = 0.86; ke =
    # (998-85)/998 = 0.915.
    #
    # Disclosed uncertainty, not resolved here: the adjusted insulin-
    # sensitivity difference's own 95% CI is -1.4 to 0.0 mg/min per
    # uU/ml -- the upper bound touches exactly zero, meaning this
    # adjusted p3 effect is not statistically robust in this cohort
    # (n=59), unlike every other disease-state citation in this
    # function. Used anyway (matching this project's own "an
    # imperfect real anchor beats an invented one" convention), but
    # flagged honestly rather than presented as equally solid.
    #
    # Secretion explicitly confirmed unaffected by the same source:
    # "We did not observe a CKD-specific deficit in insulin
    # secretion" -- secretion_scale stays at its neutral default.
    #
    # Not derived from this citation, deliberately left untouched:
    # Gb (fasting glucose) and Ib (fasting insulin) -- the source
    # reports insulin tolerance-test-derived sensitivity/clearance
    # rates, not absolute fasting concentrations for either, and its
    # own full text is not openly accessible to extract them (PMC
    # explicitly blocks XML download for this article). Inventing
    # either from the crude clamp ratios would risk double-counting
    # against ke's own already-applied reduction (fasting insulin is
    # a function of both basal secretion AND clearance). A future
    # chunk could add a real Gb/Ib anchor if one is found -- not
    # guessed here.
    return {
        "p3": 0.86,
        "ke": 0.915,
    }


def _pcos_modifiers() -> dict[str, float]:
    # Polycystic ovary syndrome (Track 1 row 1.29) -- peripheral
    # insulin resistance. Mechanistically distinct from every other
    # disease state above (not a beta-cell secretion-loss or hepatic-
    # glucose-production state -- both explicitly confirmed NULL by
    # the same sources, same "leave what the citation says is
    # unaffected" convention as CKD/GCK-MODY above).
    #
    # Three sources independently re-verified live via NCBI efetch
    # before this chunk (not trusted from an earlier scoping report
    # alone):
    #
    # **Dunaif A, Segal KR, Futterweit W, Dobrjansky A. "Profound
    # peripheral insulin resistance, independent of obesity, in
    # polycystic ovary syndrome." Diabetes. 1989;38(9):1165-74.
    # PMID:2670645.** Euglycemic clamp, obese (n=19) and nonobese
    # (n=10) PCO women vs. matched normal ovulatory controls (n=11
    # obese, n=8 nonobese). Confirmed verbatim: "PCO women have
    # significant insulin resistance that is independent of obesity,
    # changes in body composition, and impairment of glucose
    # tolerance" and "the metabolic clearance rate of insulin did
    # not differ in the four groups" (ke unaffected).
    #
    # **Lewy VD, Danadian K, Witchel SF, Arslanian S. "Early
    # metabolic abnormalities in adolescent girls with polycystic
    # ovarian syndrome." J Pediatr. 2001;138(1):38-44. PMID:11148510.**
    # n=12 obese PCOS adolescents vs. n=10 obese non-hyperandrogenic
    # controls, hyperinsulinaemic-euglycaemic clamp. Confirmed
    # verbatim: "insulin sensitivity was lower in the PCOS group
    # (1.4+/-0.2 vs 2.7+/-0.3 mg/kg/min per microU/mL, P=.002)" (p3
    # ratio 1.4/2.7=0.52) and "Fasting glucose and hepatic glucose
    # production were comparable between the 2 groups" (Gb/hepatic
    # production unaffected).
    #
    # **Ketel IJ et al. "Insulin-induced capillary recruitment is
    # impaired in both lean and obese women with PCOS." Hum Reprod.
    # 2011;26(11):3130-7. PMID:21920943.** n=36 PCOS (20 lean, 16
    # obese) vs. n=27 matched controls (14 lean, 13 obese), clamp.
    # Confirmed verbatim: obese PCOS vs. obese controls M/I "0.5...
    # versus 0.8... (P<0.01)" (ratio 0.625, SIGNIFICANT); lean PCOS
    # vs. lean controls M/I "1.5... versus 1.7... (P=0.17)" (ratio
    # 0.882, NOT significant).
    #
    # A GENUINE, UNRESOLVED literature conflict exists on whether the
    # resistance is obesity-INDEPENDENT (Dunaif's own explicit
    # conclusion) or obesity-GATED (Ketel's own lean-vs-lean
    # comparison finds no significant difference) -- disclosed here,
    # not resolved by fiat, the same discipline already applied to
    # hypothyroidism's own three-way conflict above. Anchor used:
    # Ketel's own obese-stratum ratio (0.625) -- the most recent of
    # the three citations -- disclosed explicitly as OBESE-specific,
    # not a general PCOS magnitude; this platform has no separate
    # lean/obese PCOS split to represent Ketel's own lean-arm null
    # differently, so the obese figure is applied regardless of the
    # composed profile's own BMI (a real, disclosed limitation, not
    # hidden). Lewy's own obese-cohort ratio (1.4/2.7=0.52, also
    # directly reported and statistically significant, P=.002) was
    # a real, comparably-qualified candidate too -- NOT passed over
    # for a qualitative-vs-quantitative reason, but for the same
    # "most recent primary source" tiebreaker already used elsewhere
    # in this file when multiple directly-reported numbers compete
    # (code-reviewer finding: named explicitly here rather than left
    # for a future reader to re-derive that Lewy also had a usable
    # obese-cohort number). GOLD tier -- a ratio computed directly
    # from a single primary source's own two reported numbers, no
    # unit conversion, matching this file's own recently-clarified
    # GOLD definition
    # (see parameter_seed.py's own module docstring, Section 8
    # housekeeping).
    #
    # ke/Gb/hepatic glucose production explicitly confirmed NULL by
    # Dunaif (ke) and Lewy (Gb/hepatic production) -- secretion_scale
    # also not touched: none of the three sources report a
    # secretion-specific magnitude (Lewy's own "hyperinsulinemia"
    # finding is a compensatory consequence of the resistance itself,
    # not an independently-measured secretion-capacity change, same
    # "compensatory, not independent" reasoning already applied to
    # hydrochlorothiazide/atorvastatin's own analogous findings).
    return {"p3": 0.625}


def _acromegaly_modifiers() -> dict[str, float]:
    # Chronic growth-hormone (and downstream IGF-1) excess (Track 1
    # row 1.37) -- Gb-ONLY, matching the already-accepted partial-
    # implementation precedent (GCK_MODY above, tirzepatide's own
    # Gb-only route in drug_interactions.py) rather than guessing an
    # insulin-sensitivity magnitude the source doesn't report.
    #
    # Verified live via PubMedFetcher (literature-curator agent,
    # independently re-confirmed, not trusted from an earlier
    # scoping report): **Moller N, Schmitz O, Jorgensen JO, Astrup
    # J, Bak JF, Christensen SE, Alberti KG, Weeke J. "Basal- and
    # insulin-stimulated substrate metabolism in patients with
    # active acromegaly before and after adenomectomy." J Clin
    # Endocrinol Metab. 1992;74(5):1012-9. PMID:1569148.** n=6
    # active (pre-adenomectomy) acromegaly vs. n=7 controls.
    # Confirmed verbatim: "increased plasma glucose (5.9 +/- 0.2 vs.
    # 5.2 +/- 0.2 mmol/L)... (P<0.05)" -- Gb factor below is this
    # within-study ratio (5.9/5.2=1.1346), applied to this
    # platform's own default, same convention as every other
    # disease-state preset in this module.
    #
    # Insulin-sensitivity magnitude explicitly NOT extractable, not
    # guessed at: "During the clamp the patients with active
    # acromegaly were substantially resistant to the actions of
    # insulin on both glucose and lipid metabolism" -- qualitative
    # only, no glucose-infusion-rate/M-value number anywhere in the
    # source. p3/ke/secretion_scale therefore deliberately left
    # untouched (return dict has no such key) -- an honest partial
    # implementation, not an assumed-null claim (unlike GCK_MODY/
    # CKD/PCOS above, where the source explicitly confirms a null
    # rather than simply omitting a number).
    #
    # Post-adenomectomy figure NOT used: the abstract reports only
    # the qualitative "following pituitary surgery all of these
    # metabolic abnormalities were abolished," with no extractable
    # post-op glucose number -- this preset represents ACTIVE,
    # untreated acromegaly only, matching the citation's own
    # pre-operative measurement.
    #
    # Relationship to Drug.GROWTH_HORMONE (row 1.35) decided at
    # implementation time, not left ambiguous: chronic acromegaly
    # (this, Gb-only) and acute exogenous GH infusion (that, p3-only)
    # are almost certainly the same underlying mechanism at a
    # different exposure duration, but touch different
    # PathwayParameters fields by construction -- composing both is
    # clinically paradoxical but not a double-count of any single
    # number, so no mutual exclusion was added in profile.py.
    #
    # Composition-boundary check performed before shipping (G1
    # "survives composition" criterion): this Gb-raising factor
    # (5.673 mM at this platform's own default) was checked against
    # MIN_GB_THRESHOLD_MARGIN_MM (personalizer.py) composed with the
    # platform's own strongest Gb-LOWERING mechanisms. Tirzepatide
    # at its max dose (15mg) already clamps to the floor with ZERO
    # acromegaly present (raw 2.98 mM vs. floor 4.8 mM) -- adding
    # acromegaly's own raising factor cannot make this worse (raw
    # composed value moves up, to ~3.38 mM, still clamped, same
    # already-existing clamp_warnings fires, not a new finding).
    # Metformin at its max chronic dose (2000mg) composed with
    # acromegaly lands at a raw ~4.84 mM -- just ABOVE the 4.8 mM
    # floor, not clamped, but within ~0.04 mM of it -- confirmed by
    # direct computation, not asserted; pinned as a regression test
    # below so a future change to either factor that pushes this
    # combination below the floor is caught rather than silently
    # accepted.
    return {"Gb": 5.9 / 5.2}  # -> 5.0*1.1346=5.673 mM


def _osa_modifiers() -> dict[str, float]:
    # Obstructive sleep apnoea as a CHRONIC clinical diagnosis (Track 1
    # row 1.39) -- deliberately distinct from
    # intermittent_hypoxia_modifiers() above, which represents an
    # ACUTE, wakefulness-only exposure rate at a specific measured
    # desaturation frequency (Louis & Punjabi 2009). This preset is for
    # a diagnosed OSA patient with no measured rate on hand -- the two
    # are mutually exclusive at the profile.py level (see DiseaseState
    # .OSA's own enum comment and PhysiologicalProfile.__post_init__)
    # specifically because both lower p3 via the same real mechanism
    # and would silently double-count if composed.
    #
    # Re-verified live via PubMedFetcher before writing this (not
    # trusted from this row's own earlier scoping note): **Harsch IA,
    # Schahin SP, Radespiel-Troger M, et al. "Continuous positive
    # airway pressure treatment rapidly improves insulin sensitivity
    # in patients with obstructive sleep apnea syndrome." Am J Respir
    # Crit Care Med. 2004;169(2):156-62. PMID:14512265.** n=40
    # (AHI>20), hyperinsulinemic euglycemic clamp before, 2 days after,
    # and 3 months after CPAP. Confirmed verbatim: "Insulin sensitivity
    # significantly increased after 2 days (5.75 +/- 4.20 baseline
    # versus 6.79 +/- 4.91 micromol/kg.min; p = 0.003) and remained
    # stable after 3 months of treatment" -- untreated (baseline,
    # pre-CPAP) p3 ratio = 5.75/6.79 = 0.8469.
    #
    # A GENUINE, UNRESOLVED literature conflict is disclosed, not
    # resolved by fiat (same discipline as PCOS/hypothyroidism above):
    # **Smurra M, Philip P, Taillard J, et al. "CPAP treatment does not
    # affect glucose-insulin metabolism in sleep apneic patients." Sleep
    # Med. 2001;2(3):207-13. PMID:11311683.** Re-verified live: n=6
    # OSAS clamp cohort found "no change in mean glycemia, insulin, C
    # peptide... No difference in the amount of glucose infused during
    # the duration of the clamp was noted either" after CPAP --
    # directly contradicting Harsch. Both studies' own conclusions
    # explicitly name OBESITY as a likely confound/effect-modifier
    # (Harsch: "the effect of CPAP on insulin sensitivity is smaller in
    # obese patients... suggesting that in obese individuals insulin
    # sensitivity is mainly determined by obesity"; Smurra: "Obesity,
    # when present, is the important variable") -- this platform
    # composes OSA's factor with `bmi_modifiers()` regardless (no
    # obesity-gating mechanism exists for any disease state, same
    # "disclosed, not solved" precedent as PCOS's own obesity confound
    # above), so a lean-OSA and obese-OSA profile receive the identical
    # p3 factor here despite Harsch's own data suggesting the true
    # effect is smaller in the obese case specifically.
    #
    # ke/Gb/secretion_scale not touched: Harsch's own clamp measures
    # only insulin-mediated glucose disposal (p3); no fasting glucose
    # or secretion-capacity number is reported in either source.
    # GOLD tier -- a ratio computed directly from Harsch's own two
    # reported clamp numbers, no unit conversion, matching this file's
    # established GOLD definition.
    return {"p3": 5.75 / 6.79}


def _nafld_modifiers() -> dict[str, float]:
    # Non-alcoholic fatty liver disease (Track 1 row 1.45) -- a p3
    # (peripheral insulin action) factor only.
    #
    # Verified live via PubMedFetcher (found by a fresh Track 1 scoping
    # fan-out, then independently re-fetched, not trusted from the
    # agent's report): **Bugianesi E, Gastaldelli A, Vanni E, Gambino R,
    # Cassader M, Baldi S, Ponti V, Pagano G, Ferrannini E, Rizzetto M.
    # "Insulin resistance in non-diabetic patients with non-alcoholic
    # fatty liver disease: sites and mechanisms." Diabetologia. 2005.
    # PMID:15747110.** Two-step euglycaemic clamp (1.5 and 6 pmol/min/kg)
    # + [6,6-2H2]glucose/[2H5]glycerol tracers in n=12 normolipidaemic,
    # normotensive, NON-DIABETIC patients with biopsy-proven NAFLD
    # (described as non-obese; per a reviewer's reading of the paper's
    # table, mean BMI ~26.4 -- in this platform's own "overweight" band
    # -- BMI-matched to controls at ~25.8) vs. n=6 controls. Confirmed
    # verbatim: "endogenous glucose production (basal and during the
    # clamp) was normal; however, peripheral glucose disposal was
    # markedly decreased (by 30% and 45% at the low and high insulin
    # doses, respectively, p<0.0001) at higher plasma insulin levels
    # (p=0.05)".
    #
    # ANCHOR CHOICE -- a disclosed JUDGMENT CALL the source does not
    # settle, NOT a demonstrably conservative or more-physiological
    # pick (a first draft of this comment claimed both; two independent
    # reviewers found the claim unsupported -- corrected here): 0.70 =
    # the LOW-dose step's 30% deficit; the 45% high-dose step (0.55) is
    # not used. Read from the paper's own tables by a security reviewer
    # (NOT independently re-verified by the author of this comment --
    # treat as indicative): insulin was not higher in NAFLD at the low
    # step, and normalising the disposal increment by the insulin
    # increment gives a sensitivity ratio anywhere from ~1.0 (low step)
    # to ~0.37 (high step). 0.70 is a defensible midpoint of the two
    # reported deficits, nothing more; insulin-normalised sensitivity
    # is unresolved by this source.
    #
    # EGP/Gb/ke/secretion_scale deliberately NOT touched: the abstract
    # directly measures endogenous glucose production as NORMAL (so no
    # hepatic_egp_delta double-counting -- same "leave what the
    # citation says is unaffected" discipline as CKD/GCK-MODY/PCOS), and
    # reports no magnitude for the others. DISCLOSED gap: per the same
    # reviewer's reading of the full text, it also reports a higher
    # hepatic insulin-resistance index (p=0.05) and a non-significant
    # upward shift of the EGP-suppression curve (p=0.2), and its
    # Discussion describes co-existing hepatic AND peripheral insulin
    # resistance -- real hepatic IR is NOT modelled here (no extractable
    # EGP magnitude), only the peripheral deficit.
    #
    # COMPOSITION, disclosed not blocked (same precedent as PCOS/OSA):
    # (a) Obesity: the cohort is BMI-matched to its controls, so 0.70
    # is an increment over a comparator of similar BMI -- composing with
    # bmi_modifiers() is consistent for overweight profiles, but
    # NAFLD+obese (0.70 x the obesity factor) extrapolates beyond what
    # was measured. (b) Other insulin-resistance states: the cohort
    # excluded diabetes and the metabolic syndrome, and NAFLD is itself
    # a hepatic manifestation of the metabolic syndrome -- composing it
    # with INSULIN_RESISTANT (0.70, uncited, same magnitude) or
    # TYPE_2_DIABETES (0.55; fatty liver is very common in T2D) is
    # probably PARTIAL double-counting of one underlying
    # insulin-resistant state (0.49 / 0.385), not independent effects.
    # Also legal with T1D/GDM, populations the source did not study.
    # (c) NEFA: the same abstract reports higher basal glycerol
    # appearance/lipid oxidation, less suppressed by insulin -- what
    # NEFAParameters.for_obesity() encodes. Not an automatic double-
    # count today (ScenarioDefinition.model_nefa_dynamics defaults to
    # False; runner.py builds a bare, basal-neutral NEFAPathway()), but
    # opting into NEFA dynamics AND NAFLD together describes adipose
    # insulin resistance twice.
    return {"p3": 1.0 - 0.30}


def _single_disease_modifier(state: "DiseaseState") -> dict[str, float]:
    """
    One state's modifier dict, via _disease_modifier_table(). Accepts a
    DiseaseState or its raw string value, like the former `==` chain did
    (DiseaseState is a str Enum); anything else contributes no modifier.
    """
    from ..profile import DiseaseState

    try:
        key = DiseaseState(state)
    except ValueError:
        return {}
    modifier = _disease_modifier_table().get(key)
    return modifier() if modifier is not None else {}


@functools.cache
def _disease_modifier_table() -> "dict[DiseaseState, Callable[[], dict[str, float]]]":
    """Built on first use: profile.py is only imported lazily from this module."""
    from ..profile import DiseaseState

    return {
        DiseaseState.HEALTHY: _healthy_modifiers,
        DiseaseState.INSULIN_RESISTANT: _insulin_resistant_modifiers,
        DiseaseState.TYPE_2_DIABETES: _type_2_diabetes_modifiers,
        DiseaseState.TYPE_1_DIABETES: _type_1_diabetes_modifiers,
        DiseaseState.HYPOTHYROIDISM: _hypothyroidism_modifiers,
        DiseaseState.HYPERTHYROIDISM: _hyperthyroidism_modifiers,
        DiseaseState.GESTATIONAL_DIABETES: _gestational_diabetes_modifiers,
        DiseaseState.GCK_MODY: _gck_mody_modifiers,
        DiseaseState.CHRONIC_KIDNEY_DISEASE: _chronic_kidney_disease_modifiers,
        DiseaseState.PCOS: _pcos_modifiers,
        DiseaseState.ACROMEGALY: _acromegaly_modifiers,
        DiseaseState.OSA: _osa_modifiers,
        DiseaseState.NAFLD: _nafld_modifiers,
    }
