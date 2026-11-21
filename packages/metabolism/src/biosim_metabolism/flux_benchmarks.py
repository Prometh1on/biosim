# SAVE TO: packages\metabolism\src\biosim_metabolism\flux_benchmarks.py
"""
Real, literature-anchored flux-type benchmarks for Track 4 (MASTER_PLAN
4.2 -- "the load-bearing prerequisite for the entire track").

All 6 of Biosim's registered clinical benchmarks (`biosim_validation`)
are concentration/timing measures (e.g. peak glucose in mM, time-to-
peak in minutes) -- there was previously no way to validate any Flux
Variability Analysis (FVA)-derived output at all, since FVA produces
flux RANGES (e.g. umol/kg/min of glucose production), not
concentrations. This module holds real, primary-source, tracer-measured
flux values with a defensible, documented unit-conversion path into this
platform's own mM/min convention, for future use once a real
stoichiometric hepatic network exists (MASTER_PLAN 4.3b) to run
flux_variability_analysis() against and compare (4.3c).

Literature search performed via live NCBI E-utilities across four
distinct search angles: maximal/near-maximal EGP under hypoglycemia or
glucagon stimulation; tracer-based gluconeogenesis-vs-glycogenolysis
partition; 13C-MRS hepatic glycogen turnover; indirect-calorimetry
respiratory-quotient data. No single clean "maximal EGP rate" number was
found in an abstract-extractable form -- the single most promising lead
(Matsuda/DeFronzo 2002, PMID:12200754, a graded glucagon dose-response
study up to a supraphysiological 8 ng/kg/min) is not open-access, and
its abstract reports only qualitative/statistical-significance language,
not the absolute numbers its results tables almost certainly contain.
Documented honestly rather than forced -- see CLAUDE.md's own "DONE --
Acquired a real flux-type validation benchmark" entry for the full
search trail and every candidate considered, not just the two adopted
below.

What WAS found and independently re-verified live (not trusted from a
single search pass) is a real, controlled, tracer-based BASAL/near-
basal EGP measurement with clean, directly-quotable absolute numbers in
a convertible unit -- used here as the primary flux-type benchmark
anchor, plus a second, separate, stimulated-increment reference kept
deliberately un-combined with the first (different study, different
population, no shared basal reference point -- combining them would
repeat the exact "different conditions" confound this project's own
semaglutide-benchmark work already caught and rejected once).

2026-08-12 deep-search plan, item F1: also acquired a real, fasted-
state de novo lipogenesis (DNL) benchmark -- Rosqvist et al. 2019,
PMID:30721918, a genuine 2H2O-labeling fractional-DNL measurement in
VLDL-TG palmitate. This is a FRACTIONAL contribution, not a flux rate
-- represented via a new, separate FractionalContributionBenchmark
class rather than forced through FluxBenchmark's own glucose-specific
unit conversion. Per MASTER_PLAN.md row 4.3d, this benchmark's
existence is the gate for extending the hand-curated hepatic network
with PPP/lipogenesis reactions -- acquiring it does not itself build
any such reactions (none exist yet).
"""

from __future__ import annotations

from dataclasses import dataclass

from biosim_ode.meal import DEFAULT_BODY_WEIGHT_KG, DEFAULT_VG_L_PER_KG


def umol_per_kg_min_to_mM_per_min(
    value_umol_per_kg_min: float,
    body_weight_kg: float = DEFAULT_BODY_WEIGHT_KG,
    vg_l_per_kg: float = DEFAULT_VG_L_PER_KG,
) -> float:
    """
    Convert a whole-body-weight-normalized flux rate (umol per kg body
    weight per min) into this platform's own concentration-rate
    convention (mM/min), via the SAME glucose distribution volume
    already established and consolidated elsewhere in this codebase
    (`MealInput.DEFAULT_VG_L_PER_KG` == `BERGMAN_IVGTT_VG_L_PER_KG` ==
    0.230 L/kg, Mora-Gonzalez et al. 2024, PMID:38389141). Extracted as
    a public, standalone function (rather than kept as
    `FluxBenchmark`'s own private method) so any other flux-rate value
    in this package -- not just a `FluxBenchmark` instance's own fields
    -- can be converted with the exact same formula, never a
    re-implemented copy that could silently drift.

    rate [umol/min] = value_umol_per_kg_min * body_weight_kg
    rate [mmol/min] = rate [umol/min] / 1000
    Vg   [L]        = vg_l_per_kg * body_weight_kg
    rate [mM/min]   = rate [mmol/min] / Vg [L]
    """
    rate_umol_per_min = value_umol_per_kg_min * body_weight_kg
    rate_mmol_per_min = rate_umol_per_min / 1000.0
    vg_l = vg_l_per_kg * body_weight_kg
    return rate_mmol_per_min / vg_l


@dataclass(frozen=True)
class FluxBenchmark:
    """
    A real, literature-measured whole-body glucose flux rate, with a
    defensible, documented conversion path into this platform's own
    concentration-rate convention (mM/min).
    """

    name: str
    pubmed_id: str
    value_umol_per_kg_min: float
    stderr_umol_per_kg_min: float
    population: str
    condition: str
    caveat: str = ""

    def value_mM_per_min(
        self,
        body_weight_kg: float = DEFAULT_BODY_WEIGHT_KG,
        vg_l_per_kg: float = DEFAULT_VG_L_PER_KG,
    ) -> float:
        """
        Convert this whole-body-weight-normalized flux rate (umol per kg
        body weight per min) into this platform's own concentration-rate
        convention (mM/min), via the SAME glucose distribution volume
        already established and consolidated elsewhere in this codebase
        (`MealInput.DEFAULT_VG_L_PER_KG` == `BERGMAN_IVGTT_VG_L_PER_KG`
        == 0.230 L/kg, Mora-Gonzalez et al. 2024, PMID:38389141) --
        reused directly rather than re-derived, so this conversion
        cannot silently drift out of sync with the rest of the
        platform's own distribution-volume convention.

        rate [umol/min] = value_umol_per_kg_min * body_weight_kg
        rate [mmol/min] = rate [umol/min] / 1000
        Vg   [L]        = vg_l_per_kg * body_weight_kg
        rate [mM/min]   = rate [mmol/min] / Vg [L]
        """
        return umol_per_kg_min_to_mM_per_min(self.value_umol_per_kg_min, body_weight_kg, vg_l_per_kg)

    def stderr_mM_per_min(
        self,
        body_weight_kg: float = DEFAULT_BODY_WEIGHT_KG,
        vg_l_per_kg: float = DEFAULT_VG_L_PER_KG,
    ) -> float:
        """Same conversion applied to the reported standard error."""
        return umol_per_kg_min_to_mM_per_min(self.stderr_umol_per_kg_min, body_weight_kg, vg_l_per_kg)


# Adkins A, Basu R, Persson M, Dicke B, Shah P, Vella A, Schwenk WF,
# Rizza R. "Higher insulin concentrations are required to suppress
# gluconeogenesis than glycogenolysis in nondiabetic humans." Diabetes.
# 2003;52(9):2213-2220. PMID:12941759 -- verified live via
# PubMedFetcher.fetch_by_pubmed_id() this session, abstract quoted
# directly, not assumed from a secondary search summary. Graded insulin
# infusion (0.25/0.375/0.5 mU/kg/min) at clamped euglycemia
# (~5.5 mmol/L) in healthy nondiabetic humans; EGP measured by
# [3-3H]glucose tracer. The LOWEST insulin infusion (~75 pmol/L, close
# to a physiological fasting/basal insulin level) gives the closest-to-
# basal EGP value this study reports.
ADKINS_2003_BASAL_EGP = FluxBenchmark(
    name="Adkins et al. 2003 near-basal EGP",
    pubmed_id="12941759",
    value_umol_per_kg_min=13.1,
    stderr_umol_per_kg_min=1.3,
    population="healthy nondiabetic adults",
    condition="euglycemic clamp (~5.5 mmol/L), lowest insulin infusion (~75 pmol/L, near-basal)",
    caveat=(
        "This is a near-BASAL rate (the lowest of three graded insulin "
        "infusions this study tested), not a maximal/stress-stimulated "
        "ceiling -- EGP_MAX_MM_MIN represents a maximum, so this "
        "benchmark alone can only sanity-check that a hepatic network's "
        "FVA-derived range plausibly contains this basal value, not "
        "that the range's own upper bound is correctly anchored."
    ),
)

# Drouin R, Lavoie C, Bourque J, Ducros F, Poisson D, Chiasson JL.
# "Increased hepatic glucose production response to glucagon in trained
# subjects." Am J Physiol. 1998;274(1):E23-E28. PMID:9458743 -- verified
# live via PubMedFetcher this session. Somatostatin pancreatic clamp +
# physiological glucagon infusion (1.5 ng/kg/min) in healthy men; EGP
# measured by tracer. Reports the INCREASE in HGP area-under-curve above
# each subject's own basal, not a total absolute rate -- kept
# deliberately separate from ADKINS_2003_BASAL_EGP (a different study,
# different population, no shared basal reference to add them together
# against) rather than combined into one derived "total" number.
# Untrained-subject value used (closer to a general, non-athlete
# population than this study's own trained-subject arm).
DROUIN_1998_GLUCAGON_STIMULATED_INCREMENT = FluxBenchmark(
    name="Drouin et al. 1998 glucagon-stimulated EGP increment (untrained, first hour)",
    pubmed_id="9458743",
    value_umol_per_kg_min=7.4,
    stderr_umol_per_kg_min=1.6,
    population="healthy untrained men",
    condition="somatostatin pancreatic clamp + physiological glucagon infusion (1.5 ng/kg/min), first hour",
    caveat=(
        "This is an INCREMENT above each subject's own basal EGP, not a "
        "total rate -- do not add directly to ADKINS_2003_BASAL_EGP's "
        "value (different study, different population, no shared basal "
        "reference point to combine against). Qualitative/directional "
        "corroboration only: confirms EGP can rise substantially above "
        "basal under a real physiological stimulus, supporting the "
        "existence of real headroom between EGP_BASAL_MM_MIN and "
        "EGP_MAX_MM_MIN, without directly quantifying EGP_MAX_MM_MIN "
        "itself."
    ),
)


@dataclass(frozen=True)
class FractionalContributionBenchmark:
    """
    A real, literature-measured FRACTIONAL contribution of one metabolic
    pathway to a larger pool (e.g. de novo lipogenesis's own fractional
    contribution to VLDL-triglyceride palmitate) -- a genuinely different
    physical quantity from FluxBenchmark's absolute whole-body flux rate
    (umol/kg/min), so it is NOT expressed via that class's own glucose-
    distribution-volume conversion (a fraction has no such conversion at
    all -- forcing it through that machinery would be a category error,
    not just an inconvenience).
    """

    name: str
    pubmed_id: str
    fraction: float          # 0-1, e.g. 0.065 for a 6.5% contribution
    fraction_low: float      # lower bound of the reported spread (IQR/range)
    fraction_high: float     # upper bound of the reported spread (IQR/range)
    population: str
    condition: str
    method: str
    caveat: str = ""


# Rosqvist F, McNeil CA, Pramfalk C, Parry SA, Low WS, Cornfield T,
# Fielding BA, Hodson L. "Fasting hepatic de novo lipogenesis is not
# reliably assessed using circulating fatty acid markers." Am J Clin
# Nutr. 2019 -- verified live via PubMedFetcher (full text, PMC) this
# session, not assumed from a secondary summary. Pooled cross-sectional
# analysis of 5 cohorts, n=149 (64 male/85 female), nondiabetic, no
# lipid/glucose-affecting medication. Method: 2H2O (deuterated water)
# labeling -- subjects consumed 2H2O the evening before, maintained
# plasma-water enrichment ~0.3%, then AFTER AN OVERNIGHT FAST a fasting
# blood sample was drawn and deuterium incorporation into VLDL-TG
# palmitate was measured by GC-MS.
#
# CORRECTION (found by a Track 4.3d scoping pass, computational-
# biologist agent, independently re-verified via a second live full-
# text fetch before accepting): this is NOT the "genuinely FASTED-state
# measurement" this comment previously claimed. Confirmed directly from
# the paper's own Methods section (quoted verbatim): "The evening
# before the study day (~12 h before assessment), subjects consumed
# 2H2O... over the course of the evening... Subjects were advised to
# consume a low-fat dinner the evening before measurement of DNL.
# After an overnight fast, subjects came to the Clinical Research Unit,
# and a blood sample was taken." The label is administered THE EVENING
# BEFORE, overlapping with dinner, and isotope incorporation into
# VLDL-TG palmitate accumulates continuously over that ~12h window --
# only the final BLOOD DRAW happens after an overnight fast, not the
# synthesis being measured. "Fasting hepatic DNL" describes when the
# sample was taken, not when the labeled synthesis occurred -- a real,
# consequential distinction for using this as a validation target
# against `fasted_hepatic_network.py`'s own specifically fasted-state
# scope (see that scoping pass's own conclusion: MASTER_PLAN row 4.3d
# was NOT built as originally planned because of this finding, among
# others -- see the row's own updated text).
#
# Quoted directly: "the fasting hepatic DNL was 6.5% (range 0.1-38.7%)"
# for the whole cohort (median, IQR 3.8-11.5%); a n=21-male subgroup
# measured across multiple lipid fractions gave "DNL in VLDL-TG...
# 13.7 +/- 7.9% (mean +/- SD)". The whole-cohort median (6.5%) is used
# here as the primary value -- broader-population, larger n -- with the
# male-subgroup mean noted in this benchmark's own caveat for context,
# not combined into one derived number (different sub-sample, different
# summary statistic type, matching this module's own established
# "don't merge different-condition numbers" convention for
# DROUIN_1998_GLUCAGON_STIMULATED_INCREMENT above).
ROSQVIST_2019_FASTED_HEPATIC_DNL = FractionalContributionBenchmark(
    name="Rosqvist et al. 2019 fasted hepatic de novo lipogenesis (VLDL-TG palmitate)",
    pubmed_id="30721918",
    fraction=0.065,
    fraction_low=0.038,
    fraction_high=0.115,
    population="nondiabetic adults, pooled 5-cohort cross-sectional analysis (n=149, BMI 19.5-39.5)",
    condition="fasted (overnight fast), deuterated-water (2H2O) labeling, VLDL-TG palmitate GC-MS",
    method="2H2O deuterated water labeling",
    caveat=(
        "This is a pooled, disease-free-but-not-narrowly-lean cohort "
        "(BMI up to 39.5), not a single tightly homogeneous healthy-lean "
        "sample -- the 6.5% median (IQR 3.8-11.5%) is used as the "
        "primary value here; a n=21-male subgroup in the same paper "
        "reported a higher 13.7 +/- 7.9% (mean +/- SD) for VLDL-TG DNL "
        "specifically, kept separate rather than combined (different "
        "sub-sample, different summary-statistic type). This is a "
        "FRACTIONAL contribution, not a flux rate -- do not attempt to "
        "convert it via FluxBenchmark's own umol_per_kg_min_to_mM_per_min "
        "machinery, which assumes a different physical quantity "
        "entirely. NOT a purely fasted-state measurement despite its own "
        "name: the 2H2O label is dosed the evening before, overlapping "
        "with a low-fat dinner, and isotope incorporation into VLDL-TG "
        "palmitate accumulates over that whole ~12h window -- only the "
        "final blood draw happens after an overnight fast, not the "
        "labeled synthesis itself (verified directly from the paper's "
        "own Methods section). Do not treat this as a clean validation "
        "target for a specifically fasted-state network without "
        "accounting for that. Per MASTER_PLAN.md row 4.3d, this "
        "benchmark's own existence does NOT license building PPP/DNL "
        "reactions on `fasted_hepatic_network.py` as originally scoped -- "
        "see that row's own updated text for why."
    ),
)
