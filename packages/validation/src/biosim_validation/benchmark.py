# SAVE TO: packages\validation\src\biosim_validation\benchmark.py
"""
Clinical validation benchmarks.

These are published reference values from well-replicated clinical studies.
The model is considered validated on a benchmark when its output falls within
the defined tolerance band.

A key principle: validation benchmarks are defined BEFORE running validation.
Post-hoc selection of benchmarks that the model happens to pass is circular.
All benchmarks here are pre-registered against their literature sources.

Benchmark sources:
    OGTT_75G_HEALTHY:  Van Cauter et al. 1992, Am J Physiol. PMID:1566835
                       Numerous meta-analyses of healthy adult OGTT responses.

    BERGMAN_MINIMAL:   Bergman et al. 1981, J Clin Invest. PMID:7033284
                       Original minimal model validation dataset.

Validation tolerance policy:
    TIGHT:  model must be within 10% of reference value
    NORMAL: model must be within 25% of reference value
    WIDE:   model must be within 50% of reference value (for poorly constrained outputs)
    STRICT: model must be within 20% of reference value (headline band of the
            held-out group, heldout.py; a miss gives MISS, like WIDE)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional


class Tolerance(str, Enum):
    TIGHT  = "tight"   # within 10% of reference
    NORMAL = "normal"  # within 25% of reference
    WIDE   = "wide"    # within 50% of reference
    STRICT = "strict"  # within 20% of reference (held-out group)


@dataclass(frozen=True)
class BenchmarkPoint:
    """
    A single clinical reference value with its tolerance band.

    name:         Human-readable name for reporting
    description:  What this benchmark tests
    reference:    The expected (literature) value
    unit:         Units of the reference value (for display)
    tolerance:    How close the model must be
    pubmed_id:    Primary literature source
    population:   Population the reference applies to
    """
    name:        str
    description: str
    reference:   float
    unit:        str
    tolerance:   Tolerance
    pubmed_id:   str
    population:  str = "healthy_adults"
    notes:       str = ""
    # Gold-standard grading fields (S0, 2026-11-20; docs/GOLD_STANDARD.md).
    # All optional and filled only from fetched text. They never change the
    # PASS/MISS verdict, only the separate grade (grading.py).
    spread_kind: str = "none"            # "ci_diff" | "sd" | "sem" | "figure" | "none"
    n:           Optional[int] = None
    sd:          Optional[float] = None   # between-subject SD of the reference
    sem:         Optional[float] = None   # SEM of the reference (SD = SEM*sqrt(n))
    ci_low:      Optional[float] = None   # 95% CI of a treatment difference
    ci_high:     Optional[float] = None
    sampling_interval_min: Optional[float] = None
    quote:       str = ""                 # verbatim source text for the spread
    source_location: str = ""             # e.g. "abstract", "Table 2"
    benchmark_set: str = "dev"            # "dev" (may guide changes) | "locked"

    @property
    def is_scored(self) -> bool:
        """
        False for a point whose reference value has no real literature
        anchor (pubmed_id == "UNVERIFIED") -- e.g. mtor_activation, whose
        own notes admit "no direct clinical comparator exists". Checking
        the model against an unscored point's reference is not a genuine
        validation, however the pass/fail comes out -- see
        ModelValidator._unscored_point_notes(), which flags this
        explicitly in every report rather than let it read as one of "n"
        real checks passed.
        """
        return self.pubmed_id != "UNVERIFIED"

    @property
    def tolerance_fraction(self) -> float:
        return {"tight": 0.10, "normal": 0.25, "wide": 0.50, "strict": 0.20}[self.tolerance.value]

    @property
    def lower_bound(self) -> float:
        return self.reference * (1.0 - self.tolerance_fraction)

    @property
    def upper_bound(self) -> float:
        return self.reference * (1.0 + self.tolerance_fraction)

    def check(self, model_value: float) -> "BenchmarkCheckResult":
        passed = self.lower_bound <= model_value <= self.upper_bound
        deviation = (model_value - self.reference) / self.reference * 100.0
        return BenchmarkCheckResult(
            benchmark=self,
            model_value=model_value,
            passed=passed,
            deviation_percent=deviation,
        )


@dataclass(frozen=True)
class BenchmarkCheckResult:
    """Result of checking one benchmark point."""
    benchmark:         BenchmarkPoint
    model_value:       float
    passed:            bool
    deviation_percent: float   # positive = model overshoots

    @property
    def name(self) -> str:
        return self.benchmark.name

    @property
    def status(self) -> str:
        """PASS if the model value is inside the tolerance band, MISS if not."""
        return "PASS" if self.passed else "MISS"


@dataclass
class BenchmarkDataset:
    """A named collection of benchmark points for one validation scenario."""
    name:        str
    description: str
    points:      list[BenchmarkPoint] = field(default_factory=list)

    def add(self, point: BenchmarkPoint) -> None:
        self.points.append(point)

    def __len__(self) -> int:
        return len(self.points)


# ---------------------------------------------------------------------------
# Pre-registered benchmark datasets
# ---------------------------------------------------------------------------

def ogtt_75g_healthy_adults() -> BenchmarkDataset:
    """
    75g OGTT reference values for healthy adults (non-diabetic, BMI 18.5-25,
    age 20-50, no medications affecting glucose metabolism).

    All glucose values in mM, insulin in pM, time in minutes.
    """
    ds = BenchmarkDataset(
        name="OGTT_75g_healthy",
        description="75g OGTT in healthy adults — plasma glucose and insulin dynamics",
    )

    # Fasting glucose
    ds.add(BenchmarkPoint(
        name="fasting_glucose",
        description="Fasting blood glucose before meal",
        reference=5.0, unit="mM", tolerance=Tolerance.TIGHT,
        pubmed_id="1566835",
        notes="Population mean ~4.9-5.2 mM. Tight tolerance as this is the initial condition.",
    ))

    # Peak glucose
    ds.add(BenchmarkPoint(
        name="peak_glucose",
        description="Peak blood glucose after 75g glucose load",
        reference=7.8, unit="mM", tolerance=Tolerance.NORMAL,
        pubmed_id="1566835",
        notes="Population range ~6.5-9.5 mM. Wide inter-individual variation.",
    ))

    # Time to peak glucose
    ds.add(BenchmarkPoint(
        name="time_to_peak_glucose_min",
        description="Time of peak glucose after meal (minutes)",
        reference=45.0, unit="min", tolerance=Tolerance.WIDE,
        pubmed_id="1566835",
        notes="Typically 30-60 min. Highly variable.",
    ))

    # Glucose at 120 min
    ds.add(BenchmarkPoint(
        name="glucose_at_120min",
        description="Blood glucose at 2 hours post-meal (WHO OGTT criterion)",
        reference=6.5, unit="mM", tolerance=Tolerance.NORMAL,
        # Was cited to DeFronzo 1979 (PMID:382871), the glucose-clamp methods
        # paper, which reports no OGTT value (B4, 2026-11-16). No source for
        # 6.5 mM was found, so the point is UNVERIFIED and unscored.
        pubmed_id="UNVERIFIED",
        notes="No source: 6.5 mM is a chosen healthy target below the WHO "
              "impaired-glucose-tolerance threshold (7.8 mM). "
              "WHO threshold for impaired glucose tolerance is 7.8 mM. "
              "Healthy adults should be below 7.8 mM. Target ~6-7 mM.",
    ))

    # Peak insulin
    ds.add(BenchmarkPoint(
        name="peak_insulin",
        description="Peak plasma insulin after 75g glucose load",
        reference=350.0, unit="pM", tolerance=Tolerance.WIDE,
        pubmed_id="1566835",
        notes="Wide range ~200-600 pM across healthy subjects. Wide tolerance appropriate.",
    ))

    # Time to peak insulin
    ds.add(BenchmarkPoint(
        name="time_to_peak_insulin_min",
        description="Time of peak insulin (minutes)",
        reference=45.0, unit="min", tolerance=Tolerance.WIDE,
        pubmed_id="1566835",
        notes="Insulin peaks shortly after glucose. Typically 30-60 min.",
    ))

    # mTOR at peak
    ds.add(BenchmarkPoint(
        name="mtor_activation",
        description="mTOR fractional activity at peak (post-meal)",
        reference=0.55, unit="fraction", tolerance=Tolerance.WIDE,
        pubmed_id="UNVERIFIED",
        notes="mTOR activation is estimated from insulin signalling magnitude. "
              "No direct clinical comparator exists for this derived quantity -- "
              "left explicitly unverified rather than attached to a plausible-"
              "looking but unrelated PMID (see CLAUDE.md KNOWN BUGS FIXED #9/#13). "
              "Wide tolerance is appropriate. Searched PubMed 2026-11-14: "
              "human muscle studies report phospho/total fold-changes, not "
              "an absolute 0-1 activity; the best meal time course found "
              "(Gran 2014, PMID:25302072, biopsies at 0/2/4 h) gives "
              "p70S6K only in figures and used a 31 g protein meal, so "
              "leucine, not insulin alone, drives it. No usable number.",
    ))

    return ds


def bergman_iv_glucose() -> BenchmarkDataset:
    """
    Bergman 1981 IV glucose bolus (20g in 70kg subject).
    Intravenous — faster dynamics than oral.
    """
    ds = BenchmarkDataset(
        name="Bergman_IV_glucose",
        description="20g IV glucose bolus — Bergman 1981 minimal model",
    )
    ds.add(BenchmarkPoint(
        name="glucose_return_time_min",
        description="Time for glucose to return within 10% of baseline (IV bolus)",
        reference=90.0, unit="min", tolerance=Tolerance.WIDE,
        pubmed_id="7033284",
        notes="Bergman 1981 Fig 1: glucose returns to baseline by ~90-120 min.",
        spread_kind="figure",   # a reading from a figure, not a number in the text
    ))
    return ds


def semaglutide_postprandial_glucose() -> BenchmarkDataset:
    """
    Semaglutide (1.0mg/week, steady state) vs. placebo, same meal.

    First validation benchmark in this codebase for a clinically
    significant drug pathway (GLP-1 agonists) rather than the baseline
    glucose/insulin model alone -- see CLAUDE.md "DONE -- Built the first validation benchmark for a clinically significant drug pathway".
    """
    ds = BenchmarkDataset(
        name="Semaglutide_1mg_postprandial",
        description=(
            "Semaglutide 1.0mg/week steady-state vs. placebo, standard meal -- "
            "Hjerpsted et al. 2018"
        ),
    )
    ds.add(BenchmarkPoint(
        name="glucose_iauc_5h_reduction_mmol_h_per_L",
        description=(
            "Reduction in postprandial glucose incremental AUC (0-5h) with "
            "semaglutide vs. placebo, same meal (a positive value = semaglutide "
            "lowers iAUC; magnitude of the reduction, not the signed difference)"
        ),
        reference=1.34, unit="mmol.h/L", tolerance=Tolerance.WIDE,
        pubmed_id="28941314",
        notes=(
            "Hjerpsted JB, Flint A, Brooks A, Axelsen MB, Kvist T, Blundell J. "
            "\"Semaglutide improves postprandial glucose and lipid metabolism, "
            "and delays first-hour gastric emptying in subjects with obesity.\" "
            "Diabetes Obes Metab. 2018;20(3):610-619. Randomized, double-blind, "
            "placebo-controlled crossover trial, n=30, once-weekly 1.0mg "
            "semaglutide vs. placebo: iAUC(0-5h) estimated treatment difference "
            "-1.34 mmol.h/L (glucose) vs. placebo."
        ),
        # Reduction sign: the CI [-2.42, -0.27] of the difference is a
        # reduction of 0.27 to 2.42.
        spread_kind="ci_diff", ci_low=0.27, ci_high=2.42, n=30,
        quote=("estimated treatment difference: glucose -1.34 mmol h/L "
               "[-2.42, -0.27]"),
        source_location="abstract (PubMed/Europe PMC, fetched 2026-11-20)",
    ))
    return ds


def steroid_t2d_postprandial_glucose() -> BenchmarkDataset:
    """
    Prednisone 40mg vs. no steroid, same meal, in a type 2 diabetes
    population -- Shah et al. 2020, PMID:32369480.

    See CLAUDE.md's steroid/prednisone entry (follows the semaglutide
    validation entry): real evidence shows this effect is concentrated
    in T2D (reduced beta-cell reserve), not healthy subjects (healthy
    beta cells largely compensate, n=37 study PMID:27424036 found no
    significant change at any dose up to 60mg) -- this benchmark
    deliberately targets the T2D population where a real, statistically
    significant effect was actually measured, not the healthy default.
    """
    ds = BenchmarkDataset(
        name="Prednisone_40mg_T2D_postprandial",
        description=(
            "Prednisone 40mg vs. no steroid, same meal, type 2 diabetes -- "
            "Shah et al. 2020"
        ),
    )
    ds.add(BenchmarkPoint(
        name="glucose_auc_percent_increase_t2d",
        description=(
            "Percent increase in postprandial glucose AUC with prednisone "
            "40mg vs. no steroid, same meal, in a T2D population"
        ),
        reference=19.0, unit="%", tolerance=Tolerance.WIDE,
        pubmed_id="32369480",
        notes=(
            "Shah M, Adel MM, Tahsin B, Guerra Y, Fogelfeld L. \"Effect of "
            "short-term prednisone on beta-cell function in subjects with "
            "type 2 diabetes mellitus and healthy subjects.\" PLoS One. "
            "2020. n=5 T2DM subjects, prednisone 40mg/day x3 days, 75g "
            "OGTT performed 2h after the 3rd prednisone dose (own Methods "
            "text, fetched and re-verified against this paper's full text "
            "this session, not just its abstract). "
            "Regression: this used to say 'glucose AUC increased ~50% vs. "
            "baseline (p=0.016 vs. the healthy group's non-significant "
            "~19% increase)' -- WRONG, apparently conflated with this "
            "same paper's HOMA-%B/insulinogenic-index percentile-"
            "difference numbers (50.4%/69.2% T2DM vs. 19%/32.2% healthy, "
            "a DIFFERENT metric). The paper's own Table 1 gives AUC "
            "glucose directly: T2DM 616.8 -> 734.1 mmol/L/min (+19.0%), "
            "healthy 366.5 -> 409.7 mmol/L/min (+11.8%) -- 'the increase "
            "in AUC Glc was significantly higher (p=0.016)' in T2DM than "
            "healthy, i.e. 19.0% vs. 11.8%, not 50% vs. 19%. Small sample "
            "(n=5) -- WIDE tolerance is appropriate, and this is a single "
            "study, not a meta-analysis."
        ),
    ))
    return ds


def metformin_postprandial_glucose_t2d() -> BenchmarkDataset:
    """
    Metformin (14-day chronic dosing) vs. placebo, same 4h mixed meal
    test, in a type 2 diabetes population -- Hansen et al. 2024,
    PMID:39049802.

    Same comparative-benchmark shape as validate_semaglutide()/
    validate_steroid_meal_t2d() -- this is the metformin acute-meal-
    challenge benchmark CLAUDE.md's VALIDATION STATUS section already
    flagged as the explicit next step (a directly comparable oral-meal
    AUC study, not the intraduodenal-infusion tracer studies -- e.g.
    Borg et al. 2025 -- that don't measure the same physical quantity
    this platform's MealInput/Ra(t) represents).
    """
    ds = BenchmarkDataset(
        name="Metformin_T2D_postprandial",
        description=(
            "Metformin (14-day chronic dosing) vs. placebo, same 4h mixed "
            "meal test, type 2 diabetes -- Hansen et al. 2024"
        ),
    )
    ds.add(BenchmarkPoint(
        name="glucose_iauc_4h_reduction_mmol_min_per_L_t2d",
        description=(
            "Reduction in postprandial glucose incremental AUC (0-4h) with "
            "chronic metformin vs. placebo, same meal, in a T2D population "
            "(magnitude of the reduction, not the signed difference)"
        ),
        reference=186.0, unit="mmol.min/L", tolerance=Tolerance.WIDE,
        pubmed_id="39049802", population="type_2_diabetes",
        notes=(
            "Hansen LS, Gasbjerg LS, Bronden A, et al. \"The role of "
            "glucagon-like peptide 1 in the postprandial effects of "
            "metformin in type 2 diabetes: a randomized crossover trial.\" "
            "Eur J Endocrinol. 2024;191(2):192-203. Randomized, "
            "placebo-controlled, double-blind crossover, n=15 T2D (median "
            "HbA1c 50 mmol/mol [6.7%], BMI 30.1, age 71), 14 days metformin "
            "vs. 14 days placebo, each followed by a 4h mixed meal test: "
            "postprandial glucose excursion iAUC reduced by mean-of-"
            "differences 186 mmol/L*min [95% CI 64-307] during saline "
            "infusion (the non-GLP-1R-blocked arm -- this platform doesn't "
            "model exendin(9-39) antagonism, so the saline arm is the "
            "correct comparator, not the antagonist arm's larger 268 "
            "mmol/L*min figure). Metformin dose was titrated 1000-2000mg/day "
            "over the 14 days (per ClinicalTrials.gov NCT03246451: "
            "\"2-4 x 500 mg\"); this benchmark uses 2000mg (the titration "
            "ceiling), matching this platform's existing convention "
            "elsewhere for a standard chronic metformin dose. The reported "
            "95% CI (64-307) is wider on its own than this benchmark's WIDE "
            "tolerance band around the point estimate -- WIDE is still the "
            "right existing category (this project uses a fixed 3-tier "
            "tolerance system, not per-benchmark custom bounds), and the "
            "real CI width is noted here for transparency rather than "
            "implying tighter certainty than the source has. Exact meal "
            "composition wasn't stated in the abstract; this benchmark uses "
            "this platform's standard 75g reference meal (same convention "
            "as the OGTT/semaglutide/steroid benchmarks above), not a "
            "value taken from the paper."
        ),
        spread_kind="ci_diff", ci_low=64.0, ci_high=307.0, n=15,
        quote=("postprandial plasma glucose excursions during both saline "
               "infusion (MD 186 mmol/L × min [95% CI 64-307])"),
        source_location="abstract (PubMed/Europe PMC, fetched 2026-11-20)",
    ))
    return ds


def canagliflozin_postprandial_glucose_and_insulin() -> BenchmarkDataset:
    """
    Canagliflozin 300mg vs. placebo, 20min before a 600-kcal mixed meal,
    healthy subjects -- Polidori et al. 2013, PMID:23412078.

    Same comparative-benchmark shape as validate_semaglutide()/
    validate_metformin(). This is the only same-meal drug-vs-placebo SGLT2i
    dataset independently verified for this platform (see CLAUDE.md's
    SGLT2i pathway milestone) -- but canagliflozin at this dose ALSO
    delays intestinal glucose absorption via SGLT1 co-inhibition
    (confirmed in the source abstract's own text: "canagliflozin reduces
    postprandial plasma glucose and insulin by increasing UGE [via renal
    SGLT2 inhibition] and delaying RaO, likely due to intestinal SGLT1
    inhibition"), which this platform's RenalPathway does NOT model (it
    represents only the renal glucose-clearance mechanism, not intestinal
    SGLT1/delayed gastric emptying). The model is therefore EXPECTED to
    undershoot this benchmark's magnitude -- documented here and in
    ModelValidator.validate_sglt2i(), not silently accepted as a target to
    hit. WIDE tolerance is used per this project's standard 3-tier system,
    not as an attempt to force a pass on a benchmark with a known,
    structural confound.
    """
    ds = BenchmarkDataset(
        name="Canagliflozin_300mg_postprandial",
        description=(
            "Canagliflozin 300mg vs. placebo, 20min before a 600-kcal mixed "
            "meal, healthy subjects -- Polidori et al. 2013"
        ),
    )
    ds.add(BenchmarkPoint(
        name="glucose_auc_2h_percent_reduction",
        description=(
            "Percent reduction in postprandial glucose AUC (0-2h) with "
            "canagliflozin 300mg vs. placebo, same meal (magnitude of the "
            "reduction, not the signed difference). Real value is NOT "
            "purely attributable to the renal mechanism this platform "
            "models -- see this dataset's own docstring."
        ),
        reference=35.0, unit="%", tolerance=Tolerance.WIDE,
        pubmed_id="23412078",
        notes=(
            "Polidori D, Sha S, Mudaliar S, et al. \"Canagliflozin lowers "
            "postprandial glucose and insulin by delaying intestinal "
            "glucose absorption in addition to increasing urinary glucose "
            "excretion.\" Diabetes Care. 2013;36(8):2154-61. Randomized, "
            "double-blind, crossover, n=20 healthy, canagliflozin 300mg vs. "
            "placebo 20min before a 600-kcal mixed meal: glucose AUC(0-2h) "
            "-35%, insulin AUC(0-2h) -43% (both P<0.001), UGE(0-6h) "
            "18.2+/-5.6 vs. <0.2g (P<0.001). Confounded by intestinal SGLT1 "
            "co-inhibition, not purely a renal effect -- expect this "
            "platform's renal-only model to undershoot."
        ),
    ))
    ds.add(BenchmarkPoint(
        name="insulin_auc_2h_percent_reduction",
        description=(
            "Percent reduction in postprandial insulin AUC (0-2h) with "
            "canagliflozin 300mg vs. placebo, same meal. Same SGLT1 "
            "confound caveat as glucose_auc_2h_percent_reduction above."
        ),
        reference=43.0, unit="%", tolerance=Tolerance.WIDE,
        pubmed_id="23412078",
        notes="See glucose_auc_2h_percent_reduction's notes -- same source and caveat.",
    ))
    return ds


PROTEIN_COINGESTION_NOTES = (
    "Thondre PS, Achebe I, Sampson A, et al. \"Co-ingestion of NUTRALYS(R) pea "
    "protein and a high-carbohydrate beverage influences the glycaemic, "
    "insulinaemic, glucose-dependent insulinotropic polypeptide (GIP) and "
    "glucagon-like peptide-1 (GLP-1) responses: preliminary results of a "
    "randomised controlled trial.\" Eur J Nutr. 2021. PMID:33515092 (verified "
    "live 2026-10-20). Randomised crossover, n=31 healthy, 50 g glucose alone "
    "vs. with 25 g (test 1) or 50 g (test 2) pea protein, capillary samples "
    "for 180 min: glucose iAUC180 \"(- 53%)\" test 2 and \"(- 31%)\" test 1 vs. "
    "control; insulin iAUC180 \"(+ 28%)\" test 1 and \"(+ 40%)\" test 2. "
    "Independent of the nutrient-secretion calibration sources (Krezowski "
    "1986, PMID:3538843; Nuttall 1984, PMID:6389060): different lab, protein "
    "source and cohort."
)


def pea_protein_coingestion_50g_glucose() -> BenchmarkDataset:
    """
    50 g glucose with 0 / 25 / 50 g pea protein, healthy subjects --
    Thondre et al. 2021, PMID:33515092.

    Exercises NutrientSecretionPathway (protein-stimulated insulin
    secretion), which had no registered benchmark. The model represents
    ONLY that insulin effect: protein's slowing of gastric emptying is not
    modelled, and pea protein's insulinotropic potency need not match the
    calibration studies' protein sources. (Protein's glucagon response is
    unmodelled too, but it pushes the other way: it would make the model
    over-predict glucose lowering.) The model is therefore EXPECTED to undershoot --
    most sharply for glucose, whose real reduction plausibly comes in
    large part from slower gastric emptying (an inference -- the study did
    not measure emptying). WIDE tolerance per the 3-tier system; the
    point is to make the gap measurable, not to force a pass.
    """
    ds = BenchmarkDataset(
        name="Pea_protein_coingestion_50g_glucose",
        description=(
            "50 g glucose alone vs. with 25 g or 50 g pea protein, healthy "
            "subjects, 180-min iAUCs -- Thondre et al. 2021"
        ),
    )
    for grams, glucose_ref, insulin_ref in ((25, 31.0, 28.0), (50, 53.0, 40.0)):
        ds.add(BenchmarkPoint(
            name=f"glucose_iauc_3h_percent_reduction_{grams}g_protein",
            description=(
                f"Percent reduction in glucose iAUC (0-180 min) with {grams} g "
                "pea protein vs. 50 g glucose alone (magnitude of the reduction)."
            ),
            reference=glucose_ref, unit="%", tolerance=Tolerance.WIDE,
            pubmed_id="33515092", notes=PROTEIN_COINGESTION_NOTES,
        ))
        ds.add(BenchmarkPoint(
            name=f"insulin_iauc_3h_percent_increase_{grams}g_protein",
            description=(
                f"Percent increase in insulin iAUC (0-180 min) with {grams} g "
                "pea protein vs. 50 g glucose alone."
            ),
            reference=insulin_ref, unit="%", tolerance=Tolerance.WIDE,
            pubmed_id="33515092", notes=PROTEIN_COINGESTION_NOTES,
        ))
    return ds


T2D_MEAL_INSULIN_NOTES = (
    "Basu et al. 2009, Diabetes Care, PMID:19196896: 14 T2D (A1C 6.8%, oral "
    "agents stopped 3 weeks earlier, 'relatively early' disease) vs. 11 "
    "age/weight/BMI-matched non-diabetic subjects, mixed meal (1.2 g/kg "
    "glucose within a 45/40/15 carbohydrate/fat/protein meal). Every figure "
    "was re-fetched from the full text with biosim_literature.verification. "
    "Scored as T2D/healthy RATIOS so the model's absolute insulin scale "
    "cancels -- but absolute insulin is also far off (model T2D peak ~106 pM "
    "vs. 382 pmol/l observed), so the ratios are not the whole gap. The "
    "glucose iAUC is printed by the paper as 'mmol/l over 6 h' (almost "
    "certainly mmol*min/l); only the unit-free ratio is used. DISCLOSED "
    "MISMATCHES: (1) the model arm is a flat 75 g glucose drink with no "
    "fat/protein (avoids a NutrientSecretionPathway confound), whereas the "
    "source dose is body-weight-scaled (probably ~100 g) inside a mixed "
    "meal; the peak-DELAY point is the one that is not scale-free and is "
    "meal-sensitive (the source's 62 min healthy peak reflects fat/protein "
    "slowing gastric emptying; the model's glucose drink peaks at ~39 min) -- "
    "a 100 g run changes the ratios only slightly (peak 0.23, 6 h 0.12). "
    "(2) The source cohort is more hyperglycaemic than the preset (fasting "
    "glucose 9.1 vs. 6.5 mM). In the source, fasting insulin did not differ "
    "(59 vs. 37 pmol/l), the early response was blunted and delayed, and the "
    "total 6 h insulin response did not differ; the T2D preset gives only "
    "~0.1-0.25x of the healthy insulin response on every metric. Whether that "
    "is realistic for the preset's intended population is MASTER_PLAN row "
    "1.46's OPEN question -- this benchmark makes it measurable. WIDE "
    "tolerance: NOT a target to tune secretion against."
)


def t2d_vs_nondiabetic_meal_insulin() -> BenchmarkDataset:
    """
    T2D vs. matched non-diabetic subjects after the same oral glucose load --
    Basu et al. 2009, PMID:19196896. Exercises the default
    DiseaseState.TYPE_2_DIABETES preset's postprandial insulin response,
    which had no registered benchmark and whose p3-only mechanisms are
    near-inert as a result (MASTER_PLAN row 1.46). Expected to MISS; see
    T2D_MEAL_INSULIN_NOTES.
    """
    ds = BenchmarkDataset(
        name="T2D_vs_nondiabetic_meal_insulin",
        description=(
            "T2D vs. matched non-diabetic subjects, same oral glucose load, "
            "insulin and glucose responses as T2D/healthy ratios -- Basu et al. 2009"
        ),
    )
    points = (
        ("insulin_peak_ratio_t2d_vs_nondiabetic",
         "Peak plasma insulin, T2D / non-diabetic (382 vs. 673 pmol/l).",
         382.0 / 673.0, "ratio"),
        ("insulin_iauc_2h_ratio_t2d_vs_nondiabetic",
         "Insulin increment above basal over 0-2 h, T2D / non-diabetic (21 vs. 52 nmol/l).",
         21.0 / 52.0, "ratio"),
        ("insulin_iauc_6h_ratio_t2d_vs_nondiabetic",
         "Insulin increment above basal over 0-6 h, T2D / non-diabetic (79 vs. 80 nmol/l).",
         79.0 / 80.0, "ratio"),
        ("insulin_peak_delay_min_t2d_minus_nondiabetic",
         "Time of peak insulin, T2D minus non-diabetic (156 vs. 62 min).",
         156.0 - 62.0, "min"),
        ("glucose_iauc_6h_ratio_t2d_vs_nondiabetic",
         "Glucose increment above basal over 0-6 h, T2D / non-diabetic (1,671 vs. 557).",
         1671.0 / 557.0, "ratio"),
    )
    for name, description, reference, unit in points:
        ds.add(BenchmarkPoint(
            name=name, description=description, reference=reference, unit=unit,
            tolerance=Tolerance.WIDE, pubmed_id="19196896",
            population="type_2_diabetes", notes=T2D_MEAL_INSULIN_NOTES,
        ))
    return ds


CAFFEINE_OGTT_NOTES = (
    "Graham et al. 2001, Can J Physiol Pharmacol, PMID:11478588: 18 young fit "
    "men, caffeine 5 mg/kg or placebo (double blind), 75 g dextrose 1 h later; "
    "serum insulin AUC 60% greater (C-peptide 37%), glucose AUC +24% but not "
    "significant (P = 0.20), and NO difference in pre-OGTT glucose between "
    "trials. Robinson et al. 2004, J Nutr, PMID:15465742: 12 obese men with "
    "T2D (BMI 32; fasting glucose ~6.7 mM, close to the preset's 6.5), same "
    "protocol over 3 h, oral agents withdrawn 48 h; insulin AUC 25% greater "
    "(P < 0.05), glucose also higher. Every quoted figure was re-fetched from "
    "the live abstract. Independent of the caffeine calibration source "
    "(Keijzers 2002, PMID:11815511). DISCLOSED: neither abstract says whether "
    "'AUC' is total or incremental (the model's value uses TOTAL AUC over 0-3 "
    "h of the OGTT; incremental gives the same signs) and Graham's OGTT length "
    "is not stated, so the point names carry no window; the model uses the "
    "default 70 kg profile (350 mg) rather than either cohort's weight, age or "
    "fitness, a dose beyond the uncited dose-response's 200 mg reference; "
    "caffeine's effect is a static snapshot at t = 1 h after dosing held "
    "constant for the whole run (not a time course); glucose points are not "
    "scored (Graham's +24% is not significant, Robinson reports no number). "
    "At the time of writing the model gives the WRONG SIGN for healthy "
    "(about -5%) and essentially no effect for T2D (about -1%): the p3 term "
    "alone gives +0.9% / +0.0%, and an uncited Gb factor (peak +2%, ~+1% at "
    "1 h; pubmed_id UNVERIFIED in parameter_seed.py) lowers insulin further "
    "while also raising fasting glucose, which Graham's own data contradict. "
    "Removing the Gb term would NOT close the gap (+60% / +25% observed): the "
    "insulin-raising mechanism of caffeine is not modelled. WIDE tolerance, "
    "NOT a target to tune the Gb term against."
)


def caffeine_ogtt_insulin() -> BenchmarkDataset:
    """
    Caffeine 5 mg/kg one hour before a 75 g OGTT -- insulin AUC increase in
    healthy young men (Graham 2001, PMID:11478588) and men with type 2
    diabetes (Robinson 2004, PMID:15465742). Caffeine had no registered
    benchmark. Expected to MISS (wrong sign); see CAFFEINE_OGTT_NOTES.
    """
    ds = BenchmarkDataset(
        name="Caffeine_5mgkg_OGTT_insulin",
        description=(
            "Caffeine 5 mg/kg vs. placebo 1 h before a 75 g OGTT, serum insulin "
            "AUC change -- Graham 2001 (healthy) and Robinson 2004 (T2D)"
        ),
    )
    ds.add(BenchmarkPoint(
        name="insulin_auc_percent_increase_healthy_caffeine",
        description="Percent increase in OGTT insulin AUC with caffeine vs. placebo, healthy men.",
        reference=60.0, unit="%", tolerance=Tolerance.WIDE,
        pubmed_id="11478588", population="healthy_adults", notes=CAFFEINE_OGTT_NOTES,
    ))
    ds.add(BenchmarkPoint(
        name="insulin_auc_percent_increase_t2d_caffeine",
        description="Percent increase in OGTT insulin AUC with caffeine vs. placebo, men with T2D.",
        reference=25.0, unit="%", tolerance=Tolerance.WIDE,
        pubmed_id="15465742", population="type_2_diabetes", notes=CAFFEINE_OGTT_NOTES,
    ))
    return ds
