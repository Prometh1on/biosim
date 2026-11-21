# SAVE TO: packages\data-seeds\src\biosim_seeds\population.py
"""
Population simulation engine — cohort-level analysis.

Demographic defaults in PopulationSpec (age/BMI/T2D/prediabetes) are
literature-anchored to represent a general, unselected US adult
population (18-80 years), matching this platform's existing US/NHANES-
anchored precedent elsewhere (Bergman 1981, NHANES-derived meal GI
values, DeFronzo age modifiers, etc.) -- see PopulationSpec's own
docstring for each value's source. Every one of them was previously an
uncited placeholder (age_mean=40/sd=10, bmi_mean=26/sd=4,
t2d_fraction=0.00, insulin_resistant_fraction=0.10) with no literature
anchor at all; verified live this session (via PubMedFetcher, US Census
Bureau bulk data, and CDC/NCHS reports) that every one of them was
substantially off from the real population statistics -- most strikingly
t2d_fraction=0.00, meaning the previous default cohort had literally
zero diabetic individuals, not representative of any real adult
population (real prevalence is ~12-14%).
"""

from __future__ import annotations
import dataclasses
import math
import random
from concurrent.futures import Executor, ProcessPoolExecutor
from dataclasses import dataclass, field
from typing import Optional
import numpy as np

from biosim_api import SimulationRunner, ScenarioDefinition
from biosim_api.scenario_warnings import experimental_drug_warnings
from biosim_physiology import (
    BiologicalSex, DiseaseState, FitnessLevel, PhysiologicalProfile,
)


def _run_scenario(scenario):
    """
    Module-level, picklable worker entry point for PopulationRunner's
    optional parallel execution -- same shape and same reasoning as
    biosim_expert.batch._run_scenario (not imported directly: reaching
    into a sibling package's private, underscore-prefixed function for
    a 2-line body would be a worse coupling than duplicating it, the
    same encapsulation call this session already made once this
    session for LogNormalDistribution._mu/._sigma). Must be a plain
    module-level function, not a bound method or closure, for
    ProcessPoolExecutor to pickle it across the process boundary.
    """
    return SimulationRunner().run(scenario)


@dataclass
class PopulationSpec:
    """
    Demographic/disease-state distribution parameters for synthetic
    cohort generation.

    Literature sources (verified live this session against primary
    data, not a secondary paraphrase):

    age_mean/age_sd: US Census Bureau, Population Division, "Annual
        Estimates of the Resident Population by Single Year of Age and
        Sex for the United States: April 1, 2020 to July 1, 2023"
        (vintage 2023 bulk data, nc-est2023-agesex-res.csv). Mean/SD
        computed directly from the raw single-year-of-age population
        counts, truncated to 18-80 (matching age_min/age_max exactly):
        mean=46.50, SD=17.26 (n=250,727,757). The real SD is ~70% larger
        than the previous uncited default (10.0) -- the 18-80 adult
        population is much flatter across decades than a tight cluster
        around 40 would suggest.

    bmi_mean/bmi_sd: Fryar CD, Gu Q, Afful J, Carroll MD, Ogden CL.
        "Anthropometric Reference Data for Children and Adults: United
        States, August 2021-August 2023." National Center for Health
        Statistics, Vital and Health Statistics, Series 3, No. 50,
        June 2025. PMID:41670955. Table 10 (NHANES-weighted, ages 20+):
        men mean 29.4 kg/m^2 (P5=20.9, P95=42.1), women mean 30.0 kg/m^2
        (P5=19.7, P95=44.9). Combined (population-weighted ~49.5/50.5
        M/F per the same Census data above): mean ~= 29.7 kg/m^2. NHANES
        publishes percentiles and SE-of-mean, not population SD directly
        -- SD estimated from the reported percentile spread (IQR-based,
        10-90-based, and 5-95-based estimates agree to within ~10%):
        combined SD ~= 7.0 kg/m^2 (up from the previous uncited 4.0).
        The three percentile-pair SD estimates disagreeing with each
        other (more for men than women) is itself direct empirical
        evidence of real right skew -- see BMI_LOGNORMAL_CV note below.
        Caveat stated honestly in the source's own methods: this NHANES
        cycle's exam response rate was only 25.6% (COVID-era protocol
        change), though it remains the current official NCHS reference.

    t2d_fraction/insulin_resistant_fraction: Menke A, Casagrande S,
        Geiss L, Cowie CC. "Prevalence of and Trends in Diabetes Among
        Adults in the United States, 1988-2012." JAMA. 2015;314(10):
        1021-1029. PMID:26348752. NHANES 2011-2012 (n=2,781), using the
        combined hemoglobin-A1c-or-FPG-or-2-hour-PG case definition (the
        broadest of the paper's three definitions, chosen so both
        fractions below use the SAME diagnostic criteria basis --
        mixing a narrow definition for one and a broad one for the
        other would make this single synthetic cohort internally
        incoherent): total diabetes 14.3% (95% CI 12.2-16.8%),
        prediabetes 38.0% (95% CI 34.7-41.3%). insulin_resistant_fraction
        is the correct target for the prediabetes figure specifically --
        confirmed DiseaseState.INSULIN_RESISTANT's own docstring in
        physiology/profile.py already documents it as "pre-diabetes /
        metabolic syndrome". These prevalence rates were confirmed
        still current via CDC's live National Diabetes Statistics
        Report (2023 data): ~15.2% adult diabetes, ~44.0% adult
        prediabetes -- both directionally consistent with, and higher
        than, Menke 2015's decade-old figures (diabetes prevalence has
        been rising), so Menke's own reported numbers are used as the
        directly-quotable peer-reviewed source rather than this
        session's own arithmetic derived from a live government web
        page.

    male_fraction: 0.50 (previous default) confirmed still correct --
        same Census data gives 0.4946, within 0.5 percentage points.
        No change.

    fitness_weights: NOT verified this session -- no standard
        epidemiological survey reports a population distribution across
        this platform's specific five-tier fitness classification
        (sedentary/lightly_active/moderately_active/very_active/athlete),
        unlike age/BMI/diabetes status, which are all standard, directly
        measured NHANES/Census quantities. Left as the pre-existing,
        still-uncited judgment call -- flagged here rather than silently
        presented as equally well-anchored to the values above.

    sample_genetic_variants: opt-in (default False, matching this
        project's convention for other experimental features --
        model_nefa_dynamics, exercise_bout_in_progress -- rather than
        silently changing every existing caller's cohort composition).
        When True, each individual is independently assigned a genotype
        at each of the 7 default T2D-risk variants in
        biosim_genotype.effect_db.VariantEffectDB.default() that has a
        verified population_frequency (all 7 do, as of this session --
        see effect_db.py's own module docstring for the TOPMed source),
        via Hardy-Weinberg equilibrium (assuming linkage equilibrium
        across these loci -- a standard simplification for unlinked
        GWAS hits on different genes/chromosomes) and threaded into that
        individual's simulation via PersonalizationSpec.genetic_variants
        (already-existing machinery -- GenotypePersonalizer and its
        SimulationRunner wiring were NOT built by this change, only
        newly driven by real, sampled data instead of never being
        exercised by population-level cohorts at all).
    """
    n_individuals:      int   = 100
    age_mean:           float = 46.50
    age_sd:             float = 17.26
    age_min:            float = 18.0
    age_max:            float = 80.0
    bmi_mean:           float = 29.7
    bmi_sd:             float = 7.0
    bmi_min:            float = 18.5
    bmi_max:            float = 45.0
    male_fraction:      float = 0.50
    t2d_fraction:       float = 0.143
    insulin_resistant_fraction: float = 0.380
    fitness_weights:    dict  = field(default_factory=lambda: {
        "sedentary": 0.25, "lightly_active": 0.30, "moderately_active": 0.30,
        "very_active": 0.12, "athlete": 0.03,
    })
    sample_genetic_variants: bool = False
    time_of_day_hours: Optional[float] = None
    random_seed:        int   = 42

    def __post_init__(self) -> None:
        # Found during the mandatory security-reviewer pass on the
        # demographic-anchoring fix above, see CLAUDE.md: this dataclass
        # had no validation at all, unlike every other dataclass this
        # project has hardened with fail-fast constructor validation
        # this session (DrugEvent, PhysiologicalProfile,
        # ScenarioDefinition). Confirmed reachable via the public
        # constructor: PopulationSpec(bmi_mean=0.0) raised an unhandled
        # ZeroDivisionError inside _sample_bmi()'s cv=bmi_sd/bmi_mean;
        # PopulationSpec(bmi_mean=-5.0) raised ValueError: math domain
        # error from math.log() on the negative mean -- two different,
        # confusing exception types depending on the sign of one field,
        # neither caught anywhere. Separately, t2d_fraction/
        # insulin_resistant_fraction were completely unbounded -- a sum
        # exceeding 1.0 doesn't crash, but silently makes HEALTHY
        # unreachable and produces a cohort with unintended proportions,
        # not what the caller asked for. This is a materially bigger
        # risk now than before this fix: previously 0.00+0.10=0.10 left
        # 0.90 of headroom before a caller-supplied override could
        # exceed 1.0; the new, real-prevalence-anchored defaults
        # (0.143+0.380=0.523) leave only 0.477.
        # math.isfinite first: NaN fails every comparison below (so would
        # slip through), and inf bmi_mean silently pins every individual
        # at bmi_max instead of raising.
        if not math.isfinite(self.bmi_mean) or self.bmi_mean <= 0:
            raise ValueError(f"bmi_mean must be finite and positive, got {self.bmi_mean}")
        if not math.isfinite(self.bmi_sd) or self.bmi_sd < 0:
            raise ValueError(f"bmi_sd must be finite and non-negative, got {self.bmi_sd}")
        # bmi_min/bmi_max clip the sampled BMI, which then becomes a
        # PhysiologicalProfile -- whose own constructor rejects a BMI
        # outside [MIN_PLAUSIBLE_BMI, MAX_PLAUSIBLE_BMI]. Validate the clip
        # bounds here against the same limits so a bad bound fails at spec
        # construction, not mid-cohort for whichever individual first
        # samples past it. The `not (lo <= x <= hi)` form also rejects NaN.
        lo_limit = PhysiologicalProfile.MIN_PLAUSIBLE_BMI
        hi_limit = PhysiologicalProfile.MAX_PLAUSIBLE_BMI
        for name, value in (("bmi_min", self.bmi_min), ("bmi_max", self.bmi_max)):
            if not (lo_limit <= value <= hi_limit):
                raise ValueError(
                    f"{name} must be between {lo_limit} and {hi_limit} "
                    f"(PhysiologicalProfile's plausible-BMI limits), got {value}"
                )
        if self.bmi_min > self.bmi_max:
            raise ValueError(
                f"bmi_min ({self.bmi_min}) must not exceed bmi_max ({self.bmi_max})"
            )
        if self.age_mean <= 0:
            raise ValueError(f"age_mean must be positive, got {self.age_mean}")
        if self.age_sd < 0:
            raise ValueError(f"age_sd must be non-negative, got {self.age_sd}")
        if not 0.0 <= self.male_fraction <= 1.0:
            raise ValueError(f"male_fraction must be in [0, 1], got {self.male_fraction}")
        if not 0.0 <= self.t2d_fraction <= 1.0:
            raise ValueError(f"t2d_fraction must be in [0, 1], got {self.t2d_fraction}")
        if not 0.0 <= self.insulin_resistant_fraction <= 1.0:
            raise ValueError(
                f"insulin_resistant_fraction must be in [0, 1], got "
                f"{self.insulin_resistant_fraction}"
            )
        if self.t2d_fraction + self.insulin_resistant_fraction > 1.0:
            raise ValueError(
                f"t2d_fraction ({self.t2d_fraction}) + "
                f"insulin_resistant_fraction ({self.insulin_resistant_fraction}) "
                f"= {self.t2d_fraction + self.insulin_resistant_fraction} exceeds "
                "1.0 -- these are mutually exclusive disease-state categories "
                "assigned from the same [0,1) draw, so their sum must not "
                "exceed 1.0 or HEALTHY becomes silently unreachable."
            )


@dataclass
class CohortSummary:
    n_total:               int
    n_succeeded:           int
    glucose_peak_mean_mM:  float
    glucose_peak_sd_mM:    float
    glucose_peak_p5_mM:    float
    glucose_peak_p95_mM:   float
    insulin_peak_mean_pM:  float
    mtor_peak_mean:        float
    ampk_nadir_mean:       float
    pct_above_7p8_mM:      float
    pct_above_11p1_mM:     float
    individual_profiles:   list[PhysiologicalProfile]
    individual_peaks_mM:   list[float]
    failure_reasons:       list[str] = field(default_factory=list)
    # Scenario-level disclosures that apply to every individual, e.g. the
    # "EXPERIMENTAL drug module" note when the scenario applies a drug.
    warnings:              list[str] = field(default_factory=list)

    @property
    def n_failed(self) -> int:
        return self.n_total - self.n_succeeded

    def summary_dict(self) -> dict:
        return {
            "n_total":            self.n_total,
            "n_succeeded":        self.n_succeeded,
            "n_failed":           self.n_failed,
            "failure_reasons":    list(self.failure_reasons),
            "warnings":           list(self.warnings),
            "glucose_peak_mM":    f"{self.glucose_peak_mean_mM:.2f} +/- {self.glucose_peak_sd_mM:.2f}",
            "glucose_peak_p5_95": f"{self.glucose_peak_p5_mM:.1f} - {self.glucose_peak_p95_mM:.1f}",
            "pct_above_7p8_mM":   f"{self.pct_above_7p8_mM:.1f}%",
            "pct_above_11p1_mM":  f"{self.pct_above_11p1_mM:.1f}%",
            "insulin_peak_pM":    f"{self.insulin_peak_mean_pM:.0f}",
            "mtor_peak_mean":     f"{self.mtor_peak_mean:.3f}",
        }


class PopulationRunner:
    def __init__(
        self,
        spec: Optional[PopulationSpec] = None,
        use_multiprocessing: bool = False,
        max_workers: int = 8,
        executor_cls: type[Executor] = ProcessPoolExecutor,
    ) -> None:
        """
        use_multiprocessing: opt-in (default False), matching this
            project's established convention for other experimental/
            performance features (model_nefa_dynamics,
            exercise_bout_in_progress, sample_genetic_variants) rather
            than silently changing every existing caller's behaviour --
            unlike BatchRunner (Track 3's prior milestone), which
            switched its DEFAULT to ProcessPoolExecutor, benchmarking
            directly here found the real crossover point sits above
            PopulationSpec's own default n_individuals=100: at n=100,
            8-worker ProcessPoolExecutor measured 2.326s vs. 2.222s
            sequential (SLOWER, small-n process-spawn overhead not yet
            amortized -- Windows uses "spawn," every worker re-imports
            this platform's full dependency tree); at n=400 it measured
            a genuine 1.73x speedup. Flipping the default would make the
            platform's own most common call pattern slightly worse, not
            better -- opt in explicitly for large cohorts (n_individuals
            in the hundreds or more) where the real benefit applies.
        executor_cls: same testability escape hatch as BatchRunner --
            a spawned worker process never sees a parent-process
            monkeypatch (e.g. `monkeypatch.setattr(SimulationRunner,
            "run", ...)`), so failure-injection tests pass
            ThreadPoolExecutor explicitly when both use_multiprocessing
            and monkeypatch-based failure simulation are needed
            together.

        Privacy note (security-reviewer pass, no code change needed):
            with sample_genetic_variants=True, raw Variant/rsid data
            transiently exists in up to max_workers worker processes'
            memory simultaneously rather than only the parent process's
            -- a modestly wider OS-level memory-exposure surface (e.g.
            a crash-dump scenario), not a new class of leak (the same
            raw data is already necessarily in memory to run any
            genetically-personalized simulation at all, sequential or
            not). Confirmed directly: no new disk-persistence path
            (Windows ProcessPoolExecutor ships pickled data over an
            in-memory duplex pipe, not temp files), and every exception
            path that could touch genetic data (GenotypePersonalizer
            failures in runner.py, composer.py's population-mismatch
            warnings) is already sanitized at its own source, before
            crossing back into this module's _run_all(), regardless of
            which process produced it.
        """
        self.spec = spec or PopulationSpec()
        self.use_multiprocessing = use_multiprocessing
        self.max_workers = max_workers
        self.executor_cls = executor_cls

    def _sample_bmi(self, np_rng: np.random.Generator) -> float:
        """
        Log-normal, not Normal -- real adult BMI is right-skewed (a long
        upper tail into obesity), confirmed directly this session: three
        independent SD estimates from NHANES's own reported percentiles
        (IQR-based, 10-90-based, 5-95-based) disagree with each other
        specifically because the upper tail is fatter than a symmetric
        distribution would produce (see PopulationSpec's own docstring).
        A clipped Normal would cut off that tail unrealistically; a
        log-normal (already this project's standard shape for always-
        positive, multiplicatively-uncertain quantities -- see
        biosim_uncertainty.distributions.LogNormalDistribution, same
        mean/CV parameterisation reused here directly rather than
        importing that package, to avoid a new cross-package dependency
        for one formula) represents it more faithfully while still
        matching the real reported mean and SD exactly.
        """
        cv = self.spec.bmi_sd / self.spec.bmi_mean
        sigma = math.sqrt(math.log(1.0 + cv ** 2))
        mu = math.log(self.spec.bmi_mean) - 0.5 * sigma ** 2
        return float(np_rng.lognormal(mean=mu, sigma=sigma))

    @staticmethod
    def _sample_genotype_dosage(np_rng: np.random.Generator, frequency: float) -> int:
        """
        Hardy-Weinberg equilibrium: given a risk-allele frequency q,
        P(0 copies)=(1-q)^2, P(1 copy)=2q(1-q), P(2 copies)=q^2. Standard
        population-genetics assumption for a biallelic locus not under
        recent selection/non-random mating -- the same assumption
        implicit in every GWAS allele-frequency estimate this platform
        already cites (Grant 2006, Altshuler 2000, etc. all report
        allele, not genotype, frequencies and rely on HWE to interpret
        them clinically).
        """
        return int(np_rng.choice([0, 1, 2], p=[(1 - frequency) ** 2, 2 * frequency * (1 - frequency), frequency ** 2]))

    def _sample_genetic_variants_for_individual(self, np_rng: np.random.Generator, db=None) -> list:
        """
        One independent Hardy-Weinberg draw per default T2D-risk variant
        that has a verified population_frequency (all 7 defaults, as of
        this session -- see effect_db.py). Assumes linkage equilibrium
        across loci (standard for unlinked GWAS hits on different genes/
        chromosomes, same assumption VariantEffectDB/EffectComposer
        already make when combining per-variant effects independently).

        Constructs each Variant with NO ref_allele/alt_allele/
        observed_alleles set, deliberately -- Variant.dosage_for_allele()
        falls back to "assume alt=risk" in that case (see its own
        docstring), which is exactly correct here since the sampled
        dosage IS already a risk-allele count by construction (Hardy-
        Weinberg draws directly in those terms), not an ambiguous raw
        genotype call needing separate ref/alt orientation.

        `db` lets callers (run(), which loops per-individual) build
        VariantEffectDB.default() once and reuse it -- matching the
        established convention elsewhere in this codebase
        (GenotypePersonalizer.__init__ caches it as self._effect_db
        rather than rebuilding per call). Rebuilding it here by default
        keeps this method usable standalone.
        """
        from biosim_genotype import Genotype, Variant, VariantEffectDB

        dosage_to_genotype = {
            0: Genotype.HOMOZYGOUS_REF,
            1: Genotype.HETEROZYGOUS,
            2: Genotype.HOMOZYGOUS_ALT,
        }
        db = db if db is not None else VariantEffectDB.default()
        variants = []
        for rsid in db.all_rsids():
            record = db.get(rsid)[0]
            if record.population_frequency is None:
                continue
            dosage = self._sample_genotype_dosage(np_rng, record.population_frequency)
            variants.append(Variant(rsid=rsid, genotype=dosage_to_genotype[dosage]))
        return variants

    def generate_cohort(self) -> list[PhysiologicalProfile]:
        rng    = random.Random(self.spec.random_seed)
        np_rng = np.random.default_rng(self.spec.random_seed)
        fitness_levels = list(self.spec.fitness_weights.keys())
        fitness_probs  = list(self.spec.fitness_weights.values())
        total_p = sum(fitness_probs)
        fitness_probs = [p / total_p for p in fitness_probs]
        profiles = []

        for _ in range(self.spec.n_individuals):
            age    = float(np.clip(np_rng.normal(self.spec.age_mean, self.spec.age_sd),
                                   self.spec.age_min, self.spec.age_max))
            bmi    = float(np.clip(self._sample_bmi(np_rng),
                                   self.spec.bmi_min, self.spec.bmi_max))
            weight = bmi * (1.70 ** 2)
            sex    = BiologicalSex.MALE if rng.random() < self.spec.male_fraction else BiologicalSex.FEMALE
            fkey   = rng.choices(fitness_levels, weights=fitness_probs)[0]
            disease = [DiseaseState.HEALTHY]
            r = rng.random()
            if r < self.spec.t2d_fraction:
                disease = [DiseaseState.TYPE_2_DIABETES]
            elif r < self.spec.t2d_fraction + self.spec.insulin_resistant_fraction:
                disease = [DiseaseState.INSULIN_RESISTANT]

            profiles.append(PhysiologicalProfile(
                age_years=age, biological_sex=sex,
                body_weight_kg=weight, height_cm=170.0,
                fitness_level=FitnessLevel(fkey),
                disease_states=disease,
                time_of_day_hours=self.spec.time_of_day_hours,
            ))
        return profiles

    def run(self, scenario: ScenarioDefinition) -> CohortSummary:
        """
        Run the full cohort through the given scenario.

        Cohort generation and genetic-variant sampling (both above)
        stay strictly sequential regardless of `use_multiprocessing` --
        both depend on drawing from ordered RNG streams keyed to
        `self.spec.random_seed`, and distributing those draws across
        worker processes would make the resulting cohort composition
        depend on process-scheduling timing, silently breaking
        reproducibility for the same seed. Only the actual per-
        individual ODE simulation -- the expensive, independent,
        embarrassingly-parallel step -- is optionally parallelized.
        """
        from biosim_api import PersonalizationSpec
        profiles = self.generate_cohort()
        peaks_mM, insulin_pM, mtor_peaks, ampk_nadirs = [], [], [], []
        succeeded_profiles = []
        failure_reasons = []

        # Offset-seeded, matching the pattern BayesianCalibrator's own
        # posterior-resampling RNG already uses: a dedicated RNG for this
        # independent random process (genotype sampling), not reusing
        # generate_cohort()'s RNG (which would silently desync the two
        # streams the moment one draws a different number of randoms per
        # individual than the other, e.g. skipping a disease-state draw).
        genetics_rng = (
            np.random.default_rng(self.spec.random_seed + 1)
            if self.spec.sample_genetic_variants else None
        )
        effect_db = None
        if genetics_rng is not None:
            from biosim_genotype import VariantEffectDB
            effect_db = VariantEffectDB.default()

        personalised_scenarios = []
        for profile in profiles:
            genetic_variants = (
                self._sample_genetic_variants_for_individual(genetics_rng, db=effect_db)
                if genetics_rng is not None else []
            )
            personalised_scenarios.append(dataclasses.replace(
                scenario,
                personalization=PersonalizationSpec(
                    physiological_profile=profile,
                    genetic_variants=genetic_variants,
                ),
            ))

        for profile, (output, error) in zip(profiles, self._run_all(personalised_scenarios)):
            try:
                if error is not None:
                    raise error
                if output.success:
                    m = output.expert_metrics
                    peaks_mM.append(m.glucose_peak_mM)
                    insulin_pM.append(m.insulin_peak_pM)
                    mtor_peaks.append(m.mtor_peak)
                    ampk_nadirs.append(m.ampk_nadir)
                    succeeded_profiles.append(profile)
                else:
                    failure_reasons.append(
                        f"age={profile.age_years:.0f} sex={profile.biological_sex.value}: "
                        + "; ".join(output.warnings or ["simulation reported failure, no warning detail"])
                    )
            except Exception as e:
                # Individual profiles can legitimately fail (e.g. an
                # extreme age/BMI combination the ODE solver can't
                # integrate) without invalidating the whole cohort run --
                # but per ECC error-handling rules, that must be recorded,
                # not silently discarded (see bandit B110 finding this
                # fixed: try/except/pass with no error context at all).
                failure_reasons.append(
                    f"age={profile.age_years:.0f} sex={profile.biological_sex.value}: "
                    f"{type(e).__name__}: {e}"
                )

        n_ok = len(peaks_mM)
        if n_ok == 0:
            return CohortSummary(n_total=len(profiles), n_succeeded=0,
                glucose_peak_mean_mM=0, glucose_peak_sd_mM=0,
                glucose_peak_p5_mM=0, glucose_peak_p95_mM=0,
                insulin_peak_mean_pM=0, mtor_peak_mean=0, ampk_nadir_mean=0,
                pct_above_7p8_mM=0, pct_above_11p1_mM=0,
                individual_profiles=[], individual_peaks_mM=[],
                failure_reasons=failure_reasons,
                warnings=experimental_drug_warnings(scenario))

        peaks = np.array(peaks_mM)
        return CohortSummary(
            n_total=len(profiles), n_succeeded=n_ok,
            glucose_peak_mean_mM=float(np.mean(peaks)),
            glucose_peak_sd_mM=float(np.std(peaks)),
            glucose_peak_p5_mM=float(np.percentile(peaks, 5)),
            glucose_peak_p95_mM=float(np.percentile(peaks, 95)),
            insulin_peak_mean_pM=float(np.mean(insulin_pM)),
            mtor_peak_mean=float(np.mean(mtor_peaks)),
            ampk_nadir_mean=float(np.mean(ampk_nadirs)),
            pct_above_7p8_mM=float(np.mean(peaks > 7.8) * 100),
            pct_above_11p1_mM=float(np.mean(peaks > 11.1) * 100),
            individual_profiles=succeeded_profiles,
            individual_peaks_mM=peaks_mM,
            failure_reasons=failure_reasons,
            warnings=experimental_drug_warnings(scenario),
        )

    def _run_all(self, scenarios: list) -> list:
        """
        Run each scenario, returning (output, None) on success or
        (None, exception) on failure, in the SAME order as `scenarios`
        -- the caller zips this directly against `profiles`, so order
        must be preserved exactly (unlike BatchRunner's as_completed()-
        based approach, which doesn't need order since it re-associates
        results by scenario name afterward).

        Sequential by default (self.use_multiprocessing=False) --
        see __init__'s docstring for the benchmark data behind that
        default.
        """
        if not self.use_multiprocessing:
            runner = SimulationRunner()
            results = []
            for s in scenarios:
                try:
                    results.append((runner.run(s), None))
                except Exception as e:
                    results.append((None, e))
            return results

        results = [None] * len(scenarios)
        with self.executor_cls(max_workers=self.max_workers) as executor:
            futures = [executor.submit(_run_scenario, s) for s in scenarios]
            for i, future in enumerate(futures):
                try:
                    results[i] = (future.result(), None)
                except Exception as e:
                    results[i] = (None, e)
        return results
