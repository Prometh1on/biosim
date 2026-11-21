# SAVE TO: packages\validation\src\biosim_validation\validator.py
"""
Model validator -- runs the simulation and checks outputs against benchmarks.

Validation methodology:
    1. Run the specified scenario (default: 75g OGTT, healthy population average)
    2. Extract relevant quantities from SimulationOutput
    3. Check each against its BenchmarkPoint tolerance band
    4. Produce a ValidationReport with PASS/MISS per point and an overall verdict

Verdict criteria (checked in this order):
    NO DATA: No cited (scored) point at all -- nothing was validated
    FAIL:  Any TIGHT benchmark fails, or majority of NORMAL fail
    WARN:  All TIGHT benchmarks pass; a minority of NORMAL fail
    MISS:  Every TIGHT/NORMAL point passes, but at least one cited
           (scored) point lies outside its band -- in practice a WIDE
           (+/-50%) point the model misses by 50% or more
    PASS:  Every cited point lies inside its tolerance band

Points with no verified citation (pubmed_id="UNVERIFIED") are
informational only and never change the verdict.

The validator is intentionally conservative. A model that barely passes
every benchmark is not well-validated — it just hasn't been falsified yet.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Callable, Optional
import numpy as np

from .benchmark import BenchmarkCheckResult, BenchmarkDataset, Tolerance


class ValidationVerdict(str, Enum):
    PASS = "PASS"
    MISS = "MISS"
    NO_DATA = "NO DATA"
    WARN = "WARN"
    FAIL = "FAIL"


@dataclass
class ValidationReport:
    """Full validation report for one benchmark dataset."""
    dataset_name:   str
    scenario_name:  str
    results:        list[BenchmarkCheckResult]
    verdict:        ValidationVerdict
    n_pass:         int
    n_fail:         int
    notes:          list[str] = field(default_factory=list)
    # "" for non-drug benchmarks; for a drug benchmark, says the drug module
    # is experimental (see _drug_module_label). Never changes the verdict.
    module_label:   str = ""

    @property
    def pass_rate(self) -> float:
        total = self.n_pass + self.n_fail
        return self.n_pass / total if total > 0 else 0.0

    @property
    def failed_benchmarks(self) -> list[BenchmarkCheckResult]:
        return [r for r in self.results if not r.passed]

    def summary(self) -> dict:
        return {
            "dataset":     self.dataset_name,
            "scenario":    self.scenario_name,
            "verdict":     self.verdict.value,
            "n_pass":      self.n_pass,
            "n_fail":      self.n_fail,
            "pass_rate":   round(self.pass_rate, 3),
            "failed":      [r.name for r in self.failed_benchmarks],
            "module_label": self.module_label,
            "points": [
                {
                    "name":              r.name,
                    "tolerance":         r.benchmark.tolerance.value,
                    "deviation_percent": round(r.deviation_percent, 1),
                    "status":            r.status,
                }
                for r in self.results
            ],
        }

    def as_text(self) -> str:
        lines = [
            f"VALIDATION REPORT: {self.dataset_name}",
            f"Scenario: {self.scenario_name}",
            f"Verdict: {self.verdict.value}  "
            f"({self.n_pass}/{self.n_pass+self.n_fail} points within tolerance)",
            "",
        ]
        for r in self.results:
            tol = r.benchmark.tolerance.value.upper()
            dev = f"{r.deviation_percent:+.1f}%"
            ref = f"{r.benchmark.reference} {r.benchmark.unit}"
            mod = f"{r.model_value:.3g} {r.benchmark.unit}"
            # A point with no real citation still has a true .status (used
            # elsewhere, e.g. failed_benchmarks) but displaying it as a
            # bare PASS/MISS here would read as a genuine validation check
            # -- see _unscored_point_notes(), which explains this same
            # point in the Notes section below.
            tag = r.status if r.benchmark.is_scored else f"{r.status}/INFO"
            lines.append(
                f"  [{tag}] {r.name} [{tol}]  "
                f"ref={ref}  model={mod}  deviation={dev}"
            )
        if self.module_label:
            lines += ["", self.module_label]
        if self.notes:
            lines += ["", "Notes:"] + [f"  - {n}" for n in self.notes]
        return "\n".join(lines)


def _drug_module_label(drug_name: str) -> str:
    """
    Label for a drug benchmark, from the single drug status table
    (biosim_physiology.drug_interactions.validation_status), so a reader sees
    a PASS or MISS in context. The verdict is computed separately and is
    never changed by this label.
    """
    from biosim_physiology.drug_interactions import Drug
    from biosim_physiology.drug_interactions.validation_status import (
        DRUG_VALIDATION,
        DrugValidationStatus,
    )

    record = DRUG_VALIDATION[Drug(drug_name)]
    if record.status is DrugValidationStatus.VALIDATED:
        return (
            f"Drug benchmark for a validated module ('{drug_name}', status "
            f"VALIDATED, reviewed). Caveat: {record.caveat}"
        )
    status = record.status.value
    return (
        f"Drug benchmark for an {status.lower()} module ('{drug_name}', status "
        f"{status}): an experimental module is not validated against clinical "
        "data, whatever this report's verdict."
    )


class ModelValidator:
    """
    Validates the model against clinical benchmark datasets.

    Usage:
        validator = ModelValidator()
        report = validator.validate_ogtt()
        print(report.verdict)
        print(report.as_text())
    """

    def validate(
        self,
        dataset: BenchmarkDataset,
        output,   # SimulationOutput
    ) -> ValidationReport:
        """
        Check a simulation output against a benchmark dataset.

        Args:
            dataset: BenchmarkDataset with reference values
            output:  SimulationOutput from SimulationRunner

        Returns:
            ValidationReport with per-benchmark results and overall verdict
        """
        extracted = self._extract_model_values(output)
        results = []
        notes = []

        for point in dataset.points:
            if point.name not in extracted:
                notes.append(f"Could not extract '{point.name}' from simulation output — skipped.")
                continue
            model_value = extracted[point.name]
            result = point.check(model_value)
            results.append(result)

        n_pass = sum(1 for r in results if r.passed)
        n_fail = sum(1 for r in results if not r.passed)
        verdict = self._determine_verdict(results)
        notes.extend(self._wide_failure_notes(results, verdict))
        notes.extend(self._unscored_point_notes(results))

        return ValidationReport(
            dataset_name=dataset.name,
            scenario_name=output.scenario_name,
            results=results,
            verdict=verdict,
            n_pass=n_pass,
            n_fail=n_fail,
            notes=notes,
        )

    @staticmethod
    def _unscored_point_notes(results: list[BenchmarkCheckResult]) -> list[str]:
        """
        A point with no real literature anchor (BenchmarkPoint.is_scored
        == False) still gets a PASS/MISS label and counts toward "(n/n
        points within tolerance)" in as_text() -- it never changes the
        verdict (see ModelValidator._determine_verdict()), but reading
        "7/7" without qualification understates that one of those seven is
        checking the model against a reference number nobody can actually
        stand behind, not a real clinical comparator. Make that explicit
        every time, regardless of whether it happened to pass or fail.
        """
        unscored = [r for r in results if not r.benchmark.is_scored]
        if not unscored:
            return []
        names = ", ".join(r.name for r in unscored)
        return [
            f"{len(unscored)} benchmark(s) have no real literature anchor "
            f"(pubmed_id=UNVERIFIED) and are informational only, not a "
            f"genuine validation check regardless of their pass/fail "
            f"status: {names}."
        ]

    @staticmethod
    def _wide_failure_notes(results: list[BenchmarkCheckResult], verdict: "ValidationVerdict") -> list[str]:
        """
        Name every cited point behind a MISS verdict, with its deviation,
        so the report says exactly how far the model is from the literature.
        """
        if verdict != ValidationVerdict.MISS:
            return []
        missed = [r for r in results if r.benchmark.is_scored and not r.passed]
        details = ", ".join(f"{r.name} ({r.deviation_percent:+.1f}%)" for r in missed)
        return [
            f"{len(missed)} cited benchmark point(s) missed their tolerance "
            f"band: {details}."
        ]

    def validate_ogtt(self, personalisation_spec=None) -> ValidationReport:
        """
        Convenience method: run the standard 75g OGTT validation.

        This is the primary validation scenario. All model releases
        must pass this before being marked as PRODUCTION maturity.
        """
        from biosim_api import ogtt_75g, SimulationRunner, PersonalizationSpec
        from .benchmark import ogtt_75g_healthy_adults

        spec = personalisation_spec or PersonalizationSpec()
        scenario = ogtt_75g(personalization=spec)
        output = SimulationRunner().run(scenario)
        return self.validate(ogtt_75g_healthy_adults(), output)

    def calibrate_ogtt(
        self,
        n_prior_samples: int = 1000,
        random_seed: int = 42,
    ):
        """
        Bayesian-calibrate p1/p3/KG_secretion against the 75g OGTT
        benchmark (Sampling-Importance-Resampling -- see
        biosim_uncertainty.calibration.BayesianCalibrator's module
        docstring for method and honesty caveats). Returns a
        biosim_uncertainty.CalibrationResult, not a ValidationReport --
        a genuinely different kind of output (a posterior over
        parameters, not a PASS/FAIL check).

        Chosen target and parameters (not an arbitrary pick): OGTT is
        this project's own designated primary validation scenario ("all
        model releases must pass this before PRODUCTION maturity", see
        validate_ogtt() above), with the richest set of benchmark points
        (7, several NORMAL/TIGHT). p1/p3/KG_secretion each already carry
        a real literature point-estimate anchor in PathwayParameters
        (Pacini 1998, Bergman 1981 Fig 6, Rebelato 2018 respectively --
        see pathway.py) and cleanly separate mechanistically: p1/p3
        drive glucose kinetics (peak_glucose, time_to_peak_glucose_min,
        glucose_at_120min), KG_secretion drives insulin secretion onset
        (peak_insulin, time_to_peak_insulin_min) -- each benchmark point
        is informative about a distinguishable subset of the three,
        a reasonable first calibration case. Deliberately NOT a
        comparative drug benchmark (semaglutide/metformin/SGLT2i/steroid)
        -- those are two-simulation differences, harder to attach a
        single-run likelihood to, and at least one (SGLT2i/canagliflozin)
        is documented as EXPECTED to undershoot due to an unmodeled
        mechanism -- calibrating against it would pull parameters toward
        compensating for a confound the model doesn't represent, the
        same curve-fitting-to-a-benchmark anti-pattern this project's
        own CLAUDE.md repeatedly rejects.

        Priors use SILVER-tier CV (0.30, this project's own documented
        expert-judgment default for a literature point-estimate with no
        directly-cited cross-study variance -- see distributions.py) on
        each parameter's current PathwayParameters default -- not tuned
        for this specific calibration.

        Do not treat this posterior as validated on its own -- see
        posterior_predictive_check_bergman_iv() below, which checks it
        against a DIFFERENT benchmark that contributed zero weight here.
        """
        from biosim_ode import PathwayParameters
        from biosim_ode.meal import MealInput
        from biosim_ode.gut_absorption import GutAbsorptionModel
        from biosim_uncertainty import BayesianCalibrator, CALIBRATION_METRIC_EXTRACTORS
        from biosim_uncertainty.distributions import LogNormalDistribution
        from .benchmark import ogtt_75g_healthy_adults

        defaults = PathwayParameters()
        distributions = {
            "p1": LogNormalDistribution.from_confidence_tier(defaults.p1, "silver"),
            "p3": LogNormalDistribution.from_confidence_tier(defaults.p3, "silver"),
            "KG_secretion": LogNormalDistribution.from_confidence_tier(defaults.KG_secretion, "silver"),
        }
        # Matches ogtt_75g()'s own meal exactly (biosim_api/library.py):
        # glycaemic_index=100.0 there routes to this same calibrated
        # absorption preset via MealEvent.to_meal_input() -- built
        # directly here since BayesianCalibrator/MonteCarloEngine work
        # at the bare ODE-engine level (MealInput), not through
        # ScenarioDefinition/MealEvent.
        meal = MealInput(carbs_g=75.0, absorption_model=GutAbsorptionModel.ogtt_calibrated_glucose())

        calibrator = BayesianCalibrator(random_seed=random_seed)
        return calibrator.calibrate(
            distributions=distributions,
            calibration_targets=ogtt_75g_healthy_adults().points,
            metric_extractors=CALIBRATION_METRIC_EXTRACTORS,
            meals=[meal],
            duration_s=14400,
            n_prior_samples=n_prior_samples,
        )

    def posterior_predictive_check_bergman_iv(self, calibration_result) -> ValidationReport:
        """
        Re-simulates the Bergman 1981 IV glucose tolerance test (see
        validate_bergman_iv()) using EVERY posterior parameter draw from
        `calibration_result` (e.g. from calibrate_ogtt()), and checks the
        resulting glucose_return_time_min distribution's MEDIAN against
        the real reference. Bergman IV contributed zero weight to the
        OGTT calibration's likelihood -- a posterior that systematically
        misses a benchmark it never saw is real evidence the likelihood-
        sigma choice or priors are miscalibrated, not noise. This is the
        statistical-honesty safeguard `BayesianCalibrator`'s own module
        docstring calls for: never declare a posterior correct on the
        strength of the benchmark it was fit to.

        Returns a ValidationReport (same shape as every other
        validate_*() method) so this reads exactly like a normal
        benchmark check -- scenario_name records that it's a posterior-
        predictive check, not a point-estimate run, so a caller reading
        the report text isn't misled into thinking this is the same kind
        of check as validate_bergman_iv().
        """
        import dataclasses
        import numpy as np
        from biosim_ode import GlucoseInsulinPathway, ODESystem, PathwayParameters
        from .benchmark import bergman_iv_glucose

        n_posterior = calibration_result.n_posterior_samples
        param_names = list(calibration_result.posterior_parameter_samples.keys())
        return_times: list[float] = []
        n_draw_failures = 0

        for i in range(n_posterior):
            # Regression: this loop used to have no exception handling at
            # all, unlike MonteCarloEngine.propagate()'s equivalent
            # forward-sampling loop (which wraps the same
            # pathway-construction/simulate() calls in try/except and
            # accumulates failure_reasons) -- found during the mandatory
            # security-reviewer pass, see CLAUDE.md. This matters
            # specifically because each posterior draw already succeeded
            # once under the OGTT scenario's own numerics during
            # calibration; re-simulating it here under a DIFFERENT
            # scenario (IV bolus, no meals, different y0/duration) is not
            # guaranteed to succeed just because the OGTT run did -- one
            # posterior draw raising used to crash this entire posterior-
            # predictive check (the project's own documented statistical-
            # honesty safeguard) rather than excluding just that draw.
            try:
                overrides = {
                    name: float(calibration_result.posterior_parameter_samples[name][i])
                    for name in param_names
                }
                params = dataclasses.replace(PathwayParameters(), **overrides)
                pathway = GlucoseInsulinPathway(params=params)
                y0 = pathway.bergman_ivgtt_y0(dose_g=20.0, body_weight_kg=70.0)
                sim = ODESystem(pathway=pathway).simulate(
                    duration_s=14400, meals=[], y0_override=y0, n_output_points=1000,
                )
            except Exception:
                n_draw_failures += 1
                continue
            if not sim.success:
                n_draw_failures += 1
                continue
            G = sim.species["G"]
            t_min = sim.times_s / 60.0
            Gb_uM = params.Gb * 1000.0
            below_threshold = np.where(G <= Gb_uM * 1.10)[0]
            if len(below_threshold) > 0:
                return_times.append(float(t_min[below_threshold[0]]))

        dataset = bergman_iv_glucose()
        results = []
        notes = [
            "Posterior-predictive check: every one of the below-listed "
            f"{n_posterior} posterior parameter draws was re-simulated "
            "through the Bergman IV scenario independently -- this "
            "benchmark contributed zero weight to the OGTT calibration "
            "that produced this posterior.",
        ]
        if n_draw_failures > 0:
            notes.append(
                f"{n_draw_failures}/{n_posterior} posterior draws raised an "
                "exception or failed to solve under the Bergman IV scenario "
                "and were excluded (not silently ignored) -- a genuinely "
                "unstable posterior draw, or a real bug in this check, "
                "should be visible here rather than just shrinking "
                "return_times' denominator with no explanation."
            )
        if return_times:
            median_return_time = float(np.median(return_times))
            point = next(p for p in dataset.points if p.name == "glucose_return_time_min")
            results.append(point.check(median_return_time))
            notes.append(
                f"{len(return_times)}/{n_posterior} posterior draws produced "
                "a valid return time; median used for the check above."
            )
        else:
            notes.append("No posterior draw produced a valid glucose_return_time_min.")

        n_pass = sum(1 for r in results if r.passed)
        n_fail = sum(1 for r in results if not r.passed)
        verdict = self._determine_verdict(results)
        notes.extend(self._wide_failure_notes(results, verdict))
        return ValidationReport(
            dataset_name=dataset.name,
            scenario_name="bergman_iv_posterior_predictive_check",
            results=results,
            verdict=verdict,
            n_pass=n_pass,
            n_fail=n_fail,
            notes=notes,
        )

    def validate_bergman_iv(self, dose_g: float = 20.0, body_weight_kg: float = 70.0) -> ValidationReport:
        """
        Convenience method: run the Bergman 1981 IVGTT validation.

        Unlike validate_ogtt(), this runs at the ODE-engine level directly
        (no SimulationRunner/ScenarioDefinition -- there is no IV-bolus
        event type at the scenario-API layer, see CLAUDE.md VALIDATION
        STATUS) via GlucoseInsulinPathway.bergman_ivgtt_y0() +
        ODESystem.simulate(y0_override=...). No meals, no personalisation --
        this is a narrower, ODE-model-only check than validate_ogtt().
        """
        import numpy as np
        from biosim_ode import GlucoseInsulinPathway, ODESystem
        from .benchmark import bergman_iv_glucose

        pathway = GlucoseInsulinPathway()
        y0 = pathway.bergman_ivgtt_y0(dose_g=dose_g, body_weight_kg=body_weight_kg)
        sim = ODESystem(pathway=pathway).simulate(
            duration_s=14400, meals=[], y0_override=y0, n_output_points=2000,
        )

        dataset = bergman_iv_glucose()
        results = []
        notes = []
        if sim.success:
            G = sim.species["G"]
            t_min = sim.times_s / 60.0
            Gb_uM = pathway.params.Gb * 1000.0
            below_threshold = np.where(G <= Gb_uM * 1.10)[0]
            if len(below_threshold) > 0:
                return_time_min = float(t_min[below_threshold[0]])
                point = next(p for p in dataset.points if p.name == "glucose_return_time_min")
                results.append(point.check(return_time_min))
            else:
                notes.append("Glucose never returned within 10% of baseline in the simulated window.")
        else:
            notes.append(f"ODE solver failed: {sim.message}")

        n_pass = sum(1 for r in results if r.passed)
        n_fail = sum(1 for r in results if not r.passed)
        verdict = self._determine_verdict(results)
        notes.extend(self._wide_failure_notes(results, verdict))
        notes.extend(self._unscored_point_notes(results))
        return ValidationReport(
            dataset_name=dataset.name,
            scenario_name="Bergman_IV_glucose_ODE_only",
            results=results,
            verdict=verdict,
            n_pass=n_pass,
            n_fail=n_fail,
            notes=notes,
        )

    def validate_semaglutide(self, dose_mg: float = 1.0, carbs_g: float = 75.0) -> ValidationReport:
        """
        Convenience method: compare semaglutide (steady-state) vs. placebo,
        same meal, against Hjerpsted et al. 2018 (PMID:28941314).

        This is a COMPARATIVE benchmark (drug vs. placebo iAUC difference),
        unlike validate_ogtt()/validate_bergman_iv()'s single-simulation
        absolute-value checks -- runs two full SimulationRunner scenarios
        internally with a 5h duration (the trial's iAUC window) rather than
        semaglutide_scenario()'s own default 4h duration, which is too
        short for this specific comparison.

        First validation of a clinically significant drug pathway in this
        codebase -- see CLAUDE.md "DONE -- Built the first validation
        benchmark for a clinically significant drug pathway". As of this
        writing the model OVERSHOOTS the real effect substantially
        (documented, not silently accepted) -- this method exists to
        make that precisely measurable and regression-testable, not to
        force a pass.
        """
        from biosim_api import (
            DrugEvent, MealEvent, PersonalizationSpec, ScenarioDefinition, SimulationRunner,
        )
        from .benchmark import semaglutide_postprandial_glucose

        FIVE_HOUR_DURATION_S = 18000.0
        runner = SimulationRunner()

        placebo_scenario = ScenarioDefinition(
            name="semaglutide_benchmark_placebo",
            duration_s=FIVE_HOUR_DURATION_S,
            meal_events=[MealEvent(time_s=0.0, carbs_g=carbs_g)],
            personalization=PersonalizationSpec(),
        )
        drug_scenario = ScenarioDefinition(
            name="semaglutide_benchmark_drug",
            duration_s=FIVE_HOUR_DURATION_S,
            meal_events=[MealEvent(time_s=0.0, carbs_g=carbs_g)],
            drug_events=[DrugEvent(time_s=0.0, drug="semaglutide", dose_mg=dose_mg)],
            personalization=PersonalizationSpec(),
        )

        placebo_out = runner.run(placebo_scenario)
        drug_out = runner.run(drug_scenario)

        dataset = semaglutide_postprandial_glucose()
        results = []
        notes = []
        if placebo_out.success and drug_out.success:
            reduction_mmol_h_per_L = (
                placebo_out.expert_metrics.glucose_iauc_5h_mmol_h_per_L
                - drug_out.expert_metrics.glucose_iauc_5h_mmol_h_per_L
            )
            point = next(
                p for p in dataset.points
                if p.name == "glucose_iauc_5h_reduction_mmol_h_per_L"
            )
            results.append(point.check(reduction_mmol_h_per_L))
        else:
            if not placebo_out.success:
                notes.append(f"Placebo scenario ODE solver failed: {placebo_out.warnings}")
            if not drug_out.success:
                notes.append(f"Drug scenario ODE solver failed: {drug_out.warnings}")

        n_pass = sum(1 for r in results if r.passed)
        n_fail = sum(1 for r in results if not r.passed)
        verdict = self._determine_verdict(results)
        notes.extend(self._wide_failure_notes(results, verdict))
        notes.extend(self._unscored_point_notes(results))
        return ValidationReport(
            dataset_name=dataset.name,
            module_label=_drug_module_label("semaglutide"),
            scenario_name="semaglutide_vs_placebo_5h",
            results=results,
            verdict=verdict,
            n_pass=n_pass,
            n_fail=n_fail,
            notes=notes,
        )

    def validate_steroid_meal_t2d(self, prednisone_mg: float = 40.0, carbs_g: float = 75.0) -> ValidationReport:
        """
        Convenience method: compare prednisone 40mg vs. no steroid, same
        meal, in a type 2 diabetes population, against Shah et al. 2020
        (PMID:32369480).

        Deliberately targets the T2D population, not the platform's
        healthy default: real evidence (see CLAUDE.md's steroid/prednisone
        entry, right after the semaglutide one) shows the acute
        glucose-raising effect of short-term prednisone is concentrated in
        T2D (reduced beta-cell reserve) and was NOT statistically
        significant in healthy subjects across two independent studies
        (n=37 and n=5) -- benchmarking against the healthy population
        would be checking the model against an effect real data says is
        near-null there. Same comparative-benchmark shape as
        validate_semaglutide() -- builds both arms explicitly (T2D + same
        meal, only the steroid differs) rather than reusing
        type2_diabetes_meal_scenario()/steroid_meal()'s own independent
        defaults, which could reintroduce the meal-size-mismatch class of
        confound already caught once this session.

        Dose-to-meal interval, now genuinely load-bearing (see CLAUDE.md's
        glucocorticoid PK/PD milestone entry): now that CortisolPathway's
        glycaemic effect is real PD-timing-gated (near-zero for the first
        ~2h, peaking ~5.7h post-dose), the gap between the prednisone dose
        and the OGTT is no longer a free/inert choice -- it has to match
        Shah et al. 2020's own protocol or this becomes exactly the
        "different conditions labeled placebo vs. drug" confound class the
        semaglutide benchmark hit on its first attempt. Shah's own Methods
        text (fetched and re-verified against full text this session, not
        assumed from the abstract): "Two hours later, another OGTT and
        blood samples were collected" after the third prednisone dose --
        the OGTT was performed 2h post-dose, not simultaneously with it.
        DrugEvent.hours_before_simulation=2.0 places the dose 2h before
        the meal (t=0, matching every other benchmark's meal-timing
        convention) rather than moving the meal later, so
        expert_metrics.glucose_auc_uM_s still integrates from the meal's
        own t=0.
        """
        from biosim_api import (
            DrugEvent, MealEvent, PersonalizationSpec, ScenarioDefinition, SimulationRunner,
        )
        from biosim_physiology import DiseaseState, PhysiologicalProfile
        from .benchmark import steroid_t2d_postprandial_glucose

        DURATION_S = 14400.0
        DOSE_TO_MEAL_INTERVAL_HOURS = 2.0   # Shah et al. 2020's own protocol, see docstring above
        runner = SimulationRunner()
        t2d_profile = PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES])

        no_steroid_scenario = ScenarioDefinition(
            name="steroid_benchmark_t2d_no_steroid",
            duration_s=DURATION_S,
            meal_events=[MealEvent(time_s=0.0, carbs_g=carbs_g)],
            personalization=PersonalizationSpec(physiological_profile=t2d_profile),
        )
        steroid_scenario = ScenarioDefinition(
            name="steroid_benchmark_t2d_with_steroid",
            duration_s=DURATION_S,
            meal_events=[MealEvent(time_s=0.0, carbs_g=carbs_g)],
            drug_events=[DrugEvent(
                time_s=0.0, drug="glucocorticoid", dose_mg=prednisone_mg,
                hours_before_simulation=DOSE_TO_MEAL_INTERVAL_HOURS,
            )],
            personalization=PersonalizationSpec(physiological_profile=t2d_profile),
        )

        no_steroid_out = runner.run(no_steroid_scenario)
        steroid_out = runner.run(steroid_scenario)

        dataset = steroid_t2d_postprandial_glucose()
        results = []
        notes = []
        if no_steroid_out.success and steroid_out.success:
            baseline_auc = no_steroid_out.expert_metrics.glucose_auc_uM_s
            steroid_auc = steroid_out.expert_metrics.glucose_auc_uM_s
            percent_increase = (steroid_auc - baseline_auc) / baseline_auc * 100.0
            point = next(
                p for p in dataset.points if p.name == "glucose_auc_percent_increase_t2d"
            )
            results.append(point.check(percent_increase))
        else:
            if not no_steroid_out.success:
                notes.append(f"No-steroid scenario ODE solver failed: {no_steroid_out.warnings}")
            if not steroid_out.success:
                notes.append(f"Steroid scenario ODE solver failed: {steroid_out.warnings}")

        n_pass = sum(1 for r in results if r.passed)
        n_fail = sum(1 for r in results if not r.passed)
        verdict = self._determine_verdict(results)
        notes.extend(self._wide_failure_notes(results, verdict))
        notes.extend(self._unscored_point_notes(results))
        return ValidationReport(
            dataset_name=dataset.name,
            module_label=_drug_module_label("glucocorticoid"),
            scenario_name="steroid_t2d_vs_no_steroid",
            results=results,
            verdict=verdict,
            n_pass=n_pass,
            n_fail=n_fail,
            notes=notes,
        )

    def validate_metformin(self, dose_mg: float = 2000.0, carbs_g: float = 75.0) -> ValidationReport:
        """
        Convenience method: compare chronic metformin vs. placebo, same
        4h mixed meal test, in a type 2 diabetes population, against
        Hansen et al. 2024 (PMID:39049802).

        This is the metformin acute-meal-challenge benchmark CLAUDE.md's
        VALIDATION STATUS section already flagged as the explicit next
        step: a directly comparable standard-oral-meal glucose iAUC study
        (unlike Borg et al. 2025's intraduodenal-infusion tracer
        methodology, which measures a different physical quantity than
        this platform's MealInput/Ra(t) represents). Same
        comparative-benchmark shape as validate_semaglutide()/
        validate_steroid_meal_t2d() -- builds both arms explicitly (T2D +
        same meal, only metformin differs), and uses the 4h iAUC window
        (glucose_iauc_4h_mmol_min_per_L) matching the trial's actual 4h
        mixed-meal-test duration, not the 5h window used for semaglutide
        (that mismatch is exactly the "different conditions" confound
        this codebase has caught and rejected before -- see the
        semaglutide entry's own confounded-first-attempt history).

        Metformin's existing modifiers (metformin_modifiers() in
        drug_interactions.py: p3/p1/Gb static effects) are dose-only, not
        gated behind DrugDose.is_chronic -- this model already treats
        every metformin dose as the chronic, steady-state effect the
        trial's 14-day dosing period measures, so no separate "chronic"
        wiring is needed here beyond passing the dose.
        """
        from biosim_api import (
            DrugEvent, MealEvent, PersonalizationSpec, ScenarioDefinition, SimulationRunner,
        )
        from biosim_physiology import DiseaseState, PhysiologicalProfile
        from .benchmark import metformin_postprandial_glucose_t2d

        FOUR_HOUR_DURATION_S = 14400.0
        runner = SimulationRunner()
        t2d_profile = PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES])

        placebo_scenario = ScenarioDefinition(
            name="metformin_benchmark_placebo",
            duration_s=FOUR_HOUR_DURATION_S,
            meal_events=[MealEvent(time_s=0.0, carbs_g=carbs_g)],
            personalization=PersonalizationSpec(physiological_profile=t2d_profile),
        )
        drug_scenario = ScenarioDefinition(
            name="metformin_benchmark_drug",
            duration_s=FOUR_HOUR_DURATION_S,
            meal_events=[MealEvent(time_s=0.0, carbs_g=carbs_g)],
            drug_events=[DrugEvent(time_s=0.0, drug="metformin", dose_mg=dose_mg)],
            personalization=PersonalizationSpec(physiological_profile=t2d_profile),
        )

        placebo_out = runner.run(placebo_scenario)
        drug_out = runner.run(drug_scenario)

        dataset = metformin_postprandial_glucose_t2d()
        results = []
        notes = []
        if placebo_out.success and drug_out.success:
            reduction_mmol_min_per_L = (
                placebo_out.expert_metrics.glucose_iauc_4h_mmol_min_per_L
                - drug_out.expert_metrics.glucose_iauc_4h_mmol_min_per_L
            )
            point = next(
                p for p in dataset.points
                if p.name == "glucose_iauc_4h_reduction_mmol_min_per_L_t2d"
            )
            results.append(point.check(reduction_mmol_min_per_L))
        else:
            if not placebo_out.success:
                notes.append(f"Placebo scenario ODE solver failed: {placebo_out.warnings}")
            if not drug_out.success:
                notes.append(f"Drug scenario ODE solver failed: {drug_out.warnings}")

        n_pass = sum(1 for r in results if r.passed)
        n_fail = sum(1 for r in results if not r.passed)
        verdict = self._determine_verdict(results)
        notes.extend(self._wide_failure_notes(results, verdict))
        notes.extend(self._unscored_point_notes(results))
        return ValidationReport(
            dataset_name=dataset.name,
            module_label=_drug_module_label("metformin"),
            scenario_name="metformin_vs_placebo_4h_t2d",
            results=results,
            verdict=verdict,
            n_pass=n_pass,
            n_fail=n_fail,
            notes=notes,
        )

    def validate_sglt2i(self, dose_mg: float = 300.0, carbs_g: float = 60.0) -> ValidationReport:
        """
        Convenience method: compare canagliflozin vs. placebo, same meal,
        healthy subjects, against Polidori et al. 2013 (PMID:23412078).

        Same comparative-benchmark shape as validate_semaglutide()/
        validate_metformin() -- uses a 2h AUC window (matching the trial's
        own AUC(0-2h) endpoint) for both glucose and insulin, computed
        directly from the raw simulation trajectory (this platform's
        ExpertMetrics has no dedicated 2h-window field, unlike the 4h/5h
        iAUC fields the semaglutide/metformin benchmarks use -- adding one
        for a single benchmark's use was judged not worth the schema
        churn; the 4h/5h fields exist because multiple call sites use
        them).

        This platform does not build a dedicated canagliflozin preset (per
        computational-biologist's scoping recommendation -- canagliflozin's
        Tm/GFR/threshold values are not independently verified, and its
        SGLT1 confound below means any calibration against it would fold
        a non-renal mechanism into the renal-only model's parameters). The
        drug scenario instead uses Drug.DAPAGLIFLOZIN, the closest
        verified renal-mechanism preset available -- dose_mg is accepted
        for API consistency with validate_metformin()/validate_semaglutide()
        but is not actually read by RenalClearanceParameters.
        for_dapagliflozin() (a fixed, duration-only preset -- see
        renal_pathway.py), so it has no effect on the result.

        EXPECTED TO UNDERSHOOT: canagliflozin's real effect at this dose
        is NOT purely renal -- see
        canagliflozin_postprandial_glucose_and_insulin()'s docstring for
        the confirmed SGLT1/intestinal-absorption confound. This
        benchmark exists to make that undershoot precisely measurable,
        not to force a pass -- same honesty standard as
        validate_semaglutide()'s documented overshoot.
        """
        from biosim_api import (
            DrugEvent, MealEvent, PersonalizationSpec, ScenarioDefinition, SimulationRunner,
        )
        from .benchmark import canagliflozin_postprandial_glucose_and_insulin

        TWO_HOUR_S = 7200.0
        runner = SimulationRunner()

        placebo_scenario = ScenarioDefinition(
            name="sglt2i_benchmark_placebo",
            duration_s=TWO_HOUR_S,
            meal_events=[MealEvent(time_s=0.0, carbs_g=carbs_g)],
            personalization=PersonalizationSpec(),
        )
        drug_scenario = ScenarioDefinition(
            name="sglt2i_benchmark_drug",
            duration_s=TWO_HOUR_S,
            meal_events=[MealEvent(time_s=0.0, carbs_g=carbs_g)],
            drug_events=[DrugEvent(time_s=0.0, drug="dapagliflozin", dose_mg=dose_mg)],
            personalization=PersonalizationSpec(),
        )

        placebo_out = runner.run(placebo_scenario)
        drug_out = runner.run(drug_scenario)

        dataset = canagliflozin_postprandial_glucose_and_insulin()
        results = []
        notes = []
        if placebo_out.success and drug_out.success:
            placebo_G = placebo_out.simulation.species["G"]
            drug_G = drug_out.simulation.species["G"]
            placebo_I = placebo_out.simulation.species["I"]
            drug_I = drug_out.simulation.species["I"]
            times = placebo_out.simulation.times_s

            auc_placebo_G = float(np.trapezoid(placebo_G, times))
            auc_drug_G = float(np.trapezoid(drug_G, times))
            percent_reduction_G = (auc_placebo_G - auc_drug_G) / auc_placebo_G * 100.0

            auc_placebo_I = float(np.trapezoid(placebo_I, times))
            auc_drug_I = float(np.trapezoid(drug_I, times))
            percent_reduction_I = (auc_placebo_I - auc_drug_I) / auc_placebo_I * 100.0

            g_point = next(p for p in dataset.points if p.name == "glucose_auc_2h_percent_reduction")
            i_point = next(p for p in dataset.points if p.name == "insulin_auc_2h_percent_reduction")
            results.append(g_point.check(percent_reduction_G))
            results.append(i_point.check(percent_reduction_I))
        else:
            if not placebo_out.success:
                notes.append(f"Placebo scenario ODE solver failed: {placebo_out.warnings}")
            if not drug_out.success:
                notes.append(f"Drug scenario ODE solver failed: {drug_out.warnings}")

        n_pass = sum(1 for r in results if r.passed)
        n_fail = sum(1 for r in results if not r.passed)
        verdict = self._determine_verdict(results)
        notes.extend(self._wide_failure_notes(results, verdict))
        notes.extend(self._unscored_point_notes(results))
        notes.append(
            "Model is expected to UNDERSHOOT this benchmark's magnitude -- "
            "canagliflozin's real effect includes an intestinal SGLT1 "
            "mechanism this platform's renal-only RenalPathway does not "
            "model (see canagliflozin_postprandial_glucose_and_insulin()'s "
            "docstring). A WIDE-tolerance miss here reflects that known, "
            "documented gap, not an unexplained model error."
        )
        return ValidationReport(
            dataset_name=dataset.name,
            module_label=_drug_module_label("dapagliflozin"),
            scenario_name="sglt2i_vs_placebo_2h",
            results=results,
            verdict=verdict,
            n_pass=n_pass,
            n_fail=n_fail,
            notes=notes,
        )

    def validate_protein_coingestion(self, carbs_g: float = 50.0) -> ValidationReport:
        """
        50 g glucose alone vs. with 25 g / 50 g pea protein, healthy
        subjects, against Thondre et al. 2021 (PMID:33515092) -- see
        pea_protein_coingestion_50g_glucose() for the source and why the
        model is EXPECTED to undershoot (protein-stimulated insulin
        secretion is modelled; gastric-emptying slowing is not). iAUCs over 180 min from the raw trajectory, area above each
        run's own t=0 value (negative increments ignored), matching the
        trial's iAUC180 endpoint.
        """
        from biosim_api import MealEvent, ScenarioDefinition, SimulationRunner
        from .benchmark import pea_protein_coingestion_50g_glucose

        THREE_HOURS_S = 10800.0
        runner = SimulationRunner()

        def run(protein_g: float):
            return runner.run(ScenarioDefinition(
                name=f"protein_benchmark_{protein_g:g}g",
                duration_s=THREE_HOURS_S,
                # A glucose drink: route to the calibrated OGTT absorption,
                # exactly as ogtt_75g() does.
                meal_events=[MealEvent(
                    time_s=0.0, carbs_g=carbs_g, protein_g=protein_g, glycaemic_index=100.0,
                )],
            ))

        def iauc(output, species: str) -> float:
            values = output.simulation.species[species]
            above = np.clip(values - values[0], 0.0, None)
            return float(np.trapezoid(above, output.simulation.times_s))

        dataset = pea_protein_coingestion_50g_glucose()
        points = {p.name: p for p in dataset.points}
        outputs = {grams: run(float(grams)) for grams in (0, 25, 50)}
        results = []
        notes = []
        failed = [g for g, out in outputs.items() if not out.success]
        if failed:
            notes.extend(
                f"{g} g protein scenario ODE solver failed: {outputs[g].warnings}" for g in failed
            )
        else:
            control_G, control_I = iauc(outputs[0], "G"), iauc(outputs[0], "I")
            for grams in (25, 50):
                g_reduction = (control_G - iauc(outputs[grams], "G")) / control_G * 100.0
                i_increase = (iauc(outputs[grams], "I") - control_I) / control_I * 100.0
                results.append(points[f"glucose_iauc_3h_percent_reduction_{grams}g_protein"].check(g_reduction))
                results.append(points[f"insulin_iauc_3h_percent_increase_{grams}g_protein"].check(i_increase))

        n_pass = sum(1 for r in results if r.passed)
        n_fail = sum(1 for r in results if not r.passed)
        verdict = self._determine_verdict(results)
        notes.extend(self._wide_failure_notes(results, verdict))
        notes.extend(self._unscored_point_notes(results))
        notes.append(
            "Model is expected to UNDERSHOOT -- most sharply for glucose, "
            "which protein plausibly lowers in large part by slowing gastric "
            "emptying (not measured in this study -- an inference); "
            "NutrientSecretionPathway models only protein-stimulated insulin "
            "secretion, calibrated on other protein sources. Protein's "
            "glucagon response is also unmodelled, but it works the OTHER "
            "way (it would make the model over-predict glucose lowering), "
            "so the gastric-emptying gap is if anything larger than shown. A "
            "WIDE-tolerance miss here measures that documented gap; it is not "
            "a target to tune secretion against."
        )
        return ValidationReport(
            dataset_name=dataset.name,
            scenario_name="pea_protein_vs_glucose_alone_3h",
            results=results,
            verdict=verdict,
            n_pass=n_pass,
            n_fail=n_fail,
            notes=notes,
        )

    def validate_t2d_meal_insulin(self, carbs_g: float = 75.0) -> ValidationReport:
        """
        Default T2D preset vs. a healthy profile (same age/BMI, only the
        disease state differs), same 75 g glucose drink, 6 h, against Basu et
        al. 2009 (PMID:19196896) -- see t2d_vs_nondiabetic_meal_insulin() for
        the source, the disclosed meal/cohort mismatches and why every point
        is expected to miss. Increments are measured above each run's own
        t=0 value with negative increments ignored -- the same convention
        as validate_protein_coingestion, not claimed to match Basu's own
        (probably signed) integral; no sample falls below t=0 in either arm.
        """
        from biosim_api import MealEvent, PersonalizationSpec, ScenarioDefinition, SimulationRunner
        from biosim_physiology import DiseaseState, PhysiologicalProfile
        from .benchmark import t2d_vs_nondiabetic_meal_insulin

        SIX_HOURS_S = 21600.0
        TWO_HOURS_S = 7200.0
        runner = SimulationRunner()

        def run(state: DiseaseState):
            return runner.run(ScenarioDefinition(
                name=f"t2d_insulin_benchmark_{state.value}",
                duration_s=SIX_HOURS_S,
                # A glucose drink: calibrated OGTT absorption, as ogtt_75g().
                meal_events=[MealEvent(time_s=0.0, carbs_g=carbs_g, glycaemic_index=100.0)],
                personalization=PersonalizationSpec(
                    physiological_profile=PhysiologicalProfile(disease_states=[state]),
                ),
            ))

        def iauc(output, species: str, limit_s: float) -> float:
            times = output.simulation.times_s
            values = output.simulation.species[species]
            above = np.clip(values - values[0], 0.0, None)
            mask = times <= limit_s
            return float(np.trapezoid(above[mask], times[mask]))

        def peak(output, species: str) -> tuple[float, float]:
            values = output.simulation.species[species]
            idx = int(np.argmax(values))
            return float(values[idx]), float(output.simulation.times_s[idx])

        dataset = t2d_vs_nondiabetic_meal_insulin()
        points = {p.name: p for p in dataset.points}
        outputs = {
            "healthy": run(DiseaseState.HEALTHY),
            "t2d": run(DiseaseState.TYPE_2_DIABETES),
        }
        results = []
        notes = []
        failed = [k for k, out in outputs.items() if not out.success]
        if failed:
            notes.extend(
                f"{k} arm ODE solver failed: {outputs[k].warnings}" for k in failed
            )
        else:
            h, t = outputs["healthy"], outputs["t2d"]
            (h_peak, h_t), (t_peak, t_t) = peak(h, "I"), peak(t, "I")
            measured = {
                "insulin_peak_ratio_t2d_vs_nondiabetic": t_peak / h_peak,
                "insulin_iauc_2h_ratio_t2d_vs_nondiabetic":
                    iauc(t, "I", TWO_HOURS_S) / iauc(h, "I", TWO_HOURS_S),
                "insulin_iauc_6h_ratio_t2d_vs_nondiabetic":
                    iauc(t, "I", SIX_HOURS_S) / iauc(h, "I", SIX_HOURS_S),
                "insulin_peak_delay_min_t2d_minus_nondiabetic": (t_t - h_t) / 60.0,
                "glucose_iauc_6h_ratio_t2d_vs_nondiabetic":
                    iauc(t, "G", SIX_HOURS_S) / iauc(h, "G", SIX_HOURS_S),
            }
            results.extend(points[name].check(value) for name, value in measured.items())

        n_pass = sum(1 for r in results if r.passed)
        n_fail = sum(1 for r in results if not r.passed)
        verdict = self._determine_verdict(results)
        notes.extend(self._wide_failure_notes(results, verdict))
        notes.extend(self._unscored_point_notes(results))
        if results:
            notes.append(
                "MASTER_PLAN row 1.46 (open question): this measures how the "
                "default T2D preset's postprandial insulin response compares "
                "with real, more hyperglycaemic early T2D. A WIDE-tolerance "
                "miss is a measurement, not a target to tune secretion "
                "against (changing it needs its own cited, double-counting-"
                "checked milestone)."
            )
        return ValidationReport(
            dataset_name=dataset.name,
            scenario_name="t2d_vs_healthy_75g_glucose_6h",
            results=results,
            verdict=verdict,
            n_pass=n_pass,
            n_fail=n_fail,
            notes=notes,
        )

    def validate_caffeine_ogtt(self, mg_per_kg: float = 5.0, carbs_g: float = 75.0) -> ValidationReport:
        """
        Caffeine (5 mg/kg, taken 1 h before) vs. no caffeine, same 75 g glucose
        drink, in a healthy profile and the default T2D profile, against Graham
        et al. 2001 (PMID:11478588) and Robinson et al. 2004 (PMID:15465742) --
        see caffeine_ogtt_insulin() for the source details, the disclosed
        mismatches and why both points are expected to miss with the wrong
        sign. Total insulin AUC over 3 h from the raw trajectory.
        """
        from biosim_api import (
            DrugEvent, MealEvent, PersonalizationSpec, ScenarioDefinition, SimulationRunner,
        )
        from biosim_physiology import DiseaseState, PhysiologicalProfile
        from .benchmark import caffeine_ogtt_insulin

        THREE_HOURS_S = 10800.0
        runner = SimulationRunner()

        def run(state: DiseaseState, caffeine: bool):
            profile = PhysiologicalProfile(disease_states=[state])
            drugs = [DrugEvent(
                time_s=0.0, drug="caffeine", dose_mg=mg_per_kg * profile.body_weight_kg,
                hours_before_simulation=1.0,
            )] if caffeine else []
            return runner.run(ScenarioDefinition(
                name=f"caffeine_benchmark_{state.value}_{'caffeine' if caffeine else 'placebo'}",
                duration_s=THREE_HOURS_S,
                # A glucose drink: calibrated OGTT absorption, as ogtt_75g().
                meal_events=[MealEvent(time_s=0.0, carbs_g=carbs_g, glycaemic_index=100.0)],
                drug_events=drugs,
                personalization=PersonalizationSpec(physiological_profile=profile),
            ))

        def total_auc(output, species: str) -> float:
            return float(np.trapezoid(
                output.simulation.species[species], output.simulation.times_s,
            ))

        dataset = caffeine_ogtt_insulin()
        points = {p.name: p for p in dataset.points}
        arms = {
            (state, caffeine): run(state, caffeine)
            for state in (DiseaseState.HEALTHY, DiseaseState.TYPE_2_DIABETES)
            for caffeine in (False, True)
        }
        results = []
        notes = []
        failed = [k for k, out in arms.items() if not out.success]
        if failed:
            notes.extend(
                f"{state.value}/{'caffeine' if caffeine else 'placebo'} arm ODE solver "
                f"failed: {arms[(state, caffeine)].warnings}"
                for state, caffeine in failed
            )
        else:
            for name, state in (
                ("insulin_auc_percent_increase_healthy_caffeine", DiseaseState.HEALTHY),
                ("insulin_auc_percent_increase_t2d_caffeine", DiseaseState.TYPE_2_DIABETES),
            ):
                placebo = total_auc(arms[(state, False)], "I")
                caffeine = total_auc(arms[(state, True)], "I")
                results.append(points[name].check((caffeine / placebo - 1.0) * 100.0))

        n_pass = sum(1 for r in results if r.passed)
        n_fail = sum(1 for r in results if not r.passed)
        verdict = self._determine_verdict(results)
        notes.extend(self._wide_failure_notes(results, verdict))
        notes.extend(self._unscored_point_notes(results))
        if results and any(not r.passed for r in results):
            notes.append(
                "Real caffeine RAISES OGTT insulin (+60% healthy, +25% T2D). The "
                "model's only caffeine effect is a small, cited p3 (insulin-"
                "sensitivity) reduction, which raises insulin by about 1%: the "
                "insulin-raising mechanism is not modelled. (An uncited +2% Gb "
                "term, which made the sign wrong, was removed on 2026-11-20.) "
                "A miss here is a measurement, not a target to tune against."
            )
        return ValidationReport(
            dataset_name=dataset.name,
            module_label=_drug_module_label("caffeine"),
            scenario_name="caffeine_vs_placebo_75g_ogtt_3h",
            results=results,
            verdict=verdict,
            n_pass=n_pass,
            n_fail=n_fail,
            notes=notes,
        )

    def validate_all(self) -> dict[str, ValidationReport]:
        """
        Run every pre-registered clinical validation benchmark and return
        all reports keyed by dataset name.

        Convenience entry point added after the fourth benchmark
        (`validate_ogtt`, `validate_bergman_iv`, `validate_semaglutide`,
        `validate_steroid_meal_t2d`) accumulated with no single place to
        run them all together -- useful for a CI/publication-readiness
        check or a quick "is anything broken" sweep after a change to
        shared ODE machinery. Each method already has its own default
        arguments (population, dose, meal size) matching its source
        study -- this does not override any of them.

        Does NOT introduce a new combined verdict across benchmarks --
        each dataset's PASS/WARN/FAIL stays independent (e.g. a semaglutide
        WIDE-tolerance miss should never be conflated with an OGTT TIGHT
        failure, which are different severities). Read
        `report[name].verdict` per dataset, same as calling each method
        individually.
        """
        return {
            "OGTT_75g_healthy": self.validate_ogtt(),
            "Bergman_IV_glucose": self.validate_bergman_iv(),
            "Semaglutide_1mg_postprandial": self.validate_semaglutide(),
            "Prednisone_40mg_T2D_postprandial": self.validate_steroid_meal_t2d(),
            "Metformin_T2D_postprandial": self.validate_metformin(),
            "Canagliflozin_300mg_postprandial": self.validate_sglt2i(),
            "Pea_protein_coingestion_50g_glucose": self.validate_protein_coingestion(),
            "T2D_vs_nondiabetic_meal_insulin": self.validate_t2d_meal_insulin(),
            "Caffeine_5mgkg_OGTT_insulin": self.validate_caffeine_ogtt(),
        }

    def validate_heldout(self) -> dict[str, ValidationReport]:
        """
        Held-out group (finish-plan B2): published data the model was never
        fitted to (heldout.py), run on default parameters with nothing
        adjusted. Headline band STRICT (+/-20%); the verdict uses the same
        PASS/MISS/NO DATA logic as every other report. Kept separate from
        validate_all() so the held-out result is always read on its own.
        """
        from .heldout import HELDOUT_RUNS, heldout_datasets

        return self._run_unchanged_group(
            heldout_datasets(), HELDOUT_RUNS,
            "Held-out dataset: never used to fit or calibrate the model; "
            "model run unchanged; headline band +/-20%.",
        )

    def validate_first_contact(self) -> dict[str, ValidationReport]:
        """
        H2 first-contact group (first_contact.py): three datasets never used by
        the project, pre-registered in docs/VALIDATION.md section 8 before the
        one run. Same PASS/MISS/NO DATA logic, headline band +/-20%.
        """
        from .first_contact import FIRST_CONTACT_RUNS, first_contact_datasets

        return self._run_unchanged_group(
            first_contact_datasets(), FIRST_CONTACT_RUNS,
            "First-contact dataset (H2): pre-registered, never used by the project "
            "before; model run once, unchanged; headline band +/-20%.",
        )

    def _run_unchanged_group(
        self, datasets: list[BenchmarkDataset], runs: dict[str, Callable], label: str,
    ) -> dict[str, ValidationReport]:
        from biosim_api import SimulationRunner

        runner = SimulationRunner()
        reports: dict[str, ValidationReport] = {}
        for dataset in datasets:
            values = runs[dataset.name](runner)
            # Direct lookup: a renamed key must fail loudly, never drop a point.
            results = [p.check(values[p.name]) for p in dataset.points]
            verdict = self._determine_verdict(results)
            notes = self._wide_failure_notes(results, verdict)
            reports[dataset.name] = ValidationReport(
                dataset_name=dataset.name,
                module_label=label,
                scenario_name=dataset.description,
                results=results,
                verdict=verdict,
                n_pass=sum(1 for r in results if r.passed),
                n_fail=sum(1 for r in results if not r.passed),
                notes=notes,
            )
        return reports

    def _extract_model_values(self, output) -> dict[str, float]:
        """Extract named benchmark quantities from SimulationOutput."""
        m = output.expert_metrics
        sim = output.simulation
        extracted = {}

        extracted["fasting_glucose"] = m.glucose_baseline_uM / 1000.0     # uM → mM
        extracted["peak_glucose"] = m.glucose_peak_uM / 1000.0
        extracted["time_to_peak_glucose_min"] = m.glucose_peak_time_s / 60.0
        extracted["peak_insulin"] = m.insulin_peak_uM * 1e6                # uM → pM
        extracted["time_to_peak_insulin_min"] = m.insulin_peak_time_s / 60.0
        extracted["mtor_activation"] = m.mtor_peak

        if m.time_to_glucose_return_s:
            extracted["glucose_return_time_min"] = m.time_to_glucose_return_s / 60.0

        # Glucose at 120 min
        if sim.success and len(sim.times_s) > 0:
            G_at_120 = float(np.interp(7200.0, sim.times_s, sim.species["G"]))
            extracted["glucose_at_120min"] = G_at_120 / 1000.0   # uM → mM

        return extracted

    def _determine_verdict(self, results: list[BenchmarkCheckResult]) -> ValidationVerdict:
        """
        Apply verdict logic:
            NO DATA: no cited (scored) point at all, so nothing is validated
            FAIL: any TIGHT benchmark fails, or >50% of NORMAL fail
            WARN: all TIGHT pass, a minority of NORMAL fail
            MISS: all TIGHT/NORMAL pass, but a cited WIDE or STRICT point misses
            PASS: every cited point is inside its band
        """
        scored = [r for r in results if r.benchmark.is_scored]
        if not scored:
            return ValidationVerdict.NO_DATA
        # Only scored (cited) points count: an UNVERIFIED point is
        # informational and must never change the verdict, including by
        # diluting the NORMAL-failure ratio (B4 review, 2026-11-16).
        tight_results  = [r for r in scored if r.benchmark.tolerance == Tolerance.TIGHT]
        normal_results = [r for r in scored if r.benchmark.tolerance == Tolerance.NORMAL]

        tight_failures  = [r for r in tight_results  if not r.passed]
        normal_failures = [r for r in normal_results if not r.passed]

        if tight_failures:
            return ValidationVerdict.FAIL
        if normal_results and len(normal_failures) > len(normal_results) / 2:
            return ValidationVerdict.FAIL
        if normal_failures:
            return ValidationVerdict.WARN
        if any(r.benchmark.is_scored and not r.passed for r in results):
            return ValidationVerdict.MISS
        return ValidationVerdict.PASS
