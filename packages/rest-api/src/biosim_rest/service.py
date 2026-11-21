# SAVE TO: packages\rest-api\src\biosim_rest\service.py
"""
Translation layer between REST schemas and biosim_api's internal types.

Kept out of app.py so the routes stay thin -- routes parse/validate the
HTTP boundary and call these; all ScenarioDefinition/SimulationOutput
knowledge lives here.
"""

from __future__ import annotations

from .schemas import (
    BatchResponse,
    BatchResponseItem,
    BenchmarkResultOut,
    CohortSummaryOut,
    DrugEventIn,
    ExpertMetricsOut,
    LaySummaryOut,
    MealEventIn,
    MetabolicRiskIndicatorOut,
    MetabolicRiskProfileOut,
    PopulationRequest,
    RiskIndicatorOut,
    RiskProfileOut,
    ScenarioIn,
    SimulateResponse,
    ValidationDatasetOut,
    ValidationSummaryResponse,
)


def _physiological_profile_in_to_profile(p):
    """
    Translates PhysiologicalProfileIn -> biosim_physiology.PhysiologicalProfile.
    Raises ValueError for an invalid enum string (biological_sex/
    fitness_level/disease_states) -- already caught by /simulate's
    existing `except ValueError -> 400` handler, matching DrugEvent's
    own established convention (see PhysiologicalProfileIn's docstring).
    """
    from biosim_physiology import (
        BiologicalSex, DiseaseState, FitnessLevel, MenstrualCyclePhase, PhysiologicalProfile,
    )

    return PhysiologicalProfile(
        age_years=p.age_years,
        biological_sex=BiologicalSex(p.biological_sex),
        body_weight_kg=p.body_weight_kg,
        height_cm=p.height_cm,
        body_fat_percent=p.body_fat_percent,
        fitness_level=FitnessLevel(p.fitness_level),
        disease_states=[DiseaseState(d) for d in p.disease_states],
        time_of_day_hours=p.time_of_day_hours,
        sleep_hours=p.sleep_hours,
        desaturation_events_per_hour=p.desaturation_events_per_hour,
        menstrual_cycle_phase=(
            MenstrualCyclePhase(p.menstrual_cycle_phase)
            if p.menstrual_cycle_phase is not None else None
        ),
        cold_acclimatized=p.cold_acclimatized,
    )


def scenario_in_to_definition(payload: ScenarioIn):
    from biosim_api import DrugEvent, InsulinDoseEvent, MealEvent, PersonalizationSpec, ScenarioDefinition

    personalization = PersonalizationSpec()
    if payload.physiological_profile is not None:
        personalization = PersonalizationSpec(
            physiological_profile=_physiological_profile_in_to_profile(payload.physiological_profile)
        )

    return ScenarioDefinition(
        name=payload.name,
        duration_s=payload.duration_s,
        description=payload.description,
        meal_events=[
            MealEvent(
                time_s=m.time_s, carbs_g=m.carbs_g, label=m.label,
                glycaemic_index=m.glycaemic_index, kabs=m.kabs,
                complex_carb_fraction=m.complex_carb_fraction,
                protein_g=m.protein_g,
            )
            for m in payload.meal_events
        ],
        drug_events=[
            DrugEvent(
                time_s=d.time_s, drug=d.drug, dose_mg=d.dose_mg,
                hours_before_simulation=d.hours_before_simulation,
            )
            for d in payload.drug_events
        ],
        insulin_dose_events=[
            InsulinDoseEvent(
                time_s=i.time_s, insulin_type=i.insulin_type, dose_units=i.dose_units,
                hours_before_simulation=i.hours_before_simulation,
            )
            for i in payload.insulin_dose_events
        ],
        personalization=personalization,
        run_uncertainty=payload.run_uncertainty,
        n_mc_samples=payload.n_mc_samples,
        stress_level=payload.stress_level,
        cushings_syndrome=payload.cushings_syndrome,
        model_nefa_dynamics=payload.model_nefa_dynamics,
    )


def simulation_output_to_response(output) -> SimulateResponse:
    from biosim_lay.report import assess_risk
    from biosim_lay.risk_score import compute_risk_profile

    m = output.expert_metrics
    ls = output.lay_summary
    sim = output.simulation
    risk = assess_risk(output)
    metabolic_risk = compute_risk_profile(output)

    return SimulateResponse(
        scenario_name=output.scenario_name,
        scenario_description=output.scenario_description,
        success=output.success,
        expert_metrics=ExpertMetricsOut(
            glucose_peak_mM=m.glucose_peak_mM,
            glucose_peak_time_s=m.glucose_peak_time_s,
            glucose_auc_uM_s=m.glucose_auc_uM_s,
            glucose_baseline_uM=m.glucose_baseline_uM,
            insulin_peak_pM=m.insulin_peak_pM,
            insulin_peak_time_s=m.insulin_peak_time_s,
            mtor_peak=m.mtor_peak,
            ampk_nadir=m.ampk_nadir,
            time_to_glucose_return_s=m.time_to_glucose_return_s,
            overall_watchdog_confidence=m.overall_watchdog_confidence,
            glucose_nadir_mM=m.glucose_nadir_mM,
        ),
        lay_summary=LaySummaryOut(
            times_min=ls.times_min.tolist(),
            energy_level=ls.energy_level.tolist(),
            hunger_level=ls.hunger_level.tolist(),
            glucose_narrative=ls.glucose_narrative,
            energy_narrative=ls.energy_narrative,
            confidence_note=ls.confidence_note,
            is_high_confidence=ls.is_high_confidence,
        ),
        risk_profile=RiskProfileOut(
            indicators=[
                RiskIndicatorOut(
                    label=i.label,
                    value=i.value,
                    unit=i.unit,
                    traffic_light=i.traffic_light.value,
                    description=i.description,
                    caveat=i.caveat,
                )
                for i in risk.indicators
            ],
            overall=risk.overall.value,
            summary=risk.summary,
        ),
        metabolic_risk_profile=MetabolicRiskProfileOut(
            indicators=[
                MetabolicRiskIndicatorOut(
                    name=i.name,
                    value=i.value,
                    unit=i.unit,
                    band=i.band.value,
                    interpretation=i.interpretation,
                    caveat=i.caveat,
                )
                for i in metabolic_risk.indicators
            ],
            overall_band=metabolic_risk.overall_band.value,
            headline=metabolic_risk.headline,
            disclaimer=metabolic_risk.disclaimer,
            n_elevated_or_high=metabolic_risk.n_elevated_or_high,
        ),
        times_s=sim.times_s.tolist() if output.success else [],
        glucose_uM=sim.species["G"].tolist() if output.success else [],
        insulin_uM=sim.species["I"].tolist() if output.success else [],
        warnings=output.warnings,
    )


def population_spec_from_request(payload: PopulationRequest):
    """MASTER_PLAN 3.11 -- always the PopulationRunner default
    (use_multiprocessing=False, see PopulationRequest's own docstring
    for why that flag is never exposed over REST)."""
    from biosim_seeds import PopulationSpec

    return PopulationSpec(n_individuals=payload.n_individuals, random_seed=payload.random_seed)


def cohort_summary_to_response(summary) -> CohortSummaryOut:
    """summary: biosim_seeds.population.CohortSummary."""
    return CohortSummaryOut(
        n_total=summary.n_total,
        n_succeeded=summary.n_succeeded,
        n_failed=summary.n_failed,
        glucose_peak_mean_mM=summary.glucose_peak_mean_mM,
        glucose_peak_sd_mM=summary.glucose_peak_sd_mM,
        glucose_peak_p5_mM=summary.glucose_peak_p5_mM,
        glucose_peak_p95_mM=summary.glucose_peak_p95_mM,
        insulin_peak_mean_pM=summary.insulin_peak_mean_pM,
        mtor_peak_mean=summary.mtor_peak_mean,
        ampk_nadir_mean=summary.ampk_nadir_mean,
        pct_above_7p8_mM=summary.pct_above_7p8_mM,
        pct_above_11p1_mM=summary.pct_above_11p1_mM,
        failure_reasons=summary.failure_reasons,
        warnings=summary.warnings,
    )


def batch_result_to_response(scenarios: list, result) -> BatchResponse:
    """
    scenarios: list[ScenarioDefinition] in the original request order.
    result: biosim_expert.batch.BatchResult.

    BatchResult itself only preserves order for `outputs` (re-sorted by
    scenario_name) and `failed_names`/`failure_reasons` (parallel,
    completion order) -- neither matches the client's original request
    order. Rebuilt here by looking each original scenario up by name,
    so a client always gets results back in the order it submitted them.
    """
    failure_by_name = dict(zip(result.failed_names, result.failure_reasons))
    items = []
    for s in scenarios:
        output = result.get(s.name)
        if output is not None:
            items.append(BatchResponseItem(
                name=s.name, succeeded=True,
                result=simulation_output_to_response(output),
            ))
        else:
            items.append(BatchResponseItem(
                name=s.name, succeeded=False,
                failure_reason=failure_by_name.get(s.name, "unknown failure"),
            ))
    return BatchResponse(
        n_succeeded=result.n_succeeded,
        n_failed=result.n_failed,
        results=items,
    )


def validation_reports_to_response(reports: dict) -> ValidationSummaryResponse:
    """reports: dict[str, ValidationReport] -- see ModelValidator.validate_all()."""
    return ValidationSummaryResponse(
        datasets=[
            ValidationDatasetOut(
                dataset_name=report.dataset_name,
                verdict=report.verdict.value,
                n_pass=report.n_pass,
                n_fail=report.n_fail,
                results=[
                    BenchmarkResultOut(
                        name=r.name,
                        tolerance=r.benchmark.tolerance.value,
                        is_scored=r.benchmark.is_scored,
                        reference=r.benchmark.reference,
                        unit=r.benchmark.unit,
                        model_value=r.model_value,
                        deviation_percent=r.deviation_percent,
                        passed=r.passed,
                    )
                    for r in report.results
                ],
                notes=report.notes,
                module_label=report.module_label,
            )
            for report in reports.values()
        ],
    )
