# SAVE TO: packages\data-seeds\src\biosim_seeds\quickstart.py
"""
biosim quickstart demo.

Shows the full platform capability in ~50 lines.
Run from the project root:

    python -c "from biosim_seeds.quickstart import run_demo; run_demo()"

Or import and call individual demos:

    from biosim_seeds.quickstart import demo_lay_user, demo_researcher
"""

from __future__ import annotations


def demo_lay_user() -> None:
    """
    Lay user workflow: pick foods, simulate, print plain-English report.
    """
    print("\n=== LAY USER DEMO ===")
    from biosim_lay import compose_meal, generate_report
    from biosim_api import MealEvent, ScenarioDefinition, SimulationRunner

    # Pick foods from the database
    carbs_g = compose_meal(["white_rice_200g", "cola_330ml"])
    print(f"Meal: white rice + cola = {carbs_g:.0f}g carbohydrate")

    # Simulate
    scenario = ScenarioDefinition(
        name="rice_and_cola",
        duration_s=14400,
        meal_events=[MealEvent(time_s=0.0, carbs_g=carbs_g)],
    )
    output = SimulationRunner().run(scenario)
    report = generate_report(output)

    # Plain-English output
    print(f"\nGlucose: {report.glucose_section}")
    print(f"\nEnergy:  {report.energy_section}")
    print(f"\nRisk:    {report.risk_section.split(chr(10))[0]}")


def demo_researcher() -> None:
    """
    Researcher workflow: personalised simulation with UQ, metrics, CSV export.
    """
    print("\n=== RESEARCHER DEMO ===")
    from biosim_physiology import PhysiologicalProfile, BiologicalSex, DiseaseState
    from biosim_api import SimulationRunner, PersonalizationSpec
    from biosim_expert import CSVExporter
    from biosim_seeds import morning_meal
    import dataclasses

    # T2D patient, female, 58 years old, morning meal
    profile = PhysiologicalProfile(
        age_years=58,
        biological_sex=BiologicalSex.FEMALE,
        disease_states=[DiseaseState.TYPE_2_DIABETES],
        time_of_day_hours=8.0,
    )
    scenario = morning_meal(75.0)
    scenario = dataclasses.replace(
        scenario,
        personalization=PersonalizationSpec(physiological_profile=profile),
    )

    output = SimulationRunner().run(scenario)
    m = output.expert_metrics
    print(f"Peak glucose:  {m.glucose_peak_mM:.1f} mmol/L at {m.glucose_peak_time_s/60:.0f} min")
    print(f"Peak insulin:  {m.insulin_peak_pM:.0f} pM")
    print(f"mTOR peak:     {m.mtor_peak:.2f}")
    print(f"AMPK nadir:    {m.ampk_nadir:.2f}")
    print(f"Confidence:    {m.overall_watchdog_confidence:.2f}")

    csv_str = CSVExporter().export_metrics_table([output])
    print(f"\nCSV export: {len(csv_str.splitlines())} rows")


def demo_genetic() -> None:
    """
    Genetic personalisation: TCF7L2 risk variant effect on postprandial glucose.
    """
    print("\n=== GENETIC DEMO ===")
    from biosim_genotype import Genotype, Variant, GenotypePersonalizer
    from biosim_genotype.effect_db import extended_t2d_variants
    from biosim_ode import GlucoseInsulinPathway, ODESystem, high_sugar_meal

    db = extended_t2d_variants()
    personalizer = GenotypePersonalizer(effect_db=db)

    # TCF7L2 heterozygous risk
    variants = [Variant("rs7903146", Genotype.HETEROZYGOUS)]
    result = personalizer.personalize(variants)
    print(f"Secretion scale: {result.personalized_params.secretion_scale:.3f} (vs 1.0 baseline)")
    print(f"Parameters modified: {result.n_parameters_modified}")

    # Simulate
    pathway = GlucoseInsulinPathway(params=result.personalized_params)
    system  = ODESystem(pathway=pathway)
    sim     = system.simulate_meal(high_sugar_meal(), 7200)
    peak_G, peak_t = sim.glucose_peak()
    print(f"Peak glucose (risk): {peak_G/1000:.1f} mmol/L at {peak_t/60:.0f} min")

    # Baseline comparison
    baseline_sim = ODESystem().simulate_meal(high_sugar_meal(), 7200)
    baseline_G, _ = baseline_sim.glucose_peak()
    diff = (peak_G - baseline_G) / 1000.0
    print(f"Difference from baseline: {diff:+.2f} mmol/L")


def demo_population() -> None:
    """
    Population simulation: cohort of 10 healthy adults.
    """
    print("\n=== POPULATION DEMO ===")
    from biosim_seeds import PopulationSpec, PopulationRunner
    from biosim_api import get_scenario

    spec    = PopulationSpec(n_individuals=10, random_seed=42)
    runner  = PopulationRunner(spec)
    summary = runner.run(get_scenario("standard_meal"))

    print(f"Cohort: {summary.n_total} individuals, {summary.n_succeeded} succeeded")
    print(f"Glucose peak: {summary.glucose_peak_mean_mM:.1f} +/- {summary.glucose_peak_sd_mM:.1f} mmol/L")
    print(f"Range (5-95%): {summary.glucose_peak_p5_mM:.1f} - {summary.glucose_peak_p95_mM:.1f} mmol/L")
    print(f"Above 7.8 mmol/L: {summary.pct_above_7p8_mM:.0f}%")


def demo_validation() -> None:
    """
    Clinical validation against published OGTT reference values.
    """
    print("\n=== VALIDATION DEMO ===")
    from biosim_validation import ModelValidator

    report = ModelValidator().validate_ogtt()
    print(f"Verdict: {report.verdict.value.upper()}")
    print(f"Points within tolerance: {report.n_pass}/{report.n_pass + report.n_fail}")
    for r in report.results:
        print(f"  [{r.status}] {r.name}: model={r.model_value:.2f}, ref={r.benchmark.reference:.2f} ({r.deviation_percent:+.1f}%)")


def run_demo() -> None:
    """Run all demos in sequence."""
    print("BIOSIM PLATFORM QUICKSTART")
    print("=" * 40)
    demo_lay_user()
    demo_researcher()
    demo_genetic()
    demo_population()
    demo_validation()
    print("\n" + "=" * 40)
    print("Demo complete.")


if __name__ == "__main__":
    run_demo()
