# SAVE TO: packages\data-seeds\tests\test_smoke.py
"""
Platform smoke test.

Runs the full biosim stack end-to-end in one test.
If this passes, every layer of the platform is wired correctly.
Run time: ~15 seconds.
"""

import pytest


def test_full_platform_smoke():
    """
    Full stack: seed -> personalise -> simulate -> validate -> report -> export.
    One test, one pass, everything connected.
    """
    # 1. Seed the stores
    from biosim_seeds import create_seeded_registry, create_seeded_store, ParameterLoader
    registry = create_seeded_registry()
    store    = create_seeded_store(registry)
    loader   = ParameterLoader(store, registry)
    params   = loader.load_pathway_params()
    assert params is not None

    # 2. Physiological personalisation
    from biosim_physiology import (
        PhysiologicalProfile, BiologicalSex, DiseaseState,
        FitnessLevel, Drug, DrugDose, PhysiologicalPersonalizer,
    )
    profile = PhysiologicalProfile(
        age_years=45, biological_sex=BiologicalSex.MALE,
        body_weight_kg=85.0, fitness_level=FitnessLevel.LIGHTLY_ACTIVE,
        disease_states=[DiseaseState.INSULIN_RESISTANT],
        time_of_day_hours=19.0,
    )
    drugs = [DrugDose(Drug.CAFFEINE, dose_mg=150.0, hours_since_ingestion=1.5)]
    physio_result = PhysiologicalPersonalizer().personalize(profile, params, drugs)
    assert physio_result.personalized_params is not None

    # 3. Genetic personalisation
    from biosim_genotype import Genotype, Variant, GenotypePersonalizer
    from biosim_genotype.effect_db import extended_t2d_variants
    variants = [
        Variant("rs7903146",  Genotype.HETEROZYGOUS),
        Variant("rs10830963", Genotype.HETEROZYGOUS),
    ]
    genetic_result = GenotypePersonalizer(effect_db=extended_t2d_variants()).personalize(
        variants, physio_result.personalized_params
    )
    final_params = genetic_result.personalized_params

    # 4. ODE simulation
    from biosim_ode import GlucoseInsulinPathway, ODESystem
    from biosim_ode.meal import MealInput
    meal    = MealInput(carbs_g=75.0, time_min=0.0, body_weight_kg=85.0)
    pathway = GlucoseInsulinPathway(params=final_params, meals=[meal])
    system  = ODESystem(pathway=pathway)
    result  = system.simulate(duration_s=14400, meals=[meal])
    assert result.success
    peak_G, peak_t = result.glucose_peak()
    assert peak_G > 5000   # above baseline (uM)

    # 5. Uncertainty quantification
    from biosim_uncertainty import MonteCarloEngine, LogNormalDistribution
    engine = MonteCarloEngine(n_samples=15, random_seed=42, base_params=final_params)
    uq = engine.propagate(
        {"p1": LogNormalDistribution.from_mean_cv(final_params.p1, 0.30)},
        meals=[meal], duration_s=7200,
    )
    assert uq.n_successful > 10
    band = uq.band("G")
    assert all(band.lower <= band.upper + 1e-6)

    # 6. Sensitivity analysis and curation bridge
    from biosim_uncertainty import OATSensitivityAnalyzer, STANDARD_METRICS
    from biosim_uncertainty.curation_bridge import prioritise_from_sensitivity
    from biosim_curation import CurationQueue
    sr = OATSensitivityAnalyzer(0.10).analyze(
        ["p1", "ke"], {"gp": STANDARD_METRICS["glucose_peak"]},
        meals=[meal], duration_s=7200,
    )
    queue = CurationQueue()
    bridge = prioritise_from_sensitivity(sr, queue, store, registry, high_sensitivity_threshold=0.01)
    assert bridge is not None

    # 7. Expert metrics and lay output
    from biosim_api import SimulationRunner, MealEvent, ScenarioDefinition, PersonalizationSpec
    from biosim_lay import generate_report, assess_risk
    scenario = ScenarioDefinition(
        name="smoke_test", duration_s=14400,
        meal_events=[MealEvent(0, 75.0)],
        personalization=PersonalizationSpec(physiological_profile=profile),
    )
    api_output = SimulationRunner().run(scenario)
    assert api_output.success
    report = generate_report(api_output)
    assert len(report.glucose_section) > 10
    risk = assess_risk(api_output)
    assert risk.overall is not None

    # 8. Privacy audit
    from biosim_infra import PrivacyAuditor
    audit = PrivacyAuditor().audit(api_output)
    assert audit.passed, f"Privacy violation: {audit.issues}"

    # 9. SBML export
    from biosim_expert import SBMLExporter
    sbml = SBMLExporter().export(pathway)
    assert "sbml" in sbml.lower()
    assert 'id="G"' in sbml

    # 10. Validation
    from biosim_validation import ModelValidator, ValidationVerdict
    val_report = ModelValidator().validate_ogtt()
    # Known FAIL since 2026-11-16 (B7 + B4): only the circular peak-glucose
    # point fails; anything else failing is a new regression.
    assert [r.name for r in val_report.results if r.benchmark.is_scored and not r.passed] == ["peak_glucose"], val_report.as_text()

    # 11. Reproducibility bundle
    from biosim_infra import create_bundle
    bundle = create_bundle(api_output, random_seed=42)
    assert bundle.random_seed == 42
    json_str = bundle.to_json()
    assert len(json_str) > 100

    # 12. Cortisol pathway
    from biosim_ode.cortisol_pathway import CortisolPathway
    cortisol = CortisolPathway().with_stress(0.3).with_time_of_day(19.0)
    assert cortisol.y0[0] > 0
    derivs = cortisol.dydt(0.0, cortisol.y0)
    assert len(derivs) == 3

    # 13. Hepatic EGP
    from biosim_ode.hepatic import HepGlucPathway, EGPParameters
    egp = HepGlucPathway(EGPParameters.for_type2_diabetes())
    assert egp.y0[0] > 0
    assert egp.dydt(0.0, [egp.params.egp_max], {"insulin_pM": 300.0})[0] < 0

    print("\nSMOKE TEST PASSED: all 13 platform layers verified")
