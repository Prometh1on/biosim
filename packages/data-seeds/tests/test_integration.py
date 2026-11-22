# SAVE TO: packages\data-seeds\tests\test_integration.py
"""
Comprehensive end-to-end integration tests.

These tests exercise the full platform stack from scenario definition
through to lay report, validating that all components interoperate
correctly. They are intentionally higher-level than unit tests —
each test exercises a realistic user workflow.
"""

import pytest
import numpy as np


# ---------------------------------------------------------------------------
# Workflow 1: Lay user — pick food, simulate, get report
# ---------------------------------------------------------------------------

class TestLayUserWorkflow:
    """
    A lay user selects food items, runs a simulation, and reads their report.
    No biology knowledge required.
    """

    def test_full_lay_workflow(self):
        from biosim_lay import compose_meal, get_food, generate_report, assess_risk
        from biosim_api import MealEvent, ScenarioDefinition, SimulationRunner

        # User picks foods
        carbs = compose_meal(["white_rice_200g", "cola_330ml"])
        assert carbs == pytest.approx(91.0)

        # Simulate
        scenario = ScenarioDefinition(
            name="rice_and_cola",
            duration_s=14400,
            meal_events=[MealEvent(time_s=0.0, carbs_g=carbs)],
        )
        output = SimulationRunner().run(scenario)
        assert output.success

        # Report
        report = generate_report(output)
        assert "rice_and_cola" in report.headline
        assert len(report.glucose_section) > 10
        assert report.disclaimer is not None

        # Risk
        risk = assess_risk(output)
        assert risk.overall is not None

    def test_low_vs_high_carb_risk_difference(self):
        from biosim_lay import compose_meal, assess_risk
        from biosim_api import MealEvent, ScenarioDefinition, SimulationRunner
        from biosim_lay.report import TrafficLight

        runner = SimulationRunner()

        low_carb = ScenarioDefinition(
            name="low", duration_s=7200,
            meal_events=[MealEvent(0, compose_meal(["egg_boiled", "salad_plain"]))],
        )
        high_carb = ScenarioDefinition(
            name="high", duration_s=7200,
            meal_events=[MealEvent(0, compose_meal(["white_rice_200g", "cola_330ml"]))],
        )

        out_low  = runner.run(low_carb)
        out_high = runner.run(high_carb)

        assert out_high.expert_metrics.glucose_peak_mM > out_low.expert_metrics.glucose_peak_mM


# ---------------------------------------------------------------------------
# Workflow 2: Researcher — personalised simulation with SBML export
# ---------------------------------------------------------------------------

class TestResearcherWorkflow:
    """
    A researcher runs a personalised simulation and exports to SBML.
    """

    def test_personalised_sbml_export(self):
        from biosim_physiology import PhysiologicalProfile, BiologicalSex, FitnessLevel
        from biosim_api import SimulationRunner, PersonalizationSpec
        from biosim_seeds import morning_meal
        from biosim_expert import SBMLExporter
        from biosim_ode import GlucoseInsulinPathway

        profile = PhysiologicalProfile(
            age_years=45,
            biological_sex=BiologicalSex.FEMALE,
            fitness_level=FitnessLevel.VERY_ACTIVE,
        )
        scenario = morning_meal(75.0)
        import dataclasses
        scenario = dataclasses.replace(
            scenario,
            personalization=PersonalizationSpec(physiological_profile=profile),
        )
        output = SimulationRunner().run(scenario)
        assert output.success
        assert output.physiological_personalization is not None

        sbml = SBMLExporter().export(GlucoseInsulinPathway())
        assert "sbml" in sbml.lower()
        assert "G" in sbml  # glucose species

    def test_batch_parameter_comparison(self):
        from biosim_expert import BatchRunner, CSVExporter
        from biosim_api import MealEvent, ScenarioDefinition
        import dataclasses

        base = ScenarioDefinition(name="base", meal_events=[MealEvent(0, 75.0)])
        scenarios = [
            dataclasses.replace(base, name=f"carbs_{c}g",
                meal_events=[MealEvent(0, float(c))])
            for c in [25, 50, 75]
        ]
        result = BatchRunner(max_workers=2).run_all(scenarios)
        assert result.n_succeeded == 3

        csv_str = CSVExporter().export_metrics_table(result.outputs)
        assert "glucose_peak_mM" in csv_str
        assert len(csv_str.strip().split("\n")) == 4  # header + 3 rows

    def test_validation_report(self):
        from biosim_validation import ModelValidator, ValidationVerdict
        report = ModelValidator().validate_ogtt()
        # Known FAIL since 2026-11-16 (B7 + B4): only the circular peak-glucose
        # point fails; anything else failing is a new regression.
        assert [r.name for r in report.results if r.benchmark.is_scored and not r.passed] == ["peak_glucose"], report.as_text()


# ---------------------------------------------------------------------------
# Workflow 3: Platform session — single entry point
# ---------------------------------------------------------------------------

class TestPlatformSessionWorkflow:
    @pytest.fixture(scope="class")
    @classmethod
    def session(cls):
        from biosim_seeds import PlatformSession
        return PlatformSession()

    def test_session_simulate_and_report(self, session):
        output = session.simulate("high_sugar_meal")
        assert output.success
        report = session.report(output)
        assert report is not None

    def test_session_compare(self, session):
        comparison = session.compare("high_sugar_meal", "light_snack")
        assert comparison.peak_glucose_diff_uM > 0

    def test_session_audit_clean(self, session):
        output = session.simulate("standard_meal")
        audit = session.audit_privacy(output)
        assert audit.passed

    def test_session_validate(self, session):
        from biosim_validation import ValidationVerdict
        report = session.validate()
        # Known FAIL since 2026-11-16 (B7 + B4): only the circular peak-glucose
        # point fails; anything else failing is a new regression.
        assert [r.name for r in report.results if r.benchmark.is_scored and not r.passed] == ["peak_glucose"], report.as_text()

    def test_session_sbml_export(self, session):
        sbml = session.export_sbml()
        assert isinstance(sbml, str)
        assert len(sbml) > 500

    def test_session_bundle_reproducible(self, session):
        output = session.simulate("fasted")
        bundle = session.bundle(output, random_seed=42)
        assert bundle.random_seed == 42
        json_str = bundle.to_json()
        import json
        d = json.loads(json_str)
        assert d["random_seed"] == 42


# ---------------------------------------------------------------------------
# Workflow 4: Genetic personalisation
# ---------------------------------------------------------------------------

class TestGeneticWorkflow:
    def test_high_risk_genotype_raises_glucose(self):
        from biosim_genotype import Genotype, Variant, GenotypePersonalizer
        from biosim_genotype.effect_db import extended_t2d_variants
        from biosim_ode import GlucoseInsulinPathway, ODESystem, high_sugar_meal

        db = extended_t2d_variants()
        personalizer = GenotypePersonalizer(effect_db=db)

        # Person homozygous for TCF7L2, CDKN2A, CDKAL1 risk alleles
        high_risk_variants = [
            Variant("rs7903146",  Genotype.HOMOZYGOUS_ALT),
            Variant("rs10811661", Genotype.HOMOZYGOUS_ALT),
            Variant("rs7756992",  Genotype.HOMOZYGOUS_ALT),
        ]
        result = personalizer.personalize(high_risk_variants)
        assert result.personalized_params.secretion_scale < 0.75

        # Run simulation
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        system  = ODESystem(pathway=pathway)
        sim     = system.simulate_meal(high_sugar_meal(), 7200)
        assert sim.success

        # Should produce higher peak than baseline
        baseline = ODESystem(GlucoseInsulinPathway())
        peak_baseline, _ = baseline.simulate_meal(high_sugar_meal(), 7200).glucose_peak()
        peak_risk,     _ = sim.glucose_peak()
        assert peak_risk > peak_baseline

    def test_protective_genotype_lowers_glucose(self):
        from biosim_genotype import Genotype, Variant, GenotypePersonalizer
        from biosim_genotype.effect_db import extended_t2d_variants
        from biosim_ode import GlucoseInsulinPathway, ODESystem, high_sugar_meal

        db = extended_t2d_variants()
        personalizer = GenotypePersonalizer(effect_db=db)
        variants = [Variant("rs10010131", Genotype.HOMOZYGOUS_ALT)]  # WFS1 protective
        result = personalizer.personalize(variants)
        assert result.personalized_params.secretion_scale > 1.0

    def test_privacy_preserved(self):
        from biosim_genotype import Genotype, Variant
        from biosim_api import SimulationRunner, PersonalizationSpec, high_sugar_meal_scenario
        from biosim_infra import PrivacyAuditor
        import re

        variants = [Variant("rs7903146", Genotype.HETEROZYGOUS)]
        scenario = high_sugar_meal_scenario(
            personalization=PersonalizationSpec(genetic_variants=variants)
        )
        output = SimulationRunner().run(scenario)
        audit  = PrivacyAuditor().audit(output)
        assert audit.passed
        # Double-check: no rsIDs in summary
        summary_str = str(output.summary())
        assert not re.search(r'\brs\d{4,}\b', summary_str)


# ---------------------------------------------------------------------------
# Workflow 5: Chronobiology and exercise
# ---------------------------------------------------------------------------

class TestChronoExerciseWorkflow:
    def test_morning_vs_evening_glucose(self):
        from biosim_seeds import morning_meal, evening_meal
        from biosim_api import SimulationRunner
        runner = SimulationRunner()
        out_m = runner.run(morning_meal(75.0))
        out_e = runner.run(evening_meal(75.0))
        # Evening meal should produce higher glucose peak
        assert out_e.expert_metrics.glucose_peak_mM > out_m.expert_metrics.glucose_peak_mM

    def test_exercise_reduces_postprandial_glucose(self):
        from biosim_physiology import moderate_run
        from biosim_ode import GlucoseInsulinPathway, ODESystem, high_sugar_meal, PathwayParameters
        import dataclasses

        base_params = PathwayParameters()
        e = moderate_run(45.0, hours_before=1.0)
        mods = e.post_exercise_modifiers()
        ex_params = dataclasses.replace(base_params, **{
            k: getattr(base_params, k) * v
            for k, v in mods.items() if hasattr(base_params, k)
        })

        peak_base, _ = ODESystem(GlucoseInsulinPathway(params=base_params))\
            .simulate_meal(high_sugar_meal(), 7200).glucose_peak()
        peak_ex, _   = ODESystem(GlucoseInsulinPathway(params=ex_params))\
            .simulate_meal(high_sugar_meal(), 7200).glucose_peak()

        assert peak_ex < peak_base


# ---------------------------------------------------------------------------
# Workflow 6: Uncertainty and curation loop
# ---------------------------------------------------------------------------

class TestUncertaintyCurationLoop:
    def test_sensitivity_identifies_top_params(self):
        from biosim_uncertainty import OATSensitivityAnalyzer, STANDARD_METRICS
        from biosim_ode import high_sugar_meal

        analyzer = OATSensitivityAnalyzer(perturbation_fraction=0.10)
        result = analyzer.analyze(
            parameter_names=["p1", "ke", "p3", "EC50_mTOR"],
            metrics={"glucose_peak": STANDARD_METRICS["glucose_peak"]},
            meals=[high_sugar_meal()],
            duration_s=7200,
        )
        top = result.top_parameters("glucose_peak", n=2)
        assert len(top) == 2
        # Top sensitivity parameters should have nonzero scores
        assert all(abs(score) > 0 for _, score in top)

    def test_uq_produces_confidence_bands(self):
        from biosim_uncertainty import MonteCarloEngine, LogNormalDistribution
        from biosim_ode import high_sugar_meal

        engine = MonteCarloEngine(n_samples=20, random_seed=42)
        result = engine.propagate(
            distributions={"p1": LogNormalDistribution.from_mean_cv(0.0317, 0.30)},
            meals=[high_sugar_meal()],
            duration_s=7200,
        )
        G_band = result.band("G")
        # Lower must be <= median <= upper everywhere
        assert np.all(G_band.lower <= G_band.median + 1e-6)
        assert np.all(G_band.median <= G_band.upper + 1e-6)
        assert result.n_successful > 15

    def test_curation_bridge_connects_to_queue(self):
        from biosim_uncertainty import OATSensitivityAnalyzer, STANDARD_METRICS
        from biosim_uncertainty.curation_bridge import prioritise_from_sensitivity
        from biosim_curation import CurationQueue
        from biosim_seeds import create_seeded_registry, create_seeded_store
        from biosim_ode import high_sugar_meal

        registry = create_seeded_registry()
        store    = create_seeded_store(registry)
        queue    = CurationQueue()

        sr = OATSensitivityAnalyzer(0.10).analyze(
            parameter_names=["p1", "ke"],
            metrics={"glucose_peak": STANDARD_METRICS["glucose_peak"]},
            meals=[high_sugar_meal()], duration_s=7200,
        )
        result = prioritise_from_sensitivity(
            sr, queue, store, registry,
            high_sensitivity_threshold=0.01,
        )
        assert isinstance(result.n_items_enqueued, int)


# ---------------------------------------------------------------------------
# Workflow 7: Cortisol + hepatic coupling
# ---------------------------------------------------------------------------

class TestCortisolHepatic:
    def test_cortisol_pathway_standalone(self):
        from biosim_ode.cortisol_pathway import CortisolPathway
        from scipy.integrate import solve_ivp
        p = CortisolPathway()
        sol = solve_ivp(p.dydt, [0, 120], p.y0, method="Radau")
        assert sol.success
        assert sol.y[0, -1] > 0   # cortisol still positive

    def test_hepatic_pathway_suppressed_by_insulin(self):
        from biosim_ode.hepatic import HepGlucPathway
        p = HepGlucPathway()
        # At high insulin (300 pM), EGP target should be well below max
        high_insulin_signals = {"insulin_pM": 300.0}
        derivs_high_i = p.dydt(0.0, [p.params.egp_max], signals=high_insulin_signals)
        # Derivative should be negative (EGP falling toward suppressed target)
        assert derivs_high_i[0] < 0

    def test_hepatic_t2d_params_higher_egp(self):
        from biosim_ode.hepatic import HepGlucPathway, EGPParameters
        p_healthy = HepGlucPathway()
        p_t2d     = HepGlucPathway(EGPParameters.for_type2_diabetes())
        assert p_t2d.params.egp_basal > p_healthy.params.egp_basal

    def test_hepatic_metformin_reduces_egp(self):
        # metformin_scale was retired -- for_metformin() now represents
        # its effect via a lowered egp_basal directly (same pattern as
        # T2D's elevated one), so the fixed point holds for a chronic
        # metformin preset instead of drifting away from it. See
        # test_hepatic.py's test_metformin_y0_is_genuine_fixed_point.
        from biosim_ode.hepatic import HepGlucPathway, EGPParameters
        p_control   = HepGlucPathway()
        p_metformin = HepGlucPathway(EGPParameters.for_metformin(1000.0))
        assert p_metformin.params.egp_basal < p_control.params.egp_basal

    def test_composite_glucose_insulin_cortisol(self):
        from biosim_ode.cortisol_pathway import CortisolPathway
        from biosim_ode.coupling import CompositePathway
        from biosim_ode import GlucoseInsulinPathway, ODESystem
        from scipy.integrate import solve_ivp

        gi = GlucoseInsulinPathway()
        co = CortisolPathway()
        composite = CompositePathway([gi, co])

        sol = solve_ivp(composite.dydt, [0, 60], composite.y0, method="Radau")
        assert sol.success
        assert sol.y.shape[0] == 9
