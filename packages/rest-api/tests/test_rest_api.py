# SAVE TO: packages\rest-api\tests\test_rest_api.py
"""
Tests for the Biosim REST API.

Uses FastAPI's TestClient (in-process, no real network/server) to exercise
every endpoint end-to-end through the actual SimulationRunner pipeline --
these are integration tests, not mocked route tests.
"""
import pytest
from fastapi.testclient import TestClient

from biosim_rest import create_app


@pytest.fixture(scope="module")
def client():
    return TestClient(create_app())


class TestHealth:
    def test_returns_200(self, client):
        resp = client.get("/health")
        assert resp.status_code == 200

    def test_has_status_and_version_only(self, client):
        # /health used to report a hard-coded test count that went stale
        # (2764 when the suite had 3000+). A health check reports health, not
        # a snapshot of the test suite, so the field was removed (2026-11-21).
        data = client.get("/health").json()
        assert data["status"] == "ok"
        assert isinstance(data["version"], str)
        assert set(data) == {"status", "version"}


class TestValidationEndpoint:
    """
    GET /validation exposes ModelValidator.validate_all() over HTTP --
    every registered clinical benchmark (OGTT, Bergman IVGTT, semaglutide,
    steroid/T2D) built up over this session had no way to be checked
    without direct Python access to biosim_validation.
    """

    @pytest.fixture(scope="class")
    def validation_data(self, client):
        resp = client.get("/validation")
        assert resp.status_code == 200
        return resp.json()

    def test_returns_200(self, client):
        resp = client.get("/validation")
        assert resp.status_code == 200

    def test_includes_all_nine_datasets(self, validation_data):
        names = {ds["dataset_name"] for ds in validation_data["datasets"]}
        assert names == {
            "OGTT_75g_healthy",
            "Bergman_IV_glucose",
            "Semaglutide_1mg_postprandial",
            "Prednisone_40mg_T2D_postprandial",
            "Metformin_T2D_postprandial",
            "Canagliflozin_300mg_postprandial",
            "Pea_protein_coingestion_50g_glucose",
            "T2D_vs_nondiabetic_meal_insulin",
            "Caffeine_5mgkg_OGTT_insulin",
        }

    def test_no_dataset_hard_fails(self, validation_data):
        # The in-house OGTT is a known FAIL since 2026-11-16 (see
        # test_validation.py::TestModelValidator::test_validate_ogtt_convenience).
        for ds in validation_data["datasets"]:
            if ds["dataset_name"] == "OGTT_75g_healthy":
                assert ds["verdict"] == "FAIL"
                continue
            assert ds["verdict"] != "FAIL", f"{ds['dataset_name']} hard-failed"

    def test_unscored_point_is_flagged_over_the_wire(self, validation_data):
        # mtor_activation must carry is_scored=False through to the HTTP
        # response, not just internally -- an API consumer reading raw
        # JSON needs the same "this isn't a real check" signal a Python
        # caller gets from BenchmarkPoint.is_scored.
        ogtt = next(ds for ds in validation_data["datasets"] if ds["dataset_name"] == "OGTT_75g_healthy")
        mtor = next(r for r in ogtt["results"] if r["name"] == "mtor_activation")
        assert mtor["is_scored"] is False

    def test_each_result_has_deviation_and_reference(self, validation_data):
        for ds in validation_data["datasets"]:
            for r in ds["results"]:
                assert isinstance(r["deviation_percent"], float)
                assert isinstance(r["reference"], float)
                assert r["unit"]


class TestListScenarios:
    def test_returns_200(self, client):
        resp = client.get("/scenarios")
        assert resp.status_code == 200

    def test_includes_known_scenarios(self, client):
        names = {s["name"] for s in client.get("/scenarios").json()}
        assert "high_sugar_meal" in names
        assert "fasted" in names
        assert "ogtt_75g" in names

    def test_items_have_description_and_duration(self, client):
        data = client.get("/scenarios").json()
        for item in data:
            assert isinstance(item["description"], str)
            assert item["duration_s"] > 0


class TestSimulate:
    def test_high_sugar_meal_succeeds(self, client):
        resp = client.post("/simulate", json={
            "name": "test_meal",
            "duration_s": 14400,
            "meal_events": [{"time_s": 0, "carbs_g": 75.0, "label": "test"}],
        })
        assert resp.status_code == 200
        data = resp.json()
        assert data["success"] is True
        assert data["expert_metrics"]["glucose_peak_mM"] > 5.0
        assert len(data["times_s"]) > 0
        assert len(data["glucose_uM"]) == len(data["times_s"])

    def test_fasted_stays_near_baseline(self, client):
        resp = client.post("/simulate", json={
            "name": "test_fasted",
            "duration_s": 3600,
            "meal_events": [],
        })
        data = resp.json()
        assert data["success"] is True
        assert max(data["glucose_uM"]) - min(data["glucose_uM"]) < 100.0

    def test_empty_name_returns_422(self, client):
        resp = client.post("/simulate", json={"name": "", "meal_events": []})
        assert resp.status_code == 422

    def test_negative_duration_returns_422(self, client):
        resp = client.post("/simulate", json={"name": "x", "duration_s": -1})
        assert resp.status_code == 422

    def test_stress_level_out_of_range_returns_422(self, client):
        resp = client.post("/simulate", json={"name": "x", "stress_level": 1.5})
        assert resp.status_code == 422

    def test_n_mc_samples_above_ceiling_returns_422(self, client):
        # Housekeeping-backlog fix: n_mc_samples previously had no
        # upper bound at all on this endpoint -- an unbounded,
        # uninterruptible compute-cost lever reachable by any client
        # with run_uncertainty=True (see ScenarioDefinition
        # .__post_init__'s own comment for the full reasoning and the
        # 1000 ceiling's derivation).
        resp = client.post("/simulate", json={
            "name": "x", "meal_events": [], "n_mc_samples": 1001,
        })
        assert resp.status_code == 422

    def test_n_mc_samples_at_ceiling_is_accepted(self, client):
        # Boundary check: exactly 1000 must still be accepted (only
        # values ABOVE the ceiling are rejected). Deliberately does NOT
        # set run_uncertainty=True -- actually running 1000 real Monte
        # Carlo samples would make this test slow for no added
        # discriminating value; this confirms the SCHEMA accepts the
        # boundary value, which is what this fix actually changed.
        resp = client.post("/simulate", json={
            "name": "x", "meal_events": [], "n_mc_samples": 1000,
        })
        assert resp.status_code == 200

    def test_excessive_duration_returns_422(self, client):
        # Found by the mandatory security-reviewer pass for the
        # n_mc_samples fix above -- CRITICAL, fixed in the same pass.
        # An unbounded duration_s forces the ODE solver to integrate
        # across an arbitrarily large time span, independent of
        # n_mc_samples/run_uncertainty entirely.
        resp = client.post("/simulate", json={
            "name": "x", "meal_events": [], "duration_s": 604801.0,
        })
        assert resp.status_code == 422

    def test_excessive_meal_events_returns_422(self, client):
        # Same pass: unbounded event-list length, combined with
        # MealAbsorptionPathway's per-meal ODE state allocation and the
        # Radau solver's dense finite-difference Jacobian, made a large
        # meal_events list an O(N^2)-memory DoS vector reachable with
        # no run_uncertainty involved at all.
        resp = client.post("/simulate", json={
            "name": "x",
            "meal_events": [{"time_s": i * 60, "carbs_g": 10.0} for i in range(51)],
        })
        assert resp.status_code == 422

    def test_meal_events_at_ceiling_is_accepted(self, client):
        resp = client.post("/simulate", json={
            "name": "x", "duration_s": 3600,
            "meal_events": [{"time_s": i * 60, "carbs_g": 1.0} for i in range(50)],
        })
        assert resp.status_code == 200

    def test_missing_name_returns_422(self, client):
        resp = client.post("/simulate", json={"meal_events": []})
        assert resp.status_code == 422

    def test_glp1_drug_event_reduces_glucose_peak(self, client):
        baseline = client.post("/simulate", json={
            "name": "baseline", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
        }).json()
        treated = client.post("/simulate", json={
            "name": "semaglutide",
            "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
            "drug_events": [{"time_s": 0, "drug": "semaglutide", "dose_mg": 1.0}],
        }).json()
        assert treated["expert_metrics"]["glucose_peak_mM"] < baseline["expert_metrics"]["glucose_peak_mM"]

    def test_drug_event_hours_before_simulation_actually_applies(self, client):
        # Regression: DrugEventIn had no hours_before_simulation field at
        # all, and scenario_in_to_definition() never set it -- every
        # time-dependent drug (caffeine, alcohol) submitted through the
        # REST API was silently modeled at hours_since_ingestion=0
        # regardless of the client's intended timing. For caffeine this
        # made the drug effect exactly zero (caffeine_modifiers()'s
        # concentration_fraction has a t<=0 floor), with no way to express
        # "N hours before the meal" via this API at all.
        no_drug = client.post("/simulate", json={
            "name": "nd", "meal_events": [{"time_s": 0, "carbs_g": 60.0}],
        }).json()
        caffeine_no_lead_time = client.post("/simulate", json={
            "name": "caf_now",
            "meal_events": [{"time_s": 0, "carbs_g": 60.0}],
            "drug_events": [{"time_s": 0, "drug": "caffeine", "dose_mg": 200.0}],
        }).json()
        caffeine_45min_before = client.post("/simulate", json={
            "name": "caf_lead",
            "meal_events": [{"time_s": 0, "carbs_g": 60.0}],
            "drug_events": [{
                "time_s": 0, "drug": "caffeine", "dose_mg": 200.0,
                "hours_before_simulation": 0.75,
            }],
        }).json()
        peak_nd = no_drug["expert_metrics"]["glucose_peak_mM"]
        peak_now = caffeine_no_lead_time["expert_metrics"]["glucose_peak_mM"]
        peak_lead = caffeine_45min_before["expert_metrics"]["glucose_peak_mM"]
        assert peak_now == pytest.approx(peak_nd), (
            "caffeine with zero lead time (hours_since_ingestion=0) should "
            "have no effect -- confirms the baseline comparison is meaningful"
        )
        assert peak_lead > peak_nd, (
            "caffeine 45 minutes before the meal must actually raise the "
            "glucose peak -- hours_before_simulation must reach the ODE"
        )

    def test_meal_event_glycaemic_index_actually_applies(self, client):
        # Regression: MealEventIn had no glycaemic_index/kabs fields at
        # all, and scenario_in_to_definition() never set them -- a REST
        # client could never submit a GI-adjusted meal, silently always
        # getting the default monoexponential absorption regardless of
        # the food's real glycaemic index (same shape as the
        # DrugEventIn.hours_before_simulation gap fixed above, for meals
        # instead of drugs).
        high_gi = client.post("/simulate", json={
            "name": "high_gi",
            "meal_events": [{"time_s": 0, "carbs_g": 75.0, "glycaemic_index": 100.0}],
        }).json()
        low_gi = client.post("/simulate", json={
            "name": "low_gi",
            "meal_events": [{"time_s": 0, "carbs_g": 75.0, "glycaemic_index": 30.0}],
        }).json()
        assert (
            high_gi["expert_metrics"]["glucose_peak_time_s"]
            != low_gi["expert_metrics"]["glucose_peak_time_s"]
        ), "glycaemic_index must reach the ODE's absorption model, not be silently dropped"

    def test_meal_event_complex_carb_fraction_actually_applies(self, client):
        # Same shape as the glycaemic_index gap fixed above (MASTER_PLAN
        # 1.4b(iii)): confirms the new field reaches SimulationRunner's
        # acarbose scaling through the real REST route, not just the
        # internal ScenarioDefinition/MealEvent path already tested in
        # test_simulation_api.py.
        free_sugar = client.post("/simulate", json={
            "name": "free_sugar",
            "meal_events": [{"time_s": 0, "carbs_g": 75.0, "complex_carb_fraction": 0.0}],
            "drug_events": [{"time_s": 0, "drug": "acarbose", "dose_mg": 100.0}],
        }).json()
        starch = client.post("/simulate", json={
            "name": "starch",
            "meal_events": [{"time_s": 0, "carbs_g": 75.0, "complex_carb_fraction": 1.0}],
            "drug_events": [{"time_s": 0, "drug": "acarbose", "dose_mg": 100.0}],
        }).json()
        assert (
            free_sugar["expert_metrics"]["glucose_peak_mM"]
            > starch["expert_metrics"]["glucose_peak_mM"]
        ), "complex_carb_fraction must reach the acarbose scaling, not be silently dropped"

    def test_meal_event_protein_g_actually_applies(self, client):
        # Same shape as complex_carb_fraction/glycaemic_index above --
        # MASTER_PLAN row 1.28's dietary protein co-ingestion mechanism.
        no_protein = client.post("/simulate", json={
            "name": "no_protein",
            "meal_events": [{"time_s": 0, "carbs_g": 50.0, "protein_g": 0.0}],
        }).json()
        with_protein = client.post("/simulate", json={
            "name": "with_protein",
            "meal_events": [{"time_s": 0, "carbs_g": 50.0, "protein_g": 50.0}],
        }).json()
        assert (
            with_protein["expert_metrics"]["insulin_peak_pM"]
            > no_protein["expert_metrics"]["insulin_peak_pM"]
        ), "protein_g must reach NutrientSecretionPathway, not be silently dropped"

    def test_meal_event_protein_g_out_of_range_returns_422(self, client):
        resp = client.post("/simulate", json={
            "name": "bad_protein",
            "meal_events": [{"time_s": 0, "carbs_g": 50.0, "protein_g": -1.0}],
        })
        assert resp.status_code == 422

    def test_meal_event_protein_g_above_ceiling_returns_422(self, client):
        # Mandatory security-reviewer finding, HIGH: the REST ceiling
        # (150.0) matches MAX_PLAUSIBLE_PROTEIN_G_PER_MEAL, verified
        # against the real composed ODE's insulin output before choosing
        # it -- confirm the boundary is actually enforced over REST.
        resp = client.post("/simulate", json={
            "name": "too_much_protein",
            "meal_events": [{"time_s": 0, "carbs_g": 50.0, "protein_g": 151.0}],
        })
        assert resp.status_code == 422

    def test_meal_event_protein_g_at_ceiling_is_accepted(self, client):
        resp = client.post("/simulate", json={
            "name": "max_protein",
            "meal_events": [{"time_s": 0, "carbs_g": 50.0, "protein_g": 150.0}],
        })
        assert resp.status_code == 200

    def test_meal_event_rejects_both_glycaemic_index_and_kabs(self, client):
        resp = client.post("/simulate", json={
            "name": "bad",
            "meal_events": [{
                "time_s": 0, "carbs_g": 75.0, "glycaemic_index": 80.0, "kabs": 0.05,
            }],
        })
        assert resp.status_code == 400

    def test_extreme_hours_before_simulation_returns_400_not_500(self, client):
        # Regression: an extreme (but schema-valid -- Pydantic only
        # enforces ge=0, no upper bound) hours_before_simulation for an
        # insulin_regular event overflows dose_time_min to -inf inside
        # ExogenousInsulinParameters, which raises ValueError -- that
        # exception used to propagate uncaught through SimulationRunner
        # .run(), past app.py's only except-ValueError->400 handler
        # (which only wrapped scenario_in_to_definition()), producing an
        # unhandled 500 instead of a graceful 400.
        resp = client.post("/simulate", json={
            "name": "extreme_insulin_timing",
            "drug_events": [{
                "time_s": 0, "drug": "insulin_regular", "dose_mg": 10.0,
                "hours_before_simulation": 1e307,
            }],
        })
        assert resp.status_code == 400

    def test_stress_level_raises_glucose_peak(self, client):
        no_stress = client.post("/simulate", json={
            "name": "ns", "meal_events": [{"time_s": 0, "carbs_g": 75.0}], "stress_level": 0.0,
        }).json()
        stressed = client.post("/simulate", json={
            "name": "s", "meal_events": [{"time_s": 0, "carbs_g": 75.0}], "stress_level": 0.9,
        }).json()
        assert stressed["expert_metrics"]["glucose_peak_mM"] > no_stress["expert_metrics"]["glucose_peak_mM"]

    def test_glucose_nadir_mM_present_and_plausible(self, client):
        # Regression (MASTER_PLAN 3.11b): ExpertMetrics.glucose_nadir_mM
        # existed internally (added to feed the hypoglycaemia-risk lay
        # indicator) but had no wire-schema field at all.
        resp = client.post("/simulate", json={
            "name": "nadir_check",
            "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
        }).json()
        m = resp["expert_metrics"]
        assert "glucose_nadir_mM" in m
        # The nadir can never exceed the peak, and for a real meal
        # scenario should sit at or below the fasting baseline.
        assert m["glucose_nadir_mM"] <= m["glucose_peak_mM"]
        assert m["glucose_nadir_mM"] > 0.0

    def test_risk_profile_present_with_real_indicators(self, client):
        # Regression (MASTER_PLAN 3.11b): assess_risk()'s output
        # (traffic-light indicators) was never exposed via /simulate at
        # all -- a REST client could get raw numbers or a narrative but
        # never this platform's own risk assessment.
        resp = client.post("/simulate", json={
            "name": "risk_check",
            "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
        }).json()
        assert "risk_profile" in resp
        rp = resp["risk_profile"]
        assert len(rp["indicators"]) > 0
        assert rp["overall"] in ("green", "amber", "red", "grey")
        first = rp["indicators"][0]
        assert set(first) == {"label", "value", "unit", "traffic_light", "description", "caveat"}

    def test_metabolic_risk_profile_present_with_real_indicators(self, client):
        resp = client.post("/simulate", json={
            "name": "metabolic_risk_check",
            "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
        }).json()
        assert "metabolic_risk_profile" in resp
        mrp = resp["metabolic_risk_profile"]
        assert len(mrp["indicators"]) > 0
        assert mrp["overall_band"] in ("low", "moderate", "elevated", "high")
        assert "research tool" in mrp["disclaimer"]
        first = mrp["indicators"][0]
        assert set(first) == {"name", "value", "unit", "band", "interpretation", "caveat"}

    def test_metabolic_risk_profile_is_unknown_not_healthy_when_simulation_fails(self):
        # Regression, found by the mandatory code-reviewer pass on this
        # milestone: compute_risk_profile() never checked
        # output.simulation.success at all -- for a failed simulation it
        # reported overall_band="low" and a "healthy" headline, now
        # reachable over this public REST endpoint. Fixed in
        # biosim_lay.risk_score (RiskBand.UNKNOWN) -- this test confirms
        # the fix survives the full REST translation layer, not just the
        # unit-level compute_risk_profile() call.
        from unittest.mock import patch
        from fastapi.testclient import TestClient
        from biosim_rest import create_app
        from biosim_api import SimulationRunner, ExpertMetrics, ScenarioDefinition
        import numpy as np
        from types import SimpleNamespace

        real_output = SimulationRunner().run(ScenarioDefinition(
            name="pre_failure_probe", duration_s=3600,
        ))
        failed_sim = SimpleNamespace(
            success=False, species={}, species_units={},
            times_s=np.array([0.0]), message="forced failure for this test",
        )
        failed_metrics = ExpertMetrics(
            glucose_peak_uM=0, glucose_peak_time_s=0, glucose_auc_uM_s=0,
            glucose_baseline_uM=5000.0, glucose_nadir_uM=5000.0,
            insulin_peak_uM=0, insulin_peak_time_s=0,
            mtor_peak=0, ampk_nadir=1.0,
            time_to_glucose_return_s=None,
            overall_watchdog_confidence=0.0,
            n_ode_evaluations=0,
        )
        failed_output = SimpleNamespace(
            scenario_name="forced_failure",
            scenario_description="",
            simulation=failed_sim,
            expert_metrics=failed_metrics,
            lay_summary=real_output.lay_summary,
            physiological_personalization=None,
            genetic_personalization=None,
            warnings=[],
            success=False,
        )

        isolated_client = TestClient(create_app())
        with patch.object(SimulationRunner, "run", return_value=failed_output):
            resp = isolated_client.post("/simulate", json={
                "name": "forced_failure", "meal_events": [], "duration_s": 3600,
            })
        assert resp.status_code == 200
        data = resp.json()
        assert data["success"] is False
        assert data["metabolic_risk_profile"]["overall_band"] == "unknown"
        assert "healthy" not in data["metabolic_risk_profile"]["headline"].lower()

    def test_t2d_scenario_produces_a_genuinely_elevated_risk_profile(self, client):
        # Not just "the keys exist" -- confirm the exposed risk profile
        # actually differs meaningfully for a T2D-personalized scenario
        # vs. a healthy one, through the real REST route end-to-end.
        healthy = client.post("/simulate", json={
            "name": "healthy_risk",
            "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
        }).json()
        t2d = client.post("/simulate", json={
            "name": "t2d_risk",
            "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
            "drug_events": [{"time_s": 0, "drug": "dapagliflozin", "dose_mg": 10.0}],
        }).json()
        # dapagliflozin lowers the glucose peak -- confirm the risk
        # profile's own peak-glucose value reflects the same real
        # SimulationOutput /simulate itself returns, not a stale or
        # separately-computed number.
        healthy_peak_indicator = next(
            i for i in healthy["risk_profile"]["indicators"] if i["label"] == "Peak blood sugar"
        )
        t2d_peak_indicator = next(
            i for i in t2d["risk_profile"]["indicators"] if i["label"] == "Peak blood sugar"
        )
        assert healthy_peak_indicator["value"] == pytest.approx(
            healthy["expert_metrics"]["glucose_peak_mM"], abs=0.05
        )
        assert t2d_peak_indicator["value"] == pytest.approx(
            t2d["expert_metrics"]["glucose_peak_mM"], abs=0.05
        )


class TestPhysiologicalProfileEndpoint:
    """
    MASTER_PLAN row 3.14 / process-improvement item C2: ScenarioIn had
    no field for disease state, age/sex/BMI, or fitness at all -- a REST
    client could submit drugs and meals but had no way to say WHO the
    scenario represents. Confirmed via direct grep this was previously
    entirely absent from schemas.py.
    """

    def test_t2d_profile_raises_fasting_glucose_through_the_real_rest_route(self, client):
        # The actual disease-state preset now reachable over REST, not
        # a drug-based proxy (see test_t2d_scenario_produces_a_genuinely
        # _elevated_risk_profile above, which had to use dapagliflozin as
        # a workaround before this field existed).
        healthy = client.post("/simulate", json={
            "name": "healthy_profile", "meal_events": [],
        }).json()
        t2d = client.post("/simulate", json={
            "name": "t2d_profile", "meal_events": [],
            "physiological_profile": {"disease_states": ["type_2_diabetes"]},
        }).json()
        assert t2d["expert_metrics"]["glucose_baseline_uM"] > \
            healthy["expert_metrics"]["glucose_baseline_uM"]

    def test_omitting_the_field_entirely_is_unaffected(self, client):
        # Backward-compatibility guard: a request with no
        # physiological_profile key at all must behave exactly as before
        # this field existed (personalization stays at its default,
        # empty PersonalizationSpec).
        resp = client.post("/simulate", json={
            "name": "no_profile", "meal_events": [],
        })
        assert resp.status_code == 200
        assert resp.json()["expert_metrics"]["glucose_baseline_uM"] == pytest.approx(5000.0)

    def test_full_profile_all_fields_reach_the_simulation(self, client):
        # Age/sex/BMI/fitness all compose too, not just disease_states --
        # confirm a profile combining several fields still simulates
        # successfully and produces a different result from the default.
        default = client.post("/simulate", json={
            "name": "default_profile", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
        }).json()
        custom = client.post("/simulate", json={
            "name": "custom_profile", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
            "physiological_profile": {
                "age_years": 65.0,
                "biological_sex": "female",
                "body_weight_kg": 90.0,
                "height_cm": 165.0,
                "fitness_level": "sedentary",
                "disease_states": ["insulin_resistant"],
                "time_of_day_hours": 20.0,
            },
        }).json()
        assert custom["success"] is True
        assert custom["expert_metrics"]["glucose_peak_mM"] != pytest.approx(
            default["expert_metrics"]["glucose_peak_mM"], abs=0.001
        )

    def test_sleep_hours_reaches_the_simulation_through_the_real_rest_route(self, client):
        # Fresh Track 1 scoping row 1.26 -- confirms sleep_hours composes
        # end-to-end over REST, same convention as the disease_states
        # test above.
        full_sleep = client.post("/simulate", json={
            "name": "full_sleep", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
            "physiological_profile": {"sleep_hours": 10.0},
        }).json()
        restricted = client.post("/simulate", json={
            "name": "restricted_sleep", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
            "physiological_profile": {"sleep_hours": 5.0},
        }).json()
        assert restricted["expert_metrics"]["glucose_peak_mM"] > \
            full_sleep["expert_metrics"]["glucose_peak_mM"]

    def test_out_of_range_sleep_hours_returns_422(self, client):
        resp = client.post("/simulate", json={
            "name": "bad_sleep", "meal_events": [],
            "physiological_profile": {"sleep_hours": 30.0},
        })
        assert resp.status_code == 422

    def test_generic_bound_clamp_warning_reaches_rest_output_for_t1d(self, client):
        # D1 Composition-Boundary Behavioral Audit (2026-08-18): T1D's
        # own Ib/secretion_scale factors clamp to their _PARAM_BOUNDS
        # floors (since row 1.10), previously with zero visibility
        # anywhere. Confirms the new generic clamp warning actually
        # reaches the real REST response's own `warnings` list (the
        # same field clamp_warnings has always been merged into, per
        # runner.py -- there is no separate `clamp_warnings` field on
        # SimulateResponse).
        resp = client.post("/simulate", json={
            "name": "t1d_clamp", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
            "physiological_profile": {"disease_states": ["type_1_diabetes"]},
        }).json()
        joined = " ".join(resp["warnings"])
        assert "secretion_scale" in joined
        assert "Ib" in joined

    def test_indinavir_reaches_the_simulation_through_the_real_rest_route(self, client):
        # Track 1 row 1.33.
        no_drug = client.post("/simulate", json={
            "name": "no_indinavir", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
        }).json()
        with_drug = client.post("/simulate", json={
            "name": "with_indinavir", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
            "drug_events": [{"time_s": 0, "drug": "indinavir", "dose_mg": 1200.0}],
        }).json()
        assert with_drug["expert_metrics"]["glucose_peak_mM"] > \
            no_drug["expert_metrics"]["glucose_peak_mM"]

    def test_hydrochlorothiazide_reaches_the_simulation_through_the_real_rest_route(self, client):
        # Track 1 row 1.34.
        no_drug = client.post("/simulate", json={
            "name": "no_hctz", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
        }).json()
        with_drug = client.post("/simulate", json={
            "name": "with_hctz", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
            "drug_events": [{"time_s": 0, "drug": "hydrochlorothiazide", "dose_mg": 40.0}],
        }).json()
        assert with_drug["expert_metrics"]["glucose_peak_mM"] > \
            no_drug["expert_metrics"]["glucose_peak_mM"]

    def test_growth_hormone_reaches_the_simulation_through_the_real_rest_route(self, client):
        # Track 1 row 1.35.
        no_drug = client.post("/simulate", json={
            "name": "no_gh", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
        }).json()
        with_drug = client.post("/simulate", json={
            "name": "with_gh", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
            "drug_events": [{"time_s": 0, "drug": "growth_hormone", "dose_mg": 1.0}],
        }).json()
        assert with_drug["expert_metrics"]["glucose_peak_mM"] > \
            no_drug["expert_metrics"]["glucose_peak_mM"]

    def test_two_sulfonylureas_with_a_profile_are_rejected_with_400_not_a_silent_baseline(self, client):
        # NAFLD/guard code-review MEDIUM: with a physiological_profile the
        # runner's try/except around personalisation swallowed the
        # personalizer's ValueError and returned a 200 healthy-baseline
        # simulation. The guard must fire at scenario construction.
        resp = client.post("/simulate", json={
            "name": "two_sus_profile", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
            "physiological_profile": {"disease_states": ["type_2_diabetes"]},
            "drug_events": [
                {"time_s": 0, "drug": "glibenclamide", "dose_mg": 2.5},
                {"time_s": 0, "drug": "glimepiride", "dose_mg": 1.0},
            ],
        })
        assert resp.status_code == 400
        assert "sulfonylurea" in resp.json()["detail"]

    def test_glp1_drug_with_gi_meal_is_accepted(self, client):
        # Was rejected with 400 (closed-form Ra x gastric_emptying_modifier
        # deleted ~38% of the meal's glucose). Since 2026-11-15 two-compartment
        # meals are state-tracked and mass-conserving, so this runs.
        resp = client.post("/simulate", json={
            "name": "gi_dpp4",
            "meal_events": [{"time_s": 0, "carbs_g": 75.0, "glycaemic_index": 100.0}],
            "drug_events": [{"time_s": 0, "drug": "dpp4_inhibitor", "dose_mg": 100.0}],
        })
        assert resp.status_code == 200

    def test_duplicate_hcq_events_are_rejected_with_400(self, client):
        # HCQ security-review HIGH: 20 duplicate flat-factor events used
        # to reach the p3 ceiling with only a warning.
        resp = client.post("/simulate", json={
            "name": "dup_hcq", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
            "drug_events": [
                {"time_s": 0, "drug": "hydroxychloroquine", "dose_mg": 200.0},
                {"time_s": 0, "drug": "hydroxychloroquine", "dose_mg": 200.0},
            ],
        })
        assert resp.status_code == 400
        assert "hydroxychloroquine" in resp.json()["detail"]

    def test_duplicate_dose_scaled_events_are_rejected_with_400(self, client):
        # Two 1000mg metformin events applied the full 1000mg factor twice
        # (compounding past the 2000mg ceiling) instead of one 2000mg dose.
        resp = client.post("/simulate", json={
            "name": "dup_metformin", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
            "drug_events": [
                {"time_s": 0, "drug": "metformin", "dose_mg": 1000.0},
                {"time_s": 0, "drug": "metformin", "dose_mg": 1000.0},
            ],
        })
        assert resp.status_code == 400
        assert "total dose" in resp.json()["detail"]

    def test_degenerate_height_is_rejected_not_silently_simulated_as_healthy(self, client):
        # NAFLD security-review M4: height_cm=1e-300 used to return 200 with
        # the whole profile dropped (T2D simulated as healthy) and only a
        # warning string.
        resp = client.post("/simulate", json={
            "name": "tiny_height", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
            "physiological_profile": {
                "disease_states": ["type_2_diabetes"], "height_cm": 1e-300,
            },
        })
        # Exactly 400 (a ValueError from PhysiologicalProfile construction),
        # NOT a schema 422 -- so this exercises the profile validation, and
        # the body must name the cause rather than return any simulation.
        assert resp.status_code == 400
        assert "height_cm" in resp.json()["detail"]
        assert "expert_metrics" not in resp.json()

    def test_hydroxychloroquine_reaches_the_simulation_through_the_real_rest_route(self, client):
        # Track 1 row 1.44 -- p3-RAISING, so the meal peak goes DOWN
        # (opposite direction to the p3-lowering drugs above).
        no_drug = client.post("/simulate", json={
            "name": "no_hcq", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
        }).json()
        with_drug = client.post("/simulate", json={
            "name": "with_hcq", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
            "drug_events": [{"time_s": 0, "drug": "hydroxychloroquine", "dose_mg": 400.0}],
        }).json()
        assert with_drug["expert_metrics"]["glucose_peak_mM"] <             no_drug["expert_metrics"]["glucose_peak_mM"]

    def test_two_distinct_sulfonylureas_are_rejected_with_400_not_500(self, client):
        # Glimepiride chunk's security-review finding, through the real
        # route: the personalizer's fail-fast ValueError must surface as a
        # clean 400, not an unhandled 500.
        resp = client.post("/simulate", json={
            "name": "two_sus", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
            "drug_events": [
                {"time_s": 0, "drug": "glibenclamide", "dose_mg": 2.5},
                {"time_s": 0, "drug": "glimepiride", "dose_mg": 1.0},
            ],
        })
        assert resp.status_code == 400
        assert "sulfonylurea" in resp.json()["detail"]

    def test_atorvastatin_reaches_the_simulation_through_the_real_rest_route(self, client):
        # Track 1 row 1.41.
        no_drug = client.post("/simulate", json={
            "name": "no_atorvastatin", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
        }).json()
        with_drug = client.post("/simulate", json={
            "name": "with_atorvastatin", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
            "drug_events": [{"time_s": 0, "drug": "atorvastatin", "dose_mg": 80.0}],
        }).json()
        assert with_drug["expert_metrics"]["glucose_peak_mM"] > \
            no_drug["expert_metrics"]["glucose_peak_mM"]

    def test_intermittent_hypoxia_reaches_the_simulation_through_the_real_rest_route(self, client):
        # Track 1 row 1.32 -- confirms desaturation_events_per_hour
        # composes end-to-end over REST, same convention as sleep_hours
        # above.
        no_hypoxia = client.post("/simulate", json={
            "name": "no_hypoxia", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
            "physiological_profile": {"desaturation_events_per_hour": 0.0},
        }).json()
        hypoxia = client.post("/simulate", json={
            "name": "with_hypoxia", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
            "physiological_profile": {"desaturation_events_per_hour": 24.3},
        }).json()
        assert hypoxia["expert_metrics"]["glucose_peak_mM"] > \
            no_hypoxia["expert_metrics"]["glucose_peak_mM"]

    def test_out_of_range_desaturation_events_per_hour_returns_422(self, client):
        resp = client.post("/simulate", json={
            "name": "bad_hypoxia", "meal_events": [],
            "physiological_profile": {"desaturation_events_per_hour": -1.0},
        })
        assert resp.status_code == 422

    def test_menstrual_cycle_phase_reaches_the_simulation_through_the_real_rest_route(self, client):
        # Track 1 row 1.36.
        follicular = client.post("/simulate", json={
            "name": "follicular", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
            "physiological_profile": {"menstrual_cycle_phase": "follicular"},
        }).json()
        luteal = client.post("/simulate", json={
            "name": "luteal", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
            "physiological_profile": {"menstrual_cycle_phase": "luteal"},
        }).json()
        assert luteal["expert_metrics"]["glucose_peak_mM"] > \
            follicular["expert_metrics"]["glucose_peak_mM"]

    def test_invalid_menstrual_cycle_phase_returns_400_not_500(self, client):
        resp = client.post("/simulate", json={
            "name": "bad_phase", "meal_events": [],
            "physiological_profile": {"menstrual_cycle_phase": "not_a_real_phase"},
        })
        assert resp.status_code == 400

    def test_pcos_reaches_the_simulation_through_the_real_rest_route(self, client):
        # Track 1 row 1.29.
        healthy = client.post("/simulate", json={
            "name": "healthy", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
        }).json()
        pcos = client.post("/simulate", json={
            "name": "pcos", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
            "physiological_profile": {"disease_states": ["pcos"]},
        }).json()
        assert pcos["expert_metrics"]["glucose_peak_mM"] > \
            healthy["expert_metrics"]["glucose_peak_mM"]

    def test_nafld_reaches_the_simulation_through_the_real_rest_route(self, client):
        # Track 1 row 1.45.
        healthy = client.post("/simulate", json={
            "name": "healthy", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
        }).json()
        nafld = client.post("/simulate", json={
            "name": "nafld", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
            "physiological_profile": {"disease_states": ["nafld"]},
        }).json()
        assert (
            nafld["expert_metrics"]["glucose_peak_mM"]
            > healthy["expert_metrics"]["glucose_peak_mM"]
        )

    def test_acromegaly_reaches_the_simulation_through_the_real_rest_route(self, client):
        # Track 1 row 1.37.
        healthy = client.post("/simulate", json={
            "name": "healthy", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
        }).json()
        acromegaly = client.post("/simulate", json={
            "name": "acromegaly", "meal_events": [{"time_s": 0, "carbs_g": 75.0}],
            "physiological_profile": {"disease_states": ["acromegaly"]},
        }).json()
        assert acromegaly["expert_metrics"]["glucose_peak_mM"] > \
            healthy["expert_metrics"]["glucose_peak_mM"]

    def test_invalid_disease_state_returns_400_not_500(self, client):
        resp = client.post("/simulate", json={
            "name": "bad_disease", "meal_events": [],
            "physiological_profile": {"disease_states": ["not_a_real_disease"]},
        })
        assert resp.status_code == 400

    def test_invalid_biological_sex_returns_400_not_500(self, client):
        resp = client.post("/simulate", json={
            "name": "bad_sex", "meal_events": [],
            "physiological_profile": {"biological_sex": "not_a_real_value"},
        })
        assert resp.status_code == 400

    def test_invalid_fitness_level_returns_400_not_500(self, client):
        resp = client.post("/simulate", json={
            "name": "bad_fitness", "meal_events": [],
            "physiological_profile": {"fitness_level": "not_a_real_value"},
        })
        assert resp.status_code == 400

    def test_conflicting_disease_states_returns_400_not_500(self, client):
        # PhysiologicalProfile.__post_init__'s own mutual-exclusivity
        # validation (T1D + T2D) must still fire through the REST path,
        # not just when constructed directly in Python.
        resp = client.post("/simulate", json={
            "name": "conflicting_disease", "meal_events": [],
            "physiological_profile": {
                "disease_states": ["type_1_diabetes", "type_2_diabetes"]
            },
        })
        assert resp.status_code == 400

    def test_out_of_range_age_returns_422(self, client):
        resp = client.post("/simulate", json={
            "name": "bad_age", "meal_events": [],
            "physiological_profile": {"age_years": 200.0},
        })
        assert resp.status_code == 422

    def test_excessive_disease_states_list_returns_422(self, client):
        # Only 8 real DiseaseState values exist -- bounded defensively
        # (Q2 checklist finding), matching meal_events/drug_events/
        # insulin_dose_events's own established max_length pattern.
        resp = client.post("/simulate", json={
            "name": "bad_disease_list", "meal_events": [],
            "physiological_profile": {"disease_states": ["healthy"] * 11},
        })
        assert resp.status_code == 422

    def test_excessive_body_weight_returns_422(self, client):
        # Security-reviewer finding (C2's mandatory pass): body_weight_kg
        # previously had no upper bound.
        resp = client.post("/simulate", json={
            "name": "bad_weight", "meal_events": [],
            "physiological_profile": {"body_weight_kg": 501.0},
        })
        assert resp.status_code == 422

    def test_excessive_height_returns_422(self, client):
        resp = client.post("/simulate", json={
            "name": "bad_height", "meal_events": [],
            "physiological_profile": {"height_cm": 301.0},
        })
        assert resp.status_code == 422


class TestNanInfinityValidationErrorSanitization:
    """
    MASTER_PLAN 3.13: a raw JSON `NaN`/`Infinity` literal for a float
    field is correctly rejected by Pydantic's own validation, but
    FastAPI's default RequestValidationError handler used to crash
    trying to re-serialize the rejected value back into the 422 error
    body (Starlette's JSONResponse.render() calls
    json.dumps(..., allow_nan=False)) -- turning a client mistake
    Pydantic already caught into an unhandled 500. Generic to the whole
    API, not specific to any one field -- confirmed against two
    unrelated float fields.
    """

    def test_nan_dose_mg_returns_422_not_500(self, client):
        raw = (
            b'{"name": "nan_dose_mg", "meal_events": [], "duration_s": 3600, '
            b'"drug_events": [{"time_s": 0, "drug": "metformin", "dose_mg": NaN}]}'
        )
        resp = client.post("/simulate", content=raw, headers={"Content-Type": "application/json"})
        assert resp.status_code == 422
        detail = resp.json()["detail"]
        assert detail[0]["input"] == "NaN"

    def test_nan_carbs_g_returns_422_not_500(self, client):
        raw = (
            b'{"name": "nan_carbs", "meal_events": [{"time_s": 0, "carbs_g": NaN}], '
            b'"duration_s": 3600}'
        )
        resp = client.post("/simulate", content=raw, headers={"Content-Type": "application/json"})
        assert resp.status_code == 422
        detail = resp.json()["detail"]
        assert detail[0]["input"] == "NaN"

    def test_infinity_dose_mg_still_returns_400_not_500(self, client):
        raw = (
            b'{"name": "inf_dose_mg", "meal_events": [], "duration_s": 3600, '
            b'"drug_events": [{"time_s": 0, "drug": "metformin", "dose_mg": Infinity}]}'
        )
        resp = client.post("/simulate", content=raw, headers={"Content-Type": "application/json"})
        # dose_mg=Infinity satisfies Pydantic's Field(ge=0) (inf >= 0 is
        # True), so this is rejected downstream, in DrugEvent.__post_init__
        # (math.isfinite()) -- NOT a RequestValidationError at all. Included
        # here to confirm this path was never crashing either, not just
        # the NaN one this fix specifically targets.
        assert resp.status_code == 400

    def test_ordinary_validation_error_shape_is_unaffected(self, client):
        # A ordinary (finite) out-of-range value must still produce the
        # exact same 422 shape as before this fix -- the sanitizer must
        # be a no-op for every normal value.
        resp = client.post("/simulate", json={
            "name": "ordinary_negative", "meal_events": [], "duration_s": 3600,
            "drug_events": [{"time_s": 0, "drug": "metformin", "dose_mg": -5.0}],
        })
        assert resp.status_code == 422
        detail = resp.json()["detail"]
        assert detail[0]["input"] == -5.0


class TestInsulinDoseEventEndpoint:
    """
    MASTER_PLAN 1.2b (event-type) exposed over REST: InsulinDoseEventIn
    mirrors biosim_api.scenario.InsulinDoseEvent field-for-field, the
    honestly-Unit-named alternative to DrugEventIn's dose_mg-as-Units
    hack for insulin_regular/insulin_lispro.
    """

    def test_insulin_dose_event_reaches_the_simulation(self, client):
        no_insulin = client.post("/simulate", json={
            "name": "no_insulin_via_rest", "meal_events": [], "duration_s": 36000,
        }).json()
        with_insulin = client.post("/simulate", json={
            "name": "insulin_via_rest", "meal_events": [], "duration_s": 36000,
            "insulin_dose_events": [
                {"time_s": 0, "insulin_type": "regular", "dose_units": 10.0}
            ],
        }).json()
        assert (
            with_insulin["expert_metrics"]["insulin_peak_pM"]
            > no_insulin["expert_metrics"]["insulin_peak_pM"]
        )

    def test_matches_the_legacy_drug_event_form_exactly(self, client):
        via_drug_event = client.post("/simulate", json={
            "name": "legacy_form", "meal_events": [], "duration_s": 36000,
            "drug_events": [{"time_s": 0, "drug": "insulin_regular", "dose_mg": 10.0}],
        }).json()
        via_insulin_event = client.post("/simulate", json={
            "name": "new_form", "meal_events": [], "duration_s": 36000,
            "insulin_dose_events": [
                {"time_s": 0, "insulin_type": "regular", "dose_units": 10.0}
            ],
        }).json()
        assert via_insulin_event["expert_metrics"]["insulin_peak_pM"] == pytest.approx(
            via_drug_event["expert_metrics"]["insulin_peak_pM"]
        )

    def test_invalid_insulin_type_returns_422(self, client):
        resp = client.post("/simulate", json={
            "name": "bad_insulin_type", "meal_events": [], "duration_s": 3600,
            "insulin_dose_events": [
                {"time_s": 0, "insulin_type": "glargine", "dose_units": 10.0}
            ],
        })
        assert resp.status_code == 422

    def test_combining_with_legacy_drug_event_returns_400(self, client):
        # ScenarioDefinition.__post_init__'s mutual-exclusivity check must
        # surface as a client error (400), not an unhandled 500.
        resp = client.post("/simulate", json={
            "name": "conflicting_forms", "meal_events": [], "duration_s": 3600,
            "drug_events": [{"time_s": 0, "drug": "insulin_regular", "dose_mg": 10.0}],
            "insulin_dose_events": [
                {"time_s": 0, "insulin_type": "regular", "dose_units": 10.0}
            ],
        })
        assert resp.status_code == 400


class TestSimulateCache:
    """
    /simulate is backed by biosim_infra.SimulationCache -- the same
    caching component biosim_seeds.PlatformSession already used, wired
    into this HTTP layer for the first time (previously every request,
    identical or not, re-ran the full ODE solve from scratch). Uses its
    own isolated TestClient/app (not the shared module-scoped `client`
    fixture) so its cache starts empty and isn't polluted by other tests.
    """

    @pytest.fixture
    def isolated_client(self):
        return TestClient(create_app())

    def test_identical_requests_hit_the_cache(self, isolated_client):
        from unittest.mock import patch
        from biosim_api import SimulationRunner

        payload = {
            "name": "cache_test_meal",
            "duration_s": 14400,
            "meal_events": [{"time_s": 0, "carbs_g": 60.0, "label": "cache test"}],
        }
        real_run = SimulationRunner.run
        with patch.object(SimulationRunner, "run", autospec=True, side_effect=real_run) as mock_run:
            first = isolated_client.post("/simulate", json=payload)
            second = isolated_client.post("/simulate", json=payload)

        assert first.status_code == 200
        assert second.status_code == 200
        assert mock_run.call_count == 1, (
            "identical /simulate requests should hit the cache on the second "
            "call, not re-run the ODE solver"
        )
        assert first.json() == second.json()

    def test_different_requests_do_not_share_a_cache_entry(self, isolated_client):
        from unittest.mock import patch
        from biosim_api import SimulationRunner

        real_run = SimulationRunner.run
        with patch.object(SimulationRunner, "run", autospec=True, side_effect=real_run) as mock_run:
            isolated_client.post("/simulate", json={
                "name": "a", "duration_s": 14400,
                "meal_events": [{"time_s": 0, "carbs_g": 50.0}],
            })
            isolated_client.post("/simulate", json={
                "name": "b", "duration_s": 14400,
                "meal_events": [{"time_s": 0, "carbs_g": 90.0}],
            })

        assert mock_run.call_count == 2, (
            "different scenarios must not collide in the cache or reuse "
            "each other's result"
        )


class TestCompare:
    def test_returns_200(self, client):
        resp = client.post("/compare", json={"scenario_a": "high_sugar_meal", "scenario_b": "light_snack"})
        assert resp.status_code == 200

    def test_high_sugar_higher_than_snack(self, client):
        data = client.post("/compare", json={"scenario_a": "high_sugar_meal", "scenario_b": "light_snack"}).json()
        assert data["peak_glucose_diff_mM"] > 0
        assert "higher" in data["narrative"]

    def test_unknown_scenario_returns_404(self, client):
        resp = client.post("/compare", json={"scenario_a": "not_a_real_scenario", "scenario_b": "fasted"})
        assert resp.status_code == 404

    def test_response_includes_both_full_results(self, client):
        data = client.post("/compare", json={"scenario_a": "fasted", "scenario_b": "high_sugar_meal"}).json()
        assert data["scenario_a"]["success"] is True
        assert data["scenario_b"]["success"] is True


class TestPopulation:
    """
    POST /population -- MASTER_PLAN 3.11's bounded-synchronous-size-cap
    fix. n_individuals is deliberately small in every test here (3) to
    keep the suite fast -- the cohort still runs the full ODE pipeline
    per individual, real work, not mocked.
    """

    def test_returns_200(self, client):
        resp = client.post("/population", json={
            "scenario": {
                "name": "cohort_test", "duration_s": 3600,
                "meal_events": [{"time_s": 0, "carbs_g": 60.0}],
            },
            "n_individuals": 3,
        })
        assert resp.status_code == 200

    def test_response_shape(self, client):
        data = client.post("/population", json={
            "scenario": {
                "name": "cohort_shape", "duration_s": 3600,
                "meal_events": [{"time_s": 0, "carbs_g": 60.0}],
            },
            "n_individuals": 3,
        }).json()
        assert data["n_total"] == 3
        assert data["n_succeeded"] + data["n_failed"] == 3
        assert isinstance(data["glucose_peak_mean_mM"], float)
        assert isinstance(data["failure_reasons"], list)

    def test_same_random_seed_is_reproducible(self, client):
        # Cohort generation is seeded (PopulationSpec.random_seed) --
        # the same request must produce the same cohort composition and
        # therefore the same summary statistics, not a fresh random draw
        # per HTTP call.
        payload = {
            "scenario": {
                "name": "cohort_seed", "duration_s": 3600,
                "meal_events": [{"time_s": 0, "carbs_g": 60.0}],
            },
            "n_individuals": 3, "random_seed": 7,
        }
        first = client.post("/population", json=payload).json()
        second = client.post("/population", json=payload).json()
        assert first == second

    def test_n_individuals_above_cap_returns_422(self, client):
        resp = client.post("/population", json={
            "scenario": {"name": "too_big", "meal_events": []},
            "n_individuals": 21,
        })
        assert resp.status_code == 422

    def test_n_individuals_zero_returns_422(self, client):
        resp = client.post("/population", json={
            "scenario": {"name": "too_small", "meal_events": []},
            "n_individuals": 0,
        })
        assert resp.status_code == 422

    def test_invalid_scenario_combination_returns_400_not_500(self, client):
        # Mirrors /simulate's own handling of a schema-valid-but-
        # semantically-invalid scenario (ScenarioDefinition's own
        # insulin_dose_events/DrugEvent mutual-exclusivity check).
        resp = client.post("/population", json={
            "scenario": {
                "name": "bad", "meal_events": [],
                "drug_events": [{"time_s": 0, "drug": "insulin_regular", "dose_mg": 5.0}],
                "insulin_dose_events": [
                    {"time_s": 0, "insulin_type": "regular", "dose_units": 5.0},
                ],
            },
            "n_individuals": 2,
        })
        assert resp.status_code == 400

    def test_run_uncertainty_rejected_with_400(self, client):
        # CRITICAL, mandatory code-reviewer finding: n_mc_samples has no
        # upper bound, and run_uncertainty runs a full Monte Carlo
        # ensemble PER individual -- up to n_individuals times per
        # request. Rejected outright rather than silently allowed.
        resp = client.post("/population", json={
            "scenario": {
                "name": "uq_cohort", "meal_events": [{"time_s": 0, "carbs_g": 60.0}],
                "run_uncertainty": True, "n_mc_samples": 10,
            },
            "n_individuals": 2,
        })
        assert resp.status_code == 400
        assert "run_uncertainty" in resp.json()["detail"]

    def test_n_individuals_at_cap_boundary_actually_runs(self, client):
        # The cap (20) itself was previously only exercised at 21
        # (rejected) and small values like 2-3 (accepted) -- never
        # confirmed to actually run and produce 20 real results at the
        # accepted boundary itself.
        data = client.post("/population", json={
            "scenario": {
                "name": "cohort_boundary", "duration_s": 3600,
                "meal_events": [{"time_s": 0, "carbs_g": 60.0}],
            },
            "n_individuals": 20,
        }).json()
        assert data["n_total"] == 20

    def test_never_uses_multiprocessing(self, client):
        # The actual safety property this endpoint exists to guarantee
        # (MASTER_PLAN 3.11's own resource-exhaustion finding) --
        # protected only by PopulationRunner's own default parameter
        # value until now. Asserts the real constructed instance's
        # state directly, not just how it happened to be called.
        from unittest.mock import patch

        from biosim_seeds import PopulationRunner

        real_init = PopulationRunner.__init__
        with patch.object(PopulationRunner, "__init__", autospec=True, side_effect=real_init) as mock_init:
            client.post("/population", json={
                "scenario": {
                    "name": "cohort_mp_check", "duration_s": 3600,
                    "meal_events": [{"time_s": 0, "carbs_g": 60.0}],
                },
                "n_individuals": 2,
            })

        assert mock_init.call_count == 1
        instance = mock_init.call_args.args[0]
        assert instance.use_multiprocessing is False


class TestBatch:
    """POST /batch -- MASTER_PLAN 3.11's bounded, single-threaded batch
    endpoint (BatchRunner run with max_workers=1/ThreadPoolExecutor, so
    the REST layer itself never spawns an OS process pool)."""

    def _scenario(self, name, carbs_g):
        return {
            "name": name, "duration_s": 3600,
            "meal_events": [{"time_s": 0, "carbs_g": carbs_g}],
        }

    def test_returns_200(self, client):
        resp = client.post("/batch", json={
            "scenarios": [self._scenario("b1", 40.0), self._scenario("b2", 80.0)],
        })
        assert resp.status_code == 200

    def test_response_shape_and_order_preserved(self, client):
        data = client.post("/batch", json={
            "scenarios": [self._scenario("first", 30.0), self._scenario("second", 90.0)],
        }).json()
        assert data["n_succeeded"] == 2
        assert data["n_failed"] == 0
        assert [r["name"] for r in data["results"]] == ["first", "second"]
        assert all(r["succeeded"] for r in data["results"])
        assert all(r["result"] is not None for r in data["results"])
        # A larger meal must produce a higher glucose peak -- confirms
        # each result actually corresponds to its own scenario, not a
        # shared/duplicated one.
        first_peak = data["results"][0]["result"]["expert_metrics"]["glucose_peak_mM"]
        second_peak = data["results"][1]["result"]["expert_metrics"]["glucose_peak_mM"]
        assert second_peak > first_peak

    def test_duplicate_scenario_names_return_400(self, client):
        resp = client.post("/batch", json={
            "scenarios": [self._scenario("dup", 30.0), self._scenario("dup", 90.0)],
        })
        assert resp.status_code == 400

    def test_empty_scenario_list_returns_422(self, client):
        resp = client.post("/batch", json={"scenarios": []})
        assert resp.status_code == 422

    def test_too_many_scenarios_returns_422(self, client):
        resp = client.post("/batch", json={
            "scenarios": [self._scenario(f"s{i}", 30.0) for i in range(21)],
        })
        assert resp.status_code == 422

    def test_invalid_scenario_in_batch_returns_400_not_500(self, client):
        resp = client.post("/batch", json={
            "scenarios": [
                self._scenario("ok", 30.0),
                {
                    "name": "bad", "meal_events": [],
                    "drug_events": [{"time_s": 0, "drug": "insulin_regular", "dose_mg": 5.0}],
                    "insulin_dose_events": [
                        {"time_s": 0, "insulin_type": "regular", "dose_units": 5.0},
                    ],
                },
            ],
        })
        assert resp.status_code == 400

    def test_run_uncertainty_rejected_with_400(self, client):
        # Same CRITICAL finding as TestPopulation's own version -- up to
        # 20 scenarios here, each independently able to request an
        # unbounded n_mc_samples ensemble.
        scenario = self._scenario("uq_batch", 60.0)
        scenario["run_uncertainty"] = True
        scenario["n_mc_samples"] = 10
        resp = client.post("/batch", json={"scenarios": [scenario]})
        assert resp.status_code == 400
        assert "run_uncertainty" in resp.json()["detail"]

    def test_twenty_scenarios_at_cap_boundary_actually_runs(self, client):
        data = client.post("/batch", json={
            "scenarios": [self._scenario(f"s{i}", 30.0 + i) for i in range(20)],
        }).json()
        assert data["n_succeeded"] == 20
        assert len(data["results"]) == 20

    def test_never_spawns_a_process_pool(self, client):
        # The actual safety property this endpoint exists to guarantee
        # -- protected only by app.py's hardcoded constructor arguments
        # until now. Asserts the real constructed instance's state
        # directly, so a future edit that dropped executor_cls or
        # bumped max_workers would fail this test even though it
        # wouldn't be visible from any response body.
        from concurrent.futures import ThreadPoolExecutor
        from unittest.mock import patch

        from biosim_expert import BatchRunner

        real_init = BatchRunner.__init__
        with patch.object(BatchRunner, "__init__", autospec=True, side_effect=real_init) as mock_init:
            client.post("/batch", json={
                "scenarios": [self._scenario("mp_check", 40.0)],
            })

        assert mock_init.call_count == 1
        instance = mock_init.call_args.args[0]
        assert instance.executor_cls is ThreadPoolExecutor
        assert instance.max_workers == 1


class TestUnverifiedVariantWarningInResponse:
    """The REST request schema takes no rsIDs (rule 7), but a SimulationOutput
    carrying the rs35767 note must pass it through to the response."""

    def test_note_reaches_rest_warnings(self):
        from biosim_api import PersonalizationSpec, SimulationRunner
        from biosim_api import high_sugar_meal_scenario
        from biosim_genotype import Genotype, Variant
        from biosim_genotype.effect_db import UNVERIFIED_VARIANT_NOTES
        from biosim_rest.service import simulation_output_to_response
        scenario = high_sugar_meal_scenario(
            personalization=PersonalizationSpec(
                genetic_variants=[Variant("rs35767", Genotype.HETEROZYGOUS)]
            )
        )
        response = simulation_output_to_response(SimulationRunner().run(scenario))
        assert UNVERIFIED_VARIANT_NOTES["rs35767"] in response.warnings


class TestExperimentalDrugLabel:
    def test_simulate_with_drug_returns_experimental_warning(self, client):
        from biosim_physiology.drug_interactions import Drug
        from biosim_physiology.drug_interactions.validation_status import (
            EXPERIMENTAL_DRUG_NOTES,
        )
        resp = client.post("/simulate", json={
            "name": "metformin_meal",
            "meal_events": [{"time_s": 0, "carbs_g": 50.0}],
            "drug_events": [{"time_s": 0, "drug": "metformin", "dose_mg": 1000.0}],
        })
        assert resp.status_code == 200
        assert EXPERIMENTAL_DRUG_NOTES[Drug.METFORMIN] in resp.json()["warnings"]

    def test_validation_endpoint_labels_drug_benchmarks(self, client):
        datasets = {d["dataset_name"]: d for d in client.get("/validation").json()["datasets"]}
        assert "experimental module" in datasets["Metformin_T2D_postprandial"]["module_label"]
        assert datasets["OGTT_75g_healthy"]["module_label"] == ""

    def test_population_with_drug_returns_experimental_warning(self, client):
        from biosim_physiology.drug_interactions import Drug
        from biosim_physiology.drug_interactions.validation_status import (
            EXPERIMENTAL_DRUG_NOTES,
        )
        resp = client.post("/population", json={
            "scenario": {
                "name": "metformin_pop",
                "meal_events": [{"time_s": 0, "carbs_g": 50.0}],
                "drug_events": [{"time_s": 0, "drug": "metformin", "dose_mg": 1000.0}],
                "duration_s": 3600,
            },
            "n_individuals": 2,
        })
        assert resp.status_code == 200, resp.text
        assert EXPERIMENTAL_DRUG_NOTES[Drug.METFORMIN] in resp.json()["warnings"]
