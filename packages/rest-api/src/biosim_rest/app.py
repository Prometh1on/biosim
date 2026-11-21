# SAVE TO: packages\rest-api\src\biosim_rest\app.py
"""
Biosim REST API -- exposes SimulationRunner as HTTP endpoints.

Endpoints:
    POST /simulate   -- run a scenario, return expert + lay results
    GET  /scenarios  -- list the pre-built scenario library
    POST /compare    -- run two library scenarios and return their differential
    POST /population -- run a scenario through a small synthetic cohort
                        (bounded, always sequential -- MASTER_PLAN 3.11)
    POST /batch      -- run several scenarios independently (bounded,
                        always single-threaded -- MASTER_PLAN 3.11)
    GET  /health     -- service status, version, known test count
    GET  /validation -- current status of every pre-registered clinical
                        benchmark (OGTT, Bergman IVGTT, semaglutide,
                        steroid/T2D, metformin, SGLT2i, protein co-
                        ingestion, T2D meal insulin -- see
                        ModelValidator.validate_all())

Routes are intentionally thin: parse/validate the HTTP boundary, delegate
to biosim_api and service.py, translate the result. No business logic lives
here. Routes are plain `def` (not `async def`) because SimulationRunner.run()
is pure CPU-bound work (scipy ODE solving) with no I/O -- FastAPI runs sync
def routes in a worker thread pool, which is the correct choice for
CPU-bound handlers; async def would block the event loop instead.

/simulate is backed by biosim_infra.SimulationCache (thread-safe LRU,
keyed on the scenario's own JSON representation) -- the same caching
component `biosim_seeds.PlatformSession` already uses, but until now
this HTTP layer built a fresh, uncached SimulationRunner per request:
identical repeated requests (e.g. load-testing, or many users hitting
the same pre-built library scenario) re-ran the full ODE solve from
scratch every time despite working cache infrastructure already
existing elsewhere in the codebase.
"""

from __future__ import annotations

import math
from typing import Any

from fastapi import FastAPI, HTTPException, Request
from fastapi.encoders import jsonable_encoder
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse

from .schemas import (
    BatchRequest,
    BatchResponse,
    CohortSummaryOut,
    CompareRequest,
    CompareResponse,
    ErrorResponse,
    HealthResponse,
    PopulationRequest,
    ScenarioIn,
    ScenarioListItem,
    SimulateResponse,
    ValidationSummaryResponse,
)
from .service import (
    batch_result_to_response,
    cohort_summary_to_response,
    population_spec_from_request,
    scenario_in_to_definition,
    simulation_output_to_response,
    validation_reports_to_response,
)

APP_VERSION = "0.1.0"


def _sanitize_non_finite(value: Any) -> Any:
    """
    Recursively replace NaN/Infinity/-Infinity floats with a JSON-safe
    string (MASTER_PLAN 3.13).

    Regression: a raw JSON body containing a bare `NaN`/`Infinity`
    literal for a float field (e.g. `{"dose_mg": NaN}` -- Python's own
    `json.dumps(float('nan'))` emits exactly this non-standard literal
    by default, so a real client using Python's stdlib json module can
    genuinely produce one) is correctly rejected by Pydantic's own field
    validation (e.g. `Field(ge=0)` -- NaN/Inf both fail that comparison)
    -- but FastAPI's default `RequestValidationError` handler echoes the
    rejected value back into the 422 error body's `input` field
    verbatim, and Starlette's `JSONResponse.render()` calls
    `json.dumps(..., allow_nan=False)` explicitly, which raises
    `ValueError: Out of range float values are not JSON compliant`
    instead of ever producing a response -- turning a client mistake
    that Pydantic already correctly caught into an unhandled 500. This
    is generic to every float field in this API, not specific to any
    one endpoint -- confirmed by reproducing the identical crash against
    both `DrugEventIn.dose_mg` and `InsulinDoseEventIn.dose_units`.
    """
    if isinstance(value, float) and not math.isfinite(value):
        if math.isnan(value):
            return "NaN"
        return "Infinity" if value > 0 else "-Infinity"
    if isinstance(value, dict):
        return {k: _sanitize_non_finite(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_sanitize_non_finite(v) for v in value]
    return value


def create_app() -> FastAPI:
    app = FastAPI(
        title="Biosim REST API",
        description="Postprandial glucose/insulin simulation platform.",
        version=APP_VERSION,
    )

    from biosim_infra import SimulationCache
    cache = SimulationCache()

    @app.exception_handler(RequestValidationError)
    async def _sanitized_validation_error_handler(
        request: Request, exc: RequestValidationError,
    ) -> JSONResponse:
        # Mirrors FastAPI's own default RequestValidationError handler
        # (status 422, {"detail": errors}) exactly, except the errors
        # are sanitized first (see _sanitize_non_finite's docstring) so
        # a rejected NaN/Inf input value can't crash the response itself.
        return JSONResponse(
            status_code=422,
            content={"detail": _sanitize_non_finite(jsonable_encoder(exc.errors()))},
        )

    @app.get("/health", response_model=HealthResponse)
    def health() -> HealthResponse:
        return HealthResponse(status="ok", version=APP_VERSION)

    @app.get("/validation", response_model=ValidationSummaryResponse)
    def validation() -> ValidationSummaryResponse:
        from biosim_validation import ModelValidator

        reports = ModelValidator().validate_all()
        return validation_reports_to_response(reports)

    @app.get("/scenarios", response_model=list[ScenarioListItem])
    def list_scenarios() -> list[ScenarioListItem]:
        from biosim_api import LIBRARY

        return [
            ScenarioListItem(name=name, description=s.description, duration_s=s.duration_s)
            for name, s in sorted(LIBRARY.items())
        ]

    @app.post(
        "/simulate",
        response_model=SimulateResponse,
        responses={400: {"model": ErrorResponse}},
    )
    def simulate(payload: ScenarioIn) -> SimulateResponse:
        from biosim_api import SimulationRunner

        try:
            scenario = scenario_in_to_definition(payload)
        except ValueError as e:
            raise HTTPException(status_code=400, detail=str(e))

        output = cache.get(scenario)
        if output is None:
            try:
                output = SimulationRunner().run(scenario)
            except ValueError as e:
                # A schema-valid payload can still combine into an invalid
                # simulation input downstream of scenario construction (e.g.
                # an extreme hours_before_simulation overflowing into a
                # non-finite dose_time_min for a coupled drug pathway) --
                # treat that the same as any other input-validation failure
                # rather than letting it surface as an unhandled 500.
                raise HTTPException(status_code=400, detail=str(e))
            cache.put(scenario, output)
        return simulation_output_to_response(output)

    @app.post(
        "/compare",
        response_model=CompareResponse,
        responses={404: {"model": ErrorResponse}},
    )
    def compare(payload: CompareRequest) -> CompareResponse:
        from biosim_api import ComparisonRunner, get_scenario

        try:
            scenario_a = get_scenario(payload.scenario_a)
            scenario_b = get_scenario(payload.scenario_b)
        except KeyError as e:
            raise HTTPException(status_code=404, detail=str(e).strip('"'))

        result = ComparisonRunner().compare(scenario_a, scenario_b)
        return CompareResponse(
            scenario_a=simulation_output_to_response(result.scenario_a),
            scenario_b=simulation_output_to_response(result.scenario_b),
            peak_glucose_diff_mM=result.peak_glucose_diff_uM / 1000.0,
            peak_insulin_diff_pM=result.peak_insulin_diff_uM * 1e6,
            narrative=result.narrative,
        )

    @app.post(
        "/population",
        response_model=CohortSummaryOut,
        responses={400: {"model": ErrorResponse}},
    )
    def population(payload: PopulationRequest) -> CohortSummaryOut:
        from biosim_seeds import PopulationRunner

        if payload.scenario.run_uncertainty:
            # CRITICAL (mandatory code-reviewer pass, MASTER_PLAN 3.11):
            # ScenarioIn.n_mc_samples has no upper bound (Field(ge=10)
            # only), and run_uncertainty runs a full n_mc_samples-sample
            # Monte Carlo ensemble PER individual -- multiplying an
            # already-unbounded per-scenario cost by up to
            # n_individuals (20) turns this bounded-synchronous
            # endpoint back into the unbounded-duration compute job
            # this row's own resource-exhaustion finding was written to
            # prevent. Confirmed empirically: n_individuals=3,
            # n_mc_samples=50 alone took 2.8s; schema-legal values
            # extrapolate past 300s for a single request. Rejected
            # outright rather than silently capped, since no bound on
            # n_mc_samples itself exists yet anywhere in this API (a
            # separate, pre-existing gap -- see MASTER_PLAN's
            # cross-cutting housekeeping backlog).
            raise HTTPException(
                status_code=400,
                detail="run_uncertainty is not supported for scenarios submitted to "
                       "/population -- n_mc_samples has no upper bound, and this "
                       "endpoint runs it once per individual (up to n_individuals times)",
            )

        try:
            scenario = scenario_in_to_definition(payload.scenario)
        except ValueError as e:
            raise HTTPException(status_code=400, detail=str(e))

        spec = population_spec_from_request(payload)
        # use_multiprocessing left at PopulationRunner's own default
        # (False) -- deliberately never exposed over REST, see
        # PopulationRequest's docstring.
        summary = PopulationRunner(spec=spec).run(scenario)
        return cohort_summary_to_response(summary)

    @app.post(
        "/batch",
        response_model=BatchResponse,
        responses={400: {"model": ErrorResponse}},
    )
    def batch(payload: BatchRequest) -> BatchResponse:
        from concurrent.futures import ThreadPoolExecutor

        from biosim_expert import BatchRunner

        names = [s.name for s in payload.scenarios]
        if len(set(names)) != len(names):
            # BatchResult.get()/CSVExporter both look scenarios up by
            # name, matching only the first hit -- a duplicate-named
            # pair would silently misattribute the second scenario's
            # result to the first. Reject at the boundary instead of
            # producing a misleading response.
            raise HTTPException(
                status_code=400, detail="scenario names in a /batch request must be unique",
            )

        if any(s.run_uncertainty for s in payload.scenarios):
            # Same CRITICAL finding as /population (see its own comment
            # above), same multiplier shape: up to 20 scenarios, each
            # independently able to request an unbounded n_mc_samples
            # Monte Carlo ensemble.
            raise HTTPException(
                status_code=400,
                detail="run_uncertainty is not supported for scenarios submitted to "
                       "/batch -- n_mc_samples has no upper bound, and this endpoint "
                       "runs up to 20 scenarios' worth per request",
            )

        try:
            scenarios = [scenario_in_to_definition(s) for s in payload.scenarios]
        except ValueError as e:
            raise HTTPException(status_code=400, detail=str(e))

        # max_workers=1 + ThreadPoolExecutor: never spawns an OS
        # process, unlike BatchRunner's own ProcessPoolExecutor default
        # -- a synchronous HTTP route must not let concurrent requests
        # each spawn their own process pool (MASTER_PLAN 3.11's own
        # nested-parallelism resource-exhaustion finding).
        result = BatchRunner(max_workers=1, executor_cls=ThreadPoolExecutor).run_all(scenarios)
        return batch_result_to_response(scenarios, result)

    return app
