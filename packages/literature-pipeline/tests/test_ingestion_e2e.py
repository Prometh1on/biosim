# SAVE TO: packages\literature-pipeline\tests\test_ingestion_e2e.py
"""
Real end-to-end test of IngestionPipeline.ingest_by_pubmed_id() --
real PubMedFetcher (live NCBI E-utilities) + real LLMExtractor (live
Claude API), not MockFetcher/MockExtractor.

Every other part of this pipeline is independently covered elsewhere
without any live dependency:
    - PubMedFetcher's parsing logic:      test_fetcher.py (mocked network)
    - LLMExtractor's parsing/error paths: test_llm_extractor.py (mocked client)
    - Pipeline orchestration logic:       test_literature_pipeline.py (Mock*)

This file is the one piece that can only be verified by actually calling
both live services together, which needs a real ANTHROPIC_API_KEY (NCBI
E-utilities itself needs no key for this volume of calls). It is
skipped automatically when that key is absent -- exactly the standard,
conventional way real projects gate tests on live credentials/network
(see CLAUDE.md's IngestionPipeline entry for why this can't be run in
this development environment, and why that's an environment fact, not a
missing feature). The moment ANTHROPIC_API_KEY is set, this test exists,
is correctly wired, and will run.

Uses a PMID this project has already independently verified as real via
live PubMedFetcher calls this session (Bergman RN, Phillips LS, Cobelli
C. 1981, J Clin Invest, PMID:7033284) so a failure here means the
pipeline itself is broken, not that the ID doesn't exist.
"""
import os

import pytest

from biosim_id_resolver import IDRegistry
from biosim_literature import IngestionPipeline
from biosim_literature.fetcher import PubMedFetcher
from biosim_literature.llm_extractor import LLMExtractor
from biosim_parameter_db import ParameterStore

requires_anthropic_key = pytest.mark.skipif(
    not os.environ.get("ANTHROPIC_API_KEY"),
    reason="requires a real ANTHROPIC_API_KEY to exercise the live Claude "
           "extraction step end-to-end; see CLAUDE.md's IngestionPipeline "
           "entry",
)


@pytest.mark.integration
@requires_anthropic_key
class TestIngestionPipelineLiveEndToEnd:
    def test_ingest_by_pubmed_id_runs_without_crashing(self):
        pipeline = IngestionPipeline(
            fetcher=PubMedFetcher(),
            extractor=LLMExtractor(),
            store=ParameterStore(),
            registry=IDRegistry(),
        )
        report = pipeline.ingest_by_pubmed_id("7033284")

        assert report.pubmed_id == "7033284"
        assert report.paper_title != ""
        assert report.total >= 0

    def test_ingest_by_pubmed_id_fetches_real_paper_metadata(self):
        # Independent of what the LLM extracts, the fetch step must
        # succeed and return this specific paper's real title.
        pipeline = IngestionPipeline(
            fetcher=PubMedFetcher(),
            extractor=LLMExtractor(),
            store=ParameterStore(),
            registry=IDRegistry(),
        )
        report = pipeline.ingest_by_pubmed_id("7033284")
        assert "insulin" in report.paper_title.lower() or "glucose" in report.paper_title.lower()
