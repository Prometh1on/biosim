# SAVE TO: packages\literature-pipeline\tests\test_llm_extractor.py
"""
Tests for LLMExtractor's parsing and error-handling logic.

These tests never call the real Claude API -- the anthropic.Anthropic
client is fully mocked. This validates the conversion logic (building
RawExtraction from the parsed structured-output payload, refusal handling,
API error handling) without requiring credentials or network access, and
without incurring any real API cost.
"""
from unittest.mock import MagicMock, patch

import httpx
import pytest

from biosim_literature.extraction import PaperMetadata
from biosim_literature.llm_extractor import LLMExtractor, DEFAULT_MODEL, _ExtractedParameter, _ExtractionPayload


def _paper(**overrides) -> PaperMetadata:
    defaults = dict(
        pubmed_id="12345678",
        doi="10.1234/test",
        title="Kinetics of the insulin receptor",
        abstract="The Km for insulin binding was 0.45 nM. kcat was 0.023 s-1.",
        authors=["Smith J"],
        year=2024,
        journal="J Biol Chem",
    )
    defaults.update(overrides)
    return PaperMetadata(**defaults)


def _fake_response(stop_reason="end_turn", parameters=None):
    response = MagicMock()
    response.stop_reason = stop_reason
    response.parsed_output = _ExtractionPayload(parameters=parameters or [])
    return response


class TestLLMExtractorBasics:
    def test_default_model(self):
        assert LLMExtractor()._model == DEFAULT_MODEL

    def test_name_reflects_model(self):
        assert LLMExtractor(model="claude-sonnet-5").name == "llm-extractor-claude-sonnet-5"

    def test_custom_model_and_max_tokens(self):
        e = LLMExtractor(model="claude-haiku-4-5", max_tokens=1024)
        assert e._model == "claude-haiku-4-5"
        assert e._max_tokens == 1024


class TestLLMExtractorSuccess:
    def test_successful_extraction_builds_raw_extractions(self):
        param = _ExtractedParameter(
            parameter_kind="Km", subject_name="insulin receptor",
            value=0.45, unit="nM", organism="human", extractor_confidence=0.9,
        )
        with patch("anthropic.Anthropic") as MockClient:
            MockClient.return_value.messages.parse.return_value = _fake_response(parameters=[param])
            result = LLMExtractor().extract(_paper())

        assert not result.errors
        assert len(result.extractions) == 1
        raw = result.extractions[0]
        assert raw.parameter_kind == "Km"
        assert raw.subject_name == "insulin receptor"
        assert raw.value == pytest.approx(0.45)
        assert raw.unit == "nM"
        assert raw.organism == "human"
        assert raw.extractor_confidence == pytest.approx(0.9)

    def test_empty_extraction_list_is_valid(self):
        with patch("anthropic.Anthropic") as MockClient:
            MockClient.return_value.messages.parse.return_value = _fake_response(parameters=[])
            result = LLMExtractor().extract(_paper())
        assert result.extractions == []
        assert not result.errors

    def test_uses_explicit_api_key_when_given(self):
        with patch("anthropic.Anthropic") as MockClient:
            MockClient.return_value.messages.parse.return_value = _fake_response()
            LLMExtractor(api_key="sk-ant-test").extract(_paper())
        MockClient.assert_called_once_with(api_key="sk-ant-test")

    def test_no_explicit_key_uses_env_resolution(self):
        with patch("anthropic.Anthropic") as MockClient:
            MockClient.return_value.messages.parse.return_value = _fake_response()
            LLMExtractor().extract(_paper())
        MockClient.assert_called_once_with()

    def test_includes_full_text_when_present(self):
        captured = {}

        def fake_parse(**kwargs):
            captured.update(kwargs)
            return _fake_response()

        with patch("anthropic.Anthropic") as MockClient:
            MockClient.return_value.messages.parse.side_effect = fake_parse
            LLMExtractor().extract(_paper(full_text="Full methods section text."))

        assert "Full methods section text." in captured["messages"][0]["content"]


class TestLLMExtractorRefusal:
    def test_refusal_returns_no_extractions_with_error(self):
        with patch("anthropic.Anthropic") as MockClient:
            MockClient.return_value.messages.parse.return_value = _fake_response(stop_reason="refusal")
            result = LLMExtractor().extract(_paper())
        assert result.extractions == []
        assert result.errors
        assert "refus" in result.errors[0].lower()


class TestLLMExtractorErrorHandling:
    def test_api_status_error_captured_as_extraction_error(self):
        import anthropic

        request = httpx.Request("POST", "https://api.anthropic.com/v1/messages")
        http_response = httpx.Response(status_code=429, request=request, json={"error": {"message": "rate limited"}})
        err = anthropic.APIStatusError("rate limited", response=http_response, body=None)

        with patch("anthropic.Anthropic") as MockClient:
            MockClient.return_value.messages.parse.side_effect = err
            result = LLMExtractor().extract(_paper())

        assert result.extractions == []
        assert result.errors
        assert "429" in result.errors[0]

    def test_api_connection_error_captured_as_extraction_error(self):
        import anthropic

        request = httpx.Request("POST", "https://api.anthropic.com/v1/messages")
        err = anthropic.APIConnectionError(request=request)

        with patch("anthropic.Anthropic") as MockClient:
            MockClient.return_value.messages.parse.side_effect = err
            result = LLMExtractor().extract(_paper())

        assert result.extractions == []
        assert result.errors
        assert "network error" in result.errors[0].lower()

    def test_does_not_crash_pipeline_on_error(self):
        """LLMExtractor.extract() must return a result, never raise, so
        IngestionPipeline's per-paper isolation guarantee holds."""
        import anthropic

        request = httpx.Request("POST", "https://api.anthropic.com/v1/messages")
        with patch("anthropic.Anthropic") as MockClient:
            MockClient.return_value.messages.parse.side_effect = anthropic.APIConnectionError(request=request)
            result = LLMExtractor().extract(_paper())
        assert result is not None
