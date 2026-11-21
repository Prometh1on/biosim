# SAVE TO: packages\literature-pipeline\src\biosim_literature\llm_extractor.py
"""
LLMExtractor -- calls the Claude API to extract kinetic parameters from
paper text (abstract, or abstract + full text when available).

Requires an Anthropic credential (ANTHROPIC_API_KEY, or any other source
`anthropic.Anthropic()` resolves automatically -- see the Anthropic SDK's
own credential precedence). Never exercised in the automated test suite
(no credentials in CI, and each call has a real cost) -- see
test_llm_extractor.py for the parsing/error-handling unit tests, which
mock the Claude client entirely. MockExtractor remains what the rest of
the pipeline is tested against.

Uses structured outputs (client.messages.parse with a Pydantic schema) so
the response is guaranteed to match RawExtraction's shape -- no free-text
parsing of the model's answer.
"""

from __future__ import annotations

from typing import Optional

from pydantic import BaseModel, Field

from .extraction import Extractor, ExtractionResult, PaperMetadata, RawExtraction

# Regression: was "claude-opus-4-8", a stale model ID from an earlier
# naming scheme -- never caught because this path requires a real
# ANTHROPIC_API_KEY and is deliberately never exercised in the automated
# suite or by this project's own "no paid API keys" constraint (see
# CLAUDE.md), so a wrong ID here would only surface the first time
# someone actually wires in a key. Opus tier kept deliberately (not
# downgraded to Sonnet) -- literature-to-parameter extraction feeds a
# scientific parameter database, where extraction accuracy matters more
# than latency/cost.
DEFAULT_MODEL = "claude-opus-5"

EXTRACTION_SYSTEM_PROMPT = """\
You are a biomedical literature extraction assistant for a computational \
systems biology platform. Given a paper's title and abstract (and full \
text when provided), extract every quantitative kinetic or physiological \
parameter reported for glucose/insulin metabolism, receptor kinetics, or \
related pathways (e.g. Km, kcat, EC50, IC50, Vmax, rate constants, binding \
affinities, half-lives, clearance rates).

For each parameter found, report:
- parameter_kind: a short machine-readable name for the kind of parameter \
  (e.g. "Km", "kcat", "EC50", "half_life")
- subject_name: the entity the parameter describes, as named in the paper \
  (e.g. "insulin receptor", "GLUT4") -- do not invent a canonical ID, use \
  the paper's own terminology
- value and unit: the numeric value and its unit exactly as reported
- object_name: a second entity if the parameter describes an interaction \
  (e.g. binding affinity of A for B), otherwise omit
- organism, cell_type, tissue, experimental_system: experimental context, \
  using "unknown" for organism/experimental_system if not stated
- temperature_celsius: if the measurement conditions state a temperature
- figure_or_table: the figure/table number the value came from, if cited
- raw_text: the exact sentence(s) the value was extracted from
- extractor_confidence: your own confidence (0.0-1.0) that this value was \
  extracted correctly from this specific text

Only extract values that are explicitly stated with a number and unit. Do \
not infer, estimate, or extract values from other papers you may know \
about -- extract only what this text actually states. If no such \
parameters are present, return an empty list.
"""


class _ExtractedParameter(BaseModel):
    parameter_kind: str
    subject_name: str
    value: float
    unit: str
    object_name: Optional[str] = None
    organism: str = "unknown"
    cell_type: Optional[str] = None
    tissue: Optional[str] = None
    experimental_system: str = "unknown"
    temperature_celsius: Optional[float] = None
    figure_or_table: Optional[str] = None
    raw_text: str = ""
    extractor_confidence: float = Field(default=0.7, ge=0.0, le=1.0)


class _ExtractionPayload(BaseModel):
    parameters: list[_ExtractedParameter]


class LLMExtractor(Extractor):
    """
    Extracts kinetic parameters from a paper using the Claude API.

    Usage:
        extractor = LLMExtractor()                        # ANTHROPIC_API_KEY from env
        extractor = LLMExtractor(api_key="sk-ant-...")     # explicit key
        extractor = LLMExtractor(model="claude-sonnet-5")  # cheaper model for high-volume runs
    """

    def __init__(
        self,
        model: str = DEFAULT_MODEL,
        api_key: Optional[str] = None,
        max_tokens: int = 4096,
    ) -> None:
        self._model = model
        self._api_key = api_key
        self._max_tokens = max_tokens

    @property
    def name(self) -> str:
        return f"llm-extractor-{self._model}"

    def extract(self, paper: PaperMetadata) -> ExtractionResult:
        import anthropic

        client = anthropic.Anthropic(api_key=self._api_key) if self._api_key else anthropic.Anthropic()

        text = paper.abstract
        if paper.full_text:
            text = f"{text}\n\n{paper.full_text}"

        try:
            response = client.messages.parse(
                model=self._model,
                max_tokens=self._max_tokens,
                system=EXTRACTION_SYSTEM_PROMPT,
                messages=[{
                    "role": "user",
                    "content": f"Title: {paper.title}\n\nText:\n{text}",
                }],
                output_format=_ExtractionPayload,
            )
        except anthropic.APIStatusError as e:
            return ExtractionResult(
                paper=paper, extractor_name=self.name,
                errors=[f"Claude API error ({e.status_code}): {e.message}"],
            )
        except anthropic.APIConnectionError as e:
            return ExtractionResult(
                paper=paper, extractor_name=self.name,
                errors=[f"Network error calling Claude API: {e}"],
            )

        if response.stop_reason == "refusal":
            return ExtractionResult(
                paper=paper, extractor_name=self.name,
                errors=["Extraction request was refused by the Claude API safety classifier"],
            )

        extractions = [
            RawExtraction(
                parameter_kind=p.parameter_kind,
                subject_name=p.subject_name,
                value=p.value,
                unit=p.unit,
                object_name=p.object_name,
                organism=p.organism,
                cell_type=p.cell_type,
                tissue=p.tissue,
                experimental_system=p.experimental_system,
                temperature_celsius=p.temperature_celsius,
                figure_or_table=p.figure_or_table,
                raw_text=p.raw_text,
                extractor_confidence=p.extractor_confidence,
            )
            for p in response.parsed_output.parameters
        ]
        return ExtractionResult(paper=paper, extractions=extractions, extractor_name=self.name)
