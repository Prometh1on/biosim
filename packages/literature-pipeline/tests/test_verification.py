# SAVE TO: packages\literature-pipeline\tests\test_verification.py
"""
Test suite for biosim_literature.verification.

All tests use MockFetcher -- no live network calls, matching this
project's own established convention for the offline pytest suite
(see test_fetcher.py's own module docstring).
"""
import pytest

from biosim_literature.extraction import PaperMetadata
from biosim_literature.fetcher import Fetcher, FetchError, MockFetcher
from biosim_literature.verification import ClaimCheck, CitationVerification, verify_citation


class _FailingFetcher(Fetcher):
    """Raises FetchError on every call -- used to confirm verify_citation()
    propagates a fetch failure rather than silently swallowing it."""

    def fetch_by_pubmed_id(self, pubmed_id: str) -> PaperMetadata:
        raise FetchError(f"simulated failure for {pubmed_id}")

    def fetch_by_doi(self, doi: str) -> PaperMetadata:
        raise FetchError(f"simulated failure for {doi}")


def _paper(**overrides) -> PaperMetadata:
    defaults = dict(
        pubmed_id="31695455",
        doi="10.2147/DMSO.S222282",
        title="The relationship of mitochondrial dysfunction and the "
              "development of insulin resistance in Cushing's syndrome",
        abstract=(
            "Cushing's syndrome is characterized by metabolic disturbances "
            "including insulin resistance. Mitochondrial dysfunction is one "
            "pathogenic factor in the development of insulin resistance."
        ),
        authors=["Jezkova J", "Durovcova V", "Wenchich L"],
        year=2019,
        journal="Diabetes Metab Syndr Obes",
        is_primary_source=True,
    )
    defaults.update(overrides)
    return PaperMetadata(**defaults)


class TestVerifyCitationAbstractOnly:
    def test_claim_found_in_abstract_is_reported_correctly(self):
        fetcher = MockFetcher(papers={"31695455": _paper()})
        result = verify_citation(
            "31695455", ["insulin resistance"], fetcher=fetcher,
        )
        assert result.all_claims_verified
        check = result.claim_checks[0]
        assert check.found is True
        assert check.found_in == "abstract"
        assert check.context is not None

    def test_claim_absent_everywhere_is_reported_not_found(self):
        fetcher = MockFetcher(papers={"31695455": _paper()})
        result = verify_citation(
            "31695455", ["822"], fetcher=fetcher,
        )
        assert result.all_claims_verified is False
        assert result.unverified_claims == ["822"]
        check = result.claim_checks[0]
        assert check.found is False
        assert check.found_in is None
        assert check.context is None

    def test_full_text_available_is_false_when_not_fetched(self):
        fetcher = MockFetcher(papers={"31695455": _paper()})
        result = verify_citation("31695455", ["insulin"], fetcher=fetcher)
        assert result.full_text_available is False


class TestVerifyCitationFullText:
    def test_claim_only_in_full_text_is_found_and_labeled_correctly(self):
        # Reproduces the actual real-world case this tool exists for:
        # PMID:31695455's Table 1 cortisol figures (527/822 nmol/L) are
        # absent from the abstract (confirmed live -- the abstract is
        # entirely about mitochondrial/OXPHOS gene expression) but
        # present in the full text body.
        full_text = (
            "Table 1: Plasma cortisol at 08:00 hrs (nmol/L) "
            "527+/-51 822+/-80* 286+/-71"
        )
        fetcher = MockFetcher(papers={"31695455": _paper(full_text=full_text)})
        result = verify_citation("31695455", ["822", "527"], fetcher=fetcher)
        assert result.all_claims_verified
        assert result.full_text_available is True
        by_claim = {c.claim: c for c in result.claim_checks}
        assert by_claim["822"].found_in == "full_text"
        assert by_claim["527"].found_in == "full_text"
        assert "822" in by_claim["822"].context

    def test_abstract_checked_before_full_text(self):
        # A claim present in BOTH should report "abstract" (the
        # cheaper, always-available source), not "full_text" -- the
        # search order matters and this pins it as a regression.
        fetcher = MockFetcher(papers={
            "31695455": _paper(
                abstract="insulin resistance is well documented",
                full_text="insulin resistance also appears here in full text",
            )
        })
        result = verify_citation("31695455", ["insulin resistance"], fetcher=fetcher)
        assert result.claim_checks[0].found_in == "abstract"

    def test_full_text_unavailable_claim_not_found_distinguishable(self):
        # Not-found-because-genuinely-absent and not-found-because-
        # full-text-was-unavailable must be distinguishable via
        # full_text_available, not conflated into one "not found" bucket.
        fetcher = MockFetcher(papers={"31695455": _paper()})  # no full_text
        result = verify_citation("31695455", ["822"], fetcher=fetcher)
        assert result.claim_checks[0].found is False
        assert result.full_text_available is False


class TestCitationVerificationReport:
    def test_metadata_fields_are_carried_through(self):
        fetcher = MockFetcher(papers={"31695455": _paper()})
        result = verify_citation("31695455", ["insulin"], fetcher=fetcher)
        assert result.pubmed_id == "31695455"
        assert "Cushing" in result.title
        assert result.journal == "Diabetes Metab Syndr Obes"
        assert result.year == 2019
        assert result.authors[0] == "Jezkova J"

    def test_as_text_includes_found_and_not_found_markers(self):
        fetcher = MockFetcher(papers={"31695455": _paper()})
        result = verify_citation(
            "31695455", ["insulin resistance", "822"], fetcher=fetcher,
        )
        text = result.as_text()
        assert "FOUND" in text
        assert "NOT FOUND" in text
        assert "PMID:31695455" in text

    def test_multiple_claims_all_verified_true_only_if_every_claim_found(self):
        fetcher = MockFetcher(papers={"31695455": _paper()})
        result = verify_citation(
            "31695455", ["insulin resistance", "Cushing's syndrome"], fetcher=fetcher,
        )
        assert result.all_claims_verified is True
        assert result.unverified_claims == []


class TestVerifyCitationFetchFailure:
    def test_fetch_error_propagates_not_swallowed(self):
        # An unfetchable PMID must fail loudly (raise), never silently
        # report every claim as "not found" -- that would be
        # indistinguishable from a genuinely wrong citation and defeat
        # this project's own evidence-hierarchy discipline.
        with pytest.raises(FetchError):
            verify_citation("00000000", ["anything"], fetcher=_FailingFetcher())
