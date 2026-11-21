# SAVE TO: packages\literature-pipeline\tests\test_fetcher.py
"""
Test suite for biosim_literature.fetcher.

PubMedFetcher's real network calls are never exercised by the offline
pytest suite (correctly -- CI shouldn't depend on NCBI being reachable),
but that left its actual parsing logic (_parse_pubmed_xml, the DOI-search
fallback, rate limiting) with zero automated coverage despite being the
exact code path this project's citation-verification discipline leans on
every time a PMID is checked. These tests mock urllib.request.urlopen with
realistic NCBI E-utilities response shapes so the parsing logic itself is
pinned down without any live network dependency.
"""
from unittest.mock import patch, MagicMock

import pytest

from biosim_literature.fetcher import (
    FetchError,
    Fetcher,
    MockFetcher,
    PaperMetadata,
    PubMedFetcher,
)


# ---------------------------------------------------------------------------
# MockFetcher
# ---------------------------------------------------------------------------

class TestMockFetcher:
    def test_default_paper_for_unregistered_id(self):
        paper = MockFetcher().fetch_by_pubmed_id("12345")
        assert paper.pubmed_id == "12345"
        assert paper.is_primary_source is True

    def test_returns_registered_paper_by_pubmed_id(self):
        registered = PaperMetadata(
            pubmed_id="7033284", doi="10.1172/JCI110040", title="Registered paper",
            abstract="abstract", authors=["Bergman RN"], year=1981,
            journal="J Clin Invest", is_primary_source=True,
        )
        fetcher = MockFetcher(papers={"7033284": registered})
        assert fetcher.fetch_by_pubmed_id("7033284") is registered

    def test_returns_registered_paper_by_doi(self):
        registered = PaperMetadata(
            pubmed_id="7033284", doi="10.1172/JCI110040", title="Registered paper",
            abstract="abstract", authors=["Bergman RN"], year=1981,
            journal="J Clin Invest", is_primary_source=True,
        )
        fetcher = MockFetcher(papers={"10.1172/JCI110040": registered})
        assert fetcher.fetch_by_doi("10.1172/JCI110040") is registered

    def test_unregistered_doi_falls_back_to_default(self):
        paper = MockFetcher().fetch_by_doi("10.0000/unregistered")
        assert paper.pubmed_id == "99999999"

    def test_implements_fetcher_interface(self):
        assert isinstance(MockFetcher(), Fetcher)


# ---------------------------------------------------------------------------
# PubMedFetcher -- XML parsing (mocked network)
# ---------------------------------------------------------------------------

REALISTIC_PUBMED_XML = """<?xml version="1.0"?>
<PubmedArticleSet>
  <PubmedArticle>
    <MedlineCitation>
      <Article>
        <Journal>
          <Title>The Journal of clinical investigation</Title>
        </Journal>
        <ArticleTitle>Lack of insulin-mediated glucose disposal.</ArticleTitle>
        <Abstract>
          <AbstractText>First part of the abstract.</AbstractText>
          <AbstractText>Second part of the abstract.</AbstractText>
        </Abstract>
        <AuthorList>
          <Author>
            <LastName>Bergman</LastName>
            <Initials>RN</Initials>
          </Author>
          <Author>
            <LastName>Phillips</LastName>
            <Initials>LS</Initials>
          </Author>
        </AuthorList>
      </Article>
    </MedlineCitation>
    <PubmedData>
      <History>
        <PubMedPubDate><Year>1981</Year></PubMedPubDate>
      </History>
      <ArticleIdList>
        <ArticleId IdType="pubmed">7033284</ArticleId>
        <ArticleId IdType="doi">10.1172/JCI110040</ArticleId>
      </ArticleIdList>
    </PubmedData>
  </PubmedArticle>
</PubmedArticleSet>
"""

# NCBI puts PubDate/Year under MedlineCitation/Article/Journal/JournalIssue/
# PubDate in real responses -- _parse_pubmed_xml searches `.//PubDate/Year`
# (any depth), so nest it there for a realistic fixture.
REALISTIC_PUBMED_XML = REALISTIC_PUBMED_XML.replace(
    "<Journal>\n          <Title>",
    "<Journal>\n          <JournalIssue><PubDate><Year>1981</Year></PubDate></JournalIssue>\n"
    "          <Title>",
)


def _mock_urlopen_returning(text: str):
    response = MagicMock()
    response.read.return_value = text.encode("utf-8")
    response.__enter__.return_value = response
    response.__exit__.return_value = False
    return response


class TestPubMedFetcherParsing:
    def test_fetch_by_pubmed_id_parses_all_fields(self):
        fetcher = PubMedFetcher(rate_limit_delay=0.0)
        with patch("urllib.request.urlopen", return_value=_mock_urlopen_returning(REALISTIC_PUBMED_XML)):
            paper = fetcher.fetch_by_pubmed_id("7033284")

        assert paper.pubmed_id == "7033284"
        assert paper.title == "Lack of insulin-mediated glucose disposal."
        assert paper.abstract == "First part of the abstract. Second part of the abstract."
        assert paper.authors == ["Bergman RN", "Phillips LS"]
        assert paper.year == 1981
        assert paper.journal == "The Journal of clinical investigation"
        assert paper.doi == "10.1172/JCI110040"
        assert paper.is_primary_source is True

    def test_malformed_xml_raises_fetch_error(self):
        fetcher = PubMedFetcher(rate_limit_delay=0.0)
        with patch("urllib.request.urlopen", return_value=_mock_urlopen_returning("<not><valid</xml")):
            with pytest.raises(FetchError, match="Failed to parse"):
                fetcher.fetch_by_pubmed_id("123")

    def test_no_article_in_response_raises_fetch_error(self):
        fetcher = PubMedFetcher(rate_limit_delay=0.0)
        empty = "<PubmedArticleSet></PubmedArticleSet>"
        with patch("urllib.request.urlopen", return_value=_mock_urlopen_returning(empty)):
            with pytest.raises(FetchError, match="No article found"):
                fetcher.fetch_by_pubmed_id("123")

    def test_network_failure_raises_fetch_error(self):
        fetcher = PubMedFetcher(rate_limit_delay=0.0)
        with patch("urllib.request.urlopen", side_effect=OSError("timed out")):
            with pytest.raises(FetchError, match="Failed to fetch PubMed ID"):
                fetcher.fetch_by_pubmed_id("123")

    def test_missing_optional_fields_default_sensibly(self):
        # No AuthorList, no Journal/Title, no PubDate/Year, no doi ArticleId.
        minimal_xml = """<PubmedArticleSet><PubmedArticle>
            <MedlineCitation><Article>
                <ArticleTitle>Untitled</ArticleTitle>
            </Article></MedlineCitation>
        </PubmedArticle></PubmedArticleSet>"""
        fetcher = PubMedFetcher(rate_limit_delay=0.0)
        with patch("urllib.request.urlopen", return_value=_mock_urlopen_returning(minimal_xml)):
            paper = fetcher.fetch_by_pubmed_id("999")
        assert paper.title == "Untitled"
        assert paper.abstract == ""
        assert paper.authors == []
        assert paper.year is None
        assert paper.journal is None
        assert paper.doi is None

    def test_api_key_included_in_request_params(self):
        fetcher = PubMedFetcher(api_key="test-key-123", rate_limit_delay=0.0)
        with patch("urllib.request.urlopen", return_value=_mock_urlopen_returning(REALISTIC_PUBMED_XML)) as mocked:
            fetcher.fetch_by_pubmed_id("7033284")
        called_url = mocked.call_args[0][0]
        assert "api_key=test-key-123" in called_url

    def test_api_key_from_environment_variable(self, monkeypatch):
        monkeypatch.setenv("NCBI_API_KEY", "env-key-456")
        fetcher = PubMedFetcher(rate_limit_delay=0.0)
        assert fetcher._api_key == "env-key-456"

    def test_api_key_gives_faster_rate_limit(self):
        assert PubMedFetcher(api_key="k")._rate_limit_delay == 0.1
        assert PubMedFetcher(api_key=None, rate_limit_delay=0.34)._rate_limit_delay == 0.34

    def test_abstract_with_inline_markup_is_not_truncated(self):
        # Regression: `_parse_pubmed_xml` used `el.text` per <AbstractText>,
        # which in ElementTree returns only the text up to the FIRST child
        # element -- an abstract with inline <sub>/<i> markup (extremely
        # common for subscripted units like Tm_G or italicized P-values) was
        # silently truncated right at that tag, with no error or warning.
        # This fixture is the real <Abstract> block fetched live from NCBI
        # efetch for PMID:28428225 (Al-Jobori et al., Diabetes 2017) -- the
        # exact paper where this bug was discovered (see CLAUDE.md session
        # close note): the buggy parser returned only 309 of 1400+ real
        # characters, cutting off mid-sentence at "...maximal renal glucose
        # transport (Tm" right where the first <sub> tag appears.
        inline_markup_xml = REALISTIC_PUBMED_XML.replace(
            "<Abstract>\n"
            "          <AbstractText>First part of the abstract.</AbstractText>\n"
            "          <AbstractText>Second part of the abstract.</AbstractText>\n"
            "        </Abstract>",
            "<Abstract><AbstractText>Renal glucose reabsorption was measured "
            "with the stepped hyperglycemic clamp in 15 subjects with type 2 "
            "diabetes mellitus (T2DM) and 15 without diabetes after 2 days and "
            "after more chronic (14 days) treatment with empagliflozin. "
            "Patients with T2DM had significantly greater maximal renal "
            "glucose transport (Tm<sub>G</sub>) compared with subjects "
            "without diabetes at baseline (459 &#xb1; 53 vs. 337 &#xb1; 25 "
            "mg/min; <i>P</i> &lt; 0.05). Empagliflozin treatment for 48 h "
            "reduced the Tm<sub>G</sub> in both individuals with and without "
            "diabetes by 44 &#xb1; 7 and 53 &#xb1; 6%, respectively (both "
            "<i>P</i> &lt; 0.001).</AbstractText></Abstract>",
        )
        assert inline_markup_xml != REALISTIC_PUBMED_XML  # sanity: replace matched

        fetcher = PubMedFetcher(rate_limit_delay=0.0)
        with patch("urllib.request.urlopen", return_value=_mock_urlopen_returning(inline_markup_xml)):
            paper = fetcher.fetch_by_pubmed_id("28428225")

        assert "(TmG)" in paper.abstract
        assert paper.abstract.rstrip().endswith("respectively (both P < 0.001).")
        # The buggy parser truncated at exactly this length (309 chars,
        # right after "(Tm", the point the first <sub> tag appears) --
        # the real full text (with inline markup flattened) is 554 chars.
        assert len(paper.abstract) > 400


class TestPubMedFetcherDoiSearch:
    def test_fetch_by_doi_resolves_to_pubmed_id(self):
        fetcher = PubMedFetcher(rate_limit_delay=0.0)
        search_response = _mock_urlopen_returning(
            '{"esearchresult": {"idlist": ["7033284"]}}'
        )
        efetch_response = _mock_urlopen_returning(REALISTIC_PUBMED_XML)
        with patch("urllib.request.urlopen", side_effect=[search_response, efetch_response]):
            paper = fetcher.fetch_by_doi("10.1172/JCI110040")
        assert paper.pubmed_id == "7033284"

    def test_doi_not_found_raises_fetch_error(self):
        fetcher = PubMedFetcher(rate_limit_delay=0.0)
        empty_search = _mock_urlopen_returning('{"esearchresult": {"idlist": []}}')
        with patch("urllib.request.urlopen", return_value=empty_search):
            with pytest.raises(FetchError, match="No PubMed record found"):
                fetcher.fetch_by_doi("10.9999/nonexistent")

    def test_search_network_failure_raises_fetch_error(self):
        fetcher = PubMedFetcher(rate_limit_delay=0.0)
        with patch("urllib.request.urlopen", side_effect=OSError("timed out")):
            with pytest.raises(FetchError, match="Failed to search PubMed"):
                fetcher.fetch_by_doi("10.1172/JCI110040")


class TestPubMedFetcherFullText:
    """
    Regression: fetch_full_text/EUROPE_PMC_URL were declared (constructor
    param, class constant) but nothing ever called Europe PMC -- the
    module/class docstrings promised full-text enrichment that never
    happened. Verified against the real, live Europe PMC API by hand
    before writing these (a known open-access PLOS One paper and a known
    non-open-access one -- see CLAUDE.md) -- these tests pin the parsing
    logic down with mocked responses shaped exactly like those real ones,
    same convention as TestPubMedFetcherParsing above.
    """

    EUROPE_PMC_OA_SEARCH = (
        '{"resultList": {"result": [{"pmcid": "PMC7199958", '
        '"inEPMC": "Y", "isOpenAccess": "Y"}]}}'
    )
    EUROPE_PMC_NOT_OA_SEARCH = (
        '{"resultList": {"result": [{"inEPMC": "N", "isOpenAccess": "N"}]}}'
    )
    EUROPE_PMC_NO_RESULTS = '{"resultList": {"result": []}}'
    FULLTEXT_XML = (
        '<article><body><p>Real methods section text.</p>'
        '<p>Second paragraph.</p></body></article>'
    )

    def test_fetch_full_text_false_by_default_makes_no_europe_pmc_call(self):
        fetcher = PubMedFetcher(rate_limit_delay=0.0)
        with patch(
            "urllib.request.urlopen",
            return_value=_mock_urlopen_returning(REALISTIC_PUBMED_XML),
        ) as mocked:
            paper = fetcher.fetch_by_pubmed_id("7033284")
        assert paper.full_text is None
        assert mocked.call_count == 1  # only the efetch call

    def test_fetch_full_text_true_populates_full_text_when_open_access(self):
        fetcher = PubMedFetcher(fetch_full_text=True, rate_limit_delay=0.0)
        responses = [
            _mock_urlopen_returning(REALISTIC_PUBMED_XML),
            _mock_urlopen_returning(self.EUROPE_PMC_OA_SEARCH),
            _mock_urlopen_returning(self.FULLTEXT_XML),
        ]
        with patch("urllib.request.urlopen", side_effect=responses):
            paper = fetcher.fetch_by_pubmed_id("7033284")
        assert paper.full_text == "Real methods section text. Second paragraph."
        # The abstract-only fields must still be populated correctly --
        # full text is additive, not a replacement.
        assert paper.title == "Lack of insulin-mediated glucose disposal."

    def test_full_text_stays_none_when_not_open_access(self):
        fetcher = PubMedFetcher(fetch_full_text=True, rate_limit_delay=0.0)
        responses = [
            _mock_urlopen_returning(REALISTIC_PUBMED_XML),
            _mock_urlopen_returning(self.EUROPE_PMC_NOT_OA_SEARCH),
        ]
        with patch("urllib.request.urlopen", side_effect=responses) as mocked:
            paper = fetcher.fetch_by_pubmed_id("39049802")
        assert paper.full_text is None
        assert mocked.call_count == 2  # efetch + search, no fullTextXML attempt

    def test_full_text_stays_none_when_europe_pmc_has_no_results(self):
        fetcher = PubMedFetcher(fetch_full_text=True, rate_limit_delay=0.0)
        responses = [
            _mock_urlopen_returning(REALISTIC_PUBMED_XML),
            _mock_urlopen_returning(self.EUROPE_PMC_NO_RESULTS),
        ]
        with patch("urllib.request.urlopen", side_effect=responses):
            paper = fetcher.fetch_by_pubmed_id("123")
        assert paper.full_text is None

    def test_full_text_fetch_failure_degrades_gracefully(self):
        # A network error fetching full text must NOT break the overall
        # fetch_by_pubmed_id() call -- the abstract-only result is still
        # valid and full text was always best-effort, never guaranteed.
        fetcher = PubMedFetcher(fetch_full_text=True, rate_limit_delay=0.0)
        responses = [
            _mock_urlopen_returning(REALISTIC_PUBMED_XML),
            OSError("timed out"),
        ]
        with patch("urllib.request.urlopen", side_effect=responses):
            paper = fetcher.fetch_by_pubmed_id("7033284")
        assert paper.full_text is None
        assert paper.title == "Lack of insulin-mediated glucose disposal."

    def test_malformed_fulltext_xml_degrades_gracefully(self):
        fetcher = PubMedFetcher(fetch_full_text=True, rate_limit_delay=0.0)
        responses = [
            _mock_urlopen_returning(REALISTIC_PUBMED_XML),
            _mock_urlopen_returning(self.EUROPE_PMC_OA_SEARCH),
            _mock_urlopen_returning("<not><valid</xml"),
        ]
        with patch("urllib.request.urlopen", side_effect=responses):
            paper = fetcher.fetch_by_pubmed_id("7033284")
        assert paper.full_text is None


class TestPubMedFetcherPmcFallback:
    """
    Regression found by a code review of the hydroxychloroquine chunk:
    for PMID:26197707 Europe PMC's search reports inEPMC="Y"/inPMC="Y"
    with a PMCID (an NIH author manuscript, isOpenAccess="N"), but its
    REST fullTextXML endpoint serves nothing (HTTP 404/500) -- while NCBI's
    own PMC efetch (db=pmc) returns the full ~80KB article XML for free.
    The fetcher used to give up at that point and report "no full text",
    which made verify_citation() and a hand-written code comment both
    wrongly conclude the paper had no obtainable full text. It now falls
    back to NCBI PMC efetch whenever a PMCID is known.
    """

    SEARCH_MANUSCRIPT = (
        '{"resultList": {"result": [{"pmcid": "PMC4575248", '
        '"inEPMC": "Y", "inPMC": "Y", "isOpenAccess": "N"}]}}'
    )
    PMC_EFETCH_XML = (
        '<pmc-articleset><article><body><p>Discussion text from PMC.</p>'
        '<p>Unexplained decline.</p></body></article></pmc-articleset>'
    )

    @staticmethod
    def _http_error(code):
        import urllib.error
        return urllib.error.HTTPError("https://x", code, "err", {}, None)

    def test_falls_back_to_ncbi_pmc_efetch_when_europe_pmc_fulltext_fails(self):
        fetcher = PubMedFetcher(fetch_full_text=True, rate_limit_delay=0.0)
        responses = [
            _mock_urlopen_returning(REALISTIC_PUBMED_XML),
            _mock_urlopen_returning(self.SEARCH_MANUSCRIPT),
            self._http_error(404),                       # Europe PMC fullTextXML
            _mock_urlopen_returning(self.PMC_EFETCH_XML),  # NCBI PMC efetch
        ]
        with patch("urllib.request.urlopen", side_effect=responses):
            paper = fetcher.fetch_by_pubmed_id("26197707")
        assert paper.full_text == "Discussion text from PMC. Unexplained decline."

    def test_fallback_uses_the_ncbi_pmc_database_and_the_found_pmcid(self):
        fetcher = PubMedFetcher(fetch_full_text=True, rate_limit_delay=0.0)
        responses = [
            _mock_urlopen_returning(REALISTIC_PUBMED_XML),
            _mock_urlopen_returning(self.SEARCH_MANUSCRIPT),
            self._http_error(500),
            _mock_urlopen_returning(self.PMC_EFETCH_XML),
        ]
        with patch("urllib.request.urlopen", side_effect=responses) as mocked:
            fetcher.fetch_by_pubmed_id("26197707")
        last_url = mocked.call_args_list[-1].args[0]
        assert "efetch.fcgi" in last_url
        assert "db=pmc" in last_url
        assert "PMC4575248" in last_url

    def test_full_text_stays_none_when_both_sources_fail(self):
        fetcher = PubMedFetcher(fetch_full_text=True, rate_limit_delay=0.0)
        responses = [
            _mock_urlopen_returning(REALISTIC_PUBMED_XML),
            _mock_urlopen_returning(self.SEARCH_MANUSCRIPT),
            self._http_error(404),
            self._http_error(500),
        ]
        with patch("urllib.request.urlopen", side_effect=responses):
            paper = fetcher.fetch_by_pubmed_id("26197707")
        assert paper.full_text is None
        assert paper.title == "Lack of insulin-mediated glucose disposal."

    def test_no_fallback_attempted_when_no_pmcid_is_known(self):
        fetcher = PubMedFetcher(fetch_full_text=True, rate_limit_delay=0.0)
        responses = [
            _mock_urlopen_returning(REALISTIC_PUBMED_XML),
            _mock_urlopen_returning(
                '{"resultList": {"result": [{"inEPMC": "N", "inPMC": "N"}]}}'
            ),
        ]
        with patch("urllib.request.urlopen", side_effect=responses) as mocked:
            paper = fetcher.fetch_by_pubmed_id("8458526")
        assert paper.full_text is None
        assert mocked.call_count == 2  # efetch + search only


class TestPubMedFetcherRateLimit:
    def test_rate_limit_sleeps_when_called_too_soon(self):
        fetcher = PubMedFetcher(rate_limit_delay=1.0)
        fetcher._last_request_time = 1000.0
        with patch("time.time", return_value=1000.1), patch("time.sleep") as mocked_sleep:
            fetcher._rate_limit()
        mocked_sleep.assert_called_once()
        (slept_for,), _ = mocked_sleep.call_args
        assert slept_for == pytest.approx(0.9, abs=1e-6)

    def test_rate_limit_does_not_sleep_when_enough_time_elapsed(self):
        fetcher = PubMedFetcher(rate_limit_delay=0.34)
        fetcher._last_request_time = 1000.0
        with patch("time.time", return_value=1001.0), patch("time.sleep") as mocked_sleep:
            fetcher._rate_limit()
        mocked_sleep.assert_not_called()
