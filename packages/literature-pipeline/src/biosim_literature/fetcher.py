"""
Paper fetcher — retrieves paper metadata from PubMed (NCBI E-utilities).

The Fetcher interface is abstract, with two implementations:
    MockFetcher   — returns pre-defined paper metadata; used in tests
    PubMedFetcher — calls the NCBI E-utilities API; used in production

Rate limits:
    Without API key: 3 requests/second
    With API key:    10 requests/second
    Set NCBI_API_KEY environment variable to enable higher limits.

The fetcher never stores raw paper text server-side. Full text, when
available via Europe PMC, is fetched transiently and passed to the
extractor. It is not persisted anywhere in the platform.
"""

from __future__ import annotations

import dataclasses
import os
import time
import urllib.request
import urllib.parse
import json
import defusedxml.ElementTree as ET
from abc import ABC, abstractmethod
from typing import Optional

from .extraction import PaperMetadata


class FetchError(Exception):
    """Raised when a paper cannot be fetched."""
    pass


class Fetcher(ABC):
    """Abstract base class for paper fetchers."""

    @abstractmethod
    def fetch_by_pubmed_id(self, pubmed_id: str) -> PaperMetadata:
        """Fetch paper metadata by PubMed ID."""
        ...

    @abstractmethod
    def fetch_by_doi(self, doi: str) -> PaperMetadata:
        """Fetch paper metadata by DOI."""
        ...


class MockFetcher(Fetcher):
    """
    Mock fetcher for tests and offline development.

    Returns pre-defined PaperMetadata regardless of the requested ID.
    """

    def __init__(self, papers: Optional[dict[str, PaperMetadata]] = None):
        """
        Args:
            papers: Dict mapping pubmed_id or doi to PaperMetadata.
                    If None, returns a default test paper for any ID.
        """
        self._papers = papers or {}
        self._default_paper = PaperMetadata(
            pubmed_id="99999999",
            doi="10.0000/test.2024.001",
            title="Kinetics of the insulin receptor: a quantitative study",
            abstract=(
                "We report the kinetic parameters of the human insulin receptor. "
                "The Km for insulin binding was 0.45 nM. The kcat was 0.023 s-1. "
                "Measurements were performed in primary human hepatocytes at 37°C."
            ),
            authors=["Smith J", "Jones A", "Williams B"],
            year=2024,
            journal="Journal of Biological Chemistry",
            is_primary_source=True,
        )

    def fetch_by_pubmed_id(self, pubmed_id: str) -> PaperMetadata:
        if pubmed_id in self._papers:
            return self._papers[pubmed_id]
        paper = self._default_paper
        # Return a copy with the requested ID
        return PaperMetadata(
            pubmed_id=pubmed_id,
            doi=paper.doi,
            title=paper.title,
            abstract=paper.abstract,
            authors=paper.authors,
            year=paper.year,
            journal=paper.journal,
            is_primary_source=paper.is_primary_source,
        )

    def fetch_by_doi(self, doi: str) -> PaperMetadata:
        if doi in self._papers:
            return self._papers[doi]
        return self._default_paper


class PubMedFetcher(Fetcher):
    """
    Fetches paper metadata from NCBI E-utilities API.

    Optionally enriches with full text from Europe PMC for open-access papers.

    Usage:
        fetcher = PubMedFetcher()                    # no API key, 3 req/s
        fetcher = PubMedFetcher(api_key="your_key")  # 10 req/s

        # Or set NCBI_API_KEY environment variable
    """

    BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    EUROPE_PMC_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest"

    def __init__(
        self,
        api_key: Optional[str] = None,
        fetch_full_text: bool = False,
        rate_limit_delay: float = 0.34,  # 3 req/s default
    ) -> None:
        self._api_key = api_key or os.environ.get("NCBI_API_KEY")
        self._fetch_full_text = fetch_full_text
        self._rate_limit_delay = 0.1 if self._api_key else rate_limit_delay
        self._last_request_time = 0.0

    def fetch_by_pubmed_id(self, pubmed_id: str) -> PaperMetadata:
        """Fetch paper metadata from PubMed by PMID."""
        self._rate_limit()

        params = {
            "db": "pubmed",
            "id": pubmed_id,
            "retmode": "xml",
            "rettype": "abstract",
        }
        if self._api_key:
            params["api_key"] = self._api_key

        url = f"{self.BASE_URL}/efetch.fcgi?{urllib.parse.urlencode(params)}"

        try:
            with self._urlopen_https(url) as response:
                xml_content = response.read().decode("utf-8")
        except Exception as e:
            raise FetchError(f"Failed to fetch PubMed ID {pubmed_id}: {e}")

        metadata = self._parse_pubmed_xml(xml_content, pubmed_id)

        # Regression: fetch_full_text and EUROPE_PMC_URL were declared
        # (constructor param, class constant) but nothing ever actually
        # called Europe PMC -- this module's own docstring/class docstring
        # promised full-text enrichment that never happened. Full text is
        # best-effort: most papers aren't open access, so any failure here
        # (network, not found, no OA copy) degrades to full_text=None
        # rather than failing the whole fetch -- the abstract-only result
        # above is still valid and returned either way.
        if self._fetch_full_text:
            full_text = self._fetch_europe_pmc_full_text(pubmed_id)
            if full_text is not None:
                metadata = dataclasses.replace(metadata, full_text=full_text)

        return metadata

    def fetch_by_doi(self, doi: str) -> PaperMetadata:
        """Fetch paper metadata by DOI via PubMed search."""
        self._rate_limit()

        search_params = {
            "db": "pubmed",
            "term": f"{doi}[doi]",
            "retmode": "json",
            "retmax": "1",
        }
        if self._api_key:
            search_params["api_key"] = self._api_key

        search_url = f"{self.BASE_URL}/esearch.fcgi?{urllib.parse.urlencode(search_params)}"

        try:
            with self._urlopen_https(search_url) as response:
                search_data = json.loads(response.read().decode("utf-8"))
        except Exception as e:
            raise FetchError(f"Failed to search PubMed for DOI {doi}: {e}")

        ids = search_data.get("esearchresult", {}).get("idlist", [])
        if not ids:
            raise FetchError(
                f"No PubMed record found for DOI {doi}. "
                f"The paper may not be indexed in PubMed, or the DOI may be incorrect."
            )

        return self.fetch_by_pubmed_id(ids[0])

    def _urlopen_https(self, url: str):
        """
        Open a URL, explicitly rejecting anything but https.

        `url` is always built from `self.BASE_URL` (a hardcoded https
        constant) plus urlencoded query params, so this can't actually
        fire today -- it's defense-in-depth against a future refactor
        that makes the URL constructible from less-trusted input, per a
        bandit B310 finding (unrestricted urlopen scheme).
        """
        scheme = urllib.parse.urlparse(url).scheme
        if scheme != "https":
            raise FetchError(f"Refusing to open non-https URL (scheme={scheme!r}): {url}")
        # Scheme validated above; https-only urlopen is intentional here.
        return urllib.request.urlopen(url, timeout=10)  # nosec B310

    def _rate_limit(self) -> None:
        """Enforce rate limiting between API calls."""
        now = time.time()
        elapsed = now - self._last_request_time
        if elapsed < self._rate_limit_delay:
            time.sleep(self._rate_limit_delay - elapsed)
        self._last_request_time = time.time()

    def _parse_pubmed_xml(self, xml_content: str, pubmed_id: str) -> PaperMetadata:
        """Parse PubMed XML response into a PaperMetadata instance."""
        try:
            root = ET.fromstring(xml_content)
        except ET.ParseError as e:
            raise FetchError(f"Failed to parse PubMed XML for {pubmed_id}: {e}")

        article = root.find(".//PubmedArticle")
        if article is None:
            raise FetchError(f"No article found in PubMed response for {pubmed_id}")

        # Title
        title_el = article.find(".//ArticleTitle")
        title = title_el.text or "" if title_el is not None else ""

        # Abstract
        # `el.text` alone only returns text up to the FIRST child element --
        # abstracts with inline markup (<sub>, <i>, common for subscripted
        # units and italicized statistical terms) get silently truncated
        # right at that tag. `itertext()` recursively yields all text
        # including nested elements' content, so join that instead.
        abstract_parts = article.findall(".//AbstractText")
        abstract = " ".join(
            "".join(el.itertext()) for el in abstract_parts
        ).strip()

        # Authors
        authors = []
        for author_el in article.findall(".//Author"):
            last = author_el.findtext("LastName", "")
            initials = author_el.findtext("Initials", "")
            if last:
                authors.append(f"{last} {initials}".strip())

        # Year
        year_el = article.find(".//PubDate/Year")
        year = int(year_el.text) if year_el is not None and year_el.text else None

        # Journal
        journal_el = article.find(".//Journal/Title")
        journal = journal_el.text if journal_el is not None else None

        # DOI
        doi = None
        for id_el in article.findall(".//ArticleId"):
            if id_el.get("IdType") == "doi":
                doi = id_el.text
                break

        return PaperMetadata(
            pubmed_id=pubmed_id,
            doi=doi,
            title=title,
            abstract=abstract,
            authors=authors,
            year=year,
            journal=journal,
            is_primary_source=True,
        )

    def _fetch_europe_pmc_full_text(self, pubmed_id: str) -> Optional[str]:
        """
        Best-effort: return plain-text full text, or None.

        Real calls, verified directly against the live APIs (not assumed
        from documentation):
        1. Europe PMC search?query=EXT_ID:{pmid} AND SRC:MED -- returns
           inEPMC/inPMC/pmcid for this PMID, if PubMed/MEDLINE indexed it.
        2. Europe PMC {pmcid}/fullTextXML -- the JATS XML body. Only
           served for the open-access subset: most papers are not (e.g.
           the metformin benchmark's own PMID:39049802 has none).
        3. FALLBACK, added after a code review found the gap: NCBI's own
           PMC efetch (db=pmc&id={pmcid}). Some papers -- notably NIH
           author manuscripts, e.g. PMID:26197707 / PMC4575248, which
           Europe PMC flags inEPMC="Y", inPMC="Y", isOpenAccess="N" --
           get an HTTP 404/500 from step 2 yet are freely served in full
           by NCBI. Without this, a verification pass wrongly reports
           "no obtainable full text" for a paper whose full text is public.
           Only attempted when step 1 yielded a PMCID.

        Any failure (network, not found, no OA copy, malformed XML)
        returns None rather than raising -- full text is an enrichment,
        never required for fetch_by_pubmed_id()'s core (abstract) result.
        """
        try:
            self._rate_limit()
            search_params = {
                "query": f"EXT_ID:{pubmed_id} AND SRC:MED",
                "format": "json",
            }
            search_url = f"{self.EUROPE_PMC_URL}/search?{urllib.parse.urlencode(search_params)}"
            with self._urlopen_https(search_url) as response:
                search_data = json.loads(response.read().decode("utf-8"))

            results = search_data.get("resultList", {}).get("result", [])
            if not results:
                return None
            record = results[0]
            pmcid = record.get("pmcid")
            if not pmcid or "Y" not in (record.get("inEPMC"), record.get("inPMC")):
                return None
        except Exception:
            return None

        europe_pmc_url = f"{self.EUROPE_PMC_URL}/{pmcid}/fullTextXML"
        pmc_params = {"db": "pmc", "id": pmcid, "retmode": "xml"}
        if self._api_key:
            pmc_params["api_key"] = self._api_key
        ncbi_pmc_url = f"{self.BASE_URL}/efetch.fcgi?{urllib.parse.urlencode(pmc_params)}"
        for url in (europe_pmc_url, ncbi_pmc_url):
            text = self._fetch_jats_body_text(url)
            if text is not None:
                return text
        return None

    def _fetch_jats_body_text(self, url: str) -> Optional[str]:
        """Fetch a JATS XML article from `url` and return its <body> as
        plain text, or None on any failure (network/HTTP/parse/no body)."""
        try:
            self._rate_limit()
            with self._urlopen_https(url) as response:
                xml_content = response.read().decode("utf-8")
            root = ET.fromstring(xml_content)
            body = root.find(".//body")
            if body is None:
                return None
            text = " ".join(body.itertext()).strip()
            return text or None
        except Exception:
            return None
