# SAVE TO: packages\literature-pipeline\src\biosim_literature\verification.py
"""
Citation re-verification tooling.

This project's own evidence-hierarchy rule requires every PMID cited in
source code to be checked against a REAL fetched paper, not asserted
from memory (see CLAUDE.md's KNOWN BUGS FIXED #9/#13 -- ~35 fabricated/
mismatched PMIDs were found this way). Until now that check was
performed as an ad-hoc one-off script per citation, re-derived each
time (`python -c "... PubMedFetcher().fetch_by_pubmed_id(...) ..."`).
This module makes it a single, tested, reusable call.

**Why full-text matters, not just the abstract**: a claim can be true
of a paper's own Table 1/Results body while being entirely absent from
its abstract -- an abstract-only re-verification would then wrongly
read as "not confirmed" for a citation that is, in fact, correct. This
exact gap blocked verifying the Cushing's-syndrome cortisol-ratio
citation (PMID:31695455, cortisol_pathway.py's own CUSHINGS_CORTISOL_
RATIO): its abstract discusses only mitochondrial/OXPHOS gene
expression, never mentioning cortisol at all, while its Table 1 (only
reachable via full text) directly reports the 527/822 nmol/L figures
the code cites. `verify_citation()` always attempts full text (via
PubMedFetcher's existing Europe PMC integration, `fetch_full_text=
True`) precisely so a claim like this can be confirmed rather than
reported as an apparent mismatch.

Full text is only available for open-access papers indexed in Europe
PMC -- most are not (see PubMedFetcher._fetch_europe_pmc_full_text's
own docstring). `CitationVerification.full_text_available` reports
this explicitly so a "not found" result can be read correctly: either
the claim really isn't in this paper, or full text simply wasn't
available to check beyond the abstract -- these are NOT the same
finding, and the report distinguishes them via `ClaimCheck.found_in`
staying `None` either way, but `full_text_available` telling you which
case you're in.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from .fetcher import Fetcher, PubMedFetcher

# Named so a future reader can tell -- without re-deriving it --
# how much surrounding text a ClaimCheck.context snippet carries.
DEFAULT_CONTEXT_CHARS = 120


def _normalize(text: str) -> str:
    """Collapse whitespace runs to a single space and lowercase.

    Real fetched text (XML-extracted abstracts, JATS full text) varies
    in whitespace/line-wrapping in ways a claim string written by hand
    never matches exactly -- normalizing both sides before comparing
    is what makes substring search reliable here.
    """
    return " ".join(text.split()).lower()


@dataclass(frozen=True)
class ClaimCheck:
    """Whether one claimed quote/value was found in the fetched paper."""
    claim: str
    found: bool
    found_in: Optional[str] = None  # "abstract" | "full_text" | None
    context: Optional[str] = None   # normalized surrounding snippet, for a human to eyeball


@dataclass(frozen=True)
class CitationVerification:
    """Result of re-verifying one PMID's claimed quotes/values."""
    pubmed_id: str
    title: str
    journal: Optional[str]
    year: Optional[int]
    authors: list[str] = field(default_factory=list)
    full_text_available: bool = False
    claim_checks: list[ClaimCheck] = field(default_factory=list)

    @property
    def all_claims_verified(self) -> bool:
        return all(c.found for c in self.claim_checks)

    @property
    def unverified_claims(self) -> list[str]:
        return [c.claim for c in self.claim_checks if not c.found]

    def as_text(self) -> str:
        """Human-readable report, e.g. for a CLI or a code-review comment."""
        lines = [
            f"PMID:{self.pubmed_id} -- {self.title}",
            f"  {', '.join(self.authors[:3])}{' et al.' if len(self.authors) > 3 else ''}"
            f" ({self.year}), {self.journal or 'unknown journal'}",
            f"  full text available: {self.full_text_available}",
        ]
        for check in self.claim_checks:
            status = f"FOUND (in {check.found_in})" if check.found else "NOT FOUND"
            lines.append(f"  [{status}] {check.claim!r}")
            if check.context:
                lines.append(f"      ...{check.context}...")
        return "\n".join(lines)


def verify_citation(
    pubmed_id: str,
    claims: list[str],
    fetcher: Optional[Fetcher] = None,
    context_chars: int = DEFAULT_CONTEXT_CHARS,
) -> CitationVerification:
    """
    Re-fetch `pubmed_id` and check whether each string in `claims` is
    genuinely present in the paper (abstract, then full text if the
    abstract doesn't contain it and full text is available).

    `claims` are loose substrings -- e.g. "822" and "527" for a numeric
    claim, or a short quoted phrase -- matched after whitespace
    normalization and lowercasing (see `_normalize`), not as exact
    verbatim substrings of the raw fetched text.

    `fetcher` defaults to a real `PubMedFetcher(fetch_full_text=True)`
    -- pass a `MockFetcher` in tests to avoid live network calls. If a
    fetcher without full-text support is passed (fetch_full_text=False,
    or a `MockFetcher` paper with no `full_text` set),
    `full_text_available` is correctly reported False rather than
    raising.

    Raises whatever the underlying fetcher raises (e.g. `FetchError`)
    if the PMID itself can't be fetched at all -- an unfetchable
    citation must fail loudly, not silently report every claim as
    "not found".
    """
    if fetcher is None:
        fetcher = PubMedFetcher(fetch_full_text=True)

    metadata = fetcher.fetch_by_pubmed_id(pubmed_id)

    abstract_norm = _normalize(metadata.abstract)
    full_text_norm = _normalize(metadata.full_text) if metadata.full_text is not None else None

    checks: list[ClaimCheck] = []
    for claim in claims:
        claim_norm = _normalize(claim)
        found_in: Optional[str] = None
        source: Optional[str] = None
        idx = -1

        abstract_idx = abstract_norm.find(claim_norm)
        if abstract_idx != -1:
            found_in, source, idx = "abstract", abstract_norm, abstract_idx
        elif full_text_norm is not None:
            full_text_idx = full_text_norm.find(claim_norm)
            if full_text_idx != -1:
                found_in, source, idx = "full_text", full_text_norm, full_text_idx

        context = None
        if source is not None:
            start = max(0, idx - context_chars)
            end = idx + len(claim_norm) + context_chars
            context = source[start:end]

        checks.append(ClaimCheck(
            claim=claim, found=found_in is not None, found_in=found_in, context=context,
        ))

    return CitationVerification(
        pubmed_id=pubmed_id,
        title=metadata.title,
        journal=metadata.journal,
        year=metadata.year,
        authors=list(metadata.authors),
        full_text_available=metadata.full_text is not None,
        claim_checks=checks,
    )


def _main() -> None:
    """
    CLI entry point: `python -m biosim_literature.verification <pmid> <claim> [<claim> ...]`

    Replaces the ad-hoc `python -c "..."` one-off scripts this project
    had been writing per citation with one reusable, real (network-
    calling) command.
    """
    import sys

    # Windows consoles often default stdout to cp1252, which cannot
    # encode many real author names/journal titles (e.g. "Ď", "ø") --
    # confirmed by reproducing a live UnicodeEncodeError against a real
    # PMID. Reconfigure to utf-8 with a safe fallback rather than
    # crashing on legitimate paper metadata.
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    if len(sys.argv) < 3:
        print(
            "Usage: python -m biosim_literature.verification <pmid> <claim> [<claim> ...]",
            file=sys.stderr,
        )
        raise SystemExit(2)

    pmid, claims = sys.argv[1], sys.argv[2:]
    result = verify_citation(pmid, claims)
    print(result.as_text())
    if not result.all_claims_verified:
        raise SystemExit(1)


if __name__ == "__main__":
    _main()
