# SAVE TO: packages\genotype\src\biosim_genotype\parser.py
"""
Variant input parsers.

Supported formats for v1:
    1. rsID + genotype list  — programmatic input
    2. 23andMe raw text      — tab-separated download from 23andMe

Full VCF parsing is intentionally deferred to Sprint 8.
A full genome VCF contains millions of variants; processing only
the relevant rsIDs (those in the EffectDB) is safer and faster.

PRIVACY:
    Parsers return Variant objects in memory only.
    The caller is responsible for not persisting these objects.
    Input text is not logged or stored by this module.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from .variant import Genotype, Variant


@dataclass
class ParseResult:
    """
    Result of parsing a genetic input.

    variants:         successfully parsed Variant objects
    skipped:          lines that were skipped (comments, no-calls, unknown format)
    errors:           lines that failed to parse with error messages
    source_format:    detected input format
    n_total_lines:    total input lines processed
    """
    variants:       list[Variant] = field(default_factory=list)
    skipped:        list[str]     = field(default_factory=list)
    errors:         list[str]     = field(default_factory=list)
    source_format:  str = "unknown"
    n_total_lines:  int = 0

    @property
    def n_variants(self) -> int:
        return len(self.variants)

    @property
    def success_rate(self) -> float:
        if self.n_total_lines == 0:
            return 0.0
        return self.n_variants / self.n_total_lines

    def rsids(self) -> list[str]:
        return [v.rsid for v in self.variants]


class VariantParser:
    """
    Parses genetic variant input into Variant objects.

    For privacy, this class holds no state between calls.
    Each parse call is independent and returns only in-memory objects.
    """

    def parse_rsid_list(
        self,
        rsid_genotype_pairs: list[tuple[str, str]],
    ) -> ParseResult:
        """
        Parse a list of (rsID, genotype_string) pairs.

        Args:
            rsid_genotype_pairs: e.g. [("rs7903146", "0/1"), ("rs1801282", "1/1")]
                                 Genotype strings: VCF-style (0/0, 0/1, 1/1)
                                 or two-allele codes (TT, TC, CC)

        Returns:
            ParseResult with parsed Variant objects
        """
        result = ParseResult(source_format="rsid_list", n_total_lines=len(rsid_genotype_pairs))
        for rsid, geno_str in rsid_genotype_pairs:
            rsid = rsid.strip()
            geno_str = geno_str.strip()
            if not rsid:
                result.skipped.append(f"Empty rsID")
                continue
            if geno_str in ("--", "NA", "NC", ""):
                result.skipped.append(f"{rsid}: no-call")
                continue
            try:
                genotype = Genotype.from_string(geno_str)
                # Two-letter allele codes (e.g. "CC", "CT") carry real
                # allele identity that VCF-style strings ("0/1") don't --
                # preserve it so dosage_for_allele() can count risk-allele
                # copies by direct letter match instead of guessing. See
                # Variant.dosage_for_allele()'s docstring for the
                # regression this fixes (a homozygous non-risk call like
                # "CC" was previously reported as 2 risk-allele copies,
                # the opposite of the true 0).
                observed_alleles = geno_str if len(geno_str) == 2 else None
                result.variants.append(Variant(
                    rsid=rsid, genotype=genotype, observed_alleles=observed_alleles,
                ))
            except ValueError as e:
                result.errors.append(f"{rsid}: {e}")
        return result

    def parse_23andme(self, text: str) -> ParseResult:
        """
        Parse 23andMe raw data format.

        Format (tab-separated):
            # comment lines starting with #
            rsid    chromosome    position    genotype
            rs7903146    10    112998590    CT
            rs1801282     3     12393125    CC

        Two-character genotype codes: 'CT' = heterozygous C/T
        '--' = no-call (skipped)

        Args:
            text: Full text content of the 23andMe raw data file

        Returns:
            ParseResult
        """
        result = ParseResult(source_format="23andme")
        lines = text.strip().split("\n")
        data_lines = [l for l in lines if not l.startswith("#") and l.strip()]
        result.n_total_lines = len(data_lines)

        for line in data_lines:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                result.errors.append(f"Malformed line (expected 4 cols): {line[:60]}")
                continue

            rsid = parts[0].strip()
            chrom = parts[1].strip()
            try:
                pos = int(parts[2].strip())
            except ValueError:
                pos = None
            geno_str = parts[3].strip()

            if not rsid.startswith("rs"):
                result.skipped.append(f"Non-rsID marker: {rsid}")
                continue
            if geno_str == "--":
                result.skipped.append(f"{rsid}: no-call")
                continue

            try:
                genotype = self._parse_two_allele_genotype(geno_str, rsid)
                if genotype is None:
                    result.skipped.append(f"{rsid}: ambiguous genotype '{geno_str}'")
                    continue
                result.variants.append(Variant(
                    rsid=rsid,
                    genotype=genotype,
                    chromosome=chrom,
                    position=pos,
                    # Preserve the raw allele letters -- Genotype only
                    # encodes homozygous/heterozygous SHAPE (ANY homozygous
                    # pair maps to HOMOZYGOUS_ALT), so without this,
                    # dosage_for_allele() cannot tell "CC" (0 copies of risk
                    # allele T) from "TT" (2 copies) -- see Variant
                    # .dosage_for_allele()'s docstring for the regression
                    # this fixes.
                    observed_alleles=geno_str,
                ))
            except ValueError as e:
                result.errors.append(f"{rsid}: {e}")

        return result

    def filter_to_known(
        self,
        result: ParseResult,
        known_rsids: set[str],
    ) -> ParseResult:
        """
        Return a new ParseResult containing only variants with known rsIDs.

        Use this to reduce a large input set to only variants
        that are in the EffectDB, minimising memory footprint.
        """
        filtered = ParseResult(
            source_format=result.source_format,
            n_total_lines=result.n_total_lines,
            skipped=list(result.skipped),
            errors=list(result.errors),
        )
        filtered.variants = [v for v in result.variants if v.rsid in known_rsids]
        n_irrelevant = result.n_variants - filtered.n_variants
        if n_irrelevant > 0:
            filtered.skipped.append(
                f"{n_irrelevant} variants not in effect database (ignored)"
            )
        return filtered

    def _parse_two_allele_genotype(
        self,
        geno_str: str,
        rsid: str,
    ) -> Optional[Genotype]:
        """
        Convert a two-character allele string (e.g. 'CT') to Genotype.

        Since we don't always know which allele is ref/alt in 23andMe
        format, we use homozygous/heterozygous without allele direction.
        The EffectDB handles risk allele orientation.
        """
        if len(geno_str) != 2:
            raise ValueError(f"Expected 2-char genotype, got '{geno_str}'")
        a1, a2 = geno_str[0], geno_str[1]
        if a1 == a2:
            return Genotype.HOMOZYGOUS_ALT  # Both same allele — homozygous
        else:
            return Genotype.HETEROZYGOUS
