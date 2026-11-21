# SAVE TO: packages\genotype\src\biosim_genotype\variant.py
"""
Core variant types.

PRIVACY NOTICE:
    Variant data is processed entirely in memory.
    Nothing in this module writes to disk, transmits data,
    or includes raw variant data in any log output.
    Callers must not persist Variant objects to databases
    without explicit user consent and encryption at rest.

A Variant represents a single SNP with its observed genotype.
Dosage (0, 1, 2) counts the number of non-reference alleles.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Optional


class Genotype(str, Enum):
    """
    Diploid genotype at a biallelic SNP.
    Encoded as VCF-style allele counts relative to reference.
    """
    HOMOZYGOUS_REF = "0/0"   # 0 copies of the alt allele
    HETEROZYGOUS   = "0/1"   # 1 copy  of the alt allele
    HOMOZYGOUS_ALT = "1/1"   # 2 copies of the alt allele

    @classmethod
    def from_string(cls, s: str) -> "Genotype":
        """
        Parse common genotype string formats.

        Accepts: '0/0', '0/1', '1/1', '0|0', '0|1', '1|1',
                 'AA', 'AB', 'BB', '--' (no-call)
        """
        s = s.strip().replace("|", "/")
        mapping = {
            "0/0": cls.HOMOZYGOUS_REF,
            "0/1": cls.HETEROZYGOUS,
            "1/0": cls.HETEROZYGOUS,
            "1/1": cls.HOMOZYGOUS_ALT,
        }
        if s in mapping:
            return mapping[s]
        # Two-character allele codes (23andMe format: 'AA', 'AC', 'CC', '--')
        if len(s) == 2 and s != "--":
            if s[0] == s[1]:
                # Homozygous — treat first char as the genotype indicator
                return cls.HOMOZYGOUS_ALT  # caller must interpret rel. to risk allele
            else:
                return cls.HETEROZYGOUS
        raise ValueError(
            f"Cannot parse genotype '{s}'. "
            f"Expected VCF-style (0/0, 0/1, 1/1) or two-allele code (AA, AC, CC)."
        )


_DOSAGE_MAP = {
    Genotype.HOMOZYGOUS_REF: 0,
    Genotype.HETEROZYGOUS:   1,
    Genotype.HOMOZYGOUS_ALT: 2,
}


@dataclass(frozen=True)
class Variant:
    """
    A single SNP observation for one individual.

    rsid:       dbSNP identifier, e.g. 'rs7903146'
    genotype:   observed diploid genotype
    chromosome: optional, e.g. '10'
    position:   optional, GRCh38 position
    ref_allele: optional, reference allele
    alt_allele: optional, alternate allele

    The dosage property returns the number of alt alleles (0, 1, or 2).
    Whether alt = risk allele depends on the EffectRecord in the database.
    """
    rsid:       str
    genotype:   Genotype
    chromosome: Optional[str] = None
    position:   Optional[int] = None
    ref_allele: Optional[str] = None
    alt_allele: Optional[str] = None
    # Raw 2-letter genotype call (e.g. "CC", "CT"), when the input format
    # provides actual allele letters but not ref/alt orientation (23andMe-
    # style raw data, or a two-letter code passed to parse_rsid_list()).
    # dosage_for_allele() uses this to count risk-allele copies by direct
    # letter match -- see that method's docstring for why this is needed
    # in addition to ref_allele/alt_allele.
    observed_alleles: Optional[str] = None

    def __post_init__(self) -> None:
        if not self.rsid.startswith("rs"):
            raise ValueError(
                f"rsID must start with 'rs', got '{self.rsid}'. "
                f"Only dbSNP rsIDs are supported in v1."
            )
        if self.observed_alleles is not None and len(self.observed_alleles) != 2:
            raise ValueError(
                f"observed_alleles must be exactly 2 characters, got "
                f"'{self.observed_alleles}'"
            )

    @property
    def dosage(self) -> int:
        """Number of alt alleles: 0, 1, or 2."""
        return _DOSAGE_MAP[self.genotype]

    def dosage_for_allele(self, risk_allele: str) -> int:
        """
        Count copies of a specific allele.

        Regression: for two-letter-code input (23andMe raw data, e.g. 'CC',
        'CT'), the genotype only encodes homozygous/heterozygous SHAPE
        (Genotype.from_string() maps ANY homozygous pair to HOMOZYGOUS_ALT
        regardless of which letter), with no ref_allele/alt_allele ever
        attached -- the old fallback below ("allele info not available,
        assume alt=risk") then reported dosage=2 for EVERY homozygous call,
        even when the observed letters were the actual PROTECTIVE
        homozygote (0 real risk-allele copies). Confirmed: a real-world
        'CC' call for rs7903146 (TCF7L2, real risk allele T per
        effect_db.py) reported dosage_for_allele('T')=2, the exact
        opposite of the true 0 -- silently misfires for roughly half of
        all real homozygous genotype calls (whichever letter happens NOT
        to be the risk allele), the primary real-world input format this
        parser exists to support.

        Priority order:
        1. observed_alleles (raw letters, e.g. "CC"/"CT") -- exact, counts
           by direct letter match against risk_allele. Correctly returns 0
           if risk_allele isn't among the observed letters at all.
        2. ref_allele/alt_allele -- when the caller has real VCF-style
           ref/alt orientation.
        3. Fallback (no allele info of any kind) -- standard dosage,
           assuming alt=risk. Only reachable for VCF-style '0/1'/'1/1'
           input with no observed_alleles/ref/alt set at all.
        """
        if self.observed_alleles is not None:
            return sum(1 for allele in self.observed_alleles if allele == risk_allele)
        if self.alt_allele and self.ref_allele:
            if risk_allele == self.alt_allele:
                return self.dosage
            elif risk_allele == self.ref_allele:
                return 2 - self.dosage
        # No allele info available at all — assume alt = risk (conservative)
        return self.dosage
