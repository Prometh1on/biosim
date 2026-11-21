# SAVE TO: packages\genotype\src\biosim_genotype\effect_db.py
"""
Variant effect database.

Maps rsIDs to their effects on pathway model parameters.
Seeded with well-replicated GWAS variants for the glucose/insulin axis.

EPISTEMIC WARNING:
    These genetic variants collectively explain a small fraction of
    inter-individual variance in glycaemic traits (typically <5% per variant).
    Simulation outputs should NOT be interpreted as clinical predictions.
    Effect sizes are population-specific and may not generalise across
    ancestries. See population field and warnings in PersonalizationResult.

Effect records use two composition types:
    MULTIPLICATIVE:    new_value = base * (1 + effect_per_allele * dosage)
    ADDITIVE_ABSOLUTE: new_value = base + effect_per_allele * dosage

Literature sources (verified against PubMed E-utilities -- see
KNOWN BUGS in CLAUDE.md for the citation-accuracy audit that found and
fixed the PMIDs below; several were previously fabricated/mismatched):
    Grant et al. 2006 (TCF7L2)      PMID:16415884
    Altshuler et al. 2000 (PPARG)   PMID:10973253
    Gloyn et al. 2003 (KCNJ11)      PMID:12540637
    Scott et al. 2007 (SLC30A8)     PMID:17463248
    Prokopenko et al. 2009 (MTNR1B) PMID:19060907
    Dupuis et al. 2010 (GCK)        PMID:20081858
    IGF1 rs35767: no specific single-SNP discovery paper could be
    confirmed against PubMed -- pubmed_id left unset below rather than
    asserting an unverified citation. See effect_per_allele docstring.

`population_frequency` source (all 7 default variants, added when
`biosim_seeds.population` gained genetic-variant sampling): TOPMed
(NHLBI Trans-Omics for Precision Medicine, freeze 8, n=264,690),
fetched live via NCBI's dbSNP REST API
(api.ncbi.nlm.nih.gov/variation/v0/refsnp/{rsid}) and read directly
from each response's `primary_snapshot_data.allele_annotations[*]
.frequency[*]` entries -- not a secondary aggregator's paraphrase. Used
consistently across all 7 variants (rather than mixing each GWAS
discovery paper's own study-specific case/control allele frequency,
which varies in ancestry composition and is not meant to represent
general-population prevalence) because a single, large, modern,
ancestry-diverse reference cohort is the more defensible choice for a
"general population" synthetic cohort than pooling frequencies from
studies never designed for that purpose. Extended (non-default) T2D
variants below do NOT yet have verified frequencies -- left at the
default None (not sampled) rather than guessed, a smaller, deliberately
bounded first increment; see biosim_seeds.population's own genetic-
sampling code for how it handles a record with no known frequency.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

# User-facing disclosure for variants whose effect record has
# pubmed_id="UNVERIFIED", keyed by rsID. The text itself must never contain
# the rsID or allele (CLAUDE.md rule 7). PMIDs re-fetched live 2026-11-14.
UNVERIFIED_VARIANT_NOTES: dict[str, str] = {
    "rs35767": (
        "A genetic variant in this profile changes baseline IGF-1 by an "
        "UNVERIFIED effect size: no published study reports it, and published "
        "evidence on its direction conflicts (PMID:30115944 reports the variant "
        "raises IGF-1; PMID:22689104 and PMID:21046444 found no association "
        "with circulating IGF-1). Treat the IGF-1 result as illustrative only."
    ),
}


class Population(str, Enum):
    EUR = "EUR"   # European
    EAS = "EAS"   # East Asian
    AFR = "AFR"   # African
    AMR = "AMR"   # Admixed American
    SAS = "SAS"   # South Asian
    ALL = "ALL"   # Pan-ancestry / population-agnostic


class EffectType(str, Enum):
    MULTIPLICATIVE    = "multiplicative"     # new = base * (1 + effect * dosage)
    ADDITIVE_ABSOLUTE = "additive_absolute"  # new = base + effect * dosage


@dataclass(frozen=True)
class EffectRecord:
    """
    Effect of one variant on one pathway parameter.

    parameter_field:    exact field name in PathwayParameters to modify
    effect_per_allele:  signed effect per risk allele copy
                        positive = increases parameter
                        negative = decreases parameter
    risk_allele:        which allele is the risk/effect allele
    effect_se:          uncertainty of effect_per_allele (project estimate --
                        GWAS report odds ratios, not this parameter's SE)
    confidence_tier:    'gold' (replicated >10 studies), 'silver' (3-10), 'bronze' (<3)
    population_frequency: frequency of risk_allele in the population named
                        by `population` (0.0-1.0). None if not yet
                        verified against a real source -- never a guessed
                        placeholder (see biosim_seeds.population's own
                        genetic-sampling code, the only current consumer:
                        it skips any record with population_frequency=None
                        rather than assuming a value).
    """
    rsid:             str
    gene:             str
    risk_allele:      str
    pathway_name:     str
    parameter_field:  str
    effect_per_allele: float
    effect_type:      EffectType
    population:       Population
    confidence_tier:  str
    effect_se:        float
    pubmed_id:        str
    description:      str
    population_frequency: Optional[float] = None

    def __post_init__(self) -> None:
        if self.effect_se < 0:
            raise ValueError("effect_se must be non-negative")
        if self.confidence_tier not in ("gold", "silver", "bronze"):
            raise ValueError(f"confidence_tier must be gold/silver/bronze, got '{self.confidence_tier}'")
        if self.population_frequency is not None and not 0.0 <= self.population_frequency <= 1.0:
            raise ValueError(
                f"population_frequency must be in [0, 1] or None, got "
                f"{self.population_frequency}"
            )

    def apply(self, base_value: float, dosage: int) -> float:
        """
        Apply this effect to a base parameter value.

        Args:
            base_value: The parameter's current value
            dosage: Number of risk alleles (0, 1, or 2)

        Returns:
            Modified parameter value
        """
        if dosage == 0:
            return base_value
        if self.effect_type == EffectType.MULTIPLICATIVE:
            return base_value * (1.0 + self.effect_per_allele * dosage)
        else:  # ADDITIVE_ABSOLUTE
            return base_value + self.effect_per_allele * dosage

    @property
    def uncertainty_contribution(self) -> float:
        """
        Fractional uncertainty contribution (CV equivalent).
        Used to inflate parameter uncertainty in the UQ layer.

        Regression: this used to also multiply by abs(effect_per_allele)
        after dividing by it -- those cancel exactly for any nonzero
        effect_per_allele, so this always just returned effect_se
        unchanged, regardless of effect_type. effect_se is documented as
        "standard error of effect_per_allele" -- the same units as
        effect_per_allele, which is a genuine fraction for MULTIPLICATIVE
        records but an absolute value (mM, nM) for ADDITIVE_ABSOLUTE
        records (Gb, G_threshold, IGF1b). Composed downstream as a
        dimensionless CV combined in quadrature with the tier-based CV
        (composer.py), so an absolute-unit value there was a real
        unit-consistency bug -- just low-impact today because the
        absolute effect_se values are small relative to typical
        silver-tier CVs (~0.30). The division alone (SE / effect =
        relative uncertainty of the estimate) is dimensionally correct
        for both effect types.
        """
        if abs(self.effect_per_allele) == 0:
            return 0.0
        return abs(self.effect_se / self.effect_per_allele)


class VariantEffectDB:
    """
    Database of variant-to-parameter effect mappings.

    Each rsID can have multiple effect records (different populations,
    different parameters, or different pathway contexts).
    """

    def __init__(self) -> None:
        self._records: dict[str, list[EffectRecord]] = {}

    def add(self, record: EffectRecord) -> None:
        self._records.setdefault(record.rsid, []).append(record)

    def get(self, rsid: str) -> list[EffectRecord]:
        """Return all effect records for an rsID. Empty list if not found."""
        return list(self._records.get(rsid, []))

    def has(self, rsid: str) -> bool:
        return rsid in self._records

    def all_rsids(self) -> list[str]:
        return sorted(self._records.keys())

    def count(self) -> int:
        return sum(len(v) for v in self._records.values())

    def filter_by_population(self, population: Population) -> "VariantEffectDB":
        """Return a new DB containing only records for a specific population."""
        filtered = VariantEffectDB()
        for records in self._records.values():
            for r in records:
                if r.population in (population, Population.ALL):
                    filtered.add(r)
        return filtered

    @classmethod
    def default(cls) -> "VariantEffectDB":
        """
        Return the default database seeded with glucose/insulin axis variants.
        Variants were identified in published GWAS; per-allele effect sizes
        are the project's literature-informed parameter mappings, not values
        reported by GWAS (which give disease odds ratios).
        """
        db = cls()

        # ------------------------------------------------------------------
        # TCF7L2 rs7903146 — strongest common T2D locus
        # T allele reduces Wnt signalling in beta cells → less insulin secretion
        # Effect: ~12% reduction in secretion per T allele
        # Source: Grant et al. 2006, Nature Genetics. PMID:16415884
        # ------------------------------------------------------------------
        db.add(EffectRecord(
            rsid="rs7903146",
            gene="TCF7L2",
            risk_allele="T",
            pathway_name="glucose_insulin_v1",
            parameter_field="secretion_scale",
            effect_per_allele=-0.12,
            effect_type=EffectType.MULTIPLICATIVE,
            population=Population.EUR,
            confidence_tier="gold",
            effect_se=0.015,
            pubmed_id="16415884",
            description="TCF7L2 T2D risk: reduces beta-cell insulin secretion by ~12% per T allele",
            population_frequency=0.2765,
        ))

        # ------------------------------------------------------------------
        # PPARG rs1801282 — Pro12Ala protective variant
        # G allele (Ala) → reduced PPARG activity → increased insulin sensitivity
        # Effect: ~15% increase in insulin action per G allele
        # Source: Altshuler et al. 2000, Nature Genetics. PMID:10973253
        # ------------------------------------------------------------------
        db.add(EffectRecord(
            rsid="rs1801282",
            gene="PPARG",
            risk_allele="G",
            pathway_name="glucose_insulin_v1",
            parameter_field="p3",
            effect_per_allele=0.15,
            effect_type=EffectType.MULTIPLICATIVE,
            population=Population.EUR,
            confidence_tier="gold",
            effect_se=0.025,
            pubmed_id="10973253",
            description="PPARG Pro12Ala: G allele increases insulin sensitivity ~15% per allele (protective)",
            population_frequency=0.0805,
        ))

        # ------------------------------------------------------------------
        # KCNJ11 rs5219 — Glu23Lys in ATP-sensitive K+ channel
        # T allele impairs channel closure → reduced insulin secretion
        # Effect: ~8% reduction in secretion per T allele
        # Source: Gloyn et al. 2003, Diabetes. PMID:12540637
        # ------------------------------------------------------------------
        db.add(EffectRecord(
            rsid="rs5219",
            gene="KCNJ11",
            risk_allele="T",
            pathway_name="glucose_insulin_v1",
            parameter_field="secretion_scale",
            effect_per_allele=-0.08,
            effect_type=EffectType.MULTIPLICATIVE,
            population=Population.EUR,
            confidence_tier="gold",
            effect_se=0.012,
            pubmed_id="12540637",
            description="KCNJ11 Glu23Lys: T allele reduces beta-cell insulin secretion ~8% per allele",
            population_frequency=0.2703,
        ))

        # ------------------------------------------------------------------
        # SLC30A8 rs13266634 — Arg325Trp in zinc transporter
        # T allele impairs zinc transport → less insulin crystallisation
        # Effect: ~6% reduction in secretion per T allele
        # Source: Scott et al. 2007, Science. PMID:17463248
        # ------------------------------------------------------------------
        db.add(EffectRecord(
            rsid="rs13266634",
            gene="SLC30A8",
            risk_allele="T",
            pathway_name="glucose_insulin_v1",
            parameter_field="secretion_scale",
            effect_per_allele=-0.06,
            effect_type=EffectType.MULTIPLICATIVE,
            population=Population.EUR,
            confidence_tier="silver",
            effect_se=0.014,
            pubmed_id="17463248",
            description="SLC30A8 Arg325Trp: T allele reduces beta-cell insulin secretion ~6% per allele",
            population_frequency=0.2342,
        ))

        # ------------------------------------------------------------------
        # MTNR1B rs10830963 — melatonin receptor 1B
        # G allele increases fasting glucose via circadian mechanisms
        # Effect: +0.07 mM per G allele (~1.3 mg/dL)
        # Source: Prokopenko et al. 2009, Nature Genetics. PMID:19060907
        # ------------------------------------------------------------------
        db.add(EffectRecord(
            rsid="rs10830963",
            gene="MTNR1B",
            risk_allele="G",
            pathway_name="glucose_insulin_v1",
            parameter_field="Gb",
            effect_per_allele=0.07,
            effect_type=EffectType.ADDITIVE_ABSOLUTE,
            population=Population.EUR,
            confidence_tier="gold",
            effect_se=0.008,
            pubmed_id="19060907",
            description="MTNR1B: G allele increases fasting glucose by ~0.07 mM per allele",
            population_frequency=0.2114,
        ))

        # ------------------------------------------------------------------
        # GCK rs4607517 — near glucokinase gene
        # A allele lowers glucose threshold for beta-cell activation
        # Effect: -0.09 mM shift in glucose sensing threshold per A allele
        # Source: Dupuis et al. 2010, Nature Genetics. PMID:20081858
        # ------------------------------------------------------------------
        db.add(EffectRecord(
            rsid="rs4607517",
            gene="GCK",
            risk_allele="A",
            pathway_name="glucose_insulin_v1",
            parameter_field="G_threshold",
            effect_per_allele=-0.09,
            effect_type=EffectType.ADDITIVE_ABSOLUTE,
            population=Population.EUR,
            confidence_tier="gold",
            effect_se=0.010,
            pubmed_id="20081858",
            description="GCK locus: A allele lowers glucose sensing threshold ~0.09 mM per allele",
            population_frequency=0.1514,
        ))

        # ------------------------------------------------------------------
        # IGF1 rs35767 — IGF1 gene promoter
        # Minor allele associated with higher circulating IGF-1
        # Effect: +0.5 nM per minor allele
        # Source: unverified. The previously-cited PMID:9806922 was checked
        # against PubMed E-utilities and found to be an unrelated paper
        # (Hoxa1/Krox-20 rhombomere development). No specific rs35767
        # discovery paper could be confirmed to replace it -- rather than
        # assert another unverified number, pubmed_id is left as a clear
        # placeholder and confidence downgraded to bronze until a real
        # citation is located (see CLAUDE.md known-bugs / citation audit).
        #
        # UNVERIFIED, re-checked live on PubMed 2026-11-12: no paper reports
        # an IGF-1 change per allele, so the +0.5 nM size is a placeholder.
        # The direction is contested -- Mannino 2018 (PMID:30115944) calls
        # the T allele "IGF-1-raising", while Tekcan 2012 (PMID:22689104)
        # and Burchardt 2010 (PMID:21046444) found no association with
        # circulating IGF-1.
        # ------------------------------------------------------------------
        db.add(EffectRecord(
            rsid="rs35767",
            gene="IGF1",
            risk_allele="A",
            pathway_name="glucose_insulin_v1",
            parameter_field="IGF1b",
            effect_per_allele=0.5,
            effect_type=EffectType.ADDITIVE_ABSOLUTE,
            population=Population.EUR,
            confidence_tier="bronze",
            effect_se=0.12,
            pubmed_id="UNVERIFIED",
            description=(
                "IGF1 promoter: A allele assumed to raise baseline IGF-1 "
                "(~0.5 nM per allele) -- UNVERIFIED placeholder, no source "
                "reports this effect size"
            ),
            # Allele frequency is answerable independently of the
            # unresolved discovery-paper/effect-size question above --
            # dbSNP/TOPMed reports real population allele counts for
            # this rsID regardless of which paper (if any) discovered
            # its phenotype association. Verified live (see
            # population_frequency docstring block below for the
            # source and method shared by all 7 default variants).
            population_frequency=0.2672,
        ))

        return db


def extended_t2d_variants() -> "VariantEffectDB":
    """
    Extended T2D variant database adding DIAGRAM consortium loci.

    Sources:
        Morris et al. 2012, Nature Genetics. PMID:22885922
        Voight et al. 2010, Nature Genetics. PMID:20581827
        Zeggini et al. 2008, Nature Genetics. PMID:18372903
    """
    db = VariantEffectDB.default()   # includes the original 7

    # KCNQ1 rs2237892 — voltage-gated K+ channel in beta cells
    # Strongest T2D locus in East Asian populations
    # Source: Unoki et al. 2008, Nature Genetics. PMID:18711366
    db.add(EffectRecord(
        rsid="rs2237892", gene="KCNQ1", risk_allele="T",
        pathway_name="glucose_insulin_v1",
        parameter_field="secretion_scale",
        effect_per_allele=-0.10, effect_type=EffectType.MULTIPLICATIVE,
        population=Population.EAS, confidence_tier="gold",
        effect_se=0.012, pubmed_id="18711366",
        description="KCNQ1 T2D risk (EAS): T allele reduces insulin secretion ~10% per allele",
    ))

    # CDKN2A/B rs10811661 — cell cycle regulation in beta cells
    db.add(EffectRecord(
        rsid="rs10811661", gene="CDKN2A", risk_allele="T",
        pathway_name="glucose_insulin_v1",
        parameter_field="secretion_scale",
        effect_per_allele=-0.07, effect_type=EffectType.MULTIPLICATIVE,
        population=Population.EUR, confidence_tier="gold",
        effect_se=0.011, pubmed_id="17463249",
        description="CDKN2A/B locus: T allele reduces beta-cell function ~7% per allele",
    ))

    # CDKAL1 rs7756992 — tRNA methylthiotransferase in beta cells
    db.add(EffectRecord(
        rsid="rs7756992", gene="CDKAL1", risk_allele="G",
        pathway_name="glucose_insulin_v1",
        parameter_field="secretion_scale",
        effect_per_allele=-0.08, effect_type=EffectType.MULTIPLICATIVE,
        population=Population.EUR, confidence_tier="gold",
        effect_se=0.010, pubmed_id="17460697",
        description="CDKAL1 locus: G allele reduces insulin secretion ~8% per allele",
    ))

    # IGF2BP2 rs4402960 — IGF2 mRNA binding protein
    # Also affects IGF-1 signalling pathway
    db.add(EffectRecord(
        rsid="rs4402960", gene="IGF2BP2", risk_allele="T",
        pathway_name="glucose_insulin_v1",
        parameter_field="secretion_scale",
        effect_per_allele=-0.06, effect_type=EffectType.MULTIPLICATIVE,
        population=Population.EUR, confidence_tier="gold",
        effect_se=0.009, pubmed_id="17463248",
        description="IGF2BP2 locus: T allele reduces insulin secretion ~6% per allele",
    ))

    # IGF2BP2 — secondary effect on IGF-1 production
    db.add(EffectRecord(
        rsid="rs4402960", gene="IGF2BP2", risk_allele="T",
        pathway_name="glucose_insulin_v1",
        parameter_field="IGF1b",
        effect_per_allele=-0.4, effect_type=EffectType.ADDITIVE_ABSOLUTE,
        population=Population.EUR, confidence_tier="bronze",
        effect_se=0.20, pubmed_id="17463248",
        description="IGF2BP2 locus: T allele associated with lower circulating IGF-1 (~0.4 nM per allele)",
    ))

    # ADCY5 rs11708067 — adenylate cyclase 5 in beta cells
    db.add(EffectRecord(
        rsid="rs11708067", gene="ADCY5", risk_allele="A",
        pathway_name="glucose_insulin_v1",
        parameter_field="secretion_scale",
        effect_per_allele=-0.09, effect_type=EffectType.MULTIPLICATIVE,
        population=Population.EUR, confidence_tier="silver",
        effect_se=0.016, pubmed_id="22885922",
        description="ADCY5 locus: A allele reduces insulin secretion ~9% per allele via cAMP pathway",
    ))

    # HHEX rs1111875 — homeobox transcription factor in beta cells
    db.add(EffectRecord(
        rsid="rs1111875", gene="HHEX", risk_allele="C",
        pathway_name="glucose_insulin_v1",
        parameter_field="secretion_scale",
        effect_per_allele=-0.05, effect_type=EffectType.MULTIPLICATIVE,
        population=Population.EUR, confidence_tier="gold",
        effect_se=0.008, pubmed_id="17554300",
        description="HHEX locus: C allele reduces beta-cell function ~5% per allele",
    ))

    # WFS1 rs10010131 — Wolfram syndrome protein (protective!)
    # A allele is PROTECTIVE — reduces T2D risk by improving beta-cell function
    # Source: Sandhu et al. 2007, Nature Genetics. PMID:17603484
    db.add(EffectRecord(
        rsid="rs10010131", gene="WFS1", risk_allele="A",
        pathway_name="glucose_insulin_v1",
        parameter_field="secretion_scale",
        effect_per_allele=+0.06, effect_type=EffectType.MULTIPLICATIVE,
        population=Population.EUR, confidence_tier="silver",
        effect_se=0.012, pubmed_id="17603484",
        description="WFS1 locus: A allele is PROTECTIVE — increases beta-cell function ~6% per allele",
    ))

    return db
