# SAVE TO: packages\genotype\tests\test_genotype.py
"""
Test suite for biosim_genotype.

Coverage:
    - Variant: rsID validation, dosage, allele-aware dosage
    - Genotype: parsing from VCF and 23andMe string formats
    - VariantEffectDB: seeding, lookup, filtering
    - EffectRecord: effect application, uncertainty contribution
    - VariantParser: rsID list, 23andMe text, no-call handling
    - EffectComposer: single variant, multi-variant, capping, population warnings
    - GenotypePersonalizer: full pipeline, no variants, unknown rsIDs,
                            ODE still runs after personalisation
    - PersonalizationResult: summary, meaningfulness flag, privacy
"""

import math
import pytest

from biosim_genotype import (
    ComposedEffect,
    EffectComposer,
    EffectRecord,
    EffectType,
    GenotypePersonalizer,
    Genotype,
    ParseResult,
    PersonalizationResult,
    Population,
    Variant,
    VariantEffectDB,
    VariantParser,
)
from biosim_ode import GlucoseInsulinPathway, ODESystem, PathwayParameters, high_sugar_meal


# ---------------------------------------------------------------------------
# Variant and Genotype
# ---------------------------------------------------------------------------

class TestVariant:
    def test_valid_rsid(self):
        v = Variant("rs7903146", Genotype.HETEROZYGOUS)
        assert v.rsid == "rs7903146"

    def test_invalid_rsid_raises(self):
        with pytest.raises(ValueError, match="rs"):
            Variant("7903146", Genotype.HETEROZYGOUS)

    def test_dosage_homozygous_ref(self):
        assert Variant("rs1", Genotype.HOMOZYGOUS_REF).dosage == 0

    def test_dosage_heterozygous(self):
        assert Variant("rs1", Genotype.HETEROZYGOUS).dosage == 1

    def test_dosage_homozygous_alt(self):
        assert Variant("rs1", Genotype.HOMOZYGOUS_ALT).dosage == 2

    def test_dosage_for_allele_alt_matches(self):
        v = Variant("rs1", Genotype.HETEROZYGOUS, ref_allele="C", alt_allele="T")
        assert v.dosage_for_allele("T") == 1   # alt = risk: standard dosage

    def test_dosage_for_allele_ref_matches(self):
        v = Variant("rs1", Genotype.HETEROZYGOUS, ref_allele="C", alt_allele="T")
        assert v.dosage_for_allele("C") == 1   # ref = risk: 2 - dosage = 2 - 1 = 1

    def test_dosage_for_allele_no_info(self):
        v = Variant("rs1", Genotype.HETEROZYGOUS)
        assert v.dosage_for_allele("T") == 1   # no allele info: use standard dosage

    def test_invalid_observed_alleles_length_raises(self):
        with pytest.raises(ValueError, match="observed_alleles"):
            Variant("rs1", Genotype.HOMOZYGOUS_ALT, observed_alleles="CCC")

    def test_observed_alleles_homozygous_matches_risk(self):
        # "TT" homozygous, risk_allele="T" -> both copies are the risk allele
        v = Variant("rs1", Genotype.HOMOZYGOUS_ALT, observed_alleles="TT")
        assert v.dosage_for_allele("T") == 2

    def test_observed_alleles_homozygous_does_not_match_risk(self):
        # Regression: Genotype.from_string() maps ANY homozygous 2-letter
        # code to HOMOZYGOUS_ALT regardless of which letter -- without
        # observed_alleles, the old fallback ("no allele info, assume
        # alt=risk") reported dosage=2 for a "CC" call even when the risk
        # allele is "T" (i.e. the person carries ZERO risk alleles, the
        # PROTECTIVE homozygote) -- the exact opposite of correct, and not
        # a rare edge case: it misfired for whichever homozygous letter
        # happened not to be the risk allele, roughly half of all real
        # homozygous calls in the 23andMe input format this exists for.
        v = Variant("rs1", Genotype.HOMOZYGOUS_ALT, observed_alleles="CC")
        assert v.dosage_for_allele("T") == 0

    def test_observed_alleles_heterozygous_matches_either_allele(self):
        v = Variant("rs1", Genotype.HETEROZYGOUS, observed_alleles="CT")
        assert v.dosage_for_allele("T") == 1
        assert v.dosage_for_allele("C") == 1

    def test_observed_alleles_takes_priority_over_ref_alt(self):
        # If both observed_alleles and ref/alt are somehow set, the exact
        # letter match must win -- it's the more precise source of truth.
        v = Variant(
            "rs1", Genotype.HOMOZYGOUS_ALT, observed_alleles="CC",
            ref_allele="C", alt_allele="T",
        )
        assert v.dosage_for_allele("T") == 0


class TestGenotypeParsing:
    def test_vcf_homozygous_ref(self):
        assert Genotype.from_string("0/0") == Genotype.HOMOZYGOUS_REF

    def test_vcf_heterozygous(self):
        assert Genotype.from_string("0/1") == Genotype.HETEROZYGOUS

    def test_vcf_phased(self):
        assert Genotype.from_string("0|1") == Genotype.HETEROZYGOUS

    def test_vcf_homozygous_alt(self):
        assert Genotype.from_string("1/1") == Genotype.HOMOZYGOUS_ALT

    def test_two_allele_heterozygous(self):
        assert Genotype.from_string("CT") == Genotype.HETEROZYGOUS

    def test_two_allele_homozygous(self):
        assert Genotype.from_string("CC") == Genotype.HOMOZYGOUS_ALT

    def test_invalid_raises(self):
        with pytest.raises(ValueError):
            Genotype.from_string("XYZ")


# ---------------------------------------------------------------------------
# EffectRecord
# ---------------------------------------------------------------------------

class TestEffectRecord:
    @pytest.fixture
    def tcf7l2_record(self):
        return EffectRecord(
            rsid="rs7903146", gene="TCF7L2", risk_allele="T",
            pathway_name="glucose_insulin_v1",
            parameter_field="secretion_scale",
            effect_per_allele=-0.12, effect_type=EffectType.MULTIPLICATIVE,
            population=Population.EUR, confidence_tier="gold",
            effect_se=0.015, pubmed_id="16415884",
            description="TCF7L2 T2D risk",
        )

    def test_apply_no_alleles(self, tcf7l2_record):
        assert tcf7l2_record.apply(1.0, 0) == pytest.approx(1.0)

    def test_apply_one_allele_multiplicative(self, tcf7l2_record):
        # 1.0 * (1 + (-0.12) * 1) = 0.88
        assert tcf7l2_record.apply(1.0, 1) == pytest.approx(0.88)

    def test_apply_two_alleles_multiplicative(self, tcf7l2_record):
        # 1.0 * (1 + (-0.12) * 2) = 0.76
        assert tcf7l2_record.apply(1.0, 2) == pytest.approx(0.76)

    def test_apply_additive(self):
        r = EffectRecord(
            rsid="rs10830963", gene="MTNR1B", risk_allele="G",
            pathway_name="glucose_insulin_v1",
            parameter_field="Gb",
            effect_per_allele=0.07, effect_type=EffectType.ADDITIVE_ABSOLUTE,
            population=Population.EUR, confidence_tier="gold",
            effect_se=0.008, pubmed_id="19060910", description="MTNR1B",
        )
        assert r.apply(5.0, 1) == pytest.approx(5.07)
        assert r.apply(5.0, 2) == pytest.approx(5.14)

    def test_uncertainty_contribution_multiplicative_is_relative_se(self, tcf7l2_record):
        # effect_se=0.015, effect_per_allele=-0.12 -> CV = 0.015/0.12 = 0.125
        assert tcf7l2_record.uncertainty_contribution == pytest.approx(0.125)

    def test_uncertainty_contribution_additive_is_dimensionless(self):
        # Regression: this used to also multiply by abs(effect_per_allele)
        # after dividing by it -- those cancel exactly, so it always just
        # returned effect_se unchanged (0.008 mM, an absolute value) instead
        # of a dimensionless CV, regardless of effect_type. For an
        # ADDITIVE_ABSOLUTE record, effect_se is in the SAME absolute units
        # as effect_per_allele (mM here), not already a fraction.
        r = EffectRecord(
            rsid="rs10830963", gene="MTNR1B", risk_allele="G",
            pathway_name="glucose_insulin_v1",
            parameter_field="Gb",
            effect_per_allele=0.07, effect_type=EffectType.ADDITIVE_ABSOLUTE,
            population=Population.EUR, confidence_tier="gold",
            effect_se=0.008, pubmed_id="19060910", description="MTNR1B",
        )
        # CV = 0.008 / 0.07 ≈ 0.114, NOT the bare effect_se (0.008)
        assert r.uncertainty_contribution == pytest.approx(0.008 / 0.07)
        assert r.uncertainty_contribution != pytest.approx(0.008)

    def test_negative_se_raises(self):
        with pytest.raises(ValueError):
            EffectRecord(
                rsid="rs1", gene="G", risk_allele="A",
                pathway_name="p", parameter_field="x",
                effect_per_allele=0.1, effect_type=EffectType.MULTIPLICATIVE,
                population=Population.EUR, confidence_tier="gold",
                effect_se=-0.01, pubmed_id="1", description="test",
            )

    def test_invalid_tier_raises(self):
        with pytest.raises(ValueError, match="confidence_tier"):
            EffectRecord(
                rsid="rs1", gene="G", risk_allele="A",
                pathway_name="p", parameter_field="x",
                effect_per_allele=0.1, effect_type=EffectType.MULTIPLICATIVE,
                population=Population.EUR, confidence_tier="platinum",
                effect_se=0.01, pubmed_id="1", description="test",
            )

    def test_population_frequency_defaults_to_none(self, tcf7l2_record):
        # A record with no verified frequency must default to None, not
        # a guessed placeholder -- biosim_seeds.population's genetic
        # sampling code relies on this to skip unverified variants
        # rather than assuming a value.
        assert tcf7l2_record.population_frequency is None

    def test_population_frequency_out_of_range_raises(self):
        with pytest.raises(ValueError, match="population_frequency"):
            EffectRecord(
                rsid="rs1", gene="G", risk_allele="A",
                pathway_name="p", parameter_field="x",
                effect_per_allele=0.1, effect_type=EffectType.MULTIPLICATIVE,
                population=Population.EUR, confidence_tier="gold",
                effect_se=0.01, pubmed_id="1", description="test",
                population_frequency=1.5,
            )
        with pytest.raises(ValueError, match="population_frequency"):
            EffectRecord(
                rsid="rs1", gene="G", risk_allele="A",
                pathway_name="p", parameter_field="x",
                effect_per_allele=0.1, effect_type=EffectType.MULTIPLICATIVE,
                population=Population.EUR, confidence_tier="gold",
                effect_se=0.01, pubmed_id="1", description="test",
                population_frequency=-0.1,
            )

    def test_population_frequency_boundary_values_are_fine(self):
        EffectRecord(
            rsid="rs1", gene="G", risk_allele="A",
            pathway_name="p", parameter_field="x",
            effect_per_allele=0.1, effect_type=EffectType.MULTIPLICATIVE,
            population=Population.EUR, confidence_tier="gold",
            effect_se=0.01, pubmed_id="1", description="test",
            population_frequency=0.0,
        )
        EffectRecord(
            rsid="rs1", gene="G", risk_allele="A",
            pathway_name="p", parameter_field="x",
            effect_per_allele=0.1, effect_type=EffectType.MULTIPLICATIVE,
            population=Population.EUR, confidence_tier="gold",
            effect_se=0.01, pubmed_id="1", description="test",
            population_frequency=1.0,
        )


# ---------------------------------------------------------------------------
# VariantEffectDB
# ---------------------------------------------------------------------------

class TestVariantEffectDB:
    def test_default_seeded(self):
        db = VariantEffectDB.default()
        assert db.count() == 7

    def test_has_tcf7l2(self):
        db = VariantEffectDB.default()
        assert db.has("rs7903146")

    def test_get_unknown_empty(self):
        db = VariantEffectDB.default()
        assert db.get("rs9999999") == []

    def test_all_rsids(self):
        db = VariantEffectDB.default()
        rsids = db.all_rsids()
        assert "rs7903146" in rsids
        assert "rs1801282" in rsids

    def test_all_default_variants_have_verified_population_frequency(self):
        # Added for biosim_seeds.population's genetic-variant sampling --
        # every one of the 7 default variants was verified live against
        # TOPMed (via dbSNP's REST API) this session, so none should be
        # left at the unverified-placeholder None, including rs35767
        # (IGF1) whose EFFECT SIZE citation remains unverified/bronze --
        # allele frequency is a separate, independently-answerable fact.
        db = VariantEffectDB.default()
        for rsid in db.all_rsids():
            record = db.get(rsid)[0]
            assert record.population_frequency is not None, (
                f"{rsid} ({record.gene}) has no verified population_frequency"
            )
            assert 0.0 < record.population_frequency < 1.0, (
                f"{rsid} ({record.gene}) frequency {record.population_frequency} "
                f"is not a plausible common-variant frequency"
            )

    def test_tcf7l2_frequency_matches_topmed(self):
        # Spot-check one exact value against the source fetched live
        # this session (api.ncbi.nlm.nih.gov/variation/v0/refsnp/7903146,
        # TOPMed study, risk allele T: allele_count=? total_count=264690).
        db = VariantEffectDB.default()
        record = db.get("rs7903146")[0]
        assert record.population_frequency == pytest.approx(0.2765)

    def test_filter_by_population(self):
        db = VariantEffectDB.default()
        eur_db = db.filter_by_population(Population.EUR)
        assert eur_db.count() == db.count()  # all default are EUR

    def test_filter_by_different_population(self):
        db = VariantEffectDB.default()
        eas_db = db.filter_by_population(Population.EAS)
        assert eas_db.count() == 0  # no EAS-specific records in default DB

    def test_add_and_get(self):
        db = VariantEffectDB()
        r = EffectRecord(
            rsid="rs999", gene="TEST", risk_allele="A",
            pathway_name="glucose_insulin_v1",
            parameter_field="Gb", effect_per_allele=0.1,
            effect_type=EffectType.ADDITIVE_ABSOLUTE,
            population=Population.EUR, confidence_tier="bronze",
            effect_se=0.05, pubmed_id="0", description="test",
        )
        db.add(r)
        assert db.has("rs999")
        assert len(db.get("rs999")) == 1


# ---------------------------------------------------------------------------
# VariantParser
# ---------------------------------------------------------------------------

class TestVariantParser:
    def test_parse_rsid_list(self):
        parser = VariantParser()
        result = parser.parse_rsid_list([
            ("rs7903146", "0/1"),
            ("rs1801282", "1/1"),
        ])
        assert result.n_variants == 2
        assert result.errors == []

    def test_parse_rsid_list_no_call_skipped(self):
        parser = VariantParser()
        result = parser.parse_rsid_list([("rs7903146", "--")])
        assert result.n_variants == 0
        assert len(result.skipped) == 1

    def test_parse_rsid_list_bad_genotype(self):
        parser = VariantParser()
        result = parser.parse_rsid_list([("rs7903146", "XYZ")])
        assert result.n_variants == 0
        assert len(result.errors) == 1

    def test_parse_23andme(self):
        text = (
            "# This data file generated by 23andMe\n"
            "rsid\tchromosome\tposition\tgenotype\n"
            "rs7903146\t10\t112998590\tCT\n"
            "rs1801282\t3\t12393125\tCC\n"
            "rs5219\t11\t17408630\t--\n"
        )
        parser = VariantParser()
        result = parser.parse_23andme(text)
        assert result.n_variants == 2
        assert "rs7903146" in result.rsids()
        assert "rs1801282" in result.rsids()

    def test_parse_23andme_no_call_skipped(self):
        text = "rs7903146\t10\t112998590\t--\n"
        parser = VariantParser()
        result = parser.parse_23andme(text)
        assert result.n_variants == 0
        assert len(result.skipped) >= 1

    def test_filter_to_known(self):
        parser = VariantParser()
        result = parser.parse_rsid_list([
            ("rs7903146", "0/1"),
            ("rs9999999", "1/1"),  # unknown
        ])
        filtered = parser.filter_to_known(result, {"rs7903146"})
        assert filtered.n_variants == 1
        assert filtered.variants[0].rsid == "rs7903146"

    def test_23andme_homozygous_dosage_reflects_actual_allele_not_always_two(self):
        # Regression: parse_23andme() never attached the raw allele
        # letters to the Variant, so dosage_for_allele() fell back to
        # "assume alt=risk", reporting 2 risk-allele copies for EVERY
        # homozygous call regardless of which letter was actually
        # observed. rs7903146 (TCF7L2)'s real risk allele is "T" (see
        # effect_db.py) -- a "CC" call is the PROTECTIVE homozygote
        # (0 real risk-allele copies), a "TT" call is the risk homozygote
        # (2 copies). Both must be distinguishable after parsing.
        text_cc = "rs7903146\t10\t112998590\tCC\n"
        text_tt = "rs7903146\t10\t112998590\tTT\n"
        v_cc = VariantParser().parse_23andme(text_cc).variants[0]
        v_tt = VariantParser().parse_23andme(text_tt).variants[0]
        assert v_cc.dosage_for_allele("T") == 0
        assert v_tt.dosage_for_allele("T") == 2

    def test_rsid_list_two_letter_code_dosage_reflects_actual_allele(self):
        # Same regression as above, via the parse_rsid_list() two-letter-
        # code path instead of parse_23andme().
        v_cc = VariantParser().parse_rsid_list([("rs7903146", "CC")]).variants[0]
        v_tt = VariantParser().parse_rsid_list([("rs7903146", "TT")]).variants[0]
        assert v_cc.dosage_for_allele("T") == 0
        assert v_tt.dosage_for_allele("T") == 2

    def test_rsid_list_vcf_style_still_has_no_observed_alleles(self):
        # VCF-style input ("0/1") carries no real allele letters -- must
        # NOT be misinterpreted as a 2-character allele code.
        v = VariantParser().parse_rsid_list([("rs7903146", "0/1")]).variants[0]
        assert v.observed_alleles is None

    def test_full_personalization_pipeline_distinguishes_protective_from_risk_homozygote(self):
        # End-to-end: the actual bug's real-world consequence. Before this
        # fix, both "CC" (protective) and "TT" (risk) homozygous 23andMe
        # calls for rs7903146 produced the IDENTICAL, maximal secretion_scale
        # reduction through the full GenotypePersonalizer pipeline.
        from biosim_genotype import GenotypePersonalizer
        cc_variants = VariantParser().parse_23andme("rs7903146\t10\t112998590\tCC\n").variants
        tt_variants = VariantParser().parse_23andme("rs7903146\t10\t112998590\tTT\n").variants
        result_cc = GenotypePersonalizer().personalize(cc_variants)
        result_tt = GenotypePersonalizer().personalize(tt_variants)
        assert result_cc.personalized_params.secretion_scale == pytest.approx(1.0)
        assert result_tt.personalized_params.secretion_scale < 1.0
        assert result_cc.personalized_params.secretion_scale != pytest.approx(
            result_tt.personalized_params.secretion_scale
        )

    def test_source_format_recorded(self):
        parser = VariantParser()
        r = parser.parse_rsid_list([("rs7903146", "0/1")])
        assert r.source_format == "rsid_list"


# ---------------------------------------------------------------------------
# EffectComposer
# ---------------------------------------------------------------------------

class TestEffectComposer:
    @pytest.fixture
    def db(self):
        return VariantEffectDB.default()

    @pytest.fixture
    def base_values(self):
        p = PathwayParameters()
        import dataclasses
        return dataclasses.asdict(p)

    def test_no_risk_alleles_no_effect(self, db, base_values):
        composer = EffectComposer()
        variants = [Variant("rs7903146", Genotype.HOMOZYGOUS_REF)]
        effects = composer.compose(variants, base_values, db)
        assert effects == []

    def test_tcf7l2_het_reduces_secretion(self, db, base_values):
        composer = EffectComposer()
        variants = [Variant("rs7903146", Genotype.HETEROZYGOUS)]
        effects = composer.compose(variants, base_values, db)
        sec_effect = next(e for e in effects if e.parameter_field == "secretion_scale")
        assert sec_effect.final_value < sec_effect.base_value
        assert sec_effect.percent_change == pytest.approx(-12.0, rel=0.01)

    def test_tcf7l2_hom_reduces_more(self, db, base_values):
        composer = EffectComposer()
        het = [Variant("rs7903146", Genotype.HETEROZYGOUS)]
        hom = [Variant("rs7903146", Genotype.HOMOZYGOUS_ALT)]
        eff_het = composer.compose(het, base_values, db)
        eff_hom = composer.compose(hom, base_values, db)
        het_scale = next(e.final_value for e in eff_het if e.parameter_field == "secretion_scale")
        hom_scale = next(e.final_value for e in eff_hom if e.parameter_field == "secretion_scale")
        assert hom_scale < het_scale

    def test_pparg_protective_increases_p3(self, db, base_values):
        composer = EffectComposer()
        variants = [Variant("rs1801282", Genotype.HETEROZYGOUS)]
        effects = composer.compose(variants, base_values, db)
        p3_effect = next(e for e in effects if e.parameter_field == "p3")
        assert p3_effect.final_value > p3_effect.base_value

    def test_mtnr1b_increases_gb(self, db, base_values):
        composer = EffectComposer()
        variants = [Variant("rs10830963", Genotype.HETEROZYGOUS)]
        effects = composer.compose(variants, base_values, db)
        gb_effect = next(e for e in effects if e.parameter_field == "Gb")
        assert gb_effect.final_value == pytest.approx(5.07, rel=0.001)

    def test_multiple_secretion_variants_compose(self, db, base_values):
        composer = EffectComposer()
        variants = [
            Variant("rs7903146", Genotype.HETEROZYGOUS),   # -12%
            Variant("rs5219",    Genotype.HETEROZYGOUS),   # -8%
        ]
        effects = composer.compose(variants, base_values, db)
        sec_effect = next(e for e in effects if e.parameter_field == "secretion_scale")
        # Combined: 1.0 * (1 - 0.12) * (1 - 0.08) = 0.88 * 0.92 = 0.8096
        assert sec_effect.final_value == pytest.approx(0.8096, rel=0.01)

    def test_unknown_rsid_ignored(self, db, base_values):
        composer = EffectComposer()
        variants = [Variant("rs9999999", Genotype.HETEROZYGOUS)]
        effects = composer.compose(variants, base_values, db)
        assert effects == []

    def test_population_mismatch_warning(self, db, base_values):
        composer = EffectComposer()
        variants = [Variant("rs7903146", Genotype.HETEROZYGOUS)]
        effects = composer.compose(variants, base_values, db, user_population=Population.EAS)
        sec_effect = next(e for e in effects if e.parameter_field == "secretion_scale")
        assert len(sec_effect.population_warnings) > 0

    def test_population_mismatch_warning_has_no_raw_rsid_or_gene(self, db, base_values):
        # Regression (privacy, rule 7): this used to embed the raw rsid
        # and gene name directly in the warning text -- caught by a
        # mandatory security-review of the sibling contributing_variants
        # fix, since population_warnings flows straight into
        # SimulationOutput.warnings (REST-serializable) via runner.py,
        # untouched by that fix.
        composer = EffectComposer()
        variants = [Variant("rs7903146", Genotype.HETEROZYGOUS)]
        effects = composer.compose(variants, base_values, db, user_population=Population.EAS)
        sec_effect = next(e for e in effects if e.parameter_field == "secretion_scale")
        warning = sec_effect.population_warnings[0]
        assert "rs7903146" not in warning
        assert "TCF7L2" not in warning
        assert "secretion_scale" in warning

    def test_genetic_cv_nonzero_when_applied(self, db, base_values):
        composer = EffectComposer()
        variants = [Variant("rs7903146", Genotype.HETEROZYGOUS)]
        effects = composer.compose(variants, base_values, db)
        sec = next(e for e in effects if e.parameter_field == "secretion_scale")
        assert sec.genetic_cv > 0.0

    def test_max_fold_change_capped(self, base_values):
        # Create a DB with an extreme effect
        db = VariantEffectDB()
        db.add(EffectRecord(
            rsid="rs_extreme", gene="TEST", risk_allele="A",
            pathway_name="glucose_insulin_v1",
            parameter_field="Gb",
            effect_per_allele=100.0,  # extreme
            effect_type=EffectType.MULTIPLICATIVE,
            population=Population.EUR, confidence_tier="bronze",
            effect_se=10.0, pubmed_id="0", description="extreme test",
        ))
        composer = EffectComposer(max_fold_change=3.0)
        variants = [Variant("rs_extreme", Genotype.HOMOZYGOUS_ALT)]
        effects = composer.compose(variants, base_values, db)
        effect = effects[0]
        fold = effect.final_value / effect.base_value
        assert fold <= 3.0 + 1e-9
        assert len(effect.population_warnings) > 0

    def test_max_fold_change_capped_for_decrease_too(self, base_values):
        # Regression: the cap used to be checked as
        # `abs(final/base_value) > max_fold_change` -- for a pure-decrease
        # multiplicative modifier (a product of terms each in (0,1)), that
        # ratio is always < 1.0 and can NEVER exceed max_fold_change no
        # matter how many decreasing variants stack, so extreme collapses
        # went completely uncapped and unwarned. Confirmed with this DB's
        # own real GWAS data separately: stacking all 9 real
        # secretion_scale-lowering variants homozygous collapsed
        # secretion_scale to 0.211 (79% loss) with zero warning.
        db = VariantEffectDB()
        db.add(EffectRecord(
            rsid="rs_extreme_decrease", gene="TEST", risk_allele="A",
            pathway_name="glucose_insulin_v1",
            parameter_field="secretion_scale",
            effect_per_allele=-0.45,  # large decrease, stays positive (not a sign flip)
            effect_type=EffectType.MULTIPLICATIVE,
            population=Population.EUR, confidence_tier="bronze",
            effect_se=0.05, pubmed_id="0", description="extreme decrease test",
        ))
        composer = EffectComposer(max_fold_change=3.0)
        variants = [Variant("rs_extreme_decrease", Genotype.HOMOZYGOUS_ALT)]
        effects = composer.compose(variants, base_values, db)
        effect = effects[0]
        fold = effect.final_value / effect.base_value
        assert fold >= (1.0 / 3.0) - 1e-9
        assert len(effect.population_warnings) > 0

    def test_extreme_decrease_stays_positive_without_the_old_1_percent_fallback(self, base_values):
        # An effect strong enough to flip final negative before capping
        # must still land at a small POSITIVE value (the symmetric floor,
        # base_value/max_fold_change) -- not the old, separate
        # "clamped to 1% of base value" fallback, which is now dead code
        # the fold-cap always catches first.
        db = VariantEffectDB()
        db.add(EffectRecord(
            rsid="rs_flip_negative", gene="TEST", risk_allele="A",
            pathway_name="glucose_insulin_v1",
            parameter_field="secretion_scale",
            effect_per_allele=-0.99,
            effect_type=EffectType.MULTIPLICATIVE,
            population=Population.EUR, confidence_tier="bronze",
            effect_se=0.05, pubmed_id="0", description="near-total collapse test",
        ))
        composer = EffectComposer(max_fold_change=4.0)
        variants = [Variant("rs_flip_negative", Genotype.HOMOZYGOUS_ALT)]
        effects = composer.compose(variants, base_values, db)
        effect = effects[0]
        assert effect.final_value == pytest.approx(effect.base_value / 4.0)
        assert effect.final_value > 0


# ---------------------------------------------------------------------------
# GenotypePersonalizer
# ---------------------------------------------------------------------------

class TestGenotypePersonalizer:
    def test_no_variants_returns_base_params(self):
        p = GenotypePersonalizer()
        result = p.personalize([])
        assert result.n_variants_applied == 0
        assert result.personalized_params == result.original_params

    def test_unknown_variants_no_change(self):
        p = GenotypePersonalizer()
        variants = [Variant("rs9999999", Genotype.HETEROZYGOUS)]
        result = p.personalize(variants)
        assert result.n_variants_matched == 0
        assert result.personalized_params == result.original_params

    def test_tcf7l2_risk_lowers_secretion_scale(self):
        p = GenotypePersonalizer()
        variants = [Variant("rs7903146", Genotype.HETEROZYGOUS)]
        result = p.personalize(variants)
        assert result.personalized_params.secretion_scale < 1.0
        assert result.n_parameters_modified == 1

    def test_pparg_protective_increases_p3(self):
        p = GenotypePersonalizer()
        variants = [Variant("rs1801282", Genotype.HETEROZYGOUS)]
        result = p.personalize(variants)
        base_p3 = PathwayParameters().p3
        assert result.personalized_params.p3 > base_p3

    def test_summary_structure(self):
        p = GenotypePersonalizer()
        result = p.personalize([Variant("rs7903146", Genotype.HETEROZYGOUS)])
        s = result.summary()
        assert "n_variants_input" in s
        assert "n_parameters_modified" in s
        assert "modified_parameters" in s

    def test_is_meaningfully_personalised(self):
        p = GenotypePersonalizer()
        result = p.personalize([Variant("rs7903146", Genotype.HETEROZYGOUS)])
        assert result.is_meaningfully_personalised

    def test_privacy_no_rsids_in_result(self):
        p = GenotypePersonalizer()
        result = p.personalize([Variant("rs7903146", Genotype.HETEROZYGOUS)])
        result_str = str(result.summary())
        assert "rs7903146" not in result_str

    def test_personalized_params_runs_ode(self):
        p = GenotypePersonalizer()
        variants = [
            Variant("rs7903146", Genotype.HETEROZYGOUS),
            Variant("rs1801282", Genotype.HETEROZYGOUS),
        ]
        result = p.personalize(variants)
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        system = ODESystem(pathway=pathway)
        sim = system.simulate_meal(high_sugar_meal(), duration_s=7200)
        assert sim.success
        peak_G, _ = sim.glucose_peak()
        assert peak_G > 5000

    def test_personalized_validate_y0_passes(self):
        p = GenotypePersonalizer()
        variants = [Variant("rs7903146", Genotype.HOMOZYGOUS_ALT)]
        result = p.personalize(variants)
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
        violations = pathway.validate_y0(tolerance=1e-3)
        assert violations == [], f"Personalised y0 not at steady state: {violations}"

    def test_homozygous_risk_has_larger_effect_than_heterozygous(self):
        p = GenotypePersonalizer()
        het = p.personalize([Variant("rs7903146", Genotype.HETEROZYGOUS)])
        hom = p.personalize([Variant("rs7903146", Genotype.HOMOZYGOUS_ALT)])
        assert hom.personalized_params.secretion_scale < het.personalized_params.secretion_scale

    def test_genetic_uncertainty_cv_positive_when_applied(self):
        p = GenotypePersonalizer()
        result = p.personalize([Variant("rs7903146", Genotype.HETEROZYGOUS)])
        assert result.genetic_uncertainty_cv > 0.0


class TestExtendedVariantCount:
    """Pins the documented count (CLAUDE.md, CV): 14 distinct variants,
    15 effect records, because IGF2BP2 rs4402960 has two records."""

    def test_extended_db_has_14_variants_and_15_records(self):
        from biosim_genotype.effect_db import extended_t2d_variants
        db = extended_t2d_variants()
        assert len(db.all_rsids()) == 14
        assert db.count() == 15
        assert len(db.get("rs4402960")) == 2


class TestUnverifiedIgf1VariantNote:
    """rs35767's +0.5 nM IGF-1 effect has no source; the user-facing note
    must say so, cite the conflicting evidence, and never echo the rsID."""

    def test_note_is_registered_for_rs35767_only(self):
        from biosim_genotype.effect_db import UNVERIFIED_VARIANT_NOTES
        assert set(UNVERIFIED_VARIANT_NOTES) == {"rs35767"}

    def test_note_says_unverified_and_cites_conflicting_pmids(self):
        import re
        from biosim_genotype.effect_db import UNVERIFIED_VARIANT_NOTES
        note = UNVERIFIED_VARIANT_NOTES["rs35767"]
        assert "UNVERIFIED" in note
        for pmid in ("30115944", "22689104", "21046444"):
            assert pmid in note
        assert not re.search(r"\brs\d{4,}\b", note)   # privacy rule 7
