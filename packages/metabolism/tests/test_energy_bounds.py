# SAVE TO: packages\metabolism\tests\test_energy_bounds.py
"""Tests for biosim_metabolism.energy_bounds (Track 4, MASTER_PLAN 4.3c)."""
import pytest

from biosim_metabolism.energy_bounds import (
    ATPM_HEPATIC_SPECIFIC_DERIVED_UMOL_KG_MIN,
    ATPM_WHOLE_BODY_DERIVED_UMOL_KG_MIN,
    EX_ALA_MEYER2002_UMOL_KG_MIN,
    EX_FA_WHOLE_BODY_UMOL_KG_MIN,
    EX_GLYC_MEYER2002_UMOL_KG_MIN,
    EX_LAC_MEYER2002_UMOL_KG_MIN,
    EX_O2_HEPATIC_SPECIFIC_UMOL_KG_MIN,
    EX_O2_WHOLE_BODY_UMOL_KG_MIN,
)


class TestEnergyBoundConstants:
    """
    Pins the exact literature-quoted/derived values so a future edit
    can't silently drift from the verified abstracts (Byrne et al.
    2005, PMID:15831804; Brundin & Wahren 1994, PMID:8177052; Bisschop
    et al. 2002, PMID:12067425; Rolfe & Brown 1997, PMID:9234964 --
    all independently re-verified live via PubMedFetcher this session).
    """

    def test_ex_o2_whole_body_matches_byrne_2005(self):
        # 2.6 mL O2/kg/min * (1/22.4 mol/L) * 1000 umol/mmol = 116.07...
        assert EX_O2_WHOLE_BODY_UMOL_KG_MIN == pytest.approx(116.1, abs=0.1)

    def test_ex_o2_hepatic_specific_matches_brundin_wahren_1994(self):
        # 62 mL/min / 70 kg = 0.8857 mL/kg/min -> umol/kg/min conversion
        assert EX_O2_HEPATIC_SPECIFIC_UMOL_KG_MIN == pytest.approx(39.5, abs=0.1)

    def test_ex_fa_whole_body_matches_bisschop_2002(self):
        assert EX_FA_WHOLE_BODY_UMOL_KG_MIN == pytest.approx(1.47)

    def test_atpm_whole_body_derived_matches_the_stacked_formula(self):
        # VO2 * 0.72 (Rolfe & Brown coupling fraction) * 5 (P/O=2.5 -> ATP per O2)
        assert ATPM_WHOLE_BODY_DERIVED_UMOL_KG_MIN == pytest.approx(
            EX_O2_WHOLE_BODY_UMOL_KG_MIN * 0.72 * 5.0
        )

    def test_atpm_hepatic_specific_derived_matches_the_stacked_formula(self):
        assert ATPM_HEPATIC_SPECIFIC_DERIVED_UMOL_KG_MIN == pytest.approx(
            EX_O2_HEPATIC_SPECIFIC_UMOL_KG_MIN * 0.72 * 5.0
        )

    def test_atpm_whole_body_is_larger_than_hepatic_specific(self):
        # Whole-body VO2 > splanchnic-specific VO2 -- the derived ATPM
        # estimates must preserve that same ordering.
        assert ATPM_WHOLE_BODY_DERIVED_UMOL_KG_MIN > ATPM_HEPATIC_SPECIFIC_DERIVED_UMOL_KG_MIN

    def test_carbon_precursor_bounds_match_meyer_2002_derivation(self):
        assert EX_LAC_MEYER2002_UMOL_KG_MIN == pytest.approx(4.0)
        assert EX_GLYC_MEYER2002_UMOL_KG_MIN == pytest.approx(2.5)
        assert EX_ALA_MEYER2002_UMOL_KG_MIN == pytest.approx(2.5)

    def test_hepatic_specific_o2_is_smaller_than_whole_body(self):
        # Splanchnic bed is a fraction of whole-body O2 consumption --
        # confirms the two constants weren't accidentally transposed.
        assert EX_O2_HEPATIC_SPECIFIC_UMOL_KG_MIN < EX_O2_WHOLE_BODY_UMOL_KG_MIN


class TestUmolPerKgMinToMMPerMin:
    def test_conversion_matches_flux_benchmark_convention(self):
        from biosim_metabolism.energy_bounds import EX_FA_WHOLE_BODY_UMOL_KG_MIN
        from biosim_metabolism.flux_benchmarks import (
            ADKINS_2003_BASAL_EGP,
            umol_per_kg_min_to_mM_per_min,
        )
        # Confirms the extracted, shared function reproduces
        # FluxBenchmark.value_mM_per_min()'s own result exactly --
        # RED-proofed at extraction time via the full test_flux_benchmarks.py
        # suite already passing unchanged.
        assert umol_per_kg_min_to_mM_per_min(
            ADKINS_2003_BASAL_EGP.value_umol_per_kg_min
        ) == pytest.approx(ADKINS_2003_BASAL_EGP.value_mM_per_min())
        # Sanity: also usable directly on a non-FluxBenchmark value.
        assert umol_per_kg_min_to_mM_per_min(EX_FA_WHOLE_BODY_UMOL_KG_MIN) > 0
