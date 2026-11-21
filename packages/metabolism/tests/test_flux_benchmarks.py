# SAVE TO: packages\metabolism\tests\test_flux_benchmarks.py
"""Tests for biosim_metabolism.flux_benchmarks (Track 4, MASTER_PLAN 4.2)."""
import pytest

from biosim_metabolism import (
    ADKINS_2003_BASAL_EGP,
    DROUIN_1998_GLUCAGON_STIMULATED_INCREMENT,
    ROSQVIST_2019_FASTED_HEPATIC_DNL,
    FluxBenchmark,
    FractionalContributionBenchmark,
)


class TestAdkins2003BasalEGP:
    """
    Adkins et al. 2003, PMID:12941759 -- verified live via
    PubMedFetcher this session. Pins the exact literature values so a
    future edit can't silently drift from the verified abstract quote:
    "13.1 +/- 1.3 ... micromol x kg(-1) x min(-1)" at the lowest
    (~75 pmol/L, near-basal) insulin infusion.
    """

    def test_pubmed_id(self):
        assert ADKINS_2003_BASAL_EGP.pubmed_id == "12941759"

    def test_value_matches_the_verified_abstract_quote(self):
        assert ADKINS_2003_BASAL_EGP.value_umol_per_kg_min == pytest.approx(13.1)
        assert ADKINS_2003_BASAL_EGP.stderr_umol_per_kg_min == pytest.approx(1.3)

    def test_caveat_present_and_non_empty(self):
        assert len(ADKINS_2003_BASAL_EGP.caveat) > 20

    def test_conversion_to_mM_per_min_at_default_body_weight(self):
        # Hand-computed: 13.1 umol/kg/min * 70 kg / 1000 (umol->mmol) /
        # (0.230 L/kg * 70 kg) = 0.917 mmol/min / 16.1 L = 0.05696 mM/min
        assert ADKINS_2003_BASAL_EGP.value_mM_per_min() == pytest.approx(0.056957, abs=1e-6)

    def test_stderr_conversion_to_mM_per_min(self):
        assert ADKINS_2003_BASAL_EGP.stderr_mM_per_min() == pytest.approx(0.0056522, abs=1e-6)

    def test_custom_body_weight_and_vg_scale_the_result(self):
        # A heavier body weight cancels out of the ratio (both the flux
        # AND the distribution volume scale with body_weight_kg) --
        # confirms the conversion is a genuine per-kg-normalized rate,
        # not an absolute-body-weight-dependent one.
        default = ADKINS_2003_BASAL_EGP.value_mM_per_min(body_weight_kg=70.0)
        heavier = ADKINS_2003_BASAL_EGP.value_mM_per_min(body_weight_kg=100.0)
        assert heavier == pytest.approx(default, abs=1e-9)

    def test_larger_vg_l_per_kg_lowers_the_converted_rate(self):
        # A larger distribution volume dilutes the same mass-flux into a
        # smaller concentration-rate.
        smaller_vg = ADKINS_2003_BASAL_EGP.value_mM_per_min(vg_l_per_kg=0.230)
        larger_vg = ADKINS_2003_BASAL_EGP.value_mM_per_min(vg_l_per_kg=0.286)
        assert larger_vg < smaller_vg

    def test_converted_value_is_the_same_order_of_magnitude_as_egp_basal_mm_min(self):
        # The explicit "unbelievable unit conversion -> documented
        # negative result" decision gate this row's own MASTER_PLAN
        # entry specifies -- confirms this real, tracer-measured basal
        # EGP rate, once converted through this platform's own
        # established distribution-volume convention, lands in the same
        # ballpark as (not wildly different from) the currently-uncited
        # EGP_BASAL_MM_MIN constant it will eventually help validate.
        from biosim_ode.hepatic import EGP_BASAL_MM_MIN
        converted = ADKINS_2003_BASAL_EGP.value_mM_per_min()
        ratio = converted / EGP_BASAL_MM_MIN
        assert 0.5 < ratio < 3.0, (
            f"Adkins 2003's converted basal EGP ({converted:.4f} mM/min) is "
            f"{ratio:.2f}x EGP_BASAL_MM_MIN ({EGP_BASAL_MM_MIN} mM/min) -- "
            f"outside a plausible same-order-of-magnitude range"
        )


class TestDrouin1998GlucagonStimulatedIncrement:
    """
    Drouin et al. 1998, PMID:9458743 -- verified live via PubMedFetcher
    this session. Untrained-subject, first-hour value:
    "7.4 +/- 1.6 ... mol.kg-1.min-1" (micromol, per the paper's own
    units section).
    """

    def test_pubmed_id(self):
        assert DROUIN_1998_GLUCAGON_STIMULATED_INCREMENT.pubmed_id == "9458743"

    def test_value_matches_the_verified_abstract_quote(self):
        assert DROUIN_1998_GLUCAGON_STIMULATED_INCREMENT.value_umol_per_kg_min == pytest.approx(7.4)
        assert DROUIN_1998_GLUCAGON_STIMULATED_INCREMENT.stderr_umol_per_kg_min == pytest.approx(1.6)

    def test_caveat_explicitly_warns_against_combining_with_adkins(self):
        caveat = DROUIN_1998_GLUCAGON_STIMULATED_INCREMENT.caveat.lower()
        assert "do not add directly" in caveat

    def test_conversion_to_mM_per_min(self):
        assert (
            DROUIN_1998_GLUCAGON_STIMULATED_INCREMENT.value_mM_per_min()
            == pytest.approx(0.032174, abs=1e-6)
        )

    def test_stimulated_increment_is_smaller_than_adkins_basal_rate(self):
        # Sanity check, not a scientific claim about the two studies
        # being directly comparable (they explicitly are not, per the
        # caveat above) -- just confirms both numbers are real,
        # finite, and in a sane relative order after conversion.
        assert (
            DROUIN_1998_GLUCAGON_STIMULATED_INCREMENT.value_mM_per_min()
            < ADKINS_2003_BASAL_EGP.value_mM_per_min()
        )


class TestFluxBenchmarkGenericConstruction:
    def test_construction_and_conversion_with_hand_picked_numbers(self):
        # A synthetic FluxBenchmark (not tied to any real citation) to
        # exercise the conversion formula independent of the two real
        # constants above.
        b = FluxBenchmark(
            name="synthetic", pubmed_id="0", value_umol_per_kg_min=10.0,
            stderr_umol_per_kg_min=0.0, population="test", condition="test",
        )
        # 10 umol/kg/min * 70 kg / 1000 / (0.230*70) = 0.7/16.1 = 0.043478...
        assert b.value_mM_per_min(body_weight_kg=70.0, vg_l_per_kg=0.230) == pytest.approx(
            0.043478, abs=1e-6
        )

    def test_default_caveat_is_empty_string(self):
        b = FluxBenchmark(
            name="synthetic", pubmed_id="0", value_umol_per_kg_min=1.0,
            stderr_umol_per_kg_min=0.0, population="test", condition="test",
        )
        assert b.caveat == ""


class TestRosqvist2019FastedHepaticDNL:
    """
    Rosqvist et al. 2019, PMID:30721918 -- verified live via
    PubMedFetcher (full text) this session (2026-08-12 deep-search
    plan, item F1). Pins the exact literature values so a future edit
    can't silently drift from the verified full-text quote: "the
    fasting hepatic DNL was 6.5% (range 0.1-38.7%)" (median, IQR
    3.8-11.5%).
    """

    def test_pubmed_id(self):
        assert ROSQVIST_2019_FASTED_HEPATIC_DNL.pubmed_id == "30721918"

    def test_fraction_matches_the_verified_full_text_quote(self):
        assert ROSQVIST_2019_FASTED_HEPATIC_DNL.fraction == pytest.approx(0.065)

    def test_fraction_range_matches_the_verified_iqr(self):
        assert ROSQVIST_2019_FASTED_HEPATIC_DNL.fraction_low == pytest.approx(0.038)
        assert ROSQVIST_2019_FASTED_HEPATIC_DNL.fraction_high == pytest.approx(0.115)

    def test_fraction_lies_within_its_own_reported_range(self):
        b = ROSQVIST_2019_FASTED_HEPATIC_DNL
        assert b.fraction_low <= b.fraction <= b.fraction_high

    def test_method_is_deuterated_water(self):
        assert "2H2O" in ROSQVIST_2019_FASTED_HEPATIC_DNL.method

    def test_caveat_present_and_non_empty(self):
        assert ROSQVIST_2019_FASTED_HEPATIC_DNL.caveat != ""

    def test_caveat_discloses_it_is_not_a_flux_rate(self):
        # This is the load-bearing distinction from FluxBenchmark -- a
        # future reader must not attempt to convert it via
        # umol_per_kg_min_to_mM_per_min().
        assert "not a flux rate" in ROSQVIST_2019_FASTED_HEPATIC_DNL.caveat

    def test_caveat_discloses_the_label_is_dosed_the_evening_before(self):
        # Regression: a Track 4.3d scoping pass found (and this session
        # independently re-verified live, full text) that the paper's
        # own Methods section dosed the 2H2O label "the evening before...
        # over the course of the evening" alongside "a low-fat dinner" --
        # only the final blood draw happens after an overnight fast, not
        # the labeled synthesis itself. The caveat previously (wrongly)
        # called this "a genuinely FASTED-state measurement... unlike
        # most DNL literature." Locks in the corrected disclosure so a
        # future edit can't silently revert to the overreaching claim.
        caveat = ROSQVIST_2019_FASTED_HEPATIC_DNL.caveat
        assert "NOT a purely fasted-state measurement" in caveat
        assert "evening before" in caveat

    def test_is_a_fractional_contribution_benchmark_not_a_flux_benchmark(self):
        assert isinstance(ROSQVIST_2019_FASTED_HEPATIC_DNL, FractionalContributionBenchmark)
        assert not isinstance(ROSQVIST_2019_FASTED_HEPATIC_DNL, FluxBenchmark)

    def test_has_no_value_mM_per_min_method(self):
        # Deliberate structural guard: FractionalContributionBenchmark
        # must NOT expose FluxBenchmark's own unit-conversion method --
        # a fraction has no such conversion, and accidentally inheriting
        # or duck-typing one would be a real category error.
        assert not hasattr(ROSQVIST_2019_FASTED_HEPATIC_DNL, "value_mM_per_min")


class TestFractionalContributionBenchmarkGenericConstruction:
    def test_construction_with_hand_picked_numbers(self):
        b = FractionalContributionBenchmark(
            name="synthetic", pubmed_id="0", fraction=0.10,
            fraction_low=0.05, fraction_high=0.20,
            population="test", condition="test", method="test",
        )
        assert b.fraction == pytest.approx(0.10)

    def test_default_caveat_is_empty_string(self):
        b = FractionalContributionBenchmark(
            name="synthetic", pubmed_id="0", fraction=0.10,
            fraction_low=0.05, fraction_high=0.20,
            population="test", condition="test", method="test",
        )
        assert b.caveat == ""
