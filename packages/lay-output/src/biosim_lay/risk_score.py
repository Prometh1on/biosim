# SAVE TO: packages\lay-output\src\biosim_lay\risk_score.py
"""
T2D and metabolic syndrome risk assessment from simulation outputs.

Maps postprandial simulation metrics to established epidemiological
risk indicators. This is NOT a clinical diagnostic tool — it translates
model outputs into population-level risk language for lay users.

Risk indicators implemented:
    1. Postprandial glucose excursion risk (vs WHO OGTT thresholds)
    2. Insulin resistance proxy (glucose AUC / insulin peak ratio)
    3. mTOR-AMPK balance (chronic anabolic vs catabolic state)
    4. Time-above-threshold metric (fraction of sim time with G > 7.8 mM)

Each indicator maps to a population percentile band and a plain-English
interpretation. No single indicator is diagnostic — they are presented
as a profile for informational purposes.

Reference thresholds:
    WHO 1999: 2h glucose < 7.8 mM = normal, 7.8-11.1 = IGT, > 11.1 = DM
    ADA 2023: fasting < 5.6 mM normal, 5.6-6.9 pre-DM, >= 7.0 DM
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional
import numpy as np

from .report import _HYPO_LEVEL1_MM, _HYPO_LEVEL2_MM, _secretion_appears_blunted


class RiskBand(str, Enum):
    LOW       = "low"        # < 25th percentile of population distribution
    MODERATE  = "moderate"   # 25th-75th percentile
    ELEVATED  = "elevated"   # 75th-90th percentile
    HIGH      = "high"       # > 90th percentile
    # A failed simulation, not a population percentile -- see
    # compute_risk_profile()'s own success-gating comment below. Mirrors
    # report.py's TrafficLight.GREY: never computed from real data, so it
    # must never look like a confident LOW ("healthy") result.
    UNKNOWN   = "unknown"


@dataclass(frozen=True)
class RiskIndicatorScore:
    """One risk indicator with its band, value, and interpretation."""
    name:          str
    value:         float
    unit:          str
    band:          RiskBand
    interpretation: str
    caveat:        str = ""


@dataclass
class MetabolicRiskProfile:
    """Complete metabolic risk profile from a simulation output."""
    indicators:    list[RiskIndicatorScore]
    overall_band:  RiskBand
    headline:      str
    disclaimer:    str = (
        "This is a research tool, not a medical test. "
        "These indicators are based on a mathematical model, not clinical measurements. "
        "Consult a healthcare professional for personal health assessment."
    )

    @property
    def n_elevated_or_high(self) -> int:
        return sum(
            1 for i in self.indicators
            if i.band in (RiskBand.ELEVATED, RiskBand.HIGH)
        )

    def to_dict(self) -> dict:
        # Regression: this used to omit `caveat` (per-indicator, e.g. the
        # blunted-secretion warning some indicators carry) and
        # `disclaimer` (the class-level "this is a research tool, not a
        # medical test" text) entirely -- both real fields with
        # meaningful content, silently dropped from every JSON/dict
        # consumer of this profile. Same false-reassurance bug class as
        # the identical gap just fixed in report.py's risk_section.
        return {
            "overall_band":  self.overall_band.value,
            "headline":      self.headline,
            "n_elevated":    self.n_elevated_or_high,
            "disclaimer":    self.disclaimer,
            "indicators":    [
                {
                    "name":           i.name,
                    "value":          round(i.value, 3),
                    "unit":           i.unit,
                    "band":           i.band.value,
                    "interpretation": i.interpretation,
                    "caveat":         i.caveat,
                }
                for i in self.indicators
            ],
        }


def compute_risk_profile(output) -> MetabolicRiskProfile:
    """
    Compute a metabolic risk profile from a SimulationOutput.

    Args:
        output: SimulationOutput from SimulationRunner or PlatformSession

    Returns:
        MetabolicRiskProfile with all indicators scored and an overall band
    """
    m   = output.expert_metrics
    sim = output.simulation
    indicators = []

    # Regression, found by the mandatory code-reviewer pass on the REST
    # exposure of this function (MASTER_PLAN 3.11b): unlike report.py's
    # assess_risk(), which gates every indicator on
    # overall_watchdog_confidence and turns everything TrafficLight.GREY
    # for a low-confidence/failed simulation, this function never checked
    # output.success at all. On SimulationRunner's own documented failure
    # branch (glucose_peak_mM=0, mtor_peak=0, ampk_nadir=1.0, etc.), the
    # unguarded arithmetic below produced overall_band=LOW and headline
    # "Metabolic indicators are all in a healthy range for this meal" --
    # a materially wrong, more-confident-sounding result than assess_risk()
    # gives the identical failed SimulationOutput, and (once wired into
    # biosim_rest) reachable by any REST client that ignores the separate
    # top-level `success` field. Fixed by short-circuiting before any real
    # indicator is computed, same "never let a derived-from-nothing result
    # read as a confident answer" principle as TrafficLight.GREY.
    if not sim.success:
        return MetabolicRiskProfile(
            indicators=[
                RiskIndicatorScore(
                    "Simulation status", 0.0, "", RiskBand.UNKNOWN,
                    "The simulation did not complete successfully, so no metabolic "
                    "risk indicators could be computed from it.",
                ),
            ],
            overall_band=RiskBand.UNKNOWN,
            headline="No risk assessment is available -- the underlying simulation failed.",
        )

    # --- 1. Peak glucose excursion ---
    #
    # Regression: the LOW/MODERATE bands used to be gated purely on
    # `excursion` (peak - the patient's OWN baseline) and never checked
    # the absolute peak_mM at all -- only the ELEVATED/HIGH bands did,
    # and those were only reachable once excursion >= 3.5. For a patient
    # whose own baseline is already elevated (e.g. the T2D disease-state
    # preset, ~6.5 mM), a small further rise landed the absolute peak
    # squarely in the WHO borderline/impaired-glucose-tolerance range
    # (7.8-11.1 mM) while still being called "a healthy response" --
    # reproduced directly: type2_diabetes_meal_scenario(carbs_g=25) gives
    # baseline=6.5, peak=8.4 mM, excursion=1.9 (< 2.0) -> "a healthy
    # response", while report.py's own "Peak blood sugar" indicator
    # (already gated on absolute value) correctly calls the same 8.4 mM
    # "borderline... reduced glucose tolerance" from the identical
    # SimulationOutput. Fixed by making the absolute WHO threshold the
    # primary gate (matches report.py's assess_risk()), with excursion
    # size only distinguishing LOW vs. MODERATE within the genuinely
    # below-threshold case, same as before this fix.
    peak_mM = m.glucose_peak_mM
    base_mM = m.glucose_baseline_uM / 1000.0
    excursion = peak_mM - base_mM

    if peak_mM >= 11.1:
        band = RiskBand.HIGH
        interp = (
            f"Blood sugar peak ({peak_mM:.1f} mmol/L) exceeds the WHO diabetic threshold "
            f"(11.1 mmol/L). This result is from a simulation, not a clinical test."
        )
    elif peak_mM >= 7.8:
        band = RiskBand.ELEVATED
        interp = (
            f"Blood sugar peak ({peak_mM:.1f} mmol/L) is in the borderline range. "
            f"The WHO threshold for impaired glucose tolerance is 7.8 mmol/L at 2 hours."
        )
    elif excursion >= 3.5:
        band = RiskBand.ELEVATED
        interp = f"Blood sugar peak ({peak_mM:.1f} mmol/L) is elevated but below clinical thresholds."
    elif excursion >= 2.0:
        band = RiskBand.MODERATE
        interp = f"Blood sugar peak ({peak_mM:.1f} mmol/L) is a moderate rise, typical for this meal."
    else:
        band = RiskBand.LOW
        interp = f"Blood sugar peak ({peak_mM:.1f} mmol/L) is a modest rise — a healthy response."
    indicators.append(RiskIndicatorScore("Peak glucose", peak_mM, "mmol/L", band, interp))

    # --- 2. Time above 7.8 mM (if time series available) ---
    if sim.success and hasattr(sim, "species"):
        G_mM = sim.species["G"] / 1000.0
        total_pts  = len(G_mM)
        above_pts  = int(np.sum(G_mM > 7.8))
        pct_above  = 100.0 * above_pts / max(total_pts, 1)

        if pct_above < 5:
            band   = RiskBand.LOW
            interp = f"Blood sugar stays below 7.8 mmol/L for {100-pct_above:.0f}% of the simulation."
        elif pct_above < 20:
            band   = RiskBand.MODERATE
            interp = f"Blood sugar is above 7.8 mmol/L for {pct_above:.0f}% of the simulation window."
        elif pct_above < 40:
            band   = RiskBand.ELEVATED
            interp = f"Blood sugar stays elevated (> 7.8 mmol/L) for {pct_above:.0f}% of the simulation."
        else:
            band   = RiskBand.HIGH
            interp = f"Blood sugar is elevated (> 7.8 mmol/L) for {pct_above:.0f}% of the simulation."
        indicators.append(RiskIndicatorScore("Time above 7.8 mmol/L", pct_above, "%", band, interp))

    # --- 3. Insulin resistance proxy (Matsuda-like index) ---
    # Lower insulin peak for same glucose excursion = more sensitive
    # IR proxy = glucose_AUC / insulin_peak (higher = more resistant)
    #
    # Regression: this divided by 1e6 while claiming (both in the old
    # comment and the "mM.s/pM" unit label below) to convert uM.s ->
    # mM.s -- that conversion is /1e3 (1 mM = 1000 uM); dividing by 1e6
    # actually produces M.s, a value 1000x too small. The bug was
    # invisible in every test because MetabolicRiskProfile.overall_band
    # is the WORST of all four indicators, and the other three (peak
    # glucose, time-above-threshold, mTOR/AMPK) already differentiated
    # T2D from healthy on their own -- this one silently stayed LOW for
    # every disease-state combination this platform can simulate,
    # including T2D + hypothyroidism + age 90, and told confirmed T2D
    # cases their "insulin response appears efficient."
    #
    # The threshold constants (0.8/1.5/2.5) turned out to have never
    # been validated against real simulated values either -- fixing the
    # unit bug alone pushes every scenario's value into the hundreds,
    # blowing past the old HIGH cutoff regardless of disease state. Both
    # were broken together, so both needed re-deriving, not a one-line
    # divisor swap. Recalibrated by directly running this platform's own
    # full disease-severity range (healthy meals of varying size, stress,
    # age, insulin-resistant-only, hypothyroid-only, T2D, T2D+comorbidities)
    # and choosing cutoffs that separate the two clusters actually
    # produced: healthy/mild profiles land ~175-215, confirmed T2D
    # (regardless of additional comorbidities) lands ~450-525. This is an
    # internal, simulation-behavior calibration, NOT a clinical or
    # epidemiological cutoff -- this indicator was never literature-
    # grounded to begin with (see its own caveat below), so there is no
    # real-world threshold to anchor against; replace these if a real
    # Matsuda/HOMA-IR-equivalent literature mapping is ever adopted.
    insulin_peak_pM = m.insulin_peak_pM
    glucose_auc     = m.glucose_auc_uM_s / 1e3   # uM.s -> mM.s

    if insulin_peak_pM > 0:
        ir_proxy = glucose_auc / insulin_peak_pM
        if ir_proxy < 250.0:
            band   = RiskBand.LOW
            interp = "Insulin response appears efficient relative to glucose load."
        elif ir_proxy < 400.0:
            band   = RiskBand.MODERATE
            interp = "Insulin response is in a typical range."
        elif ir_proxy < 500.0:
            band   = RiskBand.ELEVATED
            interp = "Insulin response may be disproportionate — possible early insulin resistance."
        else:
            band   = RiskBand.HIGH
            interp = "High glucose relative to insulin response — consistent with insulin resistance pattern."
        indicators.append(RiskIndicatorScore(
            "Insulin efficiency proxy", round(ir_proxy, 2), "mM·s/pM", band, interp,
            caveat="This is a model-derived proxy, not a clinical insulin resistance test.",
        ))

    # --- 4. mTOR-AMPK balance ---
    # mTOR/AMPK are driven only by insulin level in this model -- a
    # disease state that blunts insulin SECRETION (not resistance, e.g.
    # T2D's secretion_scale=0.60) lowers this ratio for a reason unrelated
    # to genuine metabolic balance, so a T2D patient can score LOW
    # ("balanced") here even at a much larger carb load than a healthy
    # comparator. See report._secretion_appears_blunted's docstring for
    # the full explanation; surfaced as an explicit caveat rather than
    # recalibrating the band thresholds themselves, which would need a
    # principled clinical basis this platform doesn't have.
    mtor_peak   = m.mtor_peak
    ampk_nadir  = m.ampk_nadir
    mtor_ampk_ratio = mtor_peak / max(ampk_nadir, 0.01)

    if mtor_ampk_ratio < 2.0:
        band   = RiskBand.LOW
        interp = "Balanced anabolic/catabolic signalling after meal."
    elif mtor_ampk_ratio < 4.0:
        band   = RiskBand.MODERATE
        interp = "Moderate shift toward anabolic (mTOR) signalling — expected after a meal."
    elif mtor_ampk_ratio < 8.0:
        band   = RiskBand.ELEVATED
        interp = "Strong mTOR activation relative to AMPK — sustained high anabolic state."
    else:
        band   = RiskBand.HIGH
        interp = (
            "Very high mTOR/AMPK ratio. Chronic high-carbohydrate meals that repeatedly "
            "produce this pattern may contribute to mTOR-related metabolic stress."
        )
    mtor_ampk_caveat = "mTOR and AMPK activity are model-estimated, not clinically measured."
    if _secretion_appears_blunted(output):
        mtor_ampk_caveat += (
            " This person's insulin secretion appears reduced (personalisation-adjusted), "
            "which independently lowers this ratio -- a LOW/MODERATE result here may reflect "
            "blunted insulin secretion rather than genuine metabolic balance."
        )
    indicators.append(RiskIndicatorScore(
        "mTOR/AMPK ratio", round(mtor_ampk_ratio, 2), "", band, interp,
        caveat=mtor_ampk_caveat,
    ))

    # --- 5. Hypoglycaemia risk (glucose nadir) ---
    # Battelino et al. 2019, PMID:31177185 (Diabetes Care) -- see
    # report.py's _HYPO_LEVEL1_MM/_HYPO_LEVEL2_MM docstring for the full
    # citation and reasoning (imported from there rather than duplicated,
    # so the two files' thresholds can't silently drift apart). Added
    # alongside report.py's own "Hypoglycaemia risk" indicator, same
    # milestone -- this platform's SGLT2i/sulfonylurea mechanisms can
    # genuinely push simulated glucose into clinical hypoglycaemia, and
    # nothing in this file's own output ever checked for it before.
    nadir_mM = m.glucose_nadir_mM
    if nadir_mM < _HYPO_LEVEL2_MM:
        band = RiskBand.HIGH
        interp = (
            f"Lowest simulated blood sugar ({nadir_mM:.1f} mmol/L) is below the clinically "
            f"significant hypoglycaemia threshold (3.0 mmol/L)."
        )
    elif nadir_mM < _HYPO_LEVEL1_MM:
        band = RiskBand.ELEVATED
        interp = (
            f"Lowest simulated blood sugar ({nadir_mM:.1f} mmol/L) is below the "
            f"hypoglycaemia alert value (3.9 mmol/L)."
        )
    else:
        band = RiskBand.LOW
        interp = f"Lowest simulated blood sugar ({nadir_mM:.1f} mmol/L) stays above the hypoglycaemia alert value."
    # Same caveat as report.py's own "Hypoglycaemia risk" indicator, found
    # by the mandatory security-reviewer pass: a LOW result here reads as
    # "no hypoglycaemia risk," but no currently-shipped scenario
    # combination reaches ELEVATED/HIGH at all, and this platform's own
    # drug-effect calibrations are independently documented elsewhere as
    # sometimes understating real clinical effect sizes.
    hypo_calibration_caveat = (
        "This reflects the simulated trajectory under this platform's current "
        "drug-effect calibrations, not a guarantee against real hypoglycaemia."
        if band == RiskBand.LOW else ""
    )
    indicators.append(RiskIndicatorScore(
        "Hypoglycaemia risk", round(nadir_mM, 2), "mmol/L", band, interp,
        caveat=hypo_calibration_caveat,
    ))

    # --- Overall band: worst of all indicators ---
    band_order = [RiskBand.LOW, RiskBand.MODERATE, RiskBand.ELEVATED, RiskBand.HIGH]
    worst = max(indicators, key=lambda i: band_order.index(i.band))
    overall = worst.band

    if overall == RiskBand.LOW:
        headline = "Metabolic indicators are all in a healthy range for this meal."
    elif overall == RiskBand.MODERATE:
        headline = "Most metabolic indicators are typical. One or more are worth monitoring."
    elif overall == RiskBand.ELEVATED:
        headline = "Some metabolic indicators are elevated. Consider lower-GI meals or more physical activity."
    else:
        headline = "Several metabolic indicators are high. This pattern, if repeated frequently, may be worth discussing with a healthcare professional."

    return MetabolicRiskProfile(
        indicators=indicators,
        overall_band=overall,
        headline=headline,
    )
