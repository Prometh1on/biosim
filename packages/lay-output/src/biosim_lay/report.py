# SAVE TO: packages\lay-output\src\biosim_lay\report.py
"""
Risk indicators and plain language report generator.

RiskIndicator provides traffic-light assessment of simulated glucose profiles.
PlainLanguageReport generates a structured, jargon-free summary for lay users.

Clinical thresholds used:
    Normal postprandial glucose:   < 7.8 mM (140 mg/dL) at 2h
    Impaired glucose tolerance:    7.8 - 11.1 mM
    Diabetic range:               >= 11.1 mM (200 mg/dL)

    Source: WHO 1999 criteria; ADA Standards of Medical Care 2023.

IMPORTANT: These thresholds are reference values for population health research.
This tool does not provide clinical diagnosis. Users should consult a
healthcare professional for personal health decisions.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional
import numpy as np


class TrafficLight(str, Enum):
    GREEN  = "green"   # within normal range
    AMBER  = "amber"   # borderline / worth monitoring
    RED    = "red"     # outside normal range
    GREY   = "grey"    # insufficient confidence to classify


@dataclass(frozen=True)
class RiskIndicator:
    """
    Traffic-light assessment of one aspect of the simulation.
    """
    label:       str           # e.g. "Peak blood sugar"
    value:       float         # The numeric value being assessed
    unit:        str           # e.g. "mmol/L"
    traffic_light: TrafficLight
    description: str           # Plain English explanation
    caveat:      str = ""      # Uncertainty caveat


@dataclass
class RiskProfile:
    """
    Complete set of risk indicators for a simulation output.
    """
    indicators:  list[RiskIndicator]
    overall:     TrafficLight
    summary:     str

    @property
    def has_red(self) -> bool:
        return any(r.traffic_light == TrafficLight.RED for r in self.indicators)

    @property
    def has_amber(self) -> bool:
        return any(r.traffic_light == TrafficLight.AMBER for r in self.indicators)

    def indicator(self, label: str) -> Optional[RiskIndicator]:
        for r in self.indicators:
            if r.label == label:
                return r
        return None


# Glucose thresholds in mM
_GLUCOSE_NORMAL_PEAK = 7.8
_GLUCOSE_IMPAIRED_PEAK = 11.1

# Fasting glucose thresholds in mM (ADA 2023, cited in this module's own
# docstring -- until this fix, never actually checked against anywhere in
# this file; only the postprandial peak thresholds above were used).
_FASTING_NORMAL_MM = 5.6
_FASTING_DIABETIC_MM = 7.0

# Battelino T, Danne T, Bergenstal RM, et al. "Clinical Targets for
# Continuous Glucose Monitoring Data Interpretation: Recommendations
# From the International Consensus on Time in Range." Diabetes Care.
# 2019;42(8):1593-1603. PMID:31177185 (verified live via NCBI efetch,
# full text): Table 2's standardized CGM metrics define "Time below
# range (TBR)" Level 1 as 54-69 mg/dL (3.0-3.8 mmol/L) and Level 2 as
# <54 mg/dL (<3.0 mmol/L) -- the standard tiered hypoglycaemia
# thresholds used throughout modern diabetes-technology/CGM research.
# Used by the "Hypoglycaemia risk" indicator below and by risk_score.py
# (imported from here to avoid a second, independently-drifting copy).
_HYPO_LEVEL1_MM = 3.9   # below this = "Level 1" hypoglycaemia alert value
_HYPO_LEVEL2_MM = 3.0   # below this = "Level 2", clinically significant

# mTOR thresholds
_MTOR_HIGH = 0.70

# AMPK thresholds
_AMPK_LOW = 0.25

# Below this, secretion_scale is considered meaningfully reduced (not
# just ordinary between-person variation) for the mTOR/AMPK caveat below.
_BLUNTED_SECRETION_SCALE = 0.9


def _uses_gi_placeholder(output) -> bool:
    """True when a meal's timing relied on the uncalibrated GI slope (the
    runner's for_gi() warning, or a gi_meal_scenario() description)."""
    from biosim_ode.gut_absorption import GI_SLOPE_PLACEHOLDER_NOTE

    texts = list(getattr(output, "warnings", None) or [])
    texts.append(getattr(output, "scenario_description", "") or "")
    return any(GI_SLOPE_PLACEHOLDER_NOTE in t for t in texts)


def _secretion_appears_blunted(output) -> bool:
    """
    True if this output's personalisation reduced secretion_scale below
    a healthy-range value (e.g. the T2D disease-state preset, secretion_scale=0.60,
    or a genetic variant with a similar effect).

    Why this matters for mTOR/AMPK specifically: in this model, mTOR and
    AMPK activity are driven ONLY by insulin level (pathway.py's
    hill_activation(I, EC50_mTOR, ...) / hill_inhibition(I, EC50_AMPK, ...)),
    never by glucose directly. A disease state that blunts insulin
    SECRETION (impaired beta-cell output, not insulin resistance) lowers
    the postprandial insulin excursion and, through this same-insulin-only
    architecture, also lowers mTOR/AMPK excursion -- so a T2D patient's
    reading here can come out "healthier-looking" than a healthy
    comparator eating the identical meal, for a reason that has nothing
    to do with genuinely balanced anabolic/catabolic signalling. Checked
    across physiological AND genetic personalisation paths (both expose
    .personalized_params.secretion_scale) since either can produce this.
    """
    for personalization in (output.physiological_personalization, output.genetic_personalization):
        if personalization is None:
            continue
        params = getattr(personalization, "personalized_params", None)
        if params is not None and params.secretion_scale < _BLUNTED_SECRETION_SCALE:
            return True
    return False


def assess_risk(output, confidence_threshold: float = 0.7) -> RiskProfile:
    """
    Compute a traffic-light risk profile from a SimulationOutput.

    Args:
        output:               SimulationOutput from SimulationRunner
        confidence_threshold: Below this watchdog score, all indicators turn grey

    Returns:
        RiskProfile with individual indicators and an overall assessment
    """
    m = output.expert_metrics
    is_confident = m.overall_watchdog_confidence >= confidence_threshold
    caveat = (
        "" if is_confident
        else "Low model confidence — treat these indicators as approximate."
    )

    indicators = []
    baseline_mM = m.glucose_baseline_uM / 1000.0

    # 0. Fasting/baseline glucose
    # Regression: this file's own docstring cites ADA fasting criteria
    # (normal <5.6, pre-diabetic 5.6-6.9, diabetic >=7.0 mM) specifically
    # to justify its clinical grounding, but until this fix nothing in
    # the file ever actually checked a value against them -- only the
    # postprandial peak thresholds (7.8/11.1) were used. Without this,
    # a patient whose own baseline is already elevated (e.g. the T2D
    # disease-state preset's Gb x1.3) had no indicator surfacing that at
    # all, and the "Recovery time" indicator below would call a fast
    # return to that elevated baseline "a healthy response."
    if not is_confident:
        gl = TrafficLight.GREY
        desc = f"Fasting blood sugar estimated at {baseline_mM:.1f} mmol/L (low confidence)."
    elif baseline_mM < _FASTING_NORMAL_MM:
        gl = TrafficLight.GREEN
        desc = f"Fasting blood sugar of {baseline_mM:.1f} mmol/L is within the normal range (<5.6 mmol/L)."
    elif baseline_mM < _FASTING_DIABETIC_MM:
        gl = TrafficLight.AMBER
        desc = (
            f"Fasting blood sugar of {baseline_mM:.1f} mmol/L is in the pre-diabetic range "
            f"(5.6-6.9 mmol/L)."
        )
    else:
        gl = TrafficLight.RED
        desc = (
            f"Fasting blood sugar of {baseline_mM:.1f} mmol/L is in the diabetic range "
            f"(>=7.0 mmol/L)."
        )
    indicators.append(RiskIndicator(
        "Fasting blood sugar", baseline_mM, "mmol/L", gl, desc, caveat
    ))
    baseline_is_elevated = baseline_mM >= _FASTING_NORMAL_MM

    # 1. Peak glucose
    peak_mM = m.glucose_peak_mM
    if not is_confident:
        gl = TrafficLight.GREY
        desc = f"Peak blood sugar estimated at {peak_mM:.1f} mmol/L (low confidence)."
    elif peak_mM < _GLUCOSE_NORMAL_PEAK:
        gl = TrafficLight.GREEN
        desc = f"Peak blood sugar of {peak_mM:.1f} mmol/L is within the normal range (<7.8 mmol/L)."
    elif peak_mM < _GLUCOSE_IMPAIRED_PEAK:
        gl = TrafficLight.AMBER
        desc = (
            f"Peak blood sugar of {peak_mM:.1f} mmol/L is in the borderline range "
            f"(7.8–11.1 mmol/L). This may indicate reduced glucose tolerance."
        )
    else:
        gl = TrafficLight.RED
        desc = (
            f"Peak blood sugar of {peak_mM:.1f} mmol/L is in the high range "
            f"(>11.1 mmol/L). This is above the normal postprandial threshold."
        )
    indicators.append(RiskIndicator(
        "Peak blood sugar", peak_mM, "mmol/L", gl, desc, caveat
    ))

    # 1.5. Hypoglycaemia risk (glucose nadir) -- see _HYPO_LEVEL1_MM/
    # _HYPO_LEVEL2_MM's own docstring for the citation. Added specifically
    # because this platform's own SGLT2i (insulin-independent renal
    # clearance) and sulfonylurea (glucose-independent secretion)
    # mechanisms can genuinely push simulated glucose into clinical
    # hypoglycaemia, and until this indicator existed nothing in this
    # file's output ever checked for it -- the same "computed a real
    # mechanism, never wired a user-facing signal" bug class already
    # caught and fixed several times in this file's own history (see
    # _secretion_appears_blunted's docstring for a prior instance).
    nadir_mM = m.glucose_nadir_mM
    if not is_confident:
        gl = TrafficLight.GREY
        desc = f"Lowest simulated blood sugar estimated at {nadir_mM:.1f} mmol/L (low confidence)."
    elif nadir_mM < _HYPO_LEVEL2_MM:
        gl = TrafficLight.RED
        desc = (
            f"Lowest simulated blood sugar of {nadir_mM:.1f} mmol/L is below the clinically "
            f"significant hypoglycaemia threshold (<3.0 mmol/L)."
        )
    elif nadir_mM < _HYPO_LEVEL1_MM:
        gl = TrafficLight.AMBER
        desc = (
            f"Lowest simulated blood sugar of {nadir_mM:.1f} mmol/L is below the hypoglycaemia "
            f"alert value (<3.9 mmol/L)."
        )
    else:
        gl = TrafficLight.GREEN
        desc = (
            f"Lowest simulated blood sugar of {nadir_mM:.1f} mmol/L stays above the "
            f"hypoglycaemia alert value (3.9 mmol/L)."
        )
    # Found by the mandatory security-reviewer pass: a GREEN result here,
    # on its own, reads as "no hypoglycaemia risk" -- but this platform's
    # own drug-effect calibrations are independently documented elsewhere
    # in CLAUDE.md as sometimes UNDERSTATING real clinical effect sizes
    # (e.g. the metformin/semaglutide benchmarks' documented deviations),
    # and no currently-shipped scenario combination reaches this
    # indicator's AMBER/RED bands at all (verified directly before
    # writing this indicator -- see the test suite's own docstring).
    # Disclosed explicitly rather than left implicit, same "hedge false
    # reassurance rather than let a correct-as-computed result be
    # over-read" pattern already used for the mTOR/AMPK and recovery-time
    # indicators above.
    hypo_calibration_caveat = (
        " This reflects the simulated trajectory under this platform's current "
        "drug-effect calibrations, not a guarantee against real hypoglycaemia."
        if gl == TrafficLight.GREEN else ""
    )
    indicators.append(RiskIndicator(
        "Hypoglycaemia risk", nadir_mM, "mmol/L", gl, desc, caveat + hypo_calibration_caveat
    ))

    # 2. Recovery time
    # Regression: "return time" measures return to the PATIENT'S OWN
    # baseline (glucose_baseline_uM = G[0]), not a population-normal
    # value. For the T2D disease-state preset (Gb x1.3, ~6.5 mM --
    # already pre-diabetic by the fasting thresholds above), a fast
    # return to that elevated baseline was unconditionally described as
    # "near-normal"/"a healthy response," the opposite of what a T2D
    # patient's report should say. Fast-return-to-an-elevated-baseline is
    # capped at AMBER and reworded, regardless of how fast it is --
    # returning to your own baseline quickly says nothing about whether
    # that baseline itself is healthy.
    if m.time_to_glucose_return_s is not None:
        return_min = m.time_to_glucose_return_s / 60.0
        if not is_confident:
            gl = TrafficLight.GREY
            desc = f"Blood sugar return time estimated at {return_min:.0f} minutes (low confidence)."
        elif return_min < 90:
            if baseline_is_elevated:
                gl = TrafficLight.AMBER
                desc = (
                    f"Blood sugar returns to your own baseline ({baseline_mM:.1f} mmol/L) "
                    f"within {return_min:.0f} minutes. Note: this baseline is itself above "
                    f"the normal fasting range (<5.6 mmol/L), so a fast return does not by "
                    f"itself indicate a healthy response."
                )
            else:
                gl = TrafficLight.GREEN
                desc = f"Blood sugar returns to near-normal within {return_min:.0f} minutes — a healthy response."
        elif return_min < 150:
            gl = TrafficLight.AMBER
            desc = f"Blood sugar takes around {return_min:.0f} minutes to return — slightly prolonged."
        else:
            gl = TrafficLight.RED
            desc = f"Blood sugar takes {return_min:.0f} minutes or more to return — a prolonged response."
        indicators.append(RiskIndicator(
            "Recovery time", return_min, "minutes", gl, desc, caveat
        ))
    elif (sim := getattr(output, "simulation", None)) is not None             and sim.success and len(sim.times_s) > 0:
        # Not back to baseline by the end of the window: the slowest case,
        # so it must be flagged, not silently dropped (reachable for a T2D
        # 75 g meal since B7, 2026-11-16).
        # Not back by W minutes means the return takes longer than W, so the
        # same thresholds as above apply: W >= 150 -> RED, 90-150 -> AMBER;
        # a window under 90 minutes is too short to judge (GREY).
        window_min = float(sim.times_s[-1]) / 60.0
        if not is_confident or window_min < 90:
            gl = TrafficLight.GREY
        elif window_min < 150:
            gl = TrafficLight.AMBER
        else:
            gl = TrafficLight.RED
        desc = (
            f"Blood sugar has not returned to your own baseline within the "
            f"{window_min:.0f} minutes simulated, so it takes longer than that."
        )
        indicators.append(RiskIndicator(
            "Recovery time", window_min, "minutes", gl, desc, caveat
        ))

    # mTOR/AMPK are driven only by insulin level in this model (see
    # _secretion_appears_blunted's docstring) -- a disease state that
    # blunts insulin SECRETION (not resistance) can make both readings
    # look "healthier" for a reason unrelated to genuine metabolic
    # balance. Surfaced as an explicit caveat rather than silently
    # recalibrating the thresholds themselves, which would need a
    # principled clinical basis this platform doesn't have (same
    # discipline as this session's other "flag the gap honestly, don't
    # force an unverified fix" decisions).
    secretion_caveat = (
        " This person's insulin secretion appears reduced (personalisation-adjusted), "
        "which independently lowers this model's mTOR/AMPK readings -- a moderate/balanced "
        "result here may reflect blunted insulin secretion rather than genuine metabolic balance."
        if _secretion_appears_blunted(output) else ""
    )

    # 3. mTOR activation (cell growth / storage signal)
    mtor = m.mtor_peak
    if not is_confident:
        gl = TrafficLight.GREY
        desc = "mTOR activation level uncertain."
    elif mtor < _MTOR_HIGH:
        gl = TrafficLight.GREEN
        desc = f"Cell growth signalling (mTOR) activation is moderate ({mtor:.0%})."
    else:
        gl = TrafficLight.AMBER
        desc = (
            f"Cell growth signalling (mTOR) activation is high ({mtor:.0%}). "
            f"This reflects strong anabolic signalling, which is normal after a large meal."
        )
    indicators.append(RiskIndicator("mTOR activation", mtor, "", gl, desc, caveat + secretion_caveat))

    # 4. AMPK suppression (energy sensor)
    ampk = m.ampk_nadir
    if not is_confident:
        gl = TrafficLight.GREY
        desc = "Energy sensor (AMPK) level uncertain."
    elif ampk > _AMPK_LOW:
        gl = TrafficLight.GREEN
        desc = f"Energy sensor (AMPK) remains active ({ampk:.0%}), indicating balanced energy status."
    else:
        gl = TrafficLight.AMBER
        desc = (
            f"Energy sensor (AMPK) is substantially suppressed ({ampk:.0%}), "
            f"indicating a strong shift toward energy storage mode."
        )
    indicators.append(RiskIndicator("AMPK suppression", ampk, "", gl, desc, caveat + secretion_caveat))

    # Overall
    if any(r.traffic_light == TrafficLight.RED for r in indicators):
        overall = TrafficLight.RED
        summary = "One or more indicators are outside the normal range. See details below."
    elif any(r.traffic_light == TrafficLight.AMBER for r in indicators):
        overall = TrafficLight.AMBER
        summary = "Some indicators are in the borderline range. This may be worth monitoring."
    elif any(r.traffic_light == TrafficLight.GREY for r in indicators):
        overall = TrafficLight.GREY
        summary = "Simulation confidence is low. Results are approximate."
    else:
        overall = TrafficLight.GREEN
        summary = "All indicators are within the normal range for this meal."

    return RiskProfile(indicators=indicators, overall=overall, summary=summary)


@dataclass
class PlainLanguageReport:
    """
    A structured plain-language report for a lay user.
    Contains no jargon. Every technical term is explained.
    """
    headline:           str
    glucose_section:    str
    energy_section:     str
    risk_section:       str
    comparison_section: Optional[str]
    disclaimer:         str
    risk_profile:       RiskProfile

    def as_text(self) -> str:
        """Return the full report as a plain text string."""
        parts = [
            "=" * 60,
            self.headline,
            "=" * 60,
            "",
            "BLOOD SUGAR",
            self.glucose_section,
            "",
            "ENERGY",
            self.energy_section,
            "",
            "HEALTH INDICATORS",
            self.risk_section,
        ]
        if self.comparison_section:
            parts += ["", "COMPARISON", self.comparison_section]
        parts += ["", "-" * 60, self.disclaimer]
        return "\n".join(parts)


def generate_report(
    output,
    comparison_output=None,
) -> PlainLanguageReport:
    """
    Generate a plain language report from a SimulationOutput.

    Args:
        output:             SimulationOutput from SimulationRunner
        comparison_output:  Optional second SimulationOutput for comparison

    Returns:
        PlainLanguageReport with all sections populated
    """
    m = output.expert_metrics
    ls = output.lay_summary
    risk = assess_risk(output)

    headline = f"Results for: {output.scenario_name}"

    # Glucose section
    peak_min = round(m.glucose_peak_time_s / 60.0)
    peak_mM = round(m.glucose_peak_mM, 1)
    baseline_mM = round(m.glucose_baseline_uM / 1000.0, 1)
    glucose_section = (
        f"After eating, your estimated blood sugar rises to around {peak_mM} mmol/L "
        f"at about {peak_min} minutes.\n"
        f"A normal response is below 7.8 mmol/L. "
    )
    if m.time_to_glucose_return_s:
        return_min = round(m.time_to_glucose_return_s / 60.0)
        # Regression: this used to say "returns toward normal" regardless
        # of what value it was actually returning to -- glucose_baseline_uM
        # is the patient's OWN starting point (G[0]), which can itself be
        # elevated (e.g. the T2D disease-state preset). "Toward normal" is
        # only true if that baseline is itself in the normal fasting range.
        if baseline_mM >= _FASTING_NORMAL_MM:
            glucose_section += (
                f"Blood sugar returns to your own baseline ({baseline_mM} mmol/L) after "
                f"around {return_min} minutes. Note: this baseline itself is above the "
                f"normal fasting range (<5.6 mmol/L)."
            )
        else:
            glucose_section += f"Blood sugar returns toward normal after around {return_min} minutes."
    else:
        glucose_section += "Blood sugar has not fully returned to baseline within the simulation window."
    if _uses_gi_placeholder(output):
        from biosim_ode.gut_absorption import GI_SLOPE_PLACEHOLDER_NOTE
        glucose_section += f"\nNote: {GI_SLOPE_PLACEHOLDER_NOTE}"

    # Energy section
    energy_section = ls.energy_narrative
    if ls.has_energy_dip:
        energy_section += (
            f"\n\nNote: A 'sugar crash' or energy dip is predicted around "
            f"{round(ls.energy_dip_min)} minutes. This happens when blood sugar "
            f"falls quickly after a high peak."
        )

    # Risk section
    #
    # Regression: RiskIndicator.caveat (e.g. the "this reading may just
    # reflect blunted insulin secretion, not metabolic health" warning
    # assess_risk() attaches to mTOR/AMPK indicators -- see CLAUDE.md's
    # "mTOR/AMPK risk indicators could score a T2D patient healthier"
    # fix) was computed and attached to every indicator, but this loop
    # only ever read .description, never .caveat -- so the exact warning
    # that fix was meant to surface never actually reached the lay report
    # text a real user reads. Same false-reassurance bug class as that
    # fix, one level further down the pipeline.
    risk_lines = [f"Overall: {risk.overall.value.upper()} — {risk.summary}", ""]
    for indicator in risk.indicators:
        light_symbol = {"green": "[OK]", "amber": "[!]", "red": "[!!]", "grey": "[?]"}
        symbol = light_symbol.get(indicator.traffic_light.value, "[?]")
        risk_lines.append(f"{symbol} {indicator.label}: {indicator.description}")
        if indicator.caveat:
            risk_lines.append(f"    Note: {indicator.caveat}")
    # The runner adds this note when an input variant's effect size is
    # unverified (biosim_genotype UNVERIFIED_VARIANT_NOTES); repeat it here.
    from biosim_genotype.effect_db import UNVERIFIED_VARIANT_NOTES
    output_warnings = getattr(output, "warnings", None) or []
    for note in UNVERIFIED_VARIANT_NOTES.values():
        if note in output_warnings:
            risk_lines.append(f"Note: {note}")
    # Same for the runner's drug-validation notes (one per drug applied):
    # EXPERIMENTAL, or VALIDATED with its caveat.
    from biosim_physiology.drug_interactions.validation_status import (
        EXPERIMENTAL_DRUG_NOTES,
        VALIDATED_DRUG_NOTES,
    )
    for note in [*EXPERIMENTAL_DRUG_NOTES.values(), *VALIDATED_DRUG_NOTES.values()]:
        if note in output_warnings:
            risk_lines.append(f"Note: {note}")
    risk_section = "\n".join(risk_lines)

    # Comparison section
    comparison_section = None
    if comparison_output is not None:
        m2 = comparison_output.expert_metrics
        diff = (m.glucose_peak_mM - m2.glucose_peak_mM)
        direction = "higher" if diff > 0 else "lower"
        comparison_section = (
            f"Compared to '{comparison_output.scenario_name}': "
            f"this scenario produces a {abs(diff):.1f} mmol/L {direction} blood sugar peak."
        )

    disclaimer = (
        "IMPORTANT: This simulation is for research and educational purposes only. "
        "It is not a medical diagnosis or clinical prediction. "
        "Individual responses vary. Consult a healthcare professional for personal health advice."
    )

    return PlainLanguageReport(
        headline=headline,
        glucose_section=glucose_section,
        energy_section=energy_section,
        risk_section=risk_section,
        comparison_section=comparison_section,
        disclaimer=disclaimer,
        risk_profile=risk,
    )
