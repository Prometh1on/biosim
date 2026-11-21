# SAVE TO: packages\physiology\src\biosim_physiology\exercise.py
"""
Exercise simulation module.

Models the metabolic effects of physical exercise on glucose/insulin dynamics.

Two time windows:
  POST-exercise: GLUT4 elevated, insulin sensitivity enhanced 20-50%
  DURING exercise: muscle contraction activates AMPK-independent glucose uptake

Literature:
    Richter & Hargreaves 2013, Physiol Rev. PMID:23899560
    Dela & Kjaer 2006, Essays Biochem. PMID:17144881
    Mikines et al. 1988, Am J Physiol. PMID:3126668 (post-exercise p3
        boost coefficient -- see post_exercise_modifiers()'s own docstring)
    Wang et al. 2013, J Appl Physiol. PMID:23493357 (post-exercise p3
        boost coefficient, corroborating)
    Hayashi et al. 2005, J Clin Endocrinol Metab. PMID:15840752;
    Farkas et al. 2022, Disabil Rehabil. PMID:33905292;
    Bordenave et al. 2008, Diabetes Metab. PMID:18448376
        (post-exercise p1_boost coefficient -- genuinely divergent
        evidence, NOT a derived point estimate; see
        post_exercise_modifiers()'s own docstring)
"""

from __future__ import annotations
import math
from dataclasses import dataclass
from enum import Enum
from typing import Optional


class ExerciseIntensity(str, Enum):
    LIGHT    = "light"
    MODERATE = "moderate"
    VIGOROUS = "vigorous"
    MAXIMAL  = "maximal"

    @property
    def percent_vo2max(self) -> float:
        return {"light": 30.0, "moderate": 50.0, "vigorous": 70.0, "maximal": 90.0}[self.value]


# MODEL ASSUMPTION (no source): extra insulin-independent glucose uptake by
# working muscle, min^-1 per unit of relative intensity. Equals the old
# during-bout increment 0.018 x 4 (B7, 2026-11-16), kept unchanged so that
# removing the resting double count does not alter exercise behaviour.
EXERCISE_GLUCOSE_UPTAKE_PER_RELATIVE_INTENSITY = 0.072


@dataclass(frozen=True)
class ExerciseBout:
    """A single exercise session."""
    intensity_percent_vo2max: float
    duration_minutes:         float
    hours_before_simulation:  float = 0.0
    label:                    str   = ""

    def __post_init__(self) -> None:
        if not (0 < self.intensity_percent_vo2max <= 100):
            raise ValueError(f"Intensity must be 0-100%, got {self.intensity_percent_vo2max}")
        if self.duration_minutes <= 0:
            raise ValueError(f"Duration must be positive, got {self.duration_minutes}")
        if self.hours_before_simulation < 0:
            raise ValueError(f"hours_before_simulation must be >= 0")

    @classmethod
    def from_intensity(
        cls, intensity: ExerciseIntensity, duration_minutes: float, hours_before: float = 0.0
    ) -> "ExerciseBout":
        return cls(
            intensity_percent_vo2max=intensity.percent_vo2max,
            duration_minutes=duration_minutes,
            hours_before_simulation=hours_before,
            label=f"{intensity.value} exercise {duration_minutes:.0f}min",
        )

    @property
    def mets(self) -> float:
        return self.intensity_percent_vo2max / 10.0

    @property
    def relative_intensity(self) -> float:
        return self.intensity_percent_vo2max / 100.0

    @property
    def duration_hours(self) -> float:
        return self.duration_minutes / 60.0

    def post_exercise_modifiers(self) -> dict[str, float]:
        """
        Parameter multipliers for post-exercise insulin sensitisation.
        Effect decays exponentially (t½ ~12h).

        `max_p3_boost`'s coefficient (2026-08-12 deep-search plan, item
        B1): anchored to two independent, live-verified real clamp
        studies rather than the prior uncited 0.50 -- Mikines KJ, Sonne
        B, Farrell PA, Tronier B, Galbo H. "Effect of physical exercise
        on sensitivity and responsiveness to insulin in humans." Am J
        Physiol. 1988. PMID:3126668 (n=7 healthy untrained men, 60min
        150W cycling, euglycaemic hyperinsulinaemic clamp immediately
        post-exercise: Km 52->43, Vmax 9.5->10.9 mg/min/kg) and Wang X
        et al. "A ~60-min brisk walk increases insulin-stimulated
        glucose disposal..." J Appl Physiol. 2013. PMID:23493357 (n=10
        nonobese older women, ~60min brisk walk, 2-stage hyperinsulinaemic
        clamp: Rd 78.4->88.0 umol/kg FFM/min).

        Both studies report BOTH a Vmax/Rd-style responsiveness shift
        AND (Mikines only) a Km sensitivity shift -- using either alone
        discards half of the same measurement. Applying one consistent
        rule to both (ratio of insulin-STIMULATED disposal above basal
        insulin, evaluated at this platform's own validated OGTT peak
        insulin, 354 pM) gives +21.8% (Mikines) and +16.0% (Wang) --
        genuinely convergent given the different sex/age/modality, not
        forced into agreement. Converting each to an implied `volume`
        requires an explicit, DISCLOSED, NOT independently literature-
        verified assumption about each study's own %VO2max (Mikines'
        150W absolute workload was assumed ~70% VO2max via the standard
        ACSM cycling equation for an average subject; Wang's own
        "brisk walk" was assumed ~50% VO2max, since neither paper states
        an exact intensity) -- giving implied coefficients of 0.311 and
        0.321 respectively, averaged to 0.31. Honest range across
        plausible intensity assumptions: 0.29-0.36. The functional form
        itself (linear in `volume`) remains unanchored -- neither study
        varied volume, so this coefficient is only directly supported
        near volume~0.5-0.7, not across its full 0-2.0 range.

        Mikines' subjects were UNTRAINED -- this coefficient represents
        the ACUTE effect for an untrained person specifically; it is
        composed multiplicatively alongside `fitness_modifiers()`'s own
        separate, CHRONIC-training p3 multiplier (a real, distinct
        mechanism, not double-counting the same effect), but the acute
        magnitude itself is not adjusted for the person's own fitness
        level (a known, undisclosed-until-now model limitation, not
        fixed this pass).

        `p1_boost`'s 0.10 coefficient and `gb_suppression`'s 0.03
        coefficient were previously uncited (fresh Track 1 scoping,
        2026-08-18) -- a literature search was run for both, live-
        verified (PubMedFetcher, before this comment was written).
        **NOT upgraded to derived point estimates for either, and this
        is a deliberate decision, not an oversight**:

        p1 (glucose effectiveness, SG in Bergman minimal-model terms) --
        real acute-exercise studies exist, but the evidence is genuinely
        DIVERGENT, not convergent like max_p3_boost's own two studies:
        Hayashi Y et al. "A single bout of exercise at higher intensity
        enhances glucose effectiveness in sedentary men." J Clin
        Endocrinol Metab. 2005;90(7):4327-33. PMID:15840752 -- n=6
        healthy men, single 20min cycle bout: a DIFFERENT two-
        compartment SG(2)* metric (dl/kg/min, not this model's SG
        min^-1) rose +237.1+/-50.5% at 70% VO2max (P<0.05), not
        significant at 50%. But TWO studies using the SAME minimal-
        model SG convention (min^-1) as this platform's own p1 are
        NULL in exactly the populations that matter most: Farkas GJ et
        al. Disabil Rehabil. 2022;44(2):343-9. PMID:33905292 -- in the
        able-bodied control group (the population closest to this
        platform's own healthy baseline), "Sg was not different at all
        timepoints (p>0.05)" (only the paraplegic group showed a
        transient +23.8% rise, gone by 24h); and Bordenave S et al.
        Diabetes Metab. 2008;34(1):50-7. PMID:18448376 -- in T2D
        patients, "S(G) remained unchanged" despite a large (+773%)
        S(I) rise in the SAME cohort. Forcing a single derived
        coefficient from Hayashi's differently-defined metric while
        ignoring two metric-matched null results would manufacture
        false precision, the same failure mode this project has
        already rejected elsewhere (KG_secretion~200, joint-fit
        p1~0.008). Direction is plausible; magnitude is not derivable
        from convergent evidence. 0.10 remains an undisclosed-precision
        placeholder, now at least backed by real (if divergent)
        citations rather than a bare guess.

        Gb (fasting glucose suppression) -- genuinely unresolved, not
        just under-searched: two independent search angles (direct
        SG/Gb terms, and fasting-glucose-recovery-timing terms) found
        no usable quantitative anchor. The one design-matched candidate
        -- Nakayama A et al. "Effect of downhill walking on next-day
        muscle damage and glucose metabolism in healthy young
        subjects." J Physiol Sci. 2019;69(4):511-9. PMID:29679309,
        single 30min bout, healthy young adults, next-day fasting
        glucose lower after the uphill trial (P<0.01) -- reports no
        extractable number in its abstract, and its full text is
        confirmed blocked (Europe PMC: "The publisher of this article
        does not allow downloading of the full text in XML form").
        Wang et al. 2013 (already cited above, PMID:23493357, same
        cohort as max_p3_boost) reports a null on basal glucose Ra (a
        flux, not the concentration Gb represents) in the same design
        window. Borer KT et al. Med Sci Sports Exerc.
        2009;41(8):1606-14. PMID:19568199 has a real number (-20.2%)
        but a mismatched design (two 2h bouts, meal-timing-dependent).
        A T2D meta-analysis (MacLeod SF et al. Diabetes Metab Res Rev.
        2013;29(8):593-603. PMID:24038928) found fasting glucose
        non-significantly lower after exercise (P=0.13). 0.03 remains
        an undisclosed-precision placeholder with NO real citation
        behind it -- an open literature gap, not a resolved one. Only
        worth revisiting if Nakayama's full text becomes accessible, or
        a new single-bout, healthy-population, quantified fasting-
        glucose study is published.
        """
        volume = min(self.relative_intensity * self.duration_hours, 2.0)
        max_p3_boost = 1.0 + 0.31 * volume
        decay = math.exp(-math.log(2) * self.hours_before_simulation / 12.0)
        p3_multiplier = 1.0 + (max_p3_boost - 1.0) * decay
        p1_boost      = 1.0 + 0.10 * volume * decay
        gb_suppression = 1.0 - 0.03 * volume * decay
        return {"p3": p3_multiplier, "p1": p1_boost, "Gb": gb_suppression}

    def during_exercise_modifiers(self) -> dict[str, float]:
        """
        Acute exercise state: extra insulin-independent glucose uptake by
        working muscle (`glucose_uptake_delta`, min^-1, ADDED to the
        clearance rate) and a secretion multiplier.

        Until 2026-11-16 (B7) the uptake was a multiplier (1 + 4 x intensity)
        on a resting `extra_g_clearance` of 0.018/min that double-counted p1.
        The resting term is now zero; the uptake keeps exactly its old
        size, 0.018 x 4 x intensity, see the constant below.
        """
        uptake = EXERCISE_GLUCOSE_UPTAKE_PER_RELATIVE_INTENSITY * self.relative_intensity
        secretion_suppress = 1.0 - 0.30 * self.relative_intensity
        return {"glucose_uptake_delta": uptake, "secretion_scale": secretion_suppress}


def light_walk(duration_minutes: float = 30.0, hours_before: float = 0.0) -> ExerciseBout:
    return ExerciseBout.from_intensity(ExerciseIntensity.LIGHT, duration_minutes, hours_before)


def moderate_run(duration_minutes: float = 45.0, hours_before: float = 0.0) -> ExerciseBout:
    return ExerciseBout.from_intensity(ExerciseIntensity.MODERATE, duration_minutes, hours_before)


def hiit_session(duration_minutes: float = 30.0, hours_before: float = 0.0) -> ExerciseBout:
    return ExerciseBout.from_intensity(ExerciseIntensity.VIGOROUS, duration_minutes, hours_before)
