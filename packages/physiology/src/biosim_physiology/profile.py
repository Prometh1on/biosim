# SAVE TO: packages\physiology\src\biosim_physiology\profile.py
"""
PhysiologicalProfile -- non-genetic individual characterisation.

This is everything about a person that isn't in their DNA but
still affects the model parameters: age, sex, body composition,
fitness, disease states, time of day, and current medications.

Design decisions:
    - BiologicalSex is distinct from gender identity. The model uses
      biological sex to parameterise physiological differences in
      insulin sensitivity and body composition. The UI layer is
      responsible for handling this sensitively.
    - time_of_day_hours drives the circadian modifier. None means
      unknown / no circadian adjustment applied.
    - disease_states is a list because comorbidities are common.
      Multiple states compose multiplicatively.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum
from typing import Optional


class BiologicalSex(str, Enum):
    MALE          = "male"
    FEMALE        = "female"
    NOT_SPECIFIED = "not_specified"


class MenstrualCyclePhase(str, Enum):
    """
    Track 1 row 1.36 -- discrete, not continuous, deliberately: Valdes &
    Elkind-Hirsch 1991 (PMID:1997519) sampled exactly these three named
    phases of the same cycle via tolbutamide-modified IVGTT + minimal
    model, not a continuous cycle-day curve, so this enum mirrors the
    source's own three anchor points rather than inventing an
    interpolated "cycle day" field with no supporting data between them.
    See modifiers.py's own menstrual_cycle_modifiers() docstring for the
    full derivation and citation.
    """
    FOLLICULAR = "follicular"  # early follicular, the paper's own reference phase
    MIDCYCLE   = "midcycle"    # NOT significantly different from follicular per
                                # the source's own conclusion -- see that
                                # function's docstring for why this is neutral,
                                # not a smaller nonzero effect.
    LUTEAL     = "luteal"      # the paper's own sole significant finding


class FitnessLevel(str, Enum):
    SEDENTARY          = "sedentary"
    LIGHTLY_ACTIVE     = "lightly_active"
    MODERATELY_ACTIVE  = "moderately_active"  # reference
    VERY_ACTIVE        = "very_active"
    ATHLETE            = "athlete"


class DiseaseState(str, Enum):
    HEALTHY            = "healthy"
    INSULIN_RESISTANT  = "insulin_resistant"   # pre-diabetes / metabolic syndrome
    TYPE_2_DIABETES    = "type_2_diabetes"
    TYPE_1_DIABETES    = "type_1_diabetes"     # autoimmune beta-cell failure --
                                                # mechanistically distinct from
                                                # TYPE_2_DIABETES/INSULIN_RESISTANT
                                                # (secretion loss, not resistance)
    HYPOTHYROIDISM     = "hypothyroidism"
    HYPERTHYROIDISM    = "hyperthyroidism"
    GESTATIONAL_DIABETES = "gestational_diabetes"  # pregnancy-specific --
                                                     # distinct from pre-
                                                     # existing T1D/T2D
                                                     # ("pregestational
                                                     # diabetes"), not
                                                     # modelled here. Unlike
                                                     # every other state in
                                                     # this enum, its
                                                     # modifiers.py preset
                                                     # folds normal
                                                     # PREGNANCY's own
                                                     # insulin-resistance
                                                     # shift together with
                                                     # GDM's additional
                                                     # decompensation into
                                                     # ONE preset -- this
                                                     # platform has no
                                                     # separate "pregnant,
                                                     # not diabetic" state to
                                                     # layer it on top of.
                                                     # See that file's own
                                                     # citation comment
                                                     # before adding a
                                                     # "pregnant" state or
                                                     # changing this one.
    GCK_MODY = "gck_mody"  # glucokinase-gene MODY (MODY2) specifically --
                            # NOT generic "MODY"; other subtypes (e.g.
                            # HNF1A/MODY3) have a clinically DIFFERENT,
                            # more severe, progressive secretion-loss
                            # phenotype (see modifiers.py's own citation
                            # comment) and are NOT represented by this
                            # value. A benign glucose-sensing THRESHOLD
                            # reset, not a resistance or secretion-loss
                            # mechanism -- insulin sensitivity is
                            # genuinely unaffected, unlike every other
                            # diabetes-related state in this enum.
    CHRONIC_KIDNEY_DISEASE = "chronic_kidney_disease"  # moderate-severe,
                            # non-diabetic CKD (mean eGFR ~38 ml/min/1.73m2)
                            # -- a renal, not a beta-cell or classically
                            # "diabetic" disease state, composable with
                            # any of the states above (a common real
                            # comorbidity, e.g. T2D+CKD) unlike the
                            # mutually-exclusive pairs enforced below.
                            # Reduces insulin SENSITIVITY and CLEARANCE
                            # (unlike every state above, which touches at
                            # most one of those two), with an explicitly
                            # confirmed NULL effect on secretion -- see
                            # modifiers.py's own citation comment.
    PCOS = "pcos"           # Polycystic ovary syndrome (Track 1 row 1.29) --
                            # peripheral insulin resistance, NOT a beta-cell
                            # secretion-loss or hepatic-glucose-production
                            # state (both explicitly confirmed NULL by the
                            # same sources -- see modifiers.py's own citation
                            # comment). Composable with any other state above
                            # (a common real comorbidity, e.g. T2D+PCOS),
                            # same category as CHRONIC_KIDNEY_DISEASE, not
                            # mutually exclusive. A genuine, UNRESOLVED
                            # literature conflict exists on whether the
                            # resistance is obesity-INDEPENDENT or obesity-
                            # GATED -- disclosed in modifiers.py, not
                            # resolved by fiat; this platform's own preset
                            # uses the obese-cohort magnitude specifically,
                            # disclosed as such.
    ACROMEGALY = "acromegaly"  # Chronic growth-hormone excess (Track 1 row
                            # 1.37) -- Gb-ONLY, same partial-implementation
                            # shape as GCK_MODY/tirzepatide (no extractable
                            # insulin-sensitivity magnitude in the source,
                            # see modifiers.py's own citation comment).
                            # Almost certainly the SAME underlying mechanism
                            # as Drug.GROWTH_HORMONE (row 1.35, acute
                            # infusion) at a chronic exposure duration --
                            # deliberately non-overlapping BY CONSTRUCTION,
                            # not by mutual exclusion: GROWTH_HORMONE
                            # touches only p3 (acute), this touches only Gb
                            # (chronic), so composing both (a person with
                            # untreated acromegaly who also receives
                            # exogenous GH, clinically paradoxical but not
                            # modelled as impossible here) does not double-
                            # count any single number. Composable with any
                            # other disease state above (e.g. acromegaly-
                            # induced diabetes is a real, recognized
                            # clinical entity) -- not added to any mutual-
                            # exclusion pair.
    OSA = "osa"             # Obstructive sleep apnoea, as a CHRONIC clinical
                            # diagnosis (Track 1 row 1.39) -- Harsch et al.
                            # 2004 (PMID:14512265), a clamp-derived p3 factor
                            # from an AHI>20 clinical cohort BEFORE CPAP
                            # treatment. Mechanistically the SAME axis as
                            # `PhysiologicalProfile.desaturation_events_per_hour`
                            # (an ACUTE, wakefulness-only intermittent-hypoxia
                            # exposure rate, Louis & Punjabi 2009,
                            # PMID:19265062) -- both represent one real
                            # person's oxyhemoglobin-desaturation burden via
                            # two independent p3-lowering routes. Composing
                            # BOTH would silently double-count (0.847*0.829
                            # ~= 0.702 on p3, from stacking one real
                            # mechanism's two representations) -- this is
                            # the first disease-state x profile-field
                            # (rather than disease-state x disease-state)
                            # mutual exclusion this platform enforces; see
                            # __post_init__ below. Composable with every
                            # other disease state (OSA-related insulin
                            # resistance is a real, common T2D/obesity
                            # comorbidity, not modelled as exclusive of
                            # anything else) -- confounds with obesity are
                            # disclosed, not resolved, matching PCOS's own
                            # precedent above.
    NAFLD = "nafld"         # Non-alcoholic fatty liver disease, biopsy-proven,
                            # non-diabetic (Track 1 row 1.45) -- Bugianesi et
                            # al. 2005 (PMID:15747110), a two-step euglycaemic
                            # clamp finding peripheral glucose disposal
                            # reduced 30%/45% at low/high insulin; measured
                            # endogenous glucose production was NORMAL, so
                            # only p3 is modelled (hepatic insulin resistance
                            # reported in the full text is NOT represented --
                            # see modifiers.py). 0.70 (the low-dose step) is a
                            # disclosed judgment call, not a conservative
                            # bound. Composable with every other state, but
                            # with INSULIN_RESISTANT/TYPE_2_DIABETES it is
                            # probably partial double-counting of one
                            # underlying metabolic-syndrome state -- disclosed
                            # in modifiers.py, not blocked (PCOS/OSA precedent).


@dataclass
class PhysiologicalProfile:
    """
    Non-genetic physiological characterisation of an individual.

    All fields have sensible population-average defaults so that
    a profile can be constructed incrementally, adding only known
    attributes.
    """
    # Generous plausibility GUARD for derived BMI, not a physiological
    # claim or a cited reference range -- it exists to reject garbage
    # (NaN/underflow/overflow/absurd pairs) loudly at construction.
    # Unannotated, so NOT dataclass fields.
    MIN_PLAUSIBLE_BMI = 5.0
    MAX_PLAUSIBLE_BMI = 250.0

    age_years:       float          = 30.0
    biological_sex:  BiologicalSex  = BiologicalSex.NOT_SPECIFIED
    body_weight_kg:  float          = 70.0
    height_cm:       float          = 170.0
    body_fat_percent: Optional[float] = None
    fitness_level:   FitnessLevel   = FitnessLevel.MODERATELY_ACTIVE
    disease_states:  list[DiseaseState] = field(
        default_factory=lambda: [DiseaseState.HEALTHY]
    )
    time_of_day_hours: Optional[float] = None  # 0–24; None = no circadian adjustment
    sleep_hours:       Optional[float] = None  # prior night's sleep duration, 0-24;
                                                # None = no sleep-restriction adjustment.
                                                # See modifiers.py's own sleep_modifiers()
                                                # docstring for the citation and derivation
                                                # -- a genuinely different axis from
                                                # time_of_day_hours (circadian PHASE) --
                                                # this is sleep DURATION/DEBT.
    desaturation_events_per_hour: Optional[float] = None  # acute intermittent
                                                # hypoxia exposure rate (oxyhemoglobin
                                                # desaturation events/hour during
                                                # wakefulness), 0+; None = no adjustment.
                                                # See modifiers.py's own
                                                # intermittent_hypoxia_modifiers()
                                                # docstring for the citation/derivation --
                                                # a genuinely different axis from
                                                # sleep_hours (this is hypoxemia during
                                                # wakefulness, not sleep duration/debt;
                                                # the source study excludes sleep loss
                                                # entirely) and from stress_level
                                                # (cortisol was directly measured
                                                # unchanged in the source study).
    menstrual_cycle_phase: Optional[MenstrualCyclePhase] = None  # None = no
                                                # adjustment (matches every other
                                                # Optional axis field's own
                                                # convention). This platform's first
                                                # female-specific PHYSIOLOGICAL axis
                                                # (distinct from the disease states
                                                # GESTATIONAL_DIABETES/PCOS, which
                                                # represent pathology, not normal
                                                # cycle physiology). See
                                                # modifiers.py's own
                                                # menstrual_cycle_modifiers()
                                                # docstring for the citation.
    cold_acclimatized: bool = False             # Track 1 row 1.40: has
                                                # completed a ~10-day CHRONIC
                                                # cold-acclimation protocol
                                                # (14-15C), Hanssen et al.
                                                # 2015 (PMID:26147760).
                                                # Deliberately narrow: NOT a
                                                # general "ambient
                                                # temperature" axis, and NOT
                                                # the same as a brief acute
                                                # cold exposure (Sellers et
                                                # al. 2021, PMID:33090911,
                                                # a SINGLE 1h shivering bout
                                                # in the SAME research
                                                # group's own follow-up work,
                                                # found the OPPOSITE
                                                # direction -- transient
                                                # post-shivering insulin
                                                # resistance, not
                                                # sensitisation). This field
                                                # represents ONLY the
                                                # chronic-acclimation regime
                                                # Hanssen measured -- see
                                                # modifiers.py's own
                                                # cold_acclimation_modifiers()
                                                # docstring for the full
                                                # derivation and the T2D-only
                                                # evidence-base disclosure.

    def __post_init__(self) -> None:
        if not (0 <= self.age_years <= 120):
            raise ValueError(f"age_years must be 0–120, got {self.age_years}")
        # NaN/inf passed the old `<= 0` checks (NaN compares False), and a
        # tiny-but-positive height (REST allowed any gt=0) underflowed
        # bmi's denominator -> ZeroDivisionError inside personalize(),
        # which the runner's catch-all swallowed, silently dropping the
        # WHOLE profile (a requested disease state simulated as healthy,
        # HTTP 200). Reject at construction instead -- REST turns this
        # ValueError into a 400. BMI is bounded to what is physically
        # conceivable (a deliberately generous guard, not a cited range),
        # so absurd weight/height pairs fail loudly rather than landing in
        # the "severely obese" or "underweight" branch as garbage.
        if not math.isfinite(self.body_weight_kg) or self.body_weight_kg <= 0:
            raise ValueError(f"body_weight_kg must be finite and positive, got {self.body_weight_kg}")
        if not math.isfinite(self.height_cm) or self.height_cm <= 0:
            raise ValueError(f"height_cm must be finite and positive, got {self.height_cm}")
        try:
            bmi = self.body_weight_kg / (self.height_cm / 100.0) ** 2
        except ZeroDivisionError:      # height so small its square underflows to 0
            bmi = float("inf")
        except OverflowError:          # height so large its square overflows
            bmi = 0.0
        if not (self.MIN_PLAUSIBLE_BMI <= bmi <= self.MAX_PLAUSIBLE_BMI):
            raise ValueError(
                f"body_weight_kg={self.body_weight_kg} / height_cm={self.height_cm} "
                f"gives an implausible BMI ({bmi:.4g}); must be between "
                f"{self.MIN_PLAUSIBLE_BMI} and {self.MAX_PLAUSIBLE_BMI}"
            )
        if self.body_fat_percent is not None:
            if not (0 < self.body_fat_percent < 100):
                raise ValueError(f"body_fat_percent must be 0–100, got {self.body_fat_percent}")
        if self.time_of_day_hours is not None:
            if not (0 <= self.time_of_day_hours < 24):
                raise ValueError(f"time_of_day_hours must be 0–24, got {self.time_of_day_hours}")
        if self.sleep_hours is not None:
            if not (0 <= self.sleep_hours <= 24):
                raise ValueError(f"sleep_hours must be 0–24, got {self.sleep_hours}")
        if self.desaturation_events_per_hour is not None:
            if not math.isfinite(self.desaturation_events_per_hour) or \
                    self.desaturation_events_per_hour < 0:
                raise ValueError(
                    "desaturation_events_per_hour must be finite and non-negative, got "
                    f"{self.desaturation_events_per_hour}"
                )
        # Found during a security-reviewer pass on the HYPERTHYROIDISM
        # addition: no field here validated disease_states at all, so a
        # duplicate entry silently double-applies that state's modifier
        # (e.g. [TYPE_2_DIABETES, TYPE_2_DIABETES] squares its Gb/p3/
        # secretion_scale factors -- a 2x-compounded, undocumented
        # severity level), and mutually-exclusive real-world conditions
        # (a patient cannot be both hypo- and hyperthyroid) could be
        # combined silently -- worse, since HYPERTHYROIDISM's static
        # modifier is deliberately {} (see modifiers.py), setting both
        # together is indistinguishable from HYPOTHYROIDISM alone at the
        # static-modifier layer while STILL composing
        # EGPParameters.for_hyperthyroidism()'s ODE pathway in
        # SimulationRunner._detect_hepatic(), producing a physiologically
        # nonsensical combination of both conditions' effects at once.
        # Fail fast here, matching DrugEvent.__post_init__'s existing
        # validate-and-reject-at-construction convention (scenario.py).
        if len(self.disease_states) != len(set(self.disease_states)):
            raise ValueError(f"disease_states contains duplicates: {self.disease_states}")
        if DiseaseState.HYPOTHYROIDISM in self.disease_states \
                and DiseaseState.HYPERTHYROIDISM in self.disease_states:
            raise ValueError(
                "disease_states cannot contain both HYPOTHYROIDISM and "
                "HYPERTHYROIDISM -- these are mutually exclusive conditions"
            )
        # TYPE_1_DIABETES (near-total secretion loss) and
        # TYPE_2_DIABETES/INSULIN_RESISTANT (resistance-driven, with
        # partially preserved -- often compensatorily elevated --
        # secretion) represent opposite beta-cell mechanisms. Combining
        # them would silently multiply TYPE_1_DIABETES's near-zero
        # secretion_scale/Ib factors by TYPE_2_DIABETES's/
        # INSULIN_RESISTANT's own factors -- a physiologically
        # nonsensical "both at once" state ("double diabetes" is a real
        # but rare clinical entity this platform does not attempt to
        # model). Same fail-fast convention as the thyroid check above.
        if DiseaseState.TYPE_1_DIABETES in self.disease_states and (
            DiseaseState.TYPE_2_DIABETES in self.disease_states
            or DiseaseState.INSULIN_RESISTANT in self.disease_states
        ):
            raise ValueError(
                "disease_states cannot contain TYPE_1_DIABETES together with "
                "TYPE_2_DIABETES or INSULIN_RESISTANT -- these represent "
                "opposite beta-cell mechanisms (secretion loss vs. "
                "resistance) this platform does not model in combination"
            )
        # GESTATIONAL_DIABETES is, by clinical definition, diabetes FIRST
        # recognized during pregnancy in a person WITHOUT pre-existing
        # diabetes -- a person already diagnosed with TYPE_1_DIABETES or
        # TYPE_2_DIABETES who becomes pregnant has "pregestational
        # diabetes" instead, a distinct clinical entity this platform
        # does not attempt to model (same "not modelling every real but
        # rare combination" scope decision as TYPE_1_DIABETES's own
        # exclusion above). INSULIN_RESISTANT is deliberately NOT
        # excluded here -- pre-existing insulin resistance/pre-diabetes
        # is a real, common GDM risk factor, not a contradictory
        # diagnosis, so composing them is unusual but not nonsensical.
        if DiseaseState.GESTATIONAL_DIABETES in self.disease_states and (
            DiseaseState.TYPE_1_DIABETES in self.disease_states
            or DiseaseState.TYPE_2_DIABETES in self.disease_states
        ):
            raise ValueError(
                "disease_states cannot contain GESTATIONAL_DIABETES together "
                "with TYPE_1_DIABETES or TYPE_2_DIABETES -- pre-existing "
                "diabetes in pregnancy is a distinct clinical entity "
                "(\"pregestational diabetes\") this platform does not model"
            )
        # GCK_MODY is a distinct monogenic diagnosis -- a person confirmed
        # to have TYPE_1_DIABETES or TYPE_2_DIABETES by definition does not
        # also carry the GCK-MODY diagnosis for the same episode (same
        # "distinct clinical entity" scope decision as GESTATIONAL_DIABETES
        # above). INSULIN_RESISTANT and GESTATIONAL_DIABETES are NOT
        # excluded -- a GCK-MODY carrier's glucose-threshold reset is an
        # independent mechanism from either, so composing them (unusual,
        # but not physiologically contradictory) is allowed.
        if DiseaseState.GCK_MODY in self.disease_states and (
            DiseaseState.TYPE_1_DIABETES in self.disease_states
            or DiseaseState.TYPE_2_DIABETES in self.disease_states
        ):
            raise ValueError(
                "disease_states cannot contain GCK_MODY together with "
                "TYPE_1_DIABETES or TYPE_2_DIABETES -- these are distinct "
                "clinical diagnoses this platform does not model in "
                "combination"
            )
        # Track 1 row 1.39: the first disease-state x profile-field mutual
        # exclusion this platform enforces (every check above is disease-
        # state-vs-disease-state). DiseaseState.OSA (a chronic clinical
        # diagnosis, Harsch 2004 clamp cohort) and
        # desaturation_events_per_hour (an acute wakefulness-only exposure
        # rate, Louis & Punjabi 2009) both lower p3 via the SAME real
        # oxyhemoglobin-desaturation mechanism in one person -- composing
        # both would silently stack two representations of it (see OSA's
        # own enum-comment derivation: 0.847*0.829~=0.702, not disclosed
        # anywhere as a combined effect). Fail fast, same convention as
        # every check above.
        #
        # Code-reviewer finding: excludes exactly 0.0, not just None --
        # intermittent_hypoxia_modifiers() itself treats 0 events/hour as
        # a neutral no-op (`if events_per_hour == 0: return {}`), so
        # composing OSA with a MEASURED rate of exactly 0.0 double-counts
        # nothing at all; rejecting it would be a false-positive
        # restriction on a legitimate input (a diagnosed OSA patient whose
        # sleep study happened to record zero desaturation events, a
        # different fact from "not measured"/None).
        if DiseaseState.OSA in self.disease_states \
                and self.desaturation_events_per_hour not in (None, 0.0):
            raise ValueError(
                "disease_states cannot contain OSA together with "
                "desaturation_events_per_hour set -- these represent the "
                "same oxyhemoglobin-desaturation mechanism via two "
                "independent p3-lowering routes and would silently "
                "double-count if composed together. Use OSA for a "
                "diagnosed chronic condition, or desaturation_events_per_hour "
                "for a specific measured acute exposure rate -- not both."
            )

    # ------------------------------------------------------------------
    # Derived properties
    # ------------------------------------------------------------------

    @property
    def bmi(self) -> float:
        """Body mass index [kg/m²]."""
        return self.body_weight_kg / (self.height_cm / 100.0) ** 2

    @property
    def bmi_category(self) -> str:
        bmi = self.bmi
        if bmi < 18.5:  return "underweight"
        if bmi < 25.0:  return "normal"
        if bmi < 30.0:  return "overweight"
        if bmi < 35.0:  return "obese"
        return "severely_obese"

    @property
    def lean_mass_kg(self) -> Optional[float]:
        if self.body_fat_percent is not None:
            return self.body_weight_kg * (1.0 - self.body_fat_percent / 100.0)
        return None

    @property
    def is_diabetic(self) -> bool:
        return DiseaseState.TYPE_2_DIABETES in self.disease_states

    @property
    def is_insulin_resistant(self) -> bool:
        return (DiseaseState.INSULIN_RESISTANT in self.disease_states
                or DiseaseState.TYPE_2_DIABETES in self.disease_states)

    def with_time(self, hour: float) -> "PhysiologicalProfile":
        """Return a new profile with time_of_day_hours set."""
        import dataclasses
        return dataclasses.replace(self, time_of_day_hours=hour)
