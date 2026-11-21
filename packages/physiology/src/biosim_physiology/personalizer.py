# SAVE TO: packages\physiology\src\biosim_physiology\personalizer.py
"""
PhysiologicalPersonalizer -- applies all non-genetic modifiers to produce
a personalised PathwayParameters instance.

Application order:
    1. Age modifiers
    2. Biological sex modifiers
    3. BMI/adiposity modifiers
    4. Fitness modifiers
    5. Disease state modifiers
    6. Circadian (time-of-day) modifiers
    7. Drug interaction modifiers

All modifiers compose multiplicatively. The resulting PathwayParameters
is clamped to physiologically valid ranges before being returned.
"""

from __future__ import annotations

import dataclasses
import math
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Optional

from biosim_ode import PathwayParameters

from .chronobiology import circadian_modifiers
from .drug_interactions import (
    DrugDose,
    SULFONYLUREA_DRUGS,
    validate_drug_combination,
    drug_modifiers,
    glucose_independent_secretion_pM_per_min,
    modifiers_for_dose,
)
from .modifiers import (
    age_modifiers,
    bmi_modifiers,
    cold_acclimation_modifiers,
    compose_all_modifiers,
    disease_modifiers,
    fitness_modifiers,
    intermittent_hypoxia_modifiers,
    menstrual_cycle_modifiers,
    sex_modifiers,
    sleep_modifiers,
)
from .profile import DiseaseState, PhysiologicalProfile

if TYPE_CHECKING:
    from .exercise import ExerciseBout


@dataclass
class AppliedModifier:
    """Records one modifier that was applied and its effect."""
    source:          str            # e.g. "age_60", "bmi_32", "metformin_500mg"
    # Regression: this used to be a single parameter/factor pair, with
    # every category's parameter hardcoded to "p3" (or "p3/secretion_scale"
    # for age, "multiple" for disease) regardless of which parameter(s)
    # that category actually modifies -- BMI/disease/metformin, for
    # instance, also move Gb/secretion_scale, and factor=net_factors
    # .get("p3", 1.0) would silently read back the neutral default for a
    # source that changed a different parameter entirely, misreporting
    # it as having no effect. A per-source category can genuinely modify
    # more than one parameter (see modifiers.py/drug_interactions.py),
    # so this is now the category's own individual modifier dict
    # (parameter name -> that source's OWN factor, not the net combined
    # factor across all sources), not a single parameter/factor pair.
    parameters:      dict[str, float]
    description:     str


@dataclass
class PhysiologicalPersonalizationResult:
    """
    Result of non-genetic physiological personalisation.
    """
    original_params:      PathwayParameters
    personalized_params:  PathwayParameters
    profile:              PhysiologicalProfile
    applied_modifiers:    list[AppliedModifier]
    drug_doses:           list[DrugDose]
    combined_factors:     dict[str, float]   # net per-parameter factor
    # Populated when a safety clamp (e.g. MIN_GB_THRESHOLD_MARGIN_MM
    # below) actually changed the composed value from what the raw
    # drug/disease modifiers alone would have produced -- added
    # 2026-08-17 (MASTER_PLAN.md 1.23) after a /verify pass found this
    # was previously silent: three tirzepatide doses (5/10/15mg, whose
    # own tirzepatide_gb_factor() is genuinely dose-dependent) produced
    # byte-identical simulated output because this clamp floored every
    # one of them at the same value, with nothing anywhere signalling
    # that the calibrated dose-response never reached the output.
    clamp_warnings:       list[str] = field(default_factory=list)

    @property
    def n_parameters_modified(self) -> int:
        return sum(1 for f in self.combined_factors.values() if abs(f - 1.0) > 0.001)

    @property
    def largest_effect_param(self) -> Optional[str]:
        if not self.combined_factors:
            return None
        return max(self.combined_factors, key=lambda k: abs(self.combined_factors[k] - 1.0))

    def summary(self) -> dict:
        return {
            "age":             self.profile.age_years,
            "sex":             self.profile.biological_sex.value,
            "bmi":             round(self.profile.bmi, 1),
            "fitness":         self.profile.fitness_level.value,
            "disease_states":  [d.value for d in self.profile.disease_states],
            "time_of_day":     self.profile.time_of_day_hours,
            "n_drugs":         len(self.drug_doses),
            "n_modified":      self.n_parameters_modified,
            "net_factors":     {k: round(v, 4) for k, v in self.combined_factors.items()
                                if abs(v - 1.0) > 0.001},
        }


# Minimum and maximum allowed values for modified parameters
_PARAM_BOUNDS: dict[str, tuple[float, float]] = {
    "Gb":              (2.0,  15.0),    # mM fasting glucose
    "Ib":              (10.0, 500.0),   # pM fasting insulin
    "p1":              (0.001, 0.20),   # glucose effectiveness
    "p2":              (0.001, 0.10),
    "p3":              (1e-8,  1e-4),
    "ke":              (0.01,  0.50),
    "secretion_scale": (0.05,  3.0),
    "G_threshold":     (2.0,   6.0),    # mM
    "IGF1b":           (1.0,   50.0),   # nM
}

# Parameters that ALSO have their own point-of-use clamp inside
# GlucoseInsulinPathway.dydt() (pathway.py's P3_EFF/GB_EFF/P1_EFF/
# SECRETION_SCALE_EFF_FLOOR/CEILING), on top of the static clamp this
# module applies before integration starts (Section 8 housekeeping,
# found during the CORTISOL_P3_SENSITIVITY security review). Whenever
# an ODE-coupled pathway (CortisolPathway, NEFAPathway, ...) ALSO
# modifies one of these via a rule-11 signal (p3_modifier, gb_modifier,
# p1_modifier, secretion_scale_modifier), the value this module clamps
# below is NOT necessarily the value the ODE solver actually integrates
# with -- the coupled signal is re-applied at every solver step, after
# this static clamp has already run. The warning messages below name
# this explicitly for these four parameters so "clamped it to X" isn't
# read as a permanent guarantee about the simulated value.
_ALSO_CLAMPED_AT_ODE_POINT_OF_USE = frozenset({"Gb", "p1", "p3", "secretion_scale"})


def _ode_coupling_caveat(param: str) -> str:
    if param not in _ALSO_CLAMPED_AT_ODE_POINT_OF_USE:
        return ""
    return (
        f" Note: {param} also has its own point-of-use clamp inside the ODE "
        "solver (pathway.py) -- if an ODE-coupled pathway (e.g. CortisolPathway, "
        "NEFAPathway) further modifies this parameter during simulation, the "
        "value actually integrated may differ from the one clamped here."
    )


class PhysiologicalPersonalizer:
    """
    Applies all non-genetic physiological modifiers to PathwayParameters.

    Usage:
        profile = PhysiologicalProfile(
            age_years=55,
            biological_sex=BiologicalSex.FEMALE,
            bmi=28.5,
            fitness_level=FitnessLevel.LIGHTLY_ACTIVE,
            disease_states=[DiseaseState.INSULIN_RESISTANT],
            time_of_day_hours=20.0,
        )
        result = PhysiologicalPersonalizer().personalize(profile)
        pathway = GlucoseInsulinPathway(params=result.personalized_params)
    """

    def personalize(
        self,
        profile: PhysiologicalProfile,
        base_params: Optional[PathwayParameters] = None,
        drug_doses: Optional[list[DrugDose]] = None,
        exercise_bout: Optional["ExerciseBout"] = None,
    ) -> PhysiologicalPersonalizationResult:
        """
        Apply physiological modifiers to produce personalised parameters.

        Args:
            profile:       Individual's physiological characteristics
            base_params:   Starting parameters (default = population average)
            drug_doses:    Current medications/substances
            exercise_bout: Prior exercise session, if any -- applies
                           post_exercise_modifiers() (p1/p3/Gb multipliers
                           for the residual insulin-sensitising effect,
                           decaying with hours_before_simulation). This is
                           layer 8 of the 9 personalisation layers listed
                           in CLAUDE.md's PERSONALISATION LAYERS section.

        Returns:
            PhysiologicalPersonalizationResult with all modifications documented
        """
        base_params = base_params or PathwayParameters()
        drug_doses  = drug_doses or []

        # Collect all modifier dicts
        all_modifier_dicts = [
            age_modifiers(profile.age_years),
            sex_modifiers(profile.biological_sex),
            bmi_modifiers(profile.bmi),
            fitness_modifiers(profile.fitness_level),
            disease_modifiers(profile.disease_states),
            circadian_modifiers(profile.time_of_day_hours),
            sleep_modifiers(profile.sleep_hours),
            intermittent_hypoxia_modifiers(profile.desaturation_events_per_hour),
            menstrual_cycle_modifiers(profile.menstrual_cycle_phase),
            cold_acclimation_modifiers(profile.cold_acclimatized),
            drug_modifiers(drug_doses),
        ]
        if exercise_bout is not None:
            all_modifier_dicts.append(exercise_bout.post_exercise_modifiers())

        # Compose into net per-parameter factors
        net_factors = compose_all_modifiers(all_modifier_dicts)

        # Apply to base params
        base_dict = dataclasses.asdict(base_params)
        overrides: dict[str, float] = {}

        # D1 Composition-Boundary Behavioral Audit (2026-08-18) found that
        # every _PARAM_BOUNDS-clamped field EXCEPT Gb/G_threshold (which
        # already have their own dedicated joint clamp+warning below, from
        # row 1.23) silently clamped with no warning at all. At the time
        # of that audit, the p3 floor was not yet breached by any then-
        # shippable mechanism combination (closest observed ~5.7x above
        # it) -- but this was a snapshot, not a permanent guarantee: a
        # mandatory security-reviewer pass on row 1.29 (PCOS) found that
        # after several more p3-reducing mechanisms landed the same
        # session (indinavir, hydrochlorothiazide, growth hormone,
        # atorvastatin, menstrual-cycle luteal phase, PCOS itself), the
        # floor genuinely IS reached by a realistic extreme stack
        # (independently re-verified: composing gestational_diabetes +
        # insulin_resistant + hypothyroidism + chronic_kidney_disease +
        # PCOS + sleep_hours=5 + desaturation_events_per_hour=24.3 +
        # menstrual_cycle_phase=LUTEAL + indinavir/hydrochlorothiazide/
        # growth_hormone/atorvastatin drives the raw composed p3 to
        # ~6.3e-9, below the 1e-8 floor). This is NOT a silent-clamp
        # incident, precisely because the generic warning below was
        # already built: it fires correctly and the calibrated-effect
        # warning reaches the real REST response. Generalized here
        # originally rather than waiting for an actual silent-clamp
        # incident, mirroring row 1.23's own Gb/G_threshold
        # convention: warn whenever the GENERIC per-parameter bound
        # actually moves a value from what the composed modifiers alone
        # would have produced. Gb/G_threshold are deliberately excluded
        # from this generic warning -- their own joint clamp below is more
        # specific (accounts for the pair's relationship, not just an
        # independent bound) and already produces its own warning; adding
        # a second, independent warning for the same field would be
        # confusing, not more informative.
        clamp_warnings: list[str] = []
        # Default (neutral) when no modifier touches Gb at all -- see the
        # "if param == 'Gb'" branch inside the loop below, which overwrites
        # this with the true raw composed value when one does.
        raw_gb_before_any_clamp = base_dict["Gb"]

        for param, factor in net_factors.items():
            if param not in base_dict:
                continue
            base_val = base_dict[param]
            if not isinstance(base_val, (int, float)):
                continue
            new_val = base_val * factor

            # Clamp to physiological bounds
            if param in _PARAM_BOUNDS:
                lo, hi = _PARAM_BOUNDS[param]
                clamped_val = max(lo, min(hi, new_val))
                # Gb/G_threshold's LOWER bound is deliberately excluded from
                # this generic warning -- the joint clamp below is more
                # specific for that direction (accounts for the pair's own
                # relationship via MIN_GB_THRESHOLD_MARGIN_MM, not just an
                # independent bound) and already produces its own warning.
                # Their UPPER bound has NO other disclosure mechanism --
                # security-reviewer finding on row 1.37 (ACROMEGALY, the
                # first mechanism to measurably narrow the margin toward
                # this ceiling, ~14.4 mM vs. the 15.0 mM bound in the
                # widest realistic stack found, not yet a breach): the
                # previous blanket `param not in ("Gb", "G_threshold")`
                # exclusion silently covered this direction too, with no
                # warning at all, on a rationale ("the more specific joint
                # clamp already warns") that was only true for the lower
                # direction. Only skip the generic warning here when the
                # LOWER bound is what engaged (new_val below lo).
                is_gb_or_threshold = param in ("Gb", "G_threshold")
                lower_bound_engaged = new_val < lo
                skip_generic_warning = is_gb_or_threshold and lower_bound_engaged
                if not skip_generic_warning and abs(clamped_val - new_val) > 1e-12:
                    clamp_warnings.append(
                        f"A drug/disease modifier's own composed effect would "
                        f"have set {param} to {new_val:.6g}, but this "
                        f"platform's own physiological bound "
                        f"({lo:.6g}-{hi:.6g}) clamped it to {clamped_val:.6g} "
                        "-- part or all of the calibrated effect on this "
                        "parameter may not be visible in this simulation's "
                        "output." + _ode_coupling_caveat(param)
                    )
                new_val = clamped_val

            if param == "Gb":
                # Captured BEFORE any clamp for the joint clamp's own
                # warning message below -- fixes a security-reviewer
                # MEDIUM finding (row 1.37): the joint clamp used to read
                # `overrides.get("Gb", ...)`, which by this point may
                # already be the GENERIC-bound-clamped value (e.g. floored
                # at 2.0), not the true raw composed value -- so an
                # extreme-enough scenario could report an understated "would
                # have set Gb to X" figure. Not reachable by any currently
                # shippable combination (requires an already-abnormally-low
                # starting Gb no shipped mechanism alone produces), but a
                # real, demonstrated flaw in the two-stage clamp design.
                raw_gb_before_any_clamp = base_val * factor

            if abs(new_val - base_val) > 1e-12:
                overrides[param] = new_val

        # Joint Gb/G_threshold safety clamp -- each parameter above is
        # bounded independently, but PathwayParameters.F_secretion_basal
        # (pathway.py) is inversely proportional to (Gb - G_threshold),
        # by design to hit an exact analytic fixed point at basal glucose.
        # No modifier here currently touches G_threshold, but several
        # ordinary, realistic ones reduce Gb (metformin, ...) -- stacking
        # a small handful of them can push Gb to within a few hundredths
        # of a mM of the default G_threshold=4.5, blowing
        # F_secretion_basal up by an order of magnitude or more (verified:
        # metformin alone at 2000mg chronic dose already shrinks this
        # margin substantially, per the D2 population x mechanism
        # coverage matrix in MASTER_PLAN.md Section 11 -- physiologically
        # backwards for an insulin SENSITISER). Alcohol is NO LONGER one
        # of the Gb-reducing modifiers this comment used to name in its
        # own worked example (a fresh Track 1 scoping fan-out, 2026-08-18,
        # found the old Gb-reduction mechanism was contradicted by both
        # available primary sources -- alcohol's real effect is now
        # modelled as a transient, ODE-coupled hepatic_egp_delta
        # contribution instead, see biosim_ode.alcohol_hepatic_pathway,
        # which does not touch this static personalization pipeline or
        # Gb at all). This clamp is a numerical artifact of the
        # fixed-point formula's own sensitivity, not a modeled biological
        # mechanism, so it's fixed here rather than by inventing a
        # clinical justification for the spike. The minimum margin below
        # is an engineering safety bound protecting this formula's
        # validity, not a literature-derived clinical value -- same
        # "honest, isolated assumption" category as HPA_FEEDBACK_GAIN and
        # the confidence-tier CV defaults.
        MIN_GB_THRESHOLD_MARGIN_MM = 0.3
        requested_gb = overrides.get("Gb", base_dict["Gb"])
        effective_g_threshold = overrides.get("G_threshold", base_dict["G_threshold"])
        if requested_gb - effective_g_threshold < MIN_GB_THRESHOLD_MARGIN_MM:
            clamped_gb = effective_g_threshold + MIN_GB_THRESHOLD_MARGIN_MM
            # Only warn if the clamp actually MOVED the value -- a request
            # landing exactly on the floor (or within float noise of it)
            # wasn't overridden, just coincidentally already there.
            if abs(clamped_gb - requested_gb) > 1e-9:
                clamp_warnings.append(
                    "A drug/disease modifier's own composed effect would have "
                    # raw_gb_before_any_clamp (the true pre-any-clamp
                    # composed value), NOT requested_gb (which may already
                    # be the generic-bound-floored intermediate) -- see
                    # this variable's own comment above the main loop.
                    f"set fasting glucose (Gb) to {raw_gb_before_any_clamp:.3f} mM, but the "
                    "MIN_GB_THRESHOLD_MARGIN_MM safety clamp raised it to "
                    f"{clamped_gb:.3f} mM (G_threshold={effective_g_threshold:.3f} mM + "
                    f"{MIN_GB_THRESHOLD_MARGIN_MM} mM margin, protecting "
                    "F_secretion_basal's fixed-point validity) -- part or all of "
                    "the calibrated effect on fasting glucose may not be visible "
                    "in this simulation's output." + _ode_coupling_caveat("Gb")
                )
            overrides["Gb"] = clamped_gb

        # First-phase ("acute") insulin response is present below a
        # fasting-glucose threshold and absent above it -- Brunzell et
        # al. 1976 (JCEM, PMID:1262429, n=66, IVGTT, verified live):
        # "The acute insulin response was present in subjects with
        # fasting glucose levels below 115 mg/dl but was absent above
        # this level." first_phase_gain (pathway.py) is this model's
        # only representation of that acute response (Mari et al. 2002,
        # PMID:12388151) -- previously a bare constant with no way for
        # any disease state to touch it, so T2D's preset scaled it
        # identically to second-phase secretion via secretion_scale
        # despite its own Gb (6.5 mM) sitting above this exact threshold.
        # Gated on the FINAL effective Gb (after the joint safety clamp
        # above), not on "is T2D present" -- this directly implements
        # Brunzell's own glucose-LEVEL finding, so it correctly leaves
        # INSULIN_RESISTANT (Gb ~5.15 mM, below threshold) untouched and
        # generalises to any modifier combination that happens to cross
        # it. Set to exactly 0.0, not an invented intermediate fraction --
        # the primary source's own finding is binary ("absent"), and
        # picking a percentage with no anchor would be the same
        # unsupported-placeholder pattern this project rejects elsewhere.
        # Safe with respect to every existing fixed-point invariant:
        # first_phase_secretion = first_phase_gain * max(dG, 0) *
        # secretion_scale is exactly zero at any basal steady state
        # (dG=0 there) for ANY gain value, so this cannot perturb
        # validate_y0() regardless of which profile triggers it.
        BRUNZELL_AIR_ABSENT_THRESHOLD_MM = 115.0 / 18.0  # ~6.39 mM
        effective_gb_final = overrides.get("Gb", base_dict["Gb"])
        if effective_gb_final >= BRUNZELL_AIR_ABSENT_THRESHOLD_MM:
            overrides["first_phase_gain"] = 0.0

        # Reject drug combinations this model can't represent faithfully
        # (two distinct sulfonylureas; duplicate flat-factor events such as
        # HCQ) -- see validate_drug_combination()'s docstring. Also enforced
        # earlier, at ScenarioDefinition construction, so REST callers get a
        # 400 rather than the runner's catch-all swallowing this.
        validate_drug_combination(drug_doses)

        # Sulfonylureas (glibenclamide) trigger insulin secretion largely
        # INDEPENDENT of ambient glucose (KATP-channel closure via SUR1
        # binding) -- a genuinely additive effect on PathwayParameters
        # .glucose_independent_secretion (base value 0.0), which the
        # multiplicative net_factors/compose_all_modifiers() pipeline above
        # cannot express (any factor times a 0.0 base stays 0.0). Applied
        # directly here, matching the same "special-cased override, bypasses
        # the generic factor loop" pattern already used for the Gb/
        # G_threshold margin clamp and Brunzell's first_phase_gain zeroing
        # above. See drug_interactions.glucose_independent_secretion_
        # pM_per_min()'s own docstring for the citation and calibration.
        #
        # Gated on dose_mg > 0 (found missing by the mandatory security-
        # reviewer pass: a dose_mg=0.0 event previously still triggered the
        # full effect, since this block only checked drug identity, not
        # dose -- every other dose-scaled modifier function in this file
        # (metformin_modifiers(), caffeine_modifiers(), etc.) already
        # guards `dose.dose_mg <= 0: return {}`; this special-cased block
        # bypasses those functions entirely so needed its own copy of the
        # same guard). Deliberately NOT gated on hours_since_ingestion/
        # hours_before_simulation -- Riefflin et al. 2015 (PMID:25297572)
        # measured a CHRONIC (14-day), steady-state effect, so (matching
        # metformin_modifiers()'s own established "always chronic, timing
        # irrelevant" precedent, see that function's docstring) this is a
        # deliberate design choice, not an unexamined gap -- tested
        # explicitly in test_physiology.py.
        #
        # Skipped entirely for DiseaseState.TYPE_1_DIABETES (found by the
        # mandatory security-reviewer pass on the T1D milestone): this
        # formula computes the effect from ke/Ib alone, with no
        # dependence on secretion_scale (deliberately, since the KATP-
        # channel mechanism is glucose-INDEPENDENT, not beta-cell-mass-
        # independent) -- so T1D's near-zero secretion_scale (beta cells
        # floored at the model's 5% engineering minimum) had no bearing
        # on whether this override applied, producing a real, nonzero,
        # ~50%-of-floor insulin rise from a drug that clinically requires
        # substantial functional beta-cell mass to work at all and is not
        # prescribed in confirmed T1D. Rather than retune the formula to
        # scale by secretion_scale (which would also change T2D's own,
        # already-calibrated-against-real-data effect size -- T2D is the
        # actual population Riefflin et al. 2015 studied), skip the
        # override for this one disease state specifically -- the
        # narrower, lower-risk fix the review itself offered as
        # sufficient ("either, not both required").
        for dose in drug_doses:
            if (
                dose.drug in SULFONYLUREA_DRUGS
                and dose.dose_mg > 0
                and DiseaseState.TYPE_1_DIABETES not in profile.disease_states
            ):
                effective_ke = overrides.get("ke", base_dict["ke"])
                effective_ib = overrides.get("Ib", base_dict["Ib"])
                overrides["glucose_independent_secretion"] = glucose_independent_secretion_pM_per_min(
                    effective_ke, effective_ib, dose.drug
                )
                break

        personalized = dataclasses.replace(base_params, **overrides)

        # Build human-readable applied modifier list
        applied = self._describe_modifiers(
            profile, drug_doses, net_factors,
            glucose_independent_secretion=overrides.get("glucose_independent_secretion", 0.0),
        )

        return PhysiologicalPersonalizationResult(
            original_params=base_params,
            personalized_params=personalized,
            profile=profile,
            applied_modifiers=applied,
            drug_doses=list(drug_doses),
            combined_factors=net_factors,
            clamp_warnings=clamp_warnings,
        )

    def _describe_modifiers(
        self,
        profile: PhysiologicalProfile,
        drug_doses: list[DrugDose],
        net_factors: dict[str, float],
        glucose_independent_secretion: float = 0.0,
    ) -> list[AppliedModifier]:
        # Each entry below calls that category's OWN modifier function
        # again to get its individual (not net-combined) parameter dict --
        # net_factors is the product across every active source, so
        # reading from it here would misattribute other sources'
        # contributions to this one. `net_factors` is still accepted (and
        # still used elsewhere) but deliberately unused for per-source
        # attribution now.
        applied = []
        from .profile import DiseaseState, BiologicalSex

        if profile.age_years != 30:
            mods = age_modifiers(profile.age_years)
            if mods:
                applied.append(AppliedModifier(
                    source=f"age_{int(profile.age_years)}",
                    parameters=mods,
                    description=f"Age {profile.age_years:.0f}: insulin sensitivity and beta-cell function adjusted",
                ))

        if profile.biological_sex != BiologicalSex.NOT_SPECIFIED:
            mods = sex_modifiers(profile.biological_sex)
            if mods:
                applied.append(AppliedModifier(
                    source=f"sex_{profile.biological_sex.value}",
                    parameters=mods,
                    description=f"Biological sex {profile.biological_sex.value}: sensitivity adjusted",
                ))

        if profile.bmi_category not in ("normal",):
            mods = bmi_modifiers(profile.bmi)
            if mods:
                applied.append(AppliedModifier(
                    source=f"bmi_{profile.bmi:.1f}_{profile.bmi_category}",
                    parameters=mods,
                    description=f"BMI {profile.bmi:.1f} ({profile.bmi_category}): insulin resistance applied",
                ))

        if DiseaseState.HEALTHY not in profile.disease_states or len(profile.disease_states) > 1:
            for ds in profile.disease_states:
                if ds != DiseaseState.HEALTHY:
                    mods = disease_modifiers([ds])
                    if mods:
                        applied.append(AppliedModifier(
                            source=f"disease_{ds.value}",
                            parameters=mods,
                            description=f"Disease state: {ds.value}",
                        ))

        if profile.time_of_day_hours is not None:
            from .chronobiology import time_description
            mods = circadian_modifiers(profile.time_of_day_hours)
            if mods:
                applied.append(AppliedModifier(
                    source=f"circadian_{profile.time_of_day_hours:.0f}h",
                    parameters=mods,
                    description=f"Time of day {profile.time_of_day_hours:.0f}:00 — {time_description(profile.time_of_day_hours)}",
                ))

        if profile.sleep_hours is not None:
            mods = sleep_modifiers(profile.sleep_hours)
            if mods:
                applied.append(AppliedModifier(
                    source=f"sleep_{profile.sleep_hours:.1f}h",
                    parameters=mods,
                    description=f"Sleep duration {profile.sleep_hours:.1f}h: insulin sensitivity adjusted",
                ))

        if profile.desaturation_events_per_hour is not None:
            mods = intermittent_hypoxia_modifiers(profile.desaturation_events_per_hour)
            if mods:
                applied.append(AppliedModifier(
                    source=f"hypoxia_{profile.desaturation_events_per_hour:.1f}events_h",
                    parameters=mods,
                    description=(
                        f"Intermittent hypoxia {profile.desaturation_events_per_hour:.1f} "
                        "desaturation events/h: insulin sensitivity and glucose "
                        "effectiveness adjusted"
                    ),
                ))

        if profile.menstrual_cycle_phase is not None:
            mods = menstrual_cycle_modifiers(profile.menstrual_cycle_phase)
            if mods:
                applied.append(AppliedModifier(
                    source=f"menstrual_cycle_{profile.menstrual_cycle_phase.value}",
                    parameters=mods,
                    description=(
                        f"Menstrual cycle phase ({profile.menstrual_cycle_phase.value}): "
                        "insulin sensitivity adjusted"
                    ),
                ))

        if profile.cold_acclimatized:
            mods = cold_acclimation_modifiers(profile.cold_acclimatized)
            if mods:
                applied.append(AppliedModifier(
                    source="cold_acclimatized",
                    parameters=mods,
                    description=(
                        "Chronic (~10-day) cold acclimation: insulin sensitivity "
                        "increased"
                    ),
                ))

        for dose in drug_doses:
            if (
                dose.drug in SULFONYLUREA_DRUGS
                and dose.dose_mg > 0
                and DiseaseState.TYPE_1_DIABETES not in profile.disease_states
            ):
                # Additive effect, not expressible via the multiplicative
                # `mods` dict at all (see _sulfonylurea_modifiers()'s own
                # docstring) -- recorded separately here so the audit trail
                # doesn't silently miss this drug's real effect the way
                # BMI/T2D/metformin's own net_factors-only reporting once
                # did (see the regression tests/entries above this class).
                # Gated on dose_mg > 0 to match the actual override block
                # above -- a zero-dose event applies no effect and must not
                # be recorded as if it did. Also gated on the same
                # TYPE_1_DIABETES exclusion as the override block itself
                # (personalizer.py above) -- when that block skips applying
                # the effect, this one falls through to modifiers_for_dose()
                # below, which returns {} for sulfonylureas unconditionally
                # (see drug_interactions.py), so no misleading "applied a
                # zero-magnitude effect" entry is recorded either.
                applied.append(AppliedModifier(
                    source=f"drug_{dose.drug.value}_{dose.dose_mg:.0f}mg",
                    parameters={"glucose_independent_secretion": glucose_independent_secretion},
                    description=(
                        f"{dose.drug.value} {dose.dose_mg:.0f} mg: +"
                        f"{glucose_independent_secretion:.2f} pM/min glucose-independent "
                        "insulin secretion"
                    ),
                ))
                continue
            mods = modifiers_for_dose(dose)
            if mods:
                applied.append(AppliedModifier(
                    source=f"drug_{dose.drug.value}_{dose.dose_mg:.0f}mg",
                    parameters=mods,
                    description=f"{dose.drug.value} {dose.dose_mg:.0f} mg, {dose.hours_since_ingestion:.1f}h ago",
                ))

        return applied
