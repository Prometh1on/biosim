# SAVE TO: packages\physiology\src\biosim_physiology\drug_interactions\classic.py
"""The original (v1) drugs: caffeine, metformin, glucocorticoid, alcohol."""

from __future__ import annotations

import math

from ._core import DrugDose

# ---------------------------------------------------------------------------
# Per-drug modifier functions
# ---------------------------------------------------------------------------

def caffeine_modifiers(dose: DrugDose) -> dict[str, float]:
    """
    Caffeine acutely reduces insulin sensitivity via adenosine antagonism
    and sympathomimetic catecholamine release.

    Peak effect at ~45 min; half-life ~5 hours.
    Effect magnitude: ~15% reduction at 200 mg; dose-dependent.

    Reference: Keijzers et al. 2002, Diabetes Care. PMID:11815511
    """
    # Pharmacokinetic: fraction of peak effect remaining
    t = dose.hours_since_ingestion
    absorption_t_half = 0.75   # hours
    elimination_t_half = 5.0   # hours

    concentration_fraction = (
        0.0 if t <= 0
        else (1.0 - math.exp(-math.log(2) / absorption_t_half * t))
             * math.exp(-math.log(2) / elimination_t_half * t)
    )

    # Max effect: 200 mg reduces p3 by ~15%, with log-dose response
    if dose.dose_mg <= 0:
        return {}
    dose_factor = math.log(1.0 + dose.dose_mg / 100.0) / math.log(3.0)
    max_reduction = 0.15 * dose_factor
    p3_factor = 1.0 - max_reduction * concentration_fraction

    # A +2% fasting-glucose (Gb) term was removed on 2026-11-20 (gold-standard
    # drug stream): it had no source (Keijzers 2002 reports only the insulin-
    # sensitivity change), and in this model a Gb change also moves the
    # beta-cell reference point (a chronic set point), so for an acute drug it
    # LOWERED insulin secretion -- the opposite of caffeine's real effect.
    return {"p3": p3_factor}


def metformin_modifiers(dose: DrugDose) -> dict[str, float]:
    """
    Metformin (chronic dosing) improves insulin sensitivity and reduces
    hepatic glucose output.

    Standard dose 500-2000 mg/day.
    Effects (see CLAUDE.md's metformin mechanism-weighting investigation
    for the full literature-verification trail):
        p3: no effect (peripheral insulin sensitivity) -- REMOVED.
            Previously +20-30%. This was flagged as a genuine, unresolved
            conflict between two placebo-controlled studies: a
            gold-standard tracer/clamp study in CHRONIC (15wk) T2D
            (PMID:8923861) found NO improvement in insulin-mediated
            glucose disposal (P=NS), while a different crossover study in
            a different population (acute, 2-day, obese IGT --
            pre-diabetic, not T2D -- PMID:8582131) found a real SI
            increase. Resolved (not by picking a side in general, but by
            checking which regime THIS function actually represents):
            `metformin_modifiers()` has no `DrugDose.is_chronic` gating
            at all -- every dose is treated as the chronic, steady-state
            effect (this function's own docstring title), and this
            platform's own metformin validation benchmark
            (`validate_metformin()`, Hansen et al. 2024, PMID:39049802)
            is a CHRONIC (14-day), T2D-population trial -- the same
            regime as PMID:8923861 (chronic, T2D), not PMID:8582131
            (acute, pre-diabetic obese IGT). Since this model only ever
            represents the chronic/T2D regime for metformin, the acute/
            pre-diabetic study is evidence about a condition this
            function doesn't claim to model -- the applicable evidence
            (PMID:8923861) says no peripheral effect. Removed rather
            than resolved by fiat in either direction.
        Gb: -14.7% (reduces fasting glucose via suppressed hepatic
            glucose production). Was an uncited -4%; corrected to match
            PMID:8923861 directly -- placebo-controlled, chronic (15wk),
            T2D, tracer-measured: basal hepatic glucose production fell
            12.9->11.0 umol/kg.min (-14.7%, P<0.01). Uses a first-order
            approximation, stated explicitly: this model's Gb is a
            steady-state glucose CONCENTRATION (mM), not the FLUX
            (umol/kg.min) the tracer study measured -- absent a
            more principled translation, a given %-reduction in hepatic
            production is assumed to produce approximately the same
            %-reduction in the resulting basal concentration (holding
            clearance kinetics fixed). Flagged as an approximation to
            revisit, not treated as an exact derivation.
        p1: no effect (glucose effectiveness) -- REMOVED. Previously +10%,
            cited to Bailey & Turner 1996 (PMID:8569826, a narrative review
            with no indexed data -- confirmed via PubMedFetcher this
            couldn't be the actual source of that number). Two independent,
            placebo-controlled studies specifically measuring Bergman
            minimal-model glucose effectiveness (S_G) found NO significant
            change with metformin: PMID:8582131 ("AUC glucose, KG and SG
            were similar in both tests", p=NS) and PMID:16968813 ("In
            contrast to metformin, pioglitazone improved... insulin-
            independent glucose disposal [glucose effectiveness (SG)]" --
            i.e. metformin did not). Both verified via PubMedFetcher.

    Reference: Bailey & Turner 1996, NEJM. PMID:8569826 (narrative review,
    does not itself support the specific p3/Gb percentages above -- see
    CLAUDE.md's mechanism-weighting investigation for the real quantitative
    sources being reconciled against it).
    """
    if dose.dose_mg <= 0:
        return {}
    # Log-linear dose response capped at 2000 mg
    effective_dose = min(dose.dose_mg, 2000.0)
    dose_factor = math.log(1.0 + effective_dose / 500.0) / math.log(5.0)
    return {
        # 0.147 = a flux->concentration approximation (see this
        # function's own docstring above), not a literature-exact value
        # -- don't read this as more precise than it is.
        "Gb": 1.0 - 0.147 * dose_factor,
    }


def glucocorticoid_modifiers(dose: DrugDose) -> dict[str, float]:
    """
    Glucocorticoids (e.g. prednisone) cause steroid-induced insulin resistance.

    Acute effect: reduces p3 ~30%, elevates Gb ~10%.
    Dose equivalent: 10 mg prednisone = reference.

    Reference: Rizza et al. 1982, J Clin Endocrinol Metab. PMID:7033265

    Note: `SimulationRunner` does NOT use this function's p3/Gb/secretion_
    scale values for a glucocorticoid drug event -- `CortisolPathway`
    (biosim_ode.cortisol_pathway) is the sole, real-PK/PD-timing-aware
    representation of all three effects for that path (see
    `_detect_cortisol()`/`glucocorticoid_secretion_scale_factor()` in
    biosim_api's runner.py, and CLAUDE.md's glucocorticoid PK/PD milestone
    entry). This function still returns real, dose-only (untimed) values
    for any OTHER caller composing DrugDose lists directly (e.g.
    `drug_modifiers()`/`compose_all_modifiers()`), not just a stub.
    """
    if dose.dose_mg <= 0:
        return {}
    # Dose relative to 10 mg prednisone equivalent
    relative_dose = dose.dose_mg / 10.0
    dose_factor = min(relative_dose, 3.0) / 3.0  # cap at 30 mg equivalent
    return {
        "p3": 1.0 - 0.30 * dose_factor,
        "Gb": 1.0 + 0.10 * dose_factor,
        "secretion_scale": 1.0 + 0.10 * dose_factor,  # compensatory
    }


def alcohol_modifiers(dose: DrugDose) -> dict[str, float]:
    """
    Intentionally returns {} -- NOT a missing case, and NOT the previous
    "uncited placeholder" state either. A fresh Track 1 scoping fan-out
    (2026-08-18) re-checked the previously-uncited Gb-reduction mechanism
    here and found it actively WRONG, contradicted by both real primary
    sources located: Siler et al. 1998 (PMID:9815011, 48g oral alcohol,
    stable-isotope tracers) explicitly states alcohol has "no effect on
    plasma glucose concentrations" despite a real, measured -12% hepatic
    glucose OUTPUT (flux) effect; Shelmet et al. 1988 (PMID:3280601, IV
    ethanol) independently confirms "no effect on glucose tolerance"
    despite real peripheral/hepatic flux effects. A Gb (fasting glucose
    CONCENTRATION set-point) modifier is therefore the wrong mechanism
    shape entirely -- alcohol's real effect is a TRANSIENT reduction in
    hepatic glucose PRODUCTION, exactly the same "flux, not concentration
    set-point" distinction already drawn for GLP-1/SGLT2i/sulfonylureas
    above, just discovered later for this drug.

    Alcohol's real effect is now modelled via a genuine ODE-coupled
    pathway instead: biosim_ode.alcohol_hepatic_pathway.
    AlcoholHepaticPathway, a stateless, transient (peak-then-decay)
    contribution to the hepatic_egp_delta rule-11 signal, calibrated
    directly to Siler's -12% figure. SimulationRunner detects
    Drug.ALCOHOL and composes it instead of applying a modifier here --
    see AlcoholHepaticPathway's own module docstring for the full
    citation detail and the reasoning for NOT using EGPParameters/
    HepGlucPathway's own chronic-preset machinery (alcohol is acute and
    self-resolving, unlike T2D/hyperthyroidism/metformin).
    """
    return {}
