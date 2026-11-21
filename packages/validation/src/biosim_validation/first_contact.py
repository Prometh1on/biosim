# SAVE TO: packages\validation\src\biosim_validation\first_contact.py
"""
H2: first-contact held-out test (pre-registered 2026-11-21, commit d804d64).

Three published datasets no part of this project had used before, run ONCE
against the model frozen at commit 85bb29e, with nothing adjusted. The
datasets, the model setup and the scoring were fixed in docs/VALIDATION.md
section 8 before any run; see there for the full method and the result.

Every reference value is from a table in open-access full text fetched from
Europe PMC on 2026-11-21. No value is read off a figure.

Deliberately NOT part of grading.locked_datasets() (which would change the
frozen manifest hash): this is a separate group, read on its own.
"""

from __future__ import annotations

from typing import Callable

import numpy as np

from .benchmark import BenchmarkDataset, BenchmarkPoint, Tolerance

PMOL_PER_MU_INSULIN = 6.0  # pmol/L per mU/L, the conversion fixed in the pre-registration
FIRST_CONTACT_PMIDS: frozenset[str] = frozenset({"20937162", "21118565", "36432554"})
LUND_TIMES_MIN = (0, 15, 30, 45, 60, 90, 120)
YOSHIDA_TIMES_MIN = (0, 30, 60, 120)

_WICKENBERG = (
    "Wickenberg et al. Effects of Curcuma longa (turmeric) on "
    "postprandial plasma glucose and insulin in healthy subjects. Nutr J 2010;9:43. "
    "PMID:20937162 (PMC2964546). 14 healthy adults (7 M / 7 F, age 29 ± 1 y, BMI 23.9 ± 2.7). "
    "Reference arm: 75 g OGTT with placebo capsules. Capillary plasma-equivalent glucose, "
    "serum insulin, 0/15/30/45/60/90/120 min. Incremental AUC 'calculated above baseline'. "
    "'All values are means ± SEM, n = 14.' "
)
_JOSIC = (
    "Josic et al. Does green tea affect "
    "postprandial glucose, insulin and satiety in healthy subjects: a randomized controlled "
    "trial. Nutr J 2010;9:63. PMID:21118565 (PMC3002911). 14 healthy adults (7 M / 7 F, "
    "age 27 ± 3 y, BMI 22.3 ± 3.4); insulin n = 13. Reference meal (with hot water): 100 g "
    "white bread + 25 g smoked turkey = 50 g carbohydrate, 12.5 g protein, 3.75 g fat. "
    "Capillary plasma-equivalent glucose (HemoCue), serum insulin, 0/15/30/45/60/90/120 min. "
    "Incremental AUC; 'All areas below the baseline were excluded'. Table 1: 'All values are "
    "means ± SEM, n = 14 for glucose, n = 13 for insulin.' "
)
_YOSHIDA = (
    "Yoshida et al. Age-Related Sex Differences in Glucose Tolerance by 75 g Oral Glucose "
    "Tolerance Test in Japanese. Nutrients 2022;14:4868. PMID:36432554 (PMC9698682). Young "
    "adults aged 22-29 without diagnosed diabetes (Table 2: 82% NGT in men, 94% in women). "
    "75 g OGTT, 0/30/60/120 min, AUC 'based on the trapezoid rule' (total AUC). Table 1: "
    "'Data presented as means ± SD.' "
)


def _point(name: str, reference: float, unit: str, pmid: str, notes: str, *,
           n: int, quote: str, source_location: str,
           sem: float | None = None, sd: float | None = None) -> BenchmarkPoint:
    spread = {"spread_kind": "sem", "sem": sem} if sem is not None else {"spread_kind": "sd", "sd": sd}
    return BenchmarkPoint(
        name=name, description=f"First contact (H2): {name}", reference=reference,
        unit=unit, tolerance=Tolerance.STRICT, pubmed_id=pmid, notes=notes,
        benchmark_set="locked", n=n, quote=quote, source_location=source_location, **spread,
    )


def wickenberg_ogtt() -> BenchmarkDataset:
    ds = BenchmarkDataset(
        name="FirstContact_Wickenberg2010_OGTT_75g_healthy",
        description="75 g OGTT (placebo arm), 14 healthy adults, incremental AUCs",
    )
    glucose_note = _WICKENBERG + (
        'Table 1, 0-120 min, Reference: "253.9 ± 32.7". The table labels the unit '
        '"mU min/L", impossible for glucose; read as mmol·min/L (pre-registered).'
    )
    ds.add(_point("glucose_iauc_0_120", 253.9, "mM*min", "20937162", glucose_note,
                  n=14, sem=32.7, quote="253.9 ± 32.7", source_location="Table 1, Reference, 0-120 min"))
    for window, mu, sem in (("0_30", 584.3, 61.2), ("0_120", 2749.5, 283.7)):
        ds.add(_point(
            f"insulin_iauc_{window}", mu * PMOL_PER_MU_INSULIN, "pM*min", "20937162",
            _WICKENBERG + f'Table 2, insulin AUC (mU min/L), Reference, {window.replace("_", "-")} '
            f'min: "{mu} ± {sem}"; converted at {PMOL_PER_MU_INSULIN} pmol per mU.',
            n=14, sem=sem * PMOL_PER_MU_INSULIN, quote=f"{mu} ± {sem}",
            source_location=f"Table 2, Reference, {window.replace('_', '-')} min"))
    return ds


def josic_bread_meal() -> BenchmarkDataset:
    ds = BenchmarkDataset(
        name="FirstContact_Josic2010_bread_meal_healthy",
        description="White bread + turkey meal (50 g carbohydrate), 14 healthy adults, incremental AUCs",
    )
    for window, ref, sem in (("0_30", 42.9, 5.4), ("0_120", 177.8, 24.7)):
        ds.add(_point(
            f"glucose_iauc_{window}", ref, "mM*min", "21118565",
            _JOSIC + f'Table 1, glucose AUC (mmol·min/l), Reference, {window.replace("_", "-")} '
            f'min: "{ref} ± {sem}".',
            n=14, sem=sem, quote=f"{ref} ± {sem}",
            source_location=f"Table 1, Reference glucose, {window.replace('_', '-')} min"))
    for window, mu, sem in (("0_30", 537.5, 60.5), ("0_120", 2902.2, 224.0)):
        ds.add(_point(
            f"insulin_iauc_{window}", mu * PMOL_PER_MU_INSULIN, "pM*min", "21118565",
            _JOSIC + f'Table 1, insulin AUC (mU·min/l), Reference, {window.replace("_", "-")} '
            f'min: "{mu} ± {sem}"; converted at {PMOL_PER_MU_INSULIN} pmol per mU.',
            n=13, sem=sem * PMOL_PER_MU_INSULIN, quote=f"{mu} ± {sem}",
            source_location=f"Table 1, Reference insulin, {window.replace('_', '-')} min"))
    return ds


def yoshida_young_ogtt() -> BenchmarkDataset:
    ds = BenchmarkDataset(
        name="FirstContact_Yoshida2022_OGTT_75g_young_adults",
        description="75 g OGTT, 505 young men and 280 young women (Japan), total AUCs",
    )
    rows = (
        ("male", 505, 768.0, 130.0, 35588.0, 20485.0),
        ("female", 280, 713.0, 137.0, 41019.0, 20112.0),
    )
    for sex, n, aucg, aucg_sd, auci, auci_sd in rows:
        group = f"Young {'Males' if sex == 'male' else 'Females'} (n = {n})"
        ds.add(_point(
            f"glucose_auc_young_{sex}", aucg, "mM*min", "36432554",
            _YOSHIDA + f'Table 1, {group}: "AUCg (mmol/L·min) {aucg:g} ± {aucg_sd:g}". '
            "Total AUC: dominated by the fasting level, so +/-20% is lenient here.",
            n=n, sd=aucg_sd, quote=f"{aucg:g} ± {aucg_sd:g}", source_location=f"Table 1, {group}, AUCg"))
        ds.add(_point(
            f"insulin_auc_young_{sex}", auci, "pM*min", "36432554",
            _YOSHIDA + f'Table 1, {group}: "AUCi (pmol/L·min) {auci:,.0f} ± {auci_sd:,.0f}".',
            n=n, sd=auci_sd, quote=f"{auci:,.0f} ± {auci_sd:,.0f}", source_location=f"Table 1, {group}, AUCi"))
    return ds


def first_contact_datasets() -> list[BenchmarkDataset]:
    return [wickenberg_ogtt(), josic_bread_meal(), yoshida_young_ogtt()]


# ---------------------------------------------------------------------------
# Model runs (fixed in the pre-registration; nothing adjusted)
# ---------------------------------------------------------------------------

def _sample(output, species: str, times_min) -> np.ndarray:
    sim = output.simulation
    values = np.interp(np.asarray(times_min, float) * 60.0, sim.times_s, sim.species[species])
    return values / 1000.0 if species == "G" else values * 1e6  # G uM->mM; I uM->pM


def _iauc_above_baseline(output, species: str, end_min: int) -> float:
    """Trapezoid of max(value - value at 0, 0) at the study's sampling times up to end_min."""
    times = [t for t in LUND_TIMES_MIN if t <= end_min]
    y = _sample(output, species, times)
    return float(np.trapezoid(np.maximum(y - y[0], 0.0), times))


def _total_auc(output, species: str) -> float:
    y = _sample(output, species, YOSHIDA_TIMES_MIN)
    return float(np.trapezoid(y, YOSHIDA_TIMES_MIN))


def _wickenberg_values(runner) -> dict[str, float]:
    from biosim_api import ogtt_75g
    out = runner.run(ogtt_75g())
    return {
        "glucose_iauc_0_120": _iauc_above_baseline(out, "G", 120),
        "insulin_iauc_0_30": _iauc_above_baseline(out, "I", 30),
        "insulin_iauc_0_120": _iauc_above_baseline(out, "I", 120),
    }


def _josic_values(runner) -> dict[str, float]:
    from biosim_api import MealEvent, ScenarioDefinition
    out = runner.run(ScenarioDefinition(
        name="first_contact_bread_meal", duration_s=14400,
        meal_events=[MealEvent(time_s=0.0, carbs_g=50.0, protein_g=12.5, label="bread + turkey")],
    ))
    return {
        f"{name}_iauc_{w}": _iauc_above_baseline(out, sp, end)
        for name, sp in (("glucose", "G"), ("insulin", "I"))
        for w, end in (("0_30", 30), ("0_120", 120))
    }


def _yoshida_values(runner) -> dict[str, float]:
    from biosim_api import PersonalizationSpec, ogtt_75g
    from biosim_physiology import PhysiologicalProfile
    from biosim_physiology.profile import BiologicalSex
    groups = {
        "male": (23.7, BiologicalSex.MALE, 65.2, 172.4),
        "female": (23.3, BiologicalSex.FEMALE, 50.8, 159.7),
    }
    values: dict[str, float] = {}
    for sex, (age, bio_sex, weight, height) in groups.items():
        out = runner.run(ogtt_75g(personalization=PersonalizationSpec(
            physiological_profile=PhysiologicalProfile(
                age_years=age, biological_sex=bio_sex, body_weight_kg=weight, height_cm=height,
            ),
        )))
        values[f"glucose_auc_young_{sex}"] = _total_auc(out, "G")
        values[f"insulin_auc_young_{sex}"] = _total_auc(out, "I")
    return values


FIRST_CONTACT_RUNS: dict[str, Callable] = {
    "FirstContact_Wickenberg2010_OGTT_75g_healthy": _wickenberg_values,
    "FirstContact_Josic2010_bread_meal_healthy": _josic_values,
    "FirstContact_Yoshida2022_OGTT_75g_young_adults": _yoshida_values,
}
