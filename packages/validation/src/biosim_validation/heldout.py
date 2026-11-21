# SAVE TO: packages\validation\src\biosim_validation\heldout.py
"""
Held-out validation datasets (finish-plan item B2, 2026-11-16).

Published human data the model was NEVER fitted or calibrated to, run against
the model exactly as it is (default parameters for the population, nothing
adjusted). Headline band: STRICT, +/-20% of the published mean. See
docs/VALIDATION.md for the method, the calibration inventory and the results.

Every reference value below was read from full text fetched from Europe PMC
on 2026-11-16; each point's `notes` quote the table row verbatim and say where
it is. Fasting values are not scored: they are the model's starting state.
"""

from __future__ import annotations

from typing import Callable

import numpy as np

from .benchmark import BenchmarkDataset, BenchmarkPoint, Tolerance

HELDOUT_STRICT_FRACTION = 0.20
MG_DL_PER_MMOL_L = 18.016  # glucose, 180.16 g/mol

# Every source used to fit, calibrate or benchmark the model before B2
# (docs/VALIDATION.md, "Calibration inventory"). A held-out source must not
# be in here; test_heldout.py also checks against every seeded parameter's
# PMID and every existing benchmark's PMID.
CALIBRATION_SOURCE_PMIDS: frozenset[str] = frozenset({
    "1566835",   # Van Cauter 1992: OGTT benchmark (absorption + first-phase gain fitted to it)
    "382871",    # DeFronzo 1979: was cited (wrongly) for the OGTT 2-h point; unscored since B4
    "7033284",   # Bergman 1981: p3, IVGTT benchmark
    "3640682",   # MINMOD: p2
    "4856884",   # Sherwin 1974: ke
    "30166558",  # Rebelato 2018: KG_secretion
    "12388151",  # Mari 2002: first-phase functional form
    "9575818",   # Pacini 1998: p1 consistency check
    "38389141",  # Mora-Gonzalez 2024: glucose distribution volume
    "651653",    # Radziuk 1978: hepatic first-pass extraction
    "19196896",  # Basu 2009: T2D p1 factor AND T2D meal-insulin benchmark
    "12643178",  # Beck-Nielsen 2002: T2D hepatic EGP preset
    "8138058",   # Orskov 1994: GLP-1 meal response
    "27686734",  # Chen 2016: acarbose kabs scale
    "3042310",   # Samad 1988: miglitol kabs scale
    "3538843",   # Krezowski 1986: protein-driven secretion
    "6389060",   # Nuttall 1984: protein-driven secretion
    "28941314", "32369480", "39049802", "23412078",  # drug benchmarks
    "33515092", "11478588", "15465742",              # protein / caffeine benchmarks
})

_NAVARRO_TABLE2 = (
    'Navarro-Masip et al. 2025, Int J Mol Sci, PMID:40564902 (PMC12193239), '
    'Table 2: "ZAG, insulin, C-peptide, and glucose values during the oral glucose '
    'tolerance test (OGTT) with 75 g of glucose, the intake of a standardized liquid '
    'nutritional supplement, and the intravenous administration of 1 mg of glucagon." '
    'Footnote: "Data are mean ± SD." n=24 healthy non-obese adults (Table 1: age 53.0 ± 16.9, '
    'BMI 24.8 ± 4.9, 12 M / 12 F). UNIT NOTE: the glucose row is labelled "(mmol/L)" but '
    'its values (94-188) are only possible in mg/dL, and Table 1 gives the same fasting '
    'value as "Fasting glucose (mg/dL) | 100.4 ± 10.8" (= Table 2 glucagon 0 min), so '
    'they are converted with /18.016. '
    'Insulin is NOT used: labelled pmol/L, but a fasting mean of 13.3 matches muU/mL, and '
    'the text does not settle which. '
)
_NAVARRO_OGTT_ROW = (
    'OGTT row (0/30/60/90/120 min): "Glucose (mmol/L) | 94.1 ± 7.6 | 167.5 ± 31.0 | '
    '187.5 ± 44.9 | 176.0 ± 53.2 | 150.4 ± 46.0".'
)
_NAVARRO_MEAL_ROW = (
    'Liquid-supplement row (0/60/90/120/180 min): "Glucose (mmol/L) | 102.4 ± 10.2 | '
    '129.4 ± 27.4 | 115.2 ± 28.5 | 104.3 ± 23.8 | 90.0 ± 10.8". Meal (Methods 4.2): '
    '"200 mL, 150 kcal; containing 6 g of protein, 5.8 g of fat, and 18.4 g of '
    'carbohydrates; Fortisip". Fat is not modelled.'
)
_ESHGHI_TABLE5 = (
    'Eshghi et al. 2017, Front Endocrinol, PMID:28744255 (PMC5504214), Table 5: '
    '"Concentrations of energy substrate and hormones during fasting and following an oral '
    'glucose tolerance test on day 2." Footnote: "Data presented as mean ± SEM." Rest '
    '(no-exercise) arm used. n=6 healthy men (age 42.5 ± 10.5) and n=6 men with T2D '
    '(age 60.5 ± 8.5, A1c 6.4%, 3.9 years since diagnosis; mild, early T2D) (Table 1). '
    'The T2D men took metformin, withheld before the study days. Protocol: 75 g glucose, 4-h OGTT, samples at '
    '15, 30, 45, 60, 90, 120, 150, 180, 210, 240 min; "The iAUC was calculated by '
    'subtracting the average of the two fasting values from the AUC" (trapezoid). '
)


def _point(name: str, reference: float, unit: str, pmid: str, notes: str,
           population: str = "healthy_adults", **spread) -> BenchmarkPoint:
    # Every held-out point is in the LOCKED test set (S0, 2026-11-20): it may
    # never be used to choose or adjust anything (docs/GOLD_STANDARD.md).
    return BenchmarkPoint(
        name=name, description=f"Held-out: {name}", reference=reference, unit=unit,
        tolerance=Tolerance.STRICT, pubmed_id=pmid, population=population, notes=notes,
        benchmark_set="locked", **spread,
    )


def navarro_ogtt() -> BenchmarkDataset:
    ds = BenchmarkDataset(
        name="HeldOut_Navarro2025_OGTT_75g_healthy",
        description="75 g OGTT, 24 healthy non-obese adults (glucose only)",
    )
    note = _NAVARRO_TABLE2 + _NAVARRO_OGTT_ROW
    for t, mgdl, sd in ((30, 167.5, 31.0), (60, 187.5, 44.9), (90, 176.0, 53.2), (120, 150.4, 46.0)):
        ds.add(_point(f"glucose_{t}min", mgdl / MG_DL_PER_MMOL_L, "mM", "40564902", note,
                      spread_kind="sd", sd=sd / MG_DL_PER_MMOL_L, n=24,
                      quote=f"{mgdl} ± {sd}", source_location="Table 2, OGTT glucose row"))
    # No separate peak / time-to-peak points: the sampled peak IS the 60-min
    # value, and with 30-min sampling a +/-20% time band only admits 60 min,
    # so both would double-count or pass trivially (B2 code review, M1).
    return ds


def navarro_liquid_meal() -> BenchmarkDataset:
    ds = BenchmarkDataset(
        name="HeldOut_Navarro2025_liquid_meal_healthy",
        description="Liquid mixed meal (18.4 g CHO, 6 g protein), same 24 adults (glucose only)",
    )
    note = _NAVARRO_TABLE2 + _NAVARRO_MEAL_ROW
    for t, mgdl, sd in ((60, 129.4, 27.4), (90, 115.2, 28.5), (120, 104.3, 23.8), (180, 90.0, 10.8)):
        ds.add(_point(f"glucose_{t}min", mgdl / MG_DL_PER_MMOL_L, "mM", "40564902", note,
                      spread_kind="sd", sd=sd / MG_DL_PER_MMOL_L, n=24,
                      quote=f"{mgdl} ± {sd}", source_location="Table 2, liquid-supplement glucose row"))
    return ds


def eshghi_ogtt() -> BenchmarkDataset:
    ds = BenchmarkDataset(
        name="HeldOut_Eshghi2017_OGTT_4h_healthy_vs_T2D",
        description="4-h 75 g OGTT, 6 healthy men vs 6 men with T2D (rest arm), iAUC",
    )
    glucose_row = (_ESHGHI_TABLE5 + 'Glucose, Rest row: "4.7 ± 0.2 | 6.7 ± 0.6 | ... | '
                   '269 ± 138 | 902 ± 193 | ..." (fasting healthy | T2D | iAUC healthy | T2D), '
                   'mmol/L x min.')
    insulin_row = (_ESHGHI_TABLE5 + 'Insulin (pmol/L), Rest row: "40.7 ± 9.2 | 32.6 ± 8.7 | '
                   '47,004 ± 14,486 | 20,144 ± 5,078 | ..." (fasting healthy | T2D | iAUC '
                   'healthy | T2D), pmol/L x min.')
    loc = "Table 5, Rest row (mean ± SEM)"
    ds.add(_point("glucose_iauc_healthy", 269.0, "mM*min", "28744255", glucose_row,
                  spread_kind="sem", sem=138.0, n=6, quote="269 ± 138", source_location=loc))
    ds.add(_point("glucose_iauc_t2d", 902.0, "mM*min", "28744255", glucose_row,
                  population="type_2_diabetes",
                  spread_kind="sem", sem=193.0, n=6, quote="902 ± 193", source_location=loc))
    ds.add(_point("insulin_iauc_healthy", 47004.0, "pM*min", "28744255", insulin_row,
                  spread_kind="sem", sem=14486.0, n=6, quote="47,004 ± 14,486", source_location=loc))
    ds.add(_point("insulin_iauc_t2d", 20144.0, "pM*min", "28744255", insulin_row,
                  population="type_2_diabetes",
                  spread_kind="sem", sem=5078.0, n=6, quote="20,144 ± 5,078", source_location=loc))
    return ds


def heldout_datasets() -> list[BenchmarkDataset]:
    return [navarro_ogtt(), navarro_liquid_meal(), eshghi_ogtt()]


# ---------------------------------------------------------------------------
# Model runs (default parameters, nothing adjusted)
# ---------------------------------------------------------------------------

def _sample_mM(output, species: str, times_min) -> np.ndarray:
    sim = output.simulation
    values = np.interp(np.asarray(times_min, float) * 60.0, sim.times_s, sim.species[species])
    return values / 1000.0 if species == "G" else values * 1e6  # G uM->mM; I uM->pM


def _iauc(output, species: str, times_min) -> float:
    """Trapezoid AUC at the study's sampling times minus fasting x duration."""
    y = _sample_mM(output, species, times_min)
    return float(np.trapezoid(y, times_min) - y[0] * (times_min[-1] - times_min[0]))


def _navarro_ogtt_values(runner) -> dict[str, float]:
    from biosim_api import ogtt_75g
    out = runner.run(ogtt_75g())
    times = [30, 60, 90, 120]
    return {f"glucose_{t}min": float(v) for t, v in zip(times, _sample_mM(out, "G", times))}


def _navarro_meal_values(runner) -> dict[str, float]:
    from biosim_api import MealEvent, ScenarioDefinition
    out = runner.run(ScenarioDefinition(
        name="heldout_fortisip", duration_s=14400,
        meal_events=[MealEvent(time_s=0.0, carbs_g=18.4, protein_g=6.0, label="Fortisip")],
    ))
    times = [60, 90, 120, 180]
    return {f"glucose_{t}min": float(v) for t, v in zip(times, _sample_mM(out, "G", times))}


def _eshghi_values(runner) -> dict[str, float]:
    from biosim_api import PersonalizationSpec, ogtt_75g
    from biosim_physiology import DiseaseState, PhysiologicalProfile
    times = [0, 15, 30, 45, 60, 90, 120, 150, 180, 210, 240]
    healthy = runner.run(ogtt_75g())
    t2d = runner.run(ogtt_75g(personalization=PersonalizationSpec(
        physiological_profile=PhysiologicalProfile(disease_states=[DiseaseState.TYPE_2_DIABETES]),
    )))
    return {
        "glucose_iauc_healthy": _iauc(healthy, "G", times),
        "glucose_iauc_t2d": _iauc(t2d, "G", times),
        "insulin_iauc_healthy": _iauc(healthy, "I", times),
        "insulin_iauc_t2d": _iauc(t2d, "I", times),
    }


HELDOUT_RUNS: dict[str, Callable] = {
    "HeldOut_Navarro2025_OGTT_75g_healthy": _navarro_ogtt_values,
    "HeldOut_Navarro2025_liquid_meal_healthy": _navarro_meal_values,
    "HeldOut_Eshghi2017_OGTT_4h_healthy_vs_T2D": _eshghi_values,
}
