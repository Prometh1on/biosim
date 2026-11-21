# SAVE TO: packages\physiology\src\biosim_physiology\modifiers\demographics.py
"""Demographic and fitness modifiers: age, sex, BMI, fitness level."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..profile import (
        BiologicalSex,
        FitnessLevel,
    )


# ---------------------------------------------------------------------------
# Age modifiers
# ---------------------------------------------------------------------------
# Reference age: 30 years. Effects compound per year beyond reference.
AGE_REFERENCE = 30.0
AGE_P3_DECLINE_PER_YEAR = 0.005     # 0.5% per year reduction in insulin sensitivity
AGE_SEC_DECLINE_PER_YEAR = 0.003    # 0.3% per year reduction in beta-cell function
AGE_P3_FLOOR = 0.30                 # never below 30% of reference
AGE_SEC_FLOOR = 0.50                # never below 50% of reference


def age_modifiers(age_years: float) -> dict[str, float]:
    """
    Age-related changes to insulin sensitivity and beta-cell function.

    Insulin sensitivity (p3) declines ~0.5% per year after 30.
    Beta-cell secretion (secretion_scale) declines ~0.3% per year.
    Both are floored to prevent physiologically impossible values.
    """
    years_over_ref = max(0.0, age_years - AGE_REFERENCE)
    p3_factor  = max(AGE_P3_FLOOR,  1.0 - AGE_P3_DECLINE_PER_YEAR  * years_over_ref)
    sec_factor = max(AGE_SEC_FLOOR, 1.0 - AGE_SEC_DECLINE_PER_YEAR * years_over_ref)
    return {"p3": p3_factor, "secretion_scale": sec_factor}


# ---------------------------------------------------------------------------
# Biological sex modifiers
# ---------------------------------------------------------------------------
# Males have ~10% higher glucose effectiveness (p1) due to lean mass.
# Pre-menopausal females have ~10% higher insulin sensitivity (p3).
# Simplified: using premenopausal female as the reference for sex=FEMALE.
# Menopausal status is not captured in v1; this is a Sprint 8 extension.

def sex_modifiers(sex: "BiologicalSex") -> dict[str, float]:
    from ..profile import BiologicalSex
    if sex == BiologicalSex.MALE:
        return {"p1": 1.10, "p3": 0.92}  # higher glucose clearance, slightly lower sensitivity
    if sex == BiologicalSex.FEMALE:
        return {"p1": 0.93, "p3": 1.10}  # lower clearance, higher sensitivity (premenopausal)
    return {}  # NOT_SPECIFIED: no adjustment


# ---------------------------------------------------------------------------
# BMI / adiposity modifiers
# ---------------------------------------------------------------------------
# Excess adipose tissue — especially visceral fat — causes insulin resistance
# via inflammatory cytokines and free fatty acid elevation.
# Effects are on p3 (insulin action) and slightly on Gb (fasting glucose).

def bmi_modifiers(bmi: float) -> dict[str, float]:
    """
    BMI-related insulin resistance and fasting glucose elevation.

    Normal (18.5–25)  : reference
    Overweight (25–30): ~20% reduction in insulin sensitivity
    Obese (30–35)     : ~35% reduction
    Severely obese >35: ~50% reduction
    """
    if bmi < 18.5:
        # Underweight: slightly increased sensitivity, reduced secretion
        return {"p3": 1.05, "secretion_scale": 0.90}
    if bmi < 25.0:
        return {}  # reference
    if bmi < 30.0:
        return {"p3": 0.80, "Gb": 1.01}  # overweight
    if bmi < 35.0:
        return {"p3": 0.65, "Gb": 1.02}  # obese
    return {"p3": 0.50, "Gb": 1.04}      # severely obese


# ---------------------------------------------------------------------------
# Fitness modifiers
# ---------------------------------------------------------------------------
# Exercise training dramatically increases GLUT4 expression and insulin
# receptor number. Athletes can have twice the insulin sensitivity of
# sedentary individuals.

_FITNESS_P3_MULTIPLIERS = {
    "sedentary":         0.72,
    "lightly_active":    0.87,
    "moderately_active": 1.00,  # reference
    "very_active":       1.28,
    "athlete":           1.55,
}

_FITNESS_SEC_MULTIPLIERS = {
    "sedentary":         0.95,
    "lightly_active":    0.98,
    "moderately_active": 1.00,
    "very_active":       1.05,
    "athlete":           1.10,
}


def fitness_modifiers(fitness: "FitnessLevel") -> dict[str, float]:
    key = fitness.value
    return {
        "p3":               _FITNESS_P3_MULTIPLIERS.get(key,  1.0),
        "secretion_scale":  _FITNESS_SEC_MULTIPLIERS.get(key, 1.0),
    }
