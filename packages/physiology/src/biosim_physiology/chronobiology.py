# SAVE TO: packages\physiology\src\biosim_physiology\chronobiology.py
"""
Circadian effects on glucose/insulin dynamics.

Insulin sensitivity varies by ~20-30% across the 24-hour day,
driven by the circadian clock in the suprachiasmatic nucleus and
peripheral clocks in liver, muscle, and adipose tissue.

Model (sinusoidal approximation):
    circadian_factor(t) = 1 + AMPLITUDE * cos(2π(t - PEAK_HOUR) / 24)

    AMPLITUDE = 0.15  (15% variation around mean)
    PEAK_HOUR = 10.0  (peak sensitivity at 10:00)

    t = 10  → factor = 1.15  (most sensitive, morning)
    t = 22  → factor = 0.85  (least sensitive, late evening)

This is applied to p3 (insulin action) and p1 (glucose effectiveness).
The effect on secretion_scale is smaller (~5% amplitude).

Literature:
    Van Cauter et al. 1992, Am J Physiol. PMID:1566835
    Scheer et al. 2009, PNAS. PMID:19255424
    Saad et al. 2012, Diabetes. PMID:22751690
"""

from __future__ import annotations

import math
from typing import Optional


CIRCADIAN_P3_AMPLITUDE    = 0.15   # 15% amplitude for insulin sensitivity
CIRCADIAN_P1_AMPLITUDE    = 0.08   # 8% amplitude for glucose effectiveness
CIRCADIAN_SEC_AMPLITUDE   = 0.05   # 5% amplitude for secretion
CIRCADIAN_PEAK_HOUR       = 10.0   # peak sensitivity at 10:00 (hours since midnight)


def circadian_factor(time_of_day_hours: float, amplitude: float, peak_hour: float = CIRCADIAN_PEAK_HOUR) -> float:
    """
    Compute the circadian multiplier for a given time of day.

    Args:
        time_of_day_hours: Current time (0.0 = midnight, 12.0 = noon, 23.99 = just before midnight)
        amplitude: Amplitude of sinusoidal variation (e.g. 0.15 = ±15%)
        peak_hour: Hour of maximum (e.g. 10.0 = 10 AM)

    Returns:
        Multiplicative factor. 1.0 at the crossover points,
        (1 + amplitude) at the peak, (1 - amplitude) at the trough.
    """
    phase = 2.0 * math.pi * (time_of_day_hours - peak_hour) / 24.0
    return 1.0 + amplitude * math.cos(phase)


def circadian_modifiers(time_of_day_hours: Optional[float]) -> dict[str, float]:
    """
    Return circadian multipliers for all affected parameters.

    Args:
        time_of_day_hours: Time of day (0–24). None = no adjustment.

    Returns:
        Dict of parameter_name → multiplier.
        Empty dict if time_of_day_hours is None.
    """
    if time_of_day_hours is None:
        return {}

    return {
        "p3": circadian_factor(
            time_of_day_hours, CIRCADIAN_P3_AMPLITUDE, CIRCADIAN_PEAK_HOUR
        ),
        "p1": circadian_factor(
            time_of_day_hours, CIRCADIAN_P1_AMPLITUDE, CIRCADIAN_PEAK_HOUR
        ),
        "secretion_scale": circadian_factor(
            time_of_day_hours, CIRCADIAN_SEC_AMPLITUDE, CIRCADIAN_PEAK_HOUR
        ),
    }


def time_description(hour: float) -> str:
    """Human-readable description of circadian state at a given hour."""
    if 6 <= hour < 12:   return "morning (good insulin sensitivity)"
    if 12 <= hour < 16:  return "midday (moderate insulin sensitivity)"
    if 16 <= hour < 20:  return "afternoon/evening (declining sensitivity)"
    if 20 <= hour < 24:  return "night (lowest insulin sensitivity)"
    return "early morning (lowest insulin sensitivity)"


def circadian_sensitivity_percent(time_of_day_hours: float) -> float:
    """
    Return insulin sensitivity as a percentage of the daily mean at this time.
    100% = daily mean, >100% = more sensitive, <100% = less sensitive.
    """
    factor = circadian_factor(time_of_day_hours, CIRCADIAN_P3_AMPLITUDE)
    return factor * 100.0
