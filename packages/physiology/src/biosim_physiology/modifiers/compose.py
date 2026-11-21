# SAVE TO: packages\physiology\src\biosim_physiology\modifiers\compose.py
"""Multiplicative composition of modifier dicts."""

from __future__ import annotations


# ---------------------------------------------------------------------------
# Composite modifier
# ---------------------------------------------------------------------------

def compose_all_modifiers(modifiers: list[dict[str, float]]) -> dict[str, float]:
    """
    Multiply all modifier dicts together to get net per-parameter factor.

    Parameters not mentioned in any modifier keep a factor of 1.0
    (i.e. they appear in the output only if at least one modifier affects them).
    """
    combined: dict[str, float] = {}
    for m in modifiers:
        for param, factor in m.items():
            combined[param] = combined.get(param, 1.0) * factor
    return combined
