from __future__ import annotations
import math

def clamp(value, lower=0.0, upper=float("inf")):
    # NaN comparisons are always False, so max(lower, min(upper, nan))
    # silently returns `upper` unchanged -- the opposite of a safe floor
    # (found during the meal-absorption mass-conservation fix's security
    # review: every pathway using clamp() as its sole defense against a
    # diverging/malformed state relied on this, incorrectly). Map NaN to
    # `lower` explicitly, matching every other pathway's intent when it
    # calls clamp(y[i], lower=0.0) to guard against a nonsensical state.
    if math.isnan(value):
        return lower
    return max(lower, min(upper, value))

def michaelis_menten(substrate, vmax, km):
    if substrate < 0: raise ValueError("substrate must be non-negative")
    if km <= 0: raise ValueError("km must be positive")
    return vmax * max(0.0, substrate) / (km + max(0.0, substrate))

def michaelis_menten_inhibited(substrate, vmax, km, inhibitor=0.0, ki=1.0):
    if substrate < 0: raise ValueError("substrate must be non-negative")
    if km <= 0: raise ValueError("km must be positive")
    if ki <= 0: raise ValueError("ki must be positive")
    denom = km * (1.0 + inhibitor / max(ki, 1e-12)) + max(0.0, substrate)
    return vmax * max(0.0, substrate) / denom

def hill_activation(x, ec50, n):
    if x < 0: raise ValueError("x must be non-negative")
    if ec50 <= 0: raise ValueError("ec50 must be positive")
    if n <= 0: raise ValueError("n must be positive")
    x = max(0.0, x)
    return (x**n) / (ec50**n + x**n)

def hill_inhibition(x, ic50, n):
    if ic50 <= 0: raise ValueError("ic50 must be positive")
    x = max(0.0, x)
    if x == 0.0: return 1.0
    return (ic50**n) / (ic50**n + x**n)

def approach_to_target(current, target, rate):
    if rate < 0: raise ValueError("rate must be non-negative")
    return rate * (target - current)

def first_order(y, k):
    return -k * max(0.0, y)

def glucose_absorption_rate(t_min, dose_mM_per_kg, peak_time_min=30.0, kabs=0.060):
    if t_min < 0: return 0.0
    if t_min <= peak_time_min:
        return kabs * dose_mM_per_kg * (t_min / peak_time_min)
    return kabs * dose_mM_per_kg * math.exp(-kabs * (t_min - peak_time_min))

def gi_adjusted_kabs(glycaemic_index, base_kabs=0.060):
    # UNVALIDATED PLACEHOLDER slope -- see GI_SLOPE_PLACEHOLDER_NOTE in
    # gut_absorption.py. Only GI=100 is calibrated anywhere in the platform.
    gi = clamp(glycaemic_index, 10.0, 100.0)
    return base_kabs * (gi / 55.0)
