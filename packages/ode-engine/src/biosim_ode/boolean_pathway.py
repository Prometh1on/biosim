# SAVE TO: packages\ode-engine\src\biosim_ode\boolean_pathway.py
"""
Boolean watchdog pathway for the glucose/insulin axis.

Implements the confidence watchdog architecture from Sprint 0:
    ODE is primary — all quantitative outputs come from the ODE.
    Boolean model runs in parallel as a consistency checker.
    Disagreement widens confidence intervals; it never blends outputs.

How the Boolean model works:
    Each species has a threshold that defines a binary HIGH/LOW state.
    Consistency rules encode expected correlations between states
    (e.g. high insulin should lead to active mTOR).
    If the ODE produces a state that violates a consistency rule,
    the confidence score decreases and uncertainty intervals widen.

This catches pathological parameter combinations that the ODE might
'solve' numerically but that are biologically implausible:
    - Insulin very high but mTOR completely inactive (broken AKT signalling)
    - Glucose very high but insulin not rising (broken beta cell response)
    - mTOR active but AMPK also active (contradictory energy states)

Thresholds are in internal ODE units (mM, pM, dimensionless).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional
import numpy as np

from .system import WatchdogResult


@dataclass(frozen=True)
class BooleanThresholds:
    """
    Threshold values for converting ODE state to Boolean (HIGH/LOW).
    All in internal ODE units.
    """
    G_high_mM:           float = 8.0     # Blood glucose HIGH if > 8 mM (~144 mg/dL)
    I_high_pM:           float = 150.0   # Insulin HIGH if > 150 pM (~25 μU/mL)
    mTOR_active:         float = 0.5     # mTOR ACTIVE if fractional activity > 0.5
    AMPK_suppressed:     float = 0.4     # AMPK SUPPRESSED if activity < 0.4
    IGF1_elevated_nM:    float = 20.0    # IGF-1 ELEVATED if > 20 nM


@dataclass
class BooleanState:
    """Boolean state of all tracked species at one time point."""
    time_min: float
    G_high: bool
    I_high: bool
    mTOR_active: bool
    AMPK_suppressed: bool
    IGF1_elevated: bool

    def as_dict(self) -> dict[str, bool]:
        return {
            "G_high":        self.G_high,
            "I_high":        self.I_high,
            "mTOR_active":   self.mTOR_active,
            "AMPK_suppressed": self.AMPK_suppressed,
            "IGF1_elevated": self.IGF1_elevated,
        }


@dataclass(frozen=True)
class ConsistencyRule:
    """
    A logical consistency rule between two Boolean states.

    If condition is True, then expected must also be True.
    Violations indicate the ODE is in a biologically implausible state.
    """
    condition: str      # Boolean state name that triggers this check
    expected: str       # Boolean state name that should be True if condition is True
    description: str    # Human-readable explanation
    severity: float     # How much a violation reduces confidence (0.0 to 1.0)


# Consistency rules for the glucose/insulin axis
GLUCOSE_INSULIN_RULES: list[ConsistencyRule] = [
    ConsistencyRule(
        condition="I_high",
        expected="mTOR_active",
        description="High insulin should activate mTOR via AKT",
        severity=0.25,
    ),
    ConsistencyRule(
        condition="I_high",
        expected="AMPK_suppressed",
        description="High insulin should suppress AMPK",
        severity=0.25,
    ),
    ConsistencyRule(
        condition="G_high",
        expected="I_high",
        description="Elevated glucose should trigger insulin secretion",
        severity=0.30,
    ),
]


class GlucoseInsulinBooleanPathway:
    """
    Boolean model of the glucose/insulin/mTOR/AMPK axis.

    Computes discrete HIGH/LOW states from ODE outputs and checks
    them against physiological consistency rules.
    """

    def __init__(
        self,
        thresholds: Optional[BooleanThresholds] = None,
        rules: Optional[list[ConsistencyRule]] = None,
    ) -> None:
        self.thresholds = thresholds or BooleanThresholds()
        self.rules = rules or GLUCOSE_INSULIN_RULES

    def compute_state(
        self,
        y_internal: list[float],
        t_min: float = 0.0,
    ) -> BooleanState:
        """
        Compute Boolean state from ODE state vector.

        Args:
            y_internal: [G_mM, I_pM, X, IGF1_nM, mTOR, AMPK]
            t_min: Current time in minutes
        """
        G, I, X, IGF1, mTOR, AMPK = y_internal[:6]
        th = self.thresholds
        return BooleanState(
            time_min=t_min,
            G_high=G > th.G_high_mM,
            I_high=I > th.I_high_pM,
            mTOR_active=mTOR > th.mTOR_active,
            AMPK_suppressed=AMPK < th.AMPK_suppressed,
            IGF1_elevated=IGF1 > th.IGF1_elevated_nM,
        )

    def check_consistency(self, state: BooleanState) -> list[ConsistencyRule]:
        """
        Return all consistency rules that are violated at this state.
        """
        state_dict = state.as_dict()
        violations = []
        for rule in self.rules:
            if state_dict.get(rule.condition, False):
                if not state_dict.get(rule.expected, False):
                    violations.append(rule)
        return violations

    def confidence_score(self, violations: list[ConsistencyRule]) -> float:
        """
        Compute confidence score from violations.
        Each violation reduces the score by its severity.
        Floor at 0.0.
        """
        penalty = sum(r.severity for r in violations)
        return max(0.0, 1.0 - penalty)


class WatchdogEngine:
    """
    Runs the Boolean watchdog alongside ODE output to produce
    confidence scores and divergence flags.

    Replaces the basic watchdog in system.py with a proper
    rule-based implementation.
    """

    def __init__(
        self,
        boolean_pathway: Optional[GlucoseInsulinBooleanPathway] = None,
    ) -> None:
        self.boolean = boolean_pathway or GlucoseInsulinBooleanPathway()

    def evaluate(
        self,
        y_internal: list[float],
        t_min: float,
    ) -> WatchdogResult:
        """
        Evaluate watchdog at a single time point.

        Args:
            y_internal: ODE state vector in internal units
            t_min: Current time in minutes

        Returns:
            WatchdogResult with confidence score and any disagreements
        """
        state = self.boolean.compute_state(y_internal, t_min)
        violations = self.boolean.check_consistency(state)
        score = self.boolean.confidence_score(violations)

        G, I, X, IGF1, mTOR, AMPK = y_internal[:6]

        return WatchdogResult(
            confidence_score=score,
            disagreements=[v.description for v in violations],
            ode_states={"G": G, "I": I, "mTOR": mTOR, "AMPK": AMPK},
            boolean_states=state.as_dict(),
        )

    def evaluate_trajectory(
        self,
        t_min_array: np.ndarray,
        y_internal_array: np.ndarray,
    ) -> list[WatchdogResult]:
        """
        Evaluate watchdog across a full simulation trajectory.

        Args:
            t_min_array: Time array in minutes, shape (n_times,)
            y_internal_array: State array in internal units, shape (n_species, n_times)

        Returns:
            List of WatchdogResult, one per time point
        """
        results = []
        for i in range(len(t_min_array)):
            y = y_internal_array[:, i].tolist()
            results.append(self.evaluate(y, float(t_min_array[i])))
        return results

    def overall_confidence(self, results: list[WatchdogResult]) -> float:
        """Mean confidence score across all time points."""
        if not results:
            return 1.0
        return float(np.mean([r.confidence_score for r in results]))

    def critical_periods(
        self,
        results: list[WatchdogResult],
        times_min: np.ndarray,
        threshold: float = 0.7,
    ) -> list[tuple[float, float]]:
        """
        Return time intervals where confidence falls below threshold.
        Returns list of (start_min, end_min) tuples.
        """
        periods = []
        in_period = False
        start = 0.0
        for i, (r, t) in enumerate(zip(results, times_min)):
            if r.confidence_score < threshold and not in_period:
                in_period = True
                start = float(t)
            elif r.confidence_score >= threshold and in_period:
                in_period = False
                periods.append((start, float(t)))
        if in_period:
            periods.append((start, float(times_min[-1])))
        return periods
