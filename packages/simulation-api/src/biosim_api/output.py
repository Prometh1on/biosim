# SAVE TO: packages\simulation-api\src\biosim_api\output.py
"""
Simulation output types.

SimulationOutput is the single result type returned by SimulationRunner.
It contains both the raw simulation data (for experts) and a lay-readable
summary (for non-expert users).

Lay summary design:
    The phenotype translation layer maps mechanistic model outputs to
    human-scale experiences. This is intentionally approximate and is
    labelled as interpretive, not derived.

    Energy level proxy (0-10):
        Rises when blood glucose is above baseline (fuel available).
        Drops sharply when glucose falls below baseline (reactive dip).
        Modulated by insulin: high insulin = slight fatigue tendency.

    Hunger proxy (0-10):
        Falls when blood sugar is high (satiated).
        Rises when blood sugar returns toward or below baseline.

    Both proxies are annotated with their uncertainty: if the underlying
    simulation has wide confidence intervals, the lay summary reflects that.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional
import numpy as np


@dataclass
class ExpertMetrics:
    """
    Quantitative summary for researcher-facing output.
    All values in platform canonical units (uM, s).
    """
    glucose_peak_uM:           float
    glucose_peak_time_s:       float
    glucose_auc_uM_s:          float
    glucose_baseline_uM:       float
    insulin_peak_uM:           float
    insulin_peak_time_s:       float
    mtor_peak:                 float
    ampk_nadir:                float
    time_to_glucose_return_s:  Optional[float]   # None if doesn't return within duration
    overall_watchdog_confidence: float
    n_ode_evaluations:         int
    # Incremental glucose AUC over the first 5h (0-18000s), clipped at 0
    # below baseline -- the standard clinical "iAUC(0-5h)" metric used in
    # real postprandial drug trials (e.g. Hjerpsted et al. 2018,
    # PMID:28941314, for semaglutide). Distinct from glucose_auc_uM_s
    # (total, unclipped AUC over the full simulated window) -- this one
    # is fixed to a clinically standard 5h window. 0.0 if the simulation
    # ran for less than 5h.
    glucose_iauc_5h_uM_s:      float = 0.0
    # Same idea as glucose_iauc_5h_uM_s, fixed to a 4h window instead --
    # the clinical standard for a "mixed meal test" endpoint (e.g. Hansen
    # et al. 2024, PMID:39049802, for metformin: a 4h mixed meal test).
    # 0.0 if the simulation ran for less than 4h. A separate field, not a
    # parameterized method, to match glucose_iauc_5h_uM_s's own existing
    # shape/convention exactly.
    glucose_iauc_4h_uM_s:      float = 0.0
    # Lowest glucose value reached anywhere in the simulated trajectory --
    # added specifically to feed a lay-output hypoglycaemia-risk indicator
    # (see biosim_lay.report/risk_score) once real glucose-independent-
    # secretion (sulfonylurea) and insulin-independent (SGLT2i) mechanisms
    # existed that could genuinely push simulated glucose low, with no
    # user-facing signal ever checking for it -- the same "computed a
    # real mechanism, never wired a user-facing signal" bug class already
    # caught and fixed several times in this codebase. 0.0 default matches
    # the failed-simulation branch's convention for every other *_uM
    # field (glucose_peak_uM=0, etc.) -- _compute_expert_metrics() always
    # sets this explicitly for a successful run.
    glucose_nadir_uM:          float = 0.0

    @property
    def glucose_excursion_uM(self) -> float:
        """Peak glucose above baseline."""
        return self.glucose_peak_uM - self.glucose_baseline_uM

    @property
    def glucose_peak_mM(self) -> float:
        return self.glucose_peak_uM / 1000.0

    @property
    def glucose_nadir_mM(self) -> float:
        return self.glucose_nadir_uM / 1000.0

    @property
    def insulin_peak_pM(self) -> float:
        return self.insulin_peak_uM * 1e6

    @property
    def glucose_iauc_5h_mmol_h_per_L(self) -> float:
        """glucose_iauc_5h_uM_s converted to mmol.h/L, the unit real
        clinical trials (e.g. Hjerpsted et al. 2018) report iAUC in."""
        return self.glucose_iauc_5h_uM_s / 1000.0 / 3600.0

    @property
    def glucose_iauc_4h_mmol_min_per_L(self) -> float:
        """glucose_iauc_4h_uM_s converted to mmol.min/L, the unit Hansen
        et al. 2024 (PMID:39049802) reports metformin's postprandial
        glucose-excursion iAUC in."""
        return self.glucose_iauc_4h_uM_s / 1000.0 / 60.0


@dataclass
class LaySummary:
    """
    Human-readable summary of a simulation.
    Intentionally approximate — labelled as interpretive output.
    """
    times_min:              np.ndarray    # time axis in minutes
    energy_level:           np.ndarray    # 0-10 proxy score
    hunger_level:           np.ndarray    # 0-10 proxy score
    energy_peak_min:        float         # time of peak energy
    energy_dip_min:         Optional[float]  # time of energy dip (if any)
    has_energy_dip:         bool
    glucose_narrative:      str           # plain English description
    energy_narrative:       str
    confidence_note:        str           # uncertainty caveat
    is_high_confidence:     bool

    @property
    def energy_at_peak(self) -> float:
        return float(np.max(self.energy_level))

    @property
    def energy_at_dip(self) -> Optional[float]:
        if self.energy_dip_min is None:
            return None
        idx = int(np.argmin(np.abs(self.times_min - self.energy_dip_min)))
        return float(self.energy_level[idx])


@dataclass
class SimulationOutput:
    """
    Complete output of a SimulationRunner.run() call.

    Contains the full simulation result, personalisation provenance,
    expert metrics, and lay summary.
    """
    scenario_name:                str
    scenario_description:         str
    simulation:                   object    # biosim_ode.SimulationResult
    expert_metrics:               ExpertMetrics
    lay_summary:                  LaySummary
    uq:                           Optional[object] = None   # UQResult if run
    physiological_personalization: Optional[object] = None
    genetic_personalization:      Optional[object] = None
    warnings:                     list[str] = field(default_factory=list)
    # The ScenarioDefinition that produced this output (biosim_api.ScenarioDefinition).
    # None only for hand-constructed SimulationOutputs (e.g. some test
    # fixtures) that don't set it -- SimulationRunner.run() always does.
    # Added so ReproducibilityBundle.create_bundle() can capture the full
    # scenario (meals/drugs/personalization), not just its name/description
    # (all output.scenario_name/scenario_description alone ever carried).
    scenario:                     Optional[object] = None

    @property
    def success(self) -> bool:
        return self.simulation.success

    @property
    def has_personalization(self) -> bool:
        return (self.physiological_personalization is not None
                or self.genetic_personalization is not None)

    def summary(self) -> dict:
        return {
            "scenario":          self.scenario_name,
            "success":           self.success,
            "glucose_peak_mM":   round(self.expert_metrics.glucose_peak_mM, 2),
            "insulin_peak_pM":   round(self.expert_metrics.insulin_peak_pM, 1),
            "mtor_peak":         round(self.expert_metrics.mtor_peak, 3),
            "ampk_nadir":        round(self.expert_metrics.ampk_nadir, 3),
            "confidence":        round(self.expert_metrics.overall_watchdog_confidence, 3),
            "has_energy_dip":    self.lay_summary.has_energy_dip,
            "energy_peak_min":   round(self.lay_summary.energy_peak_min, 1),
            "n_warnings":        len(self.warnings),
            "personalised":      self.has_personalization,
        }
