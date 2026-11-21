# SAVE TO: packages\uncertainty\src\biosim_uncertainty\__init__.py
"""
biosim_uncertainty -- uncertainty quantification and sensitivity analysis.

Public API:
    from biosim_uncertainty import (
        MonteCarloEngine, OATSensitivityAnalyzer,
        SobolSensitivityAnalyzer, SobolResult,
        UQResult, ConfidenceBand,
        LogNormalDistribution, NormalDistribution,
        UniformDistribution, FixedValue,
        STANDARD_METRICS,
        prioritise_from_sensitivity, PARAM_BIOLOGY_MAP,
    )
"""
from .calibration import (
    BayesianCalibrator, CalibrationResult, CalibrationTarget,
    CALIBRATION_METRIC_EXTRACTORS,
)
from .distributions import (
    FixedValue, LogNormalDistribution, NormalDistribution,
    ParameterDistribution, UniformDistribution,
)
from .propagation import MonteCarloEngine
from .result import ConfidenceBand, UQResult
from .sensitivity import (
    OATSensitivityAnalyzer, SensitivityResult, STANDARD_METRICS,
    glucose_auc_metric, glucose_peak_metric, insulin_peak_metric,
)
from .sobol_sensitivity import SobolSensitivityAnalyzer, SobolResult

try:
    from .curation_bridge import (
        prioritise_from_sensitivity, PARAM_BIOLOGY_MAP, BridgeResult,
    )
except ImportError:
    pass   # curation_bridge.py not yet placed

__all__ = [
    "BayesianCalibrator",
    "CALIBRATION_METRIC_EXTRACTORS",
    "CalibrationResult",
    "CalibrationTarget",
    "ConfidenceBand",
    "FixedValue",
    "LogNormalDistribution",
    "MonteCarloEngine",
    "NormalDistribution",
    "OATSensitivityAnalyzer",
    "ParameterDistribution",
    "SensitivityResult",
    "SobolResult",
    "SobolSensitivityAnalyzer",
    "STANDARD_METRICS",
    "UQResult",
    "UniformDistribution",
    "glucose_auc_metric",
    "glucose_peak_metric",
    "insulin_peak_metric",
]
__version__ = "0.2.0"
