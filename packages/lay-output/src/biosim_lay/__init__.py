# SAVE TO: packages\lay-output\src\biosim_lay\__init__.py
"""
biosim_lay -- lay user output layer.

Public API:
    from biosim_lay import (
        FoodItem, FOOD_DB, get_food, compose_meal, search_foods,
        TrafficLight, RiskIndicator, RiskProfile, assess_risk,
        PlainLanguageReport, generate_report,
        RiskBand, RiskIndicatorScore, MetabolicRiskProfile, compute_risk_profile,
    )
"""
from .food_items import (
    FOOD_DB, FoodItem, compose_meal, foods_by_gi_category, get_food,
    gi_category_for_value, search_foods,
)
from .report import (
    PlainLanguageReport, RiskIndicator, RiskProfile,
    TrafficLight, assess_risk, generate_report,
)
from .risk_score import (
    MetabolicRiskProfile, RiskBand, RiskIndicatorScore, compute_risk_profile,
)

__all__ = [
    "FOOD_DB",
    "FoodItem",
    "MetabolicRiskProfile",
    "PlainLanguageReport",
    "RiskBand",
    "RiskIndicator",
    "RiskIndicatorScore",
    "RiskProfile",
    "TrafficLight",
    "assess_risk",
    "compose_meal",
    "compute_risk_profile",
    "foods_by_gi_category",
    "generate_report",
    "get_food",
    "gi_category_for_value",
    "search_foods",
]
__version__ = "0.2.0"
