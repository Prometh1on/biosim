# SAVE TO: packages\lay-output\src\biosim_lay\food_items.py
"""
Food item database -- maps common foods to carbohydrate content.

This is the input layer for the lay user interface. Instead of asking
users to enter '75g carbohydrate', they can select 'large banana + white rice'.

Carbohydrate values are approximate population averages from:
    USDA FoodData Central (https://fdc.nal.usda.gov/)
    UK Biobank food frequency questionnaire
    McCance & Widdowson's 'The Composition of Foods' (7th ed.)

All values are total available carbohydrate in grams.
Fibre is excluded (not absorbed as glucose in this model).

Usage:
    from biosim_lay.food_items import FoodItem, FOOD_DB, compose_meal
    meal = compose_meal(["white_rice_200g", "cola_330ml"])
    # Returns carbs_g = 56 + 35 = 91g
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class FoodItem:
    """A single food item with its carbohydrate content."""
    id:          str
    name:        str
    carbs_g:     float
    portion_desc: str    # e.g. "1 medium banana", "200g cooked"
    gi_category: str     # "high" (GI>70), "medium" (55-70), "low" (<55)

    def __post_init__(self) -> None:
        if self.carbs_g < 0:
            raise ValueError(f"carbs_g cannot be negative for {self.id}")
        if self.gi_category not in ("high", "medium", "low", "negligible"):
            raise ValueError(f"gi_category must be high/medium/low/negligible, got {self.gi_category}")


def gi_category_for_value(gi: float) -> str:
    """
    Canonical implementation of this module's own stated thresholds
    (see FoodItem.gi_category's docstring: "high" (GI>70), "medium"
    (55-70), "low" (<55)) -- exists so the numeric GI values in
    biosim_seeds.gi_meals.GI_VALUES (the table actually used for
    absorption-rate math) can be checked against FOOD_DB's qualitative
    labels by a single, unambiguous rule instead of by eye. See
    test_gi_category_matches_gi_meals_numeric_values in
    packages/data-seeds/tests/test_gi_meals.py.
    """
    if gi > 70.0:
        return "high"
    if gi >= 55.0:
        return "medium"
    return "low"


# ---------------------------------------------------------------------------
# Food database (USDA FoodData Central values)
# ---------------------------------------------------------------------------

FOOD_DB: dict[str, FoodItem] = {
    # Grains and starches
    "white_rice_200g":    FoodItem("white_rice_200g",    "White rice (200g cooked)",         56.0,  "200g cooked",         "high"),
    "brown_rice_200g":    FoodItem("brown_rice_200g",    "Brown rice (200g cooked)",          44.0,  "200g cooked",         "medium"),
    "white_bread_slice":  FoodItem("white_bread_slice",  "White bread (1 slice, 30g)",        14.0,  "1 slice (30g)",       "high"),
    "wholegrain_bread":   FoodItem("wholegrain_bread",   "Wholegrain bread (1 slice, 35g)",   12.0,  "1 slice (35g)",       "low"),
    "pasta_200g":         FoodItem("pasta_200g",         "Pasta (200g cooked)",               52.0,  "200g cooked",         "medium"),
    "potato_medium":      FoodItem("potato_medium",      "Baked potato (1 medium, 150g)",     30.0,  "1 medium baked",      "high"),
    "oats_40g":           FoodItem("oats_40g",           "Porridge oats (40g dry)",           27.0,  "40g dry oats",        "medium"),

    # Fruits
    "banana_medium":      FoodItem("banana_medium",      "Banana (1 medium, 120g)",           25.0,  "1 medium banana",     "low"),
    "apple_medium":       FoodItem("apple_medium",       "Apple (1 medium, 150g)",            20.0,  "1 medium apple",      "low"),
    "orange_medium":      FoodItem("orange_medium",      "Orange (1 medium, 130g)",           15.0,  "1 medium orange",     "low"),
    "grapes_100g":        FoodItem("grapes_100g",        "Grapes (100g)",                     17.0,  "100g",                "low"),
    "mango_half":         FoodItem("mango_half",         "Mango (half, 150g)",                25.0,  "half mango (150g)",   "medium"),

    # Drinks
    "cola_330ml":         FoodItem("cola_330ml",         "Cola / fizzy drink (330ml can)",    35.0,  "330ml can",           "medium"),
    "orange_juice_200ml": FoodItem("orange_juice_200ml", "Orange juice (200ml glass)",        20.0,  "200ml glass",         "low"),
    "sports_drink_500ml": FoodItem("sports_drink_500ml", "Sports drink (500ml)",              33.0,  "500ml bottle",        "high"),
    "milk_200ml":         FoodItem("milk_200ml",         "Semi-skimmed milk (200ml)",         10.0,  "200ml glass",         "low"),
    "water_any":          FoodItem("water_any",          "Water (any amount)",                 0.0,  "any amount",          "negligible"),

    # Snacks and sweets
    "chocolate_40g":      FoodItem("chocolate_40g",      "Milk chocolate (40g, ~4 squares)",  22.0,  "40g bar",             "high"),
    "biscuit_digestive":  FoodItem("biscuit_digestive",  "Digestive biscuit (1, 14g)",        10.0,  "1 biscuit",           "medium"),
    "crisps_30g":         FoodItem("crisps_30g",         "Crisps/chips (30g packet)",         16.0,  "30g packet",          "high"),
    "cereal_bar_40g":     FoodItem("cereal_bar_40g",     "Cereal/granola bar (40g)",          27.0,  "1 bar (40g)",         "medium"),

    # Meals
    "pizza_slice":        FoodItem("pizza_slice",        "Pizza slice (1 slice, ~100g)",      28.0,  "1 slice (~100g)",     "high"),
    "burger_bun":         FoodItem("burger_bun",         "Burger bun only",                   25.0,  "1 standard bun",      "high"),
    "sandwich_white":     FoodItem("sandwich_white",     "White bread sandwich (2 slices)",   28.0,  "2 slices white bread","high"),

    # Low-carb / protein
    "egg_boiled":         FoodItem("egg_boiled",         "Boiled egg",                         0.5,  "1 large egg",         "negligible"),
    "chicken_breast":     FoodItem("chicken_breast",     "Grilled chicken breast (150g)",      0.0,  "150g cooked",         "negligible"),
    "salad_plain":        FoodItem("salad_plain",        "Plain green salad",                  2.0,  "1 bowl",              "negligible"),

    # Glucose test standard
    "glucose_75g_ogtt":   FoodItem("glucose_75g_ogtt",  "75g glucose (OGTT standard)",       75.0,  "75g dissolved glucose","high"),
}


def get_food(food_id: str) -> FoodItem:
    """Retrieve a food item by ID. Raises KeyError if not found."""
    if food_id not in FOOD_DB:
        raise KeyError(
            f"Food '{food_id}' not in database. "
            f"Available: {sorted(FOOD_DB.keys())}"
        )
    return FOOD_DB[food_id]


def compose_meal(food_ids: list[str]) -> float:
    """
    Sum the carbohydrate content of a list of food items.

    Args:
        food_ids: List of food item IDs from FOOD_DB

    Returns:
        Total carbohydrates in grams
    """
    return sum(get_food(fid).carbs_g for fid in food_ids)


def search_foods(query: str) -> list[FoodItem]:
    """
    Simple text search across food names and IDs (case-insensitive).

    Args:
        query: Search string

    Returns:
        List of matching FoodItem objects, sorted by carbs_g
    """
    q = query.lower()
    matches = [
        item for item in FOOD_DB.values()
        if q in item.name.lower() or q in item.id.lower()
    ]
    return sorted(matches, key=lambda f: f.carbs_g)


def foods_by_gi_category(category: str) -> list[FoodItem]:
    """Return all foods in a given GI category (high/medium/low/negligible)."""
    return [f for f in FOOD_DB.values() if f.gi_category == category]
