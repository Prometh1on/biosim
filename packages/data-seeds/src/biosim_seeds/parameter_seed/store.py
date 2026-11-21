# SAVE TO: packages\data-seeds\src\biosim_seeds\parameter_seed\store.py
"""Seeding entry points: populate or create a ParameterStore."""

from __future__ import annotations

from typing import Optional

from biosim_context_registry import ContextRegistry
from biosim_id_resolver import IDRegistry
from biosim_parameter_db import (
    ParameterStore,
)

from ..context_seed import create_seeded_context_registry
from .core_dynamics import (
    _seed_ampk_dynamics,
    _seed_beta_cell_secretion_km,
    _seed_glucose_dynamics,
    _seed_igf1_dynamics,
    _seed_insulin_dynamics,
    _seed_mtor_dynamics,
    _seed_nefa_pk,
    _seed_protein_secretion,
)
from .diabetes import (
    _seed_gck_mody_disease_state,
    _seed_gestational_diabetes,
    _seed_insulin_resistant,
    _seed_type_1_diabetes,
    _seed_type_2_diabetes,
)
from .diseases import (
    _seed_acromegaly,
    _seed_ckd_disease_state,
    _seed_cushings_cortisol_ratio,
    _seed_hyperthyroidism_egp,
    _seed_hypothyroidism,
    _seed_nafld,
    _seed_osa,
    _seed_pcos,
)
from .drugs_dynamic import (
    _seed_alcohol_egp,
    _seed_alpha_glucosidase_inhibitor_kabs,
    _seed_cortisol_pk,
    _seed_exogenous_insulin_pk,
    _seed_glp1_pk,
    _seed_prednisone_glucocorticoid,
    _seed_renal_sglt2i,
)
from .drugs_static import (
    _seed_atorvastatin_p3,
    _seed_caffeine_p3,
    _seed_glipizide_secretion,
    _seed_growth_hormone_p3,
    _seed_hydrochlorothiazide_p3,
    _seed_hydroxychloroquine_p3,
    _seed_indinavir_p3,
    _seed_metformin_gb,
    _seed_pioglitazone_p3,
    _seed_sulfonylurea_secretion,
    _seed_tirzepatide_gb,
)
from .lifestyle import (
    _seed_cold_acclimation,
    _seed_intermittent_hypoxia,
    _seed_menstrual_cycle,
    _seed_sleep_restriction,
)


def seed_parameter_store(
    store: ParameterStore,
    registry: IDRegistry,
    context_registry: Optional[ContextRegistry] = None,
) -> ParameterStore:
    """
    Populate a ParameterStore with literature-sourced kinetic parameters.

    Args:
        store:            ParameterStore to populate (modified in place)
        registry:         Seeded IDRegistry for MOLECULAR entity resolution
        context_registry: Seeded ContextRegistry for non-molecular (disease-
                           state/physiological-axis) subject resolution.
                           Defaults to a freshly-seeded one if not provided,
                           so every existing single-argument call site
                           (create_seeded_store(registry)) is unaffected.

    Returns:
        The same store, now populated with v1 parameters.
    """
    context_registry = context_registry or create_seeded_context_registry()
    _seed_glucose_dynamics(store, registry)
    _seed_insulin_dynamics(store, registry)
    _seed_igf1_dynamics(store, registry)
    _seed_mtor_dynamics(store, registry)
    _seed_ampk_dynamics(store, registry)
    _seed_exogenous_insulin_pk(store, registry)
    _seed_renal_sglt2i(store, registry)
    _seed_beta_cell_secretion_km(store, registry)
    _seed_cortisol_pk(store, registry)
    _seed_sulfonylurea_secretion(store, registry)
    _seed_alpha_glucosidase_inhibitor_kabs(store, registry)
    _seed_pioglitazone_p3(store, registry)
    _seed_protein_secretion(store, registry)
    _seed_ckd_disease_state(store, context_registry)
    _seed_sleep_restriction(store, context_registry)
    _seed_gck_mody_disease_state(store, context_registry)
    _seed_gestational_diabetes(store, context_registry)
    _seed_type_1_diabetes(store, context_registry)
    _seed_type_2_diabetes(store, context_registry)
    _seed_insulin_resistant(store, context_registry)
    _seed_hypothyroidism(store, context_registry)
    _seed_metformin_gb(store, registry)
    _seed_caffeine_p3(store, registry)
    _seed_tirzepatide_gb(store, registry)
    _seed_prednisone_glucocorticoid(store, registry)
    _seed_glipizide_secretion(store, registry)
    _seed_intermittent_hypoxia(store, context_registry)
    _seed_glp1_pk(store, registry)
    _seed_nefa_pk(store, registry)
    _seed_growth_hormone_p3(store, registry)
    _seed_hydrochlorothiazide_p3(store, registry)
    _seed_indinavir_p3(store, registry)
    _seed_atorvastatin_p3(store, registry)
    _seed_cold_acclimation(store, context_registry)
    _seed_menstrual_cycle(store, context_registry)
    _seed_pcos(store, context_registry)
    _seed_acromegaly(store, context_registry)
    _seed_osa(store, context_registry)
    _seed_cushings_cortisol_ratio(store, registry)
    _seed_hyperthyroidism_egp(store, context_registry)
    _seed_alcohol_egp(store, registry)
    _seed_hydroxychloroquine_p3(store, registry)
    _seed_nafld(store, context_registry)
    return store


def create_seeded_store(
    registry: IDRegistry, context_registry: Optional[ContextRegistry] = None,
) -> ParameterStore:
    """Create and return a new ParameterStore seeded with all v1 parameters."""
    store = ParameterStore()
    return seed_parameter_store(store, registry, context_registry)
