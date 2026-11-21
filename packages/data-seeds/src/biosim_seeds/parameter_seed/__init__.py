# SAVE TO: packages\data-seeds\src\biosim_seeds\parameter_seed\__init__.py
"""
Literature-sourced kinetic parameters for the v1 glucose/insulin pathway.

Every parameter is linked to its primary literature source via PubMed ID.

Confidence-tier definition (corrected -- see the housekeeping-backlog
fix below): this file previously stated tiers reflected REPLICATION
COUNT ("GOLD = replicated in >=10 independent studies"), a bar no
single parameter in this file -- including the ones already marked
GOLD -- has ever actually met, and which directly contradicted this
package's own canonical `ConfidenceTier` definition
(`biosim_parameter_db.types.ConfidenceTier`: "GOLD = Expert-validated;
SILVER = Auto-extracted, plausibility-checked; BRONZE = Inferred or
estimated"). Checked every existing tier assignment in this file
against what it ACTUALLY tracks (not the stale docstring) before
correcting the docstring to match real, established practice:
    GOLD   = a primary source's own DIRECTLY reported/measured value
             for exactly this quantity (Bergman/Sherwin/
             Pacini's foundational minimal-model kinetics; DeFronzo
             2013/Al-Jobori 2017's directly-quoted TmG figures;
             Rebelato 2018's directly-quoted GSIS EC50; Esteban 1991's
             directly-quoted cortisol production rate) -- regardless
             of how many independent studies have replicated it. This
             also covers a RATIO computed directly from a single
             primary source's own two reported numbers, with no unit
             conversion and no back-solved engineering constant
             involved (CKD's covariate-adjusted p3_factor/ke_factor;
             tirzepatide's gb_factor_at_*mg, an ANCOVA-adjusted
             between-arm delta over a same-arm baseline; pioglitazone's
             p3_multiplier_at_45mg_reference_dose, a same-arm pre/post
             ISI ratio -- re-tiered from an earlier, inconsistent
             BRONZE assignment, MASTER_PLAN Section 8 housekeeping
             item, resolved here) -- distinct from BRONZE's back-solved
             engineering constants below, which convert a clinical
             ratio into a DIFFERENT unit this platform's own machinery
             needs (e.g. a kabs scale or a pM-per-gram rate).
    SILVER = from a review/secondary-compilation source (Le Roith
             1997's IGF-1 review, De Meyts 2002's receptor-affinity
             review), once removed from a primary measurement; OR a
             directly-quoted number from a clinical consensus
             diagnostic threshold/recommendation (IADPSG 2010's fasting-
             glucose GDM cutoff) rather than a direct physiological
             measurement of an individual/cohort -- a real number, not
             an estimate, but a different evidentiary category from a
             GOLD-tier measured value (added when gestational diabetes's
             Gb_factor needed a tier one step below its own p1/p3/
             secretion_scale factors, all three of which ARE directly-
             measured GOLD values from the same disease-state entity).
    BRONZE = an indirect estimate inferred from a primary source's
             data rather than the paper's own directly-reported
             statistic (Sarbassov 2005's mTOR EC50, "estimated from
             pAKT/mTOR... dose-response curves"; Hardie 2012's AMPK
             IC50, "indirect estimate from mTOR/AMPK cross-talk
             literature"), OR a disclosed engineering back-solve/
             calibration against a real measurement (the exogenous-
             insulin PK constants, the sulfonylurea/acarbose/miglitol
             kabs-scale constants) -- not itself a directly-reported
             value, matching `types.py`'s own "inferred or estimated".

Unit convention: all values stored in ODE engine internal units
    Concentrations : mM (glucose), pM (insulin), nM (IGF-1)
    Time           : minutes
    Rate constants : min⁻¹ or min⁻² per concentration unit

Literature sources (verified against PubMed E-utilities -- see CLAUDE.md
KNOWN BUGS for the citation-accuracy audit; every PMID below was previously
fabricated/mismatched and has been corrected):
    Bergman 1981   PMID:7033284    Original minimal model parameters
    Pacini 1986    PMID:3640682    Refined p2/p3 estimates (MINMOD)
    Sherwin 1974   PMID:4856884    Insulin elimination kinetics
    Sarbassov 2005 PMID:15718470   mTOR/insulin coupling
    Hardie 2012    PMID:22436748   AMPK regulation
    Le Roith 1997  PMID:9032050    IGF-1 baseline and kinetics
    Howey 1994     PMID:8314011    Exogenous insulin (regular + lispro) PK
    DeFronzo 2013  PMID:23735727   Dapagliflozin/renal GFR/TmG (SGLT2)
    Al-Jobori 2017 PMID:28428225   Empagliflozin renal TmG (SGLT2)
    Rebelato 2018  PMID:30166558   Beta-cell GSIS EC50 (KG_secretion)
    Meikle & Tyler 1977 PMID:888843   Glucocorticoid potency equivalence
    Esteban 1991   PMID:1986026    Daily cortisol production rate
    Riefflin 2015  PMID:25297572   Sulfonylurea glucose-independent secretion
    Chen 2016      PMID:27686734   Acarbose kabs-scale calibration
    Samad 1988     PMID:3042310    Miglitol kabs-scale calibration

Package layout (split from one 3310-line module, 2026-10-20; every name
is re-exported here, so `from biosim_seeds.parameter_seed import X` keeps
working):
    store            seed_parameter_store, create_seeded_store
    _helpers         _resolve, _resolve_context
    core_dynamics    glucose/insulin/IGF-1/mTOR/AMPK, beta-cell Km, NEFA,
                     protein-stimulated secretion
    drugs_dynamic    ODE-coupled / PK drugs and meal-absorption drugs
    drugs_static     static p3/Gb/secretion-factor drugs
    diabetes         T1D, T2D, insulin resistance, GDM, GCK-MODY
    diseases         other disease states
    lifestyle        sleep, intermittent hypoxia, cold, menstrual cycle
"""

from __future__ import annotations

from .store import (  # noqa: F401 -- re-exported, see module docstring
    seed_parameter_store,
    create_seeded_store,
)
from ._helpers import (  # noqa: F401 -- re-exported, see module docstring
    _resolve,
    _resolve_context,
)
from .core_dynamics import (  # noqa: F401 -- re-exported, see module docstring
    _seed_glucose_dynamics,
    _seed_insulin_dynamics,
    _seed_igf1_dynamics,
    _seed_mtor_dynamics,
    _seed_ampk_dynamics,
    _seed_beta_cell_secretion_km,
    _seed_nefa_pk,
    _seed_protein_secretion,
)
from .drugs_dynamic import (  # noqa: F401 -- re-exported, see module docstring
    _seed_exogenous_insulin_pk,
    _seed_renal_sglt2i,
    _seed_cortisol_pk,
    _seed_glp1_pk,
    _seed_prednisone_glucocorticoid,
    _seed_alcohol_egp,
    _seed_alpha_glucosidase_inhibitor_kabs,
)
from .drugs_static import (  # noqa: F401 -- re-exported, see module docstring
    _seed_hydroxychloroquine_p3,
    _seed_growth_hormone_p3,
    _seed_hydrochlorothiazide_p3,
    _seed_indinavir_p3,
    _seed_atorvastatin_p3,
    _seed_pioglitazone_p3,
    _seed_caffeine_p3,
    _seed_metformin_gb,
    _seed_tirzepatide_gb,
    _seed_sulfonylurea_secretion,
    _seed_glipizide_secretion,
)
from .diabetes import (  # noqa: F401 -- re-exported, see module docstring
    _seed_gck_mody_disease_state,
    _seed_gestational_diabetes,
    _seed_type_1_diabetes,
    _seed_type_2_diabetes,
    _seed_insulin_resistant,
)
from .diseases import (  # noqa: F401 -- re-exported, see module docstring
    _seed_cushings_cortisol_ratio,
    _seed_hyperthyroidism_egp,
    _seed_hypothyroidism,
    _seed_nafld,
    _seed_pcos,
    _seed_acromegaly,
    _seed_osa,
    _seed_ckd_disease_state,
)
from .lifestyle import (  # noqa: F401 -- re-exported, see module docstring
    _seed_cold_acclimation,
    _seed_menstrual_cycle,
    _seed_sleep_restriction,
    _seed_intermittent_hypoxia,
)
