# SAVE TO: packages\physiology\src\biosim_physiology\modifiers\__init__.py
"""
Physiological parameter modifiers.

Each modifier is a pure function: PhysiologicalProfile → dict of
parameter_name → multiplier. Multipliers compose multiplicatively.
A multiplier of 1.0 means no change.

Literature sources (verified against PubMed E-utilities -- see CLAUDE.md
KNOWN BUGS for the citation-accuracy audit; every PMID below was previously
fabricated/mismatched and has been corrected):
    Age/insulin sensitivity:   DeFronzo 1981, Diabetes Care. PMID:7049632
    Sex differences:           Mauvais-Jarvis 2015, Biol Sex Differ. PMID:26339468
    BMI/insulin resistance:    Kahn 2006, Nature. PMID:17167471
    Fitness/insulin action:    Dela & Kjaer 2006, Essays Biochem. PMID:17144881
    T2D/beta cell loss:        Stumvoll 2005, Lancet. PMID:15823385
    Hypothyroidism/glucose:    Brenta 2011, J Thyroid Res. PMID:21941681
                               Handisurya 2008, Clin Endocrinol. PMID:18429948
                               Rochon 2003, Clin Sci. PMID:12519082
    Hyperthyroidism/glucose:   Dimitriadis & Raptis 2001, Exp Clin Endocrinol
                               Diabetes. PMID:11460573 (peripheral disposal
                               only -- the hepatic/postprandial mechanism is
                               represented via EGPParameters.for_hyperthyroidism()
                               in biosim_ode/hepatic.py, not a modifier here)

All modifiers operate on PathwayParameters fields.
Multipliers < 1.0 reduce the parameter; > 1.0 increase it.

Package layout (split from one 1197-line module, 2026-10-20; every name
is re-exported here, so `from biosim_physiology.modifiers import X` keeps
working):
    demographics   age, sex, BMI, fitness
    lifestyle      sleep, intermittent hypoxia, menstrual-cycle phase,
                   cold acclimation
    disease        disease-state modifiers
    compose        compose_all_modifiers
"""

from __future__ import annotations

from .demographics import (  # noqa: F401 -- re-exported, see module docstring
    AGE_REFERENCE,
    AGE_P3_DECLINE_PER_YEAR,
    AGE_SEC_DECLINE_PER_YEAR,
    AGE_P3_FLOOR,
    AGE_SEC_FLOOR,
    age_modifiers,
    sex_modifiers,
    bmi_modifiers,
    _FITNESS_P3_MULTIPLIERS,
    _FITNESS_SEC_MULTIPLIERS,
    fitness_modifiers,
)
from .lifestyle import (  # noqa: F401 -- re-exported, see module docstring
    SLEEP_REPLETE_HOURS,
    SLEEP_RESTRICTED_HOURS,
    SLEEP_RESTRICTED_P3_FACTOR,
    sleep_modifiers,
    HYPOXIA_REFERENCE_EVENTS_PER_HOUR,
    HYPOXIA_P3_FACTOR_AT_REFERENCE,
    HYPOXIA_P1_FACTOR_AT_REFERENCE,
    intermittent_hypoxia_modifiers,
    MENSTRUAL_CYCLE_LUTEAL_P3_FACTOR,
    menstrual_cycle_modifiers,
    COLD_ACCLIMATION_P3_FACTOR,
    cold_acclimation_modifiers,
)
from .disease import (  # noqa: F401 -- re-exported, see module docstring
    disease_modifiers,
    _single_disease_modifier,
)
from .compose import (  # noqa: F401 -- re-exported, see module docstring
    compose_all_modifiers,
)
