# SAVE TO: packages\physiology\src\biosim_physiology\drug_interactions\ode_coupled.py
"""Drugs modelled entirely by an ODE-coupled pathway (no static modifier):
SGLT2 inhibitors and exogenous insulin."""

from __future__ import annotations

from ._core import DrugDose


def _sglt2_inhibitor_modifiers(dose: DrugDose) -> dict[str, float]:
    """
    SGLT2 inhibitors (dapagliflozin, empagliflozin) intentionally return no
    PathwayParameters modifier here -- for the same reason as GLP-1 above,
    but even more strictly so: unlike metformin's hepatic coupling (which is
    genuinely COMPLEMENTARY to a separate static Gb modifier, since it is
    insulin-gated and near-zero at basal), SGLT2i's renal glucose clearance
    (biosim_ode.renal_pathway.RenalPathway) is continuously active whenever
    plasma glucose exceeds the drug-lowered threshold -- including at
    ordinary/basal glucose. A static Gb-lowering modifier alongside the ODE
    pathway would double-count exactly the way CortisolPathway/
    glucocorticoid_modifiers() once did (see CLAUDE.md KNOWN BUGS). The ODE
    pathway is the sole, complete representation; SimulationRunner detects
    these drugs and composes it instead of applying a modifier here.
    """
    return {}


def _exogenous_insulin_modifiers(dose: DrugDose) -> dict[str, float]:
    """
    Exogenous insulin (regular or the rapid-acting lispro analog --
    shared by both Drug.INSULIN_REGULAR and Drug.INSULIN_LISPRO in
    _DRUG_FUNCTIONS in registry.py) intentionally returns no PathwayParameters
    modifier here -- its effect is modelled entirely via an ODE-coupled
    pathway (biosim_ode.exogenous_insulin_pathway.ExogenousInsulinPathway,
    which itself dispatches to the correct calibration per drug),
    exactly like the GLP-1 (incretins.py) and SGLT2i drugs. SimulationRunner detects
    the specific drug and composes that pathway instead of applying a
    modifier here -- doing both would double-model the same mechanism.
    """
    return {}
