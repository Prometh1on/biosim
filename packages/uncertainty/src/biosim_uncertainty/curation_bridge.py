# SAVE TO: packages\uncertainty\src\biosim_uncertainty\curation_bridge.py
"""
Sensitivity-to-curation priority bridge.

Closes the epistemic loop: high-sensitivity parameters are automatically
escalated in the CurationQueue so experts review them first.

Flow:
    1. OAT sensitivity analysis ranks parameters by output influence
    2. This bridge maps high-sensitivity parameters to CurationQueue items
    3. Curators review highest-impact parameters first
    4. Better parameter values reduce simulation uncertainty
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class SensitivityCurationMapping:
    sensitivity_param:  str
    entity_name:        str
    parameter_kind_str: str
    description:        str


PARAM_TO_BIOLOGY: list[SensitivityCurationMapping] = [
    SensitivityCurationMapping("p1",  "D-Glucose",       "rate",
        "Glucose effectiveness - insulin-independent clearance"),
    SensitivityCurationMapping("p2",  "insulin receptor","rate",
        "X-compartment clearance - insulin action delay"),
    SensitivityCurationMapping("p3",  "insulin receptor","rate",
        "Insulin action - directly scales glucose uptake"),
    SensitivityCurationMapping("ke",  "insulin",         "rate",
        "Insulin elimination rate - plasma insulin half-life"),
    SensitivityCurationMapping("Gb",  "D-Glucose",       "abundance",
        "Basal glucose - model starting point"),
    SensitivityCurationMapping("Ib",  "insulin",         "abundance",
        "Basal insulin - sets insulin steady state"),
    SensitivityCurationMapping("EC50_mTOR", "mTOR",      "EC50",
        "mTOR activation threshold - anabolic signalling"),
    SensitivityCurationMapping("IGF1b", "IGF-1",         "abundance",
        "Basal IGF-1 level - growth factor signalling"),
]

PARAM_BIOLOGY_MAP: dict[str, SensitivityCurationMapping] = {
    m.sensitivity_param: m for m in PARAM_TO_BIOLOGY
}


@dataclass
class PrioritisationAction:
    parameter_name:    str
    sensitivity_score: float
    metric_name:       str
    new_priority:      str
    reason:            str


@dataclass
class BridgeResult:
    actions:             list[PrioritisationAction] = field(default_factory=list)
    parameters_upgraded: list[str]                  = field(default_factory=list)
    n_items_enqueued:    int = 0

    def summary(self) -> dict:
        return {
            "n_actions":           len(self.actions),
            "parameters_upgraded": self.parameters_upgraded,
            "n_items_enqueued":    self.n_items_enqueued,
        }


def prioritise_from_sensitivity(
    sensitivity_result,
    queue,
    store,
    registry,
    high_sensitivity_threshold:     float = 0.3,
    critical_sensitivity_threshold: float = 0.7,
    top_n_per_metric:               int   = 3,
) -> BridgeResult:
    """
    Enqueue or upgrade curation items based on sensitivity analysis results.
    """
    from biosim_curation import CurationItemType, CurationPriority
    from biosim_parameter_db import ParameterKind

    result = BridgeResult()

    for metric_name in sensitivity_result.metric_names:
        top_params = sensitivity_result.top_parameters(metric_name, n=top_n_per_metric)

        for param_name, score in top_params:
            abs_score = abs(score)
            if abs_score < high_sensitivity_threshold:
                continue

            if abs_score >= critical_sensitivity_threshold:
                priority_str = "CRITICAL"
                priority     = CurationPriority.CRITICAL
            else:
                priority_str = "HIGH"
                priority     = CurationPriority.HIGH

            bio = PARAM_BIOLOGY_MAP.get(param_name)
            if bio is None:
                continue

            entity_ids = registry.search_by_name(bio.entity_name)
            if not entity_ids:
                continue
            entity_id_str = str(entity_ids[0])

            try:
                kind = ParameterKind(bio.parameter_kind_str)
            except ValueError:
                kind = ParameterKind.CUSTOM

            best_param = store.best(entity_id_str, kind)
            value  = best_param.value.value if best_param else 0.0
            unit   = best_param.value.unit  if best_param else "unknown"
            source = best_param.provenance.pubmed_id if best_param else "not in store"

            item_id = f"sensitivity:{param_name}:{metric_name}"
            reason  = (
                f"Sensitivity |S_{param_name}| = {abs_score:.3f} for '{metric_name}'. "
                f"{bio.description}."
            )

            queue.enqueue(
                item_id=item_id,
                item_type=CurationItemType.PARAMETER,
                reason=reason,
                subject_name=bio.entity_name,
                parameter_kind=bio.parameter_kind_str,
                value=value,
                unit=unit,
                source_description=f"PubMed:{source}" if source else "No source",
                priority=priority,
            )

            result.actions.append(PrioritisationAction(
                parameter_name=param_name,
                sensitivity_score=abs_score,
                metric_name=metric_name,
                new_priority=priority_str,
                reason=reason,
            ))
            if param_name not in result.parameters_upgraded:
                result.parameters_upgraded.append(param_name)
            result.n_items_enqueued += 1

    return result
