# SAVE TO: packages\infrastructure\src\biosim_infra\infra.py
"""
Infrastructure layer -- simulation cache, privacy audit, versioning, reproducibility.

Sprint 14: SimulationCache
    Hash-based memoisation of simulation results.
    Identical scenarios (same name, meals, drugs, personalization, duration)
    return cached results without re-running the ODE.
    Thread-safe. LRU eviction when max_size is reached.

Sprint 13: PrivacyAudit
    Verifies that simulation outputs and personalisation results do not
    contain raw genetic data (rsIDs, allele information).
    Implements the data minimisation principle.

Sprint 15: ReproducibilityBundle and ModelVersion
    A ReproducibilityBundle captures everything needed to reproduce a
    simulation result exactly: model version, parameters, scenario,
    random seed, software versions.
    Every published result should be accompanied by its bundle.
"""

from __future__ import annotations

import hashlib
import importlib
import json
import sys
import dataclasses
from dataclasses import dataclass, field
from datetime import datetime, timezone
from collections import OrderedDict
from threading import RLock
from typing import Optional


# ---------------------------------------------------------------------------
# Sprint 14: Simulation Cache
# ---------------------------------------------------------------------------

class SimulationCache:
    """
    Thread-safe LRU cache for simulation results.

    Cache key: SHA-256 of the scenario's JSON representation.
    Same scenario + same personalisation + same duration = cache hit.

    Usage:
        cache = SimulationCache(max_size=500)
        runner = SimulationRunner()

        output = cache.get(scenario)
        if output is None:
            output = runner.run(scenario)
            cache.put(scenario, output)
    """

    def __init__(self, max_size: int = 500) -> None:
        if max_size < 1:
            raise ValueError(f"max_size must be >= 1, got {max_size}")
        self._max_size = max_size
        self._cache: OrderedDict[str, object] = OrderedDict()
        self._hits = 0
        self._misses = 0
        self._lock = RLock()

    def _key(self, scenario) -> str:
        """Compute a deterministic cache key from a scenario."""
        content = json.dumps(scenario.to_dict(), sort_keys=True)
        return hashlib.sha256(content.encode()).hexdigest()

    def get(self, scenario) -> Optional[object]:
        """Return cached result for scenario, or None if not cached."""
        with self._lock:
            k = self._key(scenario)
            if k in self._cache:
                # Move to end (most recently used)
                self._cache.move_to_end(k)
                self._hits += 1
                return self._cache[k]
            self._misses += 1
            return None

    def put(self, scenario, result) -> None:
        """Cache a simulation result. Evicts oldest entry if at capacity."""
        with self._lock:
            k = self._key(scenario)
            if k in self._cache:
                self._cache.move_to_end(k)
            self._cache[k] = result
            if len(self._cache) > self._max_size:
                self._cache.popitem(last=False)  # evict oldest

    def invalidate(self, scenario) -> bool:
        """Remove a specific scenario from the cache. Returns True if was cached."""
        with self._lock:
            k = self._key(scenario)
            if k in self._cache:
                del self._cache[k]
                return True
            return False

    def clear(self) -> None:
        """Clear all cached results."""
        with self._lock:
            self._cache.clear()
            self._hits = 0
            self._misses = 0

    @property
    def size(self) -> int:
        return len(self._cache)

    @property
    def hit_rate(self) -> float:
        total = self._hits + self._misses
        return self._hits / total if total > 0 else 0.0

    def stats(self) -> dict:
        return {
            "size":      self.size,
            "max_size":  self._max_size,
            "hits":      self._hits,
            "misses":    self._misses,
            "hit_rate":  round(self.hit_rate, 3),
        }


class CachedRunner:
    """
    SimulationRunner wrapper that uses a cache.
    Drop-in replacement for SimulationRunner for production use.
    """

    def __init__(self, cache: Optional[SimulationCache] = None, base_params=None) -> None:
        self._cache = cache or SimulationCache()
        # Same fix as SimulationRunner.__init__ (see runner.py) -- a
        # caller-supplied base (e.g. PlatformSession's literature-loaded
        # params) needs to reach the actual simulation, not just get
        # silently dropped in favour of hardcoded class defaults.
        self._base_params = base_params

    def run(self, scenario) -> object:
        result = self._cache.get(scenario)
        if result is not None:
            return result
        from biosim_api import SimulationRunner
        result = SimulationRunner(base_params=self._base_params).run(scenario)
        self._cache.put(scenario, result)
        return result

    @property
    def cache(self) -> SimulationCache:
        return self._cache


# ---------------------------------------------------------------------------
# Sprint 13: Privacy Audit
# ---------------------------------------------------------------------------

@dataclass
class PrivacyAuditResult:
    """Result of a privacy audit on a simulation output."""
    passed:     bool
    issues:     list[str]
    checked_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))

    @property
    def n_issues(self) -> int:
        return len(self.issues)

    def summary(self) -> dict:
        return {
            "passed":   self.passed,
            "n_issues": self.n_issues,
            "issues":   self.issues,
        }


class PrivacyAuditor:
    """
    Audits simulation outputs for privacy compliance.

    Checks:
    1. No raw rsIDs in simulation output summary
    2. No allele information in exported data
    3. Genetic personalisation result exposes only parameter deltas
    4. No raw variant objects in output attributes
    """

    def audit(self, output) -> PrivacyAuditResult:
        """
        Run a full privacy audit on a SimulationOutput.

        Returns:
            PrivacyAuditResult with pass/fail and any issues found
        """
        issues = []

        # Check 1: No rsIDs in simulation summary
        summary_str = json.dumps(output.summary())
        if self._contains_rsid(summary_str):
            issues.append(
                "Raw rsID found in simulation summary dict. "
                "Genetic variant identifiers must not appear in exported results."
            )

        # Check 2: Genetic personalisation result is safe
        if output.genetic_personalization is not None:
            gp_summary = json.dumps(output.genetic_personalization.summary())
            if self._contains_rsid(gp_summary):
                issues.append(
                    "Raw rsID found in genetic personalisation summary. "
                    "PersonalizationResult.summary() must contain only parameter modifiers."
                )

        # Check 3: No raw rsIDs in the ComposedEffects PersonalizationResult
        # actually carries.
        #
        # Regression: this used to check hasattr(gp, 'variants') --
        # PersonalizationResult has never had a `variants` field (see
        # personalizer.py), so this check was structurally always False,
        # a vacuous no-op. The real raw-identifier-shaped field is
        # ComposedEffect.contributing_variants
        # (PersonalizationResult.applied_effects), which this check never
        # looked at. Confirmed directly: constructing a real
        # GenotypePersonalizer().personalize(...) result and running the
        # OLD check always passed regardless of content. contributing_variants
        # is now fixed at the source to hold privacy-safe fingerprints,
        # not raw rsIDs (see composer.py's _rsid_fingerprint()) -- this
        # check is the defense-in-depth backstop for that invariant, same
        # pattern as Check 1/2's _contains_rsid() scans.
        if hasattr(output, 'genetic_personalization') and output.genetic_personalization is not None:
            gp = output.genetic_personalization
            effects_str = json.dumps(
                [list(e.contributing_variants) for e in gp.applied_effects]
            )
            if self._contains_rsid(effects_str):
                issues.append(
                    "Raw rsID found in PersonalizationResult.applied_effects "
                    "(ComposedEffect.contributing_variants). This field must "
                    "contain only privacy-safe fingerprints, never raw "
                    "variant identifiers."
                )

        # Check 4: Lay summary contains no biological identifiers
        if hasattr(output, 'lay_summary'):
            ls = output.lay_summary
            if hasattr(ls, 'glucose_narrative'):
                if self._contains_rsid(ls.glucose_narrative):
                    issues.append("rsID found in lay summary narrative.")

        # Check 5: No raw rsIDs in output.warnings.
        #
        # Regression: a mandatory security-review of Check 3 above found
        # this exact sibling leak -- ComposedEffect.population_warnings
        # (built in composer.py's population-mismatch branch) used to
        # embed the raw rsid and gene name directly in the warning text.
        # That list flows PersonalizationResult.population_warnings ->
        # (runner.py) SimulationOutput.warnings -> the REST API's
        # response verbatim (schemas.py's warnings: list[str]) --
        # Check 1's summary()-based scan never sees it because
        # SimulationOutput.summary() only exposes a warning *count*, not
        # the text. Fixed at the source in composer.py (no longer embeds
        # rsid/gene at all); this check is the defense-in-depth backstop.
        if hasattr(output, 'warnings') and output.warnings:
            warnings_str = json.dumps(list(output.warnings))
            if self._contains_rsid(warnings_str):
                issues.append(
                    "Raw rsID found in SimulationOutput.warnings. "
                    "Population-mismatch and other genetic warnings must "
                    "not embed raw variant identifiers."
                )

        return PrivacyAuditResult(passed=len(issues) == 0, issues=issues)

    def audit_personalization_result(self, result) -> PrivacyAuditResult:
        """Audit a PersonalizationResult directly."""
        issues = []
        summary_str = json.dumps(result.summary())
        if self._contains_rsid(summary_str):
            issues.append("Raw rsID in PersonalizationResult.summary()")
        # Regression: this used to check hasattr(result, 'original_variants')
        # -- PersonalizationResult has never had that field, so this was a
        # vacuous no-op (confirmed: always False, on every real result).
        # The real raw-identifier-shaped field is applied_effects[*]
        # .contributing_variants -- now fixed at the source to hold
        # privacy-safe fingerprints (composer.py's _rsid_fingerprint()),
        # this check is the defense-in-depth backstop for that invariant.
        effects_str = json.dumps(
            [list(e.contributing_variants) for e in result.applied_effects]
        )
        if self._contains_rsid(effects_str):
            issues.append(
                "Raw rsID found in applied_effects[*].contributing_variants — "
                "privacy violation."
            )
        # Same sibling leak as Check 5 in audit() -- population_warnings
        # (both the aggregate and per-effect) used to embed the raw rsid
        # and gene name; fixed at the source in composer.py, this is the
        # defense-in-depth backstop.
        warnings_str = json.dumps(
            list(result.population_warnings)
            + [w for e in result.applied_effects for w in e.population_warnings]
        )
        if self._contains_rsid(warnings_str):
            issues.append(
                "Raw rsID found in PersonalizationResult.population_warnings — "
                "privacy violation."
            )
        return PrivacyAuditResult(passed=len(issues) == 0, issues=issues)

    @staticmethod
    def _contains_rsid(text: str) -> bool:
        """Check if text contains a pattern that looks like an rsID."""
        import re
        return bool(re.search(r'\brs\d{4,}\b', text))


# ---------------------------------------------------------------------------
# Sprint 15: Versioning and Reproducibility
# ---------------------------------------------------------------------------

PLATFORM_VERSION = "0.1.0"

# Regression: this used to be its own hardcoded literal ("glucose_insulin_v1.0"),
# an independent copy of ODESystem.MODEL_VERSION rather than a reference to
# it -- the two could silently drift apart if the ODE model is ever
# versioned forward without remembering this second copy, which would
# make every future ReproducibilityBundle report the WRONG model version.
# Import the authoritative value instead of duplicating the literal.
from biosim_ode.system import ODESystem as _ODESystem
MODEL_VERSION = _ODESystem.MODEL_VERSION


@dataclass
class ModelVersion:
    """Version information for a specific model release."""
    platform_version: str
    model_id:         str
    model_version:    str
    release_date:     str
    changelog:        str = ""

    def cite(self) -> str:
        return (
            f"Biosim Platform v{self.platform_version}, "
            f"model {self.model_id} v{self.model_version}, "
            f"released {self.release_date}."
        )


CURRENT_VERSION = ModelVersion(
    platform_version=PLATFORM_VERSION,
    model_id="glucose_insulin",
    model_version="1.0.0",
    release_date="2026-01-01",
    changelog="Initial release. 6-species glucose/insulin/IGF-1/mTOR/AMPK model.",
)


@dataclass
class ReproducibilityBundle:
    """
    Everything needed to reproduce a simulation result exactly.

    Store this alongside published results. Anyone with this bundle
    and the biosim platform can reproduce the exact simulation.
    """
    platform_version:  str
    model_version:     str
    scenario_dict:     dict          # ScenarioDefinition.to_dict()
    parameter_values:  dict          # PathwayParameters as dict
    random_seed:       int
    python_version:    str
    numpy_version:     str
    scipy_version:     str
    created_at:        str

    def to_json(self) -> str:
        return json.dumps(dataclasses.asdict(self), indent=2, default=str)

    @classmethod
    def from_json(cls, json_str: str) -> "ReproducibilityBundle":
        return cls(**json.loads(json_str))

    def citation(self) -> str:
        return (
            f"Biosim Platform v{self.platform_version}, "
            f"model {self.model_version}, "
            f"created {self.created_at[:10]}. "
            f"Scenario: {self.scenario_dict.get('name', 'unknown')}."
        )


def create_bundle(output, random_seed: int = 42) -> ReproducibilityBundle:
    """
    Create a ReproducibilityBundle from a SimulationOutput.

    Args:
        output:      SimulationOutput to document
        random_seed: The random seed used for any stochastic components

    Returns:
        ReproducibilityBundle capturing full provenance

    Regression: this used to hand-build scenario_dict from just
    output.scenario_name/output.scenario_description -- SimulationOutput
    never carried the actual ScenarioDefinition, only those two strings,
    so meals/drugs/personalization were silently absent from every bundle
    despite the field comment ("ScenarioDefinition.to_dict()") and this
    class's own docstring ("Everything needed to reproduce a simulation
    result exactly") both implying otherwise. Fixed by adding
    SimulationOutput.scenario (set by SimulationRunner.run(), see
    runner.py) and using its real to_dict() here.

    One real limitation this does NOT close, documented rather than
    hidden: to_dict()'s genetic_variant_fingerprints (see scenario.py)
    are HMAC-keyed with a secret regenerated every process start,
    specifically so they can't be reversed back to a raw rsid (rule 7).
    That means a bundle for a genetically-personalized scenario is
    privacy-safe to store/publish, but its genetic dimension specifically
    can't be reconstructed byte-for-byte by a different process/session
    from the fingerprint alone -- "reproduce exactly" holds for
    meals/drugs/physiological profile, not for the exact genetic variant
    input. No fix attempted here (would need a deliberately different,
    separately-designed scheme -- see scenario.py's own docstring on this
    exact point).
    """
    import numpy as np
    import scipy

    # Extract parameter values from the simulation
    param_values = {}
    if hasattr(output, 'physiological_personalization') and output.physiological_personalization:
        params = output.physiological_personalization.personalized_params
        param_values = dataclasses.asdict(params)
    else:
        from biosim_ode import PathwayParameters
        param_values = dataclasses.asdict(PathwayParameters())

    scenario_dict = (
        output.scenario.to_dict() if getattr(output, "scenario", None) is not None
        else {"name": output.scenario_name, "description": output.scenario_description}
    )

    return ReproducibilityBundle(
        platform_version=PLATFORM_VERSION,
        model_version=MODEL_VERSION,
        scenario_dict=scenario_dict,
        parameter_values=param_values,
        random_seed=random_seed,
        python_version=sys.version.split()[0],
        numpy_version=np.__version__,
        scipy_version=scipy.__version__,
        created_at=datetime.now(timezone.utc).isoformat(),
    )
