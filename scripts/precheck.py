# SAVE TO: scripts\precheck.py
"""
Deterministic (non-agent) pre-check script for Biosim (MASTER_PLAN.md
Section 10, item C4).

Covers the ~3 bug classes this project has caught repeatedly by hand
across many milestones -- freeing the mandatory LLM code-reviewer/
security-reviewer passes for genuine judgment calls rather than
re-discovering the same mechanical pattern each time:

1. A dose-accepting function that takes a raw `dose_mg`/`dose_units:
   float` parameter (not an already-validated DrugDose/DrugEvent) with
   no `math.isfinite(...)` guard in its own body -- the exact class of
   bug found and fixed for `acarbose_kabs_scale()`/`miglitol_kabs_scale()`/
   `tirzepatide_gb_factor()` by prior security-reviewer passes (see
   CLAUDE.md's own history of each). A DrugDose/DrugEvent-typed
   parameter is exempt -- both dataclasses already validate finiteness
   in their own __post_init__ (scenario.py/drug_interactions.py).
2. A `ParameterKind.CUSTOM` registration in the `parameter_seed` package whose
   (subject, custom_name) pair collides with an earlier registration --
   `ParameterStore._index_key()` deliberately excludes `custom_kind_name`
   from a bare subject+kind comparison for `query()`/`best()`, but a
   genuinely duplicate (subject, custom_name) pair would still register
   as a real, spurious parameter conflict at seed time (the exact bug
   class fixed in MASTER_PLAN row 5.1's chunk 1) -- catching it
   statically is cheaper than running the seed and calling
   `store.unresolved_conflicts()`.
3. A function that composes a `"Gb"` (or `"G_threshold"`) modifier
   entry with no corresponding mention in MASTER_PLAN.md's Section 11
   coverage matrix -- a reminder (WARN, not a hard failure) that a new
   Gb-touching mechanism should be checked against
   MIN_GB_THRESHOLD_MARGIN_MM the way row 1.23/Section 11 already did
   for tirzepatide and metformin, not a guarantee that it will be.

This is a heuristic, not a full static-analysis tool -- false negatives
are possible (see each check's own docstring for its exact scope) and
it is not a substitute for the mandatory code-reviewer/security-reviewer
passes on anything touching dosing/safety-clamp/disease-state logic
(see CLAUDE.md's Q5 staggering rule) -- it exists to make THOSE passes
cheaper by pre-catching the mechanical, repeatedly-seen cases.

Usage:
    python scripts/precheck.py

Exit code 0 = no FAIL-level findings (WARN-level findings may still be
printed); exit code 1 = at least one FAIL-level finding.
"""
from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

DOSE_PARAM_NAMES = {"dose_mg", "dose_units"}
DOSE_TYPED_OBJECT_EXEMPT_ANNOTATIONS = {"DrugDose", "DrugEvent", "InsulinDoseEvent"}


def _display_path(path: Path) -> str:
    """Relative to REPO_ROOT when possible, else the path as-is (so
    this reports sensibly against a tmp_path fixture in tests, which
    lives outside REPO_ROOT)."""
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def _source_files(root: Path) -> list[Path]:
    """Every .py file under `root`, recursively, in a stable order."""
    return sorted(p for p in root.rglob("*.py") if "__pycache__" not in p.parts)


def _annotation_name(annotation: ast.AST | None) -> str:
    if annotation is None:
        return ""
    if isinstance(annotation, ast.Name):
        return annotation.id
    if isinstance(annotation, ast.Attribute):
        return annotation.attr
    if isinstance(annotation, ast.Constant) and isinstance(annotation.value, str):
        return annotation.value
    return ast.dump(annotation)


def _calls_isfinite_or_isnan(node: ast.AST) -> bool:
    for child in ast.walk(node):
        if isinstance(child, ast.Call):
            func = child.func
            name = func.attr if isinstance(func, ast.Attribute) else (
                func.id if isinstance(func, ast.Name) else ""
            )
            if name in ("isfinite", "isnan"):
                return True
    return False


def check_dose_functions_have_finiteness_guards(path: Path) -> list[str]:
    """
    Check 1 (see module docstring). Scope: only flags a function whose
    dose parameter is annotated `float` (or unannotated) -- a
    DrugDose/DrugEvent/InsulinDoseEvent-typed parameter is exempt, since
    those dataclasses already validate finiteness in their own
    __post_init__ before any function body ever runs.
    """
    findings = []
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in ast.walk(tree):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        for arg in node.args.args:
            if arg.arg not in DOSE_PARAM_NAMES:
                continue
            ann = _annotation_name(arg.annotation)
            if ann in DOSE_TYPED_OBJECT_EXEMPT_ANNOTATIONS:
                continue  # already validated upstream by the dataclass itself
            if not _calls_isfinite_or_isnan(node):
                findings.append(
                    f"{_display_path(path)}:{node.lineno}: "
                    f"{node.name}(..., {arg.arg}: {ann or '<untyped>'}, ...) has no "
                    "math.isfinite()/math.isnan() guard in its own body"
                )
    return findings


_CUSTOM_KIND_BLOCK_RE = re.compile(
    r"\.for_subject\((?P<subject>\w+)\)"
    r"(?:(?!\.for_subject\().)*?"
    r"custom_name=(?P<quote>[\"'])(?P<name>[^\"']+)(?P=quote)",
    re.DOTALL,
)

# A ParameterBuilder chain whose nearest preceding `store.` call is
# `store.supersede(<old_id>, ...)` (not `store.add(...)`) is a DELIBERATE
# re-registration under the same (subject, custom_name) -- the
# append-only store's own documented update mechanism (rule 5, CLAUDE.md
# ARCHITECTURE), not a duplicate-registration bug. Found as a genuine
# false positive on this exact check's first real run, against
# parameter_seed.py's own p3 supersede pair (Bergman 1981 Table 1 ->
# Fig 6 caption) -- confirmed via store.unresolved_conflicts() returning
# empty for that real pair before adding this exemption.
_STORE_CALL_RE = re.compile(r"store\.(add|supersede)\(")


def _nearest_preceding_store_call(text: str, position: int) -> str | None:
    kind = None
    for m in _STORE_CALL_RE.finditer(text, 0, position):
        kind = m.group(1)
    return kind


def check_duplicate_custom_parameter_registrations(paths: Path | list[Path]) -> list[str]:
    """
    Check 2 (see module docstring). Regex-based, not AST-based -- the
    builder pattern here is a multi-line fluent chain
    (`ParameterBuilder().for_subject(x).kind(ParameterKind.CUSTOM,
    custom_name="y")...`), and a plain textual scan for the
    (for_subject, custom_name) pair within one chain is simpler and
    robust enough for this file's own consistent style than reproducing
    the builder's call-chain semantics in AST.

    Accepts several files (parameter_seed is a package since 2026-10-20)
    and shares one `seen` table across them, so a duplicate split across
    two modules is still caught.
    """
    seen: dict[tuple[str, str], str] = {}
    findings = []
    for path in [paths] if isinstance(paths, Path) else paths:
        text = path.read_text(encoding="utf-8")
        for match in _CUSTOM_KIND_BLOCK_RE.finditer(text):
            key = (match.group("subject"), match.group("name"))
            where = f"{_display_path(path)}:{text.count(chr(10), 0, match.start()) + 1}"
            if _nearest_preceding_store_call(text, match.start()) == "supersede":
                continue  # a deliberate supersede, not a duplicate registration
            if key in seen:
                findings.append(
                    f"{where}: duplicate CUSTOM "
                    f"registration for subject={key[0]!r}, custom_name={key[1]!r} "
                    f"(first seen at {seen[key]})"
                )
            else:
                seen[key] = where
    return findings


_GB_MODIFIER_RE = re.compile(r'["\'](?:Gb|G_threshold)["\']\s*:')


def check_gb_modifiers_are_in_the_coverage_matrix(paths: list[Path]) -> list[str]:
    """
    Check 3 (see module docstring). WARN-level only -- a reminder, not a
    hard gate, since MASTER_PLAN.md's own prose format can reference a
    function name in ways this simple substring check won't recognize.
    """
    master_plan = (REPO_ROOT / "MASTER_PLAN.md").read_text(encoding="utf-8")
    warnings = []
    for path in paths:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            src = ast.get_source_segment(path.read_text(encoding="utf-8"), node) or ""
            if _GB_MODIFIER_RE.search(src) and node.name not in master_plan:
                warnings.append(
                    f"{_display_path(path)}:{node.lineno}: {node.name}() "
                    "composes a Gb/G_threshold modifier but its name doesn't appear "
                    "anywhere in MASTER_PLAN.md -- consider adding it to Section 11's "
                    "coverage matrix (is it swept against MIN_GB_THRESHOLD_MARGIN_MM?)"
                )
    return warnings


def main() -> int:
    physiology_src = REPO_ROOT / "packages" / "physiology" / "src" / "biosim_physiology"
    data_seeds_src = REPO_ROOT / "packages" / "data-seeds" / "src" / "biosim_seeds"

    fails: list[str] = []
    warns: list[str] = []

    # rglob, not glob: drug_interactions/ and modifiers/ are subpackages
    # (2026-10-20 split) and a flat glob would silently skip them.
    for py_file in _source_files(physiology_src):
        fails += check_dose_functions_have_finiteness_guards(py_file)

    seed_files = _source_files(data_seeds_src / "parameter_seed")
    if not seed_files:
        fails.append(f"{_display_path(data_seeds_src)}: no parameter_seed sources "
                     "found -- check 2 would silently not run; update precheck.py")
    fails += check_duplicate_custom_parameter_registrations(seed_files)

    warns += check_gb_modifiers_are_in_the_coverage_matrix(
        _source_files(physiology_src / "drug_interactions")
        + _source_files(physiology_src / "modifiers")
    )

    if warns:
        print("WARN (review, not blocking):")
        for w in warns:
            print(f"  - {w}")
    if fails:
        print("FAIL:")
        for f in fails:
            print(f"  - {f}")
        return 1

    print(f"precheck.py: clean ({len(warns)} warning(s), 0 failure(s))")
    return 0


if __name__ == "__main__":
    sys.exit(main())
