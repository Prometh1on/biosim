import sys
from pathlib import Path

ROOT = Path(__file__).parent

PACKAGES = [
    "id-resolver",
    "context-registry",
    "parameter-db",
    "literature-pipeline",
    "curation",
    "ode-engine",
    "uncertainty",
    "genotype",
    "physiology",
    "simulation-api",
    "expert-tools",
    "lay-output",
    "validation",
    "infrastructure",
    "data-seeds",
    "rest-api",
    "metabolism",
]

for package in PACKAGES:
    src = ROOT / "packages" / package / "src"
    if src.exists() and str(src) not in sys.path:
        sys.path.insert(0, str(src))
