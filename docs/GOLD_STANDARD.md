# Gold standard: how a benchmark result is graded

Frozen 2026-11-20 (S0), before any model change made under it. Code:
`packages/validation/src/biosim_validation/grading.py`. Tests:
`packages/validation/tests/test_grading.py`.

## Why a second grade

The validator's PASS / MISS verdict asks "is the model within a fixed band (±10%, ±25%,
±50% or ±20%) of the published number?" That ignores how certain the published number is.
The **grade** asks a different question: *is the model inside the uncertainty the paper
itself reports?* It is computed alongside the verdict and never changes it. It is not always
stricter: for small, imprecise studies the reported uncertainty is wider than the old band
(see "Read these with care" below).

## The grades

| Grade | Meaning |
|---|---|
| **SET** | The model starts at this value (fasting glucose), so it is not a prediction. |
| **NOT-ASSESSABLE** | No cited number to compare with: the reference is UNVERIFIED, or it was read off a figure (the Bergman IVGTT). |
| **CIRCULAR** | The dataset was used to fit or choose a model parameter, so agreement proves little. |
| **GOLD / SILVER / MISS** | Set by the kind of uncertainty the source reports (next table). |

| What the source reports | GOLD | SILVER | MISS |
|---|---|---|---|
| A treatment difference with a 95% CI | model inside the CI | — | outside the CI |
| A mean with a between-subject SD (or SEM, converted with SD = SEM × √n) | within 0.5 SD | within 1 SD | further |
| No spread at all | never | inside the old tolerance band | outside it |

**CIRCULAR register** (datasets used to fit or choose a parameter):
- `OGTT_75g_healthy`: the absorption constants (also the default meal) and the first-phase
  insulin gain were fitted to it.
- `T2D_vs_nondiabetic_meal_insulin`: the T2D glucose-effectiveness factor came from the same
  Basu 2009 cohort.

## Development set and locked test set

- **Development set** (`ModelValidator.validate_all()`): the 9 existing benchmarks. Later work
  may look at these when deciding what to change, but every change still needs a source.
- **Locked test set** (`ModelValidator.validate_heldout()`): the three B2 held-out datasets
  (Navarro-Masip 2025 OGTT and liquid meal; Eshghi 2017 4-h OGTT). Nothing may be chosen or
  adjusted using them. They are reported once at the end of each stage.
  **Honest caveat:** these were first run in B2 and looked at before B7, so they are
  "already seen", not untouched. No new, never-seen held-out set (H2) was added in S0.

## How the freeze is enforced

`benchmark_manifest_sha256()` hashes every benchmark value (both sets, including the spread
fields and the set labels) and the grading code. A test fails unless it equals
`frozen_manifest.FROZEN_SHA256`, and another test requires that hash in `HISTORY.md`
("Benchmark freeze: …" or "Benchmark change: …"). This is tamper-evident, not independent:
a change needs the hash and the HISTORY line updated together, which shows up in the commit. Not enforced by a test: whether a model change was *chosen* by looking at the locked
set. That rests on the written rule above and on reporting the locked set only at stage ends.

## Where the spreads come from (all fetched 2026-11-20)

| Benchmark | Spread | Source text |
|---|---|---|
| Semaglutide | 95% CI of the reduction 0.27 to 2.42 mmol·h/L, n=30 | abstract: "glucose -1.34 mmol h/L [-2.42, -0.27]" |
| Metformin | 95% CI 64 to 307 mmol/L×min, n=15 | abstract: "MD 186 mmol/L × min [95% CI 64-307]" |
| Navarro-Masip 2025 (locked) | SD per time point, n=24 | Table 2, "mean ± SD" |
| Eshghi 2017 (locked) | SEM per group, n=6 | Table 5, "mean ± SEM" |
| Prednisone (Shah 2020) | none usable (medians with IQR) | Table 1 |
| Canagliflozin, pea protein, caffeine, Basu ratios, OGTT | none reported for the benchmarked number | — |

## Results at the freeze (model at commit 4c9554c, before any S0+ change)

**Development set (24 points):** GOLD 2 (semaglutide +0.7%; metformin −31%, inside its CI),
SILVER 0, MISS 9 (prednisone, canagliflozin ×2, pea protein ×4, caffeine ×2), CIRCULAR 9,
NOT-ASSESSABLE 3, SET 1.

**Locked set (12 points):** GOLD 10, SILVER 2 (healthy 4-h glucose iAUC +95%; T2D insulin
iAUC −56%).

**Read these with care:**
- **The CI rule is generous for both drug GOLDs.** A 95% CI of a mean difference measures how
  precisely a small trial pinned the effect, and here it is wide: semaglutide's CI is −80% to
  +81% around its mean (n=30), metformin's −66% to +65% (n=15), both wider than the old ±50%
  band. "GOLD" here means "inside the trial's own uncertainty", not "close to the number".
  Do not quote "2 GOLD drug benchmarks" without this.
- **The SD rule is generous for absolute glucose.** Between-subject SDs are 2–3 mmol/L, so a
  14% miss can still be GOLD. And with n=6 the healthy glucose-area SD is so wide that a +95%
  miss is only SILVER. The rule is frozen, so this is reported rather than changed after
  seeing the results.
- **GOLD on a drug benchmark is not "validated".** The drug label rule (VALIDATED needs every
  point within ±20%, independent, *and reviewed*) is separate and stricter. Metformin is GOLD
  here at −31% and stays EXPERIMENTAL. Semaglutide meets both rules numerically; after review
  (2026-11-21) it is labelled VALIDATED, with the caveat that modelling the trial's obese
  population gives −26% to −39% (docs/VALIDATION.md).
