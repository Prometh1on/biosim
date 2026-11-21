# Parameter provenance status

Audit of 2026-11-15. No parameter value was changed. This file and the
per-parameter `verification_status` field record what the cited sources say.

**Method.** For each of the 53 cited papers, the PubMed abstract (and the Europe
PMC full text, open for 13 of them) was fetched live. Four `literature-curator`
agents classified the 77 cited store entries using only that text. Each quote was
then checked by machine to be a verbatim substring of the fetched text (all 77
passed), and every DERIVED calculation was redone by hand. The table of
results is `packages/data-seeds/src/biosim_seeds/parameter_seed/verification.py`.
Every seeded `Parameter` carries its status in `provenance.verification_status`
and its quote in `provenance.verification_quote`.
`TestProvenanceVerificationStatus` fails if any seeded value has no status.

**Statuses**
- **VERIFIED**: the fetched text states the value.
- **DERIVED**: the value follows from numbers the text states (arithmetic in the
  table). The number is checked. Whether that measurement maps onto the model
  parameter (e.g. a clamp glucose-infusion change applied as a `p3` multiplier) is
  a modelling assumption that the table's `reason` field notes.
- **PAPER-ON-TOPIC-VALUE-NOT-FOUND**: right paper, but the value isn't in the
  available text. Usually the paper is abstract-only, or the value is an
  engineering back-solve through this model.
- **WRONG-PAPER**: the paper is about a different quantity.
- **CONTRADICTED**: the paper states a different value for the same quantity.
- **UNVERIFIED**: no source is cited at all.

## Counts: seed `ParameterStore`, 80 active entries (updated 2026-11-20)

| Status | Count |
|---|---|
| VERIFIED | 10 |
| DERIVED | 34 |
| PAPER-ON-TOPIC-VALUE-NOT-FOUND | 32 |
| WRONG-PAPER | 0 (was 1, removed) |
| CONTRADICTED | 0 |
| UNVERIFIED (each labelled MODEL ASSUMPTION with a reason) | 4 (was 5: caffeine's uncited +2% fasting-glucose term was removed from the model on 2026-11-20) |

So 76 of 80 cite a source, and **44 values (10 + 34) are confirmed against the
source text**.

## B4 clean-up (2026-11-16)

- **WRONG-PAPER removed.** The `p1` replicate (0.028 min⁻¹) credited to DeFronzo,
  Tobin & Andres 1979 (PMID:382871), the glucose-clamp *methods* paper, which never
  estimates glucose effectiveness, is no longer seeded. The runtime `p1` (0.0317)
  never used it. The same paper was also cited for the in-house OGTT benchmark's
  2-hour glucose (6.5 mM); it reports no OGTT value, so that point is now
  UNVERIFIED and unscored (`benchmark.py`).
- **The 32 "value not found" entries were re-checked for open full text.** They come
  from 22 papers; Europe PMC has open full text for only one (Chen 2016, acarbose),
  which the first audit already read. Its record now quotes the Table 2 glucose values
  the 0.35 back-solve targets; the status stays "not found", because 0.35 itself is a
  model back-solve, not a number in the paper. Bergman 1981 and Sherwin 1974 are free
  on the publisher's site but only as image scans; reading them needs OCR or a person,
  and was not done. No status was upgraded.
- **The 5 UNVERIFIED entries are labelled MODEL ASSUMPTION with a reason**
  (untreated T1D fasting glucose ×2.2; the pre-diabetes preset's insulin action ×0.7,
  fasting glucose ×1.03 and secretion ×1.05; caffeine's 2% fasting-glucose rise; that last one was removed from the model on 2026-11-20).
- **The default-meal absorption constants (B3)** are labelled MODEL ASSUMPTION in
  `DEFAULT_MEAL_ABSORPTION_PROVENANCE` (`gut_absorption.py`): they reuse the
  OGTT-fitted constants; no mixed-meal source was found.

## CONTRADICTED

None found.

## VERIFIED (10)

Dapagliflozin TmG healthy 317 / T2D 420 mg/min (PMID:23735727); empagliflozin
TmG 337 / 459 mg/min (PMID:28428225); caffeine p3 reduction 0.15 (PMID:11815511;
an IV-caffeine study, so the 200 mg oral mapping is the project's choice);
fasting GLP-1 7 pM (PMID:8138058); NEFA half-life 3.7 min (PMID:1951681); NEFA
insulin IC50 68 pM (PMID:11061533); hyperthyroid basal EGP ×1.2 (PMID:3510000);
alcohol EGP reduction 0.12 at 48 g (PMID:9815011).

## DERIVED (34)

These include the IGF-1 half-life (810 min, midpoint of "12-15 h", PMID:2558477),
`KG_secretion` (10.13 − 4.5 mM, PMID:30166558), GFR 127.5, cortisol production
10 (rounded from 9.9 mg/day), and 30 disease/drug/lifestyle factors computed as
ratios or (1 − fraction) of quoted numbers. Each record's arithmetic is in
`verification.py`.

## PAPER-ON-TOPIC-VALUE-NOT-FOUND (32)

- **Core minimal model**: p1 0.0317, Gb 5.0, Ib 66 and p3 9.032e-7 (Bergman 1981,
  PMID:7033284, abstract-only; an earlier session quoted the p3 Fig 6 caption, but
  that PDF could not be re-fetched); `ke` 0.093 (Sherwin 1974, PMID:4856884; the
  abstract says the curves "are multiexponential"); `p2` 0.0123 (MINMOD software
  paper, PMID:3640682).
- **Signalling read-outs**: insulin receptor Kd 0.1 nM (PMID:12360255); mTOR EC50
  100 pM and Hill 1.5 (PMID:15718470, an mTORC2/Akt paper); AMPK IC50 200 pM
  (PMID:22436748); IGF-1 basal 13 nM (PMID:9032050, empty abstract).
- **Engineering back-solves through this model**: exogenous insulin ka and pM/U
  (×4, PMID:8314011), acarbose 0.35 and miglitol 0.1446 kabs scales, protein
  secretion 14.93 pM/g, GLP-1 meal response 4.0 pM/mM, sulfonylurea 0.5 and
  glipizide 0.9 fractions of ke·Ib, and glucocorticoid 0.30 / 0.10 / 0.10 at 30 mg.
- **Direction-only sources**: T2D factors 0.6 / 0.55 / 1.3 (Stumvoll 2005 review),
  hypothyroid 0.85 / 0.88, T1D floors 0.05 / 0.02, prednisone:cortisol 4:1
  (PMID:888843, empty abstract).

## UNVERIFIED (4)

T1D `Gb_factor` 2.2; insulin-resistant `p3_factor` 0.7, `Gb_factor` 1.03 and
`secretion_scale_factor` 1.05. (Caffeine's `gb_increase_amplitude` 0.02 was the fifth until
2026-11-20, when the term was removed from the model.)

## Model assumptions: `PathwayParameters` defaults with no store entry

These are labelled `MODEL ASSUMPTION` in `biosim_ode/pathway.py`. No source range
was found for any of them, so none is given.

| Parameter | Value | What it controls |
|---|---|---|
| `extra_g_clearance` | 0 (was 0.018 min⁻¹ until B7, 2026-11-16) | removed: it double-counted p1; exercise uptake is now its own MODEL ASSUMPTION, 0.072 min⁻¹ per unit intensity |
| `G_threshold` | 4.5 mM | glucose level above which beta cells secrete (Rebelato's rat islets started at ">5.6 mM", which this value does not use) |
| `first_phase_gain` | 121.3 pM per mM/min | size of first-phase secretion; form from Mari 2002, magnitude fitted to the OGTT benchmark |
| `IGF1_ins_sensitivity` | 0.1 | how much insulin raises IGF-1 production |
| `k_mTOR` | 0.05 min⁻¹ | how fast mTOR activity follows insulin |
| `k_AMPK` | 0.02 min⁻¹ | how fast AMPK activity follows insulin |
| `n_AMPK` | 1.0 | steepness of insulin's AMPK inhibition |
| `AMPK_max` | 1.0 | 0–1 normalisation ceiling of AMPK activity |

## Benchmarks and variant effects (not in the store)

- `mtor_activation` benchmark (0.55): UNVERIFIED. Human muscle studies report
  fold-changes, not an absolute 0–1 activity. Gran 2014 (PMID:25302072) gives
  p70S6K only in figures, after a 31 g-protein meal.
- rs35767 (IGF1) +0.5 nM/allele: UNVERIFIED, and its direction is contested
  (PMID:30115944 vs 22689104, 21046444). Runs that include it carry a warning.
- Lehmann & Deutsch 1992 (PMID:1588781): abstract-only check. Its "maximal rate of
  gastric emptying" supports the removed attribution.
