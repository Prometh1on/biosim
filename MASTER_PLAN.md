# BIOSIM MASTER PLAN

**Status: living planning document, not a memory log.** `CLAUDE.md` remains
the append-only record of what has actually been built, fixed, and
verified — read it first, always, for ground truth on current state.
This file is the forward-looking decomposition: the full shape of what
Biosim is meant to become, broken into chunks sized for one agentic
session each (roughly what a single `Agent`/orchestrating-session pass
can scope, implement, test, review, document, and commit in one sitting
— matching the actual milestones already completed, e.g. "fix
`ParameterStore.get_all_for_subject()`" or "gate first-phase insulin
secretion on Brunzell's threshold").

Update this file when a chunk's status changes (planned -> in progress
-> done -> see CLAUDE.md entry). Do not duplicate `CLAUDE.md`'s content
here — link to the relevant DONE entry instead of re-describing it.

Last synced against `CLAUDE.md` at commit `fca9e15`.

**Condensed 2026-10-20:** finished rows and Section 8 items were cut to their
first sentence (marked `…`) so this file can be read in one pass. The full
original text of every row is at `git show 177a782:MASTER_PLAN.md`; the
per-chunk details are in `HISTORY.md`.

---

## 0. Vision

Biosim simulates postprandial glucose/insulin dynamics for a specific
individual — their food, genetics, medications, disease state, exercise,
and time of day — and reports both a quantitative (expert) and
plain-language (lay) prediction, with honest uncertainty bounds and full
literature provenance on every parameter. Every number is either
measured, cited, or explicitly labeled as an isolated engineering
assumption — never fabricated, never curve-fit to a benchmark and
presented as literature-derived.

The end state this project is converging toward: a scientifically
defensible simulation core (Track 1), with real statistical honesty
about its own uncertainty (Track 2), running on infrastructure that
could support real usage at scale (Track 3), optionally deepened with
genome-scale metabolic context where that adds genuine predictive value
rather than technological prestige (Track 4) — all built free/local-only,
with a real knowledge base behind it (Track 5) and a real way for
someone other than the person who built it to actually use it
(Track 6).

**Non-negotiable constraints that apply to every chunk below** (see
`CLAUDE.md`'s ARCHITECTURE section for the full, authoritative list):
free/local-only, never fabricate a citation or parameter, verify every
agent's output yourself before trusting it, TDD with a real RED/GREEN
proof, full suite green before commit, mandatory code-reviewer pass
(+security-reviewer for anything touching drugs/dosage/physiology/
genetics/calculations), a terse `CLAUDE.md` DONE entry per chunk, commit
only when asked to (or under an already-granted continuation), never
push without being told to.

---

## 1. Track summary

| Track | Theme | Status | Depends on |
|---|---|---|---|
| 1 | Biological realism | Active, not exhausted (new findings keep surfacing) | — |
| 2 | Statistical rigor | Named candidates done; extension work available | Track 1's pathway set |
| 3 | Architecture & infrastructure | Partially blocked (Docker), several ready chunks | — |
| 4 | Substrate/tissue-resolved metabolism | Reframing **adopted by the project owner 2026-08-10**; entry criteria set, 4.1 done | A flux-type benchmark (4.2) |
| 5 | Knowledge infrastructure | **Approved by the project owner 2026-08-10** -- active track, not proposed | — |
| 6 | Platform & distribution | **Resolved by the project owner 2026-08-17: parked, not part of the roadmap** | — |

**Track order discipline** (already established, keep following it):
Track 1 -> 2 -> 3 -> 4 -> 5, skip to the next only when the current one
is *genuinely* exhausted (every named candidate done or explicitly
resolved) — a track being *blocked* (e.g. Track 3's Docker dependency,
Track 4's licensing/benchmark prerequisites) is not the same as
exhausted, and does not license skipping ahead. Track 6 stays parked
(see below) -- do not pick it up as a filler even when 1-5 are blocked.

Track 5 was **approved by the project owner on 2026-08-10** and is now a real track,
on equal footing with 1-4 in the track-order discipline above (i.e.
Track 3 exhausted/blocked -> Track 4 exhausted/blocked -> Track 5, not
"only when 1-4 are ALL blocked").

**Track 6 resolved 2026-08-17 (item G2, process-improvement Section 10):**
The project owner confirmed Biosim's audience is a research/personal tool, not
clinician-facing or an educational demo. Track 6 (platform &
distribution -- SDK packaging, a demo/UI layer, user-facing docs) is
therefore **parked, not part of the roadmap**, not merely "provisional"
-- don't resurface it as a filler chunk. This also resolves row 3.14's
own priority: `physiological_profile` REST exposure is worth doing for
Track 3's own internal completeness/testability (see below), not
because an external client needs it.

**Validation ambition resolved 2026-08-17 (item G3):** the project owner confirmed
"internally consistent" (real citations, honest negative results, the
6 existing clinical benchmarks) remains the bar -- not "citable/
reviewable by others." This means Track 3's Docker-gated DB migration
and Track 2's remaining statistical-rigor extensions should NOT be
pulled forward ahead of their own natural place in the track order;
proceed through Tracks 1-5 in the established sequence, not reprioritized
around a publication goal that doesn't exist.

---

## 2. Track 1 — Biological realism

**Exhausted for the roadmap's originally-named candidates**
(SGLT2i, thyroid crosstalk, exercise, GLP-1 fixed-point, glucocorticoid
PK/PD, metformin PK — all done or explicitly resolved). **Not exhausted
overall** — two genuinely new findings surfaced this session alone
(first-phase secretion, glucagon), and the list below names several
more that haven't been scoped yet at all.

| # | Chunk | Size | Depends on | Notes |
|---|---|---|---|---|
| 1.1 | Glucagon fasting/basal EGP effect | **Permanently blocked** (real anchor found, not usable) | — | A fresh search angle found a real primary-source anchor -- Liljenquist et al. 1977 (PMID:833282): a >50% fasting glucagon reduction produces a 75% hepatic-glucose-production reduction, insulin held constant. NOT implemented: this platform's real target scenario (semaglutide's own measured -14% fasting glucagon effect, Hjerpsted 2018) is a 3.5x-smaller perturbation than the only quantitative anchor available -- extrapolating a linear gain across that range would fabricate precision this evidence doesn't support (no dose-response curve exists at small perturbations; a targeted search for one returned zero results). See `CLAUDE.md` "DONE (negative result) — Fresh literature search for the glucagon fasting/basal-EGP gain coefficient". Only revisit if a genuine multi-point, small-perturbation dose-response study is published. |
| 1.2a | Exogenous insulin — core mechanism (regular human insulin) | **DONE** | — | New stateless `ExogenousInsulinPathway` (`exogenous_insulin_pathway.py`), a closed-form two-sequential-equal-rate-constant subcutaneous depot model (no new ODE state), additive `exogenous_insulin_delta` rule-11 signal. … |
| 1.2b (analog) | Exogenous insulin — rapid-acting analog (insulin lispro) | **DONE** | 1.2a | Added `Drug.INSULIN_LISPRO`, calibrated from the SAME within-study comparison as regular insulin (Howey et al. 1994's own arms A vs. … |
| 1.2b (event-type) | Exogenous insulin — dedicated `InsulinDoseEvent` scenario type | **DONE** | 1.2a, 1.2b (analog, done) | New `InsulinDoseEvent` (genuine `dose_units` field, not `dose_mg`-as-Units), purely additive alongside the legacy `DrugEvent` path -- `ScenarioDefinition.__post_init__` rejects combining both for the same scenario. … |
| 1.2c | Exogenous insulin — needed for a genuine T1D disease state | **DONE (verification only)** | 1.2a (done) | Verified (not built -- no new code needed): `PathwayParameters.Ib`/`secretion_scale` are already free fields, and `basal_state()`'s existing fixed-point machinery already correctly derives a near-total-beta-cell-failure basal state for any Ib, out of the box. … |
| 1.3a | Sulfonylureas — core mechanism (glibenclamide) | **DONE** | — | The original row's "reuses `secretion_scale_modifier`" premise was found WRONG during scoping (would model the drug as suppressing insulin harder during hypoglycaemia — the opposite of its real risk). … |
| 1.3b | Sulfonylureas — glimepiride/glipizide | **DONE (glipizide AND glimepiride)** | 1.3a | **UPDATE (fresh-angle retry, resolved)**: glimepiride unblocked by searching for a head-to-head RATIO to glibenclamide's already-adopted fraction instead of an absolute magnitude. … |
| 1.3c | Hypoglycaemia-risk lay-output consumer | **DONE** | 1.3a | Added `ExpertMetrics.glucose_nadir_mM` + a "Hypoglycaemia risk" indicator to both `report.py`/`risk_score.py`, citing Battelino et al. 2019 (PMID:31177185). … |
| 1.4 | Alpha-glucosidase inhibitors (acarbose) | **DONE** | — | Modelled as a direct reduction of `MealInput.kabs` (mass-conserving by construction, unlike scaling `Ra(t)`'s output), calibrated against Chen et al. 2016 (PMID:27686734)'s dual peak-height/peak-time signature. … |
| 1.4b (i) | Acarbose follow-on: GI-adjusted-meal reach | **Investigated, rejected -- structural limitation** | 1.4 | A mathematically-safe margin-scaling approach to `GutAbsorptionModel.k_absorption` was designed and calibrated against Chen et al. 2016's real target -- but even at its extreme, the achievable effect (~92% peak/1.19x delay) falls far short of the real -28.5%/1.5x effect. Not implemented (would understate the real drug effect). See `CLAUDE.md` "DONE (negative result) — Investigated extending acarbose's effect to GI-adjusted meals". |
| 1.4b (ii) | Acarbose follow-on: miglitol | **DONE** | 1.4 | Added as a second alpha-glucosidase inhibitor, calibrated against Samad et al. 1988 (PMID:3042310, BAYm 1099 = miglitol's own development code, verified live) for a -21.0% peak-height reduction on a 25g sucrose meal. … |
| 1.4b (iii) | Acarbose follow-on: substrate-composition field | **DONE** | 1.4 | Added `MealEvent.complex_carb_fraction` (default 1.0, backward-compatible), wired into `runner.py`'s acarbose scaling, the REST wire schema, and `to_dict()` (for `SimulationCache` correctness). … |
| 1.5 | Tirzepatide (dual GLP-1/GIP agonist) | **DONE (partial) -- Gb-lowering component only** | — | Deep-search plan (2026-08-12), item A6: found and independently re-verified live Furihata et al. 2022 (a real, placebo-comparable T2D trial reporting per-arm baseline and week-8 fasting-plasma-glucose values at 5/10/15mg doses), from which each arm's own … |
| 1.6 | Liraglutide benchmark | Recurring, low-cost re-check -- re-confirmed still blocked (2026-08-12) | Europe PMC full-text access to PMID:21340616 | Permanently blocked until full text becomes available (subscription-required today) — re-check periodically, don't re-search from scratch each time. Deep-search plan (2026-08-12) re-check: searched ClinicalTrials.gov (no matching registration -- likely pre-dates FDAAA 2007 mandatory registration) and every companion Flint/Kapitza paper (none analyze the same n=18 dataset); re-confirmed the only secondary source (PMC4875959/PMID:26597252) still just paraphrases the 35%/39% figures without a primary-source page. No change from the prior blocked status. |
| 1.7 | GI-timing slope (composition-dependent absorption timing) | Permanent limitation | — | Two real studies contradict each other on this; documented as genuinely unknowable with current literature, not a bug. Only revisit if a new, well-powered multi-GI-level study is published. |
| 1.8 | Renal (SGLT2i) T2D-specific calibration | **DONE** | — | Both `for_dapagliflozin()`/`for_empagliflozin()` now take `t2d: bool`, using each drug's own cited T2D-cohort TmG (420/459 mg/min vs. 317/337 healthy) and, for empagliflozin, T2D-specific reduction fractions too. … |
| 1.9 | mTOR/AMPK real energy-sensing | Large, likely superseded by Track 4 | Track 4's FBA work, if it ever produces flux-level ATP data | Currently pure insulin read-outs, already honestly disclosed as such. A real fix needs an AMP:ATP-ratio-like quantity this ODE doesn't track — don't attempt inside Track 1; revisit only if Track 4 produces something usable. |
| 1.10 | New disease states: T1D | **DONE** | — | `DiseaseState.TYPE_1_DIABETES` added, calibrated against Simonson et al. 1987 (PMID:3100365, undetectable-C-peptide IDDM cohort) for Ib/secretion_scale (clamped to the personalizer's own engineering floors), plus a disclosed, non-literature-anchored Gb … |
| 1.11 | New disease states: gestational diabetes, MODY | **DONE (both)** | — | Gestational diabetes: `DiseaseState.GESTATIONAL_DIABETES`, calibrated against Kautzky-Willer et al. 1997 (PMID:9353615) for p1/p3/secretion_scale and the IADPSG consensus (PMID:20190296) for Gb. … |
| 1.12 | Periodic citation re-audit | Recurring, ~1 chunk each time -- now formalized, see 5.5 | — | Cadence and method now formalized in `CLAUDE.md`'s own "CITATION RE-AUDIT PROTOCOL" section (MASTER_PLAN 5.5). Eight passes so far: Track 2/3 range (clean); Track 1 range (clean, 21 PMIDs); `packages/metabolism` range (11 PMIDs -- found and fixed one real mistranscribed citation, but in this file's OWN prose describing the code, not in the source code itself); Track 5 range (clean, 8 PMIDs); post-Q3 range (clean, 15 PMIDs); post-a2b633d range (clean, 14 PMIDs -- CKD/pioglitazone/sleep-restriction/NEFA-EC50/exercise-disclosure citations, effectively a complete sweep of every new citation in the range rather than a partial sample); post-c645c7b range (clean, 10 of 28 sampled -- glipizide/prednisone/glibenclamide/pioglitazone-tier/intermittent-hypoxia citations, most of that range's headline PMIDs already independently re-verified live during their own chunk's implementation, so this sample deliberately covered different ground); post-82c862f range (**1 real fix** -- PMID:33090911's author byline in `modifiers.py`'s cold-acclimation comment was mistranscribed as "Sellers AN, Ganio MS, Cheek TF" when the real paper's authors are "Sellers AJ, Pallubinsky H, Rense P..." -- PMID/title/journal/year/quote were all already correct, only the names were wrong; 30 of 31 sampled PMIDs were clean). See `HISTORY.md`'s "DONE" entries under this row's name for full detail. Next audit due per the formalized trigger (10+ commits or 15+ new citations since the marker in the most recent entry, commit `276705a`). |
| 1.13 | Deep-search plan (2026-08-12): GLP-1 endogenous secretion dose-response | **DONE** | — | Found and adopted: Orskov et al. 1994 (PMID:8138058, live-verified twice), basal 7+/-1 pM -> peak 41+/-5 pM at 90min post-breakfast. … |
| 1.14 | Deep-search plan (2026-08-12): GLP-1 receptor EC50 for the incretin effect | Blocked -- real dose-response study found, concentration numbers paywalled | — | Kjems et al. 2003 (PMID:12540611, live-verified): graded GLP-1 infusion (0/0.5/1.0/2.0 pmol/kg/min) x graded glucose infusion, n=7 T2D + 7 controls -- confirms genuine dose-dependent amplification (up to 6x T2D, 11x controls) but reports only infusion RATES, not the resulting plasma GLP-1 CONCENTRATIONS needed for a concentration-based EC50. No open-access full text (Europe PMC confirmed); two companion papers by the same group (Vilsbøll 2003, PMID:12519856; Hare 2009, PMID:19837930) checked and also don't report the needed concentrations for this specific protocol. Do not repeat this search -- only a full-text-access route or a different primary source reporting BOTH infusion rate AND resulting plasma concentration would unblock this. |
| 1.15 | Deep-search plan (2026-08-12): GLP-1 gastric-emptying dose-response | Blocked -- real dose-response studies found, concentration numbers paywalled | — | Nauck et al. 1997 (PMID:9374685) and Meier et al. 2003 (PMID:12788879), both live-verified, real graded-infusion (0.4/0.8/1.2 pmol/kg/min) studies confirming dose-dependent gastric-emptying suppression -- neither paywalled abstract reports the specific plasma-concentration-vs-suppression numbers needed to fit an IC50. One single-point anchor found (Willms 1996, PMID:8550773: ~70 pM -> near-complete suppression in T2D) -- not enough alone for a curve fit. Do not repeat this search -- only full-text access to one of the two graded-dose studies would unblock this. |
| 1.16 | Deep-search plan (2026-08-12): GLP-1 gastric-emptying/incretin interaction architecture question | Blocked, 2 search passes tried and rejected | — | Does slowed gastric emptying SUPPRESS net insulin secretion rather than adding independently to it (would mean `GLP1Pathway`'s two-independent-signals design is wrong)? Already tried twice: both originally-cited papers (PMID:24464721, PMID:39418085) turned out not to support the claim on inspection; two follow-up searches (endogenous-GLP-1 dissection studies, exendin(9-39) antagonism) didn't converge. New angle: search using the statistical framing a mediation-analysis paper would use -- "GLP-1 gastric emptying insulin secretion mediation analysis", "GLP-1 receptor agonist structural equation model glucose insulin" -- or search specifically for a systematic review/meta-analysis on this mechanism rather than another primary study. |
| 1.17 | Deep-search plan (2026-08-12): post-exercise insulin-sensitivity magnitude | **DONE** | 1.16's own exercise scoping note | `max_p3_boost`'s coefficient recalibrated from an uncited 0.50 to 0.31, anchored to Mikines et al. 1988 (PMID:3126668, Km/Vmax clamp shift) and Wang et al. 2013 (PMID:23493357, Rd clamp shift), both live-verified. … |
| 1.18 | Deep-search plan (2026-08-12): during-exercise clearance/secretion-suppression magnitude | **Investigated, rejected -- wrong shape for both mechanisms** | — | Real citations found for both (DeFronzo 1981, PMID:7033285, exercise-alone glucose-uptake increment; Wahren 1971, PMID:5129319; Jarhult & Holst 1979, PMID:230455, ~50% insulin fall during exhaustive exercise) -- but neither transfers to this model's actual mechanism shape. `extra_g_clearance_modifier` multiplies a rate constant on `(G-Gb)`, exactly zero at G=Gb -- DeFronzo's number was measured AT clamp-held euglycemia, the one condition this term contributes nothing under (confirmed by direct execution: `dG=0.0` at y0 regardless of the modifier's value). `secretion_scale_modifier` multiplies `delta_secretion`/`first_phase_secretion`, both exactly zero at basal glucose/no meal -- Jarhult measured a fall in FASTED, basal plasma insulin, a quantity this parameter structurally cannot move. Both would need new architecture (an additive rule-11 signal for clearance; treating the postprandial-increment assumption explicitly for secretion) to use these citations at all -- a milestone, not a recalibration. See `CLAUDE.md` "DONE (negative result) -- During-exercise clearance/secretion citations found, wrong shape for the model". |
| 1.19 | Deep-search plan (2026-08-12): NEFA basal-vs-obese plasma concentration | **Formally PARKED** (2 independent search passes, both counter-evidence, per Section 9's … | — | Pass 1: best-titled candidate (Opie & Walfish 1963, PMID:13940209, exact topical/MeSH match) is a pre-1975 record with no abstract and no accessible full text (subscription required) -- could not extract a single number. … |
| 1.20 | Deep-search plan (2026-08-12): NEFA plasma half-life | **DONE** | — | Re-derived from an uncited 90.0 to 3.7, anchored to Heiling et al. 1991 (PMID:1951681, live-verified twice: "the half-life of plasma [14C]FFA was 3.7 +/- 0.1 min"), independently corroborated by Hagenfeldt et al. 1975 (PMID:1213651, implied t1/2 ~0.87-1.5 min … |
| 1.21 | Deep-search plan (2026-08-12): insulin IC50 for antilipolysis | **DONE** | — | Re-derived from an uncited 80.0 to 68.0, anchored to Stumvoll et al. 2000 (PMID:11061533, live-verified twice: a real 3-step clamp study directly reporting "EC50 of... adipose tissue [lipolysis]... 68 pmol/L"). … |
| 1.22 | Deep-search plan (2026-08-12): NEFA's own insulin-resistance dose-response | **DONE (fresh-angle retry resolved it)** | — | Original pass: Gormsen et al. 2007 (PMID:17341555) and Boden et al. 1994 (PMID:8200979) both real and on-topic, but neither reports its own basal/reference FFA level alongside its dose-response curve, so no excess-above-basal EC50 could be computed without an … |
| 1.23 | Warn when a safety clamp overrides a drug/disease modifier's calibrated effect | **DONE** | — | Process-improvement item C1 (2026-08-17 self-evaluation, prompted by a `/verify` pass on 1.5/tirzepatide). … |
| Fresh Track 1 scoping (2026-08-17/18) | A dedicated `computational-biologist` fan-out, dispatched per Section 9's own rule (re-scope before declaring a track exhausted a second time) after every named candidate above was DONE/blocked/parked. Found 3 real new candidates with extractable citations (rows 1.24-1.26 below) plus one useful negative result (pramlintide/amylin: real citations found, but both routes hit the same "wrong shape for this model's architecture" failure class as row 1.18 -- either no absolute denominator to compute a fractional Ra reduction (Hinshaw 2016, PMID:26930181), or a genuine gastric-emptying LAG PHASE ("t50 values could not be calculated" per Kong 1997, PMID:9028722) that neither `MealAbsorptionPathway`'s rate-constant model nor `GutAbsorptionModel`'s biexponential model can represent without new architecture -- not pursued, logged here so it isn't re-searched). Also checked and rejected: atypical antipsychotics (Sowell 2003, PMID:14671184, null result) and beta-blockers (Ferrara 1987, PMID:2891537, null result) -- neither yields a usable anchor. Track 1 is confirmed NOT exhausted. | | |
| 1.24 | Chronic kidney disease (CKD) disease state | **DONE** | — | New `DiseaseState.CHRONIC_KIDNEY_DISEASE`, calibrated against de Boer et al. 2016 (JASN, PMID:26823551, live-verified twice) -- a moderate-severe non-diabetic CKD cohort (n=59, mean eGFR 37.6 ml/min/1.73m2) vs. n=39 healthy controls. … |
| 1.25 | Thiazolidinediones (pioglitazone) | **DONE** | — | New `Drug.PIOGLITAZONE`, calibrated against Miyazaki et al. 2002 (Diabetes Care, PMID:11874940, live-verified twice) -- a 5-arm, 26-week, diet-treated T2D RCT (n=58, placebo/7.5/15/30/45mg). … |
| 1.26 | Sleep restriction as a physiological state | **DONE** | — | New `PhysiologicalProfile.sleep_hours: Optional[float]` (None = no adjustment, same convention as `time_of_day_hours`), calibrated against Buxton et al. 2010 (Diabetes, PMID:20585000, live-verified twice) -- 12-day inpatient crossover, n=20 healthy men … |
| 1.27 | Cite or disclose `post_exercise_modifiers()`'s remaining uncited p1/Gb coefficients | **DONE (partial disclosure, no code/behavior change)** | — | Found while re-checking whether Track 1 was genuinely exhausted after 1.24-1.26: `p1_boost` (0.10) and `gb_suppression` (0.03) sat uncited right next to the already-properly-cited `max_p3_boost` in the same function. … |
| 1.28 | Dietary protein co-ingestion (new personalization axis) | **DONE** | — | Found by a dedicated re-scoping fan-out (computational-biologist) after rows 1.24-1.27 -- this platform's food model was carbohydrate-only; protein's real, direct insulinotropic effect was entirely unmodeled. … |
| 1.29 | Polycystic ovary syndrome (PCOS) disease state | **DONE** | D1 composition-boundary audit (run) | From the same fan-out as row 1.28, deliberately sequenced behind the D1 audit. … |
| 1.30 | `EGP_MAX_MM_MIN`/`EGP_SUPPRESSION_N` citation-provenance gap | Blocked -- real candidate found, ambiguous without full text | Europe PMC/institutional full-text access to PMID:6761205 | From the same fan-out as row 1.28. `hepatic.py`'s `EGP_MAX_MM_MIN`/`EGP_SUPPRESSION_N` are described in `_solve_egp_suppression_ic50()`'s own docstring as "trustworthy literature anchors" but carry no citation at all -- the same row-1.27 pattern (a constant described as anchored while being uncited), here load-bearing since `EGP_SUPPRESSION_IC50` is derived FROM `EGP_MAX_MM_MIN` and every hepatic preset re-derives through it. Candidate: Rizza RA, Mandarino LJ, Gerich JE. Diabetes. 1982;31(8 Pt 1):663-9. PMID:6761205 -- saline-control-arm clamp data (n=6 normal volunteers) giving a half-maximal-EGP-suppression insulin concentration of 20+/-3 microunits/ml. Worked both plausible readings of "half-maximal response" through this platform's own formula: one reading CORROBORATES the current derived value (predicts ~24.9 microU/ml vs. Rizza's measured 20+/-3); the other would CONFLICT (implying `egp_max` should be ~0.048 mM/min, not 0.080). The abstract alone cannot distinguish which reading Rizza's own methods intend -- needs the full methods section (not yet obtained; Europe PMC/institutional access not yet tried). Either outcome is valuable and small in scope once unblocked: a first real citation (row-1.27-style docstring-only pass) or a genuine calibration correction. Do not adopt either reading without the full text -- guessing here would repeat the exact false-precision mistake this project has rejected elsewhere. |
| 1.31 | Alcohol's glucose mechanism -- fixed a real, confirmed bug (not just a citation gap) | **DONE** | — | A fresh Track 1 scoping fan-out re-checked `drug_interactions.py`'s own explicitly-named uncited placeholder (`alcohol_modifiers()`) and found it was actively WRONG, not merely uncited: it reduced `Gb` (fasting glucose CONCENTRATION), directly contradicted by … |
| 1.32 | Intermittent hypoxia -- acute p1/p3 personalisation axis | **DONE** | — | Track 1 priority candidate #2, following alcohol (row 1.31). … |

---

## 2a. Track 1 — Exhaustive Roadmap (compiled 2026-08-18)

Full-inventory planning pass covering every open, blocked, and unscoped
item in Track 1, plus a fresh forward-scoping sweep for genuinely new
candidates (a `computational-biologist` fan-out, every PMID below
independently live-verified twice via raw NCBI `efetch` and
`PubMedFetcher` before being written here). Supersedes no existing row
-- this section organizes and sequences what rows 1.1-1.32 already
established, plus adds new rows (1.33+) for what the fresh sweep found.

### Scheduled audits (status as of 2026-10-20)

Both recurring audits ran on 2026-10-20 at commit `2748de9` (see `HISTORY.md`):
the citation re-audit was clean (6/6 PMIDs), and the D1 audit found no clamp
flattening but surfaced row 1.46 below. Next due: citation re-audit after 10 commits
or 15 new citation lines; D1 after 5 new mechanisms.

### Recommended sequencing

```
1. ~~Run the D1 audit~~ **DONE 2026-08-18** -- no flattening found;
   still implement a p3/p1/ke/secretion_scale clamp_warnings extension
   before or alongside 1.33 (see Section 10/D1's own corrected note).
2. Run the citation re-audit (overdue, cheap, independent).
3. THEN implement, in this order:
   a. A1 -- HIV protease inhibitors (indinavir)   [row 1.33]
   b. A3 -- Thiazide diuretics (HCTZ)              [row 1.34]
   c. A4 -- Growth hormone (static peak factor)    [row 1.35]
   d. Formalize + independently re-verify atorvastatin [row 1.41]
      and nicotine [row 1.42] (informally scoped in an earlier
      session, never given their own row or independent re-check)
   e. A2 -- Menstrual cycle phase                  [row 1.36]
   f. 1.29 -- PCOS (already scoped, needs the D1-informed policy
      applied, then independent re-verification of its own citations)
   g. B4 -- Acromegaly (Gb-only)                   [row 1.37]
4. Decide-then-implement (needs a design decision, not just a citation):
   B1 -- Cushing's syndrome                        [row 1.38]
5. Needs its own further scoping/full-text access before viable:
   B2 -- OSA disease state (double-counting policy needed) [row 1.39]
   B3 -- Cold acclimation (T2D-only, regime scoping needed) [row 1.40]
   Cirrhosis (Kruszynska/Bergman, PMID:8458526, full text blocked) [row 1.43]
6. Ongoing, systemic, not mechanism-specific:
   Validation-benchmark coverage gap (below)
```

### New candidates found this pass (Tier A — strong anchor, clean mapping, ready after D1)

| # | Chunk | Size | Depends on | Notes |
|---|---|---|---|---|
| 1.33 | HIV protease inhibitors (indinavir) — acute p3 reduction | **DONE** | D1 audit (run) | **Noor et al. 2002, AIDS, PMID:11964551** (independently re-verified live via NCBI efetch before implementation) -- randomized double-blind crossover, n=6 healthy HIV-seronegative men, single 1200mg oral indinavir vs. placebo, 180-min euglycemic clamp: "M/I … |
| 1.34 | Thiazide diuretics (hydrochlorothiazide) — chronic p3 reduction | **DONE** | D1 audit (run) | **Pollare, Lithell & Berne 1989, N Engl J Med, PMID:2671740** (independently re-verified live before implementation) -- randomized double-blind crossover, **n=50** (largest, best-powered study in this whole Section 2a sweep), two 4-month treatment periods … |
| 1.35 | Growth hormone / somatropin — acute p3 reduction | **DONE (route a: static peak factor)** | D1 audit (run) | **Fowelin et al. 1991, Diabetologia, PMID:1916055** (independently re-verified live before implementation) -- n=11 healthy subjects, GH infused at 3 rates x 1h, euglycemic clamp: "insulin-antagonistic effect... started after approximately 2h, was maximal … |
| 1.36 | Menstrual cycle phase — new female physiology axis (p3, not p1) | **DONE (categorical, corrected from the original 3-point-ramp scoping)** | D1 audit (run) | **Valdes & Elkind-Hirsch 1991, JCEM, PMID:1997519** (independently re-verified live before implementation) -- n=8 cycling women, tolbutamide-modified IVGTT analyzed via the Bergman minimal model directly (this platform's own p1/p3 convention exactly): SI … |

### New candidates found this pass (Tier B — real citation, needs a design decision or further scoping first)

| # | Chunk | Size | Depends on | Notes |
|---|---|---|---|---|
| 1.37 | Acromegaly — Gb-only disease state | **DONE** | D1 audit (run); relationship to row 1.35 (GH) decided | Møller N, Schmitz O, Jørgensen JO, et al. "Basal- and insulin-stimulated substrate metabolism in patients with active acromegaly before and after adenomectomy." J Clin Endocrinol Metab. 1992;74(5):1012-9. … |
| 1.38 | Cushing's syndrome — implementation-route decision needed | **DONE (2026-08-19, third attempt) -- shipped on the verified-correct foundation** | Both blockers found during earlier attempts are now DONE: Section 8's `CORTISOL_P3_SENSITIVITY` item (unblocked the p3-magnitude gap) and Section 8's CRH/Cortisol coupled fixed-point fix (unblocked the structural stability gap, third fix attempt, fully … | **UPDATE 2026-08-19, third attempt (SHIPPED)**: implemented on the now-fixed CRH/Cortisol coupled foundation. … |
| 1.39 | Obstructive sleep apnoea as a chronic disease state | **DONE (2026-08-19)** | Implemented `DiseaseState.OSA` plus this platform's FIRST disease-state-vs-profile-field mutual exclusion. | **Harsch et al. 2004, Am J Respir Crit Care Med, PMID:14512265** -- n=40 AHI&gt;20, clamp before/2 days/3 months post-CPAP: "Insulin sensitivity significantly increased after 2 days... and remained stable after 3 months" -- untreated-OSA p3 factor … |
| 1.40 | Cold acclimation (T2D-specific, sensitising) | **DONE (2026-08-19)** | Implemented `PhysiologicalProfile.cold_acclimatized: bool` -- this platform's first p3-INCREASING environmental axis, gated with a REST-visible disclosure warning rather than the "regime scoping" this row was originally parked on. | **Hanssen et al. 2015, Nat Med, PMID:26147760** -- n=8 T2D, 10 days cold acclimation (14-15C): "increased peripheral insulin sensitivity by ~43%" -- p3 factor 1.43 (5.75/6.79 rounded), T2D-only (no healthy arm), re-verified live via `PubMedFetcher` before use. … |

### Formalizing two informally-scoped candidates from an earlier session

Neither of these has ever had a MASTER_PLAN row, an independent re-verification of its own citation, or a scoping-agent report reviewed against this file's Section 9 "verify the candidate is real" process -- both are pending exactly that, not pending new literature discovery.

| # | Chunk | Size | Depends on | Notes |
|---|---|---|---|---|
| 1.41 | Atorvastatin (statins) — dose-response p3 reduction | **DONE** | D1 audit (run) | Flagged by an earlier Track 1 scoping fan-out as priority candidate #3, never independently re-verified until now. … |
| 1.42 | Nicotine | **Formally PARKED (negative result)** | — | A dedicated `computational-biologist` scoping pass (its load-bearing citation independently re-verified live before accepting the recommendation) found the healthy-subject effect is a verified NULL, not a missing citation, across THREE independent controlled … |
| 1.44 | Hydroxychloroquine (antimalarial/DMARD) -- p3-RAISING | **DONE (2026-10-19; seeds 1acc6a3)** | Fresh Track 1 scoping fan-out | **Wasko et al. 2015, Diabetologia, PMID:26197707** (verified live, plus full text via NCBI PMC) -- randomized double-blind parallel RCT, HCQ 400mg/day vs. placebo, n=17/15 non-diabetic overweight/obese adults with >=1 IR marker, 13 weeks, IVGTT + Bergman … |
| 1.45 | Non-alcoholic fatty liver disease (NAFLD) | **DONE (2026-10-19; seeds 1acc6a3)** | Fresh Track 1 scoping fan-out | **Bugianesi et al. 2005, Diabetologia, PMID:15747110** (verified live) -- two-step euglycaemic clamp, n=12 biopsy-proven NAFLD (non-diabetic, normolipidaemic, normotensive) vs. n=6: "endogenous glucose production (basal and during the clamp) was normal … |
| 1.46 | T2D preset is nearly insulin-independent, so p3 mechanisms are inert in T2D | **(a) DONE 2026-10-23 (cited T2D `p1` 0.56); (b) structural insulin gap open, deferred by the project owner 2026-10-23 (chose "cited p1 fix only")** | D1 audit 2026-10-20 | **Update 2026-10-23**: (a) shipped -- T2D `p1` factor 0.56 (Basu 2009 SG 0.019 vs. 0.034, full text re-verified live; direction corroborated by Basu 1997, PMID:9032098; seeded GOLD `p1_factor`). Measured: glucose iAUC ratio 1.11 -> 1.51 (now in band, 1/5), insulin peak 0.26 -> 0.27, 2 h iAUC 0.11 -> 0.13, 6 h 0.12 -> 0.15, delay 27 -> 32 min; metformin 24.3 -> 44.2 (ref 186), steroid 8.3 -> 7.5% (ref 19.0); all verdicts unchanged. No double-counting: `HepGlucPathway` EGP depends on insulin only. (b) a T2D fasting-insulin level and a late-phase secretion term, a model-structure change, was not taken up; the project owner chose to leave it open. See `HISTORY.md`. Original finding: Found by the 2026-10-20 D1 sweep (no clamp involved). Over a 75 g meal, insulin-dependent disposal (X*G) is 3.4% of total disposal in the default T2D preset vs 45% when healthy; simulated T2D insulin peak is 102 pM vs 541 pM healthy -- mostly the delta-secretion baseline rising with Gb (6.63 mM), plus secretion_scale 0.6. Every p3-only mechanism (indinavir, HCTZ, atorvastatin, GH, HCQ, OSA, NAFLD, menstrual phase) moves T2D glucose AUC by 0.00-0.06% vs 0.15-2.4% when healthy. Question: is a ~100 pM postprandial peak realistic for this preset's intended T2D population? **Measured 2026-10-20**: `validate_t2d_meal_insulin()` (Basu 2009, PMID:19196896, ratios vs. matched controls) -- model T2D/healthy insulin peak 0.26 (ref 0.57), 2 h iAUC 0.11 (0.40), 6 h iAUC 0.12 (0.99), peak delay 27 min (94), glucose iAUC 1.11 (3.0); 0/5 within band. Real early T2D: equal fasting insulin, blunted/delayed early response, equal 6 h total. Candidate mechanisms (agent's inference, NOT yet verified): the fixed-point secretion gain falls with raised Gb at unchanged Ib, double-counting `secretion_scale`; first-phase gain zeroed at Gb >= 6.39 mM from an IVGTT finding; p1 untouched in T2D. Any change must re-check the steroid/metformin T2D benchmarks and keep `validate_y0` clean. **Scoping diagnostic 2026-10-20** (monkeypatched `disease_modifiers`, no source change; mechanism claims above re-checked in code: `F_secretion_basal = ke*Ib*(KG+dG)/dG` really does fall ~3.2x from Gb 5.0 to 6.5 at fixed Ib, and `first_phase_gain` is zeroed at Gb >= 6.39 mM from an IVGTT finding): the one directly cited, same-cohort candidate is a T2D `p1` factor of ~0.56 (Basu 2009 Sg 0.019 vs. 0.034 min-1). Applying it alone moves the benchmark's glucose iAUC ratio 1.11 -> 1.51 (band 1.5-4.5, just inside) and the metformin benchmark's reduction 24 -> 44 (ref 186), the steroid point 8.4 -> 7.5%; insulin ratios barely move (2 h iAUC 0.11 -> 0.13). Adding `secretion_scale` 0.6 -> 1.0 on top gives 2 h iAUC 0.22, 6 h 0.24, peak 0.35, delay 32 min -- still far from 0.40/0.99/0.57/94. So no single anchored parameter closes the insulin gap: it is structural (delta-secretion has no late/potentiated component; Ib is held at 66 pM in T2D so fasting insulin cannot rise as in Basu 59 vs. 37 and Smith 86.7 vs. 52.5). Next step is a design decision, not a retune: (a) a cited `p1` T2D factor as its own milestone (check Pacini-type T2D Sg values for consistency first; composition/double-counting with the Gb elevation; re-run steroid/metformin benchmarks), and (b) whether to give the T2D preset its own fasting-insulin level and a late-phase secretion term -- a model-structure change needing the project owner's call. Do not retune secretion until a cited, double-counting-checked fix is scoped. |

### Documented negatives from this pass (checked, real search, not viable -- do not re-search without a new lead)

| Candidate | Anchor checked (live-verified) | Why rejected |
|---|---|---|
| Immunosuppressants (tacrolimus/ciclosporin) | Øzbay et al. 2012, PMID:21988494 | Direction OPPOSITE to hypothesis in the gold-standard acute healthy-volunteer data ("increased by 25%... during CsA", "13% during Tac") -- contradicts the chronic post-transplant clinical picture; no usable chronic anchor separable from transplantation/steroids. |
| Chronic pancreatitis | Ciccarelli et al. 2024, PMID:38871564 | n=146, properly stratified by OGTT status: "no differences in beta-cell function or in insulin sensitivity" vs. non-pancreatitis controls at matched glucose tolerance. No parameter delta to implement beyond the existing T2D/T1D presets. |
| Oral contraceptives | Skouby et al. 1987, PMID:3102539 | Null in the healthy population this platform is calibrated on (P=0.09, n=6); only significant in a previous-GDM subgroup, and even there "not sufficient to alter glucose tolerance." |
| Atypical antipsychotics | Hardy et al. 2011, PMID:21435142 | Significance flips direction between insulin-clamp steps for the two drugs tested (not a coherent dose-response); the source's own text attributes the effect to weight/adiposity change, not glucose directly -- would double-count `bmi_modifiers()`. Corroborates the pre-existing rejection (Sowell 2003, PMID:14671184). |
| Cirrhosis / liver disease | Kruszynska, Harry, Bergman &amp; McIntyre 1993, PMID:8458526 | **The single highest-value blocked item found this pass** -- a genuine minimal-model SI/SG IVGTT study (Bergman himself a co-author) that would map DIRECTLY onto this platform's own p1/p3 convention, but the PubMed abstract is truncated before the SI/SG values and Europe PMC confirms no open-access full text exists. Same "blocked on full text" class as rows 1.6/1.14/1.15/1.30 -- unblock route is institutional access to *Diabetologia* 1993;36(2):121-8, not further searching. See row 1.43. |
| Dietary fat co-ingestion | Gentilcore et al. 2006, PMID:16537685 | Confirms the pre-existing "blocked on missing T50" status from row 1.28's own scoping -- all P-values, zero magnitudes reported; also T2D-only and a 30-min preload rather than true co-ingestion. |
| Hypothyroid hepatic-EGP pathway (thyroid crosstalk asymmetry) | No human anchor found (2 search angles, only rodent studies) | The known asymmetry (hyperthyroidism has a real `HepGlucPathway` coupling, hypothyroidism is static-only) has no available human tracer data to close -- recommend documenting as a permanent, accepted limitation rather than continuing to search. |
| Niacin | Poynten et al. 2003, PMID:12800094 | Real, on-topic, mechanistically appealing (could route through the existing `NEFAPathway`/`NEFA_IR_EC50_UM` rather than a new modifier) but no extractable GIR magnitude or NEFA concentration in the abstract. Specific unblock: any niacin study reporting both ΔNEFA (µmol/L) and ΔGIR in one cohort. |
| Cystic fibrosis-related diabetes | Bolduc et al. 2024, PMID:39118592 | Reports hazard ratios/survival curves and Stumvoll-index quartiles only -- nothing mappable to p3/secretion_scale/Gb directly. |
| Coeliac disease | No relevant human study found | Only a mouse gliadin/high-fat-diet study surfaced. |
| Chronic/heavy alcohol as a distinct disease state | Nothing usable (rodent, or confounded with pancreatitis/withdrawal) | Row 1.31 already correctly resolved ACUTE alcohol; a chronic state would need to disentangle from cirrhosis/pancreatitis, both themselves blocked above. |
| Adipose tissue/leptin signalling as its own pathway | Only generalized-lipodystrophy metreleptin literature found | Not a general-population axis -- CLAUDE.md's own "adipose/leptin signalling" aspiration should be reframed or formally closed, not left open; the literature doesn't support what that aspiration implies. |
| Chronic altitude / non-acute hypoxia | One hit, itself still "cyclic"/intermittent (PMID:24028640) | Lands inside row 1.32's own already-disclosed Mackenzie-2012 caveat territory rather than being a genuinely distinct chronic-residence axis. |
| Beta-blockers | Re-confirmed top hit is the already-rejected Ferrara 1987 (PMID:2891537, null) | No new candidate found; existing rejection stands. |
| PCSK9 inhibitors | Zero design-matched results | Nothing to propose. |

### Systemic gap: validation-benchmark coverage

Only **9** of the platform's ~20+ literature-calibrated Track 1
mechanisms have a registered `ModelValidator` clinical benchmark
(`validate_ogtt`, `validate_bergman_iv`, `validate_semaglutide`,
`validate_steroid_meal_t2d`, `validate_metformin`, `validate_sglt2i`, and since
2026-10-20 `validate_protein_coingestion` -- Thondre 2021, which measures a large
glucose undershoot because gastric-emptying slowing is not modelled;
`validate_t2d_meal_insulin` -- Basu 2009, row 1.46; and `validate_caffeine_ogtt` --
Graham 2001 / Robinson 2004, which measures the wrong sign for healthy and ~no effect
for T2D). Caffeine therefore HAS a benchmark now; the "NO registered benchmark" list
below no longer applies to it.
Every other DONE mechanism (exogenous insulin, sulfonylureas, acarbose/
miglitol, tirzepatide, CKD, pioglitazone, sleep restriction,
alcohol, intermittent hypoxia, T1D/GDM/GCK-MODY/T2D/
insulin-resistant/hypothyroidism disease states, caffeine) has NO
registered benchmark -- meaning the D1 audit protocol (driving the real
REST server across a dose/disease range) is currently the ONLY
systematic check that a calibrated effect actually reaches simulated
output for any of them. Worth a dedicated future chunk: survey which of
these have a same-meal drug-vs-placebo (or equivalent) shape reusable
as a `validate_*()` benchmark the way `validate_sglt2i()` was built from
Polidori 2013, and prioritize the ones composing with the most other
mechanisms (pioglitazone/CKD/sleep-restriction/intermittent-hypoxia --
all `p3`-touching, all candidates for the same D1 sweep above).

---

## 3. Track 2 — Statistical rigor

**Named candidates (Bayesian calibration, population demographics,
genetic variant sampling, Sobol sensitivity) all done.** Extensions
available, several surfaced directly by this session's own work:

| # | Chunk | Size | Depends on | Notes |
|---|---|---|---|---|
| 2.1 | Wire `_compute_lay_summary()` to actually read `UQResult` | **DONE** | — | `is_high_confidence`/`confidence_note` now also downgrade when a Monte Carlo ensemble is unreliable (high sample-failure rate) or shows a wide glucose confidence band (`> 0.15` mean relative width, an engineering threshold), not just the point-estimate … |
| 2.2 | Sample extension-pathway-specific parameters in UQ | **DONE** | — | New additive `extension_distributions`/`extension_pathway_factory` pair on `propagate()`. … |
| 2.3 | `UQResult.species_pathway` provenance | **DONE** | — | `species_pathway: dict[str, str]` (bare species -> producing pathway name) + `pathway_for()`, computed by `MonteCarloEngine.propagate()` mirroring `ODESystem.simulate()`'s own prefix-stripping. … |
| 2.4 | Extend Sobol/OAT sensitivity to extension-pathway parameters | **DONE** | — | Same additive `extension_pathway_factory` pattern as 2.2. … |
| 2.5 | Incremental/streaming Bayesian re-calibration as new parameters are curated | Large, needs its own scoping pass | Track 5's curation pipeline maturity | Speculative — only worth scoping once Track 5 produces a real stream of new literature-anchored parameters to recalibrate against. |

---

## 4. Track 3 — Architecture & infrastructure

Two milestones done (BatchRunner/PopulationRunner multiprocessing), the
`ParameterStore` cross-scope bug fixed, DB migration scoped with two
open questions resolved and hygiene fixed — but genuinely blocked on
one external precondition (Docker Desktop not running in this
environment). Several *unblocked* chunks remain, found this session:

| # | Chunk | Size | Depends on | Notes |
|---|---|---|---|---|
| 3.1 | Fix `ParameterStore.query()`/`.best()` tie-order bug | **DONE** | — | `query()` now uses the same explicit `self._insertion_seq` tiebreak as `get_all_for_subject()`, replacing the old `_priority_score()`-only sort (deleted as dead code). … |
| 3.2 | Wire `CouplingSpec.validate()` into production | **DONE** | — | `CompositePathway.__init__` now calls `validate()` for every coupling, raising at construction (matching this project's fail-fast convention). … |
| 3.3 | `IDRegistry.from_mapping_files()` transitive-merge fix | **DONE** | — | `from_mapping_files()` now restores `_deprecated_index` per-file via `.update()`, mirroring `from_dict()`'s equivalent fix -- a transitive merge chain (A->B->C) no longer loses the A->B hop when loaded through this entry point. … |
| 3.4 | Neo4j invariant-4 decision | **DONE** | — | Made Python's already-deterministic insertion order an explicit `rank` property on `IDENTIFIES` edges; `cross_reference()` orders by `r.rank` (not accession string, which was silently wrong); invariant 4 now checks rank ambiguity, not multiplicity. … |
| 3.5 | Disposable Docker validation boot | 1 chunk, **blocked on Docker Desktop's engine actually starting** | Docker Desktop working locally | `docker compose up -d` -> exercise the 3 SQL functions by hand -> `docker compose down -v`. Predicted finding: `002_query_api.sql`'s functions likely need a per-function `SET search_path` given the documented `SET ROLE biosim_writer/biosim_reader` connection pattern. As of 2026-08-05: the project owner started Docker Desktop, but its engine still fails with "Docker Desktop is unable to start" (`docker info`) -- this looks like a Windows backend/virtualization issue, not just "not launched." Re-check `docker info` at the start of each session rather than assuming either state. |
| 3.6 | Backfill script — snapshot reader/validator (no SQL) | **DONE** | — | `biosim_parameter_db.snapshot_validator.validate_snapshot()` -- file-path-only input, copy-then-read, checks seq ordering + dangling supersession/conflict references. … |
| 3.7 | Backfill script — SQL/Cypher emission | 1-2 chunks | 3.5 (boot must happen first) | The actual `INSERT`/`MERGE` statement generation, informed by whatever 3.5's boot surfaces (generated-column dedup via `ON CONFLICT ... source_key`, `seq` insertion-order handling, self-FK ordering for `supersedes_parameter_id`). |
| 3.8 | Shadow-read comparison harness (migration Step 3) | Large, multi-chunk | 3.5, 3.6/3.7 | Run every test + every clinical benchmark twice (in-memory vs. read-only Postgres/Neo4j adapter), diff field-by-field. This is the actual equivalence gate — budget several chunks: build the adapter, run against the full suite, run against benchmarks, then a soak period under Step 4's dual-write. |
| 3.9 | Dual-write (migration Step 4) | Large | 3.8 | Every `add()`/`register()` call also writes to the new stores, behind a feature flag defaulting off. |
| 3.10 | Cutover (migration Step 5) | Gated, not a single chunk | 3.8, 3.9, all 7 open questions resolved | All-or-nothing checklist already specified in `infra/MIGRATION_PLAN.md`. |
| 3.11 | REST API surface growth (`/population`, `/batch` endpoints?) | **DONE** | — | Implemented the bounded-synchronous size-cap option. … |
| 3.11b | Expose `ExpertMetrics.glucose_nadir_mM` and lay-output risk indicators (`assess_risk()`/`compute_risk_profile()`) over REST | **DONE** | — | Added `glucose_nadir_mM` to `ExpertMetricsOut` and four new Pydantic models (`RiskProfileOut`/`RiskIndicatorOut`, `MetabolicRiskProfileOut`/`MetabolicRiskIndicatorOut`) exposing both functions' output on `SimulateResponse` (and therefore `/compare`, which … |
| 3.12 | SDK packaging | Gated, re-check periodically | REST API surface stabilizing (schema churn slowing down) | Checked this session: still appropriately deferred, wire schema still changing alongside biological milestones. Re-check every few sessions, not every session. |
| 3.13 | REST error-response serialization crashes on a raw NaN JSON literal (any float field, not insulin-specific) | **DONE** | — | Added `_sanitize_non_finite()` plus a custom `RequestValidationError` handler in `app.py` that mirrors FastAPI's own default handler exactly, except non-finite floats in the (rejected) input values are replaced with a JSON-safe string before serializing -- a … |
| 3.14 | Expose `physiological_profile` (disease states, age/sex/BMI, fitness) over REST | **DONE (age/sex/BMI/fitness/disease-state only)** | — | Process-improvement item C2. Added `PhysiologicalProfileIn` (schemas.py) and `_physiological_profile_in_to_profile()` (service.py), wired into `ScenarioIn`/`scenario_in_to_definition()`. … |

---

## 5. Track 4 — Substrate/tissue-resolved metabolism (reframing adopted 2026-08-10)

**Adopted by the project owner on 2026-08-10.** This session's dual scoping pass
(architect + computational-biologist) concluded the naive design (FBA
as a normal `CouplingSpec` extension pathway) is architecturally wrong
for five independent reasons, and that the track's own framing ("hybrid
ODE/FBA") should be corrected to something biology-first before any
code is written. The project owner approved the reframe as scoped: the track is now
named for its biology, not its solver, and genome-scale FBA is one
*candidate implementation* among several, not the track's definition.

**Entry criteria for any future chunk under this track** (all three
required, per the adopted reframe):
1. A named biological question the ODE paradigm structurally cannot
   answer, whose answer is concentration-**independent** (FBA produces
   fluxes, not concentrations — see the scoping's own finding that no
   real hybrid model in the literature can report metabolite
   concentrations from the FBA layer).
2. A flux/partition-type validation target (a tracer study, indirect
   calorimetry RQ, ¹³C-MRS glycogen measurement, etc.) — never a
   concentration/timing curve, which this track structurally cannot
   validate against (all 6 of Biosim's own registered benchmarks are
   concentration/timing measures).
3. A compute budget compatible with the existing test suite and the
   platform's free/local-only constraint (decision 2 in CLAUDE.md's own
   strategic roadmap) — no commercial solver, no cloud compute.

| # | Chunk | Size | Depends on | Notes |
|---|---|---|---|---|
| 4.1 | Adopt the reframing | **DONE** | — | Approved by the project owner 2026-08-10. Track renamed, entry criteria adopted verbatim as scoped (see above). See `CLAUDE.md` "DONE — Adopted the Track 4 reframing". |
| 4.2 | Acquire a real flux-type validation benchmark | **DONE** | 4.1 | No single clean "maximal EGP rate" number found in the accessible literature (the strongest lead, Matsuda/DeFronzo 2002 PMID:12200754, is paywalled) -- documented honestly, not forced. … |
| 4.3a | Stoichiometric model dataclass + `scipy.optimize.linprog`/HiGHS FVA wrapper | **DONE** | 4.1 | New `packages/metabolism` package: `StoichiometricModel` (NaN/shape/bounds-validated, defensively copied + read-only after construction) and `flux_variability_analysis()` -- deliberately objective-free (see Track 4's own scoping finding that a biomass … |
| 4.3b | Curate a small (~20-60 reaction) hepatic core network | **DONE** | 4.3a | 32-reaction, 30-metabolite FASTED-STATE-ONLY hepatic network (`fasted_hepatic_network.py`) -- glycogenolysis, gluconeogenesis from lactate/glycerol/alanine, minimal TCA, beta-oxidation, urea-cycle N-disposal. … |
| 4.3c | Run FVA, bound `hepatic.py`'s `EGP_MAX_MM_MIN` | **DONE (sensitivity table, not a point estimate)** | 4.2, 4.3a, 4.3b (all DONE) | A fresh literature search found real, usable ENERGY-side bounds (fatty acid availability, O2 consumption) that were missing -- but running the network with them applied revealed they don't end up binding: with glycogen unbounded, energy is irrelevant … |
| 4.3d | Extend the hand-curated hepatic network with PPP + de novo lipogenesis | **Scope corrected -- NOT built as originally planned** | — | A dedicated scoping pass (`computational-biologist` agent, its two most load-bearing claims independently re-verified live before accepting) found the original plan doesn't survive contact with the primary literature. **F1's own benchmark is not usable as a validation target for this network**: verified directly from Rosqvist et al. 2019's own Methods section (re-fetched independently) that the 2H2O label is dosed "the evening before... over the course of the evening" alongside "a low-fat dinner" -- isotope incorporation into VLDL-TG palmitate accumulates over that whole ~12h window, and only the final blood draw happens after an overnight fast. "Fasting hepatic DNL" describes when the SAMPLE was taken, not when the labeled synthesis occurred -- a real mismatch against this network's own specifically fasted-state scope. `ROSQVIST_2019_FASTED_HEPATIC_DNL`'s own caveat (which previously wrongly called this "a genuinely FASTED-state measurement... unlike most DNL literature") corrected to disclose this. Separately, and independently: (a) building ACC/FAS alongside the existing FAO reaction violates this network's own documented futile-cycle exclusion policy -- malonyl-CoA is simultaneously FAS's product and CPT1's inhibitor (Swanson/McGarry 1998, PMID:9794789, independently re-verified live), so both pathways running at appreciable flux is regulatorily self-contradictory, not just numerically wasteful; measured directly (scratch reconstruction) that adding them inflates FAO's currently literature-pinned FVA range 3.6-5.8x; (b) the row's own driving question ("is DNL NADPH-limited by PPP flux?") is structurally unanswerable by FVA either way -- with oxidative PPP as the sole NADPH source the answer is fixed by arithmetic alone (no LP needed); with an alternative source (e.g. malic enzyme) present, FVA reports the full feasible range and resolves nothing, since flux-partition questions need isotopomer data, not constraint-based bounds. No code written -- this is a scoping-only pass. **Replacement plan, split into two smaller, better-scoped items**: 4.3d-i (acquire a real ABSOLUTE fasted VLDL-TG/palmitate secretion-RATE benchmark, `FluxBenchmark` not `FractionalContributionBenchmark` -- needs its own fresh literature search, not yet done) and 4.3d-ii (once 4.3d-i exists: an `egp_sensitivity.py`-style bracketed table treating an imposed DNL flux as a fixed input, asking what it costs `EGP_MAX_MM_MIN`'s own ceiling -- ~8-12 reactions, FAO/FAS run as mutually-exclusive named scenarios mirroring the existing glycogen-depletion scenario convention, not simultaneously). Neither sub-item attempted yet. See `HISTORY.md` "DONE (scoping only, no code change) — Corrected the scope of PPP/DNL reaction additions (4.3d)". |
| 4.3d-i | Acquire a real ABSOLUTE fasted VLDL-TG/palmitate secretion-rate benchmark | **Blocked -- real candidate found, not usable as a DNL-specific rate without an unsupported cross-study combination** | — | Barrows & Parks 2006 (PMID:16449340) re-checked for an absolute rate: still fractional-only in the abstract, full text paywalled -- unchanged from 4.3d's own scoping pass. A fresh search found a real, independently-verified absolute VLDL-TG secretion rate: Mittendorfer, Patterson & Klein 2003 (PMID:12600845, *Am J Clin Nutr*, live-verified twice -- title/authors/exact numbers match): "[2H5]glycerol bolus + compartmental modeling, lean men 2.6+/-0.3, lean women 5.1+/-0.7, obese men 6.8+/-0.5, obese women 5.0+/-0.5 umol.L plasma-1.min-1" (basal state). **Not adopted**: this measures TOTAL VLDL-TG output, not DNL-specific flux -- deriving a DNL rate would require multiplying it by `ROSQVIST_2019_FASTED_HEPATIC_DNL`'s own 6.5% fraction, a different study, different cohort (US lean/obese men+women vs. a pooled 5-cohort European analysis), with no shared basal reference point -- exactly the cross-study-combination anti-pattern this project has already rejected once (`DROUIN_1998_GLUCAGON_STIMULATED_INCREMENT`'s own "do not add directly" caveat) and repeatedly elsewhere (the semaglutide-benchmark's own history). Additionally, this specific paper's "basal" label could not be independently confirmed as genuinely fasted/postabsorptive from ITS OWN full text (paywalled) -- only inferred from the same research group's companion papers using identical terminology. Still open -- needs either a single study reporting BOTH an absolute rate AND a DNL fraction together, or a different, DNL-specific absolute-rate paper not yet found. See `HISTORY.md`'s own entry for the full search trail (Barrows/Parks re-check, Mittendorfer/Patterson/Klein, and ~8 other candidates checked and rejected). |
| 4.3d-ii | "G6P cost of an imposed fasted DNL flux" sensitivity table | Small-medium, ~8-12 reactions | 4.3d-i | Once 4.3d-i supplies a real absolute DNL flux rate: treat it as a PINNED INPUT (never a predicted output, avoiding 4.3d's own "FVA can't resolve a partition question" finding), and ask what it costs `EGP_MAX_MM_MIN`'s own ceiling -- an `egp_sensitivity.py`-style bracketed table, matching 4.3c's own already-established pattern. Reaction set (already balanced/verified during 4.3d's scoping pass, ready to reuse): oxidative PPP (2 reactions) + a lumped non-oxidative branch (1) + ACC/FAS (2, lumped or separate) + TAG synthesis/export (2) + optionally malic enzyme (2, to avoid pre-answering the NADPH-partition question). Run FAO and FAS as MUTUALLY EXCLUSIVE named scenarios (never simultaneously), mirroring how `egp_sensitivity.py`'s existing `_SCENARIOS` already handles glycogen depletion -- do not simply add both reactions unconditionally, which would recreate 4.3d's own documented futile-cycle/malonyl-CoA-regulatory-conflict finding. |
| 4.4 | Adopt a real tissue GEM (Human-GEM-derived hepatocyte, or iHepatocytes2322) | **DONE (scoping only, deferred -- negative result)** | 4.3 (done), root `LICENSE` (done, 7.1) | Both dependencies now resolved, but a dedicated scoping pass found no genuine driving question that clears Track 4's own entry gate: every concentration-independent question a full GEM would newly answer (NADPH-limited de novo lipogenesis, serine/one-carbon … |
| 4.5 | Precomputed surrogate flux map as a smooth extension pathway | Large, multi-chunk | 4.4 | First point at which a new rule-11 signal key and full basal-invariance/fixed-point testing become necessary for this track. |
| 4.6 | Opt-in `HybridODEFBASystem` (dFBA-SOA orchestrator) | Large, multi-chunk | 4.5, and only if 4.2's benchmark demands a genuinely new prediction class | Not justified by any existing benchmark improving — only build if a real new output (lactate appearance, substrate-oxidation partitioning, ketogenesis) is actually wanted and validatable. |

---

## 6. Track 5 (approved 2026-08-10) — Knowledge infrastructure

The `literature-pipeline`/`curation`/`parameter-db`/`id-resolver`
packages already exist and have real, if narrow, seeded content (14
proteins/metabolites, a handful of literature-anchored parameters).
Growing this deliberately is a distinct effort from Track 1's
mechanism-by-mechanism additions.

| # | Chunk | Size | Depends on | Notes |
|---|---|---|---|---|
| 5.1 | Systematic parameter-breadth curation campaign | **Ongoing -- chunks 1-38 DONE (latest: hydroxychloroquine/NAFLD, `1acc6a3`)** | — | Registers each modelled constant in `ParameterStore`/`ContextRegistry` with its citation and tier, cross-checked live against the code constant (drift tests). One `HISTORY.md` entry per chunk ("Curated ... (5.1 chunk N)"); the old per-chunk log is in `git show 177a782:MASTER_PLAN.md`. Next: any modelled constant still without a seed entry. |
| 5.2 | `CurationQueue` workflow completion | **DONE (scoping only)** | — | Confirmed by direct code/grep inspection (no code changed): NO real end-to-end curator workflow exists today, only a well-built, well-tested data structure plus one real bridge function. … |
| 5.2a | Wire `ParameterStore.unresolved_conflicts()` into `CurationQueue.enqueue_conflict()` | **DONE** | 5.2 | New `biosim_parameter_db.curation_bridge.enqueue_unresolved_conflicts(store, queue, registry)`, mirroring `biosim_uncertainty.curation_bridge`'s own explicit-call (not auto-hooked into `add()`) shape. … |
| 5.2b | Wire `prioritise_from_sensitivity()` into a real orchestration entry point | **DONE** | 5.2 | Design decision made directly (opt-in method, not automatic/hooked into `simulate()`): new `PlatformSession.prioritize_curation_from_sensitivity()`, the first real caller. … |
| 5.2c | Wire `IngestionPipeline`'s `needs_curation` outcomes into `CurationQueue.enqueue()` | **DONE** | 5.2 | New `biosim_literature.curation_bridge.enqueue_needs_curation_outcomes()`, mirroring 5.2a/5.2b's own explicit-call shape. … |
| 5.2d | Build a minimal `CurationQueue` consumer (dashboard/report/CLI) | **DONE** | 5.2a/b (now done; 5.2c not needed as a prerequisite) | New `CSVExporter.export_curation_queue()` (mirrors `export_metrics_table()`'s own shape) plus `PlatformSession.export_curation_queue_csv()`. … |
| 5.3 | `LLMExtractor` real-world validation | Parked | **the project owner's explicit "keep dev free" decision** | Already correctly gated behind a real `ANTHROPIC_API_KEY` that the project owner has declined to provide, for cost reasons. Respect this — do not re-raise it as a blocker to resolve, only as an option if he changes his mind. |
| 5.4 | `IDRegistry` entity/pathway breadth | **DONE (first chunk: Track 4's metabolite set)** | — | Registered the 22 metabolites `biosim_metabolism.fasted_hepatic_network` needs (4.3b's own explicitly-deferred, non-blocking follow-on) with real ChEBI/KEGG accessions -- all 22 of each independently verified live by three separate parties this session … |
| 5.5 | Formalize the periodic citation re-audit as a named, scheduled activity | **DONE** | — | Added a "CITATION RE-AUDIT PROTOCOL" section to `CLAUDE.md`: an explicit trigger (10+ commits or 15+ new/changed citation lines since the last audit marker, whichever first) and a codified method (identify touched files since the marker commit, sample rather … |

---

## 7. Track 6 (parked 2026-08-17) — Platform & distribution

**Parked, not part of the active roadmap** (item G2, process-improvement
Section 10): the project owner confirmed Biosim's audience is a research/personal
tool. Rows 7.3-7.5 below (SDK packaging, user-facing docs, a demo/UI
layer) have no real driver under that audience and should NOT be picked
up as filler even when Tracks 1-5 are all blocked -- this is a change
from the prior "proposed, ask again later" status. 7.1/7.2 (licensing,
Python version) were genuine, already-completed housekeeping and stay
done regardless.

| # | Chunk | Size | Depends on | Notes |
|---|---|---|---|---|
| 7.1 | Root `LICENSE` file + `pyproject.toml` `license` field | **DONE** | — | the project owner chose Apache-2.0 (2026-08-10). **Correction**: this row's own earlier claim ("no root LICENSE file") was wrong -- a root `LICENSE` (MIT, commit `951c351`, from an early README/LICENSE bootstrap) already existed; nobody checked `git log`/`ls` before … |
| 7.2 | `pyproject.toml` Python version fix (`^3.11` -> actual 3.14) | **DONE** | — | Changed `python = "^3.11"` to `"^3.14"`, matching the actual installed/tested environment (confirmed `python --version` -> 3.14.5, and no `poetry.lock` exists -- this constraint was purely declarative metadata, not read by any live tooling, so the fix is a … |
| 7.3 | SDK packaging | **Parked with the rest of Track 6** -- no driver under a research/personal-tool audience | 3.12's gate | Still cross-referenced from 3.12 for visibility, but not worth pursuing given G2's resolution. |
| 7.4 | User-facing documentation (a real "how to use Biosim" guide, not just `CLAUDE.md`/`README.md`) | **Parked with the rest of Track 6** | — | `CLAUDE.md` is session-memory for Claude, not a user guide -- but the only user is the project owner himself, so this has no real driver either. |
| 7.5 | Demo/UI layer | **Parked with the rest of Track 6** -- explicitly speculative, now formally resolved as … | — | No stated need for this exists, and now won't by default -- only revisit if the project owner explicitly names a concrete use case, overriding G2's resolution. |

---

## 8. Cross-cutting housekeeping backlog

Small, independent, low-risk items surfaced across multiple scoping
passes — pick any of these up as a cheap filler chunk when nothing
larger is unblocked:

- `biosim_ode/__init__.py`'s `__all__` omits `CortisolPathway`/ `HepGlucPathway`/`NEFAPathway`/`GLP1Pathway` despite importing them **DONE** -- was actually missing 14 names, not 4 (every extension pathway's own `Parameters` class plus every `couple_with_gi`/ …
- 7.1/7.2 above -- both done.
- `ScenarioIn.n_mc_samples`/`ScenarioDefinition.n_mc_samples` have no upper bound **DONE** -- see `CLAUDE.md` "DONE — Scoped and fixed the `n_mc_samples` ceiling; mandatory security review found and fixed a WORSE, unrelated CRITICAL DoS vector on the same …
- CSV/formula-injection observation **DONE** -- see `CLAUDE.md` "DONE — Fixed CSV/formula-injection in CSVExporter and two sibling exporters". …
- Confidence-tier definition inconsistency across two files **DONE** -- see `CLAUDE.md` "DONE — Resolved the confidence-tier definition inconsistency". …
- **DONE** -- GOLD-vs-BRONZE tier inconsistency across the 5.1 curation campaign for structurally similar "ratio computed from a paper's own two reported numbers" cases: pioglitazone's `p3_multiplier_at_45mg_reference_dose` (chunk 9) was BRONZE despite a raw …
- **DONE (2026-08-19)** -- `cortisol_pathway.py`'s `p3_modifier = max(0.3, 1.0 - CORTISOL_P3_SENSITIVITY * excess_gr)` had a hard, UNREACHABLE floor far short of its own apparent design intent, found during row 1.38 (Cushing's syndrome) scoping -- a real bug in …
- **DONE (2026-08-19), part (a) -- point-of-use clamp extended to all four remaining rule-11 signals.** Found during the CORTISOL_P3_SENSITIVITY fix's own security review: the p3-specific point-of-use clamp closed only ONE instance of a systemic gap -- every …
- **DONE (2026-08-19, third attempt) -- see the full fix summary below this entry's own history.** Found during a second attempt at row 1.38 (Cushing's syndrome) -- a real, pre-existing bug in the ALREADY- SHIPPED `stress_level` feature, not something new work …
- **DONE (2026-10-20, execution plan Phase 4)** -- split the oversized modules without behaviour change: `drug_interactions`, `modifiers` and `parameter_seed` became re-exporting packages; `SimulationRunner`'s pathway assembly and reporting moved to mixins. See `HISTORY.md`.
- **DONE (2026-10-20)** -- the two functions left whole by the Phase 4 split: `SimulationRunner.run()` (~580 lines) is now a short pipeline over stage helpers (pure warning builders moved to `biosim_api/scenario_warnings.py`), and `_single_disease_modifier()` (~645 lines, a commented if-chain) is a dispatch table over one function per `DiseaseState`, each keeping its citation comments. All nine benchmarks identical; see `HISTORY.md`.
- **DONE (2026-10-20)** -- `pip-audit` advisories in `pypdf` 6.14.2, `urllib3` 2.7.0 and `pip` 26.1.1 cleared by upgrading (6.21.1 / 2.8.0 / 26.2.1); floors in `requirements.txt`. See `HISTORY.md`.

## 9. How to pick the next chunk (established process — keep using it)

**See `archive/AUTONOMOUS_LOOP_PROMPT.md`** (archived 2026-11-17) for a self-contained
prompt implementing this exact procedure, designed to run unattended
via `/loop` (dynamic self-pacing) across as many turns/sessions as the
backlog needs, without the project owner needing to re-invoke it each time.

**Search-exhaustion rule** (item D3, 2026-08-17 self-evaluation): a
blocked literature item (see the many "Permanently/genuinely blocked"
rows above -- 1.1, 1.6, 1.14, 1.15, 1.16, 1.19, 1.22 -- 1.3b/glimepiride was resolved later via a fresh-angle ratio search, which is exactly the "specific new lead" this rule exempts)
that has already had 2 independent search passes fail with genuinely
different angles is "parked" -- do not re-select it for a default
work-selection pass. Only pick it back up when a *specific* new lead
exists (a named different trial, a named different search term, a
newly-available full-text source) -- not "try searching again." This
keeps repeated blocked-item searches from crowding out Track 3/5 infra
work that's actually unblocked. Recurring re-checks explicitly exempted
from this rule (1.6, 1.12, 3.12) stay on their own stated cadence.

1. Read `CLAUDE.md`'s STRATEGIC ROADMAP and every recent DONE entry
   first — this file is not a substitute for that.
2. Respect track order (Section 1's table). Within the current track,
   prefer a chunk with no unresolved dependency over one that's blocked.
3. For anything not already fully scoped in this document, dispatch a
   small parallel fan-out (the specialist agents this project already
   has: `computational-biologist`, `literature-curator`,
   `database-architect`, `architect`, `general-purpose`) to verify the
   candidate is still real and still the best use of the next chunk —
   agent output is a hypothesis, verify every citation/code claim
   yourself before trusting it.
4. Execute with the full discipline: implement -> RED-proof via a
   scoped `git stash` -> GREEN -> full suite -> code-reviewer (+
   security-reviewer for drug/dosage/physiology/genetics/calculation
   code) -> if a review agent hits a quota outage, do not retry
   identically, perform the equivalent verification yourself and say so
   honestly -> terse `CLAUDE.md` DONE entry -> update this file's
   status for that chunk -> commit (verify `git status` first, stage
   explicitly by path, never `-A`).
5. When a track's every named chunk here is done or explicitly
   resolved, re-run a fresh scoping fan-out for that track before
   declaring it exhausted a second time — new findings keep surfacing
   (this session alone found five real bugs and one real Track 1 gap
   that four prior "exhausted" declarations had missed).

---

## 10. Process & quality improvements (2026-08-17 self-evaluation)

The project owner asked for a project self-evaluation (goals, quality, direction,
coding) plus a `/verify` pass on the most recent commit. That pass
found a real, previously-unnoticed finding (1.23 above) purely by
driving the real REST API rather than trusting the passing test suite
-- worth institutionalizing the habit, not just fixing the one
instance. Items C1/C2 got real numbered rows in their home tracks
(1.23, 3.14) above; D3 folded into Section 9's process directly. The
remainder lives here since they're either standing decisions for
The project owner, process/protocol amendments, or small cross-cutting chunks that
don't belong to one track.

**Decisions for the project owner -- resolved 2026-08-17 (both, in the same conversation):**

| # | Item | Notes |
|---|---|---|
| G2 | Force the Track 6 (platform/distribution) decision | **DONE.** the project owner: Biosim's audience is a research/personal tool. Track 6 is now formally parked (Section 7), not "proposed, ask again later" -- do not resurface SDK packaging/user-facing docs/a demo-UI layer as filler work. |
| G3 | Name the validation ambition explicitly | **DONE.** the project owner: "internally consistent" (current bar, already met) remains the standard, not "eventually citable/reviewable." Track 3's Docker-gated DB migration and Track 2's remaining statistical work stay in their normal track-order position -- do not pull them forward. |

**Protocol/process amendments (no the project owner decision needed, apply directly):**

| # | Item | Size | Notes |
|---|---|---|---|
| G1 | "Survives composition" acceptance criterion | Done | Folded into `CLAUDE.md`'s standing "What does not change" paragraph (2026-08-17) -- a milestone isn't "done" until demonstrated (a real run, not just reasoning) to reach the output for at least the calibration population, given active safety clamps/modifiers. |
| Q2 | Self-performed-equivalent-verification checklist | Done | Folded into the same `CLAUDE.md` paragraph -- named the 5 checks that were previously reconstructed ad hoc each time a review agent hit a quota outage (double-counting, fixed-point invariance, extreme-value/NaN, safety-clamp interaction, RED-proof via stash). |
| D1 | Composition-Boundary Behavioral Audit protocol | **DONE, first pass run (negative result)** | Written into `CLAUDE.md` as "COMPOSITION-BOUNDARY BEHAVIORAL AUDIT PROTOCOL", same shape as the citation re-audit protocol (5.5) -- trigger (5+ new drug/disease milestones since the last marker), method (drive the real REST server across the mechanism's full … |
| D2 | Population x mechanism coverage matrix | **DONE (first pass)** | See Section 11 below. Computed real crossover thresholds (not guessed) for every currently-known Gb-lowering mechanism against `MIN_GB_THRESHOLD_MARGIN_MM`, cross-referenced against every shipped disease preset's own Gb. … |
| Q5 | Stagger review-agent usage | **Done** | Folded into `CLAUDE.md`'s standing discipline paragraph, right after the Q2 checklist -- security-reviewer dispatch is now scoped to dosing/safety-clamp/disease-state/genetics changes specifically; lower-risk changes rely on Q2's checklist instead. … |
| D4 | Molecular-vs-context identifier boundary | **DONE** | Row 5.1's own chunk 10 note flagged sleep restriction (row 1.26) and CKD (row 1.24) as having no plausible molecular-entity subject for `IDRegistry` -- the project owner's explicit architectural guidance (2026-08-18): don't expand `IDRegistry` to absorb non-molecular … |

**Small chunks (pick up as filler, same as the housekeeping backlog above):**

| # | Item | Size | Notes |
|---|---|---|---|
| Q4 | Stop treating `KNOWN_TEST_COUNT` as mandatory-exact-every-commit | Done | Comment added in `app.py` (2026-08-17) and the `CLAUDE.md` discipline paragraph amended -- the field stays (real health-check value), but a commit is no longer blocked on getting the number exactly right. |
| Q3 | Split `CLAUDE.md` into a lean current-state doc + a separate chronological `HISTORY.md` | **DONE** | `CLAUDE.md`: 14,879 -> 1,955 lines (87% reduction). … |
| C3 | `DrugMechanism` descriptor for the growing drug roster | Medium, YAGNI-gated | 16 drugs now dispatched via a handful of parallel dicts (`_DRUG_FUNCTIONS`, `GLP1_DRUGS`, `SULFONYLUREA_DRUGS`, `ALPHA_GLUCOSIDASE_INHIBITOR_DRUGS`, `EXOGENOUS_INSULIN_DRUGS`, ...). A structured descriptor (mechanism type, known clamp interactions, citation, calibration-population range) would make D2's matrix machine-checkable and give 1.23's warning a natural lookup table. Explicitly deferred until the pattern causes real friction (one more drug class, maybe) -- don't refactor a stable system speculatively. |
| L1 | Reusable citation re-verification tool, full-text-capable | **DONE** | packages/literature-pipeline/src/biosim_literature/verification.py, new verify_citation(pubmed_id, claims, fetcher=None) -- replaces the ad-hoc python -c "... … |
| C4 | Deterministic (non-agent) pre-check script | **DONE** | `scripts/precheck.py`, run via `python scripts/precheck.py`. … |

---

## 11. Population x mechanism coverage matrix (item D2, first pass 2026-08-17)

Real, computed data (not estimated) -- every number below was produced
by direct execution against the actual codebase
(`metformin_modifiers()`, `alcohol_modifiers()`, `tirzepatide_gb_factor()`,
`disease_modifiers()`), not derived by reasoning about the formulas.
Scope of this first pass: only the confirmed joint clamp
(`MIN_GB_THRESHOLD_MARGIN_MM`, `personalizer.py`) against every
currently-known Gb-LOWERING drug mechanism. Gb-RAISING mechanisms
(every disease state's own Gb factor, glucocorticoid's +10%) are listed
for reference but don't interact with this particular clamp (a lower
bound), and the other, independent `_PARAM_BOUNDS` clamps (Ib, p1, p2,
p3, ke, secretion_scale) are NOT yet swept -- flagged as the next D1
audit pass's own scope, not computed speculatively here.

**Gb-lowering drug mechanisms, crossover threshold to clear the clamp**
(`Gb_needed = (G_threshold + 0.3) / factor`, at default `G_threshold=4.5`):

| Mechanism | Dose/condition | Gb factor | Crossover Gb needed | Warned by row 1.23's fix? |
|---|---|---|---|---|
| Metformin | 2000mg (max), chronic | 0.853 | 5.627 mM | Yes -- confirmed via direct execution (see below) |
| Tirzepatide | 5mg | 0.680 | 7.059 mM | Yes (this is what the fix was built for) |
| Tirzepatide | 10mg | 0.619 | 7.754 mM | Yes |
| Tirzepatide | 15mg | 0.596 | 8.054 mM | Yes |

**Shipped disease-state default Gb values** (base 5.0 mM x disease factor):

| Disease state | Gb factor | Effective Gb | Source function (`modifiers/disease.py`) |
|---|---|---|---|
| (healthy default) | 1.000 | 5.000 mM | `_healthy_modifiers()` |
| Insulin-resistant | 1.030 | 5.150 mM | `_insulin_resistant_modifiers()` |
| Gestational diabetes | 1.020 | 5.100 mM | `_gestational_diabetes_modifiers()` |
| Acromegaly | 1.135 | 5.673 mM | `_acromegaly_modifiers()` |
| Type 2 diabetes | 1.300 | 6.500 mM | `_type_2_diabetes_modifiers()` |
| GCK-MODY | 1.434 | 7.170 mM | `_gck_mody_modifiers()` |
| Type 1 diabetes | 2.200 | 11.000 mM | `_type_1_diabetes_modifiers()` |

(Function column and the acromegaly row added 2026-10-20 when
`_single_disease_modifier()` became a per-state dispatch table; values
re-computed via `disease_modifiers()`. Acromegaly (row 1.37) postdates this
matrix's first pass; its own composition check against metformin 2000mg
(raw ~4.84 mM, just above the 4.8 mM floor) is in that function's comment.)

**Cross-referencing the two tables, real findings**:

- **Metformin at max dose (2000mg) engages the clamp for the healthy
  default, insulin-resistant, AND gestational-diabetes populations**
  (all below the 5.627 mM crossover) -- only clears it for T2D (6.500,
  metformin's own actual calibration population, per Hansen et al.
  2024) and GCK-MODY/T1D (both well above it, though neither is a
  realistic metformin population). **Independently confirmed via
  direct execution** (not just computed) that row 1.23's
  `clamp_warnings` fix -- built and tested against tirzepatide -- is
  genuinely mechanism-agnostic and already fires correctly for all
  three flattened metformin cases, and correctly stays silent for T2D.
  No new code needed; locked in with a dedicated regression test
  (`test_clamp_warning_is_genuinely_mechanism_agnostic_not_tirzepatide_specific`).
- **Tirzepatide's flattening is universal across every currently
  clinically-realistic population**: healthy (5.0), insulin-resistant
  (5.15), and gestational diabetes (5.1) are all far below even the
  5mg crossover (7.059); T2D (6.5, tirzepatide's actual real-world
  population) is ALSO below it. Only GCK-MODY (7.17) clears the 5mg
  threshold specifically, and still falls short of 10mg/15mg's own
  higher crossover points -- but GCK-MODY is not tirzepatide's
  clinical population anyway. In practice, for every population this
  drug would actually be prescribed to, the dose-response Furihata et
  al. 2022 measured is not observable in this platform's simulated
  output at all -- it is fully absorbed by the clamp, not merely
  partially.
- **REMOVED (2026-08-18, mandatory code-reviewer finding on a fresh
  Track 1 alcohol-mechanism fix)**: the row above previously listed
  alcohol as a Gb-lowering mechanism (max ~5% reduction, barely
  clamping for the healthy default). This is now FALSE -- a scoping
  fan-out found alcohol's own Gb-reduction mechanism was contradicted
  by both available primary sources (Siler 1998 PMID:9815011, Shelmet
  1988 PMID:3280601, both live-verified: "no effect on plasma glucose
  concentrations"/"no effect on glucose tolerance"), and it was fixed
  to a transient, ODE-coupled `hepatic_egp_delta` contribution
  (`biosim_ode.alcohol_hepatic_pathway.AlcoholHepaticPathway`) instead
  -- `alcohol_modifiers()` now returns `{}` unconditionally and does
  not touch `Gb`, this clamp, or this static personalization pipeline
  at all (confirmed directly: an alcohol `DrugDose` added alongside
  metformin now produces a byte-identical `Gb` to metformin alone).
  `test_gb_threshold_margin_never_collapses_below_safety_bound`/
  `test_gb_threshold_margin_clamp_bounds_the_insulin_spike` were
  corrected in the same pass to drop the now-inert alcohol dose. See
  `HISTORY.md` "DONE — Fixed alcohol's glucose mechanism" for the full
  fix detail.

**Not yet done (explicit scope boundary for this first pass)**: sweep
the independent `_PARAM_BOUNDS` clamps (`Ib`, `p1`, `p2`, `p3`, `ke`,
`secretion_scale`) the same way -- none of these have a joint-clamp
warning mechanism at all today (only the Gb/G_threshold pair does), so
a drug whose own dose-response pushes one of these to its bound would
currently do so completely silently, with no equivalent to 1.23's fix
protecting it. This is real, unswept surface -- the natural target for
the next D1 Composition-Boundary Behavioral Audit pass.
