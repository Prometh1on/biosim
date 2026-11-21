# CLAUDE.md — Biosim Project Handoff
# Lives at the repo root.
# Claude Code reads this automatically at session start.

---


## WHAT THIS PROJECT IS

**Biosim** is a computational systems biology platform that simulates
postprandial glucose/insulin dynamics using a system of ODEs.

Users input food, genetic variants, and physiological parameters.
The platform outputs quantitative (expert) and plain-language (lay) predictions.

**The core scientific claim:** The same meal produces a measurably different
glucose response depending on time of day, prior exercise, genetic variants,
medications, and disease state. The platform quantifies this.

**Current clinical relevance:** Models GLP-1 agonists (semaglutide/Ozempic),
metformin, steroid-induced hyperglycaemia, circadian metabolic variation,
T2D genetic risk, and post-exercise insulin sensitisation.

---

## STRATEGIC ROADMAP (added 2026-07-28) — read this before picking a milestone

**See also `MASTER_PLAN.md` (added 2026-08-05)** — a full decomposition
of Tracks 1-4 below plus two proposed additional tracks (5: knowledge
infrastructure, 6: platform/distribution) into session-sized chunks,
with dependencies and sizing per chunk. This section stays the
authoritative statement of the four approved tracks and their
direction-setting decisions; `MASTER_PLAN.md` is the forward-looking
backlog derived from it, kept in sync separately (check its own "last
synced against" line before trusting it's current).

The project owner gave a standing, expanded mandate: grow Biosim continuously toward a
scientifically rigorous, eventually multi-scale/enterprise-grade platform,
not just fix isolated bugs forever. Four direction-setting decisions were
made explicitly (do not re-litigate these without asking again) before any
work under this mandate began:

1. **Track order**: start on **Track 1 (deepen biological realism)** now.
   Tracks 2-4 below are queued *after* Track 1, not simultaneous -- pick
   the next milestone from Track 1 unless the project owner says otherwise or Track 1
   is genuinely exhausted.
2. **Cost boundary**: stay **fully free/local**. Design cloud-deployment
   code if it comes up (Dockerfiles, IaC, CI configs) but never actually
   provision or pay for anything -- matches the pre-existing "no paid API
   keys" rule this project has followed all along (e.g.
   `test_ingestion_e2e.py` stays gated/skipped for the same reason).
3. **Data layer**: when Track 3 (architecture/infra) starts, go all the
   way to a **full client-server DB** (PostgreSQL for parameters/
   provenance, a graph DB for entity cross-references/pathway
   relationships) -- not the lower-risk embedded-SQLite middle step. Must
   still be **local-only** (docker-compose, no managed cloud DB) per
   decision 2, and migrated incrementally (dual-write/shadow-read, proven
   equivalent against the existing 1355+ tests and 5 clinical benchmarks
   before any cutover) -- see the `database-architect` agent below.
4. **Subagents**: authored three new custom agent types in
   `~/.claude/agents/` (not built into the base harness roster) rather
   than approximating with generic ones:
   - **`computational-biologist`** -- mechanistic/literature investigation
     for pathway and parameter questions. Read-only by design (no
     Write/Edit) -- it proposes a hypothesis, the orchestrating session
     still does the actual verify -> implement -> RED/GREEN via
     `git stash` -> full suite -> `CLAUDE.md` DONE entry -> commit
     sequence unchanged. This is deliberate: keeping it read-only
     preserves the "agent output is a hypothesis, never truth" rule this
     project runs on.
   - **`literature-curator`** -- PubMed/citation verification specialist
     (fetch via `PubMedFetcher`/live NCBI E-utilities, never assert a
     citation from memory). Also read-only.
   - **`database-architect`** -- schema/migration design for Track 3's
     parameter-db/id-resolver -> Postgres+graph-DB move. Has Write/Edit
     (needs to actually produce schema files/docker-compose/migration
     scripts), explicitly scoped to local-only zero-cost infrastructure
     and incremental, test-verified migration -- not a cutover agent.
   Use these via the `Agent` tool the same way as the existing roster
   (`architect`, `planner`, `security-reviewer`, etc.); all three still
   have the "Prompt Defense Baseline" boilerplate every other agent in
   `~/.claude/agents/` carries.

### The four tracks

**Track 1 -- Deepen biological realism (current focus).** Same proven
cadence as this session's 40+ literature-anchored fixes: new organ/
pathway modules within the existing `CompositePathway`/9-signal coupling
contract (rule 11) -- candidates not yet built: renal glucose
reabsorption/SGLT2 inhibitors (a natural next drug class alongside
GLP-1/metformin/steroids), thyroid-glucose crosstalk as a real ODE-coupled
pathway (currently only a static modifier), adipose/leptin signalling.
More PK/PD (time-varying drug concentration curves, not just static
dose-response snapshots -- caffeine/alcohol already have real
pharmacokinetics, metformin/glucocorticoids don't). More clinical
validation benchmarks (liraglutide is blocked on paywalled full text; the
GI-timing-slope question and the GLP-1 gastric-emptying/secretion-
interaction mechanism are both open literature questions, see VALIDATION
STATUS -- worth periodic re-checking, not continuous re-searching). Zero
new infrastructure, zero new cost, ships value every milestone.

**Track 2 -- Statistical rigor.** Extend `biosim_uncertainty` further:
real Bayesian parameter estimation against the clinical benchmarks
already built (rather than point-estimate + tolerance-band validation
only), literature-anchored virtual population generation (extending
`biosim_seeds.population`'s existing cohort simulation with real
demographic/genetic distribution sampling), sensitivity analysis at
scale. Builds on existing packages, no new infra, no new cost.

**Track 3 -- Architecture/infra modernization.** The full client-server
DB migration (decision 3 above, `database-architect` agent), local
multiprocessing/distributed batch execution (the existing
`BatchRunner`/`PopulationRunner` are single-process -- start with local
multi-core before anything cloud-based, per decision 2), a properly
packaged/versioned SDK once the REST API surface stabilizes.

**Track 4 -- Substrate/tissue-resolved metabolism** (renamed from
"Genome-scale metabolism (hybrid ODE/FBA)" -- reframing adopted by
The project owner 2026-08-10, see MASTER_PLAN.md Track 4 and the "DONE -- Adopted
the Track 4 reframing" entry in HISTORY.md). A different modeling paradigm
entirely (constraint-based, needs a real metabolic network
reconstruction + an LP solver) is one *candidate implementation*, not
the track's definition -- genome-scale FBA specifically requires a
named, concentration-independent biological question the ODE paradigm
can't answer, a flux-type validation target, and a free/local compute
budget (MASTER_PLAN's own adopted entry criteria) before any FBA code
is written. This remains its own dedicated architecture-planning
initiative, not a drop-in milestone -- don't start writing constraint-
based code without a separate scoping pass first (use the `architect`
agent, and likely `computational-biologist` for the biological side).

### What does not change

Every fix under any track still follows the project's existing
discipline unchanged: verify every claim yourself (agent/tool output is a
hypothesis, never truth), never fabricate a PMID/parameter/benchmark
value, TDD RED-via-`git stash`-then-GREEN, run the full suite, mandatory
code-reviewer pass (+ security-reviewer for anything touching drugs/
dosage/physiology/genetics/calculations affecting biological outputs --
if a review agent hits a session-quota outage, don't retry identically,
do the equivalent verification yourself and say so honestly), append a
terse "### DONE --" entry to **`HISTORY.md`** (not `CLAUDE.md` --
changed 2026-08-17, item Q3: `CLAUDE.md` had grown to ~14,900 lines
holding this entire chronological log, and was split into a lean
current-state doc plus `HISTORY.md`'s own append-only record; see
`HISTORY.md`'s own header for the full split rationale) and the new
test count (best-effort -- a commit isn't blocked on getting this
number exactly right; see MASTER_PLAN.md's process-improvement
section, item Q4), commit (never push without being explicitly asked
-- no git remote is configured anyway).

**"Done" additionally requires surviving composition, not just its own
calibration point** (added 2026-08-17, see MASTER_PLAN's process-
improvement section, item G1). A new drug/disease mechanism's own
isolated unit test passing is not sufficient — before marking a
milestone done, demonstrate (a real scenario run, not just reasoning
about the code) that the calibrated clinical signal actually reaches
the simulated output for at least the population the calibration data
came from, given whatever safety clamps/other active modifiers are
already in the composition. Prompted by the tirzepatide milestone: the
calibration was independently verified correct, and the effect was
still invisible in the platform's own default T2D preset because
nothing checked whether `MIN_GB_THRESHOLD_MARGIN_MM` would swallow it
first.

**When a review agent is quota-exhausted, the self-performed
equivalent verification is this fixed checklist** (named 2026-08-17,
item Q2 -- previously reconstructed ad hoc each time, now written down
so a step isn't skipped under time pressure):
1. Double-counting: does this mechanism's effect overlap with any
   already-active static modifier or coupled pathway for the same
   drug/disease state?
2. Fixed-point invariance: `validate_y0(tolerance=1e-8)` (or the
   equivalent basal-state check) still returns clean at the new
   default and at least one disease-elevated preset.
3. Extreme-value/NaN behavior: the new function doesn't blow up,
   return NaN, or silently clamp to a wrong value at 0, negative, and
   very large inputs.
4. Safety-clamp interaction: if this mechanism can move `Gb`/
   `G_threshold`/any other clamped field, check it against the real
   clamp margin at a realistic value, not just algebraically (see the
   "survives composition" criterion above -- this is the check that
   would have caught the clamp-flattening finding early).
5. RED-proof every new regression test via a scoped `git stash` on
   just the changed source file(s), confirming it fails for the
   documented reason before restoring.

**Stagger review-agent usage, don't dispatch both on every milestone
by default** (item Q5, 2026-08-17): the weekly review-agent quota has
been hit multiple times across single sessions in this project's
history. Dispatch the security-reviewer pass specifically for changes
touching dosing/safety-clamp/disease-state/genetics logic (where a
missed finding has real consequence); for lower-risk changes (docs,
standalone reference data, test-only additions, pure refactors with no
behavior change), rely on the Q2 checklist above instead of spending
an agent dispatch on it. The code-reviewer pass stays mandatory for
every code change regardless -- this only narrows when the *second*
agent is worth the quota cost.

## REPOSITORY STRUCTURE

```
Biosim/                           <- repo root, branch master, no remote
├── CLAUDE.md                     <- this file (current state + rules)
├── MASTER_PLAN.md                <- forward-looking backlog
├── HISTORY.md                    <- append-only DONE log
├── conftest.py                   <- adds every packages/*/src to sys.path
├── packages/<name>/src/<module>/ <- 17 packages, see PACKAGES below
├── infra/                        <- local docker-compose DB migration (Track 3)
├── archive/                      <- old docs kept for history (e.g. the autonomous-loop prompt)
└── scripts/                      <- precheck.py, demo.py (README one-command demo),
                                     plot_heldout.py (validation figure)
```

**Test count:** don't trust a number written here -- the latest `HISTORY.md` DONE
entry records the count at that commit, and `python -m pytest -q` gives the truth.

**Setup:**
```powershell
cd biosim
.venv\Scripts\activate
python -m pytest --tb=short -q
```

---

## ARCHITECTURE — NON-NEGOTIABLE RULES

Violating these will break tests silently in ways that are hard to debug.

### 1. ODE IS PRIMARY. BOOLEAN IS WATCHDOG ONLY.
The Boolean/hybrid model runs in parallel to check consistency.
It NEVER modifies or blends with ODE outputs.
Divergence between ODE and Boolean → widen confidence intervals only.

### 2. DELTA-SECRETION FORMULATION FOR INSULIN ODE
```python
delta_secretion = F_secretion_basal * (current_frac - basal_frac) * secretion_scale
dI = -ke * (I - Ib) + delta_secretion
```
This resolves the steady-state drift bug. Do NOT use absolute secretion.

### 3. np.trapezoid NOT np.trapz
Python 3.14 / NumPy 2.0 removed `np.trapz`. Use `np.trapezoid` everywhere.

### 4. LOG-NORMAL DISTRIBUTIONS FOR KINETIC PARAMETERS
All kinetic parameters are always-positive and multiplicatively uncertain.
Use `LogNormalDistribution.from_confidence_tier()` not NormalDistribution.

### 5. APPEND-ONLY PARAMETER DATABASE
The ParameterStore is append-only. Never delete or overwrite parameters.
Conflict detection at 2× fold threshold. GOLD/SILVER/BRONZE confidence tiers.

### 6. REACTOME IDs FOR PATHWAYS ONLY
Proteins use UniProt. Metabolites use ChEBI/HMDB. KEGG for cross-ref only.
Never assign a Reactome pathway ID to a protein entity.

### 7. PRIVACY: NO rsIDs IN OUTPUTS
PersonalizationResult must never contain raw rsIDs or allele data.
`PrivacyAuditor().audit(output)` must pass on every SimulationOutput.
Genetic personalisation exposes only parameter deltas, never raw variants.

### 8. secretion_scale FIELD IN PathwayParameters
The `secretion_scale` field (neutral default = 1.0) carries genetic and
physiological modifiers to insulin secretion. Never hardcode secretion.

### 9. FIXED GR REFERENCE FOR CORTISOL p3_modifier
```python
GR_REFERENCE = 0.15  # population-average, NOT instance.gr_occ_basal
excess_gr = max(0.0, gr_occ - GR_REFERENCE)
```

### 10. KINETICS ARGUMENT ORDER
`michaelis_menten(substrate, vmax, km)` — vmax before km.
`hill_activation(x, ec50, n)` — must raise ValueError for ec50<=0, n<=0, x<0.

### 11. EXTENSION-PATHWAY COUPLING CONTRACT
`GlucoseInsulinPathway.dydt()` reads thirteen optional `signals` keys, ALL with
neutral defaults (`signals.get(key, neutral)`), so it behaves identically
standalone or composed. (Corrected 2026-08-18: this list had drifted to only
nine keys, missing `renal_glucose_clearance_delta`/`exogenous_insulin_delta`
that were already implemented — the exact "contract drifted from
implementation" bug class `pathway.py`'s own descriptor comment warns
about, caught again while adding `nutrient_secretion_delta`.)
```
p3_modifier                  mult. on p3   (default 1.0)
gb_modifier                  mult. on Gb   (default 1.0)
insulin_resistance_modifier  mult. on p3   (default 1.0)
hepatic_egp_delta            add. to dG    (default 0.0, mM/min — EXCESS over
                              the source pathway's own basal EGP, never the
                              raw absolute value, same reasoning as rule 2)
incretin_effect               scales delta_secretion (default 0.0)
gastric_emptying_modifier      mult. on Ra   (default 1.0)
p1_modifier                    mult. on p1   (default 1.0)
extra_g_clearance_modifier      mult. on extra_g_clearance (default 1.0)
secretion_scale_modifier        mult. on secretion_scale (default 1.0)
renal_glucose_clearance_delta   add. to dG    (default 0.0, mM/min — insulin-
                                 independent SGLT2i renal glucosuria)
exogenous_insulin_delta          add. to dI    (default 0.0, pM/min — injected
                                 insulin from a subcutaneous depot)
nutrient_secretion_delta          add. to dI    (default 0.0, pM/min —
                                 protein-co-ingestion-driven endogenous
                                 insulin secretion)
exercise_glucose_uptake_delta     add. to the glucose clearance rate
                                 (default 0.0, min^-1 — working-muscle uptake
                                 during a bout; added 2026-11-16, B7)
```
Since B7 (2026-11-16) `PathwayParameters.extra_g_clearance` defaults to 0: the old
unsourced 0.018/min double-counted p1 (S_G already stands for all insulin-independent
disposal). Do not put a resting clearance term back on top of p1.
Any new coupled extension MUST use one of these keys or add a new one here
(never bypass this contract by mutating PathwayParameters in place — rule
"Immutability" in coding-style.md). Every new/changed key must keep dG=0 and
dI=0 at the analytic basal state (`GlucoseInsulinPathway().y0`) when the
extension itself is at its own basal — verify with a regression test before
merging, the way `test_basal_steady_state_preserved_with_signals_present` and
`test_composite_preserves_basal_steady_state` do.

**Important, easy to miss**: `CompositePathway.dydt()` calls every coupled
pathway's `output_signals(t, y)` with the CURRENT integration time `t` on
*every* solver step, not just `t=0` — so a coupled pathway's signals CAN be
genuinely time-varying throughout a simulation (cortisol's GR occupancy,
hepatic's EGP, NEFA, and GLP1 all already rely on this). This was
mis-stated as "doesn't exist" in an earlier version of this document (see
KNOWN BUGS FIXED #16's `fasted_run()` note) before `p1_modifier`/
`extra_g_clearance_modifier` were added specifically to make time-gated
exercise couplings (during-bout vs. post-bout) possible.

---

## PACKAGES AND WHAT THEY DO

| Package | Module | Purpose |
|---------|--------|---------|
| id-resolver | biosim_id_resolver | Maps biological entities (proteins, metabolites, pathways) to canonical IDs with cross-references (UniProt, ChEBI, HMDB, Reactome, KEGG) |
| parameter-db | biosim_parameter_db | Append-only kinetic parameter store with confidence tiers and provenance tracking |
| literature-pipeline | biosim_literature | Fetches and normalises kinetic parameters from PubMed abstracts |
| curation | biosim_curation | Priority queue for expert review of model parameters |
| ode-engine | biosim_ode | Core ODE system: Bergman minimal model (G/I/X) + three insulin-driven read-outs (IGF-1/mTOR/AMPK) + pathway extensions |
| uncertainty | biosim_uncertainty | Monte Carlo UQ, OAT sensitivity analysis, confidence bands |
| genotype | biosim_genotype | GWAS variant → parameter modifier pipeline (privacy-preserving) |
| physiology | biosim_physiology | Age/sex/BMI/fitness/disease/circadian/drug/exercise modifiers |
| simulation-api | biosim_api | Scenario definition → full simulation pipeline → SimulationOutput |
| expert-tools | biosim_expert | SBML export, batch runner, CSV export, molecular timeline, sensitivity reporter |
| lay-output | biosim_lay | Food database, traffic-light risk profile, plain-English reports, UQ lay report |
| validation | biosim_validation | Clinical benchmark validation against published OGTT reference values |
| infrastructure | biosim_infra | LRU simulation cache, privacy auditor, reproducibility bundles |
| data-seeds | biosim_seeds | Real literature parameter seeds, PlatformSession, population simulation, chrono/drug/GI/exercise scenarios |
| rest-api | biosim_rest | FastAPI HTTP layer over biosim_api: POST /simulate, GET /scenarios, POST /compare, GET /health |
| metabolism | biosim_metabolism | Track 4: stoichiometric hepatic network, FVA (HiGHS), flux benchmarks, EGP sensitivity table (offline, not wired into the ODE) |
| context-registry | biosim_context_registry | Registry for non-molecular context entities (disease states, physiological axes), separate from IDRegistry |

---

## ODE MODEL — 3 COUPLED STATES + 3 INSULIN-DRIVEN READ-OUTS

G, I and X form the coupled Bergman minimal model. IGF1, mTOR and AMPK are
insulin-driven read-outs: they read I, but nothing in dG/dI/dX reads them, so
they never feed back into glucose or insulin. Say "three insulin-driven
read-outs", not "six coupled ODEs".

```
State:  [G_mM, I_pM, X_min-1, IGF1_nM, mTOR_01, AMPK_01]
Units:  internal = mM, pM, min | canonical output = uM, s
Basal:  [5.0,   66.0,  0.0,    13.0,    0.15,    0.60]
```

**Key parameters (PathwayParameters):**
- `p1` = glucose effectiveness (0.0317 min⁻¹, Bergman 1981, confirmed against
  Pacini et al. 1998 PMID:9575818's directly-comparable normal-subject range)
- `p2` = X-compartment clearance (0.0123 min⁻¹, Pacini 1986)
- `p3` = insulin action (9.032e-7 min⁻² per pM -- corrected this session from
  an uncited 5.03e-6 to match Bergman 1981's own directly-quoted SI, Fig 6
  caption, PMID:7033284 -- see VALIDATION STATUS below)
- `ke` = insulin elimination (0.093 min⁻¹, Sherwin 1974)
- `KG_secretion` = beta-cell glucose-sensing Michaelis constant (5.63 mM,
  Rebelato et al. 2018 PMID:30166558 -- corrected this session from an
  uncited 2.0)
- `first_phase_gain` = dynamic/first-phase secretion gain (121.3 pM per
  mM/min -- fitted this session; functional form from Mari et al. 2002
  PMID:12388151)
- `secretion_scale` = genetic/physiological modifier (default 1.0)

Both pre-registered clinical validation benchmarks (`ModelValidator
.validate_ogtt()`, `.validate_bergman_iv()`) currently PASS -- see
VALIDATION STATUS below for full numbers.

---

## PATHWAY EXTENSIONS (in biosim_ode/)

| File | Class | Maturity | Couples to | Activated by |
|------|-------|----------|------------|---------------|
| cortisol_pathway.py | CortisolPathway | BETA | GlucoseInsulinPathway via p3_modifier, gb_modifier | `ScenarioDefinition.stress_level > 0` |
| hepatic.py | HepGlucPathway | BETA | GlucoseInsulinPathway via hepatic_egp_delta | `DiseaseState.TYPE_2_DIABETES` and/or `Drug.METFORMIN`, or `DiseaseState.HYPERTHYROIDISM` alone (no combined preset yet -- see `_detect_hepatic()`'s own docstring for the documented gap) |
| nefa_pathway.py | NEFAPathway | EXPERIMENTAL | GlucoseInsulinPathway via insulin_resistance_modifier | `ScenarioDefinition.model_nefa_dynamics = True` (opt-in) |
| glp1_pathway.py | GLP1Pathway | BETA | GlucoseInsulinPathway via incretin_effect, gastric_emptying_modifier | `Drug.SEMAGLUTIDE` / `LIRAGLUTIDE` / `DPP4_INHIBITOR` |
| exercise_pathway.py | ExercisePathway | EXPERIMENTAL | GlucoseInsulinPathway via p1_modifier, p3_modifier, gb_modifier, exercise_glucose_uptake_delta, secretion_scale_modifier | `ScenarioDefinition.exercise_bout_in_progress` (opt-in; distinct from `PersonalizationSpec.exercise_bout`, see KNOWN BUGS FIXED #16) |
| renal_pathway.py | RenalPathway | EXPERIMENTAL | GlucoseInsulinPathway via renal_glucose_clearance_delta (reads glucose_mM as an input signal -- the first extension whose signal is a function of GI's own state, not just time/its own dynamics) | `Drug.DAPAGLIFLOZIN` / `Drug.EMPAGLIFLOZIN` |

Four more extensions exist beyond the table, also composed by
`SimulationRunner._build_pathway()`: `MealAbsorptionPathway`
(`meal_absorption_pathway.py`, mass-conserving gastric emptying, composed only
with an active GLP-1 drug), `ExogenousInsulinPathway` (`exogenous_insulin_delta`,
insulin dose events), `NutrientSecretionPathway` (`nutrient_secretion_delta`,
meals with `protein_g > 0`) and `AlcoholHepaticPathway` (`hepatic_egp_delta`,
alcohol drug events).

All of them are wired into `SimulationRunner._build_pathway()`, which
flat-composes GI with every active extension into one `CompositePathway`
(nesting composites would silently break cross-extension signal passing,
since `CouplingSpec` injection only resolves direct siblings within a
single `CompositePathway` -- see `_build_pathway()` docstring). Each module
exposes `<name>_couplings(gi_name, ext)` (pure `CouplingSpec` list, used for
multi-extension composition) and `couple_with_gi(gi, ext)` (2-pathway
convenience wrapper, used in tests/standalone use).

Deliberate non-overlap decisions (do not "fix" these without re-reading the
reasoning -- each was chosen to avoid double-counting against the existing
static personalisation modifiers in `physiology/modifiers/` and
`physiology/chronobiology.py`):
  - **Cortisol** is NOT used for routine time-of-day meal scenarios --
    `chronobiology.circadian_modifiers()` stays authoritative for that.
    Cortisol only activates via the new `stress_level` field (acute/chronic
    psychological stress), a capability the existing system has no analogue
    for at all.
  - **Hepatic** coupling is ADDITIVE to (not a replacement for) the static
    Gb/p3/secretion_scale disease-state and metformin modifiers: the delta
    formulation makes `hepatic_egp_delta` exactly zero at each preset's own
    basal, so it contributes only the postprandial *kinetics* (impaired EGP
    suppression for T2D, enhanced suppression for metformin) that the static
    modifiers cannot represent -- verified by
    `test_t2d_hepatic_coupling_impairs_postmeal_suppression`.
  - **NEFA** is intentionally NOT triggered by BMI/obesity -- unlike
    hepatic's delta formulation, `NEFAPathway.for_obesity()`'s modifier is
    NOT neutral at its own elevated basal (by design, matching the fixed
    `NEFA_IR_REFERENCE`), so composing it for every obese profile would
    double-count `bmi_modifiers()`. It's opt-in via
    `model_nefa_dynamics=True` for scenarios where its own kinetics add real
    value (extended fasting, post-exercise NEFA mobilisation) using the
    default (basal-neutral) parameters, not the obesity preset.

`GlucoseInsulinPathway.dydt()` consumes every coupling signal key (see
rule 11 above); `ODESystem` accepts any `BasePathway`, bare or composite.

---

## PERSONALISATION LAYERS (applied in order)

```
1. Age modifiers          (physiology/modifiers/demographics.py)
2. Biological sex         (physiology/modifiers/demographics.py)
3. BMI/adiposity          (physiology/modifiers/demographics.py)
4. Fitness level          (physiology/modifiers/demographics.py)
5. Disease state          (physiology/modifiers/disease.py)
6. Circadian (time of day)(physiology/chronobiology.py)
7. Drug interactions      (physiology/drug_interactions/, one module per drug class)
8. Exercise               (physiology/exercise.py)
9. Genetic variants       (genotype/personalizer.py)
10. Sleep duration        (physiology/modifiers/lifestyle.py, sleep_modifiers())
11. Intermittent hypoxia  (physiology/modifiers/lifestyle.py, intermittent_hypoxia_modifiers())
12. Menstrual cycle phase (physiology/modifiers/lifestyle.py, menstrual_cycle_modifiers())
```

`modifiers`, `drug_interactions` and `biosim_seeds.parameter_seed` are packages
(split 2026-10-20) whose `__init__` re-exports every name, so the old import paths
still work; see each `__init__` docstring for the module map. Pathway assembly and
reporting live in `biosim_api/pathway_assembly.py` and `reporting.py`, mixed into
`SimulationRunner`.

All compose multiplicatively into PathwayParameters overrides. Layer 8
(exercise) was listed here but not actually wired into
`PhysiologicalPersonalizer`/`ScenarioDefinition` until this session --
see KNOWN BUGS FIXED #16. All listed layers are genuinely applied.

---

## GENOTYPE — SEEDED VARIANTS

**Default 7 (VariantEffectDB.default()):**
rs7903146 (TCF7L2), rs1801282 (PPARG), rs5219 (KCNJ11),
rs13266634 (SLC30A8), rs10830963 (MTNR1B), rs4607517 (GCK), rs35767 (IGF1)

rs35767 (IGF1) is **UNVERIFIED** (`pubmed_id="UNVERIFIED"`, bronze): no paper
reports its +0.5 nM effect size, and published evidence on even the direction
conflicts -- see its comment in `effect_db.py`. Any run that includes it carries
`UNVERIFIED_VARIANT_NOTES` in its warnings (runner, REST, lay report).

**Extended: 14 distinct variants, 15 effect records (extended_t2d_variants()):**
Above + rs2237892 (KCNQ1/EAS), rs10811661 (CDKN2A), rs7756992 (CDKAL1),
rs4402960 (IGF2BP2), rs11708067 (ADCY5), rs1111875 (HHEX),
rs10010131 (WFS1/protective). IGF2BP2 rs4402960 has two records (one on
`secretion_scale`, one on `IGF1b`), which is why there are 15 records but
14 variants.

Effect sizes are per-allele parameter changes chosen by this project from
each variant's reported biology; GWAS papers report disease odds ratios, not
these parameter values.

---

## DATA SEEDS — SEEDED ENTITIES

**IDRegistry:** 12 proteins (insulin, INSR, IGF-1, IGF1R, mTOR, AMPK, AKT1,
IRS-1, GCK, GLUT4, GLUT2, GH), 8 metabolites (glucose, pyruvate, lactate,
ATP, AMP, glycogen, G6P, acetyl-CoA), 4 pathways (glycolysis, insulin
signalling, mTORC1, IGF-1 signalling).

**ParameterStore:** Bergman 1981, Pacini 1986, Sherwin 1974,
Sarbassov 2005, Hardie 2012, Le Roith 1997.

---

## PLATFORM SESSION (single entry point)

```python
from biosim_seeds import PlatformSession
session = PlatformSession()
output  = session.simulate("high_sugar_meal")
report  = session.report(output)
print(report.glucose_section)
```

Full scenario names: high_sugar_meal, standard_meal, light_snack, fasted,
ogtt_75g, caffeine_and_meal. Condition/drug scenarios (biosim_seeds):
semaglutide_scenario, liraglutide_scenario, dpp4_inhibitor_scenario,
stress_meal_scenario, type2_diabetes_meal_scenario, plus the existing
caffeine/metformin/steroid/alcohol scenarios in drug_scenarios.py.

---

## KNOWN BUGS FIXED — DO NOT REINTRODUCE

One line each; the full write-ups (evidence, reasoning, tests) are in
`HISTORY.md` under "Archived from CLAUDE.md (2026-10-20)".

1. **Root file shadowing** -- never leave `types.py`, `queue.py`, `__init__.py`
   or any package module in the repo root; they shadow the stdlib.
2. **Kinetics argument order** -- `michaelis_menten(substrate, vmax, km)`.
3. **NEFA antilipolysis** uses `hill_activation`: at zero insulin NEFA must RISE.
4. **Cortisol GR reference** is a fixed module constant, never
   `instance.gr_occ_basal` (see #10 for how it is derived).
5. **`np.trapezoid`**, never `np.trapz` (removed in NumPy 2).
6. **Delta secretion** -- `dI` uses `(current_frac - basal_frac)`, not absolute secretion.
7. **Composite coupling was architecturally inert** -- `GlucoseInsulinPathway.dydt()`
   must consume every rule-11 signal; `ODESystem` reads species by name, so it runs
   any `BasePathway`, bare or composite.
8. **NEFA insulin-resistance modifier** uses NEFA in excess of `NEFA_IR_REFERENCE`,
   so it is exactly neutral at basal.
9. **~35 fabricated PMIDs** were found and corrected. Verify every citation live
   (`PubMedFetcher.fetch_by_pubmed_id()` or NCBI) before committing; never write a
   PMID from memory.
10. **Never hardcode the basal value of a Hill-relaxation state** -- derive it from
    the same formula the ODE uses (`GR_OCC_BASAL`), and integrate forward to check:
    `validate_y0()` at t=0 alone cannot catch post-t=0 drift.
11. **Hepatic `EGP_SUPPRESSION_IC50`** is derived by `_solve_egp_suppression_ic50()`
    from the anchored `EGP_BASAL_MM_MIN`/`EGP_MAX_MM_MIN`, not hardcoded.
12. **NEFA fixed point** needed a normalization (antilipolysis relative to its value
    at reference insulin), not a constant swap. Before "fixing" a fixed-point bug,
    check whether the formula's algebra even permits an exact fixed point.
13. **Citation audits grep the whole repo** -- a per-file checklist missed three
    fabricated PMIDs in `benchmark.py`.
14. **T2D EGP preset** -- Beck-Nielsen 2002 (PMID:12643178): basal EGP +12%, IC50
    derived, no maximum-velocity defect (`t2d_egp_scale` stays 1.0).
    `for_metformin()` is intentionally NOT a fixed point (that drift is drug onset).
15. **Computed-but-discarded values** -- `gi_meal_scenario()` computed `kabs` and
    dropped it; now carried by `MealEvent.kabs`. After finding one instance of a
    bug class, sweep the codebase for the same class.
16. **Exercise layer was never wired** -- now applied via
    `PersonalizationSpec.exercise_bout` (personalizer) and the time-gated
    `ExercisePathway` (`ScenarioDefinition.exercise_bout_in_progress`).
17. **IVGTT glucose volume** is `BERGMAN_IVGTT_VG_L_PER_KG = 0.230`
    (Mora-Gonzalez 2024, PMID:38389141), used only by `bergman_ivgtt_y0()`; the oral
    meal `Vg` is deliberately unchanged. The return-time gap was -45%; B7 (2026-11-16)
    removed a double-counted resting clearance (`extra_g_clearance`), leaving -20%. Do not
    tune `p1`/`p2` against it: the 90-min reference is itself an unverified figure readout.
18. **`mtor_activation` is not a real benchmark** -- `pubmed_id="UNVERIFIED"` makes
    `is_scored` False and reports show `[PASS/INFO]`; don't invent a reference value.

Later fixes (2026-08 onward) are recorded only in `HISTORY.md`.

---

## VALIDATION STATUS

Run `ModelValidator().validate_all()` for the current verdicts. Since 2026-11-12 a
cited point outside its band gives the verdict **MISS** (it used to be hidden behind
PASS whenever the band was WIDE, +/-50%). Current verdicts (after B7, 2026-11-16, which
removed a double-counted resting glucose clearance; nothing re-tuned; full before/after in
docs/VALIDATION.md section 6): **FAIL** -- `validate_ogtt` (peak glucose +33%, its only scored
NORMAL point since B4 made the wrongly-cited 2-hour point unscored; unscored points no longer
count in the verdict at all; peak insulin
+37%: its absorption and first-phase constants were fitted with the double count in place;
mTOR point unscored). **PASS** -- `validate_bergman_iv` (-20%), `validate_semaglutide`
(+0.7%, one point; VALIDATED since 2026-11-21 by the user's decision, with a caveat: modelling the trial's obese population gives -26% to -39%), `validate_metformin` (-31%, WIDE band).
**MISS** -- prednisone/T2D -68%, SGLT2i -85%/-81%,
`validate_protein_coingestion` (1/4 within band), `validate_t2d_meal_insulin` (1/5;
the insulin misses are EXPECTED, MASTER_PLAN row 1.46) and `validate_caffeine_ogtt`
(0/2; right sign since 2026-11-20 but ~+1% vs +25-60%). No drug effect is validated independently of the
OGTT-fitted default meal (semaglutide's PASS rests on it). Held-out group
(`validate_heldout()`, B2): 7 of 12 within +/-20% before B7, 10 of 12 after (seen).
Known calibration gaps since B7, pinned not re-tuned: GLP-1 meal response 50 pM on its
calibration setup (target 41), protein insulin-AUC ratio 0.22 (target 0.28).
Semaglutide is labelled **VALIDATED** with that caveat (commit 73c4328); every other drug
module is labelled **EXPERIMENTAL** (since 2026-11-16, single table in
`drug_interactions/validation_status.py`, `Drug.validation_status`; VALIDATED needs
every point of an independent cited benchmark within +/-20%). Drug runs carry an
"EXPERIMENTAL drug module" warning (semaglutide: "VALIDATED drug module" plus its caveat); drug benchmark reports carry `module_label`. A
new `Drug` value must get a row there or `test_drug_validation_status.py` fails.
SBML export encodes both meal curves (two-compartment since 2026-11-17); it still omits the
dG-driven first-phase secretion term, drugs and all coupled extensions (listed in
sbml_export.py's module docstring). Most
other calibrated mechanisms have no registered benchmark (MASTER_PLAN Section 2a,
"Systemic gap").

Caveat on the OGTT PASS: `GutAbsorptionModel.ogtt_calibrated_glucose()` was fitted to
three of the OGTT glucose points it is then checked against, so the glucose-timing
agreement is partly circular.

Key literature-anchored core parameters (each derived from independent data first,
then checked against the benchmarks -- never curve-fit to them):
- `p3 = 9.032e-7` -- Bergman 1981 (PMID:7033284), Fig 6 caption SI.
- `KG_secretion = 5.63` -- Rebelato 2018 (PMID:30166558), GSIS EC50 - G_threshold.
- `first_phase_gain = 121.3` -- functional form from Mari 2002 (PMID:12388151).
- `p1 = 0.0317` -- consistent with Pacini 1998 (PMID:9575818).
- `GutAbsorptionModel.ogtt_calibrated_glucose()` -- the calibrated OGTT absorption
  preset; `for_gi(100)` equals it exactly.

Accepted, documented limitations (not open tasks):
- **GI-to-timing slope is unknowable from current literature** (Cocate 2011 and
  Louie 2012/2013 conflict), so the GI adjustment is an **unvalidated placeholder**
  (only GI=100 is calibrated). `GI_SLOPE_PLACEHOLDER_NOTE` (`gut_absorption.py`) is
  shown in `SimulationOutput.warnings` (and REST) for any `glycaemic_index < 99.5`,
  in `gi_meal_scenario()` descriptions, and in the lay report's glucose section.
- **Meal absorption is a two-compartment first-order (Bateman) curve**, not the
  Lehmann-Deutsch model (whose gastric emptying has a maximal rate). It is the
  default for every meal without `kabs` (since 2026-11-15); a meal with no GI uses
  `GutAbsorptionModel.default_meal()` = the OGTT rate constants, a MODEL ASSUMPTION.
  Known costs: the GLP-1 meal response (calibrated on the old curve, before B7) now
  peaks ~39 pM on the default meal (Orskov 41 +/- 5) but ~50 pM on its own calibration
  setup -- documented, not re-tuned; and acarbose/miglitol (calibrated on `kabs`) delay
  the peak ~x2.1 instead of ~x1.5 on any two-compartment meal, default or GI (the run
  warns; GI meals were left unaffected until 2026-11-17).
- **No IV-bolus event at the scenario layer** -- `validate_bergman_iv()` runs at the
  ODE level via `bergman_ivgtt_y0()`.
- **A cited WIDE-tolerance miss gives the verdict MISS** (since 2026-11-12; it used
  to stay PASS). Reports name each missed point and its deviation
  (`_wide_failure_notes()`). Unverified points never change the verdict; a report
  with no cited point at all gives **NO DATA**, never PASS (since 2026-11-14).

Process lessons from the calibration work: a fit sitting at the edge of a
literature range *because* it satisfies a benchmark is still circular; check a
formula's algebra before re-deriving a constant; scanned pre-1990s PDFs' captions
extract more reliably than their tables. Full history: `HISTORY.md`.

---

## PAST WORK LOG, PART 1 (2026-07-28 through early Track 4 work)

**Moved to `HISTORY.md` on 2026-08-17** (process-improvement item Q3) --
this section used to hold the full chronological "### DONE --" log of
every completed milestone from Track 1's own kickoff through the start
of Track 4 (metabolism/FBA) work. See `HISTORY.md`'s own Part 1 for the
unabridged record. New completed-milestone entries go into `HISTORY.md`
going forward (append to its end), not here.


## AGENT TOOLING

Additive to the ECC rules in `~/.claude/rules/ecc/**`, which stay primary.
- **claude-mem** (plugin) records sessions automatically; use `mem-search` before
  re-deriving something that may already have been investigated. Treat what it
  surfaces like an agent report: cross-check against the code.
- **better-code-review-graph** (MCP) -- a structural graph for review agents.
  After a gap, run an incremental `graph update` with `base: <ref>`; full rebuilds
  are for a first index or a corrupted graph.
- **claude-code-setup** -- `claude-automation-recommender`. Its 2026-07-27 run
  rejected project-local reviewer agents (duplicates) and per-edit full-suite hooks
  (too slow; they would fight the stash-based RED/GREEN workflow).
- Project agents in `~/.claude/agents/`: `computational-biologist` and
  `literature-curator` (read-only), `database-architect` (local-only infra).

Full bootstrap notes: `HISTORY.md`, "Archived from CLAUDE.md (2026-10-20)".

---

## CITATION RE-AUDIT PROTOCOL (formalized MASTER_PLAN 5.5)

This project has run an ad hoc "periodic citation re-audit" twice
already (a Track 2/3-range pass, then the Track 1-range pass at commit
`9eb1e78` -- see that entry in HISTORY.md for its exact method). Both were
triggered by "someone noticed it's been a while," not a named rule.
This section makes the cadence and method explicit so a future session
(or the autonomous loop) can check whether one is due without
re-deriving the question from scratch.

**Trigger** (either condition, checked during Step 1 of
`archive/AUTONOMOUS_LOOP_PROMPT.md`'s own orientation, archived 2026-11-17): a re-audit is due when,
since the commit recorded in the most recent "citation re-audit" DONE
entry's own marker line —
1. **10 or more commits** have landed (`git log --oneline
   <last-audit-commit>..HEAD | wc -l`), OR
2. **15 or more new/changed `PMID:`/`pubmed_id=` occurrences** have
   been introduced (`git log <last-audit-commit>..HEAD -p -- '*.py' |
   grep -c '^+.*\(PMID:\|pubmed_id=\)'`, a cheap over-count since it
   also counts moved/reformatted lines -- treat it as a trigger
   threshold, not a precise count).

Whichever fires first. Neither condition is exact science -- both are
deliberately cheap-to-check proxies for "enough new citation surface
has accumulated that a fabrication or transcription error could be
hiding in it," matching this project's own established willingness to
use a good-enough engineering threshold rather than an unjustified
precise one (see e.g. `energy_bounds.py`'s own `0.15` mean-relative-
width UQ threshold).

**Method** (same as both prior passes, now the standard, not
reinvented per pass):
1. Identify which source files carry citations added or touched since
   the last audit (`git log <last-audit-commit>..HEAD --name-only` on
   the relevant packages, cross-referenced against `PMID:`/
   `pubmed_id=` grep hits).
2. Sample from those files -- prioritize the newest/most-touched over
   an exhaustive full-repo sweep (105+ PMIDs exist platform-wide; a
   full re-fetch every pass is disproportionate to the actual risk,
   which concentrates in recently-transcribed citations, not long-
   stable ones already audited before).
3. Live-fetch every sampled PMID via `PubMedFetcher.fetch_by_pubmed_id()`
   -- never trust a cached/remembered abstract, always re-fetch from
   real NCBI.
4. Cross-check title/author/journal/year (not just "does a paper with
   this PMID exist") against what the source file's own comment claims
   -- this is specifically the error class KNOWN BUGS FIXED #13 caught:
   a citation verified correctly during research but mistyped when
   transcribed into the actual code comment.
5. Record the result as a `### DONE (negative result, no code change) —
   Periodic citation re-audit, <range> (MASTER_PLAN 1.12)` entry (clean
   result) or a `### DONE — Fixed a mistranscribed citation, <range>`
   entry (if one is found) -- either way, **the entry's own marker
   line must record the commit hash the audit ran at**, in this exact
   form so the next orientation step's `git log` command above can
   find it automatically:

   `**Audit marker**: ran at commit <hash>, N PMIDs sampled from <range>.`

**Not a MASTER_PLAN row of its own per pass** -- row 1.12 stays the
permanent home for every future re-audit entry (matching its own
"Recurring, ~1 chunk each time" sizing), the same as this session's two
passes already used it. This section (5.5) is the protocol; 1.12 is
where each individual pass's result is logged.

---

## PAST WORK LOG, PART 2 (continuation -- Track 4/5/6 through this session)

**Moved to `HISTORY.md` on 2026-08-17** (process-improvement item Q3) --
continuation of the chronological "### DONE --" log immediately
following the citation re-audit protocol's own first real use: the rest
of Track 4 (metabolism/FBA), Track 5 (knowledge-infrastructure curation
campaign), the Track 6 parking decision, and this session's own
process/quality-improvement work (A6/F1, C1, D1/D2/Q5, C2/C4). See
`HISTORY.md`'s own Part 2 for the unabridged record.


## COMPOSITION-BOUNDARY BEHAVIORAL AUDIT PROTOCOL (formalized MASTER_PLAN 10/D1)

Same motivation and shape as the CITATION RE-AUDIT PROTOCOL above (MASTER_PLAN
5.5), for a different failure mode. That protocol catches a citation that's
wrong; this one catches a mechanism whose citation is *correct* but whose
calibrated effect never reaches the simulated output because of an
interaction with an existing safety clamp or another active modifier --
exactly what a `/verify` pass (2026-08-17) found for tirzepatide: three
literature-calibrated doses (5/10/15mg) producing byte-identical REST
output for this platform's own default T2D preset, because
`MIN_GB_THRESHOLD_MARGIN_MM` floored every one of them at the same
value. Fixed for that one case (row 1.23, `clamp_warnings` now surfaces
this), but the underlying discovery method -- actually running the
production REST API across a realistic dose/disease-state range,
instead of trusting isolated unit tests -- is worth scheduling, not
just applying once.

**Trigger**: due when 5 or more new drug/disease-mechanism milestones
(Track 1 rows adding a `Drug`/`DiseaseState` value, a new modifier
function, or a new `PathwayParameters` field touched by personalisation)
have landed since the last audit's own marker line.

**Method**:
1. Identify every drug/disease-mechanism milestone added or changed
   since the last audit (`git log <last-audit-commit>..HEAD --name-only`
   on `packages/physiology`/`packages/ode-engine`, cross-referenced
   against new `Drug`/`DiseaseState` enum values).
2. For each one, launch the real production REST server (`uvicorn
   biosim_rest:create_app --factory`, matching the `/verify` skill's
   own convention of driving the actual runtime surface rather than
   importing-and-calling) and drive `POST /simulate` across: (a) the
   mechanism's own full calibrated dose/severity range, (b) at least
   one disease-state preset the mechanism is likely to compose with
   (not just the healthy default), (c) both a fasted and a real-meal
   scenario.
3. Confirm the calibrated signal actually differentiates across the
   swept range in the JSON response -- not just "differs from a
   no-mechanism baseline," but genuinely varies dose-to-dose/severity-
   to-severity the way the underlying citation says it should.
4. If a clamp/interaction silently flattens the signal (as tirzepatide's
   did): check first whether `PhysiologicalPersonalizationResult
   .clamp_warnings` already surfaces it (row 1.23's own fix) -- if so,
   this is a known, disclosed limitation, not a new finding, just
   confirm the warning fires in the swept scenario. If the signal is
   flattened with NO warning anywhere, that is a new finding: fix it
   the same way 1.23 did (surface it, don't silently accept it), or add
   it to the population x mechanism coverage matrix (MASTER_PLAN 10/D2)
   as a documented limitation if a code fix isn't warranted yet.
5. Record the result as a `### DONE (negative result, no code change) —
   Composition-boundary behavioral audit, <range>` entry (clean result)
   or a `### DONE — Fixed a composition-boundary flattening bug, <range>`
   entry (if one is found) -- either way, the entry's own marker line
   must record the commit hash the audit ran at, in this exact form:

   `**Audit marker**: ran at commit <hash>, N mechanisms swept from <range>.`

**Not a MASTER_PLAN row of its own per pass** -- Section 10/D1 stays the
permanent home for this protocol; log each pass's result there or in a
CLAUDE.md DONE entry referencing it, matching the citation-audit
protocol's own convention.

---

## HOW TO CONTINUE A SESSION

1. Work in the repository root directly (no remote is configured).
2. Activate the venv (`.venv\Scripts\activate`) and run `python -m pytest -q`.
3. Read this file, then `MASTER_PLAN.md` Section 9 ("How to pick the next chunk")
   for the current priority order, and the last few `HISTORY.md` entries.

---


## COMMIT PATTERN

Conventional commits (`feat:`, `fix:`, `docs:`, `refactor:`, `test:`...) on
`master`. Check `git status` first and stage explicitly by path, never `git add -A`
or `git add .` (other sessions may have work in the tree). Record the test count
in the commit body and the `HISTORY.md` entry. Never push -- no remote exists.
