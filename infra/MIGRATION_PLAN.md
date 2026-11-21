# Biosim persistence-layer migration plan (Track 3)

Status: **scoping only — nothing in this plan has been wired into any real
code path yet.** `docker-compose.yml` and everything under `infra/` are
purely additive infrastructure; bringing the containers up or down has zero
effect on the existing test suite or on any simulation result. This
document describes how a *future* implementation pass would proceed, and
records the design decisions and open questions produced while scoping it.

Migration targets:

- `packages/parameter-db/src/biosim_parameter_db/` (`ParameterStore`) →
  PostgreSQL, schema in `infra/postgres/`.
- `packages/id-resolver/src/biosim_id_resolver/` (`IDRegistry`) → Neo4j
  Community Edition, schema in `infra/neo4j/`.

Both targets stay **local-only, zero-cost** (docker-compose, never a
managed/cloud database — CLAUDE.md's STRATEGIC ROADMAP decision 2), and
the migration is **incremental**: dual-write, then shadow-read, proven
equivalent against the real test suite and the platform's clinical
validation benchmarks before any cutover — never a big-bang replacement.
The existing in-memory/file-based stores and their public APIs keep
working, unchanged, at every intermediate step below.

---

## Why a graph DB, and which one

`IDRegistry`'s actual query shape is graph-shaped, not relational: resolve
an external accession to a canonical entity, walk a chain of merges to the
current live entity, and (eventually) express pathway topology as real
edges rather than name strings. Forcing that into PostgreSQL would mean
recursive CTEs standing in for graph traversal — workable, but fighting
the tool. Doing the migration in two systems adds real operational
overhead (two containers, two schemas, two backfill scripts) — accepted
here specifically because ParameterStore's actual pattern (append rows,
rank-order a query) is exactly what a relational store is for, and
inventing one hybrid schema to serve both would compromise both.

**Neo4j Community Edition**, not a cloud/managed graph offering (Neo4j
Aura is explicitly out of scope) and not a paid source-available
alternative:

- Licensing checked, not assumed: Neo4j CE server is GPLv3, the official
  `neo4j` Python driver is Apache-2.0. Self-hosting and talking to it over
  Bolt carries no copyleft obligation onto Biosim's own source — GPLv3's
  copyleft attaches on *redistribution* of the server itself (e.g.
  shipping Biosim as an appliance bundling the Neo4j binary), not on
  network calls to a separately-installed instance.
- Alternatives considered: Kùzu (embedded, no server process — a real
  contender if the container-per-service overhead ever becomes a problem,
  see open question 3), Apache AGE (Postgres extension — would collapse
  back to one database, worth a second look if the two-container overhead
  turns out not to be worth it), Memgraph and ArangoDB (both rejected as
  defaults specifically because their licences, BSL/BUSL respectively,
  are source-available but not OSI-approved open source — a real,
  material difference for a project with a standing free/open-source-only
  posture, not a technicality).

---

## What's already built (this pass)

- `docker-compose.yml` — Postgres 17 + Neo4j 5.26 Community, both bound to
  `127.0.0.1` only (not reachable from the local network), throwaway
  local-dev credentials documented as such, healthchecks, a one-shot
  `neo4j-init` service applying the Cypher schema (Neo4j has no
  `/docker-entrypoint-initdb.d` equivalent).
- `infra/postgres/001_schema.sql` — `parameter`/`literature_source`/
  `parameter_supersession`/`parameter_conflict`/`conflict_resolution`/
  `migration_audit` tables plus `ref_*` mirrors of the Python rank dicts.
  Append-only by construction: the two in-place mutations the Python store
  actually performs (`Parameter.supersede()` sets `is_superseded=True`;
  `ParameterStore.resolve_conflict()` sets `resolved=True` +
  `resolution_notes`) are re-modelled as append-only event rows, so no
  legitimate write path ever needs `UPDATE`/`DELETE` on a fact table.
- `infra/postgres/002_query_api.sql` — `v_parameter` (reconstructs a
  `Parameter` one-for-one, deriving `is_superseded` from the event table)
  and `parameter_query()`/`parameter_best()`/`parameter_all_for_subject()`/
  `parameter_stats()`/`would_conflict_with()`, each an exact translation of
  the matching `ParameterStore` method — including its exact-match (not
  wildcard) object/reaction scoping and its stable-sort tiebreak order,
  both verified against the running Python store directly rather than
  read off the source.
- `infra/postgres/003_append_only_enforcement.sql` — two layers
  (`GRANT`/`REVOKE` privileges for a `biosim_writer` NOLOGIN role, plus
  `BEFORE UPDATE/DELETE/TRUNCATE` triggers that reject even the owner and
  superusers) enforcing rule 5 ("never delete or overwrite parameters") as
  a database guarantee instead of pure discipline, with a self-test that
  fails loudly if the guard triggers aren't actually installed.
- `infra/neo4j/001_constraints.cypher` / `002_reference_queries.cypher` —
  `(:Entity)`/`(:ExternalId)`/`(:OntologySource)` node model,
  `IDENTIFIES`/`SUPERSEDED_BY` immutable edges, and translations of
  `resolve()`/`get_record()`/`cross_reference()`/`search_by_name()`.

## A real, separate finding from this scoping pass (already fixed)

Designing the graph model required reasoning precisely about what a
multi-hop merge chain (entity A merged into B, then B itself later merged
into C) needs to behave like — and that surfaced three genuine, verified
bugs in the CURRENT in-memory `IDRegistry`, independent of any migration:

1. `get_record(A)` reported "the current canonical ID is B" as final,
   even when B was itself already deprecated into C by the time A was
   looked up.
2. `_name_index` was never updated on merge — `search_by_name()` for the
   old entity's original name kept returning a dead canonical id forever.
3. `to_dict()`/`from_dict()` only round-tripped each surviving record's
   own direct deprecation list, silently losing earlier hops in a
   multi-step chain (A's mapping was nowhere in the serialized form after
   A→B→C).

All three were verified directly against the real code (not assumed from
the graph design's own reasoning) and fixed in `registry.py` this same
session, with 4 new regression tests, RED/GREEN-proofed via `git stash` —
see CLAUDE.md's own DONE entry for the fix. **This means open question 4
below is already resolved**, not merely scoped: the in-memory registry's
`search_by_name()` now returns the live, repointed canonical id, which is
Cypher option (A)'s outcome, reached by repointing rather than filtering.
`002_reference_queries.cypher`'s own inline comment describing this as
still-open is stale as of this fix and should be updated to reference it
rather than presenting it as an undecided behaviour choice.

---

## Dual-write / shadow-read rollout plan

**Step 1 — Bring the containers up, apply the schema, nothing else.**
`docker-compose up -d`. Confirm both healthchecks pass and the Postgres
append-only self-test in `003_append_only_enforcement.sql` reports "OK".
Zero application code changes. Verification: the existing full test suite
must still pass unchanged (it will — nothing imports anything from
`infra/` yet). Check `git log`/a fresh test run for the current count
rather than trusting any number written in this document — this is a
multi-session repo and the count drifts (this document itself once said
"1551" in five places, stale by the time this note was added).

**Step 2 — One-time backfill from a `ParameterStore`/`IDRegistry`
snapshot into the new stores.** A backfill script reads
`ParameterStore.save()`/`IDRegistry.to_dict()`'s existing JSON snapshot
format (already used for file persistence today) and inserts rows via
`ON CONFLICT (parameter_id) DO NOTHING` / Cypher `MERGE` — both
idempotent, so a re-run after a partial failure is safe. Populates
`migration_audit` with counts. The backfill NEVER writes back to the
existing file-based store — a failed or buggy backfill run cannot damage
the thing it's migrating from. **Must insert in exact snapshot array
order, single-threaded** — `parameter.seq` is `GENERATED ALWAYS AS
IDENTITY` (assigned by Postgres at INSERT time, cannot be supplied by the
backfill without `OVERRIDING SYSTEM VALUE`) and `parameter_best()`/
`get_all_for_subject()`'s SQL and Python tiebreaks both now depend on seq
matching insertion order exactly (see the get_all_for_subject() fix
below) — a re-run after a partial failure in a *different* order would
silently produce a different tie order than the original store, which
`ON CONFLICT DO NOTHING`'s idempotency does not protect against on its
own. Flagged for full resolution alongside open question 7.

**Step 3 — Shadow-read comparison harness (the actual equivalence
gate).** For every test in the platform's full suite AND every registered
clinical validation benchmark (`ModelValidator.validate_all()` — 6 as of
this writing, itself subject to the same "check, don't trust a written
number" rule above), run the scenario twice: once against the existing
in-memory store (today's behaviour, untouched) and once by additionally
routing every
`ParameterStore.query()`/`.best()`/`.get_all_for_subject()`/`.stats()` and
`IDRegistry.resolve()`/`.get_record()`/`.search_by_name()`/`.stats()` call
through a **read-only** adapter hitting Postgres/Neo4j instead, comparing
results field-by-field. This is a genuine gate, not a formality — it must
actually run and actually pass, not be asserted as "should be equivalent
by construction." Concretely:

- Every parameter/entity field returned by the new-store adapter must
  equal the old store's, for every call the existing test suite already
  makes (no new test scenarios need writing — the existing suite IS the
  equivalence corpus).
- **No special cases.** `parameter_all_for_subject()` was a KNOWN,
  documented exception at the time this step was first drafted (see
  open question 5, now RESOLVED below) — `ParameterStore
  .get_all_for_subject()` was fixed to match the SQL function's already-
  correct behaviour, so the comparator now has zero sanctioned
  divergences. The adapter must GROUP this function's flat row set by
  `kind` while PRESERVING row order (matching Python's
  `dict[kind -> list[Parameter]]` shape) rather than re-sorting
  client-side.
- Any other divergence found during this step is a genuine migration bug
  and blocks cutover until fixed — not "close enough."

**Step 4 — Dual-write, still reading only from the old store.** Every
`ParameterStore.add()`/`IDRegistry.register()` call also writes to
Postgres/Neo4j (application code gains a second write, behind a feature
flag defaulting off), while every read still comes from the existing
in-memory store exclusively. Run the full shadow-read comparison from
Step 3 continuously in this configuration, not just once, until confident
the two stores never drift over ordinary usage (not just at one
snapshot-and-compare instant).

**Step 5 — Cutover criteria (all required, none optional):**
1. Every test in the full suite passes with reads routed through the new
   store, not just the ones directly exercising `ParameterStore`/
   `IDRegistry`.
2. All clinical validation benchmarks (`ModelValidator.validate_all()`)
   produce numerically identical results, not just "still PASS" — a
   benchmark could still pass at a shifted value inside its tolerance
   band, which would hide a real migration bug.
3. Shadow-read has run clean (zero unexplained divergences) for a
   deliberately long soak period covering realistic usage patterns, not
   a single comparison run.
4. Every open question below is explicitly resolved (a real decision
   made and documented), not silently defaulted.
5. A rollback path exists and has itself been tested: flipping the read
   flag back to the in-memory store with zero data loss, since the old
   store keeps receiving every write throughout dual-write.

**What never changes, at every step above:** the existing
`ParameterStore`/`IDRegistry` public APIs, return types, and the file-based
snapshot persistence they already have. Nothing currently depending on
them (which is most of this platform) needs to know a migration is in
progress until Step 5 actually completes.

---

## Open questions (must be resolved before cutover, not before scoping)

0. **Provenance verification fields (added 2026-11-15).** `Provenance` now
   carries `verification_status` (a `VerificationStatus` value, default
   `UNCHECKED`) and `verification_quote`. `001_schema.sql` has no columns for
   them yet. Before any dual-write, add `verification_status text NOT NULL
   DEFAULT 'UNCHECKED'` with a CHECK on the enum values, plus
   `verification_quote text NOT NULL DEFAULT ''`. Otherwise the equivalence
   check will drop the audit result.

1. **`ConflictRecord` field parity** (referenced in `001_schema.sql`'s
   `parameter_conflict` table comment). The Python `ConflictRecord` carries
   `subject_canonical_id` + `kind` but NOT `object_canonical_id`/
   `reaction_canonical_id`, so two conflicts on the same subject+kind for
   DIFFERENT objects are indistinguishable in the record itself — the SQL
   schema exposes the scope via a join instead of inventing new columns
   the Python type doesn't have. Should `ConflictRecord` itself gain those
   fields? That's a Python-side design decision independent of this
   migration, not something the schema should force by fiat.

2. **Kùzu/embedded-graph reconsideration.** Two containers (Postgres +
   Neo4j) is real, non-trivial local operational overhead for what is
   currently a small (tens of entities) registry. Worth revisiting once
   real usage patterns at Track-3-migration scale are known — Kùzu (no
   separate server process) or Apache AGE (a Postgres extension, avoiding
   a second database engine entirely) are both real alternatives to
   re-evaluate, not commitments to keep Neo4j specifically.

3. **Neo4j Community Edition's limitations are accepted, not solved.** No
   role-based access control, no online/hot backup, single user database
   — meaningfully weaker than Postgres's DB-enforced append-only
   guarantees (`003_append_only_enforcement.sql`). `IDENTIFIES` edges are
   immutable BY APPLICATION CONVENTION in the current Cypher queries, not
   by a database-level guarantee the way Postgres's triggers provide one.
   If CE's limitations become a real operational problem, this is the
   point to revisit whether the two-database split (question 2) or a
   different graph engine is the right call, rather than trying to work
   around CE's constraints indefinitely.

4. **RESOLVED, not open** (see "a real, separate finding" above) —
   `search_by_name()`'s deprecated-id-return behaviour. Fixed at the
   Python source this session; `002_reference_queries.cypher`'s inline
   comment describing this as an open behaviour choice should be updated
   to reference the fix rather than presenting option (A) vs (B) as still
   undecided.

5. **RESOLVED, not open** — `parameter_all_for_subject()` vs.
   `ParameterStore.get_all_for_subject()`. Option (a) was adopted:
   `get_all_for_subject()` was fixed in `store.py` to accumulate across
   scopes instead of overwriting, sorted with an explicit ascending
   insertion-sequence tiebreak (a monotonic counter mirroring Postgres's
   `parameter.seq GENERATED ALWAYS AS IDENTITY`, assigned in both
   `add()` and `load()`) rather than the naive
   `sorted(key=self._priority_score, reverse=True)`, which would have
   left a second, subtler divergence in place -- NOT because
   `reverse=True` reverses tie order (Python's sort is stable regardless
   of `reverse=True`; verified directly, ties keep their original
   relative input order either way), but because the naive fix's
   accumulation loop groups by which SCOPE KEY was first created, not
   true global insertion order: interleaving inserts into different
   scopes (scope A, then scope B, then scope A again) puts both A
   entries ahead of B's in the naive version, even though B was added
   in between -- confirmed empirically by running both the naive and
   real implementations side by side against identical input. Verified
   zero production callers existed at fix time (only the method's own
   pre-existing test called it), so this carried no behaviour-change
   risk to any real caller. See
   `infra/postgres/002_query_api.sql`'s updated comment on
   `parameter_all_for_subject()` for the full resolution and the shape
   contract a shadow-read adapter must honour (group by kind, preserve
   row order, never re-sort client-side). Regression tests:
   `packages/parameter-db/tests/test_parameter_db.py::TestObjectScopedParameters`
   (`*_returns_every_scope_not_just_the_last`,
   `*_ties_break_by_global_insertion_order`,
   `*_orders_by_confidence_then_organism_then_system`,
   `*_include_superseded_filters_across_scopes`,
   `*_groups_by_kind_across_scopes`). Step 3's shadow-read comparator no
   longer needs any special case for this function.

   **A second, independent bug was found and fixed in the same pass**:
   `ParameterStore._deserialize_parameter()` never passed `created_at`
   or `Provenance.extracted_at` through to their respective constructors
   on `load()`, so both were silently regenerated via
   `default_factory=datetime.now(timezone.utc)` on every save/load
   round trip — reproduced directly. This matters specifically for this
   migration: `parameter.created_at` is `timestamptz NOT NULL` in the
   append-only, trigger-protected `001_schema.sql` schema, so a backfill
   reading via `ParameterStore.load()` would have permanently stamped
   every row with the migration's wall-clock time instead of the
   parameter's real creation time — unfixable afterward without
   deliberately disabling the append-only guard triggers. Fixed with a
   graceful `.get()`-based round trip (`datetime.fromisoformat`, falling
   back to the default for pre-fix snapshots that never wrote the
   field), matching `IDRegistry.from_dict()`'s existing degradation
   convention. Regression tests:
   `TestPersistence::test_created_at_survives_save_load_round_trip`,
   `::test_extracted_at_survives_save_load_round_trip`.

6. **Pathway/reaction topology as real graph edges.** IDRegistry's model
   above only covers entity resolution (rule 6: Reactome IDs for pathways
   only, never assigned to a protein entity) — it does not yet represent
   pathway membership or reaction topology as actual `(:Entity)`-
   `(:Pathway)` edges, which is the other half of "a graph database for
   entity cross-references AND pathway relationships" from the roadmap's
   own decision 3. Out of scope for this pass; a real design question for
   whenever pathway-level queries (e.g. "every protein in glycolysis")
   become an actual access pattern, not before.

7. **RESOLVED, not open** — who runs the backfill, and when, in a
   multi-session/concurrent-editing environment. Adopted exactly the
   rule this entry's own earlier draft leaned toward, made explicit and
   operational rather than implicit:

   - **The backfill NEVER opens a live `ParameterStore`/`IDRegistry`
     object, in-process or otherwise.** It takes a file PATH to a JSON
     snapshot as its only input, matching `migration_audit
     .source_snapshot`'s own documented contract (path/hash of a
     snapshot file, never a live object reference,
     `infra/postgres/001_schema.sql:440-444`).
   - **Copy-then-read, not read-in-place.** The backfill script copies
     the snapshot file to a read-only temp location before reading it,
     so a concurrent session's own `store.save()` call to the SAME
     path mid-backfill cannot produce a torn/partial read.
   - **The operator (not a session) decides when to run it**, at a
     point where no other session is actively writing to the specific
     snapshot file being migrated — this project has no file-locking
     mechanism and none is being added for this; the rule is procedural
     (pick a quiet moment, verify `git status` on the snapshot file
     first), matching this project's own existing convention of
     checking `git status` before every stash/commit rather than
     building automated coordination machinery for a single-operator
     project.
   - **The backfill NEVER writes back to the file-based store** (already
     stated in Step 2 above) — so even a genuinely bad read (a torn
     file, a mid-write snapshot) can only produce a bad Postgres/Neo4j
     insert, never corrupt the source of truth it's migrating from.

   **A second, related invariant surfaced while resolving this
   question, not previously named as its own open question**: the
   backfill must insert parameter rows in EXACT snapshot array order,
   single-threaded. `parameter.seq` is `bigint GENERATED ALWAYS AS
   IDENTITY` (`001_schema.sql:214`) — Postgres assigns it at INSERT
   time in whatever order rows arrive, and cannot be supplied directly
   without `OVERRIDING SYSTEM VALUE`. `parameter_best()`'s SQL and
   `ParameterStore.get_all_for_subject()`'s Python (see the DONE entry
   above) both now tiebreak on `seq ASC`/`self._insertion_seq`
   specifically because they're meant to represent the SAME concept —
   an out-of-order backfill (parallel insertion, or a re-run after a
   partial failure in a different order) would silently diverge the
   two stores' tie order for any subject+kind with genuinely tied
   parameters, with `ON CONFLICT (parameter_id) DO NOTHING` providing
   zero protection against it (idempotency guards re-insertion of the
   same row, not insertion order across separate runs). Concretely: run
   the backfill as one single-threaded pass over the snapshot's
   `parameters` array in its on-disk order, and verify post-backfill
   with a query asserting snapshot array index order equals `seq`
   order for every row before considering a backfill run complete.

8. **RESOLVED, not open** — the Neo4j invariant-4 decision
   (MASTER_PLAN Track 3.4). `002_reference_queries.cypher`'s invariant 4
   used to reject ANY active entity with two accessions from the same
   ontology source, contradicting this schema's own "no constraint
   stricter than the Python invariant" rule: `registry.py`'s own
   already-shipped merge-inheritance fix makes exactly this a
   legitimate, routine outcome (merging entity A, UniProt P00001, into
   entity B, UniProt P99999, correctly produces B.cross_references =
   {'uniprot': ['P99999', 'P00001']} -- verified directly against a
   live `IDRegistry`). Worse than a false-positive constraint: the
   PRIOR `cross_reference()` query used `ORDER BY accession ASC LIMIT 1`
   (alphabetically-first) specifically *because* invariant 4 was
   assumed to guarantee at most one candidate ever existed -- since that
   assumption is false, the query would have silently returned
   'P00001' for the scenario above, disagreeing with Python's real
   answer ('P99999', the merged entity's own restated/primary
   accession, per `ResolutionRecord.get_external_id()`'s `refs[0]`
   semantics) -- a genuine equivalence-suite-breaking divergence, not
   just an overly strict check.

   Adopted: make Python's already-deterministic insertion order an
   EXPLICIT, persisted graph property rather than fighting it. Step 2's
   write path now stamps each `IDENTIFIES` edge with `rank` (0-based
   position in Python's own `cross_references[source]` list at backfill
   time, set `ON CREATE` only, matching this document's own question 3
   -- `IDENTIFIES` edges are immutable by convention). `cross_reference()`'s
   query now orders by `r.rank` instead of the accession string,
   matching Python exactly. Invariant 4 itself now checks the REAL
   invariant `cross_reference()` depends on -- no two accessions may
   share the same rank for one (entity, source) pair (an actual
   ambiguous-primary bug) -- rather than rejecting legitimate
   multiplicity outright. The backfill script (3.7, not yet written)
   must assign `rank` from each accession's real list-index position,
   in order, for this to hold -- flagged here so that requirement isn't
   rediscovered from scratch when 3.7 is picked up.

   **Disclosed limit (code review)**: invariant 4 catches ambiguity
   (two accessions sharing one rank), not misassignment (unique-but-
   wrong-order ranks from a backfill bug -- e.g. a reversed list). Zero
   rows from invariant 4 means "unambiguous," not "correct." Only the
   equivalence suite (question/step 3.8, not yet built) can catch a
   unique-but-wrong-order bug, by comparing against Python's own real
   answer directly -- see `002_reference_queries.cypher`'s own updated
   invariant-4 comment for the full reasoning.

---

## Explicitly NOT done this pass

- **Update**: `biosim_parameter_db.store.ParameterStore` HAS since been
  modified (see the `get_all_for_subject()`/timestamp-round-trip DONE
  entry above) -- but only to fix real, independently-discovered bugs
  this migration-equivalence work surfaced, not to build any actual
  migration machinery. `biosim_id_resolver` remains completely
  untouched this pass; its own earlier IDRegistry fixes (transitive
  deprecation chains, etc.) predate this document.
- **Update (MASTER_PLAN 3.6, done)**: the pure-Python snapshot READER/
  VALIDATOR half of Step 2 now exists --
  `biosim_parameter_db.snapshot_validator.validate_snapshot()`
  (`packages/parameter-db/src/biosim_parameter_db/snapshot_validator.py`).
  It implements exactly the operational contract question 7 resolved
  (file-path-only input, copy-then-read, never writes back) and checks
  the array-order/seq-order invariant plus every supersession/conflict
  reference resolves -- but emits NO SQL/Cypher and makes no database
  connection, so it needed no Docker boot and is fully tested today. The
  actual INSERT/MERGE statement generation (row 3.7) is still not
  written -- that's the next, larger half of Step 2, gated on the
  disposable Docker validation boot below actually happening first.
- No shadow-read comparison harness exists yet (Step 3 above is a plan).
- No dual-write flag or adapter layer exists yet (Step 4).
- Pathway/reaction graph topology (open question 6) is not designed.
- The docker-compose stack has not been booted even once this pass
  (Docker Desktop's engine was not running in this environment when
  checked) -- the disposable validation boot recommended after this
  document's own hygiene fixes (README's Docker-Desktop-must-be-running
  precondition, the re-apply-changed-SQL commands, the `.env`/`.gitignore`
  fix) has not happened yet. Do this first, before writing the backfill
  script -- it is the cheapest way to catch a real defect this pass's
  own SQL review predicted but could not confirm without execution
  (function-body `search_path` resolution -- see `002_query_api.sql`'s
  and `003_append_only_enforcement.sql`'s `SET search_path` usage).

This is intentional: per this project's own standing discipline, a large,
multi-phase infrastructure migration gets scoped and reasoned about
before code is written for it, the same way every biological-pathway
milestone this session began with a design pass before implementation
(e.g. the SGLT2i renal-clearance pathway's "static vs. dynamic"
resolution, or the Sobol-sensitivity method choice) rather than starting
from whichever approach was easiest to code first.
