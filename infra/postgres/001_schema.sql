-- =====================================================================
-- Biosim parameter store — PostgreSQL schema (Track 3).
--
-- Target of migration: packages/parameter-db/src/biosim_parameter_db/
--   types.py   Parameter / ParameterValue / Provenance / ExperimentalContext
--   store.py   ParameterStore (append-only, conflict-detecting, ranked query)
--
-- This file is IDEMPOTENT: every statement is IF NOT EXISTS / ON CONFLICT
-- guarded, so it runs correctly both as a first-boot init script and when
-- re-applied by hand against a live volume.
--
-- ---------------------------------------------------------------------
-- DESIGN RATIONALE (the short version; full reasoning in
-- infra/MIGRATION_PLAN.md)
--
-- 1. APPEND-ONLY IS ENFORCED, NOT CONVENTIONAL. Rule 5 in CLAUDE.md says
--    the store is append-only, but the current Python implementation
--    performs exactly two in-place mutations:
--        ParameterStore.supersede()      -> old.is_superseded = True
--        ParameterStore.resolve_conflict() -> conflict.resolved = True + notes
--    Verified directly by running the store, not inferred from the code.
--    A naive port that made `is_superseded` a mutable BOOLEAN column would
--    make `REVOKE UPDATE` impossible. So both mutations are re-modelled as
--    APPEND-ONLY EVENT ROWS (parameter_supersession, conflict_resolution)
--    and the boolean is DERIVED in a view. With that, no code path ever
--    needs UPDATE or DELETE on a fact table, and 003_append_only_enforcement.sql
--    can revoke both outright.
--
-- 2. SELECTIVE NORMALISATION.
--    * literature_source is NORMALISED (N:1). Many parameters cite one
--      paper (Bergman 1981 backs several seeds), and "every parameter from
--      PMID X" is a real, recurring access pattern in this project — the
--      citation audits in CLAUDE.md (~35 fabricated PMIDs) are exactly that
--      query, run by hand.
--    * experimental_context is INLINED (1:1). It is an immutable value
--      object with no independent access pattern; in an append-only store no
--      update anomaly is possible, so normalising it would buy nothing and
--      cost a join on the hot path plus a NULL-in-content-hash trap.
--
-- 3. RANKED QUERY ORDER IS DATA, NOT CODE. store.py's CONFIDENCE_RANK /
--    ORGANISM_RANK / SYSTEM_RANK dicts become ref_* tables so the ORDER BY
--    in 002_query_api.sql can reproduce _priority_score() exactly.
--    !! These ranks are now DUPLICATED between store.py and this file.
--    That is precisely the "same fact in two places, nothing cross-checks
--    it" bug class CLAUDE.md has already fixed three times (Vg_L_per_kg,
--    MODEL_VERSION, gi_category). The migration plan therefore REQUIRES a
--    test that asserts the ref_* table contents equal the Python dicts,
--    in the same spirit as
--    test_model_version_matches_the_actual_ode_engine_not_a_separate_copy.
--
-- 4. NO CONSTRAINT MAY BE STRICTER THAN THE PYTHON INVARIANT. During
--    dual-write, any check the DB enforces but Python does not becomes a
--    write that succeeds in one store and fails in the other — a divergence
--    caused by the migration itself. Every CHECK below is a direct mirror
--    of an existing __post_init__ validation. Two deliberate, documented
--    exceptions are marked "DEVIATION" inline.
-- =====================================================================

CREATE SCHEMA IF NOT EXISTS biosim;
SET search_path = biosim, public;


-- ---------------------------------------------------------------------
-- Reference tables (mirrors of the Python enums and rank dicts).
-- These are the only tables in this schema that are legitimately
-- mutable: adding a new organism is a schema-evolution act performed by
-- the owner role, not an application write.
-- ---------------------------------------------------------------------

CREATE TABLE IF NOT EXISTS ref_parameter_kind (
    kind            text PRIMARY KEY,
    canonical_unit  text NOT NULL,
    description     text NOT NULL DEFAULT ''
);
COMMENT ON TABLE ref_parameter_kind IS
    'Mirror of biosim_parameter_db.types.ParameterKind. canonical_unit is '
    'documentation only: unit conversion happens at ingestion time in Python, '
    'never at query time, so this column is intentionally not enforced against '
    'parameter.unit.';

INSERT INTO ref_parameter_kind (kind, canonical_unit, description) VALUES
    ('Km',          'uM',       'Michaelis constant'),
    ('Vmax',        'uM/s',     'Maximum reaction rate'),
    ('Kd',          'uM',       'Dissociation constant'),
    ('Ki',          'uM',       'Inhibition constant'),
    ('EC50',        'uM',       'Half-maximal effective concentration'),
    ('IC50',        'uM',       'Half-maximal inhibitory concentration'),
    ('kcat',        's-1',      'Catalytic rate constant'),
    ('hill',        'dimensionless', 'Hill coefficient'),
    ('rate',        's-1',      'Generic rate constant'),
    ('fold_change', 'dimensionless', 'Expression or activity fold change'),
    ('half_life',   's',        'Molecular half-life'),
    ('abundance',   'uM',       'Basal concentration'),
    ('custom',      'varies',   'Explicitly named custom parameter')
ON CONFLICT (kind) DO NOTHING;


CREATE TABLE IF NOT EXISTS ref_confidence_tier (
    tier text PRIMARY KEY,
    rank smallint NOT NULL
);
COMMENT ON TABLE ref_confidence_tier IS
    'Mirror of biosim_parameter_db.store.CONFIDENCE_RANK. Higher rank wins in '
    'parameter_query(). Must be kept equal to the Python dict by test.';

INSERT INTO ref_confidence_tier (tier, rank) VALUES
    ('gold', 3), ('silver', 2), ('bronze', 1)
ON CONFLICT (tier) DO NOTHING;


CREATE TABLE IF NOT EXISTS ref_organism (
    organism text PRIMARY KEY,
    rank     smallint NOT NULL
);
COMMENT ON TABLE ref_organism IS
    'Mirror of biosim_parameter_db.store.ORGANISM_RANK. Note the ranks are '
    'deliberately NOT unique (rat = mouse = 3, other = unknown = 0); ties are '
    'broken downstream by parameter.seq, reproducing Python''s stable sort.';

INSERT INTO ref_organism (organism, rank) VALUES
    ('homo_sapiens', 5),
    ('rattus_norvegicus', 3),
    ('mus_musculus', 3),
    ('saccharomyces_cerevisiae', 1),
    ('escherichia_coli', 1),
    ('other', 0),
    ('unknown', 0)
ON CONFLICT (organism) DO NOTHING;


CREATE TABLE IF NOT EXISTS ref_experimental_system (
    system text PRIMARY KEY,
    rank   smallint NOT NULL
);
COMMENT ON TABLE ref_experimental_system IS
    'Mirror of biosim_parameter_db.store.SYSTEM_RANK.';

INSERT INTO ref_experimental_system (system, rank) VALUES
    ('clinical', 5), ('in_vivo', 4), ('ex_vivo', 3),
    ('in_vitro', 2), ('in_silico', 1), ('unknown', 0)
ON CONFLICT (system) DO NOTHING;


-- ---------------------------------------------------------------------
-- literature_source — the identity of a cited source, deduplicated
-- across parameters.
--
-- Mirrors the SOURCE-IDENTITY half of types.Provenance. The
-- EXTRACTION-SPECIFIC half (extracted_by, extracted_at, figure_or_table,
-- is_primary_source) stays on the parameter row, because two curators
-- extracting two different numbers from the same paper legitimately have
-- different extraction metadata for the same source.
-- ---------------------------------------------------------------------

CREATE TABLE IF NOT EXISTS literature_source (
    source_id          bigint GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    pubmed_id          text,
    doi                text,
    arxiv_id           text,
    source_description text NOT NULL DEFAULT '',
    created_at         timestamptz NOT NULL DEFAULT now(),

    -- Direct mirror of Provenance.__post_init__: a parameter cannot enter
    -- the database without a traceable origin.
    CONSTRAINT literature_source_has_a_traceable_origin CHECK (
        pubmed_id IS NOT NULL
        OR doi IS NOT NULL
        OR arxiv_id IS NOT NULL
        OR source_description <> ''
    ),

    -- Deterministic dedup key. Makes the migration idempotent (re-running
    -- it re-finds the same source row instead of inserting a duplicate)
    -- WITHOUT collapsing genuinely different citations: two parameters
    -- citing the same PMID with different source_description text remain
    -- two rows, exactly as two distinct Provenance value objects would.
    source_key text GENERATED ALWAYS AS (
        lower(coalesce(pubmed_id, ''))          || E'\x1F' ||
        lower(coalesce(doi, ''))                || E'\x1F' ||
        lower(coalesce(arxiv_id, ''))           || E'\x1F' ||
        lower(source_description)
    ) STORED
);

CREATE UNIQUE INDEX IF NOT EXISTS literature_source_key_uidx
    ON literature_source (source_key);
CREATE INDEX IF NOT EXISTS literature_source_pubmed_idx
    ON literature_source (pubmed_id) WHERE pubmed_id IS NOT NULL;

COMMENT ON COLUMN literature_source.source_key IS
    'Deduplication key using an ASCII unit-separator (0x1F) delimiter, which '
    'cannot appear in a real PMID/DOI/arXiv id. Enables ON CONFLICT DO NOTHING '
    'for an idempotent, re-runnable backfill.';


-- ---------------------------------------------------------------------
-- parameter — the fact table. Immutable rows, insert only.
-- ---------------------------------------------------------------------

CREATE TABLE IF NOT EXISTS parameter (
    -- The Python-assigned UUID is carried through as the primary key.
    -- This is what makes the backfill idempotent (ON CONFLICT DO NOTHING)
    -- and makes shadow-read comparison an exact id-for-id match rather
    -- than a fuzzy content match.
    parameter_id            uuid PRIMARY KEY,

    -- Monotonic insertion sequence. NOT decoration: ParameterStore.query()
    -- sorts with Python's STABLE sort, so equal-priority parameters are
    -- returned in insertion order and best() returns the FIRST one
    -- inserted. Verified by direct experiment. Without this tiebreaker the
    -- SQL ORDER BY would be non-deterministic among ties and best() could
    -- silently return a different parameter than the in-memory store —
    -- which, via ParameterLoader, would change a simulated result.
    seq                     bigint GENERATED ALWAYS AS IDENTITY NOT NULL UNIQUE,

    kind                    text NOT NULL REFERENCES ref_parameter_kind (kind),
    custom_kind_name        text,

    -- ParameterValue, inlined.
    value                   double precision NOT NULL,
    unit                    text NOT NULL,
    uncertainty_lower       double precision,
    uncertainty_upper       double precision,
    uncertainty_type        text NOT NULL DEFAULT 'none',

    -- Scoping triple. subject is required; object/reaction are the
    -- disambiguators added in the object-scoped/reaction-scoped fix.
    subject_canonical_id    text NOT NULL,
    object_canonical_id     text,
    reaction_canonical_id   text,

    confidence_tier         text NOT NULL DEFAULT 'silver'
                                 REFERENCES ref_confidence_tier (tier),

    -- Provenance: source identity by reference, extraction metadata inline.
    source_id               bigint NOT NULL REFERENCES literature_source (source_id),
    extracted_by            text NOT NULL DEFAULT '',
    extracted_at            timestamptz NOT NULL DEFAULT now(),
    figure_or_table         text,
    is_primary_source       boolean NOT NULL DEFAULT true,

    -- ExperimentalContext, inlined (see rationale 2 in the header).
    organism                text NOT NULL DEFAULT 'unknown'
                                 REFERENCES ref_organism (organism),
    experimental_system     text NOT NULL DEFAULT 'unknown'
                                 REFERENCES ref_experimental_system (system),
    cell_type               text,
    tissue                  text,
    temperature_celsius     double precision,
    ph                      double precision,
    context_notes           text NOT NULL DEFAULT '',

    curator_notes           text NOT NULL DEFAULT '',
    created_at              timestamptz NOT NULL DEFAULT now(),

    -- Immutable successor link, mirroring Parameter.supersedes_parameter_id.
    -- Set once at insert; the "is this parameter superseded" fact lives in
    -- parameter_supersession instead of being a mutable flag here.
    supersedes_parameter_id uuid REFERENCES parameter (parameter_id),

    -- Mirror of Parameter.__post_init__.
    CONSTRAINT parameter_custom_kind_needs_a_name CHECK (
        kind <> 'custom' OR custom_kind_name IS NOT NULL
    ),

    -- Mirror of ParameterValue.__post_init__, including its documented
    -- exception for dimensionless/fold-change values.
    CONSTRAINT parameter_value_non_negative_unless_dimensionless CHECK (
        NOT (value < 0 AND unit NOT IN ('dimensionless', 'fold_change'))
    ),

    -- DEVIATION 1 (stricter than Python, deliberate): reject the empty
    -- string as a canonical id. Python's ParameterBuilder already requires a
    -- truthy subject, but direct Parameter(...) construction does not.
    -- Needed because the NULL-safe scope columns below use a sentinel, and
    -- an empty-string id must not be confusable with "no scope". Verified
    -- by experiment that ParameterStore treats object_canonical_id='' as a
    -- DISTINCT key from None, so the sentinel (0x1F) was chosen instead of
    -- COALESCE(...,'') precisely to preserve that distinction. The
    -- equivalence suite must confirm no existing test or seed uses ''.
    CONSTRAINT parameter_ids_not_empty CHECK (
        length(subject_canonical_id) > 0
        AND (object_canonical_id IS NULL OR object_canonical_id <> E'\x1F')
        AND (reaction_canonical_id IS NULL OR reaction_canonical_id <> E'\x1F')
        AND subject_canonical_id <> E'\x1F'
    ),

    -- NULL-safe scope columns. ParameterStore's index key is the 4-tuple
    -- (subject, kind, object, reaction) where object/reaction are usually
    -- None, and query() matches them EXACTLY — passing nothing matches only
    -- subject-only parameters, it is NOT a wildcard. In SQL,
    -- `object_canonical_id = NULL` never matches, so a direct translation
    -- would silently return zero rows for the common case. These generated
    -- columns give NULL-safe equality that a btree index can serve.
    -- The sentinel MUST match the one used in 002_query_api.sql.
    object_scope            text GENERATED ALWAYS AS (
                                coalesce(object_canonical_id, E'\x1F')
                            ) STORED,
    reaction_scope          text GENERATED ALWAYS AS (
                                coalesce(reaction_canonical_id, E'\x1F')
                            ) STORED
);

COMMENT ON TABLE parameter IS
    'Append-only kinetic parameter facts. Never UPDATE or DELETE a row here: '
    'a better measurement is a new row plus a parameter_supersession row. '
    'Enforced by grants and triggers in 003_append_only_enforcement.sql.';
COMMENT ON COLUMN parameter.uncertainty_type IS
    'Documented values: 95ci, sd, sem, range, none. Deliberately NOT '
    'constrained by CHECK — the Python ParameterValue accepts any string, and '
    'a stricter DB would reject writes the in-memory store accepts, producing '
    'dual-write divergence caused by the migration itself.';

-- The exact ParameterStore._subject_kind_index key, plus the ordering
-- tiebreaker, so parameter_query() is a single index range scan.
CREATE INDEX IF NOT EXISTS parameter_scope_key_idx
    ON parameter (subject_canonical_id, kind, object_scope, reaction_scope, seq);
-- get_all_for_subject() scans by subject across all kinds.
CREATE INDEX IF NOT EXISTS parameter_subject_idx
    ON parameter (subject_canonical_id, kind, seq);
-- "every parameter from this paper" — the citation-audit access pattern.
CREATE INDEX IF NOT EXISTS parameter_source_idx
    ON parameter (source_id);
CREATE INDEX IF NOT EXISTS parameter_supersedes_idx
    ON parameter (supersedes_parameter_id)
    WHERE supersedes_parameter_id IS NOT NULL;


-- ---------------------------------------------------------------------
-- parameter_supersession — the append-only replacement for the mutable
-- Parameter.is_superseded flag.
--
-- One row means "this parameter is no longer the working value". The
-- PRIMARY KEY makes it idempotent, matching the fact that the Python flag
-- is a boolean (marking twice is a no-op). The successor is NOT recorded
-- here: Python permits two different new parameters to supersede the same
-- old one, and the authoritative per-successor link is
-- parameter.supersedes_parameter_id. Recording a single "superseded_by"
-- here would invent a first-wins fact the in-memory store does not have.
-- ---------------------------------------------------------------------

CREATE TABLE IF NOT EXISTS parameter_supersession (
    superseded_parameter_id uuid PRIMARY KEY
                                 REFERENCES parameter (parameter_id),
    recorded_at             timestamptz NOT NULL DEFAULT now(),
    recorded_by             text NOT NULL DEFAULT '',
    note                    text NOT NULL DEFAULT ''
);


-- ---------------------------------------------------------------------
-- parameter_conflict — mirrors store.ConflictRecord.
--
-- Conflicts are not errors: they are genuine scientific disagreement,
-- flagged for expert curation and never silently resolved. Both
-- parameters stay live and queryable.
-- ---------------------------------------------------------------------

CREATE TABLE IF NOT EXISTS parameter_conflict (
    -- Python-assigned UUID again, for the same idempotency reason.
    conflict_id             uuid PRIMARY KEY,
    seq                     bigint GENERATED ALWAYS AS IDENTITY NOT NULL UNIQUE,
    existing_parameter_id   uuid NOT NULL REFERENCES parameter (parameter_id),
    incoming_parameter_id   uuid NOT NULL REFERENCES parameter (parameter_id),

    -- Denormalised copies, matching ConflictRecord's own fields exactly so
    -- the record round-trips without a join. They are historical snapshots:
    -- the values AT DETECTION TIME, which is what a curator needs to see.
    subject_canonical_id    text NOT NULL,
    kind                    text NOT NULL REFERENCES ref_parameter_kind (kind),
    existing_value          double precision NOT NULL,
    incoming_value          double precision NOT NULL,
    fold_difference         double precision NOT NULL,
    detected_at             timestamptz NOT NULL DEFAULT now(),

    CONSTRAINT parameter_conflict_distinct_parameters CHECK (
        existing_parameter_id <> incoming_parameter_id
    )
);

COMMENT ON TABLE parameter_conflict IS
    'FIDELITY NOTE: ConflictRecord in Python carries subject + kind but NOT '
    'object_canonical_id / reaction_canonical_id, so two conflicts on the same '
    'subject+kind for DIFFERENT objects are indistinguishable in the record '
    'itself. Rather than add columns the Python type does not have (which would '
    'change its to_dict() shape), the scope is exposed by joining to the '
    'referenced parameter rows in v_parameter_conflict. Whether ConflictRecord '
    'should gain those fields is a Python-side decision, flagged as an open '
    'question in infra/MIGRATION_PLAN.md.';

CREATE INDEX IF NOT EXISTS parameter_conflict_existing_idx
    ON parameter_conflict (existing_parameter_id);
CREATE INDEX IF NOT EXISTS parameter_conflict_incoming_idx
    ON parameter_conflict (incoming_parameter_id);
CREATE INDEX IF NOT EXISTS parameter_conflict_fold_idx
    ON parameter_conflict (fold_difference DESC, seq);


-- ---------------------------------------------------------------------
-- conflict_resolution — the append-only replacement for the mutable
-- ConflictRecord.resolved / .resolution_notes pair.
--
-- ParameterStore.resolve_conflict() can be called more than once for the
-- same conflict with different notes, and the last call wins. That is
-- reproduced by taking the highest resolution_seq per conflict in the
-- view, rather than by overwriting a row.
-- ---------------------------------------------------------------------

CREATE TABLE IF NOT EXISTS conflict_resolution (
    resolution_seq   bigint GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    conflict_id      uuid NOT NULL REFERENCES parameter_conflict (conflict_id),
    resolution_notes text NOT NULL,
    resolved_by      text NOT NULL DEFAULT '',
    resolved_at      timestamptz NOT NULL DEFAULT now()
);

CREATE INDEX IF NOT EXISTS conflict_resolution_conflict_idx
    ON conflict_resolution (conflict_id, resolution_seq DESC);


-- ---------------------------------------------------------------------
-- migration_audit — bookkeeping for the dual-write rollout itself.
-- Records which backfill runs happened, so a re-run can be proven to be
-- a no-op rather than assumed to be one.
-- ---------------------------------------------------------------------

CREATE TABLE IF NOT EXISTS migration_audit (
    audit_seq        bigint GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    run_label        text NOT NULL,
    source_snapshot  text NOT NULL,
    parameters_seen  integer NOT NULL DEFAULT 0,
    parameters_new   integer NOT NULL DEFAULT 0,
    conflicts_seen   integer NOT NULL DEFAULT 0,
    conflicts_new    integer NOT NULL DEFAULT 0,
    started_at       timestamptz NOT NULL DEFAULT now(),
    finished_at      timestamptz,
    notes            text NOT NULL DEFAULT ''
);

COMMENT ON TABLE migration_audit IS
    'source_snapshot is the path/hash of the ParameterStore JSON snapshot the '
    'backfill read. Deliberately NOT the live store object: the backfill reads '
    'a snapshot and never writes back to the file-based store, so the existing '
    'store cannot be damaged by a migration run.';
