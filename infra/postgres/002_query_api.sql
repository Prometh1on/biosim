-- =====================================================================
-- Biosim parameter store — query API.
--
-- Purpose: keep the ranking rules, the NULL-safe scope matching and the
-- superseded/resolved derivations in ONE place, so the Python adapter, an
-- ad-hoc psql session and any future service all get byte-identical
-- semantics. During shadow-read this is what gets compared against
-- ParameterStore.query()/best()/stats(); if the logic were duplicated in
-- the Python adapter instead, a divergence could be an adapter bug rather
-- than a real data difference and would be much harder to localise.
--
-- Idempotent: CREATE OR REPLACE throughout.
-- =====================================================================

SET search_path = biosim, public;


-- ---------------------------------------------------------------------
-- v_parameter — one row per stored parameter, with is_superseded DERIVED
-- from the append-only parameter_supersession table rather than stored as
-- a mutable flag. Columns are listed explicitly (not p.*) so that adding a
-- column to the table cannot silently change this view's row type, which
-- every function below returns.
-- ---------------------------------------------------------------------

CREATE OR REPLACE VIEW v_parameter AS
SELECT
    p.parameter_id,
    p.seq,
    p.kind,
    p.custom_kind_name,
    p.value,
    p.unit,
    p.uncertainty_lower,
    p.uncertainty_upper,
    p.uncertainty_type,
    p.subject_canonical_id,
    p.object_canonical_id,
    p.reaction_canonical_id,
    p.object_scope,
    p.reaction_scope,
    p.confidence_tier,
    p.source_id,
    src.pubmed_id,
    src.doi,
    src.arxiv_id,
    src.source_description,
    p.extracted_by,
    p.extracted_at,
    p.figure_or_table,
    p.is_primary_source,
    p.organism,
    p.experimental_system,
    p.cell_type,
    p.tissue,
    p.temperature_celsius,
    p.ph,
    p.context_notes,
    p.curator_notes,
    p.created_at,
    p.supersedes_parameter_id,
    (sup.superseded_parameter_id IS NOT NULL) AS is_superseded
FROM parameter p
JOIN literature_source src
  ON src.source_id = p.source_id
LEFT JOIN parameter_supersession sup
  ON sup.superseded_parameter_id = p.parameter_id;

COMMENT ON VIEW v_parameter IS
    'Reconstructs a biosim_parameter_db.types.Parameter one-for-one, including '
    'the derived is_superseded flag. This is the only shape the adapter should '
    'read; it must never read the parameter table directly, or it will lose the '
    'supersession derivation.';


-- ---------------------------------------------------------------------
-- parameter_query — exact translation of ParameterStore.query().
--
-- Three behaviours here are easy to get wrong and were verified against
-- the running Python store rather than read off the source:
--
--  (a) OBJECT/REACTION SCOPE IS EXACT-MATCH, NOT A WILDCARD. Omitting
--      p_object_canonical_id matches ONLY subject-only parameters. A Kd
--      stored .for_object(insulin_id) is invisible unless the object is
--      passed explicitly. This is deliberate in the Python store: silently
--      mixing a subject's Kd for insulin with its Kd for IGF-1 would treat
--      two distinct physical quantities as rival measurements of one.
--      Implemented via the sentinel-backed *_scope columns because plain
--      `= NULL` matches nothing in SQL. THE SENTINEL BELOW MUST MATCH THE
--      GENERATED COLUMN EXPRESSION IN 001_schema.sql.
--      (`object_canonical_id IS NOT DISTINCT FROM p_object_canonical_id`
--      is the semantically obvious alternative but cannot use the index.)
--
--  (b) TIE ORDER IS INSERTION ORDER. Python sorts with a stable sort, so
--      equal-priority parameters keep insertion order and best() returns
--      the earliest-inserted. `ORDER BY ... , seq ASC` reproduces that.
--      Dropping the seq tiebreaker would make best() non-deterministic
--      among ties, which propagates through ParameterLoader into
--      simulated results.
--
--  (c) min_confidence FILTERS BY RANK, NOT BY EQUALITY.
-- ---------------------------------------------------------------------

CREATE OR REPLACE FUNCTION parameter_query(
    p_subject_canonical_id  text,
    p_kind                  text,
    p_min_confidence        text    DEFAULT 'bronze',
    p_organism              text    DEFAULT NULL,
    p_include_superseded    boolean DEFAULT false,
    p_object_canonical_id   text    DEFAULT NULL,
    p_reaction_canonical_id text    DEFAULT NULL
)
RETURNS SETOF v_parameter
LANGUAGE plpgsql
STABLE
AS $$
DECLARE
    v_min_rank smallint;
BEGIN
    -- Fail loud on an unknown tier/kind/organism, mirroring the Python
    -- store's KeyError on CONFIDENCE_RANK[...] and the enum's ValueError.
    -- Returning zero rows instead would look like "no such parameter",
    -- which is exactly the kind of silent-wrong-answer this project's own
    -- bug history keeps warning about.
    SELECT rank INTO v_min_rank
      FROM ref_confidence_tier WHERE tier = p_min_confidence;
    IF v_min_rank IS NULL THEN
        RAISE EXCEPTION 'Unknown confidence tier: %', p_min_confidence
            USING HINT = 'Expected one of gold, silver, bronze.';
    END IF;

    IF NOT EXISTS (SELECT 1 FROM ref_parameter_kind WHERE kind = p_kind) THEN
        RAISE EXCEPTION 'Unknown parameter kind: %', p_kind;
    END IF;

    IF p_organism IS NOT NULL
       AND NOT EXISTS (SELECT 1 FROM ref_organism WHERE organism = p_organism) THEN
        RAISE EXCEPTION 'Unknown organism: %', p_organism;
    END IF;

    RETURN QUERY
    SELECT v.*
    FROM v_parameter v
    JOIN ref_confidence_tier     ct ON ct.tier     = v.confidence_tier
    JOIN ref_organism            o  ON o.organism  = v.organism
    JOIN ref_experimental_system es ON es.system   = v.experimental_system
    WHERE v.subject_canonical_id = p_subject_canonical_id
      AND v.kind                 = p_kind
      -- (a) NULL-safe exact scope match. Sentinel = 001_schema.sql's.
      AND v.object_scope         = coalesce(p_object_canonical_id,   E'\x1F')
      AND v.reaction_scope       = coalesce(p_reaction_canonical_id, E'\x1F')
      AND (p_include_superseded OR NOT v.is_superseded)
      AND ct.rank >= v_min_rank                                        -- (c)
      AND (p_organism IS NULL OR v.organism = p_organism)
    ORDER BY ct.rank DESC, o.rank DESC, es.rank DESC, v.seq ASC;       -- (b)
END;
$$;


-- ---------------------------------------------------------------------
-- parameter_best — ParameterStore.best(). The primary lookup used by the
-- ODE engine via biosim_seeds.loader.ParameterLoader, so this is the one
-- function whose output can change a simulation result. Returns zero rows
-- where Python returns None.
-- ---------------------------------------------------------------------

CREATE OR REPLACE FUNCTION parameter_best(
    p_subject_canonical_id  text,
    p_kind                  text,
    p_min_confidence        text DEFAULT 'bronze',
    p_organism              text DEFAULT NULL,
    p_object_canonical_id   text DEFAULT NULL,
    p_reaction_canonical_id text DEFAULT NULL
)
RETURNS SETOF v_parameter
LANGUAGE sql
STABLE
AS $$
    SELECT *
    FROM parameter_query(
        p_subject_canonical_id  => p_subject_canonical_id,
        p_kind                  => p_kind,
        p_min_confidence        => p_min_confidence,
        p_organism              => p_organism,
        p_include_superseded    => false,
        p_object_canonical_id   => p_object_canonical_id,
        p_reaction_canonical_id => p_reaction_canonical_id
    )
    LIMIT 1;
$$;


-- ---------------------------------------------------------------------
-- parameter_all_for_subject — the SQL counterpart of
-- ParameterStore.get_all_for_subject().
--
-- RESOLVED (open question 5) — the divergence documented here was closed
-- by fixing PYTHON, not this function. This function's body was always
-- correct (its ORDER BY already did the right thing); the Python method
-- buckled by KIND ONLY with a plain assignment
-- (`result[kind_str] = sorted(params)`) inside a loop over the 4-tuple
-- (subject, kind, object, reaction) index, so a later index key
-- OVERWROTE the earlier bucket instead of accumulating into it --
-- reproduced directly (three Kd parameters for one subject, one
-- returned) and fixed in store.py's get_all_for_subject() by
-- accumulating across scopes, then sorting once per kind with an
-- explicit ascending insertion-sequence tiebreak
-- (packages/parameter-db/tests/test_parameter_db.py::TestObjectScopedParameters,
-- the *_returns_every_scope_not_just_the_last and
-- *_ties_break_by_global_insertion_order regression tests).
--
-- The two are now equivalent MODULO SHAPE, not byte-for-byte: this
-- function returns a flat SETOF v_parameter ordered
-- `v.kind, ct.rank DESC, o.rank DESC, es.rank DESC, v.seq ASC`; Python
-- returns `dict[kind -> list[Parameter]]`. A shadow-read adapter must
-- GROUP this function's rows by kind while PRESERVING row order, and
-- must not re-sort client-side -- Python's tiebreak is now an explicit
-- counter that `parameter.seq` (GENERATED ALWAYS AS IDENTITY) is meant
-- to mirror exactly, so a backfill that does not insert in snapshot
-- array order will silently diverge from this invariant again (see
-- MIGRATION_PLAN.md's backfill-ordering note, added alongside this
-- resolution).
--
-- Zero production callers of get_all_for_subject() existed at the time
-- of this fix (grepped the whole repo; only its own test called it), so
-- there was no live blast radius to check beyond the new regression
-- tests themselves.
-- ---------------------------------------------------------------------

CREATE OR REPLACE FUNCTION parameter_all_for_subject(
    p_subject_canonical_id text,
    p_include_superseded   boolean DEFAULT false
)
RETURNS SETOF v_parameter
LANGUAGE sql
STABLE
AS $$
    SELECT v.*
    FROM v_parameter v
    JOIN ref_confidence_tier     ct ON ct.tier    = v.confidence_tier
    JOIN ref_organism            o  ON o.organism = v.organism
    JOIN ref_experimental_system es ON es.system  = v.experimental_system
    WHERE v.subject_canonical_id = p_subject_canonical_id
      AND (p_include_superseded OR NOT v.is_superseded)
    ORDER BY v.kind, ct.rank DESC, o.rank DESC, es.rank DESC, v.seq ASC;
$$;


-- ---------------------------------------------------------------------
-- v_parameter_conflict — ConflictRecord with resolved/resolution_notes
-- derived from the latest conflict_resolution row (last call wins, as in
-- ParameterStore.resolve_conflict()). Also exposes the object/reaction
-- scope of both sides by join, which ConflictRecord itself cannot carry.
-- ---------------------------------------------------------------------

CREATE OR REPLACE VIEW v_parameter_conflict AS
SELECT
    c.conflict_id,
    c.seq,
    c.existing_parameter_id,
    c.incoming_parameter_id,
    c.subject_canonical_id,
    c.kind,
    c.existing_value,
    c.incoming_value,
    c.fold_difference,
    c.detected_at,
    (r.resolution_seq IS NOT NULL)          AS resolved,
    coalesce(r.resolution_notes, '')        AS resolution_notes,
    r.resolved_by,
    r.resolved_at,
    -- Scope of the conflicting pair, derived rather than duplicated.
    pe.object_canonical_id                  AS existing_object_canonical_id,
    pe.reaction_canonical_id                AS existing_reaction_canonical_id,
    pi.object_canonical_id                  AS incoming_object_canonical_id,
    pi.reaction_canonical_id                AS incoming_reaction_canonical_id,
    (pe.object_scope   = pi.object_scope
     AND pe.reaction_scope = pi.reaction_scope) AS same_scope
FROM parameter_conflict c
JOIN parameter pe ON pe.parameter_id = c.existing_parameter_id
JOIN parameter pi ON pi.parameter_id = c.incoming_parameter_id
LEFT JOIN LATERAL (
    SELECT cr.resolution_seq, cr.resolution_notes, cr.resolved_by, cr.resolved_at
    FROM conflict_resolution cr
    WHERE cr.conflict_id = c.conflict_id
    ORDER BY cr.resolution_seq DESC
    LIMIT 1
) r ON true;

COMMENT ON COLUMN v_parameter_conflict.same_scope IS
    'A conflict where same_scope = false would mean two parameters with '
    'different object/reaction scopes were compared as rival measurements of '
    'the same quantity. The current store cannot produce one (its index key '
    'includes the scope), so this is a standing invariant check: any row with '
    'same_scope = false indicates the scoping disambiguation has regressed.';


-- ---------------------------------------------------------------------
-- v_unresolved_conflict — ParameterStore.unresolved_conflicts(), largest
-- fold difference first. Python has no tiebreaker beyond dict insertion
-- order, so seq is used to make the order deterministic.
-- ---------------------------------------------------------------------

CREATE OR REPLACE VIEW v_unresolved_conflict AS
SELECT *
FROM v_parameter_conflict
WHERE NOT resolved
ORDER BY fold_difference DESC, seq ASC;


-- ---------------------------------------------------------------------
-- parameter_stats — ParameterStore.stats(), same keys, same shape.
-- by_confidence_tier / by_parameter_kind count ACTIVE parameters only,
-- matching the Python loop's `if p.is_superseded: continue`, and omit
-- zero-count keys, matching its increment-on-demand dict.
-- ---------------------------------------------------------------------

CREATE OR REPLACE FUNCTION parameter_stats()
RETURNS jsonb
LANGUAGE sql
STABLE
AS $$
    SELECT jsonb_build_object(
        'total_parameters',
            (SELECT count(*) FROM parameter),
        'active_parameters',
            (SELECT count(*) FROM v_parameter WHERE NOT is_superseded),
        'superseded_parameters',
            (SELECT count(*) FROM parameter_supersession),
        'unresolved_conflicts',
            (SELECT count(*) FROM v_parameter_conflict WHERE NOT resolved),
        'total_conflicts',
            (SELECT count(*) FROM parameter_conflict),
        'by_confidence_tier', coalesce(
            (SELECT jsonb_object_agg(confidence_tier, n) FROM (
                SELECT confidence_tier, count(*) AS n
                FROM v_parameter WHERE NOT is_superseded
                GROUP BY confidence_tier
            ) t), '{}'::jsonb),
        'by_parameter_kind', coalesce(
            (SELECT jsonb_object_agg(kind, n) FROM (
                SELECT kind, count(*) AS n
                FROM v_parameter WHERE NOT is_superseded
                GROUP BY kind
            ) t), '{}'::jsonb)
    );
$$;


-- ---------------------------------------------------------------------
-- Conflict detection, as a read-only helper rather than a trigger.
--
-- Why NOT a trigger: ParameterStore.load() deliberately restores conflicts
-- from the snapshot WITHOUT re-running detection. A BEFORE/AFTER INSERT
-- trigger that detected conflicts would fire during the backfill and mint
-- brand-new conflict_ids alongside the migrated ones, double-counting every
-- historical conflict and breaking stats() equivalence on the very first
-- shadow-read comparison. Detection therefore stays in the application
-- layer (where it already lives), and this function exists only so a
-- reviewer can ask "what WOULD be detected for this row" without writing.
-- ---------------------------------------------------------------------

CREATE OR REPLACE FUNCTION would_conflict_with(
    p_parameter_id uuid,
    p_fold_threshold double precision DEFAULT 2.0
)
RETURNS TABLE (
    other_parameter_id uuid,
    other_value        double precision,
    this_value         double precision,
    fold_difference    double precision
)
LANGUAGE sql
STABLE
AS $$
    WITH me AS (
        SELECT * FROM v_parameter WHERE parameter_id = p_parameter_id
    )
    SELECT
        o.parameter_id,
        o.value,
        me.value,
        greatest(o.value / me.value, me.value / o.value)
    FROM me
    JOIN v_parameter o
      ON  o.subject_canonical_id = me.subject_canonical_id
      AND o.kind                 = me.kind
      AND o.object_scope         = me.object_scope
      AND o.reaction_scope       = me.reaction_scope
      AND o.parameter_id        <> me.parameter_id
    WHERE NOT o.is_superseded
      -- Mirrors _detect_conflict()'s zero guard, which skips the pair
      -- entirely rather than dividing by zero.
      AND o.value <> 0 AND me.value <> 0
      AND greatest(o.value / me.value, me.value / o.value) >= p_fold_threshold;
$$;
