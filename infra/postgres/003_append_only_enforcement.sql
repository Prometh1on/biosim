-- =====================================================================
-- Biosim parameter store — append-only enforcement.
--
-- Answers the design question directly: YES, append-only should be a
-- database-level guarantee, not an application convention. Rule 5 in
-- CLAUDE.md ("never delete or overwrite parameters") is currently held up
-- by nothing but discipline; a single stray UPDATE from a future script or
-- a psql session would silently rewrite a literature-anchored value and
-- change simulation results with no trace. The schema in 001_schema.sql was
-- shaped specifically so that this is enforceable: because supersession and
-- conflict resolution are append-only EVENT rows rather than mutable flags,
-- no legitimate code path needs UPDATE or DELETE on any fact table at all.
--
-- Two layers, deliberately:
--   1. PRIVILEGES  — the application role is granted INSERT and SELECT and
--      nothing else. This is the primary control.
--   2. TRIGGERS    — a hard rejection of UPDATE/DELETE/TRUNCATE that also
--      applies to the table owner and to superusers, who bypass privilege
--      checks entirely. This is defence in depth, and it is honest about
--      its limit: the owner can ALTER TABLE ... DISABLE TRIGGER. That
--      friction is the point — a genuine data repair should require a
--      deliberate, visible, auditable act.
--
-- Idempotent: safe to re-apply.
--
-- CREDENTIALS: no passwords are set here and none are committed anywhere in
-- this repository. Both roles are NOLOGIN. The application connects with the
-- ordinary container user and immediately issues `SET ROLE biosim_writer;`,
-- after which privilege checks apply to the writer role (superuser bypass
-- stops applying once current_user is a non-superuser role). Stated plainly:
-- this is an ACCIDENT GUARD, not a security boundary — the session can
-- RESET ROLE. Anything beyond local development must instead give the
-- application its own LOGIN role whose password comes from the environment
-- or a secret manager at deploy time. See infra/README.md.
-- =====================================================================

SET search_path = biosim, public;


-- ---------------------------------------------------------------------
-- Roles
-- ---------------------------------------------------------------------

DO $$
BEGIN
    IF NOT EXISTS (SELECT 1 FROM pg_roles WHERE rolname = 'biosim_reader') THEN
        CREATE ROLE biosim_reader NOLOGIN;
    END IF;
    IF NOT EXISTS (SELECT 1 FROM pg_roles WHERE rolname = 'biosim_writer') THEN
        CREATE ROLE biosim_writer NOLOGIN;
    END IF;
END
$$;

COMMENT ON ROLE biosim_reader IS
    'Read-only. Use for shadow-read comparison, reporting and any analysis '
    'session. Cannot INSERT, so a comparison run can never mutate the store '
    'it is comparing against.';
COMMENT ON ROLE biosim_writer IS
    'Append-only writer: INSERT + SELECT. Deliberately has no UPDATE, DELETE '
    'or TRUNCATE on any fact table.';

GRANT USAGE ON SCHEMA biosim TO biosim_reader, biosim_writer;


-- ---------------------------------------------------------------------
-- Privileges
-- ---------------------------------------------------------------------

-- Reference tables: readable by everyone, writable by the owner only.
GRANT SELECT ON
    ref_parameter_kind,
    ref_confidence_tier,
    ref_organism,
    ref_experimental_system
TO biosim_reader, biosim_writer;

-- Fact tables: INSERT + SELECT for the writer, SELECT for the reader.
GRANT SELECT ON
    literature_source,
    parameter,
    parameter_supersession,
    parameter_conflict,
    conflict_resolution,
    migration_audit
TO biosim_reader, biosim_writer;

GRANT INSERT ON
    literature_source,
    parameter,
    parameter_supersession,
    parameter_conflict,
    conflict_resolution,
    migration_audit
TO biosim_writer;

-- migration_audit is the single exception that needs UPDATE: a backfill run
-- opens a row, then stamps finished_at/counts when it completes. Granted
-- narrowly and column-scoped rather than table-wide, and it is bookkeeping
-- about the migration, never biological data.
GRANT UPDATE (parameters_seen, parameters_new, conflicts_seen, conflicts_new,
              finished_at, notes)
    ON migration_audit TO biosim_writer;

-- Views and query functions.
GRANT SELECT ON v_parameter, v_parameter_conflict, v_unresolved_conflict
    TO biosim_reader, biosim_writer;
GRANT EXECUTE ON FUNCTION
        parameter_query(text, text, text, text, boolean, text, text),
        parameter_best(text, text, text, text, text, text),
        parameter_all_for_subject(text, boolean),
        parameter_stats(),
        would_conflict_with(uuid, double precision)
    TO biosim_reader, biosim_writer;

-- Belt and braces: revoke the destructive verbs explicitly, so that a
-- future `GRANT ALL` in some convenience script does not quietly hand them
-- back without anyone noticing.
REVOKE UPDATE, DELETE, TRUNCATE ON
    literature_source,
    parameter,
    parameter_supersession,
    parameter_conflict,
    conflict_resolution
FROM biosim_reader, biosim_writer;

REVOKE INSERT, UPDATE, DELETE, TRUNCATE ON
    ref_parameter_kind,
    ref_confidence_tier,
    ref_organism,
    ref_experimental_system
FROM biosim_reader, biosim_writer;

-- Anything created later in this schema starts read-only for the reader
-- and append-only for the writer, instead of inheriting whatever the
-- creating session happened to leave behind.
ALTER DEFAULT PRIVILEGES IN SCHEMA biosim
    GRANT SELECT ON TABLES TO biosim_reader, biosim_writer;


-- ---------------------------------------------------------------------
-- Trigger-level rejection (applies to the owner and superusers too)
-- ---------------------------------------------------------------------

CREATE OR REPLACE FUNCTION forbid_mutation()
RETURNS trigger
LANGUAGE plpgsql
AS $$
BEGIN
    RAISE EXCEPTION
        'Table %.% is append-only (Biosim rule 5): % is not permitted.',
        TG_TABLE_SCHEMA, TG_TABLE_NAME, TG_OP
        USING HINT =
            'A better measurement is a NEW parameter row plus a '
            'parameter_supersession row. A resolved conflict is a NEW '
            'conflict_resolution row. If this really is a data repair, an '
            'owner must ALTER TABLE ... DISABLE TRIGGER deliberately and '
            'record why in migration_audit.';
END;
$$;

COMMENT ON FUNCTION forbid_mutation() IS
    'Shared BEFORE UPDATE/DELETE/TRUNCATE guard for the append-only fact '
    'tables. Row-level for UPDATE/DELETE, statement-level for TRUNCATE '
    '(TRUNCATE has no row-level trigger form).';

DO $$
DECLARE
    t text;
BEGIN
    FOREACH t IN ARRAY ARRAY[
        'literature_source',
        'parameter',
        'parameter_supersession',
        'parameter_conflict',
        'conflict_resolution'
    ]
    LOOP
        EXECUTE format(
            'DROP TRIGGER IF EXISTS %I ON biosim.%I',
            t || '_no_mutation', t);
        EXECUTE format(
            'CREATE TRIGGER %I BEFORE UPDATE OR DELETE ON biosim.%I '
            'FOR EACH ROW EXECUTE FUNCTION biosim.forbid_mutation()',
            t || '_no_mutation', t);

        EXECUTE format(
            'DROP TRIGGER IF EXISTS %I ON biosim.%I',
            t || '_no_truncate', t);
        EXECUTE format(
            'CREATE TRIGGER %I BEFORE TRUNCATE ON biosim.%I '
            'FOR EACH STATEMENT EXECUTE FUNCTION biosim.forbid_mutation()',
            t || '_no_truncate', t);
    END LOOP;
END
$$;


-- ---------------------------------------------------------------------
-- Self-test. Runs on every application of this file and fails loudly if
-- the guarantee is not actually in place, rather than assuming it is
-- because the DDL above executed without error.
-- ---------------------------------------------------------------------

DO $$
DECLARE
    v_blocked boolean := false;
BEGIN
    BEGIN
        -- No rows exist yet on first init, and an UPDATE matching zero rows
        -- does not fire a row-level trigger, so assert against the catalog
        -- instead of attempting a real mutation.
        SELECT count(*) = 10 INTO v_blocked
        FROM pg_trigger tg
        JOIN pg_class c ON c.oid = tg.tgrelid
        JOIN pg_namespace n ON n.oid = c.relnamespace
        WHERE n.nspname = 'biosim'
          AND NOT tg.tgisinternal
          AND tg.tgname LIKE '%\_no\_%';
    END;

    IF NOT v_blocked THEN
        RAISE EXCEPTION
            'Append-only enforcement self-test FAILED: expected 10 guard '
            'triggers (5 tables x UPDATE/DELETE + TRUNCATE) in schema biosim.';
    END IF;

    RAISE NOTICE 'Biosim append-only enforcement: OK (guards installed, '
                 'UPDATE/DELETE/TRUNCATE revoked from biosim_writer).';
END
$$;
