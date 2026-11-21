# Biosim local persistence layer (Track 3 scoping)

**This is scoping infrastructure, not yet wired into any application
code.** Bringing these containers up or down has zero effect on the
existing test suite or on any simulation result. See
`infra/MIGRATION_PLAN.md` for the full design, open questions, and the
dual-write/shadow-read rollout plan a real implementation pass would
follow.

## Quick start

**Precondition: Docker Desktop must actually be running first** (not
just installed) -- `docker compose up` fails with a connection error to
the Docker API otherwise. Check with:

```powershell
docker info               # fails with "cannot connect to the Docker API" if Desktop isn't running
docker compose ps         # once it's running, safe to call any time
```

```powershell
docker compose up -d
docker compose ps        # both services should report "healthy"
```

First boot applies `infra/postgres/*.sql` automatically (via Postgres's
own `/docker-entrypoint-initdb.d` mechanism) and `infra/neo4j/001_constraints.cypher`
via the one-shot `neo4j-init` service. Every script is idempotent, so
re-running `docker compose up -d` against an existing data volume is
safe.

**`/docker-entrypoint-initdb.d` scripts run ONLY ONCE, on an empty
volume.** Editing a `.sql`/`.cypher` file after the first boot has no
effect on an already-initialized container -- re-apply it by hand, or
wipe and reinitialize:

```powershell
# Re-apply one changed SQL file by hand (Postgres container already running):
docker exec -i biosim-postgres psql -v ON_ERROR_STOP=1 -U biosim -d biosim -f /docker-entrypoint-initdb.d/002_query_api.sql

# Or wipe local data and let initdb.d re-run everything from scratch:
docker compose down -v
docker compose up -d
```

```powershell
docker compose down          # stop, keep data
docker compose down -v       # stop and delete all local data
```

## Credentials

Both services use throwaway local-development defaults
(`biosim_local_dev_only`), documented as such directly in
`docker-compose.yml`. These are not secrets — override them by exporting
`BIOSIM_PG_PASSWORD` / `BIOSIM_NEO4J_PASSWORD` before running
`docker compose up`, or by copying `infra/.env.example` to a `.env` file
**at the repo root** (`<repo root>/.env`, NOT `infra/.env` --
docker-compose reads `.env` from the directory containing the compose
file, which is the repo root, since `docker-compose.yml` lives there)
and editing it. `.env`/`.env.*` are gitignored (with `.env.example`
explicitly un-ignored) specifically so creating this file is safe --
never commit a populated `.env` file regardless.

## Connecting manually

```powershell
# Postgres
docker exec -it biosim-postgres psql -U biosim -d biosim

# Neo4j Browser: http://127.0.0.1:7474 (user: neo4j)
# Or cypher-shell:
docker exec -it biosim-neo4j cypher-shell -u neo4j -p biosim_local_dev_only
```

## Verifying the append-only guarantee

`infra/postgres/003_append_only_enforcement.sql` runs a self-test on
every application (including the automatic first-boot run) and raises
loudly if the guard triggers aren't actually installed. Confirm it ran
cleanly:

```powershell
docker compose logs postgres | Select-String "append-only"
```

Expect: `Biosim append-only enforcement: OK (guards installed, UPDATE/DELETE/TRUNCATE revoked from biosim_writer).`

## Cost and scope

Local-only, zero-cost, self-hosted — both images are pulled from Docker
Hub and run entirely on this machine. No account, no API key, no egress,
no managed/cloud service of any kind (matches CLAUDE.md's STRATEGIC
ROADMAP decision 2: stay fully free/local). Nothing here is a
requirement for running Biosim's existing test suite or simulations —
it's exclusively for the Track 3 migration work described in
`infra/MIGRATION_PLAN.md`.
