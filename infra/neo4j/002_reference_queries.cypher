// =====================================================================
// Biosim ID registry — reference Cypher for every IDRegistry public method,
// plus the invariant checks that stand in for the constraints Neo4j
// Community Edition cannot declare.
//
// THIS FILE IS NOT AUTO-APPLIED. docker-compose's neo4j-init service runs
// only 001_constraints.cypher. The queries here take $parameters and are
// the specification the Python adapter implements and the equivalence suite
// asserts against — keeping them in one reviewable place is the same
// reasoning as putting the Postgres ranking rules in SQL functions rather
// than in the adapter.
//
// All reads are parameterised. Never build a Cypher string by
// concatenating an accession or a name into it.
// =====================================================================


// ---------------------------------------------------------------------
// current_canonical_id — the primitive every read below is built on.
// Follows the SUPERSEDED_BY chain to its end. Fixes bug 1 from
// 001_constraints.cypher: the in-memory _deprecated_index is one hop deep,
// so get_record() on a twice-merged id hands the caller another dead id.
// ---------------------------------------------------------------------
// :param canonical_id => 'BIOSIM:METABOLITE:000001'
MATCH (start:Entity {canonical_id: $canonical_id})
OPTIONAL MATCH path = (start)-[:SUPERSEDED_BY*1..]->(end_entity:Entity)
WHERE NOT (end_entity)-[:SUPERSEDED_BY]->()
RETURN
    start.canonical_id                                AS requested,
    start.status                                      AS requested_status,
    coalesce(end_entity.canonical_id, start.canonical_id) AS current,
    length(path)                                      AS merge_hops;


// ---------------------------------------------------------------------
// resolve(ExternalID) -> CanonicalID
//
// Returns the DISTINCT set of current entities an accession reaches.
//   0 rows  -> ResolutionError
//   1 row   -> the canonical id
//   >1 rows -> AmbiguousResolutionError (a real data-integrity problem
//              requiring curation, never silently pick one)
//
// Because IDENTIFIES edges are immutable, an accession that two entities
// once claimed but which were later merged reaches ONE current entity, so
// no false ambiguity is introduced by keeping the history. A genuine
// third-party collision still reaches two, exactly as the in-memory
// registry's collision check intends.
// ---------------------------------------------------------------------
// :param source => 'uniprot'  :param accession => 'P06213'
MATCH (x:ExternalId {key: $source + ':' + toUpper($accession)})-[:IDENTIFIES]->(owner:Entity)
OPTIONAL MATCH (owner)-[:SUPERSEDED_BY*1..]->(cur:Entity)
WHERE NOT (cur)-[:SUPERSEDED_BY]->()
WITH DISTINCT coalesce(cur, owner) AS resolved
RETURN resolved.canonical_id AS canonical_id, resolved.entity_type AS entity_type
ORDER BY canonical_id;


// ---------------------------------------------------------------------
// resolve_any(accession) -> list[CanonicalID]
// Source-agnostic discovery, ordered by RESOLUTION_PRIORITY (the
// OntologySource reference nodes), de-duplicated, first occurrence wins —
// matching the `seen` set in the Python implementation.
// ---------------------------------------------------------------------
// :param accession => 'C00031'
MATCH (x:ExternalId {accession: toUpper($accession)})-[:IDENTIFIES]->(owner:Entity)
MATCH (s:OntologySource {source: x.source})
OPTIONAL MATCH (owner)-[:SUPERSEDED_BY*1..]->(cur:Entity)
WHERE NOT (cur)-[:SUPERSEDED_BY]->()
WITH coalesce(cur, owner) AS resolved, s.priority AS priority
ORDER BY priority ASC, resolved.canonical_id ASC
RETURN collect(DISTINCT resolved.canonical_id) AS canonical_ids;


// ---------------------------------------------------------------------
// get_record(CanonicalID) -> ResolutionRecord
// One query returns the record, every cross-reference, and — for a
// deprecated id — the genuinely current id rather than the next hop.
// ---------------------------------------------------------------------
// :param canonical_id => 'BIOSIM:PROTEIN:000001'
MATCH (e:Entity {canonical_id: $canonical_id})
OPTIONAL MATCH (x:ExternalId)-[xr:IDENTIFIES]->(e)
OPTIONAL MATCH (old:Entity)-[:SUPERSEDED_BY]->(e)
OPTIONAL MATCH (e)-[:SUPERSEDED_BY*1..]->(cur:Entity)
WHERE NOT (cur)-[:SUPERSEDED_BY]->()
RETURN
    e.canonical_id  AS canonical_id,
    e.entity_type   AS entity_type,
    e.name          AS name,
    e.description   AS description,
    e.status        AS status,
    collect(DISTINCT {source: x.source, accession: x.accession, rank: xr.rank}) AS cross_references,
    collect(DISTINCT old.canonical_id)                           AS deprecated_canonical_ids,
    cur.canonical_id AS superseded_by_current;


// ---------------------------------------------------------------------
// cross_reference(CanonicalID, target_source) -> accession | null
//
// RESOLVED (was: MASTER_PLAN Track 3.4's "Neo4j invariant-4 decision").
// Python returns cross_references[source][0] -- the FIRST accession in
// insertion order. A real merge legitimately leaves an active entity
// with MULTIPLE accessions for one source (verified directly: merging
// entity A, UniProt P00001, into entity B, UniProt P99999, produces
// cross_references={'uniprot': ['P99999', 'P00001']} -- B's own
// restated accession first, A's inherited one second -- and
// cross_reference() deterministically returns 'P99999', not
// alphabetically-first 'P00001'). This is NOT an error case invariant 4
// (below) should reject; it's the routine, working shape of
// registry.py's own already-shipped merge-inheritance fix (see
// CLAUDE.md "IDRegistry's merge workflow left deprecated entities' data
// permanently stale"). The PRIOR version of this query used
// `ORDER BY accession ASC LIMIT 1` (alphabetically-first) specifically
// BECAUSE "insertion order is not a property of a graph" -- but that
// choice silently disagreed with Python's real answer the moment
// invariant 4's assumption (never more than one candidate) turned out
// to be false, which it provably is. Fixed by making insertion order an
// explicit, persisted graph property: the write path (Step 2 below) now
// stamps each IDENTIFIES edge with `rank` (0-based position in Python's
// own cross_references[source] list at the time it was written), and
// this query orders by that rank instead of the accession string.
// ---------------------------------------------------------------------
// :param canonical_id => 'BIOSIM:METABOLITE:000001'  :param source => 'kegg'
MATCH (x:ExternalId {source: $source})-[r:IDENTIFIES]->(e:Entity {canonical_id: $canonical_id})
RETURN x.accession AS accession
ORDER BY r.rank ASC
LIMIT 1;


// ---------------------------------------------------------------------
// search_by_name(name) -> list[CanonicalID]
//
// RESOLVED (was open question 4 in infra/MIGRATION_PLAN.md, now fixed at
// the Python source, not just scoped): the in-memory _name_index used to
// return a DEPRECATED canonical id after a merge (never cleaned on
// merge), and get_record() on that result then raised. Fixed in
// registry.py this session by repointing (not removing) the old name's
// index entry to the new canonical id at merge time -- functionally
// option (A) below, reached by repointing rather than filtering. The
// query below still implements (A) directly against the graph model
// (active-status entities only) and still returns status so option (B)
// -- return everything, let callers filter -- remains a one-line change
// if ever needed, but this is no longer an undecided behaviour choice.
// ---------------------------------------------------------------------
// :param name => 'insulin receptor'
MATCH (e:Entity {name_lower: toLower(trim($name))})
WHERE e.status = 'active'
RETURN e.canonical_id AS canonical_id, e.status AS status
ORDER BY canonical_id;


// ---------------------------------------------------------------------
// is_registered(ExternalID) -> bool
// ---------------------------------------------------------------------
// :param source => 'hmdb'  :param accession => 'HMDB0000122'
RETURN exists {
    MATCH (:ExternalId {key: $source + ':' + toUpper($accession)})
} AS is_registered;


// ---------------------------------------------------------------------
// stats() -> dict
// total_entities counts ACTIVE entities only, matching the in-memory
// version's post-fix behaviour (merged-away entities are popped from
// _canonical_index and no longer counted).
// ---------------------------------------------------------------------
CALL {
    MATCH (e:Entity) WHERE e.status = 'active'
    RETURN count(e) AS total_entities
}
CALL {
    MATCH (x:ExternalId) RETURN count(x) AS total_xrefs
}
CALL {
    MATCH (:Entity)-[r:SUPERSEDED_BY]->(:Entity) RETURN count(r) AS deprecated_mappings
}
CALL {
    MATCH (e:Entity) WHERE e.status = 'active'
    RETURN collect(DISTINCT [e.entity_type, 1]) AS _unused_type_placeholder
}
RETURN total_entities, total_xrefs, deprecated_mappings;
// by_entity_type / by_ontology_source are separate aggregations; kept apart
// so each is a plain, readable group-by rather than one unreadable query:
//   MATCH (e:Entity) WHERE e.status='active'
//   RETURN e.entity_type AS entity_type, count(*) AS n;
//   MATCH (x:ExternalId)-[:IDENTIFIES]->(e:Entity) WHERE e.status='active'
//   RETURN x.source AS source, count(DISTINCT e) AS n;


// =====================================================================
// WRITE PATH — register() / merge
// =====================================================================

// ---------------------------------------------------------------------
// Step 1: collision pre-check, the graph equivalent of register()'s
// "Cross-reference collision" ValueError. An accession collides only if
// its CURRENT owner is an entity other than the one being registered and
// other than any entity this call is deprecating. This reproduces the
// exemption that made the merge workflow functional in the first place,
// without needing to enumerate it as a special case: the chain traversal
// already resolves a to-be-deprecated owner to its successor.
// ---------------------------------------------------------------------
// :param accessions => [{source:'hmdb', accession:'HMDB0000122'}]
// :param deprecated => ['BIOSIM:METABOLITE:000001']
UNWIND $accessions AS xref
MATCH (x:ExternalId {key: xref.source + ':' + toUpper(xref.accession)})-[:IDENTIFIES]->(owner:Entity)
OPTIONAL MATCH (owner)-[:SUPERSEDED_BY*1..]->(cur:Entity)
WHERE NOT (cur)-[:SUPERSEDED_BY]->()
WITH xref, coalesce(cur, owner) AS current_owner
WHERE NOT current_owner.canonical_id IN $deprecated
RETURN xref.source AS source, xref.accession AS accession,
       collect(DISTINCT current_owner.canonical_id) AS blocking_owners;
// Any row returned => raise the collision error. Zero rows => proceed.


// ---------------------------------------------------------------------
// Step 2: create the entity and attach cross-references. Idempotent, so a
// re-run of the backfill neither duplicates nodes nor duplicates edges.
//
// `xref.rank` (Track 3.4 fix): 0-based position of this accession in
// Python's own cross_references[source] list at backfill time -- rank 0
// is Python's "primary" (types.py's ResolutionRecord.get_external_id()
// returns refs[0]). The backfill MUST pass this exactly, in the SAME
// order the source snapshot's list holds it, for cross_reference()'s
// query above to agree with Python's real answer. `ON CREATE SET`
// (not plain `SET`) matches this project's "IDENTIFIES edges are
// immutable by application convention" rule (MIGRATION_PLAN.md open
// question 3) -- a rank, once written, is never revised by a later
// backfill run.
// ---------------------------------------------------------------------
// :param canonical_id => 'BIOSIM:METABOLITE:000002'
// :param entity_type  => 'METABOLITE'   :param numeric_id => 2
// :param name => 'D-Glucose'  :param description => ''
// :param accessions => [{source:'hmdb', accession:'HMDB0000122', rank:0}]
MERGE (e:Entity {canonical_id: $canonical_id})
ON CREATE SET
    e.entity_type = $entity_type,
    e.numeric_id  = $numeric_id,
    e.created_at  = datetime()
SET e.name        = $name,
    e.name_lower  = toLower(trim($name)),
    e.description = $description,
    e.status      = coalesce(e.status, 'active')
WITH e
UNWIND $accessions AS xref
MERGE (x:ExternalId {key: xref.source + ':' + toUpper(xref.accession)})
ON CREATE SET x.source = xref.source, x.accession = toUpper(xref.accession)
MERGE (x)-[r:IDENTIFIES]->(e)
ON CREATE SET r.rank = xref.rank
RETURN e.canonical_id AS canonical_id, count(*) AS cross_references_attached;


// ---------------------------------------------------------------------
// Step 3: the merge/deprecation edge. Note what is NOT here: no
// re-pointing of IDENTIFIES edges, no name-index cleanup, no removal of
// the old entity. Those three omissions are the three bugs from
// 001_constraints.cypher's header, and they are omissions because the
// model makes them unnecessary, not because they were forgotten.
// ---------------------------------------------------------------------
// :param new_canonical_id => 'BIOSIM:METABOLITE:000002'
// :param deprecated => ['BIOSIM:METABOLITE:000001']  :param reason => 'duplicate HMDB accession'
MATCH (new_entity:Entity {canonical_id: $new_canonical_id})
UNWIND $deprecated AS old_id
MATCH (old:Entity {canonical_id: old_id})
WHERE old.canonical_id <> new_entity.canonical_id
SET old.status = 'deprecated'
MERGE (old)-[r:SUPERSEDED_BY]->(new_entity)
ON CREATE SET r.merged_at = datetime(), r.reason = $reason
RETURN old.canonical_id AS deprecated_id, new_entity.canonical_id AS current_id;


// =====================================================================
// INVARIANT CHECKS
//
// These replace the declarative constraints Community Edition does not
// offer. Every one must return ZERO rows. The migration plan requires them
// to run as assertions in the equivalence suite (step 4) and as a periodic
// health check afterwards — an undeclarable constraint that nothing checks
// is not a constraint at all.
// =====================================================================

// 1. Required properties present (stands in for NOT NULL).
MATCH (e:Entity)
WHERE e.canonical_id IS NULL OR e.entity_type IS NULL
   OR e.numeric_id IS NULL   OR e.name IS NULL
   OR e.name_lower IS NULL   OR e.status IS NULL
RETURN e.canonical_id AS offending_entity, 'missing required property' AS problem;

// 2. canonical_id is internally consistent with entity_type + numeric_id.
//    Guards against a hand-written node whose string and fields disagree,
//    which would make resolve() and get_record() reachable by different ids.
MATCH (e:Entity)
WHERE e.canonical_id <> 'BIOSIM:' + e.entity_type + ':' +
      right('00000' + toString(e.numeric_id), 6)
RETURN e.canonical_id AS offending_entity, 'canonical_id/field mismatch' AS problem;

// 3. No SUPERSEDED_BY cycles. A cycle would make the "end of chain"
//    traversal find no terminal entity, so resolve() would silently return
//    nothing for every accession routed through it.
MATCH (e:Entity)
WHERE exists { MATCH (e)-[:SUPERSEDED_BY*1..]->(e) }
RETURN e.canonical_id AS offending_entity, 'SUPERSEDED_BY cycle' AS problem;

// 4. RESOLVED (was: MASTER_PLAN Track 3.4's "Neo4j invariant-4
//    decision"). The PRIOR version of this check rejected ANY active
//    entity with two accessions from the same ontology source --
//    contradicting this schema's own "no constraint stricter than the
//    Python invariant" rule, because registry.py's own already-shipped
//    merge-inheritance fix makes this a LEGITIMATE, routine outcome:
//    merging entity A (UniProt P00001) into entity B (UniProt P99999)
//    correctly produces B.cross_references={'uniprot': ['P99999',
//    'P00001']} -- verified directly against a live IDRegistry. The
//    prior check would have flagged every such merge as a data-
//    integrity violation, even though nothing is wrong.
//
//    The REAL invariant cross_reference() actually depends on (see that
//    query's own updated comment above) is narrower: each accession
//    must have a well-defined, UNAMBIGUOUS rank, so "the primary
//    accession for this (entity, source) pair" has exactly one answer.
//    This checks that directly -- two DISTINCT accessions at the SAME
//    rank for the same (entity, source) is the actual ambiguity to
//    catch (a genuine backfill bug: two rows both claiming rank 0), not
//    multiplicity itself, which is expected and fine.
//
//    Disclosed limit (code review, MASTER_PLAN 3.4): zero rows from this
//    check means every (entity, source) pair has an UNAMBIGUOUS primary
//    -- it does NOT mean every rank is CORRECT. A backfill bug that
//    assigns technically-unique-but-wrong-order ranks (e.g. a reversed
//    list, an off-by-one in the source enumeration) produces zero
//    invariant violations while cross_reference() still silently
//    returns the wrong accession as "primary" -- exactly the class of
//    Python/Neo4j divergence this whole redesign exists to eliminate,
//    just not catchable from the graph's own data alone. Only the
//    equivalence suite (Track 3.8, not yet built) can catch a
//    unique-but-wrong-order backfill bug, by comparing against Python's
//    own real answer directly. The backfill script (3.7, not yet
//    written) MUST assign rank from each accession's true list-index
//    position, in order -- this invariant is a backstop against
//    ambiguity, not a substitute for that correctness requirement.
MATCH (x:ExternalId)-[r:IDENTIFIES]->(e:Entity)
WHERE e.status = 'active'
WITH e, x.source AS source, r.rank AS rank, collect(DISTINCT x.accession) AS accessions
WHERE size(accessions) > 1
RETURN e.canonical_id AS offending_entity, source, rank, accessions,
       'ambiguous primary: multiple accessions share one rank for one source' AS problem;

// 5. Every accession resolves to exactly one CURRENT entity.
//    A row here is precisely the AmbiguousResolutionError condition, found
//    proactively instead of at the moment some simulation needs the lookup.
MATCH (x:ExternalId)-[:IDENTIFIES]->(owner:Entity)
OPTIONAL MATCH (owner)-[:SUPERSEDED_BY*1..]->(cur:Entity)
WHERE NOT (cur)-[:SUPERSEDED_BY]->()
WITH x, collect(DISTINCT coalesce(cur, owner).canonical_id) AS current_owners
WHERE size(current_owners) > 1
RETURN x.key AS offending_accession, current_owners, 'ambiguous resolution' AS problem;

// 6. CLAUDE.md rule 6: Reactome ids are for PATHWAYs only. Proteins use
//    UniProt, metabolites use ChEBI/HMDB. Never expressible as a CE
//    constraint; entirely expressible as a check.
MATCH (x:ExternalId {source: 'reactome'})-[:IDENTIFIES]->(e:Entity)
WHERE e.entity_type <> 'PATHWAY'
RETURN e.canonical_id AS offending_entity, e.entity_type AS entity_type,
       x.accession AS reactome_accession,
       'Reactome id assigned to a non-pathway entity (rule 6)' AS problem;

// 7. A deprecated entity must have an outgoing SUPERSEDED_BY edge, and an
//    active one must not. Catches a half-applied merge.
MATCH (e:Entity)
WHERE (e.status = 'deprecated' AND NOT (e)-[:SUPERSEDED_BY]->())
   OR (e.status = 'active'     AND (e)-[:SUPERSEDED_BY]->())
RETURN e.canonical_id AS offending_entity, e.status AS status,
       'status inconsistent with SUPERSEDED_BY edges' AS problem;
