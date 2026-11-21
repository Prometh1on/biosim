// =====================================================================
// Biosim ID registry — Neo4j graph schema (Track 3).
//
// Target of migration: packages/id-resolver/src/biosim_id_resolver/
//   types.py     CanonicalID / ExternalID / ResolutionRecord / OntologySource
//   registry.py  IDRegistry (_canonical_index, _xref_index,
//                _deprecated_index, _name_index)
//
// Applied automatically by the `neo4j-init` service in docker-compose.yml.
// Every statement is IF NOT EXISTS / MERGE, so re-running is a no-op.
//
// ---------------------------------------------------------------------
// THE MODEL
//
//   (:Entity {
//       canonical_id  : 'BIOSIM:PROTEIN:000001'   // unique, the whole string
//       entity_type   : 'PROTEIN'                 // EntityType
//       numeric_id    : 1
//       name          : 'insulin receptor'
//       name_lower    : 'insulin receptor'         // search_by_name key
//       description   : '...'
//       status        : 'active' | 'deprecated'
//   })
//
//   (:ExternalId {
//       key       : 'uniprot:P06213'               // unique = source + ':' + accession
//       source    : 'uniprot'                      // OntologySource
//       accession : 'P06213'                       // already upper-cased by ExternalID
//   })
//
//   (:ExternalId)-[:IDENTIFIES]->(:Entity)
//   (:Entity)-[:SUPERSEDED_BY {merged_at, reason}]->(:Entity)
//   (:OntologySource {source, priority})           // reference data
//
// ---------------------------------------------------------------------
// WHY THIS SHAPE — and the three real merge-workflow bugs it removes
// STRUCTURALLY rather than by careful coding.
//
// All three were reproduced against the CURRENT in-memory registry in this
// session, not inferred from reading it. Chain A -> B -> C (each merge
// deprecating the previous entity, all sharing an HMDB accession):
//
//   BUG 1 — the deprecation hint is only one hop deep.
//     get_record(A) reports "the current canonical ID is B". B is itself
//     deprecated. A caller that follows the hint once lands on another
//     dead id. _deprecated_index is a flat map and does not transitively
//     collapse.
//   BUG 2 — _name_index is never cleaned on merge.
//     search_by_name('entity a') still returns A's deprecated canonical
//     id; get_record() on that result then raises. This is the same
//     "stale index entry after merge" class CLAUDE.md already records as
//     a previously-shipped-and-fixed bug, still present via this index.
//   BUG 3 — a transitive merge silently loses history across save/load.
//     to_dict() serialises only _canonical_index, and each merge pops the
//     deprecated record out of it. After A->B then B->C, the A->B mapping
//     exists nowhere in the serialised form: reloading and calling
//     get_record(A) returns "No record found" instead of pointing at C.
//     A canonical id that once resolved becomes permanently unresolvable.
//
// How this model removes each one:
//
//   * Deprecated entities are NODES THAT ARE NEVER DELETED, marked
//     status='deprecated' and linked by SUPERSEDED_BY. There is no
//     separate index to fall out of sync (bug 2) and nothing to drop on
//     serialisation, because the graph IS the persistence (bug 3).
//   * "Current" id is a variable-length traversal to the end of the
//     SUPERSEDED_BY chain, so an arbitrarily long merge history collapses
//     in one query (bug 1). This is the one genuinely graph-shaped query
//     in the whole registry.
//   * IDENTIFIES EDGES ARE IMMUTABLE. A merge does NOT rewrite them.
//     Resolution follows IDENTIFIES and then the chain, so an accession
//     can never be left pointing at a stale owner — the failure mode the
//     in-memory registry needed an explicit re-pointing fix for simply
//     cannot occur. It also preserves which entity ORIGINALLY claimed the
//     accession, which the in-memory version destroys.
//     Semantics are preserved, not changed: two accession owners that
//     later merged collapse to one current entity (so no false
//     ambiguity), while a genuine third-party collision still yields two
//     distinct current entities and must still raise. See
//     002_reference_queries.cypher.
//
// ---------------------------------------------------------------------
// COMMUNITY EDITION LIMITS THAT SHAPE WHAT IS EXPRESSIBLE HERE
//
// Only UNIQUENESS constraints and INDEXES exist in Neo4j Community.
// Property-existence constraints (`REQUIRE x.prop IS NOT NULL`), node-key
// constraints, and relationship-property constraints are Enterprise-only.
// So "every Entity must have an entity_type" CANNOT be declared here the
// way a Postgres NOT NULL can. That gap is covered by the invariant
// queries in 002_reference_queries.cypher, which the migration plan
// requires running as an assertion in the equivalence suite. This is a
// real, concrete asymmetry between the two halves of this migration, and
// one of the reasons open question 3 in the migration plan asks whether a
// graph server is warranted for this dataset at all.
// =====================================================================


// ---------------------------------------------------------------------
// Uniqueness constraints (each also creates a backing index)
// ---------------------------------------------------------------------

CREATE CONSTRAINT entity_canonical_id_unique IF NOT EXISTS
FOR (e:Entity) REQUIRE e.canonical_id IS UNIQUE;

CREATE CONSTRAINT external_id_key_unique IF NOT EXISTS
FOR (x:ExternalId) REQUIRE x.key IS UNIQUE;

CREATE CONSTRAINT ontology_source_unique IF NOT EXISTS
FOR (s:OntologySource) REQUIRE s.source IS UNIQUE;


// ---------------------------------------------------------------------
// Indexes
// ---------------------------------------------------------------------

// search_by_name(): exact, case-insensitive match on the lower-cased name.
CREATE INDEX entity_name_lower IF NOT EXISTS
FOR (e:Entity) ON (e.name_lower);

// Filtering the active set, and stats() by entity type.
CREATE INDEX entity_status IF NOT EXISTS
FOR (e:Entity) ON (e.status);
CREATE INDEX entity_type IF NOT EXISTS
FOR (e:Entity) ON (e.entity_type);

// resolve(): (source, accession) point lookup. The unique key index above
// already serves the composite form; this supports resolve_any(), which
// searches an accession across ALL sources in priority order.
CREATE INDEX external_id_accession IF NOT EXISTS
FOR (x:ExternalId) ON (x.accession);


// ---------------------------------------------------------------------
// Reference data: ontology sources and their resolution priority.
//
// Mirror of types.RESOLUTION_PRIORITY. Stored as data for the same reason
// the confidence/organism ranks became ref_* tables on the Postgres side:
// so ordering logic lives in one place. Lower priority number wins.
// KEGG is cross-reference only and must never be a primary resolution
// source (CLAUDE.md rule 6's companion rule, types.py OntologySource).
// ---------------------------------------------------------------------

MERGE (s:OntologySource {source: 'uniprot'})
  SET s.priority = 1, s.primary_for = 'proteins',   s.cross_reference_only = false;
MERGE (s:OntologySource {source: 'chebi'})
  SET s.priority = 2, s.primary_for = 'metabolites', s.cross_reference_only = false;
MERGE (s:OntologySource {source: 'hmdb'})
  SET s.priority = 3, s.primary_for = 'metabolites', s.cross_reference_only = false;
MERGE (s:OntologySource {source: 'reactome'})
  SET s.priority = 4, s.primary_for = 'pathways',    s.cross_reference_only = false;
MERGE (s:OntologySource {source: 'kegg'})
  SET s.priority = 5, s.primary_for = 'none',        s.cross_reference_only = true;


// ---------------------------------------------------------------------
// FUTURE (sketch only — nothing below is created yet, and no Biosim code
// produces this data today).
//
// The roadmap's stated justification for a graph database is "entity
// cross-references AND PATHWAY RELATIONSHIPS". Cross-references are the
// part that exists today and, honestly, they are a 1-hop lookup that
// Postgres would serve perfectly well. Pathway topology is the part that
// genuinely wants a graph, and it does not exist in the codebase at all
// yet — Reactome ids are currently stored only as flat cross-references.
// The intended shape, so the schema above does not have to be reworked
// when it arrives:
//
//   (:Entity {entity_type:'PROTEIN'})   -[:PARTICIPATES_IN]-> (:Entity {entity_type:'PATHWAY'})
//   (:Entity {entity_type:'PROTEIN'})   -[:CATALYSES]->       (:Entity {entity_type:'REACTION'})
//   (:Entity {entity_type:'METABOLITE'})-[:SUBSTRATE_OF]->    (:Entity {entity_type:'REACTION'})
//   (:Entity {entity_type:'REACTION'})  -[:PRODUCES]->        (:Entity {entity_type:'METABOLITE'})
//   (:Entity {entity_type:'REACTION'})  -[:PART_OF]->         (:Entity {entity_type:'PATHWAY'})
//
// Those are the multi-hop questions ("which pathways connect glucose to
// mTOR", "what reactions does this protein reach within N steps") that a
// relational schema handles poorly and that would make the second
// container earn its keep. Until then, see open question 3 in
// infra/MIGRATION_PLAN.md.
// ---------------------------------------------------------------------
