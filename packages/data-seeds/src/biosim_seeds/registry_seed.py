# SAVE TO: packages\data-seeds\src\biosim_seeds\registry_seed.py
"""
Real biological entity registrations for the IDRegistry.

Every entity in the v1 glucose/insulin pathway is registered here with
its canonical cross-references from Reactome, UniProt, ChEBI, HMDB, and KEGG.

All accessions are verified against their respective databases as of 2025.
IDs are stable and versioned; changes are tracked in the changelog below.

Changelog:
    2025-01-01  Initial seed: 12 proteins, 8 metabolites, 4 pathways
    2026-08-11  +22 metabolites (MASTER_PLAN 5.4): the fasted hepatic
                network's own metabolite set (`biosim_metabolism.
                fasted_hepatic_network.METABOLITES`, 4.3b), registered
                deliberately ahead of need per Track 5's own scoping,
                not reactively alongside a new Track 1 mechanism like
                every prior addition here. 8 of the network's 30
                metabolites already had entries above under a different
                display name (D-Glucose/glc, pyruvate/pyr, lactate/lac,
                ATP/atp, AMP/amp, glycogen/glycogen,
                glucose-6-phosphate/g6p, acetyl-CoA/accoa) and are not
                re-registered. Every ChEBI and KEGG accession below was
                independently verified live this session -- ChEBI via
                direct `https://www.ebi.ac.uk/chebi/searchId.do?
                chebiId=CHEBI:nnnnn` page fetches, KEGG via the plain-
                text `https://rest.kegg.jp/get/Cnnnnn` REST API -- both
                confirmed to actually resolve to the claimed compound
                name, not merely "a page loaded." HMDB accessions WERE
                found for most of these (via a research pass) but are
                deliberately NOT included in this batch: `hmdb.ca`
                returns HTTP 403 Forbidden to every automated fetch
                attempted this session (confirmed independently, not a
                one-off tool failure), so no HMDB accession here could
                be verified to this project's own standard (a real,
                loaded page confirming the name) rather than a weaker
                search-snippet inference -- a disclosed gap, not a
                silent omission, matching this project's own repeated
                preference for an honest gap over a lower-confidence
                inclusion. Revisit if a working HMDB verification path
                (e.g. their bulk XML/SDF download) is found later.
    2026-08-11  +1 protein, SGLT2/SLC5A2 (MASTER_PLAN 5.1): the SGLT2
                inhibitor drug class's molecular target, registered as
                the subject entity for renal_pathway.py's own real
                clinical GFR/TmG constants (the parameter-curation
                campaign's own next chunk). UniProt/KEGG both
                independently verified live -- see the registration's
                own inline comment for the exact fetch URLs used.
    2026-08-11  +1 metabolite, cortisol (MASTER_PLAN 5.1, chunk 6): the
                subject entity for cortisol_pathway.py's own exogenous-
                steroid-dosing PK/potency constants. ChEBI/KEGG both
                independently verified live. HMDB0000063 was located but
                deliberately NOT included, same disclosed reason as the
                Track 4 metabolite batch above (hmdb.ca blocks automated
                fetches).
    2026-08-12  +2 metabolites, acarbose/miglitol (MASTER_PLAN 5.1,
                chunk 8): the subject entities for drug_interactions.py's
                own alpha-glucosidase-inhibitor kabs-scale constants.
                Both accessions cross-verified via two independent live
                sources each (KEGG's own compound<->drug/ChEBI conversion
                APIs, then the ChEBI page itself) -- an initial guess at
                both ChEBI IDs was wrong and caught before registering.
    2026-08-18  +1 metabolite, pioglitazone (MASTER_PLAN 5.1, chunk 9):
                the subject entity for drug_interactions.py's own
                pioglitazone p3 (whole-body insulin sensitivity)
                dose-response constant (Track 1 row 1.25). ChEBI/KEGG
                both cross-verified via the same two-independent-source
                convention as acarbose/miglitol above.
    2026-08-18  +1 metabolite, protein (MASTER_PLAN 5.1, chunk 10): the
                subject entity for nutrient_secretion_pathway.py's own
                protein-co-ingestion insulin-secretion constant (Track 1
                row 1.28). Deliberately a GENERIC macronutrient-class
                entity (CHEBI:16541/KEGG:C00017), not a specific amino
                acid, matching Krezowski et al. 1986's own "50g protein"
                test substance. A mandatory code-reviewer pass correctly
                found the initial verification (KEGG's own DBLINKS field
                plus KEGG's own conversion API) was NOT actually two
                independent sources -- both read the same KEGG record.
                Fixed by independently checking ChEBI's OWN side: CHEBI:
                16541 reciprocally cross-references KEGG:C00017 in its
                own Xrefs section (genuinely independent confirmation),
                while the alternative candidate CHEBI:36080 has NO KEGG
                cross-reference at all on its own page -- a dead end for
                this campaign's verification convention, despite its
                more plainly-matching name. CHEBI:16541 used.
    2026-08-18  +1 metabolite, metformin (MASTER_PLAN 5.1, chunk 17): the
                subject entity for drug_interactions.py's own metformin
                Gb (fasting-glucose) reduction constant. ChEBI:6801/
                KEGG:C07151 (compound)/D04966 (drug) verified via two
                genuinely independent sources: KEGG's own D04966 DBLINKS
                field lists "ChEBI: 6801"; ChEBI:6801's own page
                (backend API `database_accessions.MANUAL_X_REF`)
                reciprocally lists "KEGG COMPOUND C07151" and "KEGG DRUG
                D04966" -- a genuine two-way cross-reference, not one
                source's assertion alone (same discipline as the
                protein/pioglitazone chunks' own corrected convention).
    2026-08-18  +1 metabolite, caffeine (MASTER_PLAN 5.1, chunk 18): the
                subject entity for drug_interactions.py's own caffeine
                p3 (insulin sensitivity) dose-response constant.
                ChEBI:27732/KEGG:C07481 (compound)/D00528 (drug)
                verified via the same two-independent-source convention
                as metformin above: KEGG's own D00528 DBLINKS field
                lists "ChEBI: 27732"; ChEBI:27732's own backend API
                reciprocally lists "KEGG COMPOUND C07481"/"KEGG DRUG
                D00528" under its own MANUAL_X_REF.
    2026-08-18  +1 metabolite, tirzepatide (MASTER_PLAN 5.1, chunk 19):
                the subject entity for drug_interactions.py's own
                tirzepatide Gb dose-response constants (3 real anchor
                points, Track 1 row 1.5). ChEBI:194186/KEGG:D11360 (DRUG
                entry only -- tirzepatide is a synthetic peptide with no
                small-molecule KEGG COMPOUND record) verified via two
                genuinely independent sources: KEGG's own D11360 DBLINKS
                field lists "ChEBI: 194186"; ChEBI:194186's own backend
                API reciprocally lists "KEGG DRUG D11360" under its own
                MANUAL_X_REF, name confirmed exactly "tirzepatide".
    2026-08-18  +1 metabolite, prednisone (MASTER_PLAN 5.1, chunk 20):
                the subject entity for drug_interactions.py's own
                glucocorticoid_modifiers() p3/Gb/secretion_scale dose-
                response constants. ChEBI:8382/KEGG:C07370 (COMPOUND
                entry ONLY -- deliberately did NOT also register KEGG's
                separate D00473 DRUG entry, since unlike metformin/
                caffeine/tirzepatide, neither KEGG record cross-links to
                the other for this substance, i.e. C07370 has no "Same
                as: D00473" REMARK and D00473 has no ChEBI DBLINK at
                all -- including D00473 anyway would assert an
                unverified compound<->drug relationship this campaign's
                own discipline doesn't permit). Verified via two
                genuinely independent sources: KEGG's own C07370
                DBLINKS field lists "ChEBI: 8382"; ChEBI:8382's own
                backend API reciprocally lists "KEGG COMPOUND C07370"
                under its own MANUAL_X_REF, name confirmed exactly
                "prednisone".
    2026-08-18  +2 metabolites, glibenclamide/glipizide (MASTER_PLAN 5.1,
                chunk 21): a documentation-completeness gap, not a new
                constant -- both drugs' OWN secretion-fraction constants
                (SULFONYLUREA_SECRETION_FRACTION_OF_KE_IB, GLIPIZIDE_
                SECRETION_FRACTION_OF_KE_IB) were already registered in
                ParameterStore under the "insulin" subject (chunks 7 and
                the Track 1 row-1.3b glipizide extension), correctly --
                the quantity is fundamentally a fraction of insulin's own
                clearance kinetics (ke*Ib), not an inherent chemical
                property of either drug, so it stays there. Neither drug
                itself, however, had ever been given its own canonical
                ChEBI/KEGG-cross-referenced IDRegistry entity, unlike
                every other real drug in this campaign (pioglitazone,
                metformin, caffeine, tirzepatide, prednisone). Both
                verified via two genuinely independent sources:
                glibenclamide (KEGG's own "Same as: C07022" REMARK links
                its D00336 DRUG entry to its own C07022 COMPOUND entry,
                whose DBLINKS lists "ChEBI: 5441"; ChEBI:5441's own
                backend API reciprocally lists both "KEGG COMPOUND
                C07022" and "KEGG DRUG D00336", name="glyburide" -- the
                USAN name for the same drug INN-named glibenclamide);
                glipizide (KEGG's own D00335 DRUG entry lists "DBLINKS
                ChEBI: 5384" directly, with NO separate COMPOUND entry
                or "Same as" link -- confirmed absent, matching
                prednisone's own precedent of registering only the
                verified accession rather than asserting an unverified
                link; ChEBI:5384's own backend API reciprocally lists
                "KEGG DRUG D00335", name="glipizide" exactly).
    2026-08-19  +1 protein, GLP-1/GCG (MASTER_PLAN 5.1, chunk 23): the
                subject entity for glp1_pathway.py's own basal/meal-
                response PK constants. Registered under proglucagon's own
                real accession (P01275), not a separate GLP-1-specific
                gene product, since GLP-1 is a cleavage product of one
                precursor -- matching UniProt's own "contains" annotation
                on P01275, which explicitly lists "Glucagon-like peptide
                1(7-37)"/"(7-36)" as processed products of this single
                entry. UniProt/KEGG both independently verified live:
                P01275's own record (uniProtkbId "GLUC_HUMAN", gene
                "GCG") and KEGG's hsa:2641 (SYMBOL line "GCG, GLP-1,
                GLP1, GLP2, GRPP") both resolve to the same underlying
                proglucagon gene product from two independent databases.
    2026-08-19  +1 metabolite, NEFA (MASTER_PLAN 5.1, chunk 24): the
                subject entity for nefa_pathway.py's own basal-clearance/
                antilipolysis-potency PK constants. Deliberately a
                GENERIC "fatty acid" class entity (CHEBI:35366/
                KEGG:C00162), same category as the "protein" chunk above
                -- NEFA is a clinical aggregate measurement (total
                circulating free fatty acids), not one specific
                molecular species. Both accessions independently
                verified live: ChEBI via the OLS API (label "fatty
                acid", own cross-reference lists "kegg.compound:
                C00162"); KEGG via `rest.kegg.jp/get/C00162` (NAME
                "Fatty acid") -- a genuine two-way cross-reference.
    2026-08-19  +0 entities, growth hormone p3 constant curated
                (MASTER_PLAN 5.1, chunk 25, following Track 1 row 1.35):
                no new registry entry needed -- "growth hormone"
                (UniProt P01241) was already registered in the original
                v1 seed (force_numeric_id=12). Caught before committing:
                an initial draft attempted a duplicate registration,
                correctly rejected by `ValueError: Cross-reference
                collision` at seed time. See parameter_seed.py's
                `_seed_growth_hormone_p3()` for the actual new content
                this chunk adds (GROWTH_HORMONE_P3_FACTOR, PMID:1916055).
    2026-08-19  +1 metabolite, hydrochlorothiazide (MASTER_PLAN 5.1,
                chunk 26, following Track 1 row 1.34): the subject
                entity for drug_interactions.py's own
                THIAZIDE_P3_FACTOR_AT_TESTED_DOSE chronic p3 dose-
                response constant. ChEBI:5778/KEGG:D00340 verified via
                two genuinely independent sources: ChEBI's own
                cross-reference lists "kegg.drug:D00340"; KEGG's own
                D00340 DBLINKS field reciprocally lists "ChEBI: 5778".
    2026-08-19  +1 metabolite, indinavir (MASTER_PLAN 5.1, chunk 27,
                following Track 1 row 1.33): the subject entity for
                drug_interactions.py's own
                INDINAVIR_P3_FACTOR_AT_TESTED_DOSE acute p3 dose-response
                constant. ChEBI:44032/KEGG:C07051 verified via two
                genuinely independent sources: ChEBI's own cross-
                reference lists "kegg.compound:C07051"; KEGG's own
                C07051 DBLINKS field reciprocally lists "ChEBI: 44032".
    2026-08-19  +1 metabolite, atorvastatin (MASTER_PLAN 5.1, chunk 28,
                following Track 1 row 1.41): the subject entity for
                drug_interactions.py's own ATORVASTATIN_P3_FACTOR_AT_DOSE
                chronic p3 dose-response constants (3 real anchor
                points). ChEBI:39548/KEGG:C06834 verified via two
                genuinely independent sources: ChEBI's own cross-
                reference lists "kegg.compound:C06834"; KEGG's own
                C06834 DBLINKS field reciprocally lists "ChEBI: 39548".
    2026-08-20  +1 metabolite, alcohol/ethanol (MASTER_PLAN 5.1, chunk
                36): the subject entity for alcohol_hepatic_pathway.py's
                own ALCOHOL_EGP_REDUCTION_FRACTION_AT_TESTED_DOSE acute
                hepatic-glucose-production-suppression constant.
                ChEBI:16236/KEGG:C00469 verified via two genuinely
                independent sources: KEGG's own C00469 DBLINKS field
                lists "ChEBI: 16236"; ChEBI's own cross-reference
                reciprocally lists "kegg.compound:C00469".
    2026-10-19  +1 metabolite, hydroxychloroquine (MASTER_PLAN 5.1, chunk
                37, following Track 1 row 1.44): the subject entity for
                drug_interactions.py's own HYDROXYCHLOROQUINE_P3_FACTOR.
                ChEBI:5801/KEGG:C07043 verified via two genuinely
                independent sources: ChEBI's own cross-reference lists
                "kegg.compound:C07043" (and "kegg.drug:D08050"); KEGG's
                own C07043 AND D08050 DBLINKS fields each reciprocally
                list "ChEBI: 5801".
"""

from __future__ import annotations

from biosim_id_resolver import EntityType, IDRegistry, OntologySource


def seed_registry(registry: IDRegistry) -> IDRegistry:
    """
    Populate an IDRegistry with all v1 pathway entities.

    Args:
        registry: IDRegistry to seed (modified in place)

    Returns:
        The same registry, now populated.
    """
    _seed_proteins(registry)
    _seed_metabolites(registry)
    _seed_hepatic_network_metabolites(registry)
    _seed_pathways(registry)
    return registry


def create_seeded_registry() -> IDRegistry:
    """Create and return a new IDRegistry pre-seeded with all v1 entities."""
    return seed_registry(IDRegistry())


# ---------------------------------------------------------------------------
# Proteins
# ---------------------------------------------------------------------------

def _seed_proteins(r: IDRegistry) -> None:
    # Insulin
    r.register(
        name="insulin",
        entity_type=EntityType.PROTEIN,
        description="Pancreatic hormone regulating glucose uptake and metabolism",
        cross_references={
            OntologySource.UNIPROT:  ["P01308"],
            OntologySource.HMDB:     ["HMDB0000981"],
            OntologySource.REACTOME: ["R-HSA-264876"],
        },
        force_numeric_id=1,
    )

    # Insulin receptor (INSR)
    r.register(
        name="insulin receptor",
        entity_type=EntityType.PROTEIN,
        description="Receptor tyrosine kinase mediating cellular insulin response",
        cross_references={
            OntologySource.UNIPROT:  ["P06213"],
            OntologySource.KEGG:     ["K04527"],
            OntologySource.REACTOME: ["R-HSA-2428928"],
        },
        force_numeric_id=2,
    )

    # Insulin-like growth factor 1 (IGF-1)
    r.register(
        name="IGF-1",
        entity_type=EntityType.PROTEIN,
        description="Insulin-like growth factor 1; mediates growth hormone signalling",
        cross_references={
            OntologySource.UNIPROT:  ["P05019"],
            OntologySource.HMDB:     ["HMDB0005000"],
        },
        force_numeric_id=3,
    )

    # IGF-1 receptor (IGF1R)
    r.register(
        name="IGF-1 receptor",
        entity_type=EntityType.PROTEIN,
        description="Receptor tyrosine kinase for IGF-1; homologous to insulin receptor",
        cross_references={
            OntologySource.UNIPROT:  ["P08069"],
            OntologySource.KEGG:     ["K05069"],
        },
        force_numeric_id=4,
    )

    # mTORC1 (represented by mTOR catalytic subunit)
    r.register(
        name="mTOR",
        entity_type=EntityType.PROTEIN,
        description="Mechanistic target of rapamycin; central regulator of cell growth and metabolism",
        cross_references={
            OntologySource.UNIPROT:  ["P42345"],
            OntologySource.KEGG:     ["K07203"],
        },
        force_numeric_id=5,
    )

    # AMPK alpha-1 subunit (PRKAA1)
    r.register(
        name="AMPK",
        entity_type=EntityType.PROTEIN,
        description="AMP-activated protein kinase; master energy sensor",
        cross_references={
            OntologySource.UNIPROT:  ["Q13131"],
            OntologySource.KEGG:     ["K07198"],
            OntologySource.REACTOME: ["R-HSA-380972"],
        },
        force_numeric_id=6,
    )

    # AKT1 (Protein kinase B)
    r.register(
        name="AKT1",
        entity_type=EntityType.PROTEIN,
        description="Serine/threonine kinase AKT; mediates insulin signalling to mTOR",
        cross_references={
            OntologySource.UNIPROT:  ["P31749"],
            OntologySource.KEGG:     ["K04456"],
            OntologySource.REACTOME: ["R-HSA-198522"],
        },
        force_numeric_id=7,
    )

    # IRS-1 (Insulin receptor substrate 1)
    r.register(
        name="IRS-1",
        entity_type=EntityType.PROTEIN,
        description="Insulin receptor substrate 1; docking protein for insulin signalling",
        cross_references={
            OntologySource.UNIPROT:  ["P35568"],
            OntologySource.KEGG:     ["K16172"],
        },
        force_numeric_id=8,
    )

    # Glucokinase (GCK) — glucose sensor in beta cells
    r.register(
        name="glucokinase",
        entity_type=EntityType.PROTEIN,
        description="Hexokinase-4; glucose sensor in pancreatic beta cells and liver",
        cross_references={
            OntologySource.UNIPROT:  ["P35557"],
            OntologySource.KEGG:     ["K12407"],
            OntologySource.REACTOME: ["R-HSA-70262"],
        },
        force_numeric_id=9,
    )

    # GLUT4 (SLC2A4) — insulin-regulated glucose transporter
    r.register(
        name="GLUT4",
        entity_type=EntityType.PROTEIN,
        description="Insulin-regulated glucose transporter; translocates to membrane in response to insulin",
        cross_references={
            OntologySource.UNIPROT:  ["P14672"],
            OntologySource.KEGG:     ["K07191"],
            OntologySource.REACTOME: ["R-HSA-170796"],
        },
        force_numeric_id=10,
    )

    # GLUT2 (SLC2A2) — hepatic/beta-cell glucose transporter
    r.register(
        name="GLUT2",
        entity_type=EntityType.PROTEIN,
        description="High-Km glucose transporter in liver and beta cells",
        cross_references={
            OntologySource.UNIPROT:  ["P11168"],
            OntologySource.KEGG:     ["K07593"],
        },
        force_numeric_id=11,
    )

    # Growth hormone (GH) — stimulates IGF-1 production
    r.register(
        name="growth hormone",
        entity_type=EntityType.PROTEIN,
        description="Pituitary growth hormone; stimulates hepatic IGF-1 production",
        cross_references={
            OntologySource.UNIPROT:  ["P01241"],
            OntologySource.HMDB:     ["HMDB0062579"],
        },
        force_numeric_id=12,
    )

    # SGLT2 (SLC5A2) — renal sodium/glucose cotransporter, the SGLT2
    # inhibitor drug class's molecular target (MASTER_PLAN 5.1, Track 5
    # curation campaign -- registered ahead of the renal_pathway.py
    # parameter curation that needs a subject entity for it). Both
    # accessions independently verified live this session: UniProt via
    # https://rest.uniprot.org/uniprotkb/search?query=gene:SLC5A2 (P31639,
    # reviewed Swiss-Prot, "Sodium/glucose cotransporter 2"); KEGG via
    # https://rest.kegg.jp/get/hsa:6524 (ORTHOLOGY K14382, DRUG_TARGET
    # list explicitly includes dapagliflozin/empagliflozin/canagliflozin).
    r.register(
        name="SGLT2",
        entity_type=EntityType.PROTEIN,
        description="Sodium/glucose cotransporter 2; primary renal glucose reabsorption transporter, target of SGLT2 inhibitor drugs",
        cross_references={
            OntologySource.UNIPROT:  ["P31639"],
            OntologySource.KEGG:     ["K14382"],
        },
        force_numeric_id=13,
    )

    # GLP-1 (glucagon-like peptide 1) -- registered as the subject entity
    # for glp1_pathway.py's own basal/meal-response PK constants
    # (MASTER_PLAN 5.1, Track 5 curation campaign chunk 23). GLP-1 is a
    # proglucagon-cleavage product, not its own separate gene product --
    # registered under proglucagon's own real accession, same convention
    # UniProt itself uses (P01275's own "contains" section explicitly
    # lists "Glucagon-like peptide 1(7-37)"/"(7-36)" as processed
    # products of this one precursor). Both accessions independently
    # verified live this session: UniProt via
    # https://rest.uniprot.org/uniprotkb/P01275.json (uniProtkbId
    # "GLUC_HUMAN", recommendedName "Pro-glucagon", gene "GCG", contains
    # "Glucagon-like peptide 1" explicitly); KEGG via
    # https://rest.kegg.jp/get/hsa:2641 (SYMBOL line lists "GCG, GLP-1,
    # GLP1, GLP2, GRPP", NAME "pro-glucagon preproprotein") -- the same
    # underlying gene product confirmed from two independent databases.
    r.register(
        name="GLP-1",
        entity_type=EntityType.PROTEIN,
        description=(
            "Glucagon-like peptide 1; gut-derived incretin hormone, a "
            "proglucagon (GCG) cleavage product -- amplifies glucose-"
            "dependent insulin secretion, suppresses glucagon, slows "
            "gastric emptying"
        ),
        cross_references={
            OntologySource.UNIPROT:  ["P01275"],
            OntologySource.KEGG:     ["hsa:2641"],
        },
        force_numeric_id=14,
    )

    # Growth hormone's own subject entity for drug_interactions.py's
    # GROWTH_HORMONE_P3_FACTOR (MASTER_PLAN 5.1, Track 5 curation
    # campaign chunk 25, following Track 1 row 1.35) does NOT need a new
    # registration -- "growth hormone" (UniProt P01241) was already
    # registered in the original v1 seed above (force_numeric_id=12,
    # see `_seed_proteins()`'s own "Growth hormone (GH)" entry). Caught
    # before committing: an initial draft of this chunk attempted to
    # re-register it under a new force_numeric_id, which correctly
    # raised `ValueError: Cross-reference collision` at seed time --
    # confirming the existing entity, not a gap. `_seed_growth_hormone_p3()`
    # in parameter_seed.py resolves this ALREADY-EXISTING entity by name.


# ---------------------------------------------------------------------------
# Metabolites
# ---------------------------------------------------------------------------

def _seed_metabolites(r: IDRegistry) -> None:
    # D-Glucose
    r.register(
        name="D-Glucose",
        entity_type=EntityType.METABOLITE,
        description="Primary blood sugar; principal substrate for glycolysis",
        cross_references={
            OntologySource.CHEBI:    ["CHEBI:17234"],
            OntologySource.HMDB:     ["HMDB0000122"],
            OntologySource.KEGG:     ["C00031"],
            OntologySource.REACTOME: ["R-ALL-29926"],
        },
        force_numeric_id=1,
    )

    # Pyruvate
    r.register(
        name="pyruvate",
        entity_type=EntityType.METABOLITE,
        description="End product of glycolysis; substrate for TCA cycle and lactate",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:15361"],
            OntologySource.HMDB:  ["HMDB0000243"],
            OntologySource.KEGG:  ["C00022"],
        },
        force_numeric_id=2,
    )

    # Lactate
    r.register(
        name="lactate",
        entity_type=EntityType.METABOLITE,
        description="Anaerobic glycolysis product; gluconeogenic substrate",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:16651"],
            OntologySource.HMDB:  ["HMDB0000190"],
            OntologySource.KEGG:  ["C00186"],
        },
        force_numeric_id=3,
    )

    # ATP
    r.register(
        name="ATP",
        entity_type=EntityType.METABOLITE,
        description="Adenosine triphosphate; primary cellular energy currency",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:30616"],
            OntologySource.HMDB:  ["HMDB0000538"],
            OntologySource.KEGG:  ["C00002"],
        },
        force_numeric_id=4,
    )

    # AMP
    r.register(
        name="AMP",
        entity_type=EntityType.METABOLITE,
        description="Adenosine monophosphate; AMPK activator when AMP/ATP ratio rises",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:16027"],
            OntologySource.HMDB:  ["HMDB0000045"],
            OntologySource.KEGG:  ["C00020"],
        },
        force_numeric_id=5,
    )

    # Glycogen
    r.register(
        name="glycogen",
        entity_type=EntityType.METABOLITE,
        description="Glucose polymer; primary carbohydrate storage in liver and muscle",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:28087"],
            OntologySource.HMDB:  ["HMDB0000757"],
            OntologySource.KEGG:  ["C00182"],
        },
        force_numeric_id=6,
    )

    # Glucose-6-phosphate
    r.register(
        name="glucose-6-phosphate",
        entity_type=EntityType.METABOLITE,
        description="First committed intermediate of glycolysis",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:4170"],
            OntologySource.HMDB:  ["HMDB0001401"],
            OntologySource.KEGG:  ["C00092"],
        },
        force_numeric_id=7,
    )

    # Acetyl-CoA
    r.register(
        name="acetyl-CoA",
        entity_type=EntityType.METABOLITE,
        description="Central metabolic hub; substrate for TCA cycle and lipogenesis",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:15351"],
            OntologySource.HMDB:  ["HMDB0001206"],
            OntologySource.KEGG:  ["C00024"],
        },
        force_numeric_id=8,
    )


# ---------------------------------------------------------------------------
# Fasted hepatic network metabolites (MASTER_PLAN 5.4, Track 4's 4.3b network)
# ---------------------------------------------------------------------------

def _seed_hepatic_network_metabolites(r: IDRegistry) -> None:
    """
    The 22 metabolites `biosim_metabolism.fasted_hepatic_network` needs
    that weren't already registered above under a different name. See
    this module's own docstring for the verification method and the
    disclosed HMDB gap (found but not independently verifiable this
    session, so deliberately omitted rather than guessed).

    `force_numeric_id` continues the METABOLITE sequence from 9 (the
    last id `_seed_metabolites` above assigned was 8, acetyl-CoA).
    """
    r.register(
        name="glucose-1-phosphate",
        entity_type=EntityType.METABOLITE,
        description="Product of glycogen phosphorylase; isomerizes to glucose-6-phosphate",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:29042"],
            OntologySource.KEGG:  ["C00103"],
        },
        force_numeric_id=9,
    )
    r.register(
        name="fructose-6-phosphate",
        entity_type=EntityType.METABOLITE,
        description="Glycolytic/gluconeogenic intermediate, isomer of glucose-6-phosphate",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:15946"],
            OntologySource.KEGG:  ["C00085"],
        },
        force_numeric_id=10,
    )
    r.register(
        name="fructose-1,6-bisphosphate",
        entity_type=EntityType.METABOLITE,
        description="Aldolase substrate; splits into two triose-phosphate units",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:16905"],
            OntologySource.KEGG:  ["C00354"],
        },
        force_numeric_id=11,
    )
    r.register(
        name="glyceraldehyde-3-phosphate",
        entity_type=EntityType.METABOLITE,
        description=(
            "Triose phosphate; in this network's own lumped 'FBA_TPI' "
            "reaction this represents a combined DHAP/G3P pool, not "
            "only the stereospecific G3P species registered here"
        ),
        cross_references={
            OntologySource.CHEBI: ["CHEBI:29052"],
            OntologySource.KEGG:  ["C00118"],
        },
        force_numeric_id=12,
    )
    r.register(
        name="phosphoenolpyruvate",
        entity_type=EntityType.METABOLITE,
        description="PEPCK's product; final gluconeogenic-backbone intermediate before glycolysis reverses",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:44897"],
            OntologySource.KEGG:  ["C00074"],
        },
        force_numeric_id=13,
    )
    r.register(
        name="oxaloacetate",
        entity_type=EntityType.METABOLITE,
        # ChEBI's own primary name for CHEBI:16452 is "oxaloacetate(2-)"
        # (the anionic form); "oxaloacetate" is a listed synonym on the
        # same entry, not a different/wrong accession -- same naming
        # pattern as alpha-ketoglutarate/ADP/inorganic phosphate/oxygen
        # below, disclosed here once rather than repeated at each site.
        description="TCA-cycle/gluconeogenesis entry point; PEPCK/PC substrate-product",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:16452"],
            OntologySource.KEGG:  ["C00036"],
        },
        force_numeric_id=14,
    )
    r.register(
        name="alpha-ketoglutarate",
        entity_type=EntityType.METABOLITE,
        # CHEBI:16810's own primary name is "2-oxoglutarate(2-)"; same
        # anion/common-name naming pattern as oxaloacetate above.
        description="TCA-cycle intermediate; transamination partner for alanine/glutamate cycling",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:16810"],
            OntologySource.KEGG:  ["C00026"],
        },
        force_numeric_id=15,
    )
    r.register(
        name="glutamate",
        entity_type=EntityType.METABOLITE,
        description="ALT/GLUD1 transamination partner; nitrogen carrier toward the urea sink",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:16015"],
            OntologySource.KEGG:  ["C00025"],
        },
        force_numeric_id=16,
    )
    r.register(
        name="alanine",
        entity_type=EntityType.METABOLITE,
        description="Gluconeogenic amino-acid precursor (ALT transamination to pyruvate)",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:16977"],
            OntologySource.KEGG:  ["C00041"],
        },
        force_numeric_id=17,
    )
    r.register(
        name="ammonium",
        entity_type=EntityType.METABOLITE,
        description="GLUD1 deamination product; consumed by this network's urea-disposal sink",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:28938"],
            OntologySource.KEGG:  ["C01342"],  # NH4+ specifically, NOT C00014 (neutral NH3)
        },
        force_numeric_id=18,
    )
    r.register(
        name="glycerol",
        entity_type=EntityType.METABOLITE,
        description="Gluconeogenic precursor from lipolysis; the one precursor that funds its own conversion",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:17754"],
            OntologySource.KEGG:  ["C00116"],
        },
        force_numeric_id=19,
    )
    r.register(
        name="glycerol-3-phosphate",
        entity_type=EntityType.METABOLITE,
        description="GLYK's product; GPD1 substrate en route to the triose-phosphate pool",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:15978"],
            OntologySource.KEGG:  ["C00093"],
        },
        force_numeric_id=20,
    )
    r.register(
        name="palmitate",
        entity_type=EntityType.METABOLITE,
        description=(
            "C16:0 fatty acid; this network's FAO reaction stoichiometry "
            "(8 CoA, 7 FAD, 7 NAD consumed per firing) is specifically "
            "palmitate's own beta-oxidation stoichiometry, not a generic "
            "fatty acid"
        ),
        cross_references={
            OntologySource.CHEBI: ["CHEBI:15756"],
            OntologySource.KEGG:  ["C00249"],
        },
        force_numeric_id=21,
    )
    r.register(
        name="carbon dioxide",
        entity_type=EntityType.METABOLITE,
        description="PC/PEPCK byproduct; this network's own CO2 exchange bound (energy_bounds.py)",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:16526"],
            OntologySource.KEGG:  ["C00011"],
        },
        force_numeric_id=22,
    )
    r.register(
        name="oxygen",
        entity_type=EntityType.METABOLITE,
        # ChEBI's own primary name for CHEBI:15379 is "dioxygen" (the
        # molecular O2 species); same anion/common-name naming pattern
        # as oxaloacetate above, just IUPAC-vs-common here rather than
        # neutral-vs-anion.
        description="Terminal electron acceptor; this network's own O2 exchange bound (energy_bounds.py)",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:15379"],
            OntologySource.KEGG:  ["C00007"],
        },
        force_numeric_id=23,
    )
    r.register(
        name="ADP",
        entity_type=EntityType.METABOLITE,
        description="Adenosine diphosphate; ATP's discharged counterpart in every phosphorylation reaction",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:16761"],
            OntologySource.KEGG:  ["C00008"],
        },
        force_numeric_id=24,
    )
    r.register(
        name="inorganic phosphate",
        entity_type=EntityType.METABOLITE,
        # CHEBI:18367's own primary name is "phosphate(3-)"; same
        # anion/common-name naming pattern as oxaloacetate above.
        description="Orthophosphate; released or consumed by nearly every reaction in this network",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:18367"],
            OntologySource.KEGG:  ["C00009"],
        },
        force_numeric_id=25,
    )
    r.register(
        name="NAD+",
        entity_type=EntityType.METABOLITE,
        description="Oxidized redox cofactor; consumed by this network's dehydrogenase reactions",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:15846"],
            OntologySource.KEGG:  ["C00003"],
        },
        force_numeric_id=26,
    )
    r.register(
        name="NADH",
        entity_type=EntityType.METABOLITE,
        description="Reduced redox cofactor; produced by this network's dehydrogenase reactions",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:16908"],
            OntologySource.KEGG:  ["C00004"],
        },
        force_numeric_id=27,
    )
    r.register(
        name="FAD",
        entity_type=EntityType.METABOLITE,
        description="Oxidized redox cofactor; consumed by this network's FAO reaction",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:16238"],
            OntologySource.KEGG:  ["C00016"],
        },
        force_numeric_id=28,
    )
    r.register(
        name="FADH2",
        entity_type=EntityType.METABOLITE,
        description="Reduced redox cofactor; produced by this network's FAO reaction",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:17877"],
            OntologySource.KEGG:  ["C01352"],
        },
        force_numeric_id=29,
    )
    r.register(
        name="coenzyme A",
        entity_type=EntityType.METABOLITE,
        description="Free thiol form (CoA-SH); acyl-group carrier, distinct from acetyl-CoA",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:15346"],
            OntologySource.KEGG:  ["C00010"],
        },
        force_numeric_id=30,
    )

    # Cortisol -- registered as the subject entity for cortisol_pathway.py's
    # own exogenous-steroid-dosing PK/potency constants (MASTER_PLAN 5.1,
    # chunk 6). Both accessions independently verified live: ChEBI via
    # https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:17650 (page
    # title "cortisol (CHEBI:17650)", definition confirms "a corticosteroid
    # hormone or glucocorticoid produced by zona fasciculata of the adrenal
    # cortex"); KEGG via https://rest.kegg.jp/get/C00735 (NAME block lists
    # "Cortisol; Hydrocortisone; ..."). HMDB0000063 was located during
    # research but NOT included, same disclosed reason as the Track 4
    # metabolite batch above: hmdb.ca returns HTTP 403 to every automated
    # fetch attempted this session, so it could not be verified to this
    # project's own standard (a real, loaded page confirming the name).
    r.register(
        name="cortisol",
        entity_type=EntityType.METABOLITE,
        description="Primary glucocorticoid counter-regulatory hormone; subject for exogenous-steroid-dosing PK/potency constants",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:17650"],
            OntologySource.KEGG:  ["C00735"],
        },
        force_numeric_id=31,
    )

    # Acarbose/miglitol -- registered as the subject entities for
    # drug_interactions.py's own alpha-glucosidase-inhibitor kabs-scale
    # constants (MASTER_PLAN 5.1, chunk 8). Both accessions independently
    # verified live via TWO sources each: KEGG's own compound<->drug
    # cross-reference (https://rest.kegg.jp/get/D00216 / D00625, each
    # REMARK-linking to its compound-side C-number) AND KEGG's own
    # compound-to-ChEBI conversion API
    # (https://rest.kegg.jp/conv/chebi/cpd:C06802 / C07708), then the
    # resulting ChEBI IDs confirmed directly via
    # https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:nnnnn (page
    # titles "acarbose (CHEBI:2376)" / "Miglitol (CHEBI:6935)"). An
    # initial guess at both ChEBI IDs (2549, 6906) was WRONG (those
    # resolve to albuterol and metoprolol tartrate respectively) --
    # caught before registering anything, by cross-checking the guess
    # against the KEGG-derived ID rather than trusting either alone.
    r.register(
        name="acarbose",
        entity_type=EntityType.METABOLITE,
        description="Alpha-glucosidase inhibitor drug; subject for the kabs (meal absorption rate) reduction constants",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:2376"],
            OntologySource.KEGG:  ["C06802"],
        },
        force_numeric_id=32,
    )
    r.register(
        name="miglitol",
        entity_type=EntityType.METABOLITE,
        description="Alpha-glucosidase inhibitor drug (second in class); subject for the kabs (meal absorption rate) reduction constants",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:6935"],
            OntologySource.KEGG:  ["C07708"],
        },
        force_numeric_id=33,
    )

    # Pioglitazone -- registered as the subject entity for
    # drug_interactions.py's own pioglitazone p3 (whole-body insulin
    # sensitivity) dose-response constant (MASTER_PLAN 5.1, chunk 9;
    # Track 1 row 1.25). Both accessions independently verified live via
    # TWO sources, same convention as acarbose/miglitol above: KEGG's own
    # drug entry (https://rest.kegg.jp/get/D08378, "Pioglitazone (INN)")
    # lists "DBLINKS ChEBI: 8228" and "Same as: C07675" directly; KEGG's
    # own compound-to-ChEBI conversion API
    # (https://rest.kegg.jp/conv/chebi/cpd:C07675) independently returns
    # the identical chebi:8228; ChEBI's own page
    # (https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:8228) title
    # confirms "pioglitazone (CHEBI:8228)".
    r.register(
        name="pioglitazone",
        entity_type=EntityType.METABOLITE,
        description="Thiazolidinedione (PPAR-gamma) insulin-sensitizer drug; subject for the p3 dose-response constant",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:8228"],
            OntologySource.KEGG:  ["C07675"],
        },
        force_numeric_id=34,
    )

    # Protein (generic macronutrient/polypeptide-chain class, not a
    # specific amino acid) -- registered as the subject entity for
    # nutrient_secretion_pathway.py's own protein-co-ingestion insulin-
    # secretion constant (MASTER_PLAN 5.1, chunk 10; Track 1 row 1.28).
    # Unlike every prior chunk's SPECIFIC small-molecule drug/hormone
    # subject, this is deliberately a GENERIC macronutrient-class entity
    # (matching Krezowski et al. 1986's own "50g protein" test substance,
    # not one specific protein/amino acid).
    #
    # ChEBI:16541 vs. CHEBI:36080 -- a real, disclosed judgment call, not
    # a simple transcription-error check like acarbose/miglitol's own
    # precedent (a mandatory code-reviewer pass correctly flagged that
    # "KEGG's own DBLINKS field + KEGG's own conversion API" are NOT two
    # independent sources -- both read the same underlying KEGG record,
    # so this needed a genuinely separate, ChEBI-side check before being
    # trusted, which was then done):
    #   - CHEBI:16541's OWN page (not KEGG's) independently lists a
    #     manually-curated, RECIPROCAL cross-reference to "KEGG COMPOUND
    #     C00017" under its own Xrefs/Databases section -- this is ChEBI
    #     curators confirming the same link from their side, genuinely
    #     independent of KEGG's own DBLINKS assertion. CHEBI:16541's own
    #     definition ("A naturally occurring polypeptide synthesized at
    #     the ribosome") matches "a protein" as a listed synonym/name.
    #   - CHEBI:36080 (the initial ChEBI-search-only candidate) was
    #     checked directly and has NO KEGG cross-reference at all on its
    #     own page (confirmed absent, not merely unmentioned) -- a dead
    #     end for this campaign's own two-independent-source
    #     verification convention. Its own definition ("A biological
    #     macromolecule minimally consisting of one polypeptide chain
    #     synthesized at the ribosome") is near-synonymous with
    #     CHEBI:16541's, so nothing of chemical substance is lost by
    #     using the properly cross-referenced entry instead.
    # CHEBI:16541 is therefore the correct choice: independently
    # confirmed from BOTH KEGG's side (compound entry + conversion API)
    # AND ChEBI's own side (reciprocal Xref + matching definition/
    # synonym), while CHEBI:36080 has no comparable cross-database
    # verification path at all.
    r.register(
        name="protein",
        entity_type=EntityType.METABOLITE,
        description="Generic dietary protein/polypeptide macronutrient class; subject for the protein-co-ingestion insulin-secretion constant",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:16541"],
            OntologySource.KEGG:  ["C00017"],
        },
        force_numeric_id=35,
    )

    # Metformin -- registered as the subject entity for
    # drug_interactions.py's own metformin_modifiers() Gb (fasting-
    # glucose) reduction constant (MASTER_PLAN 5.1, chunk 17). Both
    # accessions independently verified live via TWO genuinely separate
    # sources, same convention as protein/pioglitazone above: KEGG's own
    # drug entry (https://rest.kegg.jp/get/D04966, "Metformin (USAN/INN)",
    # "Same as: C07151") lists "DBLINKS ChEBI: 6801" directly; ChEBI:6801's
    # own page (fetched via its backend JSON API, not just KEGG's
    # assertion) reciprocally lists "KEGG COMPOUND C07151" and "KEGG DRUG
    # D04966" under its own database_accessions.MANUAL_X_REF section --
    # a genuine two-way cross-reference between the two databases, not
    # one source's claim taken alone.
    r.register(
        name="metformin",
        entity_type=EntityType.METABOLITE,
        description="Biguanide antidiabetic drug; subject for the Gb (fasting-glucose) reduction constant",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:6801"],
            OntologySource.KEGG:  ["C07151"],
        },
        force_numeric_id=36,
    )

    # Caffeine -- registered as the subject entity for
    # drug_interactions.py's own caffeine_modifiers() p3 (insulin
    # sensitivity) dose-response constant (MASTER_PLAN 5.1, chunk 18).
    # Both accessions independently verified live via TWO genuinely
    # separate sources, same convention as metformin above: KEGG's own
    # drug entry (https://rest.kegg.jp/get/D00528, "Same as: C07481")
    # lists "DBLINKS ChEBI: 27732" directly; ChEBI:27732's own page
    # (backend JSON API) reciprocally lists "KEGG COMPOUND C07481" and
    # "KEGG DRUG D00528" under its own database_accessions.MANUAL_X_REF
    # section, name confirmed exactly "caffeine".
    r.register(
        name="caffeine",
        entity_type=EntityType.METABOLITE,
        description="Xanthine alkaloid stimulant; subject for the p3 (insulin sensitivity) dose-response constant",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:27732"],
            OntologySource.KEGG:  ["C07481"],
        },
        force_numeric_id=37,
    )

    # Tirzepatide -- registered as the subject entity for
    # drug_interactions.py's own tirzepatide_gb_factor() dose-response
    # constants (MASTER_PLAN 5.1, chunk 19; Track 1 row 1.5). A
    # synthetic peptide (dual GIP/GLP-1 agonist) with no KEGG COMPOUND
    # record -- only a KEGG DRUG entry exists, unlike every prior small-
    # molecule drug in this campaign. Both accessions independently
    # verified live via TWO genuinely separate sources: KEGG's own drug
    # entry (https://rest.kegg.jp/get/D11360, "Tirzepatide (JAN/USAN)")
    # lists "DBLINKS ChEBI: 194186" directly; ChEBI:194186's own page
    # (backend JSON API) reciprocally lists "KEGG DRUG D11360" under its
    # own database_accessions.MANUAL_X_REF section, name confirmed
    # exactly "tirzepatide". Registered as METABOLITE, same convention
    # every other drug in this campaign uses (this platform's
    # EntityType enum has no dedicated "peptide drug" category).
    r.register(
        name="tirzepatide",
        entity_type=EntityType.METABOLITE,
        description="Dual GIP/GLP-1 receptor agonist peptide drug; subject for the Gb dose-response constants",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:194186"],
            OntologySource.KEGG:  ["D11360"],
        },
        force_numeric_id=38,
    )

    # Prednisone -- registered as the subject entity for
    # drug_interactions.py's own glucocorticoid_modifiers() p3/Gb/
    # secretion_scale dose-response constants (MASTER_PLAN 5.1, chunk
    # 20). Both accessions independently verified live via TWO
    # genuinely separate sources: KEGG's own COMPOUND entry
    # (https://rest.kegg.jp/get/C07370) lists "DBLINKS ChEBI: 8382"
    # directly; ChEBI:8382's own page (backend JSON API) reciprocally
    # lists "KEGG COMPOUND C07370" under its own
    # database_accessions.MANUAL_X_REF section, name confirmed exactly
    # "prednisone". Deliberately did NOT also register KEGG's separate
    # D00473 DRUG entry -- unlike metformin/caffeine/tirzepatide, no
    # "Same as" cross-link exists between C07370 and D00473 in KEGG's
    # own data (checked directly), so asserting they're the same
    # substance here would go beyond what either database's own record
    # verifies.
    r.register(
        name="prednisone",
        entity_type=EntityType.METABOLITE,
        description="Synthetic glucocorticoid; subject for the p3/Gb/secretion_scale dose-response constants",
        cross_references={
            OntologySource.CHEBI: ["CHEBI:8382"],
            OntologySource.KEGG:  ["C07370"],
        },
        force_numeric_id=39,
    )

    # Glibenclamide/glipizide -- MASTER_PLAN 5.1, chunk 21. A
    # documentation-completeness gap, not a new constant: both drugs'
    # own secretion-fraction constants (SULFONYLUREA_SECRETION_FRACTION_
    # OF_KE_IB, GLIPIZIDE_SECRETION_FRACTION_OF_KE_IB) already live in
    # ParameterStore under the "insulin" subject (chunks 7 and the Track
    # 1 row-1.3b glipizide extension) -- correctly, since the quantity
    # is fundamentally a fraction of insulin's OWN clearance kinetics
    # (ke*Ib), not an inherent chemical property of either drug, so it
    # stays there (NOT moved or duplicated here). Neither drug itself
    # had ever been given its own canonical ChEBI/KEGG-cross-referenced
    # entity, unlike every other real drug in this campaign.
    #
    # Glibenclamide: KEGG's own D00336 DRUG entry ("Same as: C07022")
    # links to its own C07022 COMPOUND entry, whose DBLINKS lists
    # "ChEBI: 5441"; ChEBI:5441's own backend API reciprocally lists
    # both "KEGG COMPOUND C07022" and "KEGG DRUG D00336" -- name
    # "glyburide" (the USAN name for the same drug INN-named
    # glibenclamide, confirmed via KEGG's own D00336 entry listing both
    # names together: "Glyburide (USP); Glibenclamide (JP19/INN)").
    r.register(
        name="glibenclamide",
        entity_type=EntityType.METABOLITE,
        description=(
            "Sulfonylurea (KATP-channel-closing insulin secretagogue); "
            "its own glucose-independent-secretion constant is "
            "registered under the 'insulin' subject instead (a fraction "
            "of insulin's own clearance kinetics, not this drug's own "
            "chemical property)"
        ),
        cross_references={
            OntologySource.CHEBI: ["CHEBI:5441"],
            OntologySource.KEGG:  ["C07022"],
        },
        force_numeric_id=40,
    )

    # Glipizide: KEGG's own D00335 DRUG entry lists "DBLINKS ChEBI: 5384"
    # directly, with NO separate COMPOUND entry and NO "Same as" link --
    # confirmed absent (both via a direct compound-name search and a
    # direct fetch of D00335's own REMARK field), matching prednisone's
    # own precedent of registering only the verified accession rather
    # than asserting an unverified compound<->drug relationship.
    # ChEBI:5384's own backend API reciprocally lists "KEGG DRUG D00335",
    # name confirmed exactly "glipizide".
    r.register(
        name="glipizide",
        entity_type=EntityType.METABOLITE,
        description=(
            "Sulfonylurea (KATP-channel-closing insulin secretagogue); "
            "its own glucose-independent-secretion constant is "
            "registered under the 'insulin' subject instead (same "
            "reasoning as glibenclamide above)"
        ),
        cross_references={
            OntologySource.CHEBI: ["CHEBI:5384"],
            OntologySource.KEGG:  ["D00335"],
        },
        force_numeric_id=41,
    )

    # NEFA (non-esterified/free fatty acids) -- registered as the subject
    # entity for nefa_pathway.py's own basal-clearance/antilipolysis PK
    # constants (MASTER_PLAN 5.1, chunk 24). Deliberately a GENERIC
    # "fatty acid" class entity, same category as the "protein" chunk
    # above -- NEFA is a clinical/physiological aggregate measurement
    # (total circulating free fatty acids), not one specific molecular
    # species, matching this pathway's own "NEFA" terminology exactly.
    # Both accessions independently verified live this session: ChEBI
    # via the OLS API (https://www.ebi.ac.uk/ols4/api/ontologies/chebi/
    # terms?iri=.../CHEBI_35366 -- label "fatty acid", own
    # database_cross_reference lists "kegg.compound:C00162"); KEGG via
    # https://rest.kegg.jp/get/C00162 (NAME "Fatty acid") -- a genuine
    # two-way cross-reference between independent databases, not one
    # source's assertion alone.
    r.register(
        name="NEFA",
        entity_type=EntityType.METABOLITE,
        description=(
            "Non-esterified (free) fatty acids; circulating lipolysis "
            "product, insulin-suppressed, subject for nefa_pathway.py's "
            "own basal-clearance and antilipolysis-potency constants"
        ),
        cross_references={
            OntologySource.CHEBI: ["CHEBI:35366"],
            OntologySource.KEGG:  ["C00162"],
        },
        force_numeric_id=42,
    )

    # Hydrochlorothiazide -- registered as the subject entity for
    # drug_interactions.py's own THIAZIDE_P3_FACTOR_AT_TESTED_DOSE
    # constant (MASTER_PLAN 5.1, chunk 26, following Track 1 row 1.34).
    # Both accessions independently verified live this session, genuinely
    # reciprocal: ChEBI's OLS API (CHEBI:5778, synonym list explicitly
    # includes "hydrochlorothiazide", own cross-reference lists
    # "kegg.drug:D00340") and KEGG's own D00340 entry (`curl
    # https://rest.kegg.jp/get/D00340`, NAME "Hydrochlorothiazide (JP19/
    # USP/INN)", own DBLINKS field lists "ChEBI: 5778") -- a genuine
    # two-way cross-reference between independent databases.
    r.register(
        name="hydrochlorothiazide",
        entity_type=EntityType.METABOLITE,
        description=(
            "Thiazide diuretic; subject for its own chronic p3 "
            "(insulin action) dose-response constant"
        ),
        cross_references={
            OntologySource.CHEBI: ["CHEBI:5778"],
            OntologySource.KEGG:  ["D00340"],
        },
        force_numeric_id=43,
    )

    # Indinavir -- registered as the subject entity for
    # drug_interactions.py's own INDINAVIR_P3_FACTOR_AT_TESTED_DOSE
    # constant (MASTER_PLAN 5.1, chunk 27, following Track 1 row 1.33).
    # Existing registry checked first (grepped for "indinavir" -- no
    # prior entry, unlike chunk 25's own near-miss). Both accessions
    # independently verified live this session, genuinely reciprocal:
    # ChEBI's OLS search API (CHEBI:44032, label "indinavir" exactly,
    # own cross-reference lists "kegg.compound:C07051") and KEGG's own
    # C07051 entry (`curl https://rest.kegg.jp/get/C07051`, NAME
    # "Indinavir", own DBLINKS field lists "ChEBI: 44032") -- a genuine
    # two-way cross-reference between independent databases.
    r.register(
        name="indinavir",
        entity_type=EntityType.METABOLITE,
        description=(
            "HIV protease inhibitor; subject for its own acute p3 "
            "(insulin action) dose-response constant"
        ),
        cross_references={
            OntologySource.CHEBI: ["CHEBI:44032"],
            OntologySource.KEGG:  ["C07051"],
        },
        force_numeric_id=44,
    )

    # Atorvastatin -- registered as the subject entity for
    # drug_interactions.py's own ATORVASTATIN_P3_FACTOR_AT_DOSE dose-
    # response constants (MASTER_PLAN 5.1, chunk 28, following Track 1
    # row 1.41). Existing registry checked first (grepped for
    # "atorvastatin" -- no prior entry, before writing anything).
    # Both accessions independently verified live this session,
    # genuinely reciprocal: ChEBI's OLS search API (CHEBI:39548, label
    # "atorvastatin" exactly, own cross-reference lists
    # "kegg.compound:C06834") and KEGG's own C06834 entry (`curl
    # https://rest.kegg.jp/get/C06834`, NAME "Atorvastatin", own
    # DBLINKS field lists "ChEBI: 39548", plus its own REMARK line
    # "Same as: D07474" confirming the COMPOUND/DRUG cross-link) -- a
    # genuine two-way cross-reference between independent databases.
    r.register(
        name="atorvastatin",
        entity_type=EntityType.METABOLITE,
        description=(
            "HMG-CoA reductase inhibitor (statin); subject for its own "
            "chronic p3 (insulin action) dose-response constants"
        ),
        cross_references={
            OntologySource.CHEBI: ["CHEBI:39548"],
            OntologySource.KEGG:  ["C06834"],
        },
        force_numeric_id=45,
    )

    # Alcohol (ethanol) -- registered as the subject entity for
    # alcohol_hepatic_pathway.py's own
    # ALCOHOL_EGP_REDUCTION_FRACTION_AT_TESTED_DOSE constant (MASTER_PLAN
    # 5.1, chunk 36, following Track 1's alcohol/AlcoholHepaticPathway
    # work). Existing registry checked first (grepped for
    # "alcohol"/"ethanol" -- no prior entry, before writing anything).
    # Both accessions independently verified live this session,
    # genuinely reciprocal: KEGG's own C00469 entry (`curl
    # https://rest.kegg.jp/get/cpd:C00469`, NAME "Ethanol", own DBLINKS
    # field lists "ChEBI: 16236") and ChEBI's own OLS API record
    # (CHEBI:16236, own cross-reference lists "kegg.compound:C00469"
    # AND "kegg.drug:D00068"/"kegg.drug:D06542", consistent with the
    # KEGG compound entry's own "Same as: D00068 D02798 D04855 D06542"
    # REMARK line) -- a genuine two-way cross-reference between
    # independent databases.
    r.register(
        name="alcohol",
        entity_type=EntityType.METABOLITE,
        description=(
            "Ethanol; subject for its own acute hepatic-glucose-"
            "production (EGP) suppression constant"
        ),
        cross_references={
            OntologySource.CHEBI: ["CHEBI:16236"],
            OntologySource.KEGG:  ["C00469"],
        },
        force_numeric_id=46,
    )

    # Hydroxychloroquine -- registered as the subject entity for
    # drug_interactions.py's own HYDROXYCHLOROQUINE_P3_FACTOR (MASTER_PLAN
    # 5.1, chunk 37, following Track 1 row 1.44). Existing registry checked
    # first (grepped for "hydroxychloroquine" -- no prior entry). Both
    # accessions verified live, genuinely reciprocal: ChEBI's OLS record
    # (CHEBI:5801, label "hydroxychloroquine", own cross-reference lists
    # "kegg.compound:C07043" and "kegg.drug:D08050") and KEGG's own C07043
    # and D08050 entries (`curl https://rest.kegg.jp/get/C07043`, NAME
    # "Hydroxychloroquine"), each DBLINKS field listing "ChEBI: 5801".
    r.register(
        name="hydroxychloroquine",
        entity_type=EntityType.METABOLITE,
        description=(
            "Antimalarial/DMARD; subject for its own chronic p3 (insulin "
            "action) constant"
        ),
        cross_references={
            OntologySource.CHEBI: ["CHEBI:5801"],
            OntologySource.KEGG:  ["C07043"],
        },
        force_numeric_id=47,
    )


# ---------------------------------------------------------------------------
# Pathways
# ---------------------------------------------------------------------------

def _seed_pathways(r: IDRegistry) -> None:
    # Glycolysis
    r.register(
        name="glycolysis",
        entity_type=EntityType.PATHWAY,
        description="Cytoplasmic conversion of glucose to pyruvate",
        cross_references={
            OntologySource.REACTOME: ["R-HSA-70171"],
            OntologySource.KEGG:     ["hsa00010"],
        },
        force_numeric_id=1,
    )

    # Insulin signalling
    r.register(
        name="insulin signalling pathway",
        entity_type=EntityType.PATHWAY,
        description="Cellular response to insulin via IRS-1/PI3K/AKT axis",
        cross_references={
            OntologySource.REACTOME: ["R-HSA-74752"],
            OntologySource.KEGG:     ["hsa04910"],
        },
        force_numeric_id=2,
    )

    # mTORC1 signalling
    r.register(
        name="mTORC1 signalling",
        entity_type=EntityType.PATHWAY,
        description="Nutrient and growth factor sensing via mTOR complex 1",
        cross_references={
            OntologySource.REACTOME: ["R-HSA-165159"],
            OntologySource.KEGG:     ["hsa04150"],
        },
        force_numeric_id=3,
    )

    # IGF-1 signalling
    r.register(
        name="IGF-1 signalling pathway",
        entity_type=EntityType.PATHWAY,
        description="Growth and survival signalling via IGF-1 receptor",
        cross_references={
            OntologySource.REACTOME: ["R-HSA-9006931"],
            OntologySource.KEGG:     ["hsa04014"],
        },
        force_numeric_id=4,
    )
