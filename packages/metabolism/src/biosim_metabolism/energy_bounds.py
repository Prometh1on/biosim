# SAVE TO: packages\metabolism\src\biosim_metabolism\energy_bounds.py
"""
Literature-anchored energy-side and carbon-precursor exchange bounds
for the fasted hepatic network (Track 4, MASTER_PLAN 4.3c).

`fasted_hepatic_network.py`'s own scoping pass found that with only
carbon-precursor bounds set, FVA's answer is a trivial arithmetic sum
of those bounds -- a genuinely stoichiometry-driven answer needs the
network's energy side (fatty acid availability, O2 consumption, ATP
maintenance) bounded too, and NO literature anchor for those three
existed at the time. This module holds what a dedicated follow-up
literature search found for exactly that gap -- plus the
already-existing (less confidently anchored, explicitly flagged as
such) carbon-precursor bounds from `fasted_hepatic_network.py`'s own
docstring, gathered here as named constants so `run_egp_scenarios()`
doesn't repeat undocumented magic numbers.

All values verified live via `PubMedFetcher.fetch_by_pubmed_id()` this
session before being written down.

--- Energy side (the genuinely new literature for this chunk) ---

EX_o2 (O2 consumption, umol/kg/min):
    Byrne NM, Hills AP, Hunter GR, Weinsier RL, Schutz Y. "Metabolic
    equivalent: one size does not fit all." J Appl Physiol.
    2005;99(3):1112-1119. PMID:15831804. n=769 (642 women, 127 men,
    18-74y), resting VO2 measured by indirect calorimetry (ventilated
    hood): "2.6 +/- 0.4 ml O2 x kg(-1) x min(-1)" -- a large,
    contemporary, whole-body measurement, explicitly shown by this same
    paper to be ~35% LOWER than the classic "1 MET = 3.5 mL O2/kg/min"
    textbook convention (population here skews toward higher body fat,
    which independently lowers per-kg VO2 -- disclosed, not silently
    reconciled). Converted via the ideal gas law (22.4 L/mol at STP):
    2.6 mL/kg/min = 116.1 umol O2/kg/min.

    Brundin T, Wahren J. "Influence of protein ingestion on human
    splanchnic and whole-body oxygen consumption, blood flow, and
    blood temperature." Metabolism. 1994;43(5):626-632. PMID:8177052.
    n=10 healthy subjects, hepatic-vein catheterization (Fick
    principle) + indirect calorimetry, fasted baseline: whole-body O2
    uptake "272 +/- 11 mL/min," splanchnic (liver+gut+spleen+pancreas,
    NOT liver alone -- the closest real regional measurement available
    via this classic method, but an overestimate of pure-hepatic O2
    consumption) O2 uptake "62 +/- 5 mL/min" (~22.8% of whole-body).
    This paper reports absolute mL/min, not per-kg -- converting to a
    per-kg rate requires an assumed body weight, since the paper
    doesn't report the 10 subjects' own weight; this platform's own
    `DEFAULT_BODY_WEIGHT_KG=70.0` is used here as an ILLUSTRATIVE
    assumption, not a value from the source. 62/70 = 0.886 mL/kg/min =
    39.5 umol O2/kg/min (splanchnic-specific); 272/70 = 3.886 mL/kg/min
    = 173.5 umol O2/kg/min (whole-body, notably higher than Byrne
    2005's contemporary, much-larger-sample value -- real inter-study
    variability from a smaller, older cohort and an assumed rather than
    measured body weight, disclosed rather than silently reconciled).

EX_fa (fatty acid/palmitate availability, umol/kg/min):
    Bisschop PH, Ackermans MT, Endert E, et al. "The effect of
    carbohydrate and fat variation in euenergetic diets on
    postabsorptive free fatty acid release." Br J Nutr.
    2002;87(6):555-559. PMID:12067425. n=6 healthy men, overnight fast,
    [2,2-2H2]palmitate tracer, control (euenergetic, habitual) diet:
    "Ra palmitate ... 1.47 (se 0.15) micromol/kg per min" -- already in
    the network's own target unit, no conversion needed. Whole-body
    (systemic rate of appearance), not hepatic-specific -- no
    hepatic-specific fatty-acid-uptake number was found this session
    (a fractional-hepatic-extraction search came back empty), so this
    whole-body value is used for BOTH the whole-body and
    hepatic-specific energy scenarios below -- a real, disclosed gap
    rather than an invented proportional split.

ATPM (basal ATP maintenance turnover, umol/kg/min -- DERIVED, not a
    single directly-measured primary source):
    No primary source reporting a whole-body or hepatic basal ATP
    turnover rate directly in umol/kg/min was found this session (a
    genuine, reportable negative result for this specific sub-question
    -- the 31P-saturation-transfer literature that does exist measures
    resting SKELETAL-MUSCLE intracellular Pi->ATP flux in mM/s, a
    different tissue and unit convention entirely). Derived instead via
    a stacked, explicitly-disclosed estimate: VO2 (umol O2/kg/min) x
    0.72 (the fraction of standard-state O2 consumption coupled to ATP
    synthesis, per Rolfe DF, Brown GC. "Cellular energy utilization and
    molecular origin of standard metabolic rate in mammals." Physiol
    Rev. 1997;77(3):731-758. PMID:9234964 -- "approximately 90% of
    mammalian oxygen consumption in the standard state is
    mitochondrial, of which approximately 20% is uncoupled by the
    mitochondrial proton leak and 80% is coupled to ATP synthesis," i.e.
    0.90*0.80=0.72) x 5 (ATP per O2 molecule at a P/O ratio of 2.5,
    standard textbook biochemistry, already the same convention this
    network's own ETC_NADH/ETC_FADH2 reactions use). This is a THREE-
    SOURCE STACKED DERIVATION (a measured VO2 from one population, a
    coupling fraction from an unrelated cross-mammalian review, and a
    textbook conversion ratio), not a single validated measurement --
    presented and used as such, never as if it were a primary source.

--- Carbon precursor side (carried over from fasted_hepatic_network.py's
    own docstring, given names here rather than left as unlabeled
    scratch numbers) ---

    Meyer C, Stumvoll M, Dostou J, et al. "Renal substrate exchange and
    gluconeogenesis in normal postabsorptive humans." Am J Physiol
    Endocrinol Metab. 2002;282(2):E428-434. PMID:11788376. Systemic
    lactate gluconeogenesis "1.97 +/- 0.12 umol/kg/min ... about
    threefold greater than that from glycerol, glutamine, and alanine"
    (each individually) -- rounded up for a permissive bound (not the
    tightest possible reading) and the SAME EXPLICIT UNITS AMBIGUITY
    already flagged in fasted_hepatic_network.py's own docstring
    (glucose-equivalent vs. substrate-mass units are not disambiguated
    in the abstract) applies here too, unresolved.
"""

# --- Energy side ---
EX_O2_WHOLE_BODY_UMOL_KG_MIN = 116.1          # Byrne et al. 2005, PMID:15831804
EX_O2_HEPATIC_SPECIFIC_UMOL_KG_MIN = 39.5     # Brundin & Wahren 1994, PMID:8177052 (illustrative 70kg)
EX_FA_WHOLE_BODY_UMOL_KG_MIN = 1.47           # Bisschop et al. 2002, PMID:12067425
ATPM_WHOLE_BODY_DERIVED_UMOL_KG_MIN = EX_O2_WHOLE_BODY_UMOL_KG_MIN * 0.72 * 5.0
ATPM_HEPATIC_SPECIFIC_DERIVED_UMOL_KG_MIN = EX_O2_HEPATIC_SPECIFIC_UMOL_KG_MIN * 0.72 * 5.0

# --- Carbon precursor side (permissive rounding of Meyer et al. 2002's
# "about threefold" ratio, same explicit units-ambiguity caveat as
# fasted_hepatic_network.py's own docstring) ---
EX_LAC_MEYER2002_UMOL_KG_MIN = 4.0
EX_GLYC_MEYER2002_UMOL_KG_MIN = 2.5
EX_ALA_MEYER2002_UMOL_KG_MIN = 2.5
