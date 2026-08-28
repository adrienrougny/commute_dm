# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`commute_dm` is a research analysis pipeline (the "COMMUTE" disease-map project). It takes
CellDesigner disease maps for COVID-19 and Parkinson's disease (PD), plus BEL knowledge
graphs (AD/PD/COVID), loads them into a **Neo4j** graph DB via the `momapy_kb` stack, finds
the **interface** of proteins shared across the collections, then runs **gene-set enrichment
(GOAT)** and intersection analyses against bulk RNA-seq DE results, and renders new
CellDesigner maps of the interface subgraphs.

**The knowledge graphs are read from their `.bel` files**, not from a database dump, and
**the AD one is an ordinary CellDesigner collection.** `1_10` turns its influence-graph
projection into one CellDesigner file, next to the two disease maps' activity-flow files, and
`2_00` stores that file as `AD_KG_CD_AF`, so every analysis sees three activity-flow
CellDesigner collections — `COVID_DM_CD_AF`, `PD_DM_CD_AF`, `AD_KG_CD_AF` — and no BEL. Only
`1_10` (through `bel_projection.py` and `bel2cd.py`) knows the AD side was ever BEL.

The Python package (`src/commute_dm`) is the library; `scripts/main_analysis/*.ipynb` are the
numbered, ordered driver notebooks. Outputs land under `data/`.

## Environment & tooling

- Managed with **`uv`** (`uv.lock`, `pyproject.toml`); Python ≥3.10.
- Five core deps are **local editable installs** at `/home/rougny/code/{momapy,
  momapy_kb,fieldz_kb,pylpg,momapy_bel}` (see `[tool.uv.sources]`). These are the heart of the
  system — `momapy` (disease-map data model + CellDesigner read/write + SKIA rendering),
  `momapy_kb` (Neo4j-backed LPG sessions over momapy objects), `fieldz_kb`/`pylpg` (the LPG
  node-class machinery), `momapy_bel` (the BEL data model, its `BELReader` — which is how the
  knowledge graphs enter the pipeline — and its `BELWriter`). `momapy_bel` is a **runtime**
  dependency now, not a dev one, and it must be the version whose abundances take a tuple of
  `variants` (two AD proteins carry two `var()`s). These paths have gone missing in past checkouts (the full
  `momapy_kb.lpg.backends.neo4j` import would then fail) — if that recurs, query Neo4j directly
  with the `neo4j` driver as a fallback (DB is `bolt://localhost:7687`).
- **These editable installs freeze their metadata at install time**, and two failures follow from
  it. `momapy.io.core` finds the BEL reader through the `momapy.readers` **entry point**
  `momapy_bel` declares, so an install predating that declaration leaves an empty
  `entry_points.txt` in `momapy_bel-*.dist-info` and every `.bel` read fails with "could not find
  a suitable registered reader" — `import momapy_bel.io.bel` registers it by hand if reinstalling
  is not an option. And a `dynamic` version taken from git tags is cached by uv against the
  **mtime of `src/`** only (`commit: null, tags: null` in the dist-info's `uv_cache.json`), so a
  new tag that touches no file is invisible: pylpg kept resolving as 0.2.4 against `fieldz_kb`'s
  `pylpg[neo4j]>=0.3.0` with `v0.3.0` tagged and `uv sync` refusing to resolve at all.
  `uv sync --refresh-package <name>` is the way out; `cache-keys = [{ git = { commit = true, tags
  = true } }]` under `[tool.uv]` in the dependency's own `pyproject.toml` is the fix.
- **Interpreters: Claude uses `.venv-sandbox`, the user uses `.venv`.** When *you* (Claude) run
  anything, use `.venv-sandbox/bin/python` — it has a working interpreter where the full
  `momapy_kb` + `pylpg` stack imports cleanly. The user's own environment is `.venv` (referenced
  by `.pdm-python`); don't repoint or "fix" `.venv` to match the sandbox. Verify the interpreter
  works before running anything.
- `gea.py` imports **`rpy2`** and drives the R package **`goat`** (auto-installed from CRAN at
  runtime), so GOAT analyses need a working R. `rpy2` **is** an active dependency
  (`pyproject.toml`), pinned by `uv.lock` at **3.6.7**, which requires **R ≥ 4.5.0**: it declares
  and binds `R_getVar` / `R_ParentEnv` with no version gate (`rinterface_lib/_rinterface_capi.py`),
  both added in R 4.5. Against an older R the API-mode extension fails on `R_ParentEnv`, rpy2 falls
  back to ABI mode, and that dies with `symbol 'R_getVar' not found in library 'libR.so'`. Both
  machines are fine now: the user's R is 4.6.1, and **Claude's sandbox (a Debian 12 container) has
  R 4.6.1 too** — it used to be 4.2.2, which is why an older note said GOAT could not run there.
  The R version is the thing to check, not the venv.
  Installing `goat` in the sandbox is the one wrinkle. There is **no root and no `sudo`**, and the
  `-dev` headers `systemfonts` needs are absent, so a CRAN *source* install fails all the way up
  the chain (`systemfonts` → `ggforce`/`sass` → `ggraph`/`shiny`/`treemap` → `goat`); the runtime
  `libfontconfig.so.1` / `libfreetype.so.6` *are* present, so **Posit's binary repo installs
  cleanly**:
  `install.packages("goat", repos="https://packagemanager.posit.co/cran/__linux__/bookworm/latest",
  lib="/tmp/claudehome/.local/lib/R/library")` with `options(HTTPUserAgent=...)` set so P3M serves
  binaries. That library is `.libPaths()[1]`, which is what `gea._ensure_r_packages` writes to.
- There is no test suite (only an empty `tests/__init__.py`) and no lint/build config beyond
  hatchling packaging. "Running" the project means executing the notebooks in order.

## Credentials

Notebooks `import credentials` from `scripts/main_analysis/credentials.py` (git-ignored,
already present). It holds `NEO4J_URI/USERNAME/PASSWORD/DATABASE` and MINERVA creds. The DB
is shared mutable state — `2_0` calls `session.delete_all()` and rebuilds it.

## Pipeline (run notebooks in numeric order, from `scripts/main_analysis/`)

Notebooks are numbered `<stage>_<step>`, with **`<step>` gapped by 10** (`_00`, `_10`, `_20`, …)
so a new notebook can be inserted between two others (e.g. `4_05`) without renumbering the rest.

Every notebook starts with `%store -r`, which restores path constants saved by
`0_00` via IPython's `%store`. So **`0_00_make_paths` must be run once first** to populate the
store; afterwards any notebook can be opened standalone.

- `0_00` — define and store all `*_DATA_DIR` / `*_FILE` path constants (rooted at `data/`); every
  other notebook restores them with `%store -r`.
- `1_00_make_bel_files` — prepare BEL KG inputs.
- `1_10_make_af_files` — the **three activity-flow files**: the COVID and PD disease maps through
  `pd2af`, and the AD knowledge graph through `bel_projection.py` + `bel2cd.py`. The AD cell reads
  `data/ad_kg/bel/ad_kg.bel`, projects it, builds the map with the HGNC annotations, renumbers and
  writes one ~16 MB file to `AD_KG_CD_AF_BUILD_DIR`. No database is involved, which is why it lives
  in stage 1.
- `2_00_save_collections` — wipe the DB, then **one call per return type**: the five CellDesigner
  collections (`COVID_DM_CD`, `PD_DM_CD`, `COVID_DM_CD_AF`, `PD_DM_CD_AF`, `AD_KG_CD_AF`) as maps,
  and the four BEL knowledge graphs (`AD_KG_BEL`, `PD_KG_BEL`, `COVID_KG_BEL`, `CBM_KG_BEL`) as
  models, read straight from their `.bel` files. The two cells are symmetric — a list, then the
  call, identical but for `return_type` — and both use `with_membership_edges=True` (the
  `HAS_MEMBER_MODEL_ELEMENT` edge that `get_collections_for_nodes` relies on) and
  `integration_mode="hash"`. **No `object_key_to_node` is passed**, and that is deliberate: a call
  given none builds its own (`session.py:132`), so integration *within* a call is untouched, and
  the two calls have nothing to integrate *across* — under hash integration the key is the object
  itself, and the AD map's 26925 objects and the AD KG's 12954 have an intersection of exactly 0,
  a `GenericProtein` never being equal to a `ProteinAbundance`. Sharing a dict looked like it did
  something and did not. **`AD_KG_CD_AF` is saved in the same call as the other maps**, which is
  what makes the elements it shares with them one node — that is the integration that matters, and
  it is the within-call one.
  There is no `cypher-shell` import any more, and no annotation pass: the AD map carries its own
  HGNC annotations, written by the transformation (see `bel2cd.py`).
- `2_05_get_collections_statistics` — descriptive statistics of every collection. Its BEL cells
  **read the four stored collections back out of the database** — one query per collection
  returning its `BELModel` nodes, rebuilt by `session.execute_query_as_objects` (CBM's 59 entries
  give 59 models, merged into one), through one `node_id_to_object` shared by the four calls so a
  node two collections share gives one object. `commute_dm.bel_lpg` must be imported for the
  rebuilding to work at all. They then count on the models: node and relationship
  breakdowns are `type(element).__name__` and `type(statement).__name__` over the model, and the
  relation families are tuples of statement classes in the parameter cell. It builds its influence
  graphs by calling `bel_projection.get_bel_influence_graph_projection` per model — the same
  function `1_10` calls — and wrapping the returned `BELModel` in a `networkx.MultiDiGraph`:
  **one node per momapy_bel element, one edge per relation statement**, labelled with
  `bel_projection.get_element_bel_string(element)` and `type(element).__name__` (which is what the
  phenotype test compares against). `get_bel_isolated_node_reason_breakdown` takes the models
  rather than a session: an isolated node's reason is which statements the element does appear in,
  which the whole model answers.
  Statements and projection nodes/edges, AD / CBM / COVID / PD: 8542, 3987/4911; 3130, 1920/1900;
  3479, 2285/1690; 3128, 1881/1760. `2_00` stores exactly the `.bel` files, and the hydrated
  models give what parsing those files gives, to the number — **measured, since nothing had ever
  read a `BELModel` back out of the database before**. Hydrating is also no slower: 4.3 s for the
  four collections against about 7 s to parse the 62 files.
  **Against the cypher dumps this has now been checked term by term, not by counting**, since a
  count can agree while the content does not. The dump carries each node's BEL string, so the
  comparison is between *sets of causal edges written as BEL triples* — the dump restricted to the
  same causal relations (`INCREASES`, `DIRECTLY_INCREASES`, `DECREASES`, `DIRECTLY_DECREASES`,
  `REGULATES`, `TRANSLATED_TO`) and the same projected classes, both sides canonicalised for
  quoting, whitespace and the order of a complex's members and a protein's modifiers, which BEL
  does not fix. The result, dump vs `.bel`:
  - **AD: identical. 4910 = 4910 edges, 2788 = 2788 nodes, and the same 776 `p(HGNC:…)` symbols,
    with an empty difference in both directions.** This is the one that matters for everything
    downstream, and it settles it: the AD substrate did not change.
  - **COVID: 1690 = 1690, one edge apart** — a `composite()` naming the same member twice, which
    `frozenset[Species]` collapses (the documented homomultimer limitation).
  - **PD: 1760 = 1760, nine edges apart**, of three kinds. Three are `gmod(TestNS:TestName)`
    against `gmod(Me)`, and here **the dump is the wrong one**: the file says `gmod(Me)` (all 36 of
    its `gmod`s do) and pybel let a placeholder namespace through. Four are `pmod(Ph,Y)` against
    `pmod(Ph,Tyr)` — the *file itself* writes both spellings (`Tyr` for MAPK9/MAPK14/IGF1R/IRS1/
    INSR, `Y` for Grin2a/Grin2b), pybel canonicalises three-letter to one-letter and `momapy_bel`
    keeps what is written; same residue, and **verbatim reading splits nothing**: 0 of the 1315
    nodes are one term written both ways. The last two are `complex(p(HGNC:PARK7),p(HGNC:PARK7))`
    against `complex(p(HGNC:PARK7))` — a real loss, the same homomultimer collapse as COVID's.
  **Total content actually lost by the new path across all four KGs: 3 edges, all homomultimers**
  (2 PD, 1 COVID). Everything else that differs is rendering, or the dump being wrong.
  - **CBM: 2146 against 1900, the one real difference.** 246 of the 257 dump-only edges are
    supported *only* by four PMIDs that have no `.bel` file at all (22899644, 30342839, 36231087,
    36526429) — the dump covers 63 papers, the directory has 59, and all 59 are in the dump. The
    remaining 11-or-so pairs are curation edits, the same assertion re-namespaced
    (`a(NCIT:Remdesivir)` → `a(CHEBI:remdesivir)`). So CBM's gap is **corpus coverage, not the
    reader**: the four papers are missing from `data/cbm_kg/bel/`.
  The earlier note here compared raw dump totals (3979/4910, 1832/1760, 2255/1690, 2074/2145) and
  read PD as far off; that was an artefact of counting the dump's non-projected nodes (reactions,
  translocations, degradations, lists) rather than a content difference. Node counts in the table
  above still run *higher* than the dumps' because this projection prunes less — it drops a
  complex's members and an activity's subject, while the graph-based one also dropped anything
  reached by pybel's derived base→proteoform edges, so `p(HGNC:"ACE2",loc(MESHA:"Kidney"))` is
  again the top-level term it always was in the document — and because the projection keeps
  isolated elements, which have no edge to appear in.
- `3_00_get_interfaces` — query proteins present across collections (joined on **UniProt RDF
  annotations**), write `data/.../interface/*.json`. The COVID×AD interface is **159** two-way and
  **111** three-way, smaller than the 189/127 of the old BEL pairing because `1_10` drops the
  isolated species on purpose: those identifiers could never seed a walk, so after the drop the
  interface *is* the seedable set. These are the **full** interfaces, subunits included, which is
  the right answer to "which proteins do these maps share". Every notebook now asks for them the
  same way (`with_subunits=True`, the default), so the identifier counts are the same everywhere;
  what the flag changes is which **species** each protein stands on (see `get_interface`'s
  `with_subunits`). The archived rows record that: `species_*` is the species standing for the
  protein and `subunit_*` the protein itself whenever that species is a complex. Three-way:
  **1190 rows, 572 of them through a complex**.
- `4_00_make_goat_gene_lists` — turn raw RNA-seq DE CSVs (`data/rnaseq/`) into GOAT gene-list
  CSVs, mapping Ensembl→Entrez/symbol via the HGNC dataset.
- `4_10_make_interface_graphs` — **both pairings**, in a `PAIRINGS` list of the same shape as
  `4_20`'s: per pairing the collection names, the interface tuple, `max_levels` and the output
  directory. COVID→PD keeps the three-way interface filter (now with `AD_KG_CD_AF`); COVID→AD is
  the two-way `COVID_DM_CD_AF × AD_KG_CD_AF`. Both sides of both pairings are stored CellDesigner
  collections, so `core.load_submap_inputs` has one path and there is no BEL case anywhere.
  `max_levels` stays per pairing, and that is a property of the downstream graph, not of the code:
  `[1,2,3]` for AD against `[2,3,4,5,6]` for PD, because the AD influence graph is much denser.
  Both pairings take their interface **with subunits** (the default), which is what makes a
  protein seedable rather than what makes it unseedable: a protein stands on every species of its
  collection that **is it, or contains it**. A protein a collection holds only inside a complex
  stands on that complex; one that is both free and bound stands on both.
  Current output: **180 maps over 77 proteins** COVID→AD and **375 over 76** COVID→PD, against 95
  over 45 and 190 over 40 when the interface was taken without subunits. Everything gained is a
  protein whose wiring lives on a complex — `NTRK1`, `NTRK2`, `TLR2`, `TLR4`, `NFKB1/2`, `STAT3`,
  `SOCS1/3`, `TGFB1`, `CCL2`, `IFNG`, `NGFR`, `TRAF6`, `TRADD`, `FADD`, `APAF1`, `JUN`, `FOS`, and
  the mitochondrial respiratory chain (`COX4I1`, `COX5A`, `NDUFS1`, `MT-CO2`).
  The free protein and the complex holding it are usually **not** linked by any arrow: on the AD
  map they never are, since BEL states no complex formation — 73 of the 83 COVID→AD pairs of that
  shape have no edge between them at any hop. So standing a protein on its complexes is the only
  way those complexes are reachable, and it is not a shortcut for something the walk would find
  anyway.
  Accepted consequence: a protein that is one of a **large** complex's members gets that whole
  complex's neighbourhood attributed to it. Complexes are 3–5 members at the median but run to 49,
  which is what brings the respiratory chain in. Refusing to seed on complexes above a size would
  drop that group and almost nothing else (66 and 72 proteins instead of 76 and 77); it is one more
  setting to tune and was deliberately left out.
  The binding constraint is `MIN_N_NODES=5` on the **COVID upstream** side. The last cell reads
  every written file back, which is the assertion the identity invariants exist for.
- `4_20_make_interface_goat_analysis` — GOAT enrichment + intersection analyses of interface
  subgraphs vs. the gene lists, for **both pairings and all three modes** (`upstream`,
  `downstream`, `upstream_and_downstream`). Influences come straight from
  `submaps.load_signed_influences` — seconds, since no map is hydrated.
  Runs in Claude's sandbox now that its R is 4.6.1 and `goat` is installed there (see
  "Environment & tooling"); the note that it could not is stale.
  Its interface is taken **with subunits**, as `4_10`'s is. **The two flags are unrelated.**
  `4_20`'s own `WITH_SUBUNITS` stays `True` — it is `gea`'s widening of a selection to its species'
  subunits, what makes a subunit contribute its gene — while `get_interface`'s decides which
  **species** a protein stands on, and so what its walks start from.
  Note the two compose: a seed that is a complex contributes its members' genes through
  `WITH_SUBUNITS`, so standing a protein on a complex widens its gene set by that complex's
  members. That is intended — the complex is what the collection asserts — but it is why every
  GOAT and intersection number moves with this change.
  `4_20`'s protein list must equal `4_10`'s: same interface, same `MIN_N_NODES`, so a divergence
  means one of the two moved.
- `5_00_make_annotations_for_bel` — annotation tooling (see `annotations/` helper scripts).
- `6_00` / `6_10` — immuno-target analysis and graphs.

## Library architecture (`src/commute_dm`)

The package has no `__init__` exports; modules are imported by full path. Data flows
**Neo4j (stored AF maps) → element-id selection → CellDesignerMap / gene sets**.

### `submaps.py` — the comorbidity sub-maps
The AF collections (`COVID_DM_CD_AF`, `PD_DM_CD_AF`, from `pd2af` in `keep-reactions` mode)
already *contain* what an influence-graph round-trip used to reconstruct: an AF map has 0
reactions and its `Modulation`s point straight at `Species`, and the glyphs, arcs and
`LayoutModelMapping` are all stored. There are no derived `_*INFLUENCE` relationships in the DB,
and no `InfluenceGraph` class. The module **walks node ids and assembles stored momapy objects**:
the walk is a small edge-list query plus a BFS, the drawing material is loaded once as objects.
- **The one thing to understand before touching it**: collections are imported with
  `integration_mode="hash"`, so two content-equal elements are *one* DB node; hydrating
  everything through one shared `node_id_to_object` cache makes that one node exactly one Python
  object, across submaps and across collections. Object identity — which `LayoutModelMapping` is
  keyed on, and which the CellDesigner writer resolves aliases by — then holds for free. That is
  why there is no interning, no canonical-element table and no identity repair. Never hydrate the
  same element twice through a different cache (and never use
  `session.cypher_query_as_layout_elements`: no cache parameter → fresh objects → `alias=""`).
- `load_signed_influences(session, collection_names)` → `Influences` (four small queries, ~0.2 s):
  `upstream` / `downstream` (level-bounded BFS), `close_over_gates` (a retained gate pulls in all
  its inputs), `species_only` (drop gates — they carry no annotation and must not reach `gea`).
  Nodes are Species **and** BooleanLogicGates, so one influence through a gate is two hops.
  `SIGNED_MODULATION_CLASSES` is the single source of truth: pylpg labels a node with its class
  *and every ancestor*, so the tuple serves both the Cypher label test and `isinstance`. Matching
  `:Modulation` already excludes every `Unknown*` arc (`UnknownModulation` is a *sibling* of
  `Modulation`); the loader **raises** on a modulation class that is neither signed nor in
  `KNOWN_UNSIGNED_MODULATION_CLASSES` rather than skipping it. One `Influences` covers both
  collections — their influence graphs share no node, so per-collection loading is equivalent.
- `load_collections_as_map(session, names, node_id_to_object)` — every stored map of the
  collections, hydrated and merged by `merge_maps`. **COVID+PD AF: 37 s, 135237 objects, 0.7 GB
  peak; COVID+AD AF: 26 s, 92610 objects** (~3500 objects/s), against the ~3 min this used to take
  before pylpg's batched subgraph prefetch. Only the map path pays it; `4_20` never calls it.
- **`make_submap_from_model_elements`** assembles one sub-map out of frozen momapy objects,
  reusing the stored layout element of each model element and the stored arc of each modulation,
  with `dataclasses.replace` for the few that need a different fill or endpoints. Complex-subunit
  mapping entries are read **forward** out of the stored mapping (without them the writer refers
  to aliases it never writes); compartments are **kept** and closed over `outside` (the reader's
  `get_ordered_compartment_aliases` crashes on a dangling `outside`, not on missing compartment
  glyphs); one `LogicArcLayout` per gate input, without which the writer emits an input-less gate;
  every `id_` is renumbered except the literal `default` compartment. `extra_species` /
  `extra_influences` are the hook for elements with no stored counterpart (the synthetic central
  node; the BEL side comes in through the merged `source_map` instead) and are the only
  synthesized modulations. The modulation list is **deduplicated**, induced first so the *stored*
  arc is reused rather than a synthesised one: `model.modulations` is a frozenset while the loop
  appends one arc per list element, and `extra_influences` can hold duplicates whenever
  `node_id_to_object` is many-to-one — which it is on the BEL side.
- **The four identity invariant checks** (`check_species_template_identity`,
  `check_modification_residue_identity`, `check_subunit_layouts`, `check_compartment_identity`,
  and `check_identity_invariants` running all four) fire on every sub-map just before
  `renumber_ids`. momapy elements are equal by value with `id_` excluded, and both this module and
  the CellDesigner writer rely on *object identity* where a value-equal duplicate goes unnoticed
  until the file is read back: a duplicated template collapses in `species_templates` and the other
  species emit a `<proteinReference>` to an id never written; a duplicated `ModificationResidue`
  makes a `<modification>` name a residue the `<listOfModificationResidues>` does not declare; a
  subunit with no mapped nested glyph gets no alias; a duplicated compartment leaves a species
  pointing at a compartment id never written. They are not BEL-specific — they hold for
  stored-only data for free (hash integration made them true); the BEL side is the one that has to
  *establish* them, in `bel2cd.py`. Verified against the untouched COVID→PD pipeline first,
  before being trusted on BEL data.
- **A selected species that a selected complex also holds as a subunit is promoted to its own
  object**, at the top of `make_submap_from_model_elements` — and this is a fix to a
  **pre-existing, non-BEL defect the checks surfaced**. Hash integration makes such a species *one
  object playing both roles*, and the writer, which indexes subunits by `id()`, then emits neither
  `<species>` nor `<speciesAlias>` for the standalone occurrence while its modulations still point
  at it: the file writes cleanly and fails to read back with `KeyError`. Measured on COVID×AD
  before the fix: 7 pairs over 5 of 101 maps (`CASP1` in `NLRP3 oligomer:PYCARD:CASP1`, `TRAF2` in
  `TRAF2:ERN1:unfolded protein`) — and those 5 were **exactly** the 5 that failed to read back.
  Latent in COVID→PD too; its maps just never select such a pair. The copy is value-*equal*, which
  is what keeps the rest of the function unchanged (`layout_element_of` and `fills` are
  equality-keyed); only `model.species` membership and the modulations' endpoints care about
  identity. Two consequences worth knowing: `get_layout_element_for_model_element` takes
  `preferred_ids` so the standalone glyph is picked over the one nested in the complex, and the
  modulation list carries `(modulation, stored modulation)` pairs, because a rebuilt modulation
  must not lose the stored arc that `get_arc_for_modulation` finds by identity.
  (The 201-shared-subunits note below is the same phenomenon from the other side, and its "none is
  a modulation endpoint" qualifier does not cover this: such a node enters a selection as a gate
  input or as an interface seed.)
- `make_fitted_synthetic_layout` / `_measure` / `_wrap_label_text` live here, not in
  `bel2cd`: `core` needs the first for every sub-map's synthetic central node, in both
  pairings, and `bel2cd` is export-only. `bel2cd` imports them back; acyclic.
- `iter_subunits` / `iter_species_and_subunits` — `model.species` is top-level only, and a BEL
  sub-map keeps half its modifications and templates on complex members.
- `LayoutModelMapping`'s forward dict is **equality-keyed** (so `count_mapping_entries` is a real
  guard: two content-equal layout keys collapse and one glyph→element entry is lost) while its
  inverse is identity-keyed.
- Gotchas worth keeping in mind: `Node.position` / `width` / `height` have **no defaults**, so a
  layout built directly as a frozen object must be given geometry; and
  `pd2af…_MODULATION_CLASS_TO_LAYOUT_CLASS` has no entry for `Inhibition`, `Catalysis` or
  `PhysicalStimulation`, which `_make_modulation_arc` supplies locally.
- Output geometry: callers must run `pd2af.utils.make_auto_layout(cd_map)` — stored arc geometry
  is not reusable across a merge; it recomputes every position and segment and fits the root
  layout, which is why the assembly does not call `set_fit` itself.

### `bel2cd.py` — a BEL projection as a real activity-flow CellDesigner map
`bel2cd` imports `bel_projection`, never the reverse. It knows nothing about Neo4j: everything is
read off the elements' fields — `modifications`, `variants`, `fragment`, `location`, `members`,
`abundance` — and the only BEL string it uses is the glyph sort key. Its entry point is

    cd_map, element_to_annotations = make_cd_map_from_bel_influence_graph_projection(
        projection, hgnc_symbol_to_annotations
    )

and it gets there in three steps:
- `make_species_fields(elements)` → `{element: fields}`, where `fields` is a **plain dict** —
  `species_class`, `template_class`, `template_name`, `species_name`, `active`,
  `compartment_name`, `residue_states`, `structural_states`. It exists for one reason: a species is
  frozen and its template must already declare every residue *any* proteoform of that protein
  carries (invariant (b)), so nothing can be built until every element has been read. Hence the
  fields carry residue **names**, not residue objects, and a template **name**, not a template.
  Pure, DB-free, memoised on the element itself, and it covers the given elements **and** every
  element reached through their members. The **active set is derived here**, from the activities
  among the elements: an `Activity` takes its subject's fields whole, `ma()` code dropped, so
  `act(X)` and `X` give one species. Measured on all four KGs, every activity subject is reachable
  from a projected activity, so there is nothing to thread in from outside.
  The translation itself: `ProteinAbundance("GDE1")` → `GenericProtein("GDE1")`; `pmod(Ph,Y,18)`
  and `var("P, T, 212")` → badges (whitespace normalised, so `"P, S, 9" == "P,S,9"`);
  `var("K,670,N")` → `StructuralState("K670N")`; a `ProteinModification` **with a namespace** is an
  ontology pmod (`pmod(GO:…)`) → `StructuralState`, one without is a code → badge, so the split is
  a field, not a regex; `fragment` → `TruncatedProtein` + `TruncatedProteinTemplate` + a
  `StructuralState`; `CompositeAbundance` → `Complex` + `StructuralState("composite")`; an unnamed
  complex is named by joining its members' names. Modifiers are routed by the target class's
  *capabilities*; a class with no container for one (a `Gene`, an `Unknown`) folds it into the
  **species** name as `NAME [x|y]` — never into the **template** name, or one protein would declare
  N `<protein>` entries.
- **A `var()` can be a modification, and it is told apart from a substitution by which field holds
  the number.** BEL1-era curation writes a `pmod` inside a `var()`, so `var(P,S,396)` is
  phosphorylation of S396 — `(code, amino acid, position)` — while `var(K,670,N)` is a substitution
  — `(amino acid, position, amino acid)`. The letters cannot decide it (`P` is both proline and
  phosphorylation; `A` alanine and acetylation), only the position of the numeric field can, and
  `_VARIANT_MODIFICATION` is tried before `_VARIANT_SUBSTITUTION`.
  `VARIANT_CODE_TO_MODIFICATION_STATE` holds the spec's `bel1_migration.protein_modifications`
  codes that have a CellDesigner state (`P`/`Ph`, `A`/`Ac`, `M`/`Me`, `U`/`Ub`, `G`/`Glyco`,
  `H`/`Hy`); `F`, `R` and `S` are deliberately absent, which also keeps `S` free to be serine. A
  code outside the table is not read as a modification at all, so it falls through to the
  substitution and verbatim paths exactly as before.
  The evidence that this reading is right, measured over the four KGs: the first field of the 90
  modification-shaped variants is only ever `P`/`Ph`/`A`/`M`/`U`; the residue of a `P`/`Ph` one is
  only ever S, T or Y; `Ph` is not an amino acid at all; and **the same site is written both ways**
  — MAPT S202/S235/S262/S356/S422/T212, SNCA S129, GSK3B S9, GSK3A S21, IRS1 S312 each carry a
  `pmod(Ph,…)` node *and* a `var(P,…)` node. Before the table, the six non-`P` ones (HINFP K5/K12
  and H3F3A K14 acetylation, H3F3A K4 methylation, BACE1 L501/W280 ubiquitination) were drawn as
  the literal text `A,K,12`. Adding them changes no count — 2418/1206/4631/5, 954 annotated, 778
  UniProt, 779 subunits are unchanged — only 6 structural states become 6 modification badges.
- **Nothing is ever invented for a sub-term that says nothing.** A `frag("?")` (range unknown, no
  descriptor) still switches the class to `TruncatedProtein`, which is what says "a fragment of",
  but adds no `StructuralState("?")` on top — the placeholder was noise, and it turned up as
  `<structuralState structuralState="?">` in the output. `var("?")` and `var("p.?")` are the
  opposite case and are **kept verbatim**: there `"?"` is the value BEL states, and it is the only
  thing telling that proteoform apart from the plain protein. Same principle as the nameless
  modification residue below.
- `make_templates(species_fields)` → `{(template class, name): template}` and
  `make_compartments(species_fields, root_name)` → `({location name: compartment}, root)` — the
  shared objects, interned by value, carrying invariants (a) and (d). Nothing is interned against a
  source map: the export is standalone, and identity with the stored collections' elements is the
  *save's* job (hash integration with a seeded `object_key_to_node`). The residue objects are *not*
  tabulated: `make_species` reads them off the template it is about to hand the species, which is
  invariant (b) by construction.
- `make_species` (always a **fresh object tree**) and `make_species_layouts` (measure, then build).

The map assembly, `make_cd_map_from_bel_influence_graph_projection`:
- Species are **interned by value**: two `ma()` codes of one activity describe identically, as do
  whitespace-variant proteoforms — and so do `act(p(X))` and `p(X)`, which is the point. AD:
  3987 projected elements → 3598 species, 2457 after the isolated drop.
- **`act(X)` and `X` are one species.** The projection has no relation between the two, so keeping
  them apart severs 6889 upstream→downstream two-hop paths that three hops cannot recover; the
  price is that X carries an active border whenever BEL asserts any activity of it.
- **Species carrying no modulation are dropped.** Nothing walks or draws them, so they would stand
  stranded. `drop_isolated_species` runs **after** interning and after the modulations are built,
  so a protein isolated in its own right but wired through its activity form survives; a dropped
  species that belongs to a kept complex still appears as that complex's subunit.
- **Every species built is recorded, subunits included**, which is what the annotations need: a
  subunit is a per-occurrence object, and 16 of the interface's UniProt identifiers are carried
  only by complex members.
- `BEL_RELATION_CLASS_TO_MODULATION_CLASS` draws every causal relation, `Regulates` included: it
  is unsigned and becomes CellDesigner's unsigned `Modulation`, which `submaps`'
  `KNOWN_UNSIGNED_MODULATION_CLASSES` already expects. The sub-map walk and
  `make_submap_from_model_elements` both read signed modulations only, so these 144 arcs appear on
  this map and on `2_05`'s statistics and never on a sub-map — but they do keep 39 species alive
  that the signed relations alone would have dropped.
- **Self-loops are dropped here**, where the reason lives: pd2af's arc geometry needs two distinct
  endpoints, and interning creates species-level self-loops (an `act(X) → X` relation) that no
  element-level filter could see.
- `make_hgnc_symbol_to_annotations(hgnc_file_path)` → `{HGNC symbol: frozenset[RDFAnnotation]}`,
  and `_make_element_to_annotations` maps those onto the species. **The annotations are part of
  the transformation**, not of a database pass: the interface joins collections on UniProt
  annotations, so the map has to carry them. A protein gets **several annotations, one per
  identifier** — the identifiers are alternatives, so a symbol with two UniProt accessions means
  one or the other — each with a single resource and the `BQBiol IS` qualifier, all in one bag,
  which is what the stored CellDesigner maps hold for their own species. Four namespaces:
  `uniprot` (what the interface joins on), `ncbigene` (what GOAT matches), `hgnc.symbol` (what the
  intersection matches) and `hgnc`. Keyed by **gene symbol**, because that is the identifier the
  knowledge graph itself uses (`p(HGNC:"MAPT")`), and only HGNC-namespace protein abundances are
  annotated. An activity's subject is followed, a complex's members are not, since each is drawn
  as its own subunit species and annotated there. Accumulated **by value**, so a top-level species
  and a value-equal subunit share one entry, which is what the writer looks up. Result on AD: 961
  annotated species, 780 UniProt ids.
- `Abundance` is split on the **`namespace` field** — `a(CHEBI:…)` → `SimpleMolecule`, anything
  else → `Unknown`: MESH/CONSO/GO abundances are proteins, aggregates and cellular components, not
  molecules. `BEL_CLASS_TO_CD_CLASSES` selects the template class too, so an entry decides whether
  a species declares a `<protein>` / `<gene>` / `<rna>`.
- Accepted consequences: a fragmented gene yields two `<protein>` entries (GENERIC + TRUNCATED),
  which is idiomatic CellDesigner; a complex listing the same member twice loses the multiplicity
  to `frozenset[Species]` (`homomultimer` would be the fix, out of scope).
- No global auto-layout: `make_auto_layout` (graphviz `dot`) does not finish on a graph this size.
  Glyphs get local, non-overlapping positions; the sub-maps of `4_10` lay out their own selections.
- Current AD result: 2457 species, 1217 templates, 4775 modulations (4631 signed + 144
  `Regulates`), 5 compartments, 790 subunits, 780 UniProt ids on 961 annotated species. The
  written file round-trips: read back, it has the same counts and the same 780 UniProt ids.

**The four identity invariants** replace the old BEL-string naming invariant. Getting one wrong
produces a file that writes without error and fails to read back with `KeyError` — the measured
baseline for a naive symbol-naming attempt was 12 of 20 maps.
- **(a) one template object per `(class, name)` within the export.** It is no longer interned
  against a source map: value-equal templates become one database node at save time.
- **(b) `Modification.residue` *is* an object in its species' template's residue container** — so a
  template's residues are the union over *every* proteoform, and the lookup is re-derived from the
  interned template. Note the reader **fills in an "empty" `Modification` for every template
  residue a species does not carry** (`reader.py` ~795): that is CellDesigner's own proteoform
  semantics, not a round-trip defect, and it is why a round-trip comparison must ignore
  `state is None` modifications. A `pmod(Ph)` that names no site gives a residue with
  `name=None`, not a `"?"` — the writer omits the attribute and the reader reads `None` back, so
  the nameless residue survives the round trip. What must never be `None` is the **residue
  itself**: the writer would emit `residue=""` and the reader looks that up unguarded. Sort
  residue names through `bel2cd._residue_sort_key`, since `None` does not compare with `str`.
- **(c) a subunit object is never the same object as a top-level species** — `p(HGNC:"APP")` is
  both a projected element and a member of several complexes, so memoising species by element (the
  obvious way to get value-collapse) would make APP's own glyph disappear. Interning happens only
  at the top level, on the finished object, in
  `make_cd_map_from_bel_influence_graph_projection`.
  **This one is an in-memory property of the export and does not survive the save**: a top-level
  `p(X)` and the same `p(X)` inside a complex are value-equal, so hash integration stores them as
  one node and hydration returns one object. A `source_map` built from `AD_KG_CD_AF` therefore has
  651 subunits that *are* their top-level species, and that is expected —
  `make_submap_from_model_elements` promotes such a species to its own object per sub-map, before
  the checks run. **Never run `check_identity_invariants` on a source map**; it is a pre-write
  check for sub-maps.
- **(d) one compartment object per location, hanging off an undrawn root** (named after the
  collections the projection was loaded from). The root is **not decorative**: `pd2af.utils._build_dot_graph` (`utils.py:339-346`) only
  attaches a compartment's dot cluster to the graph when `compartment.outside is not None` — there
  is no `else` — so a drawn compartment with `outside=None` gets an orphan cluster, its species
  never reach graphviz, and they keep their throwaway positions. (A real pd2af bug, worth
  reporting upstream; worked around here, so no pd2af change is needed.) The root also keeps a BEL
  `nucleus` (`outside=BEL`, the root's name) value-distinct from a stored one
  (`outside=default`), and
  therefore the species inside them too. 18 of the 84 BEL `loc()` names already exist as stored
  compartment names, so one sub-map can show two same-labelled boxes from different sources —
  that is the price of keeping them distinct.

Layout notes:
- **Badge angles are residue-derived, not species-local.** CellDesigner stores the angle on the
  *template*, and `_writing.find_residue_angle` takes it from the **first** species using that
  template with a modification on that residue, so a species-local index would move every other
  proteoform's badge on read-back. The angle is `2π · residue.order / len(template residues)`, and
  the reader's exact forward transform is reproduced so `_writing.compute_cd_angle` inverts it
  (verified to ~1e-6 rad over the AD badges, of which the map now has 212).
- A subunit alias is emitted **only** when `mapping.get_child_layout_elements(subunit, complex)`
  returns a layout that is *also* in `complex_layout.layout_elements` — so each subunit needs both
  its own nested layout element and a mapping entry. (This refutes the old claim that "a subunit's
  sub-glyph does not survive a CellDesigner round trip".)
- **Glyph positions are throwaway but not arbitrary.** A badge is a 16×16 white square whose only
  distinguishing content is its position and one-character label, so two badges of two species at
  nearby positions are *value-equal* and the second `add_mapping` evicts the first from the
  equality-keyed forward dict. `make_species_layouts` therefore spreads the species along x so
  their `[center ± width/2]` intervals are disjoint. Measured: with pd2af's synthetic 1pt spacing,
  one AD modification silently lost its mapping entry.
- Glyphs are sized to their (wrapped) content: a leaf to its label, a complex to its stacked
  subunits plus a label band, and both grown by a badge margin when decorated, never below the
  layout class's default. This reaches the output because `make_auto_layout` *preserves* the size
  it is given (`_build_dot_graph` sets each dot node's size from the layout element,
  `_reposition_from_dot` resizes only compartments), and nested geometry survives it because only
  top-level elements become dot nodes and the subtree is translated whole.
  `core`'s synthetic central node still uses `make_fitted_synthetic_layout` (an `Unknown`'s
  default is 60×30, too small for some display names).
- `renumber_ids` and `make_auto_layout` both round-trip through `momapy.builder`, rebuilding every
  object, so **any dict keyed on pre-round-trip objects is stale afterwards** — which is why the
  nested geometry is built here rather than by calling `momapy.celldesigner.utils`'
  `set_complexes_to_fit_content` / `set_modifications_to_borders`.
- Gene/RNA modification badges are **deferred**: their residues live in `template.regions` as
  `ModificationSite`, which the writer emits differently. `TEMPLATE_RESIDUE_FIELDS` is the
  descriptor table that makes enabling them a one-line change. (AD has zero gene/RNA modifiers;
  PD's 26 `gmod()`s are the only ones.) The deferral is why `_read_species_fields` reads
  **`momapy_bel.core.ProteinModification` only** and drops the rest: `momapy_bel` reads `gmod()`
  now, so `GeneAbundance`, `RNAAbundance` and `MicroRNAAbundance` have a `modifications` field
  too, and one of those reaching `_read_protein_modification` fails on `amino_acid`.

### `bel_projection.py` — the influence-graph projection of a BEL model
**Only `bel2cd` and `2_05` import it**, it imports nothing from the package, and it has **two
public functions**:

    projection = get_bel_influence_graph_projection(bel_model)  # -> BELModel
    label = get_element_bel_string(element)

It takes a model and returns a model: no session, no collection names. The projection is a
`momapy_bel.core.BELModel` whose **relations are its statements** (`Increases`,
`DirectlyIncreases`, `TranslatedTo`, `Decreases`, `DirectlyDecreases`, `Regulates`), plus every
projected element carrying none of them, as a bare abundance, because `2_05` counts isolated
nodes.
- The projected elements are the **molecular entities and phenotypes** of the model
  (`PROJECTED_ELEMENT_CLASSES`, matched on the *exact* class since every abundance derives from
  `Abundance`; `PopulationAbundance` and the process terms are out) **minus the structural
  constituents** — an element another entity holds as a sub-term (a complex member, an activity's
  subject) and that carries no causal relation of its own, whose wiring lives on its container.
  That is what makes an isolated node mean something: a protein nothing acts on, not a protein
  that only ever appears inside a complex.
- **Elements are value-equal, so two BEL terms describing the same thing are one element.** BEL
  equality is *stricter* than the CellDesigner equality `bel2cd` goes on to apply (which
  additionally merges `act(X)` with `X`, `bp` with `path`, and whitespace variants of one
  `var()`), so nothing that merges here could have stayed apart in the map.
- **Self-loops are kept** — the statistics run on this, and a filter here would silently move a
  node's degree. `bel2cd` drops them instead, where the reason lives.
- `get_element_bel_string` is the `BELWriter` round trip (`_bel_element_to_string`), used as a
  label by `2_05` and as the sort key that keeps the glyph order — and so `renumber_ids`'
  numbering — stable: a complex's members are sorted by the writer, so the string is canonical.
  Its natural home is `momapy_bel`, as a public function next to the reader and the writer; the
  wrapper is here only because `momapy_bel` was not touched.

### `bel_lpg.py` — pylpg node classes for the BEL elements
A copy of `momapy_kb/lpg/celldesigner.py` pointed at `momapy_bel.core`, with the layout base
classes dropped since BEL has no layout. Importing it registers the 44 element classes plus
`BELModel` and `BELMap` — measured, 9 registered classes become 90 — which is what lets
`execute_query_as_objects` rebuild a stored `BELModel`. It exists because there is no BEL module
in `momapy_kb` and there will not be one; that is the whole of it.

### `queries.py` — Neo4j read helpers
- Hydrating objects from query results needs their pylpg node classes registered, which a session
  does not do on its own: it knows only what it saved and what a static type hint names, and a
  type hint only ever names the general class (`CellDesignerModel.species` says `Species`, never
  `GenericProtein`). `import momapy_kb.lpg.celldesigner` registers the CellDesigner ones (202
  classes) and `import commute_dm.bel_lpg` the BEL ones. **Both are imported for their effect
  alone** — `submaps` does the first, `2_05` the second — and neither has a visible use, so
  neither is dead. This replaces the `prewarm_session` that used to live here, written when
  `momapy_kb.lpg.celldesigner` did not import.
- `get_collections_for_nodes`, `get_subunits`, `get_identifiers`, `get_annotations`,
  `get_ids_and_context` — membership / cross-reference lookups. These match `(node:ModelElement)`,
  which BEL nodes now carry too (`2_00`), so they answer for both kinds. `get_collections_for_nodes`
  is the exception: it traverses `(:CellDesignerMap)`, so it still returns `{}` for BEL nodes.
- **`get_annotated_nodes(session, nodes, with_subunits)`** — "which nodes' annotations stand for
  this node". CellDesigner-only: a species is always its own annotated node, which is what makes
  the caller uniform, and reaches its subunits by `HAS_SUBUNIT` under `with_subunits`. The BEL
  `HAS__*` branch it used to carry is gone with the BEL branch of the pipeline — the AD KG is a
  CellDesigner collection by the time anything queries it, and the export reads its annotations by
  `(class, namespace, identifier)` instead.
- **`get_top_level_species_for_nodes(session, nodes, collection_name)`** — "which species of this
  collection are this node, or contain it". A species here is a top-level element of the
  collection's model, so a complex counts and a complex nested in another complex does not — its
  outermost complex is the answer. A node the model lists gives itself, every listed complex
  reaching it by `HAS_SUBUNIT` gives that complex, and a node can give several. **The collection
  is part of the query and must be**: under hash integration one node can be a species of one
  collection and a subunit of another, so an unscoped lookup would let a walk start in the wrong
  collection. This is what `get_interface(..., with_subunits=True)` stands its rows on.
- **`get_symbols_for_identifiers(session, identifiers)`** — UniProt accession → gene symbols, keyed
  on the **accession** and never on a node. One element's annotations sit in one bag, so the symbol
  is found next to the accession. `core.get_interface_display_names` is built on it.
- `get_nodes(session, element_ids)` — the bridge from the id-only traversal back to DB nodes, for
  the node-taking `gea` helpers.

### `core.py` — analysis orchestration
Ties queries + submaps + gea together for the interface workflows. **It knows nothing about BEL**:
every collection it sees is an ordinary CellDesigner one, the AD KG included.
**No pairing is hard-coded in the library.** `get_interface`, `load_submap_inputs` and the three
interface entry points all take the collection names as **required** arguments, and both `4_10`
and `4_20` define them per pairing in a `PAIRINGS` list in their parameter cell — COVID→PD and
COVID→AD.
Every name in an interface tuple must be **distinct**:
`get_interface` compares `size(collect(DISTINCT collection.name))` against
`size($collection_names)`, so a repeated name returns `{}` with no error.
`get_interface` takes **`with_subunits`** (default `True`), and it decides **which species a
protein stands on in each collection**, not whether the protein is in the interface at all.
`False` keeps only the proteins a collection lists among its species, each standing on itself, and
is the old behaviour. `True` stands a protein on **every species of the collection that is it, or
that contains it** — so a protein held only inside a complex stands on that complex, and a protein
that is both free and bound stands on both. **Nothing passes `False` any more**; it is kept
because it is the plain reading of the question "which proteins does this collection list".
The expansion is a second query, `queries.get_top_level_species_for_nodes`, not a condition on the
first, and it is **per collection**: under hash integration one node can be a species of one
collection and a subunit of another, so looking for the holding complex anywhere would let a walk
start in the wrong collection. It is also **per row**, not per protein, which is what reaches a
protein that is a species in one collection and a member in the other.
A row gains a fourth key, `subunit`: the protein a complex stands for, `None` when the node is the
protein itself. **Nothing branches on it** — only `3_00` reads it, to archive why a complex
represents a protein.
Note `4_20`'s own `WITH_SUBUNITS` is a different flag — `gea`'s widening of a selection to its
species' subunits — and stays `True`. The two compose: a seed that is a complex contributes its
members' genes.
`get_interface` joins collections on shared UniProt annotations, and is **CellDesigner-only** —
it keys on `(:Item)-[:HAS_KEY]->(:Protein)`, which a BEL `ProteinAbundance` does not carry, and a
stored BEL collection holds no cross-reference anyway (its annotations are BEL's own per-statement
ones). That is not a loss: the AD KG's UniProt ids reach the join on the species `bel2cd` writes
them onto. It used to cover both, back when the BEL side came from the cypher dump and its protein
nodes were labelled `:Protein`; a `[..., "AD_KG_BEL"]` triple now simply returns `{}`, which is
why `3_00` no longer asks for one.
`load_submap_inputs` returns the whole per-run setup — `(influences, source_map,
node_id_to_object)` — in one call, hydrating everything through the one `node_id_to_object` cache,
which is what makes a shared database node exactly one Python object. The gene-set analyses of
`4_20` do not call it: they walk node ids and read annotations off the DB, never hydrating a
momapy object, so they call `submaps.load_signed_influences` directly (a fraction of a second,
against the source map's half a minute and its 0.7 GB).
`make_and_write_submaps_from_interface` assembles and writes the
sub-map upstream of the upstream seeds and downstream of the downstream seeds around each
interface protein, with a synthetic central node and a fill per walk direction; then
`make_goat_analysis_from_interface`, `make_intersection_analysis_from_interface`,
`make_goat_gene_lists`.
`get_interface_display_names` names each map by the **accession**, through
`queries.get_symbols_for_identifiers`: one element's annotations sit in one bag, so the symbol is
found next to the accession. It deliberately does **not** read the symbol off the rows' nodes any
more — a complex standing for a protein carries no symbol, and the name must not depend on which
node represents the protein. It still suffixes the accession when two accessions share a symbol
(`BBC3`), and still falls back to the accession when an accession gives no single symbol
(`O14920` is annotated alongside both `IKBKB` and `CHUK`; `P09429` and `P37840` alongside four
symbols each). The old fallback to a BEL node's `name` is gone — it resolved nothing.
`_split_interface_seeds` returns `"upstream"` / `"downstream"` keyed **node ids**, so pointing
the analysis at another pair of collections is a parameter change. It **filters nothing**, and
needs to filter nothing: `get_interface` has already made every row's node a species of the
collection that row names. It used to take `node_id_to_object` and `source_map` and drop seeds
that are not species of the merged source map, and that test was **collection-blind** — it asked
whether a node is a species of *either* collection, not of the one its row names — so it was
removed rather than left doing a right test's job wrongly.
`load_submap_inputs` still returns both — the map assembly needs them. There is no seed widening any more — it existed because a BEL KG kept the
UniProt annotation on `p(HGNC:X)` and the wiring on `act(p(HGNC:X))`, and the export makes those
one species. `max_level` means hops, plainly. All three interface entry points take
`influences` and the two collection names, but only the map one takes the `source_map`.
`import commute_dm.gea` is **lazy** (`_gea()`) because `gea` imports `rpy2` at module
scope, which fails wherever R is unusable (e.g. `.venv-sandbox`).

### `gea.py` — gene-set enrichment
Builds gene sets from node annotations, writes GMT, and runs the R `goat` package via `rpy2`
(`make_goat_analysis`). **The namespace is not a detail**: GOAT gene sets are `ncbigene`, because
`goat.test_genesets` joins them against the gene lists' `gene` column, which `make_goat_gene_lists`
writes from the HGNC dataset's `entrez_id`; the intersection analysis uses `hgnc.symbol`, because
it matches the `symbol` column. `uniprot` is what the *interface* joins on and is used by neither.
`make_gene_set_from_nodes` goes through `queries.get_annotated_nodes` rather than `get_subunits`,
which is what lets a selection be widened to a species' subunits in one place.

### `utils.py` — IO/cypher helpers
Directory (re)creation, CellDesigner XML annotation injection (`add_annotation_to_file`),
`subgraph_to_cypherl` (export a query result as a CREATE/MATCH cypher dump), `merge_relationship`.

## Conventions & gotchas

- Cypher queries are mostly **f-string interpolated** (including element-id lists). Keep that
  style for consistency, but never interpolate untrusted input.
- Collection names (`"COVID_DM_CD_AF"`, `"PD_DM_CD_AF"`, `"AD_KG_CD_AF"`, …) are referenced by
  string, but **only from the notebooks** — each defines the pairing it runs on in its parameter
  cell and passes it in. The library holds no collection name at all: `2_00` names the
  collections it saves.
- `queries.get_annotations` / `get_subunits` / `get_ids_and_context` all `MATCH
  (node:ModelElement)`. Note this only gets a node's *own* annotations: an `Activity` or `Complex`
  has none. Nothing in the analysis queries a BEL node: the AD KG is a CellDesigner collection by
  then, and the transformation writes its annotations itself.
- **The knowledge graphs' `.bel` files are the source now**, and the cypher dumps stay in `data/`
  unused. `data/covid_kg/bel/covid_kg.bel` is not curated: `scripts/rebuild_covid_kg_bel.py`
  writes it out of the stored `COVID_KG_BEL` collection, which is the version every analysis has
  been run against and which the curated per-paper files in that directory only ~95% cover. Run it
  from the repository root, **while the collection is still loaded** — `2_00` wipes the database.
  It is deterministic: two runs give the same bytes.
  Its one subtlety is what counts as a top-level term. A term is a sub-term when another term
  points **at** it through a double-underscore `HAS__*` edge; it is not one because it *has*
  sub-terms, which every complex does. Getting that backwards is what made an earlier revision
  emit 6 bare terms instead of 786, silently dropping 730 standalone `complex(...)` assertions and
  750 nodes from `2_05`'s COVID row while leaving the edge count untouched. Six of the 59 CBM files did not parse and were fixed in place: multi-line `SET Support`
  values need a backslash at the end of every line of the span but the last (the files are CRLF,
  so the backslash goes *before* the carriage return), `20234358.bel` had two `association` glued
  to the following term, and `35063125.bel` wrote three hyphenated HGNC symbols unquoted. All 59
  parse afterwards, so `2_00` reads the whole directory and needs no list of exceptions.
- The "201 shared species … none is a modulation endpoint" note in `submaps.py` / `core.py` needs
  the qualifier "**in either AF collection**": 3 of the 201 are endpoints of `PD_DM_CD`
  (process-description) modulations. The two AF collections also share 289 `SpeciesTemplate` nodes
  and 1 `Compartment`, not only the 201 species. Nor does "not a modulation endpoint" mean a
  subunit cannot be *selected*: it can enter a selection as a boolean-gate input or as an
  interface seed, which is how `submaps.find_subunit_top_level_aliases` fires.
- **`pd2af.utils._build_dot_graph` drops the dot cluster of any compartment with `outside=None`**
  (`utils.py:339-346`, no `else`), so a drawn top-level compartment silently gets no box and its
  species keep their throwaway positions. Never exercised before, because the only stored root
  (`default`) is undrawn and gets no cluster at all. Worked around by the undrawn BEL root
  compartment (`bel2cd` invariant (d)), so no `pd2af` change is needed — but it is a real bug
  worth fixing upstream.
- `5_00` and `6_00`/`6_10` were **not** updated for the changeover: `5_00` will hit a new
  `StopIteration` on the exported `Unknown` species and `6_00`'s intersection becomes 5-way. Both
  were already stale or non-running.
- `6_10_make_immuno_targets_graphs.ipynb` is **stale**: it calls a
  `make_and_write_cd_maps_from_interface` that no longer exists, with seeds from the non-AF
  collections. It was already non-functional before the `submaps.py` rewrite (it rendered via the
  IG edges that no longer exist).
- Map versions in use are recorded in `data/maps/versions.txt`.
- `poubelle/` is a scratch/trash dir; ignore it.
