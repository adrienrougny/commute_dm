# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`commute_dm` is a research analysis pipeline (the "COMMUTE" disease-map project). It takes
CellDesigner disease maps for COVID-19 and Parkinson's disease (PD), plus BEL knowledge
graphs (AD/PD/COVID), loads them into a **Neo4j** graph DB via the `momapy_kb` stack, finds
the **interface** of proteins shared across the collections, then runs **gene-set enrichment
(GOAT)** and intersection analyses against bulk RNA-seq DE results, and renders new
CellDesigner maps of the interface subgraphs.

**The AD knowledge graph is an ordinary CellDesigner collection.** `2_01` exports its
influence-graph projection to one CellDesigner file and imports that file as `AD_KG_CD_AF`, so
every analysis sees three activity-flow CellDesigner collections — `COVID_DM_CD_AF`,
`PD_DM_CD_AF`, `AD_KG_CD_AF` — and no BEL. Only `2_01` (through `bel_export.py` and
`bel2cd.py`) knows the AD side was ever BEL.

The Python package (`src/commute_dm`) is the library; `scripts/main_analysis/*.ipynb` are the
numbered, ordered driver notebooks. Outputs land under `data/`.

## Environment & tooling

- Managed with **`uv`** (`uv.lock`, `pyproject.toml`); Python ≥3.10.
- Five core deps are **local editable installs** at `/home/rougny/code/{momapy,
  momapy_kb,fieldz_kb,pylpg,momapy_bel}` (see `[tool.uv.sources]`). These are the heart of the
  system — `momapy` (disease-map data model + CellDesigner read/write + SKIA rendering),
  `momapy_kb` (Neo4j-backed LPG sessions over momapy objects), `fieldz_kb`/`pylpg` (the LPG
  node-class machinery), `momapy_bel` (the BEL data model the KG is loaded into, and its
  `BELWriter`). `momapy_bel` is a **runtime** dependency now, not a dev one, and it must be the
  version whose abundances take a tuple of `variants` (two AD proteins carry two `var()`s). These paths have gone missing in past checkouts (the full
  `momapy_kb.lpg.backends.neo4j` import would then fail) — if that recurs, query Neo4j directly
  with the `neo4j` driver as a fallback (DB is `bolt://localhost:7687`).
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
  back to ABI mode, and that dies with `symbol 'R_getVar' not found in library 'libR.so'`. The
  user's machine has R 4.6.1 and is fine; **Claude's sandbox is a Debian 12 container whose R is
  4.2.2**, so GOAT cannot be run from there — the R version is the thing to check, not the venv.
  Everything upstream of the R call (gene sets, GMT, the summary assembly) is testable in the
  sandbox by shimming `sys.modules` for `rpy2*`, and the intersection analysis needs no R at all
  beyond `gea`'s module-scope import.
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
- `2_00_save_collections` — wipe the DB, import BEL KG cypher dumps (`cypher-shell`), then load
  the CellDesigner XML maps into named **Collections** (`COVID_DM_CD`, `PD_DM_CD`, plus
  `*_KG_BEL`). Uses `with_membership_edges=True` (the `HAS_MODEL_ELEMENT` edge that
  `get_collections_for_nodes` relies on) and `integration_mode="hash"`.
  Two things it does to the BEL side make the whole BEL KG legible to the generic helpers:
  every imported node is labelled **`ModelElement` as well as `BELModelElement`** (the helpers in
  `queries.py` match `(node:ModelElement)`, and without it they returned nothing for a BEL node,
  *silently*); and `add_hgnc_annotations_to_bel_kg` writes **four** namespaces on `AD_KG_BEL`'s
  HGNC-encoded proteins — `uniprot` (what the interface joins on), `ncbigene` (what GOAT matches),
  `hgnc.symbol` (what the intersection matches) and `hgnc` — all looked up in the HGNC dataset by
  the symbol the node is already keyed by (`p(HGNC:"MAPT")` → `name='MAPT'`). Coverage is the same
  for all four (1332 of 1343 HGNC-namespace AD proteins), so **the interface is unchanged** by
  the added namespaces. The encoding matches the stored maps exactly: one `Item`/`Bag` per
  protein, **one single-resource `RDFAnnotation` per identifier**, one shared `BQBiol IS`
  qualifier — a stored `BCL2` likewise carries `hgnc.symbol`, `ncbigene`, `uniprot`, `hgnc` … as
  separate annotations in one bag. The function is idempotent (it deletes the collection's
  existing annotation items first), and that is why its qualifier `MERGE` ends in
  `ORDER BY … LIMIT 1`: **`MERGE` binds every existing match, one row each**, and loading the
  CellDesigner maps leaves a *second* `BQBiol IS` node, so on a re-run against an already-loaded
  DB the MERGE yields two rows and every `CREATE` under it silently runs twice. A clean run never
  hits it (BEL is annotated before the maps are loaded); a re-run always would.
- `2_01_make_ad_kg_cd_af_collection` — **exports the AD BEL KG as an ordinary CellDesigner
  collection**, `AD_KG_CD_AF`, and is the only importer of `bel_export.py`. Ten cells, doing
  nothing but the work: build the map (~2 s), read the annotations off the BEL nodes, renumber,
  write one ~16 MB file, seed the integration cache, import. **No counts, no prints, no guard, no
  read-back check, no post-import assertions** — those were verified once and are recorded here
  instead. It runs *before* `2_05`, which is why it is `2_01`.
  Three things it decides: **`act(X)` and `X` are one species** (the projection has no edge
  between them, so keeping them apart severs 6889 upstream→downstream 2-hop paths that 3 hops
  cannot recover; the price is that X carries an active border whenever BEL asserts any activity
  of it); **species carrying no signed modulation are dropped** (1228 of 3979 projected nodes —
  nothing walks or draws them); and **some BEL terms merge**, because CellDesigner does not
  distinguish `bp` from `path`, `m` from `r`, two `ma()` codes of one activity, or whitespace
  variants of one `var()` (107 groups; none is an interface protein).
  Result: 2418 species, 1206 templates, 4631 modulations, 5 compartments, 779 subunits, 778
  UniProt ids on 954 annotated species. The written file round-trips: read back, it has the same
  counts and the same 778 UniProt ids.
  **The notebook no longer runs `check_identity_invariants` or the round-trip read-back.** Both
  still exist (`submaps.check_identity_invariants`) and are the first thing to reach for if the
  file ever fails to read back with a `KeyError` — that is the failure mode they were written
  for, and the invariants themselves are documented under `bel2cd.py` below.
  **The save spans calls, and that is the one thing to understand here.**
  `integration_mode="hash"` only ever deduplicated *within* one save call — `save_from_objects`
  built its `object_to_node` map fresh and never read the DB — so a second call would have
  duplicated every element shared with the stored maps, and duplicates break the writer
  (value-equal templates collapse in a frozenset while `<proteinReference>` resolves by object
  identity: writes clean, fails to read back). `momapy_kb`/`fieldz_kb` now expose the cache as
  `object_key_to_node`; the notebook seeds it from the DB with one query (4 s, 31852 entries) and
  passes it to `save_collections_from_file_paths`. **The seed set must be closed under descent** —
  a cache hit skips the walk over that object's descendants, so an unseeded child would get no
  `HAS_MODEL_ELEMENT` edge — and it is, being `model.descendants()`. Measured after import: 0
  elements without a membership edge; 258 templates / 94 top-level species / 22 subunits fused
  with the stored collections, 0 compartments and 0 modulations; and **no fused top-level species
  is a stored top-level species** (a stored one carries a compartment, a BEL one does not), so
  each is top-level in AD and a subunit in a COVID/PD complex — the case
  `submaps.make_submap_from_model_elements`' subunit promotion already covers.
  **The save is a plain `CREATE`, and nothing guards it any more**: running the last cell twice
  gives two copies of the collection. Drop it first if re-importing —
  `MATCH (c:Collection {name: 'AD_KG_CD_AF'})-[:HAS_ENTRY]->(e)-[:HAS_OBJ]->(m) DETACH DELETE c, e, m`
  (that leaves the model elements behind, which the next import re-integrates).
- `2_05_get_collections_statistics` — descriptive statistics of every collection. It builds its
  BEL influence graphs by calling `bel_export.get_bel_influence_graph_projection` per collection
  and wrapping the returned `BELModel` in a `networkx.MultiDiGraph` — **one node per momapy_bel
  element, one edge per relation statement** — labelled with `bel_export.get_bel_string(element)`
  and `type(element).__name__` (which is what the phenotype test compares against), and carrying
  `element.id_`, the id of the database node the element was loaded from, which is all the
  isolated-reason cell needs to go back to the DB. Node/edge counts: 3979/4910, **1832**/1760,
  2255/1690, 2074/2145. Only the PD node count moved, from 1835: momapy_bel elements are
  value-equal and its three `gmod(TestNS:"TestName")` genes merge with their plain genes, `gmod`
  having no place in a `GeneAbundance`. They are a test artefact and AD has none.
- `3_00_get_interfaces` — query proteins present across collections (joined on **UniProt RDF
  annotations**), write `data/.../interface/*.json`. The COVID×AD interface is **159** two-way and
  **111** three-way, smaller than the 189/127 of the old BEL pairing because `2_01` drops the
  isolated nodes on purpose: those identifiers could never seed a walk, so after the drop the
  interface *is* the seedable set.
- `4_00_make_goat_gene_lists` — turn raw RNA-seq DE CSVs (`data/rnaseq/`) into GOAT gene-list
  CSVs, mapping Ensembl→Entrez/symbol via the HGNC dataset.
- `4_10_make_interface_graphs` — **both pairings**, in a `PAIRINGS` list of the same shape as
  `4_20`'s: per pairing the collection names, the interface tuple, `max_levels` and the output
  directory. COVID→PD keeps the three-way interface filter (now with `AD_KG_CD_AF`); COVID→AD is
  the two-way `COVID_DM_CD_AF × AD_KG_CD_AF`. Both sides of both pairings are stored CellDesigner
  collections, so `core.load_submap_inputs` has one path and there is no BEL case anywhere.
  `max_levels` stays per pairing, and that is a property of the downstream graph, not of the code:
  `[1,2,3]` for AD against `[2,3,4,5,6]` for PD, because the AD influence graph is much denser.
  Current COVID→AD output: **113 maps over 51 proteins**, from an interface of 159 — more than the
  old BEL pairing's 107 over 47 despite the smaller interface, because merging `act(X)` with `X`
  makes the downstream side denser (median downstream selection at 3 hops goes 7 → 19). The
  binding constraint is `MIN_N_NODES=5` on the **COVID upstream** side. The last cell reads every
  written file back, which is the assertion the identity invariants exist for.
- `4_20_make_interface_goat_analysis` — GOAT enrichment + intersection analyses of interface
  subgraphs vs. the gene lists, for **both pairings and all three modes** (`upstream`,
  `downstream`, `upstream_and_downstream`). Influences come straight from
  `submaps.load_signed_influences` — seconds, since no map is hydrated.
  *Not runnable from Claude's sandbox*, whose Debian 12 R is 4.2.2 against a pinned rpy2 3.6.7
  that needs R ≥ 4.5 (see "Environment & tooling"). Nothing to do with the venv — the user's R is
  4.6.1 and the cell runs there.
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
  collections, hydrated (**~3 min, 0.2 GB**) and merged by `merge_maps`. Only the map path pays
  this; `4_20` never calls it.
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

### `bel2cd.py` — a momapy_bel element as real activity-flow CellDesigner content
`bel_export` imports `bel2cd`, never the reverse. It knows nothing about Neo4j and nothing about
BEL strings: everything is read off the elements' fields — `modifications`, `variants`,
`fragment`, `location`, `members`, `abundance`. Three steps:
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
  at the top level, on the finished object, in `bel_export`.
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
  `nucleus` (`outside=AD_KG_BEL`, the root's name) value-distinct from a stored one
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
  the 14 `HAS__GMOD` edges are a PD-only test artefact.)

### `bel_export.py` — the AD BEL KG as an exportable CellDesigner map
**Only `2_01` and `2_05` import it**, and it has **three public functions**:

    projection = get_bel_influence_graph_projection(session, collection_names)  # -> BELModel
    cd_map, element_to_annotations = make_cd_map_from_bel_influence_graph_projection(
        session, collection_names, projection
    )
    label = get_bel_string(element)

The projection is a `momapy_bel.core.BELModel`: **the relations are its statements**
(`Increases`, `DirectlyIncreases`, `TranslatedTo`, `Decreases`, `DirectlyDecreases`,
`Regulates`), and every projected element carrying none of them rides along as a bare
abundance, because `2_05` counts isolated nodes. There is no node-id table and no separate term
graph — the sub-terms are *nested in the elements* — but each element keeps the id of the node it
was loaded from in `id_`, which is excluded from equality and is what lets `2_05` ask the DB about
an element. `get_bel_string` is the `BELWriter` round trip, used as a label and as the sort key
that keeps the glyph order (and so `renumber_ids`' numbering) stable: a complex's members are
sorted by the writer, so the string is canonical.
- **Elements are value-equal, so equal nodes become one element.** Measured over the four KGs,
  the only merge is PD's three `gmod(TestNS:"TestName")` genes (§ `2_05` above). BEL equality is
  *stricter* than the CellDesigner equality the export goes on to apply (which additionally merges
  `act(X)` with `X`, `bp` with `path`, and whitespace variants of one `var()`), so nothing that
  merges here could have stayed apart in the map.
- **The loading raises rather than skips** — on a node whose labels match no known entity class,
  on an unclassified `HAS__*` edge type, on an activity without exactly one subject. Silently
  skipping is how the old influence-graph code rotted into a smaller graph than it claimed.
- **`gmod()` is dropped**: `momapy_bel`'s `GeneAbundance` has no place for it, and the only ones
  in the data are PD's three test artefacts. `var()` is *not* — `momapy_bel` takes a tuple of
  `variants`, which two AD proteins need (`p(MGI:"App",var("D,23,N"),var("E,22,Q"))`).
- Species are **interned by value** in the CD conversion: two `ma()` codes of one activity describe
  identically, as do whitespace-variant proteoforms — and so do `act(p(X))` and `p(X)`, which is
  the point (§ `2_01` above). AD: 3979 projected elements → 2418 species after the isolated drop.
- **Every species built is recorded, subunits included**, which is what the annotations need: a
  subunit is a per-occurrence object, and **16 of the interface's UniProt identifiers are carried
  only by complex members** (without them the two-way interface is 173 rather than 189).
- `_make_element_to_annotations` → `{species: frozenset[RDFAnnotation]}`, called at the end of
  `make_cd_map_from_bel_influence_graph_projection`. **One query**, keyed on
  `(element class, namespace, identifier)` — what a momapy_bel element knows about itself, and
  what `2_00` keyed the annotations on; the class is part of the key because a `g(HGNC:"APP")` must
  not pick up the protein's cross-references. An activity's subject is followed, a complex's
  members are not, since each is drawn as its own subunit species and annotated there. Accumulated
  **by value**, so a top-level species and a value-equal subunit share one entry, which is what the
  writer looks up. Result: 954 annotated species, 778 UniProt ids.
- `drop_isolated_species` runs **after** interning and after the modulations are built, so a
  protein isolated in its own right but wired through its activity form survives. A dropped species
  that belongs to a kept complex still appears as that complex's subunit.
- `REGULATES` **stays in the projection and is dropped from the map**: it is unsigned, and
  `submaps.SIGNED_MODULATION_CLASSES` excludes `Modulation`, so mapping it would give an edge that
  is *walked but never drawn* — but it is 2.6% of the AD KG's causal edges, which `2_05` counts,
  so it is `BEL_RELATION_CLASS_TO_MODULATION_CLASS` (the map's table) that leaves it out, not
  `BEL_RELATION_TYPE_TO_RELATION_CLASS` (the projection's).
- `Abundance` is split on the **`namespace` field** — `a(CHEBI:…)` → `SimpleMolecule` (350 in AD),
  anything else → `Unknown` (77): MESH/CONSO/GO abundances are proteins, aggregates and cellular
  components, not molecules. The rest of `bel2cd.BEL_CLASS_TO_CD_CLASSES` selects the template
  class too, so an entry decides whether a species declares a `<protein>` / `<gene>` / `<rna>`.
- `PROJECTABLE_NODE_CLASSES`, `BEL_NODE_CLASS_TO_ELEMENT_CLASS` and the two relation tables are
  the only vocabulary this module owns, being export decisions; the node classes and relation
  types come from `bel_vocabulary.py`.
- **Self-loops are kept in the projection** — it is what `2_05`'s statistics run on, and a filter
  there would silently move a node's degree. The CD conversion drops them instead, where the reason
  lives: pd2af's arc geometry needs two distinct endpoints, and interning creates species-level
  self-loops (an `act(X) → X` relation) that no element-level filter could see anyway.
- Accepted consequences: a fragmented gene yields two `<protein>` entries (GENERIC + TRUNCATED),
  which is idiomatic CellDesigner; a complex listing the same member twice loses the multiplicity
  to `frozenset[Species]` (`homomultimer` would be the fix, out of scope).
- No global auto-layout: `make_auto_layout` (graphviz `dot`) does not finish on a graph this size.
  Glyphs get local, non-overlapping positions; the sub-maps of `4_10` lay out their own selections.

### `bel_vocabulary.py` — the BEL vocabulary
Constants only, importing nothing from the package: the node classes
(`MOLECULAR_ENTITY_NODE_TYPES`, `PHENOTYPE_NODE_TYPES`, `REACTION_NODE_TYPES`, …), the relation
types by family (`INFLUENCE_RELATIONSHIP_TYPES` and the four others, plus the two unions) and the
`HAS__*` structural edge types (`MEMBER_EDGE_TYPES`, `MODIFIER_EDGE_TYPES`,
`OTHER_STRUCTURAL_EDGE_TYPES`). **Every reader of a BEL graph takes its vocabulary from here** —
`bel_export` and `2_05`, each binding the names it interpolates. `queries` is no longer one of
them: no BEL node is queried by the analysis any more, so its BEL branch is gone.

### `queries.py` — Neo4j read helpers
- `prewarm_session` **must be called** before hydrating CellDesigner objects from query results —
  it registers pylpg node classes for momapy subclasses the static type walk misses.
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
`get_interface` joins collections on shared UniProt annotations;
`load_submap_inputs` returns the whole per-run setup — `(influences, source_map,
node_id_to_object)` — in one call, hydrating everything through the one `node_id_to_object` cache,
which is what makes a shared database node exactly one Python object. The gene-set analyses of
`4_20` do not call it: they walk node ids and read annotations off the DB, never hydrating a
momapy object, so they call `submaps.load_signed_influences` directly (seconds, against the source
map's couple of minutes and couple of hundred megabytes).
`make_and_write_submaps_from_interface` assembles and writes the
sub-map upstream of the upstream seeds and downstream of the downstream seeds around each
interface protein, with a synthetic central node and a fill per walk direction; then
`make_goat_analysis_from_interface`, `make_intersection_analysis_from_interface`,
`make_goat_gene_lists`.
`_split_interface_seeds` returns `"upstream"` / `"downstream"` keyed **node ids**, so pointing
the analysis at another pair of collections is a parameter change; passing
`node_id_to_object` and `source_map` also drops seeds that exist only as complex subunits, which
cannot seed a walk. There is no seed widening any more — it existed because a BEL KG kept the
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
  cell and passes it in. The library holds no collection name at all: `2_01` names `AD_KG_BEL` and
  `AD_KG_CD_AF` in its parameter cell, `2_00` names the collections it imports.
- `queries.get_annotations` / `get_subunits` / `get_ids_and_context` all `MATCH
  (node:ModelElement)`. BEL nodes used to carry `:BELModelElement` and **not** `:ModelElement`, so
  these returned nothing for them, silently — that is fixed at the source, in `2_00`, which now
  sets both labels. If a BEL query comes back empty for no reason, check the label is there before
  looking anywhere else. Note this only gets a node's *own* annotations: an `Activity` or `Complex`
  has none. Nothing in the analysis queries a BEL node any more: the AD KG is a CellDesigner
  collection by then, and the export reads the BEL annotations itself.
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
