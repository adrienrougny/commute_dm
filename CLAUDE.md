# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`commute_dm` is a research analysis pipeline (the "COMMUTE" disease-map project). It takes
CellDesigner disease maps for COVID-19 and Parkinson's disease (PD), plus BEL knowledge
graphs (AD/PD/COVID), loads them into a **Neo4j** graph DB via the `momapy_kb` stack, finds
the **interface** of proteins shared across the collections, then runs **gene-set enrichment
(GOAT)** and intersection analyses against bulk RNA-seq DE results, and renders new
CellDesigner maps of the interface subgraphs.

**The AD knowledge graph is an ordinary CellDesigner collection.** `2_10` exports its
influence-graph projection to one CellDesigner file and imports that file as `AD_KG_CD_AF`, so
every analysis sees three activity-flow CellDesigner collections — `COVID_DM_CD_AF`,
`PD_DM_CD_AF`, `AD_KG_CD_AF` — and no BEL. Only `2_10` (through `bel_export.py` and
`bel_terms.py`) knows the AD side was ever BEL.

The Python package (`src/commute_dm`) is the library; `scripts/main_analysis/*.ipynb` are the
numbered, ordered driver notebooks. Outputs land under `data/`.

## Environment & tooling

- Managed with **`uv`** (`uv.lock`, `pyproject.toml`); Python ≥3.10.
- Four core deps are **local editable installs** at `/home/rougny/code/{momapy,
  momapy_kb,fieldz_kb,pylpg}` (see `[tool.uv.sources]`). These are the heart of the system —
  `momapy` (disease-map data model + CellDesigner read/write + SKIA rendering), `momapy_kb`
  (Neo4j-backed LPG sessions over momapy objects), `fieldz_kb`/`pylpg` (the LPG node-class
  machinery). These paths have gone missing in past checkouts (the full
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
- `2_10_make_ad_kg_cd_af_collection` — **exports the AD BEL KG as an ordinary CellDesigner
  collection**, `AD_KG_CD_AF`, and is the only importer of `bel_export.py`. It builds the map
  (~2 s), reads the annotations off the BEL nodes, renumbers, writes one ~16 MB file, reads it
  back and checks the four identity invariants, then imports it.
  Three things it decides, all measured in the notebook: **`act(X)` and `X` are one species**
  (the projection has no edge between them, so keeping them apart severs 6889 upstream→downstream
  2-hop paths that 3 hops cannot recover; the price is that X carries an active border whenever
  BEL asserts any activity of it); **species carrying no signed modulation are dropped** (1228 of
  3979 projected nodes — nothing walks or draws them); and **some BEL terms merge**, because
  CellDesigner does not distinguish `bp` from `path`, `m` from `r`, two `ma()` codes of one
  activity, or whitespace variants of one `var()` (107 groups; none is an interface protein).
  Result: 2418 species, 1206 templates, 4631 modulations, 5 compartments, 779 subunits, 778
  UniProt ids on 954 annotated species.
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
  It runs *after* `2_05` in numeric order, so re-run `2_05` afterwards to pick the collection up.
- `3_00_get_interfaces` — query proteins present across collections (joined on **UniProt RDF
  annotations**), write `data/.../interface/*.json`. The COVID×AD interface is **159** two-way and
  **111** three-way, smaller than the 189/127 of the old BEL pairing because `2_10` drops the
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
  *establish* them, in `bel_terms.py`. Verified against the untouched COVID→PD pipeline first,
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
  `bel_terms`: `core` needs the first for every sub-map's synthetic central node, in both
  pairings, and `bel_terms` is export-only. `bel_terms` imports them back; acyclic.
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

### `bel_terms.py` — a BEL term as real activity-flow CellDesigner content
`bel_export` imports `bel_terms`, never the reverse. **The BEL term graph is fully relational in
Neo4j** — every sub-term (`pmod`, `var`, `frag`, `loc`, each complex member, each activity
subject) is its own node reached by a typed `HAS__*` edge, and entity nodes carry `name` /
`namespace` — so **nothing here parses a BEL string**. Four steps:
- `load_bel_terms(session, names)` → `{node_id: BelTerm}` (two queries, <1 s, 5024 AD nodes,
  3178 sub-term edges). `BelTerm` is `eq=False`: the graph is recursive and identity is what the
  memo and the build context are keyed on. `BEL_TERM_NODE_CLASSES` resolves the concrete class
  (pylpg labels a node with its class *and* every ancestor, so `labels(node)[0]` is
  order-dependent); **exactly one** must match or it raises. The **double** underscore in
  `HAS__*` excludes the sparse single-underscore statement route (`HAS_MEMBERS`,
  `HAS_COMPONENTS`) and the `HAS_NODE` / `HAS_OBJ` plumbing. Edge types split into
  `MEMBER_EDGE_TYPES` (complex/composite members **and** the activity subject),
  `MODIFIER_EDGE_TYPES` and `OTHER_STRUCTURAL_EDGE_TYPES` (BEL *process* terms — translocations,
  reactions — which are never projected, enumerated only so loading raises on an unseen type).
- `describe_term(term, cache)` → `TermDescription` — pure, DB-free, unit-testable, memoised on
  `id(term)`; raises on a cycle, an unmapped class and an unroutable modifier. `p(HGNC:"GDE1")` →
  `GenericProtein("GDE1")`; `pmod(Ph,Y,18)` and `var("P, T, 212")` → `Modification` badges
  (whitespace normalised, so `"P, S, 9" == "P,S,9"`); `var("K,670,N")` → `StructuralState("K670N")`;
  ontology `pmod(GO:…)` → `StructuralState`; `frag(...)` → `TruncatedProtein` +
  `TruncatedProteinTemplate` + a `StructuralState`; `act(...)` → the **subject's own class** with
  `active=True` and the `ma()` code dropped; `composite(...)` → `Complex` +
  `StructuralState("composite")`; an unnamed complex is named by joining its members' names.
  Modifiers are routed by the target class's *capabilities*; a class with no container for one
  (a `Gene`, an `Unknown`) folds it into the **species** name as `NAME [x|y]` — never into the
  **template** name, or one protein would declare N `<protein>` entries.
- **Nothing is ever invented for a sub-term that says nothing.** `_describe_modifier` returns a
  `None` payload and the caller records no state: a `frag("?")` (range unknown, no descriptor)
  still switches the class to `TruncatedProtein`, which is what says "a fragment of", but adds no
  `StructuralState("?")` on top — the placeholder was noise, and it turned up as
  `<structuralState structuralState="?">` in the output. `var("?")` and `var("p.?")` are the
  opposite case and are **kept verbatim**: there `"?"` is the value BEL states, and it is the only
  thing telling that proteoform apart from the plain protein. Same principle as the nameless
  modification residue above.
- `make_build_context(terms, active_term_ids)` — three passes, in this order, because a species
  is frozen and its template cannot be patched afterwards: describe every term (roots **and**
  recursive members) → accumulate `{(template class, name): {residue names}}` and build one
  template per key → **re-derive** the residue lookup from that object. Nothing is interned
  against a source map any more: the export is standalone, and identity with the stored
  collections' elements is the *save's* job (hash integration with a seeded `object_key_to_node`).
  `active_term_ids` comes from `collect_activity_subject_term_ids` and is what makes `act(X)` and
  `X` describe identically — pass **every loaded term**, not just the projected roots, since a
  protein is drawn active whenever BEL asserts any activity of it.
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
  residue names through `bel_terms._residue_sort_key`, since `None` does not compare with `str`.
- **(c) a subunit object is never the same object as a top-level species** — `p(HGNC:"APP")` is
  both a projected node and a member of several complexes, so memoising species by node id (the
  obvious way to get value-collapse) would make APP's own glyph disappear. Interning happens only
  at the top level, on the finished object, in `bel_export`.
  **This one is an in-memory property of the export and does not survive the save**: a top-level
  `p(X)` and the same `p(X)` inside a complex are value-equal, so hash integration stores them as
  one node and hydration returns one object. A `source_map` built from `AD_KG_CD_AF` therefore has
  651 subunits that *are* their top-level species, and that is expected —
  `make_submap_from_model_elements` promotes such a species to its own object per sub-map, before
  the checks run. **Never run `check_identity_invariants` on a source map**; it is a pre-write
  check for sub-maps.
- **(d) one compartment object per `(collection, location)`, hanging off an undrawn per-collection
  root.** The root is **not decorative**: `pd2af.utils._build_dot_graph` (`utils.py:339-346`) only
  attaches a compartment's dot cluster to the graph when `compartment.outside is not None` — there
  is no `else` — so a drawn compartment with `outside=None` gets an orphan cluster, its species
  never reach graphviz, and they keep their throwaway positions. (A real pd2af bug, worth
  reporting upstream; worked around here, so no pd2af change is needed.) The root also keeps a BEL
  `nucleus` (`outside=AD_KG_BEL`) value-distinct from a stored one (`outside=default`), and
  therefore the species inside them too. 18 of the 84 BEL `loc()` names already exist as stored
  compartment names, so one sub-map can show two same-labelled boxes from different sources —
  that is the price of keeping them distinct.

Layout notes:
- **Badge angles are residue-derived, not species-local.** CellDesigner stores the angle on the
  *template*, and `_writing.find_residue_angle` takes it from the **first** species using that
  template with a modification on that residue, so a species-local index would move every other
  proteoform's badge on read-back. The angle is `2π · residue.order / len(template residues)`, and
  the reader's exact forward transform is reproduced so `_writing.compute_cd_angle` inverts it
  (verified to ~1e-6 rad over all 217 AD badges).
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
**Only `2_10` imports it.** It loads the KG's influence-graph projection (the one `2_05` defines
and documents — this module holds the executable copy of the queries) and
`make_bel_map(nodes, terms, edges, drop_isolated_species=True)` →
`(cd_map, species_by_node_id, node_id_species_pairs, stats)` turns it into a map of real
activity-flow content via `bel_terms`. `2_10` writes that map to a file and imports the file, after
which the KG *is* a CellDesigner collection and nothing else in the package knows it was ever BEL.
- **`species_by_node_id` is many-to-one**, because top-level species are interned by value: two
  `ma()` codes of one activity describe identically, as do whitespace-variant proteoforms — and so
  do `act(p(X))` and `p(X)`, which is the point (§ `2_10` above). Measured on AD: 3979 projected
  node ids → 3593 species before the isolated drop, 2418 after, over 2803 surviving node ids.
- **`node_id_species_pairs` is every species built, subunits included**, which is what
  `make_element_to_annotations` needs: a subunit is a per-occurrence object no `{node_id: species}`
  map can recover, and **16 of the interface's UniProt identifiers are carried only by complex
  members** (without them the two-way interface is 173 rather than 189).
- `make_element_to_annotations(session, node_id_species_pairs)` →
  `{species: frozenset[RDFAnnotation]}`. `queries.get_annotated_nodes(..., with_subunits=False)`
  says whose annotations stand for a node — an activity's subject is followed, a complex's members
  are not, since each is drawn as its own subunit species and annotated there — then one query
  rebuilds the payloads exactly as `2_00` wrote them. Accumulated **by value**, so a top-level
  species and a value-equal subunit share one entry, which is what the writer looks up.
- `drop_isolated_species` runs **after** interning and after the modulations are built, so a
  protein isolated in its own right but wired through its activity form survives. A dropped species
  that belongs to a kept complex still appears as that complex's subunit.
- `REGULATES` is dropped from both the map and the walk: it is unsigned, and
  `submaps.SIGNED_MODULATION_CLASSES` excludes `Modulation`, so mapping it would give an edge
  that is *walked but never drawn*. 2.6% of the AD KG's causal edges.
- `Abundance` is split on the **`namespace` property** (not a regex on the BEL string) —
  `a(CHEBI:…)` → `SimpleMolecule` (350 in AD), anything else → `Unknown` (77): MESH/CONSO/GO
  abundances are proteins, aggregates and cellular components, not molecules. The rest of
  `bel_terms.BEL_CLASS_TO_CD_CLASSES` selects the template class too, so an entry decides whether
  a species declares a `<protein>` / `<gene>` / `<rna>`.
- `PROJECTABLE_NODE_CLASSES` is what `load_bel_projection` validates against — a separate question
  from what a node is *drawn* as (an `Activity` is projected but has no class of its own).
- Accepted consequences: a fragmented gene yields two `<protein>` entries (GENERIC + TRUNCATED),
  which is idiomatic CellDesigner; a complex listing the same member twice loses the multiplicity
  to `frozenset[Species]` (`homomultimer` would be the fix, out of scope).
- No global auto-layout: `make_auto_layout` (graphviz `dot`) does not finish on a graph this size.
  Glyphs get local, non-overlapping positions; the sub-maps of `4_10` lay out their own selections.

### `queries.py` — Neo4j read helpers
- `prewarm_session` **must be called** before hydrating CellDesigner objects from query results —
  it registers pylpg node classes for momapy subclasses the static type walk misses.
- `get_collections_for_nodes`, `get_subunits`, `get_identifiers`, `get_annotations`,
  `get_ids_and_context` — membership / cross-reference lookups. These match `(node:ModelElement)`,
  which BEL nodes now carry too (`2_00`), so they answer for both kinds. `get_collections_for_nodes`
  is the exception: it traverses `(:CellDesignerMap)`, so it still returns `{}` for BEL nodes.
- **`get_annotated_nodes(session, nodes, with_subunits)`** — "which nodes' annotations stand for
  this node", and the one place the CellDesigner/BEL difference lives. A CellDesigner species is
  always its own annotated node and reaches its subunits by `HAS_SUBUNIT` under `with_subunits`.
  A BEL node's cross-references sit on the `:Protein` terms it merely *contains*, reached by
  `bel_terms.MEMBER_EDGE_TYPES` (`HAS__PROTEIN`, `HAS__COMPLEX`, …), and the same distinction is
  made: an **activity's subject is always followed** — `act(p(X))` is X in another form, as an
  active CellDesigner species is still that species — while a **complex's or composite's members
  are followed only under `with_subunits`**, exactly as `HAS_SUBUNIT` is. Verified on AD: every
  `Activity` reaches its subject either way, `Complex`/`Composite` only under `with_subunits`, a
  plain protein needs no traversal. The edge-type list is *repeated* here rather than imported so
  `queries` depends on none of the BEL machinery — `bel_terms` is where the classification is
  decided, and a new member edge type has to be added in both.
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
which is what makes a BEL node contribute the annotations of the entity terms it is built from.

### `utils.py` — IO/cypher helpers
Directory (re)creation, CellDesigner XML annotation injection (`add_annotation_to_file`),
`subgraph_to_cypherl` (export a query result as a CREATE/MATCH cypher dump), `merge_relationship`.

## Conventions & gotchas

- Cypher queries are mostly **f-string interpolated** (including element-id lists). Keep that
  style for consistency, but never interpolate untrusted input.
- Collection names (`"COVID_DM_CD_AF"`, `"PD_DM_CD_AF"`, `"AD_KG_CD_AF"`, …) are referenced by
  string, but **only from the notebooks** — each defines the pairing it runs on in its parameter
  cell and passes it in. The library holds no collection name at all: `2_10` names `AD_KG_BEL` and
  `AD_KG_CD_AF` in its parameter cell, `2_00` names the collections it imports.
- `queries.get_annotations` / `get_subunits` / `get_ids_and_context` all `MATCH
  (node:ModelElement)`. BEL nodes used to carry `:BELModelElement` and **not** `:ModelElement`, so
  these returned nothing for them, silently — that is fixed at the source, in `2_00`, which now
  sets both labels. If a BEL query comes back empty for no reason, check the label is there before
  looking anywhere else. Note this only gets a node's *own* annotations: an `Activity` or `Complex`
  has none, and `get_subunits` traverses `HAS_SUBUNIT`, which no BEL node has — use
  `get_annotated_nodes` for the BEL `HAS__*` routes.
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
  compartment (`bel_terms` invariant (d)), so no `pd2af` change is needed — but it is a real bug
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
