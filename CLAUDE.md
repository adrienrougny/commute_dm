# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`commute_dm` is a research analysis pipeline (the "COMMUTE" disease-map project). It takes
CellDesigner disease maps for COVID-19 and Parkinson's disease (PD), plus BEL knowledge
graphs (AD/PD/COVID), loads them into a **Neo4j** graph DB via the `momapy_kb` stack, derives
an **influence graph** from the reaction/modulation network, finds the **interface** of
proteins shared across the maps, then runs **gene-set enrichment (GOAT)** and intersection
analyses against bulk RNA-seq DE results, and renders new CellDesigner maps of the
interface subgraphs.

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
  runtime). `rpy2` is commented out in `pyproject.toml` deps — GOAT analyses need a working R.
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
- `3_00_get_interfaces` — query proteins present across collections (joined on **UniProt RDF
  annotations**), write `data/.../interface/*.json` — four files (`covid_pd`, `covid_ad`,
  `covid_pd_ad`, and the AF three-way `covid_af_pd_af_ad`, which is the one `4_10`/`4_20` use).
- `4_00_make_goat_gene_lists` — turn raw RNA-seq DE CSVs (`data/rnaseq/`) into GOAT gene-list
  CSVs, mapping Ensembl→Entrez/symbol via the HGNC dataset.
- `4_10_make_interface_graphs` — load the influence structure and the drawing material once
  (`core.load_submap_inputs`, ~3 min), then assemble and write a CellDesigner map per interface
  protein at several `MAX_LEVELS`.
- `4_20_make_interface_goat_analysis` — GOAT enrichment + intersection analyses of interface
  subgraphs vs. the gene lists.
- `5_00_make_annotations_for_bel` — annotation tooling (see `annotations/` helper scripts).
- `6_00` / `6_10` — immuno-target analysis and graphs.

## Library architecture (`src/commute_dm`)

The package has no `__init__` exports; modules are imported by full path. Data flows
**Neo4j (stored AF maps) → element-id selection → CellDesignerMap / gene sets**.

### `submaps.py` — the comorbidity sub-maps (replaced `ig.py`)
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
  node today; BEL nodes later) and are the only synthesized modulations.
- Gotchas worth keeping in mind: `Node.position` / `width` / `height` have **no defaults**, so a
  layout built directly as a frozen object must be given geometry; and
  `pd2af…_MODULATION_CLASS_TO_LAYOUT_CLASS` has no entry for `Inhibition`, `Catalysis` or
  `PhysicalStimulation`, which `_make_modulation_arc` supplies locally.
- Output geometry: callers must run `pd2af.utils.make_auto_layout(cd_map)` — stored arc geometry
  is not reusable across a merge; it recomputes every position and segment and fits the root
  layout, which is why the assembly does not call `set_fit` itself.

### `queries.py` — Neo4j read helpers
- `prewarm_session` **must be called** before hydrating CellDesigner objects from query results —
  it registers pylpg node classes for momapy subclasses the static type walk misses.
- `get_collections_for_nodes`, `get_subunits`, `get_identifiers`, `get_annotations`,
  `get_ids_and_context` — membership / cross-reference lookups. Note `get_collections_for_nodes`
  traverses `(:CellDesignerMap)`, so it returns `{}` for BEL nodes.
- `get_nodes(session, element_ids)` — the bridge from the id-only traversal back to DB nodes, for
  the node-taking `gea` helpers.

### `core.py` — analysis orchestration
Ties queries + submaps + gea together for the interface workflows. `UPSTREAM_COLLECTION_NAME` /
`DOWNSTREAM_COLLECTION_NAME` / `INTERFACE_COLLECTION_NAMES` name the AF collections the analysis
runs on. `get_interface` joins collections on shared UniProt annotations;
`load_submap_inputs` returns the whole per-run setup (`influences`, `source_map`,
`node_id_to_object`) in one call; `make_and_write_submaps_from_interface` assembles and writes the
sub-map upstream of the upstream seeds and downstream of the downstream seeds around each
interface protein, with a synthetic central node and a fill per walk direction; then
`make_goat_analysis_from_interface` / `_from_pd`, `make_intersection_analysis_from_interface`,
`make_goat_gene_lists`.
`_split_interface_seeds` returns `"upstream"` / `"downstream"` keyed **node ids**, so pointing
the analysis at another pair of collections (e.g. COVID → AD) is a parameter change; passing
`node_id_to_object` and `source_map` also drops seeds that exist only as complex subunits, which
cannot seed a walk. `max_level` means hops, plainly. All three interface entry points take
`influences` and the two collection names, but only the map one takes the `source_map`. `import commute_dm.gea` is **lazy** (`_gea()`) because `gea` imports `rpy2` at module
scope, which fails wherever R is unusable (e.g. `.venv-sandbox`).

### `gea.py` — gene-set enrichment
Builds gene sets from node annotations (`ncbigene`/`hgnc.symbol`, optionally including complex
subunits), writes GMT, and runs the R `goat` package via `rpy2` (`make_goat_analysis`).

### `utils.py` — IO/cypher helpers
Directory (re)creation, CellDesigner XML annotation injection (`add_annotation_to_file`),
`subgraph_to_cypherl` (export a query result as a CREATE/MATCH cypher dump), `merge_relationship`.

## Conventions & gotchas

- Cypher queries are mostly **f-string interpolated** (including element-id lists). Keep that
  style for consistency, but never interpolate untrusted input.
- Collection names (`"COVID_DM_CD_AF"`, `"PD_DM_CD_AF"`, `"AD_KG_BEL"`, …) are referenced by
  string across modules — the ones the analysis uses live in `core.py` as
  `UPSTREAM_COLLECTION_NAME` / `DOWNSTREAM_COLLECTION_NAME` / `INTERFACE_COLLECTION_NAMES`.
- `6_10_make_immuno_targets_graphs.ipynb` is **stale**: it calls a
  `make_and_write_cd_maps_from_interface` that no longer exists, with seeds from the non-AF
  collections. It was already non-functional before the `submaps.py` rewrite (it rendered via the
  IG edges that no longer exist).
- Map versions in use are recorded in `data/maps/versions.txt`.
- `poubelle/` is a scratch/trash dir; ignore it.
