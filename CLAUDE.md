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

Every notebook starts with `%run 0_10_load_paths.ipynb`, which restores path constants saved by
`0_00` via IPython's `%store`. So **`0_00_make_paths` must be run once first** to populate the
store; afterwards any notebook can be opened standalone.

- `0_00` / `0_10` — define and load all `*_DATA_DIR` / `*_FILE` path constants (rooted at `data/`).
- `1_00_make_bel_files` — prepare BEL KG inputs.
- `2_00_save_collections` — wipe the DB, import BEL KG cypher dumps (`cypher-shell`), then load
  the CellDesigner XML maps into named **Collections** (`COVID_DM_CD`, `PD_DM_CD`, plus
  `*_KG_BEL`). Uses `with_membership_edges=True` (the `HAS_MODEL_ELEMENT` edge that
  `get_collections_for_nodes` relies on) and `integration_mode="hash"`.
- `2_10_make_ig` — `commute_dm.ig.make_ig_in_db(session)`: materializes the influence graph as
  **direct typed relationships between Species nodes** (see below).
- `3_00_get_interfaces` — query proteins present across maps (by `hgnc.symbol` annotation),
  write `data/.../interface/*.json` (`covid_pd`, `covid_ad`, `covid_pd_ad`).
- `4_00_make_goat_gene_lists` — turn raw RNA-seq DE CSVs (`data/rnaseq/`) into GOAT gene-list
  CSVs, mapping Ensembl→Entrez/symbol via the HGNC dataset.
- `4_10_make_interface_graphs` — render CellDesigner maps of interface subgraphs at several
  `MAX_LEVELS`.
- `4_20_make_interface_goat_analysis` — GOAT enrichment + intersection analyses of interface
  subgraphs vs. the gene lists.
- `5_00_make_annotations_for_bel` — annotation tooling (see `annotations/` helper scripts).
- `6_00` / `6_10` — immuno-target analysis and graphs.

## Library architecture (`src/commute_dm`)

The package has no `__init__` exports; modules are imported by full path. Data flows
**Neo4j → InfluenceGraph → CellDesignerMap / gene sets**.

### `ig.py` — the influence graph (core abstraction)
- Influence types are module constants: `REACTANT_TO_PRODUCT`, `POSITIVE_INFLUENCE`,
  `NEGATIVE_INFLUENCE`, `NECESSARY_POSITIVE_INFLUENCE`, plus `UNCERTAIN_*`. `INFLUENCES` is the
  subset used for subgraph traversal.
- `make_ig_in_db` writes these as real Neo4j edges between Species, derived from CellDesigner
  Reactions (reactant→product, modifiers by class) and BEL Modulations; BooleanLogicGate inputs
  are expanded (NotGates excluded). Each derived edge carries `model_ids` + `reaction_ids`/
  `modulation_ids` provenance.
- `InfluenceGraph(dict)` maps `node -> set(incoming relationships)`; `get_modulators` etc. read
  edge types. `make_ig_from_nodes_and_relationships` builds one from query results.
- **`make_celldesigner_map_from_ig`** is the most intricate code in the repo and mirrors the
  external `pd2af` builder. It rebuilds a `CellDesignerMap` (model + layout + mapping) from an IG.
  The hard, load-bearing invariant — read the long docstrings before touching any of this:
  momapy dataclasses are **frozen and exclude `id_` from equality/hashing**, so species/
  templates/modulations are **interned by content** (`_register_or_reuse`, `_canon_template`,
  `_normalize_species`) to one canonical instance; then layout-element `id_`s are renumbered
  globally unique (`_assign_unique_layout_ids`) because the writer keys XML on `id_` and the
  reader indexes by it. Layout positions come from graphviz (`pydot`, `rankdir=BT`). Compartments
  are collapsed into one `DEFAULT_COMPARTMENT` to avoid the reader's compartment-ordering crash.
  `_collect_complex_subunit_mappings` pairs complex-subunit glyphs by reading the source maps'
  stored `LayoutModelMapping` **forward** (alias→species) — see the auto-memory on the
  `[0][0]` alias-pick bug.

### `queries.py` — Neo4j read helpers
- `get_subgraph` wraps `apoc.path.subgraphAll` with `mode` (`downstream`/`upstream`/`all`),
  `relationship_types`, level bounds, blacklists. This is how interface subgraphs are extracted.
- `prewarm_session` **must be called** before hydrating CellDesigner objects from query results —
  it registers pylpg node classes for momapy subclasses the static type walk misses.
- `get_collections_for_nodes`, `get_subunits`, `get_identifiers`, `get_annotations`,
  `get_ids_and_context` — membership / cross-reference lookups.

### `core.py` — analysis orchestration
Ties queries + ig + gea together for the interface workflows:
`get_interface` (proteins in ≥3 collections), `make_and_write_cd_maps_from_interface` (renders
the upstream-COVID + downstream-PD subgraph around each interface protein, with a synthetic
central node and per-collection node coloring), `make_goat_analysis_from_interface` /
`_from_pd`, `make_intersection_analysis_from_interface`, `make_goat_gene_lists`.
Note `_adjust_max_level` subtracts 1 to compensate for an old dummy-seed scheme.

### `gea.py` — gene-set enrichment
Builds gene sets from node annotations (`ncbigene`/`hgnc.symbol`, optionally including complex
subunits), writes GMT, and runs the R `goat` package via `rpy2` (`make_goat_analysis`).

### `utils.py` — IO/cypher helpers
Directory (re)creation, CellDesigner XML annotation injection (`add_annotation_to_file`),
`subgraph_to_cypherl` (export a query result as a CREATE/MATCH cypher dump), `merge_relationship`.

## Conventions & gotchas

- Cypher queries are mostly **f-string interpolated** (including element-id lists). Keep that
  style for consistency, but never interpolate untrusted input.
- Collection names (`"COVID_DM_CD"`, `"PD_DM_CD"`, `"AD_KG_BEL"`, …) and the influence-type
  constants are referenced by string/identity across modules — change them in one place only.
- Map versions in use are recorded in `data/maps/versions.txt`.
- `poubelle/` is a scratch/trash dir; ignore it.
