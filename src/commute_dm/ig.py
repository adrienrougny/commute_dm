"""Selection of sub-maps out of the **stored** activity-flow CellDesigner maps.

The DB holds activity-flow maps (`COVID_DM_CD_AF`, `PD_DM_CD_AF`) produced by
`pd2af` in `keep-reactions` mode. An AF map has 0 reactions and its
`Modulation`s point straight at `Species`, so it already *contains*, as stored
objects, everything an influence-graph round-trip used to reconstruct: species,
signed modulations, boolean logic gates, their glyphs and arcs, and the
`LayoutModelMapping` tying them together.

This module therefore **selects and assembles** stored elements instead of
deriving edges and re-synthesizing model elements and geometry:

1. :func:`load_af_index` -- one element-id-only adjacency index over the stored
   structure, built once per run.
2. :func:`select_nodes` / :func:`close_over_gates` /
   :func:`induced_modulations` -- level-bounded BFS, then the *induced* element
   set.
3. :func:`make_celldesigner_map_from_selection` -- hydrate the selected stored
   elements through one shared cache and assemble them into a map.

**The load-bearing invariant is object identity.** `LayoutModelMapping` extends
`FrozenIdentitySurjectionDict`, so its inverse is keyed by `id(value)`:
`get_mapping(model_element)` only resolves for the *same object* stored as a
mapping value, and the CellDesigner writer relies on that
(`make_celldesigner_modulation_reaction` compares with `is`, then falls back to
`get_layouts`; when that yields nothing it emits `alias=""`, which the reader
chokes on). Hence: hydrate everything through one `execute_query_as_objects`
call path sharing a single `node_id_to_object` cache, convert once to builders
with one shared `object_to_builder` cache, and express every transformation as
in-place mutation of that one builder graph. Never construct a replacement
element.

Note that `session.cypher_query_as_layout_elements` must **not** be used: it has
no `node_id_to_object` parameter and so mints fresh objects.
"""

import momapy.builder
import momapy.celldesigner
import momapy.core.layout
import momapy.core.mapping
import momapy.coloring
import momapy.geometry
import momapy.positioning


# ---------------------------------------------------------------------------
# adjacency index over the stored AF structure
# ---------------------------------------------------------------------------

# Signed / directed `Modulation` subclasses we traverse. Plain `Modulation`
# (sign unspecified) and every `Unknown*` arc are deliberately excluded. Note
# that matching the `:Modulation` label already excludes the `Unknown*` classes
# for free: `UnknownModulation` is a *sibling* of `Modulation` under
# `KnownOrUnknownModulation`, not a subclass.
KEPT_MODULATION_LABELS = (
    "PositiveInfluence",
    "NegativeInfluence",
    "Triggering",
    "Inhibition",
    "PhysicalStimulation",
    "Catalysis",
)

# Labels that carry no class information; ignored when classifying.
_STRUCTURAL_LABELS = frozenset(
    {
        "BaseNode",
        "MapElement",
        "ModelElement",
        "CellDesignerModelElement",
        "KnownOrUnknownModulation",
        "Modulation",
    }
)

# Known but deliberately-excluded concrete classes. Anything outside this set
# and `KEPT_MODULATION_LABELS` raises, rather than being silently skipped:
# silent skipping is exactly how the old `make_ig_in_db` rotted into dead code
# after momapy renamed relationships under it.
_KNOWN_EXCLUDED_LABELS = frozenset(
    {
        "UnknownModulation",
        "UnknownPositiveInfluence",
        "UnknownNegativeInfluence",
        "UnknownInhibition",
        "UnknownTriggering",
        "UnknownPhysicalStimulation",
        "UnknownCatalysis",
    }
)


class AfIndex:
    """Element-id-only traversal index over the stored AF structure.

    Nodes are `Species` **and** `BooleanLogicGate`s, so one influence through a
    gate is two hops. The index holds no momapy objects -- it is a traversal
    index over the stored modulations, not a re-typed representation of them.
    """

    def __init__(self):
        self.species = set()
        self.species_to_collections = {}
        self.gates = set()
        self.out_edges = {}
        self.in_edges = {}
        self.modulation_endpoints = {}
        self.gate_inputs = {}
        self.modulation_class = {}
        self.n_excluded_modulations = 0

    def _add_edge(self, source_id, target_id, modulation_id):
        self.out_edges.setdefault(source_id, set()).add((target_id, modulation_id))
        self.in_edges.setdefault(target_id, set()).add((source_id, modulation_id))
        self.modulation_endpoints[modulation_id] = (source_id, target_id)

    def neighbours(self, node_id, mode):
        if mode == "downstream":
            return self.out_edges.get(node_id, ())
        if mode == "upstream":
            return self.in_edges.get(node_id, ())
        raise ValueError(f"unknown mode {mode!r}")

    def get_collections(self, node_id):
        return self.species_to_collections.get(node_id, frozenset())


def _classify_modulation(labels, element_id):
    """Return the kept concrete class label, or `None` if deliberately excluded."""
    concrete = [label for label in labels if label not in _STRUCTURAL_LABELS]
    kept = [label for label in concrete if label in KEPT_MODULATION_LABELS]
    if kept:
        if len(kept) > 1:
            raise RuntimeError(
                f"modulation {element_id} matched several kept classes: {kept}"
            )
        return kept[0]
    if not concrete:
        # plain `Modulation`: inside the :Modulation label but sign unspecified
        return None
    if all(label in _KNOWN_EXCLUDED_LABELS for label in concrete):
        return None
    raise RuntimeError(
        f"unrecognised modulation class for {element_id}: {sorted(concrete)}. "
        "Add it to KEPT_MODULATION_LABELS or _KNOWN_EXCLUDED_LABELS."
    )


_SPECIES_QUERY = """
MATCH (collection:Collection)-[:HAS_ENTRY]->(:CollectionEntry)
    -[:HAS_OBJ]->(:CellDesignerMap)-[:HAS_MODEL]->(:CellDesignerModel)
    -[:HAS_MODEL_ELEMENT]->(species:Species)
WHERE collection.name IN $collection_names
RETURN elementId(species) AS species_id,
       collect(DISTINCT collection.name) AS collection_names
"""

_MODULATION_QUERY = """
MATCH (collection:Collection)-[:HAS_ENTRY]->(:CollectionEntry)
    -[:HAS_OBJ]->(:CellDesignerMap)-[:HAS_MODEL]->(:CellDesignerModel)
    -[:HAS_MODEL_ELEMENT]->(modulation:Modulation)
WHERE collection.name IN $collection_names
MATCH (modulation)-[:HAS_SOURCE]->(source)
MATCH (modulation)-[:HAS_TARGET]->(target)
RETURN DISTINCT elementId(modulation) AS modulation_id,
       labels(modulation) AS modulation_labels,
       elementId(source) AS source_id,
       elementId(target) AS target_id
"""

# Gate inputs go through a `SimpleSpeciesReference`, not straight to the
# species -- getting this path wrong is what made the old gate queries match 0
# rows.
_GATE_QUERY = """
MATCH (collection:Collection)-[:HAS_ENTRY]->(:CollectionEntry)
    -[:HAS_OBJ]->(:CellDesignerMap)-[:HAS_MODEL]->(:CellDesignerModel)
    -[:HAS_MODEL_ELEMENT]->(gate:BooleanLogicGate)
WHERE collection.name IN $collection_names
MATCH (gate)-[:HAS_INPUT]->(:SimpleSpeciesReference)
    -[:HAS_REFERRED_ELEMENT]->(species:Species)
RETURN elementId(gate) AS gate_id,
       collect(DISTINCT elementId(species)) AS input_ids
"""


def load_af_index(session, collection_names):
    """Build the adjacency index for the given AF collections (once per run)."""
    index = AfIndex()
    params = {"collection_names": list(collection_names)}

    for row in session.execute_query(_SPECIES_QUERY, params):
        index.species.add(row["species_id"])
        index.species_to_collections.setdefault(row["species_id"], set()).update(
            row["collection_names"]
        )

    for row in session.execute_query(_GATE_QUERY, params):
        index.gates.add(row["gate_id"])
        index.gate_inputs.setdefault(row["gate_id"], set()).update(row["input_ids"])
        # a species feeding a gate is an edge species -> gate (no modulation node)
        for species_id in row["input_ids"]:
            index.out_edges.setdefault(species_id, set()).add((row["gate_id"], None))
            index.in_edges.setdefault(row["gate_id"], set()).add((species_id, None))

    for row in session.execute_query(_MODULATION_QUERY, params):
        class_label = _classify_modulation(
            row["modulation_labels"], row["modulation_id"]
        )
        if class_label is None:
            index.n_excluded_modulations += 1
            continue
        index.modulation_class[row["modulation_id"]] = class_label
        index._add_edge(row["source_id"], row["target_id"], row["modulation_id"])

    return index


def select_nodes(index, seed_ids, mode, max_level):
    """Level-bounded BFS from `seed_ids` over the index. Returns element ids."""
    reached = set(seed_ids)
    frontier = set(seed_ids)
    level = 0
    while frontier and (max_level < 0 or level < max_level):
        next_frontier = set()
        for node_id in frontier:
            for neighbour_id, _ in index.neighbours(node_id, mode):
                if neighbour_id not in reached:
                    reached.add(neighbour_id)
                    next_frontier.add(neighbour_id)
        frontier = next_frontier
        level += 1
    return reached


def induced_modulations(index, node_ids):
    """Every stored modulation whose source and target are both in `node_ids`.

    The modulation set must be *induced*, not the BFS tree: that is what
    `apoc.path.subgraphAll` did, and a tree would silently lose arcs.
    """
    return {
        modulation_id
        for modulation_id, (source_id, target_id) in index.modulation_endpoints.items()
        if source_id in node_ids and target_id in node_ids
    }


def close_over_gates(index, node_ids):
    """A retained gate pulls in *all* of its input species.

    Without this the writer emits a gate with a dangling input.
    """
    closed = set(node_ids)
    for node_id in [node_id for node_id in node_ids if node_id in index.gates]:
        closed |= index.gate_inputs.get(node_id, set())
    return closed


# ---------------------------------------------------------------------------
# hydration of the selected stored elements
# ---------------------------------------------------------------------------

# One arc per modulation, taken from the stored frozenset key {arc, src, tgt}.
_MODULATION_MAPPING_QUERY = """
MATCH (:LayoutModelMapping)-[:HAS_ITEM]->(item:Item)-[:HAS_VALUE]->(modulation:Modulation),
      (item)-[:HAS_KEY]->(key:FrozenSet)-[:HAS_ITEM]->(arc:Arc)
WHERE elementId(modulation) IN $modulation_ids
RETURN modulation AS modulation, arc AS arc
"""

# One glyph per model element (species or boolean logic gate). Restricted to
# TOP-LEVEL glyphs: a glyph whose parent is another CellDesignerNode is a
# complex subunit or a state variable and cannot stand alone in the output
# layout. The choice is made HERE, in Cypher, so the subunit query below can key
# off exactly the same glyph -- one source of truth for "which alias did we
# keep". A species is an arc endpoint in several source maps via a different
# alias each time; reusing every alias would draw the same protein repeatedly.
_CHOSEN_GLYPH_MATCH = """
MATCH (:LayoutModelMapping)-[:HAS_ITEM]->(item:Item)-[:HAS_VALUE]->(element),
      (item)-[:HAS_KEY]->(key:CellDesignerNode)
WHERE elementId(element) IN $element_ids
  AND NOT EXISTS { MATCH (:CellDesignerNode)-[:HAS_LAYOUT_ELEMENT]->(key) }
WITH element, key ORDER BY key.id_, elementId(key)
WITH element, head(collect(key)) AS chosen
"""

_GLYPH_QUERY = _CHOSEN_GLYPH_MATCH + """
RETURN element AS element, chosen AS chosen
"""

# Complex-subunit glyph pairs, read FORWARD out of the stored
# LayoutModelMapping (alias -> species) for the glyphs we kept. The stored
# mapping already knows which child glyph is which subunit, so this replaces the
# old max-overlap alias heuristic -- but the entries are still needed: without
# them the writer references child-glyph aliases it never defines.
_SUBUNIT_QUERY = _CHOSEN_GLYPH_MATCH + """
MATCH (chosen)-[:HAS_LAYOUT_ELEMENT*1..]->(child:CellDesignerNode)
MATCH (:LayoutModelMapping)-[:HAS_ITEM]->(child_item:Item)-[:HAS_VALUE]->(subunit:Species),
      (child_item)-[:HAS_KEY]->(child)
RETURN child AS child, subunit AS subunit
"""

# The stored top-level glyph of every compartment the selected species live in,
# closed over `outside` -- the same "one glyph per element, chosen in Cypher"
# rule as the species glyphs above. 429 of the 431 stored compartments have one;
# the two without are the `default` roots, which are not drawn anyway.
_COMPARTMENT_GLYPH_QUERY = """
MATCH (species) WHERE elementId(species) IN $element_ids
MATCH (species)-[:HAS_COMPARTMENT]->(:Compartment)-[:HAS_OUTSIDE*0..]->(element:Compartment)
WITH DISTINCT element
MATCH (:LayoutModelMapping)-[:HAS_ITEM]->(item:Item)-[:HAS_VALUE]->(element),
      (item)-[:HAS_KEY]->(key:CellDesignerNode)
WHERE NOT EXISTS { MATCH (:CellDesignerNode)-[:HAS_LAYOUT_ELEMENT]->(key) }
WITH element, key ORDER BY key.id_, elementId(key)
WITH element, head(collect(key)) AS chosen
RETURN element AS element, chosen AS chosen
"""

_ELEMENTS_BY_ID_QUERY = """
UNWIND $element_ids AS element_id
MATCH (element) WHERE elementId(element) = element_id
RETURN element AS element
"""


def _iter_templates(species):
    """The template of `species` and of every transitive subunit.

    The writer iterates `model.species_templates`, not the species'
    back-reference, and silently drops a template class it has no branch for --
    so every template must be an explicit member of the model.
    """
    template = getattr(species, "template", None)
    if template is not None:
        yield template
    for subunit in getattr(species, "subunits", ()) or ():
        yield from _iter_templates(subunit)


def _walk_builders(root):
    """Every builder reachable from `root`, by identity."""
    seen = {}
    stack = [root]
    while stack:
        obj = stack.pop()
        if isinstance(obj, momapy.builder.Builder):
            if id(obj) in seen:
                continue
            seen[id(obj)] = obj
            for field in obj.__dataclass_fields__:
                stack.append(getattr(obj, field, None))
        elif isinstance(obj, (list, tuple, set, frozenset)):
            stack.extend(obj)
        elif isinstance(obj, dict):
            stack.extend(obj.keys())
            stack.extend(obj.values())
    return seen


class Selection:
    """The hydrated stored elements chosen for one output map."""

    def __init__(self):
        self.modulation_to_arc = {}
        self.element_to_glyph = {}
        self.subunit_entries = []
        self.dropped_no_glyph = set()


def hydrate(session, cache, element_ids, modulation_ids):
    """One hydration pass over a single shared cache -> object identity holds."""
    selection = Selection()

    if element_ids:
        params = {"element_ids": sorted(element_ids)}
        for element, glyph in session.execute_query_as_objects(
            _GLYPH_QUERY, params, node_id_to_object=cache
        ):
            selection.element_to_glyph[element] = glyph
        for child, subunit in session.execute_query_as_objects(
            _SUBUNIT_QUERY, params, node_id_to_object=cache
        ):
            selection.subunit_entries.append((child, subunit))

    if modulation_ids:
        for modulation, arc in session.execute_query_as_objects(
            _MODULATION_MAPPING_QUERY,
            {"modulation_ids": sorted(modulation_ids)},
            node_id_to_object=cache,
        ):
            # one arc per modulation: a modulation drawn in several source maps
            # has several stored frozenset keys and would otherwise be duplicated
            if modulation not in selection.modulation_to_arc:
                selection.modulation_to_arc[modulation] = arc

    # drop modulations whose endpoints have no top-level glyph
    for modulation in list(selection.modulation_to_arc):
        for endpoint in (modulation.source, modulation.target):
            if endpoint is None or endpoint not in selection.element_to_glyph:
                selection.dropped_no_glyph.add(modulation)
                del selection.modulation_to_arc[modulation]
                break

    return selection


def hydrate_elements(session, cache, element_ids):
    """Hydrate bare elements by element id through the shared cache.

    Returns the very objects :func:`hydrate` produced -- the shared cache is
    what makes them identical, and identity is what the mapping is keyed on.
    The element ids themselves are *not* recoverable from the result:
    `execute_query_as_objects` drops non-node columns, so a query cannot return
    an id alongside its object. Callers that need the correspondence therefore
    query one group at a time (see :func:`_resolve_group`) rather than relying
    on row order.
    """
    if not element_ids:
        return []
    return [
        row[0]
        for row in session.execute_query_as_objects(
            _ELEMENTS_BY_ID_QUERY,
            {"element_ids": sorted(element_ids)},
            node_id_to_object=cache,
        )
        if row
    ]


def _resolve_group(session, cache, keys):
    """Resolve a group of keys (element ids and/or objects) to elements."""
    element_ids = {key for key in keys if isinstance(key, str)}
    resolved = [key for key in keys if not isinstance(key, str)]
    resolved += hydrate_elements(session, cache, element_ids)
    return resolved


def _resolve_one(session, cache, key):
    """Resolve a single key (an element id or an object) to one element."""
    if not isinstance(key, str):
        return key
    elements = hydrate_elements(session, cache, [key])
    if len(elements) != 1:
        raise RuntimeError(f"element id {key!r} resolved to {len(elements)} elements")
    return elements[0]


# ---------------------------------------------------------------------------
# assembly
# ---------------------------------------------------------------------------


def _compartment_closure(species_builders):
    """Every compartment reachable from the retained species, following `outside`.

    Compartments are **kept**, not collapsed: collapsing them makes species that
    differ only by compartment content-equal (momapy excludes `id_` from
    eq/hash) and `model.species` is a frozenset, so it silently loses them.

    The set must be closed under `outside`, though. The reader's
    `get_ordered_compartment_aliases` topologically orders compartments by the
    `outside` chain, repeatedly looking for one with `outside is None`; if every
    remaining compartment's `outside` points at a compartment *not* in the set,
    its inner loop never breaks, `to_delete` keeps a stale value and `del`
    raises `KeyError: '<compartment id>'`. The closure terminates because a
    CellDesigner map always carries a root `default` compartment with
    `outside=None`. Compartment *glyphs* are irrelevant -- sorting an empty
    alias list is fine.
    """
    closed = {}
    stack = []
    for builder in species_builders:
        compartment = getattr(builder, "compartment", None)
        if compartment is not None:
            stack.append(compartment)
    while stack:
        compartment = stack.pop()
        if id(compartment) in closed:
            continue
        closed[id(compartment)] = compartment
        outside = getattr(compartment, "outside", None)
        if outside is not None:
            stack.append(outside)
    return closed


def _intern(builders):
    """Map each builder key to the canonical builder for its *content*.

    Returns `({key: canonical_builder}, {id(canonical): key})`. The content key
    is the frozen form of the builder, so it uses momapy's own eq/hash -- which
    exclude `id_`, and are therefore exactly the equality a `frozenset` field
    will apply. With compartments kept, this normally interns nothing away (the
    selected elements are distinct as stored); it is the guard that keeps
    `model.species` / `model.modulations` and the layout in agreement should two
    selected elements ever be content-equal, since the frozenset would otherwise
    drop one while the layout kept a glyph for it. Freezing uses a throwaway
    cache so it does not pin the builder's final frozen form.
    """
    by_content = {}
    canonical = {}
    canonical_key = {}
    for key in sorted(builders, key=lambda key: (builders[key].id_ or "")):
        builder = builders[key]
        content = momapy.builder.object_from_builder(builder, builder_to_object={})
        chosen = by_content.setdefault(content, builder)
        canonical[key] = chosen
        canonical_key.setdefault(id(chosen), key)
    return canonical, canonical_key


def _renumber_ids(map_builder):
    """Renumber every builder `id_` in the graph, except the literal `default` id.

    Ids are local to a source file, and a selection spans dozens of files: 16
    compartment ids collide over 40 nodes, and layout-element ids collide too
    (`pdme90` is a free-standing alias in one map and a complex-subunit alias in
    another). The writer keys the xml on `id_` and the reader indexes by it, so
    a collision makes the reader resolve one id to two elements. Model elements
    must be covered as well, not just the layout tree (templates and their
    modification residues collide the same way).

    The literal `default` compartment id is excluded: the writer and reader
    special-case it, and renumbering it produces `KeyError: ''`. Real
    compartments are renumbered like everything else.
    """
    counter = 0
    for builder in _walk_builders(map_builder).values():
        if getattr(builder, "id_", None) == "default":
            continue
        if "id_" not in getattr(builder, "__dataclass_fields__", {}):
            continue
        builder.id_ = f"e{counter}"
        counter += 1
    return counter


def make_celldesigner_map_from_selection(
    session,
    cache,
    element_ids,
    modulation_ids,
    color_element_ids=None,
    extra_elements=None,
    extra_influences=None,
    with_compartment_layouts=False,
):
    """Assemble a `CellDesignerMap` from the selected stored AF elements.

    Args:
        session: the momapy_kb session.
        cache: the run-wide `node_id_to_object` hydration cache. Hoisted to the
            whole run for *correctness* as much as for speed -- see the module
            docstring.
        element_ids: element ids of the selected species and boolean logic
            gates.
        modulation_ids: element ids of the selected (induced) modulations.
        color_element_ids: list of `(keys, color_name)`; each key is either an
            element id or one of the `extra_elements` model elements. The named
            `momapy.coloring` colour is applied to that element's glyph.
        extra_elements: list of `(model_element, layout_element)` pairs to add
            to the map on top of the stored selection -- the hook for elements
            that have no stored counterpart (the synthetic central node; later,
            BEL nodes, which have neither layout nor model element).
        extra_influences: list of `(source, target, modulation_class,
            arc_class)`, where source and target are either an element id or an
            extra model element. These are the only *synthesized* modulations;
            their geometry is a placeholder, recomputed by
            `pd2af.utils.make_auto_layout`.
        with_compartment_layouts: draw the compartments. Compartments are always
            in the *model*; without this the layout has no compartment glyph, so
            no box is drawn and -- because `make_auto_layout` builds one dot
            cluster per *mapped* compartment -- species of one compartment are
            not grouped either. With it, each compartment's stored top-level
            glyph is selected like any other element, and species cluster inside
            their compartment, nested by `outside`. Note that compartments are
            not merged by name: two source maps each contributing a `microglia`
            compartment give two boxes.

    Returns:
        `(cd_map, stats)`. The returned map still has the stored arc geometry:
        stored geometry is not reusable across a merge, so the caller must run
        `pd2af.utils.make_auto_layout` on it.
    """
    color_element_ids = color_element_ids or []
    extra_elements = extra_elements or []
    extra_influences = extra_influences or []

    selection = hydrate(session, cache, element_ids, modulation_ids)
    # A species drawn *only* inside complexes has no top-level glyph, so it
    # cannot stand alone in the output layout and is left out (along with its
    # arcs, see `hydrate`). Counted rather than silently dropped.
    n_elements_without_glyph = len(set(element_ids)) - len(selection.element_to_glyph)
    for model_element, layout_element in extra_elements:
        selection.element_to_glyph[model_element] = layout_element

    # ---- one shared object->builder cache over model, layout and mapping ----
    object_to_builder = {}

    def to_builder(obj):
        return momapy.builder.builder_from_object(
            obj, object_to_builder=object_to_builder
        )

    element_builders = {}
    glyph_builders = {}
    for element, glyph in selection.element_to_glyph.items():
        element_builders[id(element)] = to_builder(element)
        glyph_builders[id(element)] = to_builder(glyph)

    template_builders = {}
    for element in selection.element_to_glyph:
        for template in _iter_templates(element):
            template_builders[id(template)] = to_builder(template)

    gate_builder_cls = momapy.builder.get_or_make_builder_cls(
        momapy.celldesigner.BooleanLogicGate
    )

    compartment_builders = _compartment_closure(
        [
            builder
            for builder in element_builders.values()
            if not isinstance(builder, gate_builder_cls)
        ]
    )
    # `default` is the one id `_renumber_ids` leaves alone (the writer and reader
    # special-case it), so two root compartments carrying it would be written as
    # two `<compartment id="default">` and the reader would silently keep one --
    # losing a compartment, and resolving one root's species to the other. That
    # happens when the source maps disagree on the root's attributes, so it is
    # fixed in the maps; this only makes the breakage loud if it comes back.
    n_default_roots = sum(
        1
        for builder in compartment_builders.values()
        if getattr(builder, "id_", None) == "default"
    )
    if n_default_roots > 1:
        raise RuntimeError(
            f"{n_default_roots} distinct root compartments with id 'default' in "
            "one selection: the source maps disagree on the root compartment's "
            "attributes, so they were not merged on import. Make them identical "
            "and re-run 2_00_save_collections."
        )

    # Compartment glyphs are selected, not synthesized -- same shared cache, so
    # the hydrated compartment is the very object the species already point at.
    compartment_glyph_builders = {}
    if with_compartment_layouts and element_ids:
        for compartment, glyph in session.execute_query_as_objects(
            _COMPARTMENT_GLYPH_QUERY,
            {"element_ids": sorted(element_ids)},
            node_id_to_object=cache,
        ):
            # `compartment_builders` is keyed by the *builder*'s id, and the
            # shared `object_to_builder` cache maps this hydrated compartment to
            # the very builder the species already point at.
            compartment_b = to_builder(compartment)
            if id(compartment_b) not in compartment_builders:
                continue  # not on the retained species' `outside` chains
            compartment_glyph_builders[id(compartment_b)] = to_builder(glyph)

    # ---- in-place mutation only, from here on ----
    canonical_element, canonical_key = _intern(element_builders)
    # Use the canonical element's OWN glyph. Taking some other duplicate's glyph
    # would break the complex-subunit mapping, which is identity-based: the
    # child glyph's stored subunit is an object under *that* duplicate, not
    # under the canonical one.
    canonical_glyph = {
        id(builder): glyph_builders[canonical_key[id(builder)]]
        for builder in canonical_element.values()
    }

    # Every reference to a species must point at the canonical instance, not
    # just the ones in `model.species`: the writer resolves each reference by
    # object identity, so a stale duplicate yields an empty alias. Gate inputs
    # are the easy one to miss -- they produce `KeyError: ''` on read-back.
    builder_to_canonical = {
        id(builder): canonical_element[key] for key, builder in element_builders.items()
    }
    for builder in canonical_element.values():
        if not isinstance(builder, gate_builder_cls):
            continue
        for gate_input in builder.inputs or ():
            referred = getattr(gate_input, "referred_element", None)
            canonical = builder_to_canonical.get(id(referred))
            if canonical is not None:
                gate_input.referred_element = canonical

    # ---- influences: the stored ones, plus any synthesized extra ones ----
    influence_records = []
    for modulation, arc in selection.modulation_to_arc.items():
        influence_records.append(
            (
                to_builder(modulation),
                to_builder(arc),
                modulation.source,
                modulation.target,
            )
        )
    n_extra_influences_dropped = 0
    for source, target, modulation_class, arc_class in extra_influences:
        source_element = _resolve_one(session, cache, source)
        target_element = _resolve_one(session, cache, target)
        if (
            source_element not in selection.element_to_glyph
            or target_element not in selection.element_to_glyph
        ):
            # an endpoint with no top-level glyph is not in the output map
            n_extra_influences_dropped += 1
            continue
        influence_records.append(
            (
                momapy.builder.new_builder_object(modulation_class),
                momapy.builder.new_builder_object(arc_class),
                source_element,
                target_element,
            )
        )

    # re-point modulations at the canonical species, then intern them too: two
    # modulations that differed only by which duplicate species they touched are
    # now equal, and would collapse in `model.modulations` the same way.
    modulation_builders = {}
    for modulation_b, _, source_element, target_element in influence_records:
        modulation_b.source = canonical_element[id(source_element)]
        modulation_b.target = canonical_element[id(target_element)]
        modulation_builders[id(modulation_b)] = modulation_b
    canonical_modulation, _ = _intern(modulation_builders)

    # re-point arcs at the chosen glyph. Stored arc *geometry* is not reusable
    # across a merge, but the stored arc *object* (class, style, mapping) is;
    # `make_auto_layout` recomputes the segments.
    mapping_entries = []
    kept_arcs = {}
    kept_modulations = {}
    for modulation_b, arc_b, source_element, target_element in influence_records:
        modulation_b = canonical_modulation[id(modulation_b)]
        if id(modulation_b) in kept_arcs:
            continue  # one arc per canonical modulation
        source_glyph_b = canonical_glyph[id(canonical_element[id(source_element)])]
        target_glyph_b = canonical_glyph[id(canonical_element[id(target_element)])]
        arc_b.source = source_glyph_b
        arc_b.target = target_glyph_b
        if not arc_b.segments:
            # placeholder; `make_auto_layout` rebuilds every arc's segments
            arc_b.segments = [
                momapy.geometry.Segment(
                    source_glyph_b.center(), target_glyph_b.center()
                )
            ]
        kept_arcs[id(modulation_b)] = arc_b
        kept_modulations[id(modulation_b)] = modulation_b
        mapping_entries.append((arc_b, modulation_b, source_glyph_b, target_glyph_b))

    canonical_elements = {id(b): b for b in canonical_element.values()}
    species_builders = {
        key: builder
        for key, builder in canonical_elements.items()
        if not isinstance(builder, gate_builder_cls)
    }
    gate_builders = {
        key: builder
        for key, builder in canonical_elements.items()
        if isinstance(builder, gate_builder_cls)
    }

    # ---- assemble model and layout ----
    model_b = momapy.builder.new_builder_object(momapy.celldesigner.CellDesignerModel)
    model_b.species = frozenset(species_builders.values())
    model_b.species_templates = frozenset(template_builders.values())
    model_b.compartments = frozenset(compartment_builders.values())
    model_b.modulations = frozenset(kept_modulations.values())
    model_b.boolean_logic_gates = frozenset(gate_builders.values())

    layout_b = momapy.builder.new_builder_object(
        momapy.celldesigner.CellDesignerLayout
    )
    # Compartments first: they render an opaque interior, so a later one would
    # hide the species drawn inside it (pd2af orders them the same way).
    layout_b.layout_elements = (
        tuple(compartment_glyph_builders.values())
        + tuple(canonical_glyph.values())
        + tuple(kept_arcs.values())
    )

    n_colored = 0
    for keys, color_name in color_element_ids:
        color = getattr(momapy.coloring, color_name)
        for element in _resolve_group(session, cache, keys):
            if element not in selection.element_to_glyph:
                continue
            glyph_b = canonical_glyph.get(id(canonical_element[id(element)]))
            if glyph_b is not None:
                glyph_b.fill = color
                n_colored += 1

    map_b = momapy.builder.new_builder_object(momapy.celldesigner.CellDesignerMap)
    map_b.model = model_b
    map_b.layout = layout_b

    # ---- id renumber over the WHOLE builder graph ----
    n_renumbered = _renumber_ids(map_b)

    # ---- mapping, with anchors re-registered ----
    #
    # `LayoutModelMappingPlugin` only persists `obj.items()`, and
    # `_singleton_to_key` is not a dataclass field, so stored anchors are never
    # saved and always come back empty. `get_mapping`'s third lookup step is the
    # anchor fallback, and the writer's `frozenset_mapping=None` path scans for a
    # *non*-frozenset key -- so a frozenset-keyed modulation with no anchor
    # resolves to no layout at all. Anchors must therefore be re-registered here.
    mapping_b = momapy.core.mapping.LayoutModelMappingBuilder()
    for builder in canonical_element.values():
        mapping_b.add_mapping(canonical_glyph[id(builder)], builder)
    for key, glyph_b in compartment_glyph_builders.items():
        # `make_auto_layout` looks compartments up through the mapping to build
        # its dot clusters, so an unmapped compartment glyph would be a box that
        # nothing is placed inside.
        mapping_b.add_mapping(glyph_b, compartment_builders[key])
    n_subunit_entries = 0
    for child, subunit in selection.subunit_entries:
        child_b = object_to_builder.get(id(child))
        subunit_b = object_to_builder.get(id(subunit))
        if child_b is None or subunit_b is None:
            continue
        mapping_b.add_mapping(child_b, subunit_b)
        n_subunit_entries += 1
    for arc_b, modulation_b, source_glyph_b, target_glyph_b in mapping_entries:
        mapping_b.add_mapping(
            frozenset({arc_b, source_glyph_b, target_glyph_b}),
            modulation_b,
            anchor=arc_b,
        )
    map_b.layout_model_mapping = mapping_b

    momapy.positioning.set_fit(
        layout_b, layout_b.layout_elements, xsep=15.0, ysep=15.0
    )
    cd_map = momapy.builder.object_from_builder(map_b)
    stats = {
        "n_species": len(cd_map.model.species),
        "n_modulations": len(cd_map.model.modulations),
        "n_gates": len(cd_map.model.boolean_logic_gates),
        "n_compartments": len(cd_map.model.compartments),
        "n_compartment_layouts": len(compartment_glyph_builders),
        "n_templates": len(cd_map.model.species_templates),
        "n_renumbered": n_renumbered,
        "n_elements_without_glyph": n_elements_without_glyph,
        "n_modulations_dropped_no_glyph": len(selection.dropped_no_glyph),
        "n_extra_influences_dropped": n_extra_influences_dropped,
        "n_subunit_entries": n_subunit_entries,
        "n_interned_away": len(element_builders) - len(canonical_glyph),
        "n_colored": n_colored,
    }
    return cd_map, stats
