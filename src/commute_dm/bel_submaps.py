"""The BEL side of a comorbidity sub-map, as an *in-memory* `CellDesignerMap`.

`commute_dm.submaps` assembles a sub-map out of the elements a `source_map` stores.
A BEL knowledge graph has no CellDesigner representation at all, so this module
builds one: it loads the KG's **influence-graph projection** (the one
`2_05_get_collections_statistics` defines and documents) and turns it into a
`CellDesignerMap` of real activity-flow content -- proteins with phosphorylation
badges, complexes with nested subunit glyphs, `act(...)` as CellDesigner's active
decoration, `loc(...)` as a drawn compartment -- with the signed causal relations
as modulations.

That map is then merged into `source_map` with `commute_dm.submaps.merge_maps`, and
**everything downstream keeps working unchanged**. The reason is worth stating,
because it is what makes this module small: every assumption in `submaps` and
`core` about a selected node id is "it resolves, through `node_id_to_object`, to a
species of `source_map.model`". Keying the BEL species on their real Neo4j node
ids and merging them into the source map satisfies that assumption, so
`Influences`, `fills`, `extra_influences`, `close_over_gates`,
`get_layout_element_for_model_element` and `_split_interface_seeds`' `keep()` need
no BEL special case.

Nothing here writes to the database, and no CellDesigner file is produced for the
KG itself: the projection is rebuilt as objects in ~2 s, which is cheaper than
writing and re-reading a 17 MB XML, and far cheaper than laying one out
(`pd2af.utils.make_auto_layout` runs graphviz `dot`, which does not finish on a
3979-node graph).

**What a term becomes, and the four identity invariants that keep the result
readable, live in `commute_dm.bel_terms`.** This module is only the assembly:
which nodes are projected, how they are wired, and how the resulting map merges
into the source map.
"""

import collections

import momapy.celldesigner
import momapy.core.mapping
import momapy.geometry
import pd2af.celldesigner.building_model

import commute_dm.bel_terms
import commute_dm.submaps


# The BEL node classes the projection keeps -- the validation set for
# :func:`load_bel_projection`. What each is *drawn* as is
# `commute_dm.bel_terms.BEL_CLASS_TO_CD_CLASSES`, which is a different question:
# an `Activity` is projected but has no class of its own, taking its subject's.
#
# The interface is joined on the annotations of the BEL collection's own
# `:Protein` nodes (`commute_dm.core.get_interface`), so nothing here decides
# which proteins are in the interface or which nodes can seed a walk.

# The signed BEL causal relations. `REGULATES` is deliberately absent: it carries
# no sign, and `commute_dm.submaps.SIGNED_MODULATION_CLASSES` -- the filter for the
# modulations induced from `source_map` -- excludes the unsigned `Modulation`, so a
# `REGULATES` edge mapped to one would be *walked but never drawn*, giving sub-maps
# with invisible connections. Leaving it out of both the map and the walk keeps the
# two consistent, at the cost of 2.6% of the AD KG's causal edges.
BEL_RELATION_TO_MODULATION_CLASS = {
    "INCREASES": momapy.celldesigner.PositiveInfluence,
    "DIRECTLY_INCREASES": momapy.celldesigner.PositiveInfluence,
    "TRANSLATED_TO": momapy.celldesigner.PositiveInfluence,
    "DECREASES": momapy.celldesigner.NegativeInfluence,
    "DIRECTLY_DECREASES": momapy.celldesigner.NegativeInfluence,
}

# The node and relation vocabulary of the projection. Kept here rather than
# imported from the notebook, which is where it is *documented*: `2_05` explains
# what the projection is and why, this is the executable copy the analysis runs on.
MOLECULAR_ENTITY_NODE_TYPES = frozenset(
    {
        "Protein",
        "Complex",
        "Abundance",
        "Composite",
        "Gene",
        "Rna",
        "MicroRna",
        "Activity",
    }
)
PHENOTYPE_NODE_TYPES = frozenset({"BiologicalProcess", "Pathology"})
# What `load_bel_projection` validates a projected node's class against.
PROJECTABLE_NODE_CLASSES = MOLECULAR_ENTITY_NODE_TYPES | PHENOTYPE_NODE_TYPES
# `Population` is an entity for the structural-constituent test below, but is not
# projected: it is not a thing an influence graph should carry.
ENTITY_NODE_TYPES = MOLECULAR_ENTITY_NODE_TYPES | PHENOTYPE_NODE_TYPES | {"Population"}
CAUSAL_RELATIONSHIP_TYPES = frozenset(
    {
        "INCREASES",
        "DECREASES",
        "DIRECTLY_INCREASES",
        "DIRECTLY_DECREASES",
        "REGULATES",
        "TRANSLATED_TO",
    }
)
# Every relation type that is a BEL *statement*, so that a `HAS_*` structural edge
# can be told apart from one by exclusion -- which is how the structural-constituent
# test below works.
ALL_RELATIONSHIP_TYPES = CAUSAL_RELATIONSHIP_TYPES | {
    "ASSOCIATION",
    "POSITIVE_CORRELATION",
    "NEGATIVE_CORRELATION",
    "IS_A",
    "ORTHOLOGOUS",
    "EQUIVALENT_TO",
    "CAUSES_NO_CHANGE",
}

# The projected nodes: the molecular-entity and phenotype nodes, minus the
# structural constituents -- a node another BEL entity points to through a `HAS_*`
# edge (a complex subunit, a variant/modified base form) *and* that carries no
# causal edge of its own. In BEL such a node's causal wiring lives on its
# container. Everything else is kept, including nodes that end up isolated.
_PROJECTION_NODES_QUERY = """
MATCH (collection:Collection)-[:HAS_ENTRY]->(:CollectionEntry)
    -[:HAS_OBJ]->(model:BELModel)-[:HAS_NODE]->(node)
WHERE collection.name IN $collection_names
    AND any(label IN labels(node)
        WHERE label IN $molecular_types OR label IN $phenotype_types)
    AND NOT (
        EXISTS {
            (model)-[:HAS_NODE]->(container)-[part]->(node)
            WHERE NOT type(part) IN $all_relationship_types
                AND any(label IN labels(container) WHERE label IN $entity_types)
        }
        AND NOT EXISTS {
            (node)-[causal]-(other)
            WHERE type(causal) IN $causal_types
                AND (model)-[:HAS_NODE]->(other)
                AND any(label IN labels(other)
                    WHERE label IN $molecular_types OR label IN $phenotype_types)
        }
    )
RETURN elementId(node) AS node_id,
       node.bel AS bel,
       [label IN labels(node)
        WHERE label IN $molecular_types OR label IN $phenotype_types][0] AS node_class
"""

# The causal relations induced between projected nodes. Both endpoints are
# filtered again in Python: a relation whose endpoint was dropped as a structural
# constituent must not become an edge.
_PROJECTION_EDGES_QUERY = """
MATCH (collection:Collection)-[:HAS_ENTRY]->(:CollectionEntry)
    -[:HAS_OBJ]->(model:BELModel)
MATCH (model)-[:HAS_NODE]->(source)-[relationship]->(target)<-[:HAS_NODE]-(model)
WHERE type(relationship) IN $causal_types
    AND collection.name IN $collection_names
    AND any(label IN labels(source)
        WHERE label IN $molecular_types OR label IN $phenotype_types)
    AND any(label IN labels(target)
        WHERE label IN $molecular_types OR label IN $phenotype_types)
RETURN DISTINCT elementId(source) AS source_node_id,
       elementId(target) AS target_node_id,
       type(relationship) AS relation_type
"""

# The activity forms of a protein: `act(p(X))`, `act(p(X),ma(kin))`. Used to widen
# an interface protein's seeds -- see `load_activity_seed_expansion`.
_ACTIVITY_FORMS_QUERY = """
MATCH (:Collection)-[:HAS_ENTRY]->(:CollectionEntry)-[:HAS_OBJ]->(model:BELModel)
    -[:HAS_NODE]->(activity:Activity)
MATCH (activity)-[:HAS__PROTEIN]->(protein)
WHERE elementId(protein) IN $protein_node_ids
RETURN elementId(protein) AS protein_node_id,
       collect(DISTINCT elementId(activity)) AS activity_node_ids
"""


def load_bel_projection(session, collection_names):
    """The influence-graph projection of the given BEL collections (~0.4 s).

    Returns `(nodes, edges)`, where `nodes` is `{node_id: (bel, node_class)}` and
    `edges` is `{(source_node_id, target_node_id): [relation_type, ...]}`. Node ids
    are Neo4j `elementId`s, the same keys `commute_dm.submaps.Influences` and
    momapy_kb's `node_id_to_object` cache use.

    Raises on a node class the projection does not know, rather than dropping it:
    a new BEL class must be admitted deliberately. Silently skipping is how the
    old influence-graph code rotted into a smaller graph than it claimed.
    """
    parameters = {
        "collection_names": list(collection_names),
        "molecular_types": list(MOLECULAR_ENTITY_NODE_TYPES),
        "phenotype_types": list(PHENOTYPE_NODE_TYPES),
        "entity_types": list(ENTITY_NODE_TYPES),
        "causal_types": list(CAUSAL_RELATIONSHIP_TYPES),
        "all_relationship_types": list(ALL_RELATIONSHIP_TYPES),
    }
    nodes = {}
    for row in session.execute_query(_PROJECTION_NODES_QUERY, parameters):
        if row["node_class"] not in PROJECTABLE_NODE_CLASSES:
            raise RuntimeError(
                f"BEL node {row['node_id']} ({row['bel']!r}) has class "
                f"{row['node_class']!r}, which the projection does not admit. Add it "
                "to MOLECULAR_ENTITY_NODE_TYPES or to PHENOTYPE_NODE_TYPES."
            )
        nodes[row["node_id"]] = (row["bel"], row["node_class"])
    edges = collections.defaultdict(list)
    for row in session.execute_query(_PROJECTION_EDGES_QUERY, parameters):
        source_node_id = row["source_node_id"]
        target_node_id = row["target_node_id"]
        if source_node_id not in nodes or target_node_id not in nodes:
            continue
        # A self-loop is drawable but says nothing, and pd2af's arc geometry needs
        # two distinct endpoints. The AD KG has none; this is a guard, not a filter.
        if source_node_id == target_node_id:
            continue
        edges[(source_node_id, target_node_id)].append(row["relation_type"])
    return nodes, dict(edges)


def make_bel_map(nodes, terms, edges, existing_templates=()):
    """`(cd_map, {node_id: species}, stats)` for a projected BEL influence graph.

    The map holds exactly what `commute_dm.submaps.make_submap_from_model_elements`
    reads out of a `source_map`: the model elements, a layout element per model
    element, and the `LayoutModelMapping` tying them together. No arcs -- a
    modulation the source map does not draw gets one from
    `submaps._make_modulation_arc`.

    Two things are many-to-one, both on purpose:

    * `species_by_node_id`, because top-level species are **interned by value**.
      `act(p(X),ma(kin))` and `act(p(X),ma(pep))` describe identically once the
      `ma()` code is dropped, as do whitespace-variant proteoforms, so several
      projected node ids can stand for one species. The walk is unaffected --
      `Influences` is keyed on node ids -- only the drawing merges. Where two
      collapsed ids are reached by different walk directions, the second `fills`
      assignment in `commute_dm.core` wins; harmless while BEL is downstream only.
    * a projected node can *also* be drawn as a subunit of a selected complex, as
      a **different object** (`bel_terms` invariant (c)).

    `existing_templates` is what the BEL templates are interned against, and must
    be the source map's templates *as `make_submap_from_model_elements` re-derives
    them* -- so the source map has to exist first (see
    `commute_dm.core.load_submap_inputs`).

    The glyph positions are throwaway, as everywhere here: the caller runs
    `pd2af.utils.make_auto_layout` on the finished sub-map, which recomputes every
    position and segment.
    """
    root_terms = [terms[node_id] for node_id in sorted(nodes)]
    context = commute_dm.bel_terms.make_build_context(root_terms, existing_templates)
    species_by_node_id = {}
    interned = {}
    for node_id in sorted(nodes):
        species = commute_dm.bel_terms.make_species(terms[node_id], context)
        species_by_node_id[node_id] = interned.setdefault(species, species)
    # Sorted by node id through `species_by_node_id`, so a given knowledge graph
    # always yields the same glyph order and the same `renumber_ids` numbering.
    distinct_species = list(dict.fromkeys(species_by_node_id.values()))

    mapping = momapy.core.mapping.LayoutModelMappingBuilder()
    layout_elements = list(
        commute_dm.bel_terms.make_species_layouts(distinct_species, mapping)
    )
    for species, layout_element in zip(distinct_species, layout_elements):
        mapping.add_mapping(layout_element, species)

    # Only the compartments a *top-level* species sits in: a subunit carries
    # none, and an unused `loc()` would be an empty box. The undrawn per-collection
    # roots come along as `outside`, and get no layout and no mapping entry --
    # exactly like the stored maps' `default` root, so no box is drawn for them.
    drawn_compartments = list(
        dict.fromkeys(
            species.compartment
            for species in distinct_species
            if species.compartment is not None
        )
    )
    compartment_roots = list(
        dict.fromkeys(
            compartment.outside
            for compartment in drawn_compartments
            if compartment.outside is not None
        )
    )
    for compartment in drawn_compartments:
        compartment_layout = commute_dm.bel_terms.make_compartment_layout(compartment)
        layout_elements.append(compartment_layout)
        mapping.add_mapping(compartment_layout, compartment)

    modulations = []
    n_self_loops = 0
    for (source_node_id, target_node_id), relation_types in edges.items():
        source = species_by_node_id[source_node_id]
        target = species_by_node_id[target_node_id]
        # `load_bel_projection` already drops node-id self-loops; interning can
        # create species-level ones the node-id guard cannot see. A self-loop is
        # drawable but says nothing, and pd2af's arc geometry needs two distinct
        # endpoints.
        if source is target:
            n_self_loops += 1
            continue
        for relation_type in relation_types:
            modulation_class = BEL_RELATION_TO_MODULATION_CLASS.get(relation_type)
            if modulation_class is None:  # unsigned, see the class map
                continue
            # One modulation per class, so a pair asserted to both increase and
            # decrease keeps both: BEL sources disagree, and dropping one would
            # pick a winner. `modulations` is a frozenset, so the two relation
            # types that map to the same class collapse on their own.
            modulations.append(modulation_class(source=source, target=target))

    species = frozenset(distinct_species)
    cd_map = momapy.celldesigner.CellDesignerMap(
        model=momapy.celldesigner.CellDesignerModel(
            compartments=frozenset(drawn_compartments + compartment_roots),
            species=species,
            species_templates=frozenset(
                pd2af.celldesigner.building_model.collect_templates_from_species(species)
            ),
            boolean_logic_gates=frozenset(),
            modulations=frozenset(modulations),
        ),
        # This map is only ever a lookup structure, but `Node.position`, `width`
        # and `height` have no defaults, so they have to be given.
        layout=momapy.celldesigner.CellDesignerLayout(
            position=momapy.geometry.Point(0.0, 0.0),
            width=0.0,
            height=0.0,
            layout_elements=tuple(layout_elements),
        ),
        layout_model_mapping=mapping.build(),
    )
    stats = dict(context.stats)
    stats.update(
        {
            "n_projected_node_ids": len(nodes),
            "n_species": len(distinct_species),
            "n_species_collapsed": len(nodes) - len(distinct_species),
            "n_subunits": sum(
                len(list(commute_dm.submaps.iter_subunits(one_species))) for one_species in distinct_species
            ),
            "max_complex_depth": max(
                (_complex_depth(one_species) for one_species in distinct_species),
                default=0,
            ),
            "n_modifications": sum(
                len(getattr(one_species, "modifications", ()) or ())
                for one_species in commute_dm.submaps.iter_species_and_subunits(distinct_species)
            ),
            "n_structural_states": sum(
                len(getattr(one_species, "structural_states", ()) or ())
                for one_species in commute_dm.submaps.iter_species_and_subunits(distinct_species)
            ),
            "n_active": sum(one_species.active for one_species in distinct_species),
            "n_named_by_suffix": sum(
                1 for one_species in distinct_species if " [" in (one_species.name or "")
            ),
            "n_drawn_compartments": len(drawn_compartments),
            "n_self_loop_edges_dropped": n_self_loops,
            "n_modulations": len(cd_map.model.modulations),
            "n_mapping_entries": len(cd_map.layout_model_mapping),
        }
    )
    return cd_map, species_by_node_id, stats


def _complex_depth(species):
    subunits = getattr(species, "subunits", ()) or ()
    return 1 + max((_complex_depth(subunit) for subunit in subunits), default=0)


def make_bel_influences(nodes, edges):
    """The projection as a `commute_dm.submaps.Influences`, on the same node ids.

    The class is reused as it is. A BEL influence graph has no boolean logic gate,
    so `gate_input_node_ids` is empty -- which makes `close_over_gates` a no-op --
    and every projected node is a "species", which makes `species_only` the
    identity. Both are the right behaviour for a gate-less graph, so neither the
    class nor its callers need a BEL case.

    Only the signed relations are walked, matching what :func:`make_bel_map` draws.
    """
    influencing_node_ids = collections.defaultdict(set)
    influenced_node_ids = collections.defaultdict(set)
    for (source_node_id, target_node_id), relation_types in edges.items():
        if not any(
            relation_type in BEL_RELATION_TO_MODULATION_CLASS
            for relation_type in relation_types
        ):
            continue
        influencing_node_ids[target_node_id].add(source_node_id)
        influenced_node_ids[source_node_id].add(target_node_id)
    return commute_dm.submaps.Influences(
        influencing_node_ids=dict(influencing_node_ids),
        influenced_node_ids=dict(influenced_node_ids),
        gate_input_node_ids={},
        species_node_ids=frozenset(nodes),
    )


def merge_influences(*influences):
    """One `Influences` out of several whose node ids are disjoint.

    Disjoint is the point: a CellDesigner collection's node ids and a BEL
    collection's share nothing, so one `Influences` still covers both walk
    directions and `commute_dm.core._select_around_seeds` needs no second object.
    Should a future pair of sources share a node, that stops being true and this is
    where it would have to change.
    """

    def merge(attribute):
        merged = {}
        for one_influences in influences:
            for node_id, neighbour_node_ids in getattr(
                one_influences, attribute
            ).items():
                merged.setdefault(node_id, set()).update(neighbour_node_ids)
        return merged

    return commute_dm.submaps.Influences(
        influencing_node_ids=merge("influencing_node_ids"),
        influenced_node_ids=merge("influenced_node_ids"),
        gate_input_node_ids=merge("gate_input_node_ids"),
        species_node_ids=frozenset().union(
            *(one_influences.species_node_ids for one_influences in influences)
        ),
    )


def merge_with_source_map(source_map, bel_map):
    """`submaps.merge_maps`, plus the checks that say what may and may not fuse.

    momapy elements are equal by value, so anything a BEL element happens to
    equal in the stored maps becomes **one** object on merging. Three different
    answers, and this is where each is asserted:

    * **Species and compartments must not fuse.** A fused species would silently
      join the two influence graphs. Nothing prevents it by construction any
      more -- the names are bare HGNC symbols now -- so the guarantee is
      structural instead: a BEL species either has `compartment=None`, and every
      one of the 6342 stored AF species has a compartment, or it hangs off the
      BEL root, and every stored non-`default` compartment hangs off `default`.
      A fusion means one of those two broke.
    * **Mapping entries must not fuse.** `merge_maps` updates an equality-keyed
      dict, and the BEL map now contributes a key per subunit and per badge, so
      the collision surface is an order of magnitude larger than it was; a
      collision silently destroys a glyph -> model element entry.
    * **Templates *may* fuse** -- interning them against the source map is
      deliberate (`bel_terms` invariant (a)) and is exactly what keeps a BEL
      `MAPT` from emitting a `<proteinReference>` to an id that was never
      written. What must hold is the invariant itself: no two value-equal
      template objects are distinct objects. That is checked directly.

    The checks are complete: a collision inside any sub-map's selection would
    also exist in the union.
    """
    merged = commute_dm.submaps.merge_maps([source_map, bel_map])
    for attribute in ("species", "compartments"):
        n_merged = len(getattr(merged.model, attribute))
        n_parts = len(getattr(source_map.model, attribute)) + len(
            getattr(bel_map.model, attribute)
        )
        if n_merged != n_parts:
            raise RuntimeError(
                f"{n_parts - n_merged} {attribute} fused on merging the BEL map into "
                f"the source map ({n_merged} merged, {n_parts} apart). A BEL element "
                "is value-equal to a stored CellDesigner one, which would join the "
                "two influence graphs. A BEL species has no compartment or hangs off "
                "the BEL root compartment, so one of those guarantees broke."
            )
    n_merged_entries, n_part_entries = commute_dm.submaps.count_mapping_entries(
        [source_map, bel_map]
    )
    if n_merged_entries != n_part_entries:
        raise RuntimeError(
            f"{n_part_entries - n_merged_entries} layout-model mapping entries fused "
            f"on merging the BEL map into the source map ({n_merged_entries} merged, "
            f"{n_part_entries} apart). A BEL glyph, subunit glyph or badge is "
            "value-equal to a stored one, and one glyph -> model element entry is "
            "lost."
        )
    commute_dm.submaps.check_species_template_identity(merged)
    commute_dm.submaps.check_compartment_identity(merged)
    return merged


def load_activity_seed_expansion(session, protein_node_ids, projected_node_ids=None):
    """`{protein_node_id: {activity_node_id, ...}}` for the given BEL proteins.

    The UniProt annotations that define the interface sit on the KG's `:Protein`
    nodes, but BEL keeps a protein's causal wiring on its *activity* form, so
    seeding a walk at the annotated nodes alone reaches little: of the 169 proteins
    the COVID activity-flow maps share with the AD KG, only 86 have an annotated
    node with an outgoing causal edge, against 103 once the activity forms are
    added.

    Complexes and composites containing the protein are deliberately **not**
    included: `act(p(X))` is X in another form, whereas a complex X is a member of
    is a different entity.

    `projected_node_ids` restricts the result to nodes the projection kept, so an
    expansion cannot introduce a seed that has no species.
    """
    protein_node_ids = list(protein_node_ids)
    if not protein_node_ids:
        return {}
    projected = None if projected_node_ids is None else frozenset(projected_node_ids)
    expansion = {}
    for row in session.execute_query(
        _ACTIVITY_FORMS_QUERY, {"protein_node_ids": protein_node_ids}
    ):
        activity_node_ids = set(row["activity_node_ids"])
        if projected is not None:
            activity_node_ids &= projected
        if activity_node_ids:
            expansion[row["protein_node_id"]] = activity_node_ids
    return expansion
