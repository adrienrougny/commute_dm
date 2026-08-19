"""The AD BEL KG as an *exportable* CellDesigner map.

`2_10_make_ad_kg_cd_af_collection` is the only importer. It loads the knowledge
graph's **influence-graph projection** (the one `2_05` defines and documents),
turns it into a `CellDesignerMap` of real activity-flow content -- proteins with
phosphorylation badges, complexes with nested subunit glyphs, `act(...)` as
CellDesigner's active decoration, `loc(...)` as a drawn compartment, the signed
causal relations as modulations -- writes it to one file and imports that file as
an ordinary collection, `AD_KG_CD_AF`.

After that the knowledge graph *is* a CellDesigner collection, so nothing else in
the package knows it was ever BEL: `commute_dm.core` and `commute_dm.submaps` walk
its stored species and modulations exactly as they walk COVID's and PD's. The
in-memory merge this module used to do -- building the BEL side into an already
hydrated `source_map` on every run -- is gone with it.

Two things the export decides, both measured in `2_10`:

* **`act(X)` and `X` are one species.** The projection has no edge between the
  two, so keeping them apart severs 6889 upstream -> downstream two-hop paths that
  three hops cannot recover. The price is that X carries an active border whenever
  BEL asserts any activity of it.
* **Species carrying no signed modulation are dropped** (1228 of 3979 projected
  nodes). Nothing walks or draws them, so they would stand stranded.

**What a term becomes, and the identity invariants that keep the result readable,
live in `commute_dm.bel_terms`.** This module is only the assembly: which nodes
are projected, how they are wired, and which annotations they carry over.
"""

import collections

import momapy.celldesigner
import momapy.core.mapping
import momapy.geometry
import momapy.sbml.model
import pd2af.celldesigner.building_model

import commute_dm.bel_terms
import commute_dm.queries
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


def make_bel_map(nodes, terms, edges, drop_isolated_species=True):
    """`(cd_map, {node_id: species}, [(node_id, species)], stats)` for a projection.

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

    The third return value is the `(node_id, species)` pairs of **every** species
    built, subunits included, which is what
    :func:`make_element_to_annotations` needs: a subunit is a per-occurrence
    object no `{node_id: species}` map can recover, and it carries interface
    identifiers of its own.

    `drop_isolated_species` removes the species that carry no signed modulation
    (1228 of the AD KG's 3979 projected nodes). Nothing walks or draws them, so
    they would stand stranded in the collection. The drop happens **after** the
    species are interned and the modulations built, so a protein isolated in its
    own right but wired through its activity form survives -- `act(p(X))` and
    `p(X)` are the same species (`bel_terms.collect_activity_subject_term_ids`).
    A dropped species that is a member of a kept complex still appears as that
    complex's subunit; only its own top-level glyph goes.

    The glyph positions are throwaway, as everywhere here: the caller runs
    `pd2af.utils.make_auto_layout` on the finished sub-map, which recomputes every
    position and segment.
    """
    root_terms = [terms[node_id] for node_id in sorted(nodes)]
    # Over **every loaded term**, not just the projected roots: a protein is
    # drawn active whenever BEL asserts any activity of it, and the activity term
    # itself need not be projected.
    active_term_ids = commute_dm.bel_terms.collect_activity_subject_term_ids(
        terms.values()
    )
    context = commute_dm.bel_terms.make_build_context(root_terms, active_term_ids)
    species_by_node_id = {}
    interned = {}
    node_id_species_pairs = []
    for node_id in sorted(nodes):
        species = commute_dm.bel_terms.make_species(
            terms[node_id], context, record=node_id_species_pairs
        )
        species_by_node_id[node_id] = interned.setdefault(species, species)
    # Sorted by node id through `species_by_node_id`, so a given knowledge graph
    # always yields the same glyph order and the same `renumber_ids` numbering.
    distinct_species = list(dict.fromkeys(species_by_node_id.values()))

    # The modulations first, because whether a species is isolated is a property
    # of them.
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

    n_isolated_dropped = 0
    if drop_isolated_species:
        wired = set()
        for modulation in modulations:
            wired.add(id(modulation.source))
            wired.add(id(modulation.target))
        kept = [species for species in distinct_species if id(species) in wired]
        n_isolated_dropped = len(distinct_species) - len(kept)
        distinct_species = kept
        kept_ids = {id(species) for species in distinct_species}
        species_by_node_id = {
            node_id: species
            for node_id, species in species_by_node_id.items()
            if id(species) in kept_ids
        }
    # Only the species the map actually holds may contribute annotations, and
    # membership is by **value**: a top-level species and a value-equal subunit
    # share one `element_to_annotations` entry, which is exactly what the writer
    # looks up.
    drawn_values = set(
        commute_dm.submaps.iter_species_and_subunits(distinct_species)
    )
    node_id_species_pairs = [
        (node_id, species)
        for node_id, species in node_id_species_pairs
        if species in drawn_values
    ]

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
            "n_species_collapsed": (
                len(species_by_node_id) - len(distinct_species)
            ),
            "n_isolated_species_dropped": n_isolated_dropped,
            "n_node_ids_kept": len(species_by_node_id),
            "n_node_id_species_pairs": len(node_id_species_pairs),
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
    return cd_map, species_by_node_id, node_id_species_pairs, stats


def _complex_depth(species):
    subunits = getattr(species, "subunits", ()) or ()
    return 1 + max((_complex_depth(subunit) for subunit in subunits), default=0)


# The annotation payloads of a set of nodes, in the encoding `2_00` wrote them
# in: a per-entry `Mapping` of `Item` -> `Bag` of single-resource
# `RDFAnnotation`s, each qualified by a shared `BQBiol` node. Read back here so
# the exported file carries the same cross-references the BEL nodes carry, which
# is what puts its proteins in the interface.
_ANNOTATION_PAYLOADS_QUERY = """
UNWIND $element_ids AS element_id
MATCH (node) WHERE elementId(node) = element_id
MATCH (node)<-[:HAS_KEY]-(:Item)-[:HAS_VALUE]->(:Bag)-[:HAS_ITEM]->(a:RDFAnnotation)
MATCH (a)-[:HAS_QUALIFIER]->(q)
RETURN elementId(node) AS node_id,
       q.name AS qualifier_name,
       a.resources AS resources
"""


def make_element_to_annotations(session, node_id_species_pairs):
    """`{species: frozenset[RDFAnnotation]}` for the exported BEL species.

    Two queries. `commute_dm.queries.get_annotated_nodes` says which nodes'
    annotations stand for a given node -- an activity's subject is followed, a
    complex's members are **not** (`with_subunits=False`), because each member is
    drawn as its own subunit species and annotated in its own right. Then one
    query for the payloads, rebuilt as `momapy.sbml.model.RDFAnnotation`s: a
    faithful copy of what `2_00` wrote.

    Accumulated **by value**. A top-level species and a value-equal subunit are
    two objects but one key, which is what the CellDesigner writer's
    `element_to_annotations.get(species)` looks up -- both for a `<species>` and
    for an included species, whose RDF goes inside `<celldesigner:notes>`.
    """
    node_ids = sorted({node_id for node_id, _ in node_id_species_pairs})
    if not node_ids:
        return {}
    nodes = commute_dm.queries.get_nodes(session, node_ids)
    annotated_node_ids_by_node_id = {
        node.element_id: [part.element_id for part in parts]
        for node, parts in commute_dm.queries.get_annotated_nodes(
            session, nodes, with_subunits=False
        )
    }
    payload_element_ids = sorted(
        {
            annotated_node_id
            for annotated_node_ids in annotated_node_ids_by_node_id.values()
            for annotated_node_id in annotated_node_ids
        }
    )
    annotations_by_element_id = collections.defaultdict(set)
    for row in session.execute_query(
        _ANNOTATION_PAYLOADS_QUERY, params={"element_ids": payload_element_ids}
    ):
        annotations_by_element_id[row["node_id"]].add(
            momapy.sbml.model.RDFAnnotation(
                qualifier=momapy.sbml.model.BQBiol[row["qualifier_name"]],
                resources=frozenset(row["resources"]),
            )
        )

    element_to_annotations = collections.defaultdict(set)
    for node_id, species in node_id_species_pairs:
        for annotated_node_id in annotated_node_ids_by_node_id.get(node_id, ()):
            element_to_annotations[species] |= annotations_by_element_id.get(
                annotated_node_id, set()
            )
    return {
        species: frozenset(annotations)
        for species, annotations in element_to_annotations.items()
        if annotations
    }
