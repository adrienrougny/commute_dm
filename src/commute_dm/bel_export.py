"""The AD BEL KG as an *exportable* CellDesigner map.

Loads the knowledge graph's **influence-graph projection** as a
`momapy_bel.core.BELModel`, then turns that model into a `CellDesignerMap` of
real activity-flow content -- proteins with phosphorylation badges, complexes
with nested subunit glyphs, `act(...)` as CellDesigner's active decoration,
`loc(...)` as a drawn compartment, the signed causal relations as modulations --
and hands back a map that can be written to a file and imported as an ordinary
collection.

Once imported, the knowledge graph *is* a CellDesigner collection, so nothing
else in the package knows it was ever BEL: `commute_dm.core` and
`commute_dm.submaps` walk its stored species and modulations exactly as they walk
COVID's and PD's.

Two things the export decides:

* **`act(X)` and `X` are one species.** The projection has no relation between
  the two, so keeping them apart severs 6889 upstream -> downstream two-hop paths
  that three hops cannot recover. The price is that X carries an active border
  whenever BEL asserts any activity of it.
* **Species carrying no signed modulation are dropped** (1228 of 3979 projected
  elements). Nothing walks or draws them, so they would stand stranded.

**What an element becomes, and the identity invariants that keep the result
readable, live in `commute_dm.bel2cd`.** This module is only the loading and the
assembly: which nodes are projected, how they are wired, and which annotations
they carry over.
"""

import collections

import momapy.celldesigner
import momapy.core.mapping
import momapy.geometry
import momapy.sbml.model
import momapy_bel.core
import momapy_bel.io.bel
import pd2af.celldesigner.building_model

import commute_dm.bel2cd
import commute_dm.bel_vocabulary
import commute_dm.submaps


# The BEL node classes the projection keeps -- the validation set for
# :func:`get_bel_influence_graph_projection`. What each is *drawn* as is
# `commute_dm.bel2cd.BEL_CLASS_TO_CD_CLASSES`, which is a different question: an
# `Activity` is projected but has no class of its own, taking its subject's.
#
# The interface is joined on the annotations of the BEL collection's own
# `:Protein` nodes, so nothing here decides which proteins are in the interface
# or which nodes can seed a walk.
PROJECTABLE_NODE_CLASSES = (
    commute_dm.bel_vocabulary.MOLECULAR_ENTITY_NODE_TYPES
    | commute_dm.bel_vocabulary.PHENOTYPE_NODE_TYPES
)

# The BEL entity node classes, as the momapy_bel elements they are loaded as.
# Nothing else is ever built: a projected node is one of these, and so is every
# member of one -- BEL *process* terms (reactions, translocations) are neither.
BEL_NODE_CLASS_TO_ELEMENT_CLASS = {
    "Protein": momapy_bel.core.ProteinAbundance,
    "Gene": momapy_bel.core.GeneAbundance,
    "Rna": momapy_bel.core.RNAAbundance,
    "MicroRna": momapy_bel.core.MicroRNAAbundance,
    "Abundance": momapy_bel.core.Abundance,
    "Complex": momapy_bel.core.ComplexAbundance,
    "Composite": momapy_bel.core.CompositeAbundance,
    "Activity": momapy_bel.core.Activity,
    "BiologicalProcess": momapy_bel.core.BiologicalProcess,
    "Pathology": momapy_bel.core.Pathology,
}

# The BEL causal relation types, as the momapy_bel statements they are loaded as.
BEL_RELATION_TYPE_TO_RELATION_CLASS = {
    "INCREASES": momapy_bel.core.Increases,
    "DIRECTLY_INCREASES": momapy_bel.core.DirectlyIncreases,
    "TRANSLATED_TO": momapy_bel.core.TranslatedTo,
    "DECREASES": momapy_bel.core.Decreases,
    "DIRECTLY_DECREASES": momapy_bel.core.DirectlyDecreases,
    "REGULATES": momapy_bel.core.Regulates,
}

# The signed relations. `Regulates` is deliberately absent: it carries no sign,
# and `commute_dm.submaps.SIGNED_MODULATION_CLASSES` -- the filter for the
# modulations induced from `source_map` -- excludes the unsigned `Modulation`, so
# mapping it would give an edge that is *walked but never drawn*, and sub-maps
# with invisible connections. It stays in the projection, which is what `2_05`'s
# statistics run on, and is dropped here, at the cost of 2.6% of the AD KG's
# causal edges.
BEL_RELATION_CLASS_TO_MODULATION_CLASS = {
    momapy_bel.core.Increases: momapy.celldesigner.PositiveInfluence,
    momapy_bel.core.DirectlyIncreases: momapy.celldesigner.PositiveInfluence,
    momapy_bel.core.TranslatedTo: momapy.celldesigner.PositiveInfluence,
    momapy_bel.core.Decreases: momapy.celldesigner.NegativeInfluence,
    momapy_bel.core.DirectlyDecreases: momapy.celldesigner.NegativeInfluence,
}

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
RETURN elementId(node) AS node_id
"""

# The causal relations induced between projected nodes. Both endpoints are
# filtered again in Python: a relation whose endpoint was dropped as a structural
# constituent must not become a statement.
_PROJECTION_RELATIONSHIPS_QUERY = """
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

# Every `:BELModel` node of the collections, with the properties an element is
# built from. `default` (the activity code), `type` (the pmod type), `position`
# and `range` are Cypher reserved-ish words and must be backticked. `DISTINCT`
# because two collection entries can reach the same model.
_NODES_QUERY = """
MATCH (collection:Collection)-[:HAS_ENTRY]->(:CollectionEntry)
    -[:HAS_OBJ]->(model:BELModel)-[:HAS_NODE]->(node)
WHERE collection.name IN $collection_names
RETURN DISTINCT elementId(node) AS node_id,
       labels(node) AS node_labels,
       node.bel AS bel,
       node.name AS name,
       node.namespace AS namespace,
       node.`default` AS activity_code,
       node.hgvs AS hgvs,
       node.`type` AS modification_type,
       node.amino_acid AS amino_acid,
       toString(node.`position`) AS residue,
       node.`range` AS range,
       node.descriptor AS descriptor
"""

# The structural (sub-term) edges. The **double** underscore is deliberate: it
# excludes the sparse single-underscore statement route (`HAS_MEMBERS`,
# `HAS_COMPONENTS`, `HAS_MODIFIED_GENE`, ...) and the `HAS_NODE` / `HAS_OBJ`
# plumbing.
_EDGES_QUERY = """
MATCH (collection:Collection)-[:HAS_ENTRY]->(:CollectionEntry)
    -[:HAS_OBJ]->(model:BELModel)-[:HAS_NODE]->(container)
MATCH (container)-[part]->(node)
WHERE collection.name IN $collection_names
    AND type(part) STARTS WITH 'HAS__'
RETURN DISTINCT elementId(container) AS container_node_id,
       elementId(node) AS node_id,
       type(part) AS edge_type
"""


def get_bel_string(element):
    """The BEL string of one momapy_bel element -- its label wherever one is due.

    Canonical: a complex's members are sorted, so the string is the same from one
    run to the next, which is what makes it a usable sort key.
    """
    return momapy_bel.io.bel.BELWriter._bel_element_to_string(element)


def get_bel_influence_graph_projection(session, collection_names):
    """The influence-graph projection of the given BEL collections, as a `BELModel`.

    Its `statements` are the induced causal relations -- `Increases`,
    `DirectlyIncreases`, `TranslatedTo`, `Decreases`, `DirectlyDecreases`,
    `Regulates` -- between projected elements, plus every projected element that
    carries none of them, as a bare abundance: an isolated node is part of the
    projection and `2_05` counts it.

    The projected nodes are the molecular-entity and phenotype nodes minus the
    structural constituents -- a node another BEL entity points at through a
    `HAS__*` edge *and* that carries no causal relationship of its own, whose
    wiring lives on its container. Everything else is kept, isolated nodes
    included.

    Elements are **value-equal**, so two nodes that describe the same BEL term are
    one element. Measured over the four knowledge graphs, only PD's three
    `gmod(TestNS:"TestName")` genes merge that way -- `momapy_bel` has no place
    for a gene modification, and they are a test artefact.

    Raises on a node class or a `HAS__*` edge type this module does not know,
    rather than dropping it: silently skipping is how the old influence-graph code
    rotted into a smaller graph than it claimed.
    """
    parameters = {
        "collection_names": list(collection_names),
        "molecular_types": list(commute_dm.bel_vocabulary.MOLECULAR_ENTITY_NODE_TYPES),
        "phenotype_types": list(commute_dm.bel_vocabulary.PHENOTYPE_NODE_TYPES),
        "entity_types": list(commute_dm.bel_vocabulary.ENTITY_NODE_TYPES),
        "causal_types": list(commute_dm.bel_vocabulary.INFLUENCE_RELATIONSHIP_TYPES),
        "all_relationship_types": list(
            commute_dm.bel_vocabulary.ALL_RELATIONSHIP_TYPES
        ),
    }
    rows = {
        row["node_id"]: row
        for row in session.execute_query(_NODES_QUERY, parameters)
    }
    sub_node_ids = collections.defaultdict(lambda: collections.defaultdict(list))
    for row in session.execute_query(_EDGES_QUERY, parameters):
        edge_type = row["edge_type"]
        if edge_type not in commute_dm.bel_vocabulary.STRUCTURAL_EDGE_TYPES:
            raise RuntimeError(
                f"BEL node {row['container_node_id']} has a {edge_type} sub-term edge, "
                "which is not classified. Add it to MEMBER_EDGE_TYPES, to "
                "MODIFIER_EDGE_TYPES or to OTHER_STRUCTURAL_EDGE_TYPES in "
                "commute_dm.bel_vocabulary."
            )
        sub_node_ids[row["container_node_id"]][edge_type].append(row["node_id"])
    # Sorted so that a given knowledge graph always yields the same subunit order,
    # the same synthetic positions and the same `renumber_ids` numbering.
    for edge_type_to_node_ids in sub_node_ids.values():
        for node_ids in edge_type_to_node_ids.values():
            node_ids.sort(key=lambda node_id: (rows[node_id]["bel"] or "", node_id))

    elements = {}
    node_id_to_element = {}
    for row in session.execute_query(_PROJECTION_NODES_QUERY, parameters):
        node_id = row["node_id"]
        node_class = _resolve_node_class(node_id, rows[node_id]["node_labels"])
        if node_class not in PROJECTABLE_NODE_CLASSES:
            raise RuntimeError(
                f"BEL node {node_id} ({rows[node_id]['bel']!r}) has class "
                f"{node_class!r}, which the projection does not admit. Add it to "
                "MOLECULAR_ENTITY_NODE_TYPES or to PHENOTYPE_NODE_TYPES in "
                "commute_dm.bel_vocabulary."
            )
        node_id_to_element[node_id] = _make_element(
            node_id, rows, sub_node_ids, elements
        )

    statements = []
    wired = set()
    for row in session.execute_query(_PROJECTION_RELATIONSHIPS_QUERY, parameters):
        source = node_id_to_element.get(row["source_node_id"])
        target = node_id_to_element.get(row["target_node_id"])
        # A relation whose endpoint was dropped as a structural constituent must
        # not become a statement.
        if source is None or target is None:
            continue
        relation_class = BEL_RELATION_TYPE_TO_RELATION_CLASS[row["relation_type"]]
        statements.append(relation_class(source=source, target=target))
        wired.update((source, target))
    statements.extend(
        element for element in node_id_to_element.values() if element not in wired
    )
    return momapy_bel.core.BELModel(statements=frozenset(statements))


def _resolve_node_class(node_id, node_labels):
    """The concrete BEL class of a node.

    pylpg labels a node with its own class *and* every ancestor (`['GeneticFlow',
    'Protein', 'BELModelElement']`), so the class is resolved by intersecting the
    labels with the known ones rather than by taking `labels(node)[0]`, which is
    order-dependent. Exactly one must match, so an unexpected label combination
    fails loudly.
    """
    matched = [
        label for label in node_labels if label in BEL_NODE_CLASS_TO_ELEMENT_CLASS
    ]
    if len(matched) != 1:
        raise RuntimeError(
            f"BEL node {node_id} has labels {sorted(node_labels)}, matching "
            f"{len(matched)} known entity classes {matched}. Exactly one is "
            "expected; add the class to BEL_NODE_CLASS_TO_ELEMENT_CLASS."
        )
    return matched[0]


def _make_element(node_id, rows, sub_node_ids, elements):
    """The momapy_bel element of one BEL node, memoised, built bottom up.

    The node id rides along as the element's `id_`, which is excluded from
    equality: it is what lets `2_05` ask the database about a projected element
    without a second index. Two nodes that describe the same term become one
    element carrying one of their two ids.
    """
    element = elements.get(node_id)
    if element is not None:
        return element
    element = _build_element(node_id, rows, sub_node_ids, elements)
    elements[node_id] = element
    return element


def _build_element(node_id, rows, sub_node_ids, elements):
    row = rows[node_id]
    element_class = BEL_NODE_CLASS_TO_ELEMENT_CLASS[
        _resolve_node_class(node_id, row["node_labels"])
    ]
    parts = sub_node_ids[node_id]
    fields = {"id_": node_id}
    if element_class is momapy_bel.core.Activity:
        members = _member_node_ids(parts)
        if len(members) != 1:
            raise RuntimeError(
                f"BEL activity {node_id} ({row['bel']!r}) has {len(members)} subject "
                "sub-terms; exactly one is expected."
            )
        fields["abundance"] = _make_element(members[0], rows, sub_node_ids, elements)
        fields["molecular_activity"] = _make_molecular_activity(row)
        return element_class(**fields)
    fields["namespace"] = row["namespace"]
    fields["identifier"] = row["name"]
    if issubclass(element_class, momapy_bel.core.Abundance):
        fields["location"] = _make_location(parts, rows)
    if element_class in (
        momapy_bel.core.ComplexAbundance,
        momapy_bel.core.CompositeAbundance,
    ):
        fields["members"] = frozenset(
            _make_element(member_node_id, rows, sub_node_ids, elements)
            for member_node_id in _member_node_ids(parts)
        )
    if element_class in (
        momapy_bel.core.ProteinAbundance,
        momapy_bel.core.GeneAbundance,
        momapy_bel.core.RNAAbundance,
        momapy_bel.core.MicroRNAAbundance,
    ):
        # `gmod()` has no place in a `GeneAbundance` and is dropped: the only ones
        # in the data are PD's three `gmod(TestNS:"TestName")` test artefacts.
        fields["fusion"] = None
        fields["variants"] = tuple(
            momapy_bel.core.Variant(id_=sub_node_id, descriptor=rows[sub_node_id]["hgvs"])
            for sub_node_id in parts["HAS__VARIANT"]
        )
    if element_class is momapy_bel.core.ProteinAbundance:
        fields["fragment"] = _make_fragment(parts, rows)
        fields["modifications"] = tuple(
            momapy_bel.core.ProteinModification(
                id_=sub_node_id,
                namespace=rows[sub_node_id]["namespace"] or "",
                identifier=rows[sub_node_id]["name"]
                or rows[sub_node_id]["modification_type"],
                amino_acid=rows[sub_node_id]["amino_acid"],
                residue=rows[sub_node_id]["residue"],
            )
            for sub_node_id in parts["HAS__PMOD"]
        )
    return element_class(**fields)


def _member_node_ids(parts):
    return [
        node_id
        for edge_type in sorted(commute_dm.bel_vocabulary.MEMBER_EDGE_TYPES)
        for node_id in parts[edge_type]
    ]


def _make_molecular_activity(row):
    if row["name"]:
        return momapy_bel.core.MolecularActivity(
            namespace=row["namespace"] or "", identifier=row["name"]
        )
    if row["activity_code"]:
        return momapy_bel.core.MolecularActivity(
            namespace="", identifier=row["activity_code"]
        )
    return None


def _make_location(parts, rows):
    for node_id in parts["HAS__LOCATION"]:
        return momapy_bel.core.Location(
            id_=node_id,
            namespace=rows[node_id]["namespace"],
            identifier=rows[node_id]["name"],
        )
    return None


def _make_fragment(parts, rows):
    for node_id in parts["HAS__FRAGMENT"]:
        return momapy_bel.core.Fragment(
            id_=node_id,
            start_stop=rows[node_id]["range"],
            descriptor=rows[node_id]["descriptor"],
        )
    return None


def make_cd_map_from_bel_influence_graph_projection(
    session, collection_names, projection, drop_isolated_species=True
):
    """`(cd_map, element_to_annotations)` for an influence-graph projection.

    The session is needed for one thing only: reading the projected elements'
    cross-references back onto the species built for them, which is what puts this
    collection's proteins in the interface.

    The map holds exactly what `commute_dm.submaps.make_submap_from_model_elements`
    reads out of a `source_map`: the model elements, a layout element per model
    element, and the `LayoutModelMapping` tying them together. No arcs -- a
    modulation the source map does not draw gets one from
    `submaps._make_modulation_arc`.

    Top-level species are **interned by value**: `act(p(X),ma(kin))` and
    `act(p(X),ma(pep))` describe identically once the `ma()` code is dropped, as
    do `act(p(X))` and `p(X)`, so several elements can stand for one species.

    `drop_isolated_species` removes the species that carry no signed modulation
    (1228 of the AD KG's 3979 projected elements). Nothing walks or draws them, so
    they would stand stranded in the collection. The drop happens **after** the
    species are interned and the modulations built, so a protein isolated in its
    own right but wired through its activity form survives. A dropped species that
    is a member of a kept complex still appears as that complex's subunit; only
    its own top-level glyph goes.

    The glyph positions are throwaway, as everywhere here: the caller runs
    `pd2af.utils.make_auto_layout` on the finished sub-map, which recomputes every
    position and segment.
    """
    relation_classes = tuple(BEL_RELATION_TYPE_TO_RELATION_CLASS.values())
    relations = [
        statement
        for statement in projection.statements
        if isinstance(statement, relation_classes)
    ]
    # Sorted by BEL string, so a given knowledge graph always yields the same glyph
    # order and the same `renumber_ids` numbering.
    elements = sorted(
        {element for relation in relations for element in (relation.source, relation.target)}
        | {
            statement
            for statement in projection.statements
            if not isinstance(statement, relation_classes)
        },
        key=get_bel_string,
    )
    species_fields = commute_dm.bel2cd.make_species_fields(elements)
    templates = commute_dm.bel2cd.make_templates(species_fields)
    compartments, compartment_root = commute_dm.bel2cd.make_compartments(
        species_fields, ", ".join(collection_names)
    )
    element_to_species = {}
    interned = {}
    species_records = []
    for element in elements:
        species = commute_dm.bel2cd.make_species(
            element, species_fields, templates, compartments, record=species_records
        )
        element_to_species[element] = interned.setdefault(species, species)
    distinct_species = list(dict.fromkeys(element_to_species.values()))

    # The modulations first, because whether a species is isolated is a property
    # of them.
    modulations = []
    for relation in relations:
        modulation_class = BEL_RELATION_CLASS_TO_MODULATION_CLASS.get(type(relation))
        if modulation_class is None:  # unsigned, see the class map
            continue
        source = element_to_species[relation.source]
        target = element_to_species[relation.target]
        # A self-loop is drawable but says nothing, and pd2af's arc geometry needs
        # two distinct endpoints. Interning creates species-level ones (an
        # `act(X) -> X` relation, say) that no element-level filter could see.
        if source is target:
            continue
        # One modulation per class, so a pair asserted to both increase and
        # decrease keeps both: BEL sources disagree, and dropping one would pick a
        # winner. `modulations` is a frozenset, so two relation classes that map to
        # the same modulation class collapse on their own.
        modulations.append(modulation_class(source=source, target=target))

    if drop_isolated_species:
        wired = set()
        for modulation in modulations:
            wired.add(id(modulation.source))
            wired.add(id(modulation.target))
        distinct_species = [
            species for species in distinct_species if id(species) in wired
        ]
    # Only the species the map actually holds may contribute annotations, and
    # membership is by **value**: a top-level species and a value-equal subunit
    # share one `element_to_annotations` entry, which is exactly what the writer
    # looks up.
    drawn_values = set(
        commute_dm.submaps.iter_species_and_subunits(distinct_species)
    )
    # **Every** species built, subunits included, which is what the annotations
    # need: a subunit is a per-occurrence object, and it carries interface
    # identifiers of its own.
    element_to_species_and_subunits = collections.defaultdict(list)
    for element, species in species_records:
        if species in drawn_values:
            element_to_species_and_subunits[element].append(species)

    mapping = momapy.core.mapping.LayoutModelMappingBuilder()
    layout_elements = list(
        commute_dm.bel2cd.make_species_layouts(distinct_species, mapping)
    )
    for species, layout_element in zip(distinct_species, layout_elements):
        mapping.add_mapping(layout_element, species)

    # Only the compartments a *top-level* species sits in: a subunit carries none,
    # and an unused `loc()` would be an empty box. The undrawn root comes along as
    # `outside`, and gets no layout and no mapping entry -- exactly like the stored
    # maps' `default` root, so no box is drawn for it.
    drawn_compartments = list(
        dict.fromkeys(
            species.compartment
            for species in distinct_species
            if species.compartment is not None
        )
    )
    for compartment in drawn_compartments:
        compartment_layout = commute_dm.bel2cd.make_compartment_layout(compartment)
        layout_elements.append(compartment_layout)
        mapping.add_mapping(compartment_layout, compartment)

    species = frozenset(distinct_species)
    cd_map = momapy.celldesigner.CellDesignerMap(
        model=momapy.celldesigner.CellDesignerModel(
            compartments=frozenset(drawn_compartments + [compartment_root]),
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
    element_to_annotations = _make_element_to_annotations(
        session, collection_names, element_to_species_and_subunits
    )
    return cd_map, element_to_annotations


# The annotation payloads of a BEL collection's nodes, in the encoding `2_00`
# wrote them in: a per-entry `Mapping` of `Item` -> `Bag` of single-resource
# `RDFAnnotation`s, each qualified by a shared `BQBiol` node. Read back here so
# the exported file carries the same cross-references the BEL nodes carry, which
# is what puts its proteins in the interface. Keyed by what a momapy_bel element
# knows about itself -- its class, namespace and identifier -- since that is what
# `2_00` keyed the annotations on.
_ANNOTATION_PAYLOADS_QUERY = """
MATCH (collection:Collection)-[:HAS_ENTRY]->(:CollectionEntry)
    -[:HAS_OBJ]->(model:BELModel)-[:HAS_NODE]->(node)
WHERE collection.name IN $collection_names
MATCH (node)<-[:HAS_KEY]-(:Item)-[:HAS_VALUE]->(:Bag)-[:HAS_ITEM]->(a:RDFAnnotation)
MATCH (a)-[:HAS_QUALIFIER]->(q)
RETURN DISTINCT elementId(node) AS node_id,
       labels(node) AS node_labels,
       node.namespace AS namespace,
       node.name AS name,
       q.name AS qualifier_name,
       a.resources AS resources
"""


def _make_element_to_annotations(
    session, collection_names, element_to_species_and_subunits
):
    """`{species: frozenset[RDFAnnotation]}` for the exported BEL species.

    One query, rebuilding the payloads as `momapy.sbml.model.RDFAnnotation`s: a
    faithful copy of what `2_00` wrote. An **activity's subject is followed** --
    `act(p(X))` is X in another form -- while a **complex's members are not**,
    because each is drawn as its own subunit species and annotated in its own
    right.

    Accumulated **by value**. A top-level species and a value-equal subunit are
    two objects but one key, which is what the CellDesigner writer's
    `element_to_annotations.get(species)` looks up -- both for a `<species>` and
    for an included species, whose RDF goes inside `<celldesigner:notes>`.
    """
    annotations_by_key = collections.defaultdict(set)
    for row in session.execute_query(
        _ANNOTATION_PAYLOADS_QUERY,
        params={"collection_names": list(collection_names)},
    ):
        key = (
            BEL_NODE_CLASS_TO_ELEMENT_CLASS[
                _resolve_node_class(row["node_id"], row["node_labels"])
            ],
            row["namespace"],
            row["name"],
        )
        annotations_by_key[key].add(
            momapy.sbml.model.RDFAnnotation(
                qualifier=momapy.sbml.model.BQBiol[row["qualifier_name"]],
                resources=frozenset(row["resources"]),
            )
        )
    element_to_annotations = collections.defaultdict(set)
    for element, species_list in element_to_species_and_subunits.items():
        while isinstance(element, momapy_bel.core.Activity):
            element = element.abundance
        annotations = annotations_by_key.get(
            (type(element), element.namespace, element.identifier)
        )
        if annotations is None:
            continue
        for species in species_list:
            element_to_annotations[species] |= annotations
    return {
        species: frozenset(annotations)
        for species, annotations in element_to_annotations.items()
    }
