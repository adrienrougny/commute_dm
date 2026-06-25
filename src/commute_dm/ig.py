import copy
import dataclasses
import math

import pydot
import momapy.geometry
import momapy.core
import momapy.core.mapping
import momapy.coloring
import momapy.positioning
import momapy.celldesigner
import momapy.builder


REACTANT_TO_PRODUCT = "_REACTANT_TO_PRODUCT"
POSITIVE_INFLUENCE = "_POSITIVE_INFLUENCE"
NEGATIVE_INFLUENCE = "_NEGATIVE_INFLUENCE"
NECESSARY_POSITIVE_INFLUENCE = "_NECESSARY_POSITIVE_INFLUENCE"
UNCERTAIN_POSITIVE_INFLUENCE = "_UNCERTAIN_POSITIVE_INFLUENCE"
UNCERTAIN_NECESSARY_POSITIVE_INFLUENCE = "_UNCERTAIN_NECESSARY_POSITIVE_INFLUENCE"
UNCERTAIN_NEGATIVE_INFLUENCE = "_UNCERTAIN_NEGATIVE_INFLUENCE"

INFLUENCES = [
    REACTANT_TO_PRODUCT,
    POSITIVE_INFLUENCE,
    NEGATIVE_INFLUENCE,
    NECESSARY_POSITIVE_INFLUENCE,
]
POINTS_PER_INCH = 72


class InfluenceGraph(dict):
    def get_nodes(self):
        return set(self.keys())

    def add_node(self, node):
        self[node] = set()

    def add_relationship(self, relationship):
        self[relationship.end_node].add(relationship)

    def get_relationships(self):
        relationships = set([])
        for node in self.get_nodes():
            relationships.update(self[node])
        return relationships

    def get_stimulators(self, node):
        return [
            relationship.start_node
            for relationship in self[node]
            if relationship.type == POSITIVE_INFLUENCE
        ]

    def get_inhibitors(self, node):
        return [
            relationship.start_node
            for relationship in self[node]
            if relationship.type == NEGATIVE_INFLUENCE
        ]

    def get_necessary_stimulators(self, node):
        return [
            relationship.start_node
            for relationship in self[node]
            if relationship.type == NECESSARY_POSITIVE_INFLUENCE
            or relationship.type == REACTANT_TO_PRODUCT
        ]

    def get_modulators(self, node):
        return (
            self.get_stimulators(node)
            + self.get_inhibitors(node)
            + self.get_necessary_stimulators(node)
        )

    def remove_node(self, node):
        del self[node]
        for other_node in self.get_nodes():
            for relationship in list(self[other_node]):
                if relationship.start_node == node:
                    self[other_node].remove(relationship)


def make_ig_in_db(session):
    """Materialize the influence graph as direct relationships between
    Species nodes in the DB. For each Reaction/Modulation, derive
    REACTANT_TO_PRODUCT / POSITIVE_INFLUENCE / NEGATIVE_INFLUENCE / etc.
    edges between the participating species. Each derived relationship is
    annotated with the source `model_ids` and `reaction_ids` (or
    `modulation_ids`) so consumers can trace back to provenance.
    """
    queries = []
    queries.append(
        f"""
        MATCH
            (reaction:Reaction),
            (reaction)-[:HAS_REACTANT]->(reactant:Reactant),
            (reactant)-[:HAS_REFERRED_SPECIES]->(reactant_species:Species),
            (reaction)-[:HAS_PRODUCT]->(product:Product),
            (product)-[:HAS_REFERRED_SPECIES]->(product_species:Species)
        MERGE
            (reactant_species)-[r:{REACTANT_TO_PRODUCT}]->(product_species)
        WITH r AS r, reaction AS reaction
        MATCH
            (entry:CollectionEntry)-[:HAS_OBJ]->(:CellDesignerMap)
                -[:HAS_MODEL]->(model:CellDesignerModel)-[:HAS_REACTION]->(reaction),
            (entry)-[:HAS_SOURCE_ID_TO_MODEL_ELEMENT]->(ids:Mapping)
                -[:HAS_ITEM]->(reaction_item:Item),
            (reaction_item)-[:HAS_VALUE]->(reaction),
            (reaction_item)-[:HAS_KEY]->(reaction_id:String),
            (ids)-[:HAS_ITEM]->(model_item:Item),
            (model_item)-[:HAS_VALUE]->(model),
            (model_item)-[:HAS_KEY]->(model_id:String)
        WITH r AS r,
             collect(DISTINCT model_id.value) AS model_ids,
             collect(DISTINCT reaction_id.value) AS reaction_ids
        SET r.reaction_ids = reaction_ids,
            r.model_ids = model_ids
        RETURN r
        """
    )
    for labels, relationship_type in [
        (["PhysicalStimulator", "Catalyzer"], POSITIVE_INFLUENCE),
        (["Trigger"], NECESSARY_POSITIVE_INFLUENCE),
        (["Inhibitor"], NEGATIVE_INFLUENCE),
        (
            ["UnknownPhysicalStimulator", "UnknownCatalyzer"],
            UNCERTAIN_POSITIVE_INFLUENCE,
        ),
        (["UnknownInhibitor"], UNCERTAIN_NEGATIVE_INFLUENCE),
        (["UnknownTrigger"], UNCERTAIN_NECESSARY_POSITIVE_INFLUENCE),
    ]:
        for label in labels:
            queries.append(
                f"""
                MATCH
                    (reaction:Reaction),
                    (reaction)-[:HAS_MODIFIER]->(modulator:{label}),
                    (modulator)-[:HAS_REFERRED_SPECIES]->(modulator_species:Species),
                    (reaction)-[:HAS_PRODUCT]->(product:Product),
                    (product)-[:HAS_REFERRED_SPECIES]->(product_species:Species)
                MERGE
                    (modulator_species)-[r:{relationship_type}]->(product_species)
                WITH r AS r, reaction AS reaction
                MATCH
                    (entry:CollectionEntry)-[:HAS_OBJ]->(:CellDesignerMap)
                        -[:HAS_MODEL]->(model:CellDesignerModel)
                        -[:HAS_REACTION]->(reaction),
                    (entry)-[:HAS_SOURCE_ID_TO_MODEL_ELEMENT]->(ids:Mapping)
                        -[:HAS_ITEM]->(reaction_item:Item),
                    (reaction_item)-[:HAS_VALUE]->(reaction),
                    (reaction_item)-[:HAS_KEY]->(reaction_id:String),
                    (ids)-[:HAS_ITEM]->(model_item:Item),
                    (model_item)-[:HAS_VALUE]->(model),
                    (model_item)-[:HAS_KEY]->(model_id:String)
                WITH r AS r,
                     collect(DISTINCT model_id.value) AS model_ids,
                     collect(DISTINCT reaction_id.value) AS reaction_ids
                SET r.reaction_ids = reaction_ids,
                    r.model_ids = model_ids
                RETURN r
                """
            )
            queries.append(
                f"""
                MATCH
                    (reaction:Reaction),
                    (reaction)-[:HAS_MODIFIER]->(modulator:{label}),
                    (modulator)-[:HAS_REFERRED_SPECIES]->(modulator_gate:BooleanLogicGate),
                    (modulator_gate)-[:HAS_INPUT]->(modulator_species:Species),
                    (reaction)-[:HAS_PRODUCT]->(product:Product),
                    (product)-[:HAS_REFERRED_SPECIES]->(product_species:Species)
                WHERE NOT modulator_gate:NotGate
                MERGE
                    (modulator_species)-[r:{relationship_type}]->(product_species)
                WITH r AS r, reaction AS reaction
                MATCH
                    (entry:CollectionEntry)-[:HAS_OBJ]->(:CellDesignerMap)
                        -[:HAS_MODEL]->(model:CellDesignerModel)
                        -[:HAS_REACTION]->(reaction),
                    (entry)-[:HAS_SOURCE_ID_TO_MODEL_ELEMENT]->(ids:Mapping)
                        -[:HAS_ITEM]->(reaction_item:Item),
                    (reaction_item)-[:HAS_VALUE]->(reaction),
                    (reaction_item)-[:HAS_KEY]->(reaction_id:String),
                    (ids)-[:HAS_ITEM]->(model_item:Item),
                    (model_item)-[:HAS_VALUE]->(model),
                    (model_item)-[:HAS_KEY]->(model_id:String)
                WITH r AS r,
                     collect(DISTINCT model_id.value) AS model_ids,
                     collect(DISTINCT reaction_id.value) AS reaction_ids
                SET r.reaction_ids = reaction_ids,
                    r.model_ids = model_ids
                RETURN r
                """
            )
    for labels, relationship_type in [
        (
            ["Catalysis", "PhysicalStimulation", "PositiveInfluence"],
            POSITIVE_INFLUENCE,
        ),
        (["Triggering"], NECESSARY_POSITIVE_INFLUENCE),
        (["Inhibition", "NegativeInfluence"], NEGATIVE_INFLUENCE),
        (
            [
                "UnknownCatalysis",
                "UnknownPhysicalStimulation",
                "UnknownPositiveInfluence",
            ],
            UNCERTAIN_POSITIVE_INFLUENCE,
        ),
        (["UnknownTriggering"], UNCERTAIN_NECESSARY_POSITIVE_INFLUENCE),
        (
            ["UnknownInhibition", "UnknownNegativeInfluence"],
            UNCERTAIN_NEGATIVE_INFLUENCE,
        ),
    ]:
        for label in labels:
            queries.append(
                f"""
                MATCH
                    (modulation:{label}),
                    (modulation)-[:HAS_SOURCE]->(source_species:Species),
                    (modulation)-[:HAS_TARGET]->(target_species:Species)
                MERGE
                    (source_species)-[r:{relationship_type}]->(target_species)
                WITH r AS r, modulation AS modulation
                MATCH
                    (entry:CollectionEntry)-[:HAS_OBJ]->(:CellDesignerMap)
                        -[:HAS_MODEL]->(model:CellDesignerModel)
                        -[:HAS_MODULATION]->(modulation),
                    (entry)-[:HAS_SOURCE_ID_TO_MODEL_ELEMENT]->(ids:Mapping)
                        -[:HAS_ITEM]->(modulation_item:Item),
                    (modulation_item)-[:HAS_VALUE]->(modulation),
                    (modulation_item)-[:HAS_KEY]->(modulation_id:String),
                    (ids)-[:HAS_ITEM]->(model_item:Item),
                    (model_item)-[:HAS_VALUE]->(model),
                    (model_item)-[:HAS_KEY]->(model_id:String)
                WITH r AS r,
                     collect(DISTINCT model_id.value) AS model_ids,
                     collect(DISTINCT modulation_id.value) AS modulation_ids
                SET r.modulation_ids = modulation_ids,
                    r.model_ids = model_ids
                RETURN r
                """
            )
            queries.append(
                f"""
                MATCH
                    (modulation:{label}),
                    (modulation)-[:HAS_SOURCE]->(source:BooleanLogicGate),
                    (source)-[:HAS_INPUT]->(source_species:Species),
                    (modulation)-[:HAS_TARGET]->(target_species:Species)
                WHERE NOT source:NotGate
                MERGE
                    (source_species)-[r:{relationship_type}]->(target_species)
                WITH r AS r, modulation AS modulation
                MATCH
                    (entry:CollectionEntry)-[:HAS_OBJ]->(:CellDesignerMap)
                        -[:HAS_MODEL]->(model:CellDesignerModel)
                        -[:HAS_MODULATION]->(modulation),
                    (entry)-[:HAS_SOURCE_ID_TO_MODEL_ELEMENT]->(ids:Mapping)
                        -[:HAS_ITEM]->(modulation_item:Item),
                    (modulation_item)-[:HAS_VALUE]->(modulation),
                    (modulation_item)-[:HAS_KEY]->(modulation_id:String),
                    (ids)-[:HAS_ITEM]->(model_item:Item),
                    (model_item)-[:HAS_VALUE]->(model),
                    (model_item)-[:HAS_KEY]->(model_id:String)
                WITH r AS r,
                     collect(DISTINCT model_id.value) AS model_ids,
                     collect(DISTINCT modulation_id.value) AS modulation_ids
                SET r.modulation_ids = modulation_ids,
                    r.model_ids = model_ids
                RETURN r
                """
            )
    for query in queries:
        session.execute_query(query)


def make_ig_from_subgraph(subgraph):
    return make_ig_from_nodes_and_relationships(subgraph[0], subgraph[1])


def make_ig_from_nodes_and_relationships(nodes, relationships):
    ig = InfluenceGraph()
    for node in nodes:
        node._graph = None
        ig.add_node(node)
    for relationship in relationships:
        relationship.start_node._graph = None
        relationship.end_node._graph = None
        if (
            relationship.start_node in ig.get_nodes()
            and relationship.end_node in ig.get_nodes()
        ):
            ig.add_relationship(relationship)
    return ig


def prune_trivial_nodes(ig, recursive=True):
    ig = copy.deepcopy(ig)
    deleted = True
    while deleted:
        deleted = False
        for node in ig.get_nodes():
            if not ig.get_modulators(node):
                ig.remove_node(node)
                deleted = True
        if not recursive:
            break
    return ig


def _get_number_and_size_of_clusters(list_of_sets):
    i = 0
    while i < len(list_of_sets):
        j = i + 1
        while j < len(list_of_sets):
            if list_of_sets[i].intersection(list_of_sets[j]):
                list_of_sets[i] = list_of_sets[i].union(list_of_sets[j])
                del list_of_sets[j]
                j = i + 1
            else:
                j += 1
        i += 1
    return len(list_of_sets), [len(set_) for set_ in list_of_sets]


def get_number_and_size_of_components(ig):
    list_of_sets = []
    for node in ig.get_nodes():
        for modulator in ig.get_modulators(node):
            list_of_sets.append(set([node, modulator]))
    return _get_number_and_size_of_clusters(list_of_sets)


def _get_flatten_dot_nodes(dot_graph):
    nodes = list(dot_graph.get_nodes())
    for subgraph in dot_graph.get_subgraphs():
        nodes += _get_flatten_dot_nodes(subgraph)
    return nodes


def _get_coordinates_from_pydot_node(dot_node):
    pos = dot_node.get("pos")
    if pos is None:
        return None
    pos = pos.strip('"')
    return tuple(float(c) for c in pos.split(","))


def _translate_layout_element(layout_element, tx, ty):
    layout_element.position = layout_element.position + (tx, ty)
    for sub_layout_element in layout_element.children():
        _translate_layout_element(sub_layout_element, tx, ty)


RELATIONSHIP_TYPE_TO_CD_LAYOUT_CLASS = {
    POSITIVE_INFLUENCE: momapy.celldesigner.PositiveInfluenceLayout,
    NECESSARY_POSITIVE_INFLUENCE: momapy.celldesigner.TriggeringLayout,
    REACTANT_TO_PRODUCT: momapy.celldesigner.TriggeringLayout,
    NEGATIVE_INFLUENCE: momapy.celldesigner.InhibitionLayout,
}

RELATIONSHIP_TYPE_TO_CD_MODEL_CLASS = {
    POSITIVE_INFLUENCE: momapy.celldesigner.PositiveInfluence,
    NECESSARY_POSITIVE_INFLUENCE: momapy.celldesigner.Triggering,
    REACTANT_TO_PRODUCT: momapy.celldesigner.Triggering,
    NEGATIVE_INFLUENCE: momapy.celldesigner.NegativeInfluence,
}

# Single shared compartment for every species. We do not model compartments in
# the influence-graph view; collapsing everything into one `default`
# compartment (no `outside` parent) sidesteps the CellDesigner reader's
# compartment-ordering, which crashes on dangling `outside` references.
DEFAULT_COMPARTMENT = momapy.celldesigner.Compartment(id_="default", name="default")


def _walk_subunits(species):
    """Yield every species transitively contained in ``species.subunits``."""
    for subunit in getattr(species, "subunits", None) or ():
        yield subunit
        yield from _walk_subunits(subunit)


def _iter_templates(species):
    """Yield the template of ``species`` and of every (transitive) subunit.

    Subunit (included-species) protein templates must also be declared in the
    model's ``species_templates`` for the writer to emit them.
    """
    template = getattr(species, "template", None)
    if template is not None:
        yield template
    for subunit in getattr(species, "subunits", None) or ():
        yield from _iter_templates(subunit)


def _canon_template(template, template_canon):
    """Return a content-canonical copy of ``template`` with a globally unique
    ``id_``.

    Species reconstructed from different source maps reuse the same template
    ids (e.g. ``p_4`` names a different protein in each map). The CellDesigner
    writer keys xml ids by ``id_``, so two content-distinct templates sharing an
    ``id_`` would be emitted as one ``<protein>`` and the reader would then find
    several model elements under that id (``get_one`` raises). Interning by
    content — ``id_`` is excluded from template equality/hashing — collapses
    content-equal templates to one instance and gives content-distinct
    templates distinct ids.
    """
    existing = template_canon.get(template)
    if existing is not None:
        return existing
    rebuilt = dataclasses.replace(
        template, id_=f"species_template_{len(template_canon)}"
    )
    template_canon[template] = rebuilt
    return rebuilt


def _register_or_reuse(element, cache):
    """Intern ``element`` by content in ``cache`` (pd2af's ``register_or_reuse``).

    Frozen momapy dataclasses exclude ``id_`` from equality/hashing, so two
    content-equal elements collapse to a single Python identity end to end.
    This is what keeps ``model.species`` / ``model.modulations`` and the
    layout-to-model mapping in agreement: every modulation endpoint and every
    mapped glyph references the same canonical instance that lives in the model,
    so nothing is silently dropped by the model's ``frozenset`` fields.
    """
    existing = cache.get(element)
    if existing is not None:
        return existing
    cache[element] = element
    return element


def _normalize_species(species, top_level, species_canon, template_canon):
    """Return the content-canonical, compartment-normalized form of ``species``.

    Top-level species get the shared ``DEFAULT_COMPARTMENT``; species reached as
    a complex's subunit get ``compartment=None``. This single difference keeps a
    free-floating entity distinct from its complexed twin (so the writer emits
    the free one top-level and the complexed one only as an included species)
    without the fragile id-based "promotion" — content interning then does the
    rest. Templates are content-canonicalized too, and the result is interned in
    ``species_canon`` so content-equal species share one instance.
    """
    replacements = {
        "compartment": DEFAULT_COMPARTMENT if top_level else None
    }
    template = getattr(species, "template", None)
    if template is not None:
        replacements["template"] = _canon_template(template, template_canon)
    subunits = getattr(species, "subunits", None)
    if subunits:
        replacements["subunits"] = frozenset(
            _normalize_species(subunit, False, species_canon, template_canon)
            for subunit in subunits
        )
    normalized = dataclasses.replace(species, **replacements)
    existing = species_canon.get(normalized)
    if existing is not None:
        return existing
    # Give the canonical a globally unique ``id_``. Species ids from different
    # source maps are namespaced and so usually distinct, but interning by
    # content (``id_`` excluded) plus a fresh unique id guarantees the writer
    # never emits two content-distinct species under one id (the reader indexes
    # model elements by id and ``get_one`` raises on a collision).
    canonical = dataclasses.replace(normalized, id_=f"species_{len(species_canon)}")
    species_canon[canonical] = canonical
    return canonical


def _populate_ig_layout(
    session,
    ig,
    layout_builder,
    color_node_ids,
    label,
    extra_node_layout_elements,
    node_id_to_species=None,
    obj_cache=None,
    species_canon=None,
    template_canon=None,
):
    """Fill ``layout_builder`` with the glyphs and arcs of ``ig``.

    Glyphs are fetched from the DB (or taken from
    ``extra_node_layout_elements``), positioned with graphviz, and arcs are
    created for every influence relationship. Returns:
      - ``node_id_to_layout_element``: {node["id_"]: glyph builder} (moved
        to their final position),
      - ``arc_pairs``: list of ``(arc_builder, relationship)`` for every arc
        created, so callers can build a layout-to-model mapping, and
      - ``subunit_mappings``: list of ``(child_layout_builder,
        subunit_model_element)`` for the subunit glyphs of complex nodes, so a
        complex is drawn with its subunits. Empty unless ``node_id_to_species``
        / ``obj_cache`` / ``species_canon`` are supplied (the model-building path).
    """
    relationship_type_to_cd_class = RELATIONSHIP_TYPE_TO_CD_LAYOUT_CLASS
    arc_pairs = []
    ig_nodes = list(ig.get_nodes())
    db_nodes = [n for n in ig_nodes if n not in extra_node_layout_elements]
    element_ids = [n.element_id for n in db_nodes]
    rows = session.cypher_query_as_layout_elements(
        "UNWIND $eids AS eid MATCH (n) WHERE elementId(n) = eid RETURN n AS node",
        params={"eids": element_ids},
    )
    node_to_layout_element = {}
    for ig_node, layout_elements in zip(db_nodes, rows):
        if not layout_elements:
            raise RuntimeError(
                f"no layout for IG node id_={ig_node.get('id_')} "
                f"labels={set(ig_node.labels)}"
            )
        node_to_layout_element[ig_node] = layout_elements[0]
    for ig_node, layout_element in extra_node_layout_elements.items():
        node_to_layout_element[ig_node] = layout_element
    dot_graph = pydot.Dot(graph_type="digraph")
    id_to_dot_node = {}
    for node in ig_nodes:
        layout_element = node_to_layout_element[node]
        own_bbox = layout_element.own_bbox()
        dot_node = pydot.Node(node["id_"])
        dot_node.set("width", own_bbox.width / POINTS_PER_INCH)
        dot_node.set("height", own_bbox.height / POINTS_PER_INCH)
        dot_graph.add_node(dot_node)
        id_to_dot_node[node["id_"]] = dot_node
    for node in ig_nodes:
        for modulator in ig.get_modulators(node):
            dot_graph.add_edge(pydot.Edge(modulator["id_"], node["id_"]))
    dot_graph.set("ranksep", 1.0)
    dot_graph.set("nodesep", 0.5)
    dot_graph.set("rankdir", "BT")
    dot = dot_graph.create_dot(prog="dot").decode("utf-8")
    dot_graph = pydot.graph_from_dot_data(dot)[0]
    node_id_to_coordinates = {}
    for dot_node in _get_flatten_dot_nodes(dot_graph):
        dot_node_id = dot_node.get_name().strip('"')
        if dot_node_id in ("graph", "node", "edge"):
            continue
        coords = _get_coordinates_from_pydot_node(dot_node)
        if coords is not None:
            node_id_to_coordinates[dot_node_id] = coords
    node_id_to_layout_element_moved = {}
    for node, layout_element in node_to_layout_element.items():
        layout_element_builder = momapy.builder.builder_from_object(layout_element)
        coordinates = node_id_to_coordinates[node["id_"]]
        old_position = layout_element_builder.position
        new_position = momapy.geometry.Point.from_tuple(coordinates)
        tx, ty = new_position - old_position
        _translate_layout_element(layout_element_builder, tx, ty)
        layout_builder.layout_elements.append(layout_element_builder)
        node_id_to_layout_element_moved[node["id_"]] = layout_element_builder
    directed_pairs = {
        (r.start_node["id_"], r.end_node["id_"]) for r in ig.get_relationships()
    }
    bezier_offset = 30.0
    for relationship in ig.get_relationships():
        start_node = relationship.start_node
        end_node = relationship.end_node
        if (
            start_node["id_"] not in node_id_to_layout_element_moved
            or end_node["id_"] not in node_id_to_layout_element_moved
        ):
            continue
        start_layout_element = node_id_to_layout_element_moved[start_node["id_"]]
        end_layout_element = node_id_to_layout_element_moved[end_node["id_"]]
        arc = momapy.builder.new_builder_object(
            relationship_type_to_cd_class[relationship.type]
        )
        # Record the specific endpoint glyphs on the arc (as pd2af does). The
        # writer reads `arc.source`/`arc.target` to pick the baseReactant /
        # baseProduct aliases; without them it falls back to the modulation's
        # source/target species and picks that species' *first* glyph for both
        # ends. With content interning a single species can have several glyphs,
        # so that fallback would emit the wrong alias — and for an edge between
        # two glyphs of the same species it would emit the *same* alias for both
        # ends, a fake self-loop the reader cannot lay out.
        arc.source = start_layout_element
        arc.target = end_layout_element
        start_id = start_node["id_"]
        end_id = end_node["id_"]
        is_self_loop = start_id == end_id
        is_bidirectional = (
            not is_self_loop and (end_id, start_id) in directed_pairs
        )
        if is_self_loop:
            # Build the self-loop between two NAMED anchors (as pd2af does) so
            # the writer can match the segment endpoints to anchor names and
            # emit `<linkAnchor position="NNW/NNE">`. Border/own_angle points
            # do not match a named anchor, so the writer would fall back to the
            # "center" anchor and the reader's own_border(center) would return
            # None (zero-length line) and crash.
            start_point = start_layout_element.anchor_point("north_north_west")
            end_point = start_layout_element.anchor_point("north_north_east")
            arc.segments = [momapy.geometry.Segment(start_point, end_point)]
            layout_builder.layout_elements.append(arc)
            arc_pairs.append((arc, relationship))
            continue
        if is_bidirectional:
            start_center = start_layout_element.center()
            end_center = end_layout_element.center()
            delta_x = end_center.x - start_center.x
            delta_y = end_center.y - start_center.y
            length = math.hypot(delta_x, delta_y)
            if length == 0:
                start_point = start_layout_element.own_border(end_center)
                end_point = end_layout_element.own_border(start_center)
                if start_point is None:
                    start_point = start_layout_element.north_west()
                if end_point is None:
                    end_point = end_layout_element.north_east()
                arc.segments = [momapy.geometry.Segment(start_point, end_point)]
            else:
                normal_x = -delta_y / length
                normal_y = delta_x / length
                middle_x = (start_center.x + end_center.x) / 2
                middle_y = (start_center.y + end_center.y) / 2
                control_point = momapy.geometry.Point(
                    middle_x + bezier_offset * normal_x,
                    middle_y + bezier_offset * normal_y,
                )
                start_point = start_layout_element.own_border(control_point)
                end_point = end_layout_element.own_border(control_point)
                if start_point is None:
                    start_point = start_layout_element.north_west()
                if end_point is None:
                    end_point = end_layout_element.north_east()
                arc.segments = [
                    momapy.geometry.QuadraticBezierCurve(
                        start_point, end_point, control_point
                    )
                ]
        else:
            start_point = start_layout_element.own_border(end_layout_element.center())
            end_point = end_layout_element.own_border(start_layout_element.center())
            if start_point is None:
                start_point = start_layout_element.north_west()
            if end_point is None:
                end_point = end_layout_element.north_east()
            arc.segments = [momapy.geometry.Segment(start_point, end_point)]
        layout_builder.layout_elements.append(arc)
        arc_pairs.append((arc, relationship))
    for node_ids, color_name in color_node_ids:
        color = getattr(momapy.coloring, color_name)
        for node_id in node_ids:
            layout_element_builder = node_id_to_layout_element_moved.get(node_id)
            if layout_element_builder is not None:
                layout_element_builder.fill = color
    if label is not None:
        bbox = momapy.positioning.fit(layout_builder.layout_elements)
        text_layout = momapy.core.TextLayout(
            text=label, position=momapy.positioning.below_of(bbox.south(), 50)
        )
        layout_builder.layout_elements.append(text_layout)
    momapy.positioning.set_fit(
        layout_builder, layout_builder.layout_elements, xsep=15.0, ysep=15.0
    )
    subunit_mappings = _collect_complex_subunit_mappings(
        session,
        ig_nodes,
        extra_node_layout_elements,
        node_id_to_layout_element_moved,
        node_id_to_species,
        obj_cache,
        species_canon,
        template_canon,
    )
    return node_id_to_layout_element_moved, arc_pairs, subunit_mappings


def _collect_complex_subunit_mappings(
    session,
    ig_nodes,
    extra_node_layout_elements,
    node_id_to_layout_element,
    node_id_to_species,
    obj_cache,
    species_canon,
    template_canon,
):
    """For every complex node, pair each subunit glyph (a child of the complex's
    ``ComplexLayout``) with its subunit model element.

    The CellDesigner writer only draws a subunit when ``mapping.get_mapping``
    yields a layout that is an actual child of the complex layout (writer.py
    ``_collect_complex_aliases``). We get those pairs from the **input maps'
    own** ``LayoutModelMapping`` stored in the DB, read in the *forward*
    (alias -> model element) direction: each child alias of the complex's
    ``ComplexLayout`` has exactly one subunit species as its ``HAS_VALUE``.

    Reading forward avoids the species -> aliases ambiguity that a reverse
    lookup suffers from. A subunit species is shared across every complex and
    map it appears in (content interning), so it resolves to *many* aliases and
    there is no way to tell which one belongs to *this* complex; the child
    alias, by contrast, names its subunit unambiguously.

    The reconstructed ``complex_glyph`` lost its source-node identity, so we
    join the stored child rows to its child builders by ``id_`` (unique within
    a single source map) and, for a complex drawn in several maps, select the
    source ``ComplexLayout`` whose child set overlaps the reconstructed glyph's
    children the most. Each subunit model element is normalized as a subunit
    (``top_level=False``) and interned in the shared ``species_canon``, so it is
    the very canonical instance held in the complex's ``subunits``.

    Returns a list of ``(child_layout_builder, normalized_subunit_model)``.
    """
    subunit_mappings = []
    if (
        node_id_to_species is None
        or obj_cache is None
        or species_canon is None
        or template_canon is None
    ):
        return subunit_mappings
    for node in ig_nodes:
        if node in extra_node_layout_elements:
            continue
        species = node_id_to_species.get(node["id_"])
        if not isinstance(species, momapy.celldesigner.Complex):
            continue
        complex_glyph = node_id_to_layout_element.get(node["id_"])
        if complex_glyph is None:
            continue
        id_to_child = {}
        stack = list(getattr(complex_glyph, "layout_elements", []) or [])
        while stack:
            child = stack.pop()
            child_id = getattr(child, "id_", None)
            if child_id is not None:
                id_to_child[child_id] = child
            stack.extend(getattr(child, "layout_elements", []) or [])
        if not id_to_child:
            continue
        # Read the input maps' stored LayoutModelMapping in the forward
        # direction: for every child alias of every ComplexLayout of this
        # complex species, the subunit species it was mapped to in its source
        # map. One alias maps to exactly one subunit, so this sidesteps the
        # shared-subunit -> many-aliases ambiguity of the reverse lookup.
        # ``subunit:Species`` also matches nested complexes (CellDesigner
        # Complex is a Species subclass), so HAS_LAYOUT_ELEMENT* recursion pairs
        # subunits at every nesting level.
        mapping_rows = session.execute_query(
            "MATCH (c) WHERE elementId(c) = $eid "
            "MATCH (c)<-[:HAS_VALUE]-(:Item)-[:HAS_KEY]->(cl:ComplexLayout) "
            "MATCH (cl)-[:HAS_LAYOUT_ELEMENT*]->(child)"
            "<-[:HAS_KEY]-(:Item)-[:HAS_VALUE]->(subunit:Species) "
            "RETURN elementId(cl) AS alias_id, child.id_ AS child_id, "
            "subunit AS subunit",
            params={"eid": node.element_id},
        )
        alias_to_children = {}
        for mapping_row in mapping_rows:
            child_id = mapping_row["child_id"]
            if child_id is None:
                continue
            alias_to_children.setdefault(mapping_row["alias_id"], {})[
                child_id
            ] = mapping_row["subunit"]
        if not alias_to_children:
            continue
        # The reconstructed complex_glyph came from one source map's alias; pick
        # the alias whose child ids overlap the glyph's children the most (id_
        # is unique within a map but collides across maps, so the whole set, not
        # a single id, identifies the source alias).
        source_alias = max(
            alias_to_children,
            key=lambda alias: len(
                alias_to_children[alias].keys() & id_to_child.keys()
            ),
        )
        child_id_to_subunit_node = {
            child_id: subunit_node
            for child_id, subunit_node in alias_to_children[source_alias].items()
            if child_id in id_to_child
        }
        if not child_id_to_subunit_node:
            continue
        # Reconstruct each mapped subunit species once; normalizing as a subunit
        # interns it (via species_canon) to the same canonical instance held in
        # the complex's ``subunits``, so the writer pairs glyph and model by
        # identity.
        subunit_eids = list(
            {node_.element_id for node_ in child_id_to_subunit_node.values()}
        )
        model_rows = session.execute_query_as_objects(
            "UNWIND $eids AS eid MATCH (n) WHERE elementId(n) = eid "
            "RETURN n AS node",
            params={"eids": subunit_eids},
            node_id_to_object=obj_cache,
        )
        eid_to_subunit_model = {}
        for subunit_eid, objects in zip(subunit_eids, model_rows):
            if not objects:
                continue
            eid_to_subunit_model[subunit_eid] = _normalize_species(
                objects[0], False, species_canon, template_canon
            )
        for child_id, subunit_node in child_id_to_subunit_node.items():
            child_builder = id_to_child.get(child_id)
            subunit_model = eid_to_subunit_model.get(subunit_node.element_id)
            if child_builder is None or subunit_model is None:
                continue
            subunit_mappings.append((child_builder, subunit_model))
    return subunit_mappings


def make_map_layout_from_ig(
    session,
    ig,
    color_node_ids: list[tuple[list[str], str]] | None = None,
    label=None,
    extra_node_layout_elements=None,
):
    if color_node_ids is None:
        color_node_ids = []
    if extra_node_layout_elements is None:
        extra_node_layout_elements = {}
    layout_builder = momapy.builder.new_builder_object(
        momapy.celldesigner.CellDesignerLayout
    )
    _populate_ig_layout(
        session,
        ig,
        layout_builder,
        color_node_ids,
        label,
        extra_node_layout_elements,
    )
    return momapy.builder.object_from_builder(layout_builder)


def _make_ig_model(
    session, ig, extra_node_model_elements, obj_cache, species_canon, template_canon
):
    """Build the CellDesigner model (species + templates + modulations) of an IG.

    Mirrors pd2af's model pass: each IG node's species is reconstructed, its
    compartment normalized and its template canonicalized, then **interned by
    content** (``_register_or_reuse``). Because momapy excludes ``id_`` from
    equality, content-equal species from different source maps collapse to a
    single canonical instance — several IG nodes can therefore share one model
    species (drawn later as several aliases), and crucially every modulation
    endpoint and every mapped glyph then references the exact instance that
    lives in ``model.species`` (nothing is dropped by the model's frozenset).
    Modulations are interned the same way. Returns
    ``(node_id_to_species, model, relationship_to_modulation)``.
    """
    node_id_to_species = {}
    ig_nodes = list(ig.get_nodes())
    db_nodes = [n for n in ig_nodes if n not in extra_node_model_elements]
    element_ids = [n.element_id for n in db_nodes]
    rows = session.execute_query_as_objects(
        "UNWIND $eids AS eid MATCH (n) WHERE elementId(n) = eid RETURN n AS node",
        params={"eids": element_ids},
        node_id_to_object=obj_cache,
    )
    for ig_node, objects in zip(db_nodes, rows):
        if not objects:
            raise RuntimeError(
                f"no model element for IG node id_={ig_node.get('id_')} "
                f"labels={set(ig_node.labels)}"
            )
        node_id_to_species[ig_node["id_"]] = _normalize_species(
            objects[0], True, species_canon, template_canon
        )
    for ig_node, species in extra_node_model_elements.items():
        node_id_to_species[ig_node["id_"]] = _normalize_species(
            species, True, species_canon, template_canon
        )

    species = list(node_id_to_species.values())
    species_templates = frozenset(
        template for s in species for template in _iter_templates(s)
    )
    relationship_to_modulation = {}
    modulation_canon = {}
    for relationship in ig.get_relationships():
        modulation_class = RELATIONSHIP_TYPE_TO_CD_MODEL_CLASS.get(relationship.type)
        if modulation_class is None:
            continue
        source = node_id_to_species.get(relationship.start_node["id_"])
        target = node_id_to_species.get(relationship.end_node["id_"])
        if source is None or target is None:
            continue
        relationship_to_modulation[relationship] = _register_or_reuse(
            modulation_class(source=source, target=target), modulation_canon
        )

    model = momapy.celldesigner.CellDesignerModel(
        species=frozenset(species),
        species_templates=species_templates,
        compartments=frozenset([DEFAULT_COMPARTMENT]),
        modulations=frozenset(relationship_to_modulation.values()),
    )
    return node_id_to_species, model, relationship_to_modulation


def _assign_unique_layout_ids(layout_builder):
    """Give every layout element (and descendant) a globally unique ``id_``.

    Layout-element (alias) ids are not namespaced across the source maps, so
    glyphs from different maps collide on ``id_`` (e.g. ``pdme90`` is a
    free-standing CASP1 alias in one map and a complex-subunit CASP1 alias in
    another). The writer emits ``<speciesAlias id=...>`` from the glyph ``id_``
    and the reader indexes model elements by alias id, so a duplicate makes the
    reader resolve one alias to two species. The layout-to-model mapping is keyed
    by object identity, so renumbering ids here is safe.
    """
    counter = 0
    stack = list(layout_builder.layout_elements)
    while stack:
        element = stack.pop()
        if hasattr(element, "id_"):
            # Glyphs and arcs are builders (mutable); the only frozen layout
            # elements are extras created directly, e.g. the map-label
            # TextLayout, whose uuid id never collides — skip those.
            try:
                element.id_ = f"le{counter}"
                counter += 1
            except dataclasses.FrozenInstanceError:
                pass
        children = getattr(element, "layout_elements", None)
        if children:
            stack.extend(children)


def make_celldesigner_map_from_ig(
    session,
    ig,
    color_node_ids: list[tuple[list[str], str]] | None = None,
    label=None,
    extra_node_layout_elements=None,
    extra_node_model_elements=None,
):
    """Build a full ``CellDesignerMap`` (model + layout + mapping) from an IG.

    Mirrors pd2af's ``build_map``: a model pass (``_make_ig_model``) reconstructs
    the species and builds the ``CellDesignerModel`` and modulations; a layout
    pass (``_populate_ig_layout``) builds the glyphs, the influence arcs and the
    complex-subunit glyph pairs; then the ``LayoutModelMapping`` ties them
    together and the whole map is finalized in a single ``object_from_builder``
    pass so layout builders used as mapping keys resolve to the same frozen
    objects as in the layout. ``obj_cache`` / ``species_canon`` are shared across
    both passes so a complex's subunit and its mapped layout agree on identity.
    """
    if color_node_ids is None:
        color_node_ids = []
    if extra_node_layout_elements is None:
        extra_node_layout_elements = {}
    if extra_node_model_elements is None:
        extra_node_model_elements = {}

    obj_cache = {}
    species_canon = {}
    template_canon = {}

    node_id_to_species, model, relationship_to_modulation = _make_ig_model(
        session, ig, extra_node_model_elements, obj_cache, species_canon, template_canon
    )

    layout_builder = momapy.builder.new_builder_object(
        momapy.celldesigner.CellDesignerLayout
    )
    node_id_to_layout_element, arc_pairs, subunit_mappings = _populate_ig_layout(
        session,
        ig,
        layout_builder,
        color_node_ids,
        label,
        extra_node_layout_elements,
        node_id_to_species=node_id_to_species,
        obj_cache=obj_cache,
        species_canon=species_canon,
        template_canon=template_canon,
    )
    # Renumber layout-element ids to be globally unique (the complex-subunit
    # matching in _populate_ig_layout has already run on the original DB ids).
    _assign_unique_layout_ids(layout_builder)

    mapping_builder = momapy.core.mapping.LayoutModelMappingBuilder()
    for node_id, layout_element in node_id_to_layout_element.items():
        modeled_species = node_id_to_species.get(node_id)
        if modeled_species is not None:
            mapping_builder.add_mapping(layout_element, modeled_species)
    for child_layout_element, subunit_model in subunit_mappings:
        mapping_builder.add_mapping(child_layout_element, subunit_model)
    for arc, relationship in arc_pairs:
        modulation = relationship_to_modulation.get(relationship)
        if modulation is None:
            continue
        source_le = node_id_to_layout_element.get(relationship.start_node["id_"])
        target_le = node_id_to_layout_element.get(relationship.end_node["id_"])
        if source_le is None or target_le is None:
            continue
        mapping_builder.add_mapping(
            frozenset([arc, source_le, target_le]),
            modulation,
            anchor=arc,
        )

    map_builder = momapy.builder.new_builder_object(
        momapy.celldesigner.CellDesignerMap
    )
    map_builder.model = model
    map_builder.layout = layout_builder
    map_builder.layout_model_mapping = mapping_builder
    return momapy.builder.object_from_builder(map_builder)
