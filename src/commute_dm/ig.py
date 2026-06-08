import copy
import math

import pydot
import momapy.geometry
import momapy.core
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

    relationship_type_to_cd_class = {
        POSITIVE_INFLUENCE: momapy.celldesigner.PositiveInfluenceLayout,
        NECESSARY_POSITIVE_INFLUENCE: momapy.celldesigner.TriggeringLayout,
        REACTANT_TO_PRODUCT: momapy.celldesigner.TriggeringLayout,
        NEGATIVE_INFLUENCE: momapy.celldesigner.InhibitionLayout,
    }
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
    layout_builder = momapy.builder.new_builder_object(
        momapy.celldesigner.CellDesignerLayout
    )
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
        start_id = start_node["id_"]
        end_id = end_node["id_"]
        is_self_loop = start_id == end_id
        is_bidirectional = (
            not is_self_loop and (end_id, start_id) in directed_pairs
        )
        if is_self_loop:
            start_point = start_layout_element.own_angle(120)
            end_point = start_layout_element.own_angle(60)
            if start_point is None:
                start_point = start_layout_element.north_west()
            if end_point is None:
                end_point = start_layout_element.north_east()
            center = start_layout_element.center()
            start_dx, start_dy = start_point.x - center.x, start_point.y - center.y
            end_dx, end_dy = end_point.x - center.x, end_point.y - center.y
            start_length = math.hypot(start_dx, start_dy) or 1.0
            end_length = math.hypot(end_dx, end_dy) or 1.0
            start_control_point = momapy.geometry.Point(
                start_point.x + start_dx / start_length * bezier_offset,
                start_point.y + start_dy / start_length * bezier_offset,
            )
            end_control_point = momapy.geometry.Point(
                end_point.x + end_dx / end_length * bezier_offset,
                end_point.y + end_dy / end_length * bezier_offset,
            )
            arc.segments = [
                momapy.geometry.CubicBezierCurve(
                    start_point, end_point, start_control_point, end_control_point
                )
            ]
            layout_builder.layout_elements.append(arc)
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
    layout = momapy.builder.object_from_builder(layout_builder)
    return layout
