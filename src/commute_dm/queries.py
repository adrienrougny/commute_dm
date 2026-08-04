import fieldz_kb.lpg.core
import momapy.celldesigner


def prewarm_session(session):
    """Ensure pylpg node classes exist for every momapy type that may appear
    in saved CellDesigner maps (model elements, layout elements, glyphs,
    shapes, drawing primitives). Required before hydrating query results,
    because the session only registers node classes for types it has saved
    or for types reachable via static type hints — concrete subclasses are
    missed by the static walk.
    """
    import momapy.core.elements
    import momapy.core.layout
    import momapy.core.model
    import momapy.drawing
    import momapy.geometry
    import momapy.coloring

    ctx = session._context
    fieldz_kb.lpg.core.get_or_make_node_class_from_type(
        ctx, momapy.celldesigner.CellDesignerMap
    )

    def _all_subclasses(cls):
        seen = set()
        stack = [cls]
        while stack:
            current = stack.pop()
            for sub in current.__subclasses__():
                if sub not in seen:
                    seen.add(sub)
                    stack.append(sub)
        return seen

    bases = []
    for module in (
        momapy.core.elements,
        momapy.core.layout,
        momapy.core.model,
        momapy.drawing,
        momapy.geometry,
        momapy.coloring,
    ):
        for attr in vars(module).values():
            if isinstance(attr, type):
                bases.append(attr)

    for base in bases:
        try:
            fieldz_kb.lpg.core.get_or_make_node_class_from_type(ctx, base)
        except Exception:
            pass
        for cls in _all_subclasses(base):
            try:
                fieldz_kb.lpg.core.get_or_make_node_class_from_type(ctx, cls)
            except Exception:
                pass


def get_collections_for_nodes(session, nodes):
    """Return {node: set[collection_name]} for each input node.

    Uses the HAS_MODEL_ELEMENT membership edge emitted when collections are
    saved with `with_membership_edges=True`.
    """
    element_ids = [n.element_id for n in nodes]
    query = """
        UNWIND $eids AS eid
        MATCH (n) WHERE elementId(n) = eid
        OPTIONAL MATCH
            (n)<-[:HAS_MODEL_ELEMENT]-(model)<-[:HAS_MODEL]-(:CellDesignerMap)
                <-[:HAS_OBJ]-(:CollectionEntry)<-[:HAS_ENTRY]-(c:Collection)
        RETURN n AS node, collect(DISTINCT c.name) AS collections
    """
    result = session.execute_query(query, params={"eids": element_ids})
    return {row["node"]: set(row["collections"]) for row in result}


def get_subunits(session, nodes, recursive=True):
    node_element_ids = [node.element_id for node in nodes]
    cardinality = "*" if recursive else ""
    query = f"""
        MATCH
            (node:ModelElement)
        WHERE
            elementId(node) IN {node_element_ids}
        OPTIONAL MATCH
            (node)-[:HAS_SUBUNIT{cardinality}]->(subunit)
        RETURN node AS node, collect(subunit) AS subunits
    """
    result = session.execute_query(query)
    return [(row["node"], row["subunits"]) for row in result]


# The edges along which a node's annotations may sit on *another* node.
# `HAS_SUBUNIT` is the CellDesigner one; the rest are
# `commute_dm.bel_terms.MEMBER_EDGE_TYPES`, the `HAS__*` edges by which a BEL term
# holds another entity -- a complex's or composite's members, and an activity's
# subject. Repeated here rather than imported so this module keeps depending on
# none of the BEL machinery; `bel_terms` is where the classification is decided,
# and a new member edge type has to be added in both.
_ANNOTATED_NODE_EDGE_TYPES = (
    "HAS_SUBUNIT",
    "HAS__PROTEIN",
    "HAS__ABUNDANCE",
    "HAS__COMPLEX",
    "HAS__COMPOSITE",
    "HAS__GENE",
    "HAS__RNA",
    "HAS__MICRO_RNA",
)
# A BEL complex nests a handful of levels deep at most; the bound is a guard
# against a cycle in the term graph, not a modelling choice.
_ANNOTATED_NODE_MAX_DEPTH = 10


def get_annotated_nodes(session, nodes, with_subunits=False):
    """`[(node, [node, ...]), ...]`: the nodes whose annotations stand for each input.

    A CellDesigner species carries its own cross-references, so it is always its
    own annotated node, and its subunits are reached by `HAS_SUBUNIT` when
    `with_subunits`. A BEL node is different: the cross-references are attached
    (by `2_00`) to the `:Protein` terms, which the nodes the influence walk
    actually selects -- mostly `Activity` and `Complex` -- only *contain*. So the
    `HAS__*` member edges are followed, with the same distinction the
    CellDesigner side makes:

    * an **activity's subject is always followed**: `act(p(X))` is X in another
      form, the way an active CellDesigner species is still that species, so its
      genes are its subject's whether or not subunits are wanted;
    * a **complex's or composite's members are followed only when
      `with_subunits`**, exactly as `HAS_SUBUNIT` is.

    The input node is always included in its own list, which is what makes the
    caller uniform: a node that is itself annotated needs no traversal, and one
    that is not simply contributes nothing of its own.
    """
    node_element_ids = [node.element_id for node in nodes]
    if not node_element_ids:
        return []
    edge_types = "|".join(_ANNOTATED_NODE_EDGE_TYPES)
    query = f"""
        UNWIND $element_ids AS element_id
        MATCH (node) WHERE elementId(node) = element_id
        OPTIONAL MATCH
            path = (node)-[:{edge_types}*1..{_ANNOTATED_NODE_MAX_DEPTH}]->(part)
        WHERE
            all(inner IN nodes(path)[0..-1] WHERE
                CASE
                    WHEN inner:BELModelElement THEN inner:Activity OR $with_subunits
                    ELSE $with_subunits
                END)
        RETURN node AS node, collect(DISTINCT part) AS parts
    """
    result = session.execute_query(
        query,
        params={
            "element_ids": node_element_ids,
            "with_subunits": bool(with_subunits),
        },
    )
    return [(row["node"], [row["node"]] + row["parts"]) for row in result]


def get_identifiers(session, nodes, namespace="ncbigene"):
    formatted_result = []
    for node, annotations in get_annotations(session, nodes):
        identifiers = []
        for annotation in annotations:
            for resource in annotation["resources"]:
                if namespace in resource:
                    identifier = resource.split(":")[-1]
                    identifiers.append(identifier)
        formatted_result.append((node, identifiers))
    return formatted_result


def get_annotations(session, nodes):
    node_element_ids = [node.element_id for node in nodes]
    query = f"""
        MATCH
            (node:ModelElement)
        WHERE
            elementId(node) IN {node_element_ids}
        OPTIONAL MATCH
            (node)<-[:HAS_KEY]-(node_item:Item)-[:HAS_VALUE]->(node_bag:Bag)-[:HAS_ITEM]->(node_annotation:RDFAnnotation)
        RETURN node AS node, collect(DISTINCT node_annotation) AS annotations
    """
    result = session.execute_query(query)
    return [(row["node"], row["annotations"]) for row in result]


def get_ids_and_context(session, nodes):
    node_element_ids = [node.element_id for node in nodes]
    query = f"""
        MATCH
            (node:ModelElement)
        WHERE
            elementId(node) IN {node_element_ids}
        OPTIONAL MATCH
            (node)<-[:HAS_VALUE]-(node_item:Item)<-[:HAS_ITEM]-(ids:Mapping)<-[:HAS_SOURCE_ID_TO_MODEL_ELEMENT]-(entry:CollectionEntry)<-[:HAS_ENTRY]-(collection:Collection),
            (node_item)-[:HAS_KEY]->(node_id:String)
        RETURN node AS node,
               [x IN collect([node_id.value, entry, collection]) WHERE x[1] IS NOT NULL] AS ids_and_context
    """
    result = session.execute_query(query)
    formatted_result = []
    for row in result:
        node = row["node"]
        ids_and_context = [tuple(element) for element in row["ids_and_context"]]
        formatted_result.append((node, ids_and_context))
    return formatted_result


def get_nodes(session, element_ids):
    """Return the DB nodes for the given element ids.

    The AF walk (`commute_dm.submaps`) works on node ids only; this is the
    bridge to the node-taking helpers here and in `commute_dm.gea`.
    """
    element_ids = list(element_ids)
    if not element_ids:
        return []
    query = """
        UNWIND $element_ids AS element_id
        MATCH (node) WHERE elementId(node) = element_id
        RETURN node AS node
    """
    result = session.execute_query(query, params={"element_ids": element_ids})
    return [row["node"] for row in result]
