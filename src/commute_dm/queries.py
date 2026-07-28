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

    The AF traversal (`commute_dm.ig`) works on element ids only; this is the
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
