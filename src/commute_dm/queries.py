def get_collections_for_nodes(session, nodes):
    """Return {node: set[collection_name]} for each input node.

    Uses the HAS_MEMBER_MODEL_ELEMENT membership edge emitted when
    collections are saved with `with_membership_edges=True`.
    """
    element_ids = [n.element_id for n in nodes]
    query = """
        UNWIND $eids AS eid
        MATCH (n) WHERE elementId(n) = eid
        OPTIONAL MATCH
            (n)<-[:HAS_MEMBER_MODEL_ELEMENT]-(model)<-[:HAS_MODEL]-(:CellDesignerMap)
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


# A complex nests a handful of levels deep at most; the bound is a guard against
# a cycle, not a modelling choice.
_SUBUNIT_MAX_DEPTH = 10


def get_annotated_nodes(session, nodes, with_subunits=False):
    """`[(node, [node, ...]), ...]`: the nodes whose annotations stand for each input.

    A CellDesigner species carries its own cross-references, so it is always its
    own annotated node -- which is what makes the caller uniform -- and its
    subunits are reached by `HAS_SUBUNIT` when `with_subunits`.
    """
    node_element_ids = [node.element_id for node in nodes]
    if not with_subunits or not node_element_ids:
        return [(node, [node]) for node in nodes]
    query = f"""
        UNWIND $element_ids AS element_id
        MATCH (node) WHERE elementId(node) = element_id
        OPTIONAL MATCH (node)-[:HAS_SUBUNIT*1..{_SUBUNIT_MAX_DEPTH}]->(subunit)
        RETURN node AS node, collect(DISTINCT subunit) AS subunits
    """
    result = session.execute_query(query, params={"element_ids": node_element_ids})
    return [(row["node"], [row["node"]] + row["subunits"]) for row in result]


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
