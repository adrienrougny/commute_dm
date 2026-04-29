import os
import pathlib
import shutil
import lxml.etree
import functools
import operator

import neo4j


def flatten_list(input_list):
    def _flatten_rec(a, b):
        if isinstance(b, list):
            b = flatten_list(b)
        else:
            b = [b]
        return operator.iconcat(a, b)

    return functools.reduce(_flatten_rec, input_list, [])


def remake_dir(path):
    if pathlib.Path(path).exists():
        shutil.rmtree(path)
    os.makedirs(path, exist_ok=True)


def list_dir(path):
    files = []
    for file_name in os.listdir(path):
        if not file_name.startswith("."):
            files.append((file_name, os.path.join(path, file_name)))
    return files


def rename_file(file_name, new_extension):
    pre, ext = os.path.splitext(file_name)
    return f"{pre}.{new_extension}"


def add_annotation_to_file(
    annotation, entity_id, input_file_path, output_file_path
):

    _NSMAP = {
        "rdf": "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
        "bqbiol": "http://biomodels.net/biology-qualifiers/",
        "sbml": "http://www.sbml.org/sbml/level2/version4",
    }

    def _make_lxml_element(
        tag, namespace=None, attributes=None, text=None, nsmap=None
    ):
        if namespace is not None:
            lxml_tag = f"{{{namespace}}}{tag}"
        else:
            lxml_tag = tag
        if nsmap is None:
            nsmap = {}
        if attributes is None:
            attributes = {}
        lxml_element = lxml.etree.Element(lxml_tag, nsmap=nsmap, **attributes)
        if text is not None:
            lxml_element.text = text
        return lxml_element

    def _get_or_make_and_add_child_element(
        parent_element,
        tag,
        namespace=None,
        attributes=None,
        text=None,
        nsmap=None,
    ):
        child_element = None
        if namespace is not None:
            prefix = None
            for key, value in nsmap.items():
                if value == namespace:
                    prefix = key
                    break
            if prefix is None:
                raise ValueError(
                    f"Namespace '{namespace}' should be in the nsmap"
                )
            query = f"./{prefix}:{tag}"
        else:
            query = f"./{tag}"
        results = parent_element.xpath(query, namespaces=nsmap)
        for child_element in results:
            break
        if child_element is None:
            child_element = _make_lxml_element(
                tag=tag,
                namespace=namespace,
                attributes=attributes,
                text=text,
                nsmap=nsmap,
            )
            parent_element.append(child_element)
        return child_element

    cd_document = lxml.etree.parse(input_file_path)
    query = f"//*[@id='{entity_id}']"
    results = cd_document.xpath(query)
    cd_entity = None
    for cd_entity in results:
        break
    if cd_entity is None:
        print(
            f"Entity with ID '{entity_id}' not found in file '{input_file_path}'"
        )
    else:
        cd_annotation = _get_or_make_and_add_child_element(
            parent_element=cd_entity,
            tag="annotation",
            namespace=_NSMAP["sbml"],
            nsmap=_NSMAP,
        )
        cd_rdf = _get_or_make_and_add_child_element(
            parent_element=cd_annotation,
            tag="RDF",
            namespace=_NSMAP["rdf"],
            nsmap=_NSMAP,
        )
        cd_description = _get_or_make_and_add_child_element(
            parent_element=cd_rdf,
            tag="Description",
            namespace=_NSMAP["rdf"],
            nsmap=_NSMAP,
        )
        cd_description.attrib[f"{{{_NSMAP['rdf']}}}about"] = f"#{entity_id}"
        cd_is = _make_lxml_element("is", _NSMAP["bqbiol"], nsmap=_NSMAP)
        cd_bag = _make_lxml_element("Bag", _NSMAP["rdf"], nsmap=_NSMAP)
        for resource in annotation:
            cd_li = _make_lxml_element("li", _NSMAP["rdf"], nsmap=_NSMAP)
            cd_li.attrib[f"{{{_NSMAP['rdf']}}}resource"] = resource
            cd_bag.append(cd_li)
        cd_is.append(cd_bag)
        cd_description.append(cd_is)
    cd_document.write(
        output_file_path,
        pretty_print=True,
        xml_declaration=True,
        encoding="UTF-8",
    )
    return cd_document


def replace_label(old_label, new_label, input_file_path, output_file_path):
    with open(input_file_path) as f:
        s = f.read()
    s = s.replace(old_label, new_label)
    with open(output_file_path, "w") as f:
        f.write(s)


def subgraph_to_cypherl(
    nodes, relationships, output_file_path
):  # query must return nodes, relationships

    def _get_node_labels(node):
        return node.labels

    def _get_node_properties(node):
        properties = dict(node.items())
        properties["id"] = node.id
        return properties

    def _get_relationship_type(relationship):
        return relationship.type

    def _get_relationship_start_node(relationship):
        return relationship.start_node

    def _get_relationship_end_node(relationship):
        return relationship.end_node

    def _make_property_value_string(property_value):
        if isinstance(property_value, str):
            property_value = property_value.replace("\n", " ")
            return f'"{property_value}"'
        elif isinstance(property_value, list):
            return f"[{','.join([_make_property_value_string(list_value) for list_value in property_value])}]"
        else:
            return property_value

    def _make_properties_string(node):
        property_string = f" {{{','.join([f'{property_name}: {_make_property_value_string(property_value)}' for property_name, property_value in _get_node_properties(node).items()])}}}"
        return property_string

    def _make_node_string(
        node, variable=None, with_labels=True, with_properties=True
    ):
        if variable is not None:
            variable_string = variable
        else:
            variable_string = ""
        if with_labels:
            labels_string = f':{":".join(_get_node_labels(node))}'
        else:
            labels_string = ""
        if with_properties:
            properties_string = _make_properties_string(node)
        else:
            properties_string = ""
        return f"({variable_string}{labels_string}{properties_string})"

    def _make_create_node_string(node):
        node_string = _make_node_string(node)
        create_node_string = f"CREATE {node_string};"
        return create_node_string

    def _make_relationship_string(relationship):
        relationship_string = f"{_make_node_string(_get_relationship_start_node(relationship))}-[:{_get_relationship_type(relationship)}]->{_make_node_string(_get_relationship_end_node(relationship))}"
        return relationship_string

    def _make_merge_relationship_string(relationship):
        relationship_string = _make_relationship_string(relationship)
        create_relationship_string = f"MERGE {relationship_string};"
        return create_relationship_string

    def _make_match_create_relationship_string(relationship):
        start_node = _get_relationship_start_node(relationship)
        end_node = _get_relationship_end_node(relationship)
        start_node_variable = "u"
        end_node_variable = "v"
        match_create_relationship_string = (
            f"MATCH {_make_node_string(start_node, variable=start_node_variable, with_labels=True, with_properties=False)},"
            f"{_make_node_string(end_node, variable=end_node_variable, with_labels=True, with_properties=False)} "
            f"WHERE {start_node_variable}.id = {start_node.id} AND {end_node_variable}.id = {end_node.id} "
            f"CREATE ({start_node_variable})-[:{_get_relationship_type(relationship)}]->({end_node_variable});"
        )
        return match_create_relationship_string

    cypherl_commands = []
    for node in nodes:
        cypherl_command = _make_create_node_string(node)
        cypherl_commands.append(cypherl_command)
    for relationship in relationships:
        cypherl_command = _make_match_create_relationship_string(relationship)
        cypherl_commands.append(cypherl_command)
    with open(output_file_path, "w") as f:
        f.write("\n".join(cypherl_commands))


def merge_relationship(session, start_node, end_node, relationship_type):
    query = f"""
        MATCH (start_node), (end_node)
        WHERE elementId(start_node) = '{start_node.element_id}'
        AND elementId(end_node) = '{end_node.element_id}'
        MERGE (start_node)-[r:{relationship_type}]->(end_node)
    """
    session.execute_query(query)
