"""Rebuild a BEL document from the COVID BEL KG stored in Neo4j.

The curated BEL files under `data/covid_kg/bel/` are a repository snapshot that
covers only ~95% of the stored `COVID_KG_BEL` collection (it was imported from
an intermediate revision, and from both curation branches at once). This script
goes the other way: it reads the collection back out of the database and writes
one BEL document that reproduces it exactly.

What is emitted:

* one statement per causal, correlative or membership edge -- the subject and
  object are the endpoints' `bel` property, verbatim, so no term is rebuilt
  from fields;
* the edge's provenance as `SET` lines (Citation, Support, and the
  `annotation*` properties), grouped by citation;
* the nodes carrying no such edge as bare term statements, so the isolated
  nodes a triple-shaped export cannot carry are kept;
* the `DEFINE NAMESPACE` / `DEFINE ANNOTATION` lines, which the graph does not
  store, lifted from the curated files.

The structural `HAS__*` edges (and the derived `HAS_MODIFIED_PROTEIN`,
`HAS_VARIANT_PROTEIN`, `HAS_LOCATED_*`, `HAS_FRAGMENTED_PROTEIN` ones) are
dropped: they are pybel's reification of what the `bel` strings already say,
and emitting them would invent statements the document never made.

A term is a sub-term when another term points **at** it through a double
underscore `HAS__*` edge -- a complex member, an activity's subject, a `pmod()`.
It is not one because it *has* sub-terms of its own, which is what a complex
always does. The derived single-underscore edges run the other way, from a
protein to its own proteoform (`p(APP) -> p(APP,frag("672_711"))`), and a
proteoform is a term in its own right, so they do not make one a sub-term
either.
"""

import collections
import pathlib
import re
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent / "main_analysis"))

import credentials
import neo4j

COLLECTION_NAME = "COVID_KG_BEL"
CURATED_BEL_DIR = pathlib.Path("data/covid_kg/bel")
OUTPUT_FILE_PATH = pathlib.Path("data/covid_kg/bel/covid_kg.bel")

DOCUMENT_VALUES = [
    ("Name", "COVID-19 Knowledge Graph"),
    (
        "Description",
        "The COVID-19 BEL knowledge graph as stored in the COMMUTE Neo4j "
        "database, serialised back to BEL. Rebuilt from the graph, not from "
        "the curation repository: the two differ, the graph being the version "
        "every analysis has been run against.",
    ),
    ("Version", "1.0.0"),
    ("Authors", "COVID-19 Disease Map / BMS curation teams"),
    ("ContactInfo", "adrien.rougny@aist.go.jp"),
]

# The edge types that are BEL statements. Everything else in the graph is
# structure (see the module docstring).
RELATION_TYPE_TO_RELATION = {
    "INCREASES": "increases",
    "DECREASES": "decreases",
    "DIRECTLY_INCREASES": "directlyIncreases",
    "DIRECTLY_DECREASES": "directlyDecreases",
    "ASSOCIATION": "association",
    "POSITIVE_CORRELATION": "positiveCorrelation",
    "NEGATIVE_CORRELATION": "negativeCorrelation",
    "REGULATES": "regulates",
    "CAUSES_NO_CHANGE": "causesNoChange",
    "IS_A": "isA",
    "TRANSLATED_TO": "translatedTo",
    "HAS_MEMBERS": "hasMembers",
    "HAS_MEMBER": "hasMember",
    "HAS_COMPONENTS": "hasComponents",
    "HAS_COMPONENT": "hasComponent",
    "ORTHOLOGOUS": "orthologous",
    "EQUIVALENT_TO": "equivalentTo",
    "BIOMARKER_FOR": "biomarkerFor",
    "PROGNOSTIC_BIOMARKER_FOR": "prognosticBiomarkerFor",
    "RATE_LIMITING_STEP_OF": "rateLimitingStepOf",
    "SUB_PROCESS_OF": "subProcessOf",
    "NO_CORRELATION": "noCorrelation",
    "HAS_ACTIVITY": "hasActivity",
    "PART_OF": "partOf",
    "TRANSCRIBED_TO": "transcribedTo",
    "ANALOGOUS": "analogous",
}

# `citationType`, `citationRef` and `evidence` become Citation and Support;
# every other `annotation*` property is a generic annotation whose name is the
# property name minus the prefix.
ANNOTATION_PROPERTY_PREFIX = "annotation"

_DEFINE_PATTERN = re.compile(
    r'^DEFINE\s+(NAMESPACE|ANNOTATION)\s+(\S+)\s+AS\s+(URL|PATTERN|LIST)\s+(.*)$'
)
# A namespace prefix only ever starts a function argument, so it is anchored on
# the `(` or `,` before it -- otherwise the colon inside a quoted identifier
# (`a(CHEBI:"NADH:ubiquinone reductase")`) reads as one too.
_NAMESPACE_PATTERN = re.compile(r'[(,]\s*([A-Za-z0-9_.]+)\s*:')
_NAMESPACE_IDENTIFIER_PATTERN = re.compile(
    r'[(,]\s*([A-Za-z0-9_.]+)\s*:\s*("[^"]*"|[^,)\s]+)'
)


def get_definitions_from_curated_files(directory):
    """Collect the `DEFINE` lines of the curated BEL files.

    A name may be defined several ways across the files (two Disease Ontology
    versions, `NCBITAXON` as both a URL and a pattern). The most frequent
    definition wins, ties going to the lexicographically greatest value, which
    for these dated URLs is the most recent one.
    """
    kind_and_name_to_values = collections.defaultdict(collections.Counter)
    for file_path in sorted(directory.glob("*.bel")):
        if file_path.resolve() == OUTPUT_FILE_PATH.resolve():
            continue
        with open(file_path, encoding="utf-8") as f:
            for line in f:
                match = _DEFINE_PATTERN.match(line.strip())
                if match is not None:
                    kind, name, as_type, value = match.groups()
                    kind_and_name_to_values[(kind, name)][
                        (as_type, value.strip())
                    ] += 1
    definitions = {}
    for (kind, name), values in kind_and_name_to_values.items():
        best = max(values.items(), key=lambda item: (item[1], item[0]))
        definitions[(kind, name)] = best[0]
    return definitions


def get_statements_and_isolated_nodes(session, collection_name):
    """Return the statement edges and the top-level terms carrying none of them.

    A term with no statement edge is written as a bare term line, so that the
    document keeps the isolated nodes a triple-shaped export cannot carry. Its
    own sub-terms are not written: reading the line rebuilds them.
    """
    query = """
        MATCH (:Collection {name: $collection_name})-[:HAS_ENTRY]->()
              -[:HAS_OBJ]->(:BELModel)-[:HAS_NODE]->(source)-[relation]->(target)
        WHERE type(relation) IN $relation_types
        RETURN source.bel AS source, type(relation) AS relation_type,
               target.bel AS target, properties(relation) AS properties
    """
    result = session.run(
        query,
        collection_name=collection_name,
        relation_types=list(RELATION_TYPE_TO_RELATION),
    )
    statements = [record.data() for record in result]
    query = """
        MATCH (:Collection {name: $collection_name})-[:HAS_ENTRY]->()
              -[:HAS_OBJ]->(:BELModel)-[:HAS_NODE]->(node)
        WHERE NOT EXISTS {
            MATCH (node)-[relation]-()
            WHERE type(relation) IN $relation_types
        }
        AND NOT EXISTS {
            MATCH ()-[relation]->(node)
            WHERE type(relation) STARTS WITH 'HAS__'
        }
        RETURN DISTINCT node.bel AS bel
    """
    result = session.run(
        query,
        collection_name=collection_name,
        relation_types=list(RELATION_TYPE_TO_RELATION),
    )
    isolated_bels = sorted(record["bel"] for record in result if record["bel"])
    return statements, isolated_bels


def get_namespaces(bel_strings):
    namespaces = set()
    for bel_string in bel_strings:
        namespaces.update(_NAMESPACE_PATTERN.findall(bel_string))
    return namespaces


def get_namespace_to_identifiers(bel_strings):
    namespace_to_identifiers = collections.defaultdict(set)
    for bel_string in bel_strings:
        for namespace, identifier in _NAMESPACE_IDENTIFIER_PATTERN.findall(bel_string):
            namespace_to_identifiers[namespace].add(identifier.strip('"'))
    return namespace_to_identifiers


def make_widened_definition(definition, used_values):
    """Widen an `AS LIST` definition to the values the document actually uses.

    A `LIST` definition lifted from one curated file enumerates that file's
    values, not the graph's. Left alone it would declare a vocabulary the
    document then violates, so the values in use are merged into it.
    """
    as_type, value = definition
    if as_type != "LIST" or not used_values:
        return definition
    listed_values = set(re.findall(r'"([^"]*)"', value))
    values = sorted(listed_values | set(used_values))
    return (as_type, "{" + ", ".join(make_quoted_string(v) for v in values) + "}")


def make_quoted_string(value):
    escaped = str(value).replace('"', '\\"')
    return f'"{escaped}"'


def make_set_string(name, values):
    if len(values) > 1:
        values_string = "{" + ", ".join(make_quoted_string(v) for v in values) + "}"
    else:
        values_string = make_quoted_string(values[0])
    return f"SET {name} = {values_string}"


def make_annotation_name(property_name, annotation_definition_names):
    name = property_name[len(ANNOTATION_PROPERTY_PREFIX):]
    # The graph stores `annotationMeshanatomy` for a `MeSHAnatomy` annotation:
    # the case is lost, so it is recovered from the curated definitions.
    for definition_name in annotation_definition_names:
        if definition_name.lower() == name.lower():
            return definition_name
    return name


def make_bel_document(statements, isolated_bels, definitions):
    annotation_definition_names = [
        name for (kind, name) in definitions if kind == "ANNOTATION"
    ]
    output_lines = []
    output_lines.append("#" * 70)
    output_lines.append("# Document properties")
    output_lines.append("#" * 70)
    for document_key, value in DOCUMENT_VALUES:
        output_lines.append(make_set_string(f"DOCUMENT {document_key}", [value]))
    used_namespaces = get_namespaces(
        [statement["source"] for statement in statements]
        + [statement["target"] for statement in statements]
        + isolated_bels
    )
    namespace_to_identifiers = get_namespace_to_identifiers(
        [statement["source"] for statement in statements]
        + [statement["target"] for statement in statements]
        + isolated_bels
    )
    annotation_name_to_values = collections.defaultdict(set)
    for statement in statements:
        for property_name, values in statement["properties"].items():
            if not property_name.startswith(ANNOTATION_PROPERTY_PREFIX):
                continue
            name = make_annotation_name(property_name, annotation_definition_names)
            if not isinstance(values, list):
                values = [values]
            annotation_name_to_values[name].update(str(v) for v in values)
    used_annotation_names = set(annotation_name_to_values)
    output_lines.append("")
    output_lines.append("#" * 70)
    output_lines.append("# Definitions")
    output_lines.append("#" * 70)
    undefined_namespaces = []
    for name in sorted(used_namespaces):
        definition = definitions.get(("NAMESPACE", name))
        if definition is None:
            # A namespace the curated files never define: kept, with a
            # permissive pattern, rather than dropped with its statements.
            definition = ("PATTERN", '".*"')
            undefined_namespaces.append(name)
        as_type, value = make_widened_definition(
            definition, namespace_to_identifiers.get(name, set())
        )
        output_lines.append(f"DEFINE NAMESPACE {name} AS {as_type} {value}")
    output_lines.append("")
    undefined_annotations = []
    for name in sorted(used_annotation_names):
        definition = definitions.get(("ANNOTATION", name))
        if definition is None:
            definition = ("PATTERN", '".*"')
            undefined_annotations.append(name)
        as_type, value = make_widened_definition(
            definition, annotation_name_to_values.get(name, set())
        )
        output_lines.append(f"DEFINE ANNOTATION {name} AS {as_type} {value}")
    output_lines.append("")
    output_lines.append("#" * 70)
    output_lines.append("# Statements")
    output_lines.append("#" * 70)
    # Grouped by citation, the way a curated document is written: the citation
    # is set once and the statements it supports follow.
    citation_to_statements = collections.defaultdict(list)
    for statement in statements:
        properties = statement["properties"]
        citation = (
            properties.get("citationType", "Other"),
            properties.get("citationRef", ""),
        )
        citation_to_statements[citation].append(statement)
    for citation in sorted(citation_to_statements):
        citation_type, citation_ref = citation
        output_lines.append("")
        output_lines.append(make_set_string("Citation", [citation_type, citation_ref]))
        # The evidence is part of the key: two statements can share endpoints
        # and relation, and without it their order follows the database's row
        # order, which differs from one run to the next.
        for statement in sorted(
            citation_to_statements[citation],
            key=lambda statement: (
                statement["source"],
                statement["relation_type"],
                statement["target"],
                str(statement["properties"].get("evidence", "")),
                sorted(str(item) for item in statement["properties"].items()),
            ),
        ):
            properties = statement["properties"]
            unset_names = []
            evidence = properties.get("evidence")
            if evidence:
                # A raw newline inside a quoted value is what makes a BEL
                # document unparsable, so the evidence is kept on one line.
                evidence = " ".join(str(evidence).split())
                output_lines.append(make_set_string("Support", [evidence]))
                unset_names.append("Support")
            for property_name in sorted(properties):
                if not property_name.startswith(ANNOTATION_PROPERTY_PREFIX):
                    continue
                values = properties[property_name]
                if not isinstance(values, list):
                    values = [values]
                if not values:
                    continue
                name = make_annotation_name(property_name, annotation_definition_names)
                output_lines.append(make_set_string(name, values))
                unset_names.append(name)
            relation = RELATION_TYPE_TO_RELATION[statement["relation_type"]]
            output_lines.append(
                f"{statement['source']} {relation} {statement['target']}"
            )
            for name in unset_names:
                output_lines.append(f"UNSET {name}")
        output_lines.append("UNSET Citation")
    if isolated_bels:
        output_lines.append("")
        output_lines.append("#" * 70)
        output_lines.append("# Terms carrying no statement")
        output_lines.append("#" * 70)
        output_lines.extend(isolated_bels)
    output_lines.append("")
    return (
        "\n".join(output_lines),
        undefined_namespaces,
        undefined_annotations,
    )


def main():
    definitions = get_definitions_from_curated_files(CURATED_BEL_DIR)
    print(
        f"{len(definitions)} definitions from {CURATED_BEL_DIR}: "
        f"{sum(1 for kind, _ in definitions if kind == 'NAMESPACE')} namespaces, "
        f"{sum(1 for kind, _ in definitions if kind == 'ANNOTATION')} annotations"
    )
    driver = neo4j.GraphDatabase.driver(
        f"bolt://{credentials.NEO4J_URI}:7687",
        auth=(credentials.NEO4J_USERNAME, credentials.NEO4J_PASSWORD),
    )
    with driver.session(database=credentials.NEO4J_DATABASE) as session:
        statements, isolated_bels = get_statements_and_isolated_nodes(
            session, COLLECTION_NAME
        )
    print(f"{len(statements)} statement edges, {len(isolated_bels)} isolated terms")
    bel_document, undefined_namespaces, undefined_annotations = make_bel_document(
        statements, isolated_bels, definitions
    )
    if undefined_namespaces:
        print(f"namespaces defined as a permissive pattern: {undefined_namespaces}")
    if undefined_annotations:
        print(f"annotations defined as a permissive pattern: {undefined_annotations}")
    OUTPUT_FILE_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_FILE_PATH, "w", encoding="utf-8") as f:
        f.write(bel_document)
    print(
        f"wrote {OUTPUT_FILE_PATH} "
        f"({len(bel_document.splitlines())} lines, "
        f"{OUTPUT_FILE_PATH.stat().st_size / 1e6:.1f} MB)"
    )


if __name__ == "__main__":
    main()
