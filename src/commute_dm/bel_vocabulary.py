"""The BEL vocabulary: the node classes, relation types and structural edge types
of a BEL knowledge graph as it is stored in Neo4j.

Constants only -- no queries, no imports from the rest of the package. Everything
that reads a BEL graph reads its vocabulary from here: the influence-graph
projection and the element loader (`commute_dm.bel_export`) and the statistics
notebook.
"""

# --- Node classes -----------------------------------------------------------

MOLECULAR_ENTITY_NODE_TYPES = {
    "Protein",
    "Complex",
    "Abundance",
    "Composite",
    "Gene",
    "Rna",
    "MicroRna",
    "Activity",
}
PHENOTYPE_NODE_TYPES = {"BiologicalProcess", "Pathology"}
POPULATION_NODE_TYPES = {"Population"}
ENTITY_NODE_TYPES = (
    MOLECULAR_ENTITY_NODE_TYPES | PHENOTYPE_NODE_TYPES | POPULATION_NODE_TYPES
)
REACTION_NODE_TYPES = {
    "Reaction",
    "Translocation",
    "Degradation",
    "CellSecretion",
    "CellSurfaceExpression",
}
LOCATION_NODE_TYPES = {"Location", "ToLocation", "FromLocation"}

ALL_NODE_TYPES = ENTITY_NODE_TYPES | REACTION_NODE_TYPES | LOCATION_NODE_TYPES

# --- Relation types ---------------------------------------------------------

INFLUENCE_RELATIONSHIP_TYPES = {
    "INCREASES",
    "DECREASES",
    "DIRECTLY_INCREASES",
    "DIRECTLY_DECREASES",
    "REGULATES",
    "TRANSLATED_TO",
}
CORRELATIVE_RELATIONSHIP_TYPES = {
    "ASSOCIATION",
    "POSITIVE_CORRELATION",
    "NEGATIVE_CORRELATION",
}
HIERARCHICAL_RELATIONSHIP_TYPES = {"IS_A"}
EQUIVALENCE_RELATIONSHIP_TYPES = {"ORTHOLOGOUS", "EQUIVALENT_TO"}
NO_EFFECT_RELATIONSHIP_TYPES = {"CAUSES_NO_CHANGE"}

INFLUENCE_OR_CORRELATIVE_RELATIONSHIP_TYPES = (
    INFLUENCE_RELATIONSHIP_TYPES | CORRELATIVE_RELATIONSHIP_TYPES
)

ALL_RELATIONSHIP_TYPES = (
    INFLUENCE_RELATIONSHIP_TYPES
    | CORRELATIVE_RELATIONSHIP_TYPES
    | HIERARCHICAL_RELATIONSHIP_TYPES
    | EQUIVALENCE_RELATIONSHIP_TYPES
    | NO_EFFECT_RELATIONSHIP_TYPES
)

# --- Structural (`HAS__*`) edge types ---------------------------------------

# A complex's or composite's members, and an activity's subject: the sub-terms
# that are entities in their own right.
MEMBER_EDGE_TYPES = frozenset(
    {
        "HAS__PROTEIN",
        "HAS__ABUNDANCE",
        "HAS__COMPLEX",
        "HAS__COMPOSITE",
        "HAS__GENE",
        "HAS__RNA",
        "HAS__MICRO_RNA",
    }
)
# The sub-terms that qualify an entity rather than being one.
MODIFIER_EDGE_TYPES = frozenset(
    {
        "HAS__PMOD",
        "HAS__VARIANT",
        "HAS__FRAGMENT",
        "HAS__LOCATION",
        "HAS__GMOD",
    }
)
# `HAS__*` edges of BEL *process* terms (translocations, reactions, degradations,
# lists). None of these terms is projected, so none is ever loaded -- but they are
# enumerated so that `commute_dm.bel_export` can raise on an edge type it has never
# seen instead of silently building an element with a missing part.
OTHER_STRUCTURAL_EDGE_TYPES = frozenset(
    {
        "HAS__TO_LOCATION",
        "HAS__FROM_LOCATION",
        "HAS__REACTANTS",
        "HAS__PRODUCTS",
    }
)
STRUCTURAL_EDGE_TYPES = (
    MEMBER_EDGE_TYPES | MODIFIER_EDGE_TYPES | OTHER_STRUCTURAL_EDGE_TYPES
)
