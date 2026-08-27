"""The influence-graph projection of a BEL model.

A knowledge graph is read from its BEL files as a `momapy_bel.core.BELModel`.
This module keeps of that model what an influence graph is made of: the causal
relations and the entities they connect, plus the entities that are nothing
else's constituent. It takes a model and returns a model -- no session, no
collection names -- and imports nothing from the rest of the package.

Both `1_10`, which turns the AD knowledge graph into a CellDesigner map through
`commute_dm.bel2cd`, and `2_05`, which counts on the projection, call it.
"""

import momapy_bel.core
import momapy_bel.io.bel


# The BEL element classes the projection keeps: the molecular entities and the
# phenotypes. `PopulationAbundance` is deliberately absent, and so are the
# process terms (reactions, translocations, degradations), which are neither an
# entity nor a relation. Matched on the **exact** class, since every abundance
# derives from `Abundance`.
PROJECTED_ELEMENT_CLASSES = frozenset(
    {
        momapy_bel.core.Abundance,
        momapy_bel.core.ProteinAbundance,
        momapy_bel.core.GeneAbundance,
        momapy_bel.core.RNAAbundance,
        momapy_bel.core.MicroRNAAbundance,
        momapy_bel.core.ComplexAbundance,
        momapy_bel.core.CompositeAbundance,
        momapy_bel.core.Activity,
        momapy_bel.core.BiologicalProcess,
        momapy_bel.core.Pathology,
    }
)

# The causal relations the projection keeps. The correlative, hierarchical,
# equivalence and no-effect ones are not causal and are left out.
INFLUENCE_RELATION_CLASSES = (
    momapy_bel.core.Increases,
    momapy_bel.core.DirectlyIncreases,
    momapy_bel.core.TranslatedTo,
    momapy_bel.core.Decreases,
    momapy_bel.core.DirectlyDecreases,
    momapy_bel.core.Regulates,
)


def get_element_bel_string(element):
    """The BEL string of one momapy_bel element -- its label wherever one is due.

    Canonical: a complex's members are sorted, so the string is the same from one
    run to the next, which is what makes it a usable sort key.
    """
    return momapy_bel.io.bel.BELWriter._bel_element_to_string(element)


def get_bel_influence_graph_projection(bel_model):
    """The `BELModel` of the influence relations of `bel_model`.

    Its `statements` are the causal relations between projected elements, plus
    every projected element carrying none of them, as a bare abundance: an
    isolated node is part of the projection and `2_05` counts it.

    The projected elements are the molecular entities and the phenotypes of the
    model minus the **structural constituents** -- an element another entity
    holds as a sub-term (a complex member, an activity's subject) and that
    carries no causal relation of its own, whose wiring lives on its container.
    Everything else is kept, including an element that only appears in a
    reaction or in a correlative statement.

    Elements are **value-equal**, so two BEL terms that describe the same thing
    are one element.
    """
    elements = {}
    for element in bel_model.descendants():
        if type(element) in PROJECTED_ELEMENT_CLASSES:
            elements.setdefault(element, element)
    constituents = {
        sub_element
        for element in elements
        for sub_element in _get_sub_elements(element)
        if type(sub_element) in PROJECTED_ELEMENT_CLASSES
    }
    relations = [
        statement
        for statement in bel_model.statements
        if isinstance(statement, INFLUENCE_RELATION_CLASSES)
        and statement.source in elements
        and statement.target in elements
    ]
    causally_wired = {
        endpoint
        for relation in relations
        for endpoint in (relation.source, relation.target)
    }
    projected_elements = {
        element
        for element in elements
        if element not in constituents or element in causally_wired
    }
    statements = []
    wired = set()
    for relation in relations:
        source = elements[relation.source]
        target = elements[relation.target]
        if source not in projected_elements or target not in projected_elements:
            continue
        statements.append(type(relation)(source=source, target=target))
        wired.update((source, target))
    statements.extend(
        element for element in projected_elements if element not in wired
    )
    return momapy_bel.core.BELModel(statements=frozenset(statements))


def _get_sub_elements(element):
    if isinstance(element, momapy_bel.core.Activity):
        return [element.abundance]
    if isinstance(
        element,
        (momapy_bel.core.ComplexAbundance, momapy_bel.core.CompositeAbundance),
    ):
        return list(element.members)
    return []
